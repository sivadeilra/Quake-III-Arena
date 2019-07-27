//!
//! The rendezvous protocol allows two tracers to rendezvous and compare
//! two values, one from each tracer. The protocol has two participants, PRIMARY and SECONDARY.
//! PRIMARY and SECONDARY have different behaviors.

use lazy_static::lazy_static;
use log::{debug, error, info, trace};
use std::cell::RefCell;
use std::collections::VecDeque;
use std::fmt::Write;
use std::sync::mpsc::{self, Receiver, SyncSender};
use std::sync::{Arc, Condvar, Mutex, MutexGuard};
use std::thread::{JoinHandle, Thread};

enum MsgToPrimary {
    Value(TraceValue),
}

enum MsgToSecondary {
    Ack,
}

struct Session {
    max_history: usize,
    history: VecDeque<HistoryEntry>,
    total_values: u64,
}

struct PrimaryRole {
    rx: Receiver<MsgToPrimary>,
    tx: SyncSender<MsgToSecondary>,
    session: Mutex<Session>,
}

struct SecondaryRole {
    rx: Receiver<MsgToSecondary>,
    tx: SyncSender<MsgToPrimary>,
}

enum Role {
    Primary(PrimaryRole),
    Secondary(SecondaryRole),
}

pub struct Tracer {
    name: &'static str,
    role: Role,
}

#[derive(Clone, PartialEq)]
pub enum TraceValue {
    I32(i32),
    U32(u32),
    String(String),
    Str(&'static str),
    F32(f32),
}

const TRACE_VALUE_F32_EPSILON: f32 = 0.00001;

impl TraceValue {
    fn nearly_eq(&self, other: &Self) -> bool {
        if self == other {
            return true;
        }
        match (self, other) {
            (TraceValue::F32(a), TraceValue::F32(b)) => (a - b).abs() < TRACE_VALUE_F32_EPSILON,
            _ => false,
        }
    }
}

use core::fmt::Debug;
impl Debug for TraceValue {
    fn fmt(&self, fmt: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            TraceValue::I32(n) => write!(fmt, "{:?}", n),
            TraceValue::U32(n) => write!(fmt, "{:?}", n),
            TraceValue::Str(s) => write!(fmt, "{:?}", s),
            TraceValue::F32(f) => write!(fmt, "{:?}", f),
            TraceValue::String(s) => write!(fmt, "{:?}", s),
        }
    }
}

#[derive(Clone, Debug)]
struct HistoryEntry {
    value: TraceValue,
    // call stack?
}

impl Tracer {
    pub fn trace_value(&self, value: TraceValue) {
        match &self.role {
            Role::Primary(ref primary) => {
                trace!("primary: waiting for value from secondary");
                let msg = primary.rx.recv().unwrap();
                match msg {
                    MsgToPrimary::Value(secondary_value) => {
                        trace!("primary: received value from secondary: {:?}", value);
                        if !value.nearly_eq(&secondary_value) {
                            let session = primary.session.lock().unwrap();
                            let mut history_text = String::new();
                            for e in session.history.iter() {
                                write!(history_text, "{:?}\n", e.value).unwrap();
                            }
                            error!(
                                "DIVERGENCE DETECTED!  History:\n{}\nPRIMARY: {:?}\nSECONDARY: {:?}\n",
                                history_text, value, secondary_value
                            );
                            DebugBreak();
                            panic!(
                                "DIVERGENCE DETECTED!\n{} : {:?}\n{}: {:?}\nHistory:\n{}",
                                "PRIMARY", value, "SECONDARY", secondary_value, history_text
                            );
                        }
                        // message looks ok, release the secondary
                        trace!("primary: looks good, sending Ack to secondary");
                        {
                            let mut session = primary.session.lock().unwrap();
                            let seq_num = session.total_values;
                            session.total_values += 1;
                            trace!("primary: # {:4} good: {:?}", seq_num, value);
                            if session.history.len() == session.max_history {
                                session.history.pop_front();
                            }
                            session.history.push_back(HistoryEntry { value });
                        }
                        primary.tx.send(MsgToSecondary::Ack).unwrap();
                    }
                }
            }
            Role::Secondary(ref secondary) => {
                trace!("secondary: sending value to primary: {:?}", value);
                secondary.tx.send(MsgToPrimary::Value(value)).unwrap();
                trace!("secnodary: waiting for Ack");
                let _response = secondary.rx.recv().unwrap();
                trace!("secondary: got Ack");
            }
        }
    }

    pub fn trace_string(&self, s: String) {
        trace_value(TraceValue::String(s))
    }

    pub fn trace_str(&self, s: &'static str) {
        trace_value(TraceValue::Str(s));
    }

    pub fn trace_i32(&self, i: i32) {
        trace_value(TraceValue::I32(i))
    }
}

std::thread_local! {
    static AMBIENT_TRACER: RefCell<Option<Arc<Tracer>>> = RefCell::new(None);
}

/// This runs two closures in parallel, one in the current thread and one in
/// a worker thread. Because this code always calls join on the second thread,
/// we can safely use spawn_unchecked(). This allows the closures to contain
/// references to items on the stack.
fn run_two_threads<A, B, OutputA, OutputB>(a: A, b: B) -> (OutputA, OutputB)
where
    A: FnOnce() -> OutputA + Send,
    B: FnOnce() -> OutputB + Send,
    OutputA: Send,
    OutputB: Send,
{
    let builder = std::thread::Builder::new();
    let joiner = unsafe {
        builder.spawn_unchecked(move || {
            trace!("secondary thread is starting");
            let output = b();
            trace!("secondary thread is done");
            output
        })
    }
    .unwrap();
    let output_a = a();
    trace!("waiting for secondary thread to finish...");
    let output_b = joiner.join().unwrap();
    (output_a, output_b)
}

/// Returns (ref_output, test_output)
pub fn parallel_trace<RefCode, TestCode, Output>(
    ref_impl: RefCode,
    test_impl: TestCode,
) -> (Output, Output)
where
    RefCode: FnOnce(Arc<Tracer>) -> Output + Send,
    TestCode: FnOnce(Arc<Tracer>) -> Output + Send,
    Output: Send,
{
    fn call_impl<F: FnOnce(Arc<Tracer>) -> Output + Send, Output: Send>(
        code_impl: F,
        tracer: Arc<Tracer>,
    ) -> Output {
        let name = tracer.name;
        AMBIENT_TRACER.with(|t| *t.borrow_mut() = Some(tracer.clone()));
        trace!("{}: calling", name);
        let output = code_impl(tracer);
        AMBIENT_TRACER.with(|t| *t.borrow_mut() = None);
        trace!("{}: done", name);
        output
    }

    let secondary_name = "REF";
    let primary_name = "TEST";

    let (primary_tx, secondary_rx) = mpsc::sync_channel(1);
    let (secondary_tx, primary_rx) = mpsc::sync_channel(1);

    let primary_tracer = Tracer {
        role: Role::Primary(PrimaryRole {
            tx: primary_tx,
            rx: primary_rx,
            session: Mutex::new(Session {
                max_history: 100,
                history: VecDeque::new(),
                total_values: 0,
            }),
        }),
        name: primary_name,
    };

    let secondary_tracer = Tracer {
        role: Role::Secondary(SecondaryRole {
            tx: secondary_tx,
            rx: secondary_rx,
        }),
        name: secondary_name,
    };

    // secondary runs on a separate thread
    // primary runs on this thread

    let outputs = run_two_threads(
        move || call_impl(test_impl, Arc::new(primary_tracer)),
        move || call_impl(ref_impl, Arc::new(secondary_tracer)),
    );
    trace!("tracing is complete.");
    outputs
}

pub fn trace_value(value: TraceValue) {
    AMBIENT_TRACER.with(|t| {
        if let Some(tracer) = t.borrow().as_ref() {
            tracer.trace_value(value);
        } else {
        // ignore the call, we're not in a trace context
        }
    });
}

pub fn trace_str(s: &'static str) {
    trace_value(TraceValue::Str(s));
}

pub fn trace_string(s: String) {
    trace_value(TraceValue::String(s));
}

mod kernel32 {
    #[link(name = "kernel32")]
    extern "stdcall" {
        pub fn DebugBreak();
    }
}
fn DebugBreak() {
    unsafe {
        kernel32::DebugBreak();
    }
}

#[no_mangle]
pub extern "C" fn trace_str_raw(s: *const u8, len: usize) {
    unsafe {
        let bytes = core::slice::from_raw_parts::<'static, u8>(s, len);
        trace_str(core::str::from_utf8_unchecked(bytes))
    }
}

#[no_mangle]
pub extern "C" fn trace_string_copy_raw(s: *const u8, len: usize) {
    unsafe {
        let bytes = core::slice::from_raw_parts::<'_, u8>(s, len);
        trace_string(core::str::from_utf8_unchecked(bytes).to_string())
    }
}

#[no_mangle]
pub extern "C" fn trace_i32(i: i32) {
    trace_value(TraceValue::I32(i));
}

#[no_mangle]
pub extern "C" fn trace_f32(f: f32) {
    trace_value(TraceValue::F32(f))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn do_right_thing(t: Arc<Tracer>) {
        t.trace_i32(42);
        t.trace_str("...");
        t.trace_str("blah blah blah");
    }

    fn do_right_thing_ambient(t: Arc<Tracer>) {
        trace_i32(42);
        trace_str("...");
        trace_str("blah blah blah");
    }

    #[test]
    fn test_tracer() {
        env_logger::init();
        parallel_trace(do_right_thing, do_right_thing_ambient);
    }

}

// find a better home for this
pub fn trace_vec3(v: crate::vec3::vec3_t) {
    trace_f32(v[0]);
    trace_f32(v[1]);
    trace_f32(v[2]);
}
