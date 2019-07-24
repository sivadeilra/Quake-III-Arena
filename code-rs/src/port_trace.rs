
use std::cell::RefCell;
use std::collections::VecDeque;
use std::sync::{Arc, Condvar, Mutex, MutexGuard};
use std::thread::{JoinHandle, Thread};
use lazy_static::lazy_static;
use log::{info, debug, error};

struct Session {
    next_value: Option<TraceValue>,
    max_history: usize,
    history: VecDeque<HistoryEntry>,
    total_values: u64
}

#[derive(Clone)]
pub struct Tracer {
    session_ptr: Arc<(Mutex<Session>, Condvar)>,
    is_ref: bool,
    name: &'static str
}

#[derive(Clone, PartialEq, Eq)]
pub enum TraceValue {
    I32(i32),
    U32(u32),
    String(String),
    Str(&'static str),
}

use core::fmt::Write;
impl core::fmt::Debug for TraceValue {
    fn fmt(&self, fmt: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            TraceValue::I32(n) => write!(fmt, "{:?}", n),
            TraceValue::U32(n) => write!(fmt, "{:?}", n),
            TraceValue::Str(s) => write!(fmt, "{:?}", s),
            TraceValue::String(s) => write!(fmt, "{:?}", s),
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct HistoryEntry {
    value: TraceValue,
    // call stack?
}

/*
lazy_static::lazy_static! {

// static ref TRACE_STATE: Mutex<Option<TraceSession>> = Mutex::new(None);
static ref TRACE_GLOBALS: TraceGlobals = TraceGlobals {
    data: Mutex::new(TraceData {
        history: Vec::new(),
        next_value: None,
        test_thread_joiner: None
    }),
    condvar: Condvar::new()
};
}

*/


impl Tracer {
    pub fn trace_value(&self, value: TraceValue) {
        info!("{}: trace_value: {:?}", self.name, value);
        use std::fmt::Write;
        let mut g = self.session_ptr.0.lock().unwrap();
        let session = &mut *g;
        if session.next_value.is_some() {
            // match them up
            let other_value = session.next_value.take().unwrap();
            if value != other_value {
                let mut history_text = String::new();
                for e in session.history.iter() {
                    write!(history_text, "{:?}\n", e.value).unwrap();
                }
                error!("DIVERGENCE DETECTED!\nREF : {:?}\nTEST: {:?}\nHistory:\n{}",
                    if self.is_ref { &value } else { &other_value },
                    if self.is_ref { &other_value } else { &value },
                    history_text);
                DebugBreak();
                panic!("DIVERGENCE DETECTED!\nREF : {:?}\nTEST: {:?}\nHistory:\n{}",
                    if self.is_ref { &value } else { &other_value },
                    if self.is_ref { &other_value } else { &value },
                    history_text);
            }
            info!("{}: found existing matching value", self.name);
            drop(other_value);
            if session.history.len() == session.max_history {
                drop(session.history.pop_front());
            }
            session.history.push_back(HistoryEntry {
                value
            });
            session.total_values += 1;
            self.session_ptr.1.notify_one();
        } else {
            // waiting for peer
            debug!("{}: adding value", self.name);
            session.next_value = Some(value);
            // we contributed a value. unblock the peer, 
            while g.next_value.is_some() {
                debug!("{}: waiting for other impl to observe value", self.name);
                g = self.session_ptr.1.wait(g).unwrap();
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

std::thread_local!{
    static AMBIENT_TRACER: RefCell<Option<Tracer>> = RefCell::new(None);
}

pub fn parallel_trace<R, T>(ref_impl: R, test_impl: T) 
    where R: FnOnce(Tracer) + 'static + Send,
    T: FnOnce(Tracer) + 'static + Send
{

    fn call_impl<F: FnOnce(Tracer) + 'static + Send>(code_impl: F, tracer: Tracer) {
        let name = tracer.name;
        AMBIENT_TRACER.with(|t| *t.borrow_mut() = Some(tracer.clone()));
        info!("{}: calling", name);
        code_impl(tracer);
        AMBIENT_TRACER.with(|t| *t.borrow_mut() = None);
        info!("{}: done", name);
    }

    // ref runs on a separate thread
    // test runs on this thread

    let session_arc = Arc::new(
        (Mutex::new(Session {
            max_history: 50,
            history: VecDeque::new(),
            next_value: None,
            total_values: 0       
        }), Condvar::new()));

    let ref_session_arc = session_arc.clone();
    let test_session_arc = ref_session_arc.clone();

    let ref_thread = std::thread::spawn(move || {
        info!("ref thread is starting");
        let ref_tracer = Tracer { session_ptr: ref_session_arc, is_ref: true, name: "REF_",  };
        call_impl(ref_impl, ref_tracer);
    });

    // Now run the test code
    {
        let test_tracer = Tracer { session_ptr: test_session_arc, is_ref: false, name: "TEST" };
        call_impl(test_impl, test_tracer);
    }

    info!("waiting for ref thread to finish...");
    ref_thread.join().unwrap();

    {
        let g = session_arc.0.lock().unwrap();
        info!("tracing is complete.  total entries: {}", g.total_values);
    }
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
    unsafe { kernel32::DebugBreak(); }
}
pub extern "C" fn trace_string_raw(s: *const u8, len: usize) {
    unsafe {
        let bytes = core::slice::from_raw_parts::<u8>(s, len);
        trace_string(core::str::from_utf8_unchecked(bytes).to_string())
    }
}


#[no_mangle]
pub extern "C" fn trace_i32(i: i32) {
    trace_value(TraceValue::I32(i));
}

#[cfg(test)]
mod tests {
use super::*;

fn do_right_thing(t: Tracer) {
    t.trace_i32(42);
    t.trace_str("...");
    t.trace_str("blah blah blah");
}

fn do_right_thing_ambient(t: Tracer) {
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

