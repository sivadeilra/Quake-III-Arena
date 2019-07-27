#![allow(dead_code)]
#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_macros)]
#![allow(non_upper_case_globals)]
#![allow(unused_imports)]
#![feature(thread_spawn_unchecked)]

macro_rules! expand_axis {
    ($axis:ident) => {};
    ($axis:ident #axis $($rest:tt)*) => { expand_axis!($axis $($rest)*) };
    ($axis:ident $next:tt $($rest:tt)*) => { $next expand_axis!($axis $($rest)*) }
}

// replaces #axis with x, y, z, three times
macro_rules! for_each_axis {
    ($($t:tt)*) => {
        expand_axis!(x $($t)*);
        expand_axis!(y $($t)*);
        expand_axis!(z $($t)*)
    }
}

macro_rules! todo_type {
    ($name:ident) => {
        pub type $name = ();
    };
}

pub mod bounds;
pub mod cvar;
pub mod dbg_logger;
pub mod game;
pub mod math;
pub mod perf;
pub mod port_trace;
pub mod q_math;
pub mod q_shared;
pub mod qcommon;
pub mod qfiles;
pub mod vec3;
pub mod zerocopy;

pub mod prelude {
    pub use crate::bounds::*;
    pub use crate::cvar::*;
    pub use crate::math::*;
    pub use crate::num_utils::*;
    pub use crate::q_math::*;
    pub use crate::q_shared::*;
    pub use crate::side::*;
    pub use crate::vec3::*;
}

pub fn is_close_to(a: f32, b: f32, epsilon: f32) -> bool {
    let d = a - b;
    d.abs() <= epsilon
}

mod side {
    #[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Default, Hash)]
    pub struct Side(pub u8);
    impl Side {
        pub const FRONT: Side = Side(0);
        pub const BACK: Side = Side(1);
        pub const ON: Side = Side(2);
        pub const CROSS: Side = Side(3);
        pub fn index(self) -> usize {
            self.0 as usize
        }
    }
    pub const SIDE_FRONT: Side = Side::FRONT;
    pub const SIDE_BACK: Side = Side::BACK;
    pub const SIDE_ON: Side = Side::ON;
    pub const SIDE_CROSS: Side = Side::CROSS;
}

pub mod num_utils {
    pub fn is_odd<N>(n: N) -> bool
    where
        N: std::ops::Rem<N, Output = N>,
        N: From<u8>,
        N: PartialEq,
    {
        n % N::from(2u8) != N::from(0u8)
    }

    pub fn fmin(a: f32, b: f32) -> f32 {
        if b < a {
            b
        } else {
            a
        }
    }
    pub fn fmax(a: f32, b: f32) -> f32 {
        if b > a {
            b
        } else {
            a
        }
    }
}

fn range_len<N>(start: N, count: N) -> std::ops::Range<N>
where
    N: Copy + std::ops::Add<N, Output = N>,
{
    start..start + count
}

fn fill<T: Copy>(items: &mut [T], value: T) {
    for ii in items.iter_mut() {
        *ii = value;
    }
}

pub mod encoding {
    pub fn get_le_u32(bytes: &[u8]) -> u32 {
        (bytes[0] as u32)
            | (bytes[1] as u32) << 8
            | (bytes[2] as u32) << 16
            | (bytes[3] as u32) << 24
    }

    pub fn get_le_i32(bytes: &[u8]) -> i32 {
        get_le_u32(bytes) as i32
    }
}

#[no_mangle]
extern "C" fn rust_init_runtime() {
    log::set_logger(&crate::dbg_logger::DBG_LOGGER).unwrap();
    log::set_max_level(log::LevelFilter::Debug);
    log::debug!("this is debug");
    log::warn!("this is warn");
    log::info!("this is info");
    log::error!("this is error");
}
