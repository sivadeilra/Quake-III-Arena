#![allow(dead_code)]
#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_macros)]
#![allow(unused_parens)]

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
pub mod perf;
pub mod q_shared;
pub mod qcommon;
pub mod vec3;

pub mod prelude {
    pub use crate::bounds::*;
    pub use crate::num_utils::*;
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
}
