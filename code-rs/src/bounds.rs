use crate::prelude::*;

#[derive(Copy, Clone, Debug, Default)]
pub struct vec3_bounds {
    pub mins: vec3_t,
    pub maxs: vec3_t,
}

pub fn ClearBounds() -> vec3_bounds {
    vec3_bounds {
        mins: vec3_t::from_scalar(99999.0),
        maxs: vec3_t::from_scalar(-99999.0),
    }
}

pub fn AddPointToBounds(v: vec3_t, bounds: &mut vec3_bounds) {
    for i in 0..3 {
        let val = v[i];
        if val < bounds.mins[i] {
            bounds.mins[i] = val;
        }
        if val > bounds.maxs[i] {
            bounds.maxs[i] = val;
        }
    }
}
