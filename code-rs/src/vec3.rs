pub type vec_t = f32;

// see game/q_shared.h

#[derive(Copy, Clone, Debug, Default)]
#[repr(C)]
pub struct vec3_t {
    pub x: f32,
    pub y: f32,
    pub z: f32,
}

pub fn v3(x: f32, y: f32, z: f32) -> vec3_t {
    vec3_t { x, y, z }
}

impl vec3_t {
    pub const ORIGIN: vec3_t = vec3_t {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    };

    pub fn from_scalar(f: f32) -> vec3_t {
        vec3_t { x: f, y: f, z: f }
    }

    pub fn dot(self, other: vec3_t) -> f32 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    pub fn scale(self, s: vec_t) -> vec3_t {
        vec3_t {
            x: self.x * s,
            y: self.y * s,
            z: self.z * s,
        }
    }

    pub fn length(self) -> vec_t {
        self.dot(self).sqrt()
    }

    pub fn length2(self) -> vec_t {
        self.dot(self)
    }

    pub fn map(f: impl Fn(f32) -> f32, a: vec3_t) -> vec3_t {
        vec3_t {
            x: f(a.x),
            y: f(a.y),
            z: f(a.z),
        }
    }

    pub fn map_2(f: impl Fn(f32, f32) -> f32, a: vec3_t, b: vec3_t) -> vec3_t {
        vec3_t {
            x: f(a.x, b.x),
            y: f(a.y, b.y),
            z: f(a.z, b.z),
        }
    }
    pub fn map_3(f: impl Fn(f32, f32, f32) -> f32, a: vec3_t, b: vec3_t, c: vec3_t) -> vec3_t {
        vec3_t {
            x: f(a.x, b.x, c.x),
            y: f(a.y, b.y, c.y),
            z: f(a.z, b.z, c.z),
        }
    }

    pub fn mid(self: vec3_t, b: vec3_t) -> vec3_t {
        (self + b) * 0.5
    }

    pub fn is_close_to(self, other: vec3_t, epsilon: vec_t) -> bool {
        let d = self - other;
        d.x.abs() < epsilon && d.y.abs() < epsilon && d.z.abs() < epsilon
    }
}

use std::ops::{Add, AddAssign, Index, IndexMut, Mul, Sub, SubAssign};

impl std::ops::Neg for vec3_t {
    type Output = vec3_t;
    fn neg(self) -> Self::Output {
        vec3_t {
            x: -self.x,
            y: -self.y,
            z: -self.z,
        }
    }
}

impl Index<usize> for vec3_t {
    type Output = vec_t;
    fn index(&self, index: usize) -> &vec_t {
        match index {
            0 => &self.x,
            1 => &self.y,
            2 => &self.z,
            _ => panic!(),
        }
    }
}

impl IndexMut<usize> for vec3_t {
    fn index_mut(&mut self, index: usize) -> &mut vec_t {
        match index {
            0 => &mut self.x,
            1 => &mut self.y,
            2 => &mut self.z,
            _ => panic!(),
        }
    }
}

impl Add<vec3_t> for vec3_t {
    type Output = vec3_t;
    fn add(self, other: vec3_t) -> vec3_t {
        vec3_t {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

impl AddAssign<vec3_t> for vec3_t {
    fn add_assign(&mut self, other: vec3_t) {
        self.x += other.x;
        self.y += other.y;
        self.z += other.z;
    }
}

impl Sub<vec3_t> for vec3_t {
    type Output = vec3_t;
    fn sub(self, other: vec3_t) -> vec3_t {
        vec3_t {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}

impl SubAssign<vec3_t> for vec3_t {
    fn sub_assign(&mut self, other: vec3_t) {
        self.x -= other.x;
        self.y -= other.y;
        self.z -= other.z;
    }
}

impl Mul<vec_t> for vec3_t {
    type Output = vec3_t;
    fn mul(self, other: vec_t) -> vec3_t {
        vec3_t {
            x: self.x * other,
            y: self.y * other,
            z: self.z * other,
        }
    }
}

impl Mul<vec3_t> for vec_t {
    type Output = vec3_t;
    fn mul(self, other: vec3_t) -> vec3_t {
        vec3_t {
            x: self * other.x,
            y: self * other.y,
            z: self * other.z,
        }
    }
}

pub fn DotProduct(a: vec3_t, b: vec3_t) -> f32 {
    a.dot(b)
}

// last arg was output
pub fn VectorSubtract(a: vec3_t, b: vec3_t) -> vec3_t {
    a - b
}

// last arg was output
pub fn VectorAdd(a: vec3_t, b: vec3_t) -> vec3_t {
    a + b
}

// was (a,b) where b was output
pub fn VectorCopy(a: vec3_t) -> vec3_t {
    a
}

// last arg was output
pub fn VectorScale(v: vec3_t, s: f32) -> vec3_t {
    v.scale(s)
}

// last arg was output
pub fn VectorMA(v: vec3_t, s: f32, b: vec3_t) -> vec3_t {
    vec3_t {
        x: v.x + b.x * s,
        y: v.y + b.y * s,
        z: v.z + b.z * s,
    }
}

#[derive(Copy, Clone, Debug)]
pub struct Vec3MinMax {
    pub mins: vec3_t,
    pub maxs: vec3_t,
}

pub fn VectorSet(x: vec_t, y: vec_t, z: vec_t) -> vec3_t {
    vec3_t { x, y, z }
}

pub fn VectorNormalize_mut(v: &mut vec3_t) -> f32 {
    let length = v.length();
    if length != 0.0 {
        let ilength = 1.0 / length;
        v.x *= ilength;
        v.y *= ilength;
        v.z *= ilength;
    }
    length
}

impl vec3_t {
    pub fn normalize_with_len(&self, v: vec3_t) -> Option<(vec3_t, vec_t)> {
        let length = v.length();
        if length != 0.0 {
            Some((self.scale(1.0 / length), length))
        } else {
            None
        }
    }
    pub fn normalize(&self) -> Option<vec3_t> {
        let length = self.length();
        if length != 0.0 {
            Some(self.scale(1.0 / length))
        } else {
            None
        }
    }
}

pub fn VectorNormalize2_mut(v: &mut vec3_t) -> f32 {
    let length = v.length();
    if length != 0.0 {
        let ilength = 1.0 / length;
        *v = v.scale(ilength);
    } else {
        *v = vec3_t::ORIGIN;
    }
    length
}

pub fn CrossProduct(v1: vec3_t, v2: vec3_t) -> vec3_t {
    vec3_t {
        x: v1.y * v2.z - v1.z * v2.y,
        y: v1.z * v2.x - v1.x * v2.z,
        z: v1.x * v2.y - v1.y * v2.x,
    }
}

pub fn VectorLength(v: vec3_t) -> vec_t {
    v.dot(v).sqrt()
}

#[derive(Copy, Clone, Debug)]
pub struct vec4_t {
    pub x: f32,
    pub y: f32,
    pub z: f32,
    pub w: f32,
}

impl vec4_t {
    pub fn from_vec3_w(v: vec3_t, w: vec_t) -> vec4_t {
        vec4_t {
            x: v.x,
            y: v.y,
            z: v.z,
            w: w,
        }
    }
}
