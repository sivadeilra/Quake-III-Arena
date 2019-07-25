use crate::num_utils::{fmax, fmin};

pub type vec_t = f32;

// see game/q_shared.h

#[derive(Copy, Clone, Debug, Default, PartialEq)]
#[repr(C)]
pub struct vec3_t(pub [f32; 3]);

pub const vec3_origin: vec3_t = vec3_t([0.0, 0.0, 0.0]);

pub fn v3(x: f32, y: f32, z: f32) -> vec3_t {
    vec3_t([x, y, z])
}

impl vec3_t {
    pub const ORIGIN: vec3_t = vec3_t([0.0, 0.0, 0.0]);

    pub fn from_scalar(f: f32) -> vec3_t {
        vec3_t([f, f, f])
    }

    pub fn dot(self, other: vec3_t) -> f32 {
        self.0[0] * other.0[0] + self.0[1] * other.0[1] + self.0[2] * other.0[2]
    }

    pub fn scale(self, s: vec_t) -> vec3_t {
        vec3_t([self.0[0] * s, self.0[1] * s, self.0[2] * s])
    }

    pub fn length(self) -> vec_t {
        self.dot(self).sqrt()
    }

    pub fn length2(self) -> vec_t {
        self.dot(self)
    }

    pub fn unit(axis: u8) -> vec3_t {
        match axis {
            0 => v3(1.0, 0.0, 0.0),
            1 => v3(0.0, 1.0, 0.0),
            2 => v3(0.0, 0.0, 1.0),
            _ => panic!("invalid"),
        }
    }

    pub fn map_to_array<T>(self, f: impl Fn(f32) -> T) -> [T; 3] {
        [f(self.0[0]), f(self.0[1]), f(self.0[2])]
    }

    pub fn map(self, f: impl Fn(f32) -> f32) -> vec3_t {
        vec3_t(self.map_to_array(f))
    }

    pub fn map2(self, other: vec3_t, f: impl Fn(f32, f32) -> f32) -> vec3_t {
        vec3_t([
            f(self.0[0], other.0[0]),
            f(self.0[1], other.0[1]),
            f(self.0[2], other.0[2]),
        ])
    }

    pub fn map_3(f: impl Fn(f32, f32, f32) -> f32, a: vec3_t, b: vec3_t, c: vec3_t) -> vec3_t {
        vec3_t([
            f(a[0], b[0], c[0]),
            f(a[1], b[1], c[1]),
            f(a[2], b[2], c[2]),
        ])
    }

    pub fn mid(self: vec3_t, b: vec3_t) -> vec3_t {
        (self + b) * 0.5
    }

    pub fn is_close_to(self, other: vec3_t, epsilon: vec_t) -> bool {
        let d = self - other;
        d[0].abs() < epsilon && d[1].abs() < epsilon && d[2].abs() < epsilon
    }

    pub fn min(self, other: vec3_t) -> vec3_t {
        self.map2(other, fmin)
    }
    pub fn max(self, other: vec3_t) -> vec3_t {
        self.map2(other, fmax)
    }
}

impl std::ops::Neg for vec3_t {
    type Output = vec3_t;
    fn neg(self) -> Self::Output {
        vec3_t([-self.0[0], -self.0[1], -self.0[2]])
    }
}

impl core::ops::Index<usize> for vec3_t {
    type Output = vec_t;
    fn index(&self, index: usize) -> &vec_t {
        &self.0[index]
    }
}

impl core::ops::IndexMut<usize> for vec3_t {
    fn index_mut(&mut self, index: usize) -> &mut vec_t {
        &mut self.0[index]
    }
}

impl core::ops::Add<vec3_t> for vec3_t {
    type Output = vec3_t;
    fn add(self, other: vec3_t) -> vec3_t {
        vec3_t([self[0] + other[0], self[1] + other[1], self[2] + other[2]])
    }
}

impl core::ops::AddAssign<vec3_t> for vec3_t {
    fn add_assign(&mut self, other: vec3_t) {
        self.0[0] += other.0[0];
        self.0[1] += other.0[1];
        self.0[2] += other.0[2];
    }
}

impl core::ops::Sub<vec3_t> for vec3_t {
    type Output = vec3_t;
    fn sub(self, other: vec3_t) -> vec3_t {
        vec3_t([self[0] - other[0], self[1] - other[1], self[2] - other[2]])
    }
}

impl core::ops::SubAssign<vec3_t> for vec3_t {
    fn sub_assign(&mut self, other: vec3_t) {
        self.0[0] -= other.0[0];
        self.0[1] -= other.0[1];
        self.0[2] -= other.0[2];
    }
}

impl core::ops::Mul<vec_t> for vec3_t {
    type Output = vec3_t;
    fn mul(self, other: vec_t) -> vec3_t {
        vec3_t([self[0] * other, self[1] * other, self[2] * other])
    }
}

impl core::ops::Mul<vec3_t> for vec_t {
    type Output = vec3_t;
    fn mul(self, other: vec3_t) -> vec3_t {
        vec3_t([self * other[0], self * other[1], self * other[2]])
    }
}

impl core::ops::MulAssign<vec_t> for vec3_t {
    fn mul_assign(&mut self, other: vec_t) {
        self[0] *= other;
        self[1] *= other;
        self[2] *= other;
    }
}

impl core::convert::From<[vec_t; 3]> for vec3_t {
    fn from(v: [vec_t; 3]) -> vec3_t {
        vec3_t(v)
    }
}

impl core::convert::From<vec3_t> for [vec_t; 3] {
    fn from(v: vec3_t) -> [vec_t; 3] {
        v.0
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
    v + b * s
}

#[derive(Copy, Clone, Debug)]
pub struct Vec3MinMax {
    pub mins: vec3_t,
    pub maxs: vec3_t,
}

pub fn VectorSet(x: vec_t, y: vec_t, z: vec_t) -> vec3_t {
    vec3_t([x, y, z])
}

pub fn VectorNormalize_mut(v: &mut vec3_t) -> f32 {
    let length = v.length();
    if length != 0.0 {
        let ilength = 1.0 / length;
        *v *= ilength;
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
    vec3_t([
        v1[1] * v2[2] - v1[2] * v2[1],
        v1[2] * v2[0] - v1[0] * v2[2],
        v1[0] * v2[1] - v1[1] * v2[0],
    ])
}

pub fn VectorLength(v: vec3_t) -> vec_t {
    v.dot(v).sqrt()
}

pub fn get_sign_bits(v: vec3_t) -> u8 {
    ((v[0] < 0.0) as u8) | (((v[1] < 0.0) as u8) << 1) | (((v[2] < 0.0) as u8) << 2) //
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
            x: v[0],
            y: v[1],
            z: v[2],
            w: w,
        }
    }
}
