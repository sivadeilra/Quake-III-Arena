/*
===========================================================================
Copyright (C) 1999-2005 Id Software, Inc.

This file is part of Quake III Arena source code.

Quake III Arena source code is free software; you can redistribute it
and/or modify it under the terms of the GNU General Public License as
published by the Free Software Foundation; either version 2 of the License,
or (at your option) any later version.

Quake III Arena source code is distributed in the hope that it will be
useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Foobar; if not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
===========================================================================
*/

use crate::port_trace::*;
use crate::prelude::*;
use crate::qcommon::cm_patch::plane_t;
use log::warn;

pub type winding_t = Vec<vec3_t>;
pub type winding_slice = [vec3_t];

pub const MAX_POINTS_ON_WINDING: usize = 64;

pub const CLIP_EPSILON: f32 = 0.1;

pub const MAX_MAP_BOUNDS: f32 = 65535.0;

// you can define on_epsilon in the makefile as tighter
pub const ON_EPSILON: f32 = 0.1;

fn pw(w: &winding_t) {
    for p in w.iter() {
        println!("({:5.1}, {:5.1}, {:5.1})", p[0], p[1], p[2]);
    }
}

pub fn AllocWinding(points: usize) -> winding_t {
    Vec::with_capacity(points)
}

fn FreeWinding(w: winding_t) {
    drop(w)
}

fn RemoveColinearPoints(w: &mut winding_t) {
    let mut p: [vec3_t; MAX_POINTS_ON_WINDING] = [vec3_t::ORIGIN; MAX_POINTS_ON_WINDING];
    let mut nump: usize = 0;
    for (i, &w_i) in w.iter().enumerate() {
        let j = (i + 1) % w.len();
        let k = (i + w.len() - 1) % w.len();
        let v1 = (w[j] - w_i).normalize().unwrap_or(vec3_t::ORIGIN);
        let v2 = (w_i - w[k]).normalize().unwrap_or(vec3_t::ORIGIN);
        if v1.dot(v2) < 0.999 {
            p[nump] = w[i];
            nump += 1;
        }
    }

    if nump == w.len() {
        return;
    }

    w.clear();
    w.extend(&p[..nump]);
}

pub struct WindingPlaneResult {
    pub normal: vec3_t,
    pub dist: vec_t,
}

pub fn WindingPlane(w: &winding_slice) -> WindingPlaneResult {
    let v1 = w[1] - w[0];
    let v2 = w[2] - w[0];
    let normal = CrossProduct(v2, v1).normalize().unwrap();
    WindingPlaneResult {
        normal,
        dist: w[0].dot(normal),
    }
}

pub fn WindingArea(w: &winding_slice) -> vec_t {
    let mut total: vec_t = 0.0;
    for i in 2..w.len() {
        let d1 = w[i - 1] - w[0];
        let d2 = w[i] - w[0];
        let cross = CrossProduct(d1, d2);
        total += 0.5 * VectorLength(cross);
    }
    total
}

// was: (w, mins, maxs)
pub fn WindingBounds(w: &winding_t) -> Vec3MinMax {
    let mut mins = vec3_t([MAX_MAP_BOUNDS, MAX_MAP_BOUNDS, MAX_MAP_BOUNDS]);
    let mut maxs = vec3_t([-MAX_MAP_BOUNDS, -MAX_MAP_BOUNDS, -MAX_MAP_BOUNDS]);
    for p in w.iter() {
        mins = mins.min(*p);
        maxs = maxs.max(*p);
    }

    Vec3MinMax { mins, maxs }
}

pub fn WindingCenter(w: &winding_t) -> vec3_t {
    let mut center = vec3_t::ORIGIN;
    for p in w.iter() {
        center += *p;
    }

    let scale: f32 = 1.0 / (w.len() as f32);
    center.scale(scale)
}

pub fn BaseWindingForPlane(normal: vec3_t, dist: vec_t) -> winding_t {
    trace_str("BaseWindingForPlane");
    trace_vec3(normal);
    trace_f32(dist);

    // find the major axis
    /*
    let v_x = normal[0].abs();
    let v_y = normal[1].abs();
    let v_z = normal[2].abs();
    let vup = if v_x > v_y && v_x > v_z {
        // x is largest
        vec3_t::unit(2)
    } else if v_y > v_z {
        // y is largest
        vec3_t::unit(2)
    } else {
        // z is largest
        vec3_t::unit(0)
    };
    */
    let mut vup = vec3_t::ORIGIN;
    {
        let mut max: f32 = -MAX_MAP_BOUNDS;
        let mut x: i32 = -1;
        for i in 0..3 {
            let v = normal[i].abs();
            if v > max {
                x = i as i32;
                max = v;
            }
        }
        assert!(x != -1);
        match x {
            0 | 1 => vup[2] = 1.0,
            2 => vup[0] = 1.0,
            _ => panic!(),
        }
    }
    trace_str("vup");
    trace_vec3(vup);

    let org = normal.scale(dist);
    let v = vup.dot(normal);
    let mut vup = VectorMA(vup, -v, normal);
    vup = vup.normalize().unwrap_or(vec3_t::ORIGIN);
    let mut vright = CrossProduct(vup, normal);
    vup = vup.scale(MAX_MAP_BOUNDS);
    vright = vright.scale(MAX_MAP_BOUNDS);

    // project a really big axis aligned box onto the plane
    let mut w = Vec::with_capacity(4);
    w.push(org - vright + vup);
    w.push(org + vright + vup);
    w.push(org + vright - vup);
    w.push(org - vright - vup);
    trace_vec3(w[0]);
    trace_vec3(w[1]);
    trace_vec3(w[2]);
    trace_vec3(w[3]);
    trace_str(".");
    w
}

pub fn CopyWinding(w: &winding_t) -> winding_t {
    w.clone()
}

pub fn ReverseWinding(w: &winding_t) -> winding_t {
    w.iter().map(|&v| v).rev().collect()
}

struct ClipWindingInfo {
    dists: [f32; MAX_POINTS_ON_WINDING + 4],
    sides: [Side; MAX_POINTS_ON_WINDING + 4],
    counts: [usize; 3],
}
fn get_clip_winding_info(w: &winding_slice, plane: plane_t, epsilon: vec_t) -> ClipWindingInfo {
    let mut dists = [0.0f32; MAX_POINTS_ON_WINDING + 4];
    let mut sides = [Side(0); MAX_POINTS_ON_WINDING + 4];
    let mut counts: [usize; 3] = [0; 3];

    // determine sides for each point
    for (i, &p) in w.iter().enumerate() {
        let dot = plane.distance_to(p);
        dists[i] = dot;
        let this_side = if dot > epsilon {
            SIDE_FRONT
        } else if dot < -epsilon {
            SIDE_BACK
        } else {
            SIDE_ON
        };
        sides[i] = this_side;
        counts[this_side.index()] += 1;
    }
    sides[w.len()] = sides[0];
    dists[w.len()] = dists[0];
    ClipWindingInfo {
        dists,
        sides,
        counts,
    }
}

pub struct ClipWindingEpsilonResult {
    pub front: winding_t,
    pub back: winding_t,
}
pub fn ClipWindingEpsilon(
    in_: &winding_slice,
    plane: plane_t,
    epsilon: vec_t,
) -> ClipWindingEpsilonResult {
    type Output = ClipWindingEpsilonResult;

    let ClipWindingInfo {
        dists,
        sides,
        counts,
    } = get_clip_winding_info(in_, plane, epsilon);

    if counts[0] == 0 {
        return Output {
            front: Vec::new(),
            back: in_.to_vec(),
        };
    }
    if counts[1] == 0 {
        return Output {
            front: in_.to_vec(),
            back: Vec::new(),
        };
    }

    // cant use counts[0]+2 because
    // of fp grouping errors
    let maxpts = in_.len() + 4;

    let mut f = Vec::with_capacity(maxpts);
    let mut b = Vec::with_capacity(maxpts);

    for (i, &p1) in in_.iter().enumerate() {
        if sides[i] == SIDE_ON {
            f.push(p1);
            b.push(p1);
            continue;
        }
        if sides[i] == SIDE_FRONT {
            f.push(p1);
        }
        if sides[i] == SIDE_BACK {
            b.push(p1);
        }
        if sides[i + 1] == SIDE_ON || sides[i + 1] == sides[i] {
            continue;
        }

        // generate a split point
        let p2 = in_[(i + 1) % in_.len()];

        // avoid round off error when possible
        let dot = dists[i] / (dists[i] - dists[i + 1]);
        let mid = vec3_t::map_3(
            |normal_c, p1_c, p2_c| {
                if normal_c == 1.0 {
                    plane.dist
                } else if normal_c == -1.0 {
                    -plane.dist
                } else {
                    p1_c + dot * (p2_c - p1_c)
                }
            },
            plane.normal,
            p1,
            p2,
        );

        f.push(mid);
        b.push(mid);
    }

    if f.len() > maxpts || b.len() > maxpts {
        warn!("ClipWinding: points exceeded estimate");
    }
    if f.len() > MAX_POINTS_ON_WINDING || b.len() > MAX_POINTS_ON_WINDING {
        warn!("ClipWinding: MAX_POINTS_ON_WINDING");
    }

    Output { front: f, back: b }
}

pub fn ChopWindingInPlace(inout: &mut winding_t, plane: plane_t, epsilon: vec_t) {
    let ClipWindingInfo {
        dists,
        sides,
        counts,
    } = get_clip_winding_info(inout, plane, epsilon);

    if counts[0] == 0 {
        inout.clear();
        return;
    }
    if counts[1] == 0 {
        // inout stays the same
        return;
    }

    let in_ = inout;

    let maxpts = in_.len() + 4;
    // cant use counts[0]+2 because
    // of fp grouping errors

    let mut f = Vec::with_capacity(maxpts);

    for (i, &p1) in in_.iter().enumerate() {
        if sides[i] == SIDE_ON {
            f.push(p1);
            continue;
        }

        if sides[i] == SIDE_FRONT {
            f.push(p1);
        }

        if sides[i + 1] == SIDE_ON || sides[i + 1] == sides[i] {
            continue;
        }

        // generate a split point
        let p2 = in_[(i + 1) % in_.len()];

        // avoid round off error when possible
        let dot = dists[i] / (dists[i] - dists[i + 1]);
        let mid = vec3_t::map_3(
            |normal_c, p1_c, p2_c| {
                if normal_c == 1.0 {
                    plane.dist
                } else if normal_c == -1.0 {
                    -plane.dist
                } else {
                    p1_c + dot * (p2_c - p1_c)
                }
            },
            plane.normal,
            p1,
            p2,
        );
        f.push(mid);
    }

    if f.len() > maxpts {
        warn!("ClipWinding: points exceeded estimate");
    }
    if f.len() > MAX_POINTS_ON_WINDING {
        warn!("ClipWinding: MAX_POINTS_ON_WINDING");
    }

    *in_ = f;
}

/// Returns the fragment of in that is on the front side
/// of the cliping plane.  The original is freed.
pub fn ChopWinding(w: winding_t, plane: plane_t) -> winding_t {
    let r = ClipWindingEpsilon(&w, plane, ON_EPSILON);
    drop(w);
    drop(r.back);
    r.front
}

pub fn CheckWinding(w: &winding_t) {
    if w.len() < 3 {
        warn!("CheckWinding: {} points", w.len());
    }

    let area = WindingArea(w);
    if area < 1.0 {
        warn!("CheckWinding: {} area", area);
    }

    let WindingPlaneResult {
        normal: facenormal,
        dist: facedist,
    } = WindingPlane(w);

    for (i, &p1) in w.iter().enumerate() {
        for j in 0..3 {
            if p1[j] > MAX_MAP_BOUNDS || p1[j] < -MAX_MAP_BOUNDS {
                warn!("CheckFace: BUGUS_RANGE: {}", p1[j]);
            }
        }

        // check the point is on the face plane
        let d = p1.dot(facenormal) - facedist;
        if d < -ON_EPSILON || d > ON_EPSILON {
            warn!("CheckWinding: point off plane");
        }

        // check the edge isnt degenerate
        let p2 = w[if i + 1 == w.len() { 0 } else { i + 1 }];
        let dir = p2 - p1;

        if dir.length() < ON_EPSILON {
            warn!("CheckWinding: degenerate edge");
        }

        let edgenormal = CrossProduct(facenormal, dir).normalize().unwrap();
        let edgedist = p1.dot(edgenormal) + ON_EPSILON;

        // all other points must be on front side
        for j in 0..w.len() {
            if j == i {
                continue;
            }
            let d = w[j].dot(edgenormal);
            if d > edgedist {
                warn!("CheckWinding: non-convex");
            }
        }
    }
}

pub use winding_on_plane_side as WindingOnPlaneSide;

pub fn winding_on_plane_side(w: &winding_slice, normal: vec3_t, dist: vec_t) -> Side {
    let mut front = false;
    let mut back = false;
    for &p in w.iter() {
        let d = p.dot(normal) - dist;
        if d < -ON_EPSILON {
            if front {
                return SIDE_CROSS;
            }
            back = true;
            continue;
        }
        if d > ON_EPSILON {
            if back {
                return SIDE_CROSS;
            }
            front = true;
            continue;
        }
    }

    if back {
        SIDE_BACK
    } else if front {
        SIDE_FRONT
    } else {
        SIDE_ON
    }
}

const MAX_HULL_POINTS: usize = 128;

pub use add_winding_to_convex_hull as AddWindingToConvexHull;

/// Both w and *hull are on the same plane
pub fn add_winding_to_convex_hull(w: &winding_t, hull: &mut winding_t, normal: vec3_t) {
    let mut hull_points = [vec3_t::ORIGIN; MAX_HULL_POINTS];
    let mut new_hull_points = [vec3_t::ORIGIN; MAX_HULL_POINTS];
    let mut hull_dirs = [vec3_t::ORIGIN; MAX_HULL_POINTS];
    let mut hull_side: [bool; MAX_HULL_POINTS] = [false; MAX_HULL_POINTS];

    if hull.len() == 0 {
        hull.extend(w);
    }

    let mut num_hull_points = hull.len();
    hull_points[..num_hull_points].copy_from_slice(&hull);

    for &p in w.iter() {
        // calculate hull side vectors
        for j in 0..num_hull_points {
            let k = (j + 1) % num_hull_points;
            let dir = (hull_points[k] - hull_points[j]).normalize().unwrap();
            hull_dirs[j] = CrossProduct(normal, dir);
        }

        let mut outside = false;
        for j in 0..num_hull_points {
            let dir = p - hull_points[j];
            let d = dir.dot(hull_dirs[j]);
            if d >= ON_EPSILON {
                outside = true;
            }
            hull_side[j] = d >= -ON_EPSILON;
        }

        // if the point is effectively inside, do nothing
        if !outside {
            continue;
        }

        // find the back side to front side transition
        let mut j: usize = 0;
        while j < num_hull_points {
            if !hull_side[j % num_hull_points] && hull_side[(j + 1) % num_hull_points] {
                break;
            }
            j += 1;
        }
        if j == num_hull_points {
            continue;
        }

        // insert the point here

        new_hull_points[0] = p;
        let mut numNew: usize = 1;

        // copy over all points that aren't double fronts
        let j = (j + 1) % num_hull_points;
        for k in 0..num_hull_points {
            if hull_side[(j + k) % num_hull_points] && hull_side[(j + k + 1) % num_hull_points] {
                continue;
            }
            new_hull_points[numNew] = hull_points[(j + k + 1) % num_hull_points];
            numNew += 1;
        }

        num_hull_points = numNew;
        hull_points[..num_hull_points].copy_from_slice(&new_hull_points[..num_hull_points]);
    }

    *hull = hull_points[..num_hull_points].to_vec();
}
