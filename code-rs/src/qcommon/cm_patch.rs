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

#![allow(non_upper_case_globals)]

use crate::is_close_to;
use crate::port_trace::*;
use crate::prelude::*;
use crate::qcommon::cm_load::Error;
use crate::qcommon::cm_local::*;
use crate::qcommon::cm_polylib::*;
use crate::qcommon::cm_trace::traceWork_t;
use log::warn;
use std::mem::swap;

//#define   CULL_BBOX

/*

Issues for collision against curved surfaces:

Surface edges need to be handled differently than surface planes

Plane expansion causes raw surfaces to expand past expanded bounding box

Position test of a volume against a surface is tricky.

Position test of a point against a surface is not well defined, because the surface has no volume.


Tracing leading edge points instead of volumes?
Position test by tracing corner to corner? (8*7 traces -- ouch)

coplanar edges
triangulated patches
degenerate patches

  endcaps
  degenerate

WARNING: this may misbehave with meshes that have rows or columns that only
degenerate a few triangles.  Completely degenerate rows and columns are handled
properly.
*/

pub const MAX_FACETS: usize = 1024;
pub const MAX_PATCH_PLANES: usize = 2048;

#[derive(Copy, Clone, Debug, Default)]
pub struct plane_t {
    pub normal: vec3_t,
    pub dist: vec_t,
}
impl plane_t {
    pub fn flip(self) -> plane_t {
        plane_t {
            normal: -self.normal,
            dist: -self.dist,
        }
    }

    pub fn distance_to(&self, v: vec3_t) -> f32 {
        self.normal.dot(v) - self.dist
    }
}

#[derive(Clone, Debug)]
pub struct patchPlane_t {
    pub plane: plane_t, // [f32; 4],
    pub signbits: i32,  // signx + (signy<<1) + (signz<<2), used as lookup during collision
}

#[derive(Clone, Debug, Default)]
pub struct facet_t {
    pub surfacePlane: i32,
    pub numBorders: i32, // 3 or four + 6 axial bevels + 4 or 3 * 4 edge bevels
    pub borderPlanes: [i32; 4 + 6 + 16],
    pub borderInward: [bool; 4 + 6 + 16],
    pub borderNoAdjust: [bool; 4 + 6 + 16],
}

impl facet_t {
    pub fn iter_border_planes(&self) -> impl Iterator<Item = i32> + '_ {
        self.borderPlanes[..self.numBorders as usize]
            .iter()
            .copied()
    }
}

pub type patchCollide_s = patchCollide_t;
#[derive(Clone, Debug)]
pub struct patchCollide_t {
    pub bounds: vec3_bounds,

    // surface planes plus edge planes
    pub planes: Vec<patchPlane_t>,

    pub facets: Vec<facet_t>,
}

pub const MAX_GRID_SIZE: usize = 129;

#[derive(Clone)]
pub struct cGrid_t {
    pub width: usize,
    pub height: usize,
    pub wrapWidth: bool,
    pub wrapHeight: bool,
    // [width][height]
    pub points: [[vec3_t; MAX_GRID_SIZE]; MAX_GRID_SIZE],
}
impl cGrid_t {
    pub fn empty() -> Self {
        Self {
            width: 0,
            height: 0,
            wrapWidth: false,
            wrapHeight: false,
            points: [[vec3_t::ORIGIN; MAX_GRID_SIZE]; MAX_GRID_SIZE],
        }
    }
}

pub const SUBDIVIDE_DISTANCE: f32 = 16.0; //4 // never more than this units away from curve
pub const PLANE_TRI_EPSILON: f32 = 0.1;
pub const WRAP_POINT_EPSILON: f32 = 0.1;

type GridPlanes = [[[i32; 2]; MAX_GRID_SIZE]; MAX_GRID_SIZE];

use crate::perf::StaticCounter;

static c_totalPatchBlocks: StaticCounter = StaticCounter::new("total patch blocks");
static c_totalPatchSurfaces: StaticCounter = StaticCounter::new("total patch surfaces");
static c_totalPatchEdges: StaticCounter = StaticCounter::new("total patch edges");

/*
static const patchCollide_t *debugPatchCollide;
static const facet_t        *debugFacet;
static bool     debugBlock;
static vec3_t       debugBlockPoints[4];
*/

/*
=================
CM_ClearLevelPatches
=================
*/
pub fn CM_ClearLevelPatches() {
    /*
    debugPatchCollide = NULL;
    debugFacet = NULL;
    */
}

/*
=================
CM_SignbitsForNormal
=================
*/
fn CM_SignbitsForNormal(normal: vec3_t) -> i32 {
    get_sign_bits(normal) as i32
}

/*
=====================
CM_PlaneFromPoints

Returns false if the triangle is degenrate.
The normal will point out of the clock for clockwise ordered points
=====================
*/
fn CM_PlaneFromPoints(a: vec3_t, b: vec3_t, c: vec3_t) -> Option<plane_t> {
    trace_str("CM_PlaneFromPoints");

    let d1 = b - a;
    let d2 = c - a;
    let normal = if let Some(cp) = CrossProduct(d2, d1).normalize() {
        cp
    } else {
        trace_str("cross product is zero");
        return None;
    };
    trace_vec3(normal);
    let dist = a.dot(normal);
    trace_f32(dist);
    Some(plane_t { normal, dist })
}

/*
================================================================================

GRID SUBDIVISION

================================================================================
*/

/*
=================
CM_NeedsSubdivision

Returns true if the given quadratic curve is not flat enough for our
collision detection purposes
=================
*/
fn CM_NeedsSubdivision(a: vec3_t, b: vec3_t, c: vec3_t) -> bool {
    // calculate the linear midpoint
    let lmid = a.mid(c);

    let ab_mid = a.mid(b);
    let bc_mid = b.mid(c);

    // calculate the exact curve midpoint
    let cmid = vec3_t::mid(ab_mid, bc_mid);

    // see if the curve is far enough away from the linear mid
    let delta = cmid - lmid;
    let dist = delta.length();
    dist >= SUBDIVIDE_DISTANCE
}

/*
===============
CM_Subdivide

a, b, and c are control points.
the subdivided sequence will be: a, out1, out2, out3, c
===============
*/
fn CM_Subdivide(a: vec3_t, b: vec3_t, c: vec3_t) -> (vec3_t, vec3_t, vec3_t) {
    let out1 = vec3_t::mid(a, b);
    let out3 = vec3_t::mid(b, c);
    let out2 = vec3_t::mid(out1, out3);
    (out1, out2, out3)
}

/*
=================
CM_TransposeGrid

Swaps the rows and columns in place
=================
*/
fn CM_TransposeGrid(grid: &mut cGrid_t) {
    let width = grid.width;
    let height = grid.height;
    let points = &mut grid.points;

    if width > height {
        for i in 0..height {
            for j in i + 1..width {
                if j < height {
                    // swap the value
                    let temp = points[i][j];
                    points[i][j] = points[j][i];
                    points[j][i] = temp;
                } else {
                    // just copy
                    points[i][j] = points[j][i];
                }
            }
        }
    } else {
        for i in 0..width {
            for j in i + 1..height {
                if j < width {
                    // swap the value
                    let temp = points[j][i];
                    points[j][i] = points[i][j];
                    points[i][j] = temp;
                } else {
                    // just copy
                    points[j][i] = points[i][j];
                }
            }
        }
    }

    swap(&mut grid.width, &mut grid.height);
    swap(&mut grid.wrapWidth, &mut grid.wrapHeight);
}

/*
===================
CM_SetGridWrapWidth

If the left and right columns are exactly equal, set grid.wrapWidth true
===================
*/
fn CM_SetGridWrapWidth(grid: &mut cGrid_t) {
    let wrap_width = compute_grid_wrap_width(grid);
    grid.wrapWidth = wrap_width;
}

fn compute_grid_wrap_width(grid: &cGrid_t) -> bool {
    let height = grid.height;
    let width = grid.width;
    let points = &grid.points;
    for i in 0..height {
        let d_v = points[0][i] - points[width - 1][i];
        for j in 0..3 {
            let d = d_v[j];
            if d < -WRAP_POINT_EPSILON || d > WRAP_POINT_EPSILON {
                return false;
            }
        }
    }

    return true;
}

/*
=================
CM_SubdivideGridColumns

Adds columns as necessary to the grid until
all the aproximating points are within SUBDIVIDE_DISTANCE
from the true curve
=================
*/
fn CM_SubdivideGridColumns(grid: &mut cGrid_t) {
    let mut i: usize = 0;
    while i < grid.width - 2 {
        // grid.points[i][x] is an interpolating control point
        // grid.points[i+1][x] is an aproximating control point
        // grid.points[i+2][x] is an interpolating control point

        //
        // first see if we can collapse the aproximating collumn away
        //
        let any_needs_subdivision = (0..grid.height).any(|j| {
            CM_NeedsSubdivision(
                grid.points[i][j],
                grid.points[i + 1][j],
                grid.points[i + 2][j],
            )
        });
        if !any_needs_subdivision {
            // all of the points were close enough to the linear midpoints
            // that we can collapse the entire column away
            for j in 0..grid.height {
                // remove the column
                for k in i + 2..grid.width {
                    grid.points[k - 1][j] = grid.points[k][j];
                }
            }
            grid.width -= 1;
            // go to the next curve segment
            i += 1;
            continue;
        }

        //
        // we need to subdivide the curve
        //
        for j in 0..grid.height {
            // save the control points now
            let prev = grid.points[i][j];
            let mid = grid.points[i + 1][j];
            let next = grid.points[i + 2][j];

            // make room for two additional columns in the grid
            // columns i+1 will be replaced, column i+2 will become i+4
            // i+1, i+2, and i+3 will be generated
            for k in (i..grid.width).rev() {
                grid.points[k + 2][j] = grid.points[k][j];
            }

            // generate the subdivided points
            let (new1, new2, new3) = CM_Subdivide(prev, mid, next);
            grid.points[i + 1][j] = new1;
            grid.points[i + 2][j] = new2;
            grid.points[i + 3][j] = new3;
        }

        grid.width += 2;

        // the new aproximating point at i+1 may need to be removed
        // or subdivided farther, so don't advance i
    }
}

pub const POINT_EPSILON: f32 = 0.1;

fn near_zero(f: f32) -> bool {
    f <= POINT_EPSILON && f >= -POINT_EPSILON
}

fn CM_ComparePoints(a: vec3_t, b: vec3_t) -> bool {
    let vd = a - b;
    let nears = vd.map_to_array(near_zero);
    nears[0] && nears[1] && nears[2]
}

/// If there are any identical columns, remove them
fn CM_RemoveDegenerateColumns(grid: &mut cGrid_t) {
    let mut i: usize = 0;
    while i < grid.width - 1 {
        let row_i = &grid.points[i];
        let row_next = &grid.points[i + 1];
        if (0..grid.height).any(|j| !CM_ComparePoints(row_i[j], row_next[j])) {
            // at least one pair of points was not close to each other
            i += 1;
            continue;
        }

        for j in 0..grid.height {
            // remove the column
            for k in i + 2..grid.width {
                grid.points[k - 1][j] = grid.points[k][j];
            }
        }
        grid.width -= 1;
        // check against the next column; do not advance i
    }
}

/*
================================================================================

PATCH COLLIDE GENERATION

================================================================================
*/

const NORMAL_EPSILON: f32 = 0.0001;
const DIST_EPSILON: f32 = 0.02;

pub struct PlaneEqualOutput {
    pub is_flipped: bool,
    pub is_equal: bool,
}
pub fn CM_PlaneEqual(a: plane_t, b: plane_t) -> PlaneEqualOutput {
    trace_str("CM_PlaneEqual");
    trace_vec3(a.normal);
    trace_f32(a.dist);
    trace_vec3(b.normal);
    trace_f32(b.dist);
    fn is_close_to(a: plane_t, b: plane_t) -> bool {
        a.normal.is_close_to(b.normal, NORMAL_EPSILON) && (a.dist - b.dist).abs() < DIST_EPSILON
    }

    if is_close_to(a, b) {
        trace_str("planes are close");
        return PlaneEqualOutput {
            is_flipped: false,
            is_equal: true,
        };
    }

    let inv_b = plane_t {
        normal: -b.normal,
        dist: -b.dist,
    };
    if is_close_to(a, inv_b) {
        trace_str("planes are close (flipped)");
        return PlaneEqualOutput {
            is_flipped: true,
            is_equal: true,
        };
    }

    trace_str("planes are not equal");
    PlaneEqualOutput {
        is_flipped: false,
        is_equal: false,
    }
}

pub fn CM_SnapVector(normal: &mut vec3_t) {
    let snapped = CM_SnapVector_func(*normal);
    *normal = snapped;
}
pub fn CM_SnapVector_func(normal: vec3_t) -> vec3_t {
    for i in 0..3 {
        if (normal[i] - 1.0).abs() < NORMAL_EPSILON {
            let mut r = vec3_t::ORIGIN;
            r[i] = 1.0;
            return r;
        }
        if (normal[i] - -1.0).abs() < NORMAL_EPSILON {
            let mut r = vec3_t::ORIGIN;
            r[i] = -1.0;
            return r;
        }
    }

    normal
}

pub fn CM_FindPlane2(
    planes: &mut Vec<patchPlane_t>,
    plane: plane_t,
) -> (/*index:*/ usize, /*is_flipped:*/ bool) {
    // see if the points are close enough to an existing plane
    for (i, p) in planes.iter().enumerate() {
        let PlaneEqualOutput {
            is_flipped,
            is_equal,
        } = CM_PlaneEqual(p.plane, plane);
        if is_equal {
            return (i, is_flipped);
        }
    }

    // add a new plane
    let index = planes.len();
    add_plane(planes, plane);
    (index, false)
}

fn add_plane(planes: &mut Vec<patchPlane_t>, plane: plane_t) {
    assert!(planes.len() < MAX_PATCH_PLANES);
    trace_str("adding plane at index:");
    trace_i32(planes.len() as i32);
    trace_vec3(plane.normal);
    trace_f32(plane.dist);
    trace_str(".");
    planes.push(patchPlane_t {
        plane: plane,
        signbits: CM_SignbitsForNormal(plane.normal),
    });
}

fn CM_FindPlane(
    planes: &mut Vec<patchPlane_t>,
    p1: vec3_t,
    p2: vec3_t,
    p3: vec3_t,
) -> Result<i32, Error> {
    trace_str("CM_FindPlane");
    trace_vec3(p1);
    trace_vec3(p2);
    trace_vec3(p3);

    let new_plane = if let Some(new_plane) = CM_PlaneFromPoints(p1, p2, p3) {
        new_plane
    } else {
        trace_str("points do not form plane");
        return Ok(-1);
    };
    trace_vec3(new_plane.normal);
    trace_f32(new_plane.dist);

    // see if the points are close enough to an existing plane
    for (i, p) in planes.iter().enumerate() {
        trace_str("testing plane:");
        trace_i32(i as i32);
        trace_vec3(p.plane.normal);
        trace_f32(p.plane.dist);
        let p = p.plane;
        if new_plane.normal.dot(p.normal) < 0.0 {
            trace_str("negative dot product");
            continue; // allow backwards planes?
        }
        trace_str("checking points");
        let point_is_good =
            move |pN: vec3_t| is_close_to(pN.dot(p.normal), p.dist, PLANE_TRI_EPSILON);

        if point_is_good(p1) && point_is_good(p2) && point_is_good(p3) {
            trace_str("using existing plane:");
            trace_i32(i as i32);
            return Ok(i as i32);
        } else {
            trace_str("skipping because point is too far from plane");
        }
    }

    // add a new plane
    let index = planes.len();
    add_plane(planes, new_plane);
    Ok(index as i32)
}

fn CM_PointOnPlaneSide(planes: &[patchPlane_t], p: vec3_t, planeNum: i32) -> Side {
    if planeNum == -1 {
        return SIDE_ON;
    }

    let plane = planes[planeNum as usize].plane;
    let d = plane.distance_to(p);

    if d > PLANE_TRI_EPSILON {
        SIDE_FRONT
    } else if d < -PLANE_TRI_EPSILON {
        SIDE_BACK
    } else {
        SIDE_ON
    }
}

fn CM_GridPlane(gridPlanes: &GridPlanes, i: usize, j: usize, tri: usize) -> Result<usize, Error> {
    trace_str("CM_GridPlane");

    let p = gridPlanes[i][j][tri];
    if p != -1 {
        trace_str("returning first grid plane");
        trace_i32(p);
        return Ok(p as usize);
    }
    let p = gridPlanes[i][j][tri ^ 1];
    if p != -1 {
        trace_str("returning second grid plane");
        trace_i32(p);
        return Ok(p as usize);
    }

    // should never happen
    warn!("WARNING: CM_GridPlane unresolvable");
    Err(Error::Str("CM_GridPlane unresolvable"))
}

fn CM_EdgePlaneNum(
    planes: &mut Vec<patchPlane_t>,
    grid: &cGrid_t,
    gridPlanes: &GridPlanes,
    i: usize,
    j: usize,
    k: usize,
) -> Result<i32, Error> {
    let up;
    let p;
    let p1;
    let p2;

    trace_str("CM_EdgePlaneNum");
    trace_i32(i as i32);
    trace_i32(j as i32);
    trace_i32(k as i32);

    match k {
        0 => {
            // top border
            p1 = grid.points[i][j];
            p2 = grid.points[i + 1][j];
            p = CM_GridPlane(gridPlanes, i, j, 0)?;
            up = VectorMA(p1, 4.0, planes[p].plane.normal);
            CM_FindPlane(planes, p1, p2, up)
        }

        2 => {
            // bottom border
            p1 = grid.points[i][j + 1];
            p2 = grid.points[i + 1][j + 1];
            p = CM_GridPlane(gridPlanes, i, j, 1)?;
            up = VectorMA(p1, 4.0, planes[p].plane.normal);
            CM_FindPlane(planes, p2, p1, up)
        }

        3 => {
            // left border
            p1 = grid.points[i][j];
            p2 = grid.points[i][j + 1];
            p = CM_GridPlane(gridPlanes, i, j, 1)?;
            up = VectorMA(p1, 4.0, planes[p].plane.normal);
            CM_FindPlane(planes, p2, p1, up)
        }

        1 => {
            // right border
            p1 = grid.points[i + 1][j];
            p2 = grid.points[i + 1][j + 1];
            p = CM_GridPlane(gridPlanes, i, j, 0)?;
            up = VectorMA(p1, 4.0, planes[p].plane.normal);
            CM_FindPlane(planes, p1, p2, up)
        }

        4 => {
            // diagonal out of triangle 0
            p1 = grid.points[i + 1][j + 1];
            p2 = grid.points[i][j];
            p = CM_GridPlane(gridPlanes, i, j, 0)?;
            up = VectorMA(p1, 4.0, planes[p].plane.normal);
            CM_FindPlane(planes, p1, p2, up)
        }

        5 => {
            // diagonal out of triangle 1
            p1 = grid.points[i][j];
            p2 = grid.points[i + 1][j + 1];
            p = CM_GridPlane(gridPlanes, i, j, 1)?;
            up = VectorMA(p1, 4.0, planes[p].plane.normal);
            CM_FindPlane(planes, p1, p2, up)
        }
        _ => panic!(),
    }
    // plane_num can be -1, and that's OK, and is different from an error.
}

fn CM_SetBorderInward(
    planes: &[patchPlane_t],
    facet: &mut facet_t,
    grid: &cGrid_t,
    _gridPlanes: &GridPlanes,
    i: usize,
    j: usize,
    which: i32,
) {
    let points: [vec3_t; 4];
    let numPoints: usize;
    match which {
        -1 => {
            points = [
                grid.points[i][j],
                grid.points[i + 1][j],
                grid.points[i + 1][j + 1],
                grid.points[i][j + 1],
            ];
            numPoints = 4;
        }
        0 => {
            points = [
                grid.points[i][j],
                grid.points[i + 1][j],
                grid.points[i + 1][j + 1],
                vec3_t::ORIGIN,
            ];
            numPoints = 3;
        }
        1 => {
            points = [
                grid.points[i + 1][j + 1],
                grid.points[i][j + 1],
                grid.points[i][j],
                vec3_t::ORIGIN,
            ];
            numPoints = 3;
        }
        _ => panic!("CM_SetBorderInward: bad parameter"),
    }

    for k in 0..facet.numBorders as usize {
        let mut front = 0;
        let mut back = 0;
        for &p in points[..numPoints].iter() {
            match CM_PointOnPlaneSide(planes, p, facet.borderPlanes[k]) {
                SIDE_FRONT => front += 1,
                SIDE_BACK => back += 1,
                _ => {}
            }
        }

        if front != 0 && back == 0 {
            facet.borderInward[k] = true;
        } else if back != 0 && front == 0 {
            facet.borderInward[k] = false;
        } else if front == 0 && back == 0 {
            // flat side border
            facet.borderPlanes[k] = -1;
        } else {
            // bisecting side border
            warn!("CM_SetBorderInward: mixed plane sides");
            facet.borderInward[k] = false;
            /*
            if ( !debugBlock ) {
                debugBlock = true;
                debugBlockPoints[0] = grid.points[i][j];
                debugBlockPoints[1] = grid.points[i+1][j]);
                debugBlockPoints[2] = grid.points[i+1][j+1]);
                debugBlockPoints[3] = grid.points[i][j+1]);
            }
            */
        }
    }
}

/// If the facet isn't bounded by its borders, we screwed up.
fn CM_ValidateFacet(planes: &[patchPlane_t], facet: &facet_t) -> bool {
    if facet.surfacePlane == -1 {
        return false;
    }

    let plane = planes[facet.surfacePlane as usize].plane;
    let mut w = BaseWindingForPlane(plane.normal, plane.dist);
    for j in 0..facet.numBorders as usize {
        if facet.borderPlanes[j] == -1 {
            return false;
        }
        let mut plane = planes[facet.borderPlanes[j] as usize].plane;
        if !facet.borderInward[j] {
            plane.normal = -plane.normal;
            plane.dist = -plane.dist;
        }
        ChopWindingInPlace(&mut w, plane.normal, plane.dist, 0.1);
    }

    if w.is_empty() {
        return false; // winding was completely chopped away
    }

    // see if the facet is unreasonably large
    let bounds = WindingBounds(&w);
    for j in 0..3 {
        if bounds.maxs[j] - bounds.mins[j] > MAX_MAP_BOUNDS {
            return false; // we must be missing a plane
        }
        if bounds.mins[j] >= MAX_MAP_BOUNDS {
            return false;
        }
        if bounds.maxs[j] <= -MAX_MAP_BOUNDS {
            return false;
        }
    }
    // winding is fine
    return true;
}

fn CM_AddFacetBevels(planes: &mut Vec<patchPlane_t>, facet: &mut facet_t) {
    trace_str("CM_AddFacetBevels");
    let mut w = {
        let plane = planes[facet.surfacePlane as usize].plane;
        BaseWindingForPlane(plane.normal, plane.dist)
    };
    for j in 0..facet.numBorders as usize {
        if facet.borderPlanes[j] == facet.surfacePlane {
            continue;
        }
        let mut plane = planes[facet.borderPlanes[j] as usize].plane;
        if !facet.borderInward[j] {
            plane = plane.flip();
        }

        ChopWindingInPlace(&mut w, plane.normal, plane.dist, 0.1);
    }
    if w.is_empty() {
        trace_str("winding is empty");
        return;
    }

    let bounds = WindingBounds(&w);
    trace_str("winding bounds");
    trace_vec3(bounds.mins);
    trace_vec3(bounds.maxs);

    // add the axial planes
    trace_str("add the axial planes");
    for axis in 0..3 {
        for &dir in [-1, 1].into_iter() {
            let mut plane = plane_t {
                normal: vec3_t::ORIGIN,
                dist: if dir == 1 {
                    bounds.maxs[axis]
                } else {
                    -bounds.mins[axis]
                },
            };
            plane.normal[axis] = dir as f32;

            // if it's the surface plane
            if CM_PlaneEqual(planes[facet.surfacePlane as usize].plane, plane).is_equal {
                continue;
            }
            // see if the plane is allready present
            let any_plane_equal = facet
                .iter_border_planes()
                .any(|plane_num| CM_PlaneEqual(planes[plane_num as usize].plane, plane).is_equal);
            if !any_plane_equal {
                if facet.numBorders > 4 + 6 + 16 {
                    warn!("ERROR: too many bevels\n");
                }
                let (plane_index, flipped) = CM_FindPlane2(planes, plane);
                let border_index = facet.numBorders as usize;
                facet.borderPlanes[border_index] = plane_index as i32;
                facet.borderNoAdjust[border_index] = false;
                facet.borderInward[border_index] = flipped;
                facet.numBorders += 1;
            }
        }
    }
    trace_str("add the edge bevels");
    //
    // add the edge bevels
    //
    // test the non-axial plane edges
    for j in 0..w.len() {
        let k = (j + 1) % w.len();
        let mut vec = w[j] - w[k];
        // if it's a degenerate edge
        if VectorNormalize_mut(&mut vec) < 0.5 {
            continue;
        }
        vec = CM_SnapVector_func(vec);
        let any_axial = (0..3)
            .map(|k| vec[k])
            .any(|value| value == -1.0 || value == 1.0);
        if any_axial {
            continue; // only test non-axial edges
        }

        // try the six possible slanted axials from this edge
        for axis in 0..3 {
            for &dir in [-1.0, 1.0].into_iter() {
                // construct a plane
                let mut vec2 = vec3_t::ORIGIN;
                vec2[axis] = dir;
                let mut plane_normal = CrossProduct(vec, vec2);
                if VectorNormalize_mut(&mut plane_normal) < 0.5 {
                    continue;
                }
                let plane = plane_t {
                    normal: plane_normal,
                    dist: w[j].dot(plane_normal),
                };

                // if all the points of the facet winding are
                // behind this plane, it is a proper edge bevel
                if w.iter().any(|&p| plane.distance_to(p) > 0.1) {
                    continue;
                }

                //if it's the surface plane
                if CM_PlaneEqual(planes[facet.surfacePlane as usize].plane, plane).is_equal {
                    continue;
                }
                // see if the plane is allready present
                if !facet.borderPlanes[..facet.numBorders as usize]
                    .iter()
                    .any(|&plane_num| {
                        CM_PlaneEqual(planes[plane_num as usize].plane, plane).is_equal
                    })
                {
                    if facet.numBorders > 4 + 6 + 16 {
                        warn!("ERROR: too many bevels\n");
                    }
                    let (plane_index, flipped) = CM_FindPlane2(planes, plane);
                    facet.borderPlanes[facet.numBorders as usize] = plane_index as i32;

                    for k in 0..facet.numBorders as usize {
                        if facet.borderPlanes[facet.numBorders as usize] == facet.borderPlanes[k] {
                            warn!("WARNING: bevel plane already used\n");
                        }
                    }

                    facet.borderNoAdjust[facet.numBorders as usize] = false;
                    facet.borderInward[facet.numBorders as usize] = flipped;
                    //
                    let mut w2 = CopyWinding(&w);
                    let mut newplane =
                        planes[facet.borderPlanes[facet.numBorders as usize] as usize].plane;
                    if !facet.borderInward[facet.numBorders as usize] {
                        newplane = newplane.flip();
                    }
                    ChopWindingInPlace(&mut w2, newplane.normal, newplane.dist, 0.1);
                    if w2.is_empty() {
                        warn!("WARNING: CM_AddFacetBevels... invalid bevel\n");
                        continue;
                    }
                    //
                    facet.numBorders += 1;
                    //already got a bevel
                    //                  break;
                }
            }
        }
    }
    drop(w);

    #[cfg(not(BSPC))]
    {
        //add opposite plane
        facet.borderPlanes[facet.numBorders as usize] = facet.surfacePlane;
        facet.borderNoAdjust[facet.numBorders as usize] = false;
        facet.borderInward[facet.numBorders as usize] = true;
        facet.numBorders += 1;
    }

    trace_str("end");
}

pub type edgeName_t = usize;
pub const EN_TOP: edgeName_t = 0;
pub const EN_RIGHT: edgeName_t = 1;
pub const EN_BOTTOM: edgeName_t = 2;
pub const EN_LEFT: edgeName_t = 3;

fn CM_PatchCollideFromGrid(grid: &cGrid_t) -> Result<patchCollide_t, Error> {
    let mut gridPlanes: GridPlanes = [[[0; 2]; MAX_GRID_SIZE]; MAX_GRID_SIZE];
    let mut facets: Vec<facet_t> = Vec::new();
    let mut planes: Vec<patchPlane_t> = Vec::new();

    trace_str("CM_PatchCollideFromGrid");

    // find the planes for each triangle of the grid
    trace_str("find the planes for each triangle of the grid");
    for i in 0..grid.width as usize - 1 {
        for j in 0..grid.height as usize - 1 {
            let p1 = grid.points[i][j];
            let p2 = grid.points[i + 1][j];
            let p3 = grid.points[i + 1][j + 1];
            gridPlanes[i][j][0] = CM_FindPlane(&mut planes, p1, p2, p3)?;

            let p1 = grid.points[i + 1][j + 1];
            let p2 = grid.points[i][j + 1];
            let p3 = grid.points[i][j];
            gridPlanes[i][j][1] = CM_FindPlane(&mut planes, p1, p2, p3)?;
        }
    }

    // create the borders for each facet
    trace_str("create the borders for each facet");
    for i in 0..grid.width as usize - 1 {
        for j in 0..grid.height as usize - 1 {
            trace_str("iter: i,j");
            trace_i32(i as i32);
            trace_i32(j as i32);

            let mut borders: [i32; 4] = [0; 4];
            let mut noAdjust: [bool; 4] = [false; 4];

            borders[EN_TOP] = -1;
            if j > 0 {
                borders[EN_TOP] = gridPlanes[i][j - 1][1];
            } else if grid.wrapHeight {
                borders[EN_TOP] = gridPlanes[i][grid.height - 2][1];
            }
            noAdjust[EN_TOP] = borders[EN_TOP] == gridPlanes[i][j][0];
            if borders[EN_TOP] == -1 || noAdjust[EN_TOP] {
                borders[EN_TOP] = CM_EdgePlaneNum(&mut planes, grid, &gridPlanes, i, j, 0)?;
            }
            trace_str("borders[EN_TOP]:");
            trace_i32(borders[EN_TOP]);

            borders[EN_BOTTOM] = -1;
            if j < grid.height - 2 {
                borders[EN_BOTTOM] = gridPlanes[i][j + 1][0];
            } else if grid.wrapHeight {
                borders[EN_BOTTOM] = gridPlanes[i][0][0];
            }
            noAdjust[EN_BOTTOM] = borders[EN_BOTTOM] == gridPlanes[i][j][1];
            if borders[EN_BOTTOM] == -1 || noAdjust[EN_BOTTOM] {
                borders[EN_BOTTOM] = CM_EdgePlaneNum(&mut planes, grid, &gridPlanes, i, j, 2)?;
            }
            trace_str("borders[EN_BOTTOM]:");
            trace_i32(borders[EN_BOTTOM]);

            // EN_LEFT
            borders[EN_LEFT] = -1;
            if i > 0 {
                borders[EN_LEFT] = gridPlanes[i - 1][j][0];
            } else if grid.wrapWidth {
                borders[EN_LEFT] = gridPlanes[grid.width - 2][j][0];
            }
            noAdjust[EN_LEFT] = borders[EN_LEFT] == gridPlanes[i][j][1];
            if borders[EN_LEFT] == -1 || noAdjust[EN_LEFT] {
                borders[EN_LEFT] = CM_EdgePlaneNum(&mut planes, grid, &gridPlanes, i, j, 3)?;
            }
            trace_str("borders[EN_LEFT]:");
            trace_i32(borders[EN_LEFT]);

            // EN_RIGHT
            borders[EN_RIGHT] = -1;
            if i < grid.width - 2 {
                borders[EN_RIGHT] = gridPlanes[i + 1][j][1];
            } else if grid.wrapWidth {
                borders[EN_RIGHT] = gridPlanes[0][j][1];
            }
            noAdjust[EN_RIGHT] = borders[EN_RIGHT] == gridPlanes[i][j][0];
            if borders[EN_RIGHT] == -1 || noAdjust[EN_RIGHT] {
                borders[EN_RIGHT] = CM_EdgePlaneNum(&mut planes, grid, &gridPlanes, i, j, 1)?;
            }
            trace_str("borders[EN_RIGHT]:");
            trace_i32(borders[EN_RIGHT]);

            if gridPlanes[i][j][0] == gridPlanes[i][j][1] {
                if gridPlanes[i][j][0] == -1 {
                    trace_str("degenerate");
                    continue; // degenerate
                }
                trace_str("adding facet:");
                let mut facet = facet_t::default();
                facet.surfacePlane = gridPlanes[i][j][0];
                facet.numBorders = 4;
                facet.borderPlanes[0] = borders[EN_TOP];
                facet.borderNoAdjust[0] = noAdjust[EN_TOP];
                facet.borderPlanes[1] = borders[EN_RIGHT];
                facet.borderNoAdjust[1] = noAdjust[EN_RIGHT];
                facet.borderPlanes[2] = borders[EN_BOTTOM];
                facet.borderNoAdjust[2] = noAdjust[EN_BOTTOM];
                facet.borderPlanes[3] = borders[EN_LEFT];
                facet.borderNoAdjust[3] = noAdjust[EN_LEFT];
                CM_SetBorderInward(&planes, &mut facet, grid, &gridPlanes, i, j, -1);
                if CM_ValidateFacet(&planes, &facet) {
                    CM_AddFacetBevels(&mut planes, &mut facet);
                    facets.push(facet);
                }
            } else {
                // two seperate triangles
                trace_str("two separate triangles");
                {
                    let mut facet = facet_t::default();
                    facet.surfacePlane = gridPlanes[i][j][0];
                    facet.numBorders = 3;
                    facet.borderPlanes[0] = borders[EN_TOP];
                    facet.borderNoAdjust[0] = noAdjust[EN_TOP];
                    facet.borderPlanes[1] = borders[EN_RIGHT];
                    facet.borderNoAdjust[1] = noAdjust[EN_RIGHT];
                    facet.borderPlanes[2] = gridPlanes[i][j][1];
                    if facet.borderPlanes[2] == -1 {
                        facet.borderPlanes[2] = borders[EN_BOTTOM];
                        if facet.borderPlanes[2] == -1 {
                            facet.borderPlanes[2] =
                                CM_EdgePlaneNum(&mut planes, grid, &gridPlanes, i, j, 4)?;
                        }
                    }
                    CM_SetBorderInward(&planes, &mut facet, grid, &gridPlanes, i, j, 0);
                    if CM_ValidateFacet(&planes, &facet) {
                        CM_AddFacetBevels(&mut planes, &mut facet);
                        facets.push(facet);
                    }
                }

                {
                    let mut facet = facet_t::default();

                    facet.surfacePlane = gridPlanes[i][j][1];
                    facet.numBorders = 3;
                    facet.borderPlanes[0] = borders[EN_BOTTOM];
                    facet.borderNoAdjust[0] = noAdjust[EN_BOTTOM];
                    facet.borderPlanes[1] = borders[EN_LEFT];
                    facet.borderNoAdjust[1] = noAdjust[EN_LEFT];
                    facet.borderPlanes[2] = gridPlanes[i][j][0];
                    if facet.borderPlanes[2] == -1 {
                        facet.borderPlanes[2] = borders[EN_TOP];
                        if facet.borderPlanes[2] == -1 {
                            facet.borderPlanes[2] =
                                CM_EdgePlaneNum(&mut planes, grid, &gridPlanes, i, j, 5)?;
                        }
                    }
                    CM_SetBorderInward(&planes, &mut facet, grid, &gridPlanes, i, j, 1);
                    if CM_ValidateFacet(&planes, &facet) {
                        CM_AddFacetBevels(&mut planes, &mut facet);
                        facets.push(facet);
                    }
                }
            }
        }
    }

    trace_str("facets:");
    trace_i32(facets.len() as i32);
    for f in facets.iter() {
        trace_str("facet:");
        trace_i32(f.surfacePlane);
        trace_i32(f.numBorders);
        trace_str("...borderPlanes");
        for &b in f.borderPlanes[..f.numBorders as usize].iter() {
            trace_i32(b);
        }
        trace_str("...borderInward");
        for &inward in f.borderInward[..f.numBorders as usize].iter() {
            trace_i32(inward as i32);
        }
        trace_str("...borderNoAdjust");
        for &na in f.borderNoAdjust[..f.numBorders as usize].iter() {
            trace_i32(na as i32);
        }
    }
    trace_str("planes:");
    trace_i32(planes.len() as i32);
    for p in planes.iter() {
        trace_str("plane:");
        trace_vec3(p.plane.normal);
        trace_f32(p.plane.dist);
    }
    trace_str(".");

    // copy the results out
    Ok(patchCollide_t {
        facets,
        planes,
        bounds: ClearBounds(),
    })
}

fn trace_grid(grid: &cGrid_t) {
    trace_str("grid");
    trace_i32(grid.wrapWidth as i32);
    trace_i32(grid.wrapHeight as i32);
    for i in 0..grid.width {
        for j in 0..grid.height {
            trace_vec3(grid.points[i][j]);
        }
    }
    trace_str(".");
}
/*
===================
CM_GeneratePatchCollide

Creates an internal structure that will be used to perform
collision detection with a patch mesh.

Points is packed as concatenated rows.
===================
*/
pub fn CM_GeneratePatchCollide(
    width: usize,
    height: usize,
    points: &[vec3_t],
) -> Result<patchCollide_t, Error> {
    assert!(width > 2);
    assert!(height > 2);
    assert!(!points.is_empty());
    assert!(
        is_odd(width) && is_odd(height),
        "CM_GeneratePatchFacets: even sizes are invalid for quadratic meshes"
    );
    assert!(width <= MAX_GRID_SIZE);
    assert!(height <= MAX_GRID_SIZE);

    let mut grid: cGrid_t = cGrid_t::empty();

    // build a grid
    grid.width = width;
    grid.height = height;
    grid.wrapWidth = false;
    grid.wrapHeight = false;
    for i in 0..width {
        for j in 0..height {
            grid.points[i][j] = points[j * width + i];
        }
    }

    trace_grid(&grid);

    // subdivide the grid
    CM_SetGridWrapWidth(&mut grid);
    trace_grid(&grid);
    CM_SubdivideGridColumns(&mut grid);
    trace_grid(&grid);
    CM_RemoveDegenerateColumns(&mut grid);
    trace_grid(&grid);

    CM_TransposeGrid(&mut grid);
    trace_grid(&grid);

    CM_SetGridWrapWidth(&mut grid);
    trace_grid(&grid);
    CM_SubdivideGridColumns(&mut grid);
    trace_grid(&grid);
    CM_RemoveDegenerateColumns(&mut grid);
    trace_grid(&grid);

    // we now have a grid of points exactly on the curve
    // the aproximate surface defined by these points will be
    // collided against
    let mut bounds = ClearBounds();
    for i in 0..grid.width as usize {
        for j in 0..grid.height as usize {
            AddPointToBounds(grid.points[i][j], &mut bounds);
        }
    }

    c_totalPatchBlocks.add_usize((grid.width - 1) * (grid.height - 1));

    // generate a bsp tree for the surface
    let mut pf = CM_PatchCollideFromGrid(&grid)?;

    // expand by one unit for epsilon purposes
    bounds.mins[0] -= 1.0;
    bounds.mins[1] -= 1.0;
    bounds.mins[2] -= 1.0;

    bounds.maxs[0] += 1.0;
    bounds.maxs[1] += 1.0;
    bounds.maxs[2] += 1.0;

    pf.bounds = bounds;
    Ok(pf)
}

/*
================================================================================

TRACE TESTING

================================================================================
*/

/// Special case for point traces because the patch collide "brushes" have no volume
fn CM_TracePointThroughPatchCollide(tw: &mut traceWork_t, pc: &patchCollide_t) {
    let mut frontFacing = [false; MAX_PATCH_PLANES];
    let mut intersection = [0.0f32; MAX_PATCH_PLANES];

    /*
    #ifndef BSPC
        static cvar_t *cv;
    #endif //BSPC

    #ifndef BSPC
        if ( !cm_playerCurveClip->integer || !tw.isPoint ) {
            return;
        }
    #endif
    */

    // determine the trace's relationship to all planes
    for (i, planes) in pc.planes.iter().enumerate() {
        let offset = tw.offsets[planes.signbits as usize].dot(planes.plane.normal);
        let d1 = planes.plane.distance_to(tw.start) + offset;
        let d2 = planes.plane.distance_to(tw.end) + offset;
        frontFacing[i] = d1 > 0.0;
        intersection[i] = if d1 == d2 {
            99999.0
        } else {
            let ix = d1 / (d1 - d2);
            if ix > 0.0 {
                ix
            } else {
                99999.0
            }
        };
    }

    // see if any of the surface planes are intersected
    for facet in pc.facets.iter() {
        if !frontFacing[facet.surfacePlane as usize] {
            continue;
        }
        let intersect = intersection[facet.surfacePlane as usize];
        if intersect < 0.0 {
            continue; // surface is behind the starting point
        }
        if intersect > tw.trace.fraction {
            continue; // already hit something closer
        }
        let mut hit_facet = true;
        for j in 0..facet.numBorders as usize {
            let k = facet.borderPlanes[j] as usize;
            if frontFacing[k] ^ facet.borderInward[j] {
                if intersection[k] > intersect {
                    hit_facet = false;
                    break;
                }
            } else {
                if intersection[k] < intersect {
                    hit_facet = false;
                    break;
                }
            }
        }
        if hit_facet {
            // we hit this facet
            /*
            #ifndef BSPC
                        if (!cv) {
                            cv = Cvar_Get( "r_debugSurfaceUpdate", "1", 0 );
                        }
                        if (cv->integer) {
                            debugPatchCollide = pc;
                            debugFacet = facet;
                        }
            #endif //BSPC
            */
            let planes = &pc.planes[facet.surfacePlane as usize];

            // calculate intersection with a slight pushoff
            let offset = tw.offsets[planes.signbits as usize].dot(planes.plane.normal);
            let d1 = planes.plane.distance_to(tw.start) + offset;
            let d2 = planes.plane.distance_to(tw.end) + offset;
            tw.trace.fraction = (d1 - SURFACE_CLIP_EPSILON) / (d1 - d2);

            if tw.trace.fraction < 0.0 {
                tw.trace.fraction = 0.0;
            }

            tw.trace.plane.normal = planes.plane.normal;
            tw.trace.plane.dist = planes.plane.dist;
        }
    }
}

struct CheckFacetPlaneOutput {
    pub result: bool,
    pub hit: bool,
}
fn CM_CheckFacetPlane(
    plane: &plane_t,
    start: vec3_t,
    end: vec3_t,
    enterFrac: &mut f32,
    leaveFrac: &mut f32,
) -> CheckFacetPlaneOutput {
    type Output = CheckFacetPlaneOutput;

    let mut hit = false;

    let d1 = plane.distance_to(start);
    let d2 = plane.distance_to(end);

    // if completely in front of face, no intersection with the entire facet
    if d1 > 0.0 && (d2 >= SURFACE_CLIP_EPSILON || d2 >= d1) {
        return Output { result: false, hit };
    }

    // if it doesn't cross the plane, the plane isn't relevent
    if d1 <= 0.0 && d2 <= 0.0 {
        return Output { result: true, hit };
    }

    // crosses face
    if d1 > d2 {
        // enter
        let f = fmax(0.0, (d1 - SURFACE_CLIP_EPSILON) / (d1 - d2));
        // always favor previous plane hits and thus also the surface plane hit
        if f > *enterFrac {
            *enterFrac = f;
            hit = true;
        }
    } else {
        // leave
        let f = fmin(1.0, (d1 + SURFACE_CLIP_EPSILON) / (d1 - d2));
        if f < *leaveFrac {
            *leaveFrac = f;
        }
    }
    return Output { result: true, hit };
}

pub fn CM_TraceThroughPatchCollide(tw: &mut traceWork_t, pc: &patchCollide_t) {
    if tw.isPoint {
        CM_TracePointThroughPatchCollide(tw, pc);
        return;
    }

    let mut bestplane = plane_t::default();
    /*
    #ifndef BSPC
        static cvar_t *cv;
    #endif //BSPC
    */

    for facet in pc.facets.iter() {
        let mut enterFrac = -1.0;
        let mut leaveFrac = 1.0;
        let mut hitnum: i32 = -1;
        //
        {
            let planes = &pc.planes[facet.surfacePlane as usize];
            let mut plane = planes.plane;
            let startp;
            let endp;
            let offset;
            if tw.sphere.use_.into() {
                // adjust the plane distance apropriately for radius
                plane.dist += tw.sphere.radius;

                // find the closest point on the capsule to the plane
                let t = plane.normal.dot(tw.sphere.offset);
                if t > 0.0 {
                    startp = tw.start - tw.sphere.offset;
                    endp = tw.end - tw.sphere.offset;
                } else {
                    startp = tw.start + tw.sphere.offset;
                    endp = tw.end + tw.sphere.offset;
                }
            } else {
                offset = tw.offsets[planes.signbits as usize].dot(plane.normal);
                plane.dist -= offset;
                startp = tw.start;
                endp = tw.end;
            }

            let check = CM_CheckFacetPlane(&plane, startp, endp, &mut enterFrac, &mut leaveFrac);
            if !check.result {
                continue;
            }
            if check.hit {
                bestplane = plane;
            }
        }

        let mut early_break = false;
        for j in 0..facet.numBorders as usize {
            let planes = &pc.planes[facet.borderPlanes[j] as usize];
            let mut plane = if facet.borderInward[j] {
                planes.plane.flip()
            } else {
                planes.plane
            };
            let startp;
            let endp;
            let offset;
            if tw.sphere.use_.into() {
                // adjust the plane distance apropriately for radius
                plane.dist += tw.sphere.radius;
                // find the closest point on the capsule to the plane
                let t = plane.normal.dot(tw.sphere.offset);
                if t > 0.0 {
                    startp = tw.start - tw.sphere.offset;
                    endp = tw.end - tw.sphere.offset;
                } else {
                    startp = tw.start + tw.sphere.offset;
                    endp = tw.end + tw.sphere.offset;
                }
            } else {
                // NOTE: this works even though the plane might be flipped because the bbox is centered
                offset = tw.offsets[planes.signbits as usize].dot(plane.normal);
                plane.dist += offset.abs();
                startp = tw.start;
                endp = tw.end;
            }

            let check = CM_CheckFacetPlane(&plane, startp, endp, &mut enterFrac, &mut leaveFrac);
            if !check.result {
                early_break = true;
                break;
            }
            if check.hit {
                hitnum = j as i32;
                bestplane = plane;
            }
        }
        if early_break {
            continue;
        }
        //never clip against the back side
        if hitnum == facet.numBorders - 1 {
            continue;
        }

        if enterFrac < leaveFrac && enterFrac >= 0.0 {
            if enterFrac < tw.trace.fraction {
                if enterFrac < 0.0 {
                    enterFrac = 0.0;
                }
                /*
                #ifndef BSPC
                                if (!cv) {
                                    cv = Cvar_Get( "r_debugSurfaceUpdate", "1", 0 );
                                }
                                if (cv && cv->integer) {
                                    debugPatchCollide = pc;
                                    debugFacet = facet;
                                }
                #endif //BSPC
                */

                tw.trace.fraction = enterFrac;
                tw.trace.plane.normal = bestplane.normal;
                tw.trace.plane.dist = bestplane.dist;
            }
        }
    }
}

/*
=======================================================================

POSITION TEST

=======================================================================
*/

pub fn CM_PositionTestInPatchCollide(tw: &mut traceWork_t, pc: &patchCollide_t) -> bool {
    if tw.isPoint {
        return false;
    }
    for facet in pc.facets.iter() {
        {
            let planes = &pc.planes[facet.surfacePlane as usize];
            let startp;
            let mut plane = planes.plane;
            if tw.sphere.use_.into() {
                // adjust the plane distance apropriately for radius
                plane.dist += tw.sphere.radius;

                // find the closest point on the capsule to the plane
                let t = plane.normal.dot(tw.sphere.offset);
                startp = if t > 0.0 {
                    tw.start - tw.sphere.offset
                } else {
                    tw.start + tw.sphere.offset
                };
            } else {
                plane.dist -= tw.offsets[planes.signbits as usize].dot(plane.normal);
                startp = tw.start;
            }

            if plane.distance_to(startp) > 0.0 {
                continue;
            }
        }

        let mut early_break = false;
        for j in 0..facet.numBorders as usize {
            let planes = &pc.planes[facet.borderPlanes[j] as usize];
            let mut plane = if facet.borderInward[j] {
                planes.plane.flip()
            } else {
                planes.plane
            };
            let startp;
            if tw.sphere.use_.into() {
                // adjust the plane distance apropriately for radius
                plane.dist += tw.sphere.radius;

                // find the closest point on the capsule to the plane
                let t = plane.normal.dot(tw.sphere.offset);
                startp = if t > 0.0 {
                    tw.start - tw.sphere.offset
                } else {
                    tw.start + tw.sphere.offset
                };
            } else {
                // NOTE: this works even though the plane might be flipped because the bbox is centered
                plane.dist += tw.offsets[planes.signbits as usize].dot(plane.normal).abs();
                startp = tw.start;
            }

            if plane.distance_to(startp) > 0.0 {
                early_break = true;
                break;
            }
        }
        if early_break {
            continue;
        }
        // inside this patch facet
        return true;
    }
    return false;
}

/*

/*
=======================================================================

DEBUGGING

=======================================================================
*/


/*
==================
CM_DrawDebugSurface

Called from the renderer
==================
*/
#ifndef BSPC
void BotDrawDebugPolygons(void (*drawPoly)(int color, int numPoints, float *points), int value);
#endif

void CM_DrawDebugSurface( void (*drawPoly)(int color, int numPoints, float *points) ) {
    static cvar_t   *cv;
#ifndef BSPC
    static cvar_t   *cv2;
#endif
    const patchCollide_t    *pc;
    facet_t         *facet;
    winding_t       *w;
    int             i, j, k, n;
    int             curplanenum, planenum, curinward, inward;
    float           plane[4];
    vec3_t mins = {-15, -15, -28}, maxs = {15, 15, 28};
    //vec3_t mins = {0, 0, 0}, maxs = {0, 0, 0};
    vec3_t v1, v2;

#ifndef BSPC
    if ( !cv2 )
    {
        cv2 = Cvar_Get( "r_debugSurface", "0", 0 );
    }

    if (cv2->integer != 1)
    {
        BotDrawDebugPolygons(drawPoly, cv2->integer);
        return;
    }
#endif

    if ( !debugPatchCollide ) {
        return;
    }

#ifndef BSPC
    if ( !cv ) {
        cv = Cvar_Get( "cm_debugSize", "2", 0 );
    }
#endif
    pc = debugPatchCollide;

    for ( i = 0, facet = pc->facets ; i < pc->numFacets ; i++, facet++ ) {

        for ( k = 0 ; k < facet.numBorders + 1; k++ ) {
            //
            if (k < facet.numBorders) {
                planenum = facet.borderPlanes[k];
                inward = facet.borderInward[k];
            }
            else {
                planenum = facet.surfacePlane;
                inward = false;
                //continue;
            }

            Vector4Copy( pc->planes[ planenum ].plane, plane );

            //planenum = facet.surfacePlane;
            if ( inward ) {
                VectorSubtract( vec3_origin, plane, plane );
                plane[3] = -plane[3];
            }

            plane[3] += cv->value;
            // *
            for (n = 0; n < 3; n++)
            {
                if (plane[n] > 0) v1[n] = maxs[n];
                else v1[n] = mins[n];
            } //end for
            VectorNegate(plane, v2);
            plane[3] += fabs(DotProduct(v1, v2));
            // * /

            w = BaseWindingForPlane( plane,  plane[3] );
            for ( j = 0 ; j < facet.numBorders + 1 && w; j++ ) {
                //
                if (j < facet.numBorders) {
                    curplanenum = facet.borderPlanes[j];
                    curinward = facet.borderInward[j];
                }
                else {
                    curplanenum = facet.surfacePlane;
                    curinward = false;
                    //continue;
                }
                //
                if (curplanenum == planenum) continue;

                Vector4Copy( pc->planes[ curplanenum ].plane, plane );
                if ( !curinward ) {
                    VectorSubtract( vec3_origin, plane, plane );
                    plane[3] = -plane[3];
                }
        //          if ( !facet.borderNoAdjust[j] ) {
                    plane[3] -= cv->value;
        //          }
                for (n = 0; n < 3; n++)
                {
                    if (plane[n] > 0) v1[n] = maxs[n];
                    else v1[n] = mins[n];
                } //end for
                VectorNegate(plane, v2);
                plane[3] -= fabs(DotProduct(v1, v2));

                ChopWindingInPlace( &w, plane, plane[3], 0.1f );
            }
            if ( w ) {
                if ( facet == debugFacet ) {
                    drawPoly( 4, w.len(), w[0] );
                    //Com_Printf("blue facet has %d border planes\n", facet.numBorders);
                } else {
                    drawPoly( 1, w.len(), w[0] );
                }
                FreeWinding( w );
            }
            else
                Com_Printf("winding chopped away by border planes\n");
        }
    }

    // draw the debug block
    {
        vec3_t          v[3];

        VectorCopy( debugBlockPoints[0], v[0] );
        VectorCopy( debugBlockPoints[1], v[1] );
        VectorCopy( debugBlockPoints[2], v[2] );
        drawPoly( 2, 3, v[0] );

        VectorCopy( debugBlockPoints[2], v[0] );
        VectorCopy( debugBlockPoints[3], v[1] );
        VectorCopy( debugBlockPoints[0], v[2] );
        drawPoly( 2, 3, v[0] );
    }

#if 0
    vec3_t          v[4];

    v[0][0] = pc->bounds[1][0];
    v[0][1] = pc->bounds[1][1];
    v[0][2] = pc->bounds[1][2];

    v[1][0] = pc->bounds[1][0];
    v[1][1] = pc->bounds[0][1];
    v[1][2] = pc->bounds[1][2];

    v[2][0] = pc->bounds[0][0];
    v[2][1] = pc->bounds[0][1];
    v[2][2] = pc->bounds[1][2];

    v[3][0] = pc->bounds[0][0];
    v[3][1] = pc->bounds[1][1];
    v[3][2] = pc->bounds[1][2];

    drawPoly( 4, v[0] );
#endif
}
*/
