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

use crate::prelude::*;
use crate::qcommon::cm_patch::patchCollide_t;
// use super::cm_polylib::*;

pub const MAX_SUBMODELS: usize = 256;
pub const BOX_MODEL_HANDLE: usize = 255;
pub const CAPSULE_MODEL_HANDLE: usize = 254;

#[derive(Clone)]
pub struct cNode_t {
    pub plane: *mut cplane_t,
    pub children: [i32; 2], // negative numbers are leafs
}

#[derive(Clone, Debug)]
pub struct cLeaf_t {
    pub cluster: i32,
    pub area: i32,

    pub firstLeafBrush: i32,
    pub numLeafBrushes: i32,

    pub firstLeafSurface: i32,
    pub numLeafSurfaces: i32,
}

pub type cmodel_s = cmodel_t;

#[derive(Clone, Debug)]
pub struct cmodel_t {
    pub mins: vec3_t,
    pub maxs: vec3_t,
    pub leaf: cLeaf_t, // submodels don't reference the main tree
}

#[derive(Clone, Debug)]
pub struct cbrushside_t {
    pub plane: *mut cplane_t,
    pub surfaceFlags: i32,
    pub shaderNum: i32,
}

#[derive(Clone, Debug)]
pub struct cbrush_t {
    pub shaderNum: i32, // the shader that determined the contents
    pub contents: i32,
    pub bounds: vec3_bounds,
    pub numsides: i32,
    pub sides: *mut cbrushside_t,
    pub checkcount: i32, // to avoid repeated testings
}

#[derive(Clone, Debug)]
pub struct cPatch_t {
    pub checkcount: i32, // to avoid repeated testings
    pub surfaceFlags: i32,
    pub contents: i32,
    pub pc: *mut patchCollide_t,
}

#[derive(Clone, Debug)]
pub struct cArea_t {
    pub floodnum: i32,
    pub floodvalid: i32,
}

todo_type! {dshader_t}

pub struct clipMap_t {
    pub name: String,

    pub shaders: Vec<dshader_t>,
    pub brushsides: Vec<cbrushside_t>,
    pub planes: Vec<cplane_t>,
    pub nodes: Vec<cNode_t>,
    pub leafs: Vec<cLeaf_t>,
    pub leafbrushes: Vec<i32>,
    pub leafsurfaces: Vec<i32>,
    pub cmodels: Vec<cmodel_t>, // len was numSubModels, not numCModels
    pub brushes: Vec<cbrush_t>,

    pub numClusters: i32,
    pub clusterBytes: i32,
    pub visibility: *mut u8,
    pub vised: bool, // if false, visibility is just a single cluster of ffs

    pub entityString: String,

    pub areas: Vec<cArea_t>,
    pub areaPortals: *mut i32, // [ numAreas*numAreas ] reference counts

    pub surfaces: Vec<*mut cPatch_t>, // non-patches will be NULL

    pub floodvalid: i32,
    pub checkcount: i32, // incremented on each trace
}

// keep 1/8 unit away to keep the position valid before network snapping
// and to avoid various numeric issues
pub const SURFACE_CLIP_EPSILON: f32 = 0.125;

/*
extern  clipMap_t   cm;
extern  i32         c_pointcontents;
extern  i32         c_traces, c_brush_traces, c_patch_traces;
extern  cvar_t      *cm_noAreas;
extern  cvar_t      *cm_noCurves;
extern  cvar_t      *cm_playerCurveClip;
*/

// cm_test.c

// Used for oriented capsule collision detection
pub struct sphere_t {
    pub use_: bool,
    pub radius: f32,
    pub halfheight: f32,
    pub offset: vec3_t,
}

pub struct traceWork_t {
    pub start: vec3_t,
    pub end: vec3_t,
    pub size: [vec3_t; 2],    // size of the box being swept through the model
    pub offsets: [vec3_t; 8], // [signbits][x] = either size[0][x] or size[1][x]
    pub maxOffset: f32,       // longest corner length from origin
    pub extents: vec3_t,      // greatest of abs(size[0]) and abs(size[1])
    pub bounds: [vec3_t; 2],  // enclosing box of start and end surrounding by size
    pub modelOrigin: vec3_t,  // origin of the model tracing through
    pub contents: i32,        // ored contents of the model tracing through
    pub isPoint: bool,        // optimized case
    pub trace: trace_t,       // returned from trace call
    pub sphere: sphere_t,     // sphere for oriendted capsule collision
}

pub type leafList_s = leafList_t;
pub struct leafList_t {
    pub count: i32,
    pub maxcount: i32,
    pub overflowed: bool,
    pub list: *mut i32,
    pub bounds: vec3_bounds,
    pub lastLeaf: i32, // for overflows where each leaf can't be stored individually
    pub storeLeafs: fn(ll: *mut leafList_s, nodenum: i32),
}
