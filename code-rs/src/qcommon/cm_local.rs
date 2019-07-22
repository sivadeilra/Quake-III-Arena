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
use crate::qfiles::*;
use std::ops::Range;

pub const MAX_SUBMODELS: usize = 256;
pub const BOX_MODEL_HANDLE: clipHandle_t = 255;
pub const CAPSULE_MODEL_HANDLE: clipHandle_t = 254;

pub type clipHandle_t = i32;

#[derive(Clone)]
pub struct cNode_t {
    // these are pointers into clipMap_t::planes[]
    // pub plane: *mut cplane_t,

    // index into clipMap_t::planes[]
    pub plane_index: i32,

    pub children: [i32; 2], // negative numbers are leafs
}

impl cNode_t {
    pub fn plane<'a>(&self, cm: &'a clipMap_t) -> &'a cplane_t {
        &cm.planes[self.plane_index as usize]
    }
}

#[derive(Clone, Debug, Default)]
pub struct cLeaf_t {
    pub cluster: i32,
    pub area: i32,

    pub firstLeafBrush: i32,
    pub numLeafBrushes: i32,

    pub firstLeafSurface: i32,
    pub numLeafSurfaces: i32,
}

impl cLeaf_t {
    pub fn leaf_brushes_range(&self) -> Range<usize> {
        let start = self.firstLeafBrush as usize;
        let count = self.numLeafBrushes as usize;
        start..start + count
    }
    pub fn leaf_brushes_index<'a>(&self, cm: &'a clipMap_t) -> &'a [i32] {
        &cm.leafbrushes[self.leaf_brushes_range()]
    }

    pub fn leaf_brushes<'a>(&self, cm: &'a clipMap_t) -> impl Iterator<Item = &'a cbrush_t> {
        self.leaf_brushes_index(cm)
            .iter()
            .map(move |&i| &cm.brushes[i as usize])
    }

    pub fn leaf_surfaces_range(&self) -> Range<usize> {
        self.firstLeafSurface as usize
            ..self.firstLeafSurface as usize + self.numLeafSurfaces as usize
    }

    pub fn leaf_surfaces<'a>(&self, cm: &'a clipMap_t) -> &'a [Option<Box<cPatch_t>>] {
        &cm.surfaces[self.leaf_surfaces_range()]
    }
    pub fn leaf_surfaces_mut<'a>(
        &self,
        cm_surfaces: &'a mut [Option<Box<cPatch_t>>],
    ) -> &'a mut [Option<Box<cPatch_t>>] {
        &mut cm_surfaces[self.leaf_surfaces_range()]
    }
}

pub type cmodel_s = cmodel_t;

#[derive(Clone, Debug, Default)]
pub struct cmodel_t {
    pub mins: vec3_t,
    pub maxs: vec3_t,
    pub leaf: cLeaf_t, // submodels don't reference the main tree
}

#[derive(Clone, Debug, Default)]
pub struct cbrushside_t {
    // pub plane: *mut cplane_t,
    /// index into clipMap_t::planes
    pub plane_num: i32,

    pub surfaceFlags: i32,
    pub shaderNum: i32,
}

impl cbrushside_t {
    pub fn plane<'a>(&self, cm: &'a clipMap_t) -> &'a cplane_t {
        &cm.planes[self.plane_num as usize]
    }
}

#[derive(Clone, Debug, Default)]
pub struct cbrush_t {
    pub shaderNum: i32, // the shader that determined the contents
    pub contents: i32,
    pub bounds: vec3_bounds,
    pub numsides: i32,

    //
    // pub sides: *mut cbrushside_t,
    /// index into clipMap_t::brushsides
    pub firstSide: i32,

    pub checkcount: i32, // to avoid repeated testings
}

impl cbrush_t {
    pub fn sides<'a>(&self, cm: &'a clipMap_t) -> &'a [cbrushside_t] {
        let start = self.firstSide as usize;
        let count = self.numsides as usize;
        &cm.brushsides[start..start + count]
    }
}

#[derive(Clone, Debug)]
pub struct cPatch_t {
    pub checkcount: i32, // to avoid repeated testings
    pub surfaceFlags: i32,
    pub contents: i32,
    // pub pc: *mut patchCollide_t,
    pub pc: patchCollide_t,
}

#[derive(Clone, Debug, Default)]
pub struct cArea_t {
    pub floodnum: i32,
    pub floodvalid: i32,
}

#[derive(Clone)]
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
    pub visibility: Vec<u8>,
    pub vised: bool, // if false, visibility is just a single cluster of ffs

    pub entityString: String,

    // this gets reused / rewritten many times
    pub areas: Vec<cArea_t>,
    pub areaPortals: Vec<i32>, // [ numAreas*numAreas ] reference counts

    pub surfaces: Vec<Option<Box<cPatch_t>>>, // non-patches will be NULL

    pub floodvalid: i32,
    pub checkcount: i32, // incremented on each trace

    // these were globals
    pub box_model: cmodel_t,
    // pub box_planes: Vec<cplane_t>,
    // pub box_brush: Vec<cbrush_t>,
    /// index into Self::brushes for box_brush
    pub box_brush_index: usize,

    /// range in Self::planes where box_planes is
    pub box_planes_range: Range<usize>,
}

pub struct ClipMapVis {
    pub numClusters: i32,
    pub clusterBytes: i32,
    pub visibility: Vec<u8>,
    pub vised: bool, // if false, visibility is just a single cluster of ffs
}

impl clipMap_t {
    pub fn box_brush(&self) -> &cbrush_t {
        &self.brushes[self.box_brush_index]
    }
    pub fn box_brush_mut(&mut self) -> &mut cbrush_t {
        &mut self.brushes[self.box_brush_index]
    }

    pub fn box_planes(&self) -> &[cplane_t] {
        &self.planes[self.box_planes_range.clone()]
    }
    pub fn box_planes_mut(&mut self) -> &mut [cplane_t] {
        &mut self.planes[self.box_planes_range.clone()]
    }
}

// keep 1/8 unit away to keep the position valid before network snapping
// and to avoid various numeric issues
pub const SURFACE_CLIP_EPSILON: f32 = 0.125;

/*
extern  clipMap_t   cm;
extern  i32         c_pointcontents;
extern  i32         c_traces, c_brush_traces, c_patch_traces;
*/
pub static cm_noAreas: cvar_t = cvar_t::new("cm_noAreas");
pub static cm_noCurves: cvar_t = cvar_t::new("cm_noCurves");
pub static cm_playerCurveClip: cvar_t = cvar_t::new("cm_playerCurveClip");

// cm_test.c

// Used for oriented capsule collision detection
#[derive(Clone, Debug)]
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
    pub count: usize,
    pub overflowed: bool,
    pub bounds: vec3_bounds,
    pub lastLeaf: i32, // for overflows where each leaf can't be stored individually
}
