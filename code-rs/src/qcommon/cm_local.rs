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
use crate::range_len;
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
        range_len(self.firstLeafBrush as usize, self.numLeafBrushes as usize)
    }

    pub fn leaf_brushes<'a>(&self, cm: &'a clipMap_t) -> impl Iterator<Item = &'a cbrush_t> {
        let brushes = &cm.brushes;
        cm.leafbrushes[self.leaf_brushes_range()]
            .iter()
            .map(move |&i| &brushes[i as usize])
    }

    pub fn leaf_surfaces_range(&self) -> Range<usize> {
        range_len(
            self.firstLeafSurface as usize,
            self.numLeafSurfaces as usize,
        )
    }

    pub fn leaf_surfaces<'a>(&self, cm: &'a clipMap_t) -> &'a [Option<Box<cPatch_t>>] {
        &cm.surfaces[self.leaf_surfaces_range()]
    }
}

#[derive(Clone, Debug, Default)]
pub struct cmodel_t {
    pub mins: vec3_t,
    pub maxs: vec3_t,
    pub leaf: cLeaf_t, // submodels don't reference the main tree
}

#[derive(Clone, Debug, Default)]
pub struct cbrushside_t {
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

    /// index into clipMap_t::brushsides
    pub firstSide: i32,
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
    pub surfaceFlags: i32,
    pub contents: i32,
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

    // these were globals
    pub box_model: cmodel_t,
    // pub box_planes: Vec<cplane_t>,
    // pub box_brush: Vec<cbrush_t>,
    /// index into Self::brushes for box_brush
    pub box_brush_index: usize,

    /// range in Self::planes where box_planes is
    pub box_planes_range: Range<usize>,
}

type CheckCount = i32;
const CHECK_COUNT_MAX: CheckCount = std::i32::MAX;

#[derive(Clone)]
pub struct ClipMapClient {
    /// incremented on each trace
    checkcount: CheckCount,

    /// corresponds to clipMap_t::brushes
    brushes_checkcount: Vec<CheckCount>,

    /// corresponds to clipMap_t::surfaces
    patches_checkcount: Vec<CheckCount>,
}

impl ClipMapClient {
    pub fn new(cm: &clipMap_t) -> Self {
        Self {
            checkcount: 1,
            brushes_checkcount: vec![0; cm.brushes.len()],
            patches_checkcount: vec![0; cm.surfaces.len()],
        }
    }

    /// Increments the generation counter, which invalidates all previous checkcounts for brushes
    /// and patches. This is effectively clearing a cache.
    pub fn next_checkcount(&mut self) {
        use crate::fill;
        if self.checkcount == CHECK_COUNT_MAX {
            // Reset all checkcounts. Without this, it would be possible to accidentally use the
            // wrong results.
            fill(&mut self.brushes_checkcount, 0);
            fill(&mut self.patches_checkcount, 0);
            self.checkcount = 1;
        } else {
            self.checkcount += 1;
        }
    }

    /// Returns true if the brush has ALREADY been checked.
    /// Returns false if the brush NEEDS to be checked. In this case, the checkcount for the brush
    /// is also updated, so that the next check_brush() call will return true.
    pub fn check_brush(&mut self, brush_num: usize) -> bool {
        if self.brushes_checkcount[brush_num] == self.checkcount {
            true
        } else {
            self.brushes_checkcount[brush_num] = self.checkcount;
            false
        }
    }

    pub fn check_patch(&mut self, patch_num: usize) -> bool {
        if self.patches_checkcount[patch_num] == self.checkcount {
            true
        } else {
            self.patches_checkcount[patch_num] = self.checkcount;
            false
        }
    }
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
#[derive(Clone, Debug, Copy, Default)]
#[repr(C)]
pub struct sphere_t {
    pub use_: qboolean,
    pub radius: f32,
    pub halfheight: f32,
    pub offset: vec3_t,
}
