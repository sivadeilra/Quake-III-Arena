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
use crate::zerocopy::*;

macro_rules! define_fixed_c_string {
    (str type $name:ident = [$len:expr];) => {
        #[derive(Copy, Clone)]
        pub struct $name([c_char; $len]);
        impl $name {
            pub const MAX_LEN: usize = $len;
            pub fn as_bytes(&self) -> &[u8] {
                &self.0[..self.len()]
            }
            pub fn as_str(&self) -> Result<&str, core::str::Utf8Error> {
                core::str::from_utf8(self.as_bytes())
            }
            pub fn len(&self) -> usize {
                for (i, &b) in self.0.iter().enumerate() {
                    if b == 0 {
                        return i;
                    }
                }
                $len
            }
        }
        impl Default for $name {
            fn default() -> Self {
                Self([0; $len])
            }
        }
        impl std::fmt::Debug for $name {
            fn fmt(&self, fmt: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                if let Ok(s) = self.as_str() {
                    write!(fmt, "{:?}", s)
                } else {
                    write!(fmt, "#error")
                }
            }
        }
    };
}

define_fixed_c_string! { str type FixedString_QPATH = [MAX_QPATH]; }
define_fixed_c_string! { str type FixedString_16 = [16]; }

//
// qfiles.h: quake file formats
// This file must be identical in the quake and utils directories
//

// surface geometry should not exceed these limits
pub const SHADER_MAX_VERTEXES: i32 = 1000;
pub const SHADER_MAX_INDEXES: i32 = 6 * SHADER_MAX_VERTEXES;

// the maximum size of game relative pathnames
pub const MAX_QPATH: usize = 64;

/*
========================================================================

QVM files

========================================================================
*/

pub const VM_MAGIC: u32 = 0x12721444;

#[repr(C)]
#[derive(Clone, Copy)]
pub struct vmHeader_t {
    pub vmMagic: i32,

    pub instructionCount: i32,

    pub codeOffset: i32,
    pub codeLength: i32,

    pub dataOffset: i32,
    pub dataLength: i32,
    pub litLength: i32, // ( dataLength - litLength ) should be byteswapped on load
    pub bssLength: i32, // zero filled memory appended to datalength
}

/*
========================================================================

PCX files are used for 8 bit images

========================================================================
*/

#[repr(C)]
#[derive(Clone, Copy)]
pub struct pcx_t {
    pub manufacturer: u8,
    pub version: u8,
    pub encoding: u8,
    pub bits_per_pixel: u8,

    pub xmin: u16,
    pub ymin: u16,
    pub xmax: u16,
    pub ymax: u16,
    pub hres: u16,
    pub vres: u16,
    pub palette: [u8; 48],
    pub reserved: u8,
    pub color_planes: u8,
    pub bytes_per_line: u16,
    pub palette_type: u16,
    pub filler: [u8; 58],
    pub data: [u8; 1], // unbounded
}

/*
========================================================================

TGA files are used for 24/32 bit images

========================================================================
*/

#[repr(C)]
#[derive(Clone, Copy)]
pub struct TargaHeader {
    pub id_length: u8,
    pub colormap_type: u8,
    pub image_type: u8,
    pub colormap_index: u16,
    pub colormap_length: u16,
    pub colormap_size: u8,
    pub x_origin: u16,
    pub y_origin: u16,
    pub width: u16,
    pub height: u16,
    pub pixel_size: u8,
    pub attributes: u8,
}

/*
========================================================================

.MD3 triangle model file format

========================================================================
*/

pub const MD3_IDENT: u32 =
    (('3' as u32) << 24) + (('P' as u32) << 16) + (('D' as u32) << 8) + 'I' as u32;
pub const MD3_VERSION: i32 = 15;

// limits
pub const MD3_MAX_LODS: i32 = 3;
pub const MD3_MAX_TRIANGLES: i32 = 8192; // per surface
pub const MD3_MAX_VERTS: i32 = 4096; // per surface
pub const MD3_MAX_SHADERS: i32 = 256; // per surface
pub const MD3_MAX_FRAMES: i32 = 1024; // per model
pub const MD3_MAX_SURFACES: i32 = 32; // per model
pub const MD3_MAX_TAGS: i32 = 16; // per frame

// vertex scales
pub const MD3_XYZ_SCALE: f32 = 1.0 / 64.0;

#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct md3Frame_t {
    bounds: [vec3_t; 2],
    localOrigin: vec3_t,
    radius: f32,
    name: FixedString_16,
}

#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct md3Tag_t {
    name: FixedString_QPATH, // tag name
    origin: vec3_t,
    axis: [vec3_t; 3],
}

/*
** md3Surface_t
**
** CHUNK            SIZE
** header           sizeof( md3Surface_t )
** shaders          sizeof( md3Shader_t ) * numShaders
** triangles[0]     sizeof( md3Triangle_t ) * numTriangles
** st               sizeof( md3St_t ) * numVerts
** XyzNormals       sizeof( md3XyzNormal_t ) * numVerts * numFrames
*/
#[derive(Clone, Copy, Debug)]
#[repr(C)]
pub struct md3Surface_t {
    pub ident: i32,

    pub name: FixedString_QPATH, // polyset name

    pub flags: i32,
    pub numFrames: i32, // all surfaces in a model should have the same

    pub numShaders: i32, // all surfaces in a model should have the same
    pub numVerts: i32,

    pub numTriangles: i32,
    pub ofsTriangles: i32,

    pub ofsShaders: i32,    // offset from start of md3Surface_t
    pub ofsSt: i32,         // texture coords are common for all frames
    pub ofsXyzNormals: i32, // numVerts * numFrames

    pub ofsEnd: i32, // next surface follows
}

#[derive(Clone, Copy, Debug)]
#[repr(C)]
pub struct md3Shader_t {
    pub name: FixedString_QPATH,
    pub shaderIndex: i32, // for in-game use
}

#[derive(Clone, Copy, Debug)]
#[repr(C)]
pub struct md3Triangle_t {
    pub indexes: [i32; 3],
}

#[derive(Clone, Copy, Debug)]
#[repr(C)]
pub struct md3St_t {
    pub st: [f32; 2],
}

#[derive(Clone, Copy, Debug)]
#[repr(C)]
pub struct md3XyzNormal_t {
    pub xyz: [i16; 3],
    pub normal: i16,
}

#[derive(Clone, Copy, Debug)]
#[repr(C)]
pub struct md3Header_t {
    pub ident: i32,
    pub version: i32,

    pub name: FixedString_QPATH, // model name

    pub flags: i32,

    pub numFrames: i32,
    pub numTags: i32,
    pub numSurfaces: i32,

    pub numSkins: i32,

    pub ofsFrames: i32,   // offset for first frame
    pub ofsTags: i32,     // numFrames * numTags
    pub ofsSurfaces: i32, // first surface, others follow

    pub ofsEnd: i32, // end of file
}

/*
==============================================================================

MD4 file format

==============================================================================
*/

pub const fn make_tag_u32(a: char, b: char, c: char, d: char) -> u32 {
    ((a as u32) << 24) | ((b as u32) << 16) | ((c as u32) << 8) | (d as u32)
}

pub const MD4_IDENT: u32 = make_tag_u32('4', 'P', 'D', 'I');
pub const MD4_VERSION: i32 = 1;
pub const MD4_MAX_BONES: i32 = 128;

#[derive(Clone, Copy, Debug)]
#[repr(C)]
pub struct md4Weight_t {
    pub boneIndex: i32,  // these are indexes into the boneReferences,
    pub boneWeight: f32, // not the global per-frame bone list
    pub offset: vec3_t,
}

#[derive(Clone, Copy, Debug)]
#[repr(C)]
pub struct md4Vertex_t {
    pub normal: vec3_t,
    pub texCoords: [vec_t; 2],
    pub numWeights: i32,
    pub weights: [md4Weight_t; 1], // variable sized
}

#[derive(Clone, Copy, Debug)]
#[repr(C)]
pub struct md4Triangle_t {
    pub indexes: [i32; 3],
}

#[derive(Clone, Copy, Debug)]
#[repr(C)]
pub struct md4Surface_t {
    pub ident: i32,

    pub name: FixedString_QPATH, // polyset name
    pub shader: FixedString_QPATH,
    pub shaderIndex: i32, // for in-game use

    pub ofsHeader: i32, // this will be a negative number

    pub numVerts: i32,
    pub ofsVerts: i32,

    pub numTriangles: i32,
    pub ofsTriangles: i32,

    // Bone references are a set of ints representing all the bones
    // present in any vertex weights for this surface.  This is
    // needed because a model may have surfaces that need to be
    // drawn at different sort times, and we don't want to have
    // to re-interpolate all the bones for each surface.
    pub numBoneReferences: i32,
    pub ofsBoneReferences: i32,

    pub ofsEnd: i32, // next surface follows
}

#[derive(Clone, Copy, Debug)]
#[repr(C)]
pub struct md4Bone_t {
    // f32       matrix[3][4];
    pub matrix: [[f32; 3]; 4], // TODO: is the matrix order right??
}

#[derive(Clone, Copy, Debug)]
#[repr(C)]
pub struct md4Frame_t {
    pub bounds: [vec3_t; 2], // bounds of all surfaces of all LOD's for this frame
    pub localOrigin: vec3_t, // midpoint of bounds, used for sphere cull
    pub radius: f32,         // dist from localOrigin to corner
    pub bones: [md4Bone_t; 1], // [numBones]
}

#[derive(Clone, Copy, Debug)]
#[repr(C)]
struct md4LOD_t {
    pub numSurfaces: i32,
    pub ofsSurfaces: i32, // first surface, others follow
    pub ofsEnd: i32,      // next lod follows
}

#[derive(Clone, Copy, Debug)]
#[repr(C)]
pub struct md4Header_t {
    pub ident: i32,
    pub version: i32,

    pub name: FixedString_QPATH, // model name

    // frames and bones are shared by all levels of detail
    pub numFrames: i32,
    pub numBones: i32,
    pub ofsBoneNames: i32, // char name[ MAX_QPATH ]
    pub ofsFrames: i32,    // md4Frame_t[numFrames]

    // each level of detail has completely separate sets of surfaces
    pub numLODs: i32,
    pub ofsLODs: i32,

    pub ofsEnd: i32, // end of file
}

/*
==============================================================================

  .BSP file format

==============================================================================
*/

// little-endian "IBSP"
pub const BSP_IDENT: u32 = make_tag_u32('P', 'S', 'B', 'I');

pub const BSP_VERSION: i32 = 46;

// there shouldn't be any problem with increasing these values at the
// expense of more memory allocation in the utilities
pub const MAX_MAP_MODELS: i32 = 0x400;
pub const MAX_MAP_BRUSHES: i32 = 0x8000;
pub const MAX_MAP_ENTITIES: i32 = 0x800;
pub const MAX_MAP_ENTSTRING: i32 = 0x40000;
pub const MAX_MAP_SHADERS: i32 = 0x400;

pub const MAX_MAP_AREAS: i32 = 0x100; // MAX_MAP_AREA_BYTES in q_shared must match!
pub const MAX_MAP_FOGS: i32 = 0x100;
pub const MAX_MAP_PLANES: i32 = 0x20000;
pub const MAX_MAP_NODES: i32 = 0x20000;
pub const MAX_MAP_BRUSHSIDES: i32 = 0x20000;
pub const MAX_MAP_LEAFS: i32 = 0x20000;
pub const MAX_MAP_LEAFFACES: i32 = 0x20000;
pub const MAX_MAP_LEAFBRUSHES: i32 = 0x40000;
pub const MAX_MAP_PORTALS: i32 = 0x20000;
pub const MAX_MAP_LIGHTING: i32 = 0x800000;
pub const MAX_MAP_LIGHTGRID: i32 = 0x800000;
pub const MAX_MAP_VISIBILITY: i32 = 0x200000;

pub const MAX_MAP_DRAW_SURFS: i32 = 0x20000;
pub const MAX_MAP_DRAW_VERTS: i32 = 0x80000;
pub const MAX_MAP_DRAW_INDEXES: i32 = 0x80000;

// key / value pair sizes in the entities lump
pub const MAX_KEY: i32 = 32;
pub const MAX_VALUE: i32 = 1024;

// the editor uses these predefined yaw angles to orient entities up or down
pub const ANGLE_UP: i32 = -1;
pub const ANGLE_DOWN: i32 = -2;

pub const LIGHTMAP_WIDTH: i32 = 128;
pub const LIGHTMAP_HEIGHT: i32 = 128;

pub const MAX_WORLD_COORD: i32 = 128 * 1024;
pub const MIN_WORLD_COORD: i32 = -128 * 1024;
pub const WORLD_SIZE: i32 = MAX_WORLD_COORD - MIN_WORLD_COORD;

//=============================================================================

#[derive(Copy, Clone, Debug)]
#[repr(C)]
pub struct lump_t {
    pub fileofs: i32,
    pub filelen: i32,
}
unsafe impl FromBytes for lump_t {}

pub type lump_type_t = usize;

pub const LUMP_ENTITIES: lump_type_t = 0;
pub const LUMP_SHADERS: lump_type_t = 1;
pub const LUMP_PLANES: lump_type_t = 2;
pub const LUMP_NODES: lump_type_t = 3;
pub const LUMP_LEAFS: lump_type_t = 4;
pub const LUMP_LEAFSURFACES: lump_type_t = 5;
pub const LUMP_LEAFBRUSHES: lump_type_t = 6;
pub const LUMP_MODELS: lump_type_t = 7;
pub const LUMP_BRUSHES: lump_type_t = 8;
pub const LUMP_BRUSHSIDES: lump_type_t = 9;
pub const LUMP_DRAWVERTS: lump_type_t = 10;
pub const LUMP_DRAWINDEXES: lump_type_t = 11;
pub const LUMP_FOGS: lump_type_t = 12;
pub const LUMP_SURFACES: lump_type_t = 13;
pub const LUMP_LIGHTMAPS: lump_type_t = 14;
pub const LUMP_LIGHTGRID: lump_type_t = 15;
pub const LUMP_VISIBILITY: lump_type_t = 16;
pub const HEADER_LUMPS: usize = 17;

/// The header of the level file, at file offset 0.
#[repr(C)]
#[derive(Copy, Clone, Debug)]
pub struct dheader_t {
    pub ident: i32,
    pub version: i32,
    pub lumps: [lump_t; HEADER_LUMPS],
}
unsafe impl FromBytes for dheader_t {}

#[repr(C)]
#[derive(Copy, Clone, Debug)]
pub struct dmodel_t {
    pub mins: vec3_t,
    pub maxs: vec3_t,
    pub firstSurface: i32,
    pub numSurfaces: i32,
    pub firstBrush: i32,
    pub numBrushes: i32,
}
unsafe impl FromBytes for dmodel_t {}

#[repr(C)]
#[derive(Copy, Clone, Debug)]
pub struct dshader_t {
    pub shader: FixedString_QPATH,
    pub surfaceFlags: i32,
    pub contentFlags: i32,
}
unsafe impl FromBytes for dshader_t {}


// planes x^1 is allways the opposite of plane x

#[repr(C)]
#[derive(Copy, Clone, Debug)]
pub struct dplane_t {
    pub normal: vec3_t,
    pub dist: f32,
}
unsafe impl FromBytes for dplane_t {}

#[repr(C)]
#[derive(Copy, Clone, Debug)]
pub struct dnode_t {
    pub planeNum: i32,
    pub children: [i32; 2], // negative numbers are -(leafs+1), not nodes
    pub mins: [i32; 3],     // for frustom culling
    pub maxs: [i32; 3],
}
unsafe impl FromBytes for dnode_t {}

#[repr(C)]
#[derive(Copy, Clone, Debug)]
pub struct dleaf_t {
    pub cluster: i32, // -1 = opaque cluster (do I still store these?)
    pub area: i32,

    pub mins: [i32; 3], // for frustum culling
    pub maxs: [i32; 3],

    pub firstLeafSurface: i32,
    pub numLeafSurfaces: i32,

    pub firstLeafBrush: i32,
    pub numLeafBrushes: i32,
}
unsafe impl FromBytes for dleaf_t {}

#[repr(C)]
#[derive(Copy, Clone, Debug)]
pub struct dbrushside_t {
    pub planeNum: i32, // positive plane side faces out of the leaf
    pub shaderNum: i32,
}
unsafe impl FromBytes for dbrushside_t {}


#[repr(C)]
#[derive(Copy, Clone, Debug)]
pub struct dbrush_t {
    pub firstSide: i32,
    pub numSides: i32,
    pub shaderNum: i32, // the shader that determines the contents flags
}
unsafe impl FromBytes for dbrush_t {}

pub type c_char = u8;

#[repr(C)]
#[derive(Copy, Clone, Debug)]
pub struct dfog_t {
    pub shader: FixedString_QPATH,
    pub brushNum: i32,
    pub visibleSide: i32, // the brush side that ray tests need to clip against (-1 == none)
}
unsafe impl FromBytes for dfog_t {}

#[repr(C)]
#[derive(Copy, Clone, Debug)]
pub struct drawVert_t {
    pub xyz: vec3_t,
    pub st: [f32; 2],
    pub lightmap: [f32; 2],
    pub normal: vec3_t,
    pub color: [u8; 4],
}
unsafe impl FromBytes for drawVert_t {}

pub type mapSurfaceType_t = i32;
pub const MST_BAD: mapSurfaceType_t = 0;
pub const MST_PLANAR: mapSurfaceType_t = 1;
pub const MST_PATCH: mapSurfaceType_t = 2;
pub const MST_TRIANGLE_SOUP: mapSurfaceType_t = 3;
pub const MST_FLARE: mapSurfaceType_t = 4;

#[repr(C)]
#[derive(Copy, Clone, Debug)]
pub struct dsurface_t {
    pub shaderNum: i32,
    pub fogNum: i32,
    pub surfaceType: i32,

    pub firstVert: i32,
    pub numVerts: i32,

    pub firstIndex: i32,
    pub numIndexes: i32,

    pub lightmapNum: i32,
    pub lightmapX: i32,
    pub lightmapY: i32,
    pub lightmapWidth: i32,
    pub lightmapHeight: i32,

    pub lightmapOrigin: vec3_t,
    pub lightmapVecs: [vec3_t; 3], // for patches, [0] and [1] are lodbounds

    pub patchWidth: i32,
    pub patchHeight: i32,
}
unsafe impl FromBytes for dsurface_t {}
