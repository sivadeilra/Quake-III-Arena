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

//! Model Loading

use crate::game::surfaceflags::CONTENTS_BODY;
use crate::prelude::*;
use crate::qcommon::cm_local::*;
use crate::qcommon::cm_patch::*;
use crate::qcommon::cm_test::*;
use crate::qcommon::md4::Com_BlockChecksum;
use crate::qfiles::*;
use crate::range_len;
use crate::zerocopy::*;
use log::debug;
use std::ops::Range;

/*
#ifdef BSPC

#include "../bspc/l_qfiles.h"

pub fn SetPlaneSignbits (cplane_t *out) {
    int bits, j;

    // for fast box on planeside test
    bits = 0;
    for (j=0 ; j<3 ; j++) {
        if (out.normal[j] < 0) {
            bits |= 1<<j;
        }
    }
    out.signbits = bits;
}
#endif //BSPC
*/

pub fn get_sign_bits(v: vec3_t) -> u8 {
    ((v.x < 0.0) as u8) | (((v.y < 0.0) as u8) << 1) | (((v.z < 0.0) as u8) << 2)
}

// to allow boxes to be treated as brush models, we allocate
// some extra indexes along with those needed by the map
pub const BOX_BRUSHES: i32 = 1;
pub const BOX_SIDES: i32 = 6;
pub const BOX_LEAFS: i32 = 2;
pub const BOX_PLANES: i32 = 12;

fn LittleLong(x: i32) -> i32 {
    x
}
fn LittleLong_u32(x: u32) -> u32 {
    x
}
fn LL(x: i32) -> i32 {
    x
}
//#define LL(x) x=LittleLong(x)

fn LittleFloat(x: f32) -> f32 {
    x
}

fn LittleVec3(v: vec3_t) -> vec3_t {
    v
}

/*
clipMap_t   cm;
int         c_pointcontents;
int         c_traces, c_brush_traces, c_patch_traces;

#ifndef BSPC
cvar_t      *cm_noAreas;
cvar_t      *cm_noCurves;
cvar_t      *cm_playerCurveClip;
#endif

*/

struct FileLoader<'a> {
    data: &'a [u8],
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Error {
    BadFileRange,
    RequiredLumpEmpty,
    Str(&'static str),
    String(String),
}

impl From<&'static str> for Error {
    fn from(value: &'static str) -> Error {
        Error::Str(value)
    }
}

impl From<String> for Error {
    fn from(value: String) -> Error {
        Error::String(value)
    }
}

impl<'a> FileLoader<'a> {
    pub fn get_range(&self, range: Range<usize>) -> Result<&[u8], Error> {
        if range.start <= self.data.len() && range.end <= self.data.len() {
            Ok(&self.data[range])
        } else {
            Err(Error::BadFileRange)
        }
    }

    pub fn get_lump(&self, lump: lump_t) -> Result<&[u8], Error> {
        self.get_range(lump.fileofs as usize..lump.fileofs as usize + lump.filelen as usize)
    }

    pub fn get_lump_required(&self, lump: lump_t) -> Result<&[u8], Error> {
        if lump.filelen == 0 {
            return Err(Error::RequiredLumpEmpty);
        }
        self.get_range(lump.fileofs as usize..lump.fileofs as usize + lump.filelen as usize)
    }

    pub fn get_lump_required_as_slice<T: FromBytes>(&self, lump: lump_t) -> Result<&[T], Error> {
        Ok(slice_from_bytes::<T>(self.get_lump_required(lump)?))
    }
}

fn CMod_LoadShaders(lump: lump_t, file: &FileLoader) -> Result<Vec<dshader_t>, Error> {
    Ok(slice_from_bytes::<dshader_t>(file.get_lump_required(lump)?)
        .iter()
        .map(|dshader| dshader_t {
            shader: dshader.shader,
            contentFlags: LittleLong(dshader.contentFlags),
            surfaceFlags: LittleLong(dshader.surfaceFlags),
        })
        .collect())
}

fn CMod_LoadSubmodels(
    leafbrushes: &mut Vec<i32>,
    leafsurfaces: &mut Vec<i32>,
    lump: lump_t,
    file: &FileLoader,
) -> Result<Vec<cmodel_t>, Error> {
    let models = file.get_lump_required_as_slice::<dmodel_t>(lump)?;
    if models.len() > MAX_SUBMODELS {
        return Err(Error::Str("MAX_SUBMODELS exceeded"));
    }

    Ok(models
        .iter()
        .enumerate()
        .map(|(i, model)| {
            cmodel_t {
                // spread the mins / maxs by a pixel
                mins: model.mins - vec3_t::from_scalar(1.0),
                maxs: model.maxs + vec3_t::from_scalar(1.0),
                leaf: if i != 0 {
                    // make a "leaf" just to hold the model's brushes and surfaces
                    let numLeafBrushes = LittleLong(model.numBrushes);
                    let firstLeafBrush = leafbrushes.len() as i32;
                    leafbrushes.extend(range_len(LittleLong(model.firstBrush), numLeafBrushes));

                    let numLeafSurfaces = LittleLong(model.numSurfaces);
                    let firstLeafSurface = leafsurfaces.len() as i32;
                    leafsurfaces.extend(range_len(LittleLong(model.firstSurface), numLeafSurfaces));

                    cLeaf_t {
                        area: -1,
                        cluster: -1,
                        numLeafBrushes,
                        firstLeafBrush,
                        numLeafSurfaces,
                        firstLeafSurface,
                    }
                } else {
                    // world model doesn't need other info
                    cLeaf_t::default()
                },
            }
        })
        .collect())
}

fn CMod_LoadNodes(lump: lump_t, file: &FileLoader) -> Result<Vec<cNode_t>, Error> {
    let dnodes = file.get_lump_required_as_slice::<dnode_t>(lump)?;
    if dnodes.is_empty() {
        return Err(Error::Str("Map has no nodes"));
    }

    Ok(dnodes
        .iter()
        .map(|dnode| cNode_t {
            plane_index: LittleLong(dnode.planeNum),
            children: [LittleLong(dnode.children[0]), LittleLong(dnode.children[1])],
        })
        .collect::<Vec<cNode_t>>())
}

fn check_index<T>(index: i32, items: &[T]) -> Result<usize, Error> {
    if index < 0 {
        return Err(Error::Str("Index is negative"));
    }
    let index_u = index as usize;
    if index_u >= items.len() {
        return Err(Error::Str("Index is out of range"));
    }
    Ok(index_u)
}

fn get_item_checked<T>(index: i32, items: &[T]) -> Result<&T, Error> {
    if index < 0 {
        return Err(Error::Str("Index is negative"));
    }
    let index_u = index as usize;
    if index_u >= items.len() {
        return Err(Error::Str("Index is out of range"));
    }
    Ok(&items[index_u])
}

fn CMod_LoadBrushes(
    shaders: &[dshader_t],
    brushsides: &[cbrushside_t],
    planes: &[cplane_t],
    lump: lump_t,
    file: &FileLoader,
) -> Result<Vec<cbrush_t>, Error> {
    let bound_brush = move |firstSide: i32, numSides: i32| -> vec3_bounds {
        let b_sides = &brushsides[range_len(firstSide as usize, numSides as usize)];
        let pd = move |s: usize| -> f32 { planes[b_sides[s].plane_num as usize].dist };
        vec3_bounds {
            mins: v3(-pd(0), -pd(2), -pd(4)),
            maxs: v3(pd(1), pd(3), pd(5)),
        }
    };
    file.get_lump_required_as_slice::<dbrush_t>(lump)?
        .iter()
        .map(|dbrush| {
            let shaderNum = LittleLong(dbrush.shaderNum);
            let shader = get_item_checked(shaderNum, shaders)?;
            let firstSide = LittleLong(dbrush.firstSide);
            let numSides = LittleLong(dbrush.numSides);
            Ok(cbrush_t {
                firstSide: firstSide,
                numsides: numSides,
                shaderNum: shaderNum,
                contents: shader.contentFlags,
                bounds: bound_brush(firstSide, numSides),
            })
        })
        .collect()
}

struct LoadLeafsOutput {
    pub leafs: Vec<cLeaf_t>,
    pub num_clusters: usize,
    pub num_areas: usize,
    pub num_area_portals: usize,
}

fn CMod_LoadLeafs(lump: lump_t, file: &FileLoader) -> Result<LoadLeafsOutput, Error> {
    let mut max_area: usize = 0;
    let mut max_cluster: usize = 0;
    let leafs: Vec<cLeaf_t> = file
        .get_lump_required_as_slice::<dleaf_t>(lump)?
        .iter()
        .map(|dleaf| -> cLeaf_t {
            let cluster = LittleLong(dleaf.cluster);
            let area = LittleLong(dleaf.area);
            max_area = max_area.max(area as usize);
            max_cluster = max_cluster.max(cluster as usize);
            cLeaf_t {
                cluster,
                area,
                firstLeafBrush: LittleLong(dleaf.firstLeafBrush),
                numLeafBrushes: LittleLong(dleaf.numLeafBrushes),
                firstLeafSurface: LittleLong(dleaf.firstLeafSurface),
                numLeafSurfaces: LittleLong(dleaf.numLeafSurfaces),
            }
        })
        .collect();

    let num_clusters = max_cluster + 1;
    let num_areas = max_area + 1;

    Ok(LoadLeafsOutput {
        leafs,
        num_areas,
        num_area_portals: num_areas * num_areas,
        num_clusters,
    })
}

fn CMod_LoadPlanes(lump: lump_t, file: &FileLoader) -> Result<Vec<cplane_t>, Error> {
    let planes: &[dplane_t] = file.get_lump_required_as_slice::<dplane_t>(lump)?;
    if planes.is_empty() {
        return Err(Error::Str("Map with no planes"));
    }
    Ok(planes
        .iter()
        .map(|dplane| -> cplane_t {
            let normal = LittleVec3(dplane.normal);
            let mut bits: u8 = 0;
            for j in 0..3 {
                if normal[j] < 0.0 {
                    bits |= 1u8 << j;
                }
            }
            cplane_t {
                normal: normal,
                dist: LittleFloat(dplane.dist),
                type_: PlaneTypeForNormal(normal),
                signbits: bits,
            }
        })
        .collect())
}

fn CMod_LoadLeafBrushes(lump: lump_t, file: &FileLoader) -> Result<Vec<i32>, Error> {
    Ok(file
        .get_lump_required_as_slice::<i32>(lump)?
        .iter()
        .map(|&i| LittleLong(i))
        .collect())
}

fn CMod_LoadLeafSurfaces(lump: lump_t, file: &FileLoader) -> Result<Vec<i32>, Error> {
    Ok(file
        .get_lump_required_as_slice::<i32>(lump)?
        .iter()
        .map(|&i| LittleLong(i))
        .collect())
}

fn CMod_LoadBrushSides(
    shaders: &[dshader_t],
    lump: lump_t,
    file: &FileLoader,
) -> Result<Vec<cbrushside_t>, Error> {
    file.get_lump_required_as_slice::<dbrushside_t>(lump)?
        .iter()
        .map(|dbrushside| {
            let plane_num = LittleLong(dbrushside.planeNum);
            let shaderNum = LittleLong(dbrushside.shaderNum);
            let shader = get_item_checked(shaderNum, shaders)?;
            Ok(cbrushside_t {
                plane_num,
                shaderNum,
                surfaceFlags: shader.surfaceFlags,
            })
        })
        .collect()
}

fn CMod_LoadEntityString(lump: lump_t, file: &FileLoader) -> Result<String, Error> {
    let bytes = file.get_lump(lump)?;
    std::str::from_utf8(bytes)
        .map_err(|_| Error::Str("Invalid entity string"))
        .map(|s| s.to_string())
}

const VIS_HEADER: usize = 8;
fn CMod_LoadVisibility(
    num_clusters: usize,
    lump: lump_t,
    file: &FileLoader,
) -> Result<ClipMapVis, Error> {
    let in_ = file.get_lump(lump)?;

    if in_.is_empty() {
        let cluster_bytes = (num_clusters + 31) & !31;
        return Ok(ClipMapVis {
            vised: false, // added
            visibility: vec![255u8; cluster_bytes],
            numClusters: 0,
            clusterBytes: cluster_bytes as i32,
        });
    }

    if in_.len() < VIS_HEADER {
        return Err(Error::Str("Invalid vis lump"));
    }

    let numClusters = LittleLong(*from_bytes::<i32>(&in_[0..4]));
    let clusterBytes = LittleLong(*from_bytes::<i32>(&in_[4..8]));

    Ok(ClipMapVis {
        vised: true,
        numClusters: numClusters,
        clusterBytes: clusterBytes,
        visibility: in_[VIS_HEADER..].to_vec(),
    })
}

const MAX_PATCH_VERTS: usize = 1024;

fn CMod_LoadPatches(
    shaders: &[dshader_t],
    surfs: lump_t,
    verts: lump_t,
    file: &FileLoader,
) -> Result<Vec<Option<Box<cPatch_t>>>, Error> {
    let mut points: [vec3_t; MAX_PATCH_VERTS] = [vec3_t::ORIGIN; MAX_PATCH_VERTS];
    let dpatches = slice_from_bytes::<dsurface_t>(file.get_lump(surfs)?);
    let dv = slice_from_bytes::<drawVert_t>(file.get_lump(verts)?);

    // scan through all the surfaces, but only load patches,
    // not planar faces
    dpatches
        .iter()
        .map(|dpatch| {
            if LittleLong(dpatch.surfaceType) != MST_PATCH {
                return Ok(None); // ignore other surfaces
            }
            // FIXME: check for non-colliding patches

            // load the full drawverts onto the stack
            let width = LittleLong(dpatch.patchWidth);
            let height = LittleLong(dpatch.patchHeight);
            let c = width as usize * height as usize;
            if c > MAX_PATCH_VERTS {
                return Err(Error::Str("ParseMesh: MAX_PATCH_VERTS"));
            }

            let first_vert = LittleLong(dpatch.firstVert) as usize;
            for (p, dv_p) in points[..c]
                .iter_mut()
                .zip(dv[first_vert..first_vert + c].iter())
            {
                p[0] = LittleFloat(dv_p.xyz[0]);
                p[1] = LittleFloat(dv_p.xyz[1]);
                p[2] = LittleFloat(dv_p.xyz[2]);
            }

            // create the internal facet structure
            let pc = CM_GeneratePatchCollide(width as usize, height as usize, &points);

            let shaderNum = LittleLong(dpatch.shaderNum);
            let shader = get_item_checked(shaderNum, shaders)?;

            Ok(Some(Box::new(cPatch_t {
                contents: shader.contentFlags,
                surfaceFlags: shader.surfaceFlags,
                pc: pc,
                checkcount: 0,
            })))
        })
        .collect()
}

fn CM_LumpChecksum(lump: lump_t, file: &FileLoader) -> Result<u32, Error> {
    let data = file.get_lump(lump)?;
    Ok(LittleLong_u32(Com_BlockChecksum(data)))
}

fn CM_Checksum(header: &dheader_t, file: &FileLoader) -> Result<u32, Error> {
    let mut checksums: [u32; 16] = [0; 16];
    checksums[0] = CM_LumpChecksum(header.lumps[LUMP_SHADERS], file)?;
    checksums[1] = CM_LumpChecksum(header.lumps[LUMP_LEAFS], file)?;
    checksums[2] = CM_LumpChecksum(header.lumps[LUMP_LEAFBRUSHES], file)?;
    checksums[3] = CM_LumpChecksum(header.lumps[LUMP_LEAFSURFACES], file)?;
    checksums[4] = CM_LumpChecksum(header.lumps[LUMP_PLANES], file)?;
    checksums[5] = CM_LumpChecksum(header.lumps[LUMP_BRUSHSIDES], file)?;
    checksums[6] = CM_LumpChecksum(header.lumps[LUMP_BRUSHES], file)?;
    checksums[7] = CM_LumpChecksum(header.lumps[LUMP_MODELS], file)?;
    checksums[8] = CM_LumpChecksum(header.lumps[LUMP_NODES], file)?;
    checksums[9] = CM_LumpChecksum(header.lumps[LUMP_SURFACES], file)?;
    checksums[10] = CM_LumpChecksum(header.lumps[LUMP_DRAWVERTS], file)?;
    Ok(LittleLong_u32(Com_BlockChecksum(
        &slice_to_bytes(&checksums)[0..11 * 4],
    )))
}

/// Loads in the map and all submodels
pub fn CM_LoadMap(name: &str, clientload: bool, _checksum: &mut i32) -> Result<clipMap_t, Error> {
    assert!(!name.is_empty());

    /*
    #ifndef BSPC
        cm_noAreas = Cvar_Get ("cm_noAreas", "0", CVAR_CHEAT);
        cm_noCurves = Cvar_Get ("cm_noCurves", "0", CVAR_CHEAT);
        cm_playerCurveClip = Cvar_Get ("cm_playerCurveClip", "1", CVAR_ARCHIVE|CVAR_CHEAT );
    #endif
    */
    debug!("CM_LoadMap( {}, {} )", name, clientload);

    /*
        if ( !strcmp( cm.name, name ) && clientload ) {
            *checksum = last_checksum;
            return;
        }
    */

    // free old stuff
    CM_ClearLevelPatches();

    /*
        if ( !name[0] ) {
            cm.numLeafs = 1;
            cm.numClusters = 1;
            cm.numAreas = 1;
            cm.cmodels = Hunk_Alloc( sizeof( *cm.cmodels ), h_high );
            *checksum = 0;
            return;
        }
    */

    //
    // load the file
    //

    /*
    #ifndef BSPC
        length = FS_ReadFile( name, (void **)&buf );
    #else
        length = LoadQuakeFile((quakefile_t *) name, (void **)&buf);
    #endif
    */
    use std::io::Read;
    let mut buf: Vec<u8> = Vec::new();
    let mut f = std::fs::File::open(name).unwrap();
    f.read_to_end(&mut buf).unwrap();

    // last_checksum = LittleLong(Com_BlockChecksum(buf, length));
    // *checksum = last_checksum;

    let header: &dheader_t = from_bytes::<dheader_t>(&buf);
    /* TODO:
        for (i=0 ; i<sizeof(dheader_t)/4 ; i++) {
            ((int *)&header)[i] = LittleLong( ((int *)&header)[i]);
        }
    */

    if header.version != BSP_VERSION {
        return Err(Error::String(format!(
            "CM_LoadMap: {} has wrong version number ({} should be {})",
            name, header.version, BSP_VERSION
        )));
    }

    let file: FileLoader<'_> = FileLoader { data: &buf };

    // load into heap
    let shaders = CMod_LoadShaders(header.lumps[LUMP_SHADERS], &file)?;
    let LoadLeafsOutput {
        leafs,
        num_areas,
        num_clusters,
        num_area_portals,
    } = CMod_LoadLeafs(header.lumps[LUMP_LEAFS], &file)?;
    let mut leafbrushes = CMod_LoadLeafBrushes(header.lumps[LUMP_LEAFBRUSHES], &file)?;
    let mut leafsurfaces = CMod_LoadLeafSurfaces(header.lumps[LUMP_LEAFSURFACES], &file)?;
    let planes = CMod_LoadPlanes(header.lumps[LUMP_PLANES], &file)?;
    let brushsides = CMod_LoadBrushSides(&shaders, header.lumps[LUMP_BRUSHSIDES], &file)?;
    let brushes = CMod_LoadBrushes(
        &shaders,
        &brushsides,
        &planes,
        header.lumps[LUMP_BRUSHES],
        &file,
    )?;
    let cmodels = CMod_LoadSubmodels(
        &mut leafbrushes,
        &mut leafsurfaces,
        header.lumps[LUMP_MODELS],
        &file,
    )?;
    let nodes = CMod_LoadNodes(header.lumps[LUMP_NODES], &file)?;
    let entityString = CMod_LoadEntityString(header.lumps[LUMP_ENTITIES], &file)?;
    let vis = CMod_LoadVisibility(num_clusters, header.lumps[LUMP_VISIBILITY], &file)?;
    let surfaces = CMod_LoadPatches(
        &shaders,
        header.lumps[LUMP_SURFACES],
        header.lumps[LUMP_DRAWVERTS],
        &file,
    )?;

    let mut cm = clipMap_t {
        shaders,
        leafs,
        leafbrushes,
        leafsurfaces,
        planes,
        brushsides,
        brushes,
        cmodels,
        nodes,
        entityString,
        surfaces,

        areas: vec![Default::default(); num_areas],
        areaPortals: vec![Default::default(); num_area_portals],

        // vis
        numClusters: vis.numClusters,
        clusterBytes: vis.clusterBytes,
        visibility: vis.visibility,
        vised: vis.vised,

        //
        box_brush_index: 0,
        box_model: Default::default(),
        box_planes_range: 0..0,
        floodvalid: 0,
        name: name.to_string(),
    };

    CM_InitBoxHull(&mut cm);

    CM_FloodAreaConnections(&mut cm);

    // allow this to be cached if it is loaded by the server
    if !clientload {
        cm.name = name.to_string();
    }

    Ok(cm)
}

/*
pub fn CM_ClearMap(cm: &mut clipMap_t) {
    *cm = clipMap_t::default();
    CM_ClearLevelPatches();
}
*/

pub fn CM_ClipHandleToModel(cm: &clipMap_t, handle: clipHandle_t) -> &cmodel_t {
    assert!(handle >= 0);
    if (handle as usize) < cm.cmodels.len() {
        return &cm.cmodels[handle as usize];
    }
    if handle == BOX_MODEL_HANDLE {
        return &cm.box_model;
    }
    if (handle as usize) < MAX_SUBMODELS {
        panic!(
            "CM_ClipHandleToModel: bad handle {} < {} < {}",
            cm.cmodels.len(),
            handle,
            MAX_SUBMODELS
        );
    }
    panic!(
        "CM_ClipHandleToModel: bad handle {}",
        handle as usize + MAX_SUBMODELS
    );
}

fn CM_InlineModel(cm: &clipMap_t, index: i32) -> clipHandle_t {
    assert!(index >= 0);
    assert!((index as usize) < cm.cmodels.len());
    index
}

pub fn CM_NumClusters(cm: &clipMap_t) -> i32 {
    cm.numClusters as i32
}

pub fn CM_NumInlineModels(cm: &clipMap_t) -> i32 {
    cm.cmodels.len() as i32
}

pub fn CM_EntityString(cm: &clipMap_t) -> &str {
    &cm.entityString
}

pub fn CM_LeafCluster(cm: &clipMap_t, leafnum: i32) -> i32 {
    assert!(leafnum >= 0);
    assert!((leafnum as usize) < cm.leafs.len());
    cm.leafs[leafnum as usize].cluster
}

pub fn CM_LeafArea(cm: &clipMap_t, leafnum: i32) -> i32 {
    assert!(leafnum >= 0);
    assert!((leafnum as usize) < cm.leafs.len());
    cm.leafs[leafnum as usize].area
}

//=======================================================================

/// Set up the planes and nodes so that the six floats of a bounding box
/// can just be stored out and get a proper clipping hull structure.
pub fn CM_InitBoxHull(cm: &mut clipMap_t) {
    cm.box_brush_index = cm.brushes.len();
    cm.brushes.push(cbrush_t {
        numsides: 6,
        firstSide: cm.brushsides.len() as i32,
        contents: CONTENTS_BODY,
        bounds: vec3_bounds::default(),
        shaderNum: -1,
    });

    let box_model = cmodel_t {
        leaf: cLeaf_t {
            numLeafBrushes: 1,
            firstLeafBrush: cm.leafbrushes.len() as i32,
            ..Default::default()
        },
        maxs: vec3_t::ORIGIN,
        mins: vec3_t::ORIGIN,
    };

    cm.leafbrushes.push(cm.brushes.len() as i32);

    let firstPlane = cm.planes.len();

    // add sides
    for i in 0..6 {
        let side = i & 1;
        cm.brushsides.push(cbrushside_t {
            shaderNum: 0,
            plane_num: firstPlane as i32 + i * 2 + side,
            surfaceFlags: 0,
        });
    }

    // add planes
    for i in 0..6 {
        let axis = i >> 1;
        let normal = vec3_t::unit(axis);
        let inverse_normal = -normal;

        cm.planes.push(cplane_t {
            type_: axis,
            signbits: get_sign_bits(normal),
            normal: normal,
            dist: 0.0,
        });

        cm.planes.push(cplane_t {
            type_: 3 + axis,
            signbits: get_sign_bits(inverse_normal),
            normal: inverse_normal,
            dist: 0.0,
        });
    }

    cm.box_planes_range = range_len(firstPlane, 12);
    cm.box_model = box_model;
}

/// To keep everything totally uniform, bounding boxes are turned into small
/// BSP trees instead of being compared directly.
/// Capsules are handled differently though.
pub fn CM_TempBoxModel(
    cm: &mut clipMap_t,
    mins: vec3_t,
    maxs: vec3_t,
    capsule: bool,
) -> clipHandle_t {
    cm.box_model.mins = mins;
    cm.box_model.maxs = maxs;

    if capsule {
        return CAPSULE_MODEL_HANDLE;
    }

    let box_planes = cm.box_planes_mut();
    box_planes[0].dist = maxs[0];
    box_planes[1].dist = -maxs[0];
    box_planes[2].dist = mins[0];
    box_planes[3].dist = -mins[0];
    box_planes[4].dist = maxs[1];
    box_planes[5].dist = -maxs[1];
    box_planes[6].dist = mins[1];
    box_planes[7].dist = -mins[1];
    box_planes[8].dist = maxs[2];
    box_planes[9].dist = -maxs[2];
    box_planes[10].dist = mins[2];
    box_planes[11].dist = -mins[2];

    let box_brush = cm.box_brush_mut();
    box_brush.bounds.mins = mins;
    box_brush.bounds.maxs = maxs;

    BOX_MODEL_HANDLE
}

pub fn CM_ModelBounds(cm: &clipMap_t, model: clipHandle_t) -> vec3_bounds {
    let cmod = CM_ClipHandleToModel(cm, model);
    vec3_bounds {
        mins: cmod.mins,
        maxs: cmod.maxs,
    }
}
