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

use crate::qcommon::cm_load::CM_ClipHandleToModel;
use crate::prelude::*;
use crate::qcommon::cm_local::*;

pub fn CM_PointLeafnum_r(cm: &clipMap_t, p: vec3_t, mut num: i32) -> i32 {
    while num >= 0 {
        let node = &cm.nodes[num as usize];
        let plane = &cm.planes[node.plane_index as usize];
        let d = if plane.type_ < 3 {
            p[plane.type_ as usize] - plane.dist
        } else {
            DotProduct(plane.normal, p) - plane.dist
        };
        num = if d < 0.0 {
            node.children[1]
        } else {
            node.children[0]
        }
    }

    // c_pointcontents++;      // optimize counter

    -1 - num
}

fn CM_PointLeafnum(cm: &clipMap_t, p: vec3_t) -> i32 {
    if cm.nodes.is_empty() {
        // map not loaded
        return 0;
    }
    CM_PointLeafnum_r(cm, p, 0)
}

/*
======================================================================

LEAF LISTING

======================================================================
*/

/*
=============
CM_BoxLeafnums

Fills in a list of all the leafs touched
=============
*/
pub fn CM_BoxLeafnums_r(
    cm: &mut clipMap_t,
    ll: &mut leafList_t,
    mut nodenum: i32,
    storeLeafs: &mut impl FnMut(&mut clipMap_t, &mut leafList_t, i32),
) {
    loop {
        if nodenum < 0 {
            storeLeafs(cm, ll, nodenum);
            return;
        }

        let node = &cm.nodes[nodenum as usize];
        let children0 = node.children[0];
        let children1 = node.children[1];
        let plane = &cm.planes[node.plane_index as usize];
        let s = BoxOnPlaneSide(ll.bounds.mins, ll.bounds.maxs, plane);
        if s == 1 {
            nodenum = children0;
        } else if s == 2 {
            nodenum = children1;
        } else {
            // go down both
            CM_BoxLeafnums_r(cm, ll, children0, storeLeafs);
            nodenum = children1;
        }
    }
}

fn CM_BoxLeafnums(
    cm: &mut clipMap_t,
    mins: vec3_t,
    maxs: vec3_t,
    list: &mut [i32],
) -> (/*count*/ i32, /*lastLeaf*/ i32) {
    cm.checkcount += 1;

    let mut ll = leafList_t {
        bounds: vec3_bounds { mins, maxs },
        count: 0,
        maxcount: list.len() as i32,
        // list: list,
        lastLeaf: 0,
        overflowed: false,
    };

    CM_BoxLeafnums_r(cm, &mut ll, 0, &mut |cm, ll: &mut leafList_t, nodenum| {
        // was: CM_StoreLeafs
        let leafNum = -1 - nodenum;

        // store the lastLeaf even if the list is overflowed
        if cm.leafs[leafNum as usize].cluster != -1 {
            ll.lastLeaf = leafNum;
        }

        if ll.count >= ll.maxcount {
            ll.overflowed = true;
            return;
        }
        list[ll.count as usize] = leafNum;
        ll.count += 1;
    });

    (ll.count, ll.lastLeaf)
}

/// `list` returns a list of cbrush_t indices (in cm.leafbrushes)
fn CM_BoxBrushes(cm: &mut clipMap_t, mins: vec3_t, maxs: vec3_t, list: &mut [i32]) -> i32 {
    cm.checkcount += 1;

    let mut ll = leafList_t {
        bounds: vec3_bounds { mins, maxs },
        count: 0,
        maxcount: list.len() as i32,
        //list: (void *)list,
        //storeLeafs: CM_StoreBrushes,
        lastLeaf: 0,
        overflowed: false,
    };
    CM_BoxLeafnums_r(cm, &mut ll, 0, &mut |cm, ll: &mut leafList_t, nodenum| {
        // was: CM_StoreBrushes
        let leafnum = -1 - nodenum;
        let leaf = &cm.leafs[leafnum as usize];

        for k in 0..leaf.numLeafBrushes as usize {
            let brushnum = cm.leafbrushes[leaf.firstLeafBrush as usize + k];
            let b: &mut cbrush_t = &mut cm.brushes[brushnum as usize];
            if b.checkcount == cm.checkcount {
                continue; // already checked this brush in another leaf
            }
            b.checkcount = cm.checkcount;

            let any_bounds_misordered = (0..3).any(|i| {
                b.bounds.mins[i] >= ll.bounds.maxs[i] || b.bounds.maxs[i] <= ll.bounds.mins[i]
            });
            if any_bounds_misordered {
                continue;
            }
            if ll.count >= ll.maxcount {
                ll.overflowed = true;
                return;
            }
            list[ll.count as usize] = brushnum;
            ll.count += 1;
        }
        /*
        #if 0
            // store patches?
            for ( k = 0 ; k < leaf.numLeafSurfaces ; k++ ) {
                patch = cm.surfaces[ cm.leafsurfaces[ leaf.firstleafsurface + k ] ];
                if ( !patch ) {
                    continue;
                }
            }
        #endif
        */
    });
    ll.count
}

//====================================================================

fn CM_PointContents(cm: &clipMap_t, p: vec3_t, model: clipHandle_t) -> i32 {
    if cm.nodes.is_empty() {
        // map not loaded
        return 0;
    }

    let leaf: &cLeaf_t = if model != 0 {
        let clipm = CM_ClipHandleToModel(cm, model);
        &clipm.leaf
    } else {
        let leafnum = CM_PointLeafnum_r(cm, p, 0);
        &cm.leafs[leafnum as usize]
    };

    let mut contents: i32 = 0;
    for b in leaf.leaf_brushes(cm) {
        // see if the point is in the brush
        let mut early_break = false;
        for side in b.sides(cm).iter() {
            let side_plane = side.plane(cm);
            let d = p.dot(side_plane.normal);
            // FIXME test for Cash
            //          if ( d >= b.sides[i].plane.dist ) {
            if d > side_plane.dist {
                early_break = true;
                break;
            }
        }

        if !early_break {
            contents |= b.contents;
        }
    }

    contents
}

/*
==================
CM_TransformedPointContents

Handles offseting and rotation of the end points for moving and
rotating entities
==================
*/
fn CM_TransformedPointContents(
    cm: &clipMap_t,
    p: vec3_t,
    model: clipHandle_t,
    origin: vec3_t,
    angles: vec3_t,
) -> i32 {
    // subtract origin offset
    let mut p_l = p - origin;

    // rotate start and end into the models frame of reference
    if model != BOX_MODEL_HANDLE && (angles[0] != 0.0 || angles[1] != 0.0 || angles[2] != 0.0) {
        let (forward, right, up) = AngleVectors(angles);
        let temp = p_l;
        p_l[0] = temp.dot(forward);
        p_l[1] = -temp.dot(right);
        p_l[2] = temp.dot(up);
    }

    CM_PointContents(cm, p_l, model)
}

/*
===============================================================================

PVS

===============================================================================
*/

pub fn CM_ClusterPVS(cm: &clipMap_t, cluster: i32) -> &[u8] {
    if cluster < 0 || cluster >= cm.numClusters || !cm.vised {
        return &cm.visibility;
    }

    &cm.visibility[cluster as usize * cm.clusterBytes as usize..]
}

/*
===============================================================================

AREAPORTALS

===============================================================================
*/

fn CM_FloodArea_r(cm: &mut clipMap_t, areaNum: i32, floodnum: i32) {
    let area = &mut cm.areas[areaNum as usize];

    if area.floodvalid == cm.floodvalid {
        if area.floodnum == floodnum {
            return;
        }
        panic!("FloodArea_r: reflooded");
    }

    area.floodnum = floodnum;
    area.floodvalid = cm.floodvalid;
    let con_offset = areaNum as usize * cm.areas.len(); // index into cm.areaPortals
    for i in 0..cm.areas.len() {
        if cm.areaPortals[con_offset + i] > 0 {
            CM_FloodArea_r(cm, i as i32, floodnum);
        }
    }
}

pub fn CM_FloodAreaConnections(cm: &mut clipMap_t) {
    // all current floods are now invalid
    cm.floodvalid += 1;
    let mut floodnum = 0;

    for i in 0..cm.areas.len() {
        let area: &mut cArea_t = &mut cm.areas[i];
        if area.floodvalid == cm.floodvalid {
            continue; // already flooded into
        }
        floodnum += 1;
        CM_FloodArea_r(cm, i as i32, floodnum);
    }
}

pub fn CM_AdjustAreaPortalState(cm: &mut clipMap_t, area1: i32, area2: i32, open: bool) {
    if area1 < 0 || area2 < 0 {
        return;
    }

    let area1 = area1 as usize;
    let area2 = area2 as usize;

    assert!(area1 < cm.areas.len());
    assert!(area2 < cm.areas.len());

    if open {
        cm.areaPortals[area1 * cm.areas.len() + area2] += 1;
        cm.areaPortals[area2 * cm.areas.len() + area1] += 1;
    } else {
        cm.areaPortals[area1 * cm.areas.len() + area2] -= 1;
        cm.areaPortals[area2 * cm.areas.len() + area1] -= 1;
        if cm.areaPortals[area2 * cm.areas.len() + area1] < 0 {
            panic!("CM_AdjustAreaPortalState: negative reference count");
        }
    }

    CM_FloodAreaConnections(cm);
}

pub fn CM_AreasConnected(cm: &clipMap_t, area1: i32, area2: i32) -> bool {
    /*
    #ifndef BSPC
        if cm_noAreas.integer {
            return true;
        }
    #endif
    */

    if area1 < 0 || area2 < 0 {
        return false;
    }

    assert!((area1 as usize) < cm.areas.len());
    assert!((area2 as usize) < cm.areas.len());
    cm.areas[area1 as usize].floodnum == cm.areas[area2 as usize].floodnum
}

/*
=================
CM_WriteAreaBits

Writes a bit vector of all the areas
that are in the same flood as the area parameter
Returns the number of bytes needed to hold all the bits.

The bits are OR'd in, so you can CM_WriteAreaBits from multiple
viewpoints and get the union of all visible areas.

This is used to cull non-visible entities from snapshots
=================
*/
pub fn CM_WriteAreaBits(cm: &clipMap_t, buffer: &mut [u8], area: i32) -> i32 {
    let bytes = (cm.areas.len() + 7) >> 3;

    /*
    #ifndef BSPC
        if (cm_noAreas.integer || area == -1)
    #else
    */
    if area == -1
    /*
    #endif
    */
    {
        // for debugging, send everything
        for b in buffer[0..bytes].iter_mut() {
            *b = 255;
        }
    } else {
        let floodnum = cm.areas[area as usize].floodnum;
        for (i, area) in cm.areas.iter().enumerate() {
            if area.floodnum == floodnum {
                buffer[i >> 3] |= 1u8 << (i & 7);
            }
        }
    }

    bytes as i32
}
