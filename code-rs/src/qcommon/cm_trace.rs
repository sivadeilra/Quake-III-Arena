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

use crate::game::surfaceflags::CONTENTS_BODY;
use crate::prelude::*;
use crate::qcommon::cm_load::*;
use crate::qcommon::cm_local::*;
use crate::qcommon::cm_patch::*;
use crate::qcommon::cm_test::CM_BoxLeafnums_r;

// always use bbox vs. bbox collision and never capsule vs. bbox or vice versa
//#define ALWAYS_BBOX_VS_BBOX
// always use capsule vs. capsule collision and never capsule vs. bbox or vice versa
//#define ALWAYS_CAPSULE_VS_CAPSULE

//#define CAPSULE_DEBUG

/*
===============================================================================

BASIC MATH

===============================================================================
*/

pub type mat3x3 = [vec3_t; 3];

pub fn RotatePoint_mut(v: &mut vec3_t, matrix: &mat3x3) {
    // bk: FIXME
    let result = vec3_t {
        x: matrix[0].dot(*v),
        y: matrix[1].dot(*v),
        z: matrix[2].dot(*v),
    };
    *v = result;
}

pub fn TransposeMatrix(m: &mat3x3) -> mat3x3 {
    // bk: FIXME
    [
        v3(m[0].x, m[1].x, m[2].x),
        v3(m[0].y, m[1].y, m[2].y),
        v3(m[0].z, m[1].z, m[2].z),
    ]
}

pub fn CreateRotationMatrix(angles: vec3_t) -> mat3x3 {
    let mut m = AngleVectors(angles);
    m.1 = -m.1;
    [m.0, m.1, m.2]
}

pub fn CM_ProjectPointOntoVector(point: vec3_t, vStart: vec3_t, vDir: vec3_t) -> vec3_t {
    let pvec = point - vStart;
    // project onto the directional vector for this segment
    VectorMA(vStart, pvec.dot(vDir), vDir)
}

pub fn CM_DistanceFromLineSquared(p: vec3_t, lp1: vec3_t, lp2: vec3_t, dir: vec3_t) -> f32 {
    let proj = CM_ProjectPointOntoVector(p, lp1, dir);
    for j in 0..3 {
        if (proj[j] > lp1[j] && proj[j] > lp2[j]) || (proj[j] < lp1[j] && proj[j] < lp2[j]) {
            return (if (proj[j] - lp1[j]).abs() < (proj[j] - lp2[j]).abs() {
                p - lp1
            } else {
                p - lp2
            })
            .length2();
        }
    }
    (p - proj).length2()
}

pub fn CM_VectorDistanceSquared(p1: vec3_t, p2: vec3_t) -> f32 {
    (p2 - p1).length2()
}

pub fn SquareRootFloat(number: f32) -> f32 {
    /*
    long i;
    f32 x, y;
    const f32 f = 1.5F;

    x = number * 0.5F;
    y  = number;
    i  = * ( long * ) &y;
    i  = 0x5f3759df - ( i >> 1 );
    y  = * ( f32 * ) &i;
    y  = y * ( f - ( x * y * y ) );
    y  = y * ( f - ( x * y * y ) );
    return number * y;
    */
    number.sqrt()
}

/*
===============================================================================

POSITION TESTING

===============================================================================
*/

fn CM_TestBoxInBrush(cm: &clipMap_t, tw: &mut traceWork_t, brush: &cbrush_t) {
    if brush.numsides == 0 {
        return;
    }

    // special test for axial
    if tw.bounds[0][0] > brush.bounds.maxs[0]
        || tw.bounds[0][1] > brush.bounds.maxs[1]
        || tw.bounds[0][2] > brush.bounds.maxs[2]
        || tw.bounds[1][0] < brush.bounds.mins[0]
        || tw.bounds[1][1] < brush.bounds.mins[1]
        || tw.bounds[1][2] < brush.bounds.mins[2]
    {
        return;
    }

    if tw.sphere.use_ {
        // the first six planes are the axial planes, so we only
        // need to test the remainder
        let brush_sides = brush.sides(cm);
        for side in brush_sides[6..].iter() {
            let plane = side.plane(cm);

            // adjust the plane distance apropriately for radius
            let dist = plane.dist + tw.sphere.radius;
            // find the closest point on the capsule to the plane
            let t = plane.normal.dot(tw.sphere.offset);
            let startp = if t > 0.0 {
                tw.start - tw.sphere.offset
            } else {
                tw.start + tw.sphere.offset
            };
            let d1 = startp.dot(plane.normal) - dist;
            // if completely in front of face, no intersection
            if d1 > 0.0 {
                return;
            }
        }
    } else {
        // the first six planes are the axial planes, so we only
        // need to test the remainder
        for side in brush.sides(cm)[6..].iter() {
            let plane = side.plane(cm);

            // adjust the plane distance apropriately for mins/maxs
            let dist = plane.dist - tw.offsets[plane.signbits as usize].dot(plane.normal);
            let d1 = tw.start.dot(plane.normal) - dist;

            // if completely in front of face, no intersection
            if d1 > 0.0 {
                return;
            }
        }
    }

    // inside this brush
    tw.trace.startsolid = true;
    tw.trace.allsolid = true;
    tw.trace.fraction = 0.0;
    tw.trace.contents = brush.contents;
}

fn CM_TestInLeaf(cm: &mut clipMap_t, tw: &mut traceWork_t, leaf: &cLeaf_t) {
    // test box position against all brushes in the leaf
    for k in 0..leaf.numLeafBrushes as usize {
        let brushnum = cm.leafbrushes[leaf.firstLeafBrush as usize + k] as usize;

        let b = &mut cm.brushes[brushnum];
        if b.checkcount == cm.checkcount {
            continue; // already checked this brush in another leaf
        }
        b.checkcount = cm.checkcount;

        if (b.contents & tw.contents) == 0 {
            continue;
        }

        // TODO: Eliminate this copy
        let mut b_copy = b.clone();
        CM_TestBoxInBrush(cm, tw, &mut b_copy);
        cm.brushes[brushnum] = b_copy;

        if tw.trace.allsolid {
            return;
        }
    }

    // test against all patches
    // #ifdef BSPC
    //     if (1) {
    // #else
    if !cm_noCurves.as_bool() {
        // #endif //BSPC
        let cm_checkcount = cm.checkcount;
        for opt_patch in cm.surfaces[leaf.leaf_surfaces_range()].iter_mut() {
            if let Some(ref mut patch) = opt_patch {
                if patch.checkcount == cm_checkcount {
                    continue; // already checked this brush in another leaf
                }
                patch.checkcount = cm_checkcount;

                if (patch.contents & tw.contents) == 0 {
                    continue;
                }

                if CM_PositionTestInPatchCollide(tw, &patch.pc) {
                    tw.trace.startsolid = true;
                    tw.trace.allsolid = true;
                    tw.trace.fraction = 0.0;
                    tw.trace.contents = patch.contents;
                    return;
                }
            }
        }
    }
}

// capsule inside capsule check
fn CM_TestCapsuleInCapsule(cm: &clipMap_t, tw: &mut traceWork_t, model: clipHandle_t) {
    let vec3_bounds { mins, maxs } = CM_ModelBounds(cm, model);

    let mut top = tw.start + tw.sphere.offset;
    let bottom = tw.start - tw.sphere.offset;
    let offset = vec3_t::mid(mins, maxs);

    // let symetricSize0 = mins - offset;
    let symetricSize1 = maxs - offset;

    let halfwidth = symetricSize1[0];
    let halfheight = symetricSize1[2];
    let radius = fmin(halfwidth, halfheight);
    let offs = halfheight - radius;

    let r = Square(tw.sphere.radius + radius);
    let are_near = move |p1: vec3_t, p2: vec3_t| (p1 - p2).length2() < r;

    // check if any of the spheres overlap
    let mut p1 = offset;
    p1[2] += offs;

    if are_near(p1, top) {
        tw.trace.startsolid = true;
        tw.trace.allsolid = true;
        tw.trace.fraction = 0.0;
    }

    if are_near(p1, bottom) {
        tw.trace.startsolid = true;
        tw.trace.allsolid = true;
        tw.trace.fraction = 0.0;
    }

    let mut p2 = offset;
    p2[2] -= offs;
    if are_near(p2, top) {
        tw.trace.startsolid = true;
        tw.trace.allsolid = true;
        tw.trace.fraction = 0.0;
    }
    if are_near(p2, bottom) {
        tw.trace.startsolid = true;
        tw.trace.allsolid = true;
        tw.trace.fraction = 0.0;
    }
    // if between cylinder up and lower bounds
    if (top[2] >= p1[2] && top[2] <= p2[2]) || (bottom[2] >= p1[2] && bottom[2] <= p2[2]) {
        // 2d coordinates
        top[2] = 0.0;
        p1[2] = 0.0;
        // if the cylinders overlap
        if are_near(top, p1) {
            tw.trace.startsolid = true;
            tw.trace.allsolid = true;
            tw.trace.fraction = 0.0;
        }
    }
}

/// bounding box inside capsule check
fn CM_TestBoundingBoxInCapsule(cm: &mut clipMap_t, tw: &mut traceWork_t, model: clipHandle_t) {
    // mins maxs of the capsule
    let vec3_bounds { mins, maxs } = CM_ModelBounds(cm, model);

    // offset for capsule center
    let offset = vec3_t::mid(mins, maxs);
    // let size0 = mins - offset;
    let size1 = maxs - offset;
    tw.start -= offset;
    tw.end -= offset;

    // replace the bounding box with the capsule
    tw.sphere.use_ = true;
    tw.sphere.radius = fmin(size1[0], size1[2]);
    tw.sphere.halfheight = size1[2];
    tw.sphere.offset = VectorSet(0.0, 0.0, size1[2] - tw.sphere.radius);

    // replace the capsule with the bounding box
    let h = CM_TempBoxModel(cm, tw.size[0], tw.size[1], false);
    // calculate collision
    let cmod = CM_ClipHandleToModel(cm, h);

    // TODO: eliminate this copy
    let cmod_leaf = cmod.leaf.clone();

    CM_TestInLeaf(cm, tw, &cmod_leaf);
}

const MAX_POSITION_LEAFS: usize = 1024;

fn CM_PositionTest(cm: &mut clipMap_t, tw: &mut traceWork_t) {
    let mut leafs = [0i32; MAX_POSITION_LEAFS];

    // identify the leafs we are touching
    let mut ll = leafList_t {
        bounds: vec3_bounds {
            mins: tw.size[0] + tw.start - vec3_t::from_scalar(1.0),
            maxs: tw.size[1] + tw.start + vec3_t::from_scalar(1.0),
        },
        count: 0,
        lastLeaf: 0,
        overflowed: false,
    };

    cm.checkcount += 1;

    CM_BoxLeafnums_r(
        cm,
        &mut ll,
        0,
        &mut |cm, ll: &mut leafList_t, nodenum: i32| {
            // was CM_StoreLeafs
            let leafNum = -1 - nodenum;

            // store the lastLeaf even if the list is overflowed
            if cm.leafs[leafNum as usize].cluster != -1 {
                ll.lastLeaf = leafNum;
            }
            if ll.count < leafs.len() {
                leafs[ll.count] = leafNum;
                ll.count += 1;
            } else {
                ll.overflowed = true;
            }
        },
    );

    cm.checkcount += 1;

    // test the contents of the leafs
    for &leaf_num in leafs[..ll.count].iter() {
        // TODO: eliminate this copy
        let leaf = cm.leafs[leaf_num as usize].clone();
        CM_TestInLeaf(cm, tw, &leaf);
        if tw.trace.allsolid {
            break;
        }
    }
}

/*
===============================================================================

TRACING

===============================================================================
*/

fn CM_TraceThroughPatch(tw: &mut traceWork_t, patch: &cPatch_t) {
    // c_patch_traces.add_one();

    let oldFrac = tw.trace.fraction;

    CM_TraceThroughPatchCollide(tw, &patch.pc);

    if tw.trace.fraction < oldFrac {
        tw.trace.surfaceFlags = patch.surfaceFlags;
        tw.trace.contents = patch.contents;
    }
}

fn CM_TraceThroughBrush(cm: &clipMap_t, tw: &mut traceWork_t, brush: &cbrush_t) {
    let mut enterFrac: f32 = -1.0;
    let mut leaveFrac: f32 = 1.0;

    if brush.numsides == 0 {
        return;
    }

    let mut clipplane: Option<&cplane_t> = None;

    // c_brush_traces++;

    let mut getout = false;
    let mut startout = false;

    let mut leadside: Option<&cbrushside_t> = None;

    if tw.sphere.use_ {
        //
        // compare the trace against all planes of the brush
        // find the latest time the trace crosses a plane towards the interior
        // and the earliest time the trace crosses a plane towards the exterior
        //
        for side in brush.sides(cm).iter() {
            let plane = &cm.planes[side.plane_num as usize];

            // adjust the plane distance apropriately for radius
            let dist = plane.dist + tw.sphere.radius;

            // find the closest point on the capsule to the plane
            let t = plane.normal.dot(tw.sphere.offset);
            let startp;
            let endp;
            if t > 0.0 {
                startp = tw.start - tw.sphere.offset;
                endp = tw.end - tw.sphere.offset;
            } else {
                startp = tw.start + tw.sphere.offset;
                endp = tw.end + tw.sphere.offset;
            }

            let d1 = startp.dot(plane.normal) - dist;
            let d2 = endp.dot(plane.normal) - dist;

            if d2 > 0.0 {
                getout = true; // endpoint is not in solid
            }
            if d1 > 0.0 {
                startout = true;
            }

            // if completely in front of face, no intersection with the entire brush
            if d1 > 0.0 && (d2 >= SURFACE_CLIP_EPSILON || d2 >= d1) {
                return;
            }

            // if it doesn't cross the plane, the plane isn't relevent
            if d1 <= 0.0 && d2 <= 0.0 {
                continue;
            }

            // crosses face
            if d1 > d2 {
                // enter
                let mut f = (d1 - SURFACE_CLIP_EPSILON) / (d1 - d2);
                if f < 0.0 {
                    f = 0.0;
                }
                if f > enterFrac {
                    enterFrac = f;
                    clipplane = Some(plane);
                    leadside = Some(side);
                }
            } else {
                // leave
                let mut f = (d1 + SURFACE_CLIP_EPSILON) / (d1 - d2);
                if f > 1.0 {
                    f = 1.0;
                }
                if f < leaveFrac {
                    leaveFrac = f;
                }
            }
        }
    } else {
        //
        // compare the trace against all planes of the brush
        // find the latest time the trace crosses a plane towards the interior
        // and the earliest time the trace crosses a plane towards the exterior
        //
        for side in brush.sides(cm).iter() {
            let plane = side.plane(cm);

            // adjust the plane distance apropriately for mins/maxs
            let dist = plane.dist - DotProduct(tw.offsets[plane.signbits as usize], plane.normal);

            let d1 = DotProduct(tw.start, plane.normal) - dist;
            let d2 = DotProduct(tw.end, plane.normal) - dist;

            if d2 > 0.0 {
                getout = true; // endpoint is not in solid
            }
            if d1 > 0.0 {
                startout = true;
            }

            // if completely in front of face, no intersection with the entire brush
            if d1 > 0.0 && (d2 >= SURFACE_CLIP_EPSILON || d2 >= d1) {
                return;
            }

            // if it doesn't cross the plane, the plane isn't relevent
            if d1 <= 0.0 && d2 <= 0.0 {
                continue;
            }

            // crosses face
            if d1 > d2 {
                // enter
                let mut f = (d1 - SURFACE_CLIP_EPSILON) / (d1 - d2);
                if f < 0.0 {
                    f = 0.0;
                }
                if f > enterFrac {
                    enterFrac = f;
                    clipplane = Some(plane);
                    leadside = Some(side);
                }
            } else {
                // leave
                let mut f = (d1 + SURFACE_CLIP_EPSILON) / (d1 - d2);
                if f > 1.0 {
                    f = 1.0;
                }
                if f < leaveFrac {
                    leaveFrac = f;
                }
            }
        }
    }

    //
    // all planes have been checked, and the trace was not
    // completely outside the brush
    //
    if !startout {
        // original point was inside brush
        tw.trace.startsolid = true;
        if !getout {
            tw.trace.allsolid = true;
            tw.trace.fraction = 0.0;
            tw.trace.contents = brush.contents;
        }
        return;
    }

    if enterFrac < leaveFrac {
        if enterFrac > -1.0 && enterFrac < tw.trace.fraction {
            if enterFrac < 0.0 {
                enterFrac = 0.0;
            }
            tw.trace.fraction = enterFrac;
            tw.trace.plane = *clipplane.unwrap();
            tw.trace.surfaceFlags = leadside.unwrap().surfaceFlags;
            tw.trace.contents = brush.contents;
        }
    }
}

fn CM_TraceThroughLeaf(cm: &mut clipMap_t, tw: &mut traceWork_t, leaf: &cLeaf_t) {
    // trace line against all brushes in the leaf

    for i in leaf.leaf_brushes_range() {
        let brushnum = cm.leafbrushes[i];
        let b = &mut cm.brushes[brushnum as usize];
        if b.checkcount == cm.checkcount {
            continue; // already checked this brush in another leaf
        }
        b.checkcount = cm.checkcount;

        if (b.contents & tw.contents) == 0 {
            continue;
        }

        // TODO: eliminate this copy
        let b_copy = b.clone();
        CM_TraceThroughBrush(cm, tw, &b_copy);
        if tw.trace.fraction == 0.0 {
            return;
        }
    }

    // trace line against all patches in the leaf
    // #ifdef BSPC
    //     if (1) {
    // #else
    if true
    /* !cm_noCurves.integer */
    {
        //#endif
        for opt_patch in leaf.leaf_surfaces_mut(&mut cm.surfaces).iter_mut() {
            if let Some(ref mut patch) = opt_patch {
                if patch.checkcount == cm.checkcount {
                    continue; // already checked this patch in another leaf
                }
                patch.checkcount = cm.checkcount;

                if (patch.contents & tw.contents) == 0 {
                    continue;
                }

                CM_TraceThroughPatch(tw, patch);
                if tw.trace.fraction == 0.0 {
                    return;
                }
            }
        }
    }
}

const RADIUS_EPSILON: f32 = 1.0;

/// get the first intersection of the ray with the sphere
fn CM_TraceThroughSphere(
    tw: &mut traceWork_t,
    origin: vec3_t,
    radius: f32,
    start: vec3_t,
    end: vec3_t,
) {
    // if inside the sphere
    let dir = start - origin;
    let l1 = dir.length2();
    if l1 < Square(radius) {
        tw.trace.fraction = 0.0;
        tw.trace.startsolid = true;
        // test for allsolid
        let dir = end - origin;
        if dir.length2() < Square(radius) {
            tw.trace.allsolid = true;
        }
        return;
    }
    //
    let mut dir = end - start;
    let length = VectorNormalize_mut(&mut dir);
    //
    let l1 = CM_DistanceFromLineSquared(origin, start, end, dir);
    let v1 = end - origin;
    let l2 = v1.length2();
    // if no intersection with the sphere and the end point is at least an epsilon away
    if l1 >= Square(radius) && l2 > Square(radius + SURFACE_CLIP_EPSILON) {
        return;
    }

    //
    //  | origin - (start + t * dir) | = radius
    //  a = dir[0]^2 + dir[1]^2 + dir[2]^2;
    //  b = 2 * (dir[0] * (start[0] - origin[0]) + dir[1] * (start[1] - origin[1]) + dir[2] * (start[2] - origin[2]));
    //  c = (start[0] - origin[0])^2 + (start[1] - origin[1])^2 + (start[2] - origin[2])^2 - radius^2;
    //
    let v1 = start - origin;
    // dir is normalized so a = 1
    // let a = 1.0;//dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2];
    let b = 2.0 * (dir[0] * v1[0] + dir[1] * v1[1] + dir[2] * v1[2]);
    let c = v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]
        - (radius + RADIUS_EPSILON) * (radius + RADIUS_EPSILON);

    let d = b * b - 4.0 * c; // * a;
    if d > 0.0 {
        let sqrtd = SquareRootFloat(d);
        // = (- b + sqrtd) * 0.5f; // / (2.0f * a);
        let mut fraction = (-b - sqrtd) * 0.5; // / (2.0f * a);
                                               //
        if fraction < 0.0 {
            fraction = 0.0;
        } else {
            fraction /= length;
        }
        if fraction < tw.trace.fraction {
            tw.trace.fraction = fraction;
            let intersection = VectorMA(start, fraction, end - start);
            let dir = intersection - origin;
            let dir = dir * (1.0 / (radius + RADIUS_EPSILON));
            tw.trace.plane.normal = dir;
            tw.trace.plane.dist = dir.dot(tw.modelOrigin + intersection);
            tw.trace.contents = CONTENTS_BODY;
        }
    } else if d == 0.0 {
        //t1 = (- b ) / 2;
        // slide along the sphere
    }
    // no intersection at all
}

/// get the first intersection of the ray with the cylinder
/// the cylinder extends halfheight above and below the origin
fn CM_TraceThroughVerticalCylinder(
    tw: &mut traceWork_t,
    origin: vec3_t,
    radius: f32,
    halfheight: f32,
    start: vec3_t,
    end: vec3_t,
) {
    // 2d coordinates
    let start2d = VectorSet(start[0], start[1], 0.0);
    let end2d = VectorSet(end[0], end[1], 0.0);
    let org2d = VectorSet(origin[0], origin[1], 0.0);
    // if between lower and upper cylinder bounds
    if start[2] <= origin[2] + halfheight && start[2] >= origin[2] - halfheight {
        // if inside the cylinder
        if (start2d - org2d).length2() < Square(radius) {
            tw.trace.fraction = 0.0;
            tw.trace.startsolid = true;
            if (end2d - org2d).length2() < Square(radius) {
                tw.trace.allsolid = true;
            }
            return;
        }
    }
    //
    let mut dir = end2d - start2d;
    let length = VectorNormalize_mut(&mut dir);
    //
    // if no intersection with the cylinder and the end point is at least an epsilon away
    if CM_DistanceFromLineSquared(org2d, start2d, end2d, dir) >= Square(radius)
        && (end2d - org2d).length2() > Square(radius + SURFACE_CLIP_EPSILON)
    {
        return;
    }
    //
    //
    // (start[0] - origin[0] - t * dir[0]) ^ 2 + (start[1] - origin[1] - t * dir[1]) ^ 2 = radius ^ 2
    // (v1[0] + t * dir[0]) ^ 2 + (v1[1] + t * dir[1]) ^ 2 = radius ^ 2;
    // v1[0] ^ 2 + 2 * v1[0] * t * dir[0] + (t * dir[0]) ^ 2 +
    //                      v1[1] ^ 2 + 2 * v1[1] * t * dir[1] + (t * dir[1]) ^ 2 = radius ^ 2
    // t ^ 2 * (dir[0] ^ 2 + dir[1] ^ 2) + t * (2 * v1[0] * dir[0] + 2 * v1[1] * dir[1]) +
    //                      v1[0] ^ 2 + v1[1] ^ 2 - radius ^ 2 = 0
    //
    let v1 = start - origin;
    // dir is normalized so we can use a = 1
    // let a = 1.0;// * (dir[0] * dir[0] + dir[1] * dir[1]);
    let b = 2.0 * (v1[0] * dir[0] + v1[1] * dir[1]);
    let c = v1[0] * v1[0] + v1[1] * v1[1] - Square(radius + RADIUS_EPSILON);
    let d = b * b - 4.0 * c; // * a;
    if d > 0.0 {
        let sqrtd = SquareRootFloat(d);
        // = (- b + sqrtd) * 0.5f;// / (2.0f * a);
        let fraction = (-b - sqrtd) * 0.5; // / (2.0f * a);
                                           //
        let fraction = if fraction < 0.0 {
            0.0
        } else {
            fraction / length
        };
        if fraction < tw.trace.fraction {
            let dir = end - start;
            let intersection = VectorMA(start, fraction, dir);
            // if the intersection is between the cylinder lower and upper bound
            if intersection[2] <= origin[2] + halfheight
                && intersection[2] >= origin[2] - halfheight
            {
                //
                tw.trace.fraction = fraction;
                let mut dir = intersection - origin;
                dir[2] = 0.0;
                /*#ifdef CAPSULE_DEBUG
                    l2 = VectorLength(dir);
                    if (l2 <= radius) {
                        int bah = 1;
                    }
                #endif*/
                let scale = 1.0 / (radius + RADIUS_EPSILON);
                let dir = dir * scale;
                tw.trace.plane.normal = dir;
                tw.trace.plane.dist = dir.dot(tw.modelOrigin + intersection);
                tw.trace.contents = CONTENTS_BODY;
            }
        }
    } else if d == 0.0 {
        //t[0] = (- b ) / 2 * a;
        // slide along the cylinder
    }
    // no intersection at all
}

/// capsule vs. capsule collision (not rotated)
fn CM_TraceCapsuleThroughCapsule(cm: &clipMap_t, tw: &mut traceWork_t, model: clipHandle_t) {
    let vec3_bounds { mins, maxs } = CM_ModelBounds(cm, model);
    // test trace bounds vs. capsule bounds
    if tw.bounds[0][0] > maxs[0] + RADIUS_EPSILON
        || tw.bounds[0][1] > maxs[1] + RADIUS_EPSILON
        || tw.bounds[0][2] > maxs[2] + RADIUS_EPSILON
        || tw.bounds[1][0] < mins[0] - RADIUS_EPSILON
        || tw.bounds[1][1] < mins[1] - RADIUS_EPSILON
        || tw.bounds[1][2] < mins[2] - RADIUS_EPSILON
    {
        return;
    }
    // top origin and bottom origin of each sphere at start and end of trace
    let starttop = tw.start + tw.sphere.offset;
    let startbottom = tw.start - tw.sphere.offset;
    let endtop = tw.end + tw.sphere.offset;
    let endbottom = tw.end - tw.sphere.offset;

    // calculate top and bottom of the capsule spheres to collide with
    let offset = (mins + maxs) * 0.5;
    // let symetricSize0 = mins - offset;
    let symetricSize1 = maxs - offset;

    let halfwidth = symetricSize1[0];
    let halfheight = symetricSize1[2];
    let mut radius = if halfwidth > halfheight {
        halfheight
    } else {
        halfwidth
    };
    let offs = halfheight - radius;
    let mut top = offset;
    top[2] += offs;
    let mut bottom = offset;
    bottom[2] -= offs;
    // expand radius of spheres
    radius += tw.sphere.radius;
    // if there is horizontal movement
    if tw.start[0] != tw.end[0] || tw.start[1] != tw.end[1] {
        // height of the expanded cylinder is the height of both cylinders minus the radius of both spheres
        let h = halfheight + tw.sphere.halfheight - radius;
        // if the cylinder has a height
        if h > 0.0 {
            // test for collisions between the cylinders
            CM_TraceThroughVerticalCylinder(tw, offset, radius, h, tw.start, tw.end);
        }
    }
    // test for collision between the spheres
    CM_TraceThroughSphere(tw, top, radius, startbottom, endbottom);
    CM_TraceThroughSphere(tw, bottom, radius, starttop, endtop);
}

/*
================
CM_TraceBoundingBoxThroughCapsule

bounding box vs. capsule collision
================
*/
fn CM_TraceBoundingBoxThroughCapsule(
    cm: &mut clipMap_t,
    tw: &mut traceWork_t,
    model: clipHandle_t,
) {
    // mins maxs of the capsule
    let vec3_bounds { mins, maxs } = CM_ModelBounds(cm, model);

    // offset for capsule center
    let offset = (mins + maxs) * 0.5;
    // let size0 = mins - offset;
    let size1 = maxs - offset;
    tw.start -= offset;
    tw.end -= offset;

    // replace the bounding box with the capsule
    tw.sphere.use_ = true;
    tw.sphere.radius = if size1[0] > size1[2] {
        size1[2]
    } else {
        size1[0]
    };
    tw.sphere.halfheight = size1[2];
    tw.sphere.offset = VectorSet(0.0, 0.0, size1[2] - tw.sphere.radius);

    // replace the capsule with the bounding box
    let h = CM_TempBoxModel(cm, tw.size[0], tw.size[1], false);
    // calculate collision
    let cmod = CM_ClipHandleToModel(cm, h);

    let cmod_leaf = cmod.leaf.clone(); // TODO: remove copy
    CM_TraceThroughLeaf(cm, tw, &cmod_leaf);
}

/// Traverse all the contacted leafs from the start to the end position.
/// If the trace is a point, they will be exactly in order, but for larger
/// trace volumes it is possible to hit something in a later leaf with
/// a smaller intercept fraction.
fn CM_TraceThroughTree(
    cm: &mut clipMap_t,
    tw: &mut traceWork_t,
    num: i32,
    p1f: f32,
    p2f: f32,
    p1: vec3_t,
    p2: vec3_t,
) {
    if tw.trace.fraction <= p1f {
        return; // already hit something nearer
    }

    // if < 0, we are in a leaf node
    if num < 0 {
        let leaf_num = -1 - num;
        let leaf_copy = cm.leafs[leaf_num as usize].clone(); // TODO: remove copy
        CM_TraceThroughLeaf(cm, tw, &leaf_copy);
        return;
    }

    // find the point distances to the seperating plane
    // and the offset for the size of the box
    let node = &cm.nodes[num as usize];
    let plane = node.plane(cm);

    // adjust the plane distance apropriately for mins/maxs
    let t1;
    let t2;
    let offset;
    if plane.type_ < 3 {
        t1 = p1[plane.type_ as usize] - plane.dist;
        t2 = p2[plane.type_ as usize] - plane.dist;
        offset = tw.extents[plane.type_ as usize];
    } else {
        t1 = plane.distance_to(p1);
        t2 = plane.distance_to(p2);
        if tw.isPoint {
            offset = 0.0;
        } else {
            /*#if 0 // bk010201 - DEAD
                        // an axial brush right behind a slanted bsp plane
                        // will poke through when expanded, so adjust
                        // by sqrt(3)
                        offset = fabs(tw.extents[0]*plane.normal[0]) +
                            fabs(tw.extents[1]*plane.normal[1]) +
                            fabs(tw.extents[2]*plane.normal[2]);

                        offset *= 2;
                        offset = tw.maxOffset;
            #endif*/
            // this is silly
            offset = 2048.0;
        }
    }

    // see which sides we need to consider
    if t1 >= offset + 1.0 && t2 >= offset + 1.0 {
        CM_TraceThroughTree(cm, tw, node.children[0], p1f, p2f, p1, p2);
        return;
    }
    if t1 < -offset - 1.0 && t2 < -offset - 1.0 {
        CM_TraceThroughTree(cm, tw, node.children[1], p1f, p2f, p1, p2);
        return;
    }

    // put the crosspoint SURFACE_CLIP_EPSILON pixels on the near side
    let mut frac: f32;
    let mut frac2: f32;
    let this_child: i32;
    let that_child: i32;
    if t1 < t2 {
        let idist = 1.0 / (t1 - t2);
        this_child = node.children[1];
        that_child = node.children[0];
        frac2 = (t1 + offset + SURFACE_CLIP_EPSILON) * idist;
        frac = (t1 - offset + SURFACE_CLIP_EPSILON) * idist;
    } else if t1 > t2 {
        let idist = 1.0 / (t1 - t2);
        this_child = node.children[0];
        that_child = node.children[1];
        frac2 = (t1 - offset - SURFACE_CLIP_EPSILON) * idist;
        frac = (t1 + offset + SURFACE_CLIP_EPSILON) * idist;
    } else {
        this_child = node.children[0];
        that_child = node.children[1];
        frac = 1.0;
        frac2 = 0.0;
    }

    // move up to the node
    if frac < 0.0 {
        frac = 0.0;
    }
    if frac > 1.0 {
        frac = 1.0;
    }
    let midf = p1f + (p2f - p1f) * frac;
    let mid = p1 + frac * (p2 - p1);
    CM_TraceThroughTree(cm, tw, this_child, p1f, midf, p1, mid);

    // go past the node
    if frac2 < 0.0 {
        frac2 = 0.0;
    }
    if frac2 > 1.0 {
        frac2 = 1.0;
    }
    let midf = p1f + (p2f - p1f) * frac2;
    let mid = p1 + frac2 * (p2 - p1);
    CM_TraceThroughTree(cm, tw, that_child, midf, p2f, mid, p2);
}

//======================================================================

fn CM_Trace(
    cm: &mut clipMap_t,
    start: vec3_t,
    end: vec3_t,
    mins: Option<vec3_t>,
    maxs: Option<vec3_t>,
    model: clipHandle_t,
    origin: vec3_t,
    brushmask: i32,
    capsule: i32,
    sphere: Option<&sphere_t>,
) -> trace_t {
    let cmod = CM_ClipHandleToModel(cm, model);
    let cmod_leaf = cmod.leaf.clone(); // TODO: remove copy

    cm.checkcount += 1; // for multi-check avoidance

    // c_traces.add_one();             // for statistics, may be zeroed

    /*
    if cm.nodes.is_empty()) {
        // map not loaded, shouldn't happen
        return tw.trace;
    }
    */
    assert!(!cm.nodes.is_empty());

    // allow NULL to be passed in for 0,0,0
    let mins = mins.unwrap_or(vec3_t::ORIGIN);
    let maxs = maxs.unwrap_or(vec3_t::ORIGIN);

    // adjust mins, maxs, size[0], size[1] so that mins and maxs are always symetric, which
    // avoids some complications with plane expanding of rotated
    // bmodels
    let offset = (mins + maxs) * 0.5;
    let tw_size = [mins - offset, maxs - offset];

    let tw_sphere = if let Some(s) = sphere {
        // if a sphere is already specified
        s.clone()
    } else {
        let radius = fmin(tw_size[1][0], tw_size[1][2]);
        sphere_t {
            use_: capsule != 0,
            radius: radius,
            halfheight: tw_size[1][2],
            offset: VectorSet(0.0, 0.0, tw_size[1][2] - radius),
        }
    };

    let tw_start = start + offset;
    let tw_end = end + offset;

    // calculate bounds
    let tw_bounds = if tw_sphere.use_ {
        [
            vec3_t::min(tw_start, tw_end) - tw_sphere.offset.map(|c| c.abs() - tw_sphere.radius),
            vec3_t::max(tw_start, tw_end) + tw_sphere.offset.map(|c| c.abs() + tw_sphere.radius),
        ]
    } else {
        [
            vec3_t::min(tw_start, tw_end) + tw_size[0],
            vec3_t::max(tw_start, tw_end) + tw_size[1],
        ]
    };

    // fill in a default trace
    let mut tw = traceWork_t {
        trace: trace_t {
            // assume it goes the entire distance until shown otherwise
            fraction: 1.0,
            ..Default::default()
        },
        modelOrigin: origin,

        // set basic parms
        contents: brushmask,

        size: tw_size,
        start: tw_start,
        end: tw_end,

        sphere: tw_sphere,

        maxOffset: tw_size[1][0] + tw_size[1][1] + tw_size[1][2],

        // tw.offsets[signbits] = vector to apropriate corner from origin
        offsets: [
            v3(tw_size[0][0], tw_size[0][1], tw_size[0][2]),
            v3(tw_size[1][0], tw_size[0][1], tw_size[0][2]),
            v3(tw_size[0][0], tw_size[1][1], tw_size[0][2]),
            v3(tw_size[1][0], tw_size[1][1], tw_size[0][2]),
            v3(tw_size[0][0], tw_size[0][1], tw_size[1][2]),
            v3(tw_size[1][0], tw_size[0][1], tw_size[1][2]),
            v3(tw_size[0][0], tw_size[1][1], tw_size[1][2]),
            v3(tw_size[1][0], tw_size[1][1], tw_size[1][2]),
        ],

        bounds: tw_bounds,

        isPoint: false,
        extents: Default::default(),
    };

    //
    // check for position test special case
    //
    if start[0] == end[0] && start[1] == end[1] && start[2] == end[2] {
        if model != 0 {
            /*
            #ifdef ALWAYS_BBOX_VS_BBOX // bk010201 - FIXME - compile time flag?
                        if ( model == BOX_MODEL_HANDLE || model == CAPSULE_MODEL_HANDLE) {
                            tw.sphere.use = false;
                            CM_TestInLeaf( &tw, &cmod->leaf );
                        }
                        else
            #elif defined(ALWAYS_CAPSULE_VS_CAPSULE)
                        if ( model == BOX_MODEL_HANDLE || model == CAPSULE_MODEL_HANDLE) {
                            CM_TestCapsuleInCapsule( &tw, model );
                        }
                        else
            #endif
            */
            if model == CAPSULE_MODEL_HANDLE {
                if tw.sphere.use_ {
                    CM_TestCapsuleInCapsule(cm, &mut tw, model);
                } else {
                    CM_TestBoundingBoxInCapsule(cm, &mut tw, model);
                }
            } else {
                CM_TestInLeaf(cm, &mut tw, &cmod_leaf);
            }
        } else {
            CM_PositionTest(cm, &mut tw);
        }
    } else {
        //
        // check for point special case
        //
        if tw.size[0][0] == 0.0 && tw.size[0][1] == 0.0 && tw.size[0][2] == 0.0 {
            tw.isPoint = true;
            tw.extents = vec3_t::ORIGIN;
        } else {
            tw.isPoint = false;
            tw.extents = tw.size[1];
        }

        //
        // general sweeping through world
        //
        if model != 0 {
            /*
            #ifdef ALWAYS_BBOX_VS_BBOX
                        if ( model == BOX_MODEL_HANDLE || model == CAPSULE_MODEL_HANDLE) {
                            tw.sphere.use = false;
                            CM_TraceThroughLeaf( &tw, &cmod->leaf );
                        }
                        else
            #elif defined(ALWAYS_CAPSULE_VS_CAPSULE)
                        if ( model == BOX_MODEL_HANDLE || model == CAPSULE_MODEL_HANDLE) {
                            CM_TraceCapsuleThroughCapsule( &tw, model );
                        }
                        else
            #endif
            */
            if model == CAPSULE_MODEL_HANDLE {
                if tw.sphere.use_ {
                    CM_TraceCapsuleThroughCapsule(cm, &mut tw, model);
                } else {
                    CM_TraceBoundingBoxThroughCapsule(cm, &mut tw, model);
                }
            } else {
                CM_TraceThroughLeaf(cm, &mut tw, &cmod_leaf);
            }
        } else {
            let tw_start = tw.start;
            let tw_end = tw.end;
            CM_TraceThroughTree(cm, &mut tw, 0, 0.0, 1.0, tw_start, tw_end);
        }
    }

    // generate endpos from the original, unmodified start/end
    tw.trace.endpos = if tw.trace.fraction == 1.0 {
        end
    } else {
        start + tw.trace.fraction * (end - start)
    };

    // If allsolid is set (was entirely inside something solid), the plane is not valid.
    // If fraction == 1.0, we never hit anything, and thus the plane is not valid.
    // Otherwise, the normal on the plane should have unit length
    assert!(
        tw.trace.allsolid || tw.trace.fraction == 1.0 || (tw.trace.plane.normal).length2() > 0.9999
    );
    return tw.trace;
}

pub fn CM_BoxTrace(
    cm: &mut clipMap_t,
    start: vec3_t,
    end: vec3_t,
    mins: Option<vec3_t>,
    maxs: Option<vec3_t>,
    model: clipHandle_t,
    brushmask: i32,
    capsule: i32,
) -> trace_t {
    CM_Trace(
        cm,
        start,
        end,
        mins,
        maxs,
        model,
        vec3_t::ORIGIN,
        brushmask,
        capsule,
        None,
    )
}

/// Handles offseting and rotation of the end points for moving and
/// rotating entities
pub fn CM_TransformedBoxTrace(
    cm: &mut clipMap_t,
    start: vec3_t,
    end: vec3_t,
    mins: Option<vec3_t>,
    maxs: Option<vec3_t>,
    model: clipHandle_t,
    brushmask: i32,
    origin: vec3_t,
    angles: vec3_t,
    capsule: i32,
) -> trace_t {
    let mins = mins.unwrap_or(vec3_t::ORIGIN);
    let maxs = maxs.unwrap_or(vec3_t::ORIGIN);

    // adjust so that mins and maxs are always symetric, which
    // avoids some complications with plane expanding of rotated
    // bmodels
    let offset = (mins + maxs) * 0.5;
    let symetricSize0 = mins - offset;
    let symetricSize1 = maxs - offset;
    let start_l = start + offset;
    let end_l = end + offset;

    // subtract origin offset
    let mut start_l = start_l - origin;
    let mut end_l = end_l - origin;

    // rotate start and end into the models frame of reference
    let rotated =
        model != BOX_MODEL_HANDLE && (angles[0] != 0.0 || angles[1] != 0.0 || angles[2] != 0.0);

    let halfwidth = symetricSize1[0];
    let halfheight = symetricSize1[2];

    let sphere_radius = fmin(halfwidth, halfheight);
    let t = halfheight - sphere_radius;

    let matrix;

    let sphere_offset;
    if rotated {
        // rotation on trace line (start-end) instead of rotating the bmodel
        // NOTE: This is still incorrect for bounding boxes because the actual bounding
        //       box that is swept through the model is not rotated. We cannot rotate
        //       the bounding box or the bmodel because that would make all the brush
        //       bevels invalid.
        //       However this is correct for capsules since a capsule itself is rotated too.
        matrix = CreateRotationMatrix(angles);
        RotatePoint_mut(&mut start_l, &matrix);
        RotatePoint_mut(&mut end_l, &matrix);
        // rotated sphere offset for capsule
        sphere_offset = v3(matrix[0][2] * t, -matrix[1][2] * t, matrix[2][2] * t);
    } else {
        sphere_offset = v3(0.0, 0.0, t);
        matrix = Default::default();
    }

    let sphere = sphere_t {
        use_: capsule != 0,
        radius: sphere_radius,
        halfheight: halfheight,
        offset: sphere_offset,
    };

    // sweep the box through the model
    let mut trace = CM_Trace(
        cm,
        start_l,
        end_l,
        Some(symetricSize0),
        Some(symetricSize1),
        model,
        origin,
        brushmask,
        capsule,
        Some(&sphere),
    );

    // if the bmodel was rotated and there was a collision
    if rotated && trace.fraction != 1.0 {
        // rotation of bmodel collision plane
        let transpose = TransposeMatrix(&matrix);
        RotatePoint_mut(&mut trace.plane.normal, &transpose);
    }

    // re-calculate the end position of the trace because the trace.endpos
    // calculated by CM_Trace could be rotated and have an offset
    trace.endpos = start + trace.fraction * (end - start);
    trace
}
