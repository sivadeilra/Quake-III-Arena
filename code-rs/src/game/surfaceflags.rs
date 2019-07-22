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
//
// This file must be identical in the quake and utils directories

// contents flags are seperate bits
// a given brush can contribute multiple content bits

// these definitions also need to be in q_shared.h!

pub type ContentsFlag = i32;

pub const CONTENTS_SOLID          :ContentsFlag = 1;       // an eye is never valid in a solid
pub const CONTENTS_LAVA           :ContentsFlag = 8;
pub const CONTENTS_SLIME          :ContentsFlag = 16;
pub const CONTENTS_WATER          :ContentsFlag = 32;
pub const CONTENTS_FOG            :ContentsFlag = 64;

pub const CONTENTS_NOTTEAM1       :ContentsFlag = 0x0080;
pub const CONTENTS_NOTTEAM2       :ContentsFlag = 0x0100;
pub const CONTENTS_NOBOTCLIP      :ContentsFlag = 0x0200;

pub const CONTENTS_AREAPORTAL     :ContentsFlag = 0x8000;

pub const CONTENTS_PLAYERCLIP     :ContentsFlag = 0x10000;
pub const CONTENTS_MONSTERCLIP    :ContentsFlag = 0x20000;
//bot specific contents types
pub const CONTENTS_TELEPORTER     :ContentsFlag = 0x40000;
pub const CONTENTS_JUMPPAD        :ContentsFlag = 0x80000;
pub const CONTENTS_CLUSTERPORTAL  :ContentsFlag = 0x100000;
pub const CONTENTS_DONOTENTER     :ContentsFlag = 0x200000;
pub const CONTENTS_BOTCLIP        :ContentsFlag = 0x400000;
pub const CONTENTS_MOVER          :ContentsFlag = 0x800000;

pub const CONTENTS_ORIGIN         :ContentsFlag = 0x1000000;   // removed before bsping an entity

pub const CONTENTS_BODY           :ContentsFlag = 0x2000000;   // should never be on a brush, only in game
pub const CONTENTS_CORPSE         :ContentsFlag = 0x4000000;
pub const CONTENTS_DETAIL         :ContentsFlag = 0x8000000;   // brushes not used for the bsp
pub const CONTENTS_STRUCTURAL     :ContentsFlag = 0x10000000;  // brushes used for the bsp
pub const CONTENTS_TRANSLUCENT    :ContentsFlag = 0x20000000;  // don't consume surface fragments inside
pub const CONTENTS_TRIGGER        :ContentsFlag = 0x40000000;
pub const CONTENTS_NODROP         :ContentsFlag = 0x80000000;  // don't leave bodies or items (death fog, lava)

pub type SurfaceFlag = i32;
pub const SURF_NODAMAGE           :SurfaceFlag = 0x1;     // never give falling damage
pub const SURF_SLICK              :SurfaceFlag = 0x2;     // effects game physics
pub const SURF_SKY                :SurfaceFlag = 0x4;     // lighting from environment map
pub const SURF_LADDER             :SurfaceFlag = 0x8;
pub const SURF_NOIMPACT           :SurfaceFlag = 0x10;    // don't make missile explosions
pub const SURF_NOMARKS            :SurfaceFlag = 0x20;    // don't leave missile marks
pub const SURF_FLESH              :SurfaceFlag = 0x40;    // make flesh sounds and effects
pub const SURF_NODRAW             :SurfaceFlag = 0x80;    // don't generate a drawsurface at all
pub const SURF_HINT               :SurfaceFlag = 0x100;   // make a primary bsp splitter
pub const SURF_SKIP               :SurfaceFlag = 0x200;   // completely ignore, allowing non-closed brushes
pub const SURF_NOLIGHTMAP         :SurfaceFlag = 0x400;   // surface doesn't need a lightmap
pub const SURF_POINTLIGHT         :SurfaceFlag = 0x800;   // generate lighting info at vertexes
pub const SURF_METALSTEPS         :SurfaceFlag = 0x1000;  // clanking footsteps
pub const SURF_NOSTEPS            :SurfaceFlag = 0x2000;  // no footstep sounds
pub const SURF_NONSOLID           :SurfaceFlag = 0x4000;  // don't collide against curves with this set
pub const SURF_LIGHTFILTER        :SurfaceFlag = 0x8000;  // act as a light filter during q3map -light
pub const SURF_ALPHASHADOW        :SurfaceFlag = 0x10000; // do per-pixel light shadow casting in q3map
pub const SURF_NODLIGHT           :SurfaceFlag = 0x20000; // don't dlight even if solid (solid lava, skies)
pub const SURF_DUST               :SurfaceFlag = 0x40000; // leave a dust trail when walking on this surface
