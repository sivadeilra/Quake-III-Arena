use crate::zerocopy::slice_from_bytes;

/* UINT2 defines a two byte word */
type UINT2 = u16;

/* UINT4 defines a four byte word */
type UINT4 = u32;

/* MD4.H - header file for MD4C.C */

/* Copyright (C) 1991-2, RSA Data Security, Inc. Created 1991.

All rights reserved.

License to copy and use this software is granted provided that it is identified as the “RSA Data Security, Inc. MD4 Message-Digest Algorithm” in all material mentioning or referencing this software or this function.
License is also granted to make and use derivative works provided that such works are identified as “derived from the RSA Data Security, Inc. MD4 Message-Digest Algorithm” in all material mentioning or referencing the derived work.
RSA Data Security, Inc. makes no representations concerning either the merchantability of this software or the suitability of this software for any particular purpose. It is provided “as is” without express or implied warranty of any kind.

These notices must be retained in any copies of any part of this documentation and/or software. */

/* MD4 context. */
pub struct MD4_CTX {
    state: [u32; 4],  /* state (ABCD) */
    count: u64,       /* number of bits, modulo 2^64 (lsb first) */
    buffer: [u8; 64], /* input buffer */
}

/* MD4C.C - RSA Data Security, Inc., MD4 message-digest algorithm */
/* Copyright (C) 1990-2, RSA Data Security, Inc. All rights reserved.

License to copy and use this software is granted provided that it is identified as the
RSA Data Security, Inc. MD4 Message-Digest Algorithm
 in all material mentioning or referencing this software or this function.
License is also granted to make and use derivative works provided that such works are identified as
derived from the RSA Data Security, Inc. MD4 Message-Digest Algorithm
in all material mentioning or referencing the derived work.
RSA Data Security, Inc. makes no representations concerning either the merchantability of this software or the suitability of this software for any particular purpose. It is provided
as is without express or implied warranty of any kind.

These notices must be retained in any copies of any part of this documentation and/or software. */

/* Constants for MD4Transform routine.  */
const S11: u32 = 3;
const S12: u32 = 7;
const S13: u32 = 11;
const S14: u32 = 19;
const S21: u32 = 3;
const S22: u32 = 5;
const S23: u32 = 9;
const S24: u32 = 13;
const S31: u32 = 3;
const S32: u32 = 9;
const S33: u32 = 11;
const S34: u32 = 15;

static PADDING: [u8; 64] = [
    0x80, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0,
];

/* F, G and H are basic MD4 functions. */
fn F(x: u32, y: u32, z: u32) -> u32 {
    (x & y) | ((!x) & z)
}
fn G(x: u32, y: u32, z: u32) -> u32 {
    (x & y) | (x & z) | (y & z)
}
fn H(x: u32, y: u32, z: u32) -> u32 {
    x ^ y ^ z
}

/* ROTATE_LEFT rotates x left n bits. */
fn ROTATE_LEFT(x: u32, n: u32) -> u32 {
    (x << n) | (x >> (32 - n))
}

/* FF, GG and HH are transformations for rounds 1, 2 and 3 */
/* Rotation is separate from addition to prevent recomputation */
fn FF(a: &mut u32, b: u32, c: u32, d: u32, x: u32, s: u32) {
    *a += F(b, c, d).wrapping_add(x);
    *a = ROTATE_LEFT(*a, s);
}
fn GG(a: &mut u32, b: u32, c: u32, d: u32, x: u32, s: u32) {
    *a += G(b, c, d).wrapping_add(x).wrapping_add(0x5a827999u32);
    *a = ROTATE_LEFT(*a, s);
}
fn HH(a: &mut u32, b: u32, c: u32, d: u32, x: u32, s: u32) {
    *a += H(b, c, d).wrapping_add(x).wrapping_add(0x6ed9eba1u32);
    *a = ROTATE_LEFT(*a, s);
}

/* MD4 initialization. Begins an MD4 operation, writing a new context. */
pub fn MD4Init() -> MD4_CTX {
    /* Load magic initialization constants.*/
    MD4_CTX {
        count: 0,
        state: [0x67452301, 0xefcdab89, 0x98badcfe, 0x10325476],
        buffer: [0; 64],
    }
}

/* MD4 block update operation. Continues an MD4 message-digest operation, processing another message block, and updating the context. */
pub fn MD4Update(context: &mut MD4_CTX, input: &[u8]) {
    /* Compute number of bytes mod 64 */
    let mut index = ((context.count >> 3) & 0x3F) as usize;

    /* Update number of bits */
    context.count += (input.len() as u64) << 3;

    let partLen = 64 - index;

    let mut i: usize;

    /* Transform as many times as possible.*/
    if input.len() >= partLen {
        context.buffer[index..index + partLen].copy_from_slice(&input[..partLen]);
        MD4Transform(&mut context.state, &context.buffer);
        i = partLen;
        while i + 63 < input.len() {
            MD4Transform(&mut context.state, &input[i..i + partLen]);
            i += 64;
        }
        index = 0;
    } else {
        i = 0;
    }

    /* Buffer remaining input */
    context.buffer[index..index + input.len() - i].copy_from_slice(&input[i..i + input.len() - i]);
}

pub fn MD4Final_u32(context: &mut MD4_CTX) -> [u32; 4] {
    let digest = MD4Final_u8(context);
    let s = slice_from_bytes::<u32>(&digest);
    [s[0], s[1], s[2], s[3]]
}

/* MD4 finalization. Ends an MD4 message-digest operation, writing the the message digest and zeroizing the context. */
pub fn MD4Final_u8(context: &mut MD4_CTX) -> [u8; 16] {
    /* Save number of bits */
    let bits: [u8; 8] = context.count.to_le_bytes();

    /* Pad out to 56 mod 64.*/
    let index = ((context.count >> 3) & 0x3f) as usize;
    let padLen = if index < 56 { 56 - index } else { 120 - index };
    MD4Update(context, &PADDING[..padLen]);

    /* Append length (before padding) */
    MD4Update(context, &bits);

    /* Store state in digest */
    let mut digest = [0u8; 16];
    Encode(&mut digest, &context.state);

    digest
}

/* MD4 basic transformation. Transforms state based on block. */
fn MD4Transform(state: &mut [u32; 4], block: &[u8]) {
    assert_eq!(block.len(), 64);

    let mut a = state[0];
    let mut b = state[1];
    let mut c = state[2];
    let mut d = state[3];

    let mut x = [0u32; 16];
    Decode(&mut x, block);

    /* Round 1 */
    FF(&mut a, b, c, d, x[0], S11); // 1
    FF(&mut d, a, b, c, x[1], S12); // 2
    FF(&mut c, d, a, b, x[2], S13); // 3
    FF(&mut b, c, d, a, x[3], S14); // 4
    FF(&mut a, b, c, d, x[4], S11); // 5
    FF(&mut d, a, b, c, x[5], S12); // 6
    FF(&mut c, d, a, b, x[6], S13); // 7
    FF(&mut b, c, d, a, x[7], S14); // 8
    FF(&mut a, b, c, d, x[8], S11); // 9
    FF(&mut d, a, b, c, x[9], S12); // 10
    FF(&mut c, d, a, b, x[10], S13); // 11
    FF(&mut b, c, d, a, x[11], S14); // 12
    FF(&mut a, b, c, d, x[12], S11); // 13
    FF(&mut d, a, b, c, x[13], S12); // 14
    FF(&mut c, d, a, b, x[14], S13); // 15
    FF(&mut b, c, d, a, x[15], S14); // 16

    /* Round 2 */
    GG(&mut a, b, c, d, x[0], S21); // 17
    GG(&mut d, a, b, c, x[4], S22); // 18
    GG(&mut c, d, a, b, x[8], S23); // 19
    GG(&mut b, c, d, a, x[12], S24); // 20
    GG(&mut a, b, c, d, x[1], S21); // 21
    GG(&mut d, a, b, c, x[5], S22); // 22
    GG(&mut c, d, a, b, x[9], S23); // 23
    GG(&mut b, c, d, a, x[13], S24); // 24
    GG(&mut a, b, c, d, x[2], S21); // 25
    GG(&mut d, a, b, c, x[6], S22); // 26
    GG(&mut c, d, a, b, x[10], S23); // 27
    GG(&mut b, c, d, a, x[14], S24); // 28
    GG(&mut a, b, c, d, x[3], S21); // 29
    GG(&mut d, a, b, c, x[7], S22); // 30
    GG(&mut c, d, a, b, x[11], S23); // 31
    GG(&mut b, c, d, a, x[15], S24); // 32

    /* Round 3 */
    HH(&mut a, b, c, d, x[0], S31); // 33
    HH(&mut d, a, b, c, x[8], S32); // 34
    HH(&mut c, d, a, b, x[4], S33); // 35
    HH(&mut b, c, d, a, x[12], S34); // 36
    HH(&mut a, b, c, d, x[2], S31); // 37
    HH(&mut d, a, b, c, x[10], S32); // 38
    HH(&mut c, d, a, b, x[6], S33); // 39
    HH(&mut b, c, d, a, x[14], S34); // 40
    HH(&mut a, b, c, d, x[1], S31); // 41
    HH(&mut d, a, b, c, x[9], S32); // 42
    HH(&mut c, d, a, b, x[5], S33); // 43
    HH(&mut b, c, d, a, x[13], S34); // 44
    HH(&mut a, b, c, d, x[3], S31); // 45
    HH(&mut d, a, b, c, x[11], S32); // 46
    HH(&mut c, d, a, b, x[7], S33); // 47
    HH(&mut b, c, d, a, x[15], S34); // 48

    state[0] = state[0].wrapping_add(a);
    state[1] = state[1].wrapping_add(b);
    state[2] = state[2].wrapping_add(c);
    state[3] = state[3].wrapping_add(d);
}

/* Encodes input (UINT4) into output (unsigned char). Assumes len is a multiple of 4. */
fn Encode(output: &mut [u8], input: &[u32]) {
    for (i, in_u32) in input.iter().enumerate() {
        let b = in_u32.to_le_bytes();
        let outs = &mut output[i * 4..];
        outs[0] = b[0];
        outs[1] = b[1];
        outs[2] = b[2];
        outs[3] = b[3];
    }
}

/* Decodes input (unsigned char) into output (UINT4). Assumes len is a multiple of 4. */
fn Decode(output: &mut [u32], input: &[u8]) {
    for (i, out_u32) in output.iter_mut().enumerate() {
        let ii = &input[i * 4..];
        *out_u32 = u32::from_le_bytes([ii[0], ii[1], ii[2], ii[3]]);
    }
}

//===================================================================

pub fn Com_BlockChecksum(buffer: &[u8]) -> u32 {
    let mut ctx = MD4Init();
    MD4Update(&mut ctx, buffer);
    let digest = MD4Final_u32(&mut ctx);
    digest[0] ^ digest[1] ^ digest[2] ^ digest[3]
}

pub fn Com_BlockChecksumKey(buffer: &[u8], key: u32) -> u32 {
    let key_bytes = key.to_le_bytes();
    let mut ctx = MD4Init();
    MD4Update(&mut ctx, &key_bytes);
    MD4Update(&mut ctx, buffer);
    let digest = MD4Final_u32(&mut ctx);
    digest[0] ^ digest[1] ^ digest[2] ^ digest[3]
}
