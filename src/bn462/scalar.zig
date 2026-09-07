//! Scalar field for BN462 (the group order r).
//!
//! r = 0x240480360120023ffffffffff6ff0cf6b7d9bfca0000000000d812908ee1c201f7fffffffff6ff66fc7bf717f7c0000000002401b007e010800d
//!
//! This is a 462-bit prime, used for scalar multiplication on both G1 and G2.

const std = @import("std");
const crypto = std.crypto;
const mem = std.mem;

const NonCanonicalError = crypto.errors.NonCanonicalError;

/// Number of 64-bit limbs.
const limbs_count = 8;

/// Encoded scalar length in bytes.
pub const encoded_length = 58;

/// The scalar field modulus r.
pub const modulus: [limbs_count]u64 = .{
    0x2401b007e010800d,
    0xf717f7c000000000,
    0xfffffff6ff66fc7b,
    0x12908ee1c201f7ff,
    0xbfca0000000000d8,
    0xfffff6ff0cf6b7d9,
    0x80360120023fffff,
    0x0000000000002404,
};

/// Modulus as bytes (big-endian).
pub const modulus_bytes: [encoded_length]u8 = .{
    0x24, 0x04, 0x80, 0x36, 0x01, 0x20, 0x02, 0x3f,
    0xff, 0xff, 0xff, 0xff, 0xf6, 0xff, 0x0c, 0xf6,
    0xb7, 0xd9, 0xbf, 0xca, 0x00, 0x00, 0x00, 0x00,
    0x00, 0xd8, 0x12, 0x90, 0x8e, 0xe1, 0xc2, 0x01,
    0xf7, 0xff, 0xff, 0xff, 0xff, 0xf6, 0xff, 0x66,
    0xfc, 0x7b, 0xf7, 0x17, 0xf7, 0xc0, 0x00, 0x00,
    0x00, 0x00, 0x24, 0x01, 0xb0, 0x07, 0xe0, 0x10,
    0x80, 0x0d,
};

/// Arithmetic and fromBytes() use canonical big-endian values in [0, r - 1].
pub const CompressedScalar = [encoded_length]u8;

/// Check if bytes represent a canonical scalar (< r).
pub fn rejectNonCanonical(s: CompressedScalar) NonCanonicalError!void {
    if (crypto.timing_safe.compare(u8, &s, &modulus_bytes, .big) != .lt) {
        return error.NonCanonical;
    }
}

/// Convert big-endian bytes to little-endian limbs.
///
/// The modulus is not a whole number of limbs wide, so the leading two bytes
/// stand alone in the high limb.
fn bytesToLimbs(s: CompressedScalar) [limbs_count]u64 {
    var limbs: [limbs_count]u64 = @splat(0);
    limbs[7] = @as(u64, s[0]) << 8 | s[1];
    inline for (0..7) |i| {
        const j = 6 - i;
        const offset = 2 + i * 8;
        limbs[j] = mem.readInt(u64, s[offset..][0..8], .big);
    }
    return limbs;
}

/// Convert little-endian limbs to big-endian bytes.
fn limbsToBytes(limbs: [limbs_count]u64) CompressedScalar {
    var result: CompressedScalar = undefined;
    result[0] = @truncate(limbs[7] >> 8);
    result[1] = @truncate(limbs[7]);
    inline for (0..7) |i| {
        const j = 6 - i;
        const offset = 2 + i * 8;
        mem.writeInt(u64, result[offset..][0..8], limbs[j], .big);
    }
    return result;
}

/// Reduce a big-endian integer modulo r.
pub fn reduce(s: CompressedScalar) CompressedScalar {
    var limbs = bytesToLimbs(s);

    limbs = @import("../scalar_reduce.zig").reduce(modulus, &limbs);

    return limbsToBytes(limbs);
}

/// Add two canonical big-endian scalars modulo r.
pub fn add(a: CompressedScalar, b: CompressedScalar) CompressedScalar {
    const a_limbs = bytesToLimbs(a);
    const b_limbs = bytesToLimbs(b);

    var result: [limbs_count]u64 = undefined;
    var carry: u1 = 0;
    inline for (0..limbs_count) |i| {
        const sum1 = @addWithOverflow(a_limbs[i], b_limbs[i]);
        const sum2 = @addWithOverflow(sum1[0], carry);
        result[i] = sum2[0];
        carry = sum1[1] | sum2[1];
    }

    var borrow: u1 = 0;
    var reduced: [limbs_count]u64 = undefined;
    inline for (0..limbs_count) |i| {
        const diff1 = @subWithOverflow(result[i], modulus[i]);
        const diff2 = @subWithOverflow(diff1[0], borrow);
        reduced[i] = diff2[0];
        borrow = diff1[1] | diff2[1];
    }

    // Use the reduced value when the sum overflowed or did not go below r.
    const use_reduced = @intFromBool(carry == 1) | (1 - borrow);
    const mask: u64 = 0 -% @as(u64, use_reduced);

    inline for (0..limbs_count) |i| {
        result[i] = (mask & reduced[i]) | (~mask & result[i]);
    }

    return limbsToBytes(result);
}

/// Subtract two canonical big-endian scalars modulo r.
pub fn sub(a: CompressedScalar, b: CompressedScalar) CompressedScalar {
    const a_limbs = bytesToLimbs(a);
    const b_limbs = bytesToLimbs(b);

    var result: [limbs_count]u64 = undefined;
    var borrow: u1 = 0;
    inline for (0..limbs_count) |i| {
        const diff1 = @subWithOverflow(a_limbs[i], b_limbs[i]);
        const diff2 = @subWithOverflow(diff1[0], borrow);
        result[i] = diff2[0];
        borrow = diff1[1] | diff2[1];
    }

    const mask = 0 -% @as(u64, borrow);
    var carry: u1 = 0;
    inline for (0..limbs_count) |i| {
        const sum1 = @addWithOverflow(result[i], modulus[i] & mask);
        const sum2 = @addWithOverflow(sum1[0], carry);
        result[i] = sum2[0];
        carry = sum1[1] | sum2[1];
    }

    return limbsToBytes(result);
}

/// Negate a scalar modulo r.
pub fn neg(s: CompressedScalar) CompressedScalar {
    const zero: CompressedScalar = @splat(0);
    return sub(zero, s);
}

/// Multiply two canonical big-endian scalars modulo r.
pub fn mul(a: CompressedScalar, b: CompressedScalar) CompressedScalar {
    const a_limbs = bytesToLimbs(a);
    const b_limbs = bytesToLimbs(b);

    var t: [limbs_count * 2]u64 = @splat(0);
    inline for (0..limbs_count) |i| {
        var carry: u64 = 0;
        inline for (0..limbs_count) |j| {
            const product = @as(u128, a_limbs[i]) * b_limbs[j] + t[i + j] + carry;
            t[i + j] = @truncate(product);
            carry = @truncate(product >> 64);
        }
        t[i + limbs_count] = carry;
    }

    const result = @import("../scalar_reduce.zig").reduce(modulus, &t);

    return limbsToBytes(result);
}

/// Sample uniformly from [0, r - 1] in the requested byte order.
pub fn random(io: std.Io, comptime endian: std.builtin.Endian) CompressedScalar {
    var bytes: CompressedScalar = undefined;
    while (true) {
        io.random(&bytes);
        bytes[0] &= 0x3f;
        rejectNonCanonical(bytes) catch continue;
        return toBytes(bytes, endian);
    }
}

/// Check if scalar is zero.
pub fn isZero(s: CompressedScalar) bool {
    var acc: u8 = 0;
    for (s) |byte| acc |= byte;
    return acc == 0;
}

/// Convert from bytes with specified endianness.
pub fn fromBytes(bytes: CompressedScalar, endian: std.builtin.Endian) NonCanonicalError!CompressedScalar {
    const s = if (endian == .big) bytes else blk: {
        var swapped: CompressedScalar = undefined;
        for (bytes, 0..) |byte, i| swapped[57 - i] = byte;
        break :blk swapped;
    };
    try rejectNonCanonical(s);
    return s;
}

/// Convert to bytes with specified endianness.
pub fn toBytes(s: CompressedScalar, endian: std.builtin.Endian) CompressedScalar {
    if (endian == .big) return s;
    var swapped: CompressedScalar = undefined;
    for (s, 0..) |byte, i| swapped[57 - i] = byte;
    return swapped;
}

test "scalar arithmetic" {
    const a: CompressedScalar = blk: {
        var bytes: CompressedScalar = @splat(0);
        bytes[57] = 0x07;
        break :blk bytes;
    };
    const b: CompressedScalar = blk: {
        var bytes: CompressedScalar = @splat(0);
        bytes[57] = 0x0b;
        break :blk bytes;
    };

    // 7 + 11 = 18
    const sum = add(a, b);
    try std.testing.expectEqual(0x12, sum[57]);

    // 11 - 7 = 4
    const diff = sub(b, a);
    try std.testing.expectEqual(0x04, diff[57]);

    // 7 * 11 = 77
    const prod = mul(a, b);
    try std.testing.expectEqual(0x4d, prod[57]);
}

test "scalar reduction" {
    const s = modulus_bytes;
    const reduced = reduce(s);
    try std.testing.expect(isZero(reduced));
}

test "random scalar" {
    const io = std.testing.io;
    const s = random(io, .big);
    try std.testing.expect(!isZero(s));
    try rejectNonCanonical(s);
}
