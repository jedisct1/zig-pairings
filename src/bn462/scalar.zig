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

/// Compressed scalar type.
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
    limbs[7] = (@as(u64, s[0]) << 8) | @as(u64, s[1]);
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

/// Reduce any 58-byte value modulo r.
pub fn reduce(s: CompressedScalar) CompressedScalar {
    const limbs = bytesToLimbs(s);
    return limbsToBytes(reduceLimbs(limbs_count, &limbs));
}

/// Add two scalars modulo r.
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

    result = subtractModulusIfNecessary(result, carry);

    return limbsToBytes(result);
}

/// Subtract two scalars modulo r.
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

    // Add r back if the subtraction borrowed, without branching.
    const mask: u64 = 0 -% @as(u64, borrow);
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

/// Multiply two scalars modulo r.
pub fn mul(a: CompressedScalar, b: CompressedScalar) CompressedScalar {
    const a_limbs = bytesToLimbs(a);
    const b_limbs = bytesToLimbs(b);

    var t: [limbs_count * 2]u64 = @splat(0);
    inline for (0..limbs_count) |i| {
        var carry: u64 = 0;
        inline for (0..limbs_count) |j| {
            const product = @as(u128, a_limbs[i]) * @as(u128, b_limbs[j]) + @as(u128, t[i + j]) + @as(u128, carry);
            t[i + j] = @truncate(product);
            carry = @truncate(product >> 64);
        }
        t[i + limbs_count] = carry;
    }

    const result = reduceLimbs(limbs_count * 2, &t);

    return limbsToBytes(result);
}

/// Reduce an integer of any width modulo r, in constant time.
fn reduceLimbs(comptime n: usize, t: *const [n]u64) [limbs_count]u64 {
    var acc: [limbs_count]u64 = @splat(0);
    var i: usize = n * 64;
    while (i > 0) {
        i -= 1;
        // acc < r, so doubling it can't overflow.
        var carry: u64 = (t[i / 64] >> @intCast(i % 64)) & 1;
        inline for (0..limbs_count) |j| {
            const next = acc[j] >> 63;
            acc[j] = (acc[j] << 1) | carry;
            carry = next;
        }
        acc = subtractModulusIfNecessary(acc, 0);
    }
    return acc;
}

/// Subtract r once if the value, with its carry bit, is at least r, in constant time.
fn subtractModulusIfNecessary(limbs: [limbs_count]u64, carry: u1) [limbs_count]u64 {
    var diff: [limbs_count]u64 = undefined;
    var borrow: u1 = 0;
    inline for (0..limbs_count) |j| {
        const diff1 = @subWithOverflow(limbs[j], modulus[j]);
        const diff2 = @subWithOverflow(diff1[0], borrow);
        diff[j] = diff2[0];
        borrow = diff1[1] | diff2[1];
    }

    const keep_diff: u64 = 0 -% @as(u64, carry | (1 - borrow));
    var result: [limbs_count]u64 = undefined;
    inline for (0..limbs_count) |j| {
        result[j] = (diff[j] & keep_diff) | (limbs[j] & ~keep_diff);
    }
    return result;
}

/// Generate a uniformly random nonzero scalar.
pub fn random(io: std.Io, comptime endian: std.builtin.Endian) CompressedScalar {
    // Rejection sampling keeps the result uniform.
    const top_byte_mask: u8 = 0xff >> @clz(modulus_bytes[0]);
    while (true) {
        var s: CompressedScalar = undefined;
        io.random(&s);
        s[0] &= top_byte_mask;
        rejectNonCanonical(s) catch continue;
        if (isZero(s)) continue;
        return toBytes(s, endian);
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
    const zero: CompressedScalar = @splat(0);
    var one = zero;
    one[one.len - 1] = 1;
    var seven = zero;
    seven[seven.len - 1] = 7;
    var eleven = zero;
    eleven[eleven.len - 1] = 11;

    try std.testing.expectEqual(18, add(seven, eleven)[zero.len - 1]);
    try std.testing.expectEqual(4, sub(eleven, seven)[zero.len - 1]);
    try std.testing.expectEqual(77, mul(seven, eleven)[zero.len - 1]);

    const r_minus_one = sub(zero, one);
    try std.testing.expectEqualSlices(u8, &sub(modulus_bytes, one), &r_minus_one);
    try std.testing.expect(isZero(add(r_minus_one, one)));
    try std.testing.expectEqualSlices(u8, &one, &mul(r_minus_one, r_minus_one));

    // Expected product computed with Python.
    var a: CompressedScalar = undefined;
    for (&a, 1..) |*byte, i| byte.* = @intCast(i);
    var b: CompressedScalar = undefined;
    _ = try std.fmt.hexToBytes(&b, "1592a6ce9d209aa2adacafaecdac77cfd64db88eb1b0b3b2ba5c757c7e31b3b2a584878681a485e69b9cb32eaa888b8a9594078fd17113509d68");
    var expected: CompressedScalar = undefined;
    _ = try std.fmt.hexToBytes(&expected, "10426de032f043d48a50958b5ad85e514ffb2c948f9e96dfc7f14559a4bacf4b2ac70949b3863d12cd687d35e5547baaf65bd18c8f2b89d65989");
    try std.testing.expectEqualSlices(u8, &expected, &mul(reduce(a), b));
}

test "scalar reduction" {
    try std.testing.expect(isZero(reduce(modulus_bytes)));

    // 2^464 - 1 is about 7r. Expected value computed with Python.
    var expected: CompressedScalar = undefined;
    _ = try std.fmt.hexToBytes(&expected, "03e07e85f81ff040000000003f06a540f90bc179fffffffffa177e0c17d3b1f238000000003f042f189c3e5839bfffffffff03f42fc8df8c7fa4");
    try std.testing.expectEqualSlices(u8, &expected, &reduce(@splat(0xff)));
}

test "random scalar" {
    const io = std.testing.io;
    const s = random(io, .big);
    try std.testing.expect(!isZero(s));
    try rejectNonCanonical(s);
}
