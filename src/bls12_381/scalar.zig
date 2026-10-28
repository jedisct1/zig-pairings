//! Scalar field for BLS12-381 (the group order r).
//!
//! r = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
//!
//! This is a 255-bit prime, used for scalar multiplication on both G1 and G2.

const std = @import("std");
const crypto = std.crypto;
const mem = std.mem;

const NonCanonicalError = crypto.errors.NonCanonicalError;

/// Number of 64-bit limbs.
const limbs_count = 4;

/// Encoded scalar length in bytes.
pub const encoded_length = 32;

/// The scalar field modulus r.
pub const modulus: [limbs_count]u64 = .{
    0xffffffff00000001,
    0x53bda402fffe5bfe,
    0x3339d80809a1d805,
    0x73eda753299d7d48,
};

/// Modulus as bytes (big-endian).
pub const modulus_bytes: [encoded_length]u8 = .{
    0x73, 0xed, 0xa7, 0x53, 0x29, 0x9d, 0x7d, 0x48,
    0x33, 0x39, 0xd8, 0x08, 0x09, 0xa1, 0xd8, 0x05,
    0x53, 0xbd, 0xa4, 0x02, 0xff, 0xfe, 0x5b, 0xfe,
    0xff, 0xff, 0xff, 0xff, 0x00, 0x00, 0x00, 0x01,
};

/// Compressed scalar type.
pub const CompressedScalar = [encoded_length]u8;

/// Check if bytes represent a canonical scalar (< r).
pub fn rejectNonCanonical(s: CompressedScalar) NonCanonicalError!void {
    if (crypto.timing_safe.compare(u8, &s, &modulus_bytes, .big) != .lt) {
        return error.NonCanonical;
    }
}

/// Reduce any 32-byte value modulo r.
pub fn reduce(s: CompressedScalar) CompressedScalar {
    var input: [limbs_count]u64 = undefined;
    inline for (0..limbs_count) |i| {
        const j = limbs_count - 1 - i;
        input[i] = mem.readInt(u64, s[j * 8 ..][0..8], .big);
    }

    const limbs = reduceLimbs(limbs_count, &input);

    var result: CompressedScalar = undefined;
    inline for (0..limbs_count) |i| {
        const j = limbs_count - 1 - i;
        mem.writeInt(u64, result[j * 8 ..][0..8], limbs[i], .big);
    }
    return result;
}

/// Add two scalars modulo r.
pub fn add(a: CompressedScalar, b: CompressedScalar) CompressedScalar {
    var a_limbs: [limbs_count]u64 = undefined;
    var b_limbs: [limbs_count]u64 = undefined;

    inline for (0..limbs_count) |i| {
        const j = limbs_count - 1 - i;
        a_limbs[i] = mem.readInt(u64, a[j * 8 ..][0..8], .big);
        b_limbs[i] = mem.readInt(u64, b[j * 8 ..][0..8], .big);
    }

    var result: [limbs_count]u64 = undefined;
    var carry: u1 = 0;
    inline for (0..limbs_count) |i| {
        const sum1 = @addWithOverflow(a_limbs[i], b_limbs[i]);
        const sum2 = @addWithOverflow(sum1[0], carry);
        result[i] = sum2[0];
        carry = sum1[1] | sum2[1];
    }

    result = subtractModulusIfNecessary(result, carry);

    var out: CompressedScalar = undefined;
    inline for (0..limbs_count) |i| {
        const j = limbs_count - 1 - i;
        mem.writeInt(u64, out[j * 8 ..][0..8], result[i], .big);
    }
    return out;
}

/// Subtract two scalars modulo r.
pub fn sub(a: CompressedScalar, b: CompressedScalar) CompressedScalar {
    var a_limbs: [limbs_count]u64 = undefined;
    var b_limbs: [limbs_count]u64 = undefined;

    inline for (0..limbs_count) |i| {
        const j = limbs_count - 1 - i;
        a_limbs[i] = mem.readInt(u64, a[j * 8 ..][0..8], .big);
        b_limbs[i] = mem.readInt(u64, b[j * 8 ..][0..8], .big);
    }

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

    var out: CompressedScalar = undefined;
    inline for (0..limbs_count) |i| {
        const j = limbs_count - 1 - i;
        mem.writeInt(u64, out[j * 8 ..][0..8], result[i], .big);
    }
    return out;
}

/// Negate a scalar modulo r.
pub fn neg(s: CompressedScalar) CompressedScalar {
    const zero: CompressedScalar = @splat(0);
    return sub(zero, s);
}

/// Multiply two scalars modulo r.
pub fn mul(a: CompressedScalar, b: CompressedScalar) CompressedScalar {
    var a_limbs: [limbs_count]u64 = undefined;
    var b_limbs: [limbs_count]u64 = undefined;

    inline for (0..limbs_count) |i| {
        const j = limbs_count - 1 - i;
        a_limbs[i] = mem.readInt(u64, a[j * 8 ..][0..8], .big);
        b_limbs[i] = mem.readInt(u64, b[j * 8 ..][0..8], .big);
    }

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

    var out: CompressedScalar = undefined;
    inline for (0..limbs_count) |i| {
        const j = limbs_count - 1 - i;
        mem.writeInt(u64, out[j * 8 ..][0..8], result[i], .big);
    }
    return out;
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
    for (s) |b| acc |= b;
    return acc == 0;
}

/// Convert from bytes with specified endianness.
pub fn fromBytes(bytes: CompressedScalar, endian: std.builtin.Endian) NonCanonicalError!CompressedScalar {
    const s = if (endian == .big) bytes else blk: {
        var swapped: CompressedScalar = undefined;
        for (bytes, 0..) |b, i| swapped[31 - i] = b;
        break :blk swapped;
    };
    try rejectNonCanonical(s);
    return s;
}

/// Convert to bytes with specified endianness.
pub fn toBytes(s: CompressedScalar, endian: std.builtin.Endian) CompressedScalar {
    if (endian == .big) return s;
    var swapped: CompressedScalar = undefined;
    for (s, 0..) |b, i| swapped[31 - i] = b;
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
    _ = try std.fmt.hexToBytes(&b, "31b700537803265a7a72d7a6a006d3a561f713b3b1b257b3bdbcbfbfb9b8bbb9");
    var expected: CompressedScalar = undefined;
    _ = try std.fmt.hexToBytes(&expected, "3cb9cd63462c0b9cd32206b462a4e4005be1697a154a85f2793dcd124db92190");
    try std.testing.expectEqualSlices(u8, &expected, &mul(reduce(a), b));
}

test "scalar reduction" {
    try std.testing.expect(isZero(reduce(modulus_bytes)));

    // 2^256 - 1 is more than twice r. Expected value computed with Python.
    var expected: CompressedScalar = undefined;
    _ = try std.fmt.hexToBytes(&expected, "1824b159acc5056f998c4fefecbc4ff55884b7fa0003480200000001fffffffd");
    try std.testing.expectEqualSlices(u8, &expected, &reduce(@splat(0xff)));
}

test "random scalar" {
    const io = std.testing.io;
    const s = random(io, .big);
    try std.testing.expect(!isZero(s));
    try rejectNonCanonical(s);
}
