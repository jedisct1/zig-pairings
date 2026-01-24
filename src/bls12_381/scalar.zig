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

/// Reduce a scalar modulo r if necessary.
pub fn reduce(s: CompressedScalar) CompressedScalar {
    var limbs: [limbs_count]u64 = undefined;
    inline for (0..limbs_count) |i| {
        const j = limbs_count - 1 - i;
        limbs[i] = mem.readInt(u64, s[j * 8 ..][0..8], .big);
    }

    var borrow: u1 = 0;
    inline for (0..limbs_count) |i| {
        const result = @subWithOverflow(limbs[i], modulus[i]);
        const result2 = @subWithOverflow(result[0], borrow);
        borrow = result[1] | result2[1];
    }

    // No borrow means the value is >= r, so subtract r once.
    if (borrow == 0) {
        borrow = 0;
        inline for (0..limbs_count) |i| {
            const result = @subWithOverflow(limbs[i], modulus[i]);
            const result2 = @subWithOverflow(result[0], borrow);
            limbs[i] = result2[0];
            borrow = result[1] | result2[1];
        }
    }

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

    // A borrow means the result went negative, so add r back.
    if (borrow == 1) {
        var carry: u1 = 0;
        inline for (0..limbs_count) |i| {
            const sum1 = @addWithOverflow(result[i], modulus[i]);
            const sum2 = @addWithOverflow(sum1[0], carry);
            result[i] = sum2[0];
            carry = sum1[1] | sum2[1];
        }
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

    const result = barrettReduce(&t);

    var out: CompressedScalar = undefined;
    inline for (0..limbs_count) |i| {
        const j = limbs_count - 1 - i;
        mem.writeInt(u64, out[j * 8 ..][0..8], result[i], .big);
    }
    return out;
}

/// Reduce modulo r by repeated subtraction of r.
fn barrettReduce(t: *const [limbs_count * 2]u64) [limbs_count]u64 {
    var result: [limbs_count]u64 = undefined;
    inline for (0..limbs_count) |i| {
        result[i] = t[i];
    }

    while (true) {
        var borrow: u1 = 0;
        var reduced: [limbs_count]u64 = undefined;
        inline for (0..limbs_count) |i| {
            const diff1 = @subWithOverflow(result[i], modulus[i]);
            const diff2 = @subWithOverflow(diff1[0], borrow);
            reduced[i] = diff2[0];
            borrow = diff1[1] | diff2[1];
        }

        if (borrow == 1) {
            break;
        }
        result = reduced;
    }

    return result;
}

/// Generate a random scalar.
pub fn random(io: std.Io, comptime endian: std.builtin.Endian) CompressedScalar {
    var bytes: CompressedScalar = undefined;
    io.random(&bytes);

    // Clear the top bits so the value is below 2r; reduce() then brings it below r.
    if (endian == .big) {
        bytes[0] &= 0x73;
    } else {
        bytes[31] &= 0x73;
    }

    return reduce(if (endian == .big) bytes else blk: {
        var swapped: CompressedScalar = undefined;
        for (bytes, 0..) |b, i| swapped[31 - i] = b;
        break :blk swapped;
    });
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
    const a: CompressedScalar = .{
        0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x07,
    };
    const b: CompressedScalar = .{
        0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x0b,
    };

    // 7 + 11 = 18
    const sum = add(a, b);
    try std.testing.expectEqual(0x12, sum[31]);

    // 11 - 7 = 4
    const diff = sub(b, a);
    try std.testing.expectEqual(0x04, diff[31]);

    // 7 * 11 = 77
    const prod = mul(a, b);
    try std.testing.expectEqual(0x4d, prod[31]);
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

test "scalar deserialization is canonical" {
    // A serialized scalar must decode to a value strictly below r, so r itself
    // and anything larger is rejected while r - 1 is accepted.
    try std.testing.expectError(error.NonCanonical, fromBytes(modulus_bytes, .big));
    const all_ones: CompressedScalar = @splat(0xff);
    try std.testing.expectError(error.NonCanonical, fromBytes(all_ones, .big));

    var one: CompressedScalar = @splat(0);
    one[one.len - 1] = 1;
    const r_minus_one = sub(modulus_bytes, one);
    _ = try fromBytes(r_minus_one, .big);
}
