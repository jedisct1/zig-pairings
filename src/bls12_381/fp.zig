//! Base field Fp for BLS12-381.
//!
//! The field modulus is:
//! p = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab
//!
//! This is a 381-bit prime.
//!
//! Field elements are stored in Montgomery form using 6 x 64-bit limbs.

const std = @import("std");
const crypto = std.crypto;
const mem = std.mem;

const NonCanonicalError = crypto.errors.NonCanonicalError;
const NotSquareError = crypto.errors.NotSquareError;

/// Number of 64-bit limbs per field element.
const limbs_count = 6;

/// Field element in Montgomery form.
pub const Fp = struct {
    /// Limbs in little-endian order.
    limbs: [limbs_count]u64,

    /// The field modulus p.
    pub const modulus: [limbs_count]u64 = .{
        0xb9feffffffffaaab,
        0x1eabfffeb153ffff,
        0x6730d2a0f6b0f624,
        0x64774b84f38512bf,
        0x4b1ba7b6434bacd7,
        0x1a0111ea397fe69a,
    };

    /// Montgomery R = 2^384 mod p.
    const r: [limbs_count]u64 = .{
        0x760900000002fffd,
        0xebf4000bc40c0002,
        0x5f48985753c758ba,
        0x77ce585370525745,
        0x5c071a97a256ec6d,
        0x15f65ec3fa80e493,
    };

    /// R^2 mod p, used to convert into Montgomery form.
    const rr: [limbs_count]u64 = .{
        0xf4df1f341c341746,
        0x0a76e6a609d104f1,
        0x8de5476c4c95b6d5,
        0x67eb88a9939d83c0,
        0x9a793e85b519952d,
        0x11988fe592cae3aa,
    };

    /// -p^(-1) mod 2^64.
    const m0inv: u64 = 0x89f3fffcfffcfffd;

    /// Zero.
    pub const zero: Fp = .{ .limbs = .{ 0, 0, 0, 0, 0, 0 } };

    /// One.
    pub const one: Fp = .{ .limbs = r };

    /// Number of bytes in serialized form.
    pub const encoded_length = 48;

    /// Deserialize a field element from bytes.
    pub fn fromBytes(bytes: [encoded_length]u8, endian: std.builtin.Endian) NonCanonicalError!Fp {
        const s = if (endian == .big) bytes else orderSwap(bytes);

        var limbs: [limbs_count]u64 = undefined;
        inline for (0..limbs_count) |i| {
            const j = limbs_count - 1 - i;
            limbs[i] = mem.readInt(u64, s[j * 8 ..][0..8], .big);
        }

        if (!lessThanModulus(limbs)) {
            return error.NonCanonical;
        }

        const fp = Fp{ .limbs = limbs };
        return fp.toMontgomery();
    }

    /// Serialize a field element to bytes.
    pub fn toBytes(fe: Fp, endian: std.builtin.Endian) [encoded_length]u8 {
        const val = fe.fromMontgomery();

        var bytes: [encoded_length]u8 = undefined;
        inline for (0..limbs_count) |i| {
            const j = limbs_count - 1 - i;
            mem.writeInt(u64, bytes[j * 8 ..][0..8], val.limbs[i], .big);
        }

        return if (endian == .big) bytes else orderSwap(bytes);
    }

    /// Swap byte order.
    fn orderSwap(s: [encoded_length]u8) [encoded_length]u8 {
        var t: [encoded_length]u8 = undefined;
        for (s, 0..) |x, i| t[t.len - 1 - i] = x;
        return t;
    }

    /// Check if limbs represent a value less than the modulus.
    fn lessThanModulus(limbs: [limbs_count]u64) bool {
        var borrow: u1 = 0;
        inline for (0..limbs_count) |i| {
            const result = @subWithOverflow(limbs[i], modulus[i]);
            const result2 = @subWithOverflow(result[0], borrow);
            borrow = result[1] | result2[1];
        }
        return borrow == 1;
    }

    /// Convert to Montgomery form.
    fn toMontgomery(fe: Fp) Fp {
        return fe.mul(.{ .limbs = rr });
    }

    /// Convert from Montgomery form.
    fn fromMontgomery(fe: Fp) Fp {
        var result: [limbs_count]u64 = undefined;
        montgomeryReduce(&result, &.{
            fe.limbs[0], fe.limbs[1], fe.limbs[2],
            fe.limbs[3], fe.limbs[4], fe.limbs[5],
            0,           0,           0,
            0,           0,           0,
        });
        return .{ .limbs = result };
    }

    /// Return true if the element is zero.
    pub fn isZero(fe: Fp) bool {
        var acc: u64 = 0;
        inline for (fe.limbs) |limb| {
            acc |= limb;
        }
        return acc == 0;
    }

    /// Return true if the element is one.
    pub fn isOne(fe: Fp) bool {
        return fe.equivalent(one);
    }

    /// Return true if both elements are equivalent.
    pub fn equivalent(a: Fp, b: Fp) bool {
        return a.sub(b).isZero();
    }

    /// Return true if the element is odd.
    pub fn isOdd(fe: Fp) bool {
        const val = fe.fromMontgomery();
        return (val.limbs[0] & 1) == 1;
    }

    /// Conditionally replace the element with `other` when `choice` is 1 (constant-time).
    pub fn cMov(fe: *Fp, other: Fp, choice: u1) void {
        const mask: u64 = 0 -% @as(u64, choice);
        inline for (0..limbs_count) |i| {
            fe.limbs[i] ^= mask & (fe.limbs[i] ^ other.limbs[i]);
        }
    }

    /// Negate the field element.
    pub fn neg(a: Fp) Fp {
        return zero.sub(a);
    }

    /// Add two field elements.
    pub fn add(a: Fp, b: Fp) Fp {
        var result: [limbs_count]u64 = undefined;
        var carry: u1 = 0;

        inline for (0..limbs_count) |i| {
            const sum1 = @addWithOverflow(a.limbs[i], b.limbs[i]);
            const sum2 = @addWithOverflow(sum1[0], carry);
            result[i] = sum2[0];
            carry = sum1[1] | sum2[1];
        }

        return subtractModulusIfNecessary(result, carry);
    }

    /// Subtract two field elements.
    pub fn sub(a: Fp, b: Fp) Fp {
        var result: [limbs_count]u64 = undefined;
        var borrow: u1 = 0;

        inline for (0..limbs_count) |i| {
            const diff1 = @subWithOverflow(a.limbs[i], b.limbs[i]);
            const diff2 = @subWithOverflow(diff1[0], borrow);
            result[i] = diff2[0];
            borrow = diff1[1] | diff2[1];
        }

        // A borrow means the result went negative, so bring it back with a modulus add.
        if (borrow == 1) {
            var carry: u1 = 0;
            inline for (0..limbs_count) |i| {
                const sum1 = @addWithOverflow(result[i], modulus[i]);
                const sum2 = @addWithOverflow(sum1[0], carry);
                result[i] = sum2[0];
                carry = sum1[1] | sum2[1];
            }
        }

        return .{ .limbs = result };
    }

    /// Double a field element.
    pub fn dbl(a: Fp) Fp {
        return a.add(a);
    }

    /// Multiply two field elements.
    pub fn mul(a: Fp, b: Fp) Fp {
        var t: [limbs_count * 2]u64 = @splat(0);

        inline for (0..limbs_count) |i| {
            var carry: u64 = 0;
            inline for (0..limbs_count) |j| {
                const product = @as(u128, a.limbs[i]) * @as(u128, b.limbs[j]) + @as(u128, t[i + j]) + @as(u128, carry);
                t[i + j] = @truncate(product);
                carry = @truncate(product >> 64);
            }
            t[i + limbs_count] = carry;
        }

        var result: [limbs_count]u64 = undefined;
        montgomeryReduce(&result, &t);
        return .{ .limbs = result };
    }

    /// Square a field element.
    pub fn sq(a: Fp) Fp {
        return a.mul(a);
    }

    /// Montgomery reduction.
    fn montgomeryReduce(result: *[limbs_count]u64, t: *const [limbs_count * 2]u64) void {
        var wide: [limbs_count * 2]u64 = t.*;

        inline for (0..limbs_count) |i| {
            const k: u64 = wide[i] *% m0inv;
            var carry: u64 = 0;

            inline for (0..limbs_count) |j| {
                const product = @as(u128, k) * @as(u128, modulus[j]) + @as(u128, wide[i + j]) + @as(u128, carry);
                wide[i + j] = @truncate(product);
                carry = @truncate(product >> 64);
            }

            var idx = i + limbs_count;
            while (carry != 0 and idx < limbs_count * 2) : (idx += 1) {
                const sum = @addWithOverflow(wide[idx], carry);
                wide[idx] = sum[0];
                carry = sum[1];
            }
        }

        inline for (0..limbs_count) |i| {
            result[i] = wide[i + limbs_count];
        }

        result.* = subtractModulusIfNecessary(result.*, 0).limbs;
    }

    /// Subtract modulus if result >= modulus.
    fn subtractModulusIfNecessary(limbs: [limbs_count]u64, carry: u1) Fp {
        var result: [limbs_count]u64 = undefined;
        var borrow: u1 = 0;

        inline for (0..limbs_count) |i| {
            const diff1 = @subWithOverflow(limbs[i], modulus[i]);
            const diff2 = @subWithOverflow(diff1[0], borrow);
            result[i] = diff2[0];
            borrow = diff1[1] | diff2[1];
        }

        // Keep the original value only when the subtraction underflowed with no incoming carry.
        const keep_original = @intFromBool(carry == 0) & borrow;

        var final: [limbs_count]u64 = undefined;
        const mask: u64 = 0 -% @as(u64, keep_original);
        inline for (0..limbs_count) |i| {
            final[i] = (mask & limbs[i]) | (~mask & result[i]);
        }

        return .{ .limbs = final };
    }

    /// Multiplicative inverse, obtained from Fermat's little theorem.
    ///
    /// Returns zero when the input is zero.
    pub fn invert(a: Fp) Fp {
        return a.pow(U384, field_order - 2);
    }

    /// Check if the element is a quadratic residue.
    pub fn isSquare(x2: Fp) bool {
        // Euler's criterion: a^((p-1)/2) = 1 exactly when a is a square.
        const ls = x2.pow(U384, (field_order - 1) / 2);
        return ls.equivalent(one) or x2.isZero();
    }

    /// Square root, if one exists.
    pub fn sqrt(x2: Fp) NotSquareError!Fp {
        // p ≡ 3 (mod 4), so sqrt(a) = a^((p+1)/4).
        const result = x2.pow(U384, (field_order + 1) / 4);

        if (result.sq().equivalent(x2)) {
            return result;
        }
        return error.NotSquare;
    }

    /// Compute a^n using square-and-multiply.
    pub fn pow(a: Fp, comptime T: type, comptime n: T) Fp {
        var result = one;
        var base = a;
        var exp = n;

        while (exp != 0) {
            if (exp & 1 == 1) {
                result = result.mul(base);
            }
            base = base.sq();
            exp >>= 1;
        }

        return result;
    }

    /// Return true if the element is the lexicographically largest of {x, -x}.
    pub fn lexicographicallyLargest(fe: Fp) bool {
        const val = fe.fromMontgomery();
        const half_modulus: [limbs_count]u64 = .{
            0xdcff7fffffffd555,
            0x0f55ffff58a9ffff,
            0xb39869507b587b12,
            0xb23ba5c279c2895f,
            0x258dd3db21a5d66b,
            0x0d0088f51cbff34d,
        };

        var borrow: u1 = 0;
        inline for (0..limbs_count) |i| {
            const diff1 = @subWithOverflow(half_modulus[i], val.limbs[i]);
            const diff2 = @subWithOverflow(diff1[0], borrow);
            borrow = diff1[1] | diff2[1];
        }
        return borrow == 1;
    }

    /// Create from integer (comptime).
    pub fn fromInt(comptime x: comptime_int) Fp {
        @setEvalBranchQuota(100000);
        comptime var bytes: [encoded_length]u8 = undefined;
        comptime {
            var val = x;
            for (0..encoded_length) |i| {
                bytes[encoded_length - 1 - i] = @truncate(val);
                val >>= 8;
            }
        }
        const limbs: [limbs_count]u64 = comptime blk: {
            var l: [limbs_count]u64 = undefined;
            for (0..limbs_count) |i| {
                const j = limbs_count - 1 - i;
                l[i] = @as(u64, bytes[j * 8]) << 56 |
                    @as(u64, bytes[j * 8 + 1]) << 48 |
                    @as(u64, bytes[j * 8 + 2]) << 40 |
                    @as(u64, bytes[j * 8 + 3]) << 32 |
                    @as(u64, bytes[j * 8 + 4]) << 24 |
                    @as(u64, bytes[j * 8 + 5]) << 16 |
                    @as(u64, bytes[j * 8 + 6]) << 8 |
                    @as(u64, bytes[j * 8 + 7]);
            }
            break :blk l;
        };
        const fp = Fp{ .limbs = limbs };
        return fp.toMontgomery();
    }
};

/// 384-bit unsigned integer type for exponents.
const U384 = @Int(.unsigned, 384);

/// The field order p.
const field_order: U384 = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab;

test "fp basic arithmetic" {
    const a = Fp.fromInt(7);
    const b = Fp.fromInt(11);

    const sum = a.add(b);
    try std.testing.expect(sum.equivalent(Fp.fromInt(18)));

    const diff = b.sub(a);
    try std.testing.expect(diff.equivalent(Fp.fromInt(4)));

    const prod = a.mul(b);
    try std.testing.expect(prod.equivalent(Fp.fromInt(77)));

    const sq = a.sq();
    try std.testing.expect(sq.equivalent(Fp.fromInt(49)));
}

test "fp zero and one" {
    try std.testing.expect(Fp.zero.isZero());
    try std.testing.expect(!Fp.one.isZero());
    try std.testing.expect(Fp.one.isOne());
    try std.testing.expect(!Fp.zero.isOne());
}

test "fp negation" {
    const a = Fp.fromInt(7);
    const neg_a = a.neg();
    try std.testing.expect(a.add(neg_a).isZero());
}

test "fp inversion" {
    const a = Fp.fromInt(7);
    const inv_a = a.invert();
    try std.testing.expect(a.mul(inv_a).isOne());
}

test "fp serialization roundtrip" {
    const a = Fp.fromInt(12345678901234567890);
    const bytes = a.toBytes(.big);
    const b = try Fp.fromBytes(bytes, .big);
    try std.testing.expect(a.equivalent(b));
}
