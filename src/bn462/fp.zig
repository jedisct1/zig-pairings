//! Base field Fp for BN462.
//!
//! The field modulus is:
//! p = 0x240480360120023ffffffffff6ff0cf6b7d9bfca0000000000d812908f41c8020ffffffffff6ff66fc6ff687f640000000002401b00840138013
//!
//! This is a 462-bit prime.
//!
//! Field elements are stored in Montgomery form using 8 x 64-bit limbs.
//!
//! BN462 is a Barreto-Naehrig curve with parameter x = 2^114 + 2^101 - 2^14 - 1,
//! providing approximately 134-bit security.

const std = @import("std");
const crypto = std.crypto;
const mem = std.mem;

const NonCanonicalError = crypto.errors.NonCanonicalError;
const NotSquareError = crypto.errors.NotSquareError;

/// Number of 64-bit limbs per field element.
const limbs_count = 8;

/// Field element in Montgomery form.
pub const Fp = struct {
    /// Limbs in little-endian order.
    limbs: [limbs_count]u64,

    /// The field modulus p.
    pub const modulus: [limbs_count]u64 = .{
        0x2401b00840138013,
        0xf687f64000000000,
        0xfffffff6ff66fc6f,
        0x12908f41c8020fff,
        0xbfca0000000000d8,
        0xfffff6ff0cf6b7d9,
        0x80360120023fffff,
        0x0000000000002404,
    };

    /// Montgomery R = 2^512 mod p.
    const r: [limbs_count]u64 = .{
        0x3e11eeb41eee70a8,
        0x2eea3e13fec5661f,
        0x31530e1535b28153,
        0x68ceb6920acb7c40,
        0xfd303bfd856dff73,
        0xffaffe8680658205,
        0xfde7f1801dfec009,
        0x0000000000000fe5,
    };

    /// R^2 mod p, used to convert into Montgomery form.
    const rr: [limbs_count]u64 = .{
        0xffb1ffb6caf1880b,
        0xba49f8b9c4c1a8b2,
        0x9000c34490b9933a,
        0x4284c26b4ec54698,
        0x74c63c7da0391584,
        0x9cdcd35003bb0cd6,
        0x5763230bbc44e2af,
        0x0000000000000273,
    };

    /// -p^(-1) mod 2^64.
    const m0inv: u64 = 0xe718ce9e711bb5e5;

    /// Zero.
    pub const zero: Fp = .{ .limbs = .{ 0, 0, 0, 0, 0, 0, 0, 0 } };

    /// One.
    pub const one: Fp = .{ .limbs = r };

    /// Number of bytes in serialized form.
    pub const encoded_length = 58;

    /// Deserialize a field element from bytes.
    pub fn fromBytes(bytes: [encoded_length]u8, endian: std.builtin.Endian) NonCanonicalError!Fp {
        const s = if (endian == .big) bytes else orderSwap(bytes);

        var limbs: [limbs_count]u64 = undefined;

        limbs[7] = @as(u64, s[0]) << 8 | s[1];

        inline for (0..7) |i| {
            const j = 6 - i;
            const offset = 2 + i * 8;
            limbs[j] = mem.readInt(u64, s[offset..][0..8], .big);
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

        bytes[0] = @truncate(val.limbs[7] >> 8);
        bytes[1] = @truncate(val.limbs[7]);

        inline for (0..7) |i| {
            const j = 6 - i;
            const offset = 2 + i * 8;
            mem.writeInt(u64, bytes[offset..][0..8], val.limbs[j], .big);
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
            fe.limbs[0], fe.limbs[1], fe.limbs[2], fe.limbs[3],
            fe.limbs[4], fe.limbs[5], fe.limbs[6], fe.limbs[7],
            0,           0,           0,           0,
            0,           0,           0,           0,
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

        const mask = 0 -% @as(u64, borrow);
        var carry: u1 = 0;
        inline for (0..limbs_count) |i| {
            const sum1 = @addWithOverflow(result[i], modulus[i] & mask);
            const sum2 = @addWithOverflow(sum1[0], carry);
            result[i] = sum2[0];
            carry = sum1[1] | sum2[1];
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
                const product = @as(u128, a.limbs[i]) * b.limbs[j] + t[i + j] + carry;
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
                const product = @as(u128, k) * modulus[j] + wide[i + j] + carry;
                wide[i + j] = @truncate(product);
                carry = @truncate(product >> 64);
            }

            inline for (i + limbs_count..limbs_count * 2) |idx| {
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
        return a.pow(U512, field_order - 2);
    }

    /// Check if the element is a quadratic residue.
    pub fn isSquare(x2: Fp) bool {
        // Euler's criterion: a^((p-1)/2) = 1 exactly when a is a square.
        const ls = x2.pow(U512, (field_order - 1) / 2);
        return ls.equivalent(one) or x2.isZero();
    }

    /// Square root, if one exists.
    pub fn sqrt(x2: Fp) NotSquareError!Fp {
        // p ≡ 3 (mod 4), so sqrt(a) = a^((p+1)/4).
        const result = x2.pow(U512, (field_order + 1) / 4);

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
            0x1200d8042009c009,
            0xfb43fb2000000000,
            0xfffffffb7fb37e37,
            0x094847a0e40107ff,
            0xdfe500000000006c,
            0xfffffb7f867b5bec,
            0x401b0090011fffff,
            0x0000000000001202,
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
            l[7] = @as(u64, bytes[0]) << 8 | bytes[1];
            for (0..7) |i| {
                const j = 6 - i;
                const offset = 2 + i * 8;
                l[j] = @as(u64, bytes[offset]) << 56 |
                    @as(u64, bytes[offset + 1]) << 48 |
                    @as(u64, bytes[offset + 2]) << 40 |
                    @as(u64, bytes[offset + 3]) << 32 |
                    @as(u64, bytes[offset + 4]) << 24 |
                    @as(u64, bytes[offset + 5]) << 16 |
                    @as(u64, bytes[offset + 6]) << 8 |
                    bytes[offset + 7];
            }
            break :blk l;
        };
        const fp = Fp{ .limbs = limbs };
        return fp.toMontgomery();
    }
};

/// 512-bit unsigned integer type for exponents.
const U512 = @Int(.unsigned, 512);

/// The field order p.
const field_order: U512 = 0x240480360120023ffffffffff6ff0cf6b7d9bfca0000000000d812908f41c8020ffffffffff6ff66fc6ff687f640000000002401b00840138013;

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
