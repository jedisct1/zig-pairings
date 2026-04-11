//! Quadratic extension field Fp2 = Fp[u] / (u^2 + 1) for BN462.
//!
//! Elements are represented as a + b*u where a, b ∈ Fp and u^2 = -1.
//!
//! This is the base field for G2 curve points.

const std = @import("std");
const crypto = std.crypto;
const Fp = @import("fp.zig").Fp;

const NonCanonicalError = crypto.errors.NonCanonicalError;
const NotSquareError = crypto.errors.NotSquareError;

/// Element of the quadratic extension field Fp2.
pub const Fp2 = struct {
    /// Real component (coefficient of 1).
    c0: Fp,
    /// Imaginary component (coefficient of u).
    c1: Fp,

    /// Zero.
    pub const zero: Fp2 = .{ .c0 = Fp.zero, .c1 = Fp.zero };

    /// One.
    pub const one: Fp2 = .{ .c0 = Fp.one, .c1 = Fp.zero };

    /// Number of bytes in serialized form.
    pub const encoded_length = 116;

    /// Deserialize an Fp2 element from bytes.
    ///
    /// In big-endian the layout is c1 || c0, the higher-degree coefficient first.
    pub fn fromBytes(bytes: [encoded_length]u8, endian: std.builtin.Endian) NonCanonicalError!Fp2 {
        if (endian == .big) {
            return .{
                .c1 = try Fp.fromBytes(bytes[0..58].*, .big),
                .c0 = try Fp.fromBytes(bytes[58..116].*, .big),
            };
        } else {
            return .{
                .c0 = try Fp.fromBytes(bytes[0..58].*, .little),
                .c1 = try Fp.fromBytes(bytes[58..116].*, .little),
            };
        }
    }

    /// Serialize an Fp2 element to bytes.
    ///
    /// In big-endian the layout is c1 || c0, the higher-degree coefficient first.
    pub fn toBytes(fe: Fp2, endian: std.builtin.Endian) [encoded_length]u8 {
        var bytes: [encoded_length]u8 = undefined;
        if (endian == .big) {
            bytes[0..58].* = fe.c1.toBytes(.big);
            bytes[58..116].* = fe.c0.toBytes(.big);
        } else {
            bytes[0..58].* = fe.c0.toBytes(.little);
            bytes[58..116].* = fe.c1.toBytes(.little);
        }
        return bytes;
    }

    /// Return true if the element is zero.
    pub fn isZero(fe: Fp2) bool {
        return fe.c0.isZero() and fe.c1.isZero();
    }

    /// Return true if both elements are equivalent.
    pub fn equivalent(a: Fp2, b: Fp2) bool {
        return a.c0.equivalent(b.c0) and a.c1.equivalent(b.c1);
    }

    /// Conditionally replace the element with `other` when `choice` is 1 (constant-time).
    pub fn cMov(fe: *Fp2, other: Fp2, choice: u1) void {
        fe.c0.cMov(other.c0, choice);
        fe.c1.cMov(other.c1, choice);
    }

    /// Negate the element.
    pub fn neg(a: Fp2) Fp2 {
        return .{
            .c0 = a.c0.neg(),
            .c1 = a.c1.neg(),
        };
    }

    /// Add two Fp2 elements.
    pub fn add(a: Fp2, b: Fp2) Fp2 {
        return .{
            .c0 = a.c0.add(b.c0),
            .c1 = a.c1.add(b.c1),
        };
    }

    /// Subtract two Fp2 elements.
    pub fn sub(a: Fp2, b: Fp2) Fp2 {
        return .{
            .c0 = a.c0.sub(b.c0),
            .c1 = a.c1.sub(b.c1),
        };
    }

    /// Double an Fp2 element.
    pub fn dbl(a: Fp2) Fp2 {
        return .{
            .c0 = a.c0.dbl(),
            .c1 = a.c1.dbl(),
        };
    }

    /// Multiply two Fp2 elements, trading one multiplication for extra additions.
    pub fn mul(a: Fp2, b: Fp2) Fp2 {
        const aa = a.c0.mul(b.c0);
        const bb = a.c1.mul(b.c1);
        const sum_a = a.c0.add(a.c1);
        const sum_b = b.c0.add(b.c1);

        return .{
            .c0 = aa.sub(bb),
            .c1 = sum_a.mul(sum_b).sub(aa).sub(bb),
        };
    }

    /// Square an Fp2 element.
    pub fn sq(a: Fp2) Fp2 {
        const sum = a.c0.add(a.c1);
        const diff = a.c0.sub(a.c1);

        return .{
            .c0 = sum.mul(diff),
            .c1 = a.c0.mul(a.c1).dbl(),
        };
    }

    /// Multiply by a scalar in Fp.
    pub fn mulByFp(a: Fp2, scalar: Fp) Fp2 {
        return .{
            .c0 = a.c0.mul(scalar),
            .c1 = a.c1.mul(scalar),
        };
    }

    /// Conjugate: (a + bu) -> (a - bu).
    pub fn conjugate(a: Fp2) Fp2 {
        return .{
            .c0 = a.c0,
            .c1 = a.c1.neg(),
        };
    }

    /// Norm: N(a + bu) = a^2 + b^2.
    pub fn norm(a: Fp2) Fp {
        return a.c0.sq().add(a.c1.sq());
    }

    /// Multiplicative inverse.
    pub fn invert(a: Fp2) Fp2 {
        const norm_inv = a.norm().invert();
        return .{
            .c0 = a.c0.mul(norm_inv),
            .c1 = a.c1.neg().mul(norm_inv),
        };
    }

    /// Multiply by the tower non-residue (2 + u).
    pub fn mulByNonresidue(a: Fp2) Fp2 {
        const two = Fp.fromInt(2);
        return .{
            .c0 = a.c0.mul(two).sub(a.c1),
            .c1 = a.c0.add(a.c1.mul(two)),
        };
    }

    /// Frobenius map (p-th power).
    ///
    /// The shape of the modulus makes this the same as conjugation.
    pub fn frobeniusMap(a: Fp2, power: u8) Fp2 {
        if (power % 2 == 0) {
            return a;
        }
        return a.conjugate();
    }

    /// Check if the element is a quadratic residue in Fp2.
    pub fn isSquare(x2: Fp2) bool {
        // A value is a square in Fp2 iff its norm is a square in Fp.
        return x2.norm().isSquare() or x2.isZero();
    }

    /// Square root, if one exists.
    ///
    /// Uses the complex method, which applies because Fp2 = Fp[u]/(u^2 + 1).
    pub fn sqrt(x2: Fp2) NotSquareError!Fp2 {
        if (x2.isZero()) {
            return zero;
        }

        const c0_sq = x2.c0.sq();
        const c1_sq = x2.c1.sq();
        const n = c0_sq.add(c1_sq);

        const sqrt_n = n.sqrt() catch return error.NotSquare;

        // Of the two candidates for the real part, take whichever one is a square.
        const two = Fp.fromInt(2);
        const two_inv = two.invert();

        var delta = x2.c0.add(sqrt_n);
        if (!delta.isSquare()) {
            delta = x2.c0.sub(sqrt_n);
        }

        const delta_half = delta.mul(two_inv);

        const sqrt_delta_half = delta_half.sqrt() catch {
            // Shouldn't happen after the isSquare check; fall back to the other branch for safety.
            const other_delta = x2.c0.sub(sqrt_n).add(x2.c0.add(sqrt_n)).sub(delta);
            const other_half = other_delta.mul(two_inv);
            const other_sqrt = other_half.sqrt() catch return error.NotSquare;

            const result_c0 = other_sqrt;
            const result_c1 = x2.c1.mul(other_sqrt.dbl().invert());
            const result = Fp2{ .c0 = result_c0, .c1 = result_c1 };

            if (result.sq().equivalent(x2)) {
                return result;
            }
            return error.NotSquare;
        };

        const denom = sqrt_delta_half.dbl();
        const result_c1 = x2.c1.mul(denom.invert());

        const result = Fp2{ .c0 = sqrt_delta_half, .c1 = result_c1 };

        if (result.sq().equivalent(x2)) {
            return result;
        }

        // The other square root is the negation.
        const neg_result = result.neg();
        if (neg_result.sq().equivalent(x2)) {
            return neg_result;
        }

        return error.NotSquare;
    }

    /// Return true if the element is the lexicographically larger of {x, -x}.
    pub fn lexicographicallyLargest(fe: Fp2) bool {
        if (fe.c1.isZero()) {
            return fe.c0.lexicographicallyLargest();
        }
        return fe.c1.lexicographicallyLargest();
    }

    /// Create from integers (comptime).
    pub fn fromInts(comptime c0: comptime_int, comptime c1: comptime_int) Fp2 {
        return .{
            .c0 = Fp.fromInt(c0),
            .c1 = Fp.fromInt(c1),
        };
    }
};

test "fp2 basic arithmetic" {
    const a = Fp2.fromInts(3, 4);
    const b = Fp2.fromInts(1, 2);

    const sum = a.add(b);
    try std.testing.expect(sum.c0.equivalent(Fp.fromInt(4)));
    try std.testing.expect(sum.c1.equivalent(Fp.fromInt(6)));

    const diff = a.sub(b);
    try std.testing.expect(diff.c0.equivalent(Fp.fromInt(2)));
    try std.testing.expect(diff.c1.equivalent(Fp.fromInt(2)));
}

test "fp2 multiplication" {
    const a = Fp2.fromInts(3, 4);
    const b = Fp2.fromInts(1, 2);

    // (3 + 4u)(1 + 2u) = 3 + 10u + 8u^2 = -5 + 10u, since u^2 = -1.
    const prod = a.mul(b);
    try std.testing.expect(prod.c0.equivalent(Fp.fromInt(3).sub(Fp.fromInt(8))));
    try std.testing.expect(prod.c1.equivalent(Fp.fromInt(10)));
}

test "fp2 squaring" {
    const a = Fp2.fromInts(3, 4);

    const sq = a.sq();
    const expected = a.mul(a);
    try std.testing.expect(sq.equivalent(expected));
}

test "fp2 inversion" {
    const a = Fp2.fromInts(3, 4);
    const inv_a = a.invert();
    try std.testing.expect(a.mul(inv_a).equivalent(Fp2.one));
}

test "fp2 conjugate" {
    const a = Fp2.fromInts(3, 4);
    const conj = a.conjugate();
    try std.testing.expect(conj.c0.equivalent(a.c0));
    try std.testing.expect(conj.c1.equivalent(a.c1.neg()));
}

test "fp2 norm" {
    const a = Fp2.fromInts(3, 4);
    // norm = 3^2 + 4^2 = 25
    const n = a.norm();
    try std.testing.expect(n.equivalent(Fp.fromInt(25)));
}

test "fp2 mul by nonresidue" {
    const a = Fp2.fromInts(3, 4);
    // (3 + 4u)(2 + u) = 6 + 11u + 4u^2 = 2 + 11u, since u^2 = -1.
    const result = a.mulByNonresidue();
    try std.testing.expect(result.c0.equivalent(Fp.fromInt(2)));
    try std.testing.expect(result.c1.equivalent(Fp.fromInt(11)));
}
