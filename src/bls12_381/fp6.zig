//! Sextic extension field Fp6 = Fp2[v] / (v^3 - (1 + u)) for BLS12-381.
//!
//! Elements are represented as c0 + c1*v + c2*v^2 where c0, c1, c2 ∈ Fp2.
//!
//! The non-residue is ξ = 1 + u ∈ Fp2.

const std = @import("std");
const Fp = @import("fp.zig").Fp;
const Fp2 = @import("fp2.zig").Fp2;

/// Element of the sextic extension field Fp6.
pub const Fp6 = struct {
    c0: Fp2,
    c1: Fp2,
    c2: Fp2,

    /// Zero.
    pub const zero: Fp6 = .{ .c0 = Fp2.zero, .c1 = Fp2.zero, .c2 = Fp2.zero };

    /// One.
    pub const one: Fp6 = .{ .c0 = Fp2.one, .c1 = Fp2.zero, .c2 = Fp2.zero };

    /// Return true if the element is zero.
    pub fn isZero(fe: Fp6) bool {
        return fe.c0.isZero() and fe.c1.isZero() and fe.c2.isZero();
    }

    /// Return true if both elements are equivalent.
    pub fn equivalent(a: Fp6, b: Fp6) bool {
        return a.c0.equivalent(b.c0) and a.c1.equivalent(b.c1) and a.c2.equivalent(b.c2);
    }

    /// Conditionally replace the element with `other` when `choice` is 1 (constant-time).
    pub fn cMov(fe: *Fp6, other: Fp6, choice: u1) void {
        fe.c0.cMov(other.c0, choice);
        fe.c1.cMov(other.c1, choice);
        fe.c2.cMov(other.c2, choice);
    }

    /// Negate the element.
    pub fn neg(a: Fp6) Fp6 {
        return .{
            .c0 = a.c0.neg(),
            .c1 = a.c1.neg(),
            .c2 = a.c2.neg(),
        };
    }

    /// Add two Fp6 elements.
    pub fn add(a: Fp6, b: Fp6) Fp6 {
        return .{
            .c0 = a.c0.add(b.c0),
            .c1 = a.c1.add(b.c1),
            .c2 = a.c2.add(b.c2),
        };
    }

    /// Subtract two Fp6 elements.
    pub fn sub(a: Fp6, b: Fp6) Fp6 {
        return .{
            .c0 = a.c0.sub(b.c0),
            .c1 = a.c1.sub(b.c1),
            .c2 = a.c2.sub(b.c2),
        };
    }

    /// Double an Fp6 element.
    pub fn dbl(a: Fp6) Fp6 {
        return .{
            .c0 = a.c0.dbl(),
            .c1 = a.c1.dbl(),
            .c2 = a.c2.dbl(),
        };
    }

    /// Multiply two Fp6 elements, using three products instead of nine.
    pub fn mul(a: Fp6, b: Fp6) Fp6 {
        const a0b0 = a.c0.mul(b.c0);
        const a1b1 = a.c1.mul(b.c1);
        const a2b2 = a.c2.mul(b.c2);

        const t0 = a.c1.add(a.c2).mul(b.c1.add(b.c2)).sub(a1b1).sub(a2b2);
        const c0 = a0b0.add(t0.mulByNonresidue());

        const t1 = a.c0.add(a.c1).mul(b.c0.add(b.c1)).sub(a0b0).sub(a1b1);
        const c1 = t1.add(a2b2.mulByNonresidue());

        const t2 = a.c0.add(a.c2).mul(b.c0.add(b.c2)).sub(a0b0).sub(a2b2);
        const c2 = t2.add(a1b1);

        return .{ .c0 = c0, .c1 = c1, .c2 = c2 };
    }

    /// Square an Fp6 element.
    pub fn sq(a: Fp6) Fp6 {
        const s0 = a.c0.sq();
        const ab = a.c0.mul(a.c1);
        const s1 = ab.dbl();
        const s2 = a.c0.sub(a.c1).add(a.c2).sq();
        const bc = a.c1.mul(a.c2);
        const s3 = bc.dbl();
        const s4 = a.c2.sq();

        const c0 = s0.add(s3.mulByNonresidue());
        const c1 = s1.add(s4.mulByNonresidue());
        const c2 = s1.add(s2).add(s3).sub(s0).sub(s4);

        return .{ .c0 = c0, .c1 = c1, .c2 = c2 };
    }

    /// Multiply by a scalar in Fp2.
    pub fn mulByFp2(a: Fp6, scalar: Fp2) Fp6 {
        return .{
            .c0 = a.c0.mul(scalar),
            .c1 = a.c1.mul(scalar),
            .c2 = a.c2.mul(scalar),
        };
    }

    /// Multiply by v, which rotates the coefficients.
    pub fn mulByV(a: Fp6) Fp6 {
        return .{
            .c0 = a.c2.mulByNonresidue(),
            .c1 = a.c0,
            .c2 = a.c1,
        };
    }

    /// Sparse multiply by (c0 + c1*v).
    pub fn mulBy01(a: Fp6, c0: Fp2, c1: Fp2) Fp6 {
        const a_a = a.c0.mul(c0);
        const b_b = a.c1.mul(c1);

        const new_c0 = a.c2.mul(c1).mulByNonresidue().add(a_a);
        const new_c1 = c0.add(c1).mul(a.c0.add(a.c1)).sub(a_a).sub(b_b);
        const new_c2 = a.c2.mul(c0).add(b_b);

        return .{ .c0 = new_c0, .c1 = new_c1, .c2 = new_c2 };
    }

    /// Sparse multiply by (c1*v).
    pub fn mulBy1(a: Fp6, c1: Fp2) Fp6 {
        return .{
            .c0 = a.c2.mul(c1).mulByNonresidue(),
            .c1 = a.c0.mul(c1),
            .c2 = a.c1.mul(c1),
        };
    }

    /// Multiplicative inverse.
    ///
    /// Follows the cubic extension inversion of https://eprint.iacr.org/2010/354.pdf,
    /// which reduces the problem to a single Fp2 inversion.
    pub fn invert(a: Fp6) Fp6 {
        const c0_sq = a.c0.sq();
        const c1_sq = a.c1.sq();
        const c2_sq = a.c2.sq();

        const c0c1 = a.c0.mul(a.c1);
        const c0c2 = a.c0.mul(a.c2);
        const c1c2 = a.c1.mul(a.c2);

        const t0 = c0_sq.sub(c1c2.mulByNonresidue());
        const t1 = c2_sq.mulByNonresidue().sub(c0c1);
        const t2 = c1_sq.sub(c0c2);

        const t3 = a.c1.mul(t2).add(a.c2.mul(t1)).mulByNonresidue();

        // The norm, which is the only value that has to be inverted.
        const t4 = a.c0.mul(t0).add(t3);
        const t4_inv = t4.invert();

        return .{
            .c0 = t0.mul(t4_inv),
            .c1 = t1.mul(t4_inv),
            .c2 = t2.mul(t4_inv),
        };
    }

    /// Frobenius map (p-th power).
    pub fn frobeniusMap(a: Fp6, power: u8) Fp6 {
        const c0 = a.c0.frobeniusMap(power);
        const c1 = a.c1.frobeniusMap(power);
        const c2 = a.c2.frobeniusMap(power);

        // Raising to the p-th power leaves each coefficient scaled by a twist
        // factor that depends only on how many times the map is applied.
        return switch (power % 6) {
            0 => .{ .c0 = c0, .c1 = c1, .c2 = c2 },
            1 => .{
                .c0 = c0,
                .c1 = c1.mul(frobenius_coeffs_1[0]),
                .c2 = c2.mul(frobenius_coeffs_1[1]),
            },
            2 => .{
                .c0 = c0,
                .c1 = c1.mul(frobenius_coeffs_2[0]),
                .c2 = c2.mul(frobenius_coeffs_2[1]),
            },
            3 => .{
                .c0 = c0,
                .c1 = c1.neg(),
                .c2 = c2,
            },
            4 => .{
                .c0 = c0,
                .c1 = c1.mul(frobenius_coeffs_1[0].conjugate()),
                .c2 = c2.mul(frobenius_coeffs_1[1].conjugate()),
            },
            5 => .{
                .c0 = c0,
                .c1 = c1.mul(frobenius_coeffs_2[0].conjugate()),
                .c2 = c2.mul(frobenius_coeffs_2[1].conjugate()),
            },
            else => unreachable,
        };
    }

    /// Create from Fp2 components.
    pub fn fromFp2(c0: Fp2, c1: Fp2, c2: Fp2) Fp6 {
        return .{ .c0 = c0, .c1 = c1, .c2 = c2 };
    }
};

// Frobenius coefficients for power 1: ξ^((p - 1) / 3) and ξ^((2(p - 1)) / 3).
const frobenius_coeffs_1: [2]Fp2 = .{
    // ξ^((p - 1) / 3)
    Fp2{
        .c0 = Fp.zero,
        .c1 = Fp{
            .limbs = .{
                0xcd03c9e48671f071,
                0x5dab22461fcda5d2,
                0x587042afd3851b95,
                0x8eb60ebe01bacb9e,
                0x03f97d6e83d050d2,
                0x18f0206554638741,
            },
        },
    },
    // ξ^((2(p - 1)) / 3)
    Fp2{
        .c0 = Fp{
            .limbs = .{
                0x890dc9e4867545c3,
                0x2af322533285a5d5,
                0x50880866309b7e2c,
                0xa20d1b8c7e881024,
                0x14e4f04fe2db9068,
                0x14e56d3f1564853a,
            },
        },
        .c1 = Fp.zero,
    },
};

// Frobenius coefficients for power 2.
const frobenius_coeffs_2: [2]Fp2 = .{
    // ξ^((p^2 - 1) / 3)
    Fp2{
        .c0 = Fp{
            .limbs = .{
                0x30f1361b798a64e8,
                0xf3b8ddab7ece5a2a,
                0x16a8ca3ac61577f7,
                0xc26a2ff874fd029b,
                0x3636b76660701c6e,
                0x051ba4ab241b6160,
            },
        },
        .c1 = Fp.zero,
    },
    // ξ^((2(p^2 - 1)) / 3)
    Fp2{
        .c0 = Fp{
            .limbs = .{
                0xcd03c9e48671f071,
                0x5dab22461fcda5d2,
                0x587042afd3851b95,
                0x8eb60ebe01bacb9e,
                0x03f97d6e83d050d2,
                0x18f0206554638741,
            },
        },
        .c1 = Fp.zero,
    },
};

test "fp6 basic arithmetic" {
    const a = Fp6{
        .c0 = Fp2.fromInts(1, 2),
        .c1 = Fp2.fromInts(3, 4),
        .c2 = Fp2.fromInts(5, 6),
    };
    const b = Fp6{
        .c0 = Fp2.fromInts(7, 8),
        .c1 = Fp2.fromInts(9, 10),
        .c2 = Fp2.fromInts(11, 12),
    };

    const sum = a.add(b);
    const diff = sum.sub(b);
    try std.testing.expect(diff.equivalent(a));

    const neg_a = a.neg();
    try std.testing.expect(a.add(neg_a).isZero());
}

test "fp6 multiplication" {
    const a = Fp6{
        .c0 = Fp2.fromInts(1, 0),
        .c1 = Fp2.fromInts(0, 0),
        .c2 = Fp2.fromInts(0, 0),
    };

    const result = a.mul(Fp6.one);
    try std.testing.expect(result.equivalent(a));
}

test "fp6 squaring" {
    const a = Fp6{
        .c0 = Fp2.fromInts(1, 2),
        .c1 = Fp2.fromInts(3, 4),
        .c2 = Fp2.fromInts(5, 6),
    };

    const sq = a.sq();
    const prod = a.mul(a);
    try std.testing.expect(sq.equivalent(prod));
}

test "fp6 inversion" {
    const a = Fp6{
        .c0 = Fp2.fromInts(1, 2),
        .c1 = Fp2.fromInts(3, 4),
        .c2 = Fp2.fromInts(5, 6),
    };

    const inv_a = a.invert();
    const product = a.mul(inv_a);
    try std.testing.expect(product.equivalent(Fp6.one));
}
