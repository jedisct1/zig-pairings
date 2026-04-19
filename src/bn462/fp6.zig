//! Sextic extension field Fp6 = Fp2[v] / (v^3 - (2 + u)) for BN462.
//!
//! Elements are represented as c0 + c1*v + c2*v^2 where c0, c1, c2 ∈ Fp2.
//!
//! The non-residue is ξ = 2 + u ∈ Fp2.

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
                .c1 = c1.mul(frobenius_coeffs_3[0]),
                .c2 = c2.mul(frobenius_coeffs_3[1]),
            },
            4 => .{
                .c0 = c0,
                .c1 = c1.mul(frobenius_coeffs_4[0]),
                .c2 = c2.mul(frobenius_coeffs_4[1]),
            },
            5 => .{
                .c0 = c0,
                .c1 = c1.mul(frobenius_coeffs_5[0]),
                .c2 = c2.mul(frobenius_coeffs_5[1]),
            },
            else => unreachable,
        };
    }

    /// Create from Fp2 components.
    pub fn fromFp2(c0: Fp2, c1: Fp2, c2: Fp2) Fp6 {
        return .{ .c0 = c0, .c1 = c1, .c2 = c2 };
    }
};

// Frobenius coefficients for power 1: ξ^((p - 1) / 3) and ξ^((2(p - 1)) / 3), where ξ = 2 + u.
const frobenius_coeffs_1: [2]Fp2 = .{
    // ξ^((p - 1) / 3)
    Fp2{
        .c0 = Fp{ .limbs = .{
            0x14a56e51d493bcc9,
            0xb84358efb8503d8f,
            0x39b4a0b18436f576,
            0xe4aaf59436d2c30f,
            0xe091597b24ae8bb1,
            0xe7f3df49c0b63234,
            0xdc08285e422ff3d9,
            0x0000000000000d82,
        } },
        .c1 = Fp{ .limbs = .{
            0x65925b95a8a1f1c5,
            0xd7e63ef34893b931,
            0xf9eb470fcaf539e4,
            0x3c19efb402e70760,
            0xeb0f3bbef193cd74,
            0x64ddfca93ddb64e3,
            0xebc1f798b984f2e2,
            0x0000000000000fa4,
        } },
    },
    // ξ^((2(p - 1)) / 3)
    Fp2{
        .c0 = Fp{ .limbs = .{
            0xfb09fc06e7668e9c,
            0x9b1d75f7f522952e,
            0x3df9742f18207d9c,
            0x39eeb2d4241bf7f0,
            0x434639e8c5a35e62,
            0xcd01a43cc00cec1e,
            0x40ad3fed4d425a3c,
            0x0000000000000b64,
        } },
        .c1 = Fp{ .limbs = .{
            0x8c921cfa68cca674,
            0x7f5d2d647e455f04,
            0xab3bd277f5c827fb,
            0x0e1c4458547c6693,
            0x90e84453848de2dd,
            0x4a85ae502ca1bea1,
            0x9de40ca8a03a0dcf,
            0x00000000000009e4,
        } },
    },
};

// Frobenius coefficients for power 2: ξ^((p^2 - 1) / 3) and ξ^((2(p^2 - 1)) / 3).
const frobenius_coeffs_2: [2]Fp2 = .{
    // ξ^((p^2 - 1) / 3)
    Fp2{
        .c0 = Fp{ .limbs = .{
            0x1411710cf60202c8,
            0x66e69ee38f2f8e88,
            0x6280de9a92a4e382,
            0xb29d7462344c7831,
            0x7f965ccd9ee433ba,
            0x0087fb5a1511f6b7,
            0x318cc43c76021fef,
            0x0000000000000078,
        } },
        .c1 = Fp{ .limbs = .{
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
        } },
    },
    // ξ^((2(p^2 - 1)) / 3)
    Fp2{
        .c0 = Fp{ .limbs = .{
            0xd1de50472b230ca3,
            0x60b71948720b0b58,
            0x6c2c1347370f979a,
            0xf724644d88ea1b8e,
            0x43036734dbadcda9,
            0xffc7fd1e777f3f1c,
            0x50c14b636e3f2006,
            0x00000000000013a6,
        } },
        .c1 = Fp{ .limbs = .{
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
        } },
    },
};

// Frobenius coefficients for power 3: ξ^((p^3 - 1) / 3) and ξ^((2(p^3 - 1)) / 3).
const frobenius_coeffs_3: [2]Fp2 = .{
    // ξ^((p^3 - 1) / 3)
    Fp2{
        .c0 = Fp{ .limbs = .{
            0xa063be75a7a70c71,
            0x19924d7f199f4e8b,
            0x5d4bf52831b210ee,
            0x0b497580ad80cd70,
            0x514f20d59b9424b4,
            0xa57612d18e8c8891,
            0xaec25573dd56adda,
            0x00000000000001f9,
        } },
        .c1 = Fp{ .limbs = .{
            0x6e01d90cb5d41773,
            0x82a2cc75c74b2a5d,
            0x04ccfa89db55a622,
            0x099007c5b97ba2eb,
            0x6e55e86095aa6d3f,
            0x057381b6b4222a43,
            0x0d49773d7fd2f731,
            0x0000000000001c38,
        } },
    },
    // ξ^((2(p^3 - 1)) / 3)
    Fp2{
        .c0 = Fp{ .limbs = .{
            0xf908db11799adb2e,
            0x6a2a4dfad656a6d5,
            0xcdbccaa0d426e816,
            0x11afcf8a688608ed,
            0x651b26994c03bf78,
            0x921d9a8438c94e61,
            0xab3ca1cded9808f5,
            0x0000000000000ea0,
        } },
        .c1 = Fp{ .limbs = .{
            0xfb4a623c19a834da,
            0x2d4ad72d40b9daf1,
            0xe5aea58420d85803,
            0x03cc2c97678409e5,
            0xd623f3e7f689a1ef,
            0xdec29a1cfba4fe09,
            0x86080ed19e15ad66,
            0x0000000000001181,
        } },
    },
};

// Frobenius coefficients for power 4: ξ^((p^4 - 1) / 3) and ξ^((2(p^4 - 1)) / 3).
const frobenius_coeffs_4: [2]Fp2 = .{
    // ξ^((p^4 - 1) / 3)
    Fp2{
        .c0 = Fp{ .limbs = .{
            0xd1de50472b230ca3,
            0x60b71948720b0b58,
            0x6c2c1347370f979a,
            0xf724644d88ea1b8e,
            0x43036734dbadcda9,
            0xffc7fd1e777f3f1c,
            0x50c14b636e3f2006,
            0x00000000000013a6,
        } },
        .c1 = Fp{ .limbs = .{
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
        } },
    },
    // ξ^((2(p^4 - 1)) / 3)
    Fp2{
        .c0 = Fp{ .limbs = .{
            0x1411710cf60202c8,
            0x66e69ee38f2f8e88,
            0x6280de9a92a4e382,
            0xb29d7462344c7831,
            0x7f965ccd9ee433ba,
            0x0087fb5a1511f6b7,
            0x318cc43c76021fef,
            0x0000000000000078,
        } },
        .c1 = Fp{ .limbs = .{
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
        } },
    },
};

// Frobenius coefficients for power 5: ξ^((p^5 - 1) / 3) and ξ^((2(p^5 - 1)) / 3).
const frobenius_coeffs_5: [2]Fp2 = .{
    // ξ^((p^5 - 1) / 3)
    Fp2{
        .c0 = Fp{ .limbs = .{
            0x6ef88340c3d8b6d9,
            0x24b24fd12e1073e5,
            0x68ff6a1d497df60b,
            0x229c242ce3ae7f80,
            0x8de985af3fbd5072,
            0x729604e3bdb3fd13,
            0xf56b834de2b95e4b,
            0x0000000000001487,
        } },
        .c1 = Fp{ .limbs = .{
            0x746f2b6e21b0f6ee,
            0x9286e116f0211c71,
            0x0147be54588318d8,
            0xdf772709d3a175b4,
            0x262edbe078c1c6fc,
            0x95ae6f9e27efe08c,
            0x07609369cb2815ec,
            0x0000000000001c2c,
        } },
    },
    // ξ^((2(p^5 - 1)) / 3)
    Fp2{
        .c0 = Fp{ .limbs = .{
            0x2feed8efdf121649,
            0xf140324d3486c3fb,
            0xf449c127131f96bc,
            0xc6f20ce33b600f21,
            0x17689f7dee58e2fd,
            0xa0e0b83e14207d5a,
            0x944c1f64c7659ccd,
            0x00000000000009ff,
        } },
        .c1 = Fp{ .limbs = .{
            0x9c2530d1bd9ea4c5,
            0x49dff1ae4100c609,
            0x6f1587fae8c67c71,
            0x00a81e520c019f86,
            0x58bdc7c484e87c0c,
            0xd6b7ae91e4affb2e,
            0x5c49e5a5c3f044c9,
            0x000000000000089e,
        } },
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
