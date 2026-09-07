//! Dodecic extension field Fp12 = Fp6[w] / (w^2 - v) for BN462.
//!
//! Elements are represented as c0 + c1*w where c0, c1 ∈ Fp6.
//!
//! This is the target group GT for the pairing.

const std = @import("std");
const Fp = @import("fp.zig").Fp;
const Fp2 = @import("fp2.zig").Fp2;
const Fp6 = @import("fp6.zig").Fp6;

/// Element of the dodecic extension field Fp12.
pub const Fp12 = struct {
    c0: Fp6,
    c1: Fp6,

    /// Zero element.
    pub const zero: Fp12 = .{ .c0 = Fp6.zero, .c1 = Fp6.zero };

    /// One element.
    pub const one: Fp12 = .{ .c0 = Fp6.one, .c1 = Fp6.zero };

    /// Return true if the element is zero.
    pub fn isZero(fe: Fp12) bool {
        return fe.c0.isZero() and fe.c1.isZero();
    }

    /// Return true if the element is one.
    pub fn isOne(fe: Fp12) bool {
        return fe.equivalent(one);
    }

    /// Return true if both elements are equivalent.
    pub fn equivalent(a: Fp12, b: Fp12) bool {
        return a.c0.equivalent(b.c0) and a.c1.equivalent(b.c1);
    }

    /// Conditionally replace the element with `other` when `choice` is 1 (constant-time).
    pub fn cMov(fe: *Fp12, other: Fp12, choice: u1) void {
        fe.c0.cMov(other.c0, choice);
        fe.c1.cMov(other.c1, choice);
    }

    /// Negate the element.
    pub fn neg(a: Fp12) Fp12 {
        return .{
            .c0 = a.c0.neg(),
            .c1 = a.c1.neg(),
        };
    }

    /// Add two Fp12 elements.
    pub fn add(a: Fp12, b: Fp12) Fp12 {
        return .{
            .c0 = a.c0.add(b.c0),
            .c1 = a.c1.add(b.c1),
        };
    }

    /// Subtract two Fp12 elements.
    pub fn sub(a: Fp12, b: Fp12) Fp12 {
        return .{
            .c0 = a.c0.sub(b.c0),
            .c1 = a.c1.sub(b.c1),
        };
    }

    /// Double an Fp12 element.
    pub fn dbl(a: Fp12) Fp12 {
        return .{
            .c0 = a.c0.dbl(),
            .c1 = a.c1.dbl(),
        };
    }

    /// Multiply two Fp12 elements, trading one Fp6 product for extra additions.
    pub fn mul(a: Fp12, b: Fp12) Fp12 {
        const aa = a.c0.mul(b.c0);
        const bb = a.c1.mul(b.c1);

        const c0 = aa.add(bb.mulByV());

        const c1 = a.c0.add(a.c1).mul(b.c0.add(b.c1)).sub(aa).sub(bb);

        return .{ .c0 = c0, .c1 = c1 };
    }

    /// Square an Fp12 element.
    pub fn sq(a: Fp12) Fp12 {
        const ab = a.c0.mul(a.c1);

        const c0 = a.c0.add(a.c1).mul(a.c0.add(a.c1.mulByV())).sub(ab).sub(ab.mulByV());

        const c1 = ab.dbl();

        return .{ .c0 = c0, .c1 = c1 };
    }

    /// Compute the multiplicative inverse.
    pub fn invert(a: Fp12) Fp12 {
        const aa = a.c0.sq();
        const bb = a.c1.sq();
        const t = aa.sub(bb.mulByV());
        const t_inv = t.invert();

        return .{
            .c0 = a.c0.mul(t_inv),
            .c1 = a.c1.neg().mul(t_inv),
        };
    }

    /// Conjugate: (a + bw) -> (a - bw).
    pub fn conjugate(a: Fp12) Fp12 {
        return .{
            .c0 = a.c0,
            .c1 = a.c1.neg(),
        };
    }

    /// Frobenius map (p-th power).
    ///
    /// For power > 1, this is equivalent to applying frobeniusMap(1) repeatedly.
    pub fn frobeniusMap(a: Fp12, power: u8) Fp12 {
        if (power == 0) return a;

        var result = a;
        var i: u8 = 0;
        while (i < power % 12) : (i += 1) {
            result = result.frobeniusMapOnce();
        }
        return result;
    }

    /// Single Frobenius map (p-th power, applies once).
    fn frobeniusMapOnce(a: Fp12) Fp12 {
        const c0 = a.c0.frobeniusMap(1);
        const c1 = a.c1.frobeniusMap(1);

        const c1_scaled = c1.mulByFp2(frobenius_coeff);

        return .{ .c0 = c0, .c1 = c1_scaled };
    }

    /// Multiply by a line evaluation, whose only non-zero coefficients sit at
    /// positions 0, 1 and 4.
    pub fn mulBy014(a: Fp12, c0: Fp2, c1: Fp2, c4: Fp2) Fp12 {
        const aa = a.c0.mulBy01(c0, c1);
        const bb = a.c1.mulBy1(c4);
        const o = c1.add(c4);

        const new_c1 = a.c0.add(a.c1).mulBy01(c0, o).sub(aa).sub(bb);
        const new_c0 = aa.add(bb.mulByV());

        return .{ .c0 = new_c0, .c1 = new_c1 };
    }

    /// Multiply by a line evaluation on a D-type twist, whose only non-zero
    /// coefficients sit at positions 0, 3 and 4.
    pub fn mulBy034(a: Fp12, c0: Fp2, c3: Fp2, c4: Fp2) Fp12 {
        const aa = Fp6{ .c0 = a.c0.c0.mul(c0), .c1 = a.c0.c1.mul(c0), .c2 = a.c0.c2.mul(c0) };

        const bb = a.c1.mulBy01(c3, c4);

        const sum01 = a.c0.add(a.c1);
        const sum_c = c0.add(c3);
        const cc = sum01.mulBy01(sum_c, c4).sub(aa).sub(bb);

        const new_c0 = aa.add(bb.mulByV());

        return .{ .c0 = new_c0, .c1 = cc };
    }

    /// Squaring specialized for the cyclotomic subgroup, which is cheaper than
    /// the general one.
    ///
    /// The caller is responsible for only passing elements of that subgroup.
    pub fn cyclotomicSquare(a: Fp12) Fp12 {
        // Variables follow the z-indexing of the Granger-Scott cyclotomic squaring.
        var z0 = a.c0.c0;
        var z4 = a.c0.c1;
        var z3 = a.c0.c2;
        var z2 = a.c1.c0;
        var z1 = a.c1.c1;
        var z5 = a.c1.c2;

        const sq01 = fp4Square(z0, z1);
        var t0 = sq01[0];
        var t1 = sq01[1];

        z0 = t0.sub(z0);
        z0 = z0.add(z0).add(t0);

        z1 = t1.add(z1);
        z1 = z1.add(z1).add(t1);

        const sq23 = fp4Square(z2, z3);
        t0 = sq23[0];
        t1 = sq23[1];

        const sq45 = fp4Square(z4, z5);
        const t2 = sq45[0];
        const t3 = sq45[1];

        z4 = t0.sub(z4);
        z4 = z4.add(z4).add(t0);

        z5 = t1.add(z5);
        z5 = z5.add(z5).add(t1);

        t0 = t3.mulByNonresidue();
        z2 = t0.add(z2);
        z2 = z2.add(z2).add(t0);

        z3 = t2.sub(z3);
        z3 = z3.add(z3).add(t2);

        return .{
            .c0 = Fp6{ .c0 = z0, .c1 = z4, .c2 = z3 },
            .c1 = Fp6{ .c0 = z2, .c1 = z1, .c2 = z5 },
        };
    }

    /// Square an element of the Fp4 subfield, given and returned as its two
    /// Fp2 coefficients.
    fn fp4Square(a: Fp2, b: Fp2) struct { Fp2, Fp2 } {
        const t0 = a.sq();
        const t1 = b.sq();
        var t2 = t1.mulByNonresidue();
        const c0 = t2.add(t0);
        t2 = a.add(b);
        t2 = t2.sq();
        t2 = t2.sub(t0);
        const c1 = t2.sub(t1);
        return .{ c0, c1 };
    }

    /// Cyclotomic square n times.
    pub fn cyclotomicSquareN(a: Fp12, comptime n: comptime_int) Fp12 {
        var result = a;
        inline for (0..n) |_| {
            result = result.cyclotomicSquare();
        }
        return result;
    }

    /// Map a Miller loop result into the target group.
    ///
    /// This is the last step of a pairing.
    pub fn finalExponentiation(a: Fp12) Fp12 {
        return powByU64s(a, final_exponent[0..]);
    }

    fn powByU64s(base: Fp12, exponent: []const u64) Fp12 {
        var result = Fp12.one;
        var started = false;

        var limb_index = @as(i32, @intCast(exponent.len)) - 1;
        while (limb_index >= 0) : (limb_index -= 1) {
            const idx: usize = @intCast(limb_index);
            const limb = exponent[idx];
            var bit: i32 = 63;
            while (bit >= 0) : (bit -= 1) {
                const shift: u6 = @intCast(bit);
                const bit_set = ((limb >> shift) & 1) == 1;
                if (started) {
                    result = result.sq();
                } else if (!bit_set) {
                    continue;
                } else {
                    started = true;
                }

                if (bit_set) {
                    result = result.mul(base);
                }
            }
        }

        return result;
    }
};

const final_exponent = [_]u64{
    0xf296a63df31ad6f0,
    0xe461008d7500ad76,
    0xe2c7aeca88350a32,
    0x9ba18c0c4036572b,
    0x28267977feea9b34,
    0xc503e4a9f3c930ac,
    0x240b2c4bf1dc65bd,
    0x84c8d8bdb397757a,
    0x5251e225e754bda2,
    0xebceb75e594cfeb6,
    0xc087319cfe0c3195,
    0x7f7c9da691fe2208,
    0x494e2b1e047a2d4f,
    0x1531303cb5e551a8,
    0x68c5d5932e20d336,
    0xc23e3daa995e1047,
    0xd0d27735d8a2dbc9,
    0xc15d4df8599ef935,
    0x3534d942391e532a,
    0x4695115c998701b9,
    0xf98a9ee579a68093,
    0x43270b2906eae400,
    0x04792054256f102d,
    0x500af047b0f0a8cf,
    0x49ee1311341a9f85,
    0xdda5e39a64841c0c,
    0x531ac1e065a54380,
    0x65e3f6d28928b00a,
    0x8d8e83c0912074b4,
    0xda967915b67ecb54,
    0x9705b34acea30559,
    0xb4f25e74a99258cc,
    0xea9dcd7cdebf143b,
    0x872f0b45a1b1c39c,
    0x43ca1632456b0e08,
    0x6b8535da27a07a98,
    0x4905da5d3ed98480,
    0x1e301891de9159c4,
    0xfdeaa3a2eea44e11,
    0xbe1edc3745941ad4,
    0x0c765c06a22d6c3d,
    0x7e5a5cd302d08a8e,
    0x629b5f363448348e,
    0x3e4f21e0156a80df,
    0x32949e3933e469d0,
    0x870963bc7fa14c94,
    0x5ecdad4ed7a9439d,
    0x4bb4419b257b3687,
    0xdd0e0916f2548b87,
    0xfeaf263bb4f0081f,
    0x56aec8492bcf2075,
    0x0f114f3461581c24,
    0xd4ce5c19a78ec4fa,
    0xc57e7fff585f4b96,
    0xc888f3087be49d7d,
    0xf17549d8d175c2e0,
    0xcc34c5b0c235e796,
    0xb374fead792bbf52,
    0x7a0c9c1a3dd5bceb,
    0x7115f7302352c8c8,
    0x98cde3be7613b775,
    0x54e084a9fa68d2d8,
    0xd3fb18bf1cf9fd6f,
    0xa30c5affd8b213e0,
    0x2fd46d38ae579a86,
    0xd576873bc9b1337d,
    0x34287495990a428d,
    0xba5bfd335074618b,
    0x7b1775d782c049a9,
    0x483a31618186742b,
    0xae4878e73f4549c7,
    0x21563cc20d9226f3,
    0x5a87bb8501eee7be,
    0x6e18270f8bbef741,
    0x877b92dda477933e,
    0x1d96fee7b6170adb,
    0xcfa7a0359e8f5d0e,
    0xcd9b4523224445cc,
    0xe6f7887368ef71ad,
    0x000000000001d621,
};

// Frobenius correction for the tower basis.
const frobenius_coeff = Fp2{
    .c0 = Fp{ .limbs = .{
        0xc0696fa7c6e23251,
        0x9b58e25cf0f1df0d,
        0xc14560073aab9c79,
        0x224b6a68dfdbbac9,
        0x7d997c2a5792eb25,
        0x5c1a4d614a2d223b,
        0xf1724aaf3af407ff,
        0x00000000000020d3,
    } },
    .c1 = Fp{ .limbs = .{
        0xc419bc0bf0aeeac1,
        0x352ddf0031b8e848,
        0x27838eb3a4c48e55,
        0xaa3c3bf6ff656e5e,
        0x5753783f47b9fd4d,
        0x9d2e7092606b4575,
        0x7ad74e2b650ad865,
        0x0000000000000d59,
    } },
};

test "fp12 basic arithmetic" {
    const a = Fp12{
        .c0 = Fp6{
            .c0 = Fp2.fromInts(1, 2),
            .c1 = Fp2.fromInts(3, 4),
            .c2 = Fp2.fromInts(5, 6),
        },
        .c1 = Fp6{
            .c0 = Fp2.fromInts(7, 8),
            .c1 = Fp2.fromInts(9, 10),
            .c2 = Fp2.fromInts(11, 12),
        },
    };

    const sum = a.add(a);
    const diff = sum.sub(a);
    try std.testing.expect(diff.equivalent(a));

    const neg_a = a.neg();
    try std.testing.expect(a.add(neg_a).isZero());
}

test "fp12 multiplication" {
    const a = Fp12{
        .c0 = Fp6{
            .c0 = Fp2.fromInts(1, 2),
            .c1 = Fp2.fromInts(3, 4),
            .c2 = Fp2.fromInts(5, 6),
        },
        .c1 = Fp6{
            .c0 = Fp2.fromInts(7, 8),
            .c1 = Fp2.fromInts(9, 10),
            .c2 = Fp2.fromInts(11, 12),
        },
    };

    const result = a.mul(Fp12.one);
    try std.testing.expect(result.equivalent(a));
}

test "fp12 squaring" {
    const a = Fp12{
        .c0 = Fp6{
            .c0 = Fp2.fromInts(1, 2),
            .c1 = Fp2.fromInts(3, 4),
            .c2 = Fp2.fromInts(5, 6),
        },
        .c1 = Fp6{
            .c0 = Fp2.fromInts(7, 8),
            .c1 = Fp2.fromInts(9, 10),
            .c2 = Fp2.fromInts(11, 12),
        },
    };

    const sq = a.sq();
    const prod = a.mul(a);
    try std.testing.expect(sq.equivalent(prod));
}

test "fp12 inversion" {
    const a = Fp12{
        .c0 = Fp6{
            .c0 = Fp2.fromInts(1, 2),
            .c1 = Fp2.fromInts(3, 4),
            .c2 = Fp2.fromInts(5, 6),
        },
        .c1 = Fp6{
            .c0 = Fp2.fromInts(7, 8),
            .c1 = Fp2.fromInts(9, 10),
            .c2 = Fp2.fromInts(11, 12),
        },
    };

    const inv_a = a.invert();
    const product = a.mul(inv_a);
    try std.testing.expect(product.equivalent(Fp12.one));
}

test "fp12 conjugate basic" {
    const f = Fp12{
        .c0 = Fp6{
            .c0 = Fp2.fromInts(1, 2),
            .c1 = Fp2.fromInts(3, 4),
            .c2 = Fp2.fromInts(5, 6),
        },
        .c1 = Fp6{
            .c0 = Fp2.fromInts(7, 8),
            .c1 = Fp2.fromInts(9, 10),
            .c2 = Fp2.fromInts(11, 12),
        },
    };

    const conj = f.conjugate();

    try std.testing.expect(conj.c0.equivalent(f.c0));
    try std.testing.expect(conj.c1.equivalent(f.c1.neg()));

    const double_conj = conj.conjugate();
    try std.testing.expect(double_conj.equivalent(f));
}
