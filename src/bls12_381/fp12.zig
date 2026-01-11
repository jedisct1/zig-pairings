//! Dodecic extension field Fp12 = Fp6[w] / (w^2 - v) for BLS12-381.
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

    /// Zero.
    pub const zero: Fp12 = .{ .c0 = Fp6.zero, .c1 = Fp6.zero };

    /// One.
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

        const c1_scaled = c1.mulByFp2(frobenius_coeffs[0]);

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

    /// Squaring specialized for the cyclotomic subgroup, which is cheaper than
    /// the general one.
    ///
    /// The caller is responsible for only passing elements of that subgroup.
    ///
    /// Follows algorithm 5.5.4 of the Guide to Pairing-Based Cryptography.
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

        // Undo the z-indexing.
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
    ///
    /// The exponent is applied directly rather than through an addition chain,
    /// which keeps the code short at the cost of some speed.
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

// Frobenius coefficients: (1+u)^((p^i - 1) / 6) for i = 1..5
const frobenius_coeffs: [5]Fp2 = .{
    // (1+u)^((p - 1) / 6)
    Fp2{
        .c0 = Fp{
            .limbs = .{
                0x07089552b319d465,
                0xc6695f92b50a8313,
                0x97e83cccd117228f,
                0xa35baecab2dc29ee,
                0x1ce393ea5daace4d,
                0x08f2220fb0fb66eb,
            },
        },
        .c1 = Fp{
            .limbs = .{
                0xb2f66aad4ce5d646,
                0x5842a06bfc497cec,
                0xcf4895d42599d394,
                0xc11b9cba40a8e8d0,
                0x2e3813cbe5a0de89,
                0x110eefda88847faf,
            },
        },
    },
    // (1+u)^((p^2 - 1) / 6)
    Fp2{
        .c0 = Fp{
            .limbs = .{
                0xecfb361b798dba3a,
                0xc100ddb891865a2c,
                0x0ec08ff1232bda8e,
                0xd5c13cc6f1ca4721,
                0x47222a47bf7b5c04,
                0x0110f184e51c5f59,
            },
        },
        .c1 = Fp.zero,
    },
    // (1+u)^((p^3 - 1) / 6)
    Fp2{
        .c0 = Fp{
            .limbs = .{
                0x3e2f585da55c9ad1,
                0x4294213d86c18183,
                0x382844c88b623732,
                0x92ad2afd19103e18,
                0x1d794e4fac7cf0b9,
                0x0bd592fc7d825ec8,
            },
        },
        .c1 = Fp{
            .limbs = .{
                0x7bcfa7a25aa30fda,
                0xdc17dec12a927e7c,
                0x2f088dd86b4ebef1,
                0xd1ca2087da74d4a7,
                0x2da2596696cebc1d,
                0x0e2b7eedbbfd87d2,
            },
        },
    },
    // (1+u)^((p^4 - 1) / 6)
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
    // (1+u)^((p^5 - 1) / 6)
    Fp2{
        .c0 = Fp{
            .limbs = .{
                0x3726c30af242c66c,
                0x7c2ac1aad1b6fe70,
                0xa04007fbba4b14a2,
                0xef517c3266341429,
                0x0095ba654ed2226b,
                0x02e370eccc86f7dd,
            },
        },
        .c1 = Fp{
            .limbs = .{
                0x82d83cf50dbce43f,
                0xa2813e53df9d018f,
                0xc6f0caa53c65e181,
                0x7525cf528d50fe95,
                0x4a85ed50f4798a6b,
                0x171da0fd6cf8eebd,
            },
        },
    },
};

// (p^12 - 1) / r for BLS12-381, represented as little-endian u64 limbs.
const final_exponent = [_]u64{
    0xc0bcb9b55df57510,
    0x25f98630e68bfb24,
    0x4406fbc8fbd5f489,
    0x8e2f8491d12191a0,
    0x3e9d71650a6f8069,
    0x226c2f011d4cab80,
    0x67f67c4717489119,
    0xaf3f881bd88592d7,
    0x1a67e49eeed2161d,
    0xe5b78c7869aeb218,
    0xf6539314043f7bbc,
    0x73f62537f2701aae,
    0xaff1c910e9622d2a,
    0x6283313492caa9d4,
    0x2e2f3ec2bea83d19,
    0xa4c7e79fb02faa73,
    0x6c49637fd7961be1,
    0x08e88adce8817745,
    0x35de3f7a36399917,
    0x9c1d9f7c31759c36,
    0xfa9e13c24ea820b0,
    0x3fc56947a403577d,
    0xa4c1b6dcfc5cceb7,
    0x1bbd81367066bca6,
    0x0418a3ef0bc62775,
    0x49bf9b71a9f9e010,
    0x511291097db60b17,
    0x498345c6e5308f1c,
    0x6d8823b19dadd7c2,
    0x92004cedd556952c,
    0x4c6bec3ec03ef195,
    0x0a1fad20044ce6ad,
    0xc55d3109cd15948d,
    0x334f46c02c3f0bd0,
    0x3b5a62eb34c05739,
    0x724538411d1676a5,
    0x127a1b5ad0463434,
    0x61a474c5c85b0129,
    0x8dfc8e2886ef965e,
    0x96532fef459f1243,
    0x40ee7169cdc10412,
    0x9c40a68eb74bb22a,
    0x25118790f4684d0b,
    0x596bc293c8d4c01f,
    0x1064837f27611212,
    0x077ffb10bf24dde4,
    0xc49f570bcd2b01f3,
    0x1a0c5bf24c374693,
    0x350da5359bc73ab6,
    0xd2670d93e4d7acdd,
    0xd39099b86e1ab656,
    0x19328148978e2b0d,
    0xb113f414386b0e88,
    0x07a0dce2630d9aa4,
    0xa927e7bb93753318,
    0xe347aa68ad49466f,
    0x1c0ad0d6106feaf4,
    0xc872ee83ff3a0f0f,
    0x074e43b9a660835c,
    0xc0aadff5e9cfee9a,
    0x30698e8cc7deada9,
    0xd1073776ab353f2c,
    0x17848517badc3a43,
    0x7363baa13f8d14a9,
    0xd4977b3f7d4507d0,
    0x496a1c0a89ee0193,
    0xdcc825b7e1bda9c0,
    0x0000000002ee1db5,
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

test "fp12 mulBy014 sparse multiplication" {
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

    // Sparse element with c0, c1 at positions 0 and 1, and c4 at position 4.
    const c0 = Fp2.fromInts(13, 14);
    const c1 = Fp2.fromInts(15, 16);
    const c4 = Fp2.fromInts(17, 18);

    const sparse = Fp12{
        .c0 = Fp6{ .c0 = c0, .c1 = c1, .c2 = Fp2.zero },
        .c1 = Fp6{ .c0 = Fp2.zero, .c1 = c4, .c2 = Fp2.zero },
    };

    const result_sparse = a.mulBy014(c0, c1, c4);
    const result_full = a.mul(sparse);

    try std.testing.expect(result_sparse.equivalent(result_full));
}

test "conjugate equals frobenius_6" {
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
    const frob6_chained = f.frobeniusMap(1).frobeniusMap(1).frobeniusMap(1).frobeniusMap(1).frobeniusMap(1).frobeniusMap(1);
    const frob6_direct = f.frobeniusMap(6);

    try std.testing.expect(conj.equivalent(frob6_chained));
    try std.testing.expect(conj.equivalent(frob6_direct));
}

test "frobenius consistency" {
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

    const frob2_direct = f.frobeniusMap(2);
    const frob2_chained = f.frobeniusMap(1).frobeniusMap(1);

    try std.testing.expect(frob2_direct.equivalent(frob2_chained));

    const frob3_direct = f.frobeniusMap(3);
    const frob3_chained = f.frobeniusMap(1).frobeniusMap(1).frobeniusMap(1);
    try std.testing.expect(frob3_direct.equivalent(frob3_chained));

    const frob4_direct = f.frobeniusMap(4);
    const frob4_chained = f.frobeniusMap(2).frobeniusMap(2);
    try std.testing.expect(frob4_direct.equivalent(frob4_chained));
}

test "fp4Square correctness" {
    const a = Fp2.fromInts(3, 5);
    const b = Fp2.fromInts(7, 11);

    const result = Fp12.fp4Square(a, b);

    const a_sq = a.sq();
    const b_sq = b.sq();
    const xi_b_sq = b_sq.mulByNonresidue();
    const c0_expected = a_sq.add(xi_b_sq);
    const c1_expected = a.mul(b).dbl();

    try std.testing.expect(result[0].equivalent(c0_expected));
    try std.testing.expect(result[1].equivalent(c1_expected));
}

test "cyclotomic squaring consistency" {
    // For elements in the cyclotomic subgroup, cyclotomic squaring must equal regular squaring.
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

    // Map into the cyclotomic subgroup via the easy part of the final exponent.
    var t0 = f.conjugate();
    var t1 = f.invert();
    var t2 = t0.mul(t1);

    t1 = t2.frobeniusMap(2);
    const in_cyclo = t1.mul(t2);

    // In the cyclotomic subgroup, f * conjugate(f) = 1.
    const norm = in_cyclo.mul(in_cyclo.conjugate());
    try std.testing.expect(norm.isOne());

    // Independently recompute the expected square from the Granger-Scott formulas.
    const c0c0 = in_cyclo.c0.c0;
    const c0c1 = in_cyclo.c0.c1;
    const c0c2 = in_cyclo.c0.c2;
    const c1c0 = in_cyclo.c1.c0;
    const c1c1 = in_cyclo.c1.c1;
    const c1c2 = in_cyclo.c1.c2;

    const sq1 = Fp12.fp4Square(c0c0, c1c1);
    const t3 = sq1[0];
    const t4 = sq1[1];

    const sq2 = Fp12.fp4Square(c1c0, c0c2);
    const t5 = sq2[0];
    const t6 = sq2[1];

    const sq3 = Fp12.fp4Square(c0c1, c1c2);
    const t7 = sq3[0];
    const t8 = sq3[1];

    const t9 = t8.mulByNonresidue();

    const expected_c0_c0 = t3.sub(c0c0).dbl().add(t3);
    const expected_c0_c1 = t5.sub(c0c1).dbl().add(t5);
    const expected_c0_c2 = t7.sub(c0c2).dbl().add(t7);
    const expected_c1_c0 = t9.add(c1c0).dbl().add(t9);
    const expected_c1_c1 = t4.add(c1c1).dbl().add(t4);
    const expected_c1_c2 = t6.add(c1c2).dbl().add(t6);

    const expected_sq = Fp12{
        .c0 = Fp6{ .c0 = expected_c0_c0, .c1 = expected_c0_c1, .c2 = expected_c0_c2 },
        .c1 = Fp6{ .c0 = expected_c1_c0, .c1 = expected_c1_c1, .c2 = expected_c1_c2 },
    };

    const cyclo_sq = in_cyclo.cyclotomicSquare();
    const regular_sq = in_cyclo.sq();

    try std.testing.expect(cyclo_sq.equivalent(expected_sq));
    try std.testing.expect(expected_sq.equivalent(regular_sq));
    try std.testing.expect(cyclo_sq.equivalent(regular_sq));
}
