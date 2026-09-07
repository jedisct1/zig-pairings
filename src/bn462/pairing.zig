//! Optimal Ate pairing for BN462.
//!
//! Computes the bilinear pairing e: G1 × G2 → GT.
//!
//! For BN curves, the optimal ate pairing uses the loop parameter 6x + 2.

const std = @import("std");
const Fp = @import("fp.zig").Fp;
const Fp2 = @import("fp2.zig").Fp2;
const Fp12 = @import("fp12.zig").Fp12;
const G1 = @import("g1.zig").G1;
const G1Affine = @import("g1.zig").AffineCoordinates;
const G2 = @import("g2.zig").G2;
const G2Affine = @import("g2.zig").AffineCoordinates;

/// Element of the target group GT (a subgroup of Fp12*).
pub const GT = Fp12;

/// Input to a pairing product.
pub const Pair = struct { p: G1, q: G2 };

/// Ate loop parameter 6x + 2 = 0x1800bffffffffffffffffffffe7ffc, split into
/// two 64-bit words.
///
/// The BN parameter x = 0x4001fffffffffffffffffffffbfff is positive, so the
/// loop needs no final conjugation.
const ate_loop_count_lo: u64 = 0xfffffffffffe7ffc;
const ate_loop_count_hi: u64 = 0x001800bfffffffff;

/// Compute the optimal Ate pairing e(P, Q).
/// Requires validated subgroup inputs.
pub fn pair(p: G1, q: G2) GT {
    if (p.isIdentity() or q.isIdentity()) {
        return GT.one;
    }

    const q_affine = q.affineCoordinates();
    const p_affine = p.affineCoordinates();

    var f = millerLoop(p_affine, q_affine);
    f = f.finalExponentiation();

    return f;
}

/// Compute the multi-pairing: e(P1, Q1) * e(P2, Q2) * ...
///
/// Cheaper than pairing each term separately and multiplying, because the
/// final exponentiation is only paid for once.
pub fn multiPair(pairs: []const Pair) GT {
    var f = GT.one;

    for (pairs) |item| {
        if (item.p.isIdentity() or item.q.isIdentity()) {
            continue;
        }
        const ml = millerLoop(item.p.affineCoordinates(), item.q.affineCoordinates());
        f = f.mul(ml);
    }

    return f.finalExponentiation();
}

/// Check that the pairing product is one. Callers must validate subgroup
/// membership and enforce the protocol's identity policy.
pub fn pairingCheck(pairs: []const Pair) bool {
    const result = multiPair(pairs);
    return result.isOne();
}

/// Accumulate the line functions along the Miller loop.
///
/// The loop parameter spans more than one word, so it is walked from the top
/// of the high word down through the low one.
pub fn millerLoop(p: G1Affine, q: G2Affine) Fp12 {
    var r = G2{ .x = q.x, .y = q.y, .z = Fp2.one };
    var f = Fp12.one;

    var found_one = false;
    var b: i8 = 63;
    while (b >= 0) : (b -= 1) {
        const bit_idx: u6 = @intCast(b);
        const bit_set = ((ate_loop_count_hi >> bit_idx) & 1) == 1;

        if (!found_one) {
            if (bit_set) {
                found_one = true;
            }
            continue;
        }

        f = f.sq();

        const dbl = doublingStep(&r);
        f = ell(f, dbl, p);

        if (bit_set) {
            const add_coeffs = additionStep(&r, q);
            f = ell(f, add_coeffs, p);
        }
    }

    b = 63;
    while (b >= 0) : (b -= 1) {
        const bit_idx: u6 = @intCast(b);
        const bit_set = ((ate_loop_count_lo >> bit_idx) & 1) == 1;

        f = f.sq();

        const dbl = doublingStep(&r);
        f = ell(f, dbl, p);

        if (bit_set) {
            const add_coeffs = additionStep(&r, q);
            f = ell(f, add_coeffs, p);
        }
    }

    // BN curves need two extra line evaluations beyond the loop, against the
    // Frobenius images of Q.
    const q1 = frobeniusEndomorphism(q);
    const add1 = additionStep(&r, q1);
    f = ell(f, add1, p);

    var q2 = frobeniusEndomorphismSquare(q);
    q2 = q2.neg();

    // Nothing reads `r` after this point, so the accumulator may be clobbered.
    const add2 = additionStep(&r, q2);
    f = ell(f, add2, p);

    return f;
}

/// Evaluate a line at P and fold the result into the accumulator.
///
/// The D-type twist places the non-zero coefficients at positions 0, 3 and 4,
/// which is what makes the sparse multiplication applicable.
fn ell(f: Fp12, coeffs: LineCoeffs, p: G1Affine) Fp12 {
    const c0 = coeffs.c0.mulByFp(p.y);
    const c3 = coeffs.c1.mulByFp(p.x);
    const c4 = coeffs.c2;
    return f.mulBy034(c0, c3, c4);
}

/// Line function coefficients.
///
/// These match the ones BLS12-381 uses; only the way they are embedded into
/// Fp12 differs, because the twist is of the other kind.
const LineCoeffs = struct {
    /// Scales the y coordinate of P.
    c0: Fp2,
    /// Scales the x coordinate of P.
    c1: Fp2,
    /// Constant term.
    c2: Fp2,
};

/// Double the running point and return the line tangent to it.
///
/// Algorithm 26 from https://eprint.iacr.org/2010/354.pdf.
pub fn doublingStep(r: *G2) LineCoeffs {
    const tmp0 = r.x.sq();
    const tmp1 = r.y.sq();
    const tmp2 = tmp1.sq();
    var tmp3 = tmp1.add(r.x).sq().sub(tmp0).sub(tmp2);
    tmp3 = tmp3.dbl();
    const tmp4 = tmp0.dbl().add(tmp0);
    const tmp6 = r.x.add(tmp4);
    const tmp5 = tmp4.sq();
    const zsquared = r.z.sq();

    const new_x = tmp5.sub(tmp3).sub(tmp3);
    const new_z = r.z.add(r.y).sq().sub(tmp1).sub(zsquared);

    const tmp2_8 = tmp2.dbl().dbl().dbl();
    const new_y = tmp3.sub(new_x).mul(tmp4).sub(tmp2_8);

    r.* = G2{ .x = new_x, .y = new_y, .z = new_z };

    var c0 = new_z.mul(zsquared);
    c0 = c0.dbl();

    var c1 = tmp4.mul(zsquared);
    c1 = c1.dbl();
    c1 = c1.neg();

    var c2 = tmp6.sq().sub(tmp0).sub(tmp5);
    const tmp1_4 = tmp1.dbl().dbl();
    c2 = c2.sub(tmp1_4);

    return .{ .c0 = c0, .c1 = c1, .c2 = c2 };
}

/// Add to the running point and return the line through the two points.
///
/// Algorithm 27 from https://eprint.iacr.org/2010/354.pdf.
pub fn additionStep(r: *G2, q: G2Affine) LineCoeffs {
    const zsquared = r.z.sq();
    const ysquared = q.y.sq();
    const t0 = zsquared.mul(q.x);
    const t1 = q.y.add(r.z).sq().sub(ysquared).sub(zsquared).mul(zsquared);
    const t2 = t0.sub(r.x);
    const t3 = t2.sq();
    const t4 = t3.dbl().dbl();
    const t5 = t4.mul(t2);
    const t6 = t1.sub(r.y).sub(r.y);
    const t9 = t6.mul(q.x);
    const t7 = t4.mul(r.x);

    const new_x = t6.sq().sub(t5).sub(t7).sub(t7);
    const new_z = r.z.add(t2).sq().sub(zsquared).sub(t3);

    const t10 = q.y.add(new_z);
    const t8 = t7.sub(new_x).mul(t6);
    const t0_new = r.y.mul(t5).dbl();
    const new_y = t8.sub(t0_new);

    r.* = G2{ .x = new_x, .y = new_y, .z = new_z };

    const t10_sq = t10.sq().sub(ysquared);
    const ztsquared = new_z.sq();
    const t10_final = t10_sq.sub(ztsquared);

    const c2 = t9.dbl().sub(t10_final);
    const c0 = new_z.dbl();

    const t6_neg = t6.neg();
    const c1 = t6_neg.dbl();

    return .{ .c0 = c0, .c1 = c1, .c2 = c2 };
}

/// Frobenius endomorphism π(Q) for G2 points.
///
/// Conjugation gives the p-th power in Fp2, and each coordinate then picks up
/// a twist factor.
fn frobeniusEndomorphism(q: G2Affine) G2Affine {
    const x = q.x.conjugate().mul(frobenius_coeff_x);
    const y = q.y.conjugate().mul(frobenius_coeff_y);
    return .{ .x = x, .y = y };
}

/// Frobenius endomorphism squared π²(Q), which collapses to a scaling.
fn frobeniusEndomorphismSquare(q: G2Affine) G2Affine {
    return .{ .x = q.x.mul(frobenius_coeff_x_sq), .y = q.y.mul(frobenius_coeff_y_sq) };
}

// Twist factors for π: γ₁ = ξ^((p-1)/3) on x and γ₂ = ξ^((p-1)/2) on y,
// where ξ = 2 + u.
const frobenius_coeff_x = Fp2{
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
};

const frobenius_coeff_y = Fp2{
    .c0 = Fp{ .limbs = .{
        0x68df801330cc654d,
        0x13115ae98aa9d16d,
        0x4daac6cb113fe602,
        0x5794c299b6b85112,
        0x0794dc2088b9689b,
        0xc175794b223d3f20,
        0xad93dfb6fe40843e,
        0x0000000000000c9c,
    } },
    .c1 = Fp{ .limbs = .{
        0xd1bf00266198ca9a,
        0x2622b5d31553a2da,
        0x9b558d96227fcc04,
        0xaf2985336d70a224,
        0x0f29b8411172d136,
        0x82eaf296447a7e40,
        0x5b27bf6dfc81087d,
        0x0000000000001939,
    } },
};

// The same factors for π², where γ₁² = γ₁ * γ₁^p and γ₂² = γ₂ * γ₂^p.
const frobenius_coeff_x_sq = Fp2{
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
    .c1 = Fp.zero,
};

const frobenius_coeff_y_sq = Fp2{
    .c0 = Fp{ .limbs = .{
        0xe5efc15421250f6b,
        0xc79db82c013a99e0,
        0xceacf1e1c9b47b1c,
        0xa9c1d8afbd3693bf,
        0xc299c4027a920164,
        0x004ff8788c9135d3,
        0x824e0f9fe4413ff6,
        0x000000000000141e,
    } },
    .c1 = Fp.zero,
};

test "pairing with identity" {
    const g1 = G1.basePoint;
    const g2 = G2.basePoint;

    // e(O, Q) = 1
    const result1 = pair(G1.identityElement, g2);
    try std.testing.expect(result1.isOne());

    // e(P, O) = 1
    const result2 = pair(g1, G2.identityElement);
    try std.testing.expect(result2.isOne());
}

test "pairing non-degeneracy" {
    const g1 = G1.basePoint;
    const g2 = G2.basePoint;

    // e(G1, G2) should not be 1.
    const result = pair(g1, g2);
    try std.testing.expect(!result.isOne());
}

test "first doubling step" {
    const g2 = G2.basePoint;
    const q_affine = g2.affineCoordinates();

    var r = G2{ .x = q_affine.x, .y = q_affine.y, .z = Fp2.one };
    const dbl = doublingStep(&r);

    // R should have been updated, so z is no longer 1.
    try std.testing.expect(!r.z.equivalent(Fp2.one));

    try std.testing.expect(!dbl.c0.isZero());
    try std.testing.expect(!dbl.c1.isZero());
}

test "frobenius endomorphisms preserve curve equation" {
    const q = G2.basePoint.affineCoordinates();

    const q1 = frobeniusEndomorphism(q);
    const q2 = frobeniusEndomorphismSquare(q);
    const q1_then = frobeniusEndomorphism(q1);

    // Ensure π(Q) and π²(Q) are still on the twisted curve.
    try std.testing.expect(q1.y.sq().equivalent(q1.x.sq().mul(q1.x).add(G2.B)));
    try std.testing.expect(q2.y.sq().equivalent(q2.x.sq().mul(q2.x).add(G2.B)));

    // Ensure π²(Q) matches applying π twice.
    try std.testing.expect(q2.x.equivalent(q1_then.x));
    try std.testing.expect(q2.y.equivalent(q1_then.y));
}
