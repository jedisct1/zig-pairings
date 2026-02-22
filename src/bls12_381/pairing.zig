//! Optimal Ate pairing for BLS12-381.
//!
//! Computes the bilinear pairing e: G1 × G2 → GT.

const std = @import("std");
const Fp2 = @import("fp2.zig").Fp2;
const Fp12 = @import("fp12.zig").Fp12;
const G1 = @import("g1.zig").G1;
const G1Affine = @import("g1.zig").AffineCoordinates;
const G2 = @import("g2.zig").G2;
const G2Affine = @import("g2.zig").AffineCoordinates;

/// Element of the target group GT (a subgroup of Fp12*).
pub const GT = Fp12;

/// Compute the optimal Ate pairing e(P, Q).
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
pub fn multiPair(pairs: []const struct { p: G1, q: G2 }) GT {
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

/// Verify that e(P1, Q1) * e(P2, Q2) = 1.
///
/// This is the shape a BLS signature verification takes.
pub fn pairingCheck(pairs: []const struct { p: G1, q: G2 }) bool {
    const result = multiPair(pairs);
    return result.isOne();
}

/// Accumulate the line functions along the Miller loop.
///
/// The loop parameter is the BLS parameter t = -0xd201000000010000.
///
/// Because it is negative, the loop starts from -Q and subtracts at every set
/// bit rather than adding.
pub fn millerLoop(p: G1Affine, q: G2Affine) Fp12 {
    var r = G2{ .x = q.x, .y = q.y.neg(), .z = Fp2.one };
    const q_neg = G2Affine{ .x = q.x, .y = q.y.neg() };

    const t_abs: u64 = 0xd201000000010000;

    var f = Fp12.one;

    // The top bit is already accounted for by starting at -Q.
    var b: i8 = 62;
    while (b >= 0) : (b -= 1) {
        const bit_idx: u6 = @intCast(b);
        const bit_set = ((t_abs >> bit_idx) & 1) == 1;

        f = f.sq();

        const dbl = doublingStep(&r);
        f = ell(f, dbl, p);

        if (bit_set) {
            const add_coeffs = additionStep(&r, q_neg);
            f = ell(f, add_coeffs, p);
        }
    }

    return f;
}

/// Evaluate a line at P and fold the result into the accumulator.
fn ell(f: Fp12, coeffs: LineCoeffs, p: G1Affine) Fp12 {
    const c0_scaled = coeffs.c0.mulByFp(p.y);
    const c1_scaled = coeffs.c1.mulByFp(p.x);
    return f.mulBy014(coeffs.c2, c1_scaled, c0_scaled);
}

/// Line function coefficients.
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
fn doublingStep(r: *G2) LineCoeffs {
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
fn additionStep(r: *G2, q: G2Affine) LineCoeffs {
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

test "first ell application" {
    const g1 = G1.basePoint;
    const g2 = G2.basePoint;
    const p_affine = g1.affineCoordinates();
    const q_affine = g2.affineCoordinates();

    var r = G2{ .x = q_affine.x, .y = q_affine.y, .z = Fp2.one };

    const dbl = doublingStep(&r);

    var f = Fp12.one;
    f = ell(f, dbl, p_affine);

    try std.testing.expect(!f.isOne());
}

test "pairing bilinearity" {
    const g1 = G1.basePoint;
    const g2 = G2.basePoint;

    const result = pair(g1, g2);
    try std.testing.expect(!result.isOne());
}

test "miller loop output" {
    const g1 = G1.basePoint;
    const g2 = G2.basePoint;
    const q_affine = g2.affineCoordinates();
    const p_affine = g1.affineCoordinates();

    const f = millerLoop(p_affine, q_affine);

    // Known-answer value for the Miller loop, before final exponentiation.
    const expected_c0_c0_hex = "05194f5785436c8debf0eb2bab4c6ef3de7dc0633c85769173777b782bf897fa45025fd03e7be941123c4ee19910e62e";
    const actual_c0_c0_hex = std.fmt.bytesToHex(f.c0.c0.c0.toBytes(.big), .lower);
    try std.testing.expectEqualStrings(expected_c0_c0_hex, &actual_c0_c0_hex);

    const final_result = f.finalExponentiation();
    const expected_final_prefix = "11619b45f61edfe3";
    const final_hex = std.fmt.bytesToHex(final_result.c0.c0.c0.toBytes(.big), .lower);
    try std.testing.expectEqualStrings(expected_final_prefix, final_hex[0..16]);
}

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

test "pairing bilinearity - scalar" {
    // e(2P, Q) = e(P, Q)²
    const g1 = G1.basePoint;
    const g2 = G2.basePoint;
    const two_g1 = g1.dbl();

    const e_2p_q = pair(two_g1, g2);
    const e_p_q = pair(g1, g2);
    const e_p_q_sq = e_p_q.mul(e_p_q);

    try std.testing.expect(e_2p_q.equivalent(e_p_q_sq));
}

test "pairing negation" {
    // e(-P, Q) * e(P, Q) = 1
    const g1 = G1.basePoint;
    const g2 = G2.basePoint;
    const neg_g1 = g1.neg();

    const e1 = pair(neg_g1, g2);
    const e2 = pair(g1, g2);
    const product = e1.mul(e2);

    try std.testing.expect(product.isOne());
}
