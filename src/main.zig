const std = @import("std");
const pairing_lib = @import("pairing");
const bls12_381 = pairing_lib.bls12_381;
const pair = bls12_381.pairing.pair;
const G1 = bls12_381.G1;
const G2 = bls12_381.G2;

pub fn main() !void {
    const g1 = G1.basePoint;
    const g2 = G2.basePoint;

    std.debug.print("Computing pairing e(G1, G2)...\n", .{});
    const e_p_q = pair(g1, g2);

    const two_g1 = g1.dbl();
    const e_2p_q = pair(two_g1, g2);
    const e_p_q_sq = e_p_q.mul(e_p_q);

    std.debug.print("e(2P, Q) == e(P, Q)^2: {}\n", .{e_2p_q.equivalent(e_p_q_sq)});

    const two_g2 = g2.dbl();
    const e_p_2q = pair(g1, two_g2);
    std.debug.print("e(P, 2Q) == e(P, Q)^2: {}\n", .{e_p_2q.equivalent(e_p_q_sq)});

    const neg_g1 = g1.neg();
    const e_neg_p_q = pair(neg_g1, g2);
    const product = e_neg_p_q.mul(e_p_q);
    std.debug.print("e(-P, Q) * e(P, Q) == 1: {}\n", .{product.isOne()});
}

test "simple pairing test" {
    const g1 = G1.basePoint;
    const g2 = G2.basePoint;

    const e_p_q = pair(g1, g2);
    try std.testing.expect(!e_p_q.isOne());
}
