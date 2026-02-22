//! BLS12-381 pairing-friendly elliptic curve implementation.
//!
//! BLS12-381 is a pairing-friendly curve with embedding degree 12, designed for
//! approximately 128-bit security.
//!
//! It provides:
//! - G1: A 255-bit prime order subgroup of E(Fp) where y² = x³ + 4
//! - G2: A 255-bit prime order subgroup of E'(Fp²) where y² = x³ + 4(1+u)
//! - GT: The target group, a subgroup of Fp12*
//! - Optimal Ate pairing: e: G1 × G2 → GT
//!
//! Field tower construction:
//! - Fp: Base field (381-bit prime)
//! - Fp2 = Fp[u] / (u² + 1)
//! - Fp6 = Fp2[v] / (v³ - (1+u))
//! - Fp12 = Fp6[w] / (w² - v)
//!
//! References:
//! - https://hackmd.io/@benjaminion/bls12-381
//! - https://eips.ethereum.org/EIPS/eip-2537
//! - https://datatracker.ietf.org/doc/draft-irtf-cfrg-pairing-friendly-curves/

pub const Fp = @import("bls12_381/fp.zig").Fp;
pub const Fp2 = @import("bls12_381/fp2.zig").Fp2;
pub const Fp6 = @import("bls12_381/fp6.zig").Fp6;
pub const Fp12 = @import("bls12_381/fp12.zig").Fp12;
pub const G1 = @import("bls12_381/g1.zig").G1;
pub const G2 = @import("bls12_381/g2.zig").G2;
pub const scalar = @import("bls12_381/scalar.zig");
pub const pairing = @import("bls12_381/pairing.zig");

test {
    _ = Fp;
    _ = Fp2;
    _ = Fp6;
    _ = Fp12;
    _ = G1;
    _ = G2;
    _ = scalar;
    _ = pairing;
}
