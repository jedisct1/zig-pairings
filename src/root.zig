//! Pairing-friendly elliptic curve library for Zig.
//!
//! This library implements the BLS12-381 curve, intended for cryptographic
//! applications including BLS signatures and zero-knowledge proofs.
//!
//! Features:
//! - Field arithmetic: Fp, Fp2, Fp6, Fp12 with Montgomery multiplication
//! - Curve groups: G1 (over Fp) and G2 (over Fp2 via sextic twist)
//! - Optimal Ate pairing with Miller loop and final exponentiation
//! - Serialization following standard compressed/uncompressed formats
//! - Scalar arithmetic in the subgroup order field
//!
//! Example:
//! ```zig
//! const bls = @import("pairing").bls12_381;
//!
//! // Compute pairing of generators
//! const g1 = bls.G1.basePoint;
//! const g2 = bls.G2.basePoint;
//! const gt = bls.pairing.pair(g1, g2);
//!
//! // Identity handling
//! const e_identity = bls.pairing.pair(bls.G1.identityElement, g2);
//! std.debug.assert(e_identity.isOne());
//! ```

pub const bls12_381 = @import("bls12_381.zig");

pub const Fp = bls12_381.Fp;
pub const Fp2 = bls12_381.Fp2;
pub const Fp6 = bls12_381.Fp6;
pub const Fp12 = bls12_381.Fp12;
pub const G1 = bls12_381.G1;
pub const G2 = bls12_381.G2;
pub const scalar = bls12_381.scalar;
pub const pairing = bls12_381.pairing;

test {
    _ = bls12_381;
}
