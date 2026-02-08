//! Pairing-friendly elliptic curve library for Zig.
//!
//! This library implements the BLS12-381 curve, intended for cryptographic
//! applications including BLS signatures and zero-knowledge proofs.

pub const bls12_381 = @import("bls12_381.zig");

pub const Fp = bls12_381.Fp;
pub const Fp2 = bls12_381.Fp2;
pub const Fp6 = bls12_381.Fp6;
pub const Fp12 = bls12_381.Fp12;
pub const G1 = bls12_381.G1;
pub const G2 = bls12_381.G2;
pub const scalar = bls12_381.scalar;

test {
    _ = bls12_381;
}
