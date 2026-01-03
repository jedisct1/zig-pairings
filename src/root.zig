//! Pairing-friendly elliptic curve library for Zig.
//!
//! This library implements the BLS12-381 curve, intended for cryptographic
//! applications including BLS signatures and zero-knowledge proofs.

pub const bls12_381 = @import("bls12_381.zig");

pub const Fp = bls12_381.Fp;
pub const Fp2 = bls12_381.Fp2;

test {
    _ = bls12_381;
}
