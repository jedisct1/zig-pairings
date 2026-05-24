//! Pairing-friendly elliptic curve library for Zig.
//!
//! This library implements the BLS12-381 and BN462 curves, intended for
//! cryptographic applications including BLS signatures and zero-knowledge
//! proofs.
//!
//! The top-level aliases below default to BLS12-381.
//!
//! BN462 lives under the `bn462` namespace.
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
//!
//! BLS12-381 point encodings follow the format given by the IETF
//! specification, so they interoperate with any library that implements it.
//!
//! That format does not carry over to BN462, whose larger characteristic
//! leaves too few spare bits for the metadata, and the specification
//! deliberately leaves BN462 point encoding undefined.
//!
//! The BN462 encoding used here is therefore specific to this library.

pub const bls12_381 = @import("bls12_381.zig");
pub const bn462 = @import("bn462.zig");

// Top-level aliases default to the BLS12-381 curve.
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
    _ = bn462;
    _ = @import("interop_test.zig");
}
