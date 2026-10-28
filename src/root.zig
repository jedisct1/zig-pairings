//! Pairing-friendly elliptic curves for Zig.
//!
//! This library implements BLS12-381 and BN462, the two curves recommended at
//! the 128-bit security level by draft-irtf-cfrg-pairing-friendly-curves-14.
//! They are used for things like BLS signatures and zero-knowledge proofs.
//!
//! The top-level aliases below point to BLS12-381.
//! BN462 lives under the `bn462` namespace and has the same layout.
//!
//! Example:
//! ```zig
//! const bls = @import("pairings").bls12_381;
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
//! Points decoded from bytes or coordinates are always checked to be in their
//! group, and group operations keep them there.
//! A point built by hand from its fields has to be checked with
//! `isInSubgroup` before it is used.
//!
//! The draft leaves a few choices to the protocol using it:
//! - `fromCompressed` and `fromUncompressed` accept one form each, and
//!   `fromBytes` accepts both.
//! - Decoding accepts the identity element, and `rejectIdentity` rejects it.
//! - `scalar.fromBytes` accepts zero, and `scalar.isZero` detects it.
//!
//! Use `mul` for secret scalars, as it runs in constant time.
//! `mulPublic` is faster, but only for scalars that are not secret.

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
    _ = @import("bn462_interop_test.zig");
}
