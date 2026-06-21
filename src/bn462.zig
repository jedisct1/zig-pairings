//! BN462 pairing-friendly elliptic curve implementation.
//!
//! BN462 is a Barreto-Naehrig curve with embedding degree 12, designed for
//! approximately 134-bit security.
//!
//! It provides:
//! - G1: A prime order subgroup of E(Fp) where y² = x³ + 5
//! - G2: A prime order subgroup of E'(Fp²) where y² = x³ + (2 - u)
//! - GT: The target group, a subgroup of Fp12*
//! - Optimal Ate pairing: e: G1 × G2 → GT
//!
//! Field tower construction:
//! - Fp: Base field (462-bit prime)
//! - Fp2 = Fp[u] / (u² + 1)
//! - Fp6 = Fp2[v] / (v³ - (2 + u))
//! - Fp12 = Fp6[w] / (w² - v)
//!
//! BN parameter: x = 2^114 + 2^101 - 2^14 - 1
//!
//! The IETF specification deliberately defines no point serialization format
//! for BN462: a 462-bit characteristic needs 58 bytes, which leaves only two
//! spare bits in the leading byte where the metadata scheme needs three.
//!
//! The compressed and uncompressed encodings here are therefore specific to
//! this library.
//!
//! Field element and scalar encodings do follow the specification.
//!
//! References:
//! - https://datatracker.ietf.org/doc/draft-irtf-cfrg-pairing-friendly-curves/
//! - Barreto, P.S.L.M. and Naehrig, M. "Pairing-Friendly Elliptic Curves of Prime Order"

pub const Fp = @import("bn462/fp.zig").Fp;
pub const Fp2 = @import("bn462/fp2.zig").Fp2;
pub const Fp6 = @import("bn462/fp6.zig").Fp6;
pub const Fp12 = @import("bn462/fp12.zig").Fp12;
pub const G1 = @import("bn462/g1.zig").G1;
pub const G2 = @import("bn462/g2.zig").G2;
pub const scalar = @import("bn462/scalar.zig");
pub const pairing = @import("bn462/pairing.zig");

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
