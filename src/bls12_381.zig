//! BLS12-381 pairing-friendly elliptic curve implementation.
//!
//! BLS12-381 is a pairing-friendly curve with embedding degree 12, designed for
//! approximately 128-bit security.
//!
//! Field tower construction:
//! - Fp: Base field (381-bit prime)

pub const Fp = @import("bls12_381/fp.zig").Fp;

test {
    _ = Fp;
}
