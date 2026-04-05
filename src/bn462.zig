//! BN462 pairing-friendly elliptic curve implementation.
//!
//! BN462 is a Barreto-Naehrig curve with embedding degree 12, designed for
//! approximately 134-bit security.
//!
//! Field tower construction:
//! - Fp: Base field (462-bit prime)
//!
//! BN parameter: x = 2^114 + 2^101 - 2^14 - 1

pub const Fp = @import("bn462/fp.zig").Fp;

test {
    _ = Fp;
}
