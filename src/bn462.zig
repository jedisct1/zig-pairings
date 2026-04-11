//! BN462 pairing-friendly elliptic curve implementation.
//!
//! BN462 is a Barreto-Naehrig curve with embedding degree 12, designed for
//! approximately 134-bit security.
//!
//! Field tower construction:
//! - Fp: Base field (462-bit prime)
//! - Fp2 = Fp[u] / (u² + 1)
//!
//! BN parameter: x = 2^114 + 2^101 - 2^14 - 1

pub const Fp = @import("bn462/fp.zig").Fp;
pub const Fp2 = @import("bn462/fp2.zig").Fp2;

test {
    _ = Fp;
    _ = Fp2;
}
