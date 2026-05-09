//! BN462 pairing-friendly elliptic curve implementation.
//!
//! BN462 is a Barreto-Naehrig curve with embedding degree 12, designed for
//! approximately 134-bit security.
//!
//! Field tower construction:
//! - Fp: Base field (462-bit prime)
//! - Fp2 = Fp[u] / (u² + 1)
//! - Fp6 = Fp2[v] / (v³ - (2 + u))
//! - Fp12 = Fp6[w] / (w² - v)
//!
//! BN parameter: x = 2^114 + 2^101 - 2^14 - 1

pub const Fp = @import("bn462/fp.zig").Fp;
pub const Fp2 = @import("bn462/fp2.zig").Fp2;
pub const Fp6 = @import("bn462/fp6.zig").Fp6;
pub const Fp12 = @import("bn462/fp12.zig").Fp12;
pub const scalar = @import("bn462/scalar.zig");

test {
    _ = Fp;
    _ = Fp2;
    _ = Fp6;
    _ = Fp12;
    _ = scalar;
}
