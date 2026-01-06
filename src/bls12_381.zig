//! BLS12-381 pairing-friendly elliptic curve implementation.
//!
//! BLS12-381 is a pairing-friendly curve with embedding degree 12, designed for
//! approximately 128-bit security.
//!
//! Field tower construction:
//! - Fp: Base field (381-bit prime)
//! - Fp2 = Fp[u] / (u² + 1)
//! - Fp6 = Fp2[v] / (v³ - (1+u))

pub const Fp = @import("bls12_381/fp.zig").Fp;
pub const Fp2 = @import("bls12_381/fp2.zig").Fp2;
pub const Fp6 = @import("bls12_381/fp6.zig").Fp6;

test {
    _ = Fp;
    _ = Fp2;
    _ = Fp6;
}
