//! BLS12-381, as defined in Section 4.2.1 of
//! draft-irtf-cfrg-pairing-friendly-curves-14.
//!
//! Its security level is expected to be about 126 bits rather than 128.
//!
//! It provides:
//! - G1: the subgroup of order r of E(Fp): y^2 = x^3 + 4
//! - G2: the subgroup of order r of the twist E'(Fp2): y^2 = x^3 + 4(u + 1)
//! - GT: the target group, in Fp12
//! - The optimal Ate pairing e: G1 x G2 -> GT
//!
//! Field tower:
//! - Fp2 = Fp[u] / (u^2 + 1)
//! - Fp6 = Fp2[v] / (v^3 - (u + 1))
//! - Fp12 = Fp6[w] / (w^2 - v)
//!
//! Point and scalar encodings follow Section 5 of the draft, and the test
//! vectors of its Appendices C and D pass.
//!
//! References:
//! - https://datatracker.ietf.org/doc/draft-irtf-cfrg-pairing-friendly-curves/14/
//! - https://www.rfc-editor.org/rfc/rfc9380
//! - https://hackmd.io/@benjaminion/bls12-381

pub const Fp = @import("bls12_381/fp.zig").Fp;
pub const Fp2 = @import("bls12_381/fp2.zig").Fp2;
pub const Fp6 = @import("bls12_381/fp6.zig").Fp6;
pub const Fp12 = @import("bls12_381/fp12.zig").Fp12;
pub const G1 = @import("bls12_381/g1.zig").G1;
pub const G2 = @import("bls12_381/g2.zig").G2;
pub const scalar = @import("bls12_381/scalar.zig");
pub const pairing = @import("bls12_381/pairing.zig");

test {
    _ = Fp;
    _ = Fp2;
    _ = Fp6;
    _ = Fp12;
    _ = G1;
    _ = G2;
    _ = scalar;
    _ = pairing;
    _ = @import("bls12_381/ietf_test.zig");
}
