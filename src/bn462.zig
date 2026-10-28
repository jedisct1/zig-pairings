//! BN462, as defined in Section 4.2.2 of
//! draft-irtf-cfrg-pairing-friendly-curves-14.
//!
//! Its security level is about 134 bits.
//!
//! It provides:
//! - G1: the curve E(Fp): y^2 = x^3 + 5, which has prime order r
//! - G2: the subgroup of order r of the twist E'(Fp2): y^2 = x^3 + (2 - u)
//! - GT: the target group, in Fp12
//! - The optimal Ate pairing e: G1 x G2 -> GT
//!
//! Field tower:
//! - Fp2 = Fp[u] / (u^2 + 1)
//! - Fp6 = Fp2[v] / (v^3 - (u + 2))
//! - Fp12 = Fp6[w] / (w^2 - v)
//!
//! BN parameter: t = 2^114 + 2^101 - 2^14 - 1
//!
//! Scalars follow Section 5.2 of the draft, and the pairing matches the test
//! vector of its Appendix C.
//!
//! The draft doesn't define a point encoding for BN462.
//! The one used here is described in the G1 module.
//!
//! References:
//! - https://datatracker.ietf.org/doc/draft-irtf-cfrg-pairing-friendly-curves/14/
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
