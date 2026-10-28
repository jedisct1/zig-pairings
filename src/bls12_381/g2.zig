//! G2 group for BLS12-381.
//!
//! G2 is the subgroup of order r of the twisted curve E'(Fp2): y^2 = x^3 + 4(1 + u).
//! The specification calls it G'_2.
//!
//! Points decoded from bytes or coordinates are always checked to be in G2.

const std = @import("std");
const crypto = std.crypto;
const Fp = @import("fp.zig").Fp;
const Fp2 = @import("fp2.zig").Fp2;
const scalar = @import("scalar.zig");

const EncodingError = crypto.errors.EncodingError;
const IdentityElementError = crypto.errors.IdentityElementError;
const NonCanonicalError = crypto.errors.NonCanonicalError;
const NotSquareError = crypto.errors.NotSquareError;

/// G2 point in projective coordinates (X, Y, Z) where x = X/Z, y = Y/Z.
pub const G2 = struct {
    x: Fp2,
    y: Fp2,
    z: Fp2,

    /// Curve coefficient b' = 4(1 + u) for the twisted curve.
    pub const B = Fp2{
        .c0 = Fp.fromInt(4),
        .c1 = Fp.fromInt(4),
    };

    /// The generator point.
    pub const basePoint = G2{
        .x = Fp2{
            // x0 = 0x024aa2b2f08f0a91260805272dc51051c6e47ad4fa403b02b4510b647ae3d1770bac0326a805bbefd48056c8c121bdb8
            .c0 = Fp{
                .limbs = .{
                    0xf5f28fa202940a10,
                    0xb3f5fb2687b4961a,
                    0xa1a893b53e2ae580,
                    0x9894999d1a3caee9,
                    0x6f67b7631863366b,
                    0x058191924350bcd7,
                },
            },
            // x1 = 0x13e02b6052719f607dacd3a088274f65596bd0d09920b61ab5da61bbdc7f5049334cf11213945d57e5ac7d055d042b7e
            .c1 = Fp{
                .limbs = .{
                    0xa5a9c0759e23f606,
                    0xaaa0c59dbccd60c3,
                    0x3bb17e18e2867806,
                    0x1b1ab6cc8541b367,
                    0xc2b6ed0ef2158547,
                    0x11922a097360edf3,
                },
            },
        },
        .y = Fp2{
            // y0 = 0x0ce5d527727d6e118cc9cdc6da2e351aadfd9baa8cbdd3a76d429a695160d12c923ac9cc3baca289e193548608b82801
            .c0 = Fp{
                .limbs = .{
                    0x4c730af860494c4a,
                    0x597cfa1f5e369c5a,
                    0xe7e6856caa0a635a,
                    0xbbefb5e96e0d495f,
                    0x07d3a975f0ef25a2,
                    0x0083fd8e7e80dae5,
                },
            },
            // y1 = 0x0606c4a02ea734cc32acd2b02bc28b99cb3e287e85a763af267492ab572e99ab3f370d275cec1da1aaa9075ff05f79be
            .c1 = Fp{
                .limbs = .{
                    0xadc0fc92df64b05d,
                    0x18aa270a2b1461dc,
                    0x86adac6a3be4eba0,
                    0x79495c4ec93da33a,
                    0xe7175850a43ccaed,
                    0x0b2bc2a163de1bf2,
                },
            },
        },
        .z = Fp2.one,
    };

    /// The identity element (point at infinity).
    pub const identityElement = G2{ .x = Fp2.zero, .y = Fp2.one, .z = Fp2.zero };

    /// Compressed point encoding size.
    pub const compressed_length = 96;

    /// Uncompressed point encoding size.
    pub const uncompressed_length = 192;

    /// Create a point from affine coordinates, checking that it is in G2.
    ///
    /// `AffineCoordinates.identityElement` gives the identity element.
    pub fn fromAffineCoordinates(affine: AffineCoordinates) EncodingError!G2 {
        if (affine.x.isZero() and affine.y.isZero()) {
            return identityElement;
        }
        return fromCoordinates(affine.x, affine.y);
    }

    /// Build a point from its coordinates, checking that it is in G2.
    fn fromCoordinates(x: Fp2, y: Fp2) EncodingError!G2 {
        const p = G2{ .x = x, .y = y, .z = Fp2.one };
        if (!p.isInSubgroup()) {
            return error.InvalidEncoding;
        }
        return p;
    }

    /// Deserialize a point in either form, like octets_to_point_G2 in the specification.
    ///
    /// The compression bit must match the length, so that an uncompressed G1
    /// point, which is also 96 bytes long, is not read as a G2 point.
    ///
    /// Protocols that hash or compare encoded points should use
    /// `fromCompressed` or `fromUncompressed` instead.
    pub fn fromBytes(bytes: []const u8) (EncodingError || NotSquareError || NonCanonicalError)!G2 {
        if (bytes.len == 0) {
            return error.InvalidEncoding;
        }
        if ((bytes[0] & 0x80) != 0) {
            if (bytes.len != compressed_length) return error.InvalidEncoding;
            return fromCompressed(bytes[0..compressed_length].*);
        }
        if (bytes.len != uncompressed_length) return error.InvalidEncoding;
        return fromUncompressed(bytes[0..uncompressed_length].*);
    }

    /// Deserialize from compressed form.
    ///
    /// The encoding is the x-coordinate with flags in the top three bits,
    /// laid out the same way as for G1, with the sign given by sign_GF_p^2.
    ///
    /// The point is checked to be in G2.
    pub fn fromCompressed(bytes: [compressed_length]u8) (EncodingError || NotSquareError || NonCanonicalError)!G2 {
        const flags = bytes[0] >> 5;
        const is_compressed = (flags & 0b100) != 0;
        const is_infinity = (flags & 0b010) != 0;
        const y_sign = (flags & 0b001) != 0;

        // The sign bit only carries meaning on a compressed point that is not
        // the point at infinity.
        //
        // Any other encoding must leave it clear.
        if (y_sign and (!is_compressed or is_infinity)) {
            return error.InvalidEncoding;
        }

        if (!is_compressed) {
            return error.InvalidEncoding;
        }

        if (is_infinity) {
            var all_zero = true;
            for (bytes[1..]) |b| {
                if (b != 0) all_zero = false;
            }
            if (!all_zero or (bytes[0] & 0x1f) != 0) {
                return error.InvalidEncoding;
            }
            return identityElement;
        }

        var x_bytes = bytes;
        x_bytes[0] &= 0x1f;
        const x = try Fp2.fromBytes(x_bytes, .big);

        const x3_b = x.sq().mul(x).add(B);
        var y = try x3_b.sqrt();

        const y_is_largest = y.lexicographicallyLargest();
        if (y_sign != y_is_largest) {
            y = y.neg();
        }

        return fromCoordinates(x, y);
    }

    /// Serialize to compressed form.
    pub fn toCompressed(p: G2) [compressed_length]u8 {
        const affine = p.affineCoordinates();
        const is_identity = p.isIdentity();

        var bytes: [compressed_length]u8 = undefined;

        if (is_identity) {
            @memset(&bytes, 0);
            bytes[0] = 0b11000000;
        } else {
            bytes = affine.x.toBytes(.big);
            const y_sign: u8 = if (affine.y.lexicographicallyLargest()) 1 else 0;
            bytes[0] |= 0b10000000 | (y_sign << 5);
        }

        return bytes;
    }

    /// Deserialize from uncompressed form.
    ///
    /// The point is checked to be in G2.
    pub fn fromUncompressed(bytes: [uncompressed_length]u8) (EncodingError || NonCanonicalError)!G2 {
        const flags = bytes[0] >> 5;
        const is_compressed = (flags & 0b100) != 0;
        const is_infinity = (flags & 0b010) != 0;
        const y_sign = (flags & 0b001) != 0;

        // The sign bit only carries meaning on a compressed point that is not
        // the point at infinity.
        //
        // Any other encoding must leave it clear.
        if (y_sign and (!is_compressed or is_infinity)) {
            return error.InvalidEncoding;
        }

        if (is_compressed) {
            return error.InvalidEncoding;
        }

        if (is_infinity) {
            var all_zero = true;
            for (bytes[1..]) |b| {
                if (b != 0) all_zero = false;
            }
            if (!all_zero or (bytes[0] & 0x1f) != 0) {
                return error.InvalidEncoding;
            }
            return identityElement;
        }

        var x_bytes = bytes[0..96].*;
        x_bytes[0] &= 0x1f;
        const x = try Fp2.fromBytes(x_bytes, .big);
        const y = try Fp2.fromBytes(bytes[96..192].*, .big);

        return fromCoordinates(x, y);
    }

    /// Serialize to uncompressed form.
    pub fn toUncompressed(p: G2) [uncompressed_length]u8 {
        const affine = p.affineCoordinates();
        const is_identity = p.isIdentity();

        var bytes: [uncompressed_length]u8 = undefined;

        if (is_identity) {
            @memset(&bytes, 0);
            bytes[0] = 0b01000000;
        } else {
            bytes[0..96].* = affine.x.toBytes(.big);
            bytes[96..192].* = affine.y.toBytes(.big);
        }

        return bytes;
    }

    /// Return true if this is the identity element.
    pub fn isIdentity(p: G2) bool {
        return p.z.isZero();
    }

    /// Reject the identity element.
    ///
    /// Decoding accepts the identity, so protocols that don't want it should
    /// call this on the result.
    pub fn rejectIdentity(p: G2) IdentityElementError!void {
        if (p.isIdentity()) {
            return error.IdentityElement;
        }
    }

    /// Convert to affine coordinates.
    pub fn affineCoordinates(p: G2) AffineCoordinates {
        if (p.isIdentity()) {
            return AffineCoordinates.identityElement;
        }

        const z_inv = p.z.invert();
        return .{
            .x = p.x.mul(z_inv),
            .y = p.y.mul(z_inv),
        };
    }

    /// Return true if both points are equivalent.
    pub fn equivalent(a: G2, b: G2) bool {
        const x1z2 = a.x.mul(b.z);
        const x2z1 = b.x.mul(a.z);
        const y1z2 = a.y.mul(b.z);
        const y2z1 = b.y.mul(a.z);

        return x1z2.equivalent(x2z1) and y1z2.equivalent(y2z1);
    }

    /// Conditionally replace the point with `q` when `choice` is 1 (constant-time).
    pub fn cMov(p: *G2, q: G2, choice: u1) void {
        p.x.cMov(q.x, choice);
        p.y.cMov(q.y, choice);
        p.z.cMov(q.z, choice);
    }

    /// Swap two points when `choice` is 1 (constant-time).
    fn cSwap(p: *G2, q: *G2, choice: u1) void {
        const t = p.*;
        p.cMov(q.*, choice);
        q.cMov(t, choice);
    }

    /// Negate the point.
    pub fn neg(p: G2) G2 {
        return .{
            .x = p.x,
            .y = p.y.neg(),
            .z = p.z,
        };
    }

    /// Multiply by 3b' = 12(1 + u).
    fn mulByB3(x: Fp2) Fp2 {
        const sum = x.c0.add(x.c1);
        const diff = x.c0.sub(x.c1);
        return .{
            .c0 = diff.dbl().add(diff).dbl().dbl(),
            .c1 = sum.dbl().add(sum).dbl().dbl(),
        };
    }

    /// Double a point.
    ///
    /// Algorithm 9 from https://eprint.iacr.org/2015/1060.
    pub fn dbl(p: G2) G2 {
        var t0 = p.y.sq();
        var z3 = t0.dbl().dbl().dbl();
        var t1 = p.y.mul(p.z);
        var t2 = mulByB3(p.z.sq());
        var x3 = t2.mul(z3);
        var y3 = t0.add(t2);
        z3 = t1.mul(z3);
        t1 = t2.dbl();
        t2 = t1.add(t2);
        t0 = t0.sub(t2);
        y3 = t0.mul(y3).add(x3);
        t1 = p.x.mul(p.y);
        x3 = t0.mul(t1).dbl();

        return .{ .x = x3, .y = y3, .z = z3 };
    }

    /// Add two points.
    ///
    /// Algorithm 7 from https://eprint.iacr.org/2015/1060, which works for any
    /// two points, including equal points and the identity.
    pub fn add(p: G2, q: G2) G2 {
        var t0 = p.x.mul(q.x);
        var t1 = p.y.mul(q.y);
        var t2 = p.z.mul(q.z);
        const t3 = p.x.add(p.y).mul(q.x.add(q.y)).sub(t0.add(t1));
        const t4 = p.y.add(p.z).mul(q.y.add(q.z)).sub(t1.add(t2));
        var y3 = p.x.add(p.z).mul(q.x.add(q.z)).sub(t0.add(t2));
        t0 = t0.dbl().add(t0);
        t2 = mulByB3(t2);
        var z3 = t1.add(t2);
        t1 = t1.sub(t2);
        y3 = mulByB3(y3);
        const x3 = t3.mul(t1).sub(t4.mul(y3));
        y3 = t1.mul(z3).add(y3.mul(t0));
        z3 = z3.mul(t4).add(t0.mul(t3));

        return .{ .x = x3, .y = y3, .z = z3 };
    }

    /// Subtract two points.
    pub fn sub(p: G2, q: G2) G2 {
        return p.add(q.neg());
    }

    /// Multiply a point by a scalar, in constant time.
    ///
    /// The scalar doesn't have to be reduced, and zero gives the identity.
    pub fn mul(p: G2, scalar_bytes: [scalar.encoded_length]u8, endian: std.builtin.Endian) G2 {
        var s = scalar_bytes;
        if (endian == .little) std.mem.reverse(u8, &s);

        var r0 = identityElement;
        var r1 = p;
        for (s) |byte| {
            inline for (0..8) |k| {
                const bit: u1 = @truncate(byte >> (7 - k));
                cSwap(&r0, &r1, bit);
                r1 = r0.add(r1);
                r0 = r0.dbl();
                cSwap(&r0, &r1, bit);
            }
        }
        return r0;
    }

    /// Multiply a point by a public scalar, in variable time.
    pub fn mulPublic(p: G2, scalar_bytes: [scalar.encoded_length]u8, endian: std.builtin.Endian) G2 {
        var s = scalar_bytes;
        if (endian == .little) std.mem.reverse(u8, &s);
        const start = std.mem.findNone(u8, &s, &.{0}) orelse return identityElement;

        var q = identityElement;
        for (s[start..]) |byte| {
            inline for (0..8) |k| {
                q = q.dbl();
                if ((byte >> (7 - k)) & 1 == 1) {
                    q = q.add(p);
                }
            }
        }
        return q;
    }

    /// The untwist-Frobenius-twist endomorphism, used for cofactor clearing.
    pub fn psi(p: G2) G2 {
        return .{
            .x = p.x.conjugate().mul(psi_coeff_x),
            .y = p.y.conjugate().mul(psi_coeff_y),
            .z = p.z.conjugate(),
        };
    }

    /// Psi applied twice, which collapses to a scaling and a negation.
    pub fn psi2(p: G2) G2 {
        return .{
            .x = p.x.mul(psi2_coeff_x),
            .y = p.y.neg(),
            .z = p.z,
        };
    }

    /// Multiply by the curve parameter chi = -0xd201000000010000.
    fn mulByChi(p: G2) G2 {
        // chi is public, so the variable-time multiplication is fine.
        return p.mulPublic(chi_abs, .big).neg();
    }

    /// Clear the cofactor, mapping any point of E'(Fp2) into G2.
    ///
    /// This is clear_cofactor_bls12381_g2 from RFC 9380, Appendix G.3.
    pub fn clearCofactor(p: G2) G2 {
        const t1 = p.mulByChi();
        var t2 = p.psi();
        var t3 = p.dbl().psi2();
        t3 = t3.sub(t2);
        t2 = t1.add(t2).mulByChi();
        t3 = t3.add(t2);
        t3 = t3.sub(t1);
        return t3.sub(p);
    }

    /// Return true if the point satisfies the curve equation.
    fn isOnCurve(p: G2) bool {
        if (p.z.isZero()) {
            return p.x.isZero() and !p.y.isZero();
        }
        // Y^2 Z = X^3 + b' Z^3
        const lhs = p.y.sq().mul(p.z);
        const rhs = p.x.sq().mul(p.x).add(B.mul(p.z.sq().mul(p.z)));
        return lhs.equivalent(rhs);
    }

    /// Return true if the point is on the curve and in G2, that is if [r]Q = O.
    pub fn isInSubgroup(p: G2) bool {
        if (!p.isOnCurve()) {
            return false;
        }
        // Multiplying by chi is much cheaper than by r = chi^4 - chi^2 + 1.
        const chi2 = p.mulByChi().mulByChi();
        const chi4 = chi2.mulByChi().mulByChi();
        return chi4.sub(chi2).add(p).isIdentity();
    }

    /// Return a random point in G2.
    pub fn random(io: std.Io) G2 {
        const s = scalar.random(io, .little);
        return basePoint.mul(s, .little);
    }
};

/// G2 point in affine coordinates.
pub const AffineCoordinates = struct {
    x: Fp2,
    y: Fp2,

    pub const identityElement = AffineCoordinates{ .x = Fp2.zero, .y = Fp2.zero };

    pub fn neg(p: AffineCoordinates) AffineCoordinates {
        return .{ .x = p.x, .y = p.y.neg() };
    }
};

/// |chi| as a scalar, where chi = -0xd201000000010000 is the curve parameter.
const chi_abs: [scalar.encoded_length]u8 = blk: {
    var bytes: [scalar.encoded_length]u8 = @splat(0);
    std.mem.writeInt(u64, bytes[bytes.len - 8 ..], 0xd201000000010000, .big);
    break :blk bytes;
};

// Psi endomorphism coefficients: 1/((1+u)^((p-1)/3)) and 1/((1+u)^((p-1)/2)).
const psi_coeff_x = Fp2{
    .c0 = Fp.zero,
    .c1 = Fp{
        .limbs = .{
            0x890dc9e4867545c3,
            0x2af322533285a5d5,
            0x50880866309b7e2c,
            0xa20d1b8c7e881024,
            0x14e4f04fe2db9068,
            0x14e56d3f1564853a,
        },
    },
};

const psi_coeff_y = Fp2{
    .c0 = Fp{
        .limbs = .{
            0x3e2f585da55c9ad1,
            0x4294213d86c18183,
            0x382844c88b623732,
            0x92ad2afd19103e18,
            0x1d794e4fac7cf0b9,
            0x0bd592fc7d825ec8,
        },
    },
    .c1 = Fp{
        .limbs = .{
            0x7bcfa7a25aa30fda,
            0xdc17dec12a927e7c,
            0x2f088dd86b4ebef1,
            0xd1ca2087da74d4a7,
            0x2da2596696cebc1d,
            0x0e2b7eedbbfd87d2,
        },
    },
};

// Psi^2 coefficient for x: 1/((1+u)^((p^2-1)/3)).
const psi2_coeff_x = Fp2{
    .c0 = Fp{
        .limbs = .{
            0xcd03c9e48671f071,
            0x5dab22461fcda5d2,
            0x587042afd3851b95,
            0x8eb60ebe01bacb9e,
            0x03f97d6e83d050d2,
            0x18f0206554638741,
        },
    },
    .c1 = Fp.zero,
};

/// Return a point of E'(Fp2) that is not in G2.
fn pointOutsideSubgroup() G2 {
    var x = Fp2.one;
    while (true) : (x = x.add(Fp2.one)) {
        const y = x.sq().mul(x).add(G2.B).sqrt() catch continue;
        const p = G2{ .x = x, .y = y, .z = Fp2.one };
        if (!p.isInSubgroup()) return p;
    }
}

test "g2 group law" {
    const g = G2.basePoint;
    const id = G2.identityElement;
    try std.testing.expect(g.isInSubgroup());
    try std.testing.expect(id.isInSubgroup());
    try std.testing.expect(!g.isIdentity());

    try std.testing.expect(g.add(id).equivalent(g));
    try std.testing.expect(id.add(g).equivalent(g));
    try std.testing.expect(g.sub(g).isIdentity());
    try std.testing.expect(g.dbl().equivalent(g.add(g)));
    try std.testing.expect(id.dbl().isIdentity());

    try std.testing.expect((try G2.fromAffineCoordinates(id.affineCoordinates())).isIdentity());
}

test "g2 scalar multiplication" {
    const g = G2.basePoint;
    const zero: scalar.CompressedScalar = @splat(0);
    try std.testing.expect(g.mul(zero, .big).isIdentity());

    var one = zero;
    one[one.len - 1] = 1;
    const r_minus_one = scalar.neg(one);
    try std.testing.expect(g.mul(r_minus_one, .big).equivalent(g.neg()));
    try std.testing.expect(g.mulPublic(r_minus_one, .big).equivalent(g.neg()));
}

test "g2 points outside the subgroup" {
    var off_curve = G2.basePoint;
    off_curve.y = off_curve.y.add(off_curve.y);
    try std.testing.expect(!off_curve.isInSubgroup());

    const p = pointOutsideSubgroup();
    try std.testing.expectError(error.InvalidEncoding, G2.fromCompressed(p.toCompressed()));
    try std.testing.expectError(error.InvalidEncoding, G2.fromUncompressed(p.toUncompressed()));
    try std.testing.expectError(error.InvalidEncoding, G2.fromAffineCoordinates(p.affineCoordinates()));

    const cleared = p.clearCofactor();
    try std.testing.expect(!cleared.isIdentity());
    try std.testing.expect(cleared.isInSubgroup());
}

test "g2 random point" {
    const io = std.testing.io;
    const p = G2.random(io);
    try std.testing.expect(!p.isIdentity());
    try std.testing.expect(p.isInSubgroup());
}
