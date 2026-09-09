//! G2 group for BLS12-381.
//!
//! G2 is the group of points on the twisted curve E'(Fp2): y^2 = x^3 + 4(1 + u)
//!
//! This is a subgroup of order r of the full curve.

const std = @import("std");
const crypto = std.crypto;
const Fp = @import("fp.zig").Fp;
const Fp2 = @import("fp2.zig").Fp2;
const scalar = @import("scalar.zig");
const projective = @import("../projective.zig");

const EncodingError = crypto.errors.EncodingError;
const IdentityElementError = crypto.errors.IdentityElementError;
const NonCanonicalError = crypto.errors.NonCanonicalError;
const NotSquareError = crypto.errors.NotSquareError;

/// G2 point in projective coordinates.
pub const G2 = struct {
    x: Fp2,
    y: Fp2,
    z: Fp2,

    is_base: bool = false,

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
        .is_base = true,
    };

    /// The identity element (point at infinity).
    pub const identityElement = G2{ .x = Fp2.zero, .y = Fp2.one, .z = Fp2.zero };

    /// Compressed point encoding size.
    pub const compressed_length = 96;

    /// Uncompressed point encoding size.
    pub const uncompressed_length = 192;

    /// Create a curve point without checking subgroup membership.
    pub fn fromAffineCoordinates(affine: AffineCoordinates) EncodingError!G2 {
        const x = affine.x;
        const y = affine.y;

        const x3 = x.sq().mul(x);
        const y2 = y.sq();
        const rhs = x3.add(B);

        if (!y2.equivalent(rhs)) {
            return error.InvalidEncoding;
        }

        return G2{ .x = x, .y = y, .z = Fp2.one };
    }

    /// Deserialize from compressed form.
    ///
    /// The encoding is the x-coordinate with flags in the top three bits,
    /// laid out the same way as for G1.
    ///
    /// Checks subgroup membership and allows identity.
    /// Use rejectIdentity() if the protocol requires it.
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

        const p = G2{ .x = x, .y = y, .z = Fp2.one };
        if (!p.isInSubgroup()) return error.InvalidEncoding;
        return p;
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
    /// Checks subgroup membership and allows identity.
    /// Use rejectIdentity() if the protocol requires it.
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

        const p = try fromAffineCoordinates(.{ .x = x, .y = y });
        if (!p.isInSubgroup()) return error.InvalidEncoding;
        return p;
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
    pub fn rejectIdentity(p: G2) IdentityElementError!void {
        if (p.isIdentity()) {
            return error.IdentityElement;
        }
    }

    /// Convert to affine coordinates.
    pub fn affineCoordinates(p: G2) AffineCoordinates {
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

    /// Negate the point.
    pub fn neg(p: G2) G2 {
        return .{
            .x = p.x,
            .y = p.y.neg(),
            .z = p.z,
        };
    }

    pub fn dbl(p: G2) G2 {
        return projective.dbl(G2, p);
    }

    /// Add two points.
    pub fn add(p: G2, q: G2) G2 {
        return projective.add(G2, p, q);
    }

    /// Subtract two points.
    pub fn sub(p: G2, q: G2) G2 {
        return p.add(q.neg());
    }

    /// Add a point in affine coordinates.
    pub fn addMixed(p: G2, q: AffineCoordinates) G2 {
        var other = G2{ .x = q.x, .y = q.y, .z = Fp2.one };
        const is_identity = @intFromBool(q.x.isZero()) & @intFromBool(q.y.isZero());
        other.cMov(identityElement, is_identity);
        return p.add(other);
    }

    /// Subtract a point in affine coordinates.
    pub fn subMixed(p: G2, q: AffineCoordinates) G2 {
        return p.addMixed(q.neg());
    }

    /// Multiply by a secret scalar, rejecting an identity result.
    pub fn mul(p: G2, scalar_bytes: [32]u8, endian: std.builtin.Endian) IdentityElementError!G2 {
        const result = projective.mul(G2, p, &scalar_bytes, endian);

        try result.rejectIdentity();
        return result;
    }

    pub fn mulPublic(p: G2, scalar_bytes: [32]u8, endian: std.builtin.Endian) IdentityElementError!G2 {
        return p.mul(scalar_bytes, endian);
    }

    /// Endomorphism used for cofactor clearing.
    pub fn psi(p: G2) G2 {
        return G2{
            .x = p.x.conjugate().mul(psi_coeff_x),
            .y = p.y.conjugate().mul(psi_coeff_y),
            .z = p.z.conjugate(),
        };
    }

    pub fn psi2(p: G2) G2 {
        return G2{
            .x = p.x.mul(psi2_coeff_x),
            .y = p.y.neg(),
            .z = p.z,
        };
    }

    /// Map a curve point into G2 using the effective cofactor from RFC 9380.
    pub fn clearCofactor(p: G2) G2 {
        const t1 = p.mulByX();
        const t2 = p.psi();
        const t3 = p.dbl().psi2().sub(t2);
        return t3.add(t1.add(t2).mulByX()).sub(t1).sub(p);
    }

    /// Multiply by the signed BLS parameter.
    fn mulByX(p: G2) G2 {
        const t_abs = [8]u8{ 0xd2, 0x01, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00 };
        return projective.mul(G2, p, &t_abs, .big).neg();
    }

    /// Validate subgroup membership, allowing the identity.
    pub fn isInSubgroup(p: G2) bool {
        if (!projective.isOnCurve(G2, p)) return false;
        return projective.mul(G2, p, &scalar.modulus_bytes, .big).isIdentity();
    }

    /// Return a random point in G2.
    pub fn random(io: std.Io) G2 {
        while (true) {
            const s = scalar.random(io, .big);
            return basePoint.mul(s, .big) catch continue;
        }
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

test "g2 generator is on curve" {
    const g = G2.basePoint;
    const x3_b = g.x.sq().mul(g.x).add(G2.B);
    const y2 = g.y.sq();
    try std.testing.expect(x3_b.equivalent(y2));
}

test "g2 identity" {
    try std.testing.expect(G2.identityElement.isIdentity());
    try std.testing.expect(!G2.basePoint.isIdentity());
}

test "g2 addition" {
    const g = G2.basePoint;

    const result = g.add(G2.identityElement);
    try std.testing.expect(result.equivalent(g));

    const neg_g = g.neg();
    const zero = g.add(neg_g);
    try std.testing.expect(zero.isIdentity());
}

test "g2 doubling" {
    const g = G2.basePoint;
    const double_g = g.dbl();
    const g_plus_g = g.add(g);
    try std.testing.expect(double_g.equivalent(g_plus_g));
}

test "g2 serialization roundtrip" {
    const g = G2.basePoint;

    const compressed = g.toCompressed();
    const decompressed = try G2.fromCompressed(compressed);
    try std.testing.expect(g.equivalent(decompressed));

    const uncompressed = g.toUncompressed();
    const from_uncompressed = try G2.fromUncompressed(uncompressed);
    try std.testing.expect(g.equivalent(from_uncompressed));
}

test "g2 identity encoding roundtrip" {
    const id = G2.identityElement;

    const compressed = id.toCompressed();
    try std.testing.expectEqual(0xc0, compressed[0]);
    try std.testing.expect((try G2.fromCompressed(compressed)).isIdentity());

    const uncompressed = id.toUncompressed();
    try std.testing.expectEqual(0x40, uncompressed[0]);
    try std.testing.expect((try G2.fromUncompressed(uncompressed)).isIdentity());
}

test "g2 rejects invalid metadata bits" {
    // The sign bit must never be set on an uncompressed or point-at-infinity
    // encoding, so 0x20, 0x60, and 0xe0 in the leading byte are all invalid.

    // Compressed infinity with the sign bit set.
    var bad_compressed_inf = G2.identityElement.toCompressed();
    bad_compressed_inf[0] = 0xe0;
    try std.testing.expectError(error.InvalidEncoding, G2.fromCompressed(bad_compressed_inf));

    // Uncompressed generator with the sign bit set.
    var bad_uncompressed = G2.basePoint.toUncompressed();
    bad_uncompressed[0] |= 0x20;
    try std.testing.expectError(error.InvalidEncoding, G2.fromUncompressed(bad_uncompressed));

    // Uncompressed infinity with the sign bit set.
    var bad_uncompressed_inf = G2.identityElement.toUncompressed();
    bad_uncompressed_inf[0] = 0x60;
    try std.testing.expectError(error.InvalidEncoding, G2.fromUncompressed(bad_uncompressed_inf));
}

test "g2 random point" {
    const io = std.testing.io;
    const p = G2.random(io);
    try std.testing.expect(!p.isIdentity());

    const affine = p.affineCoordinates();
    const x3_b = affine.x.sq().mul(affine.x).add(G2.B);
    const y2 = affine.y.sq();
    try std.testing.expect(y2.equivalent(x3_b));
}
