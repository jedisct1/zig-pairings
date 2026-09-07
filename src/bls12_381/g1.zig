//! G1 group for BLS12-381.
//!
//! G1 is the group of points on the curve E(Fp): y^2 = x^3 + 4
//!
//! This is a subgroup of order r of the full curve.

const std = @import("std");
const crypto = std.crypto;
const Fp = @import("fp.zig").Fp;
const scalar = @import("scalar.zig");
const projective = @import("../projective.zig");

const EncodingError = crypto.errors.EncodingError;
const IdentityElementError = crypto.errors.IdentityElementError;
const NonCanonicalError = crypto.errors.NonCanonicalError;
const NotSquareError = crypto.errors.NotSquareError;

/// G1 point in projective coordinates (X, Y, Z) where x = X/Z, y = Y/Z.
pub const G1 = struct {
    x: Fp,
    y: Fp,
    z: Fp,

    is_base: bool = false,

    /// Curve coefficient b = 4.
    pub const B = Fp.fromInt(4);

    /// The generator point.
    pub const basePoint = G1{
        // x = 0x17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb
        .x = Fp{
            .limbs = .{
                0x5cb38790fd530c16,
                0x7817fc679976fff5,
                0x154f95c7143ba1c1,
                0xf0ae6acdf3d0e747,
                0xedce6ecc21dbf440,
                0x120177419e0bfb75,
            },
        },
        // y = 0x08b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af600db18cb2c04b3edd03cc744a2888ae40caa232946c5e7e1
        .y = Fp{
            .limbs = .{
                0xbaac93d50ce72271,
                0x8c22631a7918fd8e,
                0xdd595f13570725ce,
                0x51ac582950405194,
                0x0e1c8c3fad0059c0,
                0x0bbc3efc5008a26a,
            },
        },
        .z = Fp.one,
        .is_base = true,
    };

    /// The identity element (point at infinity).
    pub const identityElement = G1{ .x = Fp.zero, .y = Fp.one, .z = Fp.zero };

    /// Compressed point encoding size.
    pub const compressed_length = 48;

    /// Uncompressed point encoding size.
    pub const uncompressed_length = 96;

    /// Create a curve point without checking subgroup membership.
    pub fn fromAffineCoordinates(affine: AffineCoordinates) EncodingError!G1 {
        const x = affine.x;
        const y = affine.y;

        const x3 = x.sq().mul(x);
        const y2 = y.sq();
        const rhs = x3.add(B);

        if (!y2.equivalent(rhs)) {
            return error.InvalidEncoding;
        }

        return G1{ .x = x, .y = y, .z = Fp.one };
    }

    /// Create a point from serialized affine coordinates.
    pub fn fromSerializedAffineCoordinates(xs: [48]u8, ys: [48]u8, endian: std.builtin.Endian) (NonCanonicalError || EncodingError)!G1 {
        const x = try Fp.fromBytes(xs, endian);
        const y = try Fp.fromBytes(ys, endian);
        return fromAffineCoordinates(.{ .x = x, .y = y });
    }

    /// Deserialize from compressed form.
    ///
    /// The encoding is the x-coordinate with flags in the top three bits:
    /// - Bit 7: always 1, marking the compressed form
    /// - Bit 6: 1 if point at infinity
    /// - Bit 5: sign of y, 1 if y is the lexicographically largest root
    ///
    /// Allows identity and points outside the subgroup.
    /// Use isInSubgroup() and rejectIdentity() as required by the protocol.
    pub fn fromCompressed(bytes: [compressed_length]u8) (EncodingError || NotSquareError || NonCanonicalError)!G1 {
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
        const x = try Fp.fromBytes(x_bytes, .big);

        const x3_b = x.sq().mul(x).add(B);
        var y = try x3_b.sqrt();

        const y_is_largest = y.lexicographicallyLargest();
        if (y_sign != y_is_largest) {
            y = y.neg();
        }

        return G1{ .x = x, .y = y, .z = Fp.one };
    }

    /// Serialize to compressed form.
    pub fn toCompressed(p: G1) [compressed_length]u8 {
        const affine = p.affineCoordinates();
        const is_identity = p.isIdentity();

        var bytes: [compressed_length]u8 = undefined;

        if (is_identity) {
            @memset(&bytes, 0);
            bytes[0] = 0b11000000; // compressed + infinity flags
        } else {
            bytes = affine.x.toBytes(.big);
            const y_sign: u8 = if (affine.y.lexicographicallyLargest()) 1 else 0;
            bytes[0] |= 0b10000000 | (y_sign << 5);
        }

        return bytes;
    }

    /// Deserialize from uncompressed form.
    ///
    /// The encoding is the x and y coordinates back to back, with flags in the
    /// top bits of the first byte:
    /// - Bit 7: always 0, marking the uncompressed form
    /// - Bit 6: 1 if point at infinity
    ///
    /// Allows identity and points outside the subgroup.
    /// Use isInSubgroup() and rejectIdentity() as required by the protocol.
    pub fn fromUncompressed(bytes: [uncompressed_length]u8) (EncodingError || NonCanonicalError)!G1 {
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

        var x_bytes = bytes[0..48].*;
        x_bytes[0] &= 0x1f;
        const x = try Fp.fromBytes(x_bytes, .big);
        const y = try Fp.fromBytes(bytes[48..96].*, .big);

        return fromAffineCoordinates(.{ .x = x, .y = y });
    }

    /// Serialize to uncompressed form.
    pub fn toUncompressed(p: G1) [uncompressed_length]u8 {
        const affine = p.affineCoordinates();
        const is_identity = p.isIdentity();

        var bytes: [uncompressed_length]u8 = undefined;

        if (is_identity) {
            @memset(&bytes, 0);
            bytes[0] = 0b01000000; // infinity flag
        } else {
            bytes[0..48].* = affine.x.toBytes(.big);
            bytes[48..96].* = affine.y.toBytes(.big);
        }

        return bytes;
    }

    /// Return true if this is the identity element.
    pub fn isIdentity(p: G1) bool {
        return p.z.isZero();
    }

    /// Reject the identity element.
    pub fn rejectIdentity(p: G1) IdentityElementError!void {
        if (p.isIdentity()) {
            return error.IdentityElement;
        }
    }

    /// Convert to affine coordinates.
    pub fn affineCoordinates(p: G1) AffineCoordinates {
        const z_inv = p.z.invert();
        return .{
            .x = p.x.mul(z_inv),
            .y = p.y.mul(z_inv),
        };
    }

    /// Return true if both points are equivalent.
    pub fn equivalent(a: G1, b: G1) bool {
        // Cross-multiply rather than normalizing, which would cost two inversions.
        const x1z2 = a.x.mul(b.z);
        const x2z1 = b.x.mul(a.z);
        const y1z2 = a.y.mul(b.z);
        const y2z1 = b.y.mul(a.z);

        return x1z2.equivalent(x2z1) and y1z2.equivalent(y2z1);
    }

    /// Conditionally replace the point with `q` when `choice` is 1 (constant-time).
    pub fn cMov(p: *G1, q: G1, choice: u1) void {
        p.x.cMov(q.x, choice);
        p.y.cMov(q.y, choice);
        p.z.cMov(q.z, choice);
    }

    /// Negate the point.
    pub fn neg(p: G1) G1 {
        return .{
            .x = p.x,
            .y = p.y.neg(),
            .z = p.z,
        };
    }

    pub fn dbl(p: G1) G1 {
        return projective.dbl(G1, p);
    }

    /// Add two points.
    pub fn add(p: G1, q: G1) G1 {
        return projective.add(G1, p, q);
    }

    /// Subtract two points.
    pub fn sub(p: G1, q: G1) G1 {
        return p.add(q.neg());
    }

    /// Add a point in affine coordinates.
    pub fn addMixed(p: G1, q: AffineCoordinates) G1 {
        var other = G1{ .x = q.x, .y = q.y, .z = Fp.one };
        const is_identity = @intFromBool(q.x.isZero()) & @intFromBool(q.y.isZero());
        other.cMov(identityElement, is_identity);
        return p.add(other);
    }

    /// Subtract a point in affine coordinates.
    pub fn subMixed(p: G1, q: AffineCoordinates) G1 {
        return p.addMixed(q.neg());
    }

    /// Multiply by a secret scalar, rejecting an identity result.
    pub fn mul(p: G1, scalar_bytes: [32]u8, endian: std.builtin.Endian) IdentityElementError!G1 {
        const result = projective.mul(G1, p, &scalar_bytes, endian);

        try result.rejectIdentity();
        return result;
    }

    pub fn mulPublic(p: G1, scalar_bytes: [32]u8, endian: std.builtin.Endian) IdentityElementError!G1 {
        return p.mul(scalar_bytes, endian);
    }

    /// Clear the cofactor to ensure the point is in G1.
    pub fn clearCofactor(p: G1) G1 {
        // The cofactor, 0x396c8c005555e1568c00aaab0000aaab.
        const h_bytes: [32]u8 = .{
            0xab, 0xaa, 0x00, 0x00, 0xab, 0xaa, 0x00, 0x8c,
            0x56, 0xe1, 0x55, 0x55, 0x00, 0x8c, 0x6c, 0x39,
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        };

        return p.mul(h_bytes, .little) catch identityElement;
    }

    /// Validate subgroup membership, allowing the identity.
    pub fn isInSubgroup(p: G1) bool {
        if (!projective.isOnCurve(G1, p)) return false;
        return projective.mul(G1, p, &scalar.modulus_bytes, .big).isIdentity();
    }

    /// Return a random point in G1.
    pub fn random(io: std.Io) G1 {
        while (true) {
            const s = scalar.random(io, .big);
            return basePoint.mul(s, .big) catch continue;
        }
    }
};

/// G1 point in affine coordinates.
pub const AffineCoordinates = struct {
    x: Fp,
    y: Fp,

    pub const identityElement = AffineCoordinates{ .x = Fp.zero, .y = Fp.zero };

    pub fn neg(p: AffineCoordinates) AffineCoordinates {
        return .{ .x = p.x, .y = p.y.neg() };
    }
};

test "g1 generator is on curve" {
    const g = G1.basePoint;
    const x3_b = g.x.sq().mul(g.x).add(G1.B);
    const y2 = g.y.sq();
    try std.testing.expect(x3_b.equivalent(y2));
}

test "g1 identity" {
    try std.testing.expect(G1.identityElement.isIdentity());
    try std.testing.expect(!G1.basePoint.isIdentity());
}

test "g1 addition" {
    const g = G1.basePoint;

    const result = g.add(G1.identityElement);
    try std.testing.expect(result.equivalent(g));

    const neg_g = g.neg();
    const zero = g.add(neg_g);
    try std.testing.expect(zero.isIdentity());
}

test "g1 doubling" {
    const g = G1.basePoint;
    const double_g = g.dbl();
    const g_plus_g = g.add(g);
    try std.testing.expect(double_g.equivalent(g_plus_g));
}

test "g1 serialization roundtrip" {
    const g = G1.basePoint;

    const compressed = g.toCompressed();
    const decompressed = try G1.fromCompressed(compressed);
    try std.testing.expect(g.equivalent(decompressed));

    const uncompressed = g.toUncompressed();
    const from_uncompressed = try G1.fromUncompressed(uncompressed);
    try std.testing.expect(g.equivalent(from_uncompressed));
}

test "g1 identity encoding roundtrip" {
    const id = G1.identityElement;

    const compressed = id.toCompressed();
    try std.testing.expectEqual(0xc0, compressed[0]);
    try std.testing.expect((try G1.fromCompressed(compressed)).isIdentity());

    const uncompressed = id.toUncompressed();
    try std.testing.expectEqual(0x40, uncompressed[0]);
    try std.testing.expect((try G1.fromUncompressed(uncompressed)).isIdentity());
}

test "g1 rejects invalid metadata bits" {
    // The sign bit must never be set on an uncompressed or point-at-infinity
    // encoding, so 0x20, 0x60, and 0xe0 in the leading byte are all invalid.

    // Compressed infinity with the sign bit set.
    var bad_compressed_inf = G1.identityElement.toCompressed();
    bad_compressed_inf[0] = 0xe0;
    try std.testing.expectError(error.InvalidEncoding, G1.fromCompressed(bad_compressed_inf));

    // Uncompressed generator with the sign bit set.
    var bad_uncompressed = G1.basePoint.toUncompressed();
    bad_uncompressed[0] |= 0x20;
    try std.testing.expectError(error.InvalidEncoding, G1.fromUncompressed(bad_uncompressed));

    // Uncompressed infinity with the sign bit set.
    var bad_uncompressed_inf = G1.identityElement.toUncompressed();
    bad_uncompressed_inf[0] = 0x60;
    try std.testing.expectError(error.InvalidEncoding, G1.fromUncompressed(bad_uncompressed_inf));
}

test "g1 random point" {
    const io = std.testing.io;
    const p = G1.random(io);
    try std.testing.expect(!p.isIdentity());
    try std.testing.expect(p.isInSubgroup());
}
