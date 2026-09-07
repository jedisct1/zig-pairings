//! G1 group for BN462.
//!
//! G1 is the group of points on the curve E(Fp): y^2 = x^3 + 5
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

    /// Curve coefficient b = 5.
    pub const B = Fp.fromInt(5);

    /// The generator point.
    ///
    /// x = 0x21a6d67ef250191fadba34a0a30160b9ac9264b6f95f63b3edbec3cf4b2e689db1bbb4e69a416a0b1e79239c0372e5cd70113c98d91f36b6980d
    /// y = 0x0118ea0460f7f7abb82b33676a7432a490eeda842cccfa7d788c659650426e6af77df11b8ae40eb80f475432c66600622ecaa8a5734d36fb03de
    pub const basePoint = G1{
        .x = Fp{
            .limbs = .{
                0x60b833a0482110fc,
                0x2af93ed9abd13a14,
                0x6c44860453edecca,
                0x811fa9068dba1a13,
                0xe0a5dffddef391ed,
                0x2c6de500c0bdedd2,
                0x29a905e8abf39d13,
                0x0000000000000785,
            },
        },
        .y = Fp{
            .limbs = .{
                0x341e519317cb9454,
                0x9cc9f2c90251d92d,
                0x2f5528e0cf0b493e,
                0x152a8165fef04e6d,
                0xbf9e409a6384a317,
                0x02d5d0f1c2b68657,
                0x7b40e1d4cf6c2332,
                0x00000000000007be,
            },
        },
        .z = Fp.one,
        .is_base = true,
    };

    /// The identity element (point at infinity).
    pub const identityElement = G1{ .x = Fp.zero, .y = Fp.one, .z = Fp.zero };

    /// Compressed point encoding size (1 flag byte + 58 bytes x-coordinate).
    pub const compressed_length = 59;

    /// Uncompressed point encoding size (1 flag byte + 58 bytes x + 58 bytes y).
    pub const uncompressed_length = 117;

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
    pub fn fromSerializedAffineCoordinates(xs: [58]u8, ys: [58]u8, endian: std.builtin.Endian) (NonCanonicalError || EncodingError)!G1 {
        const x = try Fp.fromBytes(xs, endian);
        const y = try Fp.fromBytes(ys, endian);
        return fromAffineCoordinates(.{ .x = x, .y = y });
    }

    /// Deserialize from compressed form.
    ///
    /// The encoding is a flag byte followed by the x-coordinate:
    /// - Bit 7: always 1, marking the compressed form
    /// - Bit 6: 1 if point at infinity
    /// - Bit 5: sign of y, 1 if y is the lexicographically largest root
    /// - Bits 0-4: reserved, must be 0
    ///
    /// Allows identity and points outside the subgroup.
    /// Use isInSubgroup() and rejectIdentity() as required by the protocol.
    pub fn fromCompressed(bytes: [compressed_length]u8) (EncodingError || NotSquareError || NonCanonicalError)!G1 {
        const flags = bytes[0];
        const is_compressed = (flags & 0x80) != 0;
        const is_infinity = (flags & 0x40) != 0;
        const y_sign = (flags & 0x20) != 0;

        if (!is_compressed) {
            return error.InvalidEncoding;
        }

        // The reserved bits must be clear.
        if ((flags & 0x1f) != 0) {
            return error.InvalidEncoding;
        }

        if (is_infinity) {
            var all_zero = true;
            for (bytes[1..]) |b| {
                if (b != 0) all_zero = false;
            }
            if (!all_zero) {
                return error.InvalidEncoding;
            }
            return identityElement;
        }

        const x = try Fp.fromBytes(bytes[1..59].*, .big);

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
            bytes[0] = 0xC0; // compressed + infinity flags
        } else {
            bytes[0] = 0x80; // compressed flag
            if (affine.y.lexicographicallyLargest()) {
                bytes[0] |= 0x20; // y sign flag
            }
            bytes[1..59].* = affine.x.toBytes(.big);
        }

        return bytes;
    }

    /// Deserialize from uncompressed form.
    ///
    /// The encoding is a flag byte followed by the x and y coordinates:
    /// - Bit 7: always 0, marking the uncompressed form
    /// - Bit 6: 1 if point at infinity
    /// - Bits 0-5: reserved, must be 0
    ///
    /// Allows identity and points outside the subgroup.
    /// Use isInSubgroup() and rejectIdentity() as required by the protocol.
    pub fn fromUncompressed(bytes: [uncompressed_length]u8) (EncodingError || NonCanonicalError)!G1 {
        const flags = bytes[0];
        const is_compressed = (flags & 0x80) != 0;
        const is_infinity = (flags & 0x40) != 0;

        if (is_compressed) {
            return error.InvalidEncoding;
        }

        // The reserved bits must be clear, and there is no sign bit to carry.
        if ((flags & 0x3f) != 0) {
            return error.InvalidEncoding;
        }

        if (is_infinity) {
            var all_zero = true;
            for (bytes[1..]) |b| {
                if (b != 0) all_zero = false;
            }
            if (!all_zero) {
                return error.InvalidEncoding;
            }
            return identityElement;
        }

        const x = try Fp.fromBytes(bytes[1..59].*, .big);
        const y = try Fp.fromBytes(bytes[59..117].*, .big);

        return fromAffineCoordinates(.{ .x = x, .y = y });
    }

    /// Serialize to uncompressed form.
    pub fn toUncompressed(p: G1) [uncompressed_length]u8 {
        const affine = p.affineCoordinates();
        const is_identity = p.isIdentity();

        var bytes: [uncompressed_length]u8 = undefined;

        if (is_identity) {
            @memset(&bytes, 0);
            bytes[0] = 0x40; // infinity flag
        } else {
            bytes[0] = 0x00; // uncompressed, not infinity
            bytes[1..59].* = affine.x.toBytes(.big);
            bytes[59..117].* = affine.y.toBytes(.big);
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
    pub fn mul(p: G1, scalar_bytes: [scalar.encoded_length]u8, endian: std.builtin.Endian) IdentityElementError!G1 {
        const result = projective.mul(G1, p, &scalar_bytes, endian);

        try result.rejectIdentity();
        return result;
    }

    pub fn mulPublic(p: G1, scalar_bytes: [scalar.encoded_length]u8, endian: std.builtin.Endian) IdentityElementError!G1 {
        return p.mul(scalar_bytes, endian);
    }

    /// Validate subgroup membership, allowing the identity.
    pub fn isInSubgroup(p: G1) bool {
        if (!projective.isOnCurve(G1, p)) return false;
        return true; // BN462 G1 has cofactor 1.
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

test "g1 random point" {
    const io = std.testing.io;
    const p = G1.random(io);
    try std.testing.expect(!p.isIdentity());
    try std.testing.expect(p.isInSubgroup());
}
