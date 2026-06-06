//! G2 group for BN462.
//!
//! G2 is the group of points on the twisted curve E'(Fp2): y^2 = x^3 + (2 - u)
//!
//! This is a subgroup of order r of the full curve.
//!
//! BN462 uses a D-type twist.

const std = @import("std");
const crypto = std.crypto;
const Fp = @import("fp.zig").Fp;
const Fp2 = @import("fp2.zig").Fp2;
const scalar = @import("scalar.zig");

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

    /// Curve coefficient b' = (2 - u) for the twisted curve.
    pub const B = Fp2{
        .c0 = Fp.fromInt(2),
        .c1 = Fp.fromInt(1).neg(),
    };

    /// The generator point, as given by the IETF specification.
    pub const basePoint = G2{
        .x = Fp2{
            .c0 = Fp{
                .limbs = .{
                    0x76e30f7cc9ee8adc,
                    0xf81795d7672b1122,
                    0xfd72bb9fcfd06c09,
                    0x706e716d58081d77,
                    0x35545e784b044a95,
                    0x271c39c8cc4b5253,
                    0x626fd151b2e46f87,
                    0x0000000000002172,
                },
            },
            .c1 = Fp{
                .limbs = .{
                    0x987f5d0946aa5fb7,
                    0x0e249077d244e595,
                    0x8edb1b350eae7a44,
                    0x3d8ec0f53d225897,
                    0xdab10a2c4973fe32,
                    0xfb4e1be8dd040a3d,
                    0x86f7dab6b4a6f8ee,
                    0x000000000000006f,
                },
            },
        },
        .y = Fp2{
            .c0 = Fp{
                .limbs = .{
                    0xf9db9dd1920106f4,
                    0xa6e72663a34786ab,
                    0x24ba07e7a716c9fa,
                    0xcedaad780bacddbe,
                    0x24a70d3c8b1bad0f,
                    0x964fb960a9f2ff79,
                    0xd722823a70bdeade,
                    0x00000000000002c1,
                },
            },
            .c1 = Fp{
                .limbs = .{
                    0x6ea78171f0466a2b,
                    0xd2f85160a598e710,
                    0x7a9f91eb3c5e1cb9,
                    0x3dc63e25940e7057,
                    0x3437fc781eb949c7,
                    0xd21f8a07cf799962,
                    0x2ff68b889e358905,
                    0x0000000000001b03,
                },
            },
        },
        .z = Fp2.one,
        .is_base = true,
    };

    /// The identity element (point at infinity).
    pub const identityElement = G2{ .x = Fp2.zero, .y = Fp2.one, .z = Fp2.zero };

    /// Compressed point encoding size (1 flag byte + 116 bytes for Fp2 x-coordinate).
    pub const compressed_length = 117;

    /// Uncompressed point encoding size (1 flag byte + 116 bytes x + 116 bytes y).
    pub const uncompressed_length = 233;

    /// Create a point from affine coordinates, checking it's on the curve.
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

    /// Create a point from serialized affine coordinates.
    pub fn fromSerializedAffineCoordinates(xs: [Fp2.encoded_length]u8, ys: [Fp2.encoded_length]u8, endian: std.builtin.Endian) (NonCanonicalError || EncodingError)!G2 {
        const x = try Fp2.fromBytes(xs, endian);
        const y = try Fp2.fromBytes(ys, endian);
        return fromAffineCoordinates(.{ .x = x, .y = y });
    }

    /// Deserialize from compressed form.
    ///
    /// The encoding is a flag byte followed by the x-coordinate:
    /// - Bit 7: always 1, marking the compressed form
    /// - Bit 6: 1 if point at infinity
    /// - Bit 5: sign of y, 1 if y is the lexicographically largest root
    /// - Bits 0-4: reserved, must be 0
    pub fn fromCompressed(bytes: [compressed_length]u8) (EncodingError || NotSquareError || NonCanonicalError)!G2 {
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

        const x = try Fp2.fromBytes(bytes[1..117].*, .big);

        const x3_b = x.sq().mul(x).add(B);
        var y = try x3_b.sqrt();

        const y_is_largest = y.lexicographicallyLargest();
        if (y_sign != y_is_largest) {
            y = y.neg();
        }

        return G2{ .x = x, .y = y, .z = Fp2.one };
    }

    /// Serialize to compressed form.
    pub fn toCompressed(p: G2) [compressed_length]u8 {
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
            bytes[1..117].* = affine.x.toBytes(.big);
        }

        return bytes;
    }

    /// Deserialize from uncompressed form.
    ///
    /// The encoding is a flag byte followed by the x and y coordinates, with
    /// the flag byte laid out as for the compressed form minus the sign bit.
    pub fn fromUncompressed(bytes: [uncompressed_length]u8) (EncodingError || NonCanonicalError)!G2 {
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

        const x = try Fp2.fromBytes(bytes[1..117].*, .big);
        const y = try Fp2.fromBytes(bytes[117..233].*, .big);

        return fromAffineCoordinates(.{ .x = x, .y = y });
    }

    /// Serialize to uncompressed form.
    pub fn toUncompressed(p: G2) [uncompressed_length]u8 {
        const affine = p.affineCoordinates();
        const is_identity = p.isIdentity();

        var bytes: [uncompressed_length]u8 = undefined;

        if (is_identity) {
            @memset(&bytes, 0);
            bytes[0] = 0x40; // infinity flag
        } else {
            bytes[0] = 0x00;
            bytes[1..117].* = affine.x.toBytes(.big);
            bytes[117..233].* = affine.y.toBytes(.big);
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

    /// Negate the point.
    pub fn neg(p: G2) G2 {
        return .{
            .x = p.x,
            .y = p.y.neg(),
            .z = p.z,
        };
    }

    /// Double a point using the general addition formula.
    pub fn dbl(p: G2) G2 {
        return addGeneral(p, p);
    }

    /// General point addition that also handles P + P (doubling).
    fn addGeneral(p: G2, q: G2) G2 {
        if (p.isIdentity()) return q;
        if (q.isIdentity()) return p;

        // Working in affine coordinates keeps the special cases easy to spot.
        const a = p.affineCoordinates();
        const b = q.affineCoordinates();

        if (a.x.equivalent(b.x)) {
            if (a.y.equivalent(b.y)) {
                return dblAffine(a);
            } else {
                return identityElement;
            }
        }

        const lambda = b.y.sub(a.y).mul(b.x.sub(a.x).invert());
        const x3 = lambda.sq().sub(a.x).sub(b.x);
        const y3 = lambda.mul(a.x.sub(x3)).sub(a.y);

        return G2{ .x = x3, .y = y3, .z = Fp2.one };
    }

    /// Affine doubling formula.
    fn dblAffine(p: AffineCoordinates) G2 {
        const xx = p.x.sq();
        const lambda = xx.add(xx).add(xx).mul(p.y.dbl().invert());
        const x3 = lambda.sq().sub(p.x.dbl());
        const y3 = lambda.mul(p.x.sub(x3)).sub(p.y);

        return G2{ .x = x3, .y = y3, .z = Fp2.one };
    }

    /// Add two points.
    pub fn add(p: G2, q: G2) G2 {
        return addGeneral(p, q);
    }

    /// Subtract two points.
    pub fn sub(p: G2, q: G2) G2 {
        return p.add(q.neg());
    }

    /// Add a point in affine coordinates.
    pub fn addMixed(p: G2, q: AffineCoordinates) G2 {
        return p.add(G2{ .x = q.x, .y = q.y, .z = Fp2.one });
    }

    /// Subtract a point in affine coordinates.
    pub fn subMixed(p: G2, q: AffineCoordinates) G2 {
        return p.add(G2{ .x = q.x, .y = q.y.neg(), .z = Fp2.one });
    }

    /// Scalar multiplication.
    pub fn mul(p: G2, scalar_bytes: [scalar.encoded_length]u8, endian: std.builtin.Endian) IdentityElementError!G2 {
        const s = if (endian == .little) scalar_bytes else blk: {
            var swapped: [scalar.encoded_length]u8 = undefined;
            for (scalar_bytes, 0..) |b, i| swapped[scalar.encoded_length - 1 - i] = b;
            break :blk swapped;
        };

        var result = identityElement;
        var temp = p;

        for (s) |byte| {
            var b = byte;
            for (0..8) |_| {
                if (b & 1 == 1) {
                    result = result.add(temp);
                }
                temp = temp.dbl();
                b >>= 1;
            }
        }

        try result.rejectIdentity();
        return result;
    }

    /// Scalar multiplication with public scalar (variable time).
    pub fn mulPublic(p: G2, scalar_bytes: [scalar.encoded_length]u8, endian: std.builtin.Endian) IdentityElementError!G2 {
        return p.mul(scalar_bytes, endian);
    }

    /// Check if the point is in the correct subgroup G2.
    pub fn isInSubgroup(p: G2) bool {
        if (p.isIdentity()) return true;

        const affine = p.affineCoordinates();
        const x3_b = affine.x.sq().mul(affine.x).add(B);
        const y2 = affine.y.sq();
        return y2.equivalent(x3_b);
    }

    /// Return a random point in G2.
    pub fn random(io: std.Io) G2 {
        const s = scalar.random(io, .little);
        return basePoint.mul(s, .little) catch unreachable;
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

test "g2 random point" {
    const io = std.testing.io;
    const p = G2.random(io);
    try std.testing.expect(!p.isIdentity());
    try std.testing.expect(p.isInSubgroup());
}
