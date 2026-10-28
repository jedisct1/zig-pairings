//! G2 group for BN462.
//!
//! G2 is the subgroup of order r of the twisted curve E'(Fp2): y^2 = x^3 + (2 - u).
//! The specification calls it G'_2.
//! BN462 uses a D-type twist.
//!
//! Points decoded from bytes or coordinates are always checked to be in G2.
//! They are encoded like G1 points.

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
    };

    /// The identity element (point at infinity).
    pub const identityElement = G2{ .x = Fp2.zero, .y = Fp2.one, .z = Fp2.zero };

    /// Compressed point encoding size (1 flag byte + 116 bytes for Fp2 x-coordinate).
    pub const compressed_length = 1 + Fp2.encoded_length;

    /// Uncompressed point encoding size (1 flag byte + 116 bytes x + 116 bytes y).
    pub const uncompressed_length = 1 + 2 * Fp2.encoded_length;

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

    /// Create a point from serialized affine coordinates, checking that it is in G2.
    pub fn fromSerializedAffineCoordinates(xs: [Fp2.encoded_length]u8, ys: [Fp2.encoded_length]u8, endian: std.builtin.Endian) (NonCanonicalError || EncodingError)!G2 {
        const x = try Fp2.fromBytes(xs, endian);
        const y = try Fp2.fromBytes(ys, endian);
        return fromAffineCoordinates(.{ .x = x, .y = y });
    }

    /// Deserialize a point in either form.
    ///
    /// The compression bit must match the length, so that an uncompressed G1
    /// point, which is also 117 bytes long, is not read as a G2 point.
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
    /// The encoding is a flag byte followed by the x-coordinate:
    /// - Bit 7: always 1, marking the compressed form
    /// - Bit 6: 1 if point at infinity
    /// - Bit 5: sign of y (sign_GF_p^2), 0 for the point at infinity
    /// - Bits 0-4: reserved, must be 0
    ///
    /// The point is checked to be in G2.
    pub fn fromCompressed(bytes: [compressed_length]u8) (EncodingError || NotSquareError || NonCanonicalError)!G2 {
        const flags = bytes[0];
        const is_compressed = (flags & 0x80) != 0;
        const is_infinity = (flags & 0x40) != 0;
        const y_sign = (flags & 0x20) != 0;

        if (!is_compressed) {
            return error.InvalidEncoding;
        }

        // Reserved bits must be clear, and the point at infinity has no sign.
        if ((flags & 0x1f) != 0 or (is_infinity and y_sign)) {
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

        const x = try Fp2.fromBytes(bytes[1..compressed_length].*, .big);

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
            bytes[0] = 0xC0; // compressed + infinity flags
        } else {
            bytes[0] = 0x80; // compressed flag
            if (affine.y.lexicographicallyLargest()) {
                bytes[0] |= 0x20; // y sign flag
            }
            bytes[1..compressed_length].* = affine.x.toBytes(.big);
        }

        return bytes;
    }

    /// Deserialize from uncompressed form.
    ///
    /// The encoding is a flag byte followed by the x and y coordinates, with
    /// the flag byte laid out as for the compressed form minus the sign bit.
    ///
    /// The point is checked to be in G2.
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

        const x = try Fp2.fromBytes(bytes[1 .. 1 + Fp2.encoded_length].*, .big);
        const y = try Fp2.fromBytes(bytes[1 + Fp2.encoded_length ..].*, .big);

        return fromCoordinates(x, y);
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
            bytes[1 .. 1 + Fp2.encoded_length].* = affine.x.toBytes(.big);
            bytes[1 + Fp2.encoded_length ..].* = affine.y.toBytes(.big);
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

    /// Multiply by 3b' = 6 - 3u.
    fn mulByB3(x: Fp2) Fp2 {
        // (a + bu)(6 - 3u) = (6a + 3b) + (6b - 3a)u
        const a3 = x.c0.dbl().add(x.c0);
        const b3 = x.c1.dbl().add(x.c1);
        return .{
            .c0 = a3.dbl().add(b3),
            .c1 = b3.dbl().sub(a3),
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
        return p.isOnCurve() and p.mulPublic(scalar.modulus_bytes, .big).isIdentity();
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

/// Return a point of E'(Fp2) that is not in G2.
fn pointOutsideSubgroup() G2 {
    var x = Fp2.one;
    while (true) : (x = x.add(Fp2.one)) {
        const y = x.sq().mul(x).add(G2.B).sqrt() catch continue;
        const p = G2{ .x = x, .y = y, .z = Fp2.one };
        if (!p.isInSubgroup()) return p;
    }
}

/// Multiply by the cofactor h' given in Section 4.2.2 of the specification.
fn mulByCofactor(p: G2) G2 {
    var h: [scalar.encoded_length]u8 = undefined;
    _ = std.fmt.hexToBytes(&h, "240480360120023ffffffffff6ff0cf6b7d9bfca0000000000d812908fa1ce0227fffffffff6ff66fc63f5f7f4c0000000002401b008a0168019") catch unreachable;
    return p.mulPublic(h, .big);
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
    const affine = p.affineCoordinates();
    try std.testing.expectError(error.InvalidEncoding, G2.fromSerializedAffineCoordinates(affine.x.toBytes(.big), affine.y.toBytes(.big), .big));
    try std.testing.expectError(error.InvalidEncoding, G2.fromAffineCoordinates(affine));

    const cleared = mulByCofactor(p);
    try std.testing.expect(!cleared.isIdentity());
    try std.testing.expect(cleared.isInSubgroup());
}

test "g2 random point" {
    const io = std.testing.io;
    const p = G2.random(io);
    try std.testing.expect(!p.isIdentity());
    try std.testing.expect(p.isInSubgroup());
}
