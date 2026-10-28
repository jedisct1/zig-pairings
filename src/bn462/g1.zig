//! G1 group for BN462.
//!
//! G1 is the curve E(Fp): y^2 = x^3 + 5, which has prime order r.
//!
//! The specification doesn't define a point encoding for BN462.
//!
//! Points are encoded here with a leading flag byte holding the same three
//! flags as BLS12-381 (compressed 0x80, infinity 0x40, sign of y 0x20),
//! followed by the big-endian coordinates.
//!
//! Appendix E.4 of the specification describes this format.
//! It is not compatible with the SEC1 encoding of other BN462 libraries.

const std = @import("std");
const crypto = std.crypto;
const Fp = @import("fp.zig").Fp;
const scalar = @import("scalar.zig");

const EncodingError = crypto.errors.EncodingError;
const IdentityElementError = crypto.errors.IdentityElementError;
const NonCanonicalError = crypto.errors.NonCanonicalError;
const NotSquareError = crypto.errors.NotSquareError;

/// G1 point in projective coordinates (X, Y, Z) where x = X/Z, y = Y/Z.
pub const G1 = struct {
    x: Fp,
    y: Fp,
    z: Fp,

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
    };

    /// The identity element (point at infinity).
    pub const identityElement = G1{ .x = Fp.zero, .y = Fp.one, .z = Fp.zero };

    /// Compressed point encoding size (1 flag byte + 58 bytes x-coordinate).
    pub const compressed_length = 1 + Fp.encoded_length;

    /// Uncompressed point encoding size (1 flag byte + 58 bytes x + 58 bytes y).
    pub const uncompressed_length = 1 + 2 * Fp.encoded_length;

    /// Create a point from affine coordinates, checking that it is on the curve.
    ///
    /// `AffineCoordinates.identityElement` gives the identity element.
    pub fn fromAffineCoordinates(affine: AffineCoordinates) EncodingError!G1 {
        if (affine.x.isZero() and affine.y.isZero()) {
            return identityElement;
        }
        return fromCoordinates(affine.x, affine.y);
    }

    /// Build a point from its coordinates, checking that it is in G1.
    fn fromCoordinates(x: Fp, y: Fp) EncodingError!G1 {
        const p = G1{ .x = x, .y = y, .z = Fp.one };
        if (!p.isInSubgroup()) {
            return error.InvalidEncoding;
        }
        return p;
    }

    /// Create a point from serialized affine coordinates, checking that it is on the curve.
    pub fn fromSerializedAffineCoordinates(xs: [Fp.encoded_length]u8, ys: [Fp.encoded_length]u8, endian: std.builtin.Endian) (NonCanonicalError || EncodingError)!G1 {
        const x = try Fp.fromBytes(xs, endian);
        const y = try Fp.fromBytes(ys, endian);
        return fromAffineCoordinates(.{ .x = x, .y = y });
    }

    /// Deserialize a point in either form.
    ///
    /// The compression bit must match the length, so that a compressed G2
    /// point, which is also 117 bytes long, is not read as a G1 point.
    pub fn fromBytes(bytes: []const u8) (EncodingError || NotSquareError || NonCanonicalError)!G1 {
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
    /// - Bit 5: sign of y (sign_GF_p), 0 for the point at infinity
    /// - Bits 0-4: reserved, must be 0
    pub fn fromCompressed(bytes: [compressed_length]u8) (EncodingError || NotSquareError || NonCanonicalError)!G1 {
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

        const x = try Fp.fromBytes(bytes[1..compressed_length].*, .big);

        const x3_b = x.sq().mul(x).add(B);
        var y = try x3_b.sqrt();

        const y_is_largest = y.lexicographicallyLargest();
        if (y_sign != y_is_largest) {
            y = y.neg();
        }

        return fromCoordinates(x, y);
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
            bytes[1..compressed_length].* = affine.x.toBytes(.big);
        }

        return bytes;
    }

    /// Deserialize from uncompressed form.
    ///
    /// The encoding is a flag byte followed by the x and y coordinates:
    /// - Bit 7: always 0, marking the uncompressed form
    /// - Bit 6: 1 if point at infinity
    /// - Bits 0-5: reserved, must be 0
    pub fn fromUncompressed(bytes: [uncompressed_length]u8) (EncodingError || NonCanonicalError)!G1 {
        const flags = bytes[0];
        const is_compressed = (flags & 0x80) != 0;
        const is_infinity = (flags & 0x40) != 0;

        if (is_compressed) {
            return error.InvalidEncoding;
        }

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

        const x = try Fp.fromBytes(bytes[1 .. 1 + Fp.encoded_length].*, .big);
        const y = try Fp.fromBytes(bytes[1 + Fp.encoded_length ..].*, .big);

        return fromCoordinates(x, y);
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
            bytes[1 .. 1 + Fp.encoded_length].* = affine.x.toBytes(.big);
            bytes[1 + Fp.encoded_length ..].* = affine.y.toBytes(.big);
        }

        return bytes;
    }

    /// Return true if this is the identity element.
    pub fn isIdentity(p: G1) bool {
        return p.z.isZero();
    }

    /// Reject the identity element.
    ///
    /// Decoding accepts the identity, so protocols that don't want it should
    /// call this on the result.
    pub fn rejectIdentity(p: G1) IdentityElementError!void {
        if (p.isIdentity()) {
            return error.IdentityElement;
        }
    }

    /// Convert to affine coordinates.
    pub fn affineCoordinates(p: G1) AffineCoordinates {
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

    /// Swap two points when `choice` is 1 (constant-time).
    fn cSwap(p: *G1, q: *G1, choice: u1) void {
        const t = p.*;
        p.cMov(q.*, choice);
        q.cMov(t, choice);
    }

    /// Negate the point.
    pub fn neg(p: G1) G1 {
        return .{
            .x = p.x,
            .y = p.y.neg(),
            .z = p.z,
        };
    }

    /// Multiply by 3b = 15.
    fn mulByB3(x: Fp) Fp {
        return x.dbl().dbl().dbl().dbl().sub(x);
    }

    /// Double a point.
    ///
    /// Algorithm 9 from https://eprint.iacr.org/2015/1060.
    pub fn dbl(p: G1) G1 {
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
    pub fn add(p: G1, q: G1) G1 {
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
    pub fn sub(p: G1, q: G1) G1 {
        return p.add(q.neg());
    }

    /// Multiply a point by a scalar, in constant time.
    ///
    /// The scalar doesn't have to be reduced, and zero gives the identity.
    pub fn mul(p: G1, scalar_bytes: [scalar.encoded_length]u8, endian: std.builtin.Endian) G1 {
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
    pub fn mulPublic(p: G1, scalar_bytes: [scalar.encoded_length]u8, endian: std.builtin.Endian) G1 {
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

    /// Return true if the point is in G1.
    ///
    /// The cofactor is 1, so this only needs to check the curve equation.
    pub fn isInSubgroup(p: G1) bool {
        if (p.z.isZero()) {
            return p.x.isZero() and !p.y.isZero();
        }
        // Y^2 Z = X^3 + b Z^3
        const lhs = p.y.sq().mul(p.z);
        const rhs = p.x.sq().mul(p.x).add(B.mul(p.z.sq().mul(p.z)));
        return lhs.equivalent(rhs);
    }

    /// Return a random point in G1.
    pub fn random(io: std.Io) G1 {
        const s = scalar.random(io, .little);
        return basePoint.mul(s, .little);
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

test "g1 group law" {
    const g = G1.basePoint;
    const id = G1.identityElement;
    try std.testing.expect(g.isInSubgroup());
    try std.testing.expect(id.isInSubgroup());
    try std.testing.expect(!g.isIdentity());

    // E(Fp) has prime order r, which is why the curve check is enough.
    try std.testing.expect(g.mulPublic(scalar.modulus_bytes, .big).isIdentity());

    try std.testing.expect(g.add(id).equivalent(g));
    try std.testing.expect(id.add(g).equivalent(g));
    try std.testing.expect(g.sub(g).isIdentity());
    try std.testing.expect(g.dbl().equivalent(g.add(g)));
    try std.testing.expect(id.dbl().isIdentity());

    try std.testing.expect((try G1.fromAffineCoordinates(id.affineCoordinates())).isIdentity());
}

test "g1 scalar multiplication" {
    const g = G1.basePoint;
    const zero: scalar.CompressedScalar = @splat(0);
    try std.testing.expect(g.mul(zero, .big).isIdentity());

    var one = zero;
    one[one.len - 1] = 1;
    const r_minus_one = scalar.neg(one);
    try std.testing.expect(g.mul(r_minus_one, .big).equivalent(g.neg()));
    try std.testing.expect(g.mulPublic(r_minus_one, .big).equivalent(g.neg()));
}

test "g1 encoding round trip" {
    const points = [_]G1{ G1.basePoint, G1.basePoint.neg(), G1.identityElement };
    for (points) |p| {
        try std.testing.expect(p.equivalent(try G1.fromCompressed(p.toCompressed())));
        try std.testing.expect(p.equivalent(try G1.fromUncompressed(p.toUncompressed())));
        try std.testing.expect(p.equivalent(try G1.fromBytes(&p.toCompressed())));
        try std.testing.expect(p.equivalent(try G1.fromBytes(&p.toUncompressed())));
    }
    try std.testing.expectEqual(0xc0, G1.identityElement.toCompressed()[0]);
    try std.testing.expectEqual(0x40, G1.identityElement.toUncompressed()[0]);
}

test "g1 rejects malformed encodings" {
    // Same flag rules as BLS12-381, and the five low bits are reserved.
    var bad_identity = G1.identityElement.toCompressed();
    bad_identity[0] = 0xe0;
    try std.testing.expectError(error.InvalidEncoding, G1.fromCompressed(bad_identity));
    var bad_reserved = G1.basePoint.toCompressed();
    bad_reserved[0] |= 0x01;
    try std.testing.expectError(error.InvalidEncoding, G1.fromCompressed(bad_reserved));
    var bad_uncompressed = G1.basePoint.toUncompressed();
    bad_uncompressed[0] = 0x20;
    try std.testing.expectError(error.InvalidEncoding, G1.fromUncompressed(bad_uncompressed));
    var bad_uncompressed_identity = G1.identityElement.toUncompressed();
    bad_uncompressed_identity[0] = 0x60;
    try std.testing.expectError(error.InvalidEncoding, G1.fromUncompressed(bad_uncompressed_identity));

    // A SEC1 type byte is not accepted.
    bad_uncompressed[0] = 0x04;
    try std.testing.expectError(error.InvalidEncoding, G1.fromUncompressed(bad_uncompressed));

    try std.testing.expectError(error.InvalidEncoding, G1.fromBytes(&.{}));

    // Coordinates must be below p and on the curve.
    var non_canonical = G1.basePoint.toUncompressed();
    _ = try std.fmt.hexToBytes(non_canonical[1..59], "240480360120023ffffffffff6ff0cf6b7d9bfca0000000000d812908f41c8020ffffffffff6ff66fc6ff687f640000000002401b00840138013");
    try std.testing.expectError(error.NonCanonical, G1.fromUncompressed(non_canonical));
    var off_curve = G1.basePoint.toUncompressed();
    off_curve[116] ^= 1;
    try std.testing.expectError(error.InvalidEncoding, G1.fromUncompressed(off_curve));
    const all_zero: [G1.uncompressed_length]u8 = @splat(0);
    try std.testing.expectError(error.InvalidEncoding, G1.fromUncompressed(all_zero));
}

test "g1 random point" {
    const io = std.testing.io;
    const p = G1.random(io);
    try std.testing.expect(!p.isIdentity());
    try std.testing.expect(p.isInSubgroup());
}
