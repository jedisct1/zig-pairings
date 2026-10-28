//! G1 group for BLS12-381.
//!
//! G1 is the subgroup of order r of the curve E(Fp): y^2 = x^3 + 4
//!
//! Points decoded from bytes or coordinates are always checked to be in G1.

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
    };

    /// The identity element (point at infinity).
    pub const identityElement = G1{ .x = Fp.zero, .y = Fp.one, .z = Fp.zero };

    /// Compressed point encoding size.
    pub const compressed_length = 48;

    /// Uncompressed point encoding size.
    pub const uncompressed_length = 96;

    /// Create a point from affine coordinates, checking that it is in G1.
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

    /// Create a point from serialized affine coordinates, checking that it is in G1.
    pub fn fromSerializedAffineCoordinates(xs: [Fp.encoded_length]u8, ys: [Fp.encoded_length]u8, endian: std.builtin.Endian) (NonCanonicalError || EncodingError)!G1 {
        const x = try Fp.fromBytes(xs, endian);
        const y = try Fp.fromBytes(ys, endian);
        return fromAffineCoordinates(.{ .x = x, .y = y });
    }

    /// Deserialize a point in either form, like octets_to_point_G1 in the specification.
    ///
    /// The compression bit must match the length, so that a compressed G2
    /// point, which is also 96 bytes long, is not read as a G1 point.
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
    /// The encoding is the x-coordinate with flags in the top three bits:
    /// - Bit 7: always 1, marking the compressed form
    /// - Bit 6: 1 if point at infinity
    /// - Bit 5: sign of y, as defined by sign_GF_p
    ///
    /// The point is checked to be in G1.
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

        return fromCoordinates(x, y);
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
    /// The point is checked to be in G1.
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

        return fromCoordinates(x, y);
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

    /// Multiply by 3b = 12.
    fn mulByB3(x: Fp) Fp {
        return x.dbl().add(x).dbl().dbl();
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

    /// Multiply by the curve parameter t = -0xd201000000010000.
    fn mulByT(p: G1) G1 {
        // t is public, so the variable-time multiplication is fine.
        return p.mulPublic(t_abs, .big).neg();
    }

    /// Clear the cofactor, mapping any point of E(Fp) into G1.
    ///
    /// Like RFC 9380, this multiplies by h_eff = 1 - t rather than by the
    /// cofactor itself.
    pub fn clearCofactor(p: G1) G1 {
        return p.sub(p.mulByT());
    }

    /// Return true if the point satisfies the curve equation.
    fn isOnCurve(p: G1) bool {
        if (p.z.isZero()) {
            return p.x.isZero() and !p.y.isZero();
        }
        // Y^2 Z = X^3 + b Z^3
        const lhs = p.y.sq().mul(p.z);
        const rhs = p.x.sq().mul(p.x).add(B.mul(p.z.sq().mul(p.z)));
        return lhs.equivalent(rhs);
    }

    /// Return true if the point is on the curve and in G1, that is if [r]P = O.
    pub fn isInSubgroup(p: G1) bool {
        if (!p.isOnCurve()) {
            return false;
        }
        // r = t^4 - t^2 + 1, and multiplying by t is much cheaper than by r.
        const t2 = p.mulByT().mulByT();
        const t4 = t2.mulByT().mulByT();
        return t4.sub(t2).add(p).isIdentity();
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

/// |t| as a scalar, where t = -0xd201000000010000 is the curve parameter.
const t_abs: [scalar.encoded_length]u8 = blk: {
    var bytes: [scalar.encoded_length]u8 = @splat(0);
    std.mem.writeInt(u64, bytes[bytes.len - 8 ..], 0xd201000000010000, .big);
    break :blk bytes;
};

/// Return a point of E(Fp) that is not in G1.
fn pointOutsideSubgroup() G1 {
    var x = Fp.one;
    while (true) : (x = x.add(Fp.one)) {
        const y = x.sq().mul(x).add(G1.B).sqrt() catch continue;
        const p = G1{ .x = x, .y = y, .z = Fp.one };
        if (!p.isInSubgroup()) return p;
    }
}

test "g1 group law" {
    const g = G1.basePoint;
    const id = G1.identityElement;
    try std.testing.expect(g.isInSubgroup());
    try std.testing.expect(id.isInSubgroup());
    try std.testing.expect(!g.isIdentity());

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

test "g1 rejects malformed encodings" {
    // The sign bit is only allowed on compressed points other than the identity.
    var bad_identity = G1.identityElement.toCompressed();
    bad_identity[0] = 0xe0;
    try std.testing.expectError(error.InvalidEncoding, G1.fromCompressed(bad_identity));
    var bad_uncompressed = G1.basePoint.toUncompressed();
    bad_uncompressed[0] |= 0x20;
    try std.testing.expectError(error.InvalidEncoding, G1.fromUncompressed(bad_uncompressed));
    var bad_uncompressed_identity = G1.identityElement.toUncompressed();
    bad_uncompressed_identity[0] = 0x60;
    try std.testing.expectError(error.InvalidEncoding, G1.fromUncompressed(bad_uncompressed_identity));

    // The identity must be all zeros after the flags.
    var identity_with_payload = G1.identityElement.toCompressed();
    identity_with_payload[47] = 1;
    try std.testing.expectError(error.InvalidEncoding, G1.fromCompressed(identity_with_payload));

    // The length must match the compression bit.
    const compressed = G1.basePoint.toCompressed();
    try std.testing.expectError(error.InvalidEncoding, G1.fromBytes(&.{}));
    try std.testing.expectError(error.InvalidEncoding, G1.fromBytes(compressed[0..47]));
}

test "g1 points outside the subgroup" {
    var off_curve = G1.basePoint;
    off_curve.y = off_curve.y.add(off_curve.y);
    try std.testing.expect(!off_curve.isInSubgroup());

    const p = pointOutsideSubgroup();
    try std.testing.expectError(error.InvalidEncoding, G1.fromCompressed(p.toCompressed()));
    try std.testing.expectError(error.InvalidEncoding, G1.fromUncompressed(p.toUncompressed()));
    try std.testing.expectError(error.InvalidEncoding, G1.fromAffineCoordinates(p.affineCoordinates()));

    const cleared = p.clearCofactor();
    try std.testing.expect(!cleared.isIdentity());
    try std.testing.expect(cleared.isInSubgroup());
}

test "g1 random point" {
    const io = std.testing.io;
    const p = G1.random(io);
    try std.testing.expect(!p.isIdentity());
    try std.testing.expect(p.isInSubgroup());
}
