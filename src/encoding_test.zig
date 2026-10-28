//! Point encoding tests shared by both curves.

const std = @import("std");
const bls12_381 = @import("bls12_381.zig");
const bn462 = @import("bn462.zig");

test "point encodings" {
    try expectStrictEncoding(bls12_381.G1, bls12_381.Fp);
    try expectStrictEncoding(bls12_381.G2, bls12_381.Fp);
    try expectStrictEncoding(bn462.G1, bn462.Fp);
    try expectStrictEncoding(bn462.G2, bn462.Fp);
}

/// Check that points survive a round trip, and that each point has a single
/// encoding in each form, which Section 5.6.1 of the draft counts on.
fn expectStrictEncoding(comptime G: type, comptime Fp: type) !void {
    // BLS12-381 keeps the flags in the top three bits of the first coordinate,
    // and BN462 in a byte of their own.
    const header_len = G.compressed_length % Fp.encoded_length;
    const flags_mask: u8 = if (header_len == 0) 0xe0 else 0xff;

    try std.testing.expectEqual(0xc0, G.identityElement.toCompressed()[0]);
    try std.testing.expectEqual(0x40, G.identityElement.toUncompressed()[0]);

    // Changing the flags must make the encoding invalid, except for the sign
    // bit, which gives the negated point.
    // That leaves three valid encodings for P and for -P, and two for the
    // identity.
    var decoded: usize = 0;
    for ([_]G{ G.basePoint, G.basePoint.neg(), G.identityElement }) |p| {
        inline for (.{ p.toCompressed(), p.toUncompressed() }) |encoded| {
            for (0..256) |v| {
                const delta: u8 = @intCast(v);
                if (delta & ~flags_mask != 0) continue;
                var bytes = encoded;
                bytes[0] ^= delta;
                const q = G.fromBytes(&bytes) catch continue;
                if (delta == 0) try std.testing.expect(q.equivalent(p));
                const again = if (bytes.len == G.compressed_length) q.toCompressed() else q.toUncompressed();
                try std.testing.expectEqualSlices(u8, &bytes, &again);
                decoded += 1;
            }
        }
    }
    try std.testing.expectEqual(8, decoded);

    // A coordinate equal to p is rejected wherever it appears.
    // p is odd, so adding one to p - 1 only changes the last byte.
    var p_bytes = Fp.one.neg().toBytes(.big);
    p_bytes[p_bytes.len - 1] += 1;
    inline for (.{ G.basePoint.toCompressed(), G.basePoint.toUncompressed() }) |encoded| {
        var slot: usize = header_len;
        while (slot < encoded.len) : (slot += Fp.encoded_length) {
            var bytes = encoded;
            bytes[slot..][0..Fp.encoded_length].* = p_bytes;
            // Put the flags back if p overwrote them.
            bytes[0] |= encoded[0] & flags_mask;
            try std.testing.expectError(error.NonCanonical, G.fromBytes(&bytes));
        }
    }

    // The identity has nothing but zeros after the flags.
    var identity = G.identityElement.toCompressed();
    identity[identity.len - 1] = 1;
    try std.testing.expectError(error.InvalidEncoding, G.fromBytes(&identity));

    var off_curve = G.basePoint.toUncompressed();
    off_curve[off_curve.len - 1] ^= 1;
    try std.testing.expectError(error.InvalidEncoding, G.fromBytes(&off_curve));

    // The length has to match the compression flag.
    const compressed = G.basePoint.toCompressed();
    try std.testing.expectError(error.InvalidEncoding, G.fromBytes(compressed[0 .. compressed.len - 1]));
    try std.testing.expectError(error.InvalidEncoding, G.fromBytes(&.{}));
}
