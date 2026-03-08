//! Test vectors verifying that this implementation conforms to the IETF
//! specification for BLS12-381.

const std = @import("std");
const Fp = @import("fp.zig").Fp;
const Fp2 = @import("fp2.zig").Fp2;
const G1 = @import("g1.zig").G1;
const G2 = @import("g2.zig").G2;
const pairing = @import("pairing.zig");
const scalar = @import("scalar.zig");

/// Parse a hex string at compile time into bytes.
fn parseHex(comptime hex: []const u8) [hex.len / 2]u8 {
    var bytes: [hex.len / 2]u8 = undefined;
    _ = std.fmt.hexToBytes(&bytes, hex) catch unreachable;
    return bytes;
}

// Curve parameters, as given by the IETF specification.

// Field modulus p, hex without the 0x prefix.
const p_hex = "1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab";

// Subgroup order r.
const r_hex = "73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001";

// G1 generator x coordinate.
const g1_x_hex = "17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb";

// G1 generator y coordinate.
const g1_y_hex = "08b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af600db18cb2c04b3edd03cc744a2888ae40caa232946c5e7e1";

// G2 generator x' = x'_0 + x'_1 * u.
const g2_x0_hex = "024aa2b2f08f0a91260805272dc51051c6e47ad4fa403b02b4510b647ae3d1770bac0326a805bbefd48056c8c121bdb8";
const g2_x1_hex = "13e02b6052719f607dacd3a088274f65596bd0d09920b61ab5da61bbdc7f5049334cf11213945d57e5ac7d055d042b7e";

// G2 generator y' = y'_0 + y'_1 * u.
const g2_y0_hex = "0ce5d527727d6e118cc9cdc6da2e351aadfd9baa8cbdd3a76d429a695160d12c923ac9cc3baca289e193548608b82801";
const g2_y1_hex = "0606c4a02ea734cc32acd2b02bc28b99cb3e287e85a763af267492ab572e99ab3f370d275cec1da1aaa9075ff05f79be";

// The e_0 through e_11 components of e(BP, BP'), the pairing of the generators.
const e0_hex = "11619b45f61edfe3b47a15fac19442526ff489dcda25e59121d9931438907dfd448299a87dde3a649bdba96e84d54558";
const e1_hex = "153ce14a76a53e205ba8f275ef1137c56a566f638b52d34ba3bf3bf22f277d70f76316218c0dfd583a394b8448d2be7f";
const e2_hex = "095668fb4a02fe930ed44767834c915b283b1c6ca98c047bd4c272e9ac3f3ba6ff0b05a93e59c71fba77bce995f04692";
const e3_hex = "16deedaa683124fe7260085184d88f7d036b86f53bb5b7f1fc5e248814782065413e7d958d17960109ea006b2afdeb5f";
const e4_hex = "09c92cf02f3cd3d2f9d34bc44eee0dd50314ed44ca5d30ce6a9ec0539be7a86b121edc61839ccc908c4bdde256cd6048";
const e5_hex = "111061f398efc2a97ff825b04d21089e24fd8b93a47e41e60eae7e9b2a38d54fa4dedced0811c34ce528781ab9e929c7";
const e6_hex = "01ecfcf31c86257ab00b4709c33f1c9c4e007659dd5ffc4a735192167ce197058cfb4c94225e7f1b6c26ad9ba68f63bc";
const e7_hex = "08890726743a1f94a8193a166800b7787744a8ad8e2f9365db76863e894b7a11d83f90d873567e9d645ccf725b32d26f";
const e8_hex = "0e61c752414ca5dfd258e9606bac08daec29b3e2c57062669556954fb227d3f1260eedf25446a086b0844bcd43646c10";
const e9_hex = "0fe63f185f56dd29150fc498bbeea78969e7e783043620db33f75a05a0a2ce5c442beaff9da195ff15164c00ab66bdde";
const e10_hex = "10900338a92ed0b47af211636f7cfdec717b7ee43900eee9b5fc24f0000c5874d4801372db478987691c566a8c474978";
const e11_hex = "1454814f3085f0e6602247671bc408bbce2007201536818c901dbd4d2095dd86c1ec8b888e59611f60a301af7776be3d";

const FpBytes = [48]u8;

test "modulus matches IETF specification" {
    const expected: [48]u8 = parseHex(p_hex);

    var actual: [48]u8 = undefined;
    inline for (0..6) |i| {
        const j = 5 - i;
        std.mem.writeInt(u64, actual[j * 8 ..][0..8], Fp.modulus[i], .big);
    }

    try std.testing.expectEqualSlices(u8, &expected, &actual);
}

test "subgroup order matches IETF specification" {
    const expected: [32]u8 = parseHex(r_hex);
    try std.testing.expectEqualSlices(u8, &expected, &scalar.modulus_bytes);
}

test "G1 generator matches IETF specification" {
    const expected_x = parseHex(g1_x_hex);
    const expected_y = parseHex(g1_y_hex);

    const g1 = G1.basePoint.affineCoordinates();
    const actual_x = g1.x.toBytes(.big);
    const actual_y = g1.y.toBytes(.big);

    try std.testing.expectEqualSlices(u8, &expected_x, &actual_x);
    try std.testing.expectEqualSlices(u8, &expected_y, &actual_y);
}

test "G2 generator matches IETF specification" {
    const expected_x0 = parseHex(g2_x0_hex);
    const expected_x1 = parseHex(g2_x1_hex);
    const expected_y0 = parseHex(g2_y0_hex);
    const expected_y1 = parseHex(g2_y1_hex);

    const g2 = G2.basePoint.affineCoordinates();
    const actual_x0 = g2.x.c0.toBytes(.big);
    const actual_x1 = g2.x.c1.toBytes(.big);
    const actual_y0 = g2.y.c0.toBytes(.big);
    const actual_y1 = g2.y.c1.toBytes(.big);

    try std.testing.expectEqualSlices(u8, &expected_x0, &actual_x0);
    try std.testing.expectEqualSlices(u8, &expected_x1, &actual_x1);
    try std.testing.expectEqualSlices(u8, &expected_y0, &actual_y0);
    try std.testing.expectEqualSlices(u8, &expected_y1, &actual_y1);
}

test "curve equation b = 4 matches specification" {
    // The specification gives E: y^2 = x^3 + 4 and E': y^2 = x^3 + 4(u + 1).
    try std.testing.expect(G1.B.equivalent(Fp.fromInt(4)));
    try std.testing.expect(G2.B.c0.equivalent(Fp.fromInt(4)));
    try std.testing.expect(G2.B.c1.equivalent(Fp.fromInt(4)));
}

test "extension field tower matches specification" {
    // The specification defines:
    // GF(p^2) = GF(p)[u] / (u^2 + 1)
    // GF(p^6) = GF(p^2)[v] / (v^3 - (u + 1))
    // GF(p^12) = GF(p^6)[w] / (w^2 - v)

    // Verify u^2 = -1 in Fp2.
    const u = Fp2{ .c0 = Fp.zero, .c1 = Fp.one };
    const u_sq = u.sq();
    try std.testing.expect(u_sq.c0.equivalent(Fp.one.neg()));
    try std.testing.expect(u_sq.c1.isZero());

    // Verify the non-residue used to build Fp6 is (1 + u).
    const nonres = Fp2{ .c0 = Fp.one, .c1 = Fp.one };
    const test_val = Fp2.fromInts(3, 5);
    const by_nonres = test_val.mulByNonresidue();
    const direct = test_val.mul(nonres);
    try std.testing.expect(by_nonres.equivalent(direct));
}

test "pairing e(BP, BP') matches the IETF test vector" {
    const expected = [_]FpBytes{
        parseHex(e0_hex),
        parseHex(e1_hex),
        parseHex(e2_hex),
        parseHex(e3_hex),
        parseHex(e4_hex),
        parseHex(e5_hex),
        parseHex(e6_hex),
        parseHex(e7_hex),
        parseHex(e8_hex),
        parseHex(e9_hex),
        parseHex(e10_hex),
        parseHex(e11_hex),
    };

    const result = pairing.pair(G1.basePoint, G2.basePoint);

    const actual = [_]FpBytes{
        result.c0.c0.c0.toBytes(.big),
        result.c0.c0.c1.toBytes(.big),
        result.c0.c1.c0.toBytes(.big),
        result.c0.c1.c1.toBytes(.big),
        result.c0.c2.c0.toBytes(.big),
        result.c0.c2.c1.toBytes(.big),
        result.c1.c0.c0.toBytes(.big),
        result.c1.c0.c1.toBytes(.big),
        result.c1.c1.c0.toBytes(.big),
        result.c1.c1.c1.toBytes(.big),
        result.c1.c2.c0.toBytes(.big),
        result.c1.c2.c1.toBytes(.big),
    };

    inline for (0..12) |i| {
        try std.testing.expectEqualSlices(u8, &expected[i], &actual[i]);
    }
}

test "BLS parameter x matches specification" {
    // The specification gives t = -2^63 - 2^62 - 2^60 - 2^57 - 2^48 - 2^16.
    //
    // Its absolute value is what the Miller loop iterates over.
    const bls_x: u64 = 0xd201000000010000;

    const expected: u64 = (1 << 63) + (1 << 62) + (1 << 60) +
        (1 << 57) + (1 << 48) + (1 << 16);

    try std.testing.expectEqual(expected, bls_x);
}

test "serialization format follows the specification" {
    const g1 = G1.basePoint;
    const compressed = g1.toCompressed();

    // Compressed, and not the point at infinity.
    try std.testing.expect((compressed[0] & 0x80) != 0);
    try std.testing.expect((compressed[0] & 0x40) == 0);

    const decompressed = try G1.fromCompressed(compressed);
    try std.testing.expect(g1.equivalent(decompressed));

    const g2 = G2.basePoint;
    const g2_compressed = g2.toCompressed();

    try std.testing.expect((g2_compressed[0] & 0x80) != 0);
    try std.testing.expect((g2_compressed[0] & 0x40) == 0);

    const g2_decompressed = try G2.fromCompressed(g2_compressed);
    try std.testing.expect(g2.equivalent(g2_decompressed));
}

test "compressed base points match the reference vectors" {
    const g1_compressed = "97f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb";
    const g2_compressed =
        "93e02b6052719f607dacd3a088274f65596bd0d09920b61ab5da61bbdc7f5049" ++
        "334cf11213945d57e5ac7d055d042b7e024aa2b2f08f0a91260805272dc51051" ++
        "c6e47ad4fa403b02b4510b647ae3d1770bac0326a805bbefd48056c8c121bdb8";

    try std.testing.expectEqualSlices(u8, &parseHex(g1_compressed), &G1.basePoint.toCompressed());
    try std.testing.expectEqualSlices(u8, &parseHex(g2_compressed), &G2.basePoint.toCompressed());

    // The encodings round-trip back to the base points.
    try std.testing.expect((try G1.fromCompressed(parseHex(g1_compressed))).equivalent(G1.basePoint));
    try std.testing.expect((try G2.fromCompressed(parseHex(g2_compressed))).equivalent(G2.basePoint));
}

test "identity serialization follows the specification" {
    const id = G1.identityElement;
    const compressed = id.toCompressed();

    // Compressed, and flagged as the point at infinity.
    try std.testing.expect((compressed[0] & 0x80) != 0);
    try std.testing.expect((compressed[0] & 0x40) != 0);
    for (compressed[1..]) |b| {
        try std.testing.expect(b == 0);
    }
    try std.testing.expect((compressed[0] & 0x1f) == 0);
}

test "Fp2 serialization order (c1 || c0 for big endian)" {
    // The specification serializes x_0 + x_1*u as x_1 followed by x_0.
    const elem = Fp2.fromInts(123, 456);
    const bytes = elem.toBytes(.big);

    const c1_bytes = elem.c1.toBytes(.big);
    const c0_bytes = elem.c0.toBytes(.big);

    try std.testing.expectEqualSlices(u8, &c1_bytes, bytes[0..48]);
    try std.testing.expectEqualSlices(u8, &c0_bytes, bytes[48..96]);

    const parsed = try Fp2.fromBytes(bytes, .big);
    try std.testing.expect(parsed.equivalent(elem));
}
