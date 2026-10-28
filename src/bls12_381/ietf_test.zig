//! Known-answer tests for BLS12-381, taken from revision 14 of the
//! pairing-friendly curves draft and from RFC 9380.

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

/// Build a G1 point from hex coordinates, without any check.
fn testG1(comptime x_hex: []const u8, comptime y_hex: []const u8) !G1 {
    return .{
        .x = try Fp.fromBytes(parseHex(x_hex), .big),
        .y = try Fp.fromBytes(parseHex(y_hex), .big),
        .z = Fp.one,
    };
}

/// Build an Fp2 element from the hex of its two coefficients.
fn testFp2(comptime c0_hex: []const u8, comptime c1_hex: []const u8) !Fp2 {
    return .{
        .c0 = try Fp.fromBytes(parseHex(c0_hex), .big),
        .c1 = try Fp.fromBytes(parseHex(c1_hex), .big),
    };
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

    // The subgroup checks rely on r = chi^4 - chi^2 + 1.
    const chi: u256 = 0xd201000000010000;
    try std.testing.expectEqual(chi * chi * chi * chi - chi * chi + 1, std.mem.readInt(u256, &scalar.modulus_bytes, .big));
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

test "Appendix D.1 point vectors" {
    const g1_compressed = parseHex("97f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb");
    const g1_uncompressed = parseHex(g1_x_hex ++ g1_y_hex);
    const g2_compressed = parseHex("93e02b6052719f607dacd3a088274f65596bd0d09920b61ab5da61bbdc7f5049334cf11213945d57e5ac7d055d042b7e" ++
        "024aa2b2f08f0a91260805272dc51051c6e47ad4fa403b02b4510b647ae3d1770bac0326a805bbefd48056c8c121bdb8");
    const g2_uncompressed = parseHex(g2_x1_hex ++ g2_x0_hex ++ g2_y1_hex ++ g2_y0_hex);

    try std.testing.expectEqualSlices(u8, &g1_compressed, &G1.basePoint.toCompressed());
    try std.testing.expectEqualSlices(u8, &g1_uncompressed, &G1.basePoint.toUncompressed());
    try std.testing.expectEqualSlices(u8, &g2_compressed, &G2.basePoint.toCompressed());
    try std.testing.expectEqualSlices(u8, &g2_uncompressed, &G2.basePoint.toUncompressed());

    try std.testing.expect((try G1.fromBytes(&g1_compressed)).equivalent(G1.basePoint));
    try std.testing.expect((try G1.fromBytes(&g1_uncompressed)).equivalent(G1.basePoint));
    try std.testing.expect((try G2.fromBytes(&g2_compressed)).equivalent(G2.basePoint));
    try std.testing.expect((try G2.fromBytes(&g2_uncompressed)).equivalent(G2.basePoint));

    // A compressed G2 point and an uncompressed G1 point are both 96 bytes long.
    try std.testing.expectError(error.InvalidEncoding, G1.fromBytes(&g2_compressed));
    try std.testing.expectError(error.InvalidEncoding, G1.fromUncompressed(g2_compressed));
    try std.testing.expectError(error.InvalidEncoding, G2.fromBytes(&g1_uncompressed));
    try std.testing.expectError(error.InvalidEncoding, G2.fromCompressed(g1_uncompressed));

    // The identity is 0xc0 then zeros when compressed, and 0x40 then zeros
    // when uncompressed.
    var e_compressed: [48]u8 = @splat(0);
    e_compressed[0] = 0xc0;
    var e_uncompressed: [96]u8 = @splat(0);
    e_uncompressed[0] = 0x40;
    var e2_compressed: [96]u8 = @splat(0);
    e2_compressed[0] = 0xc0;
    var e2_uncompressed: [192]u8 = @splat(0);
    e2_uncompressed[0] = 0x40;

    try std.testing.expectEqualSlices(u8, &e_compressed, &G1.identityElement.toCompressed());
    try std.testing.expectEqualSlices(u8, &e_uncompressed, &G1.identityElement.toUncompressed());
    try std.testing.expectEqualSlices(u8, &e2_compressed, &G2.identityElement.toCompressed());
    try std.testing.expectEqualSlices(u8, &e2_uncompressed, &G2.identityElement.toUncompressed());

    try std.testing.expect((try G1.fromBytes(&e_compressed)).isIdentity());
    try std.testing.expect((try G1.fromBytes(&e_uncompressed)).isIdentity());
    try std.testing.expect((try G2.fromBytes(&e2_compressed)).isIdentity());
    try std.testing.expect((try G2.fromBytes(&e2_uncompressed)).isIdentity());
    try std.testing.expectError(error.IdentityElement, (try G1.fromBytes(&e_compressed)).rejectIdentity());
}

test "Appendix D.3 scalar vectors" {
    const zero = parseHex("0000000000000000000000000000000000000000000000000000000000000000");
    const one = parseHex("0000000000000000000000000000000000000000000000000000000000000001");
    const r_minus_one = parseHex("73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000000");
    const r = parseHex(r_hex);

    try std.testing.expect(scalar.isZero(try scalar.fromBytes(zero, .big)));
    try std.testing.expectEqualSlices(u8, &one, &try scalar.fromBytes(one, .big));
    try std.testing.expectEqualSlices(u8, &r_minus_one, &try scalar.fromBytes(r_minus_one, .big));
    try std.testing.expectError(error.NonCanonical, scalar.fromBytes(r, .big));
    try std.testing.expectError(error.NonCanonical, scalar.fromBytes(@splat(0xff), .big));

    // Scalars are big-endian: read as little-endian, the encoding of 1 is 2^248.
    try std.testing.expectEqualSlices(u8, &one, &scalar.toBytes(one, .big));
    try std.testing.expect(!std.mem.eql(u8, &one, &try scalar.fromBytes(one, .little)));
}

test "RFC 9380 G1 cofactor clearing" {
    // Q0, Q1 and P = clear_cofactor(Q0 + Q1) for msg = "" and "abc", from RFC 9380, Appendix J.9.1.
    const vectors = [_][6][]const u8{
        .{
            "11a3cce7e1d90975990066b2f2643b9540fa40d6137780df4e753a8054d07580db3b7f1f03396333d4a359d1fe3766fe",
            "0eeaf6d794e479e270da10fdaf768db4c96b650a74518fc67b04b03927754bac66f3ac720404f339ecdcc028afa091b7",
            "160003aaf1632b13396dbad518effa00fff532f604de1a7fc2082ff4cb0afa2d63b2c32da1bef2bf6c5ca62dc6b72f9c",
            "0d8bb2d14e20cf9f6036152ed386d79189415b6d015a20133acb4e019139b94e9c146aaad5817f866c95d609a361735e",
            "052926add2207b76ca4fa57a8734416c8dc95e24501772c814278700eed6d1e4e8cf62d9c09db0fac349612b759e79a1",
            "08ba738453bfed09cb546dbb0783dbb3a5f1f566ed67bb6be0e8c67e2e81a4cc68ee29813bb7994998f3eae0c9c6a265",
        },
        .{
            "125435adce8e1cbd1c803e7123f45392dc6e326d292499c2c45c5865985fd74fe8f042ecdeeec5ecac80680d04317d80",
            "0e8828948c989126595ee30e4f7c931cbd6f4570735624fd25aef2fa41d3f79cfb4b4ee7b7e55a8ce013af2a5ba20bf2",
            "11def93719829ecda3b46aa8c31fc3ac9c34b428982b898369608e4f042babee6c77ab9218aad5c87ba785481eff8ae4",
            "0007c9cef122ccf2efd233d6eb9bfc680aa276652b0661f4f820a653cec1db7ff69899f8e52b8e92b025a12c822a6ce6",
            "03567bc5ef9c690c2ab2ecdf6a96ef1c139cc0b2f284dca0a9a7943388a49a3aee664ba5379a7655d3c68900be2f6903",
            "0b9c15f3fe6e5cf4211f346271d7b01c8f3b28be689c8429c85b67af215533311f0b8dfaaa154fa6b88176c229f2885d",
        },
    };
    inline for (vectors) |v| {
        const q0 = try testG1(v[0], v[1]);
        const q1 = try testG1(v[2], v[3]);
        const expected = try testG1(v[4], v[5]);
        const p = q0.add(q1).clearCofactor();
        try std.testing.expect(p.equivalent(expected));
        try std.testing.expect(p.isInSubgroup());
    }
}

test "RFC 9380 G2 cofactor clearing" {
    // Q0, Q1 and P = clear_cofactor(Q0 + Q1) for msg = "" and "abc", from RFC 9380, Appendix J.10.1.
    // Each coordinate is given as c0, then c1.
    const vectors = [_][12][]const u8{
        .{
            "019ad3fc9c72425a998d7ab1ea0e646a1f6093444fc6965f1cad5a3195a7b1e099c050d57f45e3fa191cc6d75ed7458c",
            "171c88b0b0efb5eb2b88913a9e74fe111a4f68867b59db252ce5868af4d1254bfab77ebde5d61cd1a86fb2fe4a5a1c1d",
            "0ba10604e62bdd9eeeb4156652066167b72c8d743b050fb4c1016c31b505129374f76e03fa127d6a156213576910fef3",
            "0eb22c7a543d3d376e9716a49b72e79a89c9bfe9feee8533ed931cbb5373dde1fbcd7411d8052e02693654f71e15410a",
            "113d2b9cd4bd98aee53470b27abc658d91b47a78a51584f3d4b950677cfb8a3e99c24222c406128c91296ef6b45608be",
            "13855912321c5cb793e9d1e88f6f8d342d49c0b0dbac613ee9e17e3c0b3c97dfbb5a49cc3fb45102fdbaf65e0efe2632",
            "0fd3def0b7574a1d801be44fde617162aa2e89da47f464317d9bb5abc3a7071763ce74180883ad7ad9a723a9afafcdca",
            "056f617902b3c0d0f78a9a8cbda43a26b65f602f8786540b9469b060db7b38417915b413ca65f875c130bebfaa59790c",
            "0141ebfbdca40eb85b87142e130ab689c673cf60f1a3e98d69335266f30d9b8d4ac44c1038e9dcdd5393faf5c41fb78a",
            "05cb8437535e20ecffaef7752baddf98034139c38452458baeefab379ba13dff5bf5dd71b72418717047f5b0f37da03d",
            "0503921d7f6a12805e72940b963c0cf3471c7b2a524950ca195d11062ee75ec076daf2d4bc358c4b190c0c98064fdd92",
            "12424ac32561493f3fe3c260708a12b7c620e7be00099a974e259ddc7d1f6395c3c811cdd19f1e8dbf3e9ecfdcbab8d6",
        },
        .{
            "12b2e525281b5f4d2276954e84ac4f42cf4e13b6ac4228624e17760faf94ce5706d53f0ca1952f1c5ef75239aeed55ad",
            "05d8a724db78e570e34100c0bc4a5fa84ad5839359b40398151f37cff5a51de945c563463c9efbdda569850ee5a53e77",
            "02eacdc556d0bdb5d18d22f23dcb086dd106cad713777c7e6407943edbe0b3d1efe391eedf11e977fac55f9b94f2489c",
            "04bbe48bfd5814648d0b9e30f0717b34015d45a861425fabc1ee06fdfce36384ae2c808185e693ae97dcde118f34de41",
            "19f18cc5ec0c2f055e47c802acc3b0e40c337256a208001dde14b25afced146f37ea3d3ce16834c78175b3ed61f3c537",
            "15b0dadc256a258b4c68ea43605dffa6d312eef215c19e6474b3e101d33b661dfee43b51abbf96fee68fc6043ac56a58",
            "05e47c1781286e61c7ade887512bd9c2cb9f640d3be9cf87ea0bad24bd0ebfe946497b48a581ab6c7d4ca74b5147287f",
            "19f98db2f4a1fcdf56a9ced7b320ea9deecf57c8e59236b0dc21f6ee7229aa9705ce9ac7fe7a31c72edca0d92370c096",
            "02c2d18e033b960562aae3cab37a27ce00d80ccd5ba4b7fe0e7a210245129dbec7780ccc7954725f4168aff2787776e6",
            "139cddbccdc5e91b9623efd38c49f81a6f83f175e80b06fc374de9eb4b41dfe4ca3a230ed250fbe3a2acf73a41177fd8",
            "1787327b68159716a37440985269cf584bcb1e621d3a7202be6ea05c4cfe244aeb197642555a0645fb87bf7466b2ba48",
            "00aa65dae3c8d732d10ecd2c50f8a1baf3001578f71c694e03866e9f3d49ac1e1ce70dd94a733534f106d4cec0eddd16",
        },
    };
    inline for (vectors) |v| {
        const q0 = G2{ .x = try testFp2(v[0], v[1]), .y = try testFp2(v[2], v[3]), .z = Fp2.one };
        const q1 = G2{ .x = try testFp2(v[4], v[5]), .y = try testFp2(v[6], v[7]), .z = Fp2.one };
        const expected = G2{ .x = try testFp2(v[8], v[9]), .y = try testFp2(v[10], v[11]), .z = Fp2.one };
        const p = q0.add(q1).clearCofactor();
        try std.testing.expect(p.equivalent(expected));
        try std.testing.expect(p.isInSubgroup());
    }
}
