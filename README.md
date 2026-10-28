# Pairings for Zig

A pairing-friendly elliptic curve library for Zig. It does BLS12-381 and BN462.

You'd use it for BLS signatures, zero-knowledge proofs, or anything else that needs pairings.

It follows revision 14 of the [`draft-irtf-cfrg-pairing-friendly-curves`](https://datatracker.ietf.org/doc/draft-irtf-cfrg-pairing-friendly-curves/14/) spec, which recommends these two curves at the 128-bit security level.

## Curves

### BLS12-381

The main curve. Ethereum, Zcash, and a bunch of other systems already use it.

Its security level is about 126 bits, not 128.

- G1: points on E(Fp): y^2 = x^3 + 4
- G2: points on the sextic twist E'(Fp2): y^2 = x^3 + 4(u+1)
- GT: target group in Fp12, reached via the optimal Ate pairing

Points use the same encoding as Zcash and most other libraries, as described in Section 5 of the draft.

### BN462

A Barreto-Naehrig curve with about 134 bits of security.

It's a lot slower than BLS12-381, but gives you more margin against future attacks.

The draft doesn't define a point encoding for BN462, so this library made one up: three flags from BLS12-381 packed into a leading byte, then the big-endian coordinates. Appendix E.4 of the draft describes this format. It doesn't interoperate with the SEC1 encoding used by MIRACL and pfcurve-js.

## Usage

Add the package to your `build.zig.zon` dependencies, then:

```zig
const bls = @import("pairings").bls12_381;

const g1 = bls.G1.basePoint;
const g2 = bls.G2.basePoint;

// Compute a pairing
const gt = bls.pairing.pair(g1, g2);

// Scalar multiplication, in constant time
const s: bls.scalar.CompressedScalar = ...; // 32-byte big-endian scalar
const p = g1.mul(s, .big);

// Random points (requires std.Io)
const random_g1 = bls.G1.random(io);
const random_g2 = bls.G2.random(io);

// Serialization
const compressed = g1.toCompressed();
const restored = try bls.G1.fromCompressed(compressed);

// Decoding a point in either form
const decoded = try bls.G1.fromBytes(received_bytes);
```

## Validation

Decoding a point — from bytes or from coordinates — always checks that it's in the right group.

Group operations keep points in their group. The only way to get a point outside of it is to fill in the struct fields by hand, and `isInSubgroup()` is there for that case.

The draft leaves a few choices to the protocol using it, and so does this library:

- `fromCompressed()` and `fromUncompressed()` each accept one form, and `fromBytes()` accepts both. The draft recommends the compressed form for data sent over the wire.
- Decoding accepts the identity element. If your protocol doesn't want it, call `rejectIdentity()` on the result.
- `scalar.fromBytes()` accepts zero. If your protocol doesn't want it, check with `scalar.isZero()`.

Use `mul()` for secret scalars — it runs in constant time. `mulPublic()` is faster, but only for scalars that are not secret.
