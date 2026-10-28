# Pairings for Zig

A pairing-friendly elliptic curve library for Zig, implementing BLS12-381 and BN462.

Useful for BLS signatures, zero-knowledge proofs, etc.

This follows revision 14 of the [`draft-irtf-cfrg-pairing-friendly-curves`](https://datatracker.ietf.org/doc/draft-irtf-cfrg-pairing-friendly-curves/14/) draft, which recommends these two curves at the 128-bit security level.

## Curves

### BLS12-381

The primary curve, widely used in Ethereum, Zcash, and other systems.

Its security level is expected to be about 126 bits rather than 128.

- G1: points on E(Fp): y^2 = x^3 + 4
- G2: points on the sextic twist E'(Fp2): y^2 = x^3 + 4(u+1)
- GT: target group in Fp12, reached via the optimal Ate pairing

Points use the same encoding as Zcash and most other libraries, as described in Section 5 of the draft.

### BN462

A Barreto-Naehrig curve with a security level of about 134 bits.

It is a lot slower than BLS12-381, but has more margin against future attacks.

The draft doesn't define a point encoding for BN462.

This library puts the same three flags as BLS12-381 in a leading byte, followed by the big-endian coordinates.
Appendix E.4 of the draft describes this format.
It doesn't interoperate with the SEC1 encoding used by MIRACL and pfcurve-js.

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

Decoding a point, from bytes or from coordinates, always checks that it is in the right group.

Group operations keep points in their group.
The only way to get a point outside of it is to fill in the struct fields by hand, and `isInSubgroup()` is there for that case.

The draft leaves a few choices to the protocol using it, and so does this library:

- `fromCompressed()` and `fromUncompressed()` accept one form each, and `fromBytes()` accepts both.
  The draft recommends the compressed form for data sent over the wire.
- Decoding accepts the identity element.
  Protocols that don't want it can call `rejectIdentity()` on the result.
- `scalar.fromBytes()` accepts zero.
  Protocols that don't want it can check with `scalar.isZero()`.

Use `mul()` for secret scalars, as it runs in constant time.
`mulPublic()` is faster, but only for scalars that are not secret.
