# Pairings for Zig

A pairing-friendly elliptic curve library for Zig, implementing BLS12-381 and BN462 at the 128-bit security level.

Useful for BLS signatures, zero-knowledge proofs, etc.

This work is based on the [`draft-irtf-cfrg-pairing-friendly-curves`](https://datatracker.ietf.org/doc/draft-irtf-cfrg-pairing-friendly-curves/) draft.

## Curves

### BLS12-381

The primary curve, widely used in Ethereum 2.0, Zcash, and other systems.

- G1: Points on E(Fp): y^2 = x^3 + 4 (projective coordinates)
- G2: Points on the sextic twist E'(Fp2): y^2 = x^3 + 4(u+1) (projective coordinates)
- GT: Target group in Fp12, reached via optimal Ate pairing

### BN462

A Barreto-Naehrig curve at the 128-bit security level.

## Usage

Add the package to your `build.zig.zon` dependencies, then:

```zig
const bls = @import("pairings").bls12_381;

const g1 = bls.G1.basePoint;
const g2 = bls.G2.basePoint;

// Compute a pairing
const gt = bls.pairing.pair(g1, g2);

// Scalar multiplication
const s: bls.scalar.CompressedScalar = ...; // 32-byte little-endian scalar
const p = g1.mul(s, .little) catch unreachable;

// Random points (requires std.Io)
const random_g1 = bls.G1.random(io);
const random_g2 = bls.G2.random(io);

// Serialization
const compressed = g1.toCompressed();
const restored = try bls.G1.fromCompressed(compressed);
```
