# Benchmark Summary: Rust chia-vdf-verify vs C++ chiavdf

## Result

**C++ was 1.75x faster → now 1.12x faster.** The Rust VDF verifier closes
most of the gap with C++ chiavdf (GMP-backed) through pure Rust optimizations
using crates.io malachite-nz, with no C dependencies.

## Baseline (before optimization)

| | C++ (chiavdf + GMP) | Rust (chia-vdf-verify + malachite-nz) | Ratio |
|---|---|---|---|
| Single-threaded ms/proof | 8.67 ms | 15.19 ms | **C++ 1.75x faster** |

## After optimization

### Real-world: 50 mainnet end-of-slot proofs, single-threaded

| Run | C++ (ms/proof) | Rust (ms/proof) | Ratio |
|---|---|---|---|
| 1 | 11.78 | 13.18 | C++ 1.12x faster |
| 2 | 11.84 | 13.25 | C++ 1.12x faster |
| 3 | 11.99 | 13.33 | C++ 1.11x faster |
| 4 | 11.81 | 13.18 | C++ 1.12x faster |
| 5 | 11.76 | 13.21 | C++ 1.12x faster |

(Absolute times vary with system load; the ratio is stable at **~1.12x**.)

### Core operation: nudupl + reduce (1024-bit discriminant)

| Version | Time per op | Speedup |
|---|---|---|
| Before (main branch) | 51 µs | — |
| After (malachite branch) | 19 µs | **2.7x faster** |

## Optimizations applied

| # | Optimization | Impact |
|---|---|---|
| 1 | Port from num-bigint to malachite-nz | foundation |
| 2 | Fix PyO3 bindings, release GIL | correctness + parallel |
| 3 | Discriminant bytes API (avoid repeated decimal parse) | small |
| 4 | Extract limb words without allocation in Lehmer loop | ~5% |
| 5 | Eliminate clones, use in-place negation (`NegAssign`) | ~15% |
| 6 | O(n) byte-to-integer via direct limb construction | ~30% on decompression |
| 7 | Fused multiply-accumulate (`AddMulAssign`/`SubMulAssign`) | ~10% |
| 8 | Owned-argument GCD variants, avoid double-clones | small |
| 9 | Compiler: LTO=fat, codegen-units=1, target-cpu=native | ~5-10% |
| 10 | Optimize `fdiv_r`, BQFC decompression, refactor nudupl | small |
| 11 | GCD argument swap: return native Bézout cofactor | small |

## What's left

The remaining ~12% gap comes from GMP vs malachite-nz fundamentals:

- **GMP reuses pre-allocated buffers** via thread-local scratch (`mpz_t`);
  malachite allocates fresh `Vec`s per operation.
- **GMP has hand-tuned x86-64 assembly** for core limb operations
  (`mpn_addmul_1`, `mpn_mul_basecase`); malachite uses LLVM codegen.
- **GMP's `extended_gcd` skips unused cofactors**; malachite always derives
  both Bézout coefficients (the second via a full multiply + divide).

A malachite fork adding `extended_gcd_first_cofactor` (skip second-cofactor
derivation) and `Integer × i64` (avoid wrapper allocation) closes another
~14% on nudupl in testing, bringing Rust to approximate parity with C++.

## Environment

- 1024-bit class group discriminants, ~130M iterations per proof
- Rust stable, release profile with LTO
- malachite-nz 0.9.1 (crates.io, no fork)
- GMP: system package (used by chiavdf)
