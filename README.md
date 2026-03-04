# chia-vdf-verify

Pure-Rust Chia VDF (Verifiable Delay Function) proof verifier. No GMP, no C dependencies, no unsafe code.

This is a port of the verification path from [chiavdf](https://github.com/Chia-Network/chiavdf) (C++/GMP). Only proof **verification** is implemented — proof creation (proving) is not included.

## Why?

The existing `chiavdf` library depends on GMP (GNU Multiple Precision Arithmetic Library) via C/C++ linking, which is painful to build cross-platform — especially on Windows. This crate replaces GMP with [num-bigint](https://crates.io/crates/num-bigint) for a fully portable pure-Rust implementation.

## Performance

Verification is ~2x slower than the C/GMP version, which is acceptable for consensus validation (VDF verification is inherently fast by design — one proof per block).

| Variant | 512-bit disc, iters=100 | 1024-bit disc, depth=0 | 1024-bit disc, depth=2 |
|---------|------------------------|----------------------|----------------------|
| chiavdf C + GMP | ~1.1 ms | — | — |
| chia-vdf-verify (Rust) | ~2.3 ms | ~18 ms | ~57 ms |

Chia mainnet uses 1024-bit discriminants. Typical verification is under 60ms per proof.

## How VDF verification works

A [Verifiable Delay Function](https://en.wikipedia.org/wiki/Verifiable_delay_function) requires T sequential squarings to compute but is fast to verify. Chia uses the [Wesolowski scheme](https://eprint.iacr.org/2018/623) operating in [class groups](https://en.wikipedia.org/wiki/Ideal_class_group) of imaginary quadratic fields.

**Key concepts:**

- **Discriminant (D):** A large negative prime (1024 bits on mainnet) that defines the class group. Generated deterministically from a challenge hash via `CreateDiscriminant`.

- **Forms:** Elements of the class group, represented as binary quadratic forms (a, b, c) where b² − 4ac = D. These form a group under [composition](https://en.wikipedia.org/wiki/Binary_quadratic_form#Composition) (NUCOMP). The identity element is (1, 1, (1−D)/4).

- **The VDF computation:** Starting from the identity form x, compute y = x^(2^T) — i.e., square the form T times. This is inherently sequential.

- **The proof:** The prover also produces a proof form π. Verification checks:

  \[ \pi^B \cdot x^r = y \]

  where B = HashPrime(x ‖ y) is a 264-bit prime derived from the input/output, and r = 2^T mod B. This requires only a few group exponentiations — much faster than the T squarings.

- **Depth (n-Wesolowski):** A proof can be split into n segments, each with its own sub-proof. Depth 0 = single proof; higher depths break the proof into pieces with intermediate checkpoints. More segments means a larger proof blob but allows parallelized proving. Verification checks each segment in sequence.

## Usage

```rust
use chia_vdf_verify::discriminant::create_discriminant;
use chia_vdf_verify::verifier::check_proof_of_time_n_wesolowski;

let d = create_discriminant(seed, 1024);
let valid = check_proof_of_time_n_wesolowski(&d, &input_form, &proof_blob, iterations, depth);
```

## Testing

Run the standard test suite (fast, ~2 seconds):

```bash
cargo test
```

### Stress tests

The crate includes 110 test vectors extracted from chiavdf's `vdf.txt` — real VDF proofs with 1024-bit discriminants at depths 0 through 7. Two vectors run by default as a smoke test. To run all 110:

```bash
cargo test --release -- --ignored test_vdf_txt_all
```

This takes ~15 seconds in release mode.

### Benchmarks

Compare verification performance using Criterion:

```bash
cargo bench --bench verify
```

## Architecture

Ported from chiavdf's C++ verification path (~2,700 LOC):

| Module | Source | Purpose |
|--------|--------|---------|
| `verifier` | `verifier.h` | `VerifyWesolowskiProof`, `CheckProofOfTimeNWesolowski` |
| `proof_common` | `proof_common.h` | `FastPow`, `FastPowFormNucomp`, `GetB`, serialization |
| `nucomp` | `nucomp.h` | Class group form composition (`nucomp`, `nudupl`) |
| `reducer` | `Reducer.h` | Pulmark form reduction |
| `xgcd_partial` | `xgcd_partial.c` | Partial extended GCD (Lehmer-accelerated) |
| `bqfc` | `bqfc.c` | Compressed form serialization (BQFC format) |
| `primetest` | `primetest.h` | BPSW primality test, `HashPrime` |
| `discriminant` | `create_discriminant.h` | Discriminant generation from seed |
| `form` | `ClassGroup.h` | Quadratic form (a, b, c) with discriminant |
| `integer` | `integer_common.h` | BigInt wrapper, Lehmer extended GCD |

## License

[Apache License 2.0](LICENSE)
