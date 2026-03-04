# chia-vdf-verify

Pure-Rust Chia VDF (Verifiable Delay Function) proof verifier. No GMP, no C dependencies, no unsafe code.

This is a port of the verification path from [chiavdf](https://github.com/Chia-Network/chiavdf) (C++/GMP). Only proof **verification** is implemented — proof creation (proving) is not included.

## Why?

The existing `chiavdf` library depends on GMP (GNU Multiple Precision Arithmetic Library) via C/C++ linking, which is painful to build cross-platform — especially on Windows. This crate replaces GMP with [num-bigint](https://crates.io/crates/num-bigint) for a fully portable pure-Rust implementation.

## Performance

Verification is ~2x slower than the C/GMP version, which is acceptable for consensus validation (VDF verification is inherently fast by design — one proof per block, typically under 20ms for 1024-bit discriminants).

| Variant | Time (iters=100, 512-bit) |
|---------|--------------------------|
| chiavdf C + GMP | ~1.1 ms |
| chia-vdf-verify (Rust) | ~2.3 ms |

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
