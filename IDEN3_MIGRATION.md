# Poseidon377 → iden3/circom Poseidon Migration

This document describes the changes made to the Penumbra `poseidon377` crate to produce
hash outputs identical to the **iden3 Poseidon** implementation used in circom/Solidity
projects.

## Background

The original `poseidon377` was built for Penumbra's **decaf377** curve (BLS12-377 scalar
field). It used Penumbra's own parameter generation pipeline:

- **MDS matrix**: Cauchy construction
- **Round constants**: Merlin transcript–based derivation
- **Partial round counts**: Penumbra's own optimality criteria

The iden3 implementation targets the **BN254 scalar field** and uses a different parameter
generation pipeline:

- **MDS matrix**: Pre-computed constants distributed in `poseidon_constants.json`
- **Round constants**: Grain LFSR derivation
- **Partial round counts**: iden3's own rounding rules (generally more conservative /
  higher than Penumbra's)

Both implementations share the same core Poseidon algorithm — the same sponge
construction, the same S-box (α = 5), and the same optimized permutation structure with
full and partial rounds. Only the concrete parameters and the target field differ.

### Why iden3's constants are the "official" ones

The iden3 constants in `poseidon_constants.json` are pre-computed outputs of the official
`generate_parameters_grain.sage` script, written by the original Poseidon paper authors
(IAIK TU Graz, published in the
[hadeshash repo](https://extgit.iaik.tugraz.at/krypto/hadeshash)). The script takes the
field, S-box, width, round counts, and field modulus as inputs and deterministically
produces the round constants and MDS matrix via Grain LFSR. For the BN254 scalar field
with t=3, that invocation is:

```bash
sage generate_parameters_grain.sage 1 0 254 3 8 57 \
  0x30644e72e131a029b85045b68181585d2833e84879b9709143e1f593f0000001
```

iden3 ran this script for each supported width and bundled the results into
`poseidon_constants.json`. Since the constants come directly from the paper authors'
reference implementation, they are the canonical/official parameter set for Poseidon over
BN254. Penumbra's deviation — using Cauchy MDS and Merlin-derived round constants — is
what made their outputs incompatible.

## Changes Made

### 1. Field: decaf377 → BN254

Replaced the `decaf377` dependency with `cycles-curve-bn254` across all workspace crates.
The `cycles-curve-bn254` crate exposes the same `Fq` type and arkworks traits, so the
permutation code itself required no algorithmic changes — only import and Cargo.toml
updates.

Affected crates: `poseidon377`, `poseidon-permutation`, `poseidon-parameters`,
`poseidon-consistency`, `poseidon-paramgen`, `poseidon-tests`.

### 2. Parameters: Penumbra → iden3/Solana

Regenerated all hardcoded parameter files (`poseidon377/src/params/rate_*.rs`) using
iden3's pre-computed constants. Penumbra's generic parameter generator (`poseidon-paramgen`
crate) was not suitable for this because it derives its own MDS matrix and round constants
via Cauchy construction / Merlin transcripts — we need the *exact* Grain LFSR constants to
get matching outputs. A custom generator was written instead.

#### The generator: `poseidon377/examples/gen_params.rs`

This is a self-contained Rust binary (~840 lines) that reads iden3's
`poseidon_constants.json` and emits a complete `rate_*.rs` file to stdout. It:

1. **Loads raw constants** — parses the JSON file's `C[t-2]` array (round constants) and
   `M[t-2]` array (MDS matrix) for the requested width `t`, converting each hex string to
   a `cycles_curve_bn254::Fq` element.

2. **Computes optimized parameters** — the Poseidon permutation has an "optimized"
   variant (Appendix B of the Poseidon paper) that restructures the partial rounds for
   efficiency. The generator computes:
   - **Optimized ARC matrix**: round constants are transformed so that during partial
     rounds, only one S-box and one matrix-vector multiply are needed per round.
   - **Equivalent matrices**: the MDS matrix is decomposed into `M' * M''` form, plus
     per-round sparse matrices (`v_collection`, `w_hat_collection`) that enable the
     optimized partial-round computation.
   - **M_i**: the "initial" matrix applied before entering the partial-round loop.

3. **Verifies correctness** — runs the full permutation in both unoptimized and optimized
   modes on a fixed test vector and checks that both produce the same output. This appears
   on stderr during generation.

4. **Emits Rust source** — outputs a `const fn rate_N()` that returns a
   `PoseidonParameters<...>` struct populated with all matrices and constants in
   Montgomery-limb form (`Fq::from_montgomery_limbs([...])`) so they can be evaluated at
   compile time.

The generator includes a small dynamic-matrix library over `Fq` (multiplication, inverse,
transpose, sub-matrices, etc.) to avoid const-generic complexity.

#### Usage

```bash
# Build (--no-default-features avoids compiling the r1cs feature)
cargo build --example gen_params -p poseidon377 --no-default-features

# Generate parameters for width t (e.g. t=3 → rate_2.rs)
./target/debug/examples/gen_params 3 /path/to/poseidon_constants.json \
    > poseidon377/src/params/rate_2.rs
```

After generating a new `rate_*.rs`, the corresponding const-generic type signature in
`poseidon377/src/lib.rs` must also be updated to match the new round counts. For example,
`t=3` with R_P=57 gives `PoseidonParameters<3, 2, 9, 4, 65, 3, 195, 57>`.

If you need to regenerate `rate_7.rs` (which requires a large lib to compile), temporarily
comment out `RATE_7_PARAMS` in `lib.rs` before building `gen_params`, generate the file,
then restore the constant with the updated type signature.

#### Source of constants

The `poseidon_constants.json` file comes from iden3's `circomlibjs` package
(`circomlibjs/src/poseidon_constants.json`). It contains pre-computed round constants and
MDS matrices for widths t=2 through t=17, produced by the official
`generate_parameters_grain.sage` script as described above.

#### Parameter summary

| hash_n | Width (t) | Rate | R_F (full) | R_P (partial) | Total rounds | Round constants |
|--------|-----------|------|------------|---------------|--------------|-----------------|
| hash_1 | 2         | 1    | 8          | 56            | 64           | 128             |
| hash_2 | 3         | 2    | 8          | 57            | 65           | 195             |
| hash_3 | 4         | 3    | 8          | 56            | 64           | 256             |
| hash_4 | 5         | 4    | 8          | 60            | 68           | 340             |
| hash_5 | 6         | 5    | 8          | 60            | 68           | 408             |
| hash_6 | 7         | 6    | 8          | 63            | 71           | 497             |
| hash_7 | 8         | 7    | 8          | 64            | 72           | 576             |

All parameter sets use security level 128 and α = 5. All seven hash functions are fully
implemented.

### 3. Output word: state[1] → state[0]

Penumbra's sponge returned `state_words[1]` as the hash output. iden3 returns
`state_words[0]` (the capacity element). Changed the output index to `[0]` in:

- `poseidon-permutation/src/permutation.rs` — `n_to_1_fixed_hash()` and
  `unoptimized_n_to_1_fixed_hash()`
- `poseidon-permutation/src/r1cs.rs` — the R1CS (in-circuit) equivalent

### 4. Domain separator convention

iden3 does not use a domain separator — it leaves `state[0]` as zero and places inputs in
`state[1..]`. This implementation keeps an explicit `domain_separator` argument for
flexibility. To produce outputs matching iden3, pass `Fq::from(0)` as the domain
separator.

### 5. Old test vectors removed

The original Penumbra test vectors in `poseidon377/src/hash.rs` are commented out since
they are specific to the decaf377 parameterization. A test for `hash_2` against the
canonical iden3 test vector is included in `poseidon377/src/hash.rs`.

## Validation

### Canonical test vector (hash_2)

- **Inputs** (big-endian bytes): `[1u8; 32]` and `[2u8; 32]`, zero domain separator
- **Expected** (big-endian): `[13, 84, 225, 147, 143, 138, 140, 28, 125, 235, 94, 3, 85, 242, 99, 25, 32, 123, 132, 254, 156, 162, 206, 27, 38, 231, 53, 200, 41, 130, 25, 144]`

This test passes in `poseidon377/src/hash.rs` as `hash_2_solana_test_vector`.

### Comparison harness

A comparison harness exists in the sibling `poseidon-bench` repo:

- `poseidon377/examples/hash_cli.rs` — Rust CLI that hashes inputs and prints the result
  as a decimal string
- `poseidon-bench/circomlibjs/tools/hash_cli.mjs` — JS CLI using iden3's `circomlibjs`
- `poseidon-bench/scripts/compare_simple.mjs` — runs both CLIs with the same inputs and
  reports matches/mismatches

The `hash_2` (t=3) configuration has been verified to produce identical outputs across
both implementations.

## What's left

- R1CS circuit tests (`poseidon-tests/tests/poseidon377_r1cs.rs`) reference old Penumbra
  test vectors and need updating.
- The `poseidon-paramgen` crate (Penumbra's generic parameter generator) is unused in this
  fork — the custom `gen_params.rs` example replaces it for iden3/Solana-compatible
  generation.
