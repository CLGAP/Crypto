# zk-crypto

Rust implementations of cryptographic primitives building toward a full Groth16 zk-SNARK, from ECDSA and elliptic curve operations through R1CS, QAP, trusted setup, and the complete Groth16 protocol.

## Resources

RareSkills(J. Scholz), Abstract Algebra (D. Dummit and R. Foote), and Wikipedia.

For detailed mathematical theory and proofs, see [docs/theory.pdf](docs/theory.pdf).

## Modules

| Module | Description |
|---|---|
| `rational_sums` | Verify rational sum equations on elliptic curves via scalar-to-point mapping |
| `ecp_matrix` | Matrix-vector multiplication on elliptic curve points; homomorphism verification |
| `ecdsa` | ECDSA key generation, signing, and verification over secp256k1 |
| `bilinear_pairings` | Multi-pairing equation verification on BN254 |
| `r1cs` | R1CS verification: field-based (Hadamard) and pairing-based |
| `r1cs_to_qap` | Lagrange interpolation, QAP construction, polynomial quotient computation |
| `qap_trusted_setup` | SRS generation and Pinocchio-style prove/verify (prover sees tau) |
| `groth16::alpha_beta` | Groth16 with alpha/beta shifting (no ZK blinding) |
| `groth16::complete` | Full Groth16: gamma/delta separation, r/s blinding, public input binding |
| `utils` | Shared helpers: field conversion, modular inverse, random scalar generation |

## Usage

Run all demos:

```bash
cargo run --example demo
```

Run a specific module demo:

```bash
cargo run --example ecdsa
cargo run --example rational_sums
cargo run --example ecp_matrix
cargo run --example bilinear_pairings
cargo run --example r1cs
cargo run --example r1cs_to_qap
cargo run --example qap_trusted_setup
cargo run --example groth16_alpha_beta
cargo run --example groth16_complete
```

## Tests

```bash
cargo test
```

Proptest example:

```bash
cargo test demo_proptest -- --ignored --nocapture  
```

## Dependencies

- `ark-ec`, `ark-ff`, `ark-std`, `ark-poly` -- Arkworks core algebra
- `ark-secp256k1` -- secp256k1 curve (ECDSA, rational sums, ECP matrix)
- `ark-bn254` -- BN254 curve and bilinear pairings
- `sha2` -- SHA-256 message hashing
- `getrandom` -- OS-entropy random number generation
- `num-bigint`, `num-traits` -- Arbitrary precision integers

### Dev Dependencies

- `proptest` -- Property-based testing

## Building

```bash
cargo build
```

## Documentation

Mathematical theory and proofs are in [docs/theory.pdf](docs/theory.pdf).