//! Pedagogical Rust implementations of the cryptographic primitives behind
//! Groth16: elliptic-curve signatures (ECDSA over secp256k1), bilinear
//! pairings on BN254, R1CS and QAP constructions, trusted-setup evaluation,
//! and the complete Groth16 prover and verifier.
//!
//! Each module mirrors a section of the companion theory document
//! (`docs/theory.pdf`); constructions favour explicit, math-shaped code over
//! library shortcuts, and intermediate schemes deliberately retain documented
//! weaknesses that the later modules close.

pub mod bilinear_pairings;
pub mod ecdsa;
pub mod ecp_matrix;
pub mod groth16;
pub mod qap_trusted_setup;
pub mod r1cs;
pub mod r1cs_to_qap;
pub mod rational_sums;
pub mod utils;
