//! Shared utility functions for cryptographic operations.
//!
//! Provides:
//! - Field element conversion from `BigUint`
//! - Modular inverse computation
//! - Cryptographically secure random scalar generation
//! - Scalar-to-curve-point mapping

use ark_ec::Group;
use ark_ff::{Field, PrimeField};
use ark_secp256k1::{Fr, Projective};
use getrandom::getrandom;
use num_bigint::BigUint;

pub fn biguint_to_fr<F: PrimeField>(value: &BigUint) -> Result<F, String> {
    Ok(F::from_le_bytes_mod_order(&value.to_bytes_le()))
}

#[must_use]
pub fn mod_inverse(value: &Fr) -> Fr {
    value.inverse().unwrap()
}

#[must_use]
pub fn random_scalar() -> Fr {
    let mut random_bytes = [0u8; 32];
    getrandom(&mut random_bytes).expect("Failed to generate random bytes");
    Fr::from_le_bytes_mod_order(&random_bytes)
}

pub fn point_from_scalar(s: &Fr) -> Projective {
    let g = Projective::generator();
    g * s
}
