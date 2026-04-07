//! Shared utility functions for cryptographic operations.
//!
//! Provides:
//! - Field element conversion from BigUint
//! - Modular inverse computation
//! - Cryptographically secure random scalar generation
//! - Scalar-to-curve-point mapping

use ark_ec::Group;
use ark_ff::{Field, PrimeField};
use ark_secp256k1::{Fr, Projective};
use num_bigint::BigUint;
use num_traits::Num;
use std::time::{SystemTime, UNIX_EPOCH};
use getrandom::getrandom;
use sha2::{Sha256, Digest};


pub fn biguint_to_fr<F: PrimeField>(value: &BigUint) -> Result<F, String> {
    // Reduce modulo field modulus first
    let modulus_biguint = BigUint::from_str_radix(&F::MODULUS.to_string(), 10)
        .map_err(|_| "Failed to convert modulus to BigUint".to_string())?;
    let reduced_value = value % &modulus_biguint;
    
    // Try converting to u64 first for small values (fits in 8 bytes)
    let bytes = reduced_value.to_bytes_be();
    if bytes.len() <= 8 {
        let mut u64_bytes = [0u8; 8];
        let start_idx = 8 - bytes.len();
        u64_bytes[start_idx..].copy_from_slice(&bytes);
        let u64_val = u64::from_be_bytes(u64_bytes);
        return Ok(F::from(u64_val));
    }

    let bytes_le = reduced_value.to_bytes_le();
    Ok(F::from_le_bytes_mod_order(&bytes_le))
}

pub fn mod_inverse(value: &Fr) -> Fr{
    value.inverse().unwrap()
}

pub fn random_scalar() -> Fr {
    let mut random_bytes = [0u8; 32];
    getrandom(&mut random_bytes).expect("Failed to generate random bytes");

    let dust = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_secs();
    let dust_bytes = dust.to_be_bytes();

    let mut combined = Vec::new();
    combined.extend_from_slice(&random_bytes);
    combined.extend_from_slice(&dust_bytes);

    let hash = Sha256::digest(&combined);
    let mut key_bytes = [0u8; 32];
    key_bytes.copy_from_slice(&hash);
    Fr::from_random_bytes(&key_bytes).unwrap()
}

pub fn point_from_scalar(s: &Fr) -> Projective {
    let g = Projective::generator();
    g * s
}