//! ECDSA signature scheme over secp256k1.
//!
//! Implements key generation, signing, and verification
//! using the standard ECDSA algorithm with SHA-256 message hashing.
//!
//! See: Theory document, Section 3 (Digital Signatures)

use ark_ec::{CurveGroup, Group};
use ark_ff::{Field, PrimeField};
use ark_secp256k1::{Fr, Projective};
use ark_std::Zero;
use sha2::{Sha256, Digest};
use crate::utils::*; 

pub type PrivateKey = Fr;
pub type PublicKey = Projective;

#[derive(Debug, Clone, PartialEq)]
pub struct Signature {
    pub r: Fr,
    pub s: Fr,
    pub v: bool,
}

impl Signature {
    pub fn new(r: Fr, s: Fr, v: bool) -> Self {  
        if r.is_zero() {
            panic!("ECDSA signature r component cannot be zero: r is the x-coordinate of R = k*G and must be non-zero to ensure the private key term (r*private_key) is present in the signature equation s = k^(-1) * (h + r*private_key)");
        }
        if s.is_zero() {
            panic!("ECDSA signature s component cannot be zero: s must be invertible as verification requires computing s^(-1)");
        }
        Self { r, s, v }
    }
}

#[derive(Debug, Clone)]
pub struct KeyPair {
    pub private_key: PrivateKey,
    pub public_key: PublicKey,
}

impl KeyPair {
    pub fn generate() -> Self {
        let private_key = random_scalar();
        Self::from_private_key(private_key)
    }

    pub fn from_private_key(private_key: PrivateKey) -> Self {
        let public_key = point_from_scalar(&private_key);
        Self {
            private_key,
            public_key,
        }
    }
}

pub fn sign(private_key: &PrivateKey, message_hash: &Fr, v: bool) -> Signature {
    loop {
        // ECDSA signing nonce k: must be unique per (private_key_i, message_i) pair.
        // Could be deterministically derived by hashing (private_key + message_hash).
        // This is the signing nonce, separate from transaction nonces (which increments sequentially).
        // Present implementation using random_scalar() could theoretically reuse the same k for different messages from the same private_key, though the probability is negligible (approximately 1/2^256 per generation).
        let k = random_scalar(); //Todo: create hash function as outlined above
        
        let g = Projective::generator();
        let r_point = g * k;
        let r_affine = r_point.into_affine();
        let r_bigint = if v {
            r_affine.x.into_bigint()
        } else {
            r_affine.y.into_bigint()
        };

        let r: Fr = Fr::from_bigint(r_bigint).unwrap_or_default();

        if r.is_zero() {
            continue;
        }

        let s = mod_inverse(&k) * (message_hash + &(r * private_key));

        if s.is_zero() {
            continue;
        }

        return Signature::new(r ,s, v);
    }
}

pub fn verify(signature: &Signature, message: &[u8], public_key: &Projective) -> bool {
    let h = hash_message(message);
    let s_1 = mod_inverse(&signature.s);
    let g = Projective::generator();
    let r_1 = g * (h * s_1) + *public_key * (signature.r * s_1);
    let r_affine = r_1.into_affine();
    let r_bigint = if signature.v {
        r_affine.x.into_bigint()
    } else {
        r_affine.y.into_bigint()
    };
    let r_1_fr: Fr = Fr::from_bigint(r_bigint).unwrap_or_default();

    signature.r == r_1_fr
}

pub fn hash_message(message: &[u8]) -> Fr {
    let hash = Sha256::digest(message);
    Fr::from_random_bytes(&hash).unwrap()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sign_and_verify_with_x_coordinate() {
        let key_pair = KeyPair::generate();
        let message = b"Test message";
        let h = hash_message(message);
        let v = true;

        let signature = sign(&key_pair.private_key, &h, v);
        assert!(verify(&signature, message, &key_pair.public_key));
    }

    #[test]
    fn test_sign_and_verify_with_y_coordinate() {
        let key_pair = KeyPair::generate();
        let message = b"Test message";
        let h = hash_message(message);
        let v = false;

        let signature = sign(&key_pair.private_key, &h, v);
        assert!(verify(&signature, message, &key_pair.public_key));
    }

    #[test]
    fn test_sign_and_verify_with_wrong_message() {
        let key_pair = KeyPair::generate();
        let message1 = b"Test message 1";
        let h1 = hash_message(message1);
        let signature = sign(&key_pair.private_key, &h1, true);
        
        let message2 = b"Test message 2";
        assert!(!verify(&signature, message2, &key_pair.public_key));
    }
}