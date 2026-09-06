//! ECDSA signature scheme over secp256k1.
//!
//! Implements key generation, signing, and verification
//! using the standard ECDSA algorithm with SHA-256 message hashing.
//!
//! See: Theory document, Section 3 (Digital Signatures)

use ark_ec::{CurveGroup, Group};
use ark_ff::{BigInteger, PrimeField};
use ark_secp256k1::{Fr, Projective};
use ark_std::Zero;
use sha2::{Sha256, Digest};
use crate::utils::{random_scalar, point_from_scalar, mod_inverse}; 

pub type PrivateKey = Fr;
pub type PublicKey = Projective;

// Ethereum-standard signature: r is always x(R) mod n, and v is the recovery
// bit recording the parity of R.y, distinguishing the two points +-R that
// share the x-coordinate r; it enables public-key recovery from (r, s, v).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Signature {
    pub r: Fr,
    pub s: Fr,
    pub v: bool,
}

impl Signature {
    #[must_use]
    pub fn new(r: Fr, s: Fr, v: bool) -> Self {  
        assert!(!r.is_zero(), "ECDSA signature r component cannot be zero: r is the x-coordinate of R = k*G and must be non-zero to ensure the private key term (r*private_key) is present in the signature equation s = k^(-1) * (h + r*private_key)");
        assert!(!s.is_zero(), "ECDSA signature s component cannot be zero: s must be invertible as verification requires computing s^(-1)");
        Self { r, s, v }
    }
}

#[derive(Debug, Clone)]
pub struct KeyPair {
    pub private_key: PrivateKey,
    pub public_key: PublicKey,
}

impl KeyPair {
    #[must_use]
    pub fn generate() -> Self {
        let private_key = random_scalar();
        Self::from_private_key(private_key)
    }

    #[must_use]
    pub fn from_private_key(private_key: PrivateKey) -> Self {
        let public_key = point_from_scalar(&private_key);
        Self {
            private_key,
            public_key,
        }
    }
}

#[must_use]
pub fn sign(private_key: &PrivateKey, message_hash: &Fr) -> Signature {
    loop {
        // ECDSA signing nonce k: must be unique per (private_key_i, message_i) pair.
        // Could be deterministically derived by hashing (private_key + message_hash).
        // This is the signing nonce, separate from transaction nonces (which increments sequentially).
        // Present implementation using random_scalar() could theoretically reuse the same k for different messages from the same private_key, though the probability is negligible (approximately 1/2^256 per generation).
        let k = random_scalar(); //Todo: create hash function as outlined above

        let g = Projective::generator();
        let r_point = g * k;
        let r_affine = r_point.into_affine();
        // Ethereum standard: r = x(R) reduced mod n; v = parity of R.y
        let r: Fr = Fr::from_le_bytes_mod_order(&r_affine.x.into_bigint().to_bytes_le());
        let mut v = r_affine.y.into_bigint().is_odd();

        if r.is_zero() {
            continue;
        }

        let mut s = mod_inverse(&k) * (message_hash + &(r * private_key));

        if s.is_zero() {
            continue;
        }

        // EIP-2 low-s: canonicalize to s <= n/2. Negating s corresponds to the
        // nonce -k, i.e., R -> -R, so the recovery parity flips with it.
        if s.into_bigint() > Fr::MODULUS_MINUS_ONE_DIV_TWO {
            s = -s;
            v = !v;
        }

        return Signature::new(r ,s, v);
    }
}

#[must_use]
pub fn verify(signature: &Signature, message: &[u8], public_key: &Projective) -> bool {
    if signature.r.is_zero() || signature.s.is_zero() || public_key.is_zero() {
        return false;
    }
    let h = hash_message(message);
    let s_1 = mod_inverse(&signature.s);
    let g = Projective::generator();
    let r_1 = g * (h * s_1) + *public_key * (signature.r * s_1);
    let r_affine = r_1.into_affine();
    // r = x(R') mod n; v is not needed for verification, only for recovery
    let r_1_fr: Fr = Fr::from_le_bytes_mod_order(&r_affine.x.into_bigint().to_bytes_le());

    signature.r == r_1_fr
}

#[must_use]
pub fn hash_message(message: &[u8]) -> Fr {
    let hash = Sha256::digest(message);
    Fr::from_be_bytes_mod_order(&hash)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sign_and_verify() {
        let key_pair = KeyPair::generate();
        let message = b"Test message";
        let h = hash_message(message);

        let signature = sign(&key_pair.private_key, &h);
        assert!(verify(&signature, message, &key_pair.public_key));
    }

    #[test]
    fn test_sign_and_verify_with_wrong_message() {
        let key_pair = KeyPair::generate();
        let message1 = b"Test message 1";
        let h1 = hash_message(message1);
        let signature = sign(&key_pair.private_key, &h1);

        let message2 = b"Test message 2";
        assert!(!verify(&signature, message2, &key_pair.public_key));
    }

    #[test]
    fn test_verify_rejects_wrong_key() {
        let key_pair = KeyPair::generate();
        let other = KeyPair::generate();
        let message = b"Test message";
        let h = hash_message(message);
        let signature = sign(&key_pair.private_key, &h);
        assert!(!verify(&signature, message, &other.public_key));
    }

    #[test]
    fn test_verify_rejects_forged_zero_components() {
        let key_pair = KeyPair::generate();
        let message = b"Test message";
        // construct directly: an attacker sends raw bytes, bypassing
        // Signature::new's checks
        let zero_s = Signature { r: Fr::from(1u64), s: Fr::from(0u64), v: false };
        let zero_r = Signature { r: Fr::from(0u64), s: Fr::from(1u64), v: false };
        assert!(!verify(&zero_s, message, &key_pair.public_key));
        assert!(!verify(&zero_r, message, &key_pair.public_key));
    }

    #[test]
    fn test_low_s_enforced() {
        let key_pair = KeyPair::generate();
        let h = hash_message(b"Test message");
        let signature = sign(&key_pair.private_key, &h);
        assert!(signature.s.into_bigint() <= Fr::MODULUS_MINUS_ONE_DIV_TWO);
    }
}