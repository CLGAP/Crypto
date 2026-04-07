//! Rational sum verification on elliptic curves.
//!
//! Verifies that a/b + c/d = e/f by mapping each rational
//! to a curve point via scalar multiplication and checking
//! point addition equality.
//!
//! See: Theory document, Section 4 (Homomorphic Hiding)

use ark_ec::Group;
use ark_ff::Field;
use ark_secp256k1::Projective;
use ark_std::Zero;
use num_bigint::BigUint;
use crate::utils::biguint_to_fr;

fn quotient_point(num: BigUint, den: BigUint) -> Result<Projective, String> {
    let num_fr = biguint_to_fr::<ark_secp256k1::Fr>(&num)?;
    let den_fr = biguint_to_fr::<ark_secp256k1::Fr>(&den)?;

    if den.is_zero() || den_fr.is_zero() {
        return Err("Denominator has no multiplicative inverse modulo curve: denominator is zero.".to_string());
    }
    
    // Fermat's Little Theorem: x^ -1 % C == x^ (C-2) % C
    let den_inverse = den_fr.inverse();
    if den_inverse.is_none() {
        return Err("Denominator has no multiplicative inverse modulo curve: gcd(den, curve_order) != 1.".to_string());
    }

    let s = num_fr * den_inverse.unwrap();

    if s.is_zero() {
        return Err("Point at infinity is not a practical identity to verify rational sum".to_string());
    }

    let g = Projective::generator();
    let point = g*s;

    Ok(point)
}

pub fn verify_rational_sum(
    a_num: BigUint,
    a_den: BigUint,
    b_num: BigUint,
    b_den: BigUint,
    num: BigUint,
    den: BigUint,
) -> Result<bool, String> {
    let a = quotient_point(a_num, a_den)?;
    let b = quotient_point(b_num, b_den)?;
    let c = quotient_point(num, den)?;
    
    let sum = a + b;
    Ok(sum == c)
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_verify_rational_sum() {
        let result = verify_rational_sum(
            BigUint::from(2u64),
            BigUint::from(3u64),
            BigUint::from(5u64),
            BigUint::from(7u64),
            BigUint::from(29u64),
            BigUint::from(21u64),
        );
        assert_eq!(result, Ok(true));
    }
    #[test]
    fn test_quotient_point_zero_denominator() {
        let result = quotient_point(BigUint::from(1u64), BigUint::from(0u64));
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("zero"));
    }

    #[test]
    fn test_quotient_point_valid() {
        let result = quotient_point(BigUint::from(2u64), BigUint::from(3u64));
        assert!(result.is_ok());
        let point = result.unwrap();
        assert!(!point.is_zero()); 
    }
}