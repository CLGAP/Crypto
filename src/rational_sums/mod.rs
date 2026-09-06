//! Rational sum verification on elliptic curves.
//!
//! Verifies that a/b + c/d = e/f by mapping each rational
//! to a curve point via scalar multiplication and checking
//! point addition equality. Note "=" means equality in Fr:
//! e' = e + r*f passes despite being a different rational number.
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
    
    // multiplicative inverse in Fr (cf. theory doc: a^(n-2) = a^-1 by Fermat;
    // arkworks computes it via extended Euclid instead)
    let s = num_fr * den_fr.inverse().expect("nonzero element of a prime field is invertible");

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
    fn test_verify_rational_sum_wrong_sum() {
        // 2/3 + 5/7 != 30/21
        let result = verify_rational_sum(
            BigUint::from(2u64),
            BigUint::from(3u64),
            BigUint::from(5u64),
            BigUint::from(7u64),
            BigUint::from(30u64),
            BigUint::from(21u64),
        );
        assert_eq!(result, Ok(false));
    }

    #[test]
    fn test_verify_rational_sum_zero_numerator() {
        // 0/3 + 5/7 = 5/7; the identity point is a legitimate summand
        let result = verify_rational_sum(
            BigUint::from(0u64),
            BigUint::from(3u64),
            BigUint::from(5u64),
            BigUint::from(7u64),
            BigUint::from(5u64),
            BigUint::from(7u64),
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