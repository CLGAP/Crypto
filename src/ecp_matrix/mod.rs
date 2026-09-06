//! Elliptic curve matrix-vector multiplication verification.
//!
//! Demonstrates the homomorphic property: M * (`s_i` * G) = (M * s)_i * G,
//! i.e., matrix multiplication commutes with the scalar-to-point mapping.
//!
//! See: Theory document, Section 5 (Encrypted Polynomial Evaluation)

use ark_secp256k1::{Fr, Projective};
use ark_std::Zero;
use ark_std::rand::{Rng, rngs::StdRng, SeedableRng};
use getrandom::getrandom;
use crate::utils::{random_scalar, point_from_scalar};

fn random_matrix(n: usize, entry_bound: u64) -> Vec<Vec<u64>> {
    let mut seed = [0u8; 32];
    getrandom(&mut seed).expect("getrandom failed");
    let mut rng = StdRng::from_seed(seed);
    (0..n)
        .map(|_| {
            (0..n)
                .map(|_| rng.gen_range(0..=entry_bound))
                .collect()
        })
        .collect()
}

fn matrix_times_vector_on_ecp(
    matrix: &[Vec<u64>],
    points: &[Projective],
) -> Vec<Projective> {
    let n = matrix.len();
    let mut result = Vec::with_capacity(n);

    for i in 0..n {
        let mut acc = Projective::zero();

        for j in 0..n {
            let coeff = Fr::from(matrix[i][j]);
            acc += points[j] * coeff;
        }

        result.push(acc);
    }

    result
}

fn matrix_times_scalars_as_points(
    matrix: &[Vec<u64>],
    scalars: &[Fr],
) -> Vec<Projective> {
    let n = matrix.len();
    let mut result = Vec::with_capacity(n);

    for i in 0..n {
        let mut sum = Fr::from(0u64);

        for j in 0..n {
            let coeff = Fr::from(matrix[i][j]);
            sum += coeff * scalars[j];
        }

        result.push(point_from_scalar(&sum));
    }

    result
}

#[must_use]
pub fn verify_matrix_ec_multiplication(n: usize, entry_bound: u64) -> bool {
    let matrix = random_matrix(n, entry_bound);
    let scalars: Vec<Fr> = (0..n).map(|_| random_scalar()).collect();
    let points: Vec<Projective> = scalars.iter().map(point_from_scalar).collect();

    let lhs = matrix_times_scalars_as_points(&matrix, &scalars);
    let rhs = matrix_times_vector_on_ecp(&matrix, &points);

    lhs == rhs
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_verify_ec_multiplication_2x2() {
        assert!(verify_matrix_ec_multiplication(2, 100));
    }

    #[test]
    fn test_verify_ec_multiplication_11x11() {
        assert!(verify_matrix_ec_multiplication(11, 100));
    }
}