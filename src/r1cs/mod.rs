//! Rank-1 Constraint System (R1CS) verification.
//!
//! Implements:
//! - Field-based R1CS verification via Hadamard product (Lw o Rw = Ow)
//! - Pairing-based R1CS verification using BN254
//!
//! See: Theory document, Section 6.1, 6.2 (R1CS)

use std::ops::Mul;
use ark_ff::{Field, PrimeField};
use ark_ec::{AffineRepr, CurveGroup, Group};
use ark_ec::pairing::Pairing;
use ark_bn254::{Bn254, Fr, G1Projective, G2Projective, G1Affine, G2Affine};
use ark_std::Zero;


pub fn verify_r1cs_with_hadamard_product<F: PrimeField>(
    l: &[Vec<F>],
    r: &[Vec<F>],
    o: &[Vec<F>],
    witness: &[F],
) -> bool {
    let num_rows = l.len();
    assert_eq!(num_rows, r.len(), "L and R must have the same number of rows");
    assert_eq!(num_rows, o.len(), "L and O must have the same number of rows");

    if !rank_column_check(l, r, o) {
        false;
    }    

    let lw = matrix_vector_multiply(l, witness);
    let rw = matrix_vector_multiply(r, witness);
    let ow = matrix_vector_multiply(o, witness);
    
    let hadamard = hadamard_product(&lw, &rw);

    for i in 0..hadamard.len() {
        if hadamard[i] != ow[i] {
            return false;
        }
    }

    true
}

fn matrix_vector_multiply<F: Field>(
    matrix: &[Vec<F>],
    vector: &[F],
) -> Vec<F> {
    let num_rows = matrix.len();
    if num_rows == 0{
        return Vec::new();
    }
    let num_cols = matrix[0].len();
    assert_eq!(num_cols, vector.len(), "Matrix columns must match vector length");

    let mut result = Vec::with_capacity(num_rows);

    for i in 0..num_rows {
        assert_eq!(matrix[i].len(), num_cols, "All matrix rows must have the same number of columns. Row 0 has {} columns, but row {} has {} columns", num_cols, i, matrix[i].len());
        
        let mut sum = F::zero();
        for j in 0..num_cols {
            sum = sum + (matrix[i][j]*vector[j]);
        }
        result.push(sum);
    }
    result
}

fn hadamard_product<F: Field>(
    v1: &[F],
    v2: &[F],
) -> Vec<F> {
    let num_elements = v1.len();
    assert_eq!(num_elements, v2.len(), "Vectors must have equal length");
    let mut result = Vec::with_capacity(v1.len());

    for i in 0..num_elements {
        result.push(v1[i] * v2[i]);
    }
    result
}

pub fn verify_r1cs(
    l: &[Vec<Fr>],
    r: &[Vec<Fr>],
    o: &[Vec<Fr>],
    witness_g1: &[G1Affine],
    witness_g2: &[G2Affine],
) -> bool {
    let num_rows = l.len();
    assert_eq!(num_rows, r.len(), "L and R must have the same number of rows");
    assert_eq!(num_rows, o.len(), "L and O must have the same number of rows");

    if !rank_column_check(l, r, o) {
        false;
    }

    let l_s1 = matrix_point_multiply::<G1Projective>(l, witness_g1);
    let r_s2 = matrix_point_multiply::<G2Projective>(r, witness_g2);
    let o_s1 = matrix_point_multiply::<G1Projective>(o, witness_g1);

    let g2_generator = G2Affine::generator();

    for i in 0..num_rows {
        let lhs = Bn254::pairing( l_s1[i], r_s2[i]);
        let rhs = Bn254::pairing( o_s1[i], g2_generator);

        if lhs != rhs {
            return false
        }
    }
    true
}

fn matrix_point_multiply<C>(
    matrix: &[Vec<Fr>],
    points: &[C::Affine],
) -> Vec<C::Affine> 
where
    // CurveGroup enables C::Affine and into_group()
    // Mul<Fr, C>: mapping C x Fr -> C exists
    C: CurveGroup + Mul<Fr, Output = C>,
{
    let num_rows = matrix.len();
    if num_rows == 0 {
        return Vec::new();
    }
    let num_cols = matrix[0].len();
    assert_eq!(num_cols, points.len(), "Matrix columns must match array of points length");

    let mut result = Vec::with_capacity(num_rows);

    for i in 0..num_rows {
        assert_eq!(matrix[i].len(), num_cols, "All matrix rows must have the same number of columns. Row 0 has {} columns, but row {} has {} columns", num_cols, i, matrix[i].len());

        let mut sum = points[0].into_group() * Fr::zero();
        for j in 0..num_cols {
            sum = sum + (points[j].into_group() * matrix[i][j]);
        }
        result.push(sum.into());
    }
    result
}

pub fn create_witness_points(scalars: &[Fr]) -> (Vec<G1Affine>, Vec<G2Affine>) {
    let g1_gen = G1Projective::generator();
    let g2_gen = G2Projective::generator();
    
    let g1: Vec<G1Affine> = scalars.iter().map(|s| (g1_gen * s).into()).collect();
    let g2: Vec<G2Affine> = scalars.iter().map(|s| (g2_gen * s).into()).collect();
    (g1, g2)
}

pub fn rank_column_check<F: Field>(
    l: &[Vec<F>],
    r: &[Vec<F>],
    o: &[Vec<F>],
) -> bool {
    let num_rows = l.len();
    for i in 0..num_rows {
        let num_cols = l[i].len();
        assert_eq!(num_cols, r[i].len(), "Row{}: L and R must have the same number of columns", i);
        assert_eq!(num_cols, o[i].len(), "Row{}: L and O must have the same number of columns", i);

        for j in 0..num_cols {
            let l_nonzero = !l[i][j].is_zero();
            let r_nonzero = !r[i][j].is_zero();
            let o_nonzero = !o[i][j].is_zero();

            let count = (l_nonzero as u32) + (r_nonzero as u32) + (o_nonzero as u32);
            // check if constraints have one multiplication operation, needed to encode via bilinear pairings
            if count >= 2 {
                return false;
            }
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_std::One;

    #[test]
    fn test_verify_r1cs_field_based() {
        // Simple case: z = x * y
        let x = Fr::from(5u64);
        let y = Fr::from(7u64);
        let z = x * y;
        let l = vec![vec![Fr::zero(), Fr::zero(), Fr::one(), Fr::zero()]];
        let r = vec![vec![Fr::zero(), Fr::zero(), Fr::zero(), Fr::one()]];
        let o = vec![vec![Fr::zero(), Fr::one(), Fr::zero(), Fr::zero()]];
        let witness = vec![Fr::one(), z, x, y];
        assert!(verify_r1cs_with_hadamard_product(&l, &r, &o, &witness));

        // Complex case: x*y*z*u with intermediate variables
        let x = Fr::from(3u64);
        let y = Fr::from(5u64);
        let z = Fr::from(7u64);
        let u = Fr::from(11u64);
        let v1 = x * y;
        let v2 = z * u;
        let r = v1 * v2;
        let witness = vec![Fr::one(), r, x, y, z, u, v1, v2];
        let l = vec![
            vec![Fr::zero(), Fr::zero(), Fr::one(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero()],
            vec![Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::one(), Fr::zero(), Fr::zero(), Fr::zero()],
            vec![Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::one(), Fr::zero()],
        ];
        let r = vec![
            vec![Fr::zero(), Fr::zero(), Fr::zero(), Fr::one(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero()],
            vec![Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::one(), Fr::zero(), Fr::zero()],
            vec![Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::one()],
        ];
        let o = vec![
            vec![Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::one(), Fr::zero()],
            vec![Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::one()],
            vec![Fr::zero(), Fr::one(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero()],
        ];
        assert!(verify_r1cs_with_hadamard_product(&l, &r, &o, &witness));
    }

    #[test]
    fn test_verify_r1cs_field_based_failures() {
        // Wrong witness
        let x = Fr::from(5u64);
        let y = Fr::from(7u64);
        let z_wrong = Fr::from(100u64);
        let l = vec![vec![Fr::zero(), Fr::zero(), Fr::one(), Fr::zero()]];
        let r = vec![vec![Fr::zero(), Fr::zero(), Fr::zero(), Fr::one()]];
        let o = vec![vec![Fr::zero(), Fr::one(), Fr::zero(), Fr::zero()]];
        let witness = vec![Fr::one(), z_wrong, x, y];
        assert!(!verify_r1cs_with_hadamard_product(&l, &r, &o, &witness));

        // Invalid structure
        let l = vec![vec![Fr::one(), Fr::one()]];
        let r = vec![vec![Fr::one(), Fr::zero()]];
        let o = vec![vec![Fr::zero(), Fr::zero()]];
        let witness = vec![Fr::one(), Fr::one()];
        assert!(!verify_r1cs_with_hadamard_product(&l, &r, &o, &witness));
    }

    #[test]
    fn test_verify_r1cs_pairing_based() {
        let x = Fr::from(3u64);
        let y = Fr::from(5u64);
        let z = Fr::from(7u64);
        let u = Fr::from(11u64);
        let v1 = x * y;
        let v2 = z * u;
        let r = v1 * v2;
        let witness_scalars = vec![Fr::one(), r, x, y, z, u, v1, v2];
        let (witness_g1, witness_g2) = create_witness_points(&witness_scalars);

        let l = vec![
            vec![Fr::zero(), Fr::zero(), Fr::one(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero()],
            vec![Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::one(), Fr::zero(), Fr::zero(), Fr::zero()],
            vec![Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::one(), Fr::zero()],
        ];
        let r = vec![
            vec![Fr::zero(), Fr::zero(), Fr::zero(), Fr::one(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero()],
            vec![Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::one(), Fr::zero(), Fr::zero()],
            vec![Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::one()],
        ];
        let o = vec![
            vec![Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::one(), Fr::zero()],
            vec![Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::one()],
            vec![Fr::zero(), Fr::one(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero()],
        ];
        assert!(verify_r1cs(&l, &r, &o, &witness_g1, &witness_g2));
    }

    #[test]
    fn test_verify_r1cs_pairing_based_failures() {
        // Wrong witness
        let x = Fr::from(5u64);
        let y = Fr::from(7u64);
        let z_wrong = Fr::from(100u64);
        let (witness_g1, witness_g2) = create_witness_points(&vec![Fr::one(), z_wrong, x, y]);
        let l = vec![vec![Fr::zero(), Fr::zero(), Fr::one(), Fr::zero()]];
        let r = vec![vec![Fr::zero(), Fr::zero(), Fr::zero(), Fr::one()]];
        let o = vec![vec![Fr::zero(), Fr::one(), Fr::zero(), Fr::zero()]];
        assert!(!verify_r1cs(&l, &r, &o, &witness_g1, &witness_g2));

        // Invalid structure
        let (witness_g1, witness_g2) = create_witness_points(&vec![Fr::one(), Fr::one()]);
        let l = vec![vec![Fr::one(), Fr::one()]];
        let r = vec![vec![Fr::one(), Fr::zero()]];
        let o = vec![vec![Fr::zero(), Fr::zero()]];
        assert!(!verify_r1cs(&l, &r, &o, &witness_g1, &witness_g2));
    }
}