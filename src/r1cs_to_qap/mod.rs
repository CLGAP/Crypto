//! Polynomial-based verification primitives for zero-knowledge proofs.
//!
//! Implements:
//! - Lagrange interpolation
//! - Vector/matrix equality testing via polynomials (Schwartz-Zippel)
//! - Quadratic Arithmetic Programs (QAP)
//!
//! See: Theory document, Sections 6.3, 6.4, 6.5, 6.6

use ark_ff::{Field, PrimeField};
use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial, Polynomial};
use ark_std::{rand::Rng, Zero};

pub fn lagrange_interpolate<F: PrimeField + Field>(
    evaluation_points: &[F],
    values: &[F],
) -> DensePolynomial<F> {
    assert_eq!(evaluation_points.len(), values.len(), "Points and values must have same length");

    let mut result = DensePolynomial::zero();

    for i in 0..evaluation_points.len() {
        let basis = &lagrange_basis_polynomial(i, evaluation_points) * values[i];
        result += &basis;
    }
    //returning result alone leads to trailing 0s when interpolated polynomial has degree lower than number of points e.g., col [n, n] at [1, 2]
    DensePolynomial::from_coefficients_vec(result.coeffs().to_vec())
}

fn lagrange_basis_polynomial<F: PrimeField>(
    i: usize,
    evaluation_points: &[F],
) -> DensePolynomial<F> {
    let x_i = evaluation_points[i];
    let mut numerator = DensePolynomial::from_coefficients_vec(vec![F::one()]);
    let mut denominator = F::one();

    for j in 0..evaluation_points.len() {
        if j != i {
            let x_j = evaluation_points[j];
            assert_ne!(x_i, x_j, "Duplicate evaluation points at indices {} and {}: value = {}", i, j, x_i);
            let linear_factor = DensePolynomial::from_coefficients_vec(vec![-x_j, F::one()]);
            numerator = &numerator * &linear_factor;
            denominator = denominator * (x_i - x_j);
        }
    }
    let denom_inv = denominator.inverse().expect("Denominator must be non-zero");
    &numerator * denom_inv
}


pub fn evaluate_polynomial<F: Field>(
    poly: &DensePolynomial<F>,
    point: F,
) -> F {
    poly.evaluate(&point)
} 

pub fn verify_vector_equality<F: PrimeField, R: Rng>(
    v1: &[F],
    v2: &[F],
    evaluation_points: &[F],
    rng: &mut R,
) -> bool {
    let p1 = lagrange_interpolate(evaluation_points, v1);
    let p2 = lagrange_interpolate(evaluation_points, v2);

    let random_point = F::rand(rng);

    evaluate_polynomial(&p1,random_point) == evaluate_polynomial(&p2, random_point)
}

pub fn verify_matrix_vector_equality<F: PrimeField, R: Rng>(
    a: &[Vec<F>],
    b: &[Vec<F>],
    v: &[F],
    evaluation_points: &[F],
    rng: &mut R,
) -> bool {
    let master_poly_a = create_master_polynomial(a, v, evaluation_points);
    let master_poly_b = create_master_polynomial(b, v, evaluation_points);
    let random_point = F::rand(rng);
    
    evaluate_polynomial(&master_poly_a, random_point) == evaluate_polynomial(&master_poly_b, random_point)
}

fn create_master_polynomial<F: PrimeField>(
    matrix: &[Vec<F>],
    vector: &[F],
    evaluation_points: &[F],
) -> DensePolynomial<F> {
    let num_rows = matrix.len();
    if num_rows == 0 {
        return DensePolynomial::zero();
    }
    let num_cols = matrix[0].len();
    assert_eq!(num_cols, vector.len(), "Matrix columns must match vector length");
    assert_eq!(num_rows, evaluation_points.len(), "Matrix rows must match evaluation points length");
    for i in 1..num_rows {
        assert_eq!(matrix[i].len(), num_cols,"All matrix rows must have the same number of columns. Row 0 has {} columns, but row {} has {} columns", num_cols, i, matrix[i].len());
    }

    let column_polys = interpolate_matrix_columns(matrix, evaluation_points);
    combine_polynomials_with_witness(&column_polys, vector)
}

pub fn construct_qap<F: PrimeField>(
    l: &[Vec<F>],
    r: &[Vec<F>],
    o: &[Vec<F>],
    evaluation_points: &[F],
) -> (
    Vec<DensePolynomial<F>>,
    Vec<DensePolynomial<F>>,
    Vec<DensePolynomial<F>>,
) {
    let num_rows = l.len();
    assert_eq!(num_rows, r.len(), "L and R must have same number of rows");
    assert_eq!(num_rows, o.len(), "L and O must have same number of rows");
    assert_eq!(num_rows, evaluation_points.len(), "Number of rows must match evaluation point length");

    if num_rows == 0 {
        return (Vec::new(), Vec::new(), Vec::new())
    }
    let num_cols = l[0].len();
    assert_eq!(num_cols, r[0].len(), "L and R must have same number of columns");
    assert_eq!(num_cols, o[0].len(), "L and O must have same number of columns");
        
    for i in 1..num_rows {
        assert_eq!(l[i].len(), num_cols, "All L rows must have same number of columns");
        assert_eq!(r[i].len(), num_cols, "All R rows must have same number of columns");
        assert_eq!(o[i].len(), num_cols, "All O rows must have same number of columns");
    }

    let l_polys = interpolate_matrix_columns(l, evaluation_points);
    let r_polys = interpolate_matrix_columns(r, evaluation_points);
    let o_polys = interpolate_matrix_columns(o, evaluation_points);
    
    (l_polys, r_polys, o_polys)

}

fn interpolate_matrix_columns<F: PrimeField>(
    matrix: &[Vec<F>],
    evaluation_points: &[F],
) -> Vec<DensePolynomial<F>> {
    let num_rows = matrix.len();
    let num_cols = matrix[0].len();
    let mut column_polys = Vec::with_capacity(num_cols);
    for j in 0..num_cols {
        let column_values: Vec<F> = (0..num_rows).map(|i| matrix[i][j]).collect();
        let column_poly = lagrange_interpolate(evaluation_points, &column_values);
        column_polys.push(column_poly);
    }
    column_polys
}

pub fn target_polynomial<F: PrimeField>(
    evaluation_points: &[F],
) -> DensePolynomial<F> {
    let mut result = DensePolynomial::from_coefficients_vec(vec![F::one()]);

    for i in 0..evaluation_points.len() {
        result = &result * &DensePolynomial::from_coefficients_vec(vec![-evaluation_points[i], F::one()]);
    }
    result
}

pub fn verify_qap_satisfaction<F: PrimeField, R: Rng>(
    l_polys: &[DensePolynomial<F>],
    r_polys: &[DensePolynomial<F>],
    o_polys: &[DensePolynomial<F>],
    witness: &[F],
    target_poly: &DensePolynomial<F>,
    rng: &mut R,
) -> bool {
    let l_poly = combine_polynomials_with_witness(l_polys, witness);
    let r_poly = combine_polynomials_with_witness(r_polys, witness);
    let o_poly = combine_polynomials_with_witness(o_polys, witness);

    let h_poly = match compute_quotient_polynomial(&l_poly, &r_poly, &o_poly, target_poly) {
        Some(h) => h,
        None => return false,
    };

    let r = F::rand(rng);

    l_poly.evaluate(&r) * r_poly.evaluate(&r) - o_poly.evaluate(&r) == h_poly.evaluate(&r) * target_poly.evaluate(&r)
}

pub fn compute_quotient_polynomial<F: PrimeField>(
    l_poly: &DensePolynomial<F>,
    r_poly: &DensePolynomial<F>,
    o_poly: &DensePolynomial<F>,
    target_poly: &DensePolynomial<F>,
) -> Option<DensePolynomial<F>> {
    let numerator = &(l_poly * r_poly) - o_poly;

    let num_coeffs = numerator.coeffs();
    let denom_coeffs = target_poly.coeffs();

    // Zero polynomial is divisible by any non-zero polynomial; quotient is zero
    if num_coeffs.is_empty() || num_coeffs.iter().all(|c| c.is_zero()) {
        return Some(DensePolynomial::zero());
    }

    //non-zero divisor and quotient is poly
    if denom_coeffs.is_empty() || denom_coeffs.len() > num_coeffs.len() {
        return None;
    }

    let mut remainder = num_coeffs.to_vec();
    //+1 as  degree + 1 is coeff., need coeff
    let mut quotient_coeffs = vec![F::zero(); num_coeffs.len().saturating_sub(denom_coeffs.len())+1];

    for i in (0..=num_coeffs.len().saturating_sub(denom_coeffs.len())).rev() {
        let leading_pos = i + denom_coeffs.len() - 1;
        if leading_pos >= remainder.len() {
            continue;
        }
        //skip if no div needed
        if remainder[leading_pos].is_zero() {
            continue;
        }
        //for scaled div.
        let factor = remainder[leading_pos] / denom_coeffs[denom_coeffs.len() - 1];
        quotient_coeffs[i] = factor;

        for j in 0..denom_coeffs.len() {
            if i + j < remainder.len() {
                remainder[i + j] = remainder[i + j] - factor * denom_coeffs[j];
            }
        }
    }

    if remainder.iter().all(|&c| c.is_zero()) {
        Some(DensePolynomial::from_coefficients_vec(quotient_coeffs))
    } else {
        None
    }
}

pub fn combine_polynomials_with_witness<F: PrimeField>(
    column_polys: &[DensePolynomial<F>],
    witness: &[F],
) -> DensePolynomial<F> {
    assert_eq!(witness.len(), column_polys.len(), "Must have same amount of witness elements as interpolated polynomials");
    let mut sum = DensePolynomial::zero(); 

    for i in 0..witness.len() {
        sum = sum +  &column_polys[i] * witness[i];
    }
    sum
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bn254::Fr;
    use ark_std::{One, rand::{rngs::StdRng, SeedableRng}};

    #[test]
    fn test_lagrange_interpolate() {
        let points = vec![Fr::from(1u64), Fr::from(2u64), Fr::from(3u64)];
        let values = vec![Fr::from(2u64), Fr::from(4u64), Fr::from(6u64)];

        let poly = lagrange_interpolate(&points, &values);

        for i in 0..points.len() {
            let result = evaluate_polynomial(&poly, points[i]);
            assert_eq!(result, values[i], "Polynomial should evaluate to {} at point {}", values[i], points[i])
        }
    }

    #[test]
    fn test_verify_vector_equality() {
        let mut rng = StdRng::seed_from_u64(23);
        let points = vec![Fr::from(1u64), Fr::from(2u64), Fr::from(3u64)];
        let v1 = vec![Fr::from(2u64), Fr::from(4u64), Fr::from(6u64)];
        let v2 = vec![Fr::from(2u64), Fr::from(4u64), Fr::from(6u64)];
        assert!(verify_vector_equality(&v1, &v2, &points, &mut rng), "Equal vectors should verify");
        
        
        let v3 = vec![Fr::from(2u64), Fr::from(4u64), Fr::from(8u64)];
        assert!(!verify_vector_equality(&v1, &v3, &points, &mut rng), "Unequal vectors should not verify (with probability)");                
    }

    #[test]
    fn test_verify_matrix_vector_equality() {
        let mut rng = StdRng::seed_from_u64(42);
    
        let a = vec![
            vec![Fr::from(1u64), Fr::from(2u64)],
            vec![Fr::from(3u64), Fr::from(4u64)],
        ];
        let b = vec![
            vec![Fr::from(1u64), Fr::from(2u64)],
            vec![Fr::from(3u64), Fr::from(4u64)],
        ];
        let v = vec![Fr::from(5u64), Fr::from(6u64)];
        let points = vec![Fr::from(1u64), Fr::from(2u64)];
    
        assert!(verify_matrix_vector_equality(&a, &b, &v, &points, &mut rng), "equal matrices should verify");
    
        let c = vec![
            vec![Fr::from(9u64), Fr::from(9u64)],
            vec![Fr::from(9u64), Fr::from(9u64)],
        ];
        assert!(!verify_matrix_vector_equality(&a, &c, &v, &points, &mut rng), "different matrices should not verify");
    }

    #[test]
    fn test_construct_qap() {
        let l = vec![
            vec![Fr::zero(), Fr::zero(), Fr::from(3u64), Fr::zero(), Fr::zero(), Fr::zero()],
            vec![Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::one(), Fr::zero()],
            vec![Fr::zero(), Fr::zero(), Fr::one(), Fr::zero(), Fr::zero(), Fr::zero()],
        ];
        let r = vec![
            vec![Fr::zero(), Fr::zero(), Fr::one(), Fr::zero(), Fr::zero(), Fr::zero()],
            vec![Fr::zero(), Fr::zero(), Fr::zero(), Fr::one(), Fr::zero(), Fr::zero()],
            vec![Fr::zero(), Fr::zero(), Fr::zero(), Fr::from(5u64), Fr::zero(), Fr::zero()],
        ];
        let o = vec![
            vec![Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::one(), Fr::zero()],
            vec![Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::one()],
            vec![-Fr::from(3u64), Fr::one(), Fr::one(), Fr::from(2u64), Fr::zero(), -Fr::one()],
        ];
    
        let eval_points = vec![Fr::from(1u64), Fr::from(2u64), Fr::from(3u64)];
        let (l_polys, r_polys, o_polys) = construct_qap(&l, &r, &o, &eval_points);
    
        assert_eq!(l_polys.len(), 6);
        assert_eq!(r_polys.len(), 6);
        assert_eq!(o_polys.len(), 6);
    
        // Each basis poly should interpolate its column:
        // L column 2 = [3, 0, 1] at points [1, 2, 3]
        assert_eq!(l_polys[2].evaluate(&Fr::from(1u64)), Fr::from(3u64));
        assert_eq!(l_polys[2].evaluate(&Fr::from(2u64)), Fr::zero());
        assert_eq!(l_polys[2].evaluate(&Fr::from(3u64)), Fr::one());
    
        // R column 3 = [0, 1, 5] at points [1, 2, 3]
        assert_eq!(r_polys[3].evaluate(&Fr::from(1u64)), Fr::zero());
        assert_eq!(r_polys[3].evaluate(&Fr::from(2u64)), Fr::one());
        assert_eq!(r_polys[3].evaluate(&Fr::from(3u64)), Fr::from(5u64));
    
        // O column 0 = [0, 0, -3] at points [1, 2, 3]
        assert_eq!(o_polys[0].evaluate(&Fr::from(1u64)), Fr::zero());
        assert_eq!(o_polys[0].evaluate(&Fr::from(2u64)), Fr::zero());
        assert_eq!(o_polys[0].evaluate(&Fr::from(3u64)), -Fr::from(3u64));
    }

    #[test]
    fn test_verify_qap_satisfaction() {
        let mut rng = StdRng::seed_from_u64(42);

        let x = Fr::from(100u64);
        let y = Fr::from(100u64);
        let v1 = Fr::from(3u64) * x * x;
        let v2 = v1 * y;
        let out = Fr::from(3u64) * x * x * y
                + Fr::from(5u64) * x * y
                - x
                - Fr::from(2u64) * y 
                + Fr::from(3u64);

        let witness = vec![Fr::one(), out, x, y, v1, v2];

        let neg3 = -Fr::from(3u64);
        let neg1 = -Fr::one();

        let l = vec![
            vec![Fr::zero(), Fr::zero(), Fr::from(3u64), Fr::zero(), Fr::zero(), Fr::zero()],
            vec![Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::one(), Fr::zero()],
            vec![Fr::zero(), Fr::zero(), Fr::one(), Fr::zero(), Fr::zero(), Fr::zero()],
        ];
        let r = vec![
            vec![Fr::zero(), Fr::zero(), Fr::one(), Fr::zero(), Fr::zero(), Fr::zero()],
            vec![Fr::zero(), Fr::zero(), Fr::zero(), Fr::one(), Fr::zero(), Fr::zero()],
            vec![Fr::zero(), Fr::zero(), Fr::zero(), Fr::from(5u64), Fr::zero(), Fr::zero()],
        ];
        let o = vec![
            vec![Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::one(), Fr::zero()],
            vec![Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::one()],
            vec![neg3, Fr::one(), Fr::one(), Fr::from(2u64), Fr::zero(), neg1],
        ];

        let eval_points = vec![Fr::from(1u64), Fr::from(2u64), Fr::from(3u64)];
        let (l_polys, r_polys, o_polys) = construct_qap(&l, &r, &o, &eval_points);
        let t_poly = target_polynomial(&eval_points);

        let result = verify_qap_satisfaction(&l_polys, &r_polys, &o_polys, &witness, &t_poly, &mut rng);
        assert!(result, "valid witness should satisfy QAP");

    }
}

