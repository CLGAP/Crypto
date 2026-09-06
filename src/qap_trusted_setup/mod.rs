//! QAP evaluation on a trusted setup (structured reference string).
//!
//! See: Theory document, Section 7 (Trusted Setup and QAP Evaluation).

use ark_bn254::{Bn254, Fr, G1Affine, G1Projective, G2Affine, G2Projective};
use ark_ec::{AffineRepr, CurveGroup};
use ark_ec::pairing::Pairing;
use ark_ff::One;
use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial};
use ark_std::Zero;

use crate::r1cs_to_qap::{
    combine_polynomials_with_witness, compute_quotient_polynomial, construct_qap, evaluate_polynomial, target_polynomial
};

pub fn generate_srs<G: CurveGroup>(discrete_log: G::ScalarField, degree: usize) -> Box<[G::Affine]> {
    assert!(!discrete_log.is_zero(), "tau must not be zero");
    let g = G::generator();
    let mut power = G::ScalarField::one(); 
    let mut ret = Vec::with_capacity(degree+1);

    for _ in 0..=degree {
        ret.push((g*power).into());
        power *= discrete_log;
    };
    ret.into_boxed_slice()
}

pub fn inner_product<G: CurveGroup>(poly_coeffs: &[G::ScalarField], srs: &[G::Affine]) -> G::Affine {
    
    assert_eq!(srs.len(), poly_coeffs.len(), "Number of coefficients must equal number of constraints in QAP");
    
    let mut ret = G::zero();

    for i in 0..srs.len() {
        ret += srs[i].into_group() * poly_coeffs[i];
    }

    ret.into()
}

#[must_use]
pub fn pad_to_len(coeffs: &[Fr], len: usize) -> Vec<Fr> {
    assert!(coeffs.len() <= len, "Coefficient count {} exceeds target length {}", coeffs.len(), len);
    let mut v = coeffs.to_vec();
    v.resize(len, Fr::zero());
    v
}

#[must_use]
pub fn generate_srs_ht(tau: Fr, target_poly: &DensePolynomial<Fr>, h_degree: usize) -> Vec<G1Affine> {
    let srs = generate_srs::<G1Projective>(tau, h_degree);
    let target_poly_at_tau = evaluate_polynomial(target_poly, tau);
    assert_ne!(target_poly_at_tau, Fr::zero(), "Tau was a root of the target polynomial");
    srs.into_vec().into_iter()
        .map(|p| G1Affine::from(p.into_group() * target_poly_at_tau))
        .collect()
}

pub struct Proof {
    pub a1: G1Affine,  // A = L(\tau) G_1
    pub b2: G2Affine,  // B = R(\tau) G_2
    pub c1: G1Affine,  // C = (O(\tau) + H(\tau)t(\tau)) * G_1
}

//Illustrative: Prover would never usually see tau. Real setups sample tau,
//publish the SRS, and destroy the scalar.
#[must_use]
pub fn prove(
    l_matrix: &[Vec<Fr>],
    r_matrix: &[Vec<Fr>],
    o_matrix: &[Vec<Fr>],
    witness: &[Fr],
    tau: Fr,
) -> Option<Proof> {
    let constraints = l_matrix.len();
    assert!(constraints > 0, "Cannot prove with zero constraints");

    let eval_points: Vec<Fr> = (1..=constraints).map(|i| Fr::from(i as u64)).collect();
    let (l_polys, r_polys, o_polys) = construct_qap(l_matrix, r_matrix, o_matrix, &eval_points);

    let l_poly = combine_polynomials_with_witness(&l_polys, witness);
    let r_poly = combine_polynomials_with_witness(&r_polys, witness);
    let o_poly = combine_polynomials_with_witness(&o_polys, witness);

    let t_poly = target_polynomial(&eval_points);
    let h_poly = compute_quotient_polynomial(&l_poly, &r_poly, &o_poly, &t_poly)?;

    let srs_g1 = generate_srs::<G1Projective>(tau , constraints-1);
    let srs_g2 = generate_srs::<G2Projective>(tau , constraints-1);
    let srs_ht = generate_srs_ht(tau, &t_poly, constraints.saturating_sub(2));

    let l_coeffs = pad_to_len(l_poly.coeffs(), srs_g1.len());
    let r_coeffs = pad_to_len(r_poly.coeffs(), srs_g2.len());
    let o_coeffs = pad_to_len(o_poly.coeffs(), srs_g1.len());
    let h_coeffs = pad_to_len(h_poly.coeffs(), srs_ht.len());

    let a1 = inner_product::<G1Projective>(&l_coeffs, &srs_g1);
    let b2 = inner_product::<G2Projective>(&r_coeffs, &srs_g2);
    let o1 = inner_product::<G1Projective>(&o_coeffs, &srs_g1);    
    let ht1 = inner_product::<G1Projective>(&h_coeffs, &srs_ht);
    let c1 = G1Affine::from(o1 + ht1);

    if a1.is_zero() || b2.is_zero() || c1.is_zero() {
        return None;
    }

    Some( Proof {a1, b2, c1})
}

// Soundness caveat: this checks only e(A,B) = e(C,G2). For any x, y the
// triple (xG1, yG2, xyG1) passes with no circuit knowledge; groth16/ closes
// exactly this gap. See theory.pdf, Section 7.3 (Pinocchio-style scheme) and
// the attack scenarios in Section 8.1.
#[must_use]
pub fn verify(proof: &Proof) -> bool {
    if proof.a1.is_zero() || proof.b2.is_zero() || proof.c1.is_zero() {
        return false;
    }
    let g2 = G2Affine::generator();
    let e_ab = Bn254::pairing(proof.a1, proof.b2);
    let e_cg2 = Bn254::pairing(proof.c1, g2);
    e_ab == e_cg2
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_std::One;
    use ark_ec::Group;
    use ark_poly::Polynomial;

    #[test]
    fn test_generate_srs_g1() {
        let tau = Fr::from(3u64);
        let srs = generate_srs::<G1Projective>(tau, 2);
        assert_eq!(srs.len(), 3);
        let g = G1Projective::generator();
        assert_eq!(srs[0], G1Affine::from(g));
        assert_eq!(srs[1], G1Affine::from(g * Fr::from(3u64)));
        assert_eq!(srs[2], G1Affine::from(g * Fr::from(9u64)));
    }

    #[test]
    fn test_generate_srs_g2() {
        let tau = Fr::from(3u64);
        let srs = generate_srs::<G2Projective>(tau, 2);
        assert_eq!(srs.len(), 3);
        let g = G2Projective::generator();
        assert_eq!(srs[0], G2Affine::from(g));
        assert_eq!(srs[1], G2Affine::from(g * Fr::from(3u64)));
        assert_eq!(srs[2], G2Affine::from(g * Fr::from(9u64)));
    }

    #[test]
    fn test_inner_product_g1() {
        let tau = Fr::from(2u64);
        let srs = generate_srs::<G1Projective>(tau, 2);
        assert_eq!(srs.len(), 3);
        let coeffs = vec![Fr::from(1u64), Fr::from(0u64), Fr::from(4u64)];
        let result = inner_product::<G1Projective>(&coeffs, &srs);
        let expected = G1Affine::from(G1Projective::generator() * Fr::from(17u64));
        assert_eq!(result, expected);
    }

    #[test]
    fn test_inner_product_g2() {
        let tau = Fr::from(2u64);
        let srs = generate_srs::<G2Projective>(tau, 2);
        let coeffs = vec![Fr::from(1u64), Fr::from(0u64), Fr::from(4u64)];
        let result = inner_product::<G2Projective>(&coeffs, &srs);
        let expected = G2Affine::from(G2Projective::generator() * Fr::from(17u64));
        assert_eq!(result, expected);
    }


    #[test]
    fn test_generate_srs_ht() {
        use crate::r1cs_to_qap::target_polynomial;
        let tau = Fr::from(4u64); // not a root of t(x)
        let eval_points = vec![Fr::from(1u64), Fr::from(2u64), Fr::from(3u64)];
        let t_poly = target_polynomial(&eval_points);
        let h_degree = 1; // 3 constraints, max 3-2
        let srs = generate_srs_ht(tau, &t_poly, h_degree);
        assert_eq!(srs.len(), 2);
        let g = G1Projective::generator();
        let t_at_tau = t_poly.evaluate(&tau);
        assert_eq!(srs[0], G1Affine::from(g * t_at_tau));
        assert_eq!(srs[1], G1Affine::from(g * (tau * t_at_tau)));
    }

    #[test]
    fn test_prove_and_verify_simple() {
        // single constraint: z = x * y
        let x = Fr::from(5u64);
        let y = Fr::from(7u64);
        let z = x * y;
        let l_matrix = vec![vec![Fr::zero(), Fr::zero(), Fr::one(), Fr::zero()]];
        let r_matrix = vec![vec![Fr::zero(), Fr::zero(), Fr::zero(), Fr::one()]];
        let o_matrix = vec![vec![Fr::zero(), Fr::one(), Fr::zero(), Fr::zero()]];
        let witness = vec![Fr::one(), z, x, y];
        let tau = Fr::from(4u64); 
        let proof = prove(&l_matrix, &r_matrix, &o_matrix, &witness, tau).expect("valid witness should produce proof");
        assert!(verify(&proof));
    }

    #[test]
    fn test_prove_and_verify_r1cs() {
        // 3x^2y + 5xy − x − 2y + 3 (3 constraints, 6 witness elements)
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
        let l_matrix = vec![
            vec![Fr::zero(), Fr::zero(), Fr::from(3u64), Fr::zero(), Fr::zero(), Fr::zero()],
            vec![Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::one(), Fr::zero()],
            vec![Fr::zero(), Fr::zero(), Fr::one(), Fr::zero(), Fr::zero(), Fr::zero()],
        ];
        let r_matrix = vec![
            vec![Fr::zero(), Fr::zero(), Fr::one(), Fr::zero(), Fr::zero(), Fr::zero()],
            vec![Fr::zero(), Fr::zero(), Fr::zero(), Fr::one(), Fr::zero(), Fr::zero()],
            vec![Fr::zero(), Fr::zero(), Fr::zero(), Fr::from(5u64), Fr::zero(), Fr::zero()],
        ];
        let o_matrix = vec![
            vec![Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::one(), Fr::zero()],
            vec![Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::one()],
            vec![neg3, Fr::one(), Fr::one(), Fr::from(2u64), Fr::zero(), neg1],
        ];
        let tau = Fr::from(4u64); // not 1, 2, or 3 (roots of target poly)
        let proof = prove(&l_matrix, &r_matrix, &o_matrix, &witness, tau).expect("valid witness should produce proof");
        assert!(verify(&proof));
    }

    #[test]
    fn test_invalid_witness_rejected() {
        let x = Fr::from(5u64);
        let y = Fr::from(7u64);
        let z_wrong = Fr::from(100u64);
        let l_matrix = vec![vec![Fr::zero(), Fr::zero(), Fr::one(), Fr::zero()]];
        let r_matrix = vec![vec![Fr::zero(), Fr::zero(), Fr::zero(), Fr::one()]];
        let o_matrix = vec![vec![Fr::zero(), Fr::one(), Fr::zero(), Fr::zero()]];
        let witness = vec![Fr::one(), z_wrong, x, y];
        let tau = Fr::from(4u64);
        assert!(prove(&l_matrix, &r_matrix, &o_matrix, &witness, tau).is_none(), "invalid witness must not produce a proof");
    }

    #[test]
    fn test_verify_accepts_forgery() {
        // documents the soundness caveat on verify: for any x, y the triple
        // (xG1, yG2, xyG1) passes with no circuit knowledge
        let x = Fr::from(6u64);
        let y = Fr::from(7u64);
        let forged = Proof {
            a1: G1Affine::from(G1Projective::generator() * x),
            b2: G2Affine::from(G2Projective::generator() * y),
            c1: G1Affine::from(G1Projective::generator() * (x * y)),
        };
        assert!(verify(&forged), "verify has no circuit binding; forgery passes by design");
    }
}
