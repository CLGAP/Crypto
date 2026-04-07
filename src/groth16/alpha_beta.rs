//! Simplified Groth16 (Pinocchio-style) with alpha/beta shifting.
//!
//! Adds alpha, beta blinding to the CRS to prevent proof element forgery,
//! but without public/private input separation (gamma, delta) or
//! zero-knowledge blinding (r, s).
//!
//! See: Theory document, Section 8.1, 8.2.

use ark_bn254::{Bn254, Fr, G1Affine, G1Projective, G2Affine, G2Projective};
use ark_ec::pairing::{Pairing, PairingOutput};
use ark_ec::{AffineRepr, Group};
use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial};
use ark_std::Zero;

use crate::qap_trusted_setup::{
    generate_srs, generate_srs_ht, inner_product, pad_to_len
};
use crate::r1cs_to_qap::{
    target_polynomial, construct_qap, combine_polynomials_with_witness, compute_quotient_polynomial, evaluate_polynomial
};

pub struct CRS {
    pub srs_g1: Box<[G1Affine]>,   
    pub srs_g2: Box<[G2Affine]>,   
    pub eta: Vec<G1Affine>,      
    pub alpha_g1: G1Affine,      
    pub beta_g2: G2Affine,       
    pub psi: Vec<G1Affine>,
}

pub struct Proof {
    pub a1: G1Affine,  
    pub b2: G2Affine,  
    pub c1: G1Affine,  
}

fn compute_psi(
    l_polys: &[DensePolynomial<Fr>],
    r_polys: &[DensePolynomial<Fr>],
    o_polys: &[DensePolynomial<Fr>],
    tau: Fr,
    alpha: Fr,
    beta: Fr,
) -> Vec<G1Affine> {
    let g1 = G1Projective::generator();
    let mut psi: Vec<G1Affine> = Vec::new();
    for i in 0..l_polys.len() {
        let scalar = beta * evaluate_polynomial(&l_polys[i], tau)  + alpha * evaluate_polynomial(&r_polys[i], tau) + evaluate_polynomial(&o_polys[i], tau);
        let point = G1Affine::from(g1 * scalar);
        psi.push(point);
    }
    psi
}

pub fn trusted_setup(
    l_matrix: &[Vec<Fr>],
    r_matrix: &[Vec<Fr>],
    o_matrix: &[Vec<Fr>],
    alpha: Fr,
    beta: Fr,
    tau: Fr,
) -> CRS {
    let deg = l_matrix.len();

    let srs_g1 = generate_srs::<G1Projective>(tau, deg.saturating_sub(1));
    let srs_g2 = generate_srs::<G2Projective>(tau, deg.saturating_sub(1));
    let eval_points: Vec<Fr> = (1..=l_matrix.len()).map(|i| Fr::from(i as u64)).collect();
    let target_poly = target_polynomial(&eval_points);
    let eta = generate_srs_ht(tau, &target_poly, deg.saturating_sub(2));
    let alpha_g1 = G1Affine::from(G1Projective::generator() * alpha); 
    let beta_g2 = G2Affine::from(G2Projective::generator() * beta);
    let (l_polys, r_polys, o_polys) = construct_qap(&l_matrix, &r_matrix, &o_matrix, &eval_points);
    let psi = compute_psi(&l_polys, &r_polys, &o_polys, tau, alpha, beta);

    CRS { srs_g1, srs_g2, eta, alpha_g1, beta_g2, psi }
}

pub fn prove(
    l_matrix: &[Vec<Fr>],
    r_matrix: &[Vec<Fr>],
    o_matrix: &[Vec<Fr>],
    witness: &[Fr],
    crs: &CRS,
) -> Option<Proof> {
    let eval_points: Vec<Fr> = (1..=l_matrix.len()).map(|i| Fr::from(i as u64)).collect();
    let t_poly = target_polynomial(&eval_points);
    let (l_polys, r_polys, o_polys) = construct_qap(&l_matrix, &r_matrix, &o_matrix, &eval_points);
    let l_poly = combine_polynomials_with_witness(&l_polys, witness);
    let r_poly = combine_polynomials_with_witness(&r_polys, witness);
    let o_poly = combine_polynomials_with_witness(&o_polys, witness);
    let h_poly = compute_quotient_polynomial(&l_poly, &r_poly, &o_poly, &t_poly)?;
    let l_coeffs= pad_to_len(l_poly.coeffs(), crs.srs_g1.len());
    let r_coeffs= pad_to_len(r_poly.coeffs(), crs.srs_g2.len());
    let h_coeffs= pad_to_len(h_poly.coeffs(), crs.eta.len());

    let a1 = G1Affine::from(crs.alpha_g1.into_group() + inner_product::<G1Projective>(&l_coeffs, &crs.srs_g1));
    let b2 = G2Affine::from(crs.beta_g2.into_group() + inner_product::<G2Projective>(&r_coeffs, &crs.srs_g2));

    let mut psi_sum = G1Projective::zero();
    for i in 0..witness.len() {
        psi_sum += crs.psi[i] * witness[i];
    }

    let c1 = G1Affine::from(psi_sum + inner_product::<G1Projective>(&h_coeffs, &crs.eta).into_group());

    Some(Proof { a1, b2, c1 })
}

pub fn verify(proof: &Proof, crs: &CRS) -> bool {
    let target_identity = PairingOutput::<Bn254>::zero();
    let ab = Bn254::pairing(-proof.a1, proof.b2);
    let alphabeta = Bn254::pairing(crs.alpha_g1, crs.beta_g2);
    let cg_2 = Bn254::pairing(proof.c1, G2Affine::generator());
    target_identity == ab + alphabeta + cg_2
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_ec::Group;
    use ark_std::One;

    fn simple_circuit_z_equals_xy() -> (Vec<Vec<Fr>>, Vec<Vec<Fr>>, Vec<Vec<Fr>>) {
        use ark_std::One;
        let l = vec![vec![Fr::zero(), Fr::zero(), Fr::one(), Fr::zero()]];
        let r = vec![vec![Fr::zero(), Fr::zero(), Fr::zero(), Fr::one()]];
        let o = vec![vec![Fr::zero(), Fr::one(), Fr::zero(), Fr::zero()]];
        (l, r, o)
    }

    #[test]
    fn test_trusted_setup_crs_structure() {
        let (l, r, o) = simple_circuit_z_equals_xy();
        let n = l.len();
        let m = l[0].len();
        let crs = trusted_setup(&l, &r, &o, Fr::from(11u64), Fr::from(13u64), Fr::from(4u64));
        assert_eq!(crs.srs_g1.len(), n);
        assert_eq!(crs.srs_g2.len(), n);
        assert_eq!(crs.eta.len(), n.saturating_sub(1).max(1));
        assert_eq!(crs.psi.len(), m);
    }

    #[test]
    fn test_prove_valid_witness_returns_some() {
        let (l, r, o) = simple_circuit_z_equals_xy();
        let witness = vec![Fr::one(), Fr::from(35u64), Fr::from(5u64), Fr::from(7u64)];
        let crs = trusted_setup(&l, &r, &o, Fr::from(11u64), Fr::from(13u64), Fr::from(4u64));
        let proof = prove(&l, &r, &o, &witness, &crs).unwrap();
        assert!(!proof.a1.is_zero() && !proof.b2.is_zero() && !proof.c1.is_zero());
    }

    #[test]
    fn test_prove_invalid_witness_returns_none() {
        let (l, r, o) = simple_circuit_z_equals_xy();
        let witness = vec![Fr::one(), Fr::from(100u64), Fr::from(5u64), Fr::from(7u64)];
        let crs = trusted_setup(&l, &r, &o, Fr::from(11u64), Fr::from(13u64), Fr::from(4u64));
        assert!(prove(&l, &r, &o, &witness, &crs).is_none());
    }

    #[test]
    fn test_verify_accepts_valid_proof() {
        let (l, r, o) = simple_circuit_z_equals_xy();
        let witness = vec![Fr::one(), Fr::from(35u64), Fr::from(5u64), Fr::from(7u64)];
        let crs = trusted_setup(&l, &r, &o, Fr::from(11u64), Fr::from(13u64), Fr::from(4u64));
        let proof = prove(&l, &r, &o, &witness, &crs).unwrap();
        assert!(verify(&proof, &crs))
    }

    #[test]
    fn test_verify_rejects_forged_proof() {
        // Pinocchio-style forgery (e.g. A=G1, B=2*G2, C=2*G1).
        // assert!(!verify(&forged_proof, &crs));
        let (l, r, o) = simple_circuit_z_equals_xy();
        let crs = trusted_setup(&l, &r, &o, Fr::from(11u64), Fr::from(13u64), Fr::from(4u64));

        let forged_proof = Proof {
            a1: G1Affine::generator(),
            b2: G2Affine::from(G2Projective::generator() * Fr::from(2u64)),
            c1: G1Affine::from(G1Projective::generator() * Fr::from(2u64))
        };
        assert!(!verify(&forged_proof, &crs))
    }

    #[test]
    fn test_prove_then_verify_end_to_end() {
        // x^3 + x + 5 =35
        let l = vec![
            vec![Fr::zero(), Fr::zero(), Fr::one(), Fr::zero(), Fr::zero(), Fr::zero()],
            vec![Fr::zero(), Fr::zero(), Fr::zero(), Fr::one(), Fr::zero(), Fr::zero()],
            vec![Fr::zero(), Fr::zero(), Fr::one(), Fr::zero(), Fr::one(), Fr::zero()],
            vec![Fr::from(5u64), Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::one()],
        ];
        let r = vec![
            vec![Fr::zero(), Fr::zero(), Fr::one(), Fr::zero(), Fr::zero(), Fr::zero()],
            vec![Fr::zero(), Fr::zero(), Fr::one(), Fr::zero(), Fr::zero(), Fr::zero()],
            vec![Fr::one(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero()],
            vec![Fr::one(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero()],
        ];
        let o = vec![
            vec![Fr::zero(), Fr::zero(), Fr::zero(), Fr::one(), Fr::zero(), Fr::zero()],
            vec![Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::one(), Fr::zero()],
            vec![Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::one()],
            vec![Fr::zero(), Fr::one(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero()],
        ];
        let witness = vec![
            Fr::one(),
            Fr::from(35u64),
            Fr::from(3u64),
            Fr::from(9u64),   // var1 = x*x = 9
            Fr::from(27u64),  // var2 = var1*x = 27
            Fr::from(30u64),  // var3 = var2+x = 30
        ];
        let crs = trusted_setup(&l, &r, &o, Fr::from(11u64), Fr::from(2u64), Fr::from(7u64));
        let proof = prove(&l, &r, &o, &witness, &crs).unwrap();
        assert!(verify(&proof, &crs));
    }
}