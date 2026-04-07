//! Full Groth16 zk-SNARK protocol.
//!
//! Extends the alpha-beta scheme with:
//! - Public/private input separation via gamma and delta
//! - Zero-knowledge blinding via random r, s
//! - Public input binding in the verification equation
//!
//! See: Theory document, Section 8.3, 8.4, 8.5.

use ark_bn254::{Bn254, Fr, G1Affine, G1Projective, G2Affine, G2Projective};
use ark_ec::pairing::Pairing;
use ark_ec::{AffineRepr, Group};
use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial};
use ark_std::{rand::{rngs::StdRng, SeedableRng}, Zero, UniformRand};
use ark_ff::Field;
use getrandom::getrandom;

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
    pub beta_g1: G1Affine,      
    pub beta_g2: G2Affine,       
    pub psi: Vec<G1Affine>,
    pub gamma_g2: G2Affine,
    pub delta_g1: G1Affine,
    pub delta_g2: G2Affine,
    pub ell: usize,
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
    gamma: Fr,
    delta: Fr,
    ell: usize,
) -> Vec<G1Affine> {
    let g1 = G1Projective::generator();
    let mut psi: Vec<G1Affine> = Vec::new();
    for i in 0..l_polys.len() {
        let scalar = beta * evaluate_polynomial(&l_polys[i], tau)  + alpha * evaluate_polynomial(&r_polys[i], tau) + evaluate_polynomial(&o_polys[i], tau);
        let denom = if i < ell { gamma } else { delta };
        let point = G1Affine::from(g1 * (scalar/denom));
        psi.push(point);
    }
    psi
}

impl CRS {
    pub fn generate(
        l_matrix: &[Vec<Fr>],
        r_matrix: &[Vec<Fr>],
        o_matrix: &[Vec<Fr>],
        alpha: Fr, beta: Fr, tau: Fr,
        gamma: Fr, delta: Fr, ell: usize,
    ) -> CRS {
        assert!(
            gamma != delta,
            "gamma and delta must be distinct for public/private input separation"
        );
        let deg = l_matrix.len();
    
        let srs_g1 = generate_srs::<G1Projective>(tau, deg.saturating_sub(1));
        let srs_g2 = generate_srs::<G2Projective>(tau, deg.saturating_sub(1));
        let eval_points: Vec<Fr> = (1..=l_matrix.len()).map(|i| Fr::from(i as u64)).collect();
        let target_poly = target_polynomial(&eval_points);
        let delta_inv = delta.inverse().unwrap();
        let eta_raw = generate_srs_ht(tau, &target_poly, deg.saturating_sub(2));
        let eta = eta_raw.into_iter()
            .map(|p| G1Affine::from(p.into_group() * delta_inv)).collect();
        let alpha_g1 = G1Affine::from(G1Projective::generator() * alpha); 
        let beta_g1 = G1Affine::from(G1Projective::generator() * beta);
        let beta_g2 = G2Affine::from(G2Projective::generator() * beta);
        let gamma_g2 = G2Affine::from(G2Projective::generator() * gamma);
        let delta_g1 = G1Affine::from(G1Projective::generator() * delta);
        let delta_g2 = G2Affine::from(G2Projective::generator() * delta);
        let (l_polys, r_polys, o_polys) = construct_qap(&l_matrix, &r_matrix, &o_matrix, &eval_points);
        let psi = compute_psi(&l_polys, &r_polys, &o_polys, tau, alpha, beta, gamma, delta, ell);
    
        CRS { srs_g1, srs_g2, eta, alpha_g1, beta_g1, beta_g2, psi, gamma_g2, delta_g1, delta_g2, ell }
    }

}

impl Proof {
    pub fn new(
        l_matrix: &[Vec<Fr>],
        r_matrix: &[Vec<Fr>],
        o_matrix: &[Vec<Fr>],
        witness: &[Fr],
        crs: &CRS,
    ) -> Option<Self> {
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
        let mut seed = [0u8; 32];
        getrandom(&mut seed).expect("getrandom failed");
        let mut rng = StdRng::from_seed(seed);
        let r = Fr::rand(&mut rng);
        let s = Fr::rand(&mut rng);
        let a1 = G1Affine::from(crs.alpha_g1.into_group() + inner_product::<G1Projective>(&l_coeffs, &crs.srs_g1) + crs.delta_g1.into_group() * r);
        let b1 = crs.beta_g1.into_group() + inner_product::<G1Projective>(&r_coeffs, &crs.srs_g1) + crs.delta_g1.into_group() * s;
        let b2 = G2Affine::from(crs.beta_g2.into_group() + inner_product::<G2Projective>(&r_coeffs, &crs.srs_g2) + crs.delta_g2.into_group() * s);
    
        let mut psi_sum = G1Projective::zero();
        for i in crs.ell..witness.len() {
            psi_sum += crs.psi[i] * witness[i];
        }
    
        let c1 = G1Affine::from(psi_sum
            + inner_product::<G1Projective>(&h_coeffs, &crs.eta) 
            + a1.into_group() * s + b1 * r - crs.delta_g1.into_group() * (r * s));
    
        Some(Proof { a1, b2, c1 })
    }

    pub fn verify(&self, crs: &CRS, public_inputs: &[Fr]) -> bool {
        let lhs = Bn254::pairing(self.a1, self.b2);
        let alphabeta = Bn254::pairing(crs.alpha_g1, crs.beta_g2);
        let mut x = G1Projective::zero();
        for i in 0..crs.ell {
            x += crs.psi[i].into_group() * public_inputs[i];
        }
        let x_g1 = G1Affine::from(x);
        let xgamma = Bn254::pairing(x_g1, crs.gamma_g2);
        let cdelta = Bn254::pairing(self.c1, crs.delta_g2);
        let rhs = alphabeta + xgamma + cdelta;
        lhs == rhs
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_ec::Group;
    use ark_std::One;
    use proptest::{prelude::*, test_runner::Config};

    fn simple_circuit_z_equals_xy() -> (Vec<Vec<Fr>>, Vec<Vec<Fr>>, Vec<Vec<Fr>>) {
        use ark_std::One;
        let l = vec![vec![Fr::zero(), Fr::zero(), Fr::one(), Fr::zero()]];
        let r = vec![vec![Fr::zero(), Fr::zero(), Fr::zero(), Fr::one()]];
        let o = vec![vec![Fr::zero(), Fr::one(), Fr::zero(), Fr::zero()]];
        (l, r, o)
    }

    /// Helper: setup params for z=xy circuit (ell=2: constant 1 and output z are public).
    fn setup_z_equals_xy() -> (Vec<Vec<Fr>>, Vec<Vec<Fr>>, Vec<Vec<Fr>>, CRS) {
        let (l, r, o) = simple_circuit_z_equals_xy();
        let crs = CRS::generate(
            &l, &r, &o,
            Fr::from(11u64),
            Fr::from(13u64),
            Fr::from(4u64),
            Fr::one(),
            Fr::from(19u64),
            2,
        );
        (l, r, o, crs)
    }

    #[test]
    fn test_trusted_setup_crs_structure() {
        let (l, _, _, crs) = setup_z_equals_xy();
        let n = l.len();
        let m = l[0].len();
        assert_eq!(crs.srs_g1.len(), n);
        assert_eq!(crs.srs_g2.len(), n);
        assert_eq!(crs.eta.len(), n.saturating_sub(1).max(1));
        assert_eq!(crs.psi.len(), m);
        assert_eq!(crs.ell, 2);
        assert!(crs.ell <= m);
    }

    #[test]
    #[should_panic(expected = "gamma and delta must be distinct")]
    fn test_trusted_setup_rejects_gamma_equals_delta() {
        let (l, r, o) = simple_circuit_z_equals_xy();
        let _ = CRS::generate(
            &l, &r, &o,
            Fr::from(11u64),
            Fr::from(13u64),
            Fr::from(4u64),
            Fr::one(),
            Fr::one(), // gamma == delta; invalid
            2,
        );
    }

    #[test]
    fn test_prove_valid_witness_returns_some() {
        let (l, r, o, crs) = setup_z_equals_xy();
        // witness: [1, z, x, y] with z = x*y. 5*7 = 35.
        let witness = vec![Fr::one(), Fr::from(35u64), Fr::from(5u64), Fr::from(7u64)];
        let proof = Proof::new(&l, &r, &o, &witness, &crs).unwrap();
        assert!(!proof.a1.is_zero() && !proof.b2.is_zero() && !proof.c1.is_zero());
    }

    #[test]
    fn test_prove_witness_permutation_returns_none() {
        // Witness that fails the constraint: claim z=100 but 5*7=35.
        // An invalid witness
        let (l, r, o, crs) = setup_z_equals_xy();
        let witness = vec![Fr::one(), Fr::from(100u64), Fr::from(5u64), Fr::from(7u64)];
        assert!(Proof::new(&l, &r, &o, &witness, &crs).is_none());
    }

    #[test]
    fn test_verify_accepts_valid_proof() {
        let (l, r, o, crs) = setup_z_equals_xy();
        let witness = vec![Fr::one(), Fr::from(35u64), Fr::from(5u64), Fr::from(7u64)];
        let proof = Proof::new(&l, &r, &o, &witness, &crs).unwrap();
        let public_inputs = vec![Fr::one(), Fr::from(35u64)];
        assert!(proof.verify(&crs, &public_inputs));
    }

    #[test]
    fn test_verify_rejects_forged_proof() {
        let (_, _, _, crs) = setup_z_equals_xy();
        let public_inputs = vec![Fr::one(), Fr::from(35u64)];
        let forged_proof = Proof {
            a1: G1Affine::generator(),
            b2: G2Affine::from(G2Projective::generator() * Fr::from(2u64)),
            c1: G1Affine::from(G1Projective::generator() * Fr::from(2u64)),
        };
        assert!(!Proof::verify(&forged_proof, &crs, &public_inputs));
    }

    #[test]
    fn test_verify_rejects_wrong_public_inputs() {
        let (l, r, o, crs) = setup_z_equals_xy();
        let witness = vec![Fr::one(), Fr::from(35u64), Fr::from(5u64), Fr::from(7u64)];
        let proof = Proof::new(&l, &r, &o, &witness, &crs).unwrap();
        // Correct public inputs: [1, 35]. Wrong: [1, 36].
        assert!(!Proof::verify(&proof, &crs, &[Fr::one(), Fr::from(36u64)]));
    }

    #[test]
    fn test_verify_rejects_public_private_shift() {
        // Simulates a prover who claims output 0 but used 35 in private computation.
        // With gamma != delta, the equation cannot balance; verification fails.
        let (l, r, o, crs) = setup_z_equals_xy();
        let witness = vec![Fr::one(), Fr::from(35u64), Fr::from(5u64), Fr::from(7u64)];
        let proof = Proof::new(&l, &r, &o, &witness, &crs).unwrap();
        assert!(!Proof::verify(&proof, &crs, &[Fr::one(), Fr::zero()]));
    }

    #[test]
    fn test_proof_nondeterminism() {
        let (l, r, o, crs) = setup_z_equals_xy();
        let witness = vec![Fr::one(), Fr::from(35u64), Fr::from(5u64), Fr::from(7u64)];
        let public_inputs = vec![Fr::one(), Fr::from(35u64)];
        let proof1 = Proof::new(&l, &r, &o, &witness, &crs).unwrap();
        let proof2 = Proof::new(&l, &r, &o, &witness, &crs).unwrap();

        assert!(proof1.a1 != proof2.a1 || proof1.b2 != proof2.b2);
        assert!(Proof::verify(&proof1, &crs, &public_inputs));
        assert!(Proof::verify(&proof2, &crs, &public_inputs));
    }

    /// Helper: setup params for x^3 + x + 5 = out circuit (ell=2: constant 1 and out are public).
    fn setup_x_cubed_plus_x_plus_5() -> (Vec<Vec<Fr>>, Vec<Vec<Fr>>, Vec<Vec<Fr>>, CRS) {
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
        let crs = CRS::generate(
            &l, &r, &o,
            Fr::from(11u64), Fr::from(2u64), Fr::from(7u64),
            Fr::one(), Fr::from(19u64),
            2,
        );
        (l, r, o, crs)
    }

    #[test]
    fn test_prove_then_verify_end_to_end() {
        let (l, r, o, crs) = setup_x_cubed_plus_x_plus_5();
        let x = Fr::from(3u64);
        let x2 = x * x;
        let x3 = x2 * x;
        let x3_plus_x = x3 + x;
        let out = x3_plus_x + Fr::from(5u64);
        let witness = vec![Fr::one(), out, x, x2, x3, x3_plus_x];
        let public_inputs = vec![Fr::one(), out];
        let proof = Proof::new(&l, &r, &o, &witness, &crs).unwrap();
        assert!(Proof::verify(&proof, &crs, &public_inputs));
    }

    proptest! {
        #![proptest_config(Config { cases: 32, ..Config::default() })]

        #[test]
        fn prop_groth16_z_equals_xy(
            (x_u, y_u) in (1u64..1000, 1u64..1000)
        ) {
            let (l, r, o, crs) = setup_z_equals_xy();
            let x = Fr::from(x_u);
            let y = Fr::from(y_u);
            let z = x * y;
            let witness = vec![Fr::one(), z, x, y];
            let proof = Proof::new(&l, &r, &o, &witness, &crs).unwrap();

            let public_inputs = vec![Fr::one(), z];
            prop_assert!(proof.verify(&crs, &public_inputs));

            let wrong_public = vec![Fr::one(), z + Fr::one()];
            prop_assert!(!proof.verify(&crs, &wrong_public));
        }
    }

    proptest! {
        #![proptest_config(Config { cases: 32, ..Config::default() })]

        #[test]
        fn prop_groth16_x_cubed_plus_x_plus_5(
            x_u in (1u64..1000)
        ) {
            let (l, r, o, crs) = setup_x_cubed_plus_x_plus_5();
            let x = Fr::from(x_u);
            let x2 = x * x;
            let x3 = x2 * x;
            let x3_plus_x = x3 + x;
            let out = x3_plus_x + Fr::from(5u64);
            let witness = vec![Fr::one(), out, x, x2, x3, x3_plus_x];
            let proof = Proof::new(&l, &r, &o, &witness, &crs).unwrap();
            
            let public_inputs = vec![Fr::one(), out];
            prop_assert!(proof.verify(&crs, &public_inputs));

            let wrong_public = vec![Fr::one(), out + Fr::one()];
            prop_assert!(!proof.verify(&crs, &wrong_public));
        }
    }
}