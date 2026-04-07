use zk_crypto::r1cs_to_qap::{
    lagrange_interpolate, evaluate_polynomial, construct_qap,
    target_polynomial, verify_qap_satisfaction,
};
use ark_bn254::Fr;
use ark_std::{One, Zero, rand::{rngs::StdRng, SeedableRng}};

fn main() {
    println!("--- R1CS to QAP ---");

    // Lagrange interpolation: points (1,2), (2,4), (3,6) => linear polynomial
    let points = vec![Fr::from(1u64), Fr::from(2u64), Fr::from(3u64)];
    let values = vec![Fr::from(2u64), Fr::from(4u64), Fr::from(6u64)];
    let poly = lagrange_interpolate(&points, &values);
    println!("  Interpolation at x=2: {} (expected 4)", evaluate_polynomial(&poly, Fr::from(2u64)));

    // QAP construction from z = x * y
    let l = vec![vec![Fr::zero(), Fr::zero(), Fr::one(), Fr::zero()]];
    let r = vec![vec![Fr::zero(), Fr::zero(), Fr::zero(), Fr::one()]];
    let o = vec![vec![Fr::zero(), Fr::one(), Fr::zero(), Fr::zero()]];
    let eval_points = vec![Fr::from(1u64)];
    let (l_polys, r_polys, o_polys) = construct_qap(&l, &r, &o, &eval_points);
    println!("  QAP columns: L={}, R={}, O={}", l_polys.len(), r_polys.len(), o_polys.len());

    // QAP satisfaction check
    let t_poly = target_polynomial(&eval_points);
    let witness = vec![Fr::one(), Fr::from(35u64), Fr::from(5u64), Fr::from(7u64)];
    let mut rng = StdRng::seed_from_u64(42);
    let satisfied = verify_qap_satisfaction(&l_polys, &r_polys, &o_polys, &witness, &t_poly, &mut rng);
    println!("  QAP satisfied (x=5, y=7, z=35): {}", if satisfied { "yes" } else { "no" });
}
