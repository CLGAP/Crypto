//! Bilinear pairing verification using BN254.
//! Uses type conversion: `BigUint` to Fr
//!
//! Verifies the multi-pairing equation:
//!   0 = -e(A₁, B₂) + e(α₁, β₂) + e(X₁, γ₂) + e(C₁, δ₂)
//!
//! See: Theory document, Section 6 (Bilinear Pairings)

use crate::utils::biguint_to_fr;
use ark_bn254::{Bn254, Fr, G1Affine, G1Projective, G2Affine, G2Projective};
use ark_ec::pairing::{Pairing, PairingOutput};
use ark_ec::{CurveGroup, Group};
use ark_std::Zero;
use num_bigint::BigUint;

// Illustrative: alpha..delta are public here and the verifier derives
// alpha1 = alpha*G1 itself, so anyone can solve for a passing c. Real Groth16
// keeps these scalars secret in the trusted setup.
const ALPHA: u64 = 5;
const BETA: u64 = 6;
const GAMMA: u64 = 3;
const DELTA: u64 = 10;

#[must_use]
pub fn verify_pairing(
    a1: G1Affine,
    b2: G2Affine,
    c1: G1Affine,
    x1: &BigUint,
    x2: &BigUint,
    x3: &BigUint,
) -> bool {
    let g1 = G1Projective::generator();
    let g2 = G2Projective::generator();

    let alpha1 = (g1 * Fr::from(ALPHA)).into_affine();
    let beta2 = (g2 * Fr::from(BETA)).into_affine();
    let gamma2 = (g2 * Fr::from(GAMMA)).into_affine();
    let delta2 = (g2 * Fr::from(DELTA)).into_affine();

    // Could be distributed between parties, i.e., take in X = P_1 + P_2 + P_3 where P_i = x_i G1 for i in {1, 2, 3}
    let x_sum = x1 + x2 + x3;
    let x_sum = biguint_to_fr::<Fr>(&x_sum).unwrap();
    let x1_point = (g1 * x_sum).into_affine();

    let neg_a1 = -a1;
    let g1_points = [neg_a1, alpha1, x1_point, c1];
    let g2_points = [b2, beta2, gamma2, delta2];
    let product = Bn254::multi_pairing(g1_points.iter(), g2_points.iter());

    let identity = PairingOutput::<Bn254>::zero();

    product == identity
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_ff::{Field, PrimeField};
    use num_bigint::BigUint;
    use proptest::prelude::*;
    use proptest::test_runner::Config;

    #[test]
    fn test_verify_pairing_balanced_equation() {
        let g1 = G1Projective::generator();
        let g2 = G2Projective::generator();

        let a1 = (g1 * Fr::from(19u64)).into_affine();
        let b2 = (g2 * Fr::from(5u64)).into_affine();
        let c1 = (g1 * Fr::from(5u64)).into_affine();
        let x1 = BigUint::from(1u64);
        let x2 = BigUint::from(3u64);
        let x3 = BigUint::from(1u64);

        assert!(verify_pairing(a1, b2, c1, &x1, &x2, &x3));
    }

    #[test]
    fn test_verify_pairing_balanced_equation_with_overflow() {
        let g1 = G1Projective::generator();
        let g2 = G2Projective::generator();

        let a1 = (g1 * Fr::from(19u64)).into_affine();
        let b2 = (g2 * Fr::from(5u64)).into_affine();
        let c1 = (g1 * Fr::from(5u64)).into_affine();
        let x1 = BigUint::parse_bytes(&Fr::MODULUS.to_string().into_bytes(), 10).unwrap();
        let x2 = BigUint::from(3u64);
        let x3 = BigUint::from(2u64);

        assert!(verify_pairing(a1, b2, c1, &x1, &x2, &x3));
    }

    #[test]
    fn test_verify_pairing_fails_wrong_values() {
        let g1 = G1Projective::generator();
        let g2 = G2Projective::generator();

        let a1 = (g1 * Fr::from(1u64)).into_affine();
        let b2 = (g2 * Fr::from(1u64)).into_affine();
        let c1 = (g1 * Fr::from(1u64)).into_affine();
        let x1 = BigUint::from(1u64);
        let x2 = BigUint::from(1u64);
        let x3 = BigUint::from(1u64);

        assert!(!verify_pairing(a1, b2, c1, &x1, &x2, &x3));
    }

    proptest! {
        #![proptest_config(Config {
            cases: 32, // keep moderate: pairings are expensive
            ..Config::default()
        })]

        #[test]
        fn prop_verify_pairing_accepts_constructed_valid_inputs(
            (a_u, b_u, x1_u, x2_u, x3_u) in
                (1u64..1000, 1u64..1000, 0u64..1000, 0u64..1000, 0u64..1000)
        ) {
            let g1 = G1Projective::generator();
            let g2 = G2Projective::generator();

            let a = Fr::from(a_u);
            let b = Fr::from(b_u);
            let x1 = Fr::from(x1_u);
            let x2 = Fr::from(x2_u);
            let x3 = Fr::from(x3_u);
            let x_sum = x1 + x2 + x3;

            // Constants from module
            let alpha = Fr::from(ALPHA);
            let beta = Fr::from(BETA);
            let gamma = Fr::from(GAMMA);
            let delta = Fr::from(DELTA);

            // Solve: delta*c = a*b - alpha*beta - gamma*x_sum
            let numerator = a * b - alpha * beta - gamma * x_sum;
            let c = numerator * delta.inverse().unwrap();

            let a1 = (g1 * a).into_affine();
            let b2 = (g2 * b).into_affine();
            let c1 = (g1 * c).into_affine();

            let x1_bi = num_bigint::BigUint::from(x1_u);
            let x2_bi = num_bigint::BigUint::from(x2_u);
            let x3_bi = num_bigint::BigUint::from(x3_u);

            prop_assert!(verify_pairing(a1, b2, c1, &x1_bi, &x2_bi, &x3_bi));
        }
    }

    proptest! {
        #![proptest_config(Config { cases: 32, ..Config::default() })]

        #[test]
        fn prop_verify_pairing_rejects_tampered_input(
            (a_u, b_u, x1_u, x2_u, x3_u) in
                (1u64..1000, 1u64..1000, 0u64..1000, 0u64..1000, 0u64..1000)
        ) {
            let g1 = G1Projective::generator();
            let g2 = G2Projective::generator();

            let a = Fr::from(a_u);
            let b = Fr::from(b_u);
            let x1 = Fr::from(x1_u);
            let x2 = Fr::from(x2_u);
            let x3 = Fr::from(x3_u);

            let alpha = Fr::from(ALPHA);
            let beta = Fr::from(BETA);
            let gamma = Fr::from(GAMMA);
            let delta = Fr::from(DELTA);

            let numerator = a * b - alpha * beta - gamma * (x1 + x2 + x3);
            let c = numerator * delta.inverse().unwrap();

            let a1 = (g1 * a).into_affine();
            let b2 = (g2 * b).into_affine();
            let c1 = (g1 * c).into_affine();

            let x1_bi = num_bigint::BigUint::from(x1_u);
            let x2_bi = num_bigint::BigUint::from(x2_u);
            let x3_tampered = num_bigint::BigUint::from(x3_u + 1);

            prop_assert!(!verify_pairing(a1, b2, c1, &x1_bi, &x2_bi, &x3_tampered));
        }
    }

    proptest! {
        #![proptest_config(Config { cases: 32, ..Config::default() })]

        #[test]
        fn prop_pairing_is_bilinear(
            (a1_u, a2_u, b_u, c_u, d_u) in
                (0u64..1000, 0u64..1000, 0u64..1000, 0u64..1000, 0u64..1000)
        ) {
            let g1 = G1Projective::generator();
            let g2 = G2Projective::generator();

            let p1 = g1 * Fr::from(a1_u);
            let p2 = g1 * Fr::from(a2_u);
            let q = g2 * Fr::from(b_u);
            let r1 = g2 * Fr::from(c_u);
            let r2 = g2 * Fr::from(d_u);

            // Bilinear in first argument: e(P1 + P2, Q) = e(P1, Q) * e(P2, Q)
            // PairingOutput uses additive notation for the target-group law.
            let lhs_first = Bn254::pairing((p1 + p2).into_affine(), q.into_affine());
            let rhs_first = Bn254::pairing(p1.into_affine(), q.into_affine())
                + Bn254::pairing(p2.into_affine(), q.into_affine());
            prop_assert_eq!(lhs_first, rhs_first);

            // Bilinear in second argument: e(P, R1 + R2) = e(P, R1) * e(P, R2)
            let lhs_second = Bn254::pairing(p1.into_affine(), (r1 + r2).into_affine());
            let rhs_second = Bn254::pairing(p1.into_affine(), r1.into_affine())
                + Bn254::pairing(p1.into_affine(), r2.into_affine());
            prop_assert_eq!(lhs_second, rhs_second);
        }
    }

    proptest! {
        #![proptest_config(Config { cases: 32, ..Config::default() })]

        #[test]
        fn prop_pairing_is_nondegenerate_for_nonzero_generator_multiples(
            (a_u, b_u) in (1u64..1000, 1u64..1000)
        ) {
            let g1 = G1Projective::generator();
            let g2 = G2Projective::generator();
            let p = (g1 * Fr::from(a_u)).into_affine();
            let q = (g2 * Fr::from(b_u)).into_affine();

            let identity = PairingOutput::<Bn254>::zero();
            let value = Bn254::pairing(p, q);
            prop_assert_ne!(value, identity);
        }
    }

    proptest! {
        // Demo-only property to showcase shrinking behaviour.
        // We intentionally assert a false property for part of the domain:
        // x + y < 199000 does not always hold for x,y in 0..100000.
        #![proptest_config(proptest::test_runner::Config {
            cases: 10000,
            ..proptest::test_runner::Config::default()
        })]
        #[test]
        #[ignore = "proptest demo-only: intentionally fails so proptest can show shrinking"]
        fn demo_proptest((x, y) in (0u64..100_000, 0u64..100_000)) {
            prop_assert!(
                x + y < 199_000,
                "expected x + y < 199_000, got x={x}, y={y}, sum={}",
                x + y
            );
        }
    }
}
