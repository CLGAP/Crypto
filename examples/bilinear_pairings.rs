use zk_crypto::bilinear_pairings::verify_pairing;
use ark_bn254::{Fr, G1Projective, G2Projective};
use ark_ec::Group;
use num_bigint::BigUint;

fn main() {
    println!("--- Bilinear Pairings ---");
    let g1 = G1Projective::generator();
    let g2 = G2Projective::generator();

    let a1 = (g1 * Fr::from(19u64)).into();
    let b2 = (g2 * Fr::from(5u64)).into();
    let c1 = (g1 * Fr::from(5u64)).into();

    let result = verify_pairing(a1, b2, c1, &BigUint::from(1u64), &BigUint::from(3u64), &BigUint::from(1u64));
    println!("  Balanced equation: {}", if result { "verified" } else { "failed" });

    let wrong = (g1 * Fr::from(1u64)).into();
    let result_wrong = verify_pairing(wrong, (g2 * Fr::from(1u64)).into(), wrong, &BigUint::from(1u64), &BigUint::from(1u64), &BigUint::from(1u64));
    println!("  Wrong values (expect reject): {}", if result_wrong { "incorrectly accepted" } else { "correctly rejected" });
}
