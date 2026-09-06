use zk_crypto::ecdsa::{KeyPair, hash_message, sign, verify as ecdsa_verify};
use zk_crypto::rational_sums::verify_rational_sum;
use zk_crypto::bilinear_pairings::verify_pairing;
use zk_crypto::r1cs::{verify_r1cs, create_witness_points};
use zk_crypto::ecp_matrix::verify_matrix_ec_multiplication;
use zk_crypto::groth16::complete::{CRS, Proof};
use num_bigint::BigUint;
use ark_bn254::Fr;
use ark_std::{One, Zero};

fn main() {
    println!("--- Rational Sums ---");
    run_rational_sums();
    println!();

    println!("--- ECP Matrix ---");
    run_ecp_matrix();
    println!();

    println!("--- ECDSA ---");
    run_ecdsa();
    println!();

    println!("--- Bilinear Pairings ---");
    run_bilinear_pairings();
    println!();

    println!("--- R1CS ---");
    run_r1cs();
    println!();

    println!("--- Groth16 ---");
    run_groth16();
}

fn run_rational_sums() {
    let result = verify_rational_sum(
        BigUint::from(2u64),
        BigUint::from(3u64),
        BigUint::from(5u64),
        BigUint::from(7u64),
        BigUint::from(29u64),
        BigUint::from(21u64),
    );
    match result {
        Ok(true) => println!("  2/3 + 5/7 = 29/21: verified"),
        Ok(false) => println!("  2/3 + 5/7 = 29/21: failed"),
        Err(e) => println!("  Error: {e}"),
    }
}

fn run_ecp_matrix() {
    let result_2x2 = verify_matrix_ec_multiplication(2, 100);
    println!("  2x2 matrix: {}", if result_2x2 { "verified" } else { "failed" });

    let result_4x4 = verify_matrix_ec_multiplication(4, 100);
    println!("  4x4 matrix: {}", if result_4x4 { "verified" } else { "failed" });
}

fn run_ecdsa() {
    let key_pair = KeyPair::generate();
    let message = b"Raspberry Pi";
    let h = hash_message(message);

    println!("  Signing: \"{}\"", String::from_utf8_lossy(message));
    let signature = sign(&key_pair.private_key, &h);
    println!("  r = {:?}", signature.r);
    println!("  s = {:?}", signature.s);

    let valid = ecdsa_verify(&signature, message, &key_pair.public_key);
    println!("  Verification: {}", if valid { "valid" } else { "invalid" });
}

fn run_bilinear_pairings() {
    use ark_bn254::{G1Projective, G2Projective};
    use ark_ec::Group;

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

fn run_r1cs() {
    let x = Fr::from(5u64);
    let y = Fr::from(7u64);
    let z = x * y;
    println!("  z = x * y  (x={x}, y={y}, z={z})");

    let l = vec![vec![Fr::zero(), Fr::zero(), Fr::one(), Fr::zero()]];
    let r = vec![vec![Fr::zero(), Fr::zero(), Fr::zero(), Fr::one()]];
    let o = vec![vec![Fr::zero(), Fr::one(), Fr::zero(), Fr::zero()]];

    let (w_g1, w_g2) = create_witness_points(&[Fr::one(), z, x, y]);
    let result = verify_r1cs(&l, &r, &o, &w_g1, &w_g2);
    println!("  Valid witness: {}", if result { "verified" } else { "failed" });

    let (w_g1_bad, w_g2_bad) = create_witness_points(&[Fr::one(), Fr::from(100u64), x, y]);
    let result_bad = verify_r1cs(&l, &r, &o, &w_g1_bad, &w_g2_bad);
    println!("  Wrong witness (expect reject): {}", if result_bad { "incorrectly accepted" } else { "correctly rejected" });
}

fn run_groth16() {
    let l = vec![vec![Fr::zero(), Fr::zero(), Fr::one(), Fr::zero()]];
    let r = vec![vec![Fr::zero(), Fr::zero(), Fr::zero(), Fr::one()]];
    let o = vec![vec![Fr::zero(), Fr::one(), Fr::zero(), Fr::zero()]];

    let crs = CRS::generate(&l, &r, &o,
        Fr::from(11u64), Fr::from(13u64), Fr::from(4u64),
        Fr::one(), Fr::from(19u64), 2);

    let witness = vec![Fr::one(), Fr::from(35u64), Fr::from(5u64), Fr::from(7u64)];
    let proof = Proof::new(&l, &r, &o, &witness, &crs).unwrap();
    let result = proof.verify(&crs, &[Fr::one(), Fr::from(35u64)]);
    println!("  z = x * y  (x=5, y=7, z=35)");
    println!("  Proof with random blinding (r, s)");
    println!("  Verification: {}", if result { "PASS" } else { "FAIL" });
}
