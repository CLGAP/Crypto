use zk_crypto::r1cs::{verify_r1cs, verify_r1cs_with_hadamard_product, create_witness_points};
use ark_bn254::Fr;
use ark_std::{One, Zero};

fn main() {
    println!("--- R1CS ---");
    let x = Fr::from(5u64);
    let y = Fr::from(7u64);
    let z = x * y;
    println!("  z = x * y  (x={x}, y={y}, z={z})");

    let l = vec![vec![Fr::zero(), Fr::zero(), Fr::one(), Fr::zero()]];
    let r = vec![vec![Fr::zero(), Fr::zero(), Fr::zero(), Fr::one()]];
    let o = vec![vec![Fr::zero(), Fr::one(), Fr::zero(), Fr::zero()]];
    let witness = vec![Fr::one(), z, x, y];

    let field_result = verify_r1cs_with_hadamard_product(&l, &r, &o, &witness);
    println!("  Field-based (Hadamard): {}", if field_result { "verified" } else { "failed" });

    let (w_g1, w_g2) = create_witness_points(&witness);
    let pairing_result = verify_r1cs(&l, &r, &o, &w_g1, &w_g2);
    println!("  Pairing-based: {}", if pairing_result { "verified" } else { "failed" });

    let (w_g1_bad, w_g2_bad) = create_witness_points(&[Fr::one(), Fr::from(100u64), x, y]);
    let bad_result = verify_r1cs(&l, &r, &o, &w_g1_bad, &w_g2_bad);
    println!("  Wrong witness (expect reject): {}", if bad_result { "incorrectly accepted" } else { "correctly rejected" });
}
