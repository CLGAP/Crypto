use zk_crypto::groth16::alpha_beta::{trusted_setup, prove, verify};
use ark_bn254::Fr;
use ark_std::{One, Zero};

fn main() {
    println!("--- Groth16 Alpha-Beta ---");

    let l = vec![vec![Fr::zero(), Fr::zero(), Fr::one(), Fr::zero()]];
    let r = vec![vec![Fr::zero(), Fr::zero(), Fr::zero(), Fr::one()]];
    let o = vec![vec![Fr::zero(), Fr::one(), Fr::zero(), Fr::zero()]];

    let crs = trusted_setup(&l, &r, &o, Fr::from(11u64), Fr::from(13u64), Fr::from(4u64));
    let witness = vec![Fr::one(), Fr::from(35u64), Fr::from(5u64), Fr::from(7u64)];
    let proof = prove(&l, &r, &o, &witness, &crs).unwrap();
    let result = verify(&proof, &crs);
    println!("  z = x * y  (x=5, y=7, z=35)");
    println!("  Verification: {}", if result { "PASS" } else { "FAIL" });
}
