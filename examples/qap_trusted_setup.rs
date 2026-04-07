use zk_crypto::qap_trusted_setup::{prove, verify};
use ark_bn254::Fr;
use ark_std::{One, Zero};

fn main() {
    println!("--- QAP Trusted Setup (Pinocchio) ---");

    let x = Fr::from(5u64);
    let y = Fr::from(7u64);
    let z = x * y;
    println!("  z = x * y  (x={}, y={}, z={})", x, y, z);

    let l = vec![vec![Fr::zero(), Fr::zero(), Fr::one(), Fr::zero()]];
    let r = vec![vec![Fr::zero(), Fr::zero(), Fr::zero(), Fr::one()]];
    let o = vec![vec![Fr::zero(), Fr::one(), Fr::zero(), Fr::zero()]];
    let witness = vec![Fr::one(), z, x, y];
    let tau = Fr::from(4u64);

    let proof = prove(&l, &r, &o, &witness, tau).expect("valid witness should produce proof");
    let result = verify(&proof);
    println!("  Proof generated, verification: {}", if result { "PASS" } else { "FAIL" });
}
