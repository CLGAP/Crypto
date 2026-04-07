use zk_crypto::groth16::complete::{CRS, Proof};
use ark_bn254::Fr;
use ark_std::{One, Zero};

fn main() {
    println!("--- Groth16 Complete ---");

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
