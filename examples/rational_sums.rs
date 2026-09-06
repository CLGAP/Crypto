use zk_crypto::rational_sums::verify_rational_sum;
use num_bigint::BigUint;

fn main() {
    println!("--- Rational Sums ---");
    let result = verify_rational_sum(
        BigUint::from(2u64), BigUint::from(3u64),
        BigUint::from(5u64), BigUint::from(7u64),
        BigUint::from(29u64), BigUint::from(21u64),
    );
    match result {
        Ok(true) => println!("  2/3 + 5/7 = 29/21: verified"),
        Ok(false) => println!("  2/3 + 5/7 = 29/21: failed"),
        Err(e) => println!("  Error: {e}"),
    }
}
