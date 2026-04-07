use zk_crypto::ecp_matrix::verify_matrix_ec_multiplication;

fn main() {
    println!("--- ECP Matrix ---");
    let result_2x2 = verify_matrix_ec_multiplication(2, 100);
    println!("  2x2 matrix: {}", if result_2x2 { "verified" } else { "failed" });

    let result_4x4 = verify_matrix_ec_multiplication(4, 100);
    println!("  4x4 matrix: {}", if result_4x4 { "verified" } else { "failed" });
}
