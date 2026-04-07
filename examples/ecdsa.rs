use zk_crypto::ecdsa::{KeyPair, hash_message, sign, verify};

fn main() {
    println!("--- ECDSA ---");
    let key_pair = KeyPair::generate();
    let message = b"Raspberry Pi";
    let h = hash_message(message);

    println!("  Signing: \"{}\"", String::from_utf8_lossy(message));
    let signature = sign(&key_pair.private_key, &h, true);
    println!("  r = {:?}", signature.r);
    println!("  s = {:?}", signature.s);

    let valid = verify(&signature, message, &key_pair.public_key);
    println!("  Verification: {}", if valid { "valid" } else { "invalid" });
}
