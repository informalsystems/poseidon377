use ark_ed_on_bls12_381::{Fq, FqConfig};  // Replace with your field
use ark_ff::MontConfig;
use poseidon_paramgen::v1;

fn main() {
    let params = v1::generate::<Fq>(
        128,                  // 128-bit security
        3,                    // t=3 (2-to-1 hash)
        FqConfig::MODULUS,
        true,
    );
    
    println!("Generated Poseidon parameters:");
    println!("  Security level: {}", params.M);
    println!("  Width (t): {}", params.t);
    println!("  Alpha: {:?}", params.alpha);
    println!("  Full rounds: {}", params.rounds.full());
    println!("  Partial rounds: {}", params.rounds.partial());
}