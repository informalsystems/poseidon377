use core::str::FromStr;

use poseidon377::{hash_1, Fq};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 3 {
        eprintln!("usage: hash_cli <domain_separator_decimal> <value_decimal>");
        std::process::exit(2);
    }

    let domain_separator = match Fq::from_str(&args[1]) {
        Ok(v) => v,
        Err(_) => {
            eprintln!("invalid domain separator: {}", args[1]);
            std::process::exit(2);
        }
    };

    let value = match Fq::from_str(&args[2]) {
        Ok(v) => v,
        Err(_) => {
            eprintln!("invalid value: {}", args[2]);
            std::process::exit(2);
        }
    };

    let out = hash_1(&domain_separator, value);
    println!("{out:?}");
}
