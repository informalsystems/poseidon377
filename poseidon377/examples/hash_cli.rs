use core::str::FromStr;

use poseidon377::{hash_2, Fq};

fn fq_to_decimal(fq: Fq) -> String {
    let bytes = fq.to_bytes_le();
    let mut digits = vec![0u8];
    for &byte in bytes.iter().rev() {
        let mut carry = byte as u32;
        for digit in digits.iter_mut() {
            let val = (*digit as u32) * 256 + carry;
            *digit = (val % 10) as u8;
            carry = val / 10;
        }
        while carry > 0 {
            digits.push((carry % 10) as u8);
            carry /= 10;
        }
    }
    let s: String = digits.iter().rev().map(|d| (b'0' + d) as char).collect();
    let s = s.trim_start_matches('0');
    if s.is_empty() { "0".to_string() } else { s.to_string() }
}

fn parse_arg(name: &str, s: &str) -> Fq {
    Fq::from_str(s).unwrap_or_else(|_| {
        eprintln!("invalid {}: {}", name, s);
        std::process::exit(2);
    })
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 4 {
        eprintln!("usage: hash_cli <domain_decimal> <value1_decimal> <value2_decimal>");
        std::process::exit(2);
    }

    let domain = parse_arg("domain", &args[1]);
    let v1 = parse_arg("value1", &args[2]);
    let v2 = parse_arg("value2", &args[3]);

    let out = hash_2(&domain, (v1, v2));
    println!("{}", fq_to_decimal(out));
}
