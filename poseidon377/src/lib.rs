//! An instantiation of Poseidon for the BLS12-377 scalar field.
#![cfg_attr(not(feature = "std"), no_std)]

mod hash;
mod params;

pub use hash::{hash_1, hash_2, hash_3, hash_4, hash_5, hash_6, hash_7};

/// Parameters for the rate-1 instance of Poseidon.
pub const RATE_1_PARAMS: PoseidonParameters<2, 1, 4, 1, 64, 2, 128, 56> = params::rate_1::rate_1();

/// Parameters for the rate-2 instance of Poseidon.
pub const RATE_2_PARAMS: PoseidonParameters<3, 2, 9, 4, 65, 3, 195, 57> = params::rate_2::rate_2();

/// Parameters for the rate-3 instance of Poseidon.
pub const RATE_3_PARAMS: PoseidonParameters<4, 3, 16, 9, 64, 4, 256, 56> = params::rate_3::rate_3();

/// Parameters for the rate-4 instance of Poseidon.
pub const RATE_4_PARAMS: PoseidonParameters<5, 4, 25, 16, 68, 5, 340, 60> =
    params::rate_4::rate_4();

/// Parameters for the rate-5 instance of Poseidon.
pub const RATE_5_PARAMS: PoseidonParameters<6, 5, 36, 25, 68, 6, 408, 60> =
    params::rate_5::rate_5();

/// Parameters for the rate-6 instance of Poseidon.
pub const RATE_6_PARAMS: PoseidonParameters<7, 6, 49, 36, 71, 7, 497, 63> =
    params::rate_6::rate_6();

/// Parameters for the rate-7 instance of Poseidon.
pub const RATE_7_PARAMS: PoseidonParameters<8, 7, 64, 49, 39, 8, 312, 31> =
    params::rate_7::rate_7();

pub use cycles_curve_bn254::Fq;
pub use poseidon_parameters::v1::PoseidonParameters;
pub use poseidon_permutation::Instance;

#[cfg(feature = "r1cs")]
pub mod r1cs;
