use crate::{Fq, Instance};

/// Hash a single [`Fq`] element with the provided `domain_separator`.
pub fn hash_1(domain_separator: &Fq, value: Fq) -> Fq {
    let params = &crate::RATE_1_PARAMS;
    let mut state = Instance::new(&params);
    state.n_to_1_fixed_hash(&[*domain_separator, value])
}

/// Hash two [`Fq`] elements with the provided `domain_separator`.
pub fn hash_2(domain_separator: &Fq, value: (Fq, Fq)) -> Fq {
    let params = &crate::RATE_2_PARAMS;
    let mut state = Instance::new(params);
    state.n_to_1_fixed_hash(&[*domain_separator, value.0, value.1])
}

/// Hash three [`Fq`] elements with the provided `domain_separator`.
pub fn hash_3(domain_separator: &Fq, value: (Fq, Fq, Fq)) -> Fq {
    let params = &crate::RATE_3_PARAMS;
    let mut state = Instance::new(params);
    state.n_to_1_fixed_hash(&[*domain_separator, value.0, value.1, value.2])
}

/// Hash four [`Fq`] elements with the provided `domain_separator`.
pub fn hash_4(domain_separator: &Fq, value: (Fq, Fq, Fq, Fq)) -> Fq {
    let params = &crate::RATE_4_PARAMS;
    let mut state = Instance::new(params);
    state.n_to_1_fixed_hash(&[*domain_separator, value.0, value.1, value.2, value.3])
}

/// Hash five [`Fq`] elements with the provided `domain_separator`.
pub fn hash_5(domain_separator: &Fq, value: (Fq, Fq, Fq, Fq, Fq)) -> Fq {
    let params = &crate::RATE_5_PARAMS;
    let mut state = Instance::new(params);
    state.n_to_1_fixed_hash(&[
        *domain_separator,
        value.0, value.1, value.2, value.3, value.4,
    ])
}

/// Hash six [`Fq`] elements with the provided `domain_separator`.
pub fn hash_6(domain_separator: &Fq, value: (Fq, Fq, Fq, Fq, Fq, Fq)) -> Fq {
    let params = &crate::RATE_6_PARAMS;
    let mut state = Instance::new(params);
    state.n_to_1_fixed_hash(&[
        *domain_separator,
        value.0, value.1, value.2, value.3, value.4, value.5,
    ])
}

/// Hash seven [`Fq`] elements with the provided `domain_separator`.
pub fn hash_7(domain_separator: &Fq, value: (Fq, Fq, Fq, Fq, Fq, Fq, Fq)) -> Fq {
    let params = &crate::RATE_7_PARAMS;
    let mut state = Instance::new(params);
    state.n_to_1_fixed_hash(&[
        *domain_separator,
        value.0, value.1, value.2, value.3, value.4, value.5, value.6,
    ])
}

#[cfg(test)]
mod test {
    use core::str::FromStr;
    use super::*;

    /// All test vectors computed with iden3's circomlibjs poseidon_opt.js
    /// using domain_separator = 0 and inputs = [1, 2, ..., N].
    fn d() -> Fq { Fq::from(0u64) }

    #[test]
    fn test_hash_1_iden3() {
        let expected = Fq::from_str(
            "18586133768512220936620570745912940619677854269274689475585506675881198879027",
        ).unwrap();
        assert_eq!(hash_1(&d(), Fq::from(1u64)), expected);
    }

    #[test]
    fn test_hash_2_iden3() {
        let expected = Fq::from_str(
            "7853200120776062878684798364095072458815029376092732009249414926327459813530",
        ).unwrap();
        assert_eq!(hash_2(&d(), (Fq::from(1u64), Fq::from(2u64))), expected);
    }

    #[test]
    fn test_hash_3_iden3() {
        let expected = Fq::from_str(
            "6542985608222806190361240322586112750744169038454362455181422643027100751666",
        ).unwrap();
        assert_eq!(
            hash_3(&d(), (Fq::from(1u64), Fq::from(2u64), Fq::from(3u64))),
            expected,
        );
    }

    #[test]
    fn test_hash_4_iden3() {
        let expected = Fq::from_str(
            "18821383157269793795438455681495246036402687001665670618754263018637548127333",
        ).unwrap();
        assert_eq!(
            hash_4(&d(), (Fq::from(1u64), Fq::from(2u64), Fq::from(3u64), Fq::from(4u64))),
            expected,
        );
    }

    #[test]
    fn test_hash_5_iden3() {
        let expected = Fq::from_str(
            "6183221330272524995739186171720101788151706631170188140075976616310159254464",
        ).unwrap();
        assert_eq!(
            hash_5(&d(), (
                Fq::from(1u64), Fq::from(2u64), Fq::from(3u64), Fq::from(4u64), Fq::from(5u64),
            )),
            expected,
        );
    }

    #[test]
    fn test_hash_6_iden3() {
        let expected = Fq::from_str(
            "20400040500897583745843009878988256314335038853985262692600694741116813247201",
        ).unwrap();
        assert_eq!(
            hash_6(&d(), (
                Fq::from(1u64), Fq::from(2u64), Fq::from(3u64),
                Fq::from(4u64), Fq::from(5u64), Fq::from(6u64),
            )),
            expected,
        );
    }

    #[test]
    fn test_hash_7_iden3() {
        let expected = Fq::from_str(
            "12748163991115452309045839028154629052133952896122405799815156419278439301912",
        ).unwrap();
        assert_eq!(
            hash_7(&d(), (
                Fq::from(1u64), Fq::from(2u64), Fq::from(3u64), Fq::from(4u64),
                Fq::from(5u64), Fq::from(6u64), Fq::from(7u64),
            )),
            expected,
        );
    }
}

// #[cfg(test)]
// mod test_penumbra {
//     use core::str::FromStr;

//     use super::*;

//     #[test]
//     fn rate_1() {
//         let domain_sep = Fq::from_le_bytes_mod_order(b"Penumbra_TestVec");

//         let input = Fq::from_str(
//             "7553885614632219548127688026174585776320152166623257619763178041781456016062",
//         )
//         .unwrap();
//         let output = hash_1(&domain_sep, input);

//         let expected_output = Fq::from_str(
//             "2337838243217876174544784248400816541933405738836087430664765452605435675740",
//         )
//         .unwrap();

//         assert_eq!(output, expected_output);
//     }

//     #[test]
//     fn rate_2() {
//         let domain_sep = Fq::from_le_bytes_mod_order(b"Penumbra_TestVec");

//         let input = (
//             Fq::from_str(
//                 "7553885614632219548127688026174585776320152166623257619763178041781456016062",
//             )
//             .unwrap(),
//             Fq::from_str(
//                 "2337838243217876174544784248400816541933405738836087430664765452605435675740",
//             )
//             .unwrap(),
//         );

//         let output = hash_2(&domain_sep, input);

//         let expected_output = Fq::from_str(
//             "4318449279293553393006719276941638490334729643330833590842693275258805886300",
//         )
//         .unwrap();

//         assert_eq!(output, expected_output);
//     }

//     #[test]
//     fn rate_3() {
//         let domain_sep = Fq::from_le_bytes_mod_order(b"Penumbra_TestVec");

//         let input = (
//             Fq::from_str(
//                 "7553885614632219548127688026174585776320152166623257619763178041781456016062",
//             )
//             .unwrap(),
//             Fq::from_str(
//                 "2337838243217876174544784248400816541933405738836087430664765452605435675740",
//             )
//             .unwrap(),
//             Fq::from_str(
//                 "4318449279293553393006719276941638490334729643330833590842693275258805886300",
//             )
//             .unwrap(),
//         );

//         let output = hash_3(&domain_sep, input);

//         let expected_output = Fq::from_str(
//             "2884734248868891876687246055367204388444877057000108043377667455104051576315",
//         )
//         .unwrap();

//         assert_eq!(output, expected_output);
//     }

//     #[test]
//     fn rate_4() {
//         let domain_sep = Fq::from_le_bytes_mod_order(b"Penumbra_TestVec");

//         let input = (
//             Fq::from_str(
//                 "7553885614632219548127688026174585776320152166623257619763178041781456016062",
//             )
//             .unwrap(),
//             Fq::from_str(
//                 "2337838243217876174544784248400816541933405738836087430664765452605435675740",
//             )
//             .unwrap(),
//             Fq::from_str(
//                 "4318449279293553393006719276941638490334729643330833590842693275258805886300",
//             )
//             .unwrap(),
//             Fq::from_str(
//                 "2884734248868891876687246055367204388444877057000108043377667455104051576315",
//             )
//             .unwrap(),
//         );

//         let output = hash_4(&domain_sep, input);

//         let expected_output = Fq::from_str(
//             "5235431038142849831913898188189800916077016298531443239266169457588889298166",
//         )
//         .unwrap();

//         assert_eq!(output, expected_output);
//     }

//     #[test]
//     fn rate_5() {
//         let domain_sep = Fq::from_le_bytes_mod_order(b"Penumbra_TestVec");

//         let input = (
//             Fq::from_str(
//                 "7553885614632219548127688026174585776320152166623257619763178041781456016062",
//             )
//             .unwrap(),
//             Fq::from_str(
//                 "2337838243217876174544784248400816541933405738836087430664765452605435675740",
//             )
//             .unwrap(),
//             Fq::from_str(
//                 "4318449279293553393006719276941638490334729643330833590842693275258805886300",
//             )
//             .unwrap(),
//             Fq::from_str(
//                 "2884734248868891876687246055367204388444877057000108043377667455104051576315",
//             )
//             .unwrap(),
//             Fq::from_str(
//                 "5235431038142849831913898188189800916077016298531443239266169457588889298166",
//             )
//             .unwrap(),
//         );

//         let output = hash_5(&domain_sep, input);

//         let expected_output = Fq::from_str(
//             "66948599770858083122195578203282720327054804952637730715402418442993895152",
//         )
//         .unwrap();

//         assert_eq!(output, expected_output);
//     }

//     #[test]
//     fn rate_6() {
//         let domain_sep = Fq::from_le_bytes_mod_order(b"Penumbra_TestVec");

//         let input = (
//             Fq::from_str(
//                 "7553885614632219548127688026174585776320152166623257619763178041781456016062",
//             )
//             .unwrap(),
//             Fq::from_str(
//                 "2337838243217876174544784248400816541933405738836087430664765452605435675740",
//             )
//             .unwrap(),
//             Fq::from_str(
//                 "4318449279293553393006719276941638490334729643330833590842693275258805886300",
//             )
//             .unwrap(),
//             Fq::from_str(
//                 "2884734248868891876687246055367204388444877057000108043377667455104051576315",
//             )
//             .unwrap(),
//             Fq::from_str(
//                 "5235431038142849831913898188189800916077016298531443239266169457588889298166",
//             )
//             .unwrap(),
//             Fq::from_str(
//                 "66948599770858083122195578203282720327054804952637730715402418442993895152",
//             )
//             .unwrap(),
//         );

//         let output = hash_6(&domain_sep, input);

//         let expected_output = Fq::from_str(
//             "6797655301930638258044003960605211404784492298673033525596396177265014216269",
//         )
//         .unwrap();

//         assert_eq!(output, expected_output);
//     }
// }
