use cycles_curve_bn254::Fq;
use std::fmt::Write as FmtWrite;
use serde_json::Value;

// ---------------------------------------------------------------------------
// Tiny dynamic-matrix library over Fq (avoids const-generic pain in codegen)
// ---------------------------------------------------------------------------

#[derive(Clone, Debug)]
struct Mat {
    rows: usize,
    cols: usize,
    data: Vec<Fq>,
}

impl Mat {
    fn new(rows: usize, cols: usize, data: Vec<Fq>) -> Self {
        assert_eq!(data.len(), rows * cols);
        Self { rows, cols, data }
    }

    fn zeros(rows: usize, cols: usize) -> Self {
        Self::new(rows, cols, vec![Fq::from(0u64); rows * cols])
    }

    fn identity(n: usize) -> Self {
        let mut m = Self::zeros(n, n);
        for i in 0..n {
            m.set(i, i, Fq::from(1u64));
        }
        m
    }

    fn get(&self, i: usize, j: usize) -> Fq {
        self.data[i * self.cols + j]
    }

    fn set(&mut self, i: usize, j: usize, v: Fq) {
        self.data[i * self.cols + j] = v;
    }

    fn transpose(&self) -> Self {
        let mut out = Self::zeros(self.cols, self.rows);
        for i in 0..self.rows {
            for j in 0..self.cols {
                out.set(j, i, self.get(i, j));
            }
        }
        out
    }

    fn mul(&self, rhs: &Self) -> Self {
        assert_eq!(self.cols, rhs.rows);
        let mut out = Self::zeros(self.rows, rhs.cols);
        for i in 0..self.rows {
            for j in 0..rhs.cols {
                let mut s = Fq::from(0u64);
                for k in 0..self.cols {
                    s += self.get(i, k) * rhs.get(k, j);
                }
                out.set(i, j, s);
            }
        }
        out
    }

    fn sub_matrix(&self, skip_row: usize, skip_col: usize) -> Self {
        let mut elems = Vec::with_capacity((self.rows - 1) * (self.cols - 1));
        for i in 0..self.rows {
            if i == skip_row {
                continue;
            }
            for j in 0..self.cols {
                if j == skip_col {
                    continue;
                }
                elems.push(self.get(i, j));
            }
        }
        Self::new(self.rows - 1, self.cols - 1, elems)
    }

    fn determinant(&self) -> Fq {
        assert_eq!(self.rows, self.cols);
        let n = self.rows;
        if n == 1 {
            return self.get(0, 0);
        }
        if n == 2 {
            return self.get(0, 0) * self.get(1, 1) - self.get(0, 1) * self.get(1, 0);
        }
        let mut det = Fq::from(0u64);
        let mut sign = Fq::from(1u64);
        let neg = Fq::from(0u64) - Fq::from(1u64);
        for j in 0..n {
            det += sign * self.get(0, j) * self.sub_matrix(0, j).determinant();
            sign = sign * neg;
        }
        det
    }

    fn inverse(&self) -> Self {
        assert_eq!(self.rows, self.cols);
        let n = self.rows;
        let det = self.determinant();
        assert!(det != Fq::from(0u64), "matrix is singular");
        let det_inv = Fq::from(1u64) / det;

        let mut adj = Self::zeros(n, n);
        let neg = Fq::from(0u64) - Fq::from(1u64);
        for i in 0..n {
            for j in 0..n {
                let minor_det = self.sub_matrix(i, j).determinant();
                let sign = if (i + j) % 2 == 0 {
                    Fq::from(1u64)
                } else {
                    neg
                };
                adj.set(j, i, sign * minor_det * det_inv);
            }
        }
        adj
    }

    fn hat(&self) -> Self {
        assert_eq!(self.rows, self.cols);
        self.sub_matrix(0, 0)
    }

    fn v_row(&self) -> Self {
        let mut elems = Vec::with_capacity(self.cols - 1);
        for j in 1..self.cols {
            elems.push(self.get(0, j));
        }
        Self::new(1, self.cols - 1, elems)
    }

    fn w_col(&self) -> Self {
        let mut elems = Vec::with_capacity(self.rows - 1);
        for i in 1..self.rows {
            elems.push(self.get(i, 0));
        }
        Self::new(self.rows - 1, 1, elems)
    }

    fn row_vec(&self, r: usize) -> Self {
        let mut elems = Vec::with_capacity(self.cols);
        for j in 0..self.cols {
            elems.push(self.get(r, j));
        }
        Self::new(1, self.cols, elems)
    }

    fn set_row(&mut self, r: usize, vals: &[Fq]) {
        assert_eq!(vals.len(), self.cols);
        for j in 0..self.cols {
            self.set(r, j, vals[j]);
        }
    }
}

// ---------------------------------------------------------------------------
// Montgomery limb extraction via transmute
// ---------------------------------------------------------------------------

fn montgomery_limbs(fq: Fq) -> [u64; 4] {
    unsafe { std::mem::transmute(fq) }
}

fn format_fq(fq: Fq) -> String {
    let limbs = montgomery_limbs(fq);
    format!(
        "Fq::from_montgomery_limbs([\n                {}, {}, {}, {},\n            ])",
        limbs[0], limbs[1], limbs[2], limbs[3]
    )
}

// ---------------------------------------------------------------------------
// Poseidon parameter generation
// ---------------------------------------------------------------------------

const T: usize = 3;
const R_F: usize = 8;
const R_P: usize = 57;
const TOTAL_ROUNDS: usize = R_F + R_P;

fn parse_hex(s: &str) -> Fq {
    let s = s.strip_prefix("0x").unwrap_or(s);
    let padded = format!("{:0>64}", s);
    let limbs: [u64; 4] = [
        u64::from_str_radix(&padded[48..64], 16).unwrap(),
        u64::from_str_radix(&padded[32..48], 16).unwrap(),
        u64::from_str_radix(&padded[16..32], 16).unwrap(),
        u64::from_str_radix(&padded[0..16], 16).unwrap(),
    ];
    let b64 = Fq::from(1u64 << 32) * Fq::from(1u64 << 32);
    Fq::from(limbs[0]) + Fq::from(limbs[1]) * b64 + Fq::from(limbs[2]) * b64 * b64 + Fq::from(limbs[3]) * b64 * b64 * b64
}

fn load_iden3_constants(path: &str) -> (Mat, Mat) {
    let json_str = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("failed to read {}: {}", path, e));
    let root: Value = serde_json::from_str(&json_str).expect("invalid JSON");

    let c_arr = root["C"][1].as_array().expect("C[1] must be an array");
    assert_eq!(c_arr.len(), TOTAL_ROUNDS * T, "expected {} round constants, got {}", TOTAL_ROUNDS * T, c_arr.len());
    let arc_elems: Vec<Fq> = c_arr.iter()
        .map(|v| parse_hex(v.as_str().expect("C entry must be string")))
        .collect();

    let m_arr = root["M"][1].as_array().expect("M[1] must be an array");
    assert_eq!(m_arr.len(), T, "expected {} MDS rows, got {}", T, m_arr.len());
    let mut mds_elems = Vec::with_capacity(T * T);
    for row in m_arr {
        let row_arr = row.as_array().expect("M row must be an array");
        assert_eq!(row_arr.len(), T);
        for val in row_arr {
            mds_elems.push(parse_hex(val.as_str().expect("M entry must be string")));
        }
    }

    (Mat::new(T, T, mds_elems), Mat::new(TOTAL_ROUNDS, T, arc_elems))
}

// Appendix B optimized round constants
fn optimized_arc(arc: &Mat, mds: &Mat) -> Mat {
    let mut ct = arc.clone();
    let r_f = R_F / 2;
    let mds_t = mds.transpose();
    let mds_inv = mds_t.inverse();

    for r in ((r_f)..(TOTAL_ROUNDS - 1 - r_f)).rev() {
        let row_rp1 = ct.row_vec(r + 1);
        let inv_cip1 = row_rp1.mul(&mds_inv);
        assert_eq!(inv_cip1.cols, T);

        for j in 1..T {
            let cur = ct.get(r, j);
            ct.set(r, j, cur + inv_cip1.get(0, j));
        }

        let mut new_row = vec![Fq::from(0u64); T];
        new_row[0] = inv_cip1.get(0, 0);
        ct.set_row(r + 1, &new_row);
    }
    ct
}

fn prime_matrix(m_hat: &Mat) -> Mat {
    let dim = m_hat.rows + 1;
    let mut out = Mat::zeros(dim, dim);
    out.set(0, 0, Fq::from(1u64));
    for i in 1..dim {
        for j in 1..dim {
            out.set(i, j, m_hat.get(i - 1, j - 1));
        }
    }
    out
}

fn doubleprime_matrix(m_hat_inv: &Mat, w: &Mat, v: &Mat, m00: Fq) -> Mat {
    let dim = m_hat_inv.rows + 1;
    let w_hat = m_hat_inv.mul(w);
    let mut out = Mat::zeros(dim, dim);

    for i in 0..dim {
        for j in 0..dim {
            if i == 0 && j == 0 {
                out.set(i, j, m00);
            } else if i == 0 {
                out.set(i, j, v.get(0, j - 1));
            } else if j == 0 {
                out.set(i, j, w_hat.get(i - 1, 0));
            } else {
                out.set(
                    i,
                    j,
                    if i == j {
                        Fq::from(1u64)
                    } else {
                        Fq::from(0u64)
                    },
                );
            }
        }
    }
    out
}

fn calc_equivalent_matrices(mds: &Mat) -> (Mat, Vec<Mat>, Vec<Mat>) {
    let m_t = mds.transpose();
    let mut m_mul = m_t.clone();
    let mut m_i = prime_matrix(&m_mul.hat());
    let mut v_collection = Vec::with_capacity(R_P);
    let mut w_hat_collection = Vec::with_capacity(R_P);

    for _ in (0..R_P).rev() {
        let m_hat = m_mul.hat();
        let w = m_mul.w_col();
        let v = m_mul.v_row();
        v_collection.push(v);
        let m_hat_inv = m_hat.inverse();
        let w_hat = m_hat_inv.mul(&w);
        w_hat_collection.push(w_hat);
        m_i = prime_matrix(&m_hat);
        m_mul = Mat::new(T, T, m_t.mul(&m_i).data);
    }

    (m_i.transpose(), v_collection, w_hat_collection)
}

// ---------------------------------------------------------------------------
// Code emitter
// ---------------------------------------------------------------------------

fn emit_fq_array(out: &mut String, elements: &[Fq], indent: &str) {
    for (i, fq) in elements.iter().enumerate() {
        let limbs = montgomery_limbs(*fq);
        write!(
            out,
            "{indent}Fq::from_montgomery_limbs([\n{indent}    {}, {}, {}, {},\n{indent}]),\n",
            limbs[0], limbs[1], limbs[2], limbs[3]
        )
        .unwrap();
        let _ = i;
    }
}

fn emit_matrix_1xN(out: &mut String, m: &Mat, type_str: &str) {
    write!(out, "{type_str}::new_from_known([\n").unwrap();
    emit_fq_array(out, &m.data, "                    ");
    write!(out, "                ])").unwrap();
}

fn emit_square_matrix(out: &mut String, m: &Mat, type_str: &str) {
    write!(out, "{type_str}::new_from_known([\n").unwrap();
    emit_fq_array(out, &m.data, "                ");
    write!(out, "            ])").unwrap();
}

fn main() {
    assert_eq!(
        std::mem::size_of::<Fq>(),
        32,
        "Fq must be 32 bytes for transmute"
    );

    let constants_path = std::env::args().nth(1).unwrap_or_else(|| {
        eprintln!("Usage: gen_params <path-to-poseidon_constants.json>");
        std::process::exit(1);
    });
    let (mds, arc) = load_iden3_constants(&constants_path);
    let opt_arc = optimized_arc(&arc, &mds);

    let m_hat = mds.hat();
    let m_hat_inv = m_hat.inverse();
    let v = mds.v_row();
    let w = mds.w_col();
    let m_prime = prime_matrix(&m_hat);
    let m00 = mds.get(0, 0);
    let m_doubleprime = doubleprime_matrix(&m_hat_inv, &w, &v, m00);
    let m_inverse = mds.inverse();

    let (m_i, v_collection, w_hat_collection) = calc_equivalent_matrices(&mds);

    let mut s = String::with_capacity(1024 * 256);

    write!(
        s,
        r#"use cycles_curve_bn254::Fq;
use poseidon_parameters::v1::{{
    Alpha, ArcMatrix, Matrix, MdsMatrix, OptimizedArcMatrix, OptimizedMdsMatrices,
    PoseidonParameters, RoundNumbers, SquareMatrix,
}};

/// Parameters for the rate-2 instance of Poseidon.
pub const fn rate_2() -> PoseidonParameters<3, 2, 9, 4, {TOTAL_ROUNDS}, 3, {TOTAL_ELEMENTS}, {R_P}> {{
    PoseidonParameters {{
        M: 128,
        arc: ArcMatrix::<{TOTAL_ROUNDS}, 3, {TOTAL_ELEMENTS}>::new_from_known([
"#,
        TOTAL_ROUNDS = TOTAL_ROUNDS,
        TOTAL_ELEMENTS = TOTAL_ROUNDS * T,
        R_P = R_P
    )
    .unwrap();

    emit_fq_array(&mut s, &arc.data, "            ");
    write!(s, "        ]),\n").unwrap();

    // MDS
    write!(s, "        mds: MdsMatrix::<3, 2, 9, 4>::new_from_known([\n").unwrap();
    emit_fq_array(&mut s, &mds.data, "            ");
    write!(s, "        ]),\n").unwrap();

    // Alpha and rounds
    write!(
        s,
        "        alpha: Alpha::Exponent(5),\n        rounds: RoundNumbers {{ r_P: {R_P}, r_F: {R_F} }},\n"
    )
    .unwrap();

    // Optimized MDS
    write!(s, "        optimized_mds: OptimizedMdsMatrices {{\n").unwrap();

    // M_hat
    write!(s, "            M_hat: ").unwrap();
    emit_square_matrix(&mut s, &m_hat, "SquareMatrix::<2, 4>");
    write!(s, ",\n").unwrap();

    // v
    write!(s, "            v: ").unwrap();
    emit_matrix_1xN(&mut s, &v, "Matrix::<1, 2, 2>");
    write!(s, ",\n").unwrap();

    // w
    write!(s, "            w: ").unwrap();
    emit_matrix_1xN(&mut s, &w, "Matrix::<2, 1, 2>");
    write!(s, ",\n").unwrap();

    // M_prime
    write!(s, "            M_prime: ").unwrap();
    emit_square_matrix(&mut s, &m_prime, "SquareMatrix::<3, 9>");
    write!(s, ",\n").unwrap();

    // M_doubleprime
    write!(s, "            M_doubleprime: ").unwrap();
    emit_square_matrix(&mut s, &m_doubleprime, "SquareMatrix::<3, 9>");
    write!(s, ",\n").unwrap();

    // M_inverse
    write!(s, "            M_inverse: ").unwrap();
    emit_square_matrix(&mut s, &m_inverse, "SquareMatrix::<3, 9>");
    write!(s, ",\n").unwrap();

    // M_hat_inverse
    write!(s, "            M_hat_inverse: ").unwrap();
    emit_square_matrix(&mut s, &m_hat_inv, "SquareMatrix::<2, 4>");
    write!(s, ",\n").unwrap();

    // M_00
    write!(s, "            M_00: {},\n", format_fq(m00)).unwrap();

    // M_i
    write!(s, "            M_i: ").unwrap();
    emit_square_matrix(&mut s, &m_i, "Matrix::<3, 3, 9>");
    write!(s, ",\n").unwrap();

    // v_collection
    write!(s, "            v_collection: [\n").unwrap();
    for v_mat in &v_collection {
        write!(s, "                ").unwrap();
        emit_matrix_1xN(&mut s, v_mat, "Matrix::<1, 2, 2>");
        write!(s, ",\n").unwrap();
    }
    write!(s, "            ],\n").unwrap();

    // w_hat_collection
    write!(s, "            w_hat_collection: [\n").unwrap();
    for w_mat in &w_hat_collection {
        write!(s, "                ").unwrap();
        emit_matrix_1xN(&mut s, w_mat, "Matrix::<2, 1, 2>");
        write!(s, ",\n").unwrap();
    }
    write!(s, "            ],\n").unwrap();

    write!(s, "        }},\n").unwrap();

    // Optimized ARC
    write!(
        s,
        "        optimized_arc: OptimizedArcMatrix::<{TOTAL_ROUNDS}, 3, {TOTAL_ELEMENTS}>::new_from_known([\n",
        TOTAL_ROUNDS = TOTAL_ROUNDS,
        TOTAL_ELEMENTS = TOTAL_ROUNDS * T,
    )
    .unwrap();
    emit_fq_array(&mut s, &opt_arc.data, "            ");
    write!(s, "        ]),\n").unwrap();

    write!(s, "    }}\n}}\n").unwrap();

    print!("{s}");
}
