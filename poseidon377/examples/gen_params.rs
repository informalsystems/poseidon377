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
        if n == 0 {
            return Fq::from(1u64);
        }
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

const R_F: usize = 8;
const N_ROUNDS_P: [usize; 16] = [56, 57, 56, 60, 60, 63, 64, 63, 60, 66, 60, 65, 70, 60, 64, 68];

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

fn load_iden3_constants(path: &str, t: usize, total_rounds: usize) -> (Mat, Mat) {
    let json_str = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("failed to read {}: {}", path, e));
    let root: Value = serde_json::from_str(&json_str).expect("invalid JSON");

    let idx = t - 2;
    let c_arr = root["C"][idx].as_array().unwrap_or_else(|| panic!("C[{}] must be an array", idx));
    assert_eq!(c_arr.len(), total_rounds * t, "expected {} round constants, got {}", total_rounds * t, c_arr.len());
    let arc_elems: Vec<Fq> = c_arr.iter()
        .map(|v| parse_hex(v.as_str().expect("C entry must be string")))
        .collect();

    let m_arr = root["M"][idx].as_array().unwrap_or_else(|| panic!("M[{}] must be an array", idx));
    assert_eq!(m_arr.len(), t, "expected {} MDS rows, got {}", t, m_arr.len());
    let mut mds_elems = Vec::with_capacity(t * t);
    for row in m_arr {
        let row_arr = row.as_array().expect("M row must be an array");
        assert_eq!(row_arr.len(), t);
        for val in row_arr {
            mds_elems.push(parse_hex(val.as_str().expect("M entry must be string")));
        }
    }

    (Mat::new(t, t, mds_elems), Mat::new(total_rounds, t, arc_elems))
}

// Appendix B optimized round constants
fn optimized_arc(arc: &Mat, mds: &Mat, t: usize, total_rounds: usize) -> Mat {
    let mut ct = arc.clone();
    let r_f = R_F / 2;
    let mds_t = mds.transpose();
    let mds_inv = mds_t.inverse();

    for r in ((r_f)..(total_rounds - 1 - r_f)).rev() {
        let row_rp1 = ct.row_vec(r + 1);
        let inv_cip1 = row_rp1.mul(&mds_inv);
        assert_eq!(inv_cip1.cols, t);

        for j in 1..t {
            let cur = ct.get(r, j);
            ct.set(r, j, cur + inv_cip1.get(0, j));
        }

        let mut new_row = vec![Fq::from(0u64); t];
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

fn calc_equivalent_matrices(mds: &Mat, t: usize, r_p: usize) -> (Mat, Vec<Mat>, Vec<Mat>) {
    let m_t = mds.transpose();
    let mut m_mul = m_t.clone();
    let mut m_i = prime_matrix(&m_mul.hat());
    let mut v_collection = Vec::with_capacity(r_p);
    let mut w_hat_collection = Vec::with_capacity(r_p);

    for _ in (0..r_p).rev() {
        let m_hat = m_mul.hat();
        let w = m_mul.w_col();
        let v = m_mul.v_row();
        v_collection.push(v);
        let m_hat_inv = m_hat.inverse();
        let w_hat = m_hat_inv.mul(&w);
        w_hat_collection.push(w_hat);
        m_i = prime_matrix(&m_hat);
        m_mul = Mat::new(t, t, m_t.mul(&m_i).data);
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

fn fq_to_hex(fq: Fq) -> String {
    let bytes = fq.to_bytes_le();
    let mut hex = String::from("0x");
    for b in bytes.iter().rev() {
        hex.push_str(&format!("{:02x}", b));
    }
    hex
}

fn verify_optimization(
    t: usize,
    r_p: usize,
    arc: &Mat,
    mds: &Mat,
    opt_arc: &Mat,
    m_i: &Mat,
    m00: Fq,
    v_collection: &[Mat],
    w_hat_collection: &[Mat],
) {
    let total_rounds = R_F + r_p;
    let r_f = R_F / 2;

    let mut state_init = vec![Fq::from(0u64); t];
    if t >= 2 { state_init[1] = Fq::from(1u64); }

    // --- Unoptimized permutation ---
    let mut state_u = state_init.clone();
    let rc = &arc.data;
    let mut rc_idx = 0;

    for _ in 0..r_f {
        for i in 0..t { state_u[i] += rc[rc_idx]; rc_idx += 1; }
        for i in 0..t { state_u[i] = state_u[i].power([5u64]); }
        let prev = state_u.clone();
        for i in 0..t {
            state_u[i] = (0..t).map(|j| mds.get(i, j) * prev[j]).sum();
        }
    }
    for _ in 0..r_p {
        for i in 0..t { state_u[i] += rc[rc_idx]; rc_idx += 1; }
        state_u[0] = state_u[0].power([5u64]);
        let prev = state_u.clone();
        for i in 0..t {
            state_u[i] = (0..t).map(|j| mds.get(i, j) * prev[j]).sum();
        }
    }
    for _ in 0..r_f {
        for i in 0..t { state_u[i] += rc[rc_idx]; rc_idx += 1; }
        for i in 0..t { state_u[i] = state_u[i].power([5u64]); }
        let prev = state_u.clone();
        for i in 0..t {
            state_u[i] = (0..t).map(|j| mds.get(i, j) * prev[j]).sum();
        }
    }

    // --- Optimized permutation ---
    let mut state_o = state_init.clone();

    // First R_f full rounds
    for r in 0..r_f {
        for i in 0..t { state_o[i] += opt_arc.get(r, i); }
        for i in 0..t { state_o[i] = state_o[i].power([5u64]); }
        let prev = state_o.clone();
        for i in 0..t {
            state_o[i] = (0..t).map(|j| mds.get(i, j) * prev[j]).sum();
        }
    }

    // First partial round ARC (full) + M_i multiply
    let mut rcc = r_f;
    for i in 0..t { state_o[i] += opt_arc.get(rcc, i); }
    {
        let prev = state_o.clone();
        for i in 0..t {
            state_o[i] = (0..t).map(|j| m_i.get(i, j) * prev[j]).sum();
        }
    }

    // Partial rounds (R_P - 1 iterations)
    for r in 0..(r_p - 1) {
        state_o[0] = state_o[0].power([5u64]);
        rcc += 1;
        state_o[0] += opt_arc.get(rcc, 0);
        // sparse_mat_mul(r_p - r - 1)
        let sp_idx = r_p - r - 1;
        let prev = state_o.clone();
        let mut add_row = vec![Fq::from(0u64); t - 1];
        for i in 0..(t - 1) {
            add_row[i] = v_collection[sp_idx].data[i] * prev[0] + prev[i + 1];
        }
        state_o[0] = m00 * prev[0]
            + (0..(t - 1))
                .map(|i| w_hat_collection[sp_idx].data[i] * prev[i + 1])
                .sum::<Fq>();
        for i in 0..(t - 1) {
            state_o[i + 1] = add_row[i];
        }
    }

    // Last partial round
    state_o[0] = state_o[0].power([5u64]);
    {
        let sp_idx = 0;
        let prev = state_o.clone();
        let mut add_row = vec![Fq::from(0u64); t - 1];
        for i in 0..(t - 1) {
            add_row[i] = v_collection[sp_idx].data[i] * prev[0] + prev[i + 1];
        }
        state_o[0] = m00 * prev[0]
            + (0..(t - 1))
                .map(|i| w_hat_collection[sp_idx].data[i] * prev[i + 1])
                .sum::<Fq>();
        for i in 0..(t - 1) {
            state_o[i + 1] = add_row[i];
        }
    }
    rcc += 1;

    // Final R_f full rounds
    for _ in 0..r_f {
        for i in 0..t { state_o[i] += opt_arc.get(rcc, i); }
        for i in 0..t { state_o[i] = state_o[i].power([5u64]); }
        let prev = state_o.clone();
        for i in 0..t {
            state_o[i] = (0..t).map(|j| mds.get(i, j) * prev[j]).sum();
        }
        rcc += 1;
    }

    eprintln!("  Unoptimized[0]: {}", fq_to_hex(state_u[0]));
    eprintln!("  Optimized[0]:   {}", fq_to_hex(state_o[0]));
    if state_u[0] == state_o[0] {
        eprintln!("  Permutation MATCH ✓");
    } else {
        eprintln!("  Permutation MISMATCH ✗ — optimization is buggy for t={}", t);

        bisect_divergence(t, r_p, arc, mds, opt_arc, m_i, m00, v_collection, w_hat_collection);
    }

    // Always check matrix decomposition
    eprintln!("\n  Verifying M_i * sparse_product == MDS^R_P ...");
    let total_rounds = R_F + r_p;
    let mut mds_power = Mat::identity(t);
    for _ in 0..r_p {
        mds_power = mds.mul(&mds_power);
    }

    let mut sparse_product = Mat::identity(t);
    for k in (0..r_p).rev() {
        let mut sp_mat = Mat::identity(t);
        sp_mat.set(0, 0, m00);
        for i in 0..(t - 1) {
            sp_mat.set(0, i + 1, w_hat_collection[k].data[i]);
            sp_mat.set(i + 1, 0, v_collection[k].data[i]);
        }
        sparse_product = sp_mat.mul(&sparse_product);
    }

    let mi_sparse = m_i.mul(&sparse_product);
    let match_product = (0..t*t).all(|idx| mi_sparse.data[idx] == mds_power.data[idx]);
    eprintln!("  M_i * sparse_product == MDS^R_P: {}", if match_product { "YES ✓" } else { "NO ✗" });
}

fn bisect_divergence(
    t: usize,
    r_p: usize,
    arc: &Mat,
    mds: &Mat,
    opt_arc: &Mat,
    m_i: &Mat,
    m00: Fq,
    v_collection: &[Mat],
    w_hat_collection: &[Mat],
) {
    let r_f = R_F / 2;

    let mut state_init = vec![Fq::from(0u64); t];
    if t >= 2 { state_init[1] = Fq::from(1u64); }

    // Run unoptimized step-by-step, saving state after each round
    let mut state_u = state_init.clone();
    let rc = &arc.data;
    let mut rc_idx = 0;
    let mut unopt_states: Vec<Vec<Fq>> = Vec::new();

    for _ in 0..r_f {
        for i in 0..t { state_u[i] += rc[rc_idx]; rc_idx += 1; }
        for i in 0..t { state_u[i] = state_u[i].power([5u64]); }
        let prev = state_u.clone();
        for i in 0..t {
            state_u[i] = (0..t).map(|j| mds.get(i, j) * prev[j]).sum();
        }
        unopt_states.push(state_u.clone());
    }
    for _ in 0..r_p {
        for i in 0..t { state_u[i] += rc[rc_idx]; rc_idx += 1; }
        state_u[0] = state_u[0].power([5u64]);
        let prev = state_u.clone();
        for i in 0..t {
            state_u[i] = (0..t).map(|j| mds.get(i, j) * prev[j]).sum();
        }
        unopt_states.push(state_u.clone());
    }
    for _ in 0..r_f {
        for i in 0..t { state_u[i] += rc[rc_idx]; rc_idx += 1; }
        for i in 0..t { state_u[i] = state_u[i].power([5u64]); }
        let prev = state_u.clone();
        for i in 0..t {
            state_u[i] = (0..t).map(|j| mds.get(i, j) * prev[j]).sum();
        }
        unopt_states.push(state_u.clone());
    }

    // Run optimized step-by-step
    let mut state_o = state_init.clone();

    for r in 0..r_f {
        for i in 0..t { state_o[i] += opt_arc.get(r, i); }
        for i in 0..t { state_o[i] = state_o[i].power([5u64]); }
        let prev = state_o.clone();
        for i in 0..t {
            state_o[i] = (0..t).map(|j| mds.get(i, j) * prev[j]).sum();
        }
        let round_num = r;
        if state_o != unopt_states[round_num] {
            eprintln!("  Divergence at full round {} (round_num={})", r, round_num);
            eprintln!("    Unopt: {:?}", unopt_states[round_num].iter().map(|x| fq_to_hex(*x)).collect::<Vec<_>>());
            eprintln!("    Opt:   {:?}", state_o.iter().map(|x| fq_to_hex(*x)).collect::<Vec<_>>());
            return;
        }
    }
    eprintln!("  First {} full rounds match.", r_f);

    // Now the tricky part: the optimized partial rounds transform the state differently
    // (ARC absorption + sparse matrices), so intermediate states won't match round-by-round.
    // Instead, the states should match AFTER ALL partial rounds complete.
    //
    // Let's check after all partial rounds:
    let mut rcc = r_f;
    for i in 0..t { state_o[i] += opt_arc.get(rcc, i); }
    {
        let prev = state_o.clone();
        for i in 0..t {
            state_o[i] = (0..t).map(|j| m_i.get(i, j) * prev[j]).sum();
        }
    }
    for r in 0..(r_p - 1) {
        state_o[0] = state_o[0].power([5u64]);
        rcc += 1;
        state_o[0] += opt_arc.get(rcc, 0);
        let sp_idx = r_p - r - 1;
        let prev = state_o.clone();
        let mut add_row = vec![Fq::from(0u64); t - 1];
        for i in 0..(t - 1) {
            add_row[i] = v_collection[sp_idx].data[i] * prev[0] + prev[i + 1];
        }
        state_o[0] = m00 * prev[0]
            + (0..(t - 1))
                .map(|i| w_hat_collection[sp_idx].data[i] * prev[i + 1])
                .sum::<Fq>();
        for i in 0..(t - 1) {
            state_o[i + 1] = add_row[i];
        }
    }
    state_o[0] = state_o[0].power([5u64]);
    {
        let sp_idx = 0;
        let prev = state_o.clone();
        let mut add_row = vec![Fq::from(0u64); t - 1];
        for i in 0..(t - 1) {
            add_row[i] = v_collection[sp_idx].data[i] * prev[0] + prev[i + 1];
        }
        state_o[0] = m00 * prev[0]
            + (0..(t - 1))
                .map(|i| w_hat_collection[sp_idx].data[i] * prev[i + 1])
                .sum::<Fq>();
        for i in 0..(t - 1) {
            state_o[i + 1] = add_row[i];
        }
    }
    rcc += 1;

    let after_partial_idx = r_f + r_p - 1;
    if state_o == unopt_states[after_partial_idx] {
        eprintln!("  After all partial rounds: MATCH");
    } else {
        eprintln!("  After all partial rounds: MISMATCH");
        eprintln!("    Unopt: {:?}", unopt_states[after_partial_idx].iter().map(|x| fq_to_hex(*x)).collect::<Vec<_>>());
        eprintln!("    Opt:   {:?}", state_o.iter().map(|x| fq_to_hex(*x)).collect::<Vec<_>>());
        eprintln!("  Bug is in the partial round optimization (ARC or sparse matrices).");

        // Also verify: does the M_i * product-of-sparse-matrices equal M^R_P?
        eprintln!("\n  Verifying M_i * sparse_product == MDS^R_P ...");
        let mut mds_power = Mat::identity(t);
        for _ in 0..r_p {
            mds_power = mds.mul(&mds_power);
        }

        let mut sparse_product = Mat::identity(t);
        for k in (0..r_p).rev() {
            let mut sp_mat = Mat::identity(t);
            sp_mat.set(0, 0, m00);
            for i in 0..(t - 1) {
                sp_mat.set(0, i + 1, w_hat_collection[k].data[i]);
                sp_mat.set(i + 1, 0, v_collection[k].data[i]);
            }
            sparse_product = sp_mat.mul(&sparse_product);
        }

        let mi_sparse = m_i.mul(&sparse_product);

        let match_product = (0..t*t).all(|idx| mi_sparse.data[idx] == mds_power.data[idx]);
        if match_product {
            eprintln!("  M_i * sparse_product == MDS^R_P: YES ✓");
            eprintln!("  Matrix decomposition is correct; bug is in ARC optimization.");
        } else {
            eprintln!("  M_i * sparse_product == MDS^R_P: NO ✗");
            eprintln!("  Matrix decomposition is WRONG.");

            eprintln!("  MDS^R_P:");
            for i in 0..t {
                for j in 0..t {
                    eprintln!("    [{},{}] = {}", i, j, fq_to_hex(mds_power.get(i, j)));
                }
            }
            eprintln!("  M_i * sparse_product:");
            for i in 0..t {
                for j in 0..t {
                    eprintln!("    [{},{}] = {}", i, j, fq_to_hex(mi_sparse.get(i, j)));
                }
            }
        }
    }
}

fn main() {
    assert_eq!(
        std::mem::size_of::<Fq>(),
        32,
        "Fq must be 32 bytes for transmute"
    );

    let args: Vec<String> = std::env::args().collect();
    if args.len() != 3 {
        eprintln!("Usage: gen_params <t> <path-to-poseidon_constants.json>");
        std::process::exit(1);
    }
    let t: usize = args[1].parse().expect("t must be an integer");
    assert!(t >= 2 && t <= 17, "t must be in [2, 17]");
    let constants_path = &args[2];

    let r_p = N_ROUNDS_P[t - 2];
    let total_rounds = R_F + r_p;
    let rate = t - 1;
    let total_elements = total_rounds * t;
    let mds_elements = t * t;
    let sm1 = t - 1;
    let sm1_elements = sm1 * sm1;

    eprintln!("t={}, rate={}, R_F={}, R_P={}, total_rounds={}, total_elements={}", t, rate, R_F, r_p, total_rounds, total_elements);

    let (mds, arc) = load_iden3_constants(constants_path, t, total_rounds);
    let opt_arc = optimized_arc(&arc, &mds, t, total_rounds);

    let m_hat = mds.hat();
    let m_hat_inv = m_hat.inverse();
    let v = mds.v_row();
    let w = mds.w_col();
    let m_prime = prime_matrix(&m_hat);
    let m00 = mds.get(0, 0);
    let m_doubleprime = doubleprime_matrix(&m_hat_inv, &w, &v, m00);
    let m_inverse = mds.inverse();

    let (m_i, v_collection, w_hat_collection) = calc_equivalent_matrices(&mds, t, r_p);

    eprintln!("\nVerifying optimization for t={}...", t);
    verify_optimization(t, r_p, &arc, &mds, &opt_arc, &m_i, m00, &v_collection, &w_hat_collection);

    if t == 2 {
        eprintln!("\n  M_00 = {}", fq_to_hex(m00));
        eprintln!("  Dumping first 3 and last 3 sparse entries (iden3 S order = reversed):");
        for idx in [r_p - 1, r_p - 2, r_p - 3, 2usize, 1, 0] {
            eprintln!("    v_collection[{}] = {}", idx, fq_to_hex(v_collection[idx].data[0]));
            eprintln!("    w_hat_collection[{}] = {}", idx, fq_to_hex(w_hat_collection[idx].data[0]));
        }
    }

    let mut s = String::with_capacity(1024 * 256);

    write!(
        s,
        r#"use cycles_curve_bn254::Fq;
use poseidon_parameters::v1::{{
    Alpha, ArcMatrix, Matrix, MdsMatrix, OptimizedArcMatrix, OptimizedMdsMatrices,
    PoseidonParameters, RoundNumbers, SquareMatrix,
}};

/// Parameters for the rate-{RATE} instance of Poseidon.
pub const fn rate_{RATE}() -> PoseidonParameters<{T}, {SM1}, {MDS_ELEMS}, {SM1_ELEMS}, {TOTAL_ROUNDS}, {T2}, {TOTAL_ELEMENTS}, {R_P}> {{
    PoseidonParameters {{
        M: 128,
        arc: ArcMatrix::<{TOTAL_ROUNDS}, {T2}, {TOTAL_ELEMENTS}>::new_from_known([
"#,
        RATE = rate,
        T = t,
        SM1 = sm1,
        MDS_ELEMS = mds_elements,
        SM1_ELEMS = sm1_elements,
        TOTAL_ROUNDS = total_rounds,
        T2 = t,
        TOTAL_ELEMENTS = total_elements,
        R_P = r_p,
    )
    .unwrap();

    emit_fq_array(&mut s, &arc.data, "            ");
    write!(s, "        ]),\n").unwrap();

    write!(s, "        mds: MdsMatrix::<{}, {}, {}, {}>::new_from_known([\n", t, sm1, mds_elements, sm1_elements).unwrap();
    emit_fq_array(&mut s, &mds.data, "            ");
    write!(s, "        ]),\n").unwrap();

    write!(
        s,
        "        alpha: Alpha::Exponent(5),\n        rounds: RoundNumbers {{ r_P: {}, r_F: {} }},\n",
        r_p, R_F
    )
    .unwrap();

    write!(s, "        optimized_mds: OptimizedMdsMatrices {{\n").unwrap();

    write!(s, "            M_hat: ").unwrap();
    emit_square_matrix(&mut s, &m_hat, &format!("SquareMatrix::<{}, {}>", sm1, sm1_elements));
    write!(s, ",\n").unwrap();

    write!(s, "            v: ").unwrap();
    emit_matrix_1xN(&mut s, &v, &format!("Matrix::<1, {}, {}>", sm1, sm1));
    write!(s, ",\n").unwrap();

    write!(s, "            w: ").unwrap();
    emit_matrix_1xN(&mut s, &w, &format!("Matrix::<{}, 1, {}>", sm1, sm1));
    write!(s, ",\n").unwrap();

    write!(s, "            M_prime: ").unwrap();
    emit_square_matrix(&mut s, &m_prime, &format!("SquareMatrix::<{}, {}>", t, mds_elements));
    write!(s, ",\n").unwrap();

    write!(s, "            M_doubleprime: ").unwrap();
    emit_square_matrix(&mut s, &m_doubleprime, &format!("SquareMatrix::<{}, {}>", t, mds_elements));
    write!(s, ",\n").unwrap();

    write!(s, "            M_inverse: ").unwrap();
    emit_square_matrix(&mut s, &m_inverse, &format!("SquareMatrix::<{}, {}>", t, mds_elements));
    write!(s, ",\n").unwrap();

    write!(s, "            M_hat_inverse: ").unwrap();
    emit_square_matrix(&mut s, &m_hat_inv, &format!("SquareMatrix::<{}, {}>", sm1, sm1_elements));
    write!(s, ",\n").unwrap();

    write!(s, "            M_00: {},\n", format_fq(m00)).unwrap();

    write!(s, "            M_i: ").unwrap();
    emit_square_matrix(&mut s, &m_i, &format!("Matrix::<{}, {}, {}>", t, t, mds_elements));
    write!(s, ",\n").unwrap();

    write!(s, "            v_collection: [\n").unwrap();
    for v_mat in &v_collection {
        write!(s, "                ").unwrap();
        emit_matrix_1xN(&mut s, v_mat, &format!("Matrix::<1, {}, {}>", sm1, sm1));
        write!(s, ",\n").unwrap();
    }
    write!(s, "            ],\n").unwrap();

    write!(s, "            w_hat_collection: [\n").unwrap();
    for w_mat in &w_hat_collection {
        write!(s, "                ").unwrap();
        emit_matrix_1xN(&mut s, w_mat, &format!("Matrix::<{}, 1, {}>", sm1, sm1));
        write!(s, ",\n").unwrap();
    }
    write!(s, "            ],\n").unwrap();

    write!(s, "        }},\n").unwrap();

    write!(
        s,
        "        optimized_arc: OptimizedArcMatrix::<{}, {}, {}>::new_from_known([\n",
        total_rounds, t, total_elements,
    )
    .unwrap();
    emit_fq_array(&mut s, &opt_arc.data, "            ");
    write!(s, "        ]),\n").unwrap();

    write!(s, "    }}\n}}\n").unwrap();

    print!("{s}");
}
