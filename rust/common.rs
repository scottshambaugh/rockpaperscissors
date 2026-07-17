// Shared routines for the census tools. Each binary stays a standalone
// `rustc -O rust/<tool>.rs` build (no cargo): a `mod common;` line pulls this
// file in from the same directory. Everything here was validated inside the
// tool it originated in (noted per section) and then hoisted verbatim, so the
// binaries keep running the exact vetted code -- the counts in rust/ci_test.sh
// are the regression net for any edit to this file.
//
// Two graph representations coexist deliberately:
//   * beats/arc bitmasks (u16 for n<=16, or the Arc = [u64; MAXN] rows of the
//     balanced/regular family): the game as arcs, for combinatorial tests;
//   * f64 / i64 skew matrices: the same game as the payoff matrix M with
//     M[i][j] = +1 iff i beats j, for the linear-algebra tests.
#![allow(dead_code)]

// ================= digraph6 I/O (from cm_extend.rs) =================

// digraph6 decode: '&' + (63+p) + row-major arc bits. adj[i] bit j <=> i beats j.
pub fn decode(rec: &[u8], p: usize, out: &mut [[f64; 8]; 8]) {
    for row in out.iter_mut() {
        *row = [0f64; 8];
    }
    let psq = p * p;
    let payload = &rec[2..];
    let mut k = 0usize;
    'd: for &byte in payload {
        let mut bits = ((byte - 63) as u32) << 26;
        for _ in 0..6 {
            if bits & 0x8000_0000 != 0 {
                let (i, j) = (k / p, k % p);
                out[i][j] = 1.0;
                out[j][i] = -1.0;
            }
            bits <<= 1;
            k += 1;
            if k == psq {
                break 'd;
            }
        }
    }
}

// Encode an n-vertex oriented game (beats[i] bitmask: i beats j) to digraph6.
pub fn encode(beats: &[u16], n: usize, buf: &mut Vec<u8>) {
    buf.clear();
    buf.push(b'&');
    buf.push((63 + n) as u8);
    let nsq = n * n;
    let mut acc = 0u32;
    let mut nb = 0u32;
    let mut k = 0usize;
    while k < nsq {
        let (i, j) = (k / n, k % n);
        acc = (acc << 1) | ((beats[i] >> j) & 1) as u32;
        nb += 1;
        if nb == 6 {
            buf.push(63 + acc as u8);
            acc = 0;
            nb = 0;
        }
        k += 1;
    }
    if nb > 0 {
        buf.push(63 + (acc << (6 - nb)) as u8);
    }
    buf.push(b'\n');
}

// ================= exact integer Pfaffian (from cm_filter.rs) =================

// Pfaffian of the skew submatrix of m indexed by the set bits of `mask`,
// expanding along the lowest set index: Pf = sum_t (-1)^(t+1) a_{i0,it} Pf(rest).
// |Pf| of a k x k {-1,0,1} skew matrix is at most k^(k/4) (Hadamard), i.e.
// <= 4096 for k = 12, so i64 has room to spare through any n this runs at.
pub fn pf(m: &[[i64; 16]; 16], mask: u16) -> i64 {
    if mask.count_ones() == 4 {
        // closed form kills two recursion levels: Pf = a01*a23 - a02*a13 + a03*a12
        let i = mask.trailing_zeros() as usize;
        let mut r = mask & (mask - 1);
        let j = r.trailing_zeros() as usize;
        r &= r - 1;
        let k = r.trailing_zeros() as usize;
        let l = (r & (r - 1)).trailing_zeros() as usize;
        return m[i][j] * m[k][l] - m[i][k] * m[j][l] + m[i][l] * m[j][k];
    }
    if mask == 0 {
        return 1;
    }
    let i0 = mask.trailing_zeros() as usize;
    let rest = mask & !(1u16 << i0);
    if rest == 0 {
        return 0; // odd-sized minor: Pfaffian undefined -> identically zero
    }
    let mut acc = 0i64;
    let mut sign = 1i64;
    let mut r = rest;
    while r != 0 {
        let j = r.trailing_zeros() as usize;
        r &= r - 1;
        let a = m[i0][j];
        if a != 0 {
            acc += sign * a * pf(m, rest & !(1u16 << j));
        }
        sign = -sign;
    }
    acc
}

// ============ extension-method linear algebra (from cm_extend.rs) ============

const TOL: f64 = 1e-7;

// Gauss-Jordan inverse of a p x p matrix (p <= 8) with partial pivoting.
// Returns None if singular. Also usable to detect the singular parents.
pub fn inverse(a: &[[f64; 8]; 8], p: usize) -> Option<[[f64; 8]; 8]> {
    let mut m = *a;
    let mut inv = [[0f64; 8]; 8];
    for i in 0..p {
        inv[i][i] = 1.0;
    }
    for col in 0..p {
        // pivot: largest |m[.][col]| at or below `col`
        let mut piv = col;
        for r in (col + 1)..p {
            if m[r][col].abs() > m[piv][col].abs() {
                piv = r;
            }
        }
        if m[piv][col].abs() < 1e-9 {
            return None; // singular
        }
        m.swap(col, piv);
        inv.swap(col, piv);
        let d = 1.0 / m[col][col]; // one reciprocal, then multiply (divisions are slow)
        for k in 0..p {
            m[col][k] *= d;
            inv[col][k] *= d;
        }
        for r in 0..p {
            if r != col {
                let f = m[r][col];
                if f != 0.0 {
                    // m's columns <= col are already reduced (0/identity) and
                    // never read again -- only eliminate the future columns
                    for k in (col + 1)..p {
                        m[r][k] -= f * m[col][k];
                    }
                    for k in 0..p {
                        inv[r][k] -= f * inv[col][k];
                    }
                }
            }
        }
    }
    Some(inv)
}

// Bound-pruned DFS over r in {-1,0,1}^p collecting those with (W r)_i < -TOL for
// all i (equivalently v' = -W r > TOL, a completely mixed extension). `s` is the
// running W r over the coordinates fixed so far; `asum[k][i] = sum_{j>=k}|W[i][j]|`
// bounds the most-negative each component can still reach. If for any i the best
// case s[i] - asum[k][i] cannot dip below -TOL, no completion works -> prune the
// whole subtree. Because Wr<0 is p simultaneous strict constraints, almost every
// prefix is pruned within the first few coordinates, turning the 3^p leaf scan
// into a walk of a tiny live subtree.
// `wc` holds W's columns in the DFS visitation order (highest L1 first, so the
// bound tightens fastest); `s` is updated in place and backtracked (no per-node
// copy); r[k] is the coefficient of the k-th ORDERED column.
#[allow(clippy::too_many_arguments)]
pub fn dfs(
    k: usize,
    p: usize,
    s: &mut [f64; 8],
    r: &mut [i32; 8],
    wc: &[[f64; 8]; 8],
    asum: &[[f64; 8]; 9],
    out: &mut Vec<[i32; 8]>,
) {
    for i in 0..p {
        if s[i] - asum[k][i] >= -TOL {
            return; // component i can never reach < -TOL
        }
    }
    if k == p {
        out.push(*r); // prune above guarantees every s[i] < -TOL here
        return;
    }
    // val = 0 first (cheapest, no update); then +/-1 with in-place update+undo
    r[k] = 0;
    dfs(k + 1, p, s, r, wc, asum, out);
    for &val in &[-1i32, 1] {
        r[k] = val;
        let f = val as f64;
        for i in 0..p {
            s[i] += f * wc[i][k];
        }
        dfs(k + 1, p, s, r, wc, asum, out);
        for i in 0..p {
            s[i] -= f * wc[i][k];
        }
    }
    r[k] = 0;
}


// Const-generic clone of `dfs` (P = number of parent vertices known at compile
// time) so LLVM can fully unroll the prune/update loops. Measured benefit at
// n=9 is a few percent (within thermal noise on the census box) -- kept because
// it is free; the dynamic `dfs` remains for other callers.
#[allow(clippy::too_many_arguments)]
pub fn dfs_g<const P: usize>(
    k: usize,
    s: &mut [f64; 8],
    r: &mut [i32; 8],
    wc: &[[f64; 8]; 8],
    asum: &[[f64; 8]; 9],
    out: &mut Vec<[i32; 8]>,
) {
    for i in 0..P {
        if s[i] - asum[k][i] >= -TOL {
            return;
        }
    }
    if k == P {
        out.push(*r);
        return;
    }
    r[k] = 0;
    dfs_g::<P>(k + 1, s, r, wc, asum, out);
    for &val in &[-1i32, 1] {
        r[k] = val;
        let f = val as f64;
        for i in 0..P {
            s[i] += f * wc[i][k];
        }
        dfs_g::<P>(k + 1, s, r, wc, asum, out);
        for i in 0..P {
            s[i] -= f * wc[i][k];
        }
    }
    r[k] = 0;
}

// ====== canonical-deletion signature prefilter (from cm_extend.rs) ======

// Is vertex p signature-maximal among the vertices `eligible` admits? Used to
// emit each extension-built child from only (approximately) one of its parents:
// emit iff the added vertex p is maximal, under a cheap isomorphism-invariant
// signature (one round of degree refinement). Some eligible vertex is maximal,
// so at least one parent reconstructs the child with its added vertex maximal
// -- no class is ever lost (sound; p itself must satisfy `eligible`). Ties
// (symmetric games with several maximal vertices) are mopped up by the final
// `sort -u`; rigid games -- almost all at n=9 -- have a unique maximum, so the
// filter is essentially exact.
//
// Signature of vertex v: (outdeg, indeg, sorted out-neighbour (outdeg,indeg),
// sorted in-neighbour (outdeg,indeg)), compared lexicographically.
pub fn sig_maximal(beats: &[u16], n: usize, p: usize, eligible: impl Fn(usize) -> bool) -> bool {
    let mut od = [0u8; 16];
    let mut id = [0u8; 16];
    for v in 0..n {
        od[v] = beats[v].count_ones() as u8;
    }
    for v in 0..n {
        let mut inn = 0u8;
        for u in 0..n {
            if beats[u] & (1 << v) != 0 {
                inn += 1;
            }
        }
        id[v] = inn;
    }
    // Fill sorted-desc neighbour degree-signatures into stack arrays (no heap).
    // dsig(u) = (od[u]<<8)|id[u]; returns (out_len, in_len).
    let fill = |v: usize, outs: &mut [u16; 16], ins: &mut [u16; 16]| -> (usize, usize) {
        let (mut no, mut ni) = (0usize, 0usize);
        for u in 0..n {
            let du = ((od[u] as u16) << 8) | id[u] as u16;
            if beats[v] & (1 << u) != 0 {
                outs[no] = du;
                no += 1;
            }
            if beats[u] & (1 << v) != 0 {
                ins[ni] = du;
                ni += 1;
            }
        }
        outs[..no].sort_unstable_by(|a, b| b.cmp(a));
        ins[..ni].sort_unstable_by(|a, b| b.cmp(a));
        (no, ni)
    };
    // key(v) > key(p) ?  compare (od, id, out-sig, in-sig) lexicographically.
    let mut po = [0u16; 16];
    let mut pi = [0u16; 16];
    let (pno, pni) = fill(p, &mut po, &mut pi);
    let mut vo = [0u16; 16];
    let mut vi = [0u16; 16];
    for v in 0..n {
        if v == p || !eligible(v) {
            continue;
        }
        // lexicographic: od, id, then the two sorted signature slices
        let ord = od[v]
            .cmp(&od[p])
            .then(id[v].cmp(&id[p]))
            .then_with(|| {
                let (vno, vni) = fill(v, &mut vo, &mut vi);
                vo[..vno].cmp(&po[..pno]).then_with(|| vi[..vni].cmp(&pi[..pni]))
            });
        if ord == std::cmp::Ordering::Greater {
            return false;
        }
    }
    true
}

// All n deletions of a completely mixed game are (nonsingular) parents, so for
// the CM stream every vertex competes.
pub fn added_is_maximal(beats: &[u16], n: usize, p: usize) -> bool {
    sig_maximal(beats, n, p, |_| true)
}


// Full-signature comparison of vertices a vs b: (od, id, sorted-desc out- and
// in-neighbour (od,id) lists), the same key sig_maximal uses. The caller passes
// precomputed degree arrays (od/id of the SAME graph `beats` describes) so
// repeated comparisons against one pivot don't recompute them.
pub fn sig_cmp_with(beats: &[u16], n: usize, od: &[u8; 16], id: &[u8; 16], a: usize, b: usize) -> std::cmp::Ordering {
    let ord = od[a].cmp(&od[b]).then(id[a].cmp(&id[b]));
    if ord != std::cmp::Ordering::Equal {
        return ord;
    }
    let fill = |v: usize, outs: &mut [u16; 16], ins: &mut [u16; 16]| -> (usize, usize) {
        let (mut no, mut ni) = (0usize, 0usize);
        for u in 0..n {
            let du = ((od[u] as u16) << 8) | id[u] as u16;
            if beats[v] & (1 << u) != 0 {
                outs[no] = du;
                no += 1;
            }
            if beats[u] & (1 << v) != 0 {
                ins[ni] = du;
                ni += 1;
            }
        }
        outs[..no].sort_unstable_by(|x, y| y.cmp(x));
        ins[..ni].sort_unstable_by(|x, y| y.cmp(x));
        (no, ni)
    };
    let mut ao = [0u16; 16];
    let mut ai = [0u16; 16];
    let (ano, ani) = fill(a, &mut ao, &mut ai);
    let mut bo = [0u16; 16];
    let mut bi = [0u16; 16];
    let (bno, bni) = fill(b, &mut bo, &mut bi);
    ao[..ano].cmp(&bo[..bno]).then_with(|| ai[..ani].cmp(&bi[..bni]))
}

// ================= fully-mixed Phase-1 LP (from inc_fast.rs) =================

const EPS: f64 = 1e-7;

// ker(M) meets the open positive orthant?  x = 1 + z, z >= 0, M z = -M*ones.
// Phase-1 simplex, ALLOCATION-FREE: with n <= 15 the tableau is at most 15 x 31,
// so it lives in fixed stack arrays (the LP can be called billions of times, so
// a per-call heap alloc would dominate). A = M (n x n), rhs b = -M*ones;
// feasible => a strictly positive kernel vector x = 1 + z exists.
const LP_N: usize = 16;
const LP_COL: usize = 2 * LP_N + 1;
pub fn fully_mixed(out: &[u16], n: usize) -> bool {
    let k = n; // z-columns
    let ncol = k + n + 1; // [ A | I(artificials) | b ]
    let mut t = [[0f64; LP_COL]; LP_N];
    for i in 0..n {
        let mut rs = 0f64;
        let mut row = [0f64; LP_N];
        for j in 0..n {
            let v = if out[i] & (1 << j) != 0 {
                1.0
            } else if out[j] & (1 << i) != 0 {
                -1.0
            } else {
                0.0
            };
            row[j] = v;
            rs += v;
        }
        let b = -rs;
        let s = if b < 0.0 { -1.0 } else { 1.0 };
        for j in 0..n {
            t[i][j] = s * row[j];
        }
        t[i][k + i] = 1.0; // artificial
        t[i][k + n] = s * b;
    }
    let mut basis = [0usize; LP_N];
    for i in 0..n {
        basis[i] = k + i;
    }
    let mut obj = [0f64; LP_COL];
    for j in 0..k {
        obj[j] = -(0..n).map(|i| t[i][j]).sum::<f64>();
    }
    obj[k + n] = -(0..n).map(|i| t[i][k + n]).sum::<f64>();
    loop {
        let mut piv_col = usize::MAX;
        for (j, &o) in obj.iter().enumerate().take(k + n) {
            if o < -EPS {
                piv_col = j;
                break;
            }
        }
        if piv_col == usize::MAX {
            break;
        }
        let mut piv_row = usize::MAX;
        let mut best = f64::INFINITY;
        for i in 0..n {
            if t[i][piv_col] > EPS {
                let r = t[i][k + n] / t[i][piv_col];
                if r < best - EPS || (r < best + EPS && (piv_row == usize::MAX || basis[i] < basis[piv_row])) {
                    best = r;
                    piv_row = i;
                }
            }
        }
        if piv_row == usize::MAX {
            break;
        }
        let pv = t[piv_row][piv_col];
        for j in 0..ncol {
            t[piv_row][j] /= pv;
        }
        for i in 0..n {
            if i != piv_row {
                let f = t[i][piv_col];
                if f.abs() > EPS {
                    for j in 0..ncol {
                        t[i][j] -= f * t[piv_row][j];
                    }
                }
            }
        }
        let f = obj[piv_col];
        if f.abs() > EPS {
            for j in 0..ncol {
                obj[j] -= f * t[piv_row][j];
            }
        }
        basis[piv_row] = piv_col;
    }
    (-obj[k + n]).abs() < 1e-6
}

// ========== exact kernel/cone machinery (inclusive strata) ==========

// integer kernel basis of skew m (n x n) if nullity == d, else None.
// fraction-free RREF; basis rows are scaled to integers.
pub fn kernel_basis_exact(m: &[[i64; 16]; 16], n: usize, d: usize) -> Option<Vec<[i64; 16]>> {
    let mut a = [[0i128; 16]; 16];
    for i in 0..n {
        for j in 0..n {
            a[i][j] = m[i][j] as i128;
        }
    }
    let mut piv_col = [usize::MAX; 16];
    let mut is_piv = [false; 16];
    let mut row = 0usize;
    for col in 0..n {
        if row >= n {
            break;
        }
        let mut pr = usize::MAX;
        for r in row..n {
            if a[r][col] != 0 {
                pr = r;
                break;
            }
        }
        if pr == usize::MAX {
            continue;
        }
        a.swap(row, pr);
        for r in 0..n {
            if r != row && a[r][col] != 0 {
                let (num, den) = (a[r][col], a[row][col]);
                for c in 0..n {
                    a[r][c] = a[r][c] * den - num * a[row][c];
                }
                // keep numbers small: divide row by gcd
                let mut g = 0i128;
                for c in 0..n {
                    g = gcd_i128(g, a[r][c].abs());
                }
                if g > 1 {
                    for c in 0..n {
                        a[r][c] /= g;
                    }
                }
            }
        }
        is_piv[col] = true;
        piv_col[row] = col;
        row += 1;
    }
    if n - row != d {
        return None;
    }
    let mut basis = Vec::with_capacity(d);
    for col in 0..n {
        if !is_piv[col] {
            let mut v = [0i128; 16];
            // free var col = 1 (scaled): v[col] = prod of pivots; v[pc] = -a[r][col]*...
            // simple rational construction then clear denominators:
            // x[col] = 1; for pivot rows: x[pc] = -a[r][col]/a[r][pc]
            // scale by lcm of pivots
            let mut l: i128 = 1;
            for r in 0..row {
                let pc = piv_col[r];
                let pv = a[r][pc];
                l = l / gcd_i128(l.abs(), pv.abs()).max(1) * pv;
            }
            let l = l.abs().max(1);
            v[col] = l;
            for r in 0..row {
                let pc = piv_col[r];
                v[pc] = -a[r][col] * (l / a[r][pc]);
            }
            let mut g = 0i128;
            for c in 0..n {
                g = gcd_i128(g, v[c].abs());
            }
            if g > 1 {
                for c in 0..n {
                    v[c] /= g;
                }
            }
            let mut out = [0i64; 16];
            for c in 0..n {
                assert!(v[c].abs() < (1i128 << 62));
                out[c] = v[c] as i64;
            }
            basis.push(out);
        }
    }
    Some(basis)
}

fn gcd_i128(a: i128, b: i128) -> i128 {
    if b == 0 { a } else { gcd_i128(b, a % b) }
}

// does the cone { y = B^T lam : y >= 0 } contain a nonzero point? B = basis
// rows (d x n), full row rank => lambda-cone pointed => every extreme ray is
// tight on >= d-1 independent constraints: enumerate (d-1)-subsets of the n
// constraint normals (columns of B), take an integer nullspace vector of the
// subset, and test both signs against all constraints. Exact integer.
pub fn cone_has_nonneg(basis: &[[i64; 16]], n: usize, d: usize) -> bool {
    let cols: Vec<[i64; 8]> = (0..n)
        .map(|j| {
            let mut c = [0i64; 8];
            for (i, b) in basis.iter().enumerate() {
                c[i] = b[j];
            }
            c
        })
        .collect();
    // enumerate (d-1)-subsets
    let mut idx = vec![0usize; d.saturating_sub(1)];
    fn rec(
        start: usize,
        k: usize,
        idx: &mut Vec<usize>,
        pos: usize,
        n: usize,
        d: usize,
        cols: &[[i64; 8]],
    ) -> bool {
        if pos == k {
            // nullspace vector of the chosen (d-1) columns (each a d-vector):
            // lam with lam . cols[i] = 0 for chosen i. Build (d-1) x d system,
            // take generalized cross product via cofactor expansion.
            let mut mtx = [[0i128; 8]; 8];
            for (r, &ci) in idx.iter().enumerate() {
                for c in 0..d {
                    mtx[r][c] = cols[ci][c] as i128;
                }
            }
            let mut lam = [0i128; 8];
            for c in 0..d {
                // cofactor: delete column c, det of (d-1)x(d-1), sign (-1)^c
                let mut sub = [[0i128; 8]; 8];
                for r in 0..(d - 1) {
                    let mut cc = 0;
                    for c2 in 0..d {
                        if c2 == c {
                            continue;
                        }
                        sub[r][cc] = mtx[r][c2];
                        cc += 1;
                    }
                }
                let dt = det_small(&sub, d - 1);
                lam[c] = if c % 2 == 0 { dt } else { -dt };
            }
            if lam[..d].iter().all(|&x| x == 0) {
                return false; // degenerate subset
            }
            // test both signs
            'sgn: for sflip in [1i128, -1] {
                for col in cols.iter() {
                    let mut dot = 0i128;
                    for c in 0..d {
                        dot += sflip * lam[c] * col[c] as i128;
                    }
                    if dot < 0 {
                        continue 'sgn;
                    }
                }
                return true;
            }
            return false;
        }
        for i in start..n {
            idx[pos] = i;
            if rec(i + 1, k, idx, pos + 1, n, d, cols) {
                return true;
            }
        }
        false
    }
    if d == 1 {
        // 1-dim kernel: nonneg iff the single basis vector is one-signed
        let pos = basis[0][..n].iter().any(|&x| x > 0);
        let neg = basis[0][..n].iter().any(|&x| x < 0);
        return !(pos && neg);
    }
    rec(0, d - 1, &mut idx, 0, n, d, &cols)
}

fn det_small(a: &[[i128; 8]; 8], m: usize) -> i128 {
    if m == 0 {
        return 1;
    }
    if m == 1 {
        return a[0][0];
    }
    if m == 2 {
        return a[0][0] * a[1][1] - a[0][1] * a[1][0];
    }
    // Laplace along first row (m <= 4 here)
    let mut d = 0i128;
    for c in 0..m {
        if a[0][c] == 0 {
            continue;
        }
        let mut sub = [[0i128; 8]; 8];
        for r in 1..m {
            let mut cc = 0;
            for c2 in 0..m {
                if c2 == c {
                    continue;
                }
                sub[r - 1][cc] = a[r][c2];
                cc += 1;
            }
        }
        let s = if c % 2 == 0 { 1 } else { -1 };
        d += s * a[0][c] * det_small(&sub, m - 1);
    }
    d
}


// ========== arc-bitmask game tests (from balanced.rs / regular.rs) ==========

pub const MAXN: usize = 16;
pub type Arc = [u64; MAXN];

pub fn full_mask(n: usize) -> u64 {
    if n >= 64 { u64::MAX } else { (1u64 << n) - 1 }
}

pub fn rel(arc: &Arc, x: usize, s: usize) -> i32 {
    if (arc[x] >> s) & 1 == 1 { 1 } else if (arc[s] >> x) & 1 == 1 { -1 } else { 0 }
}

pub fn paradoxical(arc: &Arc, n: usize) -> bool {
    let mut beaten = 0u64;
    for i in 0..n {
        if arc[i] == 0 { return false; } // no win
        beaten |= arc[i];
    }
    beaten == full_mask(n) // everyone beaten at least once
}

pub fn connected(arc: &Arc, n: usize) -> bool {
    let mut adj = [0u64; MAXN];
    for i in 0..n { adj[i] = arc[i]; }
    for i in 0..n {
        let mut r = arc[i];
        while r != 0 {
            let j = r.trailing_zeros() as usize;
            r &= r - 1;
            adj[j] |= 1u64 << i;
        }
    }
    let mut visited = 1u64;
    let mut frontier = 1u64;
    while frontier != 0 {
        let mut next = 0u64;
        let mut f = frontier;
        while f != 0 {
            let vtx = f.trailing_zeros() as usize;
            f &= f - 1;
            next |= adj[vtx];
        }
        next &= !visited;
        visited |= next;
        frontier = next;
    }
    visited == full_mask(n)
}

pub fn twin_free(arc: &Arc, n: usize) -> bool {
    let mut lose = [0u64; MAXN];
    for i in 0..n {
        let mut r = arc[i];
        while r != 0 {
            let j = r.trailing_zeros() as usize;
            r &= r - 1;
            lose[j] |= 1u64 << i;
        }
    }
    for i in 0..n {
        for j in (i + 1)..n {
            if arc[i] == arc[j] && lose[i] == lose[j] { return false; }
        }
    }
    true
}

pub fn is_prime(arc: &Arc, n: usize) -> bool {
    if n < 3 { return false; }
    // For each seed pair {a,b}, grow the smallest module containing it: any vertex x
    // that distinguishes the pair (relates differently to a vs. some member) must lie
    // in every module containing the pair, so force it in and repeat. Each external x
    // is compared against the fixed anchor a, so a vertex once consistent stays so
    // until a newly added member splits it -> each vertex enters the stack at most
    // once, giving O(n^2) per pair. If a closure stops short of V it is a proper
    // nontrivial (size >= 2) module, so the game is not prime. All pairs reaching V
    // => no proper module => prime. Every nontrivial proper module contains a pair
    // whose closure stays within it, so nothing is missed.
    let full = full_mask(n);
    let mut stack = [0usize; MAXN];
    for a in 0..n {
        for b in (a + 1)..n {
            let mut in_s = (1u64 << a) | (1u64 << b);
            let mut sp = 0usize;
            stack[sp] = b; sp += 1; // a is the anchor/reference
            while sp > 0 {
                sp -= 1;
                let s = stack[sp];
                for x in 0..n {
                    if (in_s >> x) & 1 == 1 { continue; }
                    if rel(arc, x, s) != rel(arc, x, a) {
                        in_s |= 1u64 << x;
                        stack[sp] = x; sp += 1;
                    }
                }
                if in_s == full { break; }
            }
            if in_s != full { return false; } // proper nontrivial module
        }
    }
    true
}
