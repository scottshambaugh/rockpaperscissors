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


// Gauss-Jordan inverse of a p x p matrix (p <= 8) with partial pivoting.
// Returns None if singular. Also usable to detect the singular parents.

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


// Const-generic clone of `dfs` (P = number of parent vertices known at compile
// time) so LLVM can fully unroll the prune/update loops. Measured benefit at
// n=9 is a few percent (within thermal noise on the census box) -- kept because
// it is free; the dynamic `dfs` remains for other callers.
#[allow(clippy::too_many_arguments)]

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


// ker(M) meets the open positive orthant?  x = 1 + z, z >= 0, M z = -M*ones.
// Phase-1 simplex, ALLOCATION-FREE: with n <= 15 the tableau is at most 15 x 31,
// so it lives in fixed stack arrays (the LP can be called billions of times, so
// a per-call heap alloc would dominate). A = M (n x n), rhs b = -M*ones;
// feasible => a strictly positive kernel vector x = 1 + z exists.

// ========== exact kernel/cone machinery (inclusive strata) ==========

// integer kernel basis of skew m (n x n) if nullity == d, else None.
// fraction-free RREF; basis rows are scaled to integers.
pub fn kernel_basis_exact(m: &[[i64; 16]; 16], n: usize, d: usize) -> Option<Vec<[i64; 16]>> {
    let b = kernel_basis_any(m, n);
    if b.len() == d { Some(b) } else { None }
}

// kernel basis at whatever the nullity is (possibly empty), one RREF pass
pub fn kernel_basis_any(m: &[[i64; 16]; 16], n: usize) -> Vec<[i64; 16]> {
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
    let d = n - row;
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
    basis
}


// checked-i64 mirror of kernel_basis_any + the nullity-d filter: identical
// arithmetic in a narrower type (entries here are tiny), Err(()) on any
// overflow so the caller can fall back to the i128 version.
pub fn kernel_basis_exact_fast(m: &[[i64; 16]; 16], n: usize, d: usize) -> Option<Vec<[i64; 16]>> {
    match kernel_basis_any_i64(m, n) {
        Ok(b) => {
            if b.len() == d {
                Some(b)
            } else {
                None
            }
        }
        Err(()) => kernel_basis_exact(m, n, d),
    }
}

#[allow(dead_code)]
pub fn gcd_i64(a: i64, b: i64) -> i64 {
    if b == 0 { a } else { gcd_i64(b, a % b) }
}

pub fn kernel_basis_any_fast(m: &[[i64; 16]; 16], n: usize) -> Vec<[i64; 16]> {
    match kernel_basis_any_i64(m, n) {
        Ok(b) => b,
        Err(()) => kernel_basis_any(m, n),
    }
}

fn kernel_basis_any_i64(m: &[[i64; 16]; 16], n: usize) -> Result<Vec<[i64; 16]>, ()> {
    let mut a = [[0i64; 16]; 16];
    for i in 0..n {
        a[i][..n].copy_from_slice(&m[i][..n]);
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
                    a[r][c] = a[r][c].checked_mul(den).ok_or(())?.checked_sub(num.checked_mul(a[row][c]).ok_or(())?).ok_or(())?;
                }
                let mut g = 0i64;
                for c in 0..n {
                    g = gcd_i64(g, a[r][c].abs());
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
    let d = n - row;
    let mut basis = Vec::with_capacity(d);
    for col in 0..n {
        if !is_piv[col] {
            let mut v = [0i64; 16];
            let mut l: i64 = 1;
            for r in 0..row {
                let pc = piv_col[r];
                let pv = a[r][pc];
                l = (l / gcd_i64(l.abs(), pv.abs()).max(1)).checked_mul(pv).ok_or(())?;
            }
            let l = l.abs().max(1);
            v[col] = l;
            for r in 0..row {
                let pc = piv_col[r];
                v[pc] = a[r][col].checked_mul(l / a[r][pc]).ok_or(())?.checked_neg().ok_or(())?;
            }
            let mut g = 0i64;
            for c in 0..n {
                g = gcd_i64(g, v[c].abs());
            }
            if g > 1 {
                for c in 0..n {
                    v[c] /= g;
                }
            }
            basis.push(v);
        }
    }
    Ok(basis)
}

pub fn gcd_i128(a: i128, b: i128) -> i128 {
    if b == 0 { a } else { gcd_i128(b, a % b) }
}

// rps_autsize (nauty group order) comes back as f64. An f64 represents integers
// exactly only below 2^53; above that `as u128` silently rounds -- the same
// float-truncation class of bug that once corrupted an awk sum. |Aut| divides
// n!, so it is exact for every n we can enumerate, but GUARD it at runtime so a
// future run panics loudly instead of returning a rounded count. 9.0e15 < 2^53.
#[allow(dead_code)]
pub fn autsize_u128(x: f64) -> u128 {
    assert!(
        x >= 0.0 && x.fract() == 0.0 && x < 9.0e15,
        "rps_autsize returned {x}, not an exact non-negative integer below 2^53"
    );
    x as u128
}

// Checked narrowing casts. `overflow-checks=on` guards +,-,* but NOT `as`
// truncation, so every narrowing cast that carries a mathematical value goes
// through these -- a future (larger n) run panics loudly instead of silently
// wrapping. Bounds hold today by Hadamard (|det|<=~4096, |adj|<=~907) but are
// NOT assumed; they are verified at runtime. Edition-agnostic (no TryFrom).
#[allow(dead_code)]
pub fn narrow_i64(x: i128) -> i64 {
    assert!(x >= i64::MIN as i128 && x <= i64::MAX as i128, "value {x} overflows i64");
    x as i64
}
#[allow(dead_code)]
pub fn narrow_i32(x: i64) -> i32 {
    assert!(x >= i32::MIN as i64 && x <= i32::MAX as i64, "value {x} overflows i32");
    x as i32
}

// Build the {-1,0,1} skew matrix from an out-adjacency bitmask (beats[i] bit j
// set iff i beats j). Lets the exact feasibility test (has_positive_kernel)
// stand in for the float LP (fully_mixed) at every call site.
#[allow(dead_code)]
pub fn skew_from_beats(beats: &[u16], n: usize) -> [[i64; 16]; 16] {
    let mut m = [[0i64; 16]; 16];
    for i in 0..n {
        let mut w = beats[i];
        while w != 0 {
            let j = w.trailing_zeros() as usize;
            w &= w - 1;
            m[i][j] = 1;
            m[j][i] = -1;
        }
    }
    m
}

// Exact i64 mirror of the float `dfs`: collect every r in {-1,0,1}^p with
// (wc . r)_i < 0 for ALL i, via the same suffix-bound prune. Used for the
// completely-mixed extension: with M' nonsingular skew (det = Pf^2 > 0),
// -M'^-1 r > 0 componentwise  <=>  sgn(det)*(adj M' . r)_i < 0 for all i, so
// wc = sgn(det)*adj columns (reordered) and this enumerates the CM children
// exactly -- no f64, no tolerance. Entries are Hadamard-bounded so the running
// sums fit i64 (and every + - * is overflow-checked in release).
#[allow(dead_code)]
#[allow(clippy::too_many_arguments)]
pub fn dfs_neg(
    k: usize,
    p: usize,
    s: &mut [i64; 16],
    r: &mut [i32; 16],
    wc: &[[i64; 16]; 16],
    asum: &[[i64; 16]; 17],
    out: &mut Vec<[i32; 16]>,
) {
    for i in 0..p {
        if s[i] - asum[k][i] >= 0 {
            return; // component i can never reach < 0
        }
    }
    if k == p {
        out.push(*r);
        return;
    }
    r[k] = 0;
    dfs_neg(k + 1, p, s, r, wc, asum, out);
    for &val in &[-1i32, 1] {
        r[k] = val;
        for i in 0..p {
            s[i] += val as i64 * wc[i][k];
        }
        dfs_neg(k + 1, p, s, r, wc, asum, out);
        for i in 0..p {
            s[i] -= val as i64 * wc[i][k];
        }
    }
    r[k] = 0;
}

// Exact perpendicularity DFS: collect every r in {-1,0,1}^p with (rows_i . r)==0
// for ALL i in 0..d (d integer equality rows = r _|_ each kernel-basis vector).
// Prune when |s[i]| > asum[i][k] (remaining columns can't return component i to
// 0). Used for the odd-n nullity>=3 extension (child nullity d+1 iff r _|_ ker).
#[allow(dead_code)]
#[allow(clippy::too_many_arguments)]
pub fn dfs_perp_eq(
    k: usize,
    p: usize,
    d: usize,
    s: &mut [i64; 16],
    r: &mut [i32; 16],
    rows: &[[i64; 16]; 16],
    asum: &[[i64; 16]; 17],
    out: &mut Vec<[i32; 16]>,
) {
    for i in 0..d {
        if s[i].abs() > asum[i][k] {
            return; // component i can never return to 0
        }
    }
    if k == p {
        out.push(*r); // prune above guarantees s[i]==0 for all i here
        return;
    }
    r[k] = 0;
    dfs_perp_eq(k + 1, p, d, s, r, rows, asum, out);
    for &val in &[-1i32, 1] {
        r[k] = val;
        for i in 0..d {
            s[i] += val as i64 * rows[i][k];
        }
        dfs_perp_eq(k + 1, p, d, s, r, rows, asum, out);
        for i in 0..d {
            s[i] -= val as i64 * rows[i][k];
        }
    }
    r[k] = 0;
}

// Exact integer rank of the top-left n x n block (fraction-free Gaussian
// elimination in i128 with per-row gcd reduction, so intermediates stay small).
#[allow(dead_code)]
pub fn exact_rank(m: &[[i64; 16]; 16], n: usize) -> usize {
    let mut a = [[0i128; 16]; 16];
    for i in 0..n {
        for j in 0..n {
            a[i][j] = m[i][j] as i128;
        }
    }
    let mut rank = 0usize;
    let mut prow = 0usize;
    for col in 0..n {
        let mut piv = usize::MAX;
        for r in prow..n {
            if a[r][col] != 0 {
                piv = r;
                break;
            }
        }
        if piv == usize::MAX {
            continue;
        }
        a.swap(prow, piv);
        for r in (prow + 1)..n {
            if a[r][col] != 0 {
                let (num, den) = (a[r][col], a[prow][col]);
                for c in col..n {
                    a[r][c] = a[r][c] * den - num * a[prow][c];
                }
                let mut g = 0i128;
                for c in col..n {
                    g = gcd_i128(g, a[r][c].abs());
                }
                if g > 1 {
                    for c in col..n {
                        a[r][c] /= g;
                    }
                }
            }
        }
        rank += 1;
        prow += 1;
    }
    rank
}

// does the cone { y = B^T lam : y >= 0 } contain a nonzero point? B = basis
// rows (d x n), full row rank => lambda-cone pointed => every extreme ray is
// tight on >= d-1 independent constraints: enumerate (d-1)-subsets of the n
// constraint normals (columns of B), take an integer nullspace vector of the
// subset, and test both signs against all constraints. Exact integer.
pub fn cone_has_nonneg(basis: &[[i64; 16]], n: usize, d: usize) -> bool {
    let cols: Vec<[i64; 16]> = (0..n)
        .map(|j| {
            let mut c = [0i64; 16];
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
        cols: &[[i64; 16]],
    ) -> bool {
        if pos == k {
            // nullspace vector of the chosen (d-1) columns (each a d-vector):
            // lam with lam . cols[i] = 0 for chosen i. Build (d-1) x d system,
            // take generalized cross product via cofactor expansion.
            let mut mtx = [[0i128; 16]; 16];
            for (r, &ci) in idx.iter().enumerate() {
                for c in 0..d {
                    mtx[r][c] = cols[ci][c] as i128;
                }
            }
            let mut lam = [0i128; 16];
            for c in 0..d {
                // cofactor: delete column c, det of (d-1)x(d-1), sign (-1)^c
                let mut sub = [[0i128; 16]; 16];
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
    if d == 2 {
        // 2-dim specialization: extreme rays are perpendicular to single
        // constraint normals -- lam = rot90(col). Same enumeration as rec.
        for a in cols.iter() {
            let (l0, l1) = (-(a[1] as i128), a[0] as i128);
            if l0 == 0 && l1 == 0 {
                continue;
            }
            'sg2: for sf in [1i128, -1] {
                for c in cols.iter() {
                    if sf * (l0 * c[0] as i128 + l1 * c[1] as i128) < 0 {
                        continue 'sg2;
                    }
                }
                return true;
            }
        }
        return false;
    }
    if d == 3 {
        // 3-dim specialization: an extreme ray of the (pointed) lam-cone is
        // orthogonal to 2 constraint normals -- exactly the cross product of
        // that column pair. Same enumeration as rec, no matrix scratch.
        for i in 0..n {
            let a = &cols[i];
            for b in cols.iter().skip(i + 1) {
                let l0 = a[1] as i128 * b[2] as i128 - a[2] as i128 * b[1] as i128;
                let l1 = a[2] as i128 * b[0] as i128 - a[0] as i128 * b[2] as i128;
                let l2 = a[0] as i128 * b[1] as i128 - a[1] as i128 * b[0] as i128;
                if l0 == 0 && l1 == 0 && l2 == 0 {
                    continue;
                }
                'sg: for sf in [1i128, -1] {
                    for c in cols.iter() {
                        if sf * (l0 * c[0] as i128 + l1 * c[1] as i128 + l2 * c[2] as i128) < 0 {
                            continue 'sg;
                        }
                    }
                    return true;
                }
            }
        }
        return false;
    }
    // quick accept: lam = +-e_i is a nonzero feasible point iff coordinate i
    // of every constraint column is one-signed (and not all zero) -- a genuine
    // cone point, so this can only anticipate rec's answer, never change it
    for i in 0..d {
        let mut anypos = false;
        let mut anyneg = false;
        for c in cols.iter() {
            anypos |= c[i] > 0;
            anyneg |= c[i] < 0;
        }
        if anypos != anyneg {
            return true;
        }
    }
    rec(0, d - 1, &mut idx, 0, n, d, &cols)
}

fn det_small(a: &[[i128; 16]; 16], m: usize) -> i128 {
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
        let mut sub = [[0i128; 16]; 16];
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

// ---- shared machinery for the inclusive-census family (inc10/inc4/inc_count/
// inc_strata): none of it is specific to any one stratum engine ----

pub fn factorial(n: u64) -> u128 {
    (1..=n as u128).product::<u128>().max(1)
}

pub fn lcm_to(n: u64) -> u128 {
    fn gcd(a: u64, b: u64) -> u64 {
        if b == 0 { a } else { gcd(b, a % b) }
    }
    (1..=n).fold(1u64, |l, x| l / gcd(l, x) * x) as u128
}

// fraction-free Gauss-Jordan adjugate + determinant (Bareiss divisions exact).
// adj satisfies b0 * adj = det * I; None if singular within the leading m x m
// block. `verify` re-checks that identity at runtime -- keep it on wherever the
// call is per-parent rather than per-grandparent.
pub fn adjugate_ff<const N: usize>(
    b0: &[[i128; N]; N],
    m: usize,
    verify: bool,
) -> Option<([[i128; N]; N], i128)> {
    let mut a = *b0;
    let mut aug = [[0i128; N]; N];
    for i in 0..m {
        aug[i][i] = 1;
    }
    let mut denom: i128 = 1;
    let mut sign: i128 = 1;
    for k in 0..m {
        if a[k][k] == 0 {
            let mut piv = usize::MAX;
            for r in (k + 1)..m {
                if a[r][k] != 0 {
                    piv = r;
                    break;
                }
            }
            if piv == usize::MAX {
                return None;
            }
            a.swap(k, piv);
            aug.swap(k, piv);
            sign = -sign;
        }
        let pk = a[k][k];
        for r in 0..m {
            if r == k {
                continue;
            }
            let f = a[r][k];
            for c in 0..m {
                if c >= k {
                    a[r][c] = (a[r][c] * pk - f * a[k][c]) / denom;
                }
                aug[r][c] = (aug[r][c] * pk - f * aug[k][c]) / denom;
            }
            a[r][k] = 0;
        }
        denom = pk;
    }
    let det = sign * a[m - 1][m - 1];
    let mut adj = [[0i128; N]; N];
    for i in 0..m {
        for j in 0..m {
            adj[i][j] = sign * aug[i][j];
        }
    }
    if verify {
        for i in 0..m {
            for j in 0..m {
                let mut s = 0i128;
                for k in 0..m {
                    s += b0[i][k] * adj[k][j];
                }
                assert!(s == if i == j { det } else { 0 }, "adjugate verify failed");
            }
        }
    }
    Some((adj, det))
}


// i64 fraction-free adjugate: same algorithm as adjugate_ff, all arithmetic
// in checked i64 (FF intermediates are minors of the input, Hadamard-bounded
// well inside i64 for the census matrices). Err(()) on any overflow -- the
// caller falls back to the i128 version, so this is purely a fast path.
pub fn adjugate_ff_i64(b0w: &[[i128; 16]; 16], m: usize, verify: bool) -> Result<Option<([[i128; 16]; 16], i128)>, ()> {
    let mut a = [[0i64; 16]; 16];
    for i in 0..m {
        for j in 0..m {
            if b0w[i][j].abs() > i64::MAX as i128 {
                return Err(());
            }
            a[i][j] = b0w[i][j] as i64;
        }
    }
    let b0 = a;
    let mut aug = [[0i64; 16]; 16];
    for i in 0..m {
        aug[i][i] = 1;
    }
    let mut denom: i64 = 1;
    let mut sign: i64 = 1;
    for k in 0..m {
        if a[k][k] == 0 {
            let mut piv = usize::MAX;
            for r in (k + 1)..m {
                if a[r][k] != 0 {
                    piv = r;
                    break;
                }
            }
            if piv == usize::MAX {
                return Ok(None);
            }
            a.swap(k, piv);
            aug.swap(k, piv);
            sign = -sign;
        }
        let pk = a[k][k];
        for r in 0..m {
            if r == k {
                continue;
            }
            let f = a[r][k];
            for c in 0..m {
                if c >= k {
                    let t = a[r][c].checked_mul(pk).ok_or(())?.checked_sub(f.checked_mul(a[k][c]).ok_or(())?).ok_or(())?;
                    a[r][c] = t / denom;
                }
                let t = aug[r][c].checked_mul(pk).ok_or(())?.checked_sub(f.checked_mul(aug[k][c]).ok_or(())?).ok_or(())?;
                aug[r][c] = t / denom;
            }
            a[r][k] = 0;
        }
        denom = pk;
    }
    let det = sign * a[m - 1][m - 1];
    let mut adj = [[0i128; 16]; 16];
    for i in 0..m {
        for j in 0..m {
            adj[i][j] = (sign * aug[i][j]) as i128;
        }
    }
    if verify {
        for i in 0..m {
            for j in 0..m {
                let mut s = 0i128;
                for k in 0..m {
                    s += b0[i][k] as i128 * adj[k][j];
                }
                assert!(s == if i == j { det as i128 } else { 0 }, "adjugate verify failed");
            }
        }
    }
    Ok(Some((adj, det as i128)))
}

// vertex signature inside an n-vertex beats-bitmask game: (od, id, out-
// neighbour digest, in-neighbour digest), packed comparable. An iso-invariant;
// digest collisions only weaken rigidity certificates (callers fall back to
// canon), never correctness.
pub fn vertex_sigs(beats: &[u16; 16], n: usize, sig: &mut [u64; 16]) {
    let mut od = [0u8; 16];
    let mut id = [0u8; 16];
    for v in 0..n {
        od[v] = beats[v].count_ones() as u8;
    }
    for v in 0..n {
        let mut c = 0u8;
        for u in 0..n {
            if beats[u] & (1 << v) != 0 {
                c += 1;
            }
        }
        id[v] = c;
    }
    for v in 0..n {
        // order-invariant neighbour digest: sum and sum-of-squares of
        // neighbour (od,id) codes, separated for out and in
        let mut so = 0u32;
        let mut sq = 0u32;
        let mut si = 0u32;
        let mut sqi = 0u32;
        for u in 0..n {
            let du = ((od[u] as u32) << 5) | id[u] as u32;
            if beats[v] & (1 << u) != 0 {
                so += du;
                sq += du * du;
            }
            if beats[u] & (1 << v) != 0 {
                si += du;
                sqi += du * du;
            }
        }
        sig[v] = ((od[v] as u64) << 56)
            | ((id[v] as u64) << 48)
            | ((so as u64 & 0xFFF) << 36)
            | ((sq as u64 & 0xFFF) << 24)
            | ((si as u64 & 0xFFF) << 12)
            | (sqi as u64 & 0xFFF);
    }
}

// weak connectivity of a beats-bitmask game
pub fn connected_beats(beats: &[u16; 16], n: usize) -> bool {
    let full: u16 = ((1u32 << n) - 1) as u16;
    let mut inn = [0u16; 16];
    for i in 0..n {
        let mut w = beats[i];
        while w != 0 {
            let j = w.trailing_zeros() as usize;
            w &= w - 1;
            inn[j] |= 1 << i;
        }
    }
    let mut seen = 1u16;
    let mut fr = 1u16;
    while fr != 0 {
        let mut nf = 0u16;
        let mut f = fr;
        while f != 0 {
            let v = f.trailing_zeros() as usize;
            f &= f - 1;
            nf |= (beats[v] | inn[v]) & !seen;
        }
        seen |= nf;
        fr = nf;
    }
    seen == full
}

// paradoxical (every vertex has a win and a loss) + weakly connected
pub fn paradox_connected_beats(beats: &[u16; 16], n: usize) -> bool {
    let mut inn = [0u16; 16];
    for i in 0..n {
        let mut w = beats[i];
        while w != 0 {
            let j = w.trailing_zeros() as usize;
            w &= w - 1;
            inn[j] |= 1 << i;
        }
    }
    for i in 0..n {
        if beats[i] == 0 || inn[i] == 0 {
            return false;
        }
    }
    connected_beats(beats, n)
}

// positive dependencies (support <= d+1) of the p fixed d-columns of B:
// subsets S with a positive vector alpha, sum alpha_i col_i = 0. Returned as
// (indices, alpha) with integer alpha. Minimal supports have rank |S|-1.
pub fn positive_dependencies(cols: &[[i64; 16]], p: usize, d: usize) -> Vec<(Vec<usize>, Vec<i128>)> {
    if d == 3 {
        return positive_dependencies_d3(cols, p);
    }
    let mut deps = Vec::new();
    // zero-column bitmap once (size-1 deps, and the skip test inside subsets)
    let mut zero = [false; 16];
    for i in 0..p {
        zero[i] = cols[i][..d].iter().all(|&x| x == 0);
        if zero[i] {
            deps.push((vec![i], vec![1]));
        }
    }
    // sizes 2..=d+1: alpha from cofactors of the (k-1) x d matrix of the others.
    // All scratch lives in fixed stack arrays; heap is touched only when a
    // dependency is actually found (rare), never per candidate subset.
    let mut mtx = [[0i128; 16]; 16];
    for k in 2..=(d + 1).min(p) {
        let m = k - 1;
        if m > d {
            continue;
        }
        let mut idx = vec![0usize; k];
        subsets(0, p, 0, k, &mut idx, &mut |set: &[usize]| {
            for &i in set {
                if zero[i] {
                    return;
                }
            }
            // solve sum alpha_i col_i = 0: k unknowns, d equations. minimal =>
            // rank k-1. alpha_i = (-1)^i det(matrix without row i) using any
            // (k-1)-subset of coordinates with full rank -- take all C(d, k-1)
            // coordinate subsets until a nonzero cofactor vector appears.
            let mut cidx = [0usize; 16];
            let mut alpha = [0i128; 16];
            let mut found = false;
            comb_scan(d, m, &mut cidx, &mut |cset: &[usize]| -> bool {
                let mut nonzero = false;
                for i in 0..k {
                    // det of matrix rows = set minus i, cols = cset; mtx is
                    // reused scratch -- every read cell is written first
                    let mut rr = 0;
                    for (t, &s) in set.iter().enumerate() {
                        if t == i {
                            continue;
                        }
                        for (cc, &co) in cset.iter().enumerate() {
                            mtx[rr][cc] = cols[s][co] as i128;
                        }
                        rr += 1;
                    }
                    let dt = det_n(&mtx, m);
                    alpha[i] = if i % 2 == 0 { dt } else { -dt };
                    if dt != 0 {
                        nonzero = true;
                    }
                }
                if !nonzero {
                    return false;
                }
                // verify dependency and positivity (either global sign)
                for c in 0..d {
                    let mut s = 0i128;
                    for i in 0..k {
                        s += alpha[i] * cols[set[i]][c] as i128;
                    }
                    if s != 0 {
                        return false;
                    }
                }
                let pos = alpha[..k].iter().all(|&x| x > 0);
                let neg = alpha[..k].iter().all(|&x| x < 0);
                if neg {
                    for x in alpha[..k].iter_mut() {
                        *x = -*x;
                    }
                }
                pos || neg
            }, &mut found);
            if found {
                deps.push((set.to_vec(), alpha[..k].to_vec()));
            }
        });
    }
    deps
}


// d = 3 specialization of positive_dependencies: identical output (same
// subsets, same order, same alpha construction) with the 2x2 and 3x3 basis
// minors precomputed once and all subset enumeration as direct loops.
fn positive_dependencies_d3(cols: &[[i64; 16]], p: usize) -> Vec<(Vec<usize>, Vec<i128>)> {
    let d = 3usize;
    let mut deps = Vec::new();
    let mut zero = [false; 16];
    for i in 0..p {
        zero[i] = cols[i][..d].iter().all(|&x| x == 0);
        if zero[i] {
            deps.push((vec![i], vec![1]));
        }
    }
    let pidx = |i: usize, j: usize| i * 16 + j;
    let rp_rows = [(0usize, 1usize), (0, 2), (1, 2)];
    let mut m2 = [[0i128; 256]; 3];
    for (rp, &(r0, r1)) in rp_rows.iter().enumerate() {
        for i in 0..p {
            for j in (i + 1)..p {
                m2[rp][pidx(i, j)] = cols[i][r0] as i128 * cols[j][r1] as i128
                    - cols[i][r1] as i128 * cols[j][r0] as i128;
            }
        }
    }
    let sd2 = |rp: usize, i: usize, j: usize| -> i128 {
        if i < j { m2[rp][pidx(i, j)] } else if i > j { -m2[rp][pidx(j, i)] } else { 0 }
    };
    // size 2
    for i in 0..p {
        if zero[i] {
            continue;
        }
        'sj: for j in (i + 1)..p {
            if zero[j] {
                continue;
            }
            for c in 0..d {
                let a0 = cols[j][c] as i128;
                let a1 = -(cols[i][c] as i128);
                if a0 == 0 && a1 == 0 {
                    continue;
                }
                let mut okdep = true;
                for c2 in 0..d {
                    if a0 * cols[i][c2] as i128 + a1 * cols[j][c2] as i128 != 0 {
                        okdep = false;
                        break;
                    }
                }
                if okdep {
                    let pos = a0 > 0 && a1 > 0;
                    let neg = a0 < 0 && a1 < 0;
                    if pos {
                        deps.push((vec![i, j], vec![a0, a1]));
                    } else if neg {
                        deps.push((vec![i, j], vec![-a0, -a1]));
                    }
                }
                continue 'sj; // first nonzero cofactor vector decides
            }
        }
    }
    // size 3
    for i in 0..p {
        if zero[i] {
            continue;
        }
        for j in (i + 1)..p {
            if zero[j] {
                continue;
            }
            'sl: for l in (j + 1)..p {
                if zero[l] {
                    continue;
                }
                for rp in 0..3 {
                    let a0 = sd2(rp, j, l);
                    let a1 = -sd2(rp, i, l);
                    let a2 = sd2(rp, i, j);
                    if a0 == 0 && a1 == 0 && a2 == 0 {
                        continue;
                    }
                    let mut okdep = true;
                    for c in 0..d {
                        if a0 * cols[i][c] as i128 + a1 * cols[j][c] as i128 + a2 * cols[l][c] as i128 != 0 {
                            okdep = false;
                            break;
                        }
                    }
                    if okdep {
                        let pos = a0 > 0 && a1 > 0 && a2 > 0;
                        let neg = a0 < 0 && a1 < 0 && a2 < 0;
                        if pos {
                            deps.push((vec![i, j, l], vec![a0, a1, a2]));
                        } else if neg {
                            deps.push((vec![i, j, l], vec![-a0, -a1, -a2]));
                        }
                    }
                    continue 'sl;
                }
            }
        }
    }
    // size 4: single coordinate triple; 3x3 dets via first-row m2 expansion
    let sd3 = |i: usize, j: usize, l: usize| -> i128 {
        cols[i][0] as i128 * sd2(2, j, l) - cols[i][1] as i128 * sd2(1, j, l)
            + cols[i][2] as i128 * sd2(0, j, l)
    };
    for i in 0..p {
        if zero[i] {
            continue;
        }
        for j in (i + 1)..p {
            if zero[j] {
                continue;
            }
            for l in (j + 1)..p {
                if zero[l] {
                    continue;
                }
                for q in (l + 1)..p {
                    if zero[q] {
                        continue;
                    }
                    let a0 = sd3(j, l, q);
                    let a1 = -sd3(i, l, q);
                    let a2 = sd3(i, j, q);
                    let a3 = -sd3(i, j, l);
                    if a0 == 0 && a1 == 0 && a2 == 0 && a3 == 0 {
                        continue;
                    }
                    let mut okdep = true;
                    for c in 0..d {
                        if a0 * cols[i][c] as i128 + a1 * cols[j][c] as i128 + a2 * cols[l][c] as i128 + a3 * cols[q][c] as i128 != 0 {
                            okdep = false;
                            break;
                        }
                    }
                    if !okdep {
                        continue;
                    }
                    let pos = a0 > 0 && a1 > 0 && a2 > 0 && a3 > 0;
                    let neg = a0 < 0 && a1 < 0 && a2 < 0 && a3 < 0;
                    if pos {
                        deps.push((vec![i, j, l, q], vec![a0, a1, a2, a3]));
                    } else if neg {
                        deps.push((vec![i, j, l, q], vec![-a0, -a1, -a2, -a3]));
                    }
                }
            }
        }
    }
    deps
}
// iterate k-combinations of 0..n into idx, calling f until it returns true
// (result latched into *found); allocation-free companion to `subsets`
fn comb_scan(n: usize, k: usize, idx: &mut [usize; 16], f: &mut dyn FnMut(&[usize]) -> bool, found: &mut bool) {
    fn rec(start: usize, n: usize, pos: usize, k: usize, idx: &mut [usize; 16], f: &mut dyn FnMut(&[usize]) -> bool, found: &mut bool) {
        if *found {
            return;
        }
        if pos == k {
            *found = f(&idx[..k]);
            return;
        }
        for i in start..n {
            idx[pos] = i;
            rec(i + 1, n, pos + 1, k, idx, f, found);
            if *found {
                return;
            }
        }
    }
    rec(0, n, 0, k, idx, f, found);
}

pub fn subsets(start: usize, n: usize, pos: usize, k: usize, idx: &mut Vec<usize>, f: &mut dyn FnMut(&[usize])) {
    if pos == k {
        f(&idx[..k]);
        return;
    }
    for i in start..n {
        idx[pos] = i;
        subsets(i + 1, n, pos + 1, k, idx, f);
    }
}

pub fn det_n(a: &[[i128; 16]; 16], m: usize) -> i128 {
    match m {
        0 => 1,
        1 => a[0][0],
        2 => a[0][0] * a[1][1] - a[0][1] * a[1][0],
        3 => {
            a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
                - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
                + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
        }
        _ => {
            let mut d = 0i128;
            for c in 0..m {
                if a[0][c] == 0 {
                    continue;
                }
                let mut sub = [[0i128; 16]; 16];
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
                d += s * a[0][c] * det_n(&sub, m - 1);
            }
            d
        }
    }
}


// exists strictly positive w with Mw = 0, exact: kernel basis then Gordan
// (no nonneg nonzero dependency among the kernel-basis columns). The row
// sign precheck (every row both-signed or zero) is a fast necessary filter.
pub fn has_positive_kernel(m: &[[i64; 16]; 16], n: usize) -> bool {
    for i in 0..n {
        let mut pos = false;
        let mut neg = false;
        for j in 0..n {
            pos |= m[i][j] > 0;
            neg |= m[i][j] < 0;
        }
        if pos != neg {
            return false;
        }
    }
    let basis = kernel_basis_any_fast(m, n);
    if basis.is_empty() {
        return false;
    }
    let d = basis.len();
    let cols: Vec<[i64; 16]> = (0..n)
        .map(|j| {
            let mut col = [0i64; 16];
            for (t, b) in basis.iter().enumerate() {
                col[t] = b[j];
            }
            col
        })
        .collect();
    positive_dependencies(&cols, n, d).is_empty()
}

// ---- shared bound-pruned DFS core (structure lifted from inc10.rs's
// child_dfs: L1-descending column ordering, suffix-abs-sum bounds, last-level
// determined coordinate, optional paradox-forced coordinate masks with the
// common-case split). inc10.rs keeps its own i32 specialization (hot path,
// engine-specific mid-DFS prune); inc4.rs and f3x.rs run on this one. ----

pub const DFS_RCAP: usize = 160;

// visitation order: columns sorted by descending combined L1 over the first
// `nrows` rows, so the bounds tighten fastest. Writes ord (position ->
// original column) and the permuted rows.
pub fn order_columns_l1(
    rows: &[[i64; 16]; DFS_RCAP],
    nrows: usize,
    p: usize,
    ord: &mut [usize; 16],
    prows: &mut [[i64; 16]; DFS_RCAP],
) {
    let mut l1 = [0i64; 16];
    for (c, l) in l1.iter_mut().enumerate().take(p) {
        for row in rows.iter().take(nrows) {
            *l += row[c].abs();
        }
    }
    for (t, o) in ord.iter_mut().enumerate().take(p) {
        *o = t;
    }
    ord[..p].sort_unstable_by(|&a, &b| l1[b].cmp(&l1[a]));
    for i in 0..nrows {
        for c in 0..p {
            prows[i][c] = rows[i][ord[c]];
        }
    }
}

// asum[i][k] = sum of |row i| over columns k.. (suffix bound seeds)
pub fn suffix_abs_sums(
    prows: &[[i64; 16]; DFS_RCAP],
    nrows: usize,
    p: usize,
    asum: &mut [[i64; 17]; DFS_RCAP],
) {
    for i in 0..nrows {
        asum[i][p] = 0;
        for k in (0..p).rev() {
            asum[i][k] = asum[i][k + 1] + prows[i][k].abs();
        }
    }
}

// bound-pruned DFS over r in {-1,0,+1}^p. Rows 0..ne are EQUALITIES (dot r
// must be exactly 0 at a leaf; pruned by |s| <= suffix bound), rows ne..nt
// are STRICT (dot r < 0 at a leaf; pruned by reachability). fplus/fminus are
// forced-coordinate masks in ORDERED positions (bit k set: r[k] must be
// +1 resp. -1). The last level is solved directly from the first equality
// row with a nonzero final coefficient. The leaf callback receives r in
// ordered coordinates (invert with ord[] as needed).
#[allow(clippy::too_many_arguments)]
pub fn dfs_es(
    k: usize,
    p: usize,
    ne: usize,
    nt: usize,
    s: &mut [i64; DFS_RCAP],
    r: &mut [i32; 16],
    rows: &[[i64; 16]; DFS_RCAP],
    asum: &[[i64; 17]; DFS_RCAP],
    fplus: u16,
    fminus: u16,
    leaf: &mut impl FnMut(&[i32; 16]),
) {
    for i in 0..ne {
        if s[i].abs() > asum[i][k] {
            return;
        }
    }
    for i in ne..nt {
        if s[i] - asum[i][k] >= 0 {
            return;
        }
    }
    if k == p {
        leaf(r);
        return;
    }
    if k + 1 == p {
        let fp = (fplus >> k) & 1 != 0;
        let fm = (fminus >> k) & 1 != 0;
        let mut e0 = usize::MAX;
        for e in 0..ne {
            if rows[e][k] != 0 {
                e0 = e;
                break;
            }
        }
        if e0 != usize::MAX {
            // determined final coordinate
            let c0 = rows[e0][k];
            if s[e0] % c0 != 0 {
                return;
            }
            let v = -s[e0] / c0;
            if !(-1..=1).contains(&v) || (fp && v != 1) || (fm && v != -1) {
                return;
            }
            let val = v as i32;
            r[k] = val;
            if val != 0 {
                for i in 0..nt {
                    s[i] += val as i64 * rows[i][k];
                }
            }
            let mut ok = true;
            for i in 0..ne {
                if s[i] != 0 {
                    ok = false;
                    break;
                }
            }
            if ok {
                for i in ne..nt {
                    if s[i] >= 0 {
                        ok = false;
                        break;
                    }
                }
            }
            if ok {
                leaf(r);
            }
            if val != 0 {
                for i in 0..nt {
                    s[i] -= val as i64 * rows[i][k];
                }
            }
            r[k] = 0;
            return;
        }
        // no equality touches this column: equality sums are already final
        for i in 0..ne {
            if s[i] != 0 {
                return;
            }
        }
        for val in [0i32, -1, 1] {
            if (fp && val != 1) || (fm && val != -1) {
                continue;
            }
            r[k] = val;
            if val != 0 {
                for i in 0..nt {
                    s[i] += val as i64 * rows[i][k];
                }
            }
            let mut ok = true;
            for i in ne..nt {
                if s[i] >= 0 {
                    ok = false;
                    break;
                }
            }
            if ok {
                leaf(r);
            }
            if val != 0 {
                for i in 0..nt {
                    s[i] -= val as i64 * rows[i][k];
                }
            }
        }
        r[k] = 0;
        return;
    }
    if fplus | fminus == 0 {
        // common case: no forced coordinates anywhere below
        for val in [0i32, -1, 1] {
            r[k] = val;
            if val != 0 {
                for i in 0..nt {
                    s[i] += val as i64 * rows[i][k];
                }
            }
            dfs_es(k + 1, p, ne, nt, s, r, rows, asum, 0, 0, leaf);
            if val != 0 {
                for i in 0..nt {
                    s[i] -= val as i64 * rows[i][k];
                }
            }
        }
        r[k] = 0;
        return;
    }
    for val in [0i32, -1, 1] {
        if (fplus >> k) & 1 != 0 && val != 1 {
            continue;
        }
        if (fminus >> k) & 1 != 0 && val != -1 {
            continue;
        }
        r[k] = val;
        if val != 0 {
            for i in 0..nt {
                s[i] += val as i64 * rows[i][k];
            }
        }
        dfs_es(k + 1, p, ne, nt, s, r, rows, asum, fplus, fminus, leaf);
        if val != 0 {
            for i in 0..nt {
                s[i] -= val as i64 * rows[i][k];
            }
        }
    }
    r[k] = 0;
}
