// Enumerate all INCLUSIVE games on n (odd) vertices by EXTENDING the (n-1)-vertex
// oriented parents -- avoiding the ~A001174(n) directg-n scan entirely, the same
// trick that made cm_extend fast, now generalized past the completely-mixed
// (nullity-1) case. Reads digraph6 (n-1)-vertex oriented graphs on stdin, writes
// digraph6 of every inclusive n-vertex child; pipe through `nauty-labelg | sort -u
// | wc -l` to count iso classes.
//
//   rustc -O rust/inc_extend.rs -o /tmp/incx
//   nauty-geng 6 | nauty-directg -o | /tmp/incx 7 | nauty-labelg | sort -u | wc -l   # 10525
//
// The extension lemma (validated 2002/2002 on real games, and provable): a parent
// M' with kernel K of dim d, extended by a new-vertex vector r, yields a child of
// nullity  d+1 if r _|_ K,  else d-1.  Two cases:
//
//   * nonsingular parent (d=0): child nullity 1 for every r, kernel (-M'^-1 r, 1);
//     inclusive iff -M'^-1 r > 0 -- the completely-mixed cone, found by the same
//     bound-pruned DFS + canonical-deletion prefilter as cm_extend.
//   * singular parent (d>=2): the nullity-(d+1) children are exactly those with
//     r _|_ K (K r = 0, d linear constraints -- a hard prune on the 3^(n-1) r-scan).
//     Every nullity-(d+1) game has a nullity-d parent (a rank-preserving deletion),
//     so extending every even-nullity parent by r _|_ K reaches all nullity>=3
//     games; a Phase-1 LP then keeps the ones whose >=3-dim kernel meets the
//     positive orthant. These are emitted WITHOUT the maximal prefilter (its
//     canonical parent may sit at a different nullity level), so the final
//     `sort -u` does their dedup.
//
// inverse()/dfs()/added_is_maximal()/encode()/decode()/fully_mixed() live in
// common.rs (hoisted verbatim from cm_extend.rs / inc_fast.rs, so the CM half is
// byte-for-byte the validated tool); nullspace() and the singular branch are the
// new part.

use std::env;
use std::io::{self, BufWriter, Read, Write};

mod common;
use common::{added_is_maximal, decode, dfs, dfs_g, encode, fully_mixed, inverse, sig_maximal};

// Bound-pruned DFS collecting r in {-1,0,1}^p with N r = 0 (r _|_ K), N = kernel
// basis (d rows, columns reordered so the DFS assigns high-L1 columns first).
// `s[i]` is the running i-th constraint N_i . r; `asum[i][k] = sum_{j>=k}|Nc[i][j]|`
// bounds how far it can still move, so if |s[i]| > asum[i][k] the constraint can no
// longer reach 0 -- prune. Turns the 3^p flat r-scan into a walk of a small live
// subtree, exactly analogous to the completely-mixed cone-DFS.
#[allow(clippy::too_many_arguments)]
fn dfs_perp(
    k: usize,
    p: usize,
    d: usize,
    s: &mut [f64; 8],
    r: &mut [i32; 8],
    nc: &[[f64; 8]; 8],
    asum: &[[f64; 9]; 8],
    out: &mut Vec<[i32; 8]>,
) {
    for i in 0..d {
        if s[i].abs() > asum[i][k] + 1e-6 {
            return; // constraint i can no longer reach 0
        }
    }
    if k == p {
        out.push(*r);
        return;
    }
    for val in [-1i32, 0, 1] {
        r[k] = val;
        if val != 0 {
            let f = val as f64;
            for i in 0..d {
                s[i] += f * nc[i][k];
            }
        }
        dfs_perp(k + 1, p, d, s, r, nc, asum, out);
        if val != 0 {
            let f = val as f64;
            for i in 0..d {
                s[i] -= f * nc[i][k];
            }
        }
    }
    r[k] = 0;
}

// ---- kernel basis via RREF: rows of `basis[..d]` span ker(a) (p x p) ----
fn nullspace(a: &[[f64; 8]; 8], p: usize) -> (usize, [[f64; 8]; 8]) {
    let mut m = *a;
    let mut pivot_col = [usize::MAX; 8];
    let mut is_pivot = [false; 8];
    let mut row = 0usize;
    for col in 0..p {
        if row >= p {
            break;
        }
        let mut piv = row;
        for r in row..p {
            if m[r][col].abs() > m[piv][col].abs() {
                piv = r;
            }
        }
        if m[piv][col].abs() < 1e-9 {
            continue; // free column
        }
        m.swap(row, piv);
        let d = 1.0 / m[row][col];
        for k in 0..p {
            m[row][k] *= d;
        }
        for r in 0..p {
            if r != row {
                let f = m[r][col];
                if f != 0.0 {
                    for k in 0..p {
                        m[r][k] -= f * m[row][k];
                    }
                }
            }
        }
        is_pivot[col] = true;
        pivot_col[row] = col;
        row += 1;
    }
    let rank = row;
    let mut basis = [[0f64; 8]; 8];
    let mut bi = 0usize;
    for col in 0..p {
        if !is_pivot[col] {
            basis[bi][col] = 1.0;
            for (r, &pc) in pivot_col.iter().enumerate().take(rank) {
                basis[bi][pc] = -m[r][col];
            }
            bi += 1;
        }
    }
    (p - rank, basis)
}

// rank of the skew matrix `m` (n x n) with row/col `skip` deleted (size n-1),
// via float Gaussian elimination counting pivots.
fn rank_deleting(m: &[[f64; 9]; 9], n: usize, skip: usize) -> usize {
    let sz = n - 1;
    let mut a = [[0f64; 9]; 9];
    let mut ri = 0;
    for i in 0..n {
        if i == skip {
            continue;
        }
        let mut ci = 0;
        for j in 0..n {
            if j == skip {
                continue;
            }
            a[ri][ci] = m[i][j];
            ci += 1;
        }
        ri += 1;
    }
    let mut rank = 0usize;
    let mut prow = 0usize;
    for col in 0..sz {
        let mut piv = usize::MAX;
        let mut best = 1e-9;
        for r in prow..sz {
            if a[r][col].abs() > best {
                best = a[r][col].abs();
                piv = r;
            }
        }
        if piv == usize::MAX {
            continue;
        }
        a.swap(prow, piv);
        for r in (prow + 1)..sz {
            let f = a[r][col] / a[prow][col];
            if f != 0.0 {
                for c in col..sz {
                    a[r][c] -= f * a[prow][c];
                }
            }
        }
        rank += 1;
        prow += 1;
    }
    rank
}

// Canonical-deletion prefilter for the nullity>=3 emissions. A nullity-(d+1) child
// is built by adding vertex p to a nullity-d parent (p is always a nullity-d
// deletion, i.e. nullity(M_p)=d). Every nullity-(d+1) game has such a nullity-d
// deletion, so its signature-maximal one is a well-defined canonical parent; we
// emit only when p IS that vertex. This drops the emission from
// answer*redundancy down to ~answer, which the sparse parents (2.4 emissions each)
// would otherwise blow up. Sound: p in the candidate set, so we reject only when a
// strictly bigger nullity-d vertex exists; ties fall through and `sort -u` cleans
// them. Signature = the degree-refinement key used by added_is_maximal.
fn nullity3_canonical(beats: &[u16], n: usize, p: usize, d: usize) -> bool {
    let mut m = [[0f64; 9]; 9];
    for i in 0..n {
        let mut w = beats[i];
        while w != 0 {
            let j = w.trailing_zeros() as usize;
            w &= w - 1;
            m[i][j] = 1.0;
            m[j][i] = -1.0;
        }
    }
    // p competes only against the other nullity-d deletions (see doc above);
    // p itself qualifies by construction.
    sig_maximal(beats, n, p, |v| (n - 1) - rank_deleting(&m, n, v) == d)
}

// ---- kernel-coordinate fully-mixed test (the LP replacement) ----
//
// For a singular parent M' with kernel basis B (d rows) and extension vector q
// (q _|_ K), the child kernel is exactly { (B^T lam + t*w, t) } where M'w = -q:
// the last row q.x = 0 holds automatically (q _|_ K kills B^T lam, and
// q.w = -(M'w).w = 0 by skewness). A strictly positive kernel vector must have
// t > 0 (its last coordinate IS t), so scaling t = 1:
//
//     child fully mixed  <=>  exists lam in R^d :  B^T lam + w > 0.
//
// That is a d-dimensional feasibility problem -- d = 2 almost always -- against
// p = 8 half-planes whose NORMALS a_j = (B[0][j], B[1][j]) are fixed per parent;
// only the offset w changes with the candidate q. By Motzkin transposition the
// system { a_j . lam > -w_j } is infeasible iff some alpha >= 0 (not all 0) has
// sum alpha_j a_j = 0 and sum alpha_j w_j <= 0, and by Caratheodory a minimal
// such alpha is supported on <= 3 normals in R^2. So we precompute, once per
// parent, every positive dependency of <= 3 normals (zero normals; antiparallel
// pairs; positively-spanning triples), and each candidate is then one O(p^2)
// solve for w plus a dot product per dependency -- vs. the full n x (2n+1)
// phase-1 simplex tableau this replaces (~10x the flops). d >= 4 parents (rare)
// still take the LP fallback.
//
// w comes from the bordered system [[M', B^T],[B, 0]] (w, mu) = (-q, 0), whose
// matrix is nonsingular exactly because B spans ker(M'); its inverse is computed
// once per parent and w = -(upper-left p x p block) q per candidate.

// Gauss-Jordan inverse of the (p+d) x (p+d) bordered matrix (p+d <= 12).
fn bordered_inverse(mp: &[[f64; 8]; 8], nb: &[[f64; 8]; 8], p: usize, d: usize) -> Option<[[f64; 12]; 12]> {
    let m = p + d;
    let mut a = [[0f64; 12]; 12];
    let mut inv = [[0f64; 12]; 12];
    for i in 0..p {
        for j in 0..p {
            a[i][j] = mp[i][j];
        }
        for t in 0..d {
            a[i][p + t] = nb[t][i];
            a[p + t][i] = nb[t][i];
        }
    }
    for i in 0..m {
        inv[i][i] = 1.0;
    }
    for col in 0..m {
        let mut piv = col;
        for r in (col + 1)..m {
            if a[r][col].abs() > a[piv][col].abs() {
                piv = r;
            }
        }
        if a[piv][col].abs() < 1e-9 {
            return None; // numerically singular: caller falls back to the LP
        }
        a.swap(col, piv);
        inv.swap(col, piv);
        let dv = 1.0 / a[col][col];
        for k in 0..m {
            a[col][k] *= dv;
            inv[col][k] *= dv;
        }
        for r in 0..m {
            if r != col {
                let f = a[r][col];
                if f != 0.0 {
                    for k in (col + 1)..m {
                        a[r][k] -= f * a[col][k];
                    }
                    for k in 0..m {
                        inv[r][k] -= f * inv[col][k];
                    }
                }
            }
        }
    }
    Some(inv)
}

// A dependency (i, j, k, ai, aj, ak): a positive combination of <= 3 of the
// fixed 2-D normals that sums to zero. j == usize::MAX marks a zero normal
// (constraint reduces to w_i > 0); k == usize::MAX marks an antiparallel pair.
// The candidate is feasible iff EVERY dependency has sum alpha_j w_j > 0.
type Dep = (usize, usize, usize, f64, f64, f64);

fn helly2_deps(nb: &[[f64; 8]; 8], p: usize, deps: &mut Vec<Dep>) {
    deps.clear();
    let (ax, ay) = (&nb[0], &nb[1]);
    let norm = |j: usize| (ax[j] * ax[j] + ay[j] * ay[j]).sqrt();
    const ZT: f64 = 1e-9;
    for j in 0..p {
        if norm(j) < ZT {
            deps.push((j, usize::MAX, 0, 1.0, 0.0, 0.0));
        }
    }
    for i in 0..p {
        let ni = norm(i);
        if ni < ZT {
            continue;
        }
        for j in (i + 1)..p {
            let nj = norm(j);
            if nj < ZT {
                continue;
            }
            let cross = ax[i] * ay[j] - ay[i] * ax[j];
            if cross.abs() < 1e-7 * ni * nj {
                // parallel: a genuine dependency only if opposed (alpha = norms)
                if ax[i] * ax[j] + ay[i] * ay[j] < 0.0 {
                    deps.push((i, j, usize::MAX, nj, ni, 0.0));
                }
                continue; // triples over a parallel pair reduce to pair/single deps
            }
            for k in (j + 1)..p {
                if norm(k) < ZT {
                    continue;
                }
                // solve ai*a_i + aj*a_j = -a_k (Cramer); positive pair => 0 is in
                // the open positive hull of the triple
                let bi = (ay[k] * ax[j] - ax[k] * ay[j]) / cross;
                let bj = (ay[i] * ax[k] - ax[i] * ay[k]) / cross;
                if bi > 1e-9 && bj > 1e-9 {
                    deps.push((i, j, k, bi, bj, 1.0));
                }
            }
        }
    }
}

// Fused DFS for the d=2 singular parents: walk q in {-1,0,1}^p pruning on BOTH
// the kernel equalities (nc rows 0..ne: N q = 0, the nullity condition) and the
// fully-mixed inequalities (rows ne..nt: g q < -TOL). The inequalities are the
// Helly dependencies pushed through w = -Wq: dependency alpha is satisfied iff
// sum alpha_j w_j > TOL <=> (alpha^T W) q < -TOL, one fixed row g = alpha^T W per
// dependency. So the tree only ever reaches leaves that are already nullity-3
// AND fully mixed -- the ~100x-larger set of merely-perpendicular q is never
// enumerated, which is where the flat dfs_perp + per-candidate test spent its
// time. Equality rows prune two-sided (|s| can no longer reach 0), inequality
// rows one-sided (s - asum >= -TOL can no longer dip below), exactly the cm
// cone-DFS bound.
const MAXROWS: usize = 96;
#[allow(clippy::too_many_arguments)]
fn dfs_fused(
    k: usize,
    p: usize,
    ne: usize,
    nt: usize,
    s: &mut [f64; MAXROWS],
    r: &mut [i32; 8],
    rows: &[[f64; 8]; MAXROWS],
    asum: &[[f64; 9]; MAXROWS],
    out: &mut Vec<[i32; 8]>,
) {
    for i in 0..ne {
        if s[i].abs() > asum[i][k] + 1e-6 {
            return; // equality row i can no longer reach 0
        }
    }
    for i in ne..nt {
        if s[i] - asum[i][k] >= -1e-7 {
            return; // inequality row i can no longer dip below -TOL
        }
    }
    if k == p {
        out.push(*r);
        return;
    }
    for val in [0i32, -1, 1] {
        r[k] = val;
        if val != 0 {
            let f = val as f64;
            for i in 0..nt {
                s[i] += f * rows[i][k];
            }
        }
        dfs_fused(k + 1, p, ne, nt, s, r, rows, asum, out);
        if val != 0 {
            let f = val as f64;
            for i in 0..nt {
                s[i] -= f * rows[i][k];
            }
        }
    }
    r[k] = 0;
}

// paradoxical (every vertex has a win and a loss) + connected on decisive edges.
// Required for INCLUSIVE: unlike nullity 1, a nullity>=3 game can be fully-mixed
// yet disconnected (an all-tie/isolated vertex contributes e_i to the kernel), so
// fully_mixed alone doesn't imply it -- the singular branch must check explicitly.
fn paradox_connected(beats: &[u16], n: usize) -> bool {
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
            return false; // a vertex with no win or no loss
        }
    }
    let full: u16 = ((1u32 << n) - 1) as u16;
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

fn main() {
    let n: usize = env::args().nth(1).and_then(|s| s.parse().ok()).expect("usage: incx n");
    assert!(n >= 3 && n % 2 == 1 && n <= 9, "n must be odd, 3..=9 (8-wide fixed arrays)");
    let p = n - 1;
    let preclen = 2 + (p * p + 5) / 6 + 1;

    // optional 2nd arg: write the nullity>=3 children to this file instead of
    // stdout. CM (nullity 1) and nullity>=3 children are disjoint iso classes,
    // so the two streams dedup independently: the CM side must reproduce the
    // cm_extend census exactly (n=9 -> 583591020), a built-in full-run checksum.
    let mut hiw = env::args().nth(2).map(|pth| {
        BufWriter::with_capacity(1 << 20, std::fs::File::create(pth).expect("create hi-out file"))
    });
    let mut stdin = io::stdin().lock();
    let mut inbuf = vec![0u8; 1 << 20];
    let mut have = 0usize;
    let out = io::stdout();
    let mut w = BufWriter::with_capacity(1 << 20, out.lock());
    let mut child = Vec::with_capacity(2 + (n * n + 5) / 6 + 1);

    let mut mp = [[0f64; 8]; 8];
    let mut r = [0i32; 8];
    let mut valid: Vec<[i32; 8]> = Vec::with_capacity(256);
    let mut deps: Vec<Dep> = Vec::with_capacity(64);
    // fused-DFS constraint rows + suffix bounds, reused across parents (only the
    // first nt rows are ever written/read for a given parent)
    let mut rows = [[0f64; 8]; MAXROWS];
    let mut fasum = [[0f64; 9]; MAXROWS];
    let (mut parents, mut emit_cm, mut emit_hi) = (0u64, 0u64, 0u64);
    let (mut sing_d2, mut sing_fb) = (0u64, 0u64);

    loop {
        let got = stdin.read(&mut inbuf[have..]).unwrap();
        if got == 0 {
            break;
        }
        have += got;
        let nrec = have / preclen;
        for ri in 0..nrec {
            let rec = &inbuf[ri * preclen..(ri + 1) * preclen];
            if rec[0] != b'&' {
                continue;
            }
            parents += 1;
            decode(rec, p, &mut mp);
            // parent arc bitmask (shared by every child of this parent)
            let mut base = [0u16; 16];
            for i in 0..p {
                for j in 0..p {
                    if mp[i][j] > 0.5 {
                        base[i] |= 1 << j;
                    }
                }
            }

            match inverse(&mp, p) {
                Some(inv) => {
                    // nonsingular: completely-mixed children via the cone-DFS (cm_extend)
                    let mut ord = [0usize; 8];
                    for (k, o) in ord[..p].iter_mut().enumerate() {
                        *o = k;
                    }
                    let mut l1 = [0f64; 8];
                    for k in 0..p {
                        let mut sum = 0f64;
                        for i in 0..p {
                            sum += inv[i][k].abs();
                        }
                        l1[k] = sum;
                    }
                    ord[..p].sort_unstable_by(|&a, &b| l1[b].partial_cmp(&l1[a]).unwrap());
                    let mut wc = [[0f64; 8]; 8];
                    for k in 0..p {
                        for i in 0..p {
                            wc[i][k] = inv[i][ord[k]];
                        }
                    }
                    let mut asum = [[0f64; 8]; 9];
                    for k in (0..p).rev() {
                        for i in 0..p {
                            asum[k][i] = asum[k + 1][i] + wc[i][k].abs();
                        }
                    }
                    valid.clear();
                    let mut s0 = [0f64; 8];
                    match p {
                        8 => dfs_g::<8>(0, &mut s0, &mut r, &wc, &asum, &mut valid),
                        6 => dfs_g::<6>(0, &mut s0, &mut r, &wc, &asum, &mut valid),
                        4 => dfs_g::<4>(0, &mut s0, &mut r, &wc, &asum, &mut valid),
                        _ => dfs(0, p, &mut s0, &mut r, &wc, &asum, &mut valid),
                    }
                    if valid.is_empty() {
                        continue;
                    }
                    // parent degrees, hoisted out of the candidate loop: the child's
                    // degrees are od[c] + [c beats new], id[c] + [new beats c], so
                    // the (od, id) first stage of the maximality signature costs O(p)
                    // per candidate instead of the O(n^2) recomputation inside
                    // added_is_maximal -- which then only runs on exact (od, id) ties.
                    let mut pod = [0u8; 8];
                    let mut pid = [0u8; 8];
                    for v in 0..p {
                        pod[v] = base[v].count_ones() as u8;
                        let mut cnt = 0u8;
                        for u in 0..p {
                            if base[u] & (1 << v) != 0 {
                                cnt += 1;
                            }
                        }
                        pid[v] = cnt;
                    }
                    for rv in valid.iter() {
                        let (mut odp, mut idp) = (0u8, 0u8);
                        for k in 0..p {
                            if rv[k] < 0 {
                                odp += 1;
                            } else if rv[k] > 0 {
                                idp += 1;
                            }
                        }
                        let mut reject = false;
                        let mut tie = false;
                        for k in 0..p {
                            let c = ord[k];
                            let odv = pod[c] + (rv[k] > 0) as u8;
                            let idv = pid[c] + (rv[k] < 0) as u8;
                            if odv > odp || (odv == odp && idv > idp) {
                                reject = true; // a strictly larger vertex exists
                                break;
                            }
                            if odv == odp && idv == idp {
                                tie = true; // needs the full signature comparison
                            }
                        }
                        if reject {
                            continue;
                        }
                        let mut beats = base;
                        for k in 0..p {
                            let c = ord[k];
                            if rv[k] > 0 {
                                beats[c] |= 1 << p;
                            } else if rv[k] < 0 {
                                beats[p] |= 1 << c;
                            }
                        }
                        if tie && !added_is_maximal(&beats, n, p) {
                            continue;
                        }
                        encode(&beats[..n], n, &mut child);
                        w.write_all(&child).unwrap();
                        emit_cm += 1;
                    }
                }
                None => {
                    // singular: nullity>=3 children via r _|_ K, kept fully mixed.
                    let (d, nb) = nullspace(&mp, p);
                    if d == 0 {
                        continue; // numerical: treat as nonsingular-but-uninvertible, skip
                    }
                    // reorder kernel-basis columns by decreasing L1 for stronger pruning
                    let mut ord = [0usize; 8];
                    for (k, o) in ord[..p].iter_mut().enumerate() {
                        *o = k;
                    }
                    let mut l1 = [0f64; 8];
                    for k in 0..p {
                        let mut sum = 0f64;
                        for bv in nb.iter().take(d) {
                            sum += bv[k].abs();
                        }
                        l1[k] = sum;
                    }
                    ord[..p].sort_unstable_by(|&a, &b| l1[b].partial_cmp(&l1[a]).unwrap());

                    // d=2 (the common case): fused DFS over kernel equalities AND
                    // the Helly fully-mixed inequalities -- every leaf is already
                    // nullity-3 and fully mixed, so no per-candidate LP at all.
                    let mut fused = false;
                    if d == 2 {
                        if let Some(binv) = bordered_inverse(&mp, &nb, p, d) {
                            helly2_deps(&nb, p, &mut deps);
                            if 2 + deps.len() <= MAXROWS {
                                fused = true;
                                sing_d2 += 1;
                                let nt = 2 + deps.len();
                                for k in 0..p {
                                    for i in 0..2 {
                                        rows[i][k] = nb[i][ord[k]];
                                    }
                                }
                                for (t, &(i, j, kk, ai, aj, ak)) in deps.iter().enumerate() {
                                    // g = alpha^T W in original columns (W = upper-left
                                    // p x p block of the bordered inverse), reordered
                                    let mut g = [0f64; 8];
                                    for (c, gc) in g.iter_mut().enumerate().take(p) {
                                        *gc = if j == usize::MAX {
                                            binv[i][c]
                                        } else if kk == usize::MAX {
                                            ai * binv[i][c] + aj * binv[j][c]
                                        } else {
                                            ai * binv[i][c] + aj * binv[j][c] + ak * binv[kk][c]
                                        };
                                    }
                                    for k in 0..p {
                                        rows[2 + t][k] = g[ord[k]];
                                    }
                                }
                                for i in 0..nt {
                                    for k in (0..p).rev() {
                                        fasum[i][k] = fasum[i][k + 1] + rows[i][k].abs();
                                    }
                                }
                                valid.clear();
                                let mut s0 = [0f64; MAXROWS];
                                dfs_fused(0, p, 2, nt, &mut s0, &mut r, &rows, &fasum, &mut valid);
                                for rv in valid.iter() {
                                    let mut beats = base;
                                    for k in 0..p {
                                        let c = ord[k];
                                        if rv[k] > 0 {
                                            beats[c] |= 1 << p;
                                        } else if rv[k] < 0 {
                                            beats[p] |= 1 << c;
                                        }
                                    }
                                    if !paradox_connected(&beats[..n], n) {
                                        continue;
                                    }
                                    // canonical-deletion prefilter, kernel-coordinate
                                    // form: deleting v from a skew M drops the nullity
                                    // by 1 iff some kernel vector has a nonzero v-th
                                    // coordinate (if all kernel vectors vanish at v,
                                    // then e_v is in im(M) and its preimage w also has
                                    // w_v = w.(Mw) = 0, so the nullity RISES by 1). The
                                    // child kernel basis is {(B_1,0),(B_2,0),(w,1)} with
                                    // w = -W q already available -- so the nullity-d
                                    // deletions are just the v where (B_1,B_2,w)_v != 0,
                                    // no rank computations at all. The added vertex p is
                                    // always eligible (its coordinate is t = 1).
                                    let mut wv = [0f64; 8];
                                    for k in 0..p {
                                        if rv[k] != 0 {
                                            let f = rv[k] as f64;
                                            let c = ord[k];
                                            for i in 0..p {
                                                wv[i] -= f * binv[i][c];
                                            }
                                        }
                                    }
                                    let elig = |v: usize| {
                                        nb[0][v].abs() > 1e-7
                                            || nb[1][v].abs() > 1e-7
                                            || wv[v].abs() > 1e-7
                                    };
                                    if sig_maximal(&beats[..n], n, p, elig) {
                                        encode(&beats[..n], n, &mut child);
                                        match hiw.as_mut() {
                                            Some(f) => f.write_all(&child).unwrap(),
                                            None => w.write_all(&child).unwrap(),
                                        }
                                        emit_hi += 1;
                                    }
                                }
                            }
                        }
                    }
                    if fused {
                        continue;
                    }
                    sing_fb += 1;

                    // fallback (d >= 4, or a numerically-singular bordered matrix):
                    // flat perpendicularity DFS + per-candidate LP, as validated.
                    let mut nc = [[0f64; 8]; 8];
                    for k in 0..p {
                        for i in 0..d {
                            nc[i][k] = nb[i][ord[k]];
                        }
                    }
                    let mut asum = [[0f64; 9]; 8];
                    for k in (0..p).rev() {
                        for i in 0..d {
                            asum[i][k] = asum[i][k + 1] + nc[i][k].abs();
                        }
                    }
                    valid.clear();
                    let mut s0 = [0f64; 8];
                    dfs_perp(0, p, d, &mut s0, &mut r, &nc, &asum, &mut valid);
                    for rv in valid.iter() {
                        let mut beats = base;
                        for k in 0..p {
                            let c = ord[k];
                            if rv[k] > 0 {
                                beats[c] |= 1 << p;
                            } else if rv[k] < 0 {
                                beats[p] |= 1 << c;
                            }
                        }
                        if paradox_connected(&beats[..n], n)
                            && fully_mixed(&beats[..n], n)
                            && nullity3_canonical(&beats[..n], n, p, d)
                        {
                            encode(&beats[..n], n, &mut child);
                            match hiw.as_mut() {
                                Some(f) => f.write_all(&child).unwrap(),
                                None => w.write_all(&child).unwrap(),
                            }
                            emit_hi += 1;
                        }
                    }
                }
            }
        }
        let rem = have - nrec * preclen;
        inbuf.copy_within(nrec * preclen..have, 0);
        have = rem;
    }
    w.flush().unwrap();
    if let Some(f) = hiw.as_mut() {
        f.flush().unwrap();
    }
    eprintln!(
        "incx n={}: parents={} (d2-fused={} fallback={}) emitted cm(nullity1)={} nullity>=3={} (pipe: labelg|sort -u|wc -l)",
        n, parents, sing_d2, sing_fb, emit_cm, emit_hi
    );
}
