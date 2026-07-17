// Nullity-4 stratum of inclusive(even n): labeled count via the same
// deletion-multiplicity identity as inc10.rs, one nullity level up.
//
// Parents: (n-1)-vertex games with nullity EXACTLY 3 and a nonnegative kernel
// vector (family F3', emitted by `inc_strata f3-emit 3`). For a child C
// (nullity 4, inclusive), the valid deletions are the vertices v whose kernel
// slice K_C /\ {x_v = 0} (3-dim) contains a nonzero nonnegative vector; the
// facet structure of O = K /\ Delta guarantees at least one. Hence
//
//   L_4 = n * sum over labeled F3' parents P, valid r:  1 / z4(child)
//
// Per parent, exact integer arithmetic throughout:
//   * kernel basis B (3 x p) via fraction-free RREF (common::kernel_basis_exact)
//   * bordered adjugate of [[M', B^T],[B, 0]] gives w = -(ADJ r)/D
//   * child inclusive <=> exists lam: B^T lam + w > 0; by Motzkin the failures
//     are witnessed by positive dependencies of <= 4 kernel columns -- all
//     enumerated once per parent, each linear in r, FUSED into the r-DFS as
//     strict rows alongside the three equality rows (r _|_ each basis row)
//   * per leaf: child paradox+connected, then z4 = #{v: slice-cone nonzero}
//     via extreme-ray enumeration in the equality-projected 3-var cone.
//
// Anchor (must be exact): n=8 with parents from directg-7 f3-emit:
//   L_4 = 880,869,360.
//
//   rustc -O rust/inc4.rs -o /tmp/inc4 -C link-args="shim.o -lnauty"
//   nauty-geng 7 | nauty-directg -o | /tmp/incs 7 f3-emit 3 | /tmp/inc4 8
use std::env;
use std::io::{self, Read};
use std::os::raw::c_int;

mod common;
use common::{cone_has_nonneg, kernel_basis_exact};

extern "C" {
    fn rps_autsize(arc: *const u64, n: c_int) -> f64;
}

fn factorial(n: u64) -> u128 {
    (1..=n as u128).product::<u128>().max(1)
}

fn lcm_to(n: u64) -> u128 {
    fn gcd(a: u64, b: u64) -> u64 {
        if b == 0 { a } else { gcd(b, a % b) }
    }
    (1..=n).fold(1u64, |l, x| l / gcd(l, x) * x) as u128
}

fn gcd128(a: i128, b: i128) -> i128 {
    if b == 0 { a } else { gcd128(b, a % b) }
}

// fraction-free Gauss-Jordan adjugate for m <= 12 (i128; bordered matrices)
fn adjugate12(b0: &[[i128; 12]; 12], m: usize) -> ([[i128; 12]; 12], i128) {
    let mut a = *b0;
    let mut aug = [[0i128; 12]; 12];
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
            assert!(piv != usize::MAX, "bordered singular");
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
    let mut adj = [[0i128; 12]; 12];
    for i in 0..m {
        for j in 0..m {
            adj[i][j] = sign * aug[i][j];
        }
    }
    // verify B*adj = det*I (cheap relative to the DFS)
    for i in 0..m {
        for j in 0..m {
            let mut s = 0i128;
            for k in 0..m {
                s += b0[i][k] * adj[k][j];
            }
            assert!(s == if i == j { det } else { 0 }, "adjugate verify failed");
        }
    }
    (adj, det)
}

// positive dependencies (support <= d+1) of the p fixed d-columns of B:
// subsets S with a positive vector alpha, sum alpha_i col_i = 0. Returned as
// (indices, alpha) with integer alpha. Minimal supports have rank |S|-1.
fn positive_dependencies(cols: &[[i64; 4]], p: usize, d: usize) -> Vec<(Vec<usize>, Vec<i128>)> {
    let mut deps = Vec::new();
    // size 1: zero columns
    for i in 0..p {
        if cols[i][..d].iter().all(|&x| x == 0) {
            deps.push((vec![i], vec![1]));
        }
    }
    // sizes 2..=d+1: alpha from cofactors of the (k-1) x d matrix of the others
    let sizes: Vec<usize> = (2..=(d + 1)).collect();
    for &k in &sizes {
        let mut idx = vec![0usize; k];
        subsets(0, p, 0, k, &mut idx, &mut |set: &[usize]| {
            // skip if any member is a zero column (covered by size-1)
            for &i in set {
                if cols[i][..d].iter().all(|&x| x == 0) {
                    return;
                }
            }
            // solve sum alpha_i col_i = 0: k unknowns, d equations. minimal =>
            // rank k-1. alpha_i = (-1)^i det(matrix without row i) using any
            // (k-1)-subset of coordinates with full rank -- take all C(d, k-1)
            // coordinate subsets until a nonzero cofactor vector appears.
            let m = k - 1;
            if m > d {
                return;
            }
            let coords: Vec<usize> = (0..d).collect();
            let mut cidx = vec![0usize; m];
            let mut found: Option<Vec<i128>> = None;
            subsets(0, coords.len(), 0, m, &mut cidx, &mut |cset: &[usize]| {
                if found.is_some() {
                    return;
                }
                let mut alpha = vec![0i128; k];
                let mut nonzero = false;
                for i in 0..k {
                    // det of matrix rows = set minus i, cols = cset
                    let mut mtx = [[0i128; 8]; 8];
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
                    return;
                }
                // verify dependency and positivity (either global sign)
                for c in 0..d {
                    let mut s = 0i128;
                    for i in 0..k {
                        s += alpha[i] * cols[set[i]][c] as i128;
                    }
                    if s != 0 {
                        return;
                    }
                }
                let pos = alpha.iter().all(|&x| x > 0);
                let neg = alpha.iter().all(|&x| x < 0);
                if pos {
                    found = Some(alpha);
                } else if neg {
                    found = Some(alpha.iter().map(|&x| -x).collect());
                }
            });
            if let Some(alpha) = found {
                deps.push((set.to_vec(), alpha));
            }
        });
    }
    deps
}

fn subsets(start: usize, n: usize, pos: usize, k: usize, idx: &mut Vec<usize>, f: &mut dyn FnMut(&[usize])) {
    if pos == k {
        f(&idx[..k]);
        return;
    }
    for i in start..n {
        idx[pos] = i;
        subsets(i + 1, n, pos + 1, k, idx, f);
    }
}

fn det_n(a: &[[i128; 8]; 8], m: usize) -> i128 {
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
                d += s * a[0][c] * det_n(&sub, m - 1);
            }
            d
        }
    }
}

// bound-pruned DFS: d equality rows + dependency strict rows, i128 sums
#[allow(clippy::too_many_arguments)]
fn dfs4(
    k: usize,
    p: usize,
    ne: usize,
    nt: usize,
    s: &mut [i128; 40],
    r: &mut [i32; 10],
    rows: &[[i128; 10]; 40],
    asum: &[[i128; 11]; 40],
    out: &mut Vec<[i32; 10]>,
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
        out.push(*r);
        return;
    }
    for val in [0i32, -1, 1] {
        r[k] = val;
        if val != 0 {
            let f = val as i128;
            for i in 0..nt {
                s[i] += f * rows[i][k];
            }
        }
        dfs4(k + 1, p, ne, nt, s, r, rows, asum, out);
        if val != 0 {
            let f = val as i128;
            for i in 0..nt {
                s[i] -= f * rows[i][k];
            }
        }
    }
    r[k] = 0;
}

fn main() {
    let n: usize = env::args().nth(1).and_then(|s| s.parse().ok()).expect("usage: inc4 n < f3-parents");
    assert!(n % 2 == 0 && (6..=10).contains(&n));
    let p = n - 1;
    let d = 3usize;
    let psq = p * p;
    let pbytes = (psq + 5) / 6;
    let reclen = 2 + pbytes + 1;
    let lcm = lcm_to(2 * n as u64);
    let pfact = factorial(p as u64);

    let mut stdin = io::stdin().lock();
    let mut buf = vec![0u8; 1 << 20];
    let mut have = 0usize;
    let (mut parents, mut leaves) = (0u64, 0u64);
    let mut sum: u128 = 0;
    let mut valid: Vec<[i32; 10]> = Vec::with_capacity(64);

    loop {
        let got = stdin.read(&mut buf[have..]).unwrap();
        if got == 0 {
            break;
        }
        have += got;
        let nrec = have / reclen;
        for ri in 0..nrec {
            let rec = &buf[ri * reclen..(ri + 1) * reclen];
            assert!(rec[0] == b'&' && rec[1] as usize == 63 + p);
            parents += 1;
            let mut pb = [0u16; 16];
            let payload = &rec[2..2 + pbytes];
            let mut kk = 0usize;
            'dec: for &byte in payload {
                let mut bits = ((byte - 63) as u32) << 26;
                for _ in 0..6 {
                    if bits & 0x8000_0000 != 0 {
                        pb[kk / p] |= 1 << (kk % p);
                    }
                    bits <<= 1;
                    kk += 1;
                    if kk == psq {
                        break 'dec;
                    }
                }
            }
            let mut m = [[0i64; 16]; 16];
            for i in 0..p {
                let mut w = pb[i];
                while w != 0 {
                    let j = w.trailing_zeros() as usize;
                    w &= w - 1;
                    m[i][j] = 1;
                    m[j][i] = -1;
                }
            }
            let basis = kernel_basis_exact(&m, p, d).expect("parent not nullity-3");
            // bordered adjugate
            let mut b0 = [[0i128; 12]; 12];
            for i in 0..p {
                for j in 0..p {
                    b0[i][j] = m[i][j] as i128;
                }
                for t in 0..d {
                    b0[i][p + t] = basis[t][i] as i128;
                    b0[p + t][i] = basis[t][i] as i128;
                }
            }
            let (adj, det) = adjugate12(&b0, p + d);
            let sgn: i128 = if det > 0 { 1 } else { -1 };
            // dependency list of the kernel columns
            let cols: Vec<[i64; 4]> = (0..p)
                .map(|j| {
                    let mut c = [0i64; 4];
                    for t in 0..d {
                        c[t] = basis[t][j];
                    }
                    c
                })
                .collect();
            let deps = positive_dependencies(&cols, p, d);
            // DFS rows: d equalities (basis rows), then one strict row per
            // dependency: sum_i alpha_i w_i > 0 <=> sum_i alpha_i sgn (ADJ_i . r) < 0
            let mut rows = [[0i128; 10]; 40];
            for t in 0..d {
                for j in 0..p {
                    rows[t][j] = basis[t][j] as i128;
                }
            }
            let mut nt = d;
            for (set, alpha) in &deps {
                assert!(nt < 40, "too many dependencies");
                for j in 0..p {
                    let mut g = 0i128;
                    for (t, &i) in set.iter().enumerate() {
                        g += alpha[t] * sgn * adj[i][j];
                    }
                    rows[nt][j] = g;
                }
                nt += 1;
            }
            let mut asum = [[0i128; 11]; 40];
            for i in 0..nt {
                for k in (0..p).rev() {
                    asum[i][k] = asum[i][k + 1] + rows[i][k].abs();
                }
            }
            valid.clear();
            let mut s0 = [0i128; 40];
            let mut r0 = [0i32; 10];
            dfs4(0, p, d, nt, &mut s0, &mut r0, &rows, &asum, &mut valid);
            if valid.is_empty() {
                continue;
            }
            let mut arc64 = [0u64; 16];
            for i in 0..p {
                arc64[i] = pb[i] as u64;
            }
            let aut = unsafe { rps_autsize(arc64.as_ptr(), p as c_int) } as u128;
            let wp = pfact / aut;
            for rv in valid.iter() {
                // child bitmasks
                let mut cb = pb;
                for i in 0..p {
                    if rv[i] > 0 {
                        cb[i] |= 1 << p;
                    }
                }
                let mut nr = 0u16;
                for i in 0..p {
                    if rv[i] < 0 {
                        nr |= 1 << i;
                    }
                }
                cb[p] = nr;
                if !paradox_connected(&cb, n) {
                    continue;
                }
                // w numerators: wnum[i] = -(sgn * ADJ_i . r); w_i = wnum[i]/|det|
                let mut wnum = [0i128; 10];
                for i in 0..p {
                    let mut a = 0i128;
                    for j in 0..p {
                        a += adj[i][j] * rv[j] as i128;
                    }
                    wnum[i] = -sgn * a;
                }
                // z4: count vertices v with a nonzero nonneg vector in the slice
                let mut z4 = 0u64;
                // v = new vertex: slice = {t=0} -> parent's nonneg kernel: always valid
                z4 += 1;
                for v in 0..p {
                    // equality e.(lam, t) = 0 with e = (B1v, B2v, B3v, w_v);
                    // scale w_v to integers: use (|det| * Bv..., wnum[v])
                    let e = [
                        basis[0][v] as i128 * det.abs(),
                        basis[1][v] as i128 * det.abs(),
                        basis[2][v] as i128 * det.abs(),
                        wnum[v],
                    ];
                    if slice_has_nonneg(&e, &basis, &wnum, det.abs(), p) {
                        z4 += 1;
                    }
                }
                leaves += 1;
                sum += wp * (lcm / z4 as u128);
            }
        }
        let rem = have - nrec * reclen;
        buf.copy_within(nrec * reclen..have, 0);
        have = rem;
    }
    let total = (n as u128) * sum;
    assert!(total % lcm == 0, "1/z4 weights not integral");
    println!(
        "n={}: parents={} leaves={} L_nullity4_labeled={}",
        n,
        parents,
        leaves,
        total / lcm
    );
}

// does the 4-var cone { (lam,t): B^T lam + t w >= 0 (p rows), t >= 0 }
// intersected with { e.(lam,t) = 0 } contain a nonzero point? Exact: project
// onto an integer basis of e's nullspace (3-dim), then extreme-ray enumeration
// via common::cone_has_nonneg on the transformed constraint columns.
fn slice_has_nonneg(e: &[i128; 4], basis: &[[i64; 16]], wnum: &[i128; 10], absd: i128, p: usize) -> bool {
    // integer nullspace basis of e (3 vectors in 4-space)
    let mut nb: Vec<[i128; 4]> = Vec::with_capacity(3);
    if e.iter().all(|&x| x == 0) {
        // slice = whole cone; parent membership in F3' does NOT directly apply
        // (that was the t=0 slice); the full cone contains (0,0,0,1)->w... the
        // cone contains nonzero points iff ... just use the 4 standard axes
        nb.push([1, 0, 0, 0]);
        nb.push([0, 1, 0, 0]);
        nb.push([0, 0, 1, 0]);
        // 4th direction dropped: use 3-dim subcone test on lam-space + t=free?
        // simplest sound choice: the vector (lam,t)=(0,0,0,1) gives y = w which
        // needs w >= 0 -- check directly; else fall through to lam-space test
        // (t=0): B^T lam >= 0 nonzero, i.e. the parent cone: true by family.
        return true;
    }
    let j = (0..4).find(|&i| e[i] != 0).unwrap();
    for k in 0..4 {
        if k == j {
            continue;
        }
        let mut v = [0i128; 4];
        v[k] = e[j];
        v[j] = -e[k];
        let g = gcd128(v.iter().map(|x| x.abs()).fold(0, gcd128), 0).max(1);
        for x in v.iter_mut() {
            *x /= g;
        }
        nb.push(v);
    }
    // transformed constraint columns: for each constraint row c (in 4-space),
    // column_t = c . nb[t]. Constraints: p rows (B1i, B2i, B3i, w_i-scaled) and
    // (0,0,0,1) for t >= 0. Scale consistently: row_i = (|D| B1i, |D| B2i,
    // |D| B3i, wnum_i); t-row = (0,0,0,1).
    let mut cols: Vec<[i64; 16]> = Vec::with_capacity(p + 1);
    for i in 0..p {
        let row = [
            basis[0][i] as i128 * absd,
            basis[1][i] as i128 * absd,
            basis[2][i] as i128 * absd,
            wnum[i],
        ];
        let mut col = [0i64; 16];
        for (t, b) in nb.iter().enumerate() {
            let mut dot = 0i128;
            for c in 0..4 {
                dot += row[c] * b[c];
            }
            // reduce magnitude
            col[t] = clampdown(dot);
        }
        cols.push(col);
    }
    {
        let trow = [0i128, 0, 0, 1];
        let mut col = [0i64; 16];
        for (t, b) in nb.iter().enumerate() {
            let mut dot = 0i128;
            for c in 0..4 {
                dot += trow[c] * b[c];
            }
            col[t] = clampdown(dot);
        }
        cols.push(col);
    }
    // reuse cone_has_nonneg by presenting the transposed data as a 3 x (p+1)
    // "basis": cone_has_nonneg(basis=d rows over n cols) tests { y = B^T lam }
    // -- exactly our transformed cone with n = p+1 constraints
    let mut tb: Vec<[i64; 16]> = vec![[0i64; 16]; 3];
    for (ci, col) in cols.iter().enumerate() {
        for t in 0..3 {
            tb[t][ci] = col[t];
        }
    }
    cone_has_nonneg(&tb, p + 1, 3)
}

fn clampdown(x: i128) -> i64 {
    assert!(x.abs() < (1i128 << 62), "slice projection overflow");
    x as i64
}

fn paradox_connected(beats: &[u16; 16], n: usize) -> bool {
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
    for i in 0..n {
        if beats[i] == 0 || inn[i] == 0 {
            return false;
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
