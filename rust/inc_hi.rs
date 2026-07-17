// Nullity-4 (and, with arg D=6, nullity-6) stratum of inclusive(even n):
// labeled count via the same deletion-multiplicity identity as inc10.rs, one
// (resp. two) nullity levels up. Everything below is written for kernel
// dimension d = D-1 as a runtime value; only array capacities are fixed.
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
//   rustc -O rust/inc_hi.rs -o /tmp/inc_hi -C link-args="shim.o -lnauty"
//   nauty-geng 7 | nauty-directg -o | /tmp/incs 7 f3-emit 3 | /tmp/inc_hi 8
// Nullity-6 anchor: parents f3-emit 5, L_6 = 210,882:
//   nauty-geng 7 | nauty-directg -o | /tmp/incs 7 f3-emit 5 | /tmp/inc_hi 8 6
use std::env;
use std::io::{self, Read};
use std::os::raw::c_int;

mod common;
use common::{adjugate_ff, cone_has_nonneg, dfs_es, factorial, gcd_i128, kernel_basis_exact, lcm_to, order_columns_l1, paradox_connected_beats, positive_dependencies, suffix_abs_sums, DFS_RCAP};

extern "C" {
    fn rps_autsize(arc: *const u64, n: c_int) -> f64;
}

// bordered matrices here are always nonsingular; verify stays on (per-parent)
fn adjugate16(b0: &[[i128; 16]; 16], m: usize) -> ([[i128; 16]; 16], i128) {
    adjugate_ff(b0, m, true).expect("bordered singular")
}

fn main() {
    let n: usize = env::args().nth(1).and_then(|s| s.parse().ok()).expect("usage: inc_hi n [D=4|6] < f-parents");
    let dd: usize = env::args().nth(2).and_then(|s| s.parse().ok()).unwrap_or(4);
    assert!(n % 2 == 0 && (6..=10).contains(&n));
    assert!(dd == 4 || dd == 6, "stratum D must be 4 or 6");
    let p = n - 1;
    let d = dd - 1;
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
            // parent degrees, for the child degree stage at leaves
            let mut pod = [0u8; 16];
            let mut pid = [0u8; 16];
            for i in 0..p {
                pod[i] = pb[i].count_ones() as u8;
                let mut c = 0u8;
                for u in 0..p {
                    if pb[u] & (1 << i) != 0 {
                        c += 1;
                    }
                }
                pid[i] = c;
            }
            // bordered adjugate
            let mut b0 = [[0i128; 16]; 16];
            for i in 0..p {
                for j in 0..p {
                    b0[i][j] = m[i][j] as i128;
                }
                for t in 0..d {
                    b0[i][p + t] = basis[t][i] as i128;
                    b0[p + t][i] = basis[t][i] as i128;
                }
            }
            let (adj, det) = adjugate16(&b0, p + d);
            let sgn: i128 = if det > 0 { 1 } else { -1 };
            // dependency list of the kernel columns
            let cols: Vec<[i64; 10]> = (0..p)
                .map(|j| {
                    let mut c = [0i64; 10];
                    for t in 0..d {
                        c[t] = basis[t][j];
                    }
                    c
                })
                .collect();
            let deps = positive_dependencies(&cols, p, d);
            // DFS rows: d equalities (basis rows), then one strict row per
            // dependency: sum_i alpha_i w_i > 0 <=> sum_i alpha_i sgn (ADJ_i . r) < 0.
            // i64 is ample: measured max |entry| at n=8 is 680,400; the guard
            // aborts loudly if n=10 bordered adjugates ever outgrow it.
            let mut rows = [[0i64; 16]; DFS_RCAP];
            for t in 0..d {
                for j in 0..p {
                    rows[t][j] = basis[t][j];
                }
            }
            let mut nt = d;
            for (set, alpha) in &deps {
                assert!(nt < DFS_RCAP, "too many dependencies");
                for j in 0..p {
                    let mut g = 0i128;
                    for (t, &i) in set.iter().enumerate() {
                        g += alpha[t] * sgn * adj[i][j];
                    }
                    assert!(g.abs() < (1i128 << 55), "dependency row overflow");
                    rows[nt][j] = g as i64;
                }
                nt += 1;
            }
            // paradox-forced coordinates: a parent vertex with no win must
            // beat the new vertex (r=+1); no loss => lose to it (r=-1)
            let mut pinn = [0u16; 16];
            for i in 0..p {
                let mut w = pb[i];
                while w != 0 {
                    let j = w.trailing_zeros() as usize;
                    w &= w - 1;
                    pinn[j] |= 1 << i;
                }
            }
            let mut ord = [0usize; 16];
            let mut prows = [[0i64; 16]; DFS_RCAP];
            order_columns_l1(&rows, nt, p, &mut ord, &mut prows);
            let mut asum = [[0i64; 17]; DFS_RCAP];
            suffix_abs_sums(&prows, nt, p, &mut asum);
            let mut fplus = 0u16;
            let mut fminus = 0u16;
            let mut impossible = false;
            for k in 0..p {
                let c = ord[k];
                let nw = pb[c] == 0;
                let nl = pinn[c] == 0;
                if nw && nl {
                    impossible = true;
                    break;
                }
                if nw {
                    fplus |= 1 << k;
                }
                if nl {
                    fminus |= 1 << k;
                }
            }
            if impossible {
                continue;
            }
            valid.clear();
            let mut s0 = [0i64; DFS_RCAP];
            let mut r0 = [0i32; 16];
            dfs_es(0, p, d, nt, &mut s0, &mut r0, &prows, &asum, fplus, fminus, &mut |rv: &[i32; 16]| {
                let mut orig = [0i32; 10];
                for k in 0..p {
                    orig[ord[k]] = rv[k];
                }
                valid.push(orig);
            });
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
                if !paradox_connected_beats(&cb, n) {
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
                // two-sided rule (the stratum-2 trick generalized): instead
                // of 1/z_D over all valid deletions, accept iff the new vertex
                // is the sig-argmax of the child's valid-deletion set V(C),
                // fractional 1/T on sig ties. The new vertex is ALWAYS in V(C)
                // (the parent is F_d' by construction), so slice-cone tests
                // run only for vertices that degree-dominate or degree-tie the
                // new vertex -- typically 0-2 per leaf instead of p.
                let slice_ok = |v: usize| -> bool {
                    let mut e = [0i128; 6];
                    for t in 0..d {
                        e[t] = basis[t][v] as i128 * det.abs();
                    }
                    e[d] = wnum[v];
                    slice_has_nonneg(&e, &basis, &wnum, det.abs(), p, d)
                };
                let nod = nr.count_ones() as u8;
                let mut nid = 0u8;
                for i in 0..p {
                    if rv[i] > 0 {
                        nid += 1;
                    }
                }
                let mut beaten = false;
                let mut degree_ties = 0u16;
                for v in 0..p {
                    let cod = pod[v] + (rv[v] > 0) as u8;
                    let cid = pid[v] + (rv[v] < 0) as u8;
                    if cod > nod || (cod == nod && cid > nid) {
                        if slice_ok(v) {
                            beaten = true;
                            break;
                        }
                    } else if cod == nod && cid == nid {
                        degree_ties |= 1 << v;
                    }
                }
                if beaten {
                    continue;
                }
                let mut t1 = 1u64;
                if degree_ties != 0 {
                    // full signature comparison for degree-tied competitors
                    let mut csig = [0u64; 16];
                    common::vertex_sigs(&cb, n, &mut csig);
                    let ns = csig[p];
                    let mut reject = false;
                    let mut mm = degree_ties;
                    while mm != 0 {
                        let v = mm.trailing_zeros() as usize;
                        mm &= mm - 1;
                        if csig[v] > ns && slice_ok(v) {
                            reject = true;
                            break;
                        }
                        if csig[v] == ns && slice_ok(v) {
                            t1 += 1;
                        }
                    }
                    if reject {
                        continue;
                    }
                }
                leaves += 1;
                sum += wp * (lcm / t1 as u128);
            }
        }
        let rem = have - nrec * reclen;
        buf.copy_within(nrec * reclen..have, 0);
        have = rem;
    }
    let total = (n as u128) * sum;
    assert!(total % lcm == 0, "1/z weights not integral");
    println!(
        "n={}: parents={} leaves={} L_nullity{}_labeled={}",
        n,
        parents,
        leaves,
        dd,
        total / lcm
    );
}

// does the (d+1)-var cone { (lam,t): B^T lam + t w >= 0 (p rows), t >= 0 }
// intersected with { e.(lam,t) = 0 } contain a nonzero point? Exact: project
// onto an integer basis of e's nullspace (d-dim), then extreme-ray enumeration
// via common::cone_has_nonneg on the transformed constraint columns.
fn slice_has_nonneg(e: &[i128; 6], basis: &[[i64; 16]], wnum: &[i128; 10], absd: i128, p: usize, d: usize) -> bool {
    let dv = d + 1; // ambient (lam, t) dimension
    let mut nb: Vec<[i128; 6]> = Vec::with_capacity(d);
    if e[..dv].iter().all(|&x| x == 0) {
        // slice = whole cone; parent membership in F3' does NOT directly apply
        // (that was the t=0 slice); the full cone contains (0,0,0,1)->w... the
        // cone contains nonzero points iff ... just use the 4 standard axes
        // slice = whole cone; the t=0 sub-slice is the parent cone which
        // contains a nonzero nonneg vector by family membership
        return true;
    }
    let j = (0..dv).find(|&i| e[i] != 0).unwrap();
    for k in 0..dv {
        if k == j {
            continue;
        }
        let mut v = [0i128; 6];
        v[k] = e[j];
        v[j] = -e[k];
        let g = gcd_i128(v.iter().map(|x| x.abs()).fold(0, gcd_i128), 0).max(1);
        for x in v.iter_mut() {
            *x /= g;
        }
        nb.push(v);
    }
    // transformed constraint columns: for each constraint row c (in (d+1)-
    // space), column_t = c . nb[t]. Constraints: p rows (B_1i..B_di, w_i-
    // scaled) and the t >= 0 row. Scale consistently: row_i =
    // (|D| B_1i, .., |D| B_di, wnum_i); t-row = (0,..,0,1).
    let mut cols: Vec<[i64; 16]> = Vec::with_capacity(p + 1);
    for i in 0..p {
        let mut row = [0i128; 6];
        for t in 0..d {
            row[t] = basis[t][i] as i128 * absd;
        }
        row[d] = wnum[i];
        let mut col = [0i64; 16];
        for (t, b) in nb.iter().enumerate() {
            let mut dot = 0i128;
            for c in 0..dv {
                dot += row[c] * b[c];
            }
            // reduce magnitude
            col[t] = clampdown(dot);
        }
        cols.push(col);
    }
    {
        let mut trow = [0i128; 6];
        trow[d] = 1;
        let mut col = [0i64; 16];
        for (t, b) in nb.iter().enumerate() {
            let mut dot = 0i128;
            for c in 0..dv {
                dot += trow[c] * b[c];
            }
            col[t] = clampdown(dot);
        }
        cols.push(col);
    }
    // reuse cone_has_nonneg by presenting the transposed data as a d x (p+1)
    // "basis": cone_has_nonneg(basis=d rows over n cols) tests { y = B^T lam }
    // -- exactly our transformed cone with n = p+1 constraints
    let mut tb: Vec<[i64; 16]> = vec![[0i64; 16]; d];
    for (ci, col) in cols.iter().enumerate() {
        for t in 0..d {
            tb[t][ci] = col[t];
        }
    }
    cone_has_nonneg(&tb, p + 1, d)
}

fn clampdown(x: i128) -> i64 {
    assert!(x.abs() < (1i128 << 62), "slice projection overflow");
    x as i64
}

