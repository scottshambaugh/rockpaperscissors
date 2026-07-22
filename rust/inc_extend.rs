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
// EXACT integer throughout (no f64): the CM (d=0) branch mirrors cm_extend
// (adjugate_ff_i64 + common::dfs_neg; -M'^-1 r > 0 <=> sgn(det)*(adj r)_i < 0),
// and the singular (d>=2) branch uses common::kernel_basis_any + the exact
// perpendicularity DFS common::dfs_perp_eq (r _|_ ker), keeping each child that
// is paradox+connected, has a strictly positive kernel (has_positive_kernel),
// and is canonical (nullity3_canonical via exact_rank). Anchors: cm(5)=7,
// cm(7)=7268, hi(7)=3257, inclusive(7)=10525.

use std::env;
use std::io::{self, BufWriter, Read, Write};

mod common;
use common::{added_is_maximal, encode, has_positive_kernel, sig_maximal};

// Bound-pruned DFS collecting r in {-1,0,1}^p with N r = 0 (r _|_ K), N = kernel
// basis (d rows, columns reordered so the DFS assigns high-L1 columns first).
// `s[i]` is the running i-th constraint N_i . r; `asum[i][k] = sum_{j>=k}|Nc[i][j]|`
// bounds how far it can still move, so if |s[i]| > asum[i][k] the constraint can no
// longer reach 0 -- prune. Turns the 3^p flat r-scan into a walk of a small live
// subtree, exactly analogous to the completely-mixed cone-DFS.
#[allow(clippy::too_many_arguments)]

// ---- kernel basis via RREF: rows of `basis[..d]` span ker(a) (p x p) ----

// rank of the skew matrix `m` (n x n) with row/col `skip` deleted (size n-1),
// via float Gaussian elimination counting pivots.

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
    let m = common::skew_from_beats(beats, n);
    // p competes only against the other nullity-d deletions (see doc above);
    // p itself qualifies by construction. Eligibility of v = "nullity of M with
    // vertex v deleted equals d" -- exact integer rank (no f64 tolerance).
    sig_maximal(beats, n, p, |v| {
        let mut mv = [[0i64; 16]; 16];
        let mut ii = 0usize;
        for i in 0..n {
            if i == v {
                continue;
            }
            let mut jj = 0usize;
            for j in 0..n {
                if j == v {
                    continue;
                }
                mv[ii][jj] = m[i][j];
                jj += 1;
            }
            ii += 1;
        }
        (n - 1) - common::exact_rank(&mv, n - 1) == d
    })
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

// A dependency (i, j, k, ai, aj, ak): a positive combination of <= 3 of the
// fixed 2-D normals that sums to zero. j == usize::MAX marks a zero normal
// (constraint reduces to w_i > 0); k == usize::MAX marks an antiparallel pair.
// The candidate is feasible iff EVERY dependency has sum alpha_j w_j > 0.


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
#[allow(clippy::too_many_arguments)]

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
    let pbytes = (p * p + 5) / 6;
    let preclen = 2 + pbytes + 1;

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

    let mut r = [0i32; 16];
    let mut valid: Vec<[i32; 16]> = Vec::with_capacity(256);
    let (mut parents, mut emit_cm, mut emit_hi) = (0u64, 0u64, 0u64);
    let mut sing = 0u64;

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
            // unpack the p-vertex digraph6 record into out-adjacency masks
            let mut base = [0u16; 16];
            let payload = &rec[2..2 + pbytes];
            let mut ukk = 0usize;
            'dec: for &byte in payload {
                let mut bits = ((byte - 63) as u32) << 26;
                for _ in 0..6 {
                    if bits & 0x8000_0000 != 0 {
                        base[ukk / p] |= 1 << (ukk % p);
                    }
                    bits <<= 1;
                    ukk += 1;
                    if ukk == p * p {
                        break 'dec;
                    }
                }
            }
            // integer skew M' (i64 for the kernel, i128 for the adjugate)
            let mi = common::skew_from_beats(&base, p);
            let mut mp = [[0i128; 16]; 16];
            for i in 0..p {
                for j in 0..p {
                    mp[i][j] = mi[i][j] as i128;
                }
            }
            let adjopt = match common::adjugate_ff_i64(&mp, p, false) {
                Ok(x) => x,
                Err(()) => common::adjugate_ff(&mp, p, false),
            };
            match adjopt {
                Some((adj, det)) => {
                    // nonsingular: completely-mixed (nullity-1) children, exact.
                    // -M'^-1 r > 0  <=>  sgn(det)*(adj . r)_i < 0 for all i.
                    let sd: i128 = if det > 0 { 1 } else { -1 };
                    let mut ord = [0usize; 16];
                    for (k, o) in ord[..p].iter_mut().enumerate() {
                        *o = k;
                    }
                    let mut l1 = [0i128; 16];
                    for k in 0..p {
                        let mut sum = 0i128;
                        for i in 0..p {
                            sum += adj[i][k].abs();
                        }
                        l1[k] = sum;
                    }
                    ord[..p].sort_unstable_by(|&a, &b| l1[b].cmp(&l1[a]));
                    let mut wc = [[0i64; 16]; 16];
                    for k in 0..p {
                        for i in 0..p {
                            wc[i][k] = common::narrow_i64(sd * adj[i][ord[k]]);
                        }
                    }
                    let mut asum = [[0i64; 16]; 17];
                    for k in (0..p).rev() {
                        for i in 0..p {
                            asum[k][i] = asum[k + 1][i] + wc[i][k].abs();
                        }
                    }
                    valid.clear();
                    let mut s0 = [0i64; 16];
                    common::dfs_neg(0, p, &mut s0, &mut r, &wc, &asum, &mut valid);
                    if valid.is_empty() {
                        continue;
                    }
                    // parent degrees, hoisted out of the candidate loop: the child's
                    // degrees are od[c] + [c beats new], id[c] + [new beats c], so
                    // the (od, id) first stage of the maximality signature costs O(p)
                    // per candidate instead of the O(n^2) recomputation inside
                    // added_is_maximal -- which then only runs on exact (od, id) ties.
                    let mut pod = [0u8; 16];
                    let mut pid = [0u8; 16];
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
                    // singular: nullity-(d+1) children via r _|_ ker(M'), kept
                    // fully mixed. Exact: integer kernel basis, exact
                    // perpendicularity DFS, exact fully-mixed + canonical tests.
                    // (The float d==2 bordered/Helly fast path is dropped; the
                    // exact perpendicularity DFS + per-candidate cone test gives
                    // the identical children -- see README.)
                    let bvecs = common::kernel_basis_any(&mi, p);
                    let d = bvecs.len();
                    if d == 0 {
                        continue;
                    }
                    sing += 1;
                    let mut ord = [0usize; 16];
                    for (k, o) in ord[..p].iter_mut().enumerate() {
                        *o = k;
                    }
                    let mut l1 = [0i64; 16];
                    for k in 0..p {
                        let mut sum = 0i64;
                        for bv in bvecs.iter() {
                            sum += bv[k].abs();
                        }
                        l1[k] = sum;
                    }
                    ord[..p].sort_unstable_by(|&a, &b| l1[b].cmp(&l1[a]));
                    let mut rows = [[0i64; 16]; 16];
                    for (i, bv) in bvecs.iter().enumerate() {
                        for k in 0..p {
                            rows[i][k] = bv[ord[k]];
                        }
                    }
                    let mut asum = [[0i64; 16]; 17];
                    for k in (0..p).rev() {
                        for i in 0..d {
                            asum[i][k] = asum[i][k + 1] + rows[i][k].abs();
                        }
                    }
                    valid.clear();
                    let mut s0 = [0i64; 16];
                    common::dfs_perp_eq(0, p, d, &mut s0, &mut r, &rows, &asum, &mut valid);
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
                            && has_positive_kernel(&common::skew_from_beats(&beats[..n], n), n)
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
        "incx n={}: parents={} singular={} emitted cm(nullity1)={} nullity>=3={} (pipe: labelg|sort -u|wc -l)",
        n, parents, sing, emit_cm, emit_hi
    );
}
