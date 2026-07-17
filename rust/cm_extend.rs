// Enumerate completely mixed games on n (odd) vertices by EXTENDING the
// (n-1)-vertex oriented games, instead of filtering all A001174(n) oriented
// n-graphs. Reads digraph6 (n-1)-vertex oriented graphs (from `nauty-directg
// -o`) on stdin and writes digraph6 of every completely mixed n-vertex child;
// pipe through `nauty-labelg | sort -u | wc -l` to count iso classes.
//
//   rustc -O rust/cm_extend.rs -o /tmp/cmx
//   nauty-geng 8 | nauty-directg -o | /tmp/cmx 9 | nauty-labelg | sort -u | wc -l
//
// The lemma (see rust/README.md): for a game M on n vertices, delete vertex n-1
// to get the (n-1)x(n-1) skew M'. If M' is nonsingular then the extended M has
// nullity EXACTLY 1 for every new-vertex vector r in {-1,0,1}^{n-1}, with kernel
// v = (-M'^{-1} r, 1) -- the row-(n-1) equation r^T v' = 0 holds automatically
// because r^T M'^{-1} r = 0 for skew M'. Hence
//
//     M is completely mixed  <=>  -M'^{-1} r > 0 componentwise,
//
// and that single strict-positivity condition already forces paradoxical +
// connected (a zero/one-signed-failing coordinate is exactly a dropped strategy
// or a disconnecting even block). So we never build a 9-graph that isn't CM:
// each nonsingular (n-1)-parent contributes only its finitely many CM children.
//
// Since M'^{-1} = adj(M')/det with det = Pf(M')^2 > 0, the sign of v'_j is the
// sign of -(M'^{-1} r)_j; we compute M'^{-1} once in f64 (exact enough: a nonzero
// v'_j has magnitude >= 1/det >= 1e-7, far above the 1e-13 round-off) and reject
// r on the first non-positive coordinate. Validated to reproduce the exact
// Pfaffian census: n=5 -> 7, n=7 -> 7268.

use std::env;
use std::io::{self, Read, Write, BufWriter};

const TOL: f64 = 1e-7;

// Gauss-Jordan inverse of a p x p matrix (p <= 8) with partial pivoting.
// Returns None if singular. Also usable to detect the singular parents.
fn inverse(a: &[[f64; 8]; 8], p: usize) -> Option<[[f64; 8]; 8]> {
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
fn dfs(
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

// digraph6 decode: '&' + (63+p) + row-major arc bits. adj[i] bit j <=> i beats j.
fn decode(rec: &[u8], p: usize, out: &mut [[f64; 8]; 8]) {
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

// Canonical-deletion prefilter. Each n-vertex CM game M arises once per
// deletion-parent M-v (all n deletions of a CM game are nonsingular, hence all
// are parents) -- a ~n-fold redundancy the downstream `labelg | sort -u` must
// otherwise absorb. We cut it at the source: emit a child only when its added
// vertex (index p) is a MAXIMAL vertex under a cheap isomorphism-invariant
// signature (one round of degree refinement). Since some vertex is maximal, at
// least one of M's parents reconstructs it with the added vertex maximal, so no
// class is ever lost (sound); the only survivors past the filter are the
// ~n-times-fewer maximal copies, and the rare ties (symmetric games with several
// maximal vertices) are mopped up by the final `sort -u`. Rigid games -- almost
// all at n=9 -- have a unique maximum, so the filter is essentially exact.
//
// Signature of vertex v: (outdeg, indeg, sorted out-neighbour (outdeg,indeg),
// sorted in-neighbour (outdeg,indeg)), compared lexicographically. Returns true
// iff vertex p is >= every other vertex.
fn added_is_maximal(beats: &[u16], n: usize, p: usize) -> bool {
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
        if v == p {
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

// Encode an n-vertex oriented game (beats[i] bitmask: i beats j) to digraph6.
fn encode(beats: &[u16], n: usize, buf: &mut Vec<u8>) {
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

fn main() {
    let n: usize = env::args().nth(1).and_then(|s| s.parse().ok()).expect("usage: cmx n");
    assert!(n >= 3 && n % 2 == 1 && n <= 16, "n must be odd, 3..=15");
    let p = n - 1; // parent (even) size
    let preclen = 2 + (p * p + 5) / 6 + 1; // digraph6 record length for p-vertex

    let mut stdin = io::stdin().lock();
    let mut buf = vec![0u8; 1 << 20];
    let mut have = 0usize;
    let out = io::stdout();
    let mut w = BufWriter::with_capacity(1 << 20, out.lock());
    let mut child = Vec::with_capacity(2 + (n * n + 5) / 6 + 1);

    let mut mp = [[0f64; 8]; 8];
    let mut r = [0i32; 8];
    let mut valid: Vec<[i32; 8]> = Vec::with_capacity(256);
    let (mut parents, mut nonsing, mut emitted) = (0u64, 0u64, 0u64);

    loop {
        let got = stdin.read(&mut buf[have..]).unwrap();
        if got == 0 {
            break;
        }
        have += got;
        let nrec = have / preclen;
        for ri in 0..nrec {
            let rec = &buf[ri * preclen..(ri + 1) * preclen];
            if rec[0] != b'&' {
                continue;
            }
            parents += 1;
            decode(rec, p, &mut mp);
            let inv = match inverse(&mp, p) {
                Some(x) => x,
                None => continue, // singular parent: no completely mixed child
            };
            nonsing += 1;
            // order columns by descending L1 norm so the DFS bound tightens
            // fastest (more pruning), then build the reordered column matrix wc
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
            // suffix absolute column-sums of wc for the DFS bound
            let mut asum = [[0f64; 8]; 9];
            for k in (0..p).rev() {
                for i in 0..p {
                    asum[k][i] = asum[k + 1][i] + wc[i][k].abs();
                }
            }
            // collect the (few) valid r (in ordered coords) via bound-pruned DFS
            valid.clear();
            let mut s0 = [0f64; 8];
            dfs(0, p, &mut s0, &mut r, &wc, &asum, &mut valid);
            if valid.is_empty() {
                continue;
            }
            // parent arcs are the same for every child; build them once
            let mut base = [0u16; 16];
            for i in 0..p {
                for j in 0..p {
                    if mp[i][j] > 0.5 {
                        base[i] |= 1 << j;
                    }
                }
            }
            for rv in valid.iter() {
                let mut beats = base;
                for k in 0..p {
                    // rv is in reordered coords: position k is original column ord[k].
                    // coeff = M[c][p]: +1 => c beats new vertex p; -1 => p beats c
                    let c = ord[k];
                    if rv[k] > 0 {
                        beats[c] |= 1 << p;
                    } else if rv[k] < 0 {
                        beats[p] |= 1 << c;
                    }
                }
                // canonical-deletion prefilter: keep only maximal-added-vertex copies
                if !added_is_maximal(&beats, n, p) {
                    continue;
                }
                encode(&beats[..n], n, &mut child);
                w.write_all(&child).unwrap();
                emitted += 1;
            }
        }
        let rem = have - nrec * preclen;
        buf.copy_within(nrec * preclen..have, 0);
        have = rem;
    }
    w.flush().unwrap();
    eprintln!(
        "cmx n={}: parents={} nonsingular={} cm_children_emitted={} (pipe through labelg|sort -u|wc -l)",
        n, parents, nonsing, emitted
    );
}
