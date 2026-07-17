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

mod common;
use common::{added_is_maximal, decode, dfs, encode, inverse};

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
