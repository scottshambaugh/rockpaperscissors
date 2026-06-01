// Isomorph-free count of tournaments and of S2 ("two-paradox" / Erdos-Schutte)
// tournaments on n vertices, as an independent verifier / extender of the Python
// `search_two_paradox`.
//
//   rustc -O rust/s2_count.rs -o /tmp/s2 && /tmp/s2 10
//
// Method: canonical augmentation. Build one representative per isomorphism class
// of k-vertex tournaments, then form every (k+1)-extension of every rep and
// deduplicate by a nauty-style canonical code (color refinement + individualize-
// refine, min code over leaves). The *unfiltered* class count is printed at every
// level and must reproduce A000568 (1,1,2,4,12,56,456,6880,191536,...) -- that is
// the built-in correctness check. S2 = every pair of vertices has a common
// dominator (a third vertex beating both).
//
// Tournament: adj[i] is a bitmask, bit j set => i beats j. Skew by construction.

use std::collections::HashSet;
use std::env;

const MAXN: usize = 12; // vertices 0..11 -> u16 rows, C(12,2)=66 bits -> u128 code

#[derive(Clone, Copy)]
struct Tour {
    n: usize,
    adj: [u16; MAXN],
}

impl Tour {
    #[inline]
    fn beats(&self, i: usize, j: usize) -> bool {
        self.adj[i] & (1 << j) != 0
    }
}

// Coarsest equitable coloring (1-WL): refine by (color, sorted out-neighbour
// colors, sorted in-neighbour colors) until stable. Re-ranked to 0..c each round,
// so a finer input partition only splits further (never merges).
fn refine(t: &Tour, init: &[u32]) -> Vec<u32> {
    let n = t.n;
    let mut col = init.to_vec();
    loop {
        let mut sig: Vec<(u32, Vec<u32>, Vec<u32>)> = Vec::with_capacity(n);
        for i in 0..n {
            let (mut out_c, mut in_c) = (Vec::new(), Vec::new());
            for j in 0..n {
                if i == j {
                    continue;
                }
                if t.beats(i, j) {
                    out_c.push(col[j]);
                } else {
                    in_c.push(col[j]);
                }
            }
            out_c.sort_unstable();
            in_c.sort_unstable();
            sig.push((col[i], out_c, in_c));
        }
        let mut order: Vec<usize> = (0..n).collect();
        order.sort_by(|&a, &b| sig[a].cmp(&sig[b]));
        let mut newcol = vec![0u32; n];
        let mut c = 0u32;
        for idx in 0..n {
            if idx > 0 && sig[order[idx]] != sig[order[idx - 1]] {
                c += 1;
            }
            newcol[order[idx]] = c;
        }
        if newcol == col {
            return col;
        }
        col = newcol;
    }
}

fn code_of_discrete(t: &Tour, col: &[u32]) -> u128 {
    let n = t.n;
    let mut vert_at = [0usize; MAXN];
    for v in 0..n {
        vert_at[col[v] as usize] = v; // discrete: color == final position
    }
    let mut code: u128 = 0;
    let mut bit = 0;
    for a in 0..n {
        for b in (a + 1)..n {
            if t.beats(vert_at[a], vert_at[b]) {
                code |= 1u128 << bit;
            }
            bit += 1;
        }
    }
    code
}

fn canon_rec(t: &Tour, col: Vec<u32>, best: &mut u128) {
    let n = t.n;
    let maxc = *col.iter().max().unwrap();
    if maxc as usize == n - 1 {
        // discrete: every color is a singleton, the ranking is a full labeling
        let code = code_of_discrete(t, &col);
        if code < *best {
            *best = code;
        }
        return;
    }
    let mut freq = vec![0u32; (maxc + 1) as usize];
    for &c in &col {
        freq[c as usize] += 1;
    }
    let target = (0..=maxc).find(|&c| freq[c as usize] > 1).unwrap();
    // branch: individualize each vertex of the first non-singleton cell
    for v in 0..n {
        if col[v] != target {
            continue;
        }
        let mut nc = col.clone();
        for u in 0..n {
            if nc[u] > target || (nc[u] == target && u != v) {
                nc[u] += 1; // push the rest of the cell (and everything above) up by 1
            }
        }
        let refined = refine(t, &nc);
        canon_rec(t, refined, best);
    }
}

fn canon(t: &Tour) -> u128 {
    let n = t.n;
    let mask = (1u16 << n) - 1;
    let init: Vec<u32> = (0..n).map(|i| (t.adj[i] & mask).count_ones()).collect();
    let col = refine(t, &init);
    let mut best = u128::MAX;
    canon_rec(t, col, &mut best);
    best
}

fn extend(parent: &Tour, ext: u16) -> Tour {
    let n = parent.n;
    let mut t = *parent;
    t.n = n + 1;
    let newv = n;
    t.adj[newv] = 0;
    for i in 0..n {
        if ext & (1 << i) != 0 {
            t.adj[newv] |= 1 << i; // new vertex beats i
        } else {
            t.adj[i] |= 1 << newv; // i beats new vertex
        }
    }
    t
}

// S2 / two-paradox: every pair {i,j} has a common dominator k (k beats both).
fn is_s2(t: &Tour) -> bool {
    let n = t.n;
    let mut beaters = [0u16; MAXN]; // beaters[i] = {k : k beats i}
    for k in 0..n {
        let mut r = t.adj[k];
        while r != 0 {
            let i = r.trailing_zeros() as usize;
            beaters[i] |= 1 << k;
            r &= r - 1;
        }
    }
    for i in 0..n {
        for j in (i + 1)..n {
            if beaters[i] & beaters[j] == 0 {
                return false;
            }
        }
    }
    true
}

fn main() {
    let n: usize = env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(9);
    assert!(n >= 1 && n <= MAXN, "n must be 1..={}", MAXN);

    // A000568 for the built-in checksum.
    let a000568 = [1u64, 1, 2, 4, 12, 56, 456, 6880, 191536, 9733056, 903753248];

    let mut reps: Vec<Tour> = vec![Tour { n: 1, adj: [0; MAXN] }];
    println!(" n |   tournaments (A000568 check) | S2 (two-paradox)");
    println!("---+-------------------------------+-----------------");
    for k in 2..=n {
        let last = k == n;
        let mut seen: HashSet<u128> = HashSet::new();
        let mut next: Vec<Tour> = Vec::new();
        let mut s2: u64 = 0;
        for p in &reps {
            for ext in 0..(1u32 << (k - 1)) {
                let t = extend(p, ext as u16);
                if seen.insert(canon(&t)) {
                    if is_s2(&t) {
                        s2 += 1;
                    }
                    if !last {
                        next.push(t);
                    }
                }
            }
        }
        let total = seen.len() as u64;
        let check = a000568
            .get(k - 1)
            .map(|&e| if e == total { "ok".into() } else { format!("MISMATCH {}", e) })
            .unwrap_or_else(|| "?".into());
        println!("{:2} | {:>12}  [{:>10}] | {}", k, total, check, s2);
        reps = next;
    }
}
