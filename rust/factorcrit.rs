// Factor-critical filter for graph6 streams: pass through only graphs G on n
// (odd) vertices such that G - v has a perfect matching for EVERY vertex v.
//
//   rustc -O rust/factorcrit.rs -o /tmp/fc
//   nauty-geng 9 2>/dev/null | /tmp/fc 9 | nauty-directg -o 2>/dev/null | /tmp/cmf 9
//
// Why this is exact-lossless in front of the completely-mixed filter: a game is
// completely mixed only if all n principal Pfaffians Pf(M-i) are nonzero, and a
// skew Pfaffian is a signed sum over the perfect matchings of its support graph
// -- so Pf(M-i) != 0 forces G-i to have a perfect matching. Needing that for all
// i is precisely G being factor-critical. Factor-critical also implies connected
// and min-degree >= 2 (a degree-1 vertex u's neighbor v isolates u in G-v), so
// this single test subsumes geng's -c -d2 and prunes strictly more, at the cheap
// undirected level -- before directg's 2^(edges) orientation blow-up.
//
// Perfect-matching test on the <=15 remaining vertices: subset-DP over a bitmask
// (dp[S] = some vertex set S can be perfectly matched using edges within S),
// pinned to always match the lowest free vertex to cut the branching.

use std::env;
use std::io::{self, Read, Write, BufWriter};

// Does the induced subgraph on vertex set `mask` (even popcount) have a perfect
// matching? adj[v] = neighbor bitmask over all n vertices.
fn has_pm(adj: &[u16], mask: u16) -> bool {
    // memoized recursion: match the lowest set vertex to each available neighbor.
    // n<=15 so a full 2^n table is at most 32768 bytes -- allocate per call is
    // wasteful, so use an explicit stack-free recursion with a small cache.
    fn rec(adj: &[u16], mask: u16, memo: &mut [i8]) -> bool {
        if mask == 0 {
            return true;
        }
        let c = memo[mask as usize];
        if c >= 0 {
            return c != 0;
        }
        let v = mask.trailing_zeros() as usize;
        let rest = mask & !(1u16 << v);
        let mut cand = adj[v] & rest;
        let mut ok = false;
        while cand != 0 {
            let u = cand.trailing_zeros();
            cand &= cand - 1;
            if rec(adj, rest & !(1u16 << u), memo) {
                ok = true;
                break;
            }
        }
        memo[mask as usize] = ok as i8;
        ok
    }
    let mut memo = vec![-1i8; 1usize << 16];
    rec(adj, mask, &mut memo)
}

fn factor_critical(adj: &[u16], n: usize) -> bool {
    let full: u16 = ((1u32 << n) - 1) as u16;
    for v in 0..n {
        if !has_pm(adj, full & !(1u16 << v)) {
            return false;
        }
    }
    true
}

// Decode one graph6 line (n <= 62, no header byte here since n<63) into adj.
fn decode_g6(line: &[u8], n: usize, adj: &mut [u16]) {
    for a in adj.iter_mut() {
        *a = 0;
    }
    // graph6 for n<63: first byte is n+63, then upper-triangle bits column-major:
    // bit for edge (i,j), i<j, ordered j=1: (0,1); j=2: (0,2),(1,2); ...
    let data = &line[1..];
    let mut bitpos = 0usize;
    for j in 1..n {
        for i in 0..j {
            let byte = data[bitpos / 6];
            let bit = 5 - (bitpos % 6);
            if (byte - 63) & (1 << bit) != 0 {
                adj[i] |= 1 << j;
                adj[j] |= 1 << i;
            }
            bitpos += 1;
        }
    }
}

fn main() {
    let n: usize = env::args().nth(1).and_then(|s| s.parse().ok()).expect("usage: fc n");
    assert!(n >= 3 && n <= 15 && n % 2 == 1, "n must be odd, 3..=15");
    let count_only = env::args().nth(2).as_deref() == Some("-u");
    let mut input = String::new();
    io::stdin().read_to_string(&mut input).unwrap();
    let out = io::stdout();
    let mut w = BufWriter::new(out.lock());
    let mut adj = vec![0u16; n];
    let (mut total, mut kept) = (0u64, 0u64);
    for line in input.lines() {
        let b = line.as_bytes();
        if b.is_empty() {
            continue;
        }
        // skip any header/aux lines geng might emit (they start with '>')
        if b[0] == b'>' {
            continue;
        }
        debug_assert_eq!(b[0] as usize, n + 63, "graph6 order mismatch");
        total += 1;
        decode_g6(b, n, &mut adj);
        if factor_critical(&adj, n) {
            kept += 1;
            if !count_only {
                w.write_all(b).unwrap();
                w.write_all(b"\n").unwrap();
            }
        }
    }
    w.flush().unwrap();
    eprintln!("fc n={}: {} graphs in, {} factor-critical", n, total, kept);
}
