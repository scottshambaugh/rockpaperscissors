// Within(C) = T(C) = # tournaments on V(C) such that every C-edge {x,y} has NO common
// in-neighbor (equivalently every out-neighborhood is a C-independent set).
// C = caterpillar given by spine leaf-counts: spine vertices s_0..s_{m-1} form a path,
// spine vertex i additionally carries leaf[i] leaves. Vertices numbered 0..k-1.
// Usage: within_brute L0 L1 L2 ...   (the leaf counts; number of args = spine length)
use std::env;

fn main() {
    let leaves: Vec<usize> = env::args().skip(1).map(|s| s.parse().unwrap()).collect();
    let m = leaves.len();
    // assign vertex ids: spine first (0..m), then leaves grouped by spine vertex
    let mut edges: Vec<(usize, usize)> = Vec::new();
    // spine path edges
    for i in 0..m.saturating_sub(1) { edges.push((i, i + 1)); }
    // leaves
    let mut next = m;
    let mut leafverts: Vec<Vec<usize>> = vec![Vec::new(); m];
    for i in 0..m {
        for _ in 0..leaves[i] { edges.push((i, next)); leafverts[i].push(next); next += 1; }
    }
    let k = next;
    if k > 26 { eprintln!("k={} too large for brute", k); return; }
    // pair index for the complete graph on k vertices
    let mut pa = Vec::new();
    for a in 0..k { for b in (a + 1)..k { pa.push((a, b)); } }
    let np = pa.len();
    let full: u32 = (1u32 << k) - 1;

    let mut count: u64 = 0;
    let total: u64 = 1u64 << np;
    for t in 0u64..total {
        // build beats
        let mut beats = [0u32; 26];
        for (j, &(a, b)) in pa.iter().enumerate() {
            if (t >> j) & 1 == 1 { beats[a] |= 1 << b; } else { beats[b] |= 1 << a; }
        }
        let mut bb = [0u32; 26];
        for i in 0..k { bb[i] = (full ^ (1 << i)) & !beats[i]; }
        // every C-edge: no common in-neighbor
        let mut ok = true;
        for &(x, y) in &edges {
            if bb[x] & bb[y] != 0 { ok = false; break; }
        }
        if ok { count += 1; }
    }
    // describe
    let desc: Vec<String> = leaves.iter().map(|x| x.to_string()).collect();
    println!("caterpillar leaves=[{}] k={} Within={}", desc.join(","), k, count);
}
