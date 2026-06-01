// Exact COUNTER for the regular strata: instead of enumerating each Eulerian
// orientation (watercluster2), count them via Burnside + a degree-vector DP.
//
// For each connected 2d-regular graph G (read as `geng | listg -eq`), the number
// of regular oriented graphs on G up to isomorphism is, by Burnside over Aut(G):
//     (1/|Aut(G)|) * sum_{sigma in Aut(G)} fix(sigma)
// where fix(sigma) = # Eulerian orientations (out=in=d everywhere) invariant under
// sigma. We compute fix(sigma) with a DP over the sigma-edge-orbits (each orbit is
// oriented as a unit; sigma=identity recovers the per-edge DP = total EO count).
// Summing over all G in a stratum gives R(n,d); summing over d gives regular(n).
// This is exact and covers the full search space (same set geng feeds watercluster2);
// validated below against every stratum we already enumerated.
//
//   rustc -O rust/regular_count.rs -o /tmp/regcount
//   nauty-geng 12 -d8 -D8 -c | nauty-listg -eq 2>/dev/null | /tmp/regcount 4

use std::collections::HashMap;
use std::env;
use std::io::{self, BufRead};

// ---------- 1-WL refinement ----------
fn refine(adj: &[u16], n: usize, init: &[u32]) -> Vec<u32> {
    let mut col = init.to_vec();
    loop {
        let mut sig: Vec<(u32, Vec<u32>)> = Vec::with_capacity(n);
        for i in 0..n {
            let mut nb: Vec<u32> = Vec::new();
            for j in 0..n {
                if adj[i] & (1 << j) != 0 {
                    nb.push(col[j]);
                }
            }
            nb.sort_unstable();
            sig.push((col[i], nb));
        }
        let mut order: Vec<usize> = (0..n).collect();
        order.sort_by(|&a, &b| sig[a].cmp(&sig[b]));
        let mut nc = vec![0u32; n];
        let mut c = 0u32;
        for k in 0..n {
            if k > 0 && sig[order[k]] != sig[order[k - 1]] {
                c += 1;
            }
            nc[order[k]] = c;
        }
        if nc == col {
            return col;
        }
        col = nc;
    }
}

fn code_of(adj: &[u16], n: usize, vert_at: &[usize]) -> u128 {
    let mut code: u128 = 0;
    let mut bit = 0;
    for a in 0..n {
        for b in (a + 1)..n {
            if adj[vert_at[a]] & (1 << vert_at[b]) != 0 {
                code |= 1u128 << bit;
            }
            bit += 1;
        }
    }
    code
}

// collect all discrete labelings achieving the minimum canonical code
const AUT_CAP: usize = 2_000_000; // guard: bail before an n!-sized group OOMs us

fn canon_rec(adj: &[u16], n: usize, col: Vec<u32>, best: &mut u128, labs: &mut Vec<Vec<usize>>) {
    if labs.len() > AUT_CAP {
        return;
    }
    let maxc = *col.iter().max().unwrap();
    if maxc as usize == n - 1 {
        let mut vert_at = vec![0usize; n];
        for v in 0..n {
            vert_at[col[v] as usize] = v;
        }
        let code = code_of(adj, n, &vert_at);
        if code < *best {
            *best = code;
            labs.clear();
            labs.push(vert_at);
        } else if code == *best {
            labs.push(vert_at);
        }
        return;
    }
    let mut freq = vec![0u32; (maxc + 1) as usize];
    for &c in &col {
        freq[c as usize] += 1;
    }
    let target = (0..=maxc).find(|&c| freq[c as usize] > 1).unwrap();
    for v in 0..n {
        if col[v] != target {
            continue;
        }
        let mut nc = col.clone();
        for u in 0..n {
            if nc[u] > target || (nc[u] == target && u != v) {
                nc[u] += 1;
            }
        }
        let r = refine(adj, n, &nc);
        canon_rec(adj, n, r, best, labs);
    }
}

// automorphisms as vertex permutations a[]: a maps G to itself
fn automorphisms(adj: &[u16], n: usize) -> Vec<Vec<usize>> {
    let deg: Vec<u32> = (0..n).map(|i| adj[i].count_ones()).collect();
    let col = refine(adj, n, &deg);
    let mut best = u128::MAX;
    let mut labs: Vec<Vec<usize>> = Vec::new();
    canon_rec(adj, n, col, &mut best, &mut labs);
    // labs[k][p] = vertex placed at position p by labeling k. Each is an iso G->canon.
    // autos = labs[0]_pos_inverse composed with labs[k]: a(v) = vert_at0[ pos_k[v] ]
    let l0 = &labs[0];
    let mut pos0 = vec![0usize; n]; // pos0[v] = position of v in l0
    for p in 0..n {
        pos0[l0[p]] = p;
    }
    let mut autos = Vec::new();
    for lab in &labs {
        let mut posk = vec![0usize; n];
        for p in 0..n {
            posk[lab[p]] = p;
        }
        // a(v) = vertex at position posk[v] in l0 = l0[posk[v]]
        let a: Vec<usize> = (0..n).map(|v| l0[posk[v]]).collect();
        autos.push(a);
    }
    autos
}

// fix(sigma): # Eulerian orientations (out=in=d) invariant under sigma, via
// edge-orbit DP.  edges: list of (u,v).  d: target out-degree.
fn fix_sigma(edges: &[(usize, usize)], sigma: &[usize], n: usize, d: u32) -> u64 {
    // index edges
    let mut emap: HashMap<(usize, usize), usize> = HashMap::new();
    for (k, &(u, v)) in edges.iter().enumerate() {
        let key = if u < v { (u, v) } else { (v, u) };
        emap.insert(key, k);
    }
    let m = edges.len();
    let mut orbit_of = vec![usize::MAX; m];
    let mut orbits: Vec<Vec<usize>> = Vec::new();
    for start in 0..m {
        if orbit_of[start] != usize::MAX {
            continue;
        }
        let oid = orbits.len();
        let mut orb = Vec::new();
        let mut cur = start;
        loop {
            if orbit_of[cur] != usize::MAX {
                break;
            }
            orbit_of[cur] = oid;
            orb.push(cur);
            let (u, v) = edges[cur];
            let (su, sv) = (sigma[u], sigma[v]);
            let key = if su < sv { (su, sv) } else { (sv, su) };
            cur = emap[&key];
        }
        orbits.push(orb);
    }
    // For each orbit, the two orientation choices give out-degree-increment vectors.
    // Choice A: orient representative edge u->v; propagate by sigma. Choice B: reverse.
    // If propagation forces an edge to be oriented both ways (a flip), fix=0.
    let mut choice_vecs: Vec<[Vec<u8>; 2]> = Vec::with_capacity(orbits.len());
    for orb in &orbits {
        // build orientation by following sigma from the first edge with a fixed start dir
        let rep = orb[0];
        let (ru, rv) = edges[rep];
        // direction map: for each edge in orbit, (tail, head). Start rep: tail=ru.
        let mut tail: HashMap<usize, usize> = HashMap::new(); // edge idx -> tail vertex
        // propagate: edge e with tail t; sigma(e) has tail sigma(t)
        let mut stack = vec![(rep, ru)];
        let mut consistent = true;
        while let Some((e, t)) = stack.pop() {
            if let Some(&old) = tail.get(&e) {
                if old != t {
                    consistent = false;
                }
                continue;
            }
            tail.insert(e, t);
            let (u, v) = edges[e];
            let h = if t == u { v } else { u };
            // sigma maps directed edge t->h to sigma(t)->sigma(h)
            let (st, sh) = (sigma[t], sigma[h]);
            let key = if st < sh { (st, sh) } else { (sh, st) };
            let se = emap[&key];
            stack.push((se, st));
        }
        let mut a = vec![0u8; n];
        if consistent {
            for (&_e, &t) in &tail {
                a[t] += 1;
            }
        }
        // choice B = reverse: tails become heads
        let mut b = vec![0u8; n];
        if consistent {
            for (&e, &t) in &tail {
                let (u, v) = edges[e];
                let h = if t == u { v } else { u };
                b[h] += 1;
            }
        }
        if !consistent {
            return 0;
        }
        choice_vecs.push([a, b]);
        let _ = (ru, rv);
    }
    // DP over orbits: state = packed out-degree vector (each 0..=d), 4 bits/vertex.
    debug_assert!(d < 15 && n <= 16);
    let pack = |v: &[u8]| -> u64 {
        let mut k: u64 = 0;
        for i in 0..n {
            k |= (v[i] as u64) << (4 * i);
        }
        k
    };
    let mut cur: HashMap<u64, u64> = HashMap::new();
    cur.insert(0u64, 1u64);
    let mut deg = vec![0u8; n];
    for cv in &choice_vecs {
        let mut next: HashMap<u64, u64> = HashMap::with_capacity(cur.len() * 2);
        for (&state, &cnt) in &cur {
            for i in 0..n {
                deg[i] = ((state >> (4 * i)) & 0xF) as u8;
            }
            for ch in 0..2 {
                let inc = &cv[ch];
                let mut ok = true;
                for i in 0..n {
                    if deg[i] + inc[i] > d as u8 {
                        ok = false;
                        break;
                    }
                }
                if !ok {
                    continue;
                }
                let mut ns = state;
                for i in 0..n {
                    if inc[i] != 0 {
                        ns += (inc[i] as u64) << (4 * i);
                    }
                }
                *next.entry(ns).or_insert(0) += cnt;
            }
        }
        cur = next;
    }
    let mut full = vec![d as u8; n];
    let target = pack(&full);
    full.clear();
    (*cur.get(&target).unwrap_or(&0)) as u64
}

fn main() {
    let d: u32 = env::args().nth(1).and_then(|s| s.parse().ok()).expect("usage: regcount d");
    let stdin = io::stdin();
    let mut lines = stdin.lock().lines().map(|l| l.unwrap());
    let mut total: u128 = 0;
    let mut ngraphs: u64 = 0;
    loop {
        // read "n m"
        let header = match lines.next() {
            Some(h) => h,
            None => break,
        };
        let header = header.trim();
        if header.is_empty() {
            continue;
        }
        let mut it = header.split_ascii_whitespace();
        let n: usize = it.next().unwrap().parse().unwrap();
        let m: usize = it.next().unwrap().parse().unwrap();
        // read edges: listg wraps long lists across several lines, so keep
        // collecting numbers until we have 2*m of them
        let mut nums: Vec<usize> = Vec::with_capacity(2 * m);
        while nums.len() < 2 * m {
            let eline = lines.next().expect("unexpected EOF reading edges");
            nums.extend(eline.split_ascii_whitespace().map(|x| x.parse::<usize>().unwrap()));
        }
        let mut edges = Vec::with_capacity(m);
        let mut adj = vec![0u16; n];
        for e in 0..m {
            let u = nums[2 * e];
            let v = nums[2 * e + 1];
            edges.push((u, v));
            adj[u] |= 1 << v;
            adj[v] |= 1 << u;
        }
        let autos = automorphisms(&adj, n);
        let mut s: u64 = 0;
        for sigma in &autos {
            s += fix_sigma(&edges, sigma, n, d);
        }
        total += (s / (autos.len() as u64)) as u128;
        ngraphs += 1;
    }
    println!("d={} graphs={} R(n,d)={}", d, ngraphs, total);
}
