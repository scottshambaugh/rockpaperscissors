// Regular-game census (total / twin-free / modular-prime) by isomorph-free
// canonical augmentation, stratified by the common degree d. A regular game has
// every vertex out-degree = in-degree = d (so 2d decisive edges, n-1-2d ties).
// Same nauty canonical-parent engine as balanced.rs; only the prune and leaf test
// change: for a target d, a size-s partial game keeps a vertex iff o_i<=d, in_i<=d
// and (d-o_i)+(d-in_i) <= n-s (each future vertex can add at most one win or one
// loss to a given vertex). That band is closed under canonical parent, so growth
// stays isomorph-free; at the leaf it forces o_i==in_i==d exactly. Strata are
// disjoint (a game is d-regular for one d), so we run d = 1..(n-1)/2 and sum.
// Single-threaded; parallelism by process sharding (default libnauty not thread-safe).
//
//   build: gcc -O2 -c rust/balanced_shim.c -I/usr/include/x86_64-linux-gnu/nauty -o /tmp/bshim.o
//          rustc -O rust/regular.rs -o /tmp/regular -C link-args="/tmp/bshim.o -lnauty"
//   run  : /tmp/regular N SPLIT NSHARDS SHARD   -> PARTIAL total/twin_free/prime
//
// arc[i] has bit (1<<j) set iff i beats j (M[i,j]==+1); ties = neither bit.

use std::collections::HashSet;
use std::os::raw::c_int;

mod common;
use common::{connected, is_prime, paradoxical, twin_free, Arc, MAXN};

extern "C" {
    fn rps_canon(arc: *const u64, n: c_int, canong: *mut u64, lab: *mut c_int, orbits: *mut c_int);
}

#[derive(Clone)]
struct Game {
    arc: Arc,
    out: [i32; MAXN],
    ind: [i32; MAXN],
    k: usize,
    d: i32,
}



fn canon(arc: &Arc, n: usize) -> ([u64; MAXN], [i32; MAXN], [i32; MAXN]) {
    let mut canong = [0u64; MAXN];
    let mut lab = [0i32; MAXN];
    let mut orbits = [0i32; MAXN];
    unsafe {
        rps_canon(arc.as_ptr(), n as c_int, canong.as_mut_ptr(), lab.as_mut_ptr(), orbits.as_mut_ptr());
    }
    (canong, lab, orbits)
}

// One level of isomorph-free augmentation toward the d-regular stratum.
fn children(g: &Game, n: usize) -> Vec<Game> {
    let k = g.k;
    let d = g.d;
    let rem = (n - (k + 1)) as i32; // future vertices after this augmentation
    let total = 3usize.pow(k as u32);
    let mut seen: HashSet<[u64; MAXN]> = HashSet::new();
    let mut res = Vec::new();
    let mut v = [0i32; MAXN];
    for code in 0..total {
        let mut c = code;
        for i in 0..k {
            v[i] = (c % 3) as i32 - 1;
            c /= 3;
        }
        // child degrees + degree-band feasibility for target d
        let mut out = g.out;
        let mut ind = g.ind;
        let mut new_out = 0i32; // new node k's out-degree
        let mut new_in = 0i32;  // new node k's in-degree
        for i in 0..k {
            match v[i] {
                1 => { ind[i] += 1; new_out += 1; }   // new node k beats i
                -1 => { out[i] += 1; new_in += 1; }   // i beats new node k
                _ => {}
            }
        }
        out[k] = new_out;
        ind[k] = new_in;
        let mut ok = true;
        for i in 0..=k {
            if out[i] > d || ind[i] > d || (d - out[i]) + (d - ind[i]) > rem { ok = false; break; }
        }
        if !ok { continue; }
        // augment adjacency
        let mut arc = g.arc;
        let mut newbits = 0u64;
        for i in 0..k {
            if v[i] == -1 { arc[i] |= 1u64 << k; } else if v[i] == 1 { newbits |= 1u64 << i; }
        }
        arc[k] = newbits;
        // canonical-parent accept
        let (canong, lab, orbits) = canon(&arc, k + 1);
        let cc = lab[k] as usize;
        if cc != k && orbits[cc] != orbits[k] { continue; }
        let mut key = [0u64; MAXN];
        key[..k + 1].copy_from_slice(&canong[..k + 1]);
        if !seen.insert(key) { continue; }
        res.push(Game { arc, out, ind, k: k + 1, d });
    }
    res
}





fn descend(g: &Game, n: usize) -> (u64, u64, u64) {
    if g.k == n {
        if paradoxical(&g.arc, n) && connected(&g.arc, n) {
            return (1, twin_free(&g.arc, n) as u64, is_prime(&g.arc, n) as u64);
        }
        return (0, 0, 0);
    }
    let (mut t, mut f, mut p) = (0u64, 0u64, 0u64);
    for ch in children(g, n) {
        let (a, b, c) = descend(&ch, n);
        t += a; f += b; p += c;
    }
    (t, f, p)
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let n: usize = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(11);
    let split: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(6);
    let nshards: usize = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(1);
    let shard: usize = args.get(4).and_then(|s| s.parse().ok()).unwrap_or(0);

    // combined frontier across all d strata, each root tagged with its d
    let mut roots: Vec<Game> = Vec::new();
    let dmax = ((n - 1) / 2) as i32;
    for d in 1..=dmax {
        let g0 = Game { arc: [0u64; MAXN], out: [0i32; MAXN], ind: [0i32; MAXN], k: 1, d };
        let mut frontier = vec![g0];
        for _ in 1..split {
            let mut nxt = Vec::new();
            for g in &frontier { nxt.extend(children(g, n)); }
            frontier = nxt;
        }
        roots.extend(frontier);
    }
    eprintln!("[regular] n={} split={} shard={}/{} dmax={} roots={}", n, split, shard, nshards, dmax, roots.len());

    let (mut t, mut f, mut p) = (0u64, 0u64, 0u64);
    for (i, g) in roots.iter().enumerate() {
        if i % nshards != shard { continue; }
        let (a, b, c) = descend(g, n);
        t += a; f += b; p += c;
        if shard == 0 && i % (nshards * 64) == 0 {
            eprintln!("  shard0 root {}/{} (d={}) running total={}", i, roots.len(), g.d, t);
        }
    }
    println!("PARTIAL n={} shard={}/{} total={} twin_free={} prime={}", n, shard, nshards, t, f, p);
}
