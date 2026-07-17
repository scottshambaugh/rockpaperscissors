// Balanced-game census (total / twin-free / modular-prime) by isomorph-free
// canonical augmentation, with nauty (densenauty, via balanced_shim.c) for the
// canonical-parent test. Mirrors the validated Python _generate_balanced/descend:
// every balanced n-class has a unique canonical-parent chain, so we grow one rep
// per class and tally the leaves. Single-threaded; parallelism is by process
// sharding (idx % nshards) since the default libnauty is not thread-safe.
//
//   build: gcc -O2 -c rust/balanced_shim.c -I/usr/include/x86_64-linux-gnu/nauty -o /tmp/bshim.o
//          rustc -O rust/balanced.rs -o /tmp/balanced -C link-args="/tmp/bshim.o -lnauty"
//   run  : /tmp/balanced N SPLIT NSHARDS SHARD   -> prints PARTIAL total/tf/prime
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
    rowsum: [i32; MAXN],
    k: usize,
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

// One level of isomorph-free augmentation: accepted, per-parent-deduped children.
fn children(g: &Game, n: usize) -> Vec<Game> {
    let k = g.k;
    let lim = (n - (k + 1)) as i32;
    let total = 3usize.pow(k as u32);
    let mut seen: HashSet<[u64; MAXN]> = HashSet::new();
    let mut out = Vec::new();
    let mut v = [0i32; MAXN];
    for code in 0..total {
        let mut c = code;
        for i in 0..k {
            v[i] = (c % 3) as i32 - 1;
            c /= 3;
        }
        // balance feasibility (closed under canonical parent)
        let mut ok = true;
        for i in 0..k {
            let r = g.rowsum[i] - v[i];
            if r < -lim || r > lim { ok = false; break; }
        }
        if !ok { continue; }
        let s: i32 = (0..k).map(|i| v[i]).sum();
        if s < -lim || s > lim { continue; }
        // augment
        let mut arc = g.arc;
        let mut newbits = 0u64;
        for i in 0..k {
            if v[i] == -1 { arc[i] |= 1u64 << k; } else if v[i] == 1 { newbits |= 1u64 << i; }
        }
        arc[k] = newbits;
        // canonical-parent accept: added node k is in the canonical-last orbit
        let (canong, lab, orbits) = canon(&arc, k + 1);
        let cc = lab[k] as usize;
        if cc != k && orbits[cc] != orbits[k] { continue; }
        // dedup isomorphic children from this same parent
        let mut key = [0u64; MAXN];
        key[..k + 1].copy_from_slice(&canong[..k + 1]);
        if !seen.insert(key) { continue; }
        let mut rowsum = g.rowsum;
        for i in 0..k { rowsum[i] -= v[i]; }
        rowsum[k] = s;
        out.push(Game { arc, rowsum, k: k + 1 });
    }
    out
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
    let n: usize = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(9);
    let split: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(6);
    let nshards: usize = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(1);
    let shard: usize = args.get(4).and_then(|s| s.parse().ok()).unwrap_or(0);

    let g0 = Game { arc: [0u64; MAXN], rowsum: [0i32; MAXN], k: 1 };
    let mut frontier = vec![g0];
    for _ in 1..split {
        let mut nxt = Vec::new();
        for g in &frontier { nxt.extend(children(g, n)); }
        frontier = nxt;
    }
    eprintln!("[balanced] n={} split={} shard={}/{} roots={}", n, split, shard, nshards, frontier.len());

    let (mut t, mut f, mut p) = (0u64, 0u64, 0u64);
    for (i, g) in frontier.iter().enumerate() {
        if i % nshards != shard { continue; }
        let (a, b, c) = descend(g, n);
        t += a; f += b; p += c;
        if shard == 0 && i % (nshards * 64) == 0 {
            eprintln!("  shard0 root {}/{}  running total={}", i, frontier.len(), t);
        }
    }
    println!("PARTIAL n={} shard={}/{} total={} twin_free={} prime={}", n, shard, nshards, t, f, p);
}
