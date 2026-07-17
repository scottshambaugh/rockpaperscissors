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
use common::{connected, is_prime, paradoxical, sig_cmp_with, twin_free, Arc, MAXN};

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
// Same three-stage engine as balanced.rs: (1) DFS over the per-coordinate
// FEASIBLE extension values only -- +1 (new beats i) needs ind[i] < d, -1 needs
// out[i] < d, 0 needs vertex i's remaining deficit to fit in the remaining
// vertices, and the new vertex's own deficit gives a suffix lower bound on the
// number of nonzeros -- (2) reject candidates whose added vertex is not
// (od,id)-maximal, then not signature-maximal (no nauty call for either), and
// (3) canon only on exact signature ties, accepting iff the added vertex sits
// in the orbit of the canonically-last signature-maximal vertex.
fn children(g: &Game, n: usize) -> Vec<Game> {
    let k = g.k;
    let d = g.d;
    let rem = (n - (k + 1)) as i32;
    let mut seen: HashSet<[u64; MAXN]> = HashSet::new();
    let mut res = Vec::new();
    let mut allow = [[false; 3]; MAXN]; // index a+1 for a in {-1,0,1}
    for i in 0..k {
        let deficit = (d - g.out[i]) + (d - g.ind[i]);
        allow[i][0] = g.out[i] < d;
        allow[i][1] = deficit <= rem;
        allow[i][2] = g.ind[i] < d;
        if !allow[i].iter().any(|&b| b) {
            return res; // vertex i can no longer reach out = in = d
        }
    }
    // suffix count of coordinates that can still contribute a nonzero, for the
    // new vertex's deficit bound: p1 + m1 must reach at least 2d - rem
    let mut sufnz = [0i32; MAXN + 1];
    for i in (0..k).rev() {
        sufnz[i] = sufnz[i + 1] + ((allow[i][0] || allow[i][2]) as i32);
    }
    let need = 2 * d - rem;
    let mut v = [0i32; MAXN];
    let mut va = [0usize; MAXN];
    let mut p1 = [0i32; MAXN + 1]; // running count of +1 (new vertex's out-degree)
    let mut m1 = [0i32; MAXN + 1]; // running count of -1 (new vertex's in-degree)
    let mut stack_i = 0usize;
    loop {
        if stack_i == k {
            let new_out = p1[k];
            let new_in = m1[k];
            // stage 1: (od,id) maximality of the added vertex, from parent degrees
            let odk = new_out as u8;
            let idk = new_in as u8;
            let mut cod = [0u8; 16];
            let mut cid = [0u8; 16];
            let mut maximal = true;
            for i in 0..k {
                cod[i] = (g.out[i] + (v[i] == -1) as i32) as u8;
                cid[i] = (g.ind[i] + (v[i] == 1) as i32) as u8;
                if cod[i] > odk || (cod[i] == odk && cid[i] > idk) {
                    maximal = false;
                    break;
                }
            }
            if maximal {
                cod[k] = odk;
                cid[k] = idk;
                // stage 2: full-signature maximality among (od,id) ties
                let mut beats = [0u16; 16];
                for i in 0..k {
                    let mut row = (g.arc[i] & ((1u64 << k) - 1)) as u16;
                    if v[i] == -1 { row |= 1 << k; }
                    beats[i] = row;
                }
                let mut nb16 = 0u16;
                for i in 0..k {
                    if v[i] == 1 { nb16 |= 1 << i; }
                }
                beats[k] = nb16;
                let mut eqmask = 0u32;
                let mut sig_max = true;
                for i in 0..k {
                    if cod[i] == odk && cid[i] == idk {
                        match sig_cmp_with(&beats[..k + 1], k + 1, &cod, &cid, i, k) {
                            std::cmp::Ordering::Greater => { sig_max = false; break; }
                            std::cmp::Ordering::Equal => { eqmask |= 1 << i; }
                            std::cmp::Ordering::Less => {}
                        }
                    }
                }
                if sig_max {
                    // stage 3: canon on ties; accept iff k is in the designated orbit
                    let mut arc = g.arc;
                    let mut newbits = 0u64;
                    for i in 0..k {
                        if v[i] == -1 { arc[i] |= 1u64 << k; } else if v[i] == 1 { newbits |= 1u64 << i; }
                    }
                    arc[k] = newbits;
                    let (canong, lab, orbits) = canon(&arc, k + 1);
                    let mut pos = [0i32; MAXN];
                    for (p, &vv) in lab.iter().enumerate().take(k + 1) {
                        pos[vv as usize] = p as i32;
                    }
                    let mut w = k;
                    let mut m = eqmask;
                    while m != 0 {
                        let i = m.trailing_zeros() as usize;
                        m &= m - 1;
                        if pos[i] > pos[w] {
                            w = i;
                        }
                    }
                    if orbits[w] == orbits[k] {
                        let mut key = [0u64; MAXN];
                        key[..k + 1].copy_from_slice(&canong[..k + 1]);
                        if seen.insert(key) {
                            let mut out = g.out;
                            let mut ind = g.ind;
                            for i in 0..k {
                                if v[i] == -1 { out[i] += 1; } else if v[i] == 1 { ind[i] += 1; }
                            }
                            out[k] = new_out;
                            ind[k] = new_in;
                            res.push(Game { arc, out, ind, k: k + 1, d });
                        }
                    }
                }
            }
            if stack_i == 0 { break; }
            stack_i -= 1;
            continue;
        }
        let i = stack_i;
        let mut advanced = false;
        while va[i] < 3 {
            let ai = va[i];
            va[i] += 1;
            if !allow[i][ai] { continue; }
            let a = ai as i32 - 1;
            let (np, nm) = (p1[i] + (a == 1) as i32, m1[i] + (a == -1) as i32);
            if np > d || nm > d { continue; }
            if np + nm + sufnz[i + 1] < need { continue; }
            v[i] = a;
            p1[i + 1] = np;
            m1[i + 1] = nm;
            stack_i += 1;
            va[stack_i] = 0;
            advanced = true;
            break;
        }
        if !advanced {
            if i == 0 { break; }
            stack_i -= 1;
        }
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
