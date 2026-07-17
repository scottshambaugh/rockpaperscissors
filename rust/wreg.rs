// WEIGHTED-regular twin-free core enumerator, for the regular (twin-free) and
// [prime] bracket complement counts (same complement method as wbal.rs, regular
// flavor). A non-twin-free regular n-game collapses to a twin-free core C with
// multiplicities m (weights); the blow-up is d-regular iff every core vertex has
// WEIGHTED out-sum = weighted in-sum = d (the tie condition is then automatic).
// Enumerates (C, m) up to color-preserving iso (weights = colors) for one weight
// vector sorted descending; leaf requires connected + twin-free (paradox is
// implied by d >= 1). Optional quotient mode: vertex 0 gets its own target
// d0 = d - d_module for the prime-gap substitution construction.
//
//   build: gcc -O2 -c rust/balanced_shim.c -I$NAUTY_INC -o /tmp/bshim.o
//          rustc -O rust/wreg.rs -o /tmp/wreg -C link-args="/tmp/bshim.o -lnauty"
//   run  : /tmp/wreg 2,1,1,1,1,1,1,1,1,1,1 D [split nshards shard]        # cores
//          /tmp/wreg t,1,...,1 D d0 q                                     # quotients
//
// arc[i] has bit (1<<j) set iff i beats j (M[i,j]==+1); ties = neither bit.

use std::collections::HashSet;
use std::os::raw::c_int;

mod common;
use common::{connected, is_prime, paradoxical, sig_cmp_with, twin_free, Arc, MAXN};

extern "C" {
    fn rps_canon_colored(
        arc: *const u64,
        n: c_int,
        col: *const c_int,
        canong: *mut u64,
        lab: *mut c_int,
        orbits: *mut c_int,
    );
}

static mut WEIGHTS: [i32; MAXN] = [1; MAXN];
static mut D0: i32 = -1; // quotient mode: vertex 0's separate degree target
static mut QMODE: bool = false;

fn target(v: usize, d: i32) -> i32 {
    unsafe {
        if v == 0 && D0 >= 0 {
            D0
        } else {
            d
        }
    }
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
        rps_canon_colored(
            arc.as_ptr(),
            n as c_int,
            WEIGHTS.as_ptr(),
            canong.as_mut_ptr(),
            lab.as_mut_ptr(),
            orbits.as_mut_ptr(),
        );
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
    let (wk, rem) = unsafe {
        let wk = WEIGHTS[k];
        let mut s = 0i32;
        for t in (k + 1)..n {
            s += WEIGHTS[t];
        }
        (wk, s)
    };
    let mut seen: HashSet<[u64; MAXN]> = HashSet::new();
    let mut res = Vec::new();
    let mut allow = [[false; 3]; MAXN]; // index a+1 for a in {-1,0,1}
    for i in 0..k {
        let ti = target(i, d);
        let deficit = (ti - g.out[i]) + (ti - g.ind[i]);
        // every branch must leave the child's deficit within the remaining
        // weight (with weights the root can start OVER budget, so unlike the
        // unweighted enumerator this cannot be left to induction from v=0)
        allow[i][0] = g.out[i] + wk <= ti && deficit - wk <= rem;
        allow[i][1] = deficit <= rem;
        allow[i][2] = g.ind[i] + wk <= ti && deficit - wk <= rem;
        if !allow[i].iter().any(|&b| b) {
            return res; // vertex i can no longer reach its target
        }
    }
    // suffix weight sum of coordinates that can still contribute a nonzero, for
    // the new vertex's deficit bound: p1 + m1 must reach at least 2*tk - rem
    let wi: Vec<i32> = (0..k).map(|i| unsafe { WEIGHTS[i] }).collect();
    let mut sufnz = [0i32; MAXN + 1];
    for i in (0..k).rev() {
        sufnz[i] = sufnz[i + 1] + if allow[i][0] || allow[i][2] { wi[i] } else { 0 };
    }
    let tk = target(k, d);
    let need = 2 * tk - rem;
    let mut v = [0i32; MAXN];
    let mut va = [0usize; MAXN];
    let mut p1 = [0i32; MAXN + 1]; // running WEIGHTED sum of +1 (new vertex's out)
    let mut m1 = [0i32; MAXN + 1]; // running WEIGHTED sum of -1 (new vertex's in)
    let mut stack_i = 0usize;
    loop {
        if stack_i == k {
            let new_out = p1[k];
            let new_in = m1[k];
            // stage 1: weighted-(od,id) maximality among SAME-COLOR vertices
            let odk = new_out as u8;
            let idk = new_in as u8;
            let mut cod = [0u8; 16];
            let mut cid = [0u8; 16];
            let mut maximal = true;
            for i in 0..k {
                cod[i] = (g.out[i] + if v[i] == -1 { wk } else { 0 }) as u8;
                cid[i] = (g.ind[i] + if v[i] == 1 { wk } else { 0 }) as u8;
                if wi[i] != wk {
                    continue;
                }
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
                    if wi[i] != wk {
                        continue;
                    }
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
                                if v[i] == -1 { out[i] += wk; } else if v[i] == 1 { ind[i] += wk; }
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
            let (np, nm) = (
                p1[i] + if a == 1 { wi[i] } else { 0 },
                m1[i] + if a == -1 { wi[i] } else { 0 },
            );
            if np > tk || nm > tk { continue; }
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

fn descend(g: &Game, n: usize, out: &mut Vec<u8>) -> u64 {
    if g.k == n {
        let pass = connected(&g.arc, n) && twin_free(&g.arc, n);
        if pass && unsafe { QMODE } {
            let mut beats = [0u16; 16];
            for i in 0..n {
                let mut r = g.arc[i];
                while r != 0 {
                    let j = r.trailing_zeros() as usize;
                    r &= r - 1;
                    beats[i] |= 1 << j;
                }
            }
            common::encode(&beats[..n], n, out);
            use std::io::Write;
            std::io::stdout().write_all(out).unwrap();
        }
        return pass as u64;
    }
    let mut t = 0u64;
    for ch in children(g, n) {
        t += descend(&ch, n, out);
    }
    t
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let wspec = args.get(1).expect("usage: wreg w1,w2,... D [split nshards shard | D0 q]");
    let weights: Vec<i32> = wspec.split(',').map(|s| s.parse().unwrap()).collect();
    let n = weights.len();
    assert!(n >= 2 && n <= 16);
    for i in 1..n {
        assert!(weights[i] <= weights[i - 1], "weights must be sorted descending");
    }
    unsafe {
        for (i, &w) in weights.iter().enumerate() {
            WEIGHTS[i] = w;
        }
    }
    let d: i32 = args.get(2).and_then(|s| s.parse().ok()).expect("degree target D");
    // quotient mode: wreg w D d0 q
    if args.get(4).map(|s| s == "q").unwrap_or(false) {
        unsafe {
            D0 = args.get(3).and_then(|s| s.parse().ok()).expect("d0");
            QMODE = true;
        }
    }
    let (split, nshards, shard) = if unsafe { QMODE } {
        (5usize, 1usize, 0usize)
    } else {
        (
            args.get(3).and_then(|s| s.parse().ok()).unwrap_or(5),
            args.get(4).and_then(|s| s.parse().ok()).unwrap_or(1),
            args.get(5).and_then(|s| s.parse().ok()).unwrap_or(0),
        )
    };

    let g0 = Game { arc: [0u64; MAXN], out: [0i32; MAXN], ind: [0i32; MAXN], k: 1, d };
    let mut frontier = vec![g0];
    for _ in 1..split.min(n) {
        let mut nxt = Vec::new();
        for g in &frontier { nxt.extend(children(g, n)); }
        frontier = nxt;
    }
    eprintln!("[wreg] w={} d={} split={} shard={}/{} roots={}", wspec, d, split, shard, nshards, frontier.len());

    let mut t = 0u64;
    let mut outbuf: Vec<u8> = Vec::new();
    for (i, g) in frontier.iter().enumerate() {
        if i % nshards != shard { continue; }
        t += descend(g, n, &mut outbuf);
    }
    eprintln!("WREG w={} d={} shard={}/{} cores={}", wspec, d, shard, nshards, t);
    if !unsafe { QMODE } {
        println!("WREG w={} d={} shard={}/{} cores={}", wspec, d, shard, nshards, t);
    }
}
