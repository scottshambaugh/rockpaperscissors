// WEIGHTED-balanced twin-free core enumerator, for the balanced (twin-free)
// bracket complement count. A non-twin-free balanced n-game collapses (maximal
// tie-twin classes -> one vertex each) to a unique twin-free CORE C with
// multiplicities m_v >= 1, sum m = n, some m >= 2; the blow-up is balanced iff
// C is m-WEIGHT-balanced: sum_u m_u rel(v, u) = 0 for every v. This enumerates
// (C, m) pairs up to color-preserving isomorphism (weights = vertex colors,
// nauty rps_canon_colored) for ONE weight vector, sorted descending so every
// prefix is a valid parent family; sum over the weight partitions of n gives
// the non-twin-free count, and  tf(n) = balanced(n) - non_tf(n).
// Leaf conditions: weighted balance (by construction) + paradoxical +
// connected + twin-free, all on the core (each is equivalent to the same
// property of the blow-up).
//
//   build: gcc -O2 -c rust/balanced_shim.c -I$NAUTY_INC -o /tmp/bshim.o
//          rustc -O rust/wbal.rs -o /tmp/wbal -C link-args="/tmp/bshim.o -lnauty"
//   run  : /tmp/wbal 2,1,1,1,1,1,1,1,1 [SPLIT NSHARDS SHARD]
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

// weight vector, set once in main before any enumeration
static mut WEIGHTS: [i32; MAXN] = [1; MAXN];

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

// One level of isomorph-free augmentation: accepted, per-parent-deduped children.
// The extension vectors v in {-1,0,1}^k are enumerated by a DFS over the
// per-coordinate FEASIBLE sets (|rowsum[i] - v[i]| <= lim, with suffix bounds on
// the achievable total s in [-lim, lim]) instead of scanning all 3^k codes --
// at the last level lim = 0 pins every coordinate (at most ONE candidate), where
// the flat scan burned 3^(n-1) iterations per parent to find it.
fn children(g: &Game, n: usize) -> Vec<Game> {
    let k = g.k;
    // remaining weight after this augmentation = how much any weighted rowsum
    // can still move; the added vertex has weight w_k
    let (wk, lim) = unsafe {
        let wk = WEIGHTS[k];
        let mut s = 0i32;
        for t in (k + 1)..n {
            s += WEIGHTS[t];
        }
        (wk, s)
    };
    // parent degrees, for the (od,id) prefilter below
    let mut pod = [0u8; MAXN];
    let mut pid = [0u8; MAXN];
    for i in 0..k {
        pod[i] = g.arc[i].count_ones() as u8;
        let mut r = g.arc[i];
        while r != 0 {
            let j = r.trailing_zeros() as usize;
            r &= r - 1;
            pid[j] += 1;
        }
    }
    let mut seen: HashSet<[u64; MAXN]> = HashSet::new();
    let mut out = Vec::new();
    // allowed values per coordinate, and suffix min/max of the achievable sum
    let mut allow = [[false; 3]; MAXN]; // index a+1 for a in {-1,0,1}
    let mut smin = [0i32; MAXN + 1];
    let mut smax = [0i32; MAXN + 1];
    // v[i] = -1 (i beats k): R(i) += wk, R(k) -= w_i
    // v[i] = +1 (k beats i): R(i) -= wk, R(k) += w_i
    for i in 0..k {
        for (ai, a) in [-1i32, 0, 1].iter().enumerate() {
            let r = g.rowsum[i] - a * wk; // R(i) after choosing v[i] = a... sign below
            let _ = r;
            let eff = match ai {
                0 => wk,  // v = -1
                1 => 0,
                _ => -wk, // v = +1
            };
            let r2 = g.rowsum[i] + eff;
            allow[i][ai] = -lim <= r2 && r2 <= lim;
        }
        if !allow[i].iter().any(|&b| b) {
            return out; // some weighted row can no longer reach 0
        }
    }
    // suffix bounds on the ADDED vertex's weighted rowsum R(k) = sum of
    // (+w_i for v=+1, -w_i for v=-1); coordinate i contributes in [-w_i, +w_i]
    let wi: Vec<i32> = (0..k).map(|i| unsafe { WEIGHTS[i] }).collect();
    for i in (0..k).rev() {
        let lo = if allow[i][0] { -wi[i] } else { 0 }.min(if allow[i][2] { wi[i] } else { 0 });
        let hi = if allow[i][2] { wi[i] } else { 0 }.max(if allow[i][0] { -wi[i] } else { 0 });
        smin[i] = smin[i + 1] + lo;
        smax[i] = smax[i + 1] + hi;
    }
    let mut v = [0i32; MAXN];
    let mut stack_i = 0usize;
    // iterative DFS over coordinates; va[i] = next value-index (0..3) to try
    let mut va = [0usize; MAXN];
    let mut part = [0i32; MAXN + 1]; // partial sum before coordinate i
    loop {
        if stack_i == k {
            let s = part[k];
            // (s is within [-lim, lim] by the suffix pruning)
            //
            // Canonical-parent rule: accept iff the added vertex k lies in the
            // DESIGNATED orbit of the child -- the orbit of the canonically-last
            // vertex among those with lexicographically-maximal (outdeg, indeg).
            // Any isomorphism-invariant orbit selector gives exactly one accepted
            // (parent, orbit) per child class (deleting any vertex of a feasible
            // partial game stays feasible), and this one admits a ~free prefilter:
            // the child's degrees derive from the parent's plus v, so candidates
            // whose added vertex is NOT (od,id)-maximal are rejected with no
            // nauty call at all -- which is most of them.
            let odk = (0..k).filter(|&i| v[i] == 1).count() as u8;
            let idk = (0..k).filter(|&i| v[i] == -1).count() as u8;
            let mut maximal = true;
            for i in 0..k {
                if unsafe { WEIGHTS[i] } != wk {
                    continue; // different color class: never competes
                }
                let oi = pod[i] + (v[i] == -1) as u8;
                let ii = pid[i] + (v[i] == 1) as u8;
                if oi > odk || (oi == odk && ii > idk) {
                    maximal = false;
                    break;
                }
            }
            if !maximal {
                if stack_i == 0 { break; }
                stack_i -= 1;
                continue;
            }
            // Stage 2: among (od,id)-ties, compare the FULL degree-refinement
            // signature (sig_cmp_with) -- the designated orbit is now "orbit of
            // the canonically-last vertex among signature-maximal vertices".
            // Rejects most stage-1 survivors without nauty; only exact signature
            // ties reach the canon call.
            let mut beats = [0u16; 16];
            let mut cod = [0u8; 16];
            let mut cid = [0u8; 16];
            for i in 0..k {
                let mut row = (g.arc[i] & ((1u64 << k) - 1)) as u16;
                if v[i] == -1 { row |= 1 << k; }
                beats[i] = row;
                cod[i] = pod[i] + (v[i] == -1) as u8;
                cid[i] = pid[i] + (v[i] == 1) as u8;
            }
            let mut newbits16 = 0u16;
            for i in 0..k {
                if v[i] == 1 { newbits16 |= 1 << i; }
            }
            beats[k] = newbits16;
            cod[k] = odk;
            cid[k] = idk;
            let mut eqmask = 0u32; // (od,id)+signature ties with k (same color only)
            let mut sig_max = true;
            for i in 0..k {
                if unsafe { WEIGHTS[i] } != wk {
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
            if !sig_max {
                if stack_i == 0 { break; }
                stack_i -= 1;
                continue;
            }
            let mut arc = g.arc;
            let mut newbits = 0u64;
            for i in 0..k {
                if v[i] == -1 { arc[i] |= 1u64 << k; } else if v[i] == 1 { newbits |= 1u64 << i; }
            }
            arc[k] = newbits;
            let (canong, lab, orbits) = canon(&arc, k + 1);
            // w = the signature-maximal vertex with the largest canonical position
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
            let accept = orbits[w] == orbits[k];
            if accept {
                // dedup isomorphic children from this same parent
                let mut key = [0u64; MAXN];
                key[..k + 1].copy_from_slice(&canong[..k + 1]);
                if seen.insert(key) {
                    let mut rowsum = g.rowsum;
                    for i in 0..k { rowsum[i] -= v[i] * wk; }
                    rowsum[k] = s;
                    out.push(Game { arc, rowsum, k: k + 1 });
                }
            }
            // backtrack
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
            let ps = part[i] + a * wi[i]; // contribution to R(k): +w_i for v=+1
            // suffix bound: can R(k) still land in [-lim, lim]?
            if ps + smin[i + 1] > lim || ps + smax[i + 1] < -lim { continue; }
            v[i] = a;
            part[i + 1] = ps;
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
    out
}





static mut QMODE: bool = false;

// paradox for the quotient use: the special vertex 0 may get its wins/losses
// from INSIDE the module it stands for, so only vertices >= 1 are required to
// have a win and a loss here; the substitution driver re-checks the composite.
fn paradox_skip0(arc: &Arc, n: usize) -> bool {
    let mut beaten = 0u64;
    for i in 0..n {
        if i > 0 && arc[i] == 0 {
            return false;
        }
        beaten |= arc[i];
    }
    (beaten | 1) == full_mask_local(n)
}
fn full_mask_local(n: usize) -> u64 {
    if n >= 64 { u64::MAX } else { (1u64 << n) - 1 }
}

fn descend(g: &Game, n: usize, out: &mut Vec<u8>) -> u64 {
    if g.k == n {
        let qmode = unsafe { QMODE };
        let pass = if qmode {
            paradox_skip0(&g.arc, n) && connected(&g.arc, n) && twin_free(&g.arc, n)
        } else {
            // cores must be twin-free (else the collapse was not maximal);
            // paradox/connected of the blow-up == paradox/connected of the core
            paradoxical(&g.arc, n) && connected(&g.arc, n) && twin_free(&g.arc, n)
        };
        if pass && qmode {
            // emit digraph6 (vertex 0 = special, enumeration order preserved)
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
    let wspec = args.get(1).expect("usage: wbal w1,w2,... [split nshards shard]");
    let weights: Vec<i32> = wspec.split(',').map(|s| s.parse().unwrap()).collect();
    let n = weights.len();
    assert!(n >= 2 && n <= 16);
    // descending order => every prefix's weight multiset is well-defined,
    // making prefix families valid canonical parents
    for i in 1..n {
        assert!(weights[i] <= weights[i - 1], "weights must be sorted descending");
    }
    unsafe {
        for (i, &w) in weights.iter().enumerate() {
            WEIGHTS[i] = w;
        }
    }
    let split: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(5);
    let nshards: usize = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(1);
    let shard: usize = args.get(4).and_then(|s| s.parse().ok()).unwrap_or(0);
    if args.get(5).map(|s| s == "q").unwrap_or(false) {
        unsafe { QMODE = true; }
    }

    let g0 = Game { arc: [0u64; MAXN], rowsum: [0i32; MAXN], k: 1 };
    let mut frontier = vec![g0];
    for _ in 1..split.min(n) {
        let mut nxt = Vec::new();
        for g in &frontier { nxt.extend(children(g, n)); }
        frontier = nxt;
    }
    eprintln!("[wbal] w={} split={} shard={}/{} roots={}", wspec, split, shard, nshards, frontier.len());

    let mut t = 0u64;
    let mut outbuf: Vec<u8> = Vec::new();
    for (i, g) in frontier.iter().enumerate() {
        if i % nshards != shard { continue; }
        t += descend(g, n, &mut outbuf);
    }
    eprintln!("WBAL w={} shard={}/{} cores={}", wspec, shard, nshards, t);
    if !unsafe { QMODE } {
        println!("WBAL w={} shard={}/{} cores={}", wspec, shard, nshards, t);
    }
}
