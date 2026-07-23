// S2 (two-paradox, gamma>=3) tournament ISO count via Burnside over cycle types.
//   S2_iso(n) = (1/n!) sum_{lambda all-odd} (n!/z_lambda) * S(lambda)
//   S(lambda) = # of pi-invariant labeled tournaments (pi of type lambda) that are S2
//             = # orientations of the E(lambda) edge-orbits giving gamma>=3.
// For each odd partition with E(lambda) <= THRESH we compute S(lambda) by enumerating
// 2^E orientations. Terms with E > THRESH are reported as the residual "wall".
//
// Usage: s2_burnside N [THRESH]   (default THRESH=24)

use std::collections::HashMap;

fn gcd(a: u64, b: u64) -> u64 { if b == 0 { a } else { gcd(b, a % b) } }
fn factorial(n: u64) -> u128 { (1..=n).map(|x| x as u128).product::<u128>().max(1) }
fn z_lambda(parts: &[u64]) -> u128 {
    let mut mult = HashMap::new();
    for &p in parts { *mult.entry(p).or_insert(0u64) += 1; }
    let mut z: u128 = 1;
    for (&k, &m) in &mult { z *= (k as u128).pow(m as u32); z *= factorial(m); }
    z
}
fn partitions(n: u64, max: u64, cur: &mut Vec<u64>, out: &mut Vec<Vec<u64>>) {
    if n == 0 { out.push(cur.clone()); return; }
    let hi = max.min(n);
    let mut p = hi;
    while p >= 1 { cur.push(p); partitions(n - p, p, cur, out); cur.pop(); p -= 1; }
}

// build permutation of given cycle type on 0..n-1
fn build_perm(parts: &[u64], n: usize) -> Vec<usize> {
    let mut perm = vec![0usize; n];
    let mut v = 0usize;
    for &l in parts {
        let l = l as usize;
        for i in 0..l { perm[v + i] = v + (i + 1) % l; }
        v += l;
    }
    perm
}

// compute edge orbits: returns Vec of orbits, each orbit a Vec of (a,b) with a<b
fn edge_orbits(perm: &[usize], n: usize) -> Vec<Vec<(usize, usize)>> {
    let idx = |a: usize, b: usize| -> (usize, usize) { if a < b { (a, b) } else { (b, a) } };
    let mut seen = vec![false; n * n];
    let mut orbits = Vec::new();
    for a in 0..n {
        for b in (a + 1)..n {
            if seen[a * n + b] { continue; }
            let mut orb = Vec::new();
            let (mut ca, mut cb) = (a, b);
            loop {
                let (x, y) = idx(ca, cb);
                if seen[x * n + y] { break; }
                seen[x * n + y] = true;
                orb.push((x, y));
                let (na, nb) = (perm[ca], perm[cb]);
                ca = na; cb = nb;
            }
            orbits.push(orb);
        }
    }
    orbits
}

// test S2 (gamma>=3): every pair has a common dominator (some w beats both)
fn is_s2(beats: &[u32], n: usize, full: u32) -> bool {
    let mut beatenby = [0u32; 16];
    for i in 0..n { beatenby[i] = (full ^ (1 << i)) & !beats[i]; }
    for u in 0..n {
        for v in (u + 1)..n {
            if beatenby[u] & beatenby[v] == 0 { return false; }
        }
    }
    true
}

// count S(lambda) by enumerating 2^(#orbits) orientation assignments.
// orientation bit for orbit r: 0 => rep (a<b) oriented a beats b ; 1 => b beats a.
fn count_s_lambda(orbits: &[Vec<(usize, usize)>], n: usize) -> u128 {
    let full: u32 = if n == 32 { u32::MAX } else { (1u32 << n) - 1 };
    let m = orbits.len();
    // For each orbit, precompute two edge-sets? We rebuild beats per assignment.
    // To orient consistently: within an orbit, the "rep direction" propagates.
    // We stored orbit as pairs (x,y) x<y reached by iterating perm from the rep.
    // Orientation must be consistent under perm: pick rep = orbit[0]=(a0,b0).
    // If bit=0: a0 beats b0. The k-th pair was reached by applying perm k times to
    // (a0,b0) as an ORDERED pair (ca,cb) starting (a0,b0). So the ordered direction
    // at step k is (perm^k a0) beats (perm^k b0). We must recover ordered orientation.
    // Rebuild via ordered traversal.
    let mut ordered: Vec<Vec<(usize, usize)>> = Vec::with_capacity(m); // (winner,loser) if bit=0
    // reconstruct ordered traversal per orbit
    for orb in orbits {
        let (a0, b0) = orb[0];
        let mut seq = Vec::with_capacity(orb.len());
        let (mut ca, mut cb) = (a0, b0);
        for _ in 0..orb.len() {
            seq.push((ca, cb)); // bit=0 => ca beats cb
            ca = PERM_APPLY.with(|p| p.borrow()[ca]);
            cb = PERM_APPLY.with(|p| p.borrow()[cb]);
        }
        ordered.push(seq);
    }
    let mut count: u128 = 0;
    let mut beats = vec![0u32; n];
    for assign in 0u64..(1u64 << m) {
        for x in beats.iter_mut() { *x = 0; }
        for r in 0..m {
            let bit = (assign >> r) & 1;
            for &(ca, cb) in &ordered[r] {
                if bit == 0 { beats[ca] |= 1 << cb; } else { beats[cb] |= 1 << ca; }
            }
        }
        if is_s2(&beats, n, full) { count += 1; }
    }
    count
}

use std::cell::RefCell;
thread_local!(static PERM_APPLY: RefCell<Vec<usize>> = RefCell::new(Vec::new()));

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let n: usize = args[1].parse().unwrap();
    let thresh: usize = args.get(2).map(|s| s.parse().unwrap()).unwrap_or(24);

    let mut parts_list = Vec::new();
    partitions(n as u64, n as u64, &mut Vec::new(), &mut parts_list);
    let nfac = factorial(n as u64);

    let mut computed_total: u128 = 0; // sum over computed lambda of (n!/z)*S(lambda)
    let mut residual: Vec<(usize, u128, Vec<u64>)> = Vec::new(); // (E, coeff, parts)
    let mut all_cheap = true;

    for parts in &parts_list {
        if parts.iter().any(|&l| l % 2 == 0) { continue; }
        let perm = build_perm(parts, n);
        let orbits = edge_orbits(&perm, n);
        let e = orbits.len();
        let coeff = nfac / z_lambda(parts);
        if e <= thresh {
            PERM_APPLY.with(|p| *p.borrow_mut() = perm.clone());
            let s = count_s_lambda(&orbits, n);
            computed_total += coeff * s;
        } else {
            all_cheap = false;
            residual.push((e, coeff, parts.clone()));
        }
    }

    residual.sort();
    if all_cheap {
        let iso = computed_total / nfac;
        let rem = computed_total % nfac;
        println!("n={} THRESH={}  ALL terms computed.  S2_iso = {}  (exact_div={})",
                 n, thresh, iso, rem == 0);
    } else {
        println!("n={} THRESH={}  computed 12-cheap partial sum (n! * ...) = {}",
                 n, thresh, computed_total);
        println!("  residual (uncomputed) terms, E > {}:", thresh);
        for (e, coeff, parts) in &residual {
            println!("    E={:2}  coeff(n!/z)={:<20}  lambda={:?}", e, coeff, parts);
        }
        println!("  --> S2_iso(n) = (computed_partial + sum_residual coeff*S(lambda)) / {}!", n);
    }
}
