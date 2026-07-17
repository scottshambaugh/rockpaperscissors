// regular(n) by PURE COUNTING (no game enumeration) -- Rust port of
// burnside_regular.py, which validated against every known anchor (regular
// 3..12, A007079, A096368, and the regular_count strata) before producing
// regular(13) = 12050109962241.
//
// Method: iso classes of connected d-regular oriented games =
//   inverse Euler transform over m of  iso_all(m, d)
//   where iso_all(m, d) = (1/m!) sum over cycle types lambda of
//         (m!/z_lambda) * fix(lambda)  (Burnside),
// and fix(lambda) counts sigma-invariant labeled d-regular oriented graphs by
// an exchangeable-state DP: the fixed points of sigma form an exchangeable
// block (state = multiset of their partial (out, in) tallies), sigma's
// nontrivial cycles are super-vertices whose pair-orbits contribute in
// lcm-scaled degree jumps, and the antipodal orbit of an even cycle is forced
// to tie (orienting it would need both directions of one pair).
//
//   rustc -O rust/burnside_regular.rs -o /tmp/burn
//   /tmp/burn 13           # all strata + totals (validates vs known values)
//   /tmp/burn 14 4         # single stratum d=4 (for process-parallel runs)
//
// u128 throughout: the largest intermediates (m! * labeled count before the
// Burnside division) stay far below 2^127 for n <= 15.

use std::collections::HashMap;
use std::env;

const MAXN: usize = 16;

fn comb_table() -> [[u128; MAXN + 1]; MAXN + 1] {
    let mut c = [[0u128; MAXN + 1]; MAXN + 1];
    for i in 0..=MAXN {
        c[i][0] = 1;
        for j in 1..=i {
            c[i][j] = c[i - 1][j - 1] + if j <= i - 1 { c[i - 1][j] } else { 0 };
        }
    }
    c
}

fn factorial(n: usize) -> u128 {
    (1..=n as u128).product::<u128>().max(1)
}

// ---- state encodings (memo keys) ----
// fixed-point multiset: sorted (o,i) pairs, one byte each (o*16+i)
fn enc_f(fs: &[(u8, u8)]) -> Vec<u8> {
    fs.iter().map(|&(o, i)| o * 16 + i).collect()
}
// cycle state: sorted (len, o, i) triples
fn enc_c(cs: &[(u8, u8, u8)]) -> Vec<u8> {
    let mut v = Vec::with_capacity(cs.len() * 3);
    for &(l, o, i) in cs {
        v.push(l);
        v.push(o);
        v.push(i);
    }
    v
}

struct Ctx {
    d: u8,
    comb: [[u128; MAXN + 1]; MAXN + 1],
    f_memo: HashMap<Vec<u8>, u128>,
    g_memo: HashMap<Vec<u8>, u128>,
}

impl Ctx {
    fn new(d: u8) -> Self {
        Ctx { d, comb: comb_table(), f_memo: HashMap::new(), g_memo: HashMap::new() }
    }

    // ---- stage-1 DP: complete the pairs among an exchangeable multiset ----
    fn f(&mut self, state: &mut Vec<(u8, u8)>) -> u128 {
        state.sort_unstable();
        let key = enc_f(state);
        if let Some(&v) = self.f_memo.get(&key) {
            return v;
        }
        let val = self.f_inner(state);
        self.f_memo.insert(key, val);
        val
    }

    fn f_inner(&mut self, state: &Vec<(u8, u8)>) -> u128 {
        if state.is_empty() {
            return 1;
        }
        let (o0, i0) = state[0];
        let d = self.d;
        if o0 > d || i0 > d {
            return 0;
        }
        let (a, b) = ((d - o0) as usize, (d - i0) as usize);
        let rest = &state[1..];
        if rest.is_empty() {
            return if a == 0 && b == 0 { 1 } else { 0 };
        }
        // group rest by class
        let mut cls: Vec<((u8, u8), usize)> = Vec::new();
        for &t in rest {
            match cls.last_mut() {
                Some((c, n)) if *c == t => *n += 1,
                _ => cls.push((t, 1)),
            }
        }
        let mut total = 0u128;
        // distribute a outs (target gains an in) and b ins (target gains an out)
        let mut newrest: Vec<(u8, u8)> = Vec::with_capacity(rest.len());
        self.dist_f(&cls, 0, a, b, 1, &mut newrest, &mut total);
        total
    }

    fn dist_f(
        &mut self,
        cls: &[((u8, u8), usize)],
        ci: usize,
        ra: usize,
        rb: usize,
        acc: u128,
        newrest: &mut Vec<(u8, u8)>,
        total: &mut u128,
    ) {
        if ci == cls.len() {
            if ra == 0 && rb == 0 {
                let mut nr = newrest.clone();
                *total += acc * self.f(&mut nr);
            }
            return;
        }
        let ((o, i), nc) = cls[ci];
        let d = self.d;
        for x in 0..=ra.min(nc) {
            if x > 0 && i + 1 > d {
                break;
            }
            for y in 0..=rb.min(nc - x) {
                if y > 0 && o + 1 > d {
                    break;
                }
                let ways = self.comb[nc][x] * self.comb[nc - x][y];
                let base = newrest.len();
                for _ in 0..x {
                    newrest.push((o, i + 1));
                }
                for _ in 0..y {
                    newrest.push((o + 1, i));
                }
                for _ in 0..(nc - x - y) {
                    newrest.push((o, i));
                }
                self.dist_f(cls, ci + 1, ra - x, rb - y, acc * ways, newrest, total);
                newrest.truncate(base);
            }
        }
    }

    // ---- sigma-invariant count: cycles + fixed points ----
    fn g(&mut self, cstate: &mut Vec<(u8, u8, u8)>, fstate: &mut Vec<(u8, u8)>) -> u128 {
        cstate.sort_unstable();
        fstate.sort_unstable();
        if cstate.is_empty() {
            return self.f(fstate);
        }
        let mut key = enc_c(cstate);
        key.push(255);
        key.extend(enc_f(fstate));
        if let Some(&v) = self.g_memo.get(&key) {
            return v;
        }
        let (l, o, i) = cstate[0];
        let later: Vec<(u8, u8, u8)> = cstate[1..].to_vec();
        let d = self.d;
        let slots = ((l - 1) / 2) as usize;
        let mut total = 0u128;
        for t in 0..=slots {
            let (o1, i1) = (o + t as u8, i + t as u8);
            if o1 > d || i1 > d {
                break;
            }
            let ways = self.comb[slots][t] * (1u128 << t);
            total += ways * self.walk(l, o1, i1, &later, 0, fstate);
        }
        self.g_memo.insert(key, total);
        total
    }

    // decide orbits between the current cycle and later[idx..]
    fn walk(
        &mut self,
        l: u8,
        o: u8,
        i: u8,
        later: &Vec<(u8, u8, u8)>,
        idx: usize,
        fstate: &mut Vec<(u8, u8)>,
    ) -> u128 {
        if idx == later.len() {
            return self.close(l, o, i, later, fstate);
        }
        let d = self.d;
        let (l2, o2, i2) = later[idx];
        let g2 = gcd(l as usize, l2 as usize);
        let up_o = (l2 as usize / g2) as u8; // current-cycle per-vertex gain per orbit
        let up2 = (l as usize / g2) as u8; // later-cycle per-vertex gain per orbit
        let mut total = 0u128;
        for fw in 0..=g2 {
            let oo = o + (fw as u8) * up_o;
            if oo > d {
                break;
            }
            let ii2f = i2 + (fw as u8) * up2;
            if ii2f > d {
                break;
            }
            for bw in 0..=(g2 - fw) {
                let ii = i + (bw as u8) * up_o;
                if ii > d {
                    break;
                }
                let oo2 = o2 + (bw as u8) * up2;
                if oo2 > d {
                    break;
                }
                let ways = self.comb[g2][fw] * self.comb[g2 - fw][bw];
                let mut nl = later.clone();
                nl[idx] = (l2, oo2, ii2f);
                total += ways * self.walk(l, oo, ii, &nl, idx + 1, fstate);
            }
        }
        total
    }

    // distribute the cycle's remaining deficits over fixed points, then recurse
    fn close(
        &mut self,
        l: u8,
        o: u8,
        i: u8,
        later: &Vec<(u8, u8, u8)>,
        fstate: &mut Vec<(u8, u8)>,
    ) -> u128 {
        let d = self.d;
        if o > d || i > d {
            return 0;
        }
        let (need_o, need_i) = ((d - o) as usize, (d - i) as usize);
        let mut cls: Vec<((u8, u8), usize)> = Vec::new();
        fstate.sort_unstable();
        for &t in fstate.iter() {
            match cls.last_mut() {
                Some((c, n)) if *c == t => *n += 1,
                _ => cls.push((t, 1)),
            }
        }
        let mut total = 0u128;
        let mut newf: Vec<(u8, u8)> = Vec::with_capacity(fstate.len());
        self.dist_close(&cls, 0, need_o, need_i, 1, l, later, &mut newf, &mut total);
        total
    }

    #[allow(clippy::too_many_arguments)]
    fn dist_close(
        &mut self,
        cls: &[((u8, u8), usize)],
        ci: usize,
        ra: usize,
        rb: usize,
        acc: u128,
        l: u8,
        later: &Vec<(u8, u8, u8)>,
        newf: &mut Vec<(u8, u8)>,
        total: &mut u128,
    ) {
        if ci == cls.len() {
            if ra == 0 && rb == 0 {
                let mut nc = later.clone();
                let mut nf = newf.clone();
                *total += acc * self.g(&mut nc, &mut nf);
            }
            return;
        }
        let ((o, i), n_c) = cls[ci];
        let d = self.d;
        let xmax = if i + l <= d { ra.min(n_c) } else { 0 };
        for x in 0..=xmax {
            let ymax = if o + l <= d { rb.min(n_c - x) } else { 0 };
            for y in 0..=ymax {
                let ways = self.comb[n_c][x] * self.comb[n_c - x][y];
                let base = newf.len();
                for _ in 0..x {
                    newf.push((o, i + l)); // cycle -> fixed: fixed gains l ins
                }
                for _ in 0..y {
                    newf.push((o + l, i)); // fixed -> cycle: fixed gains l outs
                }
                for _ in 0..(n_c - x - y) {
                    newf.push((o, i));
                }
                self.dist_close(cls, ci + 1, ra - x, rb - y, acc * ways, l, later, newf, total);
                newf.truncate(base);
            }
        }
    }
}

fn gcd(a: usize, b: usize) -> usize {
    if b == 0 { a } else { gcd(b, a % b) }
}

fn partitions(n: usize, mx: usize, cur: &mut Vec<usize>, out: &mut Vec<Vec<usize>>) {
    if n == 0 {
        out.push(cur.clone());
        return;
    }
    for k in (1..=n.min(mx)).rev() {
        cur.push(k);
        partitions(n - k, k, cur, out);
        cur.pop();
    }
}

fn z_lambda(lam: &[usize]) -> u128 {
    let mut z = 1u128;
    let mut idx = 0;
    while idx < lam.len() {
        let l = lam[idx];
        let mut a = 0usize;
        while idx < lam.len() && lam[idx] == l {
            a += 1;
            idx += 1;
        }
        z *= (l as u128).pow(a as u32) * factorial(a);
    }
    z
}

fn iso_all(m: usize, ctx: &mut Ctx) -> u128 {
    let mut parts = Vec::new();
    partitions(m, m, &mut Vec::new(), &mut parts);
    let mfact = factorial(m);
    let mut total = 0u128;
    for lam in &parts {
        let cycles: Vec<u8> = lam.iter().filter(|&&l| l > 1).map(|&l| l as u8).collect();
        let nfixed = lam.iter().filter(|&&l| l == 1).count();
        let mut cstate: Vec<(u8, u8, u8)> = cycles.iter().map(|&l| (l, 0, 0)).collect();
        let mut fstate: Vec<(u8, u8)> = vec![(0, 0); nfixed];
        let fx = ctx.g(&mut cstate, &mut fstate);
        total += fx * (mfact / z_lambda(lam));
    }
    assert!(total % mfact == 0, "Burnside sum not divisible by m!");
    total / mfact
}

fn connected_row(nmax: usize, d: u8) -> Vec<u128> {
    let mut ctx = Ctx::new(d);
    let mut a = vec![0u128; nmax + 1];
    for m in 1..=nmax {
        a[m] = iso_all(m, &mut ctx);
        eprintln!("  d={} m={} iso_all={} (f-memo {}, g-memo {})", d, m, a[m], ctx.f_memo.len(), ctx.g_memo.len());
    }
    // inverse Euler transform (values can exceed i128 only far past n=15)
    let mut b = vec![0i128; nmax + 1];
    let mut c = vec![0i128; nmax + 1];
    for n in 1..=nmax {
        let mut bn = (n as i128) * (a[n] as i128);
        for k in 1..n {
            bn -= b[k] * (a[n - k] as i128);
        }
        b[n] = bn;
        let s: i128 = (1..n).filter(|dd| n % dd == 0).map(|dd| (dd as i128) * c[dd]).sum();
        assert!((b[n] - s) % (n as i128) == 0, "inverse Euler divisibility failed");
        c[n] = (b[n] - s) / (n as i128);
        assert!(c[n] >= 0);
    }
    c.iter().map(|&x| x as u128).collect()
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let n: usize = args.get(1).and_then(|s| s.parse().ok()).expect("usage: burn N [d]");
    assert!(n <= 15, "u128/i128 headroom validated only to n=15");
    let dmax = ((n - 1) / 2) as u8;
    let dsel: Option<u8> = args.get(2).and_then(|s| s.parse().ok());
    let mut grand = vec![0u128; n + 1];
    for d in 1..=dmax {
        if let Some(ds) = dsel {
            if d != ds {
                continue;
            }
        }
        let c = connected_row(n, d);
        let row: Vec<String> = c[1..].iter().map(|x| x.to_string()).collect();
        println!("d={}: connected iso by n: [{}]", d, row.join(", "));
        for m in 1..=n {
            grand[m] += c[m];
        }
    }
    if dsel.is_none() {
        let row: Vec<String> = grand[1..].iter().map(|x| x.to_string()).collect();
        println!("regular(n) totals: [{}]", row.join(", "));
    }
}
