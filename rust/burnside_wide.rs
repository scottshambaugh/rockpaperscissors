// WIDE (256-bit) build of burnside_regular.rs for n = 16, 17 -- identical
// algorithm, exact U256 arithmetic (u128 overflows past n = 15), plus a 2x
// memo saving from arc-reversal symmetry: reversing every arc maps a valid
// completion to a valid completion with (out, in) swapped everywhere, so f/g
// keys are canonicalized against their mirror image.
// Validated by reproducing the u128 build's n = 13..15 outputs exactly.
//
// Original header: regular(n) by PURE COUNTING -- Rust port of
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

// ---- minimal exact U256 (little-endian u64 limbs) ----
#[derive(Clone, Copy, PartialEq, Eq)]
struct W([u64; 4]);

impl W {
    const ZERO: W = W([0; 4]);
    const ONE: W = W([1, 0, 0, 0]);
    fn from_u64(x: u64) -> W {
        W([x, 0, 0, 0])
    }
    fn is_zero(&self) -> bool {
        self.0 == [0; 4]
    }
    fn add(&self, o: &W) -> W {
        let mut r = [0u64; 4];
        let mut carry = 0u128;
        for i in 0..4 {
            let s = self.0[i] as u128 + o.0[i] as u128 + carry;
            r[i] = s as u64;
            carry = s >> 64;
        }
        assert!(carry == 0, "U256 add overflow");
        W(r)
    }
    fn sub(&self, o: &W) -> W {
        let mut r = [0u64; 4];
        let mut borrow = 0i128;
        for i in 0..4 {
            let s = self.0[i] as i128 - o.0[i] as i128 - borrow;
            if s < 0 {
                r[i] = (s + (1i128 << 64)) as u64;
                borrow = 1;
            } else {
                r[i] = s as u64;
                borrow = 0;
            }
        }
        assert!(borrow == 0, "U256 sub underflow");
        W(r)
    }
    fn mul(&self, o: &W) -> W {
        let mut acc = [0u128; 8];
        for i in 0..4 {
            if self.0[i] == 0 {
                continue;
            }
            for j in 0..4 {
                if o.0[j] == 0 {
                    continue;
                }
                let p = self.0[i] as u128 * o.0[j] as u128;
                let k = i + j;
                acc[k] += p & 0xFFFF_FFFF_FFFF_FFFF;
                acc[k + 1] += p >> 64;
            }
        }
        let mut r = [0u64; 4];
        let mut carry = 0u128;
        for k in 0..8 {
            let s = acc[k] + carry;
            if k < 4 {
                r[k] = s as u64;
            } else {
                assert!(s as u64 == 0, "U256 mul overflow");
            }
            carry = s >> 64;
        }
        assert!(carry == 0, "U256 mul overflow (carry)");
        W(r)
    }
    fn mul_u64(&self, x: u64) -> W {
        self.mul(&W::from_u64(x))
    }
    // divide by a u64, assert exact when required by caller via remainder
    fn divmod_u64(&self, x: u64) -> (W, u64) {
        assert!(x != 0);
        let mut r = [0u64; 4];
        let mut rem = 0u128;
        for i in (0..4).rev() {
            let cur = (rem << 64) | self.0[i] as u128;
            r[i] = (cur / x as u128) as u64;
            rem = cur % x as u128;
        }
        (W(r), rem as u64)
    }
    fn cmp_w(&self, o: &W) -> std::cmp::Ordering {
        for i in (0..4).rev() {
            match self.0[i].cmp(&o.0[i]) {
                std::cmp::Ordering::Equal => continue,
                x => return x,
            }
        }
        std::cmp::Ordering::Equal
    }
    fn to_string(&self) -> String {
        if self.is_zero() {
            return "0".into();
        }
        let mut parts: Vec<u64> = Vec::new();
        let mut cur = *self;
        while !cur.is_zero() {
            let (q, rem) = cur.divmod_u64(1_000_000_000_000_000_000);
            parts.push(rem);
            cur = q;
        }
        let mut s = parts.pop().unwrap().to_string();
        while let Some(p) = parts.pop() {
            s.push_str(&format!("{:018}", p));
        }
        s
    }
}

const MAXN: usize = 16;

fn comb_table() -> [[u64; MAXN + 1]; MAXN + 1] {
    let mut c = [[0u64; MAXN + 1]; MAXN + 1];
    for i in 0..=MAXN {
        c[i][0] = 1;
        for j in 1..=i {
            c[i][j] = c[i - 1][j - 1] + if j <= i - 1 { c[i - 1][j] } else { 0 };
        }
    }
    c
}

fn factorial(n: usize) -> u64 {
    (1..=n as u64).product::<u64>().max(1)
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
    comb: [[u64; MAXN + 1]; MAXN + 1],
    f_memo: HashMap<Vec<u8>, W>,
    g_memo: HashMap<Vec<u8>, W>,
}

impl Ctx {
    fn new(d: u8) -> Self {
        Ctx { d, comb: comb_table(), f_memo: HashMap::new(), g_memo: HashMap::new() }
    }

    // ---- stage-1 DP: complete the pairs among an exchangeable multiset ----
    fn f(&mut self, state: &mut Vec<(u8, u8)>) -> W {
        state.sort_unstable();
        let mut mirror: Vec<(u8, u8)> = state.iter().map(|&(o, i)| (i, o)).collect();
        mirror.sort_unstable();
        let key = enc_f(state).min(enc_f(&mirror));
        if let Some(&v) = self.f_memo.get(&key) {
            return v;
        }
        let val = self.f_inner(state);
        self.f_memo.insert(key, val);
        val
    }

    fn f_inner(&mut self, state: &Vec<(u8, u8)>) -> W {
        if state.is_empty() {
            return W::ONE;
        }
        let (o0, i0) = state[0];
        let d = self.d;
        if o0 > d || i0 > d {
            return W::ZERO;
        }
        let (a, b) = ((d - o0) as usize, (d - i0) as usize);
        let rest = &state[1..];
        if rest.is_empty() {
            return if a == 0 && b == 0 { W::ONE } else { W::ZERO };
        }
        // group rest by class
        let mut cls: Vec<((u8, u8), usize)> = Vec::new();
        for &t in rest {
            match cls.last_mut() {
                Some((c, n)) if *c == t => *n += 1,
                _ => cls.push((t, 1)),
            }
        }
        let mut total = W::ZERO;
        // distribute a outs (target gains an in) and b ins (target gains an out)
        let mut newrest: Vec<(u8, u8)> = Vec::with_capacity(rest.len());
        self.dist_f(&cls, 0, a, b, W::ONE, &mut newrest, &mut total);
        total
    }

    fn dist_f(
        &mut self,
        cls: &[((u8, u8), usize)],
        ci: usize,
        ra: usize,
        rb: usize,
        acc: W,
        newrest: &mut Vec<(u8, u8)>,
        total: &mut W,
    ) {
        if ci == cls.len() {
            if ra == 0 && rb == 0 {
                let mut nr = newrest.clone();
                let fv = self.f(&mut nr);
                *total = total.add(&acc.mul(&fv));
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
                let ways = self.comb[nc][x] as u128 * self.comb[nc - x][y] as u128;
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
                let nacc = acc.mul(&W([ways as u64, (ways >> 64) as u64, 0, 0]));
                self.dist_f(cls, ci + 1, ra - x, rb - y, nacc, newrest, total);
                newrest.truncate(base);
            }
        }
    }

    // ---- sigma-invariant count: cycles + fixed points ----
    fn g(&mut self, cstate: &mut Vec<(u8, u8, u8)>, fstate: &mut Vec<(u8, u8)>) -> W {
        cstate.sort_unstable();
        fstate.sort_unstable();
        if cstate.is_empty() {
            return self.f(fstate);
        }
        let mut cm: Vec<(u8, u8, u8)> = cstate.iter().map(|&(l, o, i)| (l, i, o)).collect();
        cm.sort_unstable();
        let mut fm: Vec<(u8, u8)> = fstate.iter().map(|&(o, i)| (i, o)).collect();
        fm.sort_unstable();
        let mut key = enc_c(cstate);
        key.push(255);
        key.extend(enc_f(fstate));
        let mut key2 = enc_c(&cm);
        key2.push(255);
        key2.extend(enc_f(&fm));
        let key = key.min(key2);
        if let Some(&v) = self.g_memo.get(&key) {
            return v;
        }
        let (l, o, i) = cstate[0];
        let later: Vec<(u8, u8, u8)> = cstate[1..].to_vec();
        let d = self.d;
        let slots = ((l - 1) / 2) as usize;
        let mut total = W::ZERO;
        for t in 0..=slots {
            let (o1, i1) = (o + t as u8, i + t as u8);
            if o1 > d || i1 > d {
                break;
            }
            let ways = W::from_u64(self.comb[slots][t]).mul_u64(1u64 << t);
            let wv = self.walk(l, o1, i1, &later, 0, fstate);
            total = total.add(&ways.mul(&wv));
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
    ) -> W {
        if idx == later.len() {
            return self.close(l, o, i, later, fstate);
        }
        let d = self.d;
        let (l2, o2, i2) = later[idx];
        let g2 = gcd(l as usize, l2 as usize);
        let up_o = (l2 as usize / g2) as u8; // current-cycle per-vertex gain per orbit
        let up2 = (l as usize / g2) as u8; // later-cycle per-vertex gain per orbit
        let mut total = W::ZERO;
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
                let ways = W::from_u64(self.comb[g2][fw]).mul_u64(self.comb[g2 - fw][bw]);
                let mut nl = later.clone();
                nl[idx] = (l2, oo2, ii2f);
                let wv = self.walk(l, oo, ii, &nl, idx + 1, fstate);
                total = total.add(&ways.mul(&wv));
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
    ) -> W {
        let d = self.d;
        if o > d || i > d {
            return W::ZERO;
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
        let mut total = W::ZERO;
        let mut newf: Vec<(u8, u8)> = Vec::with_capacity(fstate.len());
        self.dist_close(&cls, 0, need_o, need_i, W::ONE, l, later, &mut newf, &mut total);
        total
    }

    #[allow(clippy::too_many_arguments)]
    fn dist_close(
        &mut self,
        cls: &[((u8, u8), usize)],
        ci: usize,
        ra: usize,
        rb: usize,
        acc: W,
        l: u8,
        later: &Vec<(u8, u8, u8)>,
        newf: &mut Vec<(u8, u8)>,
        total: &mut W,
    ) {
        if ci == cls.len() {
            if ra == 0 && rb == 0 {
                let mut nc = later.clone();
                let mut nf = newf.clone();
                let gv = self.g(&mut nc, &mut nf);
                *total = total.add(&acc.mul(&gv));
            }
            return;
        }
        let ((o, i), n_c) = cls[ci];
        let d = self.d;
        let xmax = if i + l <= d { ra.min(n_c) } else { 0 };
        for x in 0..=xmax {
            let ymax = if o + l <= d { rb.min(n_c - x) } else { 0 };
            for y in 0..=ymax {
                let ways = W::from_u64(self.comb[n_c][x]).mul_u64(self.comb[n_c - x][y]);
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
                self.dist_close(cls, ci + 1, ra - x, rb - y, acc.mul(&ways), l, later, newf, total);
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

fn z_lambda(lam: &[usize]) -> u64 {
    let mut z = 1u64;
    let mut idx = 0;
    while idx < lam.len() {
        let l = lam[idx];
        let mut a = 0usize;
        while idx < lam.len() && lam[idx] == l {
            a += 1;
            idx += 1;
        }
        z *= (l as u64).pow(a as u32) * factorial(a);
    }
    z
}

// Burnside integrality gate + modular fingerprint (see burnside_regular.rs /
// README): a residual prime factor > m cannot come from an m-smooth
// multiplicity, so it flags a wrong per-cycle Fix or an overflowed accumulator.
fn assert_wide_integral(r: u64, m: usize) {
    if r == 0 {
        return;
    }
    let (mut x, mut p, mut f) = (r, 2u64, Vec::new());
    while p * p <= x {
        while x % p == 0 {
            f.push(p);
            x /= p;
        }
        p += 1;
    }
    if x > 1 {
        f.push(x);
    }
    eprintln!("BURNSIDE-WIDE GATE FAILED (m={m}): residual = {r} = {f:?}");
    if let Some(&big) = f.iter().find(|&&q| q as usize > m) {
        eprintln!("  prime factor {big} > m={m}: not from any m-smooth multiplicity -- a per-cycle Fix or accumulator is wrong.");
    }
    panic!("Burnside sum not divisible by m!");
}

fn iso_all(m: usize, ctx: &mut Ctx) -> W {
    let mut parts = Vec::new();
    partitions(m, m, &mut Vec::new(), &mut parts);
    let mfact = factorial(m);
    let mut total = W::ZERO;
    for lam in &parts {
        let cycles: Vec<u8> = lam.iter().filter(|&&l| l > 1).map(|&l| l as u8).collect();
        let nfixed = lam.iter().filter(|&&l| l == 1).count();
        let mut cstate: Vec<(u8, u8, u8)> = cycles.iter().map(|&l| (l, 0, 0)).collect();
        let mut fstate: Vec<(u8, u8)> = vec![(0, 0); nfixed];
        let fx = ctx.g(&mut cstate, &mut fstate);
        total = total.add(&fx.mul_u64(mfact / z_lambda(lam)));
    }
    let (q, r) = total.divmod_u64(mfact);
    assert_wide_integral(r, m);
    q
}

fn connected_row(nmax: usize, d: u8) -> Vec<W> {
    let mut ctx = Ctx::new(d);
    let mut a = vec![W::ZERO; nmax + 1];
    for m in 1..=nmax {
        a[m] = iso_all(m, &mut ctx);
        eprintln!("  d={} m={} iso_all={} (f-memo {}, g-memo {})", d, m, a[m].to_string(), ctx.f_memo.len(), ctx.g_memo.len());
    }
    // inverse Euler transform, unsigned-safe: subtrahends never exceed the
    // minuend because b(n) = sum_{dd|n} dd*c(dd) >= 0 and c(n) >= 0
    let mut b = vec![W::ZERO; nmax + 1];
    let mut c = vec![W::ZERO; nmax + 1];
    for n in 1..=nmax {
        let mut sub = W::ZERO;
        for k in 1..n {
            sub = sub.add(&b[k].mul(&a[n - k]));
        }
        let lead = a[n].mul_u64(n as u64);
        assert!(sub.cmp_w(&lead) != std::cmp::Ordering::Greater, "Euler subtrahend exceeds minuend");
        b[n] = lead.sub(&sub);
        let mut s = W::ZERO;
        for dd in (1..n).filter(|dd| n % dd == 0) {
            s = s.add(&c[dd].mul_u64(dd as u64));
        }
        assert!(s.cmp_w(&b[n]) != std::cmp::Ordering::Greater);
        let (q, r) = b[n].sub(&s).divmod_u64(n as u64);
        assert!(r == 0, "inverse Euler divisibility failed");
        c[n] = q;
    }
    c
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let n: usize = args.get(1).and_then(|s| s.parse().ok()).expect("usage: burn N [d]");
    assert!(n <= 17, "U256 headroom checked to n=17");
    let dmax = ((n - 1) / 2) as u8;
    let dsel: Option<u8> = args.get(2).and_then(|s| s.parse().ok());
    let mut grand = vec![W::ZERO; n + 1];
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
            grand[m] = grand[m].add(&c[m]);
        }
    }
    if dsel.is_none() {
        let row: Vec<String> = grand[1..].iter().map(|x| x.to_string()).collect();
        println!("regular(n) totals: [{}]", row.join(", "));
    }
}
