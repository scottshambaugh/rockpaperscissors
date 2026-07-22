// Fix_sigma for cycle types (2^k, 1^f) by a k-marked-tuple class sweep --
// replaces the 3^(2*C(k,2) + ...) raw brute of sigma_fix, which is infeasible
// at n=10 for k=3 (3^24) and slow for k=4 (3^21).
//
// The quotient of a (2^k,1^f)-invariant game on n = 2k+f vertices has
// m = k+f supervertices: k doubles (the 2-cycles) and f singles. Entries:
//   * double-double (i<j): a g=2 bundle, TWO trits (t1, t2) -- class 0 links
//     a_i-a_j, b_i-b_j; class 1 links a_i-b_j, b_i-a_j.
//   * everything else: one plain trit.
// Convention: a structure corresponds to exactly one labeled plain m-game
// with ZERO in the k*(k-1)/2 double-double slots, plus the raw bundle trits.
//
// Conditions per raw assignment:
//   * balance: exists w>0 with M~ w = 0, where 2*M~ has the double-double
//     entries s_ij = t1+t2 and all other entries doubled. Depends on the
//     s-vector only -> memoized per (class, tuple).
//   * lifted paradox: a double's two members see the SAME sign multiset:
//     plain row of the double + the incident bundle trits (negated on the
//     j-side). Singles: their plain row (arcs to doubles duplicate).
//   * lifted connectivity: explicit BFS on the (2k+f)-vertex lifted graph
//     (bundle class edges as above; single-double arcs edge to both members).
//
// Labeled counting via classes: NUM = sum over classes (m!/|Aut|) * #(ordered
// pairwise-tie k-tuples passing), Fix = NUM / (m falling-factorial k).
//
// Anchors (all against sigma_fix brute):
//   k=3: geng 3 -> (2,2,2) n=6:      fix = 66
//        geng 5 -> (2,2,2,1,1) n=8:  fix = 32310
//   k=4: geng 4 -> (2,2,2,2) n=8:    fix = 20298
// Production: k=3 over geng 7 (2.1M classes) -> Fix_(2,2,2,1^4);
//             k=4 over geng 6 -> Fix_(2,2,2,2,1^2).
//
//   rustc -O rust/sigma_ktuple.rs -o /tmp/ktup -C link-args="shim.o -lnauty"
//   nauty-geng 5 | nauty-directg -o | /tmp/ktup 5 3
use std::env;
use std::io::{self, Read};
use std::os::raw::c_int;

mod common;
use common::{factorial, has_positive_kernel};

extern "C" {
    fn rps_autsize(arc: *const u64, n: c_int) -> f64;
}

const MAXB: usize = 6; // C(k,2) for k <= 4

// patterns (t1,t2) with t1+t2 = s, for s+2 in 0..5
const PATS: [&[(i8, i8)]; 5] = [
    &[(-1, -1)],
    &[(-1, 0), (0, -1)],
    &[(0, 0), (1, -1), (-1, 1)],
    &[(1, 0), (0, 1)],
    &[(1, 1)],
];

struct Pre {
    nb: usize,                   // number of bundles = C(k,2)
    bpairs: Vec<(usize, usize)>, // bundle -> (tuple positions i<j)
}

fn build_pre(k: usize) -> Pre {
    let mut bpairs = Vec::new();
    for i in 0..k {
        for j in (i + 1)..k {
            bpairs.push((i, j));
        }
    }
    let nb = bpairs.len();
    Pre { nb, bpairs }
}

fn main() {
    let m: usize = env::args().nth(1).and_then(|s| s.parse().ok()).expect("usage: ktup m k < classes.d6");
    let k: usize = env::args().nth(2).and_then(|s| s.parse().ok()).expect("usage: ktup m k");
    assert!((3..=4).contains(&k) && k <= m && m <= 8);
    let f = m - k;
    let pre = build_pre(k);
    let nlift = 2 * k + f;
    let msq = m * m;
    let mbytes = (msq + 5) / 6;
    let reclen = 2 + mbytes + 1;
    let mfact = factorial(m as u64);
    let kfact = (1..=k as u64).product::<u64>();

    let mut stdin = io::stdin().lock();
    let mut buf = vec![0u8; 1 << 20];
    let mut have = 0usize;
    let (mut classes, mut tuples_tested) = (0u64, 0u64);
    let mut num: u128 = 0;

    loop {
        let got = stdin.read(&mut buf[have..]).unwrap();
        if got == 0 {
            break;
        }
        have += got;
        let nrec = have / reclen;
        for ri in 0..nrec {
            let rec = &buf[ri * reclen..(ri + 1) * reclen];
            assert!(rec[0] == b'&' && rec[1] as usize == 63 + m, "misaligned digraph6");
            classes += 1;
            let mut beats = [0u16; 16];
            let payload = &rec[2..2 + mbytes];
            let mut kk = 0usize;
            'dec: for &byte in payload {
                let mut bits = ((byte - 63) as u32) << 26;
                for _ in 0..6 {
                    if bits & 0x8000_0000 != 0 {
                        beats[kk / m] |= 1 << (kk % m);
                    }
                    bits <<= 1;
                    kk += 1;
                    if kk == msq {
                        break 'dec;
                    }
                }
            }
            let mut inn = [0u16; 16];
            for i in 0..m {
                let mut w = beats[i];
                while w != 0 {
                    let j = w.trailing_zeros() as usize;
                    w &= w - 1;
                    inn[j] |= 1 << i;
                }
            }
            // plain row sign presence per vertex
            let mut rowsg = [0u8; 16];
            for i in 0..m {
                rowsg[i] = ((beats[i] != 0) as u8) | (((inn[i] != 0) as u8) << 1);
            }
            // vertices whose plain row lacks a sign: must be marked (their
            // bundles can repair); if more than k such, no tuple qualifies
            let mut defect = 0u16;
            for i in 0..m {
                if rowsg[i] != 3 {
                    defect |= 1 << i;
                }
            }
            if (defect.count_ones() as usize) > k {
                continue;
            }
            // enumerate k-SETS (ascending tuples): the per-ordering count is
            // invariant under permuting the tuple (bundle/class relabeling is
            // a bijection of the assignment space), so ordered = k! * sets
            let mut tup = [0usize; 4];
            let mut count_class = 0u64;
            enum_tuples(m, k, &beats, &inn, defect, &mut tup, 0, 0, &mut |tupv: &[usize; 4]| {
                tuples_tested += 1;
                count_class += tuple_count(
                    m, k, nlift, tupv, &beats, &rowsg, &pre,
                );
            });
            if count_class > 0 {
                let mut arc64 = [0u64; 16];
                for i in 0..m {
                    arc64[i] = beats[i] as u64;
                }
                let aut = common::autsize_u128(unsafe { rps_autsize(arc64.as_ptr(), m as c_int) });
                num += (mfact / aut) * (count_class * kfact) as u128;
            }
        }
        let rem = have - nrec * reclen;
        buf.copy_within(nrec * reclen..have, 0);
        have = rem;
    }
    let mut ff: u128 = 1;
    for t in 0..k {
        ff *= (m - t) as u128;
    }
    println!(
        "m={} k={} classes={} tuples={} rawnum={}",
        m, k, classes, tuples_tested, num
    );
    if num % ff == 0 {
        println!("Fix_(2^{},1^{}) = {}", k, f, num / ff);
    } else {
        println!("NOT INTEGRAL by {} -- full class stream required", ff);
    }
}

// ordered k-tuples of distinct vertices, pairwise tied, covering `defect`
fn enum_tuples(
    m: usize,
    k: usize,
    beats: &[u16; 16],
    inn: &[u16; 16],
    defect: u16,
    tup: &mut [usize; 4],
    pos: usize,
    used: u16,
    f: &mut dyn FnMut(&[usize; 4]),
) {
    if pos == k {
        if defect & !used == 0 {
            f(tup);
        }
        return;
    }
    let start = if pos == 0 { 0 } else { tup[pos - 1] + 1 };
    for v in start..m {
        if used & (1 << v) != 0 {
            continue;
        }
        // pairwise tie with all previous tuple members
        let mut ok = true;
        for &p in tup.iter().take(pos) {
            if beats[p] & (1 << v) != 0 || inn[p] & (1 << v) != 0 {
                ok = false;
                break;
            }
        }
        if !ok {
            continue;
        }
        tup[pos] = v;
        enum_tuples(m, k, beats, inn, defect, tup, pos + 1, used | (1 << v), f);
    }
}

#[allow(clippy::too_many_arguments)]
fn tuple_count(
    m: usize,
    k: usize,
    nlift: usize,
    tup: &[usize; 4],
    beats: &[u16; 16],
    rowsg: &[u8; 16],
    pre: &Pre,
) -> u64 {
    // lifted vertex layout: doubles: tuple position i -> lifted 2i, 2i+1;
    // singles: the m-k non-tuple vertices in order -> lifted 2k..
    let mut is_tup = [usize::MAX; 16];
    for (i, &v) in tup.iter().enumerate().take(k) {
        is_tup[v] = i;
    }
    let mut single_of = [usize::MAX; 16];
    let mut si = 2 * k;
    for v in 0..m {
        if is_tup[v] == usize::MAX {
            single_of[v] = si;
            si += 1;
        }
    }
    // base lifted adjacency (undirected), without bundle edges
    let mut inn = [0u16; 16];
    for i in 0..m {
        let mut w = beats[i];
        while w != 0 {
            let j = w.trailing_zeros() as usize;
            w &= w - 1;
            inn[j] |= 1 << i;
        }
    }
    let mut base = [0u16; 12];
    let mut edge = |a: usize, b: usize, base: &mut [u16; 12]| {
        base[a] |= 1 << b;
        base[b] |= 1 << a;
    };
    for u in 0..m {
        for v in (u + 1)..m {
            if beats[u] & (1 << v) == 0 && beats[v] & (1 << u) == 0 {
                continue; // tie: no edge
            }
            match (is_tup[u], is_tup[v]) {
                (usize::MAX, usize::MAX) => edge(single_of[u], single_of[v], &mut base),
                (i, usize::MAX) if i != usize::MAX => {
                    edge(2 * i, single_of[v], &mut base);
                    edge(2 * i + 1, single_of[v], &mut base);
                }
                (usize::MAX, j) => {
                    edge(2 * j, single_of[u], &mut base);
                    edge(2 * j + 1, single_of[u], &mut base);
                }
                _ => {}
            }
        }
    }
    let mut psign = [0u8; 4];
    for i in 0..k {
        psign[i] = rowsg[tup[i]];
    }
    // scaled quotient matrix base (2*rel); double-double slots filled per s
    let mut m2 = [[0i64; 16]; 16];
    for u in 0..m {
        let mut w = beats[u];
        while w != 0 {
            let v = w.trailing_zeros() as usize;
            w &= w - 1;
            m2[u][v] = 2;
            m2[v][u] = -2;
        }
    }
    let full_mask: u16 = ((1u32 << m) - 1) as u16;
    let nkeys = 5usize.pow(pre.nb as u32);
    let mut count = 0u64;
    let mut svals = [0i64; MAXB];
    for key in 0..nkeys {
        let mut kk = key;
        for b in 0..pre.nb {
            let sv = (kk % 5) as i64 - 2;
            kk /= 5;
            svals[b] = sv;
            let (i, j) = pre.bpairs[b];
            m2[tup[i]][tup[j]] = sv;
            m2[tup[j]][tup[i]] = -sv;
        }
        // balance: even m needs singularity first (Pfaffian gate); odd m is
        // always singular with the nullity-1 kernel vector given by the
        // Pfaffian cofactors v_i = (-1)^i Pf(M with row/col i deleted) --
        // strictly positive kernel iff all cofactors nonzero with alternating
        // signs matching. RREF+Gordan only for the rare nullity>=3 case.
        let bal = if m % 2 == 0 {
            common::pf(&m2, full_mask) == 0 && has_positive_kernel(&m2, m)
        } else {
            let mut allzero = true;
            let mut pos = 0usize;
            let mut neg = 0usize;
            for i in 0..m {
                let c = common::pf(&m2, full_mask & !(1 << i));
                let v = if i % 2 == 0 { c } else { -c };
                if v != 0 {
                    allzero = false;
                }
                if v > 0 {
                    pos += 1;
                } else if v < 0 {
                    neg += 1;
                }
            }
            if allzero {
                has_positive_kernel(&m2, m)
            } else {
                pos == m || neg == m
            }
        };
        if !bal {
            continue;
        }
        // enumerate the bundle patterns compatible with this s-vector:
        // paradox (per-double sign multisets) + lifted connectivity
        let mut pidx = [0usize; MAXB];
        'combo: loop {
            // evaluate this pattern combo
            let mut sg = psign;
            let mut adj = base;
            for b in 0..pre.nb {
                let (t1, t2) = PATS[(svals[b] + 2) as usize][pidx[b]];
                let (i, j) = pre.bpairs[b];
                for t in [t1, t2] {
                    if t > 0 {
                        sg[i] |= 1;
                        sg[j] |= 2;
                    } else if t < 0 {
                        sg[i] |= 2;
                        sg[j] |= 1;
                    }
                }
                if t1 != 0 {
                    adj[2 * i] |= 1 << (2 * j);
                    adj[2 * j] |= 1 << (2 * i);
                    adj[2 * i + 1] |= 1 << (2 * j + 1);
                    adj[2 * j + 1] |= 1 << (2 * i + 1);
                }
                if t2 != 0 {
                    adj[2 * i] |= 1 << (2 * j + 1);
                    adj[2 * j + 1] |= 1 << (2 * i);
                    adj[2 * i + 1] |= 1 << (2 * j);
                    adj[2 * j] |= 1 << (2 * i + 1);
                }
            }
            let mut pok = true;
            for &x in sg.iter().take(k) {
                if x != 3 {
                    pok = false;
                    break;
                }
            }
            if pok {
                let full: u16 = ((1u32 << nlift) - 1) as u16;
                let mut seen: u16 = 1;
                let mut fr: u16 = 1;
                while fr != 0 {
                    let mut nf = 0u16;
                    let mut w = fr;
                    while w != 0 {
                        let v = w.trailing_zeros() as usize;
                        w &= w - 1;
                        nf |= adj[v] & !seen;
                    }
                    seen |= nf;
                    fr = nf;
                }
                if seen == full {
                    count += 1;
                }
            }
            // next combo
            let mut b = 0;
            loop {
                if b == pre.nb {
                    break 'combo;
                }
                pidx[b] += 1;
                if pidx[b] < PATS[(svals[b] + 2) as usize].len() {
                    break;
                }
                pidx[b] = 0;
                b += 1;
            }
        }
        // restore tie slots
        for b in 0..pre.nb {
            let (i, j) = pre.bpairs[b];
            m2[tup[i]][tup[j]] = 0;
            m2[tup[j]][tup[i]] = 0;
        }
    }
    count
}
