// inclusive(even n) WITHOUT enumerating children: the deletion-multiplicity
// identity. Every labeled inclusive n-game C of nullity 2 has equilibrium-
// segment endpoints p, q (disjoint zero sets, covering all coordinates), and
// its valid deletions -- vertices whose removal leaves a nullity-1 game with a
// NONNEGATIVE kernel ("semi-CM" parent) -- are EXACTLY the endpoint zeros:
// z(C) = |Z_p| + |Z_q|. (Deleting a doubly-positive vertex leaves the mixed-
// sign slice generator q_v p - p_v q: not semi-CM.) Therefore
//
//   L_2 = #labeled nullity-2 inclusive n-games
//       = n * sum over labeled semi-CM parents P, valid r:  1 / z(child(P,r))
//
// with no isomorphism dedup at all: each labeled child is built z times and
// contributes z * (1/z) = 1. Parents stream in as digraph6 classes; each class
// contributes (n-1)!/|Aut(P)| labeled copies (|Aut| via nauty).
//
// Per parent (exact integer arithmetic throughout):
//   * kernel vector v_i = (-1)^i Pf(M'_-i) (integer, nonneg after sign flip);
//   * bordered matrix B = [[M', v],[v^T, 0]] is nonsingular; D = det B and the
//     adjugate block give w = -(ADJ r)/D with M'w = -r, so every per-child
//     quantity is linear in r;
//   * child inclusive <=> w_i > 0 for all i in Z(v)  (then lambda*v + w > 0
//     for large lambda) -- strict integer inequalities, linear in r;
//   * a bound-pruned DFS over r in {-1,0,1}^(n-1) walks only candidates with
//     v.r = 0 (nullity-2 condition) AND the Z(v) positivity rows;
//   * per leaf: child paradox+connected (bitmask), endpoint 2 = lambda* v + w
//     with lambda* = max_i(-w_i/v_i) over supp(v) (exact fraction comparisons),
//     z = (|Z(v)| + 1) + |argmax|, accumulate (n-1)!/|Aut| * LCM/z.
//
// Output: SUM (integer, scaled by LCM = lcm(1..2n)); L_2 = n * SUM / LCM.
// The n=8 instantiation must reproduce the directg-8 ground truth
// (labeled nullity-2 inclusive count) exactly before n=10 is trusted.
//
//   rustc -O rust/inc_count.rs -o /tmp/incc -C link-args="shim.o -lnauty"
//   nauty-geng 7 | nauty-directg -o | /tmp/incs 7 semi-emit | /tmp/incc 8
use std::env;
use std::io::{self, Read};
use std::os::raw::c_int;

mod common;
use common::{adjugate_ff, gcd_i128, paradox_connected_beats};

extern "C" {
    fn rps_autsize(arc: *const u64, n: c_int) -> f64;
}

fn factorial(n: usize) -> u128 {
    common::factorial(n as u64)
}

fn lcm_to(n: u64) -> u64 {
    fn gcd(a: u64, b: u64) -> u64 {
        if b == 0 { a } else { gcd(b, a % b) }
    }
    (1..=n).fold(1u64, |l, x| l / gcd(l, x) * x)
}

// Fraction-free (Bareiss) Gauss-Jordan on [B | I]: returns (adj(B), det B)
// with adj(B) = det * B^-1, all integer, exact divisions only -- no gcds.
fn adjugate(b0: &[[i128; 11]; 11], m: usize) -> ([[i128; 11]; 11], i128) {
    adjugate_ff(b0, m, true).expect("bordered singular")
}


// DFS over r in {-1,0,1}^p: eq row 0 (exact 0), ineq rows ne..nt (< 0),
// tracking-only rows nt..nl (endpoint data, no pruning). Leaves push (r, s).
#[allow(clippy::too_many_arguments)]
fn dfs_r2(
    k: usize,
    p: usize,
    ne: usize,
    nt: usize,
    nl: usize,
    s: &mut [i128; 14],
    r: &mut [i32; 10],
    rows: &[[i128; 10]; 14],
    asum: &[[i128; 11]; 14],
    out: &mut Vec<([i32; 10], [i128; 14])>,
) {
    for i in 0..ne {
        if s[i].abs() > asum[i][k] {
            return;
        }
    }
    for i in ne..nt {
        if s[i] - asum[i][k] >= 0 {
            return;
        }
    }
    if k == p {
        out.push((*r, *s));
        return;
    }
    for val in [0i32, -1, 1] {
        r[k] = val;
        if val != 0 {
            let f = val as i128;
            for i in 0..nl {
                s[i] += f * rows[i][k];
            }
        }
        dfs_r2(k + 1, p, ne, nt, nl, s, r, rows, asum, out);
        if val != 0 {
            let f = val as i128;
            for i in 0..nl {
                s[i] -= f * rows[i][k];
            }
        }
    }
    r[k] = 0;
}

fn main() {
    let n: usize = env::args().nth(1).and_then(|s| s.parse().ok()).expect("usage: incc n < parents");
    assert!(n % 2 == 0 && (6..=10).contains(&n));
    let p = n - 1; // parent size (odd)
    let psq = p * p;
    let pbytes = (psq + 5) / 6;
    let reclen = 2 + pbytes + 1;
    let full: u16 = (1u16 << p) - 1;
    let lcm = lcm_to(2 * n as u64) as u128;
    let pfact = factorial(p);

    let mut stdin = io::stdin().lock();
    let mut buf = vec![0u8; 1 << 20];
    let mut have = 0usize;
    let (mut parents, mut skipped, mut leaves) = (0u64, 0u64, 0u64);
    let mut sum: u128 = 0; // sum of (p!/|Aut|) * lcm/z
    let mut valid: Vec<([i32; 10], [i128; 14])> = Vec::with_capacity(64);

    loop {
        let got = stdin.read(&mut buf[have..]).unwrap();
        if got == 0 {
            break;
        }
        have += got;
        let nrec = have / reclen;
        for ri in 0..nrec {
            let rec = &buf[ri * reclen..(ri + 1) * reclen];
            assert!(rec[0] == b'&' && rec[1] as usize == 63 + p, "misaligned digraph6");
            parents += 1;
            // decode parent
            let mut beats = [0u16; 16];
            let payload = &rec[2..2 + pbytes];
            let mut kk = 0usize;
            'dec: for &byte in payload {
                let mut bits = ((byte - 63) as u32) << 26;
                for _ in 0..6 {
                    if bits & 0x8000_0000 != 0 {
                        beats[kk / p] |= 1 << (kk % p);
                    }
                    bits <<= 1;
                    kk += 1;
                    if kk == psq {
                        break 'dec;
                    }
                }
            }
            let mut m = [[0i64; 16]; 16];
            for i in 0..p {
                let mut w = beats[i];
                while w != 0 {
                    let j = w.trailing_zeros() as usize;
                    w &= w - 1;
                    m[i][j] = 1;
                    m[j][i] = -1;
                }
            }
            // integer kernel vector via Pfaffian cofactors, sign-normalized
            let mut v = [0i128; 10];
            let mut pos = false;
            let mut neg = false;
            let mut nonzero = 0;
            for i in 0..p {
                let pfv = common::pf(&m, full & !(1u16 << i));
                let vi = if i % 2 == 0 { pfv } else { -pfv } as i128;
                v[i] = vi;
                if vi > 0 {
                    pos = true;
                } else if vi < 0 {
                    neg = true;
                }
                if vi != 0 {
                    nonzero += 1;
                }
            }
            if nonzero == 0 || (pos && neg) {
                skipped += 1; // not a semi-CM parent
                continue;
            }
            if neg {
                for vi in v.iter_mut() {
                    *vi = -*vi;
                }
            }
            // |Aut(P)|
            let mut arc64 = [0u64; 16];
            for i in 0..p {
                arc64[i] = beats[i] as u64;
            }
            let aut = unsafe { rps_autsize(arc64.as_ptr(), p as c_int) } as u128;
            let wp = pfact / aut;
            // bordered matrix + adjugate block rows for i in 0..p
            let mut b0 = [[0i128; 11]; 11];
            for i in 0..p {
                for j in 0..p {
                    b0[i][j] = m[i][j] as i128;
                }
                b0[i][p] = v[i];
                b0[p][i] = v[i];
            }
            // ADJ columns: adjcol[j] = x with B x = det * e_j; we need rows of
            // the top-left block: w = -(1/D) * (ADJ r) with ADJ[i][j] = adj(B)[i][j]
            // = (column-j solve)[i]. Build ADJ[i][j] for i,j < p.
            let (adjm, det) = adjugate(&b0, p + 1);
            let mut adjb = [[0i128; 10]; 10];
            for i in 0..p {
                for j in 0..p {
                    adjb[i][j] = adjm[i][j];
                }
            }
            let sgn: i128 = if det > 0 { 1 } else { -1 };
            // DFS rows: eq = v; ineq rows for i in Z(v): sgn * (ADJ[i] . r) < 0
            let mut rows = [[0i128; 10]; 14];
            let mut ne = 1usize;
            for j in 0..p {
                rows[0][j] = v[j];
            }
            let mut zrows = ne;
            let mut zset = [false; 10];
            for i in 0..p {
                if v[i] == 0 {
                    zset[i] = true;
                    for j in 0..p {
                        rows[zrows][j] = sgn * adjb[i][j];
                    }
                    zrows += 1;
                }
            }
            let nt = zrows;
            ne = 1;
            let mut asum = [[0i128; 11]; 14];
            for i in 0..nt {
                for k in (0..p).rev() {
                    asum[i][k] = asum[i][k + 1] + rows[i][k].abs();
                }
            }
            // supp rows: track a_i = ADJ[i].r incrementally (no pruning) so
            // leaves read endpoint data without dot products
            let mut nsupp = 0usize;
            let mut suppv = [0usize; 10];
            let mut all_rows = rows;
            for i in 0..p {
                if v[i] != 0 {
                    for j in 0..p {
                        all_rows[nt + nsupp][j] = adjb[i][j];
                    }
                    suppv[nsupp] = i;
                    nsupp += 1;
                }
            }
            valid.clear();
            let mut s0 = [0i128; 14];
            let mut r0 = [0i32; 10];
            dfs_r2(0, p, ne, nt, nt + nsupp, &mut s0, &mut r0, &all_rows, &asum, &mut valid);
            if valid.is_empty() {
                continue;
            }
            for (rv, sv) in valid.iter() {
                // build child bitmasks: r_i = M[i][new]: +1 => i beats new
                let mut cb = [0u16; 16];
                for i in 0..p {
                    cb[i] = beats[i];
                    if rv[i] > 0 {
                        cb[i] |= 1 << p;
                    }
                }
                let mut newrow = 0u16;
                for i in 0..p {
                    if rv[i] < 0 {
                        newrow |= 1 << i;
                    }
                }
                cb[p] = newrow;
                if !paradox_connected_beats(&cb, n) {
                    continue;
                }
                // endpoint 2: lambda* = max over supp(v) of a_i/(D v_i) with
                // a_i = ADJ[i].r  (w_i = -a_i/D; lambda* = max(-w_i/v_i))
                let mut best_n: i128 = 0;
                let mut best_d: i128 = 0; // fraction best_n/best_d, d>0
                let mut argmax = 0u32;
                let mut first = true;
                for t in 0..nsupp {
                    let i = suppv[t];
                    let ai: i128 = sv[nt + t];
                    // fraction f_i = ai / (det * v[i]); normalize denom > 0
                    let (mut fn_, mut fd) = (ai, det * v[i]);
                    if fd < 0 {
                        fn_ = -fn_;
                        fd = -fd;
                    }
                    if first {
                        best_n = fn_;
                        best_d = fd;
                        argmax = 1;
                        first = false;
                    } else {
                        // compare fn_/fd vs best_n/best_d
                        let lhs = fn_ * best_d;
                        let rhs = best_n * fd;
                        if lhs > rhs {
                            best_n = fn_;
                            best_d = fd;
                            argmax = 1;
                        } else if lhs == rhs {
                            argmax += 1;
                        }
                    }
                }
                let z1 = v[..p].iter().filter(|&&x| x == 0).count() as u128 + 1;
                let z = z1 + argmax as u128;
                leaves += 1;
                sum += wp * (lcm / z);
            }
        }
        let rem = have - nrec * reclen;
        buf.copy_within(nrec * reclen..have, 0);
        have = rem;
    }
    let total = (n as u128) * sum;
    assert!(total % lcm == 0, "1/z weights did not resolve to an integer");
    println!(
        "n={}: parents={} (skipped {}) leaves={} L_nullity2_labeled={}",
        n,
        parents,
        skipped,
        leaves,
        total / lcm
    );
}

