// F_d' parent-stream generator by extension: reads a q-vertex digraph6
// grandparent stream, emits all (q+1)-vertex games with nullity exactly d and
// a nonnegative kernel vector (the parent family of the nullity-(d+1) stratum
// engine inc4.rs). directg cannot reach 9 vertices, so the n=10 parent
// streams (q=8, d=3 and d=5) are built this way instead.
//
// Theory: deleting a vertex of a nullity-d game changes nullity to d-1 iff
// some kernel vector is nonzero there, and supp(kernel) is never empty, so
// EVERY F_d' class arises by adding a vertex to a q-vertex nullity-(d-1)
// grandparent G with column c satisfying c _|_ ker(G) (then the child kernel
// is exactly {(x0 t + k, t)} of dimension d, where M x0 = -c). Emission is
// NOT deduplicated here -- pipe through `nauty-labelg | sort -u`:
//
//   nauty-geng 6 | nauty-directg -o | f3x 6 3 | nauty-labelg | sort -u \
//     | wc -l   -> 31257 (the inc_strata f3-emit 3 stream at n=8)
//
// Child membership tests (exact): nullity d automatic from c _|_ ker;
// nonneg kernel vector automatic when G itself has one, else a d-column
// cone test on the child kernel basis {(k_i, 0), (x0hat, L)}; paradoxical +
// connected on the child bitmasks.
use std::env;
use std::io::{self, Read, Write as IoWrite};

mod common;
use common::{cone_has_nonneg, encode, kernel_basis_any};

// fraction-free full reduction of [M | I]: returns (A, T, pivots) with
// T*M = A in reduced form (pivot columns cleared elsewhere), rows gcd-reduced
#[allow(clippy::type_complexity)]
fn rref_aug(m: &[[i64; 16]; 16], q: usize) -> ([[i128; 16]; 16], [[i128; 16]; 16], Vec<(usize, usize)>) {
    let mut a = [[0i128; 16]; 16];
    let mut t = [[0i128; 16]; 16];
    for i in 0..q {
        for j in 0..q {
            a[i][j] = m[i][j] as i128;
        }
        t[i][i] = 1;
    }
    let mut piv: Vec<(usize, usize)> = Vec::new(); // (row, col)
    let mut row = 0usize;
    for col in 0..q {
        if row >= q {
            break;
        }
        let mut pr = usize::MAX;
        for r in row..q {
            if a[r][col] != 0 {
                pr = r;
                break;
            }
        }
        if pr == usize::MAX {
            continue;
        }
        a.swap(row, pr);
        t.swap(row, pr);
        for r in 0..q {
            if r != row && a[r][col] != 0 {
                let (num, den) = (a[r][col], a[row][col]);
                for cc in 0..q {
                    a[r][cc] = a[r][cc] * den - num * a[row][cc];
                    t[r][cc] = t[r][cc] * den - num * t[row][cc];
                }
                let mut g = 0i128;
                for cc in 0..q {
                    g = common::gcd_i128(g, a[r][cc].abs());
                    g = common::gcd_i128(g, t[r][cc].abs());
                }
                if g > 1 {
                    for cc in 0..q {
                        a[r][cc] /= g;
                        t[r][cc] /= g;
                    }
                }
            }
        }
        piv.push((row, col));
        row += 1;
    }
    (a, t, piv)
}

// constrained cone DFS: enumerate c in {-1,0,1}^q with c . k_i = 0 for each
// kernel basis vector (suffix-bound pruned)
#[allow(clippy::too_many_arguments)]
fn c_dfs(
    k: usize,
    q: usize,
    ne: usize,
    s: &mut [i64; 4],
    c: &mut [i32; 16],
    rows: &[[i64; 16]; 4],
    asum: &[[i64; 17]; 4],
    f: &mut dyn FnMut(&[i32; 16]),
) {
    for i in 0..ne {
        if s[i].abs() > asum[i][k] {
            return;
        }
    }
    if k == q {
        f(c);
        return;
    }
    for val in [0i32, -1, 1] {
        c[k] = val;
        if val != 0 {
            let fv = val as i64;
            for i in 0..ne {
                s[i] += fv * rows[i][k];
            }
        }
        c_dfs(k + 1, q, ne, s, c, rows, asum, f);
        if val != 0 {
            let fv = val as i64;
            for i in 0..ne {
                s[i] -= fv * rows[i][k];
            }
        }
    }
    c[k] = 0;
}

fn main() {
    let q: usize = env::args().nth(1).and_then(|s| s.parse().ok()).expect("usage: f3x q d < grandparents");
    let d: usize = env::args().nth(2).and_then(|s| s.parse().ok()).expect("usage: f3x q d");
    assert!(d == 3 || d == 5);
    let gd = d - 1; // grandparent nullity
    let qsq = q * q;
    let qbytes = (qsq + 5) / 6;
    let reclen = 2 + qbytes + 1;
    let p = q + 1; // emitted parent size

    let mut stdin = io::stdin().lock();
    let stdout = io::stdout().lock();
    let mut out = io::BufWriter::new(stdout);
    let mut buf = vec![0u8; 1 << 20];
    let mut obuf: Vec<u8> = Vec::with_capacity(64);
    let mut have = 0usize;
    let (mut gseen, mut gmatch, mut emitted) = (0u64, 0u64, 0u64);

    loop {
        let got = stdin.read(&mut buf[have..]).unwrap();
        if got == 0 {
            break;
        }
        have += got;
        let nrec = have / reclen;
        for ri in 0..nrec {
            let rec = &buf[ri * reclen..(ri + 1) * reclen];
            assert!(rec[0] == b'&' && rec[1] as usize == 63 + q, "misaligned digraph6");
            gseen += 1;
            let mut gb = [0u16; 16];
            let payload = &rec[2..2 + qbytes];
            let mut kk = 0usize;
            'dec: for &byte in payload {
                let mut bits = ((byte - 63) as u32) << 26;
                for _ in 0..6 {
                    if bits & 0x8000_0000 != 0 {
                        gb[kk / q] |= 1 << (kk % q);
                    }
                    bits <<= 1;
                    kk += 1;
                    if kk == qsq {
                        break 'dec;
                    }
                }
            }
            let mut m = [[0i64; 16]; 16];
            for i in 0..q {
                let mut w = gb[i];
                while w != 0 {
                    let j = w.trailing_zeros() as usize;
                    w &= w - 1;
                    m[i][j] = 1;
                    m[j][i] = -1;
                }
            }
            let ker = kernel_basis_any(&m, q);
            if ker.len() != gd {
                continue;
            }
            gmatch += 1;
            // does G itself have a nonneg kernel vector? then every child does
            let g_semi = cone_has_nonneg(&ker, q, gd);
            // particular-solution machinery: T*M = A reduced, pivots
            let (a, t, piv) = rref_aug(&m, q);
            // scale L = lcm |pivot values|
            let mut lval: i128 = 1;
            for &(r, c) in &piv {
                let pv = a[r][c].abs();
                lval = lval / common::gcd_i128(lval, pv).max(1) * pv;
            }
            // equality rows: c . ker_i = 0
            let mut rows = [[0i64; 16]; 4];
            for (i, kv) in ker.iter().enumerate() {
                rows[i][..q].copy_from_slice(&kv[..q]);
            }
            let mut asum = [[0i64; 17]; 4];
            for i in 0..gd {
                for k in (0..q).rev() {
                    asum[i][k] = asum[i][k + 1] + rows[i][k].abs();
                }
            }
            let mut s0 = [0i64; 4];
            let mut c0 = [0i32; 16];
            c_dfs(0, q, gd, &mut s0, &mut c0, &rows, &asum, &mut |cv: &[i32; 16]| {
                // child bitmasks: cv[i] > 0 => i beats new; cv[i] < 0 => new beats i
                let mut cb = gb;
                let mut newrow = 0u16;
                let mut plus = 0u16;
                for i in 0..q {
                    if cv[i] > 0 {
                        cb[i] |= 1 << q;
                        plus |= 1 << i;
                    } else if cv[i] < 0 {
                        newrow |= 1 << i;
                    }
                }
                cb[q] = newrow;
                let _ = plus;
                // NOTE: the family is nullity-d + nonneg kernel ONLY -- parents
                // need not be paradoxical or connected (the stratum engine's
                // leaves test the CHILD for those)
                // nonneg kernel vector of the child
                if !g_semi {
                    // particular solution x0hat: M x0hat = -L c
                    let mut x0 = [0i128; 16];
                    for &(r, pc) in &piv {
                        // b = T*(-c) row r, scaled by L/A[r][pc]
                        let mut b = 0i128;
                        for j in 0..q {
                            b -= t[r][j] * cv[j] as i128;
                        }
                        x0[pc] = b * (lval / a[r][pc]);
                    }
                    // child kernel basis: (ker_i, 0) and (x0hat, L)
                    let mut basis: Vec<[i64; 16]> = Vec::with_capacity(d);
                    for kv in ker.iter() {
                        let mut v = [0i64; 16];
                        v[..q].copy_from_slice(&kv[..q]);
                        basis.push(v);
                    }
                    let mut v = [0i64; 16];
                    let mut mx: i128 = lval.abs();
                    for x in x0.iter().take(q) {
                        mx = mx.max(x.abs());
                    }
                    // gcd-reduce to keep i64
                    let mut g = lval.abs();
                    for x in x0.iter().take(q) {
                        g = common::gcd_i128(g, x.abs());
                    }
                    let g = g.max(1);
                    assert!(mx / g < (1i128 << 62), "x0 overflow");
                    for i in 0..q {
                        v[i] = (x0[i] / g) as i64;
                    }
                    v[q] = (lval / g) as i64;
                    basis.push(v);
                    if !cone_has_nonneg(&basis, p, d) {
                        return;
                    }
                }
                emitted += 1;
                encode(&cb, p, &mut obuf);
                out.write_all(&obuf).unwrap();
            });
        }
        let rem = have - nrec * reclen;
        buf.copy_within(nrec * reclen..have, 0);
        have = rem;
    }
    out.flush().unwrap();
    eprintln!("f3x: q={} d={} grandparents={} nullity{}={} emitted={}", q, d, gseen, gd, gmatch, emitted);
}
