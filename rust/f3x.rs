// F_d' parent-stream generator by extension: reads a q-vertex digraph6
// grandparent stream, emits all (q+1)-vertex games with nullity exactly d and
// a nonnegative kernel vector (the parent family of the nullity-(d+1) stratum
// engine inc_hi.rs). directg cannot reach 9 vertices, so the n=10 parent
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
// With the third arg `u` (unique mode) each class is emitted EXACTLY ONCE and
// no labelg/sort pass is needed: a construction (G, c) is accepted iff the
// new vertex is the canonical choice within supp(kernel(child)) -- the
// sig-argmax, with rigidity fast path and canon fallback (canonical-max orbit
// on ties, per-grandparent dedup of Aut(G)-equivalent columns). The deleted
// argmax vertex determines the grandparent class, so acceptance across
// different G cannot collide. Same architecture as inc10.rs's parent
// acceptance. Requires linking the nauty shim.
//   nauty-geng 6 | nauty-directg -o | f3x 6 3 u | wc -l  -> 31257
//
// Child membership tests (exact): nullity d automatic from c _|_ ker;
// nonneg kernel vector automatic when G itself has one, else a d-column
// cone test on the child kernel basis {(k_i, 0), (x0hat, L)}; paradoxical +
// connected on the child bitmasks.
use std::collections::HashSet;
use std::env;
use std::io::{self, Read, Write as IoWrite};
use std::os::raw::c_int;

mod common;
use common::{cone_has_nonneg, dfs_es, encode, kernel_basis_any, order_columns_l1, suffix_abs_sums, vertex_sigs, DFS_RCAP};

extern "C" {
    fn rps_canon(arc: *const u64, n: c_int, canong: *mut u64, lab: *mut c_int, orbits: *mut c_int);
}

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

fn main() {
    let q: usize = env::args().nth(1).and_then(|s| s.parse().ok()).expect("usage: f3x q d < grandparents");
    let d: usize = env::args().nth(2).and_then(|s| s.parse().ok()).expect("usage: f3x q d [u]");
    assert!(d == 3 || d == 5);
    let unique = env::args().nth(3).as_deref() == Some("u");
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
    let mut canon_seen: HashSet<[u64; 16]> = HashSet::new();

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
            let mut gsig = [0u64; 16];
            let mut g_rigid = true;
            let mut god = [0u8; 16];
            let mut gid = [0u8; 16];
            if unique {
                for i in 0..q {
                    god[i] = gb[i].count_ones() as u8;
                    let mut c = 0u8;
                    for u in 0..q {
                        if gb[u] & (1 << i) != 0 {
                            c += 1;
                        }
                    }
                    gid[i] = c;
                }
                vertex_sigs(&gb, q, &mut gsig);
                for i in 0..q {
                    for j in (i + 1)..q {
                        if gsig[i] == gsig[j] {
                            g_rigid = false;
                        }
                    }
                }
                canon_seen.clear();
            }
            // particular-solution machinery: T*M = A reduced, pivots
            let (a, t, piv) = rref_aug(&m, q);
            // scale L = lcm |pivot values|
            let mut lval: i128 = 1;
            for &(r, c) in &piv {
                let pv = a[r][c].abs();
                lval = lval / common::gcd_i128(lval, pv).max(1) * pv;
            }
            // equality rows: c . ker_i = 0 (shared bound-pruned DFS core;
            // the leaf receives ordered coordinates, mapped back before use)
            let mut rows = [[0i64; 16]; DFS_RCAP];
            for (i, kv) in ker.iter().enumerate() {
                rows[i][..q].copy_from_slice(&kv[..q]);
            }
            let mut ord = [0usize; 16];
            let mut prows = [[0i64; 16]; DFS_RCAP];
            order_columns_l1(&rows, gd, q, &mut ord, &mut prows);
            let mut asum = [[0i64; 17]; DFS_RCAP];
            suffix_abs_sums(&prows, gd, q, &mut asum);
            let mut s0 = [0i64; DFS_RCAP];
            let mut c0 = [0i32; 16];
            dfs_es(0, q, gd, gd, &mut s0, &mut c0, &prows, &asum, 0, 0, &mut |rv: &[i32; 16]| {
                let mut cvarr = [0i32; 16];
                for k in 0..q {
                    cvarr[ord[k]] = rv[k];
                }
                let cv = &cvarr;
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
                // particular solution x0hat: M x0hat = -L c (needed for the
                // nonneg test when G is not semi, and for supp in unique mode)
                let mut x0 = [0i128; 16];
                if !g_semi || unique {
                    for &(r, pc) in &piv {
                        // b = T*(-c) row r, scaled by L/A[r][pc]
                        let mut b = 0i128;
                        for j in 0..q {
                            b -= t[r][j] * cv[j] as i128;
                        }
                        x0[pc] = b * (lval / a[r][pc]);
                    }
                }
                if unique {
                    // supp(kernel(child)): v < q in supp iff some of
                    // (ker_i[v], x0[v]) nonzero; the new vertex always is
                    let mut supp: u16 = 1 << q;
                    for v in 0..q {
                        let mut nz = x0[v] != 0;
                        for kv in ker.iter() {
                            nz |= kv[v] != 0;
                        }
                        if nz {
                            supp |= 1 << v;
                        }
                    }
                    // stage 0: degree-only reject (degrees are the leading
                    // signature key): if any supp vertex lex-beats the new
                    // vertex on child (od,id), the new vertex cannot be the
                    // sig-argmax -- no signature or canon work needed
                    let nod = newrow.count_ones() as u8;
                    let nid = plus.count_ones() as u8;
                    let mut mm = supp & !(1 << q);
                    while mm != 0 {
                        let v = mm.trailing_zeros() as usize;
                        mm &= mm - 1;
                        let od = god[v] + ((plus >> v) & 1) as u8;
                        let id = gid[v] + ((newrow >> v) & 1) as u8;
                        if od > nod || (od == nod && id > nid) {
                            return;
                        }
                    }
                    // acceptance: new vertex is the canonical choice in supp
                    let mut csig = [0u64; 16];
                    vertex_sigs(&cb, p, &mut csig);
                    let mut maxsig = 0u64;
                    let mut mm = supp;
                    while mm != 0 {
                        let i = mm.trailing_zeros() as usize;
                        mm &= mm - 1;
                        if csig[i] > maxsig {
                            maxsig = csig[i];
                        }
                    }
                    if csig[q] != maxsig {
                        return;
                    }
                    let mut ties = 0u32;
                    let mut mm = supp;
                    while mm != 0 {
                        let i = mm.trailing_zeros() as usize;
                        mm &= mm - 1;
                        if csig[i] == maxsig {
                            ties += 1;
                        }
                    }
                    let mut c_rigid = true;
                    for i in 0..p {
                        for j in (i + 1)..p {
                            if csig[i] == csig[j] {
                                c_rigid = false;
                            }
                        }
                    }
                    if !(ties == 1 && g_rigid && c_rigid) {
                        // canon path: orbit-canonical tie choice + per-G dedup
                        let mut arc64 = [0u64; 16];
                        for i in 0..p {
                            arc64[i] = cb[i] as u64;
                        }
                        let mut canong = [0u64; 16];
                        let mut lab = [0i32; 16];
                        let mut orbits = [0i32; 16];
                        unsafe {
                            rps_canon(arc64.as_ptr(), p as c_int, canong.as_mut_ptr(), lab.as_mut_ptr(), orbits.as_mut_ptr());
                        }
                        if ties > 1 {
                            let mut pos = [0i32; 16];
                            for (qq, &vv) in lab.iter().enumerate().take(p) {
                                pos[vv as usize] = qq as i32;
                            }
                            let mut best = usize::MAX;
                            let mut mm = supp;
                            while mm != 0 {
                                let i = mm.trailing_zeros() as usize;
                                mm &= mm - 1;
                                if csig[i] == maxsig && (best == usize::MAX || pos[i] > pos[best]) {
                                    best = i;
                                }
                            }
                            if orbits[q] != orbits[best] {
                                return;
                            }
                        }
                        let mut key = [0u64; 16];
                        key[..p].copy_from_slice(&canong[..p]);
                        if !canon_seen.insert(key) {
                            return;
                        }
                    }
                }
                // nonneg kernel vector of the child -- in unique mode this
                // runs only on acceptance survivors (~200x fewer cone tests)
                if !g_semi {
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
