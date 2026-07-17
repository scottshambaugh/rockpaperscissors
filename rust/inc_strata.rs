// Ground-truth classifier for the inclusive(even n) counting project.
// Reads a digraph6 stream of n-vertex oriented games (nauty-directg -o: one rep
// per iso class) and reports, for the paradoxical+connected+fully-mixed
// (= inclusive) games, the breakdown by kernel nullity, plus the LABELED count
// (sum of n!/|Aut|, |Aut| via the nauty shim's orbit... here via canonical
// refinement: we use the automorphism count from repeated canon calls is
// expensive, so |Aut| is computed by brute orbit-stabilizer over canon labels
// -- for validation-scale streams only).
//
// Also counts, for ODD n streams, the semi-CM games (nullity 1 with NONNEGATIVE
// kernel: all Pfaffian cofactors zero-or-one-signed, not all zero) -- the
// parent family of the extension-counting identity.
//
//   rustc -O rust/inc_strata.rs -o /tmp/incs
//   nauty-geng 8 | nauty-directg -o | /tmp/incs 8      # even: strata
//   nauty-geng 7 | nauty-directg -o | /tmp/incs 7 semi # odd: semi-CM count
use std::env;
use std::io::{self, Read};

mod common;
use common::{cone_has_nonneg, fully_mixed, kernel_basis_exact, pf};
use std::os::raw::c_int;
extern "C" {
    fn rps_autsize(arc: *const u64, n: c_int) -> f64;
}

fn main() {
    let n: usize = env::args().nth(1).and_then(|s| s.parse().ok()).expect("usage: incs n [semi]");
    let semi_mode = env::args().nth(2).as_deref().map(|s| s.starts_with("semi")).unwrap_or(false);
    let emit = env::args().nth(2).as_deref() == Some("semi-emit");
    // f3-emit D: emit odd-n games with nullity exactly D and a nonneg kernel vec
    let f3_mode = env::args().nth(2).as_deref() == Some("f3-emit");
    let f3_d: usize = env::args().nth(3).and_then(|s| s.parse().ok()).unwrap_or(3);
    let nsq = n * n;
    let nbytes = (nsq + 5) / 6;
    let reclen = 2 + nbytes + 1;
    let mut stdin = io::stdin().lock();
    let mut buf = vec![0u8; 1 << 20];
    let mut have = 0usize;
    let full: u16 = (1u16 << n) - 1;
    let (mut total, mut semi, mut semi_strict) = (0u64, 0u64, 0u64);
    // inclusive counts by nullity (index = nullity): classes and labeled sums
    let mut incl = [0u64; 12];
    let mut lab_sum = [0u128; 12];
    let nfact = common::factorial(n as u64);
    loop {
        let got = stdin.read(&mut buf[have..]).unwrap();
        if got == 0 {
            break;
        }
        have += got;
        let nrec = have / reclen;
        for r in 0..nrec {
            let rec = &buf[r * reclen..(r + 1) * reclen];
            assert!(rec[0] == b'&' && rec[1] as usize == 63 + n, "misaligned digraph6");
            total += 1;
            let mut out = [0u16; 16];
            let mut inn = [0u16; 16];
            let payload = &rec[2..2 + nbytes];
            let mut k = 0usize;
            'dec: for &byte in payload {
                let mut bits = ((byte - 63) as u32) << 26;
                for _ in 0..6 {
                    if bits & 0x8000_0000 != 0 {
                        let (i, j) = (k / n, k % n);
                        out[i] |= 1 << j;
                        inn[j] |= 1 << i;
                    }
                    bits <<= 1;
                    k += 1;
                    if k == nsq {
                        break 'dec;
                    }
                }
            }
            let mut m = [[0i64; 16]; 16];
            for i in 0..n {
                let mut w = out[i];
                while w != 0 {
                    let j = w.trailing_zeros() as usize;
                    w &= w - 1;
                    m[i][j] = 1;
                    m[j][i] = -1;
                }
            }
            if f3_mode {
                assert!(n % 2 == 1);
                // nullity via rank; then integer kernel basis; then pointed-cone
                // nonneg test by extreme-ray enumeration
                if let Some(basis) = kernel_basis_exact(&m, n, f3_d) {
                    if cone_has_nonneg(&basis, n, f3_d) {
                        semi += 1;
                        use std::io::Write;
                        std::io::stdout().write_all(rec).unwrap();
                    }
                }
                continue;
            }
            if semi_mode {
                // odd n: nullity-1 nonneg kernel <=> cofactor vector v_i =
                // (-1)^i Pf(M_-i) is one-signed allowing zeros, not all zero
                assert!(n % 2 == 1);
                let mut pos = false;
                let mut neg = false;
                let mut nonzero = 0u32;
                for i in 0..n {
                    let p = pf(&m, full & !(1u16 << i));
                    if p != 0 {
                        nonzero += 1;
                        let v = if i % 2 == 0 { p } else { -p };
                        if v > 0 {
                            pos = true;
                        } else {
                            neg = true;
                        }
                    }
                }
                if nonzero > 0 && !(pos && neg) {
                    semi += 1;
                    if nonzero == n as u32 {
                        semi_strict += 1; // strictly positive kernel = CM
                    }
                    if emit {
                        use std::io::Write;
                        std::io::stdout().write_all(rec).unwrap();
                    }
                }
                continue;
            }
            // even n: paradox + connected + fully-mixed, stratified by nullity
            let mut ok = true;
            for i in 0..n {
                if out[i] == 0 || inn[i] == 0 {
                    ok = false;
                    break;
                }
            }
            if !ok {
                continue;
            }
            let mut seen = 1u16;
            let mut fr = 1u16;
            while fr != 0 {
                let mut nf = 0u16;
                let mut f = fr;
                while f != 0 {
                    let v = f.trailing_zeros() as usize;
                    f &= f - 1;
                    nf |= (out[v] | inn[v]) & !seen;
                }
                seen |= nf;
                fr = nf;
            }
            if seen != full {
                continue;
            }
            if !fully_mixed(&out, n) {
                continue;
            }
            // nullity via integer-preserving Gaussian rank over rationals (f64
            // is fine here: entries {-1,0,1}, n <= 10, validated vs Pfaffian)
            let mut a = [[0f64; 16]; 16];
            for i in 0..n {
                for j in 0..n {
                    a[i][j] = m[i][j] as f64;
                }
            }
            let mut rank = 0usize;
            let mut prow = 0usize;
            for col in 0..n {
                let mut piv = usize::MAX;
                let mut best = 1e-9;
                for rr in prow..n {
                    if a[rr][col].abs() > best {
                        best = a[rr][col].abs();
                        piv = rr;
                    }
                }
                if piv == usize::MAX {
                    continue;
                }
                a.swap(prow, piv);
                for rr in (prow + 1)..n {
                    let f = a[rr][col] / a[prow][col];
                    if f != 0.0 {
                        for c in col..n {
                            a[rr][c] -= f * a[prow][c];
                        }
                    }
                }
                rank += 1;
                prow += 1;
            }
            let d = n - rank;
            incl[d] += 1;
            let mut arc64 = [0u64; 16];
            for i in 0..n {
                arc64[i] = out[i] as u64;
            }
            let aut = unsafe { rps_autsize(arc64.as_ptr(), n as c_int) };
            let aut_u = aut as u128;
            assert!(aut > 0.0 && (nfact % aut_u == 0), "n!/|Aut| not integer");
            lab_sum[d] += nfact / aut_u;
        }
        let rem = have - nrec * reclen;
        buf.copy_within(nrec * reclen..have, 0);
        have = rem;
    }
    if f3_mode {
        eprintln!("f3-emit n={} d={}: emitted={}", n, f3_d, semi);
    } else if semi_mode {
        if emit {
            eprintln!("n={}: scanned={} semi_cm={} (strict_cm={})", n, total, semi, semi_strict);
        } else {
            println!("n={}: scanned={} semi_cm={} (strict_cm={})", n, total, semi, semi_strict);
        }
    } else {
        let tot: u64 = incl.iter().sum();
        let strata: Vec<String> = (0..12)
            .filter(|&d| incl[d] > 0)
            .map(|d| format!("nullity{}={} (labeled {})", d, incl[d], lab_sum[d]))
            .collect();
        println!("n={}: scanned={} inclusive={} [{}]", n, total, tot, strata.join(" "));
    }
}
