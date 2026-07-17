// Fast-path INCLUSIVE filter (paradoxical + connected + a fully-mixed equilibrium
// exists) for odd n. The fully-mixed question is normally an LP on every candidate
// -- which is why the earlier inclusive enumeration stopped at n=8. But for a skew
// game the answer is governed by the nullity, which the Pfaffian cofactors give
// essentially for free:
//
//   * some Pf(M_-i) != 0  =>  rank n-1, nullity 1: the kernel is the single vector
//       v_i = (-1)^i Pf(M_-i), so a fully-mixed equilibrium exists iff v is strictly
//       one-signed -- exactly the completely-mixed test. No LP.
//   * all Pf(M_-i) == 0   =>  rank <= n-3, nullity >= 3 (n odd): the kernel is a
//       >=3-dimensional subspace; only here do we run the Phase-1 LP to test whether
//       it meets the open positive orthant.
//
// Nullity 1 is generic, so the LP fires on a tiny fraction and the average cost is
// ~the Pfaffian test rather than the LP. inclusive = CM (nullity 1) + {nullity>=3
// inclusive}; the CM sub-count re-confirms the extension-method census (n=9 ->
// 583591020). pf() and fully_mixed() live in common.rs
// (validated in cm_filter.rs / here), so this is a thin combiner over vetted code.
//
//   rustc -O rust/inc_fast.rs -o /tmp/incf
//   nauty-geng 7 2>/dev/null | nauty-directg -o 2>/dev/null | /tmp/incf 7   # 10525

use std::env;
use std::io::{self, Read};

mod common;
use common::{fully_mixed, pf};

fn main() {
    let n: usize = env::args().nth(1).and_then(|s| s.parse().ok()).expect("usage: incf n");
    assert!(n >= 3 && n % 2 == 1 && n <= 15, "n must be odd, 3..=15");
    let nsq = n * n;
    let nbytes = (nsq + 5) / 6;
    let reclen = 2 + nbytes + 1;
    let mut stdin = io::stdin().lock();
    let mut buf = vec![0u8; 1 << 20];
    let mut have = 0usize;
    let (mut total, mut inc, mut cm, mut hi_null, mut lp_calls) = (0u64, 0u64, 0u64, 0u64, 0u64);
    let full: u16 = (1u16 << n) - 1;
    loop {
        let got = stdin.read(&mut buf[have..]).unwrap();
        if got == 0 {
            break;
        }
        have += got;
        let nrec = have / reclen;
        for r in 0..nrec {
            let rec = &buf[r * reclen..(r + 1) * reclen];
            assert!(
                rec[0] == b'&' && rec[1] as usize == 63 + n && rec[reclen - 1] == b'\n',
                "misaligned digraph6 record"
            );
            total += 1;
            let mut out = [0u16; 16];
            let mut inn = [0u16; 16];
            let payload = &rec[2..2 + nbytes];
            let mut k = 0usize;
            'decode: for &byte in payload {
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
                        break 'decode;
                    }
                }
            }
            // paradox prefilter: every vertex needs a win and a loss
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
            // connected on the decisive (out|in) graph
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
            if seen.count_ones() as usize != n {
                continue;
            }
            // integer skew matrix for the Pfaffian cofactors
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
            // Classify by nullity via the n principal Pfaffians v_i = (-1)^i Pf(M_-i),
            // with early exit: once a nonzero cofactor coexists with a zero one or a
            // sign clash, the game is nullity 1 but its kernel isn't strictly one-
            // signed -- not inclusive -- so we stop without computing the rest. Only
            // completely-mixed games (all cofactors nonzero, alternating) run the full
            // n; only all-zero games (nullity >= 3) fall through to the LP.
            let mut seen_zero = false;
            let mut seen_nonzero = false;
            let mut sign_clash = false;
            let mut want = 0i64;
            let mut i = 0usize;
            while i < n {
                let p = pf(&m, full & !(1u16 << i));
                if p == 0 {
                    seen_zero = true;
                } else {
                    seen_nonzero = true;
                    let sgn = if i % 2 == 0 { p.signum() } else { -p.signum() };
                    if want == 0 {
                        want = sgn;
                    } else if sgn != want {
                        sign_clash = true;
                    }
                }
                if seen_nonzero && (seen_zero || sign_clash) {
                    break; // nullity 1, kernel not strictly positive -> not inclusive
                }
                i += 1;
            }
            if i < n {
                // early-exited: nullity 1, not inclusive
            } else if seen_nonzero {
                // nullity 1, all cofactors nonzero + alternating -> completely mixed
                inc += 1;
                cm += 1;
            } else {
                // all cofactors zero -> nullity >= 3 -> the >=3-dim kernel needs the LP
                hi_null += 1;
                lp_calls += 1;
                if fully_mixed(&out, n) {
                    inc += 1;
                }
            }
        }
        let rem = have - nrec * reclen;
        buf.copy_within(nrec * reclen..have, 0);
        have = rem;
    }
    assert_eq!(have, 0, "trailing partial record");
    println!(
        "n={}: candidates={} inclusive={} (cm/nullity1={} nullity>=3_inclusive={}) lp_calls={}",
        n,
        total,
        inc,
        cm,
        inc - cm,
        lp_calls
    );
    eprintln!("  [nullity>=3 games seen: {}]", hi_null);
}
