// Count completely mixed games at odd size n: read digraph6 from
// `nauty-directg -o` (oriented graphs, one representative per iso class), keep
// the games whose EVERY Nash equilibrium is fully mixed (Kaplansky's completely
// mixed games). nauty does the isomorph-free generation; this only filters.
//
//   rustc -O rust/cm_filter.rs -o /tmp/cmf
//   nauty-geng 7 2>/dev/null | nauty-directg -o 2>/dev/null | /tmp/cmf 7
//
// Test (exact integer arithmetic, no LP, no floats): for odd skew-symmetric M
// the vector v with v_i = (-1)^i Pf(M with row+col i deleted) satisfies Mv = 0
// identically, and spans ker(M) exactly when rank = n-1. So the game is
// completely mixed  <=>  all n principal Pfaffians are nonzero and v is strictly
// one-signed (Kaplansky 1995). Positivity of v subsumes everything else: a
// vertex with no loss (or no win) would force (Mv)_i != 0, and a decisively
// disconnected game puts zeros in v -- so paradoxical + connected + fully-mixed
// + unique all fall out of this single test. A cheap paradox bitmask prefilter
// still runs first because it rejects for ~free. Early exit: Pfaffians are
// computed one vertex at a time and abandoned on the first zero or sign clash.
//
// Entries: |Pf| of a k x k {-1,0,1} skew matrix is at most sqrt(det) <=
// k^(k/4) (Hadamard), i.e. <= 4096 for k = 12, so i64 has room to spare
// through any n this will ever run at.
//
// Validated against the Python census: n=5 -> 7, n=7 -> 7268 (and the n=8
// candidate total 575016219 = A001174(8) as a completeness checksum; the CM
// count at even n must be 0 -- the parity theorem -- which doubles as a
// correctness check of the sign logic).

use std::env;
use std::io::{self, Read, Write};

mod common;
use common::pf;

// Completely mixed test: all n Pfaffian cofactors nonzero, alternating signs
// (so that v_i = (-1)^i Pf_i is one-signed). Early exit on first failure.
fn completely_mixed(m: &[[i64; 16]; 16], n: usize) -> bool {
    let full: u16 = (1u16 << n) - 1;
    let mut want = 0i64; // sign of v_0, fixed by the first cofactor
    for i in 0..n {
        let p = pf(m, full & !(1u16 << i));
        if p == 0 {
            return false;
        }
        let v = if i % 2 == 0 { p } else { -p };
        if i == 0 {
            want = v.signum();
        } else if v.signum() != want {
            return false;
        }
    }
    true
}

fn main() {
    let n: usize = env::args().nth(1).and_then(|s| s.parse().ok()).expect("usage: cmf n [emit]");
    assert!(n >= 3 && n <= 15, "n out of range");
    // "emit": re-print the digraph6 of each completely mixed game (for piping to
    // nauty-labelg to cross-check completeness against another enumeration).
    let emit = env::args().nth(2).as_deref() == Some("emit");
    let nsq = n * n;
    let nbytes = (nsq + 5) / 6; // digraph6 payload bytes after '&' + N(n)
    let reclen = 2 + nbytes + 1; // '&' + size char + payload + '\n'
    let mut stdin = io::stdin().lock();
    let mut out_w = io::BufWriter::new(io::stdout().lock());
    let mut buf = vec![0u8; 1 << 20];
    let mut have = 0usize;
    let (mut total, mut cm) = (0u64, 0u64);
    let stderr = io::stderr();
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
                "misaligned digraph6 record (unexpected header line?)"
            );
            total += 1;
            if total % 100_000_000 == 0 {
                eprint!("\r[cmf n={}] {}e8 scanned, {} completely mixed", n, total / 100_000_000, cm);
                stderr.lock().flush().unwrap();
            }
            // decode row-major adjacency bits into out/in masks
            let mut out = [0u16; 16];
            let mut inn = [0u16; 16];
            let payload = &rec[2..2 + nbytes];
            let mut k = 0usize; // bit index i*n + j
            'decode: for byte in payload {
                let mut bits = (byte - 63) as u32;
                bits <<= 26; // top 6 bits of a u32
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
            if completely_mixed(&m, n) {
                cm += 1;
                if emit {
                    out_w.write_all(rec).unwrap();
                }
            }
        }
        // move any partial record to the front
        let rem = have - nrec * reclen;
        buf.copy_within(nrec * reclen..have, 0);
        have = rem;
    }
    assert_eq!(have, 0, "trailing partial record");
    out_w.flush().unwrap();
    drop(out_w); // release the stdout lock the BufWriter holds before printing below
    // In emit mode the digraph6 records own stdout, so the summary goes to stderr;
    // otherwise it stays on stdout (where the count-only callers/tests read it).
    let summary = format!("n={}: candidates={} completely_mixed={}", n, total, cm);
    if emit {
        eprintln!("{summary}");
    } else {
        println!("{summary}");
    }
}
