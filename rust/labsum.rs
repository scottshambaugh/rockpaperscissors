// Labeled-count summation over a digraph6 class stream: reads n-vertex
// digraph6 records (one per iso class) and reports count and
// sum of n!/|Aut| (the labeled count), via the nauty autsize shim.
//
//   rustc -O rust/labsum.rs -o /tmp/labsum -C link-args="shim.o -lnauty"
//   sort -m -u out/c_*.d6 | /tmp/labsum 9
use std::env;
use std::io::{self, Read};
use std::os::raw::c_int;

mod common;
use common::factorial;

extern "C" {
    fn rps_autsize(arc: *const u64, n: c_int) -> f64;
}

fn main() {
    let n: usize = env::args().nth(1).and_then(|s| s.parse().ok()).expect("usage: labsum n < classes.d6");
    let nsq = n * n;
    let nbytes = (nsq + 5) / 6;
    let reclen = 2 + nbytes + 1;
    let nfact = factorial(n as u64);

    let mut stdin = io::stdin().lock();
    let mut buf = vec![0u8; 1 << 20];
    let mut have = 0usize;
    let mut classes = 0u64;
    let mut labeled: u128 = 0;

    loop {
        let got = stdin.read(&mut buf[have..]).unwrap();
        if got == 0 {
            break;
        }
        have += got;
        let nrec = have / reclen;
        for ri in 0..nrec {
            let rec = &buf[ri * reclen..(ri + 1) * reclen];
            assert!(rec[0] == b'&' && rec[1] as usize == 63 + n, "misaligned digraph6");
            classes += 1;
            let mut beats = [0u64; 16];
            let payload = &rec[2..2 + nbytes];
            let mut kk = 0usize;
            'dec: for &byte in payload {
                let mut bits = ((byte - 63) as u32) << 26;
                for _ in 0..6 {
                    if bits & 0x8000_0000 != 0 {
                        beats[kk / n] |= 1 << (kk % n);
                    }
                    bits <<= 1;
                    kk += 1;
                    if kk == nsq {
                        break 'dec;
                    }
                }
            }
            let aut = unsafe { rps_autsize(beats.as_ptr(), n as c_int) } as u128;
            debug_assert!(nfact % aut == 0);
            labeled += nfact / aut;
        }
        let rem = have - nrec * reclen;
        buf.copy_within(nrec * reclen..have, 0);
        have = rem;
    }
    println!("n={} classes={} labeled={}", n, classes, labeled);
}
