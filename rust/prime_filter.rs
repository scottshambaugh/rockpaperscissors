// Count modular-prime games in a digraph6 stream. Reads oriented n-graphs
// (`nauty-directg -o` output, or `nauty-labelg` canonical forms) on stdin and
// reports how many are modular-prime -- no nontrivial module (a set S every
// outside move relates to identically). Used to get the *prime* subcount of an
// already-filtered family, e.g. the completely mixed games:
//
//   rustc -O rust/prime_filter.rs -o /tmp/prime
//   cat ext/c_*.d6 | /tmp/prime 9        # prime subcount of the n=9 CM games
//
// The is_prime module-closure test (O(n^2) per pair, validated in balanced.rs
// against the Python is_prime and the balanced/regular [prime] census columns)
// lives in common.rs, so this stays a thin reader over the same vetted primality
// logic, not a re-implementation.

use std::env;
use std::io::{self, Read};

mod common;
use common::{is_prime, twin_free, MAXN};

fn main() {
    let n: usize = env::args().nth(1).and_then(|s| s.parse().ok()).expect("usage: prime n");
    assert!((3..=15).contains(&n), "n out of range");
    let nsq = n * n;
    let nbytes = (nsq + 5) / 6;
    let reclen = 2 + nbytes + 1; // '&' + N(n) + payload + '\n'
    let mut stdin = io::stdin().lock();
    let mut buf = vec![0u8; 1 << 20];
    let mut have = 0usize;
    let (mut total, mut prime, mut tf) = (0u64, 0u64, 0u64);
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
            // decode row-major arc bits: arc[i] bit j <=> i beats j
            let mut arc = [0u64; MAXN];
            let payload = &rec[2..2 + nbytes];
            let mut k = 0usize;
            'decode: for &byte in payload {
                let mut bits = ((byte - 63) as u32) << 26;
                for _ in 0..6 {
                    if bits & 0x8000_0000 != 0 {
                        arc[k / n] |= 1u64 << (k % n);
                    }
                    bits <<= 1;
                    k += 1;
                    if k == nsq {
                        break 'decode;
                    }
                }
            }
            if is_prime(&arc, n) {
                prime += 1;
            }
            if twin_free(&arc, n) {
                tf += 1;
            }
        }
        let rem = have - nrec * reclen;
        buf.copy_within(nrec * reclen..have, 0);
        have = rem;
    }
    assert_eq!(have, 0, "trailing partial record");
    // "total=... prime=..." stays a contiguous substring (ci_test.sh greps it)
    println!("n={}: total={} prime={} twin_free={}", n, total, prime, tf);
}
