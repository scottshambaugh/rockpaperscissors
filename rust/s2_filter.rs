// S2 ("two-paradox" / Erdos-Schutte) predicate filter for nauty's gentourng.
//
//   rustc -O rust/s2_filter.rs -o /tmp/s2filter
//   nauty-gentourng <n> 2>/dev/null | /tmp/s2filter <n>
//
// Relies on nauty (gentourng) for the isomorph-free generation -- the hard,
// trustworthy part -- and only applies the S2 test here. gentourng's default
// output is the upper triangle row-by-row in ascii, char at (i<j) = '1' iff i
// beats j. Counts total lines (must equal A000568(n)) and the S2 subset
// (every pair of vertices has a common dominator).
//
// Supports gentourng's res/mod splitting: run several shards and sum the S2
// counts (and totals) yourself.

use std::env;
use std::io::{self, BufRead, BufReader};

fn main() {
    let n: usize = env::args().nth(1).and_then(|s| s.parse().ok()).expect("usage: s2filter n");
    assert!(n <= 16);
    let mut total: u64 = 0;
    let mut s2: u64 = 0;
    let stdin = io::stdin();
    let mut buf = Vec::with_capacity(128);
    let mut h = BufReader::with_capacity(1 << 20, stdin.lock());
    loop {
        buf.clear();
        if h.read_until(b'\n', &mut buf).expect("read") == 0 {
            break;
        }
        // strip trailing newline / whitespace
        let mut len = buf.len();
        while len > 0 && (buf[len - 1] == b'\n' || buf[len - 1] == b'\r') {
            len -= 1;
        }
        if len == 0 || buf[0] == b'>' {
            continue; // skip blank lines and any ">A"/">Z" header
        }
        // parse upper triangle -> beaters[v] = bitmask of vertices that beat v
        let mut beaters = [0u16; 16];
        let mut idx = 0usize;
        for i in 0..n {
            for j in (i + 1)..n {
                let c = buf[idx];
                idx += 1;
                if c == b'1' {
                    beaters[j] |= 1u16 << i; // i beats j
                } else {
                    beaters[i] |= 1u16 << j; // j beats i
                }
            }
        }
        total += 1;
        if total % 20_000_000 == 0 {
            eprintln!("[s2 n={}] {} tournaments scanned, {} S2 so far", n, total, s2);
        }
        let mut ok = true;
        'outer: for i in 0..n {
            for j in (i + 1)..n {
                if beaters[i] & beaters[j] == 0 {
                    ok = false;
                    break 'outer;
                }
            }
        }
        if ok {
            s2 += 1;
        }
    }
    println!("n={}: total={} S2(two-paradox)={}", n, total, s2);
}
