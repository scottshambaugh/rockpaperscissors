// Quotient filter for the balanced prime-gap count. A twin-free non-prime
// balanced n-game is a substitution G = H[M -> special vertex]: a module M of
// size t (3 <= t, and t-1 <= k-1 where k = n+1-t) whose internal row sums are
// all 0 (they must be equal, and they sum to 0), substituted at a special
// vertex i* of a k-vertex quotient H satisfying
//     rowsum_H(i*) = 0   and   rowsum_H(v) = -(t-1) * rel(v, i*)  for v != i*.
// This tool scans a digraph6 stream of k-vertex oriented games and, for every
// game and EVERY valid choice of special vertex, emits the game relabeled with
// the special vertex first (the substitution driver's convention).
//
//   rustc -O rust/qfilter.rs -o /tmp/qf
//   nauty-geng 8 | nauty-directg -o | /tmp/qf 8 3   # quotients for t=3, k=8
use std::env;
use std::io::{self, Read, Write};

mod common;
use common::encode;

fn main() {
    let k: usize = env::args().nth(1).and_then(|s| s.parse().ok()).expect("usage: qf k t");
    let t: i32 = env::args().nth(2).and_then(|s| s.parse().ok()).expect("usage: qf k t");
    let nsq = k * k;
    let nbytes = (nsq + 5) / 6;
    let reclen = 2 + nbytes + 1;
    let mut stdin = io::stdin().lock();
    let mut w = io::BufWriter::new(io::stdout().lock());
    let mut buf = vec![0u8; 1 << 20];
    let mut have = 0usize;
    let (mut total, mut emitted) = (0u64, 0u64);
    let mut child = Vec::with_capacity(reclen);
    loop {
        let got = stdin.read(&mut buf[have..]).unwrap();
        if got == 0 {
            break;
        }
        have += got;
        let nrec = have / reclen;
        for r in 0..nrec {
            let rec = &buf[r * reclen..(r + 1) * reclen];
            assert!(rec[0] == b'&' && rec[1] as usize == 63 + k, "misaligned digraph6");
            total += 1;
            let mut beats = [0u16; 16];
            let payload = &rec[2..2 + nbytes];
            let mut bit = 0usize;
            'dec: for &byte in payload {
                let mut bits = ((byte - 63) as u32) << 26;
                for _ in 0..6 {
                    if bits & 0x8000_0000 != 0 {
                        beats[bit / k] |= 1 << (bit % k);
                    }
                    bits <<= 1;
                    bit += 1;
                    if bit == nsq {
                        break 'dec;
                    }
                }
            }
            let mut rowsum = [0i32; 16];
            for i in 0..k {
                for j in 0..k {
                    if beats[i] & (1 << j) != 0 {
                        rowsum[i] += 1;
                        rowsum[j] -= 1;
                    }
                }
            }
            for v0 in 0..k {
                if rowsum[v0] != 0 {
                    continue;
                }
                let mut ok = true;
                for v in 0..k {
                    if v == v0 {
                        continue;
                    }
                    let rel = if beats[v] & (1 << v0) != 0 {
                        1
                    } else if beats[v0] & (1 << v) != 0 {
                        -1
                    } else {
                        0
                    };
                    if rowsum[v] != -(t - 1) * rel {
                        ok = false;
                        break;
                    }
                }
                if !ok {
                    continue;
                }
                // relabel: v0 -> position 0, others keep relative order
                let mut map = [0usize; 16];
                map[v0] = 0;
                let mut pos = 1;
                for v in 0..k {
                    if v != v0 {
                        map[v] = pos;
                        pos += 1;
                    }
                }
                let mut nb = [0u16; 16];
                for i in 0..k {
                    for j in 0..k {
                        if beats[i] & (1 << j) != 0 {
                            nb[map[i]] |= 1 << map[j];
                        }
                    }
                }
                encode(&nb[..k], k, &mut child);
                w.write_all(&child).unwrap();
                emitted += 1;
            }
        }
        let rem = have - nrec * reclen;
        buf.copy_within(nrec * reclen..have, 0);
        have = rem;
    }
    w.flush().unwrap();
    eprintln!("qf k={} t={}: scanned={} quotient-emissions={}", k, t, total, emitted);
}
