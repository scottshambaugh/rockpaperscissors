// Read gentourng default ascii (upper-triangle row-by-row) on stdin.
// For each tournament: test S2 (two-paradox = every pair has a common dominator),
// and among S2 ones test primality (no nontrivial module).
// Usage: nauty-gentourng N | s2_prime N
use std::io::{self, BufRead, Write};

fn main() {
    let n: usize = std::env::args().nth(1).unwrap().parse().unwrap();
    let full: u32 = if n == 32 { u32::MAX } else { (1u32 << n) - 1 };

    let mut total: u64 = 0;
    let mut s2: u64 = 0;
    let mut s2_prime: u64 = 0;

    // beats[i] = bitmask of vertices i beats; beatenby[i] = mask of vertices that beat i
    let mut beats = [0u32; 32];
    let mut beatenby = [0u32; 32];

    let stdin = io::stdin();
    let mut lock = stdin.lock();
    let mut line = Vec::with_capacity(64);
    loop {
        line.clear();
        let read = lock.read_until(b'\n', &mut line).unwrap();
        if read == 0 { break; }
        let b: &[u8] = if line.last() == Some(&b'\n') { &line[..line.len() - 1] } else { &line[..] };
        if b.len() != n * (n - 1) / 2 { continue; }
        for x in beats.iter_mut().take(n) { *x = 0; }
        let mut k = 0usize;
        for i in 0..n {
            for j in (i + 1)..n {
                // convention: '1' at (i,j) => i beats j
                if b[k] == b'1' { beats[i] |= 1 << j; } else { beats[j] |= 1 << i; }
                k += 1;
            }
        }
        for i in 0..n { beatenby[i] = (full ^ (1 << i)) & !beats[i]; }
        total += 1;

        // S2: every pair {u,v} has a common dominator (some w beats both)
        let mut ok = true;
        'outer: for u in 0..n {
            for v in (u + 1)..n {
                if beatenby[u] & beatenby[v] == 0 { ok = false; break 'outer; }
            }
        }
        if !ok { continue; }
        s2 += 1;

        // primality: for each pair grow smallest module; if any proper => not prime
        if is_prime(&beats, n, full) { s2_prime += 1; }
    }

    let mut out = io::stdout();
    writeln!(out, "n={} total={} S2={} S2_prime={} S2_nonprime={}",
             n, total, s2, s2_prime, s2 - s2_prime).unwrap();
}

// A module M: every vertex outside M relates uniformly to all of M.
// Prime = no module with 2 <= |M| <= n-1.
fn is_prime(beats: &[u32; 32], n: usize, full: u32) -> bool {
    for a in 0..n {
        for b_ in (a + 1)..n {
            let mut m: u32 = (1 << a) | (1 << b_);
            loop {
                let mut added = 0u32;
                let outside = full & !m;
                let mut rest = outside;
                while rest != 0 {
                    let x = rest.trailing_zeros() as usize;
                    rest &= rest - 1;
                    let bm = beats[x] & m;      // members of M that x beats
                    if bm != 0 && bm != m { added |= 1 << x; } // x splits M
                }
                if added == 0 { break; }
                m |= added;
            }
            if (m & full).count_ones() < n as u32 { return false; } // proper module >=2
        }
    }
    true
}
