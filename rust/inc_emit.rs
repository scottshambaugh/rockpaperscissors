// Size the inclusive filter: read `nauty-directg -o -T` (oriented graphs as
// "nv ne from to from to ..."), count how many are paradoxical and how many are
// paradoxical+connected -- the set the (expensive) fully-mixed test would run on.
//   nauty-geng 7 | nauty-directg -o -T | /tmp/incsizer 7
use std::io::{self, BufRead, BufReader};

fn main() {
    let stdin = io::stdin();
    let mut h = BufReader::with_capacity(1 << 20, stdin.lock());
    let (mut total, mut parad, mut pc) = (0u64, 0u64, 0u64);
    let mut line = String::new();
    loop {
        line.clear();
        if h.read_line(&mut line).unwrap() == 0 {
            break;
        }
        let mut it = line.split_ascii_whitespace().map(|s| s.parse::<usize>().unwrap());
        let nv = match it.next() {
            Some(v) => v,
            None => continue,
        };
        let ne = it.next().unwrap();
        // adj[i] bit j => i beats j ; und[i] bit j => decisive edge (either way)
        let mut adj = vec![0u16; nv];
        let mut und = vec![0u16; nv];
        for _ in 0..ne {
            let a = it.next().unwrap();
            let b = it.next().unwrap();
            adj[a] |= 1u16 << b;
            und[a] |= 1u16 << b;
            und[b] |= 1u16 << a;
        }
        total += 1;
        // paradoxical: every vertex has >=1 win and >=1 loss
        let mut win = vec![0u16; nv];
        for i in 0..nv {
            win[i] = adj[i];
        }
        let mut loss = vec![0u16; nv];
        for i in 0..nv {
            for j in 0..nv {
                if adj[j] & (1 << i) != 0 {
                    loss[i] |= 1 << j;
                }
            }
        }
        let parad_ok = (0..nv).all(|i| win[i] != 0 && loss[i] != 0);
        if !parad_ok {
            continue;
        }
        parad += 1;
        // connected on decisive edges (BFS from 0)
        let mut seen = 1u16; // bit0
        let mut frontier = 1u16;
        while frontier != 0 {
            let mut nf = 0u16;
            let mut f = frontier;
            while f != 0 {
                let v = f.trailing_zeros() as usize;
                f &= f - 1;
                nf |= und[v] & !seen;
            }
            seen |= nf;
            frontier = nf;
        }
        if seen.count_ones() as usize == nv {
            let mut code = String::with_capacity(nv*(nv-1)/2);
            for i in 0..nv { for j in (i+1)..nv {
                code.push(if adj[i] & (1<<j) != 0 { '1' } else if adj[j] & (1<<i) != 0 { '2' } else { '0' });
            }}
            println!("{}", code);
            pc += 1;
        }
    }
    eprintln!("total={} paradoxical={} paradoxical_connected={}", total, parad, pc);
}
