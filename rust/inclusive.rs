// Count inclusive games at size n: read `nauty-directg -o -T` (oriented graphs as
// "nv ne from to ..."), keep the paradoxical + connected ones with a fully-mixed
// equilibrium. nauty does the isomorph-free generation; this only filters, fast
// enough to run on the ~10^8 candidates at n=8.
//
//   rustc -O rust/inclusive.rs -o /tmp/inc
//   nauty-geng 8 2>/dev/null | nauty-directg -o -T 2>/dev/null | /tmp/inc 8
//
// Fully-mixed test: a symmetric (skew M) zero-sum game has a fully-supported
// optimal strategy iff ker(M) meets the open positive orthant, i.e. there is
// x > 0 with Mx = 0. Scaling, that is feasibility of {x >= 1, Mx = 0}; with
// x = 1 + z (z >= 0) it becomes {Mz = -M*ones, z >= 0}, a Phase-1 LP. M has
// entries in {-1,0,1} so float + a tolerance is safe (validated against the
// Python has_fully_mixed for all n <= 7).

use std::env;
use std::io::{self, BufRead, BufReader};

const EPS: f64 = 1e-7;

// Phase-1 simplex feasibility of { A z = b, z >= 0 }, A is m x k. Returns true
// if feasible. Bland's rule (anti-cycling). Small dense tableau.
fn feasible(a: &[Vec<f64>], b: &[f64], m: usize, k: usize) -> bool {
    // make b >= 0 by negating rows
    let mut t = vec![vec![0f64; k + m + 1]; m]; // [ A | I(artificials) | b ]
    for i in 0..m {
        let s = if b[i] < 0.0 { -1.0 } else { 1.0 };
        for j in 0..k {
            t[i][j] = s * a[i][j];
        }
        t[i][k + i] = 1.0; // artificial
        t[i][k + m] = s * b[i];
    }
    // basis = artificials
    let mut basis: Vec<usize> = (0..m).map(|i| k + i).collect();
    // objective: minimize sum of artificials. reduced costs = -sum of constraint rows
    // (for artificial columns cost 1, others 0). We track the objective row.
    let ncol = k + m + 1;
    let mut obj = vec![0f64; ncol];
    // reduced costs for the all-artificial basis: z-columns and rhs are -(col sum);
    // artificial columns are basic, so their reduced cost is 0.
    for j in 0..k {
        obj[j] = -(0..m).map(|i| t[i][j]).sum::<f64>();
    }
    obj[k + m] = -(0..m).map(|i| t[i][k + m]).sum::<f64>();
    loop {
        // entering: first column (Bland) with negative reduced cost among z + artificials
        let mut piv_col = usize::MAX;
        for j in 0..(k + m) {
            if obj[j] < -EPS {
                piv_col = j;
                break;
            }
        }
        if piv_col == usize::MAX {
            break; // optimal
        }
        // ratio test
        let mut piv_row = usize::MAX;
        let mut best = f64::INFINITY;
        for i in 0..m {
            if t[i][piv_col] > EPS {
                let r = t[i][k + m] / t[i][piv_col];
                if r < best - EPS || (r < best + EPS && (piv_row == usize::MAX || basis[i] < basis[piv_row])) {
                    best = r;
                    piv_row = i;
                }
            }
        }
        if piv_row == usize::MAX {
            break; // unbounded (won't happen in phase-1 min)
        }
        // pivot
        let pv = t[piv_row][piv_col];
        for j in 0..ncol {
            t[piv_row][j] /= pv;
        }
        for i in 0..m {
            if i != piv_row {
                let f = t[i][piv_col];
                if f.abs() > EPS {
                    for j in 0..ncol {
                        t[i][j] -= f * t[piv_row][j];
                    }
                }
            }
        }
        let f = obj[piv_col];
        if f.abs() > EPS {
            for j in 0..ncol {
                obj[j] -= f * t[piv_row][j];
            }
        }
        basis[piv_row] = piv_col;
    }
    // objective value = -obj[rhs]; feasible iff sum of artificials ~ 0
    (-obj[k + m]).abs() < 1e-6
}

fn fully_mixed(adj: &[u16], n: usize) -> bool {
    // M[i][j] = +1 if i beats j, -1 if j beats i, 0 tie.  b = -M*ones = -rowsum.
    let mut mat = vec![vec![0f64; n]; n];
    let mut b = vec![0f64; n];
    for i in 0..n {
        let mut rs = 0f64;
        for j in 0..n {
            let v = if adj[i] & (1 << j) != 0 {
                1.0
            } else if adj[j] & (1 << i) != 0 {
                -1.0
            } else {
                0.0
            };
            mat[i][j] = v;
            rs += v;
        }
        b[i] = -rs;
    }
    feasible(&mat, &b, n, n)
}

fn main() {
    let n: usize = env::args().nth(1).and_then(|s| s.parse().ok()).expect("usage: inc n");
    let stdin = io::stdin();
    let mut h = BufReader::with_capacity(1 << 20, stdin.lock());
    let (mut total, mut inc) = (0u64, 0u64);
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
        let mut adj = vec![0u16; nv];
        let mut und = vec![0u16; nv];
        for _ in 0..ne {
            let a = it.next().unwrap();
            let bb = it.next().unwrap();
            adj[a] |= 1u16 << bb;
            und[a] |= 1u16 << bb;
            und[bb] |= 1u16 << a;
        }
        total += 1;
        if total % 10_000_000 == 0 {
            eprintln!("[inc n={}] {} candidates scanned, {} inclusive so far", n, total, inc);
        }
        // paradoxical: every vertex has a win and a loss
        let mut ok = true;
        for i in 0..nv {
            let win = adj[i] != 0;
            let mut loss = false;
            for j in 0..nv {
                if adj[j] & (1 << i) != 0 {
                    loss = true;
                    break;
                }
            }
            if !win || !loss {
                ok = false;
                break;
            }
        }
        if !ok {
            continue;
        }
        // connected on decisive edges
        let mut seen = 1u16;
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
        if seen.count_ones() as usize != nv {
            continue;
        }
        if fully_mixed(&adj, nv) {
            inc += 1;
        }
    }
    println!("n={}: candidates={} inclusive={}", n, total, inc);
}
