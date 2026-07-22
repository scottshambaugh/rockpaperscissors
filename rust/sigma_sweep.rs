// Marked-pair statistics over an n-vertex digraph6 class stream, for the
// (2,2,1^{n-4}) sigma-correction type whose quotient (m = n-2 supervertices,
// one g=2 bundle between the two 2-cycles) is too large to brute-enumerate at
// n=10. Decomposing that type by bundle profile (t1,t2):
//
//   Fix_(2,2,1^{f}) on f+4 vertices = 2*J1 + J0 + 4*H1 + 2*K0b
//
// where, over LABELED plain games on n = f+2 vertices with marked ordered
// pair (x, y):
//   J1  = #[rel(x,y)=+1, game inclusive]          (profiles (1,1) and (-1,-1))
//   J0  = #[rel(x,y)=0,  game inclusive]          (profile (0,0))
//   K0b = #[rel(x,y)=0, balance-pos-kernel, paradox except x,y, connected
//          with the extra lifted edge x-y]        (profiles (1,-1), (-1,1))
//   H1  = #[rel(x,y)=0 in the plain matrix, balance on M with the (x,y)
//          entry replaced by +1/2 (tested on 2M with a +-1 entry), x-row has
//          a loss, y-row has a win, others fully paradoxical, and the LIFTED
//          10-vertex graph (x,y split into 2-cycles, single bundle class
//          x1-y1, x2-y2) connected]               (profiles (1,0),(0,1),(-1,0),(0,-1))
//
// Labeled counts come from classes: each class C contributes
// (n!/|Aut(C)|) * #qualifying ordered pairs / (n(n-1)); the division is done
// once at the end (asserted exact). The same pass reports L_inc(n) as a
// cross-check against the strata engines.
//
// Anchor: n=6 stream (21480 classes) must give 2*J1+J0+4*H1+2*K0b = 339138
// = sigma_fix 2,2,1,1,1,1 (the brute enumeration of the same type at n=8).
use std::env;
use std::io::{self, Read};
use std::os::raw::c_int;

mod common;
use common::{connected_beats, factorial, has_positive_kernel, paradox_connected_beats};

extern "C" {
    fn rps_autsize(arc: *const u64, n: c_int) -> f64;
}

// connectivity of the plain graph plus one extra undirected edge x-y
fn connected_plus_edge(beats: &[u16; 16], n: usize, x: usize, y: usize) -> bool {
    let mut b = *beats;
    b[x] |= 1 << y;
    connected_beats(&b, n)
}

// connectivity of the lifted graph: x -> {x, n}, y -> {y, n+1}; copies share
// the original neighbourhoods; bundle edges (x,y) and (n, n+1) only
fn connected_lifted(beats: &[u16; 16], inn: &[u16; 16], n: usize, x: usize, y: usize) -> bool {
    let mut adj = [0u16; 16];
    for u in 0..n {
        if u == x || u == y {
            continue;
        }
        let nb = (beats[u] | inn[u]) & !(1 << x) & !(1 << y);
        adj[u] |= nb;
    }
    // undirected closure of the single-vertex part
    for u in 0..n {
        if u == x || u == y {
            continue;
        }
        let mut w = adj[u];
        while w != 0 {
            let v = w.trailing_zeros() as usize;
            w &= w - 1;
            adj[v] |= 1 << u;
        }
    }
    let xn = (beats[x] | inn[x]) & !(1 << y); // x's plain neighbours (rel(x,y)=0 anyway)
    let yn = (beats[y] | inn[y]) & !(1 << x);
    let x2 = n;
    let y2 = n + 1;
    let mut w = xn;
    while w != 0 {
        let v = w.trailing_zeros() as usize;
        w &= w - 1;
        adj[x] |= 1 << v;
        adj[v] |= 1 << x;
        adj[x2] |= 1 << v;
        adj[v] |= 1 << x2;
    }
    let mut w = yn;
    while w != 0 {
        let v = w.trailing_zeros() as usize;
        w &= w - 1;
        adj[y] |= 1 << v;
        adj[v] |= 1 << y;
        adj[y2] |= 1 << v;
        adj[v] |= 1 << y2;
    }
    // bundle class: x1-y1, x2-y2
    adj[x] |= 1 << y;
    adj[y] |= 1 << x;
    adj[x2] |= 1 << y2;
    adj[y2] |= 1 << x2;
    // BFS over n+2 vertices
    let full: u16 = ((1u32 << (n + 2)) - 1) as u16;
    let mut seen: u16 = 1 << x;
    let mut fr: u16 = 1 << x;
    while fr != 0 {
        let mut nf = 0u16;
        let mut f = fr;
        while f != 0 {
            let v = f.trailing_zeros() as usize;
            f &= f - 1;
            nf |= adj[v] & !seen;
        }
        seen |= nf;
        fr = nf;
    }
    seen == full
}

fn main() {
    let n: usize = env::args().nth(1).and_then(|s| s.parse().ok()).expect("usage: sigma_sweep n < classes.d6");
    assert!(n % 2 == 0);
    let nsq = n * n;
    let nbytes = (nsq + 5) / 6;
    let reclen = 2 + nbytes + 1;
    let nfact = factorial(n as u64);

    let mut stdin = io::stdin().lock();
    let mut buf = vec![0u8; 1 << 20];
    let mut have = 0usize;
    let mut classes = 0u64;
    let mut linc: u128 = 0;
    let (mut nj1, mut nj0, mut nk0b, mut nh1): (u128, u128, u128, u128) = (0, 0, 0, 0);

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
            let mut beats = [0u16; 16];
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
            let mut inn = [0u16; 16];
            for i in 0..n {
                let mut w = beats[i];
                while w != 0 {
                    let j = w.trailing_zeros() as usize;
                    w &= w - 1;
                    inn[j] |= 1 << i;
                }
            }
            // one-signed vertices (have arcs of exactly one direction)
            let mut onesigned = [false; 16];
            let mut ones = 0usize;
            let mut onelist = [0usize; 16];
            for i in 0..n {
                if (beats[i] == 0) != (inn[i] == 0) {
                    onesigned[i] = true;
                    onelist[ones] = i;
                    ones += 1;
                }
            }
            if ones > 2 {
                continue; // no statistic can hold: balance needs rows covered by {x,y}
            }
            let mut m = [[0i64; 16]; 16];
            for i in 0..n {
                let mut w = beats[i];
                while w != 0 {
                    let j = w.trailing_zeros() as usize;
                    w &= w - 1;
                    m[i][j] = 1;
                    m[j][i] = -1;
                }
            }
            let plain_balance = ones == 0 && has_positive_kernel(&m, n);
            let inclusive = plain_balance && paradox_connected_beats(&beats, n);
            let mut aut: u128 = 0; // lazy
            let mut getaut = |beats: &[u16; 16]| -> u128 {
                let mut arc64 = [0u64; 16];
                for i in 0..n {
                    arc64[i] = beats[i] as u64;
                }
                common::autsize_u128(unsafe { rps_autsize(arc64.as_ptr(), n as c_int) })
            };
            let (mut c_j1, mut c_j0, mut c_k0b, mut c_h1) = (0u64, 0u64, 0u64, 0u64);
            if inclusive {
                // J stats: pair counts by rel value
                let mut arcs = 0u64;
                for i in 0..n {
                    arcs += beats[i].count_ones() as u64;
                }
                c_j1 = arcs; // ordered pairs with rel = +1
                c_j0 = (n * (n - 1)) as u64 - 2 * arcs; // ordered ties
            }
            if plain_balance {
                // K0b: ordered tie pairs (x,y), paradox except x,y (automatic:
                // balance forced every row both-signed-or-zero; zero rows are
                // possible!), connected with the extra edge
                for x in 0..n {
                    for y in 0..n {
                        if x == y || beats[x] & (1 << y) != 0 || beats[y] & (1 << x) != 0 {
                            continue;
                        }
                        // vertices other than x,y with zero rows fail paradox
                        let mut ok = true;
                        for u in 0..n {
                            if u != x && u != y && beats[u] == 0 && inn[u] == 0 {
                                ok = false;
                                break;
                            }
                        }
                        if ok && connected_plus_edge(&beats, n, x, y) {
                            c_k0b += 1;
                        }
                    }
                }
            }
            // H1: candidate pairs must cover the one-signed set
            {
                for x in 0..n {
                    for y in 0..n {
                        if x == y || beats[x] & (1 << y) != 0 || beats[y] & (1 << x) != 0 {
                            continue;
                        }
                        // one-signed vertices must be within {x,y}
                        let mut ok = true;
                        for t in 0..ones {
                            let u = onelist[t];
                            if u != x && u != y {
                                ok = false;
                                break;
                            }
                        }
                        if !ok {
                            continue;
                        }
                        // paradox pattern: x needs a loss, y needs a win,
                        // others need both; zero-row others fail
                        if inn[x] == 0 || beats[y] == 0 {
                            continue;
                        }
                        let mut pok = true;
                        for u in 0..n {
                            if u != x && u != y && (beats[u] == 0 || inn[u] == 0) {
                                pok = false;
                                break;
                            }
                        }
                        if !pok {
                            continue;
                        }
                        // balance on 2M with (x,y) entry = +1
                        let mut m2 = [[0i64; 16]; 16];
                        for i in 0..n {
                            for j in 0..n {
                                m2[i][j] = 2 * m[i][j];
                            }
                        }
                        m2[x][y] = 1;
                        m2[y][x] = -1;
                        if !has_positive_kernel(&m2, n) {
                            continue;
                        }
                        if connected_lifted(&beats, &inn, n, x, y) {
                            c_h1 += 1;
                        }
                    }
                }
            }
            if c_j1 + c_j0 + c_k0b + c_h1 > 0 || inclusive {
                if aut == 0 {
                    aut = getaut(&beats);
                }
                let w = nfact / aut;
                if inclusive {
                    linc += w;
                }
                nj1 += w * c_j1 as u128;
                nj0 += w * c_j0 as u128;
                nk0b += w * c_k0b as u128;
                nh1 += w * c_h1 as u128;
            }
        }
        let rem = have - nrec * reclen;
        buf.copy_within(nrec * reclen..have, 0);
        have = rem;
    }
    let pairs = (n * (n - 1)) as u128;
    println!(
        "n={} classes={} L_inc={} rawJ1={} rawJ0={} rawK0b={} rawH1={}",
        n, classes, linc, nj1, nj0, nk0b, nh1
    );
    for (name, v) in [("J1", nj1), ("J0", nj0), ("K0b", nk0b), ("H1", nh1)] {
        if v % pairs == 0 {
            println!("{} = {}", name, v / pairs);
        } else {
            println!("{} = {}/{} NOT INTEGRAL (sum over full stream required)", name, v, pairs);
        }
    }
    let fix = (2 * nj1 + nj0 + 4 * nh1 + 2 * nk0b) / pairs;
    println!("Fix_(2,2,1^{}) = {}", n - 2, fix);
}
