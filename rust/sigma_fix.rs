// Fix_sigma counts for the Burnside sigma-corrections of the inclusive
// census: the number of labeled n-vertex games invariant under a fixed
// permutation sigma of a given cycle type that are inclusive (paradoxical +
// connected + exists w > 0 with Mw = 0).
//
// A sigma-invariant game is an assignment of one trit to each <sigma>-orbit
// of ordered vertex pairs (with rel(y,x) = -rel(x,y) folded in):
//   * within cycle i (length l): classes k = 1..floor((l-1)/2); the antipodal
//     class of an even cycle pairs (x, sigma^{l/2} x) with its own reverse and
//     is a FORCED TIE.
//   * between cycles i < j: g = gcd(li, lj) classes by offset congruence
//     (b - a) mod g; each class is one orbit (orbit size lcm = li*lj/g).
//
// Balance: M sigma-invariant and exists w>0 with Mw=0 iff (averaging over
// <sigma>) an invariant witness exists iff the m x m quotient row-sum matrix
//   R[i][j] = sum over y in cycle j of rel(x0, y) = (lj/g) * (sum of bundle
//   trits),  R[i][i] = 0 always (each within class cancels its mirror)
// has a strictly positive kernel vector. Strict positivity is exact via
// Gordan: no nonneg nonzero dependency among the kernel-basis columns
// (common::positive_dependencies, supports <= d+1 by Caratheodory).
//
// Paradox and connectivity are NOT quotient-local (offset classes matter:
// e.g. two 2-cycles joined by a single class split into two components while
// the quotient graph looks connected), so both are tested exactly on the
// LIFTED n-vertex game. Within-class trits never enter R, so the bundle DFS
// tests balance once per bundle assignment and the within enumeration only
// re-tests the lifted paradox + connectivity.
//
// The cycle type 1,1,...,1 (sigma = id) enumerates ALL labeled m-vertex games
// and thus computes the plain labeled inclusive count L_inc(m) with the same
// code path.
//
//   rustc -O rust/sigma_fix.rs -o /tmp/sfix
//   /tmp/sfix 2,2,1,1,1,1     -> Fix for that type at n = 8
//   /tmp/sfix 1,1,1,1,1       -> L_inc(5)
//   /tmp/sfix TYPE s OF       -> shard s of OF (OF = 3^k): the first k bundle
//                                slots are preset from s's base-3 digits;
//                                fixes sum over shards
use std::env;

mod common;
use common::{connected_beats, has_positive_kernel};

fn gcd(a: usize, b: usize) -> usize {
    if b == 0 { a } else { gcd(b, a % b) }
}

struct Bundle {
    i: usize,
    j: usize,
    wij: i64, // lj / g: R[i][j] contribution per trit
    wji: i64, // li / g: -R[j][i] contribution per trit
    plus: [u16; 16],
    minus: [u16; 16],
}

struct Ctx {
    n: usize,
    m: usize,
    row_check: Vec<Option<usize>>,
    withins: Vec<([u16; 16], [u16; 16])>, // (plus, minus) arc masks
    fix: u64,
    bundle_leaves: u64,
    balance_pass: u64,
}

// paradoxical + connected on the lifted game
fn lifted_ok(beats: &[u16; 16], n: usize) -> bool {
    let mut inn = [0u16; 16];
    for i in 0..n {
        let mut w = beats[i];
        while w != 0 {
            let j = w.trailing_zeros() as usize;
            w &= w - 1;
            inn[j] |= 1 << i;
        }
    }
    for i in 0..n {
        if beats[i] == 0 || inn[i] == 0 {
            return false;
        }
    }
    connected_beats(beats, n)
}

// exists strictly positive w with Rw = 0, exact
fn balance_ok(r: &[[i64; 8]; 8], m: usize) -> bool {
    let mut mm = [[0i64; 16]; 16];
    for i in 0..m {
        for j in 0..m {
            mm[i][j] = r[i][j];
        }
    }
    has_positive_kernel(&mm, m)
}

fn within_dfs(k: usize, beats: [u16; 16], ctx: &mut Ctx) {
    if k == ctx.withins.len() {
        if lifted_ok(&beats, ctx.n) {
            ctx.fix += 1;
        }
        return;
    }
    let (plus, minus) = ctx.withins[k];
    within_dfs(k + 1, beats, ctx);
    let mut b = beats;
    for x in 0..ctx.n {
        b[x] |= plus[x];
    }
    within_dfs(k + 1, b, ctx);
    let mut b = beats;
    for x in 0..ctx.n {
        b[x] |= minus[x];
    }
    within_dfs(k + 1, b, ctx);
}

fn bundle_dfs(k: usize, bundles: &[Bundle], beats: [u16; 16], r: &mut [[i64; 8]; 8], ctx: &mut Ctx) {
    // row-boundary prune: R rows with min-index < bundles[k].i are complete;
    // a completed row that has entries of exactly one sign kills the subtree
    // (no strictly positive kernel can balance it)
    if let Some(row) = ctx.row_check[k] {
        let mut pos = false;
        let mut neg = false;
        for j in 0..ctx.m {
            pos |= r[row][j] > 0;
            neg |= r[row][j] < 0;
        }
        if pos != neg {
            return;
        }
    }
    if k == bundles.len() {
        ctx.bundle_leaves += 1;
        if !balance_ok(r, ctx.m) {
            return;
        }
        ctx.balance_pass += 1;
        within_dfs(0, beats, ctx);
        return;
    }
    let bl = &bundles[k];
    // t = 0
    bundle_dfs(k + 1, bundles, beats, r, ctx);
    // t = +1
    let mut b = beats;
    for x in 0..ctx.n {
        b[x] |= bl.plus[x];
    }
    r[bl.i][bl.j] += bl.wij;
    r[bl.j][bl.i] -= bl.wji;
    bundle_dfs(k + 1, bundles, b, r, ctx);
    r[bl.i][bl.j] -= bl.wij;
    r[bl.j][bl.i] += bl.wji;
    // t = -1
    let mut b = beats;
    for x in 0..ctx.n {
        b[x] |= bl.minus[x];
    }
    r[bl.i][bl.j] -= bl.wij;
    r[bl.j][bl.i] += bl.wji;
    bundle_dfs(k + 1, bundles, b, r, ctx);
    r[bl.i][bl.j] += bl.wij;
    r[bl.j][bl.i] -= bl.wji;
}

fn main() {
    let arg = env::args().nth(1).expect("usage: sigma_fix l1,l2,...");
    let ls: Vec<usize> = arg.split(',').map(|x| x.parse().expect("bad cycle length")).collect();
    let n: usize = ls.iter().sum();
    let m = ls.len();
    assert!(n <= 12 && m <= 7, "brute-mode limits (m=8,9 types use reductions)");
    let mut base = vec![0usize; m];
    for i in 1..m {
        base[i] = base[i - 1] + ls[i - 1];
    }
    let mut bundles: Vec<Bundle> = Vec::new();
    for i in 0..m {
        for j in (i + 1)..m {
            let g = gcd(ls[i], ls[j]);
            for c in 0..g {
                let mut plus = [0u16; 16];
                let mut minus = [0u16; 16];
                for a in 0..ls[i] {
                    for b in 0..ls[j] {
                        if (b + g - (a % g)) % g == c {
                            let x = base[i] + a;
                            let y = base[j] + b;
                            plus[x] |= 1 << y; // t=+1: x beats y
                            minus[y] |= 1 << x; // t=-1: y beats x
                        }
                    }
                }
                bundles.push(Bundle {
                    i,
                    j,
                    wij: (ls[j] / g) as i64,
                    wji: (ls[i] / g) as i64,
                    plus,
                    minus,
                });
            }
        }
    }
    let mut withins: Vec<([u16; 16], [u16; 16])> = Vec::new();
    for i in 0..m {
        let l = ls[i];
        if l < 3 {
            continue;
        }
        for k in 1..=((l - 1) / 2) {
            let mut plus = [0u16; 16];
            let mut minus = [0u16; 16];
            for a in 0..l {
                let x = base[i] + a;
                let y = base[i] + (a + k) % l;
                plus[x] |= 1 << y;
                minus[y] |= 1 << x;
            }
            withins.push((plus, minus));
        }
    }
    // row_check[k] = Some(row) if R row `row` is fully determined once slots
    // 0..k are assigned (bundles are built sorted by (i, j))
    let mut row_check: Vec<Option<usize>> = vec![None; bundles.len() + 1];
    for i in 0..m {
        // last slot index touching row i as min-index
        let last = bundles.iter().rposition(|b| b.i == i);
        if let Some(l) = last {
            row_check[l + 1] = Some(i);
        }
    }
    let mut ctx = Ctx {
        n,
        m,
        row_check,
        withins,
        fix: 0,
        bundle_leaves: 0,
        balance_pass: 0,
    };
    let mut r = [[0i64; 8]; 8];
    let mut beats0 = [0u16; 16];
    let mut start = 0usize;
    let shard_lbl;
    if let (Some(s), Some(of)) = (
        env::args().nth(2).and_then(|x| x.parse::<usize>().ok()),
        env::args().nth(3).and_then(|x| x.parse::<usize>().ok()),
    ) {
        let mut k = 0usize;
        let mut p = 1usize;
        while p < of {
            p *= 3;
            k += 1;
        }
        assert!(p == of && s < of && k <= bundles.len(), "shard OF must be a power of 3");
        let mut sv = s;
        for t in 0..k {
            let digit = sv % 3;
            sv /= 3;
            let val: i32 = [0, 1, -1][digit];
            let bl = &bundles[t];
            if val == 1 {
                for x in 0..n {
                    beats0[x] |= bl.plus[x];
                }
                r[bl.i][bl.j] += bl.wij;
                r[bl.j][bl.i] -= bl.wji;
            } else if val == -1 {
                for x in 0..n {
                    beats0[x] |= bl.minus[x];
                }
                r[bl.i][bl.j] -= bl.wij;
                r[bl.j][bl.i] += bl.wji;
            }
        }
        start = k;
        shard_lbl = format!(" shard={}/{}", s, of);
    } else {
        shard_lbl = String::new();
    }
    bundle_dfs(start, &bundles, beats0, &mut r, &mut ctx);
    println!(
        "type={} n={}{} bundle_leaves={} balance_pass={} fix={}",
        arg, n, shard_lbl, ctx.bundle_leaves, ctx.balance_pass, ctx.fix
    );
}
