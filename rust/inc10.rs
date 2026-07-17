// FUSED inclusive(even n) nullity-2 labeled counter -- the 6-hour-budget
// version of inc_count.rs. One pass over (n-2)-vertex GRANDPARENT classes:
//
//   grandparent M'' (nonsingular)  --cone DFS over c-->  semi-CM parent P
//   --child DFS over r-->  weighted leaves (1/z), no child ever materialized.
//
// The two structural wins over inc_count.rs:
//  * NO per-parent linear algebra: for P = M'' extended by c, the parent kernel
//    is v = (-sgn * adj'' c, |det''|) and the child particular solution is
//    w = (-adj'' r1 / det'', 0) (skewness kills the coupling term; the
//    consistency row IS the r . v = 0 condition). adj'' is computed ONCE per
//    grandparent; its entries are 7x7 sign-matrix minors <= ~907, so all DFS
//    rows, dots and endpoint fractions fit i64 with big margins.
//  * NO nauty in the hot path: a parent class is accepted from exactly one
//    (G, c) construction -- the one whose added vertex has the strictly
//    maximal vertex-signature among supp(v) -- whenever both G and P are
//    rigidity-certified (all vertex signatures distinct). Degenerate cases
//    (~10%) fall back to full canon + a per-grandparent dedup set.
//
// Anchors (must be EXACT before any n=10 run):
//   n=6 grandparents(4-vtx):  L2 = 126900
//   n=8 grandparents(6-vtx):  L2 = 45897886776
//
//   rustc -O rust/inc10.rs -o /tmp/inc10 -C link-args="shim.o -lnauty"
//   nauty-geng 6 | nauty-directg -o | /tmp/inc10 8
use std::collections::HashSet;
use std::env;
use std::io::{self, Read};
use std::os::raw::c_int;

mod common;

extern "C" {
    fn rps_canon(arc: *const u64, n: c_int, canong: *mut u64, lab: *mut c_int, orbits: *mut c_int);
    fn rps_autsize(arc: *const u64, n: c_int) -> f64;
}

static mut NODES: u64 = 0;
static mut CONE_LEAVES: u64 = 0;
static mut SIG_CALLS: u64 = 0;

use common::{adjugate_ff, connected_beats, factorial, lcm_to, vertex_sigs};

// child-extension DFS state sizes: eq row + up to 8 strict rows + 9 tracking
const MAXR: usize = 20;

#[allow(clippy::too_many_arguments)]
fn child_dfs(
    k: usize,
    p: usize,
    nstrict: usize,
    nl: usize,
    s: &mut [i32; MAXR],
    r: &mut [i32; 10],
    rows: &[[i32; 16]; MAXR],
    asum: &[[i32; 11]; MAXR],
    fplus: u16,
    fminus: u16,
    nod: u8,
    nid: u8,
    ctx: &mut LeafCtx,
) {
    // row 0: equality (v . r = 0); rows 1..1+nstrict: strictly negative
    if s[0].abs() > asum[0][k] {
        return;
    }
    for i in 1..(1 + nstrict) {
        if s[i] - asum[i][k] >= 0 {
            return;
        }
    }
    if k == p {
        leaf(ctx, s, r);
        return;
    }
    if nstrict != 0 {
        // degree-domination prune: acceptance needs the new vertex to be the
        // (od,id)-lex argmax over Z(v) u {new} in the child. Its od can reach
        // at most nod + remaining; each Z-vertex's od is bounded below by its
        // parent od (+1 if its already-set coordinate is +1). If some Z-vertex
        // wins in every completion, no leaf below can be accepted.
        let cap = nod + (p - k) as u8;
        if cap < ctx.zmaxpod {
            return;
        }
        // per-vertex codmin <= zpod+1 <= zmaxpod+1, so when cap exceeds that
        // no prune or od-tie can fire and the loop is pure overhead
        if cap <= ctx.zmaxpod + 1 {
        for t in 0..nstrict {
            let posz = ctx.zpos[t] as usize;
            let (czp, czm) = if posz < k {
                let rv = r[posz];
                ((rv > 0) as u8, (rv < 0) as u8)
            } else {
                (0u8, 1u8) // od-tie scenario forces all remaining to -1
            };
            let codmin = ctx.zpod[t] + czp;
            if cap < codmin || (cap == codmin && nid < ctx.zpid[t] + czm) {
                return;
            }
        }
        }
    }
    unsafe { NODES += 1; }
    if k + 1 == p {
        // last level: the equality row DETERMINES the final coordinate --
        // s[0] + val*rows[0][k] = 0 has at most one admissible val
        let c0 = rows[0][k];
        let fp = (fplus >> k) & 1 != 0;
        let fm = (fminus >> k) & 1 != 0;
        let val = if c0 != 0 {
            if s[0] % c0 != 0 {
                r[k] = 0;
                return;
            }
            let v = -s[0] / c0;
            if v < -1 || v > 1 || (fp && v != 1) || (fm && v != -1) {
                r[k] = 0;
                return;
            }
            v
        } else {
            if s[0] != 0 {
                r[k] = 0;
                return;
            }
            // free coordinate: all three values close the equality; check each
            for val in [0i32, -1, 1] {
                if (fp && val != 1) || (fm && val != -1) {
                    continue;
                }
                r[k] = val;
                if val != 0 {
                    for i in 0..nl {
                        s[i] += val * rows[i][k];
                    }
                }
                let mut ok = true;
                for i in 1..(1 + nstrict) {
                    if s[i] >= 0 {
                        ok = false;
                        break;
                    }
                }
                if ok {
                    leaf(ctx, s, r);
                }
                if val != 0 {
                    for i in 0..nl {
                        s[i] -= val * rows[i][k];
                    }
                }
            }
            r[k] = 0;
            return;
        };
        r[k] = val;
        if val != 0 {
            for i in 0..nl {
                s[i] += val * rows[i][k];
            }
        }
        let mut ok = true;
        for i in 1..(1 + nstrict) {
            if s[i] >= 0 {
                ok = false;
                break;
            }
        }
        if ok {
            leaf(ctx, s, r);
        }
        if val != 0 {
            for i in 0..nl {
                s[i] -= val * rows[i][k];
            }
        }
        r[k] = 0;
        return;
    }
    if fplus | fminus == 0 {
        // common case: no forced coordinates anywhere in this parent
        for val in [0i32, -1, 1] {
            r[k] = val;
            if val != 0 {
                for i in 0..nl {
                    s[i] += val * rows[i][k];
                }
            }
            child_dfs(k + 1, p, nstrict, nl, s, r, rows, asum, 0, 0, nod + (val < 0) as u8, nid + (val > 0) as u8, ctx);
            if val != 0 {
                for i in 0..nl {
                    s[i] -= val * rows[i][k];
                }
            }
        }
        r[k] = 0;
        return;
    }
    for val in [0i32, -1, 1] {
        if (fplus >> k) & 1 != 0 && val != 1 {
            continue; // paradox-forced: this vertex must beat the new one
        }
        if (fminus >> k) & 1 != 0 && val != -1 {
            continue; // paradox-forced: the new vertex must beat this one
        }
        r[k] = val;
        if val != 0 {
            for i in 0..nl {
                s[i] += val * rows[i][k];
            }
        }
        child_dfs(k + 1, p, nstrict, nl, s, r, rows, asum, fplus, fminus, nod + (val < 0) as u8, nid + (val > 0) as u8, ctx);
        if val != 0 {
            for i in 0..nl {
                s[i] -= val * rows[i][k];
            }
        }
    }
    r[k] = 0;
}

struct LeafCtx {
    p: usize,             // parent size = n-1
    nz: usize,            // |Z(v)|: 0 => CM parent, trivial half-weight accept
    zverts: [usize; 10],  // original indices of the Z(v) vertices
    wp: u128,             // (n-1)!/|Aut(P)|
    lcm: u128,
    nowin: u16,
    noloss: u16,
    parent_disconnected: bool,
    pbeats: [u16; 16],
    ord: [usize; 10],
    // parent degrees (child degrees = these + r contribution)
    pod: [u8; 16],
    pid: [u8; 16],
    // Z(v) data in DFS-ordered coordinates, for the mid-DFS domination prune
    zpos: [u8; 8],
    zpod: [u8; 8],
    zpid: [u8; 8],
    zmaxpod: u8,
    sum: u128,
    leaves: u64,
    accepted: u64,
}

fn leaf(ctx: &mut LeafCtx, _s: &[i32; MAXR], r: &[i32; 10]) {
    let p = ctx.p;
    let mut plus = 0u16; // vertices that beat the new vertex
    let mut minus = 0u16; // vertices the new vertex beats
    for i in 0..p {
        let c = ctx.ord[i];
        if r[i] > 0 {
            plus |= 1 << c;
        } else if r[i] < 0 {
            minus |= 1 << c;
        }
    }
    if minus == 0 || plus == 0 {
        return;
    }
    if ctx.nowin & !plus != 0 || ctx.noloss & !minus != 0 {
        return;
    }
    if ctx.parent_disconnected {
        let mut beats = ctx.pbeats;
        for i in 0..p {
            if r[i] > 0 {
                beats[ctx.ord[i]] |= 1 << p;
            }
        }
        beats[p] = minus;
        if !connected_beats(&beats, p + 1) {
            return;
        }
    }
    ctx.leaves += 1;
    if ctx.nz == 0 {
        // CM parent: Z1 = {new}, always the argmax -- weight 1/2
        ctx.accepted += 1;
        ctx.sum += ctx.wp * (ctx.lcm / 2);
        return;
    }
    // two-sided rule: accept iff the new vertex attains the max of an
    // iso-invariant signature over Z1 = Z(v) u {new}, weight 1/(2*T1).
    // Child degrees: parent degrees + r contribution; new vertex = counts.
    let nod = minus.count_ones() as u8; // new vertex out-degree
    let nid = plus.count_ones() as u8;
    // child (od,id) for any vertex i: pod/pid + (i beats new)/(new beats i)
    let cod = |i: usize| ctx.pod[i] + ((plus >> i) & 1) as u8;
    let cid = |i: usize| ctx.pid[i] + ((minus >> i) & 1) as u8;
    // signature digest of a child vertex (matches vertex_sigs semantics):
    // (od, id, out-neighbour digest, in-neighbour digest). Neighbour digests
    // need each neighbour's child (od,id): O(p) per competitor, competitors <= 4.
    let sig_of = |v: usize| -> u64 {
        let (ovd, ivd) = if v == p { (nod, nid) } else { (cod(v), cid(v)) };
        let vbeats: u16 = if v == p { minus } else { ctx.pbeats[v] | (((plus >> v) & 1) << p) };
        let vlosers: u16 = if v == p {
            plus
        } else {
            // who beats v in the child: parents who beat v, plus new if new beats v
            let mut l = 0u16;
            for u in 0..p {
                if ctx.pbeats[u] & (1 << v) != 0 {
                    l |= 1 << u;
                }
            }
            l | (((minus >> v) & 1) << p)
        };
        let (mut so, mut sq, mut si, mut sqi) = (0u32, 0u32, 0u32, 0u32);
        for u in 0..=p {
            let du = if u == p {
                ((nod as u32) << 5) | nid as u32
            } else {
                ((cod(u) as u32) << 5) | cid(u) as u32
            };
            if vbeats & (1 << u) != 0 {
                so += du;
                sq += du * du;
            }
            if vlosers & (1 << u) != 0 {
                si += du;
                sqi += du * du;
            }
        }
        ((ovd as u64) << 56)
            | ((ivd as u64) << 48)
            | ((so as u64 & 0xFFF) << 36)
            | ((sq as u64 & 0xFFF) << 24)
            | ((si as u64 & 0xFFF) << 12)
            | (sqi as u64 & 0xFFF)
    };
    // cheap first stage: (od,id) domination
    let newkey = ((nod as u16) << 8) | nid as u16;
    let mut maybe_tied = 0u16;
    for t in 0..ctx.nz {
        let zv = ctx.zverts[t];
        let k = ((cod(zv) as u16) << 8) | cid(zv) as u16;
        if k > newkey {
            return; // beaten on degrees: new not argmax
        }
        if k == newkey {
            maybe_tied |= 1 << zv;
        }
    }
    let mut t1 = 1u64;
    if maybe_tied != 0 {
        let ns = sig_of(p);
        let mut m = maybe_tied;
        while m != 0 {
            let zv = m.trailing_zeros() as usize;
            m &= m - 1;
            let zs = sig_of(zv);
            if zs > ns {
                return;
            }
            if zs == ns {
                t1 += 1;
            }
        }
    }
    ctx.leaves += 0;
    ctx.accepted += 1;
    ctx.sum += ctx.wp * (ctx.lcm / (2 * t1 as u128));
}

// cone DFS over c in {-1,0,1}^g: parent kernel (-sgn adj'' c) must be >= 0
#[allow(clippy::too_many_arguments)]
fn cone_dfs(
    k: usize,
    g: usize,
    s: &mut [i64; 8],
    c: &mut [i32; 8],
    rows: &[[i64; 8]; 8], // rows[i] = -sgn*adj''[i][..]; need s_i >= 0 at leaf
    asum: &[[i64; 9]; 8],
    out: &mut Vec<([i32; 8], [i64; 8])>,
) {
    for i in 0..g {
        if s[i] + asum[i][k] < 0 {
            return; // can no longer reach >= 0
        }
    }
    if k == g {
        out.push((*c, *s));
        return;
    }
    for val in [0i32, -1, 1] {
        c[k] = val;
        if val != 0 {
            let f = val as i64;
            for i in 0..g {
                s[i] += f * rows[i][k];
            }
        }
        cone_dfs(k + 1, g, s, c, rows, asum, out);
        if val != 0 {
            let f = val as i64;
            for i in 0..g {
                s[i] -= f * rows[i][k];
            }
        }
    }
    c[k] = 0;
}

fn main() {
    let n: usize = env::args().nth(1).and_then(|s| s.parse().ok()).expect("usage: inc10 n < grandparents");
    assert!(n % 2 == 0 && (6..=10).contains(&n));
    let g = n - 2; // grandparent size
    let p = n - 1; // parent size
    let gsq = g * g;
    let gbytes = (gsq + 5) / 6;
    let reclen = 2 + gbytes + 1;
    let lcm = lcm_to(2 * n as u64);
    let pfact = factorial(p as u64);

    let mut stdin = io::stdin().lock();
    let mut buf = vec![0u8; 1 << 20];
    let mut have = 0usize;
    let (mut gseen, mut gsing, mut parents_fast, mut parents_canon) = (0u64, 0u64, 0u64, 0u64);
    let mut cone: Vec<([i32; 8], [i64; 8])> = Vec::with_capacity(256);
    let mut ctx = LeafCtx {
        p,
        nz: 0,
        zverts: [0; 10],
        wp: 0,
        lcm,
        nowin: 0,
        noloss: 0,
        parent_disconnected: false,
        pbeats: [0; 16],
        ord: [0; 10],
        pod: [0; 16],
        pid: [0; 16],
        zpos: [0; 8],
        zpod: [0; 8],
        zpid: [0; 8],
        zmaxpod: 0,
        sum: 0,
        leaves: 0,
        accepted: 0,
    };
    let mut canon_seen: HashSet<[u64; 16]> = HashSet::new();

    loop {
        let got = stdin.read(&mut buf[have..]).unwrap();
        if got == 0 {
            break;
        }
        have += got;
        let nrec = have / reclen;
        for ri in 0..nrec {
            let rec = &buf[ri * reclen..(ri + 1) * reclen];
            assert!(rec[0] == b'&' && rec[1] as usize == 63 + g, "misaligned digraph6");
            gseen += 1;
            let mut gb = [0u16; 16];
            let payload = &rec[2..2 + gbytes];
            let mut kk = 0usize;
            'dec: for &byte in payload {
                let mut bits = ((byte - 63) as u32) << 26;
                for _ in 0..6 {
                    if bits & 0x8000_0000 != 0 {
                        gb[kk / g] |= 1 << (kk % g);
                    }
                    bits <<= 1;
                    kk += 1;
                    if kk == gsq {
                        break 'dec;
                    }
                }
            }
            let mut gm = [[0i128; 8]; 8];
            for i in 0..g {
                let mut w = gb[i];
                while w != 0 {
                    let j = w.trailing_zeros() as usize;
                    w &= w - 1;
                    gm[i][j] = 1;
                    gm[j][i] = -1;
                }
            }
            // entries are minors of a {-1,0,1} skew matrix: |adj| <= 907,
            // |det| <= 4096, so i64 holds them with room to spare
            let (adj128, det128) = match adjugate_ff(&gm, g, false) {
                Some(x) => x,
                None => {
                    gsing += 1;
                    continue;
                }
            };
            let det = det128 as i64;
            let mut adj = [[0i64; 8]; 8];
            for i in 0..g {
                for j in 0..g {
                    debug_assert!(adj128[i][j].abs() < (1 << 40));
                    adj[i][j] = adj128[i][j] as i64;
                }
            }
            let sgn: i64 = if det > 0 { 1 } else { -1 };
            let absd = det.abs();
            // grandparent rigidity certificate + Aut(G) handling via sigs
            let mut gsig = [0u64; 16];
            vertex_sigs(&gb, g, &mut gsig);
            let mut g_rigid = true;
            for i in 0..g {
                for j in (i + 1)..g {
                    if gsig[i] == gsig[j] {
                        g_rigid = false;
                    }
                }
            }
            // cone DFS: kernel entries kv_i = (-sgn * adj c)_i >= 0
            let mut crows = [[0i64; 8]; 8];
            for i in 0..g {
                for j in 0..g {
                    crows[i][j] = -sgn * adj[i][j];
                }
            }
            let mut casum = [[0i64; 9]; 8];
            for i in 0..g {
                for k in (0..g).rev() {
                    casum[i][k] = casum[i][k + 1] + crows[i][k].abs();
                }
            }
            cone.clear();
            let mut cs = [0i64; 8];
            let mut cc = [0i32; 8];
            cone_dfs(0, g, &mut cs, &mut cc, &crows, &casum, &mut cone);
            canon_seen.clear();
            // grandparent degrees, for the staged degree-only argmax reject
            let mut god = [0u8; 8];
            let mut gid = [0u8; 8];
            for i in 0..g {
                god[i] = gb[i].count_ones() as u8;
                let mut cint = 0u8;
                for u in 0..g {
                    if gb[u] & (1 << i) != 0 {
                        cint += 1;
                    }
                }
                gid[i] = cint;
            }
            unsafe { CONE_LEAVES += cone.len() as u64; }
            let mut rows_buf = [[0i32; 16]; MAXR];
            let mut prows_buf = [[0i32; 16]; MAXR];
            let mut asum_buf = [[0i32; 11]; MAXR];
            for (cv, kv) in cone.iter() {
                // parent P = G + vertex g with column c: c_i = M'[i][g]
                // kernel v = (kv_0..kv_{g-1}, |det|) all >= 0, last > 0
                let mut pb = gb;
                let mut newrow = 0u16;
                for i in 0..g {
                    if cv[i] > 0 {
                        pb[i] |= 1 << g; // i beats new parent vertex
                    } else if cv[i] < 0 {
                        newrow |= 1 << i;
                    }
                }
                pb[g] = newrow;
                // stage 0: degree-only reject. new vertex degrees from c; parent
                // vertex degrees = grandparent degrees + c contribution. If any
                // supp vertex lex-beats new on (od, id), new cannot be the
                // sig-argmax (degrees are the signature's leading key).
                let mut nod = 0u8;
                let mut nid = 0u8;
                for i in 0..g {
                    if cv[i] < 0 {
                        nod += 1;
                    } else if cv[i] > 0 {
                        nid += 1;
                    }
                }
                let mut beaten = false;
                for i in 0..g {
                    if kv[i] == 0 {
                        continue;
                    }
                    let oi = god[i] + (cv[i] > 0) as u8;
                    let ii = gid[i] + (cv[i] < 0) as u8;
                    if oi > nod || (oi == nod && ii > nid) {
                        beaten = true;
                        break;
                    }
                }
                if beaten {
                    continue;
                }
                unsafe { SIG_CALLS += 1; }
                // signatures of the parent, rigidity + argmax acceptance
                let mut psig = [0u64; 16];
                vertex_sigs(&pb, p, &mut psig);
                let mut p_rigid = true;
                for i in 0..p {
                    for j in (i + 1)..p {
                        if psig[i] == psig[j] {
                            p_rigid = false;
                        }
                    }
                }
                // supp of v (kernel), the valid-deletion set
                let mut supp_mask = 1u16 << g; // last vertex always in supp
                for i in 0..g {
                    if kv[i] != 0 {
                        supp_mask |= 1 << i;
                    }
                }
                let wp: u128;
                // acceptance rule is a function of P alone:
                //   S = sig-argmax set within supp(v);
                //   |S| = 1: accept iff new == that vertex;
                //   |S| > 1: accept iff new is in the canonical-max orbit in S.
                // Dedup of Aut(G)-equivalent c's (same P, same G, different c):
                // needed whenever G is non-rigid or |S| > 1 -> canon(P) set.
                let mut maxsig = 0u64;
                let mut m = supp_mask;
                while m != 0 {
                    let i = m.trailing_zeros() as usize;
                    m &= m - 1;
                    if psig[i] > maxsig {
                        maxsig = psig[i];
                    }
                }
                if psig[g] != maxsig {
                    continue; // new vertex not in the argmax set: never accepted
                }
                let mut ties = 0u32;
                let mut m = supp_mask;
                while m != 0 {
                    let i = m.trailing_zeros() as usize;
                    m &= m - 1;
                    if psig[i] == maxsig {
                        ties += 1;
                    }
                }
                if ties == 1 && g_rigid && p_rigid {
                    // unique argmax = new, G rigid (one c per class-construction),
                    // P rigid (|Aut| = 1): the pure fast path
                    wp = pfact;
                    parents_fast += 1;
                } else {
                    // canon path: resolves sig ties, dedups equivalent c's, and
                    // supplies |Aut(P)|
                    let mut arc64 = [0u64; 16];
                    for i in 0..p {
                        arc64[i] = pb[i] as u64;
                    }
                    let mut canong = [0u64; 16];
                    let mut lab = [0i32; 16];
                    let mut orbits = [0i32; 16];
                    unsafe {
                        rps_canon(arc64.as_ptr(), p as c_int, canong.as_mut_ptr(), lab.as_mut_ptr(), orbits.as_mut_ptr());
                    }
                    if ties > 1 {
                        // choose the canonical-max orbit within the tied set
                        let mut pos = [0i32; 16];
                        for (q, &vv) in lab.iter().enumerate().take(p) {
                            pos[vv as usize] = q as i32;
                        }
                        let mut best = usize::MAX;
                        let mut m = supp_mask;
                        while m != 0 {
                            let i = m.trailing_zeros() as usize;
                            m &= m - 1;
                            if psig[i] == maxsig && (best == usize::MAX || pos[i] > pos[best]) {
                                best = i;
                            }
                        }
                        if orbits[g] != orbits[best] {
                            continue;
                        }
                    }
                    let mut key = [0u64; 16];
                    key[..p].copy_from_slice(&canong[..p]);
                    if !canon_seen.insert(key) {
                        continue;
                    }
                    let aut = unsafe { rps_autsize(arc64.as_ptr(), p as c_int) } as u128;
                    wp = pfact / aut;
                    parents_canon += 1;
                }
                // ---- child DFS over r in {-1,0,1}^p ----
                // rows: 0 = equality v.r (v = (kv, |det|));
                //       1..1+nz: strict rows for i in Z(v): (sgn adj)[i].r1 < 0
                //       then tracking rows for supp(v) minus last vertex
                let rows = &mut rows_buf; // cols 0..=g rewritten; DFS reads only those
                for i in 0..g {
                    rows[0][i] = kv[i] as i32;
                }
                rows[0][g] = absd as i32;
                let mut nz = 0usize;
                for i in 0..g {
                    if kv[i] == 0 {
                        for j in 0..g {
                            rows[1 + nz][j] = (sgn * adj[i][j]) as i32;
                        }
                        nz += 1;
                    }
                }
                let nl = 1 + nz;
                // order columns by descending combined L1 of the pruning rows
                // so bounds tighten fastest (leaf semantics need the inverse
                // permutation for the r bitmasks)
                let mut ord = [0usize; 10];
                for (t, o) in ord[..p].iter_mut().enumerate() {
                    *o = t;
                }
                let mut l1 = [0i32; 10];
                for c in 0..p {
                    let mut sum = 0i32;
                    for i in 0..(1 + nz) {
                        sum += rows[i][c].abs();
                    }
                    l1[c] = sum;
                }
                ord[..p].sort_unstable_by(|&a, &b| l1[b].cmp(&l1[a]));
                let prows = &mut prows_buf;
                for i in 0..nl {
                    for c in 0..p {
                        prows[i][c] = rows[i][ord[c]];
                    }
                }
                let rows = &*prows;
                let asum = &mut asum_buf;
                for i in 0..(1 + nz) {
                    asum[i][p] = 0; // reused buffer: reset the suffix seed
                    for k in (0..p).rev() {
                        asum[i][k] = asum[i][k + 1] + rows[i][k].abs();
                    }
                }
                // paradox masks of the parent
                let mut inn = [0u16; 16];
                for i in 0..p {
                    let mut w = pb[i];
                    while w != 0 {
                        let j = w.trailing_zeros() as usize;
                        w &= w - 1;
                        inn[j] |= 1 << i;
                    }
                }
                let mut nowin = 0u16;
                let mut noloss = 0u16;
                for i in 0..p {
                    if pb[i] == 0 {
                        nowin |= 1 << i;
                    }
                    if inn[i] == 0 {
                        noloss |= 1 << i;
                    }
                }
                ctx.nz = nz;
                {
                    let mut t = 0usize;
                    for i in 0..g {
                        if kv[i] == 0 {
                            ctx.zverts[t] = i;
                            t += 1;
                        }
                    }
                }
                ctx.wp = wp;
                ctx.nowin = nowin;
                ctx.noloss = noloss;
                ctx.parent_disconnected = !connected_beats(&pb, p);
                ctx.pbeats = pb;
                ctx.ord = ord;
                {
                    // parent degrees for child-degree derivation at leaves
                    for i in 0..p {
                        ctx.pod[i] = pb[i].count_ones() as u8;
                        let mut c = 0u8;
                        for u in 0..p {
                            if pb[u] & (1 << i) != 0 {
                                c += 1;
                            }
                        }
                        ctx.pid[i] = c;
                    }
                }
                {
                    // Z(v) data in ordered coordinates for the mid-DFS prune
                    let mut inv = [0u8; 16];
                    for (k2, &o) in ord[..p].iter().enumerate() {
                        inv[o] = k2 as u8;
                    }
                    ctx.zmaxpod = 0;
                    for t in 0..nz {
                        let zv = ctx.zverts[t];
                        ctx.zpos[t] = inv[zv];
                        ctx.zpod[t] = ctx.pod[zv];
                        ctx.zpid[t] = ctx.pid[zv];
                        if ctx.pod[zv] > ctx.zmaxpod {
                            ctx.zmaxpod = ctx.pod[zv];
                        }
                    }
                }
                // paradox-forced coordinates (in ordered positions): a parent
                // vertex with no win must beat the new vertex; no loss => lose
                let mut fplus = 0u16;
                let mut fminus = 0u16;
                let mut impossible = false;
                for k in 0..p {
                    let c = ord[k];
                    if c < g {
                        let nw = nowin & (1 << c) != 0;
                        let nl2 = noloss & (1 << c) != 0;
                        if nw && nl2 {
                            impossible = true;
                            break;
                        }
                        if nw {
                            fplus |= 1 << k;
                        }
                        if nl2 {
                            fminus |= 1 << k;
                        }
                    } else {
                        // c == g is the parent's own added vertex; it has arcs
                        // from the cone construction -- masks cover it too
                        let nw = nowin & (1 << c) != 0;
                        let nl2 = noloss & (1 << c) != 0;
                        if nw && nl2 {
                            impossible = true;
                            break;
                        }
                        if nw {
                            fplus |= 1 << k;
                        }
                        if nl2 {
                            fminus |= 1 << k;
                        }
                    }
                }
                if impossible {
                    continue;
                }
                let mut s0 = [0i32; MAXR];
                let mut r0 = [0i32; 10];
                child_dfs(0, p, nz, 1 + nz, &mut s0, &mut r0, rows, asum, fplus, fminus, 0, 0, &mut ctx);
            }
        }
        let rem = have - nrec * reclen;
        buf.copy_within(nrec * reclen..have, 0);
        have = rem;
    }
    let total = (n as u128) * ctx.sum;
    assert!(total % lcm == 0, "1/z weights did not resolve to an integer");
    unsafe {
        eprintln!("PROFILE nodes={} cone_leaves={} sig_calls={}", NODES, CONE_LEAVES, SIG_CALLS);
    }
    println!(
        "n={}: grandparents={} (singular {}) parents_fast={} parents_canon={} leaves={} accepted={} L_nullity2_labeled={}",
        n,
        gseen,
        gsing,
        parents_fast,
        parents_canon,
        ctx.leaves,
        ctx.accepted,
        total / lcm
    );
}
