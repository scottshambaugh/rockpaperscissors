#!/usr/bin/env python3
"""
s2_iso.py -- two-paradox count up to isomorphism via Burnside, per conjugacy class.

The number of S2 ("two-paradox" / Erdos-Schutte) tournaments up to isomorphism is

    two_paradox(n) = (1/n!) sum_{sigma in S_n} fix(sigma)
                   = sum_{lambda |- n, all parts odd} S(lambda) / z_lambda ,

where S(lambda) = # S2 tournaments fixed by a permutation of cycle type lambda, and
z_lambda = prod_i (i^{m_i} m_i!).  Only all-odd cycle types contribute: every
automorphism of a tournament has odd order, so fix(sigma)=0 unless sigma has odd order.

This is an INDEPENDENT third method for the README two-paradox column (the other two
are nauty `nauty-gentourng | s2_filter` and the from-scratch `s2_count`), and the only
one that has any chance at n>=13 where isomorph-free generation of all tournaments is
hopeless -- because Burnside replaces "generate 10^30 tournaments" by "brute a handful
of highly-symmetric fixed-point counts S(lambda)".

Per-class ingredients, all reproducible:
  * identity lambda = [1^n]:  S = L(n) = # labeled S2 tournaments (the gamma=3 column of
    the README domination table).  For small n it is just `s2_slam n 1 1 ... 1`; for large
    n supply it with --Ln (from that table / the poly pipeline).
  * single 5-cycle classes [5,1^m]:  S([5,1^m]) = 2*S([3,1^m])  (validated reduction,
    e.g. S([5,1^6])=1280=2*640), which lowers their edge-orbit count into brute range.
  * every other class:  brute with `rust/s2_slam` (nested-prune S(lambda) counter).

A class whose edge-orbit count E exceeds the brute ceiling is printed as a BLOCKER and
NOT guessed; the partial Burnside sum and the exact missing terms are reported instead.

Usage:
    rustc -O -C target-cpu=native rust/s2_slam.rs -o /tmp/s2slam
    python3 rust/s2_iso.py 10 --slam /tmp/s2slam                 # reproduces 29816
    python3 rust/s2_iso.py 13 --slam /tmp/s2slam --Ln 10059739307720354796544 --ceil 42
"""
import sys, subprocess, argparse
from math import factorial
from collections import Counter

def perm_of(parts, n):
    p = [0]*n; v = 0
    for l in parts:
        for i in range(l): p[v+i] = v + (i+1) % l
        v += l
    return p

def num_orbits(parts, n):
    perm = perm_of(parts, n); seen = set(); o = 0
    for a in range(n):
        for b in range(a+1, n):
            if (a, b) in seen: continue
            o += 1; ca, cb = a, b
            while True:
                x, y = min(ca, cb), max(ca, cb)
                if (x, y) in seen: break
                seen.add((x, y)); ca, cb = perm[ca], perm[cb]
    return o

def zlam(parts):
    c = Counter(parts); z = 1
    for k, m in c.items(): z *= (k**m) * factorial(m)
    return z

def odd_partitions(n):
    out = []
    def rec(rem, mx, cur):
        if rem == 0: out.append(cur[:]); return
        p = min(mx, rem)
        while p >= 1:
            if p % 2 == 1:
                cur.append(p); rec(rem-p, p, cur); cur.pop()
            p -= 1
    rec(n, n, []); return out

def slam(slampath, n, parts):
    r = subprocess.run([slampath, str(n)] + [str(x) for x in parts],
                       capture_output=True, text=True)
    if r.returncode != 0 or not r.stdout.strip():
        raise RuntimeError(f"s2_slam failed for {parts}: {r.stderr[:200]}")
    return int(r.stdout.strip())

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("n", type=int)
    ap.add_argument("--slam", default="/tmp/s2slam", help="path to compiled rust/s2_slam")
    ap.add_argument("--Ln", type=int, default=None, help="L(n)=# labeled S2 tournaments (identity term)")
    ap.add_argument("--ceil", type=int, default=40, help="max edge-orbits E to brute directly")
    ap.add_argument("--cache", default=None, help="optional json cache of {lambda_str: S}")
    args = ap.parse_args()
    n = args.n; nf = factorial(n)

    cache = {}
    if args.cache:
        import json, os
        if os.path.exists(args.cache): cache = json.load(open(args.cache))

    def get_S(parts):
        key = ",".join(map(str, parts))
        if key in cache: return cache[key], "cache"
        if parts == [1]*n:
            if args.Ln is not None: return args.Ln, "Ln(identity)"
            return slam(args.slam, n, parts), "brute(identity)"
        # single 5-cycle reduction: S([5,1^m]) = 2*S([3,1^m])
        c = Counter(parts)
        if c.get(5, 0) == 1 and set(parts) <= {5, 1}:
            m = c.get(1, 0)
            three = [3] + [1]*m
            if num_orbits(three, m+3) <= args.ceil:
                return 2*slam(args.slam, m+3, three), "2*S([3,1^%d])" % m
        if num_orbits(parts, n) <= args.ceil:
            return slam(args.slam, n, parts), "brute"
        return None, "BLOCKER(E=%d)" % num_orbits(parts, n)

    total = 0            # will hold n! * two_paradox  (== sum coeff*S)
    blockers = []
    print(f"# two_paradox({n}) via Burnside over odd conjugacy classes")
    for parts in odd_partitions(n):
        coeff = nf // zlam(parts)          # n!/z_lambda
        S, how = get_S(parts)
        E = num_orbits(parts, n)
        if S is None:
            blockers.append((E, coeff, parts));
            print(f"  lambda={parts} E={E} coeff={coeff}  S=?? [{how}]")
            continue
        total += coeff * S
        print(f"  lambda={parts} E={E} coeff={coeff}  S={S}  [{how}]")
    print("-"*60)
    if blockers:
        print(f"two_paradox({n}) NOT closed. Missing {len(blockers)} class(es):")
        for E, coeff, parts in sorted(blockers):
            print(f"   BLOCKER lambda={parts}  E={E}  coeff={coeff}")
        print(f"partial (n! * two_paradox, feasible+identity) = {total}")
    else:
        assert total % nf == 0, "Burnside sum not divisible by n! -- inputs inconsistent!"
        print(f"two_paradox({n}) = {total // nf}   (exact: total %% {n}! == 0)")

if __name__ == "__main__":
    main()
