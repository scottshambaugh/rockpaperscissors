#!/usr/bin/env python3
"""Burnside assembly for the inclusive census.

inclusive(n) = (1/n!) * sum over sigma in S_n of Fix_sigma(labeled inclusive)
             = (L_inc(n) + sum over nontrivial cycle types of #sigma * Fix) / n!

Fix strategies per cycle type (all validated against brute at n <= 8):
  * m = 1 or non-coprime lengths: brute quotient enumeration (rust/sigma_fix.rs,
    exact; shardable as `sfix TYPE s OF`).
  * m >= 2 with pairwise-coprime cycle lengths: Fix = 3^(sum floor((l-1)/2))
    * L_inc(m). (Balance forces every quotient row both-signed-or-zero; a zero
    row is an isolated cycle => disconnected for m >= 2; so the quotient must
    be a plain labeled inclusive m-game and all within-assignments are free.)
    Verified: (3,1^5)=3*L_inc(6), (4,1^4)=(3,2,1^3)=3*L_inc(5),
    (5,1^3)=9*L_inc(4), (5,2,1)=(6,1,1)=(4,3,1)=9*L_inc(3), and the blind
    n=8 gate prediction Fix_(2,1^6)=L_inc(7).
  * (2,2,1^{n-4}): marked-pair sweep formula 2*J1+J0+4*H1+2*K0b
    (rust/sigma_sweep.rs; anchored: reproduces brute 339138 at n=8).

usage: sigma_gate.py n --linc 3=2,4=42,5=978,... [--fix TYPE=VALUE ...] [--plan]
  --linc  labeled inclusive counts for the coprime rule
  --fix   precomputed Fix values for heavy types (from sharded sfix runs or
          the sigma_sweep formula); overrides any strategy
  --plan  print the commands needed for missing heavy types instead of failing
"""
import subprocess
import sys
from math import factorial, gcd

SFIX = "/tmp/claude-1000/sfix"
BRUTE_SLOT_LIMIT = 15  # 3^15 ~ 14M bundle leaves: inline-quick

def partitions(n, mx=None):
    if mx is None:
        mx = n
    if n == 0:
        yield []
        return
    for k in range(min(n, mx), 0, -1):
        for rest in partitions(n - k, k):
            yield [k] + rest

def nperm(lam, n):
    denom = 1
    for l in set(lam):
        a = lam.count(l)
        denom *= (l ** a) * factorial(a)
    return factorial(n) // denom

def coprime(lam):
    ls = [l for l in lam if l > 1]
    for i in range(len(ls)):
        for j in range(i + 1, len(ls)):
            if gcd(ls[i], ls[j]) > 1:
                return False
    return True

def bundle_slots(lam):
    s = 0
    for i in range(len(lam)):
        for j in range(i + 1, len(lam)):
            s += gcd(lam[i], lam[j])
    return s

def within_exp(lam):
    return sum((l - 1) // 2 for l in lam)

def main():
    n = int(sys.argv[1])
    linc = {}
    fixed = {}
    plan = "--plan" in sys.argv
    for i, a in enumerate(sys.argv):
        if a == "--linc":
            for kv in sys.argv[i + 1].split(","):
                k, v = kv.split("=")
                linc[int(k)] = int(v)
        if a == "--fix":
            k, v = sys.argv[i + 1].split("=")
            fixed[k] = int(v)
    corrections = 0
    missing = []
    for lam in partitions(n):
        if lam == [1] * n:
            continue
        arg = ",".join(map(str, lam))
        m = len(lam)
        cnt = nperm(lam, n)
        if arg in fixed:
            fix, how = fixed[arg], "given"
        elif m >= 2 and coprime(lam):
            if m not in linc:
                missing.append((arg, f"needs L_inc({m})"))
                continue
            fix, how = 3 ** within_exp(lam) * linc[m], f"coprime 3^{within_exp(lam)}*L_inc({m})"
        elif bundle_slots(lam) <= BRUTE_SLOT_LIMIT:
            out = subprocess.run([SFIX, arg], capture_output=True, text=True).stdout
            fix, how = int(out.split("fix=")[1]), "brute"
        else:
            missing.append((arg, f"heavy brute ({bundle_slots(lam)} bundle slots): "
                                 f"run sharded `sfix {arg} S OF` or sweep, pass --fix {arg}=V"))
            continue
        corrections += cnt * fix
        print(f"  type {arg:24s} #sigma={cnt:9d}  fix={fix}  [{how}]")
    if missing:
        print("MISSING heavy types:")
        for arg, why in missing:
            print(f"  {arg}: {why}")
        if not plan:
            sys.exit(1)
        return
    fact = factorial(n)
    print(f"n={n}: corrections={corrections}")
    if n in linc:
        total = linc[n] + corrections
        print(f"n={n}: identity={linc[n]} total={total}")
        if total % fact == 0:
            print(f"n={n}: inclusive({n}) = {total // fact}")
        else:
            print(f"n={n}: NOT DIVISIBLE by {n}! -- something is wrong")

if __name__ == "__main__":
    main()
