"""Tally total / twin-free / modular-prime over games piped in or generated.

Scratch helper for the balanced/inclusive narrowings; the heavy cases now run in
`rust/balanced.rs` / `rust/regular.rs`, but this stays handy as a pure-Python
cross-check. Modes: `code` (upper-triangular trit string per line), `T`
(watercluster2 edge list), `balanced` (generate via search_balanced_gen).
"""

import sys
from itertools import combinations

import numpy as np

from rpsfair import is_prime, twin_free


def tally(games):
    total = tf = pr = 0
    for m in games:
        total += 1
        tf += twin_free(m)
        pr += is_prime(m)
    return total, tf, pr


def from_code(line, n):
    m = np.zeros((n, n), np.int8)
    for c, (i, j) in zip(line.strip(), combinations(range(n), 2), strict=False):
        v = 1 if c == "1" else (-1 if c == "2" else 0)
        m[i, j] = v
        m[j, i] = -v
    return m


def from_t(parts):
    nv, ne = int(parts[0]), int(parts[1])
    m = np.zeros((nv, nv), np.int8)
    for k in range(ne):
        a, b = int(parts[2 + 2 * k]), int(parts[3 + 2 * k])
        m[a, b] = 1
        m[b, a] = -1
    return m


def main():
    mode, n = sys.argv[1], int(sys.argv[2])
    if mode == "code":
        t, tf, pr = tally(from_code(ln, n) for ln in sys.stdin if ln.strip())
    elif mode == "T":
        t, tf, pr = tally(from_t(ln.split()) for ln in sys.stdin if ln.strip())
    elif mode == "balanced":
        from rpsfair.generate import search_balanced_gen

        t, tf, pr = tally(m for m, _ in search_balanced_gen(n))
    else:
        raise SystemExit(f"unknown mode {mode!r}")
    print(f"{mode} n={n}: total={t} twin_free={tf} prime={pr}")


if __name__ == "__main__":
    main()
