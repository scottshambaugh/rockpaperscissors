import sys, numpy as np
from itertools import combinations
from rpsfair import twin_free, is_prime

def tally(Ms):
    t = tf = pr = 0
    for M in Ms:
        t += 1
        tf += twin_free(M)
        pr += is_prime(M)
    return t, tf, pr

def from_code(line, n):
    M = np.zeros((n, n), np.int8)
    for c, (i, j) in zip(line.strip(), combinations(range(n), 2)):
        v = 1 if c == '1' else (-1 if c == '2' else 0)
        M[i, j] = v; M[j, i] = -v
    return M

def from_T(parts):
    nv, ne = int(parts[0]), int(parts[1])
    M = np.zeros((nv, nv), np.int8)
    for k in range(ne):
        a, b = int(parts[2 + 2*k]), int(parts[3 + 2*k])
        M[a, b] = 1; M[b, a] = -1
    return M

mode, n = sys.argv[1], int(sys.argv[2])
if mode == 'code':
    t, tf, pr = tally(from_code(l, n) for l in sys.stdin if l.strip())
elif mode == 'T':
    t, tf, pr = tally(from_T(l.split()) for l in sys.stdin if l.strip())
elif mode == 'balanced':
    from rpsfair.generate import search_balanced_gen
    t, tf, pr = tally(M for M, _ in search_balanced_gen(n))
print(f'{mode} n={n}: total={t} twin_free={tf} prime={pr}')
