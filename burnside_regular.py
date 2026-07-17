# regular(n) WITHOUT enumeration: iso-class counts of connected d-regular
# oriented games summed over d, via Burnside over cycle types + an
# exchangeable-state DP for the sigma-invariant labeled counts, then an
# inverse Euler transform for connectivity at the iso level.
#
# fix(sigma) structure for an oriented graph invariant under sigma:
#  - vertices in one sigma-cycle share (out, in);
#  - within a cycle of length L, chord-orbits s = 1..floor((L-1)/2) each admit
#    tie / forward / backward, adding (1,1) per vertex when oriented; the
#    antipodal orbit (s = L/2, L even) is FORCED to tie (orienting it would
#    require both directions of the same pair);
#  - between cycles of lengths a, b: g = gcd(a,b) orbits; orienting one sends
#    lcm/a = b/g outs per a-side vertex and lcm/b = a/g ins per b-side vertex
#    (or the reverse); ties add nothing;
#  - fixed points interact with a cycle of length L via one orbit: the fixed
#    vertex gains L ins (or L outs), each cycle vertex gains 1 out (or 1 in);
#  - fixed-fixed pairs are ordinary pairs: the stage-1 exchangeable DP.
# Every vertex must end at out = in = d.
from functools import lru_cache
from math import comb, gcd, factorial
import sys

sys.setrecursionlimit(100000)


def make_f(d):
    """Stage-1 DP: f(state) = # completions of the pairs among 'remaining'
    vertices, state = sorted tuple of their (o, i) tallies."""

    @lru_cache(maxsize=None)
    def f(state):
        if not state:
            return 1
        u = state[0]
        rest = list(state[1:])
        a, b = d - u[0], d - u[1]
        if a < 0 or b < 0:
            return 0
        if not rest:
            return 1 if (a == 0 and b == 0) else 0
        classes = {}
        for t in rest:
            classes[t] = classes.get(t, 0) + 1
        cls = sorted(classes)
        cnt = [classes[c] for c in cls]
        k = len(cls)
        total = 0

        def rec(ci, ra, rb, acc, newrest):
            nonlocal total
            if ci == k:
                if ra == 0 and rb == 0:
                    total += acc * f(tuple(sorted(newrest)))
                return
            c = cls[ci]
            n_c = cnt[ci]
            for x in range(min(ra, n_c) + 1):
                for y in range(min(rb, n_c - x) + 1):
                    ways = comb(n_c, x) * comb(n_c - x, y)
                    nr = (newrest + [(c[0], c[1] + 1)] * x
                          + [(c[0] + 1, c[1])] * y + [c] * (n_c - x - y))
                    rec(ci + 1, ra - x, rb - y, acc * ways, nr)

        rec(0, a, b, 1, [])
        return total

    return f


def fix_count(cycles, nfixed, d, f):
    """Clean implementation. Items: nontrivial cycles (list of lengths) and
    nfixed fixed points. Process cycles one at a time; a cycle decides its
    within-chords, its orbits to every LATER cycle, and its orbits to fixed
    points, then must close at (d, d). Fixed-fixed pairs close with f()."""

    @lru_cache(maxsize=None)
    def g(cstate, fstate):
        # cstate: sorted tuple of (len, o, i) of unprocessed cycles
        # fstate: sorted tuple of fixed-point (o, i)
        if not cstate:
            return f(fstate)
        L, o, i = cstate[0]
        later = tuple(cstate[1:])
        total = 0
        slots = (L - 1) // 2
        for t in range(slots + 1):
            o1, i1 = o + t, i + t
            if o1 > d or i1 > d:
                continue
            total += comb(slots, t) * (2 ** t) * walk(L, o1, i1, later, 0, fstate)
        return total

    @lru_cache(maxsize=None)
    def walk(L, o, i, later, idx, fstate):
        if idx == len(later):
            return close(L, o, i, later, fstate)
        L2, o2, i2 = later[idx]
        g2 = gcd(L, L2)
        up_o = L2 // g2   # out gain per current-cycle vertex per fw orbit
        up2 = L // g2     # gain per later-cycle vertex per orbit
        total = 0
        for fw in range(g2 + 1):
            oo = o + fw * up_o
            if oo > d:
                break
            ii2f = i2 + fw * up2
            if ii2f > d:
                break
            for bw in range(g2 - fw + 1):
                ii = i + bw * up_o
                if ii > d:
                    break
                oo2 = o2 + bw * up2
                if oo2 > d:
                    break
                ways = comb(g2, fw) * comb(g2 - fw, bw)
                nl = later[:idx] + ((L2, oo2, ii2f),) + later[idx + 1:]
                total += ways * walk(L, oo, ii, nl, idx + 1, fstate)
        return total

    @lru_cache(maxsize=None)
    def close(L, o, i, later, fstate):
        # distribute the cycle's remaining (d-o) outs / (d-i) ins over fixed
        # points, then recurse into g() for the remaining cycles.
        need_o, need_i = d - o, d - i
        if need_o < 0 or need_i < 0:
            return 0
        classes = {}
        for t in fstate:
            classes[t] = classes.get(t, 0) + 1
        cls = sorted(classes)
        cnt = [classes[c] for c in cls]
        k = len(cls)
        total = 0

        def rec(ci, ra, rb, acc, newf):
            nonlocal total
            if ci == k:
                if ra == 0 and rb == 0:
                    total += acc * g(tuple(sorted(later)), tuple(sorted(newf)))
                return
            c = cls[ci]
            n_c = cnt[ci]
            xmax = min(ra, n_c) if c[1] + L <= d else 0
            for x in range(xmax + 1):
                ymax = min(rb, n_c - x) if c[0] + L <= d else 0
                for y in range(ymax + 1):
                    ways = comb(n_c, x) * comb(n_c - x, y)
                    nf = (newf + [(c[0], c[1] + L)] * x
                          + [(c[0] + L, c[1])] * y + [c] * (n_c - x - y))
                    rec(ci + 1, ra - x, rb - y, acc * ways, nf)

        rec(0, need_o, need_i, 1, [])
        return total

    cstate = tuple(sorted(((L, 0, 0) for L in cycles)))
    fstate = tuple([(0, 0)] * nfixed)
    return g(cstate, fstate)


def partitions(n, mx=None):
    if mx is None:
        mx = n
    if n == 0:
        yield ()
        return
    for k in range(min(n, mx), 0, -1):
        for rest in partitions(n - k, k):
            yield (k,) + rest


def z_lambda(lam):
    from collections import Counter
    z = 1
    for l, a in Counter(lam).items():
        z *= (l ** a) * factorial(a)
    return z


def iso_all(m, d, fcache={}):
    """Iso classes of (possibly disconnected) d-regular oriented graphs on m
    vertices: Burnside over cycle types."""
    if d < 0:
        return 0
    key = d
    if key not in fcache:
        fcache[key] = make_f(d)
    f = fcache[key]
    total = 0
    for lam in partitions(m):
        cycles = [l for l in lam if l > 1]
        nfixed = sum(1 for l in lam if l == 1)
        fx = fix_count(tuple(cycles), nfixed, d, f)
        total += fx * (factorial(m) // z_lambda(lam))
    assert total % factorial(m) == 0
    return total // factorial(m)


def connected_iso(N, d):
    """Inverse Euler transform: c(n) from a(n) = iso_all(n, d), n = 1..N."""
    a = [0] * (N + 1)
    for m in range(1, N + 1):
        a[m] = iso_all(m, d)
    # Euler transform inverse via the standard b-sequence:
    # b(n) = n*a(n) - sum_{k=1}^{n-1} b(k) a(n-k);  b(n) = sum_{dd|n} dd*c(dd)
    b = [0] * (N + 1)
    c = [0] * (N + 1)
    for n in range(1, N + 1):
        b[n] = n * a[n] - sum(b[k] * a[n - k] for k in range(1, n))
        s = sum(dd * c[dd] for dd in range(1, n) if n % dd == 0)
        assert (b[n] - s) % n == 0
        c[n] = (b[n] - s) // n
    return c


if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 9
    dmax = (N - 1) // 2
    grand = [0] * (N + 1)
    for d in range(1, dmax + 1):
        c = connected_iso(N, d)
        print(f"d={d}: connected iso by n: {c[1:]}", flush=True)
        for n in range(1, N + 1):
            grand[n] += c[n]
    print("regular(n) totals:", grand[1:])
    print("expect (n=3..11): 1 1 2 5 13 82 2016 154831 21171976")
