"""Isomorph-free generation of games via canonical augmentation (McKay).

Generates exactly one representative per isomorphism class of skew {-1,0,+1}
games by adding a node at a time and keeping an augmentation only when the new
node is the canonically-distinguished one. The recursion is depth-first and
lazy, so memory stays flat -- unlike brute enumerate-then-dedup (which iterates
all 3^C(n,2) labelings) and unlike the row-by-row searches (which store a whole
level and blow up by n=8). Filtering for a property (paradoxical, connected,
inclusive, ...) is done at the target size only, since none of those are
monotone under adding a node.
"""

from itertools import product

import numpy as np
from tqdm import tqdm

# Number of isomorphism classes of n-node skew {-1,0,+1} games (so the inclusive
# generator can show a real progress total + ETA when augmenting the (n-1) classes).
_ISO_COUNTS = {1: 1, 2: 2, 3: 7, 4: 42, 5: 582, 6: 21480}


def _normalize(seq):
    rank = {v: r for r, v in enumerate(sorted(set(seq)))}
    return [rank[s] for s in seq]


def _refine(M, col):
    """Equitable color refinement from an initial coloring (1-WL, respecting col).

    A vertex's new color is its old color plus the sorted multiset of (relation,
    neighbor color); since the old color is the first key, colors only ever split,
    never merge -- the result is the coarsest equitable refinement of `col`.
    """
    n = len(M)
    col = list(col)
    while True:
        sig = [
            (col[i], tuple(sorted((int(M[i, j]), col[j]) for j in range(n) if j != i)))
            for i in range(n)
        ]
        new = _normalize(sig)
        if new == col:
            return col
        col = new


def _wl_colors(M):
    """Coarsest equitable coloring, seeded by the WTL profile (wins/ties/losses)."""
    n = len(M)
    seed = _normalize(
        [(int((M[i] == 1).sum()), int((M[i] == 0).sum()) - 1, int((M[i] == -1).sum())) for i in range(n)]
    )
    return _refine(M, seed)


def _nauty_search(M):
    """Refine-and-individualize search (nauty-style).

    Refine to an equitable coloring; while a non-singleton color class remains,
    branch by individualizing each of its vertices (split it off, re-refine) and
    recurse. Each leaf is a discrete coloring = a labeling. Returns
    (canonical_key, orderings) where canonical_key is the lex-min flattened matrix
    over the leaves and `orderings` are *all* leaf labelings achieving it (one per
    automorphism). Visits ~|Aut| labelings rather than the full n!, so it is fast
    even for highly symmetric (e.g. regular) games where plain refinement collapses
    the whole graph into one color class.
    """
    M = np.asarray(M)
    n = len(M)
    best = [None]
    orderings = []

    def rec(col):
        col = _refine(M, col)
        cells = {}
        for i, c in enumerate(col):
            cells.setdefault(c, []).append(i)
        target = next((c for c in sorted(cells) if len(cells[c]) > 1), None)
        if target is None:  # discrete -> a labeling (position i -> vertex order[i])
            order = sorted(range(n), key=lambda i: col[i])
            flat = tuple(int(x) for x in M[np.ix_(order, order)].reshape(-1))
            if best[0] is None or flat < best[0]:
                best[0] = flat
                orderings.clear()
                orderings.append(order)
            elif flat == best[0]:
                orderings.append(order)
            return
        for v in cells[target]:  # individualize each vertex of the target cell
            rec(_normalize([(col[i], 0 if i == v else 1) for i in range(n)]))

    rec(_wl_colors(M))
    return best[0], orderings


def canonical_key_nauty(M):
    """Canonical key (lex-min flattened matrix) via refine-and-individualize."""
    return _nauty_search(M)[0]


def canonical_key_fast(M):
    """Iso-invariant canonical key (refine-and-individualize); see canonical_key_nauty.

    Different value than structure.canonical_key (all-n! lex-min) but still equal
    iff isomorphic, so it deduplicates correctly with one key per class (bounded
    memory) and stays fast even for symmetric games -- unlike the orbit-hash dedup
    that OOMs, or plain within-class permutation that hits n! on regular games.
    """
    return canonical_key_nauty(M)


def _canon_and_autos(M):
    """(canonical-labeling perm, automorphism perms) via refine-and-individualize.

    The canonical labeling is any leaf achieving the canonical key; the full
    automorphism group is recovered from *all* such leaves: if labelings o and o0
    both yield the canonical matrix then a = o o0^{-1} (a[o0[i]] = o[i]) fixes M.
    Fast even on symmetric games, unlike the n! within-class sweep.
    """
    n = len(M)
    _, orderings = _nauty_search(M)
    o0 = orderings[0]
    cperm = np.array(o0, dtype=np.int64)
    autos = np.empty((len(orderings), n), dtype=np.int64)
    for k, o in enumerate(orderings):
        for i in range(n):
            autos[k, o0[i]] = o[i]
    return cperm, autos


def _orbit_label(n, autos):
    """Union-find orbit id per node under the automorphism perms."""
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a in autos:
        for i in range(n):
            ri, rj = find(i), find(int(a[i]))
            if ri != rj:
                parent[ri] = rj
    return [find(i) for i in range(n)]


def _aug_orbit_reps(autos, n1, values=(-1, 0, 1)):
    """One augmentation vector per Aut(Y)-orbit of values^(n1)."""
    alist = [a.tolist() for a in autos]
    seen, reps = set(), []
    for tup in product(values, repeat=n1):
        v = np.array(tup, dtype=np.int8)
        key = min(tuple(int(x) for x in v[a]) for a in alist)
        if key not in seen:
            seen.add(key)
            reps.append(v)
    return reps


def _augment(Y, v):
    n = len(Y) + 1
    Z = np.zeros((n, n), dtype=np.int8)
    Z[: n - 1, : n - 1] = Y
    Z[n - 1, : n - 1] = v
    Z[: n - 1, n - 1] = -v
    return Z


def _accept(Z):
    """Canonical-augmentation test: the added node n-1 is the canonical last node."""
    n = len(Z)
    cperm, autos = _canon_and_autos(Z)
    c = int(cperm[n - 1])  # node placed last by the canonical labeling
    if c == n - 1:
        return True
    lab = _orbit_label(n, autos)
    return lab[n - 1] == lab[c]


def generate_up_to_iso(n, prune=None):
    """Yield one matrix per isomorphism class of n-node skew {-1,0,+1} games.

    Optional `prune(Z, k)` (k = current size) returns False to discard a partial
    game that cannot reach the target. It MUST be closed under canonical parent --
    if Z passes at size k then Z minus its canonical node passes at size k-1 -- so
    that no reachable class is lost. Used to grow only constraint-feasible classes
    (e.g. balance-feasible), keeping the generation isomorph-free with no duplicate
    leaves.
    """
    if n <= 1:
        yield np.zeros((1, 1), dtype=np.int8)
        return
    for Y in generate_up_to_iso(n - 1, prune):
        _, autos_y = _canon_and_autos(Y)
        for v in _aug_orbit_reps(autos_y, n - 1):
            Z = _augment(Y, v)
            if prune is not None and not prune(Z, n):
                continue
            if _accept(Z):
                yield Z


def count_iso_classes(n):
    """Number of isomorphism classes of n-node skew {-1,0,+1} games (closed form).

    Burnside count of `3^C(n,2)` labelings modulo the S_n relabeling action: for a
    permutation sigma the number of fixed games is `3^g`, where g is the number of
    sigma-orbits on unordered pairs that carry a consistent orientation. An orbit
    that maps some edge to its own reverse is forced to a tie (factor 1, not 3); a
    signed union-find over ordered pairs detects exactly those. Summing over cycle
    types (weighted by class size) is instant for any n -- far past where
    `generate_up_to_iso` can enumerate. Matches 1, 2, 7, 42, 582, 21480 for n <= 6;
    gives 2_142_288 / 575_016_219 / 415_939_243_032 at n = 7 / 8 / 9.
    """
    from collections import Counter
    from math import factorial

    def partitions(m, hi):
        if m == 0:
            yield ()
            return
        for k in range(min(m, hi), 0, -1):
            for rest in partitions(m - k, k):
                yield (k, *rest)

    def consistent_orbits(part):
        sigma, base = [], 0
        for length in part:  # representative permutation with these cycle lengths
            sigma += [base + (i + 1) % length for i in range(length)]
            base += length
        visited, good = set(), 0
        for i in range(n):
            for j in range(n):
                if i == j or (i, j) in visited:
                    continue
                local = {(i, j): 1}
                visited.add((i, j))
                stack, ok = [((i, j), 1)], True
                while stack:
                    (a, b), s = stack.pop()
                    for (c, d), ns in (((b, a), -s), ((sigma[a], sigma[b]), s)):
                        if (c, d) in local:
                            ok = ok and local[(c, d)] == ns
                        else:
                            local[(c, d)] = ns
                            visited.add((c, d))
                            stack.append(((c, d), ns))
                good += ok
        return good

    total = 0
    for part in partitions(n, n):
        denom = 1
        for k, mult in Counter(part).items():
            denom *= (k**mult) * factorial(mult)
        total += (factorial(n) // denom) * 3 ** consistent_orbits(part)
    return total // factorial(n)


def generate_tournaments(n):
    """Yield one matrix per isomorphism class of n-node tournaments.

    Same canonical augmentation as `generate_up_to_iso`, but the new node's
    edges are all decisive ({-1,+1}, no ties), so the leaf count is the number
    of tournaments up to iso (A000568: ... 456, 6880, 191536 at n=7, 8, 9).
    Tie-freeness is closed under removing the canonical node, so the generation
    stays isomorph-free.
    """
    if n <= 1:
        yield np.zeros((1, 1), dtype=np.int8)
        return
    for Y in generate_tournaments(n - 1):
        _, autos = _canon_and_autos(Y)
        for v in _aug_orbit_reps(autos, n - 1, values=(-1, 1)):
            Z = _augment(Y, v)
            if _accept(Z):
                yield Z


def search_balanced_gen(n):
    """Balanced games via isomorph-free canonical augmentation.

    Grows only *balance-feasible* classes: a partial k-node game is kept only if
    every |row sum| <= n-k (each of the n-k future entries can shift a sum by at
    most 1). That bound is closed under canonical parent, so generation stays
    isomorph-free; at the leaf (k=n) it forces row sums == 0, i.e. exactly the
    balanced classes -- no millions of duplicate labelings, bounded memory.
    """
    from .cache import cached
    from .structure import connected, paradoxical

    def prune(Z, k):
        return bool((np.abs(Z.sum(axis=1)) <= (n - k)).all())

    def go():
        uniform = np.ones(n) / n
        out = []
        bar = tqdm(generate_up_to_iso(n, prune), desc=f"balanced n={n}", unit="cls", leave=False)
        for M in bar:  # leaf prune already forces row sums == 0
            if paradoxical(M) and connected(M):
                out.append((M.copy(), uniform))
                bar.set_postfix_str(f"{len(out)} kept")
        bar.close()
        return out

    return cached(f"balanced_n{n}", go)


from .structure import paradoxical_batch as _paradoxical_batch  # noqa: E402


def _connected_batch(Zs):
    """Bool mask: the decisive-edge graph (ignoring ties) is connected, batched."""
    n = Zs.shape[1]
    reach = (Zs != 0) | np.eye(n, dtype=bool)
    r = reach.astype(np.int16)
    for _ in range(max(1, (n - 1).bit_length())):  # transitive closure by squaring
        r = (np.matmul(r, r) > 0).astype(np.int16)
    return (r > 0).all(axis=(1, 2))


def _has_fully_mixed_batch(Zs, tol=1e-8):
    """Bool mask: a fully-mixed (strictly positive kernel) equilibrium exists.

    Batched SVD handles the common nullity-1 case (just sign-check the lone kernel
    vector); the rare nullity>=2 case falls back to the LP existence test.
    """
    from .equilibrium import has_fully_mixed

    A = Zs.astype(float)
    _, S, Vt = np.linalg.svd(A)
    thr = tol * np.maximum(S[:, 0], 1.0)
    nul = (thr[:, None] > S).sum(axis=1)
    out = np.zeros(len(Zs), dtype=bool)
    for b in range(len(Zs)):
        if nul[b] == 1:
            v = Vt[b, -1]
            out[b] = bool((v > tol).all() or (v < -tol).all())
        elif nul[b] >= 2:
            out[b] = has_fully_mixed(Zs[b])[0]
    return out


def search_inclusive_gen(n):
    """Inclusive games at size n via canonical-augmentation generation + filter.

    Scales past the 3^C(n,2) brute-force wall. We generate every (n-1)-node class
    once (streaming canonical augmentation, flat memory), then augment each by a
    node and keep the paradoxical, connected, fully-mixed results -- the
    *expensive* canonical key is computed only on those survivors (to dedup the
    few hundred inclusive games), not on the millions of dead leaves. Augmentation
    vectors with no win or no loss for the new node are skipped up front, since
    such a game cannot be paradoxical. Returns [(M, maxmin_equilibrium)].
    """
    from .cache import cached
    from .equilibrium import has_fully_mixed, maxmin_equilibrium
    from .structure import canonical_key, connected, paradoxical

    def go():
        if n <= 2:
            return [(M, maxmin_equilibrium(M)) for M in generate_up_to_iso(n)
                    if paradoxical(M) and connected(M) and has_fully_mixed(M)[0]]
        # all new-node vectors that could be paradoxical (need a win and a loss)
        allv = np.array(
            [t for t in product((-1, 0, 1), repeat=n - 1) if 1 in t and -1 in t], dtype=np.int8
        )
        B = len(allv)
        seen, out = set(), []
        bar = tqdm(
            generate_up_to_iso(n - 1), total=_ISO_COUNTS.get(n - 1),
            desc=f"inclusive n={n}", unit="cls", leave=False,
        )
        for Y in bar:
            bar.set_postfix_str(f"{len(out)} found")
            # build every augmentation of Y at once; isomorphic duplicates from one
            # parent are fine -- the canonical-key dedup of survivors removes them.
            Zs = np.zeros((B, n, n), dtype=np.int8)
            Zs[:, : n - 1, : n - 1] = Y
            Zs[:, n - 1, : n - 1] = allv
            Zs[:, : n - 1, n - 1] = -allv
            Zs = Zs[_paradoxical_batch(Zs)]
            if len(Zs):
                Zs = Zs[_connected_batch(Zs)]
            if len(Zs):
                Zs = Zs[_has_fully_mixed_batch(Zs)]
            for Z in Zs:
                key = canonical_key(Z)
                if key not in seen:
                    seen.add(key)
                    out.append((Z, maxmin_equilibrium(Z)))
        bar.close()
        return out

    return cached(f"inclusive_n{n}", go)


# A000568: number of tournaments on n nodes up to iso -- the progress total for
# the two-paradox census, which streams every tournament class exactly once.
_TOURNAMENT_COUNTS = {1: 1, 2: 1, 3: 2, 4: 4, 5: 12, 6: 56, 7: 456, 8: 6880, 9: 191536}


def search_two_paradox(n):
    """Two-paradox (P2 / Erdos-Schutte) games at size n, up to isomorphism.

    P2 -- every *pair* of strategies has a common strict beater -- forces a
    tie-free tournament (a tie offers no beater), so the authoritative count is
    a filter over *all* tournament classes, enumerated once each by canonical
    augmentation (`generate_tournaments`). Crucially this is NOT a refinement of
    the regular/balanced/inclusive fairness ladder: the n=8 two-paradox games
    are none of those (no even-n tournament is regular or balanced, and they
    lack a fully-mixed equilibrium), so post-filtering a fairness tier
    undercounts -- e.g. filtering `search_regular` reports 0 at n=8 and 5 at
    n=9 versus the true 4 and 221. Returns [(M, maxmin_equilibrium)].
    """
    from .cache import cached
    from .equilibrium import maxmin_equilibrium
    from .structure import connected, k_paradoxical, paradoxical

    def go():
        out = []
        bar = tqdm(generate_tournaments(n), total=_TOURNAMENT_COUNTS.get(n),
                   desc=f"two-paradox n={n}", unit="cls", leave=False)
        for M in bar:
            if k_paradoxical(M, 2) and paradoxical(M) and connected(M):
                out.append((M.copy(), maxmin_equilibrium(M)))
                bar.set_postfix_str(f"{len(out)} found")
        bar.close()
        return out

    return cached(f"two_paradox_n{n}", go)
