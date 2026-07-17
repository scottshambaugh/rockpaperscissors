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
import pynauty
from tqdm import tqdm

from .structure import pynauty_graph

# Number of isomorphism classes of n-node skew {-1,0,+1} games (so the inclusive
# generator can show a real progress total + ETA when augmenting the (n-1) classes).
_ISO_COUNTS = {1: 1, 2: 2, 3: 7, 4: 42, 5: 582, 6: 21480}

# Known final counts per tier (OEIS), used purely as a tqdm `total` so a search can
# show a true % / ETA toward the answer it is about to (re)derive. A pruned
# canonical augmentation has no cheap a-priori count of its own work, but the number
# of *kept* games climbs monotonically toward these, which is a usable progress
# signal. A308239 (balanced); inclusive is this project's own census.
_BALANCED_COUNTS = {
    3: 1,
    4: 1,
    5: 4,
    6: 16,
    7: 175,
    8: 5274,
    9: 434017,
    10: 90658149,
    11: 48825116761,
    12: 68579602126387,
}
_INCLUSIVE_COUNTS = {3: 1, 4: 3, 5: 15, 6: 222, 7: 10525, 8: 1198013}


def canonical_key_fast(M):
    """Iso-invariant canonical key via nauty -- equal iff isomorphic, so it
    deduplicates with one key per class in bounded memory. Same value as
    structure.canonical_key (both nauty); kept as a separate name for the
    streaming searches that dedup with it.
    """
    return pynauty.certificate(pynauty_graph(M))


def _canon_and_autos(M):
    """(canonical-labeling perm, automorphism perms) via nauty.

    `pynauty.canon_label` gives the canonical order (position p -> vertex
    `cperm[p]`, matching the old refine-and-individualize convention), and the
    automorphism group is nauty's generating set expanded to the full group.
    """
    from .metrics import _expand_group

    n = len(M)
    g = pynauty_graph(M)
    cperm = np.array(pynauty.canon_label(g), dtype=np.int64)
    gens = pynauty.autgrp(g)[0]
    autos = np.array(_expand_group(gens, n), dtype=np.int64)
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


def _aug_orbit_reps(autos, n1, values=(-1, 0, 1), keep=None):
    """One augmentation vector per Aut(Y)-orbit of values^(n1).

    `keep(tup)` (optional) drops infeasible augmentation vectors up front -- before
    the orbit-dedup and the per-candidate canonical test downstream. This is only a
    valid speedup when the predicate is constant on each Aut(Y)-orbit (so a whole
    orbit is kept or dropped together, never split); the balance-feasibility prune
    qualifies because permuting `tup` by an automorphism yields an isomorphic
    augmentation with the same row-sum multiset. Building the int8 vector is deferred
    to surviving reps, and orbit keys use plain tuples (no per-tuple numpy alloc).
    """
    alist = [a.tolist() for a in autos]
    seen, reps = set(), []
    for tup in product(values, repeat=n1):
        if keep is not None and not keep(tup):
            continue
        key = min(tuple(tup[a[i]] for i in range(n1)) for a in alist)
        if key not in seen:
            seen.add(key)
            reps.append(np.array(tup, dtype=np.int8))
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


def _generate_balanced(n):
    """`generate_up_to_iso` specialized to the balance-feasible subtree.

    Output is *identical* to `generate_up_to_iso(n, prune)` with the balance prune
    `|row sum| <= n-k`, but the infeasible augmentation vectors are rejected up
    front (inside `_aug_orbit_reps` via `keep`) rather than built, orbit-deduped,
    augmented and canonical-tested only to be thrown away. The new node's row sum is
    `sum(v)` and each existing node `i` shifts to `ysum[i] - v[i]`, so feasibility is
    a few integer comparisons on plain ints -- no `np.sum` per candidate, which
    profiling showed was the dominant cost (numpy dispatch overhead on tiny arrays,
    paid billions of times). Parent row sums are computed once per parent.
    """

    def rec(k):
        if k <= 1:
            yield np.zeros((1, 1), dtype=np.int8)
            return
        lim = n - k
        for Y in rec(k - 1):
            ysum = Y.sum(axis=1).tolist()  # once per parent, not per candidate

            def keep(tup, ysum=ysum, lim=lim):
                s = 0
                for yi, x in zip(ysum, tup, strict=False):
                    if not -lim <= yi - x <= lim:  # existing node stays feasible
                        return False
                    s += x
                return -lim <= s <= lim  # new node feasible

            _, autos_y = _canon_and_autos(Y)
            for v in _aug_orbit_reps(autos_y, k - 1, keep=keep):
                Z = _augment(Y, v)
                if _accept(Z):
                    yield Z

    yield from rec(n)


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

    def go():
        uniform = np.ones(n) / n
        out = []
        # Total = the known balanced count (A308239): the pruned augmentation has no
        # cheap a-priori work estimate, but kept games climb toward this, so the bar
        # shows a true % / ETA. tqdm writes to stderr -- a run harness that wants the
        # bar must NOT send stderr to /dev/null.
        bar = tqdm(total=_BALANCED_COUNTS.get(n), desc=f"balanced n={n}", unit="game", leave=False)
        for M in _generate_balanced(n):  # leaves already have all row sums == 0
            if paradoxical(M) and connected(M):
                out.append((M.copy(), uniform))
                bar.update(1)
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
            return [
                (M, maxmin_equilibrium(M))
                for M in generate_up_to_iso(n)
                if paradoxical(M) and connected(M) and has_fully_mixed(M)[0]
            ]
        # all new-node vectors that could be paradoxical (need a win and a loss)
        allv = np.array(
            [t for t in product((-1, 0, 1), repeat=n - 1) if 1 in t and -1 in t], dtype=np.int8
        )
        B = len(allv)
        seen, out = set(), []
        bar = tqdm(
            generate_up_to_iso(n - 1),
            total=_ISO_COUNTS.get(n - 1),
            desc=f"inclusive n={n}",
            unit="cls",
            leave=False,
        )
        for Y in bar:
            bar.set_postfix_str(f"{len(out)}/{_INCLUSIVE_COUNTS.get(n, '?')} found")
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
        bar = tqdm(
            generate_tournaments(n),
            total=_TOURNAMENT_COUNTS.get(n),
            desc=f"two-paradox n={n}",
            unit="cls",
            leave=False,
        )
        for M in bar:
            if k_paradoxical(M, 2) and paradoxical(M) and connected(M):
                out.append((M.copy(), maxmin_equilibrium(M)))
                bar.set_postfix_str(f"{len(out)} found")
        bar.close()
        return out

    return cached(f"two_paradox_n{n}", go)
