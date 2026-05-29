"""Modular decomposition of fair games (skew {-1,0,+1} 'win/lose/tie' structures).

A *module* is a set S of moves that every outside move relates to identically (same
sign in M for all of S). Modules generalize tie-twins -- a twin pair is exactly a
size-2 all-tie module. A game with no nontrivial module (only singletons and the
whole set) is **modular-prime**: a genuine irreducible building block. Every other
game is a prime/degenerate *quotient* with sub-games substituted into its moves
(G = H[M_1, ..., M_k]), and the number of extreme equilibria factorizes over the
decomposition tree:

    n_eq(node) = sum over vertices w of O(quotient) of  prod_{i in supp(w)} n_eq(child_i)

with a leaf contributing 1. (All-tie quotient -> sum of children; transitive 'order'
quotient -> the dominant child; prime quotient -> the general mixed sum.) This is the
exact generalization of the twin 2^k blow-up, and it shows the equilibrium-count
inflation is entirely a decomposition phenomenon.
"""

from itertools import combinations

import numpy as np

from .structure import canonical_key


def _cyclic(n, dists):
    M = np.zeros((n, n), dtype=np.int8)
    for i in range(n):
        for d in dists:
            j = (i + d) % n
            M[i, j], M[j, i] = 1, -1
    return M


def _cop():
    M = np.zeros((4, 4), dtype=np.int8)
    for w, l in [(3, 0), (0, 2), (0, 1), (1, 2), (2, 3)]:
        M[w, l], M[l, w] = 1, -1
    return M


# canonical_key materializes all n! relabelings, so it is only affordable for
# small modules. All-tie blocks (automorphism group S_k) would blow it up, so they
# are recognized structurally instead and never reach canonical_key.
_MAX_CANON = 8
_KNOWN = None


def _known_games():
    """Registry {canonical_key: name} of recognizable small games (built once)."""
    global _KNOWN
    if _KNOWN is None:
        reg = {}
        reg[canonical_key(_cyclic(3, [1]))] = "RPS"
        reg[canonical_key(_cyclic(5, [1, 2]))] = "RPSLS"
        reg[canonical_key(_cop())] = "cops"
        reg[canonical_key(_cyclic(7, [1, 2, 4]))] = "Paley7"
        for k in range(4, _MAX_CANON + 1):  # rings (ring3 == RPS, already named)
            reg.setdefault(canonical_key(_cyclic(k, [1])), f"ring{k}")
        _KNOWN = reg
    return _KNOWN


def named_subgame(M):
    """Short name for game M if it is a recognized one (RPS, cop, ring4, tie2, ...), else None.

    All-tie blocks are detected structurally (any size); other games are matched by
    canonical_key, which is only computed for small modules (canonical_key is O(n!)).
    """
    M = np.asarray(M)
    k = len(M)
    if k >= 2 and not M.any():  # all-tie block -- skip canonical_key (S_k symmetry)
        return f"tie{k}"
    if k > _MAX_CANON:
        return None
    return _known_games().get(canonical_key(M))


def is_module(M, S):
    """Is S subseteq V a module: does every move outside S see all of S identically?"""
    n = len(M)
    S = list(S)
    if len(S) <= 1 or len(S) == n:
        return True
    inside = set(S)
    for x in range(n):
        if x in inside:
            continue
        col = M[x, S[0]]
        if any(M[x, s] != col for s in S[1:]):
            return False
    return True


def is_prime(M):
    """True iff the game is modular-prime: n>=3 and no proper nontrivial module."""
    n = len(M)
    if n < 3:
        return False
    for size in range(2, n):
        for S in combinations(range(n), size):
            if is_module(M, S):
                return False
    return True


def _all_modules(M):
    n = len(M)
    return [frozenset(S) for size in range(1, n + 1) for S in combinations(range(n), size) if is_module(M, S)]


def _overlap(a, b):
    return bool(a & b) and not (a <= b) and not (b <= a)


def _quotient_type(Q):
    k = len(Q)
    off = Q[~np.eye(k, dtype=bool)]
    if np.all(off == 0):
        return "tie"
    if np.all(off != 0):  # tournament quotient; transitive => 'order'
        order = all(
            not (Q[i, j] == 1 and Q[j, l] == 1 and Q[i, l] != 1)
            for i in range(k)
            for j in range(k)
            for l in range(k)
        )
        if order:
            return "order"
    return "prime"


def modular_decomposition(M):
    """Strong-module decomposition tree.

    Returns nested dicts: {'members': frozenset, 'type': 'leaf'|'tie'|'order'|'prime',
    'children': [...]}. Internal nodes carry the type of their quotient (the relation
    among their child modules).
    """
    M = np.asarray(M)
    n = len(M)
    mods = _all_modules(M)
    strong = {A for A in mods if all(not _overlap(A, B) for B in mods)}

    def build(S):
        inside = [T for T in strong if T < S]
        children = [T for T in inside if not any(T < U for U in inside)]
        if not children:
            return {"members": S, "type": "leaf", "children": []}
        reps = [min(c) for c in children]
        Q = M[np.ix_(reps, reps)]
        return {
            "members": S,
            "type": _quotient_type(Q),
            "children": [build(c) for c in sorted(children, key=min)],
        }

    return build(frozenset(range(n)))


def neq_tree(M, node=None):
    """Number of extreme equilibria computed via the decomposition tree.

    Equals num_equilibria(M); provided as an independent check of the factorization
    law and as an O(small) route once the tree is known.
    """
    from .equilibrium import equilibrium_vertices

    M = np.asarray(M)
    if node is None:
        node = modular_decomposition(M)
    if node["type"] == "leaf":
        return 1
    cn = [neq_tree(M, c) for c in node["children"]]
    reps = [min(c["members"]) for c in node["children"]]
    Q = M[np.ix_(reps, reps)]
    total = 0
    for w in equilibrium_vertices(Q):
        prod = 1
        for i, wi in enumerate(w):
            if wi > 1e-9:
                prod *= cn[i]
        total += prod
    return total


def tree_summary(node, names=None):
    """One-line bracket string of the decomposition tree, e.g. prime[cop, cop, cop]."""

    def lab(members):
        if names is not None:
            return "/".join(names[i] for i in sorted(members))
        return ",".join(str(i) for i in sorted(members))

    def rec(nd):
        if nd["type"] == "leaf":
            return lab(nd["members"])
        return f"{nd['type']}[" + ", ".join(rec(c) for c in nd["children"]) + "]"

    return rec(node)
