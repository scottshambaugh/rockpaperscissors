"""Structural metrics on fair RPS games."""

from itertools import permutations

import networkx as nx
import numpy as np


def num_roles(M, xs, tol=0.01):
    """Distinct (profile, equilibrium-rate) classes across nodes.

    Two nodes share a role iff same (W, T, L) profile AND same equilibrium
    play rate (bucketed by `tol`). `num_roles == n` means every strategy is
    structurally unique; `1` means all strategies are interchangeable.
    """
    roles = set()
    for i, r in enumerate(M):
        prof = (int((r == 1).sum()), int((r == 0).sum()) - 1, int((r == -1).sum()))
        roles.add((prof, round(xs[i] / tol)))
    return len(roles)


def automorphisms(M):
    """All node relabelings that leave M unchanged (the automorphism group)."""
    n = len(M)
    return [p for p in permutations(range(n)) if np.array_equal(M[np.ix_(p, p)], M)]


def aut_size(M):
    """Size of the automorphism group |Aut(M)|."""
    return len(automorphisms(M))


def node_orbits(M, autos=None):
    """Partition of nodes into automorphism orbits (mutually swappable nodes).

    Two nodes share an orbit iff some automorphism maps one to the other, i.e.
    they can be swapped without changing the game. Pass `autos` (from
    `automorphisms`) to avoid recomputing the group.
    """
    n = len(M)
    if autos is None:
        autos = automorphisms(M)
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for p in autos:
        for i in range(n):
            ri, rj = find(i), find(p[i])
            if ri != rj:
                parent[ri] = rj
    groups = {}
    for i in range(n):
        groups.setdefault(find(i), []).append(i)
    return list(groups.values())


def num_orbits(M, autos=None):
    """Number of automorphism orbits = count of structurally distinct strategies.

    Parameterization-free replacement for `num_roles`: it depends only on the
    game's symmetry, not on which canonical equilibrium is chosen. `num_orbits
    == n` means no two strategies are interchangeable (every one is structurally
    unique); `== 1` means all are interchangeable (RPS, regular games).
    """
    return len(node_orbits(M, autos))


def num_cuts(M):
    """Articulation points in the decisive subgraph. >0 means the game decomposes."""
    n = len(M)
    G = nx.Graph()
    G.add_nodes_from(range(n))
    for i in range(n):
        for j in range(i + 1, n):
            if M[i, j] != 0:
                G.add_edge(i, j)
    return len(list(nx.articulation_points(G)))


def gini(xs):
    """Gini coefficient of the equilibrium distribution. 0 = uniform play."""
    n = len(xs)
    s = np.sort(xs)
    c = np.cumsum(s)
    return (n + 1 - 2 * np.sum(c) / c[-1]) / n


def tie_fraction(M):
    """Fraction of all n*n matchups that are ties, including the diagonal.

    Each move ties itself, so the n diagonal cells always count as ties and
    contribute 1/n; a tournament (no off-diagonal ties) has tie fraction 1/n.
    """
    n = len(M)
    return int((M == 0).sum()) / (n * n)
