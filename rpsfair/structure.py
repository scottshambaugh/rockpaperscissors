"""Structural predicates and isomorphism handling for skew-symmetric matrices."""

import hashlib
from itertools import permutations

import numpy as np
import pynauty

_PERMS_CACHE = {}


def pynauty_graph(M):
    """A skew {-1,0,+1} game as a pynauty directed graph: arc i->j iff i beats j.

    Ties are non-edges, so each game is an oriented graph and nauty's
    directed-graph isomorphism/automorphism machinery applies directly.
    """
    M = np.asarray(M)
    n = len(M)
    adj = {i: np.flatnonzero(M[i] == 1).tolist() for i in range(n)}
    return pynauty.Graph(n, directed=True, adjacency_dict=adj)


def _perms(n):
    if n not in _PERMS_CACHE:
        _PERMS_CACHE[n] = np.array(list(permutations(range(n))))
    return _PERMS_CACHE[n]


def profile_of(M):
    """Sorted multiset of (wins, ties, losses) per node."""
    return tuple(
        sorted((int((r == 1).sum()), int((r == 0).sum()) - 1, int((r == -1).sum())) for r in M)
    )


def paradoxical(M):
    """Every node has at least one win and one loss."""
    return bool((M == 1).any(axis=1).all() and (M == -1).any(axis=1).all())


def paradoxical_batch(Zs):
    """Bool mask over a (B, n, n) stack: each node has at least one win and loss."""
    return (Zs == 1).any(axis=2).all(axis=1) & (Zs == -1).any(axis=2).all(axis=1)


def regular(M):
    """Every node has the same (W, T, L) profile. Strictly stronger than balanced."""
    W = (M == 1).sum(axis=1)
    L = (M == -1).sum(axis=1)
    return bool((W[0] == W).all() and (L[0] == L).all())


def reduce_twins(M):
    """Merge tie-twin nodes into one, returning the twin-free (fully reduced) core.

    Two nodes are twins iff they tie and share the same parents (beaters)
    and children (beaten) -- equivalently, identical rows in M. A twin pair
    is a duplicated strategy: redundant. By skew-symmetry identical rows
    imply identical columns, so merging never creates new twins; the core
    is simply the distinct rows. Returns (core_M, multiplicity) where
    multiplicity[k] counts original nodes collapsed into core node k.
    The merged node's equilibrium mass equals the sum of its members'.
    """
    M = np.asarray(M)
    seen = {}
    reps = []
    mult = []
    for i in range(len(M)):
        key = M[i].tobytes()
        if key in seen:
            mult[seen[key]] += 1
        else:
            seen[key] = len(reps)
            reps.append(i)
            mult.append(1)
    return M[np.ix_(reps, reps)], mult


def twin_free(M):
    """Twin-free: no two nodes are tie-twins, so the game can't be reduced by
    merging (all rows of M are distinct). A twin-free game is not a smaller game
    with a duplicated strategy."""
    M = np.asarray(M)
    return len({M[i].tobytes() for i in range(len(M))}) == len(M)


def k_paradoxical(M, k):
    """Every k-subset of strategies has a common strict beater.

    A strategy C strictly beats {A_1, ..., A_k} iff M[C, A_i] == +1 for
    every A_i (ties don't count). k=1 is the baseline paradoxical
    property (every strategy has a beater); k=2 is the classical
    Erdos-Schutte "two-paradox" / "paradoxical tournament" condition.
    """
    from itertools import combinations

    n = len(M)
    if k < 1 or k > n - 1:
        return False
    nodes = range(n)
    for subset in combinations(nodes, k):
        beaten = set(subset)
        if not any(all(M[c, a] == 1 for a in subset) for c in nodes if c not in beaten):
            return False
    return True


def connected(M):
    """Decisive (non-tie) subgraph is connected. BFS from node 0."""
    n = len(M)
    A = M != 0
    visited = np.zeros(n, dtype=bool)
    visited[0] = True
    frontier = [0]
    while frontier:
        nxt = []
        for v in frontier:
            for u in np.where(A[v] & ~visited)[0]:
                visited[u] = True
                nxt.append(int(u))
        frontier = nxt
    return bool(visited.all())


def canonical_key(M):
    """Canonical hashable key of M under node relabeling (nauty certificate).

    Two games get the same key iff they are isomorphic. Computed by nauty
    (via pynauty) on the directed-graph view of M -- exact and fast.
    """
    return pynauty.certificate(pynauty_graph(M))


_HASH_BYTES = 12  # 96-bit digest -> collisions negligible even at ~1e9 entries


def matrix_hash(M):
    """96-bit BLAKE2b digest of M's bytes, as a Python int orbit key.

    An int key costs ~half the memory of an equivalent `bytes` object in a
    set (no per-object bytes header), which is what keeps the large-n dedup
    set in RAM. Run-stable (unlike the salted built-in `hash`).
    """
    d = hashlib.blake2b(np.ascontiguousarray(M).tobytes(), digest_size=_HASH_BYTES).digest()
    return int.from_bytes(d, "little")


def orbit_hashes(M):
    """Set of integer hashes (see `matrix_hash`), one per node-relabeling of M.

    Memory-bounded dedup for the brute `_enumerate` search: a leaf is a
    duplicate iff `matrix_hash(leaf) in seen`, and `seen` only ever holds
    compact ints rather than full n*n matrices.
    """
    perms = _perms(len(M))
    ii = perms[:, :, None]
    jj = perms[:, None, :]
    orbit = np.ascontiguousarray(M[ii, jj])
    bl = hashlib.blake2b
    return {
        int.from_bytes(bl(o.tobytes(), digest_size=_HASH_BYTES).digest(), "little") for o in orbit
    }


def canonicalize(mats):
    """Deduplicate a list of matrices up to node relabeling."""
    out = {}
    for M in mats:
        out.setdefault(canonical_key(M), M)
    return list(out.values())
