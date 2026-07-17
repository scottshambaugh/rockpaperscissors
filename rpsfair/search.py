"""Enumerate fair RPS structures (balanced and inclusive)."""

import time
from itertools import permutations

import numpy as np
from tqdm import tqdm

from .cache import cached
from .equilibrium import has_fully_mixed, is_completely_mixed, maxmin_equilibrium
from .structure import connected, matrix_hash, orbit_hashes, paradoxical_batch

BATCH_SIZE = 8192


def _edges(n):
    return [(i, j) for i in range(n) for j in range(i + 1, n)]


def _vals_batch(start, end, v0, n_edges):
    """Generate vals[start:end] via radix-3 expansion of the trailing edges.

    vals[k, 0] = v0 (fixed by the M[0,1] >= 0 symmetry break).
    """
    bs = end - start
    out = np.empty((bs, n_edges), dtype=np.int8)
    out[:, 0] = v0
    idx = np.arange(start, end, dtype=np.int64)
    for k in range(n_edges - 1):
        out[:, k + 1] = (idx % 3).astype(np.int8) - 1
        idx //= 3
    return out


def _build_batch(vals, edge_i, edge_j, n):
    """Build a (B, n, n) batch of skew-symmetric matrices from edge values."""
    bs = vals.shape[0]
    M = np.zeros((bs, n, n), dtype=np.int8)
    M[:, edge_i, edge_j] = vals
    M[:, edge_j, edge_i] = -vals
    return M


def _enumerate(n, kind, predicate, batch_filter=None):
    """Batched brute-force enumerate skew-symmetric matrices.

    Pipeline per batch:
      1. radix-3 generate vals
      2. fancy-index build M_batch
      3. vectorized paradoxical check
      4. optional `batch_filter(M_batch) -> bool_mask` (e.g. row-sum-zero)
      5. per-survivor: connected check, orbit-bytes dedup, `predicate(M)`

    Symmetry break: every iso class has a rep with M[0,1] >= 0, so we
    skip the vals[0] == -1 branch (one third of the search).
    """
    edges = _edges(n)
    n_edges = len(edges)
    edge_i = np.array([e[0] for e in edges])
    edge_j = np.array([e[1] for e in edges])

    total_per_v0 = 3 ** (n_edges - 1)
    total = 2 * total_per_v0

    seen = set()
    results = []
    surv_paradox = surv_conn = surv_new = kept = 0
    t0 = time.perf_counter()
    pbar = tqdm(total=total, desc=f"{kind} n={n}", unit="cand", leave=False)

    for v0 in (0, 1):
        for batch_start in range(0, total_per_v0, BATCH_SIZE):
            batch_end = min(batch_start + BATCH_SIZE, total_per_v0)
            vals = _vals_batch(batch_start, batch_end, v0, n_edges)
            M_batch = _build_batch(vals, edge_i, edge_j, n)

            mask = paradoxical_batch(M_batch)
            if batch_filter is not None:
                mask &= batch_filter(M_batch)
            surv_paradox += int(mask.sum())

            for idx in np.where(mask)[0]:
                M = M_batch[idx]
                if matrix_hash(M) in seen:
                    continue
                surv_new += 1
                seen |= orbit_hashes(M)
                if not connected(M):
                    continue
                surv_conn += 1
                kept_item = predicate(M)
                if kept_item is not None:
                    kept += 1
                    results.append(kept_item)

            pbar.update(batch_end - batch_start)

    pbar.close()
    dt = time.perf_counter() - t0
    print(
        f"  {kind} n={n}: {total} candidates -> "
        f"{surv_paradox} paradox/extra -> {surv_new} new iso -> "
        f"{surv_conn} connected -> {kept} kept   ({dt:.2f}s)"
    )
    return results


def search_balanced(n):
    """Zero row sums + paradoxical + connected, up to isomorphism."""

    def go():
        uniform = np.ones(n) / n

        def predicate(M):
            return (M.copy(), uniform)

        def row_sum_zero(M_batch):
            return (M_batch.sum(axis=2) == 0).all(axis=1)

        return _enumerate(n, "balanced", predicate, batch_filter=row_sum_zero)

    return cached(f"balanced_n{n}", go)


def search_regular(n):
    """Regular structures (every node has identical W, T, L), up to isomorphism.

    Row-by-row with the exact per-row (W, T, L) constraint and a node-0 canonical
    pin, deduplicated by the nauty canonical key (bounded memory, scales past the
    orbit-hash dedup that OOMed at n=9). Thin wrapper over `search_regular_stream`.
    """
    return search_regular_stream(n)


def _msperm(p, q, z):
    """Distinct arrangements of p ones, q minus-ones, z zeros."""
    return set(permutations((1,) * p + (-1,) * q + (0,) * z))


def _stream_search(n, cache_name, desc, runs):
    """Generic streaming row-by-row enumeration with canonical-key dedup.

    Builds M one row at a time; each `expand(row, M)` in `runs` yields the
    candidate fill arrays for M[row, row+1:] (or None to prune the branch).
    A single canonical key per iso class is stored instead of the full orbit,
    so memory stays bounded and it scales past where the orbit-hash dedup OOMs.
    Multiple `runs` (e.g. one per win-count W) feed a shared dedup set.
    """
    from .generate import canonical_key_fast

    def go():
        M = np.zeros((n, n), dtype=np.int8)
        uniform = np.ones(n) / n
        seen = set()
        results = []
        leaves = new_iso = new_connected = 0
        t0 = time.perf_counter()
        pbar = tqdm(desc=f"{desc} n={n}", unit="leaf", leave=False)

        def rec(row, expand):
            nonlocal leaves, new_iso, new_connected
            if row == n:
                leaves += 1
                pbar.update(1)
                key = canonical_key_fast(M)
                if key in seen:
                    return
                seen.add(key)
                new_iso += 1
                if connected(M):
                    new_connected += 1
                    pbar.set_postfix_str(f"{new_connected} kept")
                    results.append((M.copy(), uniform))
                return
            arrs = expand(row, M)
            if arrs is None:
                return
            for arr in arrs:
                for j, v in enumerate(arr):
                    col = row + 1 + j
                    M[row, col] = v
                    M[col, row] = -v
                rec(row + 1, expand)

        for expand in runs:
            M[:] = 0
            rec(0, expand)
        pbar.close()
        print(
            f"  {desc}_stream n={n}: {leaves} leaves -> {new_iso} iso -> "
            f"{new_connected} connected   ({time.perf_counter() - t0:.1f}s)"
        )
        return results

    return cached(cache_name, go)


def search_balanced_stream(n):
    """Streaming balanced search: row-by-row build with the zero-row-sum
    constraint, dedup by one canonical key per class instead of storing every
    orbit. Memory stays bounded, so it scales past n=7 where the orbit-hash
    dedup OOMs.
    """

    def expand(row, M):
        forced = M[row, :row]
        fp = int((forced == 1).sum())
        fn_ = int((forced == -1).sum())
        target = -int(forced.sum())
        free = n - 1 - row
        if abs(target) > free:
            return None
        arrs = []
        for p in range(free + 1):
            q = p - target
            if q < 0 or p + q > free:
                continue
            if fp + p == 0 or fn_ + q == 0:  # node would never win / never lose
                continue
            arrs.extend(_msperm(p, q, free - p - q))
        return arrs

    return _stream_search(n, f"balanced_n{n}", "balanced", [expand])


def search_regular_stream(n):
    """Streaming regular search: row-by-row with the exact per-row (W, T, L)
    constraint and a node-0 canonical pin, one run per win-count W, dedup by
    canonical key so it stays in bounded memory past n=8.
    """

    def make_expand(W):
        def expand(row, M):
            forced = M[row, :row]
            fp = int((forced == 1).sum())
            fn_ = int((forced == -1).sum())
            rp = W - fp
            rn = W - fn_
            rz = (n - 1 - row) - rp - rn
            if rp < 0 or rn < 0 or rz < 0:
                return None
            if row == 0:  # pin node 0 to break the within-W relabeling symmetry
                return [tuple([1] * rp + [0] * rz + [-1] * rn)]
            return _msperm(rp, rn, rz)

        return expand

    runs = [make_expand(W) for W in range(1, (n - 1) // 2 + 1)]
    return _stream_search(n, f"regular_n{n}", "regular", runs)


def search_inclusive(n):
    """Paradoxical + connected + fully-mixed equilibrium, up to isomorphism.

    The cheap `has_fully_mixed` existence test gates membership; the stored
    equilibrium is the canonical leximin (max-min) point of O (computed only for
    kept structures), so all downstream odds/metrics are well-defined even when
    the equilibrium set is a higher-dimensional family (even n).
    """

    def go():
        def predicate(M):
            ok, _ = has_fully_mixed(M)
            if not ok:
                return None
            return (M.copy(), maxmin_equilibrium(M))

        return _enumerate(n, "inclusive", predicate)

    return cached(f"inclusive_n{n}", go)


def search_completely_mixed(n):
    """Completely mixed games (Kaplansky 1945): EVERY Nash equilibrium is fully
    mixed -- equivalently, the equilibrium is unique and fully mixed, so every
    strategy is *required* (played in every equilibrium), not merely playable.

    This is the strict top of the equilibrium axis: completely-mixed ⊂ inclusive.
    It is NOT comparable with balanced/regular the way those are with each other:
    a balanced or regular game at even n is never completely mixed (rank(M) is
    even, so even n forces a continuum of equilibria -- Kaplansky's parity
    observation), hence the immediate [] below; at odd n regular games can still
    fail (e.g. two of the twelve twin-free regular n=7 games have 4 extreme
    equilibria each).

    Implemented as a filter over `search_inclusive` (the unique equilibrium of a
    completely mixed game is fully mixed, so every completely mixed game is
    inclusive), with the cheap nullity-1 + one-signed-kernel test gating
    membership. The stored xs (leximin) is the unique equilibrium itself. No
    separate cache: the inclusive cache does the heavy lifting. Completely mixed
    implies twin-free automatically (a tie-twin pair i, j puts e_i - e_j in
    ker(M), breaking nullity 1).
    """
    if n % 2 == 0:
        return []
    return [(M, xs) for M, xs in search_inclusive(n) if is_completely_mixed(M)]


def search_balanced_fast(n):
    """Row-by-row pruned enumeration of balanced structures.

    Thin wrapper over `search_balanced_stream` (canonical-key dedup, bounded
    memory) -- the older orbit-hash dedup OOMed at n=8.
    """
    return search_balanced_stream(n)
