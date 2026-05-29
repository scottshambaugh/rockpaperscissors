"""Enumerate fair RPS structures (balanced and inclusive)."""

import time
from itertools import permutations

import numpy as np
from tqdm import tqdm

from .cache import cached
from .equilibrium import has_fully_mixed, maxmin_equilibrium
from .structure import connected, matrix_hash, orbit_hashes

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


def _paradoxical_batch(M_batch):
    has_w = (M_batch == 1).any(axis=2).all(axis=1)
    has_l = (M_batch == -1).any(axis=2).all(axis=1)
    return has_w & has_l


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

            mask = _paradoxical_batch(M_batch)
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
    """Enumerate regular structures (every node has identical W, T, L).

    Each row is forced to have exactly W ones, T zeros, L=W minus-ones —
    a much stronger per-row constraint than balanced's row-sum-zero. We
    iterate over valid W and recurse row-by-row with arithmetic pruning,
    so most balanced-but-non-regular matrices are never built.

    Symmetry break: node 0's row is pinned to the single canonical pattern
    (+1 * W, 0 * T, -1 * W) instead of all its multiset permutations. Every
    isomorphism class has a labeling with node 0 sorted this way (just permute
    the other nodes), so this is exact, and it slashes the leaf count by ~n²
    (n=8: 2.56M -> 14k leaves). Remaining duplicates are removed by the
    orbit-bytes dedup at the leaf.
    """

    def go():
        results = []
        seen = set()
        uniform = np.ones(n) / n
        leaves = 0
        new_iso = 0
        new_connected = 0
        M = np.zeros((n, n), dtype=np.int8)
        t0 = time.perf_counter()

        def msperm(p, q, z):
            return set(permutations((1,) * p + (-1,) * q + (0,) * z))

        def rec(row, W):
            nonlocal leaves, new_iso, new_connected
            if row == n:
                leaves += 1
                h = matrix_hash(M)
                if h in seen:
                    return
                new_iso += 1
                seen.update(orbit_hashes(M))
                if connected(M):
                    new_connected += 1
                    results.append((M.copy(), uniform))
                return
            forced = M[row, :row]
            fp = int((forced == 1).sum())
            fn_ = int((forced == -1).sum())
            rp = W - fp
            rn = W - fn_
            free = n - 1 - row
            rz = free - rp - rn
            if rp < 0 or rn < 0 or rz < 0:
                return
            # pin node 0 to its canonical sorted pattern (exact, see docstring)
            row0_pat = (tuple([1] * rp + [0] * rz + [-1] * rn),)
            arrs = row0_pat if row == 0 else msperm(rp, rn, rz)
            for arr in arrs:
                for j, v in enumerate(arr):
                    col = row + 1 + j
                    M[row, col] = v
                    M[col, row] = -v
                rec(row + 1, W)

        # Valid (W, T, L) with W=L>=1 (paradoxical) and W+T+L = n-1
        for W in range(1, (n - 1) // 2 + 1):
            M[:] = 0
            rec(0, W)

        dt = time.perf_counter() - t0
        print(
            f"  regular n={n}: {leaves} leaves -> "
            f"{new_iso} new iso -> {new_connected} connected   ({dt:.2f}s)"
        )
        return results

    return cached(f"regular_n{n}", go)


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


def search_balanced_fast(n):
    """Row-by-row pruned enumeration of balanced structures. Good to n=7."""

    def go():
        M = np.zeros((n, n), dtype=np.int8)
        uniform = np.ones(n) / n
        seen = set()
        results = []
        # Tracking
        leaves = 0
        new_iso = 0
        new_connected = 0
        t0 = time.perf_counter()

        def msperm(p, q, z):
            return set(permutations((1,) * p + (-1,) * q + (0,) * z))

        def rec(row):
            nonlocal leaves, new_iso, new_connected
            if row == n:
                leaves += 1
                h = matrix_hash(M)
                if h in seen:
                    return
                new_iso += 1
                seen.update(orbit_hashes(M))
                if connected(M):
                    new_connected += 1
                    results.append((M.copy(), uniform))
                return
            forced = M[row, :row]
            fp = int((forced == 1).sum())
            fn_ = int((forced == -1).sum())
            target = -int(forced.sum())
            free = n - 1 - row
            if abs(target) > free:
                return
            for p in range(free + 1):
                q = p - target
                if q < 0 or p + q > free:
                    continue
                if fp + p == 0 or fn_ + q == 0:
                    continue
                z = free - p - q
                for arr in msperm(p, q, z):
                    for j, v in enumerate(arr):
                        col = row + 1 + j
                        M[row, col] = v
                        M[col, row] = -v
                    rec(row + 1)

        rec(0)
        dt = time.perf_counter() - t0
        print(
            f"  balanced_fast n={n}: {leaves} leaves -> "
            f"{new_iso} new iso -> {new_connected} connected   ({dt:.2f}s)"
        )
        return results

    return cached(f"balanced_n{n}", go)
