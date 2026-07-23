"""Text rendering of fair RPS structures."""

import numpy as np

SYM = {1: "+", -1: "-", 0: "0"}


def letter_labels(n):
    """Return ['A', 'B', ..., n-th letter] - opt-in alternative to WTL."""
    return [chr(ord("A") + i) for i in range(n)]


def wtl_labels(M):
    """Per-node 'W·T·L' profile strings - the default label scheme."""
    M = np.asarray(M)
    return [f"{int((r == 1).sum())}·{int((r == 0).sum()) - 1}·{int((r == -1).sum())}" for r in M]


def pretty(M, labels=None):
    """Return a labeled upper-triangular text grid for M.

    Cell M[i, j] (i < j) shows '+' if row i beats col j, '-' if row i loses,
    '0' if they tie. Lower triangle and diagonal are blank. Default labels
    are the WTL profile per node; pass a custom list for game-specific names.
    """
    M = np.asarray(M)
    n = len(M)
    if labels is None:
        labels = wtl_labels(M)
    if len(labels) != n:
        raise ValueError(f"got {len(labels)} labels for n={n}")
    w = max(len(s) for s in labels)

    def cell(s):
        return f"{s:>{w}}"

    header = " " * (w + 3) + "  ".join(cell(s) for s in labels)
    sep = " " * (w + 2) + "+" + "-" * ((w + 2) * n + 1)
    rows = []
    for i in range(n):
        parts = [cell(SYM[int(M[i, j])]) if j > i else cell("") for j in range(n)]
        rows.append(f"{cell(labels[i])} | " + "  ".join(parts))
    return "\n".join([header, sep, *rows])


def show(M, labels=None):
    """Print pretty(M, labels) to stdout."""
    print(pretty(M, labels))
