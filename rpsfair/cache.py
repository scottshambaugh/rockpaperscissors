"""JSON-backed disk cache with human-readable filenames."""

import json
import os

import numpy as np

CACHE_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "cache")
os.makedirs(CACHE_DIR, exist_ok=True)


def _path(name):
    return os.path.join(CACHE_DIR, f"{name}.json")


def load(name):
    """Load a list of (M, xs) from cache/<name>.json, or None if absent."""
    path = _path(name)
    if not os.path.exists(path):
        return None
    with open(path) as f:
        data = json.load(f)
    return [
        (np.array(s["M"], dtype=np.int8), np.array(s["xs"], dtype=float))
        for s in data["structures"]
    ]


def save(name, items):
    """Persist a list of (M, xs) to cache/<name>.json.

    Inner arrays (matrix rows, xs) are kept on a single line each for
    readability.
    """
    lines = [
        "{",
        f'  "name": {json.dumps(name)},',
        f'  "count": {len(items)},',
        '  "structures": [',
    ]
    for idx, (M, xs) in enumerate(items):
        suffix = "," if idx < len(items) - 1 else ""
        rows = M.astype(int).tolist()
        xs_compact = [round(float(x), 6) for x in xs]
        lines.append("    {")
        lines.append('      "M": [')
        for ri, row in enumerate(rows):
            comma = "," if ri < len(rows) - 1 else ""
            lines.append(f"        {json.dumps(row)}{comma}")
        lines.append("      ],")
        lines.append(f'      "xs": {json.dumps(xs_compact)}')
        lines.append(f"    }}{suffix}")
    lines.append("  ]")
    lines.append("}")
    with open(_path(name), "w") as f:
        f.write("\n".join(lines) + "\n")


def cached(name, fn):
    """Return cached value if present, else compute, persist, and return."""
    hit = load(name)
    if hit is not None:
        return hit
    result = fn()
    save(name, result)
    return result
