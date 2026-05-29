"""Render fair RPS structures as graph diagrams."""

import os
from itertools import permutations

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize
from matplotlib.patches import Circle, FancyArrowPatch

# Best-effort: register Symbola if installed so emoji glyphs in labels render
# instead of showing as missing-glyph boxes. Color-emoji fonts (NotoColorEmoji)
# aren't supported by matplotlib's FreeType build.
for _path in (
    "/usr/share/fonts/truetype/ancient-scripts/Symbola_hint.ttf",
    "/usr/share/fonts/truetype/ttf-ancient-scripts/Symbola_hint.ttf",
):
    if os.path.exists(_path):
        try:
            mpl.font_manager.fontManager.addfont(_path)
            mpl.rcParams["font.family"] = ["DejaVu Sans", "Symbola"]
        except Exception:
            pass
        break

WIN = "#222222"
TIE = "#bbbbbb"

# Global play-rate color scale: 0% (cold) → 50% (hot). Same mapping for every n,
# so colors are comparable across plots.
CMAP = plt.get_cmap("viridis")
NORM = Normalize(vmin=0.0, vmax=0.5)


def _luminance(rgb):
    return 0.299 * rgb[0] + 0.587 * rgb[1] + 0.114 * rgb[2]


# Beyond this n, (n-1)!/2 exact orderings get expensive, so fall back to
# deterministic local search.
_EXACT_MAX_N = 9


def _crossing_arrays(M):
    """Precompute disjoint edge pairs (in node space) and per-pair weight class.

    Two chords cross only if their endpoints are disjoint; sharing an endpoint
    never crosses. Returns numpy arrays of the four endpoints per pair plus
    boolean masks for the solid-solid / mixed / tie-tie weight classes.
    """
    n = len(M)
    edges = [(i, j, M[i, j] != 0) for i in range(n) for j in range(i + 1, n)]
    u1, v1, u2, v2, solid_solid, tie_tie = [], [], [], [], [], []
    for a in range(len(edges)):
        i1, j1, s1 = edges[a]
        for b in range(a + 1, len(edges)):
            i2, j2, s2 = edges[b]
            if i2 in (i1, j1) or j2 in (i1, j1):
                continue
            u1.append(i1)
            v1.append(j1)
            u2.append(i2)
            v2.append(j2)
            solid_solid.append(s1 and s2)
            tie_tie.append((not s1) and (not s2))
    ss = np.array(solid_solid, dtype=bool)
    tt = np.array(tie_tie, dtype=bool)
    return (
        np.array(u1),
        np.array(v1),
        np.array(u2),
        np.array(v2),
        ss,
        tt,
        ~ss & ~tt,
    )


def _crossing_key(pos, arrays):
    """(solid-solid, mixed, tie-tie) crossing counts for a slot assignment."""
    u1, v1, u2, v2, ss, tt, mixed = arrays
    lo = np.minimum(pos[u1], pos[v1])
    hi = np.maximum(pos[u1], pos[v1])
    c, d = pos[u2], pos[v2]
    cross = ((lo < c) & (c < hi)) != ((lo < d) & (d < hi))
    return (int((cross & ss).sum()), int((cross & mixed).sum()), int((cross & tt).sum()))


def best_layout(M):
    """Cyclic node order minimizing crossings: solid-solid > solid-tie > tie-tie.

    Crossings on a circle are invariant under rotation and reflection, so we
    pin node 0 to slot 0 and break the reflection -- searching (n-1)!/2
    orderings instead of n! (a 2n-fold reduction). Exact through n=9; for
    larger n a deterministic local search (best-improvement node swaps from
    the identity order) replaces the exhaustive sweep.
    """
    n = len(M)
    if n <= 3:
        return tuple(range(n))
    arrays = _crossing_arrays(M)

    if n <= _EXACT_MAX_N:
        best_key, best = None, None
        pos = np.empty(n, dtype=np.int64)
        for tail in permutations(range(1, n)):
            if tail[0] > tail[-1]:  # reflection symmetry break
                continue
            pos[0] = 0
            pos[list(tail)] = np.arange(1, n)
            key = _crossing_key(pos, arrays)
            if best_key is None or key < best_key:
                best_key, best = key, (0, *tail)
        return best

    # Heuristic: best-improvement swaps until no pair-swap reduces crossings.
    order = list(range(n))
    pos = np.argsort(order).astype(np.int64)
    cur = _crossing_key(pos, arrays)
    improved = True
    while improved:
        improved = False
        for a in range(n):
            for b in range(a + 1, n):
                order[a], order[b] = order[b], order[a]
                pos = np.argsort(order).astype(np.int64)
                key = _crossing_key(pos, arrays)
                if key < cur:
                    cur, improved = key, True
                else:
                    order[a], order[b] = order[b], order[a]
    return tuple(order)


def draw(
    ax,
    M,
    xs,
    node_labels=None,
    title=None,
    radius=0.22,
    optimize_layout=True,
    show_percentages=True,
):
    """Draw one structure onto a matplotlib axis.

    Node fill encodes equilibrium play rate via the global CMAP/NORM
    (0% → 50%). Default node labels are the WTL profile tuple; pass
    `node_labels` to override. `radius` controls node size; text offsets
    and arrow gaps scale with it. Set `optimize_layout=False` to keep
    nodes in M's index order (skip the crossings-minimization search).
    Set `show_percentages=False` to suppress the per-node `xs[i]` label
    (useful for balanced/regular where every node is uniform).
    """
    n = len(M)
    perm = best_layout(M) if optimize_layout else tuple(range(n))
    ang = np.linspace(np.pi / 2, np.pi / 2 + 2 * np.pi, n, endpoint=False)
    pos = np.zeros((n, 2))
    for i, node in enumerate(perm):
        pos[node] = (np.cos(ang[i]), np.sin(ang[i]))

    scale = radius / 0.22
    arrow_off = radius + 0.02
    font = 9 * scale

    profiles = [(int((r == 1).sum()), int((r == 0).sum()) - 1, int((r == -1).sum())) for r in M]
    for i in range(n):
        for j in range(i + 1, n):
            if M[i, j] == 0:
                ax.plot(
                    [pos[i, 0], pos[j, 0]],
                    [pos[i, 1], pos[j, 1]],
                    color=TIE,
                    ls="--",
                    lw=1.2,
                    zorder=1,
                    alpha=0.6,
                )
    for i in range(n):
        for j in range(i + 1, n):
            if M[i, j] == 0:
                continue
            src, dst = (pos[i], pos[j]) if M[i, j] == 1 else (pos[j], pos[i])
            d = np.hypot(*(dst - src))
            s = src + (dst - src) / d * arrow_off
            e = dst - (dst - src) / d * arrow_off
            ax.add_patch(
                FancyArrowPatch(
                    s, e, arrowstyle="-|>", mutation_scale=13, color=WIN, lw=1.4, zorder=2
                )
            )
    for i, (x, y) in enumerate(pos):
        fill = CMAP(NORM(xs[i]))
        text_color = "white" if _luminance(fill) < 0.5 else "black"
        ax.add_patch(Circle((x, y), radius, facecolor=fill, edgecolor="black", lw=1.3, zorder=3))
        lab = node_labels[i] if node_labels else "{}·{}·{}".format(*profiles[i])
        lab_y = y if not show_percentages else y + 0.05 * scale
        ax.text(
            x,
            lab_y,
            lab,
            ha="center",
            va="center",
            fontsize=font,
            fontweight="bold",
            color=text_color,
            zorder=4,
        )
        if show_percentages:
            ax.text(
                x,
                y - 0.08 * scale,
                f"{xs[i] * 100:.0f}%",
                ha="center",
                va="center",
                fontsize=font,
                fontweight="bold",
                color=text_color,
                zorder=4,
            )
    if title:
        ax.set_title(title, fontsize=10, pad=6)
    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-1.5, 1.5)
    ax.set_aspect("equal")
    ax.axis("off")


def add_colorbar(fig, ax=None, cax=None, orientation="horizontal", label="equilibrium play rate"):
    """Attach the global play-rate colorbar (0%-50%) to a figure."""
    sm = plt.cm.ScalarMappable(cmap=CMAP, norm=NORM)
    cbar = fig.colorbar(sm, ax=ax, cax=cax, orientation=orientation)
    cbar.set_label(label)
    ticks = np.linspace(0.0, 0.5, 6)
    cbar.set_ticks(ticks)
    cbar.set_ticklabels([f"{int(t * 100)}%" for t in ticks])
    return cbar


def grid(
    items,
    ncols=4,
    titles=None,
    suptitle=None,
    path=None,
    labels=None,
    optimize_layout=True,
    show_percentages=True,
):
    """items: list of (M, xs). Render a grid with a shared colorbar.

    labels: optional list of label-lists (one length-n list per item).
        Use None to fall back to WTL profile labels on every node.
    optimize_layout / show_percentages: forwarded to draw().
    """
    n_items = len(items)
    ncols = min(ncols, max(1, n_items))
    nrows = (n_items + ncols - 1) // ncols
    fig_h = 3.8 * nrows + 0.4
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.8 * ncols, fig_h))
    axes = np.atleast_1d(axes).flatten()
    for i, (M, xs) in enumerate(items):
        draw(
            axes[i],
            M,
            xs,
            title=titles[i] if titles else None,
            node_labels=labels[i] if labels else None,
            optimize_layout=optimize_layout,
            show_percentages=show_percentages,
        )
    for ax in axes[n_items:]:
        ax.axis("off")
    if suptitle:
        fig.suptitle(suptitle, fontsize=12, y=0.99)
    plt.tight_layout(rect=[0, 0.05, 1, 1])
    if show_percentages:
        # Fixed physical colorbar size (inches), independent of grid height
        bar_h_in, bar_bottom_in, bar_w_in = 0.16, 0.32, 3.0
        fig_w = 3.8 * ncols
        cbar_ax = fig.add_axes(
            [
                0.5 - (bar_w_in / fig_w) / 2,
                bar_bottom_in / fig_h,
                bar_w_in / fig_w,
                bar_h_in / fig_h,
            ]
        )
        add_colorbar(fig, cax=cbar_ax)
    if path:
        plt.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    return fig
