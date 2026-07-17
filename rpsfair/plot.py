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


_MODBOX = {"tie": "#e08214", "order": "#3690c0", "prime": "#999999"}


_MOD_R = 0.22  # node radius (matches draw()'s default)
# inches of panel width per data-unit of view radius; chosen so a node renders at
# ~draw()'s physical size, keeping the fixed 9pt label legible at any content extent.
_MOD_IN_PER_UNIT = 1.25


def _modular_layout(M, R=_MOD_R, gap=0.30):
    """Bottom-up modular layout. Returns (pos, boxes, qedges, center, half).

    Each module is laid out in its own local frame and recentered on its content
    bounding box, so a module's children sit centered within its disk at every level
    (a big sub-module no longer biases the cluster). `half` is the content half-extent
    (natural square view is center ± half), exposed so callers can share a scale.
    """
    import math

    from .modular import modular_decomposition, named_subgame

    M = np.asarray(M)

    def label_of(node):
        reps = [min(c["members"]) for c in node["children"]]
        return named_subgame(M[np.ix_(reps, reps)]) or ""

    def layout(node):
        # returns (pos, boxes, qedges, radius), all centered on the content bbox
        if node["type"] == "leaf":
            return {min(node["members"]): np.zeros(2)}, [], [], R
        kids = [layout(c) for c in node["children"]]
        crad = [k[3] for k in kids]
        reps = [min(c["members"]) for c in node["children"]]
        k = len(kids)
        if k == 1:
            off = [np.zeros(2)]
        else:
            ring = max(
                (crad[i] + crad[(i + 1) % k] + gap) / (2 * math.sin(math.pi / k)) for i in range(k)
            )
            angs = np.linspace(math.pi / 2, math.pi / 2 + 2 * math.pi, k, endpoint=False)
            off = [ring * np.array([math.cos(a), math.sin(a)]) for a in angs]
        pos, boxes, qedges = {}, [], []
        for (cpos, cbox, cq, cr), o, child in zip(kids, off, node["children"], strict=False):
            for idx, p in cpos.items():
                pos[idx] = p + o
            boxes += [(bc + o, br, bt, bl) for bc, br, bt, bl in cbox]
            qedges += [(qa + o, ra, qb + o, rb, sg) for qa, ra, qb, rb, sg in cq]
            if child["type"] != "leaf":  # draw a disk for each internal child module
                boxes.append((o, cr, child["type"], label_of(child)))
        for ia in range(k):
            for ib in range(ia + 1, k):
                qedges.append((off[ia], crad[ia], off[ib], crad[ib], int(M[reps[ia], reps[ib]])))
        # recenter this module's content on its bounding box
        xe = [c for p in pos.values() for c in (p[0] - R, p[0] + R)] + [
            c for bc, br, *_ in boxes for c in (bc[0] - br, bc[0] + br)
        ]
        ye = [c for p in pos.values() for c in (p[1] - R, p[1] + R)] + [
            c for bc, br, *_ in boxes for c in (bc[1] - br, bc[1] + br)
        ]
        sh = np.array([(min(xe) + max(xe)) / 2, (min(ye) + max(ye)) / 2])
        pos = {i: p - sh for i, p in pos.items()}
        boxes = [(bc - sh, br, bt, bl) for bc, br, bt, bl in boxes]
        qedges = [(qa - sh, ra, qb - sh, rb, sg) for qa, ra, qb, rb, sg in qedges]
        rad = max(
            [float(np.hypot(*p)) + R for p in pos.values()]
            + [float(np.hypot(*bc)) + br for bc, br, *_ in boxes]
        )
        return pos, boxes, qedges, rad

    tree = modular_decomposition(M)
    pos, boxes, qedges, rad = layout(tree)  # root content already centered at origin
    return pos, boxes, qedges, (0.0, 0.0), rad + 0.12


def plot_modular(
    M,
    xs=None,
    ax=None,
    node_labels=None,
    title=None,
    path=None,
    dpi=140,
    view_half=None,
    view_center=None,
):
    """Draw a game with nodes laid out by their modular decomposition.

    Each module is a cluster; nested modules are nested translucent disks colored by
    quotient type (tie = orange, order = blue, prime = grey) and named by their
    quotient game, so the substitution structure G = H[M_1,...] is visible. Node fill
    is the equilibrium play rate (defaults to max-min); node size/font match draw().

    Pass `ax` to compose into a larger figure. Pass `view_half`/`view_center` to force
    a shared view scale across panels (so nodes render at the same physical size); the
    natural values come from `_modular_layout`.
    """
    from .equilibrium import maxmin_equilibrium

    M = np.asarray(M)
    n = len(M)
    if xs is None:
        xs = maxmin_equilibrium(M)
    R = _MOD_R
    pos, boxes, qedges, center, half = _modular_layout(M, R=R)
    if view_half is None:
        view_half = half
    if view_center is None:
        view_center = center

    own_fig = ax is None
    if own_fig:
        # size the figure so the node radius is ~draw()'s physical size at this scale
        s = max(4.5, 2 * view_half * _MOD_IN_PER_UNIT)
        fig, ax = plt.subplots(figsize=(s, s))

    for cen, radius, typ, label in sorted(boxes, key=lambda b: -b[1]):  # largest disk first
        col = _MODBOX.get(typ, "#999999")
        ax.add_patch(
            Circle(
                cen, radius, facecolor=col, alpha=0.15, edgecolor=col, lw=1.6, ls="--", zorder=0.5
            )
        )
        if label:
            ax.text(
                cen[0] - 0.71 * radius,
                cen[1] + 0.71 * radius,
                label,
                ha="center",
                va="center",
                fontsize=8.5,
                color=col,
                fontweight="bold",
                zorder=0.6,
            )

    for ca, ra, cb, rb, sign in qedges:  # one edge per sibling-module pair
        u = cb - ca
        d = np.hypot(*u)
        if d < 1e-9:
            continue
        u = u / d
        s0, e0 = ca + u * (ra + 0.05), cb - u * (rb + 0.05)
        if sign == 0:
            ax.plot([s0[0], e0[0]], [s0[1], e0[1]], color=TIE, ls="--", lw=1.3, zorder=1, alpha=0.7)
        else:
            src, dst = (s0, e0) if sign == 1 else (e0, s0)
            ax.add_patch(
                FancyArrowPatch(
                    src, dst, arrowstyle="-|>", mutation_scale=13, color=WIN, lw=1.5, zorder=2
                )
            )

    # WTL profile + play rate, rendered exactly as draw() (radius R, font 9, same offsets)
    profiles = [(int((r == 1).sum()), int((r == 0).sum()) - 1, int((r == -1).sum())) for r in M]
    for i in range(n):
        x, y = pos[i]
        fill = CMAP(NORM(xs[i]))
        tc = "white" if _luminance(fill) < 0.5 else "black"
        ax.add_patch(Circle((x, y), R, facecolor=fill, edgecolor="black", lw=1.3, zorder=3))
        lab = node_labels[i] if node_labels else "{}·{}·{}".format(*profiles[i])
        ax.text(
            x,
            y + 0.05,
            lab,
            ha="center",
            va="center",
            fontsize=9,
            fontweight="bold",
            color=tc,
            zorder=4,
        )
        ax.text(
            x,
            y - 0.08,
            f"{xs[i] * 100:.0f}%",
            ha="center",
            va="center",
            fontsize=9,
            fontweight="bold",
            color=tc,
            zorder=4,
        )
    if title:
        ax.set_title(title, fontsize=11, pad=6)
    ax.set_xlim(view_center[0] - view_half, view_center[0] + view_half)
    ax.set_ylim(view_center[1] - view_half, view_center[1] + view_half)
    ax.set_aspect("equal")
    ax.axis("off")
    if own_fig:
        if path:
            plt.savefig(path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        return fig
    return ax


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
    row_labels=None,
    dpi=140,
):
    """items: list of (M, xs). Render a grid with a shared colorbar.

    labels: optional list of label-lists (one length-n list per item).
        Use None to fall back to WTL profile labels on every node.
    A None entry in `items` leaves that cell blank -- pad an n-group to a full
    row to keep multi-n gardens grouped by n.
    row_labels: optional {row_index: text} -- a bold label in the left margin of
        that row (e.g. {0: "n=3", 1: "n=4"} to mark each group).
    optimize_layout / show_percentages: forwarded to draw().
    """
    n_items = len(items)
    ncols = min(ncols, max(1, n_items))
    nrows = (n_items + ncols - 1) // ncols
    # Fixed-inch margins so tall grids don't accrue proportional whitespace and
    # the suptitle never overlaps the top row.
    PANEL = 3.8
    top_in = 0.55 if suptitle else 0.12
    bottom_in = 0.70 if show_percentages else 0.12
    left_in = 1.05 if row_labels else 0.0
    fig_w = PANEL * ncols + left_in
    fig_h = PANEL * nrows + top_in + bottom_in
    fig, axes = plt.subplots(nrows, ncols, figsize=(fig_w, fig_h))
    axes = np.atleast_1d(axes).flatten()
    for i, item in enumerate(items):
        if item is None:
            axes[i].axis("off")
            continue
        M, xs = item
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
    rect = [left_in / fig_w, bottom_in / fig_h, 1.0, 1.0 - top_in / fig_h]
    plt.tight_layout(rect=rect)
    if suptitle:
        fig.suptitle(suptitle, fontsize=14, y=1.0 - 0.18 * top_in / fig_h, va="top")
    if row_labels:
        for row, label in row_labels.items():
            pos = axes[row * ncols].get_position()
            fig.text(
                0.5 * left_in / fig_w,
                (pos.y0 + pos.y1) / 2,
                label,
                rotation=0,
                va="center",
                ha="center",
                fontsize=22,
                fontweight="bold",
            )
    if show_percentages:
        # Fixed physical colorbar size (inches), independent of grid height
        bar_h_in, bar_bottom_in, bar_w_in = 0.16, 0.30, 3.0
        cbar_ax = fig.add_axes(
            [
                left_in / fig_w + (1 - left_in / fig_w) / 2 - (bar_w_in / fig_w) / 2,
                bar_bottom_in / fig_h,
                bar_w_in / fig_w,
                bar_h_in / fig_h,
            ]
        )
        add_colorbar(fig, cax=cbar_ax)
    if path:
        plt.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return fig


def equilibria_grid(M, path=None, ncols=4, labels=None, optimize_layout=True, dpi=140):
    """Plot every extreme equilibrium (vertex of O) of game M, plus the max-min point.

    Each panel is the same game graph -- same node layout -- with the nodes filled
    by that solution's play rates, so the support shift across the vertices (and
    where the leximin max-min point sits among them) is visible at a glance. The
    final panel is the canonical max-min equilibrium. The suptitle carries the
    game's metrics; each panel is titled by that solution's Gini coefficient.
    """
    from .equilibrium import equilibrium_vertices, maxmin_equilibrium
    from .metrics import gini, num_cuts, num_orbits, tie_fraction

    verts = equilibrium_vertices(M)
    mm = maxmin_equilibrium(M)
    items = [(M, v) for v in verts] + [(M, mm)]
    titles = [f"gini={gini(v):.2f}" for v in verts] + [f"max-min  gini={gini(mm):.2f}"]
    labels_arg = [labels] * len(items) if labels is not None else None
    parts = [
        f"n={len(M)}",
        f"orbits={num_orbits(M)}",
        f"n_eq={len(verts)}",
        f"ties={tie_fraction(M):.2f}",
        f"gini={gini(mm):.2f}",
        f"cuts={num_cuts(M)}",
    ]
    return grid(
        items,
        ncols=ncols,
        titles=titles,
        labels=labels_arg,
        suptitle="  ".join(parts),
        path=path,
        optimize_layout=optimize_layout,
        show_percentages=True,
        dpi=dpi,
    )
