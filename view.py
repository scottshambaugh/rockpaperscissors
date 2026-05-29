"""view.py — inspect a single fair RPS structure.

Examples:
    uv run view.py --n 3
    uv run view.py --n 3 --labels Rock,Paper,Scissors
    uv run view.py --n 5 --kind balanced --index 2
    uv run view.py --n 4 --index 1 --save my.png
"""

import argparse
import os

import matplotlib.pyplot as plt

from rpsfair import add_colorbar, draw, pretty, search_balanced, search_inclusive

PLOTS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "plots")


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--n", type=int, required=True, help="number of strategies")
    p.add_argument("--kind", choices=["balanced", "inclusive"], default="inclusive")
    p.add_argument(
        "--index", type=int, default=0, help="0-based index into the search results (default: 0)"
    )
    p.add_argument(
        "--labels",
        default=None,
        help="Comma-separated label list, e.g. 'Rock,Paper,Scissors'. Defaults to WTL profiles.",
    )
    p.add_argument(
        "--save",
        default=None,
        help="Path to save the figure. Default: plots/single_<kind>_n<n>_<index>.png",
    )
    p.add_argument(
        "--no-plot", action="store_true", help="Skip the figure (still prints the text grid)"
    )
    p.add_argument(
        "--no-optimize-layout",
        action="store_true",
        help="Skip the minimize-edge-crossings reordering; place nodes in matrix index order",
    )
    p.add_argument(
        "--no-percentages",
        action="store_true",
        help="Suppress per-node equilibrium percentage labels",
    )
    args = p.parse_args()

    search = {"balanced": search_balanced, "inclusive": search_inclusive}[args.kind]
    items = search(args.n)
    if not 0 <= args.index < len(items):
        raise SystemExit(
            f"index {args.index} out of range — {args.kind}_n{args.n} has {len(items)} structures"
        )
    M, xs = items[args.index]

    labels = args.labels.split(",") if args.labels else None
    if labels is not None and len(labels) != args.n:
        raise SystemExit(f"got {len(labels)} labels for n={args.n}")

    print(f"{args.kind} n={args.n} index={args.index}")
    print(pretty(M, labels))
    print("\nequilibrium:")
    rendered = labels if labels else [f"node {i}" for i in range(args.n)]
    for lab, x in zip(rendered, xs, strict=True):
        print(f"  {lab:>8s}  {float(x) * 100:5.1f}%")

    if args.no_plot:
        return

    out = args.save or os.path.join(PLOTS_DIR, f"single_{args.kind}_n{args.n}_{args.index}.png")
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    fig, ax = plt.subplots(figsize=(5.2, 5.6))
    draw(
        ax,
        M,
        xs,
        node_labels=labels,
        title=f"{args.kind} n={args.n} #{args.index}",
        radius=0.34,
        optimize_layout=not args.no_optimize_layout,
        show_percentages=not args.no_percentages,
    )
    if not args.no_percentages:
        fig.subplots_adjust(bottom=0.12)
        cbar_ax = fig.add_axes([0.25, 0.06, 0.5, 0.025])
        add_colorbar(fig, cax=cbar_ax)
    plt.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
