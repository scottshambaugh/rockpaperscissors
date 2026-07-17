"""run.py — enumerate, rank, and visualize fair RPS structures across n."""

import argparse
import os
import time

from rpsfair import (
    automorphisms,
    equilibrium_dim,
    gini,
    grid,
    is_prime,
    num_cuts,
    num_equilibria,
    num_orbits,
    pretty,
    reduce_twins,
    search_balanced,
    search_completely_mixed,
    search_inclusive,
    search_regular,
    search_two_paradox,
    tie_fraction,
    twin_free,
)

PLOTS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "plots")
os.makedirs(PLOTS_DIR, exist_ok=True)

NS = [3, 4, 5, 6]
# inclusive equilibrium is non-uniform; balanced/regular are always uniform 1/n
PERCENTAGES = {"regular": False, "balanced": False, "inclusive": True, "completely_mixed": True}


def rank(structures, sort=True):
    """Annotate structures with metrics and (optionally) sort: orbits desc, cuts asc, ties asc, gini asc.

    All connected structures are kept (matching the published counts); cut
    vertices are reported as the `cuts` metric rather than filtered out.
    """
    ranked = []
    for M, xs in structures:
        autos = automorphisms(M)  # compute the group once
        ranked.append(
            {
                "M": M,
                "xs": xs,
                # structurally distinct strategies (automorphism orbits) -- the
                # parameterization-free replacement for the old "roles" metric
                "orbits": num_orbits(M, autos),
                "aut": len(autos),
                "gini": gini(xs),
                "ties": tie_fraction(M),
                "cuts": num_cuts(M),  # articulation points (0 = robustly connected)
                "core": len(reduce_twins(M)[0]),
                # modular-prime: no nontrivial module (stronger than twin-free)
                "prime": is_prime(M),
                # Dimension of the full equilibrium polytope, including games
                # whose equilibria lie only on the simplex boundary.
                "eqdim": equilibrium_dim(M),
                # number of extreme equilibria (vertices of O): 1 = unique solution
                "nverts": num_equilibria(M),
            }
        )
    if sort:
        # orbits desc (fewest swappable nodes first), then cuts asc, ties asc, gini asc
        ranked.sort(key=lambda r: (-r["orbits"], r["cuts"], r["ties"], r["gini"]))
    return ranked


def report_and_plot(
    kind,
    n,
    structures,
    sort=True,
    optimize_layout=True,
    twin_free_only=False,
    prime_only=False,
    ncols=4,
    dpi=140,
):
    """Print ranking + top text grid, save the multi-panel plot.

    twin_free_only: keep only twin-free structures (no tie-twin duplicates).
    prime_only: keep only modular-prime structures (no nontrivial module; stronger
        than twin-free). Tags the plot/file as the prime set; takes precedence.
    ncols: panels per row in the saved grid (e.g. 8 for the wide n=6 gardens).
    """
    ranked = rank(structures, sort=sort)
    if prime_only:
        ranked = [r for r in ranked if r["prime"]]
    elif twin_free_only:
        ranked = [r for r in ranked if r["core"] == n]
    label = "ranked" if sort else "unranked"
    suffix = "_prime" if prime_only else ("_twinfree" if twin_free_only else "")
    descr = "prime " if prime_only else ("twin-free " if twin_free_only else "")
    n_family = sum(1 for r in ranked if r["eqdim"] > 0)
    n_cut = sum(1 for r in ranked if r["cuts"] > 0)
    print(f"\nn={n} {kind} {descr}({label}): {len(ranked)} structures")
    print(f"    equilibrium: {len(ranked) - n_family} unique, {n_family} with a continuous family")
    print(f"    {n_cut} have a cut vertex (articulation point)")
    for i, r in enumerate(ranked, 1):
        core_tag = "twin-free" if r["core"] == n else f"core={r['core']}"
        eq_tag = "eq=unique" if r["nverts"] == 1 else f"eq={r['nverts']} vertices ({r['eqdim']}D)"
        print(
            f"  #{i:<2d} orbits={r['orbits']} |Aut|={r['aut']} "
            f"ties={r['ties']:.2f} gini={r['gini']:.2f} cuts={r['cuts']} {core_tag} {eq_tag}"
        )
    if ranked:
        print("  Top structure (#1):")
        for line in pretty(ranked[0]["M"]).splitlines():
            print(f"    {line}")
    path = os.path.join(PLOTS_DIR, f"n{n}_{kind}{suffix}.png")

    def title(i, r):
        # gini varies only when the equilibrium is non-uniform (inclusive)
        # order follows the sort key: orbits, then ties, then gini
        parts = [
            f"#{i + 1}",
            f"orbits={r['orbits']}",
            f"n_eq={r['nverts']}",
            f"ties={r['ties']:.2f}",
        ]
        if PERCENTAGES[kind]:
            parts.append(f"gini={r['gini']:.2f}")
        if r["cuts"]:
            parts.append(f"cuts={r['cuts']}")
        return " ".join(parts)

    titles = [title(i, r) for i, r in enumerate(ranked)]
    grid(
        [(r["M"], r["xs"]) for r in ranked],
        ncols=ncols,
        titles=titles,
        suptitle=f"n={n} {kind} {descr}structures, {label}",
        path=path,
        optimize_layout=optimize_layout,
        show_percentages=PERCENTAGES[kind],
        dpi=dpi,
    )
    print(f"  wrote {path}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument(
        "--no-rank",
        action="store_true",
        help="Preserve enumeration order in plots (skip the orbits/cuts/ties/gini sort)",
    )
    ap.add_argument(
        "--no-optimize-layout",
        action="store_true",
        help="Place nodes in matrix index order (skip crossings minimization)",
    )
    args = ap.parse_args()

    # each cell shows "total(twin-free)" -- twin-free = no tie-twin duplicates
    print("n   two-paradox      regular     balanced      inclusive   compl-mixed")
    print("-   -----------   ----------   ----------   ------------   -----------")
    for n in NS:
        t0 = time.perf_counter()
        reg = search_regular(n)
        bal = search_balanced(n)
        inc = search_inclusive(n)
        cm = search_completely_mixed(n)  # unique fully-mixed eq; 0 at even n
        p2 = search_two_paradox(n)  # authoritative: a filter over all tournaments

        def cell(structs):
            n_tf = sum(twin_free(M) for M, _ in structs)
            return f"{len(structs)}({n_tf})"

        dt = time.perf_counter() - t0
        print(
            f"{n:1d}   {cell(p2):>11}   "
            f"{cell(reg):>10}   {cell(bal):>10}   {cell(inc):>12}   {cell(cm):>11}   ({dt:.2f}s)"
        )

    for n in NS:
        # nested subsets: regular ⊂ balanced ⊂ inclusive
        report_and_plot(
            "regular",
            n,
            search_regular(n),
            sort=not args.no_rank,
            optimize_layout=not args.no_optimize_layout,
        )
        report_and_plot(
            "balanced",
            n,
            search_balanced(n),
            sort=not args.no_rank,
            optimize_layout=not args.no_optimize_layout,
        )
        report_and_plot(
            "inclusive",
            n,
            search_inclusive(n),
            sort=not args.no_rank,
            optimize_layout=not args.no_optimize_layout,
        )
        # completely mixed ⊂ inclusive; empty at even n (parity theorem)
        cm = search_completely_mixed(n)
        if cm:
            report_and_plot(
                "completely_mixed",
                n,
                cm,
                sort=not args.no_rank,
                optimize_layout=not args.no_optimize_layout,
            )


if __name__ == "__main__":
    main()
