# Fair generalizations of Rock-Paper-Scissors

Enumerate and visualize "fair" generalizations of Rock-Paper-Scissors on `n` strategies, represented as `n×n` skew-symmetric matrices with entries in `{-1, 0, +1}`. Three nested notions of fairness:

- **Regular** — every node has the **same** `(W, T, L)` profile. Strongest.
- **Balanced** — every row sums to zero (uniform equilibrium); profiles may differ.
- **Inclusive** — unique optimal mixed strategy is fully supported (every strategy gets positive play). Weakest.

`regular ⊂ balanced ⊂ inclusive`. All three require the game to be paradoxical (every strategy wins and loses at least once) and decisively connected. Structures are counted up to graph isomorphism.

Metrics computed per structure: `orbits` (structurally distinct strategies = automorphism orbits; two moves share an orbit iff some relabeling swaps them without changing the game), `|Aut|` (automorphism group size), articulation cuts, Gini of the equilibrium, tie fraction.

## Counts up to isomorphism

Each fairness cell is `total (twin-free) [prime]` — the total count up to isomorphism, the twin-free subcount in parentheses, and the modular-prime subcount in brackets. The `candidates` column is the raw brute-force search space, `3^(n choose 2)` (one of `{-1, 0, +1}` per upper-triangle edge), before any filtering or isomorphism reduction.

| n | candidates = 3^C(n,2) | two-paradox  | regular        | balanced         | inclusive         |
|---|-----------------------|--------------|----------------|------------------|-------------------|
| 3 | 27                    | 0 (0) [0]    | 1 (1) [1]      | 1 (1) [1]        | 1 (1) [1]         |
| 4 | 729                   | 0 (0) [0]    | 1 (1) [1]      | 1 (1) [1]        | 3 (2) [2]         |
| 5 | 59,049                | 0 (0) [0]    | 2 (2) [2]      | 4 (3) [3]        | 15 (8) [7]        |
| 6 | 14,348,907            | 0 (0) [0]    | 5 (4) [4]      | 16 (13) [13]     | 222 (177) [169]   |
| 7 | 10,460,353,203        | 1 (1) [1]    | 13 (12) [12]   | 175 (152) [152]  | —                 |
| 8 | 22,876,792,454,961    | 0 (0) [0]    | 82 (76) [76]   | —                | —                 |

Nesting (left to right, strict subsets): `two-paradox ⊂ regular ⊂ balanced ⊂ inclusive`. All games we enumerate are P₁ (every strategy has a beater — that's the `paradoxical()` baseline filter); the smallest *P₂* (two-paradox) appears at `n=7` (the Paley tournament `Q₇`).

- **`two-paradox`** = `P₂`: every *pair* of strategies has a common strict beater — a third strategy that beats both (ties don't count). Currently computed as a post-filter on `search_regular`; for larger `n` we could push the constraint into the row-by-row enumeration.
- **twin-free**: no two nodes are tie-twins. Two nodes are *twins* if they tie and share the same parents and children (equivalently, identical rows in `M`); merging a twin pair gives an `(n-1)`-node game, so a structure with twins is just a smaller game with a duplicated strategy. The twin-free subcount is the "genuinely new at this `n`" sequence — e.g. Brick+Boulder+Paper+Scissors reduces to RPS, and the 5-element Water/Fire/Clay/Sand/Grass game reduces to the 4-strategy Cop/K-9/Perp/Witness game.
- **prime** (`[…]`): *modular-prime* — no nontrivial **module**. A module is a set of moves that every outside move relates to identically (beats all, loses to all, or ties all); a twin pair is just a size-2 all-tie module. A non-prime game is a smaller **quotient** game with sub-games substituted into its moves (`G = H[M₁,…,M_k]`), so prime games are the genuine irreducible atoms — a strictly stronger notion than twin-free (`prime ⊆ twin-free`). The gap first appears at `n=5` inclusive: 8 twin-free but only 7 prime, the one extra being RPS with a move blown up into another RPS — twin-free (no duplicated rows) yet decomposable. The number of extreme equilibria factorizes over the decomposition tree (`num_equilibria` ↔ `neq_tree`), which is why the `n_eq` blow-up is a property of *reducible* games: `twin-free ⟹ n_eq ≤ n` is **false** (e.g. `RPSLS[cop×5]` is twin-free with `n_eq = 2⁵ = 32 > 20`), but the refined `prime ⟹ n_eq ≤ n` holds across everything enumerated.

The inclusive column for `n=7` is currently out of reach via brute force (`3^21 ≈ 10B` candidates). Balanced and regular use row-by-row pruning and scale to `n=7+`.

## Number of equilibria

The symmetric Nash equilibria of a game form a convex polytope `O = { p ∈ Δ : M p ≤ 0 }` (a strategy is a best response to itself iff no move beats it on average). The number of *solutions* is the number of **vertices (extreme equilibria)** of `O`: `1` = a unique equilibrium, `2` = a line segment (the Cop game), more = a polygon/polytope of equilibria whose interior is the fully-mixed family. `num_equilibria(M)` enumerates them exactly; `equilibrium_dim(M)` gives `dim O`.

Distribution of the solution count over the **twin-free** games at each `n` (`games` is the twin-free count):

| set       | n | games | #equilibria (vertices of `O`) |
|-----------|---|-------|-------------------------------|
| regular   | 3 | 1     | `{1:1}` |
| regular   | 4 | 1     | `{2:1}` |
| regular   | 5 | 2     | `{1:2}` |
| regular   | 6 | 4     | `{2:3, 5:1}` |
| regular   | 7 | 12    | `{1:10, 4:2}` |
| regular   | 8 | 76    | `{2:68, 4:3, 5:3, 6:1, 7:1}` |
| balanced  | 3 | 1     | `{1:1}` |
| balanced  | 4 | 1     | `{2:1}` |
| balanced  | 5 | 3     | `{1:3}` |
| balanced  | 6 | 13    | `{2:12, 5:1}` |
| balanced  | 7 | 152   | `{1:102, 3:10, 4:40}` |
| inclusive | 3 | 1     | `{1:1}` |
| inclusive | 4 | 2     | `{2:2}` |
| inclusive | 5 | 8     | `{1:7, 3:1}` |
| inclusive | 6 | 177   | `{2:176, 5:1}` |

A cell `{k:m}` reads "*m* twin-free games have *k* extreme equilibria". The **parity theorem** is visible directly: `rank(M)` is always even, so `dim O ≡ n−1 (mod 2)` and a unique equilibrium (a single vertex) is possible only at **odd `n`**. Every even-`n` game is non-unique — its equilibria always form a continuum (segment, polygon, …), never an isolated point.

In general the equilibrium condition is `M p ≤ 0` (every pure reply scores `≤ 0`), *not* `M p = 0`: only the *fully-mixed* equilibria are forced into `ker(M)`, while a partial-support equilibrium can sit on a face where an unplayed move's reply is strictly negative — and such a vertex is **not** a kernel vector. For an arbitrary game you therefore intersect `M p ≤ 0` with the simplex and enumerate vertices (what `equilibrium_vertices` does), rather than just taking `ker(M)`. For *this* category, though, the two coincide: across all enumerated games (1,241 vertices over 542 games) every single equilibrium vertex satisfies `M v = 0` to machine precision, so `O = ker(M) ∩ Δ` exactly and the null-space recipe is sufficient. This is a property of fair (paradoxical, connected) games, not of skew-symmetric games at large.

## Install

```
uv sync
```

## Run

```
uv run run.py
```

For each `n` in `NS` (default `[3, 4, 5, 6]` in `run.py`): prints the counts table, then for each fairness kind ranks the structures (decomposable filtered), prints the labeled upper-triangular grid of the top structure, and renders the ranking to `plots/n{n}_{regular,balanced,inclusive}.png`.

## Test

```
uv run pytest
```

## Lint / format

```
uv run ruff check .       # report issues
uv run ruff check --fix . # auto-fix
uv run ruff format .      # apply formatter
```

## Inspect a single structure

`view.py` is a CLI for one structure — prints the labeled upper-triangular grid and saves a plot.

```
uv run view.py --n 3                                  # default: inclusive, index 0
uv run view.py --n 3 --labels Rock,Paper,Scissors     # custom node names
uv run view.py --n 5 --kind balanced --index 2
uv run view.py --n 4 --index 1 --save my.png
uv run view.py --n 4 --index 1 --no-plot              # text only
```

## Labeling nodes

By default, nodes are labeled with their WTL profile tuple (e.g. `1·0·2` = 1 win, 0 ties, 2 losses). Override with a custom list to give game-specific names:

```python
from rpsfair import pretty, search_inclusive

M, xs = search_inclusive(3)[0]
print(pretty(M, labels=["Rock", "Paper", "Scissors"]))
```

`pretty(M, labels)` returns a labeled upper-triangular text grid using `+` (row beats col), `-` (row loses), `0` (tie). `show(M, labels)` prints it directly. `letter_labels(n)` yields `['A', 'B', ...]` if you want letter labels.

## Color scale

Plots use a global colormap (`viridis`) normalized to `0%`–`50%` equilibrium play rate, so fills are directly comparable across plots and across different `n`. `grid()` and `view.py` include a horizontal colorbar; `add_colorbar(fig, ...)` is exposed for custom figures.

## Caching

Search results are written to `cache/<name>.json` with human-readable filenames (`regular_n5.json`, `balanced_n5.json`, `inclusive_n5.json`, …). The filename does **not** encode the source — if you edit a search function, wipe the cache (`rm -rf cache/`) before re-running.

## Layout

```
rpsfair/
  cache.py        JSON disk cache
  structure.py    paradoxical / connected / regular / canonical_key / orbit_bytes
  equilibrium.py  equilibrium, has_fully_mixed (SVD null-space + LP)
  search.py       search_regular / search_balanced / search_balanced_fast / search_inclusive
  metrics.py      num_orbits, aut_size, num_cuts, gini, tie_fraction
  plot.py         draw, grid, best_layout, add_colorbar (matplotlib)
  display.py      pretty, show, wtl_labels, letter_labels (text)
run.py            driver — loops over n in NS, emits balanced/regular/inclusive plots
view.py           CLI for inspecting a single structure
tests/            pytest sanity tests for known counts
```
