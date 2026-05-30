## See *[The Secret Garden of Rock-Paper-Scissors](https://theshamblog.com/the-secret-garden-of-rock-paper-scissors/)* for context, details, and results

This repo was pretty much entirely vibe coded useing Claude Code. I understood and could verify the initial implementation on low n, but was hands off for the optimizations that were necessary for searching higher n, beyond verifying that the algorithms still reproduced the values we already found. All this to say, take the results here with a grain of salt and I think it would be an interesting exercise to implement this yourself and figure out optimizations for faster game searches.

----

# Fair generalizations of Rock-Paper-Scissors

Enumerate and visualize "fair" generalizations of Rock-Paper-Scissors on `n` strategies, represented as `n×n` skew-symmetric matrices with entries in `{-1, 0, +1}`. Three nested notions of fairness:

- **Regular** — every node has the **same** `(W, T, L)` profile. Strongest.
- **Balanced** — every row sums to zero (uniform equilibrium); profiles may differ.
- **Inclusive** — unique optimal mixed strategy is fully supported (every strategy gets positive play). Weakest.

`regular ⊂ balanced ⊂ inclusive`. All three require the game to be paradoxical (every strategy wins and loses at least once) and decisively connected. Structures are counted up to graph isomorphism.

Metrics computed per structure: `orbits` (structurally distinct strategies = automorphism orbits; two moves share an orbit iff some relabeling swaps them without changing the game), `|Aut|` (automorphism group size), articulation cuts, Gini of the equilibrium, tie fraction (over all `n²` matchups including the always-tie diagonal, so a tournament with no draws is `1/n`).

## Counts up to isomorphism

Each fairness cell is `total (twin-free) [prime]` — the total count up to isomorphism, the twin-free subcount in parentheses, and the modular-prime subcount in brackets. The `candidates` column is the raw brute-force search space, `3^(n choose 2)` (one of `{-1, 0, +1}` per upper-triangle edge), before any filtering. The `iso classes` column collapses that by relabeling: the number of distinct games up to isomorphism (`candidates` modulo the `S_n` action, computed in closed form by Burnside over the cycle types — `count_iso_classes(n)`), which is the real size of the space the fair games sit inside. It matches the enumerated sequence `1, 2, 7, 42, 582, 21480, …` for small `n` and stays comparatively tame (≈ `4.2 × 10¹¹` at `n=9`) where the `candidates` column has already passed `10¹⁷`. It is *not* simply `3^C(n,2) / n!`: that ratio is only the identity-permutation term of the Burnside sum (and isn't even an integer), so it undercounts by the automorphism correction — the true value is slightly larger, by a factor `1.032 / 1.014 / 1.006` at `n = 7 / 8 / 9` (the correction vanishes as games become rigid).

| n | candidates = 3^C(n,2)   | iso classes     | two-paradox     | regular            | balanced           | inclusive           |
|---|-------------------------|-----------------|-----------------|--------------------|--------------------|---------------------|
| 3 | 27                      | 7               | 0 (0) [0]       | 1 (1) [1]          | 1 (1) [1]          | 1 (1) [1]           |
| 4 | 729                     | 42              | 0 (0) [0]       | 1 (1) [1]          | 1 (1) [1]          | 3 (2) [2]           |
| 5 | 59,049                  | 582             | 0 (0) [0]       | 2 (2) [2]          | 4 (3) [3]          | 15 (8) [7]          |
| 6 | 14,348,907              | 21,480          | 0 (0) [0]       | 5 (4) [4]          | 16 (13) [13]       | 222 (177) [169]     |
| 7 | 10,460,353,203          | 2,142,288       | 1 (1) [1]       | 13 (12) [12]       | 175 (152) [152]    | 10525 (9401) [9350] |
| 8 | 22,876,792,454,961      | 575,016,219     | 4 (4) [3]       | 82 (76) [76]       | 5274 (4921) [4917] | —                   |
| 9 | 150,094,635,296,999,121 | 415,939,243,032 | 221 (221) [197] | 2016 (1973) [1972] | —                  | —                   |

**OEIS cross-references.** Two of these columns are known sequences, which both validates the code and extends the table:

- `iso classes` is **[A001174](https://oeis.org/A001174)** — "oriented graphs (complete digraphs) on `n` unlabeled nodes" — exactly, since a skew `{-1,0,+1}` game *is* an oriented graph (each pair is a tie/non-edge or one of two orientations). A001174 even lists the same Burnside-over-partitions formula `count_iso_classes` uses, and the asymptotic `~ 3^(n(n-1)/2)/n!` quoted above.
- `balanced` is **[A308239](https://oeis.org/A308239)** — "connected Eulerian oriented graphs" — a definitional identity (zero row sums ⇔ in-degree = out-degree ⇔ Eulerian; connected + Eulerian + `n≥2` forces paradoxical). It is offset 0, so `a(n)` is our `n`-vertex count directly (`a(8)=5274`), which **predicts the un-enumerated `balanced n=9 = 434017`** (and `n=10 = 90658149`).
- `two-paradox` rests on the tournament count **[A000568](https://oeis.org/A000568)** (`…456, 6880, 191536`), verified in the tests. Its *minimum* matches **[A362137](https://oeis.org/A362137)** `a(2)=7` (smallest 2-paradoxical tournament = Paley `Q₇`) — but that sequence is smallest *sizes*, not counts.
- `regular`, `inclusive`, and the two-paradox *counts* (`1, 4, 221`) return no OEIS match (exact 5-to-7-term searches) — they appear to be new.

The **fairness ladder** is `regular ⊂ balanced ⊂ inclusive` — strict subsets, each a stronger fairness condition. All games we enumerate are P₁ (every strategy has a beater — the `paradoxical()` baseline) and decisively connected. **Two-paradox is a separate axis, not a rung on that ladder**: it strengthens the *paradox* condition (P₁ → P₂), not the fairness condition, and is contained in none of the fairness tiers (see below).

- **`two-paradox`** = `P₂`: every *pair* of strategies has a common strict beater — a third strategy that beats both (ties don't count; `k_paradoxical(M, 2)`). A tie offers no beater, so a P₂ game is necessarily a **tie-free tournament**. The count is therefore an authoritative filter over *all* tournament classes — `search_two_paradox` enumerates every tournament up to iso via canonical augmentation (`generate_tournaments`; the leaf counts are A000568 = …456, 6880, 191536 at `n=7,8,9`) and keeps the P₂, paradoxical, connected ones — **not** a post-filter on a fairness tier. That distinction is the whole point: because P₂ is orthogonal to regular/balanced/inclusive, post-filtering a fairness tier *undercounts*. At `n=8` all 4 two-paradox games are non-regular, non-balanced, *and* non-inclusive (no fully-mixed equilibrium) — no fairness tier contains a single one of them — so the old `search_regular` post-filter reported `0` instead of `4` (no even-`n` tournament is regular: 7 games can't split into equal wins and losses), and `5` instead of `221` at `n=9` (it caught only the *regular* P₂ tournaments). The smallest P₂ game is the Paley tournament `Q₇` (`n=7`, the classical Erdős–Schütte minimum). Tournaments have no tie-twins, so the twin-free subcount always equals the total; the prime subcount can be strictly smaller (e.g. at `n=9`, 197 of the 221 are modular-prime).
- **twin-free**: no two nodes are tie-twins. Two nodes are *twins* if they tie and share the same parents and children (equivalently, identical rows in `M`); merging a twin pair gives an `(n-1)`-node game, so a structure with twins is just a smaller game with a duplicated strategy. The twin-free subcount is the "genuinely new at this `n`" sequence — e.g. Brick+Boulder+Paper+Scissors reduces to RPS, and the 5-element Water/Fire/Clay/Sand/Grass game reduces to the 4-strategy Cop/K-9/Perp/Witness game.
- **prime** (`[…]`): *modular-prime* — no nontrivial **module**. A module is a set of moves that every outside move relates to identically (beats all, loses to all, or ties all); a twin pair is just a size-2 all-tie module. A non-prime game is a smaller **quotient** game with sub-games substituted into its moves (`G = H[M₁,…,M_k]`), so prime games are the genuine irreducible atoms — a strictly stronger notion than twin-free (`prime ⊆ twin-free`). The gap first appears at `n=5` inclusive: 8 twin-free but only 7 prime, the one extra being RPS with a move blown up into another RPS — twin-free (no duplicated rows) yet decomposable. The number of extreme equilibria factorizes over the decomposition tree (`num_equilibria` ↔ `neq_tree`), which is why the `n_eq` blow-up is a property of *reducible* games: `twin-free ⟹ n_eq ≤ n` is **false** (e.g. `RPSLS[cop×5]` is twin-free with `n_eq = 2⁵ = 32 > 20`). The refined `prime ⟹ n_eq ≤ n` survives through `n=8`, but it is **also false**: the first counterexample appears at `n=9` — a *prime*, rigid (`|Aut|=1`, zero nontrivial modules) regular game with profile `(3,2,3)` and `n_eq = 11 > 9`. So the `n_eq ≤ n` bound is not implied by primeness either; the inflation just needs a rich enough indecomposable game, which first occurs at `n=9`.

Brute force (`search_inclusive`, `search_balanced`) enumerates all `3^C(n,2)` labelings, so it stops at `n=6`. Beyond that we use **isomorph-free generation** (`rpsfair/generate.py`): canonical augmentation (add one node at a time, keep only the canonically-distinguished extension) with a nauty-style refine-and-individualize canonical form, so each isomorphism class is built exactly once with bounded memory. The inclusive `n=7` count above (`10525 (9401) [9350]`) was computed this way — `search_inclusive_gen` augments the 21 480 six-node classes and keeps the paradoxical, connected, fully-mixed ones. Constraint-feasibility pruning gives `search_balanced_gen` (canonical augmentation, keeping only balance-feasible partial games) for the balanced tier at higher `n`. The regular tier instead uses the streaming row-by-row `search_regular_stream` (the augmentation overshoot-prune is too loose to be worthwhile for regularity, whereas the exact per-row `(W, T, L)` constraint prunes hard); `search_two_paradox` enumerates tournaments the same canonical-augmentation way.

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
  search.py       brute + streaming searches: search_{regular,balanced,inclusive}, *_stream
  generate.py     isomorph-free generation: generate_up_to_iso / generate_tournaments,
                  nauty canonical key, search_{inclusive,balanced}_gen, search_two_paradox
  metrics.py      num_orbits, aut_size, num_cuts, gini, tie_fraction
  plot.py         draw, grid, best_layout, add_colorbar (matplotlib)
  display.py      pretty, show, wtl_labels, letter_labels (text)
run.py            driver — loops over n in NS, emits balanced/regular/inclusive plots
view.py           CLI for inspecting a single structure
tests/            pytest sanity tests for known counts
```
