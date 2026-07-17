## See *[The Secret Garden of Rock-Paper-Scissors](https://theshamblog.com/the-secret-garden-of-rock-paper-scissors/)* for context, details, and results

This repo was pretty much entirely vibe coded useing Claude Code. I understood and could verify the initial implementation on low n, but was hands off for the optimizations that were necessary for searching higher n, beyond verifying that the algorithms still reproduced the values we already found. All this to say, take the results here with a grain of salt and I think it would be an interesting exercise to implement this yourself and figure out optimizations for faster game searches.

----

# Fair generalizations of Rock-Paper-Scissors

## Contents

- [Overview](#overview)
- [Fairness definitions](#fairness-definitions)
- [Counts up to isomorphism](#counts-up-to-isomorphism)
- [How the counts were computed](#how-the-counts-were-computed)
- [OEIS cross-references](#oeis-cross-references)
- [Equilibrium structure](#equilibrium-structure)
- [Usage](#usage)

## Overview

Enumerate and visualize "fair" generalizations of Rock-Paper-Scissors on `n` strategies, represented as `n×n` skew-symmetric matrices with entries in `{-1, 0, +1}`.

Metrics computed per structure: `orbits` (structurally distinct strategies = automorphism orbits; two moves share an orbit iff some relabeling swaps them without changing the game), `|Aut|` (automorphism group size), articulation cuts, Gini of the equilibrium, tie fraction (over all `n²` matchups including the always-tie diagonal, so a tournament with no draws is `1/n`).

## Fairness definitions

### The fairness ladder

Three nested notions of fairness:

- **Regular** — every node has the **same** `(W, T, L)` profile. Strongest.
- **Balanced** — every row sums to zero (uniform equilibrium); profiles may differ.
- **Inclusive** — *some* equilibrium is fully supported (every strategy is playable). Weakest.

`regular ⊂ balanced ⊂ inclusive`. All three require the game to be paradoxical (every strategy wins and loses at least once) and decisively connected. Structures are counted up to graph isomorphism.

### Completely mixed

A fourth notion strengthens inclusive along the *equilibrium* axis rather than the profile axis:

- **Completely mixed** (Kaplansky's term) — *every* equilibrium is fully supported, so every strategy is not just playable but **required** (played in every equilibrium). This forces the equilibrium to be unique and `n` to be odd; `completely-mixed ⊂ inclusive`, but it is incomparable with regular/balanced (no even-`n` game qualifies, and some odd-`n` regular games fail it).

### The three axes

The **fairness ladder** is `regular ⊂ balanced ⊂ inclusive` — strict subsets, each a stronger fairness condition. All games we enumerate are P₁ (every strategy has a beater — the `paradoxical()` baseline) and decisively connected.

**Two-paradox is a separate axis, not a rung on that ladder**: it strengthens the *paradox* condition (P₁ → P₂), not the fairness condition, and is contained in none of the fairness tiers (see below).

**Completely mixed is a third axis**: it strengthens the *equilibrium* condition (some fully-mixed equilibrium → every equilibrium fully mixed), sits strictly inside inclusive, and cuts across regular/balanced rather than containing or being contained by them.

### Subcount notions: twin-free and prime

Each fairness cell in the table below carries two subcounts. Their definitions:

- **twin-free**: no two nodes are tie-twins. Two nodes are *twins* if they tie and share the same parents and children (equivalently, identical rows in `M`); merging a twin pair gives an `(n-1)`-node game, so a structure with twins is just a smaller game with a duplicated strategy. The twin-free subcount is the "genuinely new at this `n`" sequence — e.g. Brick+Boulder+Paper+Scissors reduces to RPS, and the 5-element Water/Fire/Clay/Sand/Grass game reduces to the 4-strategy Cop/K-9/Perp/Witness game.
- **prime** (`[…]`): *modular-prime* — no nontrivial **module**. A module is a set of moves that every outside move relates to identically (beats all, loses to all, or ties all); a twin pair is just a size-2 all-tie module. A non-prime game is a smaller **quotient** game with sub-games substituted into its moves (`G = H[M₁,…,M_k]`), so prime games are the genuine irreducible atoms — a strictly stronger notion than twin-free (`prime ⊆ twin-free`). The gap first appears at `n=5` inclusive: 8 twin-free but only 7 prime, the one extra being RPS with a move blown up into another RPS — twin-free (no duplicated rows) yet decomposable. The number of extreme equilibria factorizes over the decomposition tree (`num_equilibria` ↔ `neq_tree`), which is why the `n_eq` blow-up is a property of *reducible* games: `twin-free ⟹ n_eq ≤ n` is **false** (e.g. `RPSLS[cop×5]` is twin-free with `n_eq = 2⁵ = 32 > 20`). The refined `prime ⟹ n_eq ≤ n` survives through `n=8`, but it is **also false**: the first counterexample appears at `n=9` — a *prime*, rigid (`|Aut|=1`, zero nontrivial modules) regular game with profile `(3,2,3)` and `n_eq = 11 > 9`. So the `n_eq ≤ n` bound is not implied by primeness either; the inflation just needs a rich enough indecomposable game, which first occurs at `n=9`.

## Counts up to isomorphism

Each fairness cell is `total (twin-free) [prime]` — the total count up to isomorphism, the twin-free subcount in parentheses, and the modular-prime subcount in brackets. **Past the brute-force range the `(twin-free) [prime]` breakdown isn't computed, so those cells show the total only** (e.g. `n ≥ 12` balanced). `—` means not reached.

The `candidates` column is the raw brute-force search space, `3^(n choose 2)` (one of `{-1, 0, +1}` per upper-triangle edge), before any filtering. The `iso classes` column collapses that by relabeling: the number of distinct games up to isomorphism (`candidates` modulo the `S_n` action, computed in closed form by Burnside over the cycle types — `count_iso_classes(n)`), the real size of the space the fair games sit inside. It matches `1, 2, 7, 42, 582, 21480, …` and stays comparatively tame (≈ `4.2 × 10¹¹` at `n=9`) where `candidates` has already passed `10¹⁷`. It is *not* simply `3^C(n,2) / n!`: that ratio is only the identity-permutation term of the Burnside sum (and isn't even an integer), so it undercounts by the automorphism correction — slightly, by `1.032 / 1.014 / 1.006` at `n = 7 / 8 / 9` (vanishing as games become rigid).

The `two-paradox` column is the **raw S₂** tournament count (`A000568`-restricted; the paradoxical+connected sub-count is the slightly smaller `1, 4, 221` at `n=7,8,9` — see the OEIS notes). Counts past the brute-force wall come from nauty + the tools in [`nauty/`](nauty/) and [`rust/`](rust/), each anchored to a known checksum (`A001174` / `A000568` / `A308239`).

| n  | candidates = 3^C(n,2)            | iso classes             | two-paradox | regular                  | balanced                 | inclusive                   | completely mixed     |
|----|----------------------------------|-------------------------|-------------|--------------------------|--------------------------|-----------------------------|----------------------|
| 3  | 27                               | 7                       | 0           | 1 (1) [1]                | 1 (1) [1]                | 1 (1) [1]                   | 1 (1) [1]            |
| 4  | 729                              | 42                      | 0           | 1 (1) [1]                | 1 (1) [1]                | 3 (2) [2]                   | 0                    |
| 5  | 59049                            | 582                     | 0           | 2 (2) [2]                | 4 (3) [3]                | 15 (8) [7]                  | 7 (7) [6]            |
| 6  | 14348907                         | 21480                   | 0           | 5 (4) [4]                | 16 (13) [13]             | 222 (177) [169]             | 0                    |
| 7  | 10460353203                      | 2142288                 | 1           | 13 (12) [12]             | 175 (152) [152]          | 10525 (9401) [9350]         | 7268 (7268) [7240]   |
| 8  | 22876792454961                   | 575016219               | 5           | 82 (76) [76]             | 5274 (4921) [4917]       | 1198013 (1128896) [1127592] | 0                    |
| 9  | 150094635296999121               | 415939243032            | 226         | 2016 (1973) [1972]       | 434017 (420498) [420437] | 654809211 (645515521) [645436603] | 583591020 (583591020) [583533123] |
| 10 | 2954312706550833698643           | 816007449011040         | 29816       | 154831 (153529) [153529] | 90658149 (89434381) [89432057] | —                           | 0                    |
| 11 | 174449211009120179071170507      | 4374406209970747314     | 6959159     | 21171976 (21091698) [21091694] | 48825116761 (48546325240) [48546116371] | —                           | —                    |
| 12 | 30903154382632612361920641803529 | 64539836938720749739356 | 2629321652  | 9348286118 (9332883819) [9332883596] | 68579602126387           | —                           | 0                    |
| 13 | 16423203268260658146231467800709255289 | 2637796735571225009053373136 | — | 12050109962241 | 254236892409251687 | — | — |
| 14 | 26183890704263137277674192438430182020124347 | 300365896158980530053498490893399 | — | 26133517897247816 | 2512196983110257299898 | — | 0 |
| 15 | 125236737537878753441860054533045969266612127846243 | 95772686348981154949990837579207581003 | — | 193789800287451697002 | 66738846701458593981832314 | — | — |
| 16 | 1797010299914431210413179829509605039731475627537851106401 | 85888418249161165461289338402279023247144414 | — | 4477753123613209191571206 | 4802164152196181961357340585741 | — | 0 |

## How the counts were computed

### Isomorph-free generation

Brute force (`search_inclusive`, `search_balanced`) enumerates all `3^C(n,2)` labelings, so it stops at `n=6`. Beyond that we use **isomorph-free generation** (`rpsfair/generate.py`): canonical augmentation (add one node at a time, keep only the canonically-distinguished extension) with a nauty-style refine-and-individualize canonical form, so each isomorphism class is built exactly once with bounded memory.

The inclusive `n=7` count above (`10525 (9401) [9350]`) was computed this way — `search_inclusive_gen` augments the 21 480 six-node classes and keeps the paradoxical, connected, fully-mixed ones. Constraint-feasibility pruning gives `search_balanced_gen` (canonical augmentation, keeping only balance-feasible partial games) for the balanced tier at higher `n`. The regular tier instead uses the streaming row-by-row `search_regular_stream` (the augmentation overshoot-prune is too loose to be worthwhile for regularity, whereas the exact per-row `(W, T, L)` constraint prunes hard); `search_two_paradox` enumerates tournaments the same canonical-augmentation way.

### regular

The `n = 13..16` regular values come from the counting method in `burnside_regular.py` / `rust/burnside_regular.rs` — **pure counting**, no game enumeration at all: labeled σ-invariant d-regular counts via an exchangeable-state DP (processing vertices one at a time, the unprocessed ones are interchangeable, so the state collapses to the multiset of their partial `(out, in)` tallies; σ's vertex cycles enter as super-vertices whose pair-orbits contribute in `lcm`-sized degree jumps, with antipodal orbits of even cycles forced to tie), Burnside over the 101 cycle types of `S₁₃`, and an inverse Euler transform for connectivity (components of a `d`-regular game are `d`-regular).

It reproduces every known value `n = 3..12` — including `regular(12) = 9348286118` and the per-stratum anchors `A096368` (1223 regular tournaments at `n=11`, 1495297 at `n=13`) — before emitting `regular(13) = 12050109962241` (strata `d=1..6`: 1, 7371469, 278601715904, 11545150361450, 226349018120, 1495297). The Rust port reproduces the Python `n=13` output exactly and runs it in 1.4 s (Python took ~90 min), giving `regular(14) = 26133517897247816` in 12 s, `regular(15) = 193789800287451697002` in 114 s, and `regular(16) = 4477753123613209191571206` in 36 min; the wide build uses checked `U256` arithmetic. The `n=15` top stratum independently reproduces `A096368`'s 18400989629 regular tournaments.

Within the brute-force-plus-nauty range, `regular` reaches `n=12` via an exact **counter** (`rust/regular_count.rs`: Burnside over each graph's automorphisms + a degree-vector DP that *counts* Eulerian orientations instead of enumerating them — its dense `d=4` stratum, 7353314011, lands in 17 min vs a multi-day enumeration).

### balanced

`balanced` is closed (`A308239`) and `iso classes` is closed-form (`A001174`); `balanced(13..16)` are A308239's b-file values (its `n = 10..12` terms match our enumerations).

### inclusive

`inclusive` reaches `n=9` by the same extension trick as completely mixed, generalized to singular parents (`rust/inc_extend.rs`): a parent with kernel `K` extended by `r ⊥ K` gains a kernel dimension, so the nullity-3 inclusive games are built from the singular 8-vertex parents — with the fully-mixed condition reduced to kernel coordinates (an ∃λ: `Bᵀλ + w > 0` feasibility whose Helly/Motzkin certificates are *linear in `r`*, hence fused straight into the bound-pruned DFS as extra inequality rows) — while the nullity-1 side is exactly the `cm_extend` cone (its 583591020 dedup total doubling as the run's checksum); a filter-everything `n=9` scan would instead need the LP on all `A001174(9) ≈ 4.2 × 10¹¹` oriented graphs.

`inclusive(10)` stays out of reach for a different reason: the answer itself is likely `~10¹¹–10¹²` classes (multi-TB before dedup), and its nullity-1 parents are generic 9-vertex games, so no small parent family exists.

### completely mixed

**`completely mixed`** (Kaplansky 1945): every Nash equilibrium is fully mixed — the standard game-theory name for exactly this property. Kaplansky proved such a game has a *unique* optimal strategy, and for skew-symmetric games it forces odd `n` (even `n` gives an even rank, hence a fully-mixed *continuum* whose boundary points drop strategies — the parity theorem again, which is why the even-`n` cells are `0` by theorem, no search needed).

The characterisation is Kaplansky's Theorem 5 (1945): a value-zero game is completely mixed iff its matrix is square of rank `n−1` with all cofactors nonzero and of one sign — exactly `dim ker(M) = 1` with a strictly one-signed kernel vector, which is what `is_completely_mixed` tests via a single SVD. Kaplansky already writes out the skew case explicitly there — the parity theorem ("an `n×n` symmetric game can be completely mixed only if `n` is odd"), the `n=3` game (which he calls "a sort of generalization of Stone, Paper and Scissors"), and the five sign-expressions for `n=5` that are precisely the principal Pfaffians `Pf(M₋ᵢ)` — so `rust/cm_filter.rs` is running his own 1945 formula; the clean general-`n` "principal Pfaffians alternate in sign" statement is the 1995 sequel ([Ann. of Math. 46 (1945) 474–479](https://www.jstor.org/stable/1969164)).

Reaching `n=9` exploits that the membership test is exact integer arithmetic: for odd skew-symmetric `M` the vector `v_i = (−1)^i Pf(M₋ᵢ)` of principal-minor Pfaffians always satisfies `Mv = 0`, and the game is completely mixed iff all `n` Pfaffians are nonzero with alternating signs (Kaplansky 1995) — positivity of `v` subsuming paradox, connectivity, existence *and* uniqueness in one test. Two independent routes compute `n=9` and agree (see [`rust/`](rust/)).

- **(1) Direct filter** `rust/cm_filter.rs` runs that test over the whole `nauty-directg -o` stream (~4M candidates/s/core; a ~1-day 4-core sweep).
- **(2) Extension method** `rust/cm_extend.rs`, far faster: it *builds* the CM games from the `A001174(8) = 575016219` eight-vertex parents rather than filtering all `~4×10¹¹` nine-graphs — deleting a vertex from a CM game leaves a **nonsingular** `M'`, and conversely *any* `r ∈ {−1,0,1}⁸` extends it to nullity exactly 1 with kernel `(−M'⁻¹r, 1)`, so `M` is completely mixed iff `−M'⁻¹r > 0`, and that single condition already forces paradoxical + connected. A bound-pruned DFS over `r`, a canonical-deletion prefilter (emit only when the added vertex is signature-maximal, cutting the `~n`-fold deletion redundancy to `1.002×`), and `nauty-labelg | sort -u` for the final dedup bring it down to minutes; `directg` generation becomes the bottleneck.

Twins always break the property (a tie-twin pair `i, j` puts `e_i − e_j` in `ker(M)`), so **completely mixed ⟹ twin-free** and the twin-free subcount equals the total; the modular-prime subcounts are `1, 6, 7240, 583533123` at `n = 3, 5, 7, 9` — only `57897` of the `583591020` at `n=9` are non-prime, so `99.99%` are irreducible atoms, the prime fraction climbing toward `1` as games grow (`86% → 99.6% → 99.99%`). The `n=9` prime subcount was got by streaming the enumerated CM games (`rust/prime_filter.rs`, reusing `balanced.rs`'s `is_prime`) round-robin to parallel filters.

Cross-axis intersections at odd `n`: regular ∩ completely-mixed = `1/1, 2/2, 10/13, 1381/2016` and balanced ∩ completely-mixed = `1/1, 3/4, 102/175, 312652/434017` of the tier at `n = 3, 5, 7, 9` (matching the `#equilibria` table below: the `{1: k}` cells are the completely mixed ones). Getting *really* strict — adding rigidity (`|Aut| = 1`, every vertex structurally unique) — leaves `0, 3, 7146` at `n = 3, 5, 7`, and `3, 7143` of those (n = 5, 7) are also prime: fully fair in the equilibrium sense, yet with zero symmetry and no modular shortcut. The odd-`n` count sequence is `1, 7, 7268, 583591020` for `n = 3, 5, 7, 9` — the `n=9` term explodes because completely-mixedness stops being rare once games are large (`583591020` is ≈`0.14%` of `A001174(9)`).

*Validation.* Both routes anchor at `n = 3, 5, 7` to the Python census and at `n=8` to the parity-theorem `0` (with the `A001174(8)` candidate checksum); the `n=9` value `583591020` was cross-checked both directions — every extension-emitted game re-passes the integer-Pfaffian filter (no false positives over 60.6M sampled) and an independent `directg`-9 slice enumeration is fully contained in the extension output (no misses). See [`rust/`](rust/).

### two-paradox

**`two-paradox`** = classical tournament `P₂`: the game is a **tie-free tournament**, and every *pair* of strategies has a common strict beater — a third strategy that beats both (`k_paradoxical(M, 2)`). Tie-freeness is part of the classical definition, not a consequence of the common-beater condition alone. The count is therefore an authoritative filter over *all* tournament classes — `search_two_paradox` enumerates every tournament up to iso via canonical augmentation (`generate_tournaments`; the leaf counts are A000568 = …456, 6880, 191536 at `n=7,8,9`) and keeps the P₂, paradoxical, connected ones — **not** a post-filter on a fairness tier.

That distinction is the whole point: because P₂ is orthogonal to regular/balanced/inclusive, post-filtering a fairness tier *undercounts*. At `n=8` all 4 two-paradox games are non-regular, non-balanced, *and* non-inclusive (no fully-mixed equilibrium) — no fairness tier contains a single one of them — so the old `search_regular` post-filter reported `0` instead of `4` (no even-`n` tournament is regular: 7 games can't split into equal wins and losses), and `5` instead of `221` at `n=9` (it caught only the *regular* P₂ tournaments). The smallest P₂ game is the Paley tournament `Q₇` (`n=7`, the classical Erdős–Schütte minimum). Tournaments have no tie-twins, so the twin-free subcount always equals the total; the prime subcount can be strictly smaller (e.g. at `n=9`, 197 of the 221 are modular-prime).

`two-paradox` reaches `n=12` (8 `gentourng` shards summing to `A000568(12) = 154108311168`, ~12h wall).

## OEIS cross-references

Several columns connect to known sequences, which both validates the code and extends the table:

- `iso classes` is **[A001174](https://oeis.org/A001174)** — "oriented graphs (complete digraphs) on `n` unlabeled nodes" — exactly, since a skew `{-1,0,+1}` game *is* an oriented graph (each pair is a tie/non-edge or one of two orientations). A001174 even lists the same Burnside-over-partitions formula `count_iso_classes` uses, and the asymptotic `~ 3^(n(n-1)/2)/n!` quoted above.
- `balanced` is **[A308239](https://oeis.org/A308239)** — "connected Eulerian oriented graphs" — a definitional identity (zero row sums ⇔ in-degree = out-degree ⇔ Eulerian; connected + Eulerian + `n≥2` forces paradoxical). It is offset 0, so `a(n)` is our `n`-vertex count directly (`a(8)=5274`), which **predicts the un-enumerated `balanced n=9 = 434017`** (and `n=10 = 90658149`).
- `two-paradox` rests on the tournament count **[A000568](https://oeis.org/A000568)** (`…456, 6880, 191536`), verified in the tests; its *minimum* matches **[A362137](https://oeis.org/A362137)** `a(2)=7` (smallest 2-paradoxical tournament = Paley `Q₇`) — but A362137 is smallest *sizes*, not counts. The **count** of S₂ (two-paradox / Erdős–Schütte) tournaments per `n` is **not in OEIS**: it is `1, 5, 226, 29816, 6959159, 2629321652` at `n = 7..12` (the raw classical count, shown in the table; restricting to paradoxical+connected gives the slightly smaller `1, 4, 221` at `n=7,8,9`). Verified two independent ways — nauty (`nauty-gentourng | rust/s2_filter`) and a from-scratch generator (`rust/s2_count`), which agree through n=9 and both reproduce A000568 (the `n = 10,11,12` terms via nauty, each with the A000568 total as a completeness checksum — `A000568(12) = 154108311168`, matched exactly); see [`rust/`](rust/). This is the strongest OEIS-submission candidate: a count sequence that slots directly beside the famous Erdős–Schütte minimum.
- `regular` (= connected oriented graphs with out-degree = in-degree = `d` at every vertex, summed over `d`) returns no OEIS match, but **stratifies** into known sequences by the common win-degree `d`: `regular(n) = Σ_{d=1}^{⌊(n-1)/2⌋} R(n,d)` with `R(n,1) = 1` (the directed `n`-cycle), `R(n,2) = ` **[A219894](https://oeis.org/A219894)** ("directed 2-regular graphs without 2-cycles", `1,4,9,55,453,…`), and the top stratum `T=0` (odd `n`) `= ` **[A096368](https://oeis.org/A096368)** regular tournaments (`1,1,3,15,…`). The *sum* over `d` is the new part (and has no closed form — it contains the regular-tournament enumeration, which itself has none). Extended past the Python search with nauty (`nauty/regular.sh`, via `geng | watercluster2`, validated against `search_regular` for n ≤ 9): the sequence is `1, 1, 2, 5, 13, 82, 2016, 154831, 21171976, 9348286118` for n = 3..12 (the n=11 strata are `d=1..5: 1, 47594, 16018987, 5104171, 1223`, the last being A096368's regular tournaments on 11 nodes). For the dense high-`n` strata the `watercluster2` *enumeration* becomes hopeless (one 8-regular graph on 12 nodes has ~10⁹ Eulerian orientations, days to list), so those are computed by an exact **counter** instead — `rust/regular_count.rs` sums `(1/|Aut(G)|) Σ_σ fix(σ)` over each graph, with `fix(σ)` a degree-vector DP that counts σ-invariant Eulerian orientations (the ice-type / six-vertex partition function) without listing them. Validated against `search_regular` on every stratum through n=11 (e.g. the n=11 `d=4` value 5104171, which took `watercluster2` 80 min, the counter gets in 15 s).
- `inclusive` (`1, 3, 15, 222, 10525, 1198013, 654809211` for n = 3..9) returns no OEIS match — it too appears to be new. Independently verified *and extended* via nauty: `nauty-geng n | nauty-directg -o` enumerates all oriented graphs (= A001174), then `rust/inclusive` keeps the paradoxical, connected ones with a fully-mixed equilibrium (a Phase-1 LP feasibility test). This reproduces `search_inclusive` (`3, 15, 222, 10525` for n = 4..7) through a different generation path, then reaches **`n=8 = 1198013`** (575016219 candidates = A001174(8), checksum-matched). The `n=9` term comes from the **extension method** (`rust/inc_extend.rs`, see below and [`rust/`](rust/)): `inclusive(9) = 654809211 = 583591020 + 71218191`, the completely mixed games plus the nullity-3 fully-mixed ones — the two nullity classes are built and deduplicated separately, and the CM side reproducing the independent `cm_extend` census exactly is the run's completeness checksum.

## Equilibrium structure

### Why full support forces uniqueness

The equivalence of "completely mixed" with "unique equilibrium" is two lines — if `p > 0` spans `ker(M)` and `q` is any equilibrium (`Mq ≤ 0`), then `pᵀ(Mq) = -(Mp)ᵀq = 0` forces `Mq = 0`, so `q ∈ ker(M)` and `q = p`.

**Geometrically, this is why full support *forces* uniqueness** — the counterintuitive core of the notion. The symmetric equilibria form a convex polytope, and for these games the fully-mixed ones fill `ker(M) ∩ Δ`; so two distinct fully-mixed equilibria would span a segment lying inside the flat `ker(M)`, and sliding along that segment until a coordinate first hits `0` reaches an equilibrium that *drops* a strategy. Hence more than one equilibrium always forces a zero-support one, and `dim ker(M) = 1` — a single kernel ray, meeting the simplex in exactly one point — is the *only* way every equilibrium can stay interior.

The elemental game is this failure made visible: a 4-corner equilibrium polytope whose corners zero out Water/Fire/Clay/Grass, so it is inclusive but not completely mixed (`equilibria_grid` renders it). The equivalence is verified independently against `num_equilibria == 1` on all 10766 enumerated inclusive games (n ≤ 7).

**Required vs. playable.** Membership means every strategy is **required** (played in *every* equilibrium, `required_strategies`), not merely *playable* (played in *some* equilibrium — which is all inclusive promises: the union of equilibrium supports is everything). The gap is real and can be total: one twin-free inclusive `n=5` game has 3 extreme equilibria whose required set is **empty** — every single strategy is droppable. The 5-move elemental game shows the subtler version: its Grass/Clay twin pair makes each of the two individually droppable (only their *sum* is pinned), so it is inclusive but not completely mixed.

### Number of equilibria

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

A cell `{k:m}` reads "*m* twin-free games have *k* extreme equilibria". Every category in this table is **inclusive**, so it has a fully-mixed equilibrium and therefore `O = ker(M) ∩ Δ`. Hence `dim O = dim ker(M) − 1 ≡ n−1 (mod 2)`: within these categories a unique equilibrium is possible only at **odd `n`**, and every even-`n` game has a continuum of equilibria. This parity conclusion does not extend to arbitrary games: an even-order game with only boundary equilibria can have a unique equilibrium.

In general the equilibrium condition is `M p ≤ 0` (every pure reply scores `≤ 0`), *not* `M p = 0`. A fully-mixed equilibrium is itself forced into `ker(M)`; without one, a partial-support equilibrium can sit on a face where an unplayed move's reply is strictly negative, and such a vertex is **not** a kernel vector. For an arbitrary game you therefore intersect `M p ≤ 0` with the simplex and enumerate vertices (what `equilibrium_vertices` does), rather than just taking `ker(M)`. If a fully-mixed equilibrium exists, however, its positivity forces every equilibrium into the kernel, so `O = ker(M) ∩ Δ`; that theorem applies to all inclusive games in this table. Paradoxicality and connectivity alone do not imply it.

## Usage

### Install

```
uv sync
```

### Run

```
uv run run.py
```

For each `n` in `NS` (default `[3, 4, 5, 6]` in `run.py`): prints the counts table, then for each fairness kind ranks the structures (decomposable filtered), prints the labeled upper-triangular grid of the top structure, and renders the ranking to `plots/n{n}_{regular,balanced,inclusive}.png`.

### Test

```
uv run pytest
```

### Lint / format

```
uv run ruff check .       # report issues
uv run ruff check --fix . # auto-fix
uv run ruff format .      # apply formatter
```

### Inspect a single structure

`view.py` is a CLI for one structure — prints the labeled upper-triangular grid and saves a plot.

```
uv run view.py --n 3                                  # default: inclusive, index 0
uv run view.py --n 3 --labels Rock,Paper,Scissors     # custom node names
uv run view.py --n 5 --kind balanced --index 2
uv run view.py --n 4 --index 1 --save my.png
uv run view.py --n 4 --index 1 --no-plot              # text only
```

### Labeling nodes

By default, nodes are labeled with their WTL profile tuple (e.g. `1·0·2` = 1 win, 0 ties, 2 losses). Override with a custom list to give game-specific names:

```python
from rpsfair import pretty, search_inclusive

M, xs = search_inclusive(3)[0]
print(pretty(M, labels=["Rock", "Paper", "Scissors"]))
```

`pretty(M, labels)` returns a labeled upper-triangular text grid using `+` (row beats col), `-` (row loses), `0` (tie). `show(M, labels)` prints it directly. `letter_labels(n)` yields `['A', 'B', ...]` if you want letter labels.

### Color scale

Plots use a global colormap (`viridis`) normalized to `0%`–`50%` equilibrium play rate, so fills are directly comparable across plots and across different `n`. `grid()` and `view.py` include a horizontal colorbar; `add_colorbar(fig, ...)` is exposed for custom figures.

### Caching

Search results are written to `cache/<name>.json` with human-readable filenames (`regular_n5.json`, `balanced_n5.json`, `inclusive_n5.json`, …). The filename does **not** encode the source — if you edit a search function, wipe the cache (`rm -rf cache/`) before re-running.

### Long enumerations

The searches can run for minutes to hours. **Always wrap a long enumeration in a [`tqdm`](https://github.com/tqdm/tqdm) progress bar** so a run is never a silent black box — `search_balanced_gen`, `search_inclusive_gen`, and `search_two_paradox` already do. Two rules that make the bar actually useful:

- **Give it a `total` whenever one is knowable.** A pruned canonical augmentation has no cheap a-priori estimate of its own work, but the *final* count per tier is known (the OEIS / census tables in `generate.py`: `_BALANCED_COUNTS`, `_INCLUSIVE_COUNTS`, `_TOURNAMENT_COUNTS`, `_ISO_COUNTS`). Counting *kept* results toward that known total yields a true percentage and ETA. With no total, tqdm still shows count + rate + elapsed — which beats nothing — but you get no ETA.
- **tqdm writes to stderr, so don't discard stderr.** A harness that runs a search as `... > out.txt 2>/dev/null` throws the bar away — the work still runs, you just go blind (this is exactly how an early balanced n=9 run lost its bar). Keep `2>&1`, or send the bar to a separate file, when you want to watch progress.

For the native Rust counters (`rust/`), the equivalent is the `eprint!("\r…")` progress line — and it must be followed by `io::stderr().flush()` or it won't appear when stderr is redirected to a file.

### Layout

```
rpsfair/
  cache.py        JSON disk cache
  structure.py    paradoxical / connected / regular / canonical_key / orbit_bytes
  equilibrium.py  equilibrium, has_fully_mixed (SVD null-space + LP),
                  is_completely_mixed, required_strategies
  search.py       brute + streaming searches: search_{regular,balanced,inclusive}, *_stream,
                  search_completely_mixed (nullity-1 filter over inclusive)
  generate.py     isomorph-free generation: generate_up_to_iso / generate_tournaments,
                  nauty canonical key, search_{inclusive,balanced}_gen, search_two_paradox
  metrics.py      num_orbits, aut_size, num_cuts, gini, tie_fraction
  plot.py         draw, grid, best_layout, add_colorbar (matplotlib)
  display.py      pretty, show, wtl_labels, letter_labels (text)
run.py            driver — loops over n in NS, emits balanced/regular/inclusive plots
view.py           CLI for inspecting a single structure
tests/            pytest sanity tests for known counts
```
