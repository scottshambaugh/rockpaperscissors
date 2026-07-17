# Rust census tools

## Install requirements

On Debian/Ubuntu (including WSL):

```sh
sudo apt install nauty libnauty-dev gcc rustc
```

- **`nauty`** — the command-line tools every pipeline here uses: `nauty-geng`,
  `nauty-directg`, `nauty-labelg`, `nauty-gentourng`, `nauty-listg`,
  `nauty-watercluster2`.
- **`libnauty-dev`** — `nauty.h` + `libnauty`, needed only to build the two
  FFI binaries (`balanced.rs` / `regular.rs`, which call densenauty through
  `balanced_shim.c`). Every other tool is plain `rustc -O <file>.rs`, no cargo.
  (Note the package is `libnauty-dev`; there is no `libnauty2-dev`.)

`rust/ci_test.sh` builds everything and checks all tools against the known
small-`n` census counts — run it after any toolchain or code change.

# Fast two-paradox (S₂) tournament counting

Independent, fast enumeration of **S₂ / "two-paradox" / Erdős–Schütte tournaments**
(tournaments in which every pair of vertices has a common dominator), used to
verify and extend the Python `search_two_paradox` counts.

## Primary path — nauty (recommended)

Rely on **nauty** for the hard part (isomorph-free generation); we only supply the
S₂ predicate. `nauty-gentourng` generates one representative per isomorphism class
of tournaments; its default output is the upper triangle row-by-row in ascii
(`1` at position `i<j` ⟺ `i` beats `j`).

```sh
rustc -O s2_filter.rs -o /tmp/s2filter
nauty-gentourng 10 2>/dev/null | /tmp/s2filter 10        # total must equal A000568(n)
```

For large `n`, split with gentourng's `res/mod` across cores and sum the shards:

```sh
for r in 0 1 2 3; do
  nauty-gentourng 11 $r/4 2>/dev/null | /tmp/s2filter 11 &
done; wait                                               # sum the per-shard S2 counts
```

`nauty-gentourng -u n` counts all tournaments (no S₂) in seconds — the A000568
checksum.

## Secondary path — self-contained generator (no nauty)

`s2_count.rs` is a from-scratch canonical-augmentation generator (color refinement
+ individualize-refine canonical form) that needs no external tools. It independently
reproduced `A000568` and the S₂ counts `1, 5, 226` for n ≤ 9, which is what
cross-validated the method — but it is far slower than nauty and stores one code
per class, so it caps at n = 10. Kept as an independent check, not the workhorse.

```sh
rustc -O s2_count.rs -o /tmp/s2 && /tmp/s2 9
```

## Results (verified, both paths agree for n ≤ 9)

| n | tournaments (A000568) | S₂ tournaments |
|---|----------------------:|---------------:|
| 7 | 456                   | 1              |
| 8 | 6880                  | 5              |
| 9 | 191536                | 226            |
| 10 | 9733056              | 29816          |
| 11 | 903753248            | 6959159        |
| 12 | 154108311168         | 2629321652     |

The S₂ count is the raw classical sequence (no extra filter). The repo's
`two-paradox` table column is the slightly smaller `S₂ ∩ paradoxical ∩ connected`
sub-count (`1, 4, 221` at n = 7,8,9).

## Regular oriented graphs — `regular_count.rs` (the fast counter)

For the `regular` strata, enumerating Eulerian orientations one at a time
(`nauty/regular.sh` via `watercluster2`) becomes hopeless for dense graphs — a
single 8-regular graph on 12 nodes has ~10⁹ Eulerian orientations (days to list).
`regular_count.rs` instead **counts** them. For each connected 2d-regular graph G
(`geng | listg -eq`), it computes the number of regular oriented graphs up to iso
as `(1/|Aut(G)|) Σ_{σ∈Aut(G)} fix(σ)`, where `fix(σ)` is the number of σ-invariant
Eulerian orientations — a degree-vector DP over the σ-edge-orbits (no enumeration).
This is the ice-type / six-vertex partition function, #P-complete in general but
fast for n≤12. `Aut(G)` is found by refinement+individualization.

```sh
rustc -O regular_count.rs -o /tmp/regcount
for d in $(seq 1 5); do nauty-geng 12 -d$((2*d)) -D$((2*d)) -c | nauty-listg -eq 2>/dev/null | /tmp/regcount $d; done
```

Validated against `search_regular` on every stratum through n=11 (e.g. n=11 d=4 =
5,104,171: `watercluster2` 80 min, counter 15 s). Gives `regular(12) = 9,348,286,118`
(d=4 stratum 7,353,314,011 in 17 min, 0.6 GB). Caveat: the Burnside loop iterates
`|Aut(G)|`, so it OOMs on complete graphs Kₙ (|Aut|=n!) — fine for n=12's strata
(max |Aut| ≈ 82,944 at K(4,4,4)) but the odd-n top stratum Kₙ must use A096368.


## Completely mixed games — `cm_filter.rs`

Counts the **completely mixed** games (Kaplansky: every Nash equilibrium fully
mixed) at odd `n` by filtering `nauty-geng -c -d2 | nauty-directg -o` (both geng
restrictions are implied by complete mixedness, so nothing is lost). The test is
exact integer arithmetic, no LP and no floats: for odd skew-symmetric `M` the
vector `v_i = (−1)^i Pf(M₋ᵢ)` of principal-minor Pfaffians satisfies `Mv = 0`
identically, and the game is completely mixed iff all `n` Pfaffians are nonzero
with alternating signs (Kaplansky 1995). Positivity of `v` subsumes paradox,
connectivity, fully-mixed existence and uniqueness in one check; Pfaffians are
computed recursively (closed-form 4×4 base case) with early exit on the first
zero or sign clash — ~4M candidates/s/core at n=9.

```sh
rustc -O cm_filter.rs -o /tmp/cmf
nauty-geng 7 2>/dev/null | nauty-directg -o 2>/dev/null | /tmp/cmf 7   # 7268
```

**Factor-critical prefilter (`factorcrit.rs`).** A completely mixed game needs all
`n` principal Pfaffians `Pf(M−i)` nonzero, and a skew Pfaffian is a signed sum over
the perfect matchings of its support graph — so `Pf(M−i) ≠ 0` forces `G−i` to have
a perfect matching. Needing that for every `i` is exactly `G` being **factor-critical**
(which also implies connected and min-degree ≥ 2). `factorcrit.rs` filters the
`graph6` stream to factor-critical graphs *before* directg's `2^edges` orientation
blow-up — lossless (re-verified: `n=7` still `7268`). Honest caveat: it only trims
the sparse tail (~10% of candidates), because the orientation count is dominated by
dense graphs, nearly all of which are factor-critical. CM's real selectivity is the
*global* one-signed-kernel condition, which no local graph property captures — so
there is no large generation-side prune, and n=9 stays a ~`4×10¹¹` enumeration made
feasible only by the cheap exact filter.

```sh
rustc -O factorcrit.rs -o /tmp/fc
nauty-geng 9 2>/dev/null | /tmp/fc 9 | nauty-directg -o 2>/dev/null | /tmp/cmf 9
```

Anchors: `n = 3, 5, 7` → `1, 7, 7268` (matches the Python `search_completely_mixed`
census), `n=8` → `0` (the parity theorem — doubles as a sign-logic check) with
candidate total `575016219 = A001174(8)` as the completeness checksum. For `n=9`,
shard geng (`res/32`) and sum; the restricted candidate total plus the complement
(`comm -23 <(nauty-geng 9 | sort) <(nauty-geng 9 -c -d2 | sort) | nauty-directg -o -u`)
must reproduce `A001174(9) = 415939243032`.

## Completely mixed games, the fast way — `cm_extend.rs` (the extension method)

`cm_filter.rs` still has to look at every one of the `~4×10¹¹` oriented 9-graphs.
`cm_extend.rs` avoids that by **building** the completely mixed games from the
`A001174(8) = 575016219` eight-vertex parents instead of filtering the nine-vertex
children. The lemma: delete a vertex from a CM game `M` and the remaining `8×8`
skew `M'` is **nonsingular** (`Pf(M')≠0`); conversely, for *any* nonsingular `M'`
and *any* new-vertex vector `r ∈ {−1,0,1}⁸`, the extended matrix has nullity
**exactly 1** with kernel `v = (−M'⁻¹r, 1)` — the row equation `rᵀv'=0` holds
automatically because `rᵀM'⁻¹r = 0` for skew `M'`. Hence

    M is completely mixed  ⇔  −M'⁻¹ r > 0 componentwise,

and that single strict-positivity condition already forces paradoxical + connected
(a zero or sign-flipped coordinate is exactly a dropped strategy or a disconnecting
even block). So we never construct a non-CM 9-graph.

Pipeline: `nauty-geng 8 | nauty-directg -o | cm_extend 9 | nauty-labelg | sort -u | wc -l`
(directg emits the 8-vertex parents; `cm_extend` emits every CM child as digraph6;
`labelg` canonicalises; `sort -u` counts iso classes). Performance came from
tackling bottlenecks in turn:

- **Dedup ~9× → 1.03×.** Each CM game is generated once per deletion-parent (all `n`
  deletions are nonsingular). A *canonical-deletion prefilter* emits a child only
  when its added vertex is maximal under a cheap degree-refinement signature —
  sound (some vertex is always maximal, so no class is lost), and the residual
  ties are mopped up by `sort -u`.
- **r-scan `3⁸` → tiny.** A bound-pruned DFS over `r` (track the partial `M'⁻¹r`,
  cut a branch once a coordinate's best case can't reach `<0`) with columns ordered
  by descending magnitude and in-place backtracking. Worst-case slice `103s → 3.3s`.
- **inverse.** Reciprocal-multiply instead of per-pivot divisions; skip
  already-reduced columns.

After that, `directg` generation is the bottleneck (~5 min single-thread, shards
`res/N` across cores). Validated to reproduce `n=5 → 7` and `n=7 → 7268` — an
*independent* algorithm from `cm_filter.rs`, so agreement cross-checks both
(`rust/ci_test.sh`).

```sh
rustc -O cm_extend.rs -o /tmp/cmx
nauty-geng 6 2>/dev/null | nauty-directg -o 2>/dev/null | /tmp/cmx 7 \
  | nauty-labelg 2>/dev/null | sort -u | wc -l          # 7268
```

## Modular-prime subcount — `prime_filter.rs`

Counts how many games in a digraph6 stream are **modular-prime** (no nontrivial
module — the `[prime]` bracket of the census columns). It reuses `balanced.rs`'s
vetted `is_prime` (the O(n²)-per-pair module-closure test) as a thin reader, so it
can prime-filter any already-selected family. Used for the completely-mixed
`[prime]` subcounts, including streaming the `n=9` CM games:

```sh
rustc -O prime_filter.rs -o /tmp/prime
nauty-geng 6 | nauty-directg -o | /tmp/cmx 7 | nauty-labelg | sort -u | /tmp/prime 7   # total=7268 prime=7240
sort -m -u ext/c_*.d6 | /tmp/prime 9                                                   # n=9 CM prime subcount
```

CM prime subcounts: `1, 6, 7240` at `n = 3, 5, 7` (matching the Python census).
For `n=9`, round-robin the merge-deduped `583591020` games to parallel `prime`
instances (`split -n r/4 --filter`) — ~4× faster, streaming, no extra disk.

## Shared module — `common.rs`

The routines duplicated across tools live in `common.rs`, pulled in by a
`mod common;` line (plain `rustc` still works — no cargo). Contents, each hoisted
verbatim from the tool where it was validated: digraph6 `decode`/`encode` and the
extension-method `inverse`/`dfs`/`added_is_maximal` (from `cm_extend.rs`), the
integer Pfaffian `pf` (from `cm_filter.rs`), the allocation-free Phase-1 LP
`fully_mixed` (from `inc_fast.rs`), and the `Arc`-bitmask tests `paradoxical`/
`connected`/`twin_free`/`is_prime` (from `balanced.rs`). `sig_maximal` generalizes
the canonical-deletion signature test with an eligibility predicate (used both for
the CM all-vertices case and `inc_extend`'s nullity-restricted case). Any edit
here is regression-checked by `ci_test.sh` — every consumer's counts must hold.

## Inclusive games, the fast filter — `inc_fast.rs`

Filters a digraph6 stream for **inclusive** games (paradoxical + connected + a
fully-mixed equilibrium exists) at odd `n`, classifying by nullity via the `n`
Pfaffian cofactors before ever touching an LP: nullity 1 ⇒ inclusive iff the
kernel `v_i = (−1)^i Pf(M₋ᵢ)` is strictly one-signed (= completely mixed, no LP);
all cofactors zero ⇒ nullity ≥ 3, only then run the Phase-1 LP (~0.6% of
candidates). Validated `n=5 → 15`, `n=7 → 10525`. This makes a directg-9 scan
(~10¹¹ candidates) feasible in ~a day of CPU — but `inc_extend.rs` below makes
that scan unnecessary.

```sh
rustc -O inc_fast.rs -o /tmp/incf
nauty-geng 7 2>/dev/null | nauty-directg -o 2>/dev/null | /tmp/incf 7   # 10525
```

## Inclusive games by extension — `inc_extend.rs`

Extends `cm_extend`'s build-don't-filter trick to the whole inclusive family.
The extension lemma, now for any parent: extending an 8-vertex parent `M'` with
kernel `K` (dim `d`) by `r ∈ {−1,0,1}⁸` gives child nullity `d+1` if `r ⊥ K`,
else `d−1`. Every nullity-`k` game has a rank-preserving deletion, so all
inclusive 9-games are reached from 8-vertex parents two ways:

- **nonsingular parent (`d=0`)**: nullity-1 children; inclusive = completely
  mixed = the `cm_extend` cone-DFS, byte-for-byte (its emissions must dedup to
  exactly `CM(9) = 583591020` — a built-in full-run checksum).
- **singular parent (`d=2` almost always)**: the nullity-3 children have `r ⊥ K`
  (2 linear equalities). Fully-mixedness collapses to **kernel coordinates**: the
  child kernel is `{(Bᵀλ + t·w, t)}` with `M'w = −r` (via a bordered-matrix
  inverse computed once per parent), so *child fully mixed ⇔ ∃λ ∈ R²:
  Bᵀλ + w > 0* — 8 half-planes whose normals are fixed per parent. By Motzkin/
  Carathéodory that system is infeasible exactly when a positive dependency of
  ≤ 3 normals has `Σαᵢwᵢ ≤ 0`, and each dependency is linear in `r` — so the
  fully-mixed conditions become extra *inequality rows in the same bound-pruned
  DFS* as the perpendicularity equalities (`dfs_fused`). The tree only ever
  reaches `r` that are already nullity-3 **and** fully mixed; the flat
  `3⁸`-with-LP scan this replaces was ~25× slower on singular-heavy input.
  Survivors get `paradox_connected` (a nullity-3 game can be fully mixed yet
  disconnected — an isolated vertex contributes `eᵢ` to the kernel), and a
  canonical-deletion prefilter restricted to nullity-`d` deletions
  (`sig_maximal`), which cuts emission redundancy ~4× so disk stays ~answer-sized.
  Rare `d ≥ 4` parents fall back to the flat DFS + LP.

Validated `n=5 → 15`, `n=7 → 10525` with every intermediate count identical
across the LP, per-candidate-Helly, and fused-DFS implementations. An optional
second argument splits the two streams (`inc_extend 9 hi.d6`: nullity-3 children
to the file, CM children to stdout) — the classes are disjoint, so they dedup
independently and `inclusive(9) = 583591020 + |dedup(hi)|`.

```sh
rustc -O inc_extend.rs -o /tmp/incx
nauty-geng 6 2>/dev/null | nauty-directg -o 2>/dev/null | /tmp/incx 7 \
  | nauty-labelg 2>/dev/null | sort -u | wc -l          # 10525
```

Measured n=9 cost (uniform samples per geng edge band): ~26–40 µs/parent → ~4–5
core-hours over all `575016219` parents, ~1.5 h wall on 4 cores + dedup; ~11 GB
of shard files (CM stream ~10 GB, nullity-3 stream ~1 GB).

## Counting without enumeration — `burnside_regular.rs` / `burnside_wide.rs`

`regular(n)` by pure counting (no games ever built): σ-invariant labeled
d-regular counts via an exchangeable-state DP, Burnside over cycle types, and
an inverse Euler transform for connectivity. Validated against every known
census value (`n = 3..12`) and the A096368 tournament strata before extending
the column: `regular(13..16) = 12050109962241, 26133517897247816,
193789800287451697002, 4477753123613209191571206` (1.4 s / 12 s / 114 s /
36 min). `burnside_wide.rs` is the 256-bit build (u128 overflows past n=15)
with a 2× memo saving from arc-reversal symmetry; both are CI-checked. The
Python prototype (`../burnside_regular.py`) is the validated reference.

## The bracket complement method — `wbal.rs` + `qfilter.rs`

The `(twin-free) [prime]` brackets for a total known in closed form (balanced =
A308239) don't need the ~98.5% of games that are twin-free and prime — only the
structured exceptions:

- **non-twin-free** balanced games are blow-ups of smaller twin-free cores with
  multiplicities `m` (weights): balanced ⇔ the core is `m`-weight-balanced
  (`Σ m_u rel(v,u) = 0`). `wbal.rs` enumerates weighted cores per weight
  partition (colored nauty canon — `rps_canon_colored`); summing gives the
  non-tf count. Anchors: 1, 3, 23, 353, 13519, 1223768 at `n = 5..10`, all exact.
- **twin-free non-prime** games need a module of size ≥ 3 (all module row sums
  are equal and sum to 0, so a 2-module with an arc is impossible → every
  2-module is a tie-twin), size ≤ 6 (`t−1 ≤ k−1` for a decisive external), and
  the quotient family is exactly the weight-`(t,1,…,1)`-balanced family
  (`wbal … q` emit mode, or `qfilter.rs` over a directg stream). Substituting
  the (tiny) rowsum-0 module families and deduping with `labelg` counts the
  prime gap: 61 / 2324 at `n = 9, 10` (both quotient routes agree), 208869 at
  `n = 11`.

Result: `balanced(11) = 48825116761 (48546325240) [48546116371]` in ~70 min
total, vs ~5 days for direct enumeration.
