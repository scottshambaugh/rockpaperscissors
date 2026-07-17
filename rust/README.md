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
