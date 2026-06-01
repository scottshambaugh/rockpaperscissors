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

