# High-n enumeration via nauty

For sizes past where the self-contained Python (`rpsfair`) is practical, we lean on
**nauty/Traces** (McKay & Piperno) - the standard isomorph-free generators - and
supply only the property predicate. nauty does the hard, trustworthy part
(generating one representative per isomorphism class); every pipeline below
re-derives a known total (A000568 / A096368 strata) as a built-in checksum.

Install (Debian/Ubuntu): `sudo apt install nauty` (binaries are `nauty-*`).

## Two-paradox (S₂ / Erdős–Schütte) tournaments - `two_paradox.sh`

```sh
nauty/two_paradox.sh 11        # gentourng | rust/s2_filter, sharded over cores
```

`nauty-gentourng` generates tournaments up to iso; `rust/s2_filter` keeps those
where every pair of vertices has a common dominator. The per-shard *total* sums to
A000568(n) - the completeness checksum. Do **not** add `gentourng -c`: S₂ does not
imply strong connectivity (5 of the 226 at n=9 are not strongly connected).

| n  | 7 | 8 | 9   | 10    | 11      | 12         |
|----|---|---|-----|-------|---------|------------|
| S₂ | 1 | 5 | 226 | 29816 | 6959159 | 2629321652 |

(Also reproduced by the self-contained `rust/s2_count.rs` for n ≤ 9. The n=12 run
is 8 shards summing to `A000568(12) = 154108311168`, ~12 h wall on 4 cores.)

## Regular oriented graphs - `regular.sh`

```sh
nauty/regular.sh 10
```

`regular(n) = Σ_d` (connected oriented graphs, every vertex out-deg = in-deg = d).
For each `d`: `nauty-geng n {2d}-regular -c | nauty-watercluster2 i$d o$d S` - the
decisive edges are a connected 2d-regular graph, and `i d o d` forces the Eulerian
orientation out=in=d. `watercluster2` emits one representative per iso class.
Validated against the Python `search_regular`: `2, 5, 13, 82, 2016` at n = 5..9
(matching per-degree, e.g. n=9 = `d1:1 d2:453 d3:1547 d4:15`), then extended to
`regular(10) = 154831`. The top stratum (`d=(n-1)/2`, odd n) is the regular
tournaments, [A096368](https://oeis.org/A096368).

## Inclusive oriented games - candidate generation

`nauty-geng n | nauty-directg -o` generates every oriented graph up to iso
(= A001174). `rust/inclusive` keeps the paradoxical + connected ones and selects
those with a fully-mixed equilibrium (ker(M) meets the positive orthant - a
Phase-1 LP feasibility test, done in-filter so it scales to n=8's ~10^8 candidates):

```sh
rustc -O rust/inclusive.rs -o /tmp/inc
nauty-geng 8 2>/dev/null | nauty-directg -o -T 2>/dev/null | /tmp/inc 8   # -> 1198013
```

Reproduces `search_inclusive` = `3, 15, 222, 10525` for n = 4..7 (a different
generation path), then reaches **n=8 = 1198013** (575,016,219 candidates =
A001174(8), checksum-matched; sharded over cores with `geng res/mod`).

(`rust/inc_emit` is the earlier variant that emits survivor codes for an external
Python `has_fully_mixed` check - fine up to n=7, too slow for n=8.)
