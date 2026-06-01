# Converting the hand-rolled iso machinery to nauty

`rpsfair` reimplements, by hand, several things nauty does natively and faster.
Each skew `{-1,0,+1}` game `M` *is* an oriented graph (arc `i→j` iff `M[i,j]==1`;
ties are non-edges), so nauty's **directed-graph** mode applies directly.

**Validated representation** (matches our code exactly — `|Aut|` and canonical
certificate agree on RPS/RPSLS/COP/Paley₇/ring₆ and all 582 n=5 classes):

```python
import pynauty
def to_pynauty(M):
    n = len(M)
    adj = {i: [j for j in range(n) if M[i, j] == 1] for i in range(n)}
    return pynauty.Graph(n, directed=True, adjacency_dict=adj)
# canonical key:        pynauty.certificate(to_pynauty(M))
# automorphism group:   pynauty.autgrp(to_pynauty(M))  -> (gens, grpsize, ...)
```

## What converts

| our function (file) | does | nauty equivalent |
|---|---|---|
| `canonical_key` (structure.py) | brute `n!` canonical form | `pynauty.certificate` (in-process) / `dreadnaut` `d`+canon |
| `canonical_key_fast`, `_nauty_search`, `_canon_and_autos`, `_refine`, `_wl_colors` (generate.py) | our refine+individualize canonical form / autos | `pynauty.certificate` + `pynauty.autgrp` |
| `automorphisms`, `aut_size` (metrics.py) | automorphism group / its size | `pynauty.autgrp` (returns generators + group size) |
| `num_orbits`, `node_orbits` (metrics.py) | vertex orbits under Aut | `pynauty.autgrp` returns orbits directly |
| `orbit_bytes`, `orbit_hashes` (structure.py) | relabeling-orbit dedup keys | subsumed by `certificate` |
| `generate_tournaments` (generate.py) | iso-free tournaments | `nauty-gentourng` (already used in `nauty/`) |
| `generate_up_to_iso` (generate.py) | iso-free oriented graphs | `nauty-geng \| nauty-directg -o` (already used for inclusive) |
| `rust/s2_count.rs`, `rust/regular_count.rs` Aut code | refine+individualize | could call nauty, but self-contained on purpose |

## What does NOT convert (keep as-is)

- `modular.py` (`is_module`, `modular_decomposition`, `named_subgame`, `neq_tree`) — modular
  decomposition is a different algorithm, not an isomorphism/automorphism task.
- The brute `canonical_key` is worth keeping as the *tested reference oracle* even if
  the hot path moves to pynauty.

## Dependency note

`pynauty` (pip) is a C extension wrapping nauty; it would add a build dependency
(and CI would need it). The CLI tools (`dreadnaut`, `geng`, `gentourng`,
`directg`, `watercluster2`) need no Python dependency but cost a process spawn per
call — fine for batch/offline work (the `nauty/` scripts), not per-game in a loop.
Recommendation: keep the pure-Python `rpsfair` self-contained (no hard nauty dep)
for the tested n≤9 results; use nauty (CLI in `nauty/`, pynauty optionally) for the
high-n pushes, exactly as we already do.

## Docs

Rendered man pages for the tools we use are in this directory (`dreadnaut.txt`,
`nauty-geng.txt`, `nauty-gentourng.txt`, `nauty-directg.txt`, `nauty-listg.txt`,
`nauty-labelg.txt`, `nauty-amtog.txt`, `nauty-watercluster2.txt`,
`nauty-converseg.txt`, `nauty-countg.txt`). The Debian `dreadnaut` man page is only
a stub — the full dreadnaut command reference (the `d` digraph command, `x`
canonical labeling, automorphism output) and library docs are in the **nauty &
Traces User's Guide**: <https://pallini.di.uniroma1.it/> (McKay & Piperno).
