# The iso machinery runs on nauty

The package's isomorphism/automorphism machinery is **done** in nauty (via the
`pynauty` C extension) rather than hand-rolled. Each skew `{-1,0,+1}` game `M` *is*
an oriented graph (arc `i→j` iff `M[i,j]==1`; ties are non-edges), so nauty's
**directed-graph** mode applies directly. The representation (`structure.pynauty_graph`):

```python
import pynauty
def pynauty_graph(M):
    n = len(M)
    adj = {i: [j for j in range(n) if M[i, j] == 1] for i in range(n)}
    return pynauty.Graph(n, directed=True, adjacency_dict=adj)
# canonical key:        pynauty.certificate(g)
# automorphism group:   pynauty.autgrp(g)  -> (gens, grpsize1, grpsize2, orbits, numorbits)
# canonical labeling:   pynauty.canon_label(g)
```

## What was converted (now nauty)

| function (file) | now uses |
|---|---|
| `canonical_key` (structure.py) | `pynauty.certificate` |
| `canonical_key_fast`, `_canon_and_autos` (generate.py) | `pynauty.certificate` / `pynauty.canon_label` + `autgrp` |
| `automorphisms`, `aut_size` (metrics.py) | `pynauty.autgrp` (generators expanded; exact group order) |
| `num_orbits`, `node_orbits` (metrics.py) | `pynauty.autgrp` orbits |
| `generate_up_to_iso`, `generate_tournaments` (generate.py) | canonical augmentation keyed on nauty canon/autos |

Validated: `|Aut|` and canonical certificate match a brute n! oracle on
RPS/RPSLS/COP/Paley₇/ring₆ and all 582 n=5 classes; every enumeration count
(ISO, search, generation through n=6/n=8) is unchanged - and faster (n=6
generation ~17 s vs ~40–65 s before). Removed as dead/redundant: the refine-and-
individualize helpers (`_refine`, `_wl_colors`, `_normalize`, `_nauty_search`,
`canonical_key_nauty`) and the superseded `orbit_bytes`.

## What stays hand-rolled (on purpose)

- `modular.py` (modular decomposition) - a different algorithm, not iso/automorphism.
- `structure.orbit_hashes` / `matrix_hash` - the cheap per-candidate dedup for the
  brute `_enumerate` (one hash/candidate beats one certificate/candidate over the
  millions of n≤6 labelings).
- `rust/s2_count.rs`, `rust/regular_count.rs` - self-contained Rust (their refine+
  individualize Aut code is local by design, validated against `search_regular`).

## Dependency

`pynauty` is now a hard dependency (in `pyproject.toml` / `uv.lock`); it bundles and
builds nauty. The CLI tools (`gentourng`, `geng`, `directg`, `watercluster2`) are
still used directly by the `nauty/` and `rust/` high-n scripts.

## Docs

Rendered man pages for the tools we use are in this directory (`dreadnaut.txt`,
`nauty-geng.txt`, `nauty-gentourng.txt`, `nauty-directg.txt`, `nauty-listg.txt`,
`nauty-labelg.txt`, `nauty-amtog.txt`, `nauty-watercluster2.txt`,
`nauty-converseg.txt`, `nauty-countg.txt`). The Debian `dreadnaut` man page is only
a stub - the full dreadnaut command reference (the `d` digraph command, `x`
canonical labeling, automorphism output) and library docs are in the **nauty &
Traces User's Guide**: <https://pallini.di.uniroma1.it/> (McKay & Piperno).
