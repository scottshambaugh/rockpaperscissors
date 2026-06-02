/* Thin C shim over nauty's densenauty so Rust doesn't have to replicate the
 * optionblk struct / WORDSIZE macros. One call returns the canonical form (for
 * dedup), the canonical labeling lab[] (lab[i] = vertex at canonical position i)
 * and the automorphism orbits[] (orbits[v] = orbit representative of v).
 *
 * Input `arc` is MY adjacency convention: arc[i] has bit (1<<j) set iff i->j
 * (i beats j). The shim re-emits arcs in nauty's BITT convention via ADDONEARC0,
 * so nauty computes correct digraph automorphisms. canong is returned verbatim
 * (nauty's bit convention); Rust only hashes it for dedup, so the convention is
 * irrelevant as long as it is consistent across isomorphic graphs.
 */
#define MAXN 64
#include "nauty.h"
#include <stdint.h>

void rps_canon(const uint64_t *arc, int n,
               uint64_t *canong, int *lab, int *orbits)
{
    int m = 1, i, j;
    graph g[MAXN];
    int ptn[MAXN];
    statsblk stats;
    DEFAULTOPTIONS_GRAPH(options);
    options.digraph = TRUE;   /* treat arcs as directed (asymmetric refinement) */
    options.getcanon = TRUE;

    for (i = 0; i < n; i++) {
        g[i] = 0;          /* EMPTYGRAPH for m==1 */
        lab[i] = i;
        ptn[i] = 1;
    }
    ptn[n - 1] = 0;
    for (i = 0; i < n; i++) {
        uint64_t row = arc[i];
        while (row) {
            j = __builtin_ctzll(row);
            row &= row - 1;
            ADDONEARC0(g, i, j, m);
        }
    }
    densenauty(g, lab, ptn, orbits, &options, &stats, m, n, (graph *)canong);
}
