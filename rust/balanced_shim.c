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

/* Colored variant: `col[v]` is an arbitrary small integer per vertex; the
 * canonical form / labeling / orbits are computed under COLOR-PRESERVING
 * isomorphism (lab is initialized sorted by color, ptn marks class ends,
 * defaultptn = FALSE). Used by the weighted-balanced enumerations, where a
 * vertex's blow-up multiplicity is part of the structure. */
void rps_canon_colored(const uint64_t *arc, int n, const int *col,
                       uint64_t *canong, int *lab, int *orbits)
{
    int m = 1, i, j;
    graph g[MAXN];
    int ptn[MAXN];
    statsblk stats;
    DEFAULTOPTIONS_GRAPH(options);
    options.digraph = TRUE;
    options.getcanon = TRUE;
    options.defaultptn = FALSE;

    /* lab: vertices sorted by color (stable); ptn: 0 at each class end */
    int pos = 0;
    int cmin = col[0], cmax = col[0];
    for (i = 1; i < n; i++) {
        if (col[i] < cmin) cmin = col[i];
        if (col[i] > cmax) cmax = col[i];
    }
    int c;
    for (c = cmax; c >= cmin; c--) {   /* descending: heavier class first */
        int start = pos;
        for (i = 0; i < n; i++)
            if (col[i] == c) lab[pos++] = i;
        if (pos > start) {
            for (i = start; i < pos - 1; i++) ptn[i] = 1;
            ptn[pos - 1] = 0;
        }
    }
    for (i = 0; i < n; i++) g[i] = 0;
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


/* Automorphism group size of the digraph (for labeled-count ground truths):
 * returns stats.grpsize1 * 10^grpsize2 as a double (exact for |Aut| < 2^53). */
double rps_autsize(const uint64_t *arc, int n)
{
    int m = 1, i, j;
    graph g[MAXN];
    int lab[MAXN], ptn[MAXN], orbits[MAXN];
    statsblk stats;
    DEFAULTOPTIONS_GRAPH(options);
    options.digraph = TRUE;
    for (i = 0; i < n; i++) { g[i] = 0; lab[i] = i; ptn[i] = 1; }
    ptn[n - 1] = 0;
    for (i = 0; i < n; i++) {
        uint64_t row = arc[i];
        while (row) {
            j = __builtin_ctzll(row);
            row &= row - 1;
            ADDONEARC0(g, i, j, m);
        }
    }
    densenauty(g, lab, ptn, orbits, &options, &stats, m, n, NULL);
    double s = stats.grpsize1;
    for (i = 0; i < stats.grpsize2; i++) s *= 10.0;
    return s;
}

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
