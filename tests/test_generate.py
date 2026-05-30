"""Tests for rpsfair.generate: isomorph-free generation, the nauty-style canonical
form / automorphisms, the batched filters, and the constraint-pruned searches.

The *_gen / *_stream searches cache to disk, so the `no_cache` fixture forces a
fresh computation -- otherwise we'd only be checking the cache, not the logic.
"""

import numpy as np
import pytest
from toy import COP, PALEY7, RPS, RPSLS, ring

from rpsfair import automorphisms, canonical_key, connected, paradoxical
from rpsfair.equilibrium import has_fully_mixed
from rpsfair.generate import (
    _canon_and_autos,
    _connected_batch,
    _has_fully_mixed_batch,
    _paradoxical_batch,
    canonical_key_fast,
    generate_up_to_iso,
    search_balanced_gen,
    search_inclusive_gen,
    search_regular_gen,
)
from rpsfair.search import search_balanced_stream, search_regular_stream

# number of isomorphism classes of n-node skew {-1,0,+1} games
ISO = {1: 1, 2: 2, 3: 7, 4: 42, 5: 582}


@pytest.fixture
def no_cache(monkeypatch):
    import rpsfair.cache as C

    monkeypatch.setattr(C, "load", lambda name: None)
    monkeypatch.setattr(C, "save", lambda name, items: None)


# ---------------- isomorph-free generation ----------------
@pytest.mark.parametrize("n", [2, 3, 4, 5])
def test_generate_count(n):
    g = list(generate_up_to_iso(n))
    assert len(g) == ISO[n]
    # one representative per class: no two generated games are isomorphic
    assert len({canonical_key(M) for M in g}) == ISO[n]


def test_generate_skew_and_zero_diagonal():
    for M in generate_up_to_iso(5):
        assert np.array_equal(M, -M.T)
        assert (np.diag(M) == 0).all()


def test_generate_matches_brute_iso_count():
    # brute-force every labeling of n=4, dedup by the true canonical key
    from itertools import product

    iu = np.triu_indices(4, 1)
    keys = set()
    for vals in product((-1, 0, 1), repeat=len(iu[0])):
        M = np.zeros((4, 4), np.int8)
        M[iu] = vals
        M.T[iu] = [-x for x in vals]
        keys.add(canonical_key(M))
    assert len(keys) == ISO[4] == len(list(generate_up_to_iso(4)))


# ---------------- nauty canonical form + automorphisms ----------------
@pytest.mark.parametrize("M", [RPS, RPSLS, COP, PALEY7, ring(6)])
def test_canonical_key_fast_relabel_invariant(M):
    rng = np.random.default_rng(1)
    base = canonical_key_fast(M)
    for _ in range(8):
        p = rng.permutation(len(M))
        assert canonical_key_fast(M[np.ix_(p, p)]) == base


def test_canonical_key_fast_distinguishes_classes():
    g = list(generate_up_to_iso(5))
    assert len({canonical_key_fast(M) for M in g}) == len(g)


@pytest.mark.parametrize(
    "M", [RPS, RPSLS, COP, PALEY7, ring(6), np.zeros((4, 4), np.int8), np.zeros((1, 1), np.int8)]
)
def test_canon_autos_match_brute(M):
    _, autos = _canon_and_autos(M)
    nauty = {tuple(int(x) for x in a) for a in autos}
    assert nauty == {tuple(p) for p in automorphisms(M)}


def test_canon_autos_match_brute_all_n5():
    for M in generate_up_to_iso(5):
        _, autos = _canon_and_autos(M)
        assert {tuple(int(x) for x in a) for a in autos} == {tuple(p) for p in automorphisms(M)}


def test_canon_perm_yields_canonical_matrix():
    # M relabeled by the canonical perm equals the canonical-key matrix
    for M in list(generate_up_to_iso(4))[:20]:
        cperm, _ = _canon_and_autos(M)
        relabeled = tuple(int(x) for x in M[np.ix_(cperm, cperm)].reshape(-1))
        assert relabeled == canonical_key_fast(M)


# ---------------- batched filters match the scalar versions ----------------
def test_batch_filters_match_scalar():
    g = list(generate_up_to_iso(5))
    Zs = np.array(g)
    assert (_paradoxical_batch(Zs) == np.array([paradoxical(M) for M in g])).all()
    assert (_connected_batch(Zs) == np.array([connected(M) for M in g])).all()
    assert (_has_fully_mixed_batch(Zs) == np.array([has_fully_mixed(M)[0] for M in g])).all()


# ---------------- constraint-pruned searches reproduce known counts ----------------
@pytest.mark.parametrize("n,want", [(3, 1), (4, 3), (5, 15)])
def test_search_inclusive_gen(n, want, no_cache):
    assert len(search_inclusive_gen(n)) == want


@pytest.mark.parametrize("n,want", [(5, 4), (6, 16)])
def test_search_balanced_gen(n, want, no_cache):
    g = search_balanced_gen(n)
    assert len(g) == want
    assert all((M.sum(1) == 0).all() and paradoxical(M) and connected(M) for M, _ in g)


@pytest.mark.parametrize("n,want", [(5, 2), (6, 5)])
def test_search_regular_gen(n, want, no_cache):
    from rpsfair.structure import profile_of

    g = search_regular_gen(n)
    assert len(g) == want
    assert all(len(set(profile_of(M))) == 1 and connected(M) for M, _ in g)  # all profiles equal


@pytest.mark.parametrize("n,want", [(5, 4), (6, 16)])
def test_search_balanced_stream(n, want, no_cache):
    assert len(search_balanced_stream(n)) == want


@pytest.mark.parametrize("n,want", [(6, 5), (7, 13)])
def test_search_regular_stream(n, want, no_cache):
    assert len(search_regular_stream(n)) == want


# ---------------- slow: larger validations ----------------
@pytest.mark.slow
def test_inclusive_gen_matches_known_n6(no_cache):
    # the n=6 inclusive count via generation equals the established 222
    assert len(search_inclusive_gen(6)) == 222


@pytest.mark.slow
def test_generate_count_n6():
    assert sum(1 for _ in generate_up_to_iso(6)) == 21480


@pytest.mark.slow
def test_search_regular_gen_n8(no_cache):
    assert len(search_regular_gen(8)) == 82
