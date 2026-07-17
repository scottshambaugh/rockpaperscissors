"""Validation of the enumeration against known counts, plus structural-invariant
property tests on every game a search returns.

The fast row-by-row searches (search_regular, search_balanced_fast) are used for
counts so the suite is cache-independent; the brute-force search_balanced is only
exercised on small n where it terminates quickly. Big enumerations are `slow`.
"""

import numpy as np
import pytest

from rpsfair import (
    canonical_key,
    connected,
    has_fully_mixed,
    is_completely_mixed,
    k_paradoxical,
    num_equilibria,
    paradoxical,
    profile_of,
    regular,
    search_balanced,
    search_balanced_fast,
    search_completely_mixed,
    search_inclusive,
    search_regular,
    twin_free,
)

# Known counts: total and twin-free, indexed by n. (Established by this project.)
REGULAR = {3: (1, 1), 4: (1, 1), 5: (2, 2), 6: (5, 4), 7: (13, 12), 8: (82, 76)}
BALANCED = {3: (1, 1), 4: (1, 1), 5: (4, 3), 6: (16, 13), 7: (175, 152)}
INCLUSIVE = {3: (1, 1), 4: (3, 2), 5: (15, 8), 6: (222, 177)}
# completely mixed (Kaplansky): every equilibrium fully mixed = unique + fully
# mixed. Zero at every even n (parity theorem); twin-free always equals total.
COMPLETELY_MIXED = {3: 1, 4: 0, 5: 7, 6: 0}


def counts(structs):
    return len(structs), sum(twin_free(M) for M, _ in structs)


# ---------------- known-count validation ----------------
@pytest.mark.parametrize("n", [3, 4, 5, 6, 7])
def test_regular_counts(n):
    assert counts(search_regular(n)) == REGULAR[n]


@pytest.mark.parametrize("n", [3, 4, 5, 6])
def test_balanced_counts(n):
    assert counts(search_balanced_fast(n)) == BALANCED[n]


@pytest.mark.parametrize("n", [3, 4, 5])
def test_balanced_brute_matches_fast(n):
    # the brute-force search_balanced agrees with the pruned row-by-row search
    assert len(search_balanced(n)) == len(search_balanced_fast(n)) == BALANCED[n][0]


@pytest.mark.parametrize("n", [3, 4, 5])
def test_inclusive_counts(n):
    assert counts(search_inclusive(n)) == INCLUSIVE[n]


@pytest.mark.parametrize("n", [3, 4, 5])
def test_completely_mixed_counts(n):
    assert len(search_completely_mixed(n)) == COMPLETELY_MIXED[n]


def test_completely_mixed_even_n_is_empty_without_search():
    # Kaplansky parity: an even-order skew-symmetric game is never completely
    # mixed, so even n short-circuits (n=100 would be unenumerable otherwise)
    assert search_completely_mixed(100) == []


@pytest.mark.parametrize("n", [3, 4, 5])
def test_completely_mixed_is_unique_equilibrium(n):
    # the cheap nullity-1 + one-signed-kernel test agrees with the independent
    # vertex enumeration: completely mixed <=> exactly one extreme equilibrium
    for M, _ in search_inclusive(n):
        assert is_completely_mixed(M) == (num_equilibria(M) == 1)


@pytest.mark.parametrize("n", [3, 5])
def test_completely_mixed_is_twin_free(n):
    # a tie-twin pair puts e_i - e_j in ker(M), so twins break nullity 1
    for M, _ in search_completely_mixed(n):
        assert twin_free(M)


def test_two_paradox_first_appears_at_n7():
    # the smallest two-paradox (P2) game is the n=7 Paley tournament: exactly one
    for n in [3, 4, 5, 6]:
        assert sum(k_paradoxical(M, 2) for M, _ in search_regular(n)) == 0
    assert sum(k_paradoxical(M, 2) for M, _ in search_regular(7)) == 1


@pytest.mark.slow
@pytest.mark.parametrize("n", [8])
def test_regular_counts_big(n):
    assert counts(search_regular(n)) == REGULAR[n]


@pytest.mark.slow
def test_balanced_count_n7():
    assert counts(search_balanced_fast(7)) == BALANCED[7]


@pytest.mark.slow
def test_inclusive_count_n6():
    assert counts(search_inclusive(6)) == INCLUSIVE[6]


@pytest.mark.slow
def test_completely_mixed_count_n6():
    assert len(search_completely_mixed(6)) == COMPLETELY_MIXED[6]


# ---------------- structural invariants of returned games ----------------
@pytest.mark.parametrize("n", [3, 4, 5, 6])
def test_regular_games_are_regular_connected_paradoxical(n):
    for M, _ in search_regular(n):
        assert regular(M) and connected(M) and paradoxical(M)


@pytest.mark.parametrize("n", [3, 4, 5, 6])
def test_balanced_games_have_zero_row_sums(n):
    for M, _ in search_balanced_fast(n):
        assert (M.sum(axis=1) == 0).all()
        assert connected(M) and paradoxical(M)


@pytest.mark.parametrize("n", [3, 4, 5])
def test_inclusive_games_are_fully_mixed(n):
    for M, _ in search_inclusive(n):
        assert has_fully_mixed(M)[0]
        assert connected(M) and paradoxical(M)


@pytest.mark.parametrize(
    "search,n",
    [(search_regular, 6), (search_balanced_fast, 6), (search_inclusive, 5)],
)
def test_results_are_pairwise_non_isomorphic(search, n):
    keys = [canonical_key(M) for M, _ in search(n)]
    assert len(keys) == len(set(keys))


def test_skew_symmetric_and_zero_diagonal():
    for M, _ in search_inclusive(5):
        assert np.array_equal(M, -M.T)
        assert (np.diag(M) == 0).all()


@pytest.mark.parametrize(
    "search,n",
    [(search_regular, 6), (search_balanced_fast, 6), (search_inclusive, 5)],
)
def test_every_equilibrium_vertex_is_in_the_kernel(search, n):
    # For fair (paradoxical, connected) games O = ker(M) ∩ Δ: every extreme
    # equilibrium satisfies M v = 0 exactly, not merely M v <= 0. (Not true for
    # skew-symmetric games at large -- a dominated move breaks it -- so this is a
    # property of the category, and the null-space recipe is sufficient here.)
    from rpsfair import equilibrium_vertices

    for M, _ in search(n):
        A = M.astype(float)
        for v in equilibrium_vertices(M):
            assert np.allclose(A @ v, 0.0, atol=1e-9)


def test_nested_subsets_regular_in_balanced_in_inclusive():
    # every regular game is balanced; every balanced game is inclusive (up to iso)
    for n in [4, 5]:
        reg = {canonical_key(M) for M, _ in search_regular(n)}
        bal = {canonical_key(M) for M, _ in search_balanced_fast(n)}
        inc = {canonical_key(M) for M, _ in search_inclusive(n)}
        assert reg <= bal <= inc
