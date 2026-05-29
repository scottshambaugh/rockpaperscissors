"""Unit tests for rpsfair.structure predicates and isomorphism helpers."""

import numpy as np
import pytest
from toy import BRICK, COP, DOMINATED, ELEM, PALEY7, RPS, RPSLS, game, ring

from rpsfair import (
    canonical_key,
    canonicalize,
    connected,
    matrix_hash,
    orbit_bytes,
    orbit_hashes,
    paradoxical,
    profile_of,
    regular,
    twin_free,
)


@pytest.mark.parametrize(
    "M,expected",
    [
        (RPS, True),
        (RPSLS, True),
        (COP, True),
        (BRICK, True),
        (ELEM, True),
        (PALEY7, True),
        (DOMINATED, False),
    ],  # dominated 4th move has no win
)
def test_paradoxical(M, expected):
    assert paradoxical(M) is expected


def test_paradoxical_needs_win_and_loss():
    assert paradoxical(RPS)
    # a move with all losses, or all wins, breaks it
    M = game(3, [(0, 1), (0, 2)])  # node 0 never loses
    assert not paradoxical(M)


@pytest.mark.parametrize("M", [RPS, RPSLS, COP, BRICK, ELEM, PALEY7, ring(6)])
def test_connected(M):
    assert connected(M)


def test_disconnected():
    # two independent RPS triangles, no edge between -> not connected
    M = np.zeros((6, 6), np.int8)
    for a, b in [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3)]:
        M[a, b] = 1
        M[b, a] = -1
    assert not connected(M)


@pytest.mark.parametrize(
    "M,expected",
    [
        (RPS, True),
        (RPSLS, True),
        (ring(5), True),
        (PALEY7, True),
        (COP, False),
        (BRICK, False),
        (ELEM, False),
    ],
)
def test_regular(M, expected):
    assert regular(M) is expected


@pytest.mark.parametrize(
    "M,expected",
    [
        (RPS, True),
        (RPSLS, True),
        (COP, True),
        (PALEY7, True),
        (BRICK, False),
        (ELEM, False),
    ],  # Brick/Boulder and Clay/Grass are twins
)
def test_twin_free(M, expected):
    assert twin_free(M) is expected


def test_profile_of():
    assert profile_of(RPS) == ((1, 0, 1), (1, 0, 1), (1, 0, 1))
    assert profile_of(COP) == ((1, 0, 2), (1, 1, 1), (1, 1, 1), (2, 0, 1))
    assert profile_of(ring(5)) == tuple([(1, 2, 1)] * 5)


# --- isomorphism: canonical_key, orbit hashing, canonicalize ---
def _relabel(M, perm):
    perm = np.asarray(perm)
    return M[np.ix_(perm, perm)]


@pytest.mark.parametrize("M", [RPS, COP, BRICK, ELEM, RPSLS])
def test_canonical_key_relabel_invariant(M):
    rng = np.random.default_rng(0)
    key = canonical_key(M)
    for _ in range(5):
        perm = rng.permutation(len(M))
        assert canonical_key(_relabel(M, perm)) == key


def test_canonical_key_distinguishes_nonisomorphic():
    assert canonical_key(COP) != canonical_key(BRICK)  # same profiles, different games
    assert canonical_key(RPS) != canonical_key(ring(4))


def test_orbit_bytes_contains_self_and_relabels():
    orbit = orbit_bytes(COP)
    assert COP.tobytes() in orbit
    perm = [3, 1, 0, 2]
    assert _relabel(COP, perm).tobytes() in orbit


def test_matrix_hash_orbit_hashes_consistency():
    # every relabeling's hash is in the orbit-hash set; the set is the matrix_hash of orbit members
    assert matrix_hash(COP) in orbit_hashes(COP)
    assert matrix_hash(_relabel(COP, [3, 2, 1, 0])) in orbit_hashes(COP)
    # distinct games -> disjoint (with overwhelming probability) hash sets
    assert matrix_hash(RPS) not in orbit_hashes(ring(4))


def test_canonicalize_dedups_relabelings():
    rng = np.random.default_rng(1)
    copies = [_relabel(COP, rng.permutation(4)) for _ in range(8)] + [RPS.copy()]
    # COP copies collapse to 1, RPS is a different size/game -> 2 classes
    out = canonicalize([c for c in copies if len(c) == 4])
    assert len(out) == 1
