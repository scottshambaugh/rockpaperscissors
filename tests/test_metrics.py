"""Unit tests for rpsfair.metrics (orbits, |Aut|, cuts, gini, tie fraction)."""

import numpy as np
import pytest
from toy import BRICK, COP, ELEM, PALEY7, RPS, RPSLS, game, ring

from rpsfair import (
    aut_size,
    automorphisms,
    gini,
    node_orbits,
    num_cuts,
    num_orbits,
    tie_fraction,
)


# --- automorphisms / |Aut| ---
@pytest.mark.parametrize(
    "M,size",
    [
        (RPS, 3),
        (RPSLS, 5),
        (ring(4), 4),
        (ring(6), 6),
        (PALEY7, 21),
        (COP, 1),
        (BRICK, 2),
        (ELEM, 2),
    ],
)
def test_aut_size(M, size):
    assert aut_size(M) == size


@pytest.mark.parametrize("M", [RPS, COP, BRICK, ELEM])
def test_identity_is_an_automorphism(M):
    autos = automorphisms(M)
    assert tuple(range(len(M))) in autos
    assert len(autos) == aut_size(M)


# --- orbits ---
@pytest.mark.parametrize(
    "M,orbits",
    [
        (RPS, 1),
        (RPSLS, 1),
        (ring(5), 1),
        (PALEY7, 1),  # circulant -> 1 orbit
        (COP, 4),
        (BRICK, 3),
        (ELEM, 4),
    ],
)
def test_num_orbits(M, orbits):
    assert num_orbits(M) == orbits


def test_orbit_partition_brick_and_elem():
    # Brick & Boulder share an orbit; Clay & Grass share an orbit
    assert sorted(sorted(o) for o in node_orbits(BRICK)) == [[0, 1], [2], [3]]
    assert sorted(sorted(o) for o in node_orbits(ELEM)) == [[0], [1], [2, 4], [3]]


def test_orbits_equal_n_iff_aut_trivial():
    assert num_orbits(COP) == len(COP) and aut_size(COP) == 1
    assert num_orbits(RPS) == 1 and aut_size(RPS) > 1


# --- num_cuts (articulation points of the decisive subgraph) ---
def test_num_cuts_cycle_is_zero():
    assert num_cuts(RPS) == 0
    assert num_cuts(ring(6)) == 0


def test_num_cuts_path_has_articulation():
    # decisive edges 0->1->2 form a path; node 1 is a cut vertex
    M = game(3, [(0, 1), (1, 2)])
    assert num_cuts(M) == 1


# --- tie_fraction ---
def test_tie_fraction_tournament_is_zero():
    assert tie_fraction(RPS) == 0.0
    assert tie_fraction(RPSLS) == 0.0
    assert tie_fraction(PALEY7) == 0.0


@pytest.mark.parametrize("n", [4, 5, 6, 7])
def test_tie_fraction_ring_formula(n):
    # ring: each node ties n-3 others -> tie fraction (n-3)/(n-1)
    assert tie_fraction(ring(n)) == pytest.approx((n - 3) / (n - 1))


def test_tie_fraction_cop():
    assert tie_fraction(COP) == pytest.approx(1 / 6)  # one tie pair of six


# --- gini ---
def test_gini_uniform_is_zero():
    assert gini(np.ones(5) / 5) == pytest.approx(0.0)
    assert gini(np.array([0.25, 0.25, 0.25, 0.25])) == pytest.approx(0.0)


def test_gini_cop_equilibrium():
    assert gini(np.array([0.2, 0.2, 0.2, 0.4])) == pytest.approx(0.15)


def test_gini_increases_with_concentration():
    assert gini(np.array([0.1, 0.2, 0.3, 0.4])) > gini(np.array([0.22, 0.24, 0.26, 0.28]))
