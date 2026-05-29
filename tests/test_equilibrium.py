"""Unit tests for rpsfair.equilibrium (Nash equilibria of antisymmetric games)."""

import numpy as np
import pytest
from toy import BRICK, COP, DOMINATED, ELEM, PALEY7, RPS, RPSLS, ring

from rpsfair import (
    equilibrium_dim,
    equilibrium_info,
    equilibrium_vertices,
    has_fully_mixed,
    kernel_dim,
    max_entropy_equilibrium,
    maxmin_equilibrium,
    num_equilibria,
)

ALL = [RPS, RPSLS, ring(4), ring(5), ring(6), COP, BRICK, ELEM, PALEY7]


def is_equilibrium(M, p, tol=1e-6):
    """p is a symmetric-zero-sum equilibrium iff p is a distribution and M p <= 0."""
    p = np.asarray(p, float)
    return abs(p.sum() - 1) < tol and (p >= -tol).all() and (M.astype(float) @ p <= tol).all()


# --- kernel_dim and the parity theorem ---
@pytest.mark.parametrize(
    "M,expected",
    [
        (RPS, 1),
        (RPSLS, 1),
        (ring(5), 1),
        (ELEM, 3),
        (PALEY7, 1),
        (ring(4), 2),
        (ring(6), 2),
        (COP, 2),
    ],
)
def test_kernel_dim(M, expected):
    assert kernel_dim(M) == expected


@pytest.mark.parametrize("M", ALL)
def test_kernel_dim_parity_matches_n(M):
    # rank is always even, so dim ker(A) has the same parity as n
    assert kernel_dim(M) % 2 == len(M) % 2


# --- has_fully_mixed: existence test ---
@pytest.mark.parametrize("M", [RPS, RPSLS, ring(4), COP, BRICK, ELEM, PALEY7])
def test_fully_mixed_exists(M):
    ok, witness = has_fully_mixed(M)
    assert ok
    assert (np.asarray(witness) > 0).all()
    assert is_equilibrium(M, witness)


def test_dominated_has_no_fully_mixed():
    ok, witness = has_fully_mixed(DOMINATED)
    assert not ok and witness is None
    assert kernel_dim(DOMINATED) == 0


# --- maxmin (leximin) canonical equilibrium: exact values on toy models ---
def test_maxmin_uniform_when_symmetric():
    assert np.allclose(maxmin_equilibrium(RPS), [1 / 3] * 3)
    assert np.allclose(maxmin_equilibrium(RPSLS), [0.2] * 5)
    assert np.allclose(maxmin_equilibrium(ring(6)), [1 / 6] * 6)


def test_maxmin_cop_is_clean_2_2_2_4():
    assert np.allclose(maxmin_equilibrium(COP), [0.2, 0.2, 0.2, 0.4], atol=1e-6)


def test_maxmin_brick():
    assert np.allclose(maxmin_equilibrium(BRICK), [1 / 6, 1 / 6, 1 / 3, 1 / 3], atol=1e-6)


@pytest.mark.parametrize("M", ALL)
def test_maxmin_is_valid_equilibrium(M):
    assert is_equilibrium(M, maxmin_equilibrium(M))


def test_maxmin_balances_independent_modules_at_each_scale():
    # RPS[ tie2, cop, . ]: a twin pair and a cop module are independent free
    # directions of different scales. Leximin must even the twin pair (1/6,1/6),
    # not let one copy absorb the slack -- the multi-scale balancing fix.
    def beat(M, w, l):
        M[w, l], M[l, w] = 1, -1

    cop = np.zeros((4, 4), np.int8)
    for w, l in [(3, 0), (0, 2), (0, 1), (1, 2), (2, 3)]:
        beat(cop, w, l)
    H = np.array([[0, 1, -1], [-1, 0, 1], [1, -1, 0]], np.int8)  # RPS quotient
    blocks = [np.zeros((2, 2), np.int8), cop, np.zeros((1, 1), np.int8)]
    sizes = [2, 4, 1]
    off = np.cumsum([0, *sizes])
    G = np.zeros((7, 7), np.int8)
    for i, b in enumerate(blocks):
        G[off[i] : off[i + 1], off[i] : off[i + 1]] = b
        for j in range(3):
            if i != j:
                G[off[i] : off[i + 1], off[j] : off[j + 1]] = H[i, j]
    p = maxmin_equilibrium(G)
    assert is_equilibrium(G, p)
    assert np.allclose(p, [1 / 6, 1 / 6, 1 / 15, 1 / 15, 1 / 15, 2 / 15, 1 / 3], atol=1e-5)


# --- max-entropy alternative ---
def test_max_entropy_cop_is_structure_weighted():
    p = max_entropy_equilibrium(COP)
    # NOT the symmetric 0.2/0.2/0.2/0.4; structure-weighted ~ (.216,.177,.216,.392)
    assert np.allclose(p, [0.2156, 0.1766, 0.2156, 0.3922], atol=2e-3)
    assert is_equilibrium(COP, p)
    assert not np.allclose(p, [0.2, 0.2, 0.2, 0.4], atol=1e-2)


def test_max_entropy_uniform_when_unique():
    # when O is a single point, both canonical choices give it
    assert np.allclose(max_entropy_equilibrium(RPS), [1 / 3] * 3)
    assert np.allclose(max_entropy_equilibrium(RPSLS), [0.2] * 5)


# --- counting equilibria: vertices of the Nash polytope O ---
@pytest.mark.parametrize(
    "M,nverts",
    [
        (RPS, 1),  # unique
        (RPSLS, 1),  # unique
        (PALEY7, 1),  # unique
        (ring(5), 1),  # unique (odd n)
        (COP, 2),  # 1-D segment, 2 endpoints
        (BRICK, 2),  # twin pair -> segment
        (ring(6), 2),  # even n -> segment
        (ELEM, 4),  # 2-D family, 4 corners
    ],
)
def test_num_equilibria(M, nverts):
    assert num_equilibria(M) == nverts


@pytest.mark.parametrize("M", ALL)
def test_equilibrium_vertices_are_valid_and_distinct(M):
    V = equilibrium_vertices(M)
    assert len(V) >= 1
    for p in V:
        assert is_equilibrium(M, p)
    # vertices are pairwise distinct
    for i in range(len(V)):
        for j in range(i + 1, len(V)):
            assert not np.allclose(V[i], V[j], atol=1e-6)


@pytest.mark.parametrize("M", ALL)
def test_unique_iff_single_vertex(M):
    # exactly one extreme equilibrium <=> dim O == 0
    assert (num_equilibria(M) == 1) == (equilibrium_dim(M) == 0)


@pytest.mark.parametrize("M", ALL)
def test_equilibrium_dim_matches_kernel_for_this_family(M):
    # for these games O is exactly the fully-mixed kernel family: dim O = nullity - 1
    assert equilibrium_dim(M) == kernel_dim(M) - 1


@pytest.mark.parametrize("M", ALL)
def test_even_n_is_never_unique(M):
    # parity theorem: a unique equilibrium (single vertex) is possible only at odd n
    if len(M) % 2 == 0:
        assert num_equilibria(M) >= 2


# --- equilibrium_info summary ---
def test_equilibrium_info_unique_odd():
    info = equilibrium_info(RPS)
    assert info["nullity"] == 1
    assert info["fully_mixed"] is True
    assert info["family_dim"] == 0
    assert info["unique"] is True


def test_equilibrium_info_family_even():
    info = equilibrium_info(COP)
    assert info["nullity"] == 2
    assert info["family_dim"] == 1  # 1-D segment of equilibria
    assert info["unique"] is False


def test_equilibrium_info_high_dim_family_odd_n():
    # ELEM is n=5 but kernel dim 3 -> a 2-D family despite odd n
    info = equilibrium_info(ELEM)
    assert info["family_dim"] == 2
    assert info["unique"] is False
