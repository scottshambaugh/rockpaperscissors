"""Unit tests for rpsfair.equilibrium (Nash equilibria of antisymmetric games)."""

import numpy as np
import pytest
from toy import BRICK, COP, DOMINATED, ELEM, PALEY7, RPS, RPSLS, ring

from rpsfair import (
    equilibrium_info,
    has_fully_mixed,
    kernel_dim,
    max_entropy_equilibrium,
    maxmin_equilibrium,
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
