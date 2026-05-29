"""Unit tests for tie-twin reduction (reduce_twins) and the known reductions."""

import numpy as np
from toy import BRICK, COP, ELEM, RPS, RPSLS, ring

from rpsfair import canonical_key, maxmin_equilibrium, reduce_twins, twin_free


def test_reduce_twins_noop_on_twin_free():
    for M in [RPS, RPSLS, COP, ring(5)]:
        core, mult = reduce_twins(M)
        assert len(core) == len(M)
        assert mult == [1] * len(M)
        assert canonical_key(core) == canonical_key(M)


def test_brick_reduces_to_rps():
    core, mult = reduce_twins(BRICK)
    assert len(core) == 3
    assert sorted(mult) == [1, 1, 2]  # Brick+Boulder merge
    assert canonical_key(core) == canonical_key(RPS)
    assert not twin_free(BRICK)


def test_elemental_reduces_to_cop():
    core, mult = reduce_twins(ELEM)
    assert len(core) == 4
    assert sorted(mult) == [1, 1, 1, 2]  # Clay+Grass merge -> Witness
    assert canonical_key(core) == canonical_key(COP)
    assert not twin_free(ELEM)


def test_merged_equilibrium_mass_adds_up():
    # Witness (Clay+Grass merged) plays the sum of Clay's and Grass's mass.
    # ELEM is balanced -> uniform 1/5 each; merged node -> 2/5; cop maxmin has a 0.4 node.
    assert np.isclose(maxmin_equilibrium(ELEM).sum(), 1.0)
    cop_mm = maxmin_equilibrium(COP)
    assert np.isclose(max(cop_mm), 0.4)  # the witness (merged) node
