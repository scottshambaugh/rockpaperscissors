"""Coverage for modular decomposition and its equilibrium factorization."""

import numpy as np
from toy import BRICK, COP, ELEM, RPS, RPSLS

from rpsfair import modular_decomposition, named_subgame, neq_tree, num_equilibria, tree_summary


def test_named_subgames_and_all_tie_blocks():
    assert named_subgame(RPS) == "RPS"
    assert named_subgame(RPSLS) == "RPSLS"
    assert named_subgame(COP) == "cops"
    assert named_subgame(np.zeros((9, 9), dtype=np.int8)) == "tie9"


def test_twin_module_appears_in_decomposition_tree():
    tree = modular_decomposition(BRICK)
    assert tree["members"] == frozenset(range(4))
    assert tree_summary(tree) == "prime[tie[0, 1], 2, 3]"
    twin = tree["children"][0]
    assert twin["type"] == "tie"
    assert twin["members"] == frozenset({0, 1})


def test_elemental_game_has_expected_twin_module():
    tree = modular_decomposition(ELEM)
    modules = {child["members"] for child in tree["children"] if child["type"] == "tie"}
    assert modules == {frozenset({2, 4})}


def test_decomposition_equilibrium_factorization_matches_direct_count():
    for M in (RPS, COP, BRICK, ELEM, RPSLS):
        assert neq_tree(M) == num_equilibria(M)
