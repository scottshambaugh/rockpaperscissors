"""Unit tests for rpsfair.display (text rendering of games)."""

import numpy as np
import pytest
from toy import COP, RPS

from rpsfair import letter_labels, pretty, wtl_labels


def test_letter_labels():
    assert letter_labels(3) == ["A", "B", "C"]
    assert letter_labels(5)[-1] == "E"


def test_wtl_labels():
    assert wtl_labels(RPS) == ["1·0·1"] * 3  # 1·0·1
    # per-node order Cop, K-9, Perp, Witness (not the sorted profile multiset)
    assert wtl_labels(COP) == ["2·0·1", "1·1·1", "1·0·2", "1·1·1"]


def test_pretty_symbols_and_shape():
    out = pretty(RPS, labels=["Rock", "Paper", "Scissors"])
    assert "Rock" in out and "Paper" in out and "Scissors" in out
    # upper-triangular: contains +, -, and no stray digits from the matrix
    assert "+" in out and "-" in out
    # one header line + separator + n rows
    assert len(out.splitlines()) == 2 + len(RPS)


def test_pretty_default_uses_wtl():
    out = pretty(RPS)
    assert "1·0·1" in out


def test_pretty_label_count_mismatch_raises():
    with pytest.raises(ValueError):
        pretty(RPS, labels=["only", "two"])


def test_pretty_cell_values():
    # 2x2-ish check: RPS row 0 beats 1 (+), loses to 2 (-)
    M = np.array([[0, 1, -1], [-1, 0, 1], [1, -1, 0]], dtype=np.int8)
    out = pretty(M, labels=["a", "b", "c"])
    # first data row 'a' should show '+' (vs b) then '-' (vs c).
    # data rows contain '|'; that excludes the header and the '+---' separator.
    a_row = next(ln for ln in out.splitlines() if "|" in ln and ln.lstrip().startswith("a"))
    assert a_row.index("+") < a_row.index("-")
