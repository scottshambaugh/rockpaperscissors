"""Cache round-trip, hit/miss, and search-key isolation tests."""

import numpy as np

import rpsfair.cache as cache
import rpsfair.generate as generate_module
import rpsfair.search as search_module
from rpsfair.generate import search_balanced_gen, search_inclusive_gen
from rpsfair.search import search_balanced, search_balanced_stream, search_inclusive


def test_cache_round_trip(tmp_path, monkeypatch):
    monkeypatch.setattr(cache, "CACHE_DIR", str(tmp_path))
    M = np.array([[0, 1], [-1, 0]], dtype=np.int8)
    xs = np.array([1 / 3, 2 / 3])
    cache.save("sample", [(M, xs)])
    loaded = cache.load("sample")
    assert len(loaded) == 1
    assert np.array_equal(loaded[0][0], M)
    assert np.allclose(loaded[0][1], xs, atol=1e-6)


def test_cached_computes_once_then_hits(tmp_path, monkeypatch):
    monkeypatch.setattr(cache, "CACHE_DIR", str(tmp_path))
    calls = 0

    def compute():
        nonlocal calls
        calls += 1
        return [(np.zeros((1, 1), dtype=np.int8), np.ones(1))]

    assert len(cache.cached("once", compute)) == 1
    assert len(cache.cached("once", compute)) == 1
    assert calls == 1


def test_search_implementations_use_isolated_cache_keys(monkeypatch):
    names = []

    def record(name, _fn):
        names.append(name)
        return []

    monkeypatch.setattr(search_module, "cached", record)
    monkeypatch.setattr(generate_module, "cached", record, raising=False)
    # The generator functions import cached locally, so intercept the source
    # module as well; the search functions use their module-level binding.
    monkeypatch.setattr(cache, "cached", record)
    search_balanced(3)
    search_balanced_stream(3)
    search_balanced_gen(3)
    search_inclusive(3)
    search_inclusive_gen(3)
    assert names == [
        "balanced_brute_n3",
        "balanced_stream_n3",
        "balanced_gen_n3",
        "inclusive_brute_n3",
        "inclusive_gen_n3",
    ]
