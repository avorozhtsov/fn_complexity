import json
import math

from fn_complexity import ExchangeRateCache, RATE_CACHE_ALGORITHM


def test_seeded_rate_is_reused_without_computation(tmp_path):
    cache = ExchangeRateCache(tmp_path / "rates.json")
    assert cache.seed((2, 2), (3, 1), 0.630929753571457)

    assert cache.get((2, 2), (3, 1)) == 0.630929753571457
    assert cache.hits == 1
    assert cache.misses == 0


def test_cache_round_trip_encodes_infinity_as_strict_json(tmp_path):
    path = tmp_path / "rates.json"
    cache = ExchangeRateCache(path)
    cache.seed((2,), (1,), math.inf)
    cache.save()

    payload = json.loads(path.read_text(encoding="utf-8"))
    assert payload["algorithm"] == RATE_CACHE_ALGORITHM
    assert list(payload["rates"].values()) == ["inf"]

    loaded = ExchangeRateCache(path)
    assert math.isinf(loaded.get((2,), (1,)))
    assert loaded.hits == 1
    assert loaded.misses == 0
