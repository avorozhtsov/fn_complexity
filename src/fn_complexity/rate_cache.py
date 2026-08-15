"""Persistent cache for numerically computed exchange rates."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Iterable

from .core import exchange_rate_result, normalize_signature

ALGORITHM = "exchange-rate-v1-grid2048-tol1e-13"
DEFAULT_CACHE_PATH = Path(__file__).with_name("exchange_rates_cache.json")


def _signature_text(signature: Iterable[int]) -> str:
    return ",".join(map(str, normalize_signature(signature)))


def _cache_key(implementer: Iterable[int], implemented: Iterable[int]) -> str:
    return f"g={_signature_text(implementer)}|f={_signature_text(implemented)}"


def _parse_label(value: str) -> tuple[int, ...]:
    if not (value.startswith("{") and value.endswith("}")):
        raise ValueError(f"not a signature label: {value!r}")
    return normalize_signature(
        int(part) for part in value[1:-1].split(",") if part
    )


class ExchangeRateCache:
    """Read-through JSON cache for ``C(implementer -> implemented)``.

    The cache is tied to the numerical algorithm and its accuracy settings.
    Infinite rates are encoded as the JSON string ``"inf"`` rather than a
    non-standard JSON number.
    """

    def __init__(self, path: Path | str = DEFAULT_CACHE_PATH) -> None:
        self.path = Path(path)
        self.hits = 0
        self.misses = 0
        self._rates: dict[str, float] = {}
        if self.path.exists():
            payload = json.loads(self.path.read_text(encoding="utf-8"))
            if payload.get("algorithm") != ALGORITHM:
                raise ValueError(
                    f"cache algorithm mismatch in {self.path}: "
                    f"{payload.get('algorithm')!r}"
                )
            for key, value in payload.get("rates", {}).items():
                self._rates[key] = (
                    math.inf if value == "inf" else float(value)
                )

    def __len__(self) -> int:
        return len(self._rates)

    def seed(
        self,
        implementer: Iterable[int],
        implemented: Iterable[int],
        rate: float,
    ) -> bool:
        """Add an already-computed value, without replacing cached data."""

        key = _cache_key(implementer, implemented)
        if key in self._rates:
            return False
        self._rates[key] = float(rate)
        return True

    def seed_from_matrix_csv(self, path: Path | str) -> int:
        """Import rates from a previously generated matrix CSV.

        Signature columns are recognized by their ``{...}`` headers, so this
        remains compatible if metadata columns are added to the matrix.
        """

        csv_path = Path(path)
        if not csv_path.exists():
            return 0
        added = 0
        with csv_path.open(encoding="utf-8", newline="") as stream:
            reader = csv.reader(stream)
            header = next(reader)
            signature_columns = [
                (index, _parse_label(value))
                for index, value in enumerate(header)
                if value.startswith("{") and value.endswith("}")
            ]
            for row in reader:
                implementer = _parse_label(row[0])
                for index, implemented in signature_columns:
                    value = row[index]
                    rate = math.inf if value in {"inf", "∞"} else float(value)
                    added += self.seed(implementer, implemented, rate)
        return added

    def get(
        self,
        implementer: Iterable[int],
        implemented: Iterable[int],
    ) -> float:
        """Return a cached rate, computing it only when absent."""

        implementer_tuple = normalize_signature(implementer)
        implemented_tuple = normalize_signature(implemented)
        key = _cache_key(implementer_tuple, implemented_tuple)
        if key in self._rates:
            self.hits += 1
            return self._rates[key]
        rate = exchange_rate_result(
            implementer_tuple,
            implemented_tuple,
            grid_size=2048,
            tolerance=1e-13,
        ).rate
        self._rates[key] = rate
        self.misses += 1
        return rate

    def save(self) -> None:
        """Atomically persist the cache in deterministic key order."""

        self.path.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "version": 1,
            "algorithm": ALGORITHM,
            "rates": {
                key: "inf" if math.isinf(value) else value
                for key, value in sorted(self._rates.items())
            },
        }
        temporary_path = self.path.with_suffix(self.path.suffix + ".tmp")
        temporary_path.write_text(
            json.dumps(payload, indent=2, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        temporary_path.replace(self.path)
