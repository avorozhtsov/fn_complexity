#!/usr/bin/env python3
"""Extend the cached list of zeta zero ordinates to 2400 and cache as .npy."""

from __future__ import annotations

from pathlib import Path

import mpmath
import numpy as np

DIRECTORY = Path(__file__).resolve().parent
SOURCE = DIRECTORY / "zeta_zeros_1200.npy"
TARGET = DIRECTORY / "zeta_zeros_2400.npy"
TOTAL = 2400


def main() -> int:
    known = list(np.load(SOURCE))
    for index in range(len(known) + 1, TOTAL + 1):
        known.append(float(mpmath.zetazero(index).imag))
        if index % 100 == 0:
            print(f"{index} zeros, gamma = {known[-1]:.3f}", flush=True)
    np.save(TARGET, np.array(known))
    print(f"wrote {TARGET.name}: {len(known)} zeros, T = {known[-1]:.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
