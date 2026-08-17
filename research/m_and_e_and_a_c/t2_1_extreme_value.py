#!/usr/bin/env python3
"""T2.1, part 2 -- the deficit ``2g - max_c(-a_c)/sqrt(q)`` as an extreme value.

If the normalised traces ``T_c = -a_c / sqrt(q)`` equidistribute for the Haar
measure of ``USp(2g)`` (large monodromy, Katz--Sarnak) then

    m / sqrt(q) = max_c T_c

is the maximum of about ``q`` samples from the ``USp(2g)`` trace law, and the
deficit ``2g - m/sqrt(q)`` is a pure extreme-value quantity.

Everything below is computed from the *exact* lower-tail law.  Writing
``x_i = cos theta_i``, the Weyl integration formula for ``USp(2g)`` reads

    dmu = (2^{g^2} / (g! pi^g)) prod_{i<j}(x_i - x_j)^2 prod_i sqrt(1 - x_i^2) dx

on ``[-1,1]^g``, and ``T = 2 sum_i x_i``.  Substituting ``1 - x_i = (eps/2) v_i``
turns ``{2g - T < eps}`` into the standard simplex ``{v > 0, sum v_i < 1}``:

    P(2g - T < eps) = (2^{g^2} / (g! pi^g)) (eps/2)^{d/2} J_g(eps),   d = 2g^2 + g
    J_g(eps) = int_simplex prod_{i<j}(v_j - v_i)^2 prod_i sqrt(v_i)
                           prod_i sqrt(2 - (eps/2) v_i) [ (eps/2) v_i <= 2 ] dv

This is exact for every eps, not just asymptotically; the whole ``eps``
dependence outside the prefactor sits in a bounded analytic factor.  Letting
``eps -> 0`` gives the edge law

    P(2g - T < eps) ~ K_g eps^{d/2},   d = dim USp(2g) = 2g^2 + g,

so with ``n = q`` independent samples the deficit is Weibull with shape ``d/2``:

    E[2g - max T] = Gamma(1 + 2/d) (q K_g)^{-2/d}.

``J_g`` is evaluated by Monte Carlo on the simplex (bounded integrand, so the
estimator is well behaved); the module cross-checks it against a closed form at
``g = 1`` and against rejection sampling of the Weyl measure at ``g = 2``.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np

OUTPUT_DIRECTORY = Path(__file__).resolve().parent
RNG = np.random.default_rng(20260816)
SAMPLES = 4_000_000


def weyl_constant(genus: int) -> float:
    return 2.0 ** (genus**2) / (math.factorial(genus) * math.pi**genus)


def simplex_sample(genus: int, count: int) -> np.ndarray:
    """Uniform on {v > 0, sum v_i < 1} in R^genus (volume 1/genus!)."""

    exponentials = RNG.exponential(size=(count, genus + 1))
    return exponentials[:, :genus] / exponentials.sum(axis=1, keepdims=True)


class TailLaw:
    """Exact lower-tail law of the USp(2g) trace, P(2g - T < eps).

    ``J_g(eps)`` is estimated once on a grid of ``eps`` with common random
    numbers and then interpolated; it is smooth, bounded and slowly varying, so
    the interpolation error is far below the Monte Carlo error.
    """

    def __init__(self, genus: int, count: int = SAMPLES, knots: int = 129) -> None:
        self.genus = genus
        self.dimension = 2 * genus**2 + genus
        v = simplex_sample(genus, count)
        vandermonde = np.ones(count)
        for i in range(genus):
            for j in range(i + 1, genus):
                vandermonde *= (v[:, j] - v[:, i]) ** 2
        base = vandermonde * np.sqrt(v).prod(axis=1)
        self.prefactor = weyl_constant(genus) / math.factorial(genus)

        self.knots = np.linspace(0.0, 4.0 * genus, knots)
        values, errors = [], []
        for epsilon in self.knots:
            residual = 2.0 - (epsilon / 2.0) * v
            weight = np.where(
                (residual > 0.0).all(axis=1),
                np.sqrt(np.clip(residual, 0.0, None)).prod(axis=1),
                0.0,
            )
            sample = base * weight
            values.append(float(sample.mean()))
            errors.append(float(sample.std() / math.sqrt(count)))
        self.integral = np.array(values)
        self.integral_error = np.array(errors)

    def _interpolate(self, epsilon: float) -> float:
        return float(np.interp(epsilon, self.knots, self.integral))

    def cdf(self, epsilon: float) -> float:
        """P(2g - T < epsilon)."""

        if epsilon <= 0.0:
            return 0.0
        if epsilon >= 4.0 * self.genus:
            return 1.0
        value = (
            self.prefactor
            * (epsilon / 2.0) ** (self.dimension / 2.0)
            * self._interpolate(epsilon)
        )
        return min(1.0, value)

    def cdf_relative_error(self, epsilon: float) -> float:
        error = float(np.interp(epsilon, self.knots, self.integral_error))
        return error / max(self._interpolate(epsilon), 1e-300)

    @property
    def edge_constant(self) -> float:
        """K_g with P(2g - T < eps) ~ K_g eps^{d/2} as eps -> 0."""

        return self.prefactor * 2.0 ** (-self.dimension / 2.0) * float(self.integral[0])

    def expected_deficit(self, samples: int, grid: int = 40000) -> float:
        """E[min over `samples` draws of (2g - T)], assuming independence."""

        # Concentrate the quadrature around where the survival function turns.
        upper = min(4.0 * self.genus, 8.0 * self.median_deficit(samples))
        points = np.linspace(0.0, upper, grid + 1)
        cdf = np.clip(
            self.prefactor
            * (points / 2.0) ** (self.dimension / 2.0)
            * np.interp(points, self.knots, self.integral),
            0.0,
            1.0,
        )
        survival = (1.0 - cdf) ** samples
        return float(np.trapezoid(survival, points))

    def median_deficit(self, samples: int) -> float:
        """Median of the same minimum, from (1 - F(t))^n = 1/2."""

        target = 1.0 - 0.5 ** (1.0 / samples)
        low, high = 0.0, 4.0 * self.genus
        for _ in range(80):
            middle = 0.5 * (low + high)
            if self.cdf(middle) < target:
                low = middle
            else:
                high = middle
        return 0.5 * (low + high)


# --------------------------------------------------------------------------
# cross-checks on the tail law
# --------------------------------------------------------------------------
def check_genus_one(law: TailLaw) -> None:
    """g = 1: P(2 - 2 cos t < eps) = (2/pi)(theta0 - sin(2 theta0)/2), closed form."""

    print("  g=1 tail law vs closed form")
    for epsilon in (0.05, 0.5, 1.0, 2.0, 3.5):
        theta = math.acos(max(-1.0, 1.0 - epsilon / 2.0))
        exact = (theta - math.sin(2 * theta) / 2.0) / math.pi
        print(f"    eps={epsilon:>5.2f}  formula={law.cdf(epsilon):.8f}  "
              f"exact={exact:.8f}  rel err={abs(law.cdf(epsilon) / exact - 1):.2e}")
    print(f"    K_1 = {law.edge_constant:.8f}   2/(3 pi) = {2 / (3 * math.pi):.8f}")


def check_genus_two(law: TailLaw, count: int = 60_000_000, chunk: int = 4_000_000) -> None:
    """g = 2: rejection sampling of the Weyl measure on [0,pi]^2."""

    print("  g=2 tail law vs rejection sampling of the Weyl measure")
    accepted = []
    for _ in range(count // chunk):
        theta = RNG.uniform(0.0, math.pi, size=(chunk, 2))
        cosines = np.cos(theta)
        density = (cosines[:, 0] - cosines[:, 1]) ** 2 * (
            (1.0 - cosines[:, 0] ** 2) * (1.0 - cosines[:, 1] ** 2)
        )
        keep = RNG.uniform(0.0, 4.0, size=chunk) < density
        accepted.append(2.0 * cosines[keep].sum(axis=1))
    traces = np.concatenate(accepted)
    print(f"    accepted {traces.size} of {count} samples "
          f"(expected rate {math.factorial(2) / 2 ** (2 * 4 - 2):.4f})")
    for epsilon in (0.5, 1.0, 2.0, 3.0, 4.0):
        empirical = float((4.0 - traces < epsilon).mean())
        print(f"    eps={epsilon:>5.2f}  formula={law.cdf(epsilon):.6f}  "
              f"sampled={empirical:.6f}")
    print(f"    moments: E[T^2]={float((traces**2).mean()):.4f} (exact 1), "
          f"E[T^4]={float((traces**4).mean()):.4f} (exact 3)")


def check_arithmetic_tail(laws: dict[int, "TailLaw"]) -> None:
    """Compare the empirical tail of the arithmetic traces with the USp(2g) law.

    This is the equidistribution input to the whole extreme-value picture,
    tested directly rather than only through the maximum.
    """

    from t2_1_genus_scaling import FAMILIES, traces_hyperelliptic

    print("  empirical tail of -a_c/sqrt(q) vs the USp(2g) law (q = 1000003)")
    q = 1000003
    for genus, names in ((1, ("E1", "E1b", "E1c", "E1d")),
                         (2, ("H2", "H2b", "H2c", "H2d")),
                         (3, ("H3", "H3b", "H3c", "H3d"))):
        pooled = np.concatenate(
            [-traces_hyperelliptic(FAMILIES[name][1], q) for name in names]
        ) / math.sqrt(q)
        law = laws[genus]
        print(f"    g={genus}  ({pooled.size} traces)")
        for epsilon in (0.25, 0.5, 1.0, 2.0, 3.0)[: 2 + genus]:
            empirical = float((2 * genus - pooled < epsilon).mean())
            exact = law.cdf(epsilon)
            print(f"      eps={epsilon:>5.2f}  arithmetic={empirical:.7f}  "
                  f"USp(2g)={exact:.7f}  ratio={empirical / max(exact, 1e-300):.3f}")


def gaussian_max(samples: int) -> float:
    """Classical location of the maximum of n iid standard normals."""

    logarithm = math.log(samples)
    scale = math.sqrt(2.0 * logarithm)
    return scale - math.log(4.0 * math.pi * logarithm) / (2.0 * scale)


# --------------------------------------------------------------------------
def main() -> int:
    observed = load_observed()

    print("=" * 78)
    print("cross-checks on the exact USp(2g) lower-tail law")
    print("=" * 78)
    laws = {g: TailLaw(g) for g in (1, 2, 3, 4)}
    check_genus_one(laws[1])
    check_genus_two(laws[2])
    check_arithmetic_tail(laws)

    print()
    print("=" * 78)
    print("edge constants and predicted exponents")
    print("=" * 78)
    print(f"{'g':>3} {'d=dim USp(2g)':>14} {'edge exponent d/2':>18} "
          f"{'K_g':>14} {'deficit exponent':>17} {'crossover log q':>16}")
    for g, law in laws.items():
        print(f"{g:>3} {law.dimension:>14} {law.dimension / 2:>18.1f} "
              f"{law.edge_constant:>14.6e} {-2 / law.dimension:>17.4f} "
              f"{2 * g * g:>16}")

    print()
    print("=" * 78)
    print("observed vs predicted deficit  2g - max_c(-a_c)/sqrt(q)")
    print("=" * 78)
    rows: list[list] = []
    for g in (1, 2, 3, 4):
        law = laws[g]
        print(f"\ngenus {g} (2g = {2 * g}, families averaged: "
              f"{', '.join(observed[g]['families'])})")
        print(f"{'q':>8} {'obs deficit':>12} {'spread':>9} {'predicted E':>12} "
              f"{'predicted med':>14} {'edge asympt':>12} {'gauss max':>10} "
              f"{'obs mu':>8} {'pred mu':>8} {'obs abs':>9} {'pred abs':>9} "
              f"{'int gap':>8}")
        for q, values in sorted(observed[g]["deficits"].items()):
            mean = float(np.mean(values))
            spread = float(np.std(values))
            predicted = law.expected_deficit(q)
            median = law.median_deficit(q)
            edge = math.gamma(1 + 2 / law.dimension) * (q * law.edge_constant) ** (
                -2 / law.dimension
            )
            gauss = gaussian_max(q)
            # the same comparison in units of a_c, where a_c is an integer
            root = math.sqrt(q)
            integer_bound = math.floor(2 * g * root)
            gap = float(np.mean([integer_bound - m for m in observed[g]["maxima"][q]]))
            print(f"{q:>8} {mean:>12.4f} {spread:>9.4f} {predicted:>12.4f} "
                  f"{median:>14.4f} {edge:>12.4f} {gauss:>10.4f} "
                  f"{2 * g - mean:>8.4f} {2 * g - predicted:>8.4f} "
                  f"{mean * root:>9.2f} {predicted * root:>9.2f} {gap:>8.2f}")
            rows.append([g, q, mean, spread, predicted, median, edge, gauss,
                         mean * root, predicted * root, gap])

    print()
    print("=" * 78)
    print("fitted exponent of the deficit against q (log-log least squares)")
    print("=" * 78)
    print(f"{'g':>3} {'range of q':>22} {'fitted slope':>13} "
          f"{'edge -2/d':>10} {'gaussian slope':>15}")
    for g in (1, 2, 3, 4):
        law = laws[g]
        items = sorted(observed[g]["deficits"].items())
        big = [(q, float(np.mean(v))) for q, v in items if q >= 4000]
        xs = np.log([q for q, _ in big])
        ys = np.log([d for _, d in big])
        slope = float(np.polyfit(xs, ys, 1)[0])
        # local slope of the gaussian-regime prediction over the same window
        low, high = big[0][0], big[-1][0]
        gauss_slope = (
            math.log(2 * g - gaussian_max(high)) - math.log(2 * g - gaussian_max(low))
        ) / (math.log(high) - math.log(low)) if gaussian_max(high) < 2 * g else float("nan")
        print(f"{g:>3} {f'{low}..{high}':>22} {slope:>13.4f} "
              f"{-2 / law.dimension:>10.4f} {gauss_slope:>15.4f}")

    path = OUTPUT_DIRECTORY / "t2_1_extreme_value.csv"
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["genus", "q", "observed_deficit_mean", "observed_deficit_std",
                         "predicted_mean", "predicted_median", "edge_asymptotic",
                         "gaussian_max_location", "observed_deficit_in_a_units",
                         "predicted_deficit_in_a_units", "gap_to_integer_weil_bound"])
        writer.writerows(rows)
    print(f"\nwritten to {path.name}")
    return 0


def load_observed() -> dict[int, dict]:
    """Read the deficits produced by t2_1_genus_scaling.py."""

    path = OUTPUT_DIRECTORY / "t2_1_genus_scaling.csv"
    result: dict[int, dict] = {}
    with path.open(encoding="utf-8") as stream:
        for record in csv.DictReader(stream):
            genus = int(record["genus"])
            if genus == 0 or record["family"] == "XY":
                continue
            bucket = result.setdefault(
                genus, {"deficits": {}, "maxima": {}, "families": []}
            )
            if record["family"] not in bucket["families"]:
                bucket["families"].append(record["family"])
            bucket["deficits"].setdefault(int(record["q"]), []).append(
                float(record["deficit"])
            )
            bucket["maxima"].setdefault(int(record["q"]), []).append(int(record["m"]))
    return result


if __name__ == "__main__":
    raise SystemExit(main())
