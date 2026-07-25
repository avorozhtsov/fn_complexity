"""Exchange rates between signatures of finite onto maps.

A signature is the multiset of non-zero fiber sizes of an onto finite map.
If ``implemented`` and ``implementer`` are signatures, ``implemented`` is
implemented by ``implementer`` precisely when its decreasing rearrangement is
componentwise no larger than an equally long prefix of the decreasing
rearrangement of ``implementer``.

Tensor (Cartesian) powers multiply fiber sizes.  The asymptotic exchange rate
is defined operationally by

    C(implementer | implemented)
      = lim_{n -> infinity} max{k : implemented**k <= implementer**n} / n.

The power-sum spectrum gives the equivalent computational characterization:

    C(implementer | implemented)
      = inf_{beta >= 0} log(sum(implementer_i**beta)) /
                          log(sum(implemented_i**beta)).

Here ``beta = 1 / T`` is inverse physical temperature.  The value at
``beta=infinity`` (``T=0``) is interpreted as
log(max(target)) / log(max(source)).
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
import math
from typing import Iterable, Mapping, Sequence

Signature = tuple[int, ...]
CompressedSignature = Mapping[int, int]


@dataclass(frozen=True)
class ExchangeRateResult:
    """A rate together with the limiting inverse temperature ``beta``."""

    rate: float
    beta: float

    @property
    def attained_at_infinity(self) -> bool:
        """Whether the limiting point is ``beta=infinity`` (``T=0``)."""

        return math.isinf(self.beta)

    @property
    def temperature(self) -> float:
        """Return the physical temperature ``T=1/beta``."""

        if self.beta == 0.0:
            return math.inf
        if math.isinf(self.beta):
            return 0.0
        return 1.0 / self.beta


@dataclass(frozen=True)
class ThermodynamicPoint:
    """A point on the Gibbs entropy-energy curve (natural-log units)."""

    temperature: float
    energy: float
    entropy: float
    probabilities: tuple[float, ...]


def normalize_signature(signature: Iterable[int]) -> Signature:
    """Validate and return a signature in decreasing order.

    Booleans and non-integral values are rejected: fiber cardinalities must be
    positive integers, and an onto map has at least one fiber.
    """

    values = tuple(signature)
    if not values:
        raise ValueError("a signature must contain at least one fiber")
    for value in values:
        if isinstance(value, bool) or not isinstance(value, int):
            raise TypeError("fiber sizes must be integers")
        if value <= 0:
            raise ValueError("fiber sizes must be positive")
    return tuple(sorted(values, reverse=True))


def _log_partition(signature: Sequence[int], beta: float) -> float:
    """Compute ``log(sum(n_i**beta))`` without overflow."""

    if beta < 0 or math.isnan(beta):
        raise ValueError("beta must be non-negative")
    if math.isinf(beta):
        return math.inf if signature[0] > 1 else math.log(len(signature))
    logs = [math.log(value) for value in signature]
    largest = beta * logs[0]
    return largest + math.log(
        sum(math.exp(beta * log_value - largest) for log_value in logs)
    )


def gibbs_point(signature: Iterable[int], temperature: float) -> ThermodynamicPoint:
    """Return ``(E(T), H(T))`` for energies ``epsilon_i = -log(n_i)``.

    Probabilities obey ``p_i proportional to exp(-epsilon_i / T)``.  At
    ``T=0`` they are uniform on the largest fibers (ground states); at
    ``T=infinity`` they are uniform on all fibers.
    """

    a = normalize_signature(signature)
    if temperature < 0 or math.isnan(temperature):
        raise ValueError("temperature must be non-negative")
    if math.isinf(temperature):
        probabilities = tuple(1.0 / len(a) for _ in a)
    elif temperature == 0:
        ground_count = sum(value == a[0] for value in a)
        probabilities = tuple(
            1.0 / ground_count if value == a[0] else 0.0 for value in a
        )
    else:
        inverse_temperature = 1.0 / temperature
        log_weights = [inverse_temperature * math.log(value) for value in a]
        largest = max(log_weights)
        weights = [math.exp(log_weight - largest) for log_weight in log_weights]
        total = sum(weights)
        probabilities = tuple(weight / total for weight in weights)

    energies = tuple(-math.log(value) for value in a)
    energy = sum(probability * level for probability, level in zip(probabilities, energies))
    entropy = -sum(
        probability * math.log(probability)
        for probability in probabilities
        if probability > 0
    )
    return ThermodynamicPoint(temperature, energy, entropy, probabilities)


def gibbs_curve(
    signature: Iterable[int], temperatures: Iterable[float]
) -> tuple[ThermodynamicPoint, ...]:
    """Evaluate the entropy-energy curve at the requested temperatures."""

    a = normalize_signature(signature)
    return tuple(gibbs_point(a, temperature) for temperature in temperatures)


def _ratio(implemented: Signature, implementer: Signature, beta: float) -> float:
    denominator = _log_partition(implemented, beta)
    numerator = _log_partition(implementer, beta)
    if denominator == 0.0:
        return math.inf
    return numerator / denominator


def _infinity_ratio(implemented: Signature, implementer: Signature) -> float:
    implemented_log_max = math.log(implemented[0])
    if implemented_log_max == 0.0:
        return math.inf
    return math.log(implementer[0]) / implemented_log_max


def _search_horizon(implemented: Signature, implementer: Signature) -> float:
    """Choose beta where every non-maximal Gibbs weight is negligible."""

    gaps: list[float] = []
    for signature in (implemented, implementer):
        top = math.log(signature[0])
        gaps.extend(top - math.log(value) for value in signature if value != signature[0])
    if not gaps:
        return 64.0
    # exp(-36) is below 3e-16.  Cap only to avoid an unreasonable grid for
    # nearly equal huge Python integers; stable evaluation remains available.
    return min(36.0 / min(gaps), 1.0e7)


def _golden_minimum(function, left: float, right: float, tolerance: float) -> tuple[float, float]:
    """Dependency-free bounded minimization on a bracketed basin."""

    inverse_phi = (math.sqrt(5.0) - 1.0) / 2.0
    x1 = right - inverse_phi * (right - left)
    x2 = left + inverse_phi * (right - left)
    y1, y2 = function(x1), function(x2)
    for _ in range(160):
        if right - left <= tolerance * (1.0 + abs(x1) + abs(x2)):
            break
        if y1 <= y2:
            right, x2, y2 = x2, x1, y1
            x1 = right - inverse_phi * (right - left)
            y1 = function(x1)
        else:
            left, x1, y1 = x1, x2, y2
            x2 = left + inverse_phi * (right - left)
            y2 = function(x2)
    return (x1, y1) if y1 <= y2 else (x2, y2)


def exchange_rate_result(
    implementer: Iterable[int],
    implemented: Iterable[int],
    *,
    grid_size: int = 1024,
    tolerance: float = 1e-12,
) -> ExchangeRateResult:
    """Compute the operational rate ``C(implementer | implemented)``.

    This is the limit of ``k_max(n) / n``, where ``implemented**k`` must be
    implemented by ``implementer**n``.  The implementation evaluates its
    equivalent power-sum infimum characterization.

    ``beta`` is the non-negative exponent in the power-sum formula and equals
    inverse physical temperature ``1/T`` for energies ``-log(n_i)``.  It is
    infinity (``T=0``) when the largest-fiber constraint is limiting.

    The one-fiber identity signature ``(1,)`` has zero complexity and can be
    produced at arbitrarily large rate from every non-empty implementer.
    """

    a = normalize_signature(implemented)
    b = normalize_signature(implementer)
    if grid_size < 16:
        raise ValueError("grid_size must be at least 16")
    if tolerance <= 0:
        raise ValueError("tolerance must be positive")

    if a == (1,):
        return ExchangeRateResult(math.inf, 0.0)

    # If all source fibers have size one, the only resource is their number.
    if a[0] == 1:
        rate = math.log(len(b)) / math.log(len(a))
        return ExchangeRateResult(rate, 0.0)

    horizon = _search_horizon(a, b)
    # A geometric grid resolves both near-zero and large-beta behavior, including
    # signatures whose distinct fiber sizes are extremely close.
    smallest_positive = min(1.0e-10, horizon / grid_size)
    log_span = math.log(horizon / smallest_positive)
    points = [0.0] + [
        smallest_positive * math.exp(log_span * index / (grid_size - 1))
        for index in range(grid_size)
    ]
    values = [_ratio(a, b, point) for point in points]

    candidates: list[tuple[float, float]] = [
        (0.0, values[0]),
        (math.inf, _infinity_ratio(a, b)),
    ]
    for index in range(1, len(points) - 1):
        if values[index] <= values[index - 1] and values[index] <= values[index + 1]:
            candidates.append(
                _golden_minimum(
                    lambda beta: _ratio(a, b, beta),
                    points[index - 1],
                    points[index + 1],
                    tolerance,
                )
            )

    beta, rate = min(candidates, key=lambda item: item[1])
    if rate < 0 and rate > -tolerance:
        rate = 0.0
    return ExchangeRateResult(rate, beta)


def exchange_rate(
    implementer: Iterable[int],
    implemented: Iterable[int],
    **kwargs,
) -> float:
    """Return only ``C(implementer | implemented)``."""

    return exchange_rate_result(implementer, implemented, **kwargs).rate


def power_signature(signature: Iterable[int], exponent: int) -> Counter[int]:
    """Return a tensor power as ``fiber_size -> multiplicity``.

    The compressed representation permits exact checks substantially beyond
    what expanding all ``len(signature)**exponent`` fibers would allow.
    """

    base = normalize_signature(signature)
    if isinstance(exponent, bool) or not isinstance(exponent, int):
        raise TypeError("exponent must be an integer")
    if exponent < 0:
        raise ValueError("exponent must be non-negative")

    return _power_signature_accumulated(base, exponent)


def _power_signature_accumulated(signature: Signature, exponent: int) -> Counter[int]:
    result: Counter[int] = Counter({1: 1})
    base_counts = Counter(signature)
    for _ in range(exponent):
        next_result: Counter[int] = Counter()
        for old_size, old_count in result.items():
            for base_size, base_count in base_counts.items():
                next_result[old_size * base_size] += old_count * base_count
        result = next_result
    return result


def implements(implemented: CompressedSignature, implementer: CompressedSignature) -> bool:
    """Check exact implementation for compressed signatures.

    Each mapping associates a fiber size with the number of such fibers.
    """

    if any(size <= 0 or count <= 0 for size, count in implemented.items()):
        raise ValueError("implemented sizes and multiplicities must be positive")
    if any(size <= 0 or count <= 0 for size, count in implementer.items()):
        raise ValueError("implementer sizes and multiplicities must be positive")

    target_levels = sorted(implementer.items(), reverse=True)
    target_index = 0
    available = 0
    for source_size, source_count in sorted(implemented.items(), reverse=True):
        while target_index < len(target_levels) and target_levels[target_index][0] >= source_size:
            available += target_levels[target_index][1]
            target_index += 1
        if available < source_count:
            return False
        available -= source_count
    return True


def k_max(implementer: Iterable[int], implemented: Iterable[int], n: int) -> int:
    """Compute max k such that ``implemented**k`` is implemented by ``implementer**n``.

    This verifier is intended for modest powers.  Feasibility is monotone in
    k (restrict an implementation to a slice with one fixed input coordinate),
    so the largest feasible value is found by binary search.
    """

    a = normalize_signature(implemented)
    b = normalize_signature(implementer)
    if isinstance(n, bool) or not isinstance(n, int):
        raise TypeError("n must be an integer")
    if n < 0:
        raise ValueError("n must be non-negative")
    if a == (1,):
        raise ValueError("k_max is unbounded for implemented signature (1,)")

    rate_upper = min(
        n * math.log(len(b)) / math.log(len(a)) if len(a) > 1 else math.inf,
        (
            n * math.log(b[0]) / math.log(a[0])
            if a[0] > 1
            else math.inf
        ),
    )
    upper = max(0, math.floor(rate_upper + 1e-12))
    # A uniform implementer has len(b)**n equal fibers, each of size b[0]**n.
    # The count and largest-fiber bounds are then jointly sufficient.
    if all(size == b[0] for size in b):
        return upper

    target_power = _power_signature_accumulated(b, n)
    feasible, infeasible = 0, upper + 1
    while infeasible - feasible > 1:
        candidate = (feasible + infeasible) // 2
        if implements(_power_signature_accumulated(a, candidate), target_power):
            feasible = candidate
        else:
            infeasible = candidate
    return feasible
