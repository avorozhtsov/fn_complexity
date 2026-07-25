import math
import unittest
from fractions import Fraction

from fn_complexity import (
    continued_fraction_convergents,
    enumerate_hyperbolas,
    first_uncovered_greedy_cover,
    k_max,
    maximum_coverage_curve_per_slope,
    minimum_curve_cover,
    point_curve_map,
)


class HyperbolaAnalysisTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.values = [
            (n, k_max((2, 2), (3, 1), n))
            for n in range(1, 101)
        ]
        cls.rate = math.log(2) / math.log(3)
        cls.convergents = continued_fraction_convergents(
            cls.rate, maximum_denominator=100
        )
        cls.curves = enumerate_hyperbolas(cls.values, cls.convergents)

    def test_expected_convergents(self):
        self.assertEqual(
            self.convergents,
            (
                Fraction(1, 1),
                Fraction(1, 2),
                Fraction(2, 3),
                Fraction(5, 8),
                Fraction(12, 19),
                Fraction(41, 65),
                Fraction(53, 84),
            ),
        )

    def test_catalog_and_mapping(self):
        self.assertEqual(len(self.curves), 259)
        mapping = point_curve_map(self.values, self.curves)
        self.assertTrue(all(len(mapping[n]) == 7 for n in range(1, 101)))

    def test_first_point_greedy_is_not_minimum(self):
        greedy = first_uncovered_greedy_cover(self.values, self.curves)
        minimum = minimum_curve_cover(self.values, self.curves)
        self.assertEqual(len(greedy), 15)
        self.assertEqual(len(minimum.curves), 12)
        self.assertTrue(
            all(curve.slope == Fraction(5, 8) for curve in minimum.curves)
        )
        covered = {n for curve in minimum.curves for n in curve.points}
        self.assertEqual(covered, set(range(1, 101)))

    def test_one_maximum_curve_per_convergent(self):
        selected = maximum_coverage_curve_per_slope(self.curves)
        self.assertEqual(
            [(curve.slope, curve.offset, len(curve.points)) for curve in selected],
            [
                (Fraction(1), Fraction(2), 3),
                (Fraction(1, 2), Fraction(-1, 2), 4),
                (Fraction(2, 3), Fraction(5, 3), 10),
                (Fraction(5, 8), Fraction(1, 4), 13),
                (Fraction(12, 19), Fraction(3, 19), 6),
                (Fraction(41, 65), Fraction(2, 65), 2),
                (Fraction(53, 84), Fraction(1, 21), 2),
            ],
        )


if __name__ == "__main__":
    unittest.main()
