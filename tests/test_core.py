import math
import unittest

from fn_complexity import (
    exchange_rate,
    exchange_rate_result,
    gibbs_point,
    implements,
    k_max,
    normalize_signature,
    power_signature,
)


class ExchangeRateTests(unittest.TestCase):
    def test_three_one_implements_two_two(self):
        result = exchange_rate_result((3, 1), (2, 2))
        self.assertAlmostEqual(result.rate, 0.965384441732, places=10)
        self.assertAlmostEqual(result.beta, 0.403679910968, places=7)
        self.assertAlmostEqual(result.temperature, 2.477210266, places=7)

    def test_two_two_implements_three_one(self):
        result = exchange_rate_result((2, 2), (3, 1))
        self.assertAlmostEqual(result.rate, math.log(2) / math.log(3), places=11)
        self.assertTrue(result.attained_at_infinity)

    def test_gibbs_endpoints(self):
        cold = gibbs_point((3, 1), 0.0)
        hot = gibbs_point((3, 1), math.inf)
        self.assertEqual(cold.probabilities, (1.0, 0.0))
        self.assertAlmostEqual(cold.energy, -math.log(3))
        self.assertEqual(cold.entropy, 0.0)
        self.assertEqual(hot.probabilities, (0.5, 0.5))
        self.assertAlmostEqual(hot.entropy, math.log(2))

    def test_all_ones_and_identity_edges(self):
        self.assertTrue(math.isinf(exchange_rate((1,), (1,))))
        self.assertTrue(math.isinf(exchange_rate((3, 1), (1,))))
        self.assertEqual(exchange_rate((1,), (2,)), 0.0)
        self.assertAlmostEqual(exchange_rate((1, 1, 1), (1, 1)), math.log(3) / math.log(2))
        self.assertEqual(exchange_rate((1,), (1, 1)), 0.0)

    def test_validation(self):
        with self.assertRaises(ValueError):
            normalize_signature(())
        with self.assertRaises(ValueError):
            normalize_signature((2, 0))
        with self.assertRaises(TypeError):
            normalize_signature((2, 1.0))
        with self.assertRaises(TypeError):
            normalize_signature((True,))


class ExactImplementationTests(unittest.TestCase):
    def test_power_signature_accumulates_coincident_products(self):
        self.assertEqual(power_signature((2, 1), 2), {4: 1, 2: 2, 1: 1})

    def test_implementation_order(self):
        self.assertTrue(implements({3: 1, 1: 1}, {4: 1, 2: 1}))
        self.assertFalse(implements({3: 1, 2: 1}, {4: 1, 1: 1}))

    def test_first_rate_has_exact_floor_formula(self):
        for n in range(1, 25):
            self.assertEqual(k_max((2, 2), (3, 1), n), math.floor(n * math.log(2) / math.log(3)))

    def test_reverse_rate_converges(self):
        ratios = [k_max((3, 1), (2, 2), n) / n for n in (20, 40, 80, 120)]
        self.assertLess(abs(ratios[-1] - exchange_rate((3, 1), (2, 2))), 0.04)
        self.assertGreater(ratios[-1], 0.90)


if __name__ == "__main__":
    unittest.main()
