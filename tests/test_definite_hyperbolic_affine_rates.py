import importlib.util
import sys
from pathlib import Path
import unittest


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "analysis"))

from definite_hyperbolic_affine_rates import (  # noqa: E402
    finite_field_search,
    verify_certificate,
)


class ThreeToTwoCertificateTests(unittest.TestCase):
    def test_certificate_is_an_identity(self):
        self.assertTrue(verify_certificate())


class FiniteFieldSearchTests(unittest.TestCase):
    def test_two_split_blocks_from_three_anisotropic_copies(self):
        search = finite_field_search("definite", "hyperbolic", 2)
        self.assertTrue(search.found)
        self.assertTrue(search.target_span_is_pullback_free)
        self.assertEqual(search.pullbacks, 861)

    def test_no_two_anisotropic_blocks_from_three_split_copies(self):
        search = finite_field_search("hyperbolic", "definite", 2)
        self.assertFalse(search.found)
        self.assertTrue(search.target_span_is_pullback_free)
        self.assertEqual(search.pullbacks, 1641)
        self.assertLess(search.best_difference_rank, search.k)


@unittest.skipIf(
    importlib.util.find_spec("scipy") is None, "scipy is an analysis extra"
)
class RealScanTests(unittest.TestCase):
    def test_one_hyperbolic_block_needs_two_definite_copies(self):
        from definite_hyperbolic_affine_rates import real_probe

        self.assertFalse(real_probe("definite", "hyperbolic", 1, 1, 20, 1).feasible)
        self.assertTrue(real_probe("definite", "hyperbolic", 1, 2, 20, 2).feasible)


if __name__ == "__main__":
    unittest.main()
