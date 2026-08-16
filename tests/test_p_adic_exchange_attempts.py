import sys
from pathlib import Path
import unittest


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "analysis"))

from p_adic_exchange_attempts import (  # noqa: E402
    integral_three_to_two_mod_four,
    residue_exchange_rate,
    residue_signature,
)


class PAdicResidueRateTests(unittest.TestCase):
    def test_small_residue_signatures(self):
        self.assertEqual(residue_signature(2, 1), (8, 4, 4))
        self.assertEqual(residue_signature(2, 2), (4, 4, 4, 4))

    def test_known_mod_four_rates(self):
        q1 = residue_signature(4, 1)
        q2 = residue_signature(4, 2)
        self.assertAlmostEqual(
            residue_exchange_rate(q1, q2).rate,
            0.954242509439,
            places=11,
        )
        self.assertAlmostEqual(
            residue_exchange_rate(q2, q1).rate,
            0.991176112003,
            places=11,
        )

    def test_no_standard_integral_three_to_two_identity_mod_four(self):
        forward = integral_three_to_two_mod_four(1, 2)
        reverse = integral_three_to_two_mod_four(2, 1)
        self.assertFalse(forward.found)
        self.assertFalse(reverse.found)
        self.assertEqual(forward.pullback_forms, 136)
        self.assertEqual(reverse.pullback_forms, 256)


if __name__ == "__main__":
    unittest.main()
