import subprocess
import tempfile
import unittest
from pathlib import Path
from xml.etree import ElementTree

from fn_complexity.homogeneous_tensor_maps import (
    TENSOR_CASES,
    all_subspaces,
    compute_tensor_poset,
    invertible_matrices,
    minimum_orbit_counts_for_case6,
    rank_counts_for_case6,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]
CLI = PROJECT_ROOT / "cli" / "homogeneous_tensor_poset_cli"
IMAGE_DIRECTORY = PROJECT_ROOT / "paper_finite_fields_maps" / "images"


class HomogeneousTensorPosetTests(unittest.TestCase):
    def test_linear_groups_and_grassmannian_size(self):
        self.assertEqual(len(invertible_matrices(2)), 48)
        self.assertEqual(len(invertible_matrices(3)), 11_232)
        self.assertEqual(len(all_subspaces(6, 3)), 45_256)

    def test_exact_poset_counts_and_partitions(self):
        expected = {
            1: (5, 5),
            2: (7, 8),
            3: (50, 210),
            4: (26, 66),
            5: (19, 32),
        }
        for case_number, (class_count, cover_count) in expected.items():
            poset = compute_tensor_poset(case_number)
            self.assertEqual(len(poset.orbits), class_count)
            self.assertEqual(len(poset.covers), cover_count)
            self.assertEqual(
                sum(orbit.size for orbit in poset.orbits),
                TENSOR_CASES[case_number].tensor_count,
            )

    def test_case6_rank_counts_and_class_lower_bound(self):
        counts = rank_counts_for_case6()
        self.assertEqual(
            counts,
            {
                0: 1,
                1: 767_624,
                2: 45_325_126_704,
                3: 205_845_806_200_320,
            },
        )
        self.assertEqual(sum(counts.values()), 3**30)
        minima = minimum_orbit_counts_for_case6()
        exact_rank1 = len(compute_tensor_poset(4).orbits) - 1
        self.assertEqual(minima[3] + minima[2] + exact_rank1 + 1, 1_632_040)

    def test_checked_in_svgs_are_well_formed_and_complete(self):
        expected = {
            1: ("quadratic-ternary-form", 5, 5),
            2: ("quadratic-p1-map", 7, 8),
            3: ("quadratic-p2-map", 50, 210),
            4: ("cubic-ternary-form", 26, 66),
            5: ("cubic-p1-map", 19, 32),
            6: ("cubic-p2-map", 4, 3),
        }
        for case_number, (name, node_count, edge_count) in expected.items():
            path = (
                IMAGE_DIRECTORY
                / f"homogeneous-tensor-poset-q3-case{case_number}-{name}.svg"
            )
            ElementTree.parse(path)
            svg = path.read_text(encoding="utf-8")
            self.assertEqual(svg.count('class="node"'), node_count)
            self.assertEqual(svg.count('class="edge"'), edge_count)
        case6 = (
            IMAGE_DIRECTORY
            / "homogeneous-tensor-poset-q3-case6-cubic-p2-map.svg"
        ).read_text(encoding="utf-8")
        self.assertIn("not a Hasse diagram", case6)
        self.assertIn("1,632,040", case6)

    def test_cli_generates_selected_cases(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            result = subprocess.run(
                [
                    str(CLI),
                    "1",
                    "6",
                    "--output-dir",
                    temporary_directory,
                ],
                cwd=PROJECT_ROOT,
                text=True,
                capture_output=True,
                check=True,
            )
            self.assertIn("1\t729\t5\t5", result.stdout)
            self.assertIn("6\t205891132094649\t>=1632040\tN/A", result.stdout)


if __name__ == "__main__":
    unittest.main()
