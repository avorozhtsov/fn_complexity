import subprocess
import tempfile
import unittest
from pathlib import Path
from xml.etree import ElementTree

from fn_complexity.finite_field_maps import (
    class_sizes,
    quadratic_map_count,
    quadratic_map_covers,
    quadratic_map_classes,
)
from fn_complexity.cubic_field_maps import (
    CUBIC_Q8_QUADRATIC_CLASSES,
    CUBIC_Q8_QUADRATIC_COVERS,
    cubic_q3_classes,
    cubic_q3_covers,
    cubic_q3_map_count,
    cubic_q8_generated_class_count,
    cubic_q8_map_count,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]
CLI = PROJECT_ROOT / "cli" / "finite_field_map_poset_cli"
CUBIC_CLI = PROJECT_ROOT / "cli" / "cubic_map_poset_cli"


class QuadraticMapPosetTests(unittest.TestCase):
    def test_q3_class_sizes(self):
        self.assertEqual(
            class_sizes(3),
            {
                "constant": 3,
                "linear": 24,
                "rank1": 72,
                "parabolic": 144,
                "split": 324,
                "anisotropic": 162,
            },
        )

    def test_q5_class_sizes(self):
        self.assertEqual(
            class_sizes(5),
            {
                "constant": 5,
                "linear": 120,
                "rank1": 600,
                "parabolic": 2400,
                "split": 7500,
                "anisotropic": 5000,
            },
        )

    def test_q9_class_sizes(self):
        self.assertEqual(
            class_sizes(9),
            {
                "constant": 9,
                "linear": 720,
                "rank1": 6480,
                "parabolic": 51840,
                "split": 262440,
                "anisotropic": 209952,
            },
        )

    def test_q11_and_q13_class_sizes(self):
        self.assertEqual(
            class_sizes(11),
            {
                "constant": 11,
                "linear": 1320,
                "rank1": 14520,
                "parabolic": 145200,
                "split": 878460,
                "anisotropic": 732050,
            },
        )

    def test_q25_class_sizes(self):
        self.assertEqual(
            class_sizes(25),
            {
                "constant": 25,
                "linear": 15600,
                "rank1": 390000,
                "parabolic": 9360000,
                "split": 121875000,
                "anisotropic": 112500000,
            },
        )
        self.assertEqual(
            class_sizes(13),
            {
                "constant": 13,
                "linear": 2184,
                "rank1": 28392,
                "parabolic": 340704,
                "split": 2399124,
                "anisotropic": 2056392,
            },
        )

    def test_q2_classes_are_boolean_function_classes(self):
        self.assertEqual(
            class_sizes(2),
            {"constant": 2, "linear": 6, "singleton": 8},
        )
        self.assertEqual(quadratic_map_count(2), 16)

    def test_even_characteristic_class_sizes(self):
        self.assertEqual(
            class_sizes(4),
            {
                "constant": 4,
                "linear": 60,
                "rank1": 60,
                "separable": 180,
                "parabolic": 720,
                "split": 1920,
                "anisotropic": 1152,
            },
        )
        self.assertEqual(
            class_sizes(8),
            {
                "constant": 8,
                "linear": 504,
                "rank1": 504,
                "separable": 3528,
                "parabolic": 28224,
                "split": 129024,
                "anisotropic": 100352,
            },
        )
        self.assertEqual(
            class_sizes(16),
            {
                "constant": 16,
                "linear": 4080,
                "rank1": 4080,
                "separable": 61200,
                "parabolic": 979200,
                "split": 8355840,
                "anisotropic": 7372800,
            },
        )

    def test_classes_partition_all_coefficient_vectors(self):
        for q in (3, 4, 5, 7, 8, 9, 11, 13, 16, 25):
            self.assertEqual(sum(item.size for item in quadratic_map_classes(q)), q**6)

    def test_odd_hasse_covers(self):
        self.assertEqual(
            set(quadratic_map_covers(3)),
            {
                ("parabolic", "linear"),
                ("parabolic", "rank1"),
                ("split", "linear"),
                ("split", "rank1"),
                ("anisotropic", "rank1"),
                ("linear", "constant"),
                ("rank1", "constant"),
            },
        )

    def test_even_hasse_covers(self):
        self.assertEqual(
            set(quadratic_map_covers(4)),
            {
                ("parabolic", "linear"),
                ("parabolic", "rank1"),
                ("parabolic", "separable"),
                ("split", "linear"),
                ("split", "rank1"),
                ("split", "separable"),
                ("anisotropic", "rank1"),
                ("anisotropic", "separable"),
                ("linear", "constant"),
                ("rank1", "constant"),
                ("separable", "constant"),
            },
        )

    def test_binary_hasse_covers(self):
        self.assertEqual(
            quadratic_map_covers(2),
            (("singleton", "linear"), ("linear", "constant")),
        )

    def test_cli_generates_requested_svgs(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            result = subprocess.run(
                [
                    str(CLI),
                    "2",
                    "4",
                    "8",
                    "16",
                    "--output-dir",
                    temporary_directory,
                ],
                cwd=PROJECT_ROOT,
                text=True,
                capture_output=True,
                check=True,
            )
            self.assertIn("2\t16\t3", result.stdout)
            self.assertIn("4\t4096\t7", result.stdout)
            self.assertIn("8\t262144\t7", result.stdout)
            self.assertIn("16\t16777216\t7", result.stdout)
            cases = (
                (2, "|class| = 8", 3, 2),
                (4, "|class| = 1,920", 7, 11),
                (8, "|class| = 129,024", 7, 11),
                (16, "|class| = 8,355,840", 7, 11),
            )
            for q, expected_size, node_count, edge_count in cases:
                svg = (
                    Path(temporary_directory) / f"quadratic-map-poset-q{q}.svg"
                ).read_text(encoding="utf-8")
                self.assertIn("<svg", svg)
                field = "F" + str(q).translate(
                    str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
                )
                self.assertIn(f"Quadratic maps {field}²", svg)
                self.assertIn(expected_size, svg)
                self.assertEqual(svg.count('class="node '), node_count)
                self.assertEqual(svg.count('class="edge"'), edge_count)

    def test_cli_rejects_non_prime_power(self):
        result = subprocess.run(
            [str(CLI), "6"],
            cwd=PROJECT_ROOT,
            text=True,
            capture_output=True,
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("prime power", result.stderr)

    def test_manual_real_and_complex_diagrams(self):
        cases = {
            "R": (6, 7, "real orbit dimension = 6"),
            "C": (5, 6, "complex orbit dimension = 6"),
        }
        for field, (node_count, edge_count, dimension_label) in cases.items():
            path = (
                PROJECT_ROOT
                / "paper_finite_fields_maps"
                / "images"
                / f"quadratic-map-poset-{field}.svg"
            )
            ElementTree.parse(path)
            svg = path.read_text(encoding="utf-8")
            self.assertEqual(svg.count('class="node '), node_count)
            self.assertEqual(svg.count('class="edge"'), edge_count)
            self.assertIn("|class| = 𝔠", svg)
            self.assertIn(dimension_label, svg)

    def test_manual_two_adic_diagram(self):
        path = (
            PROJECT_ROOT
            / "paper_finite_fields_maps"
            / "images"
            / "quadratic-map-poset-Q2.svg"
        )
        ElementTree.parse(path)
        svg = path.read_text(encoding="utf-8")
        self.assertEqual(svg.count('class="node '), 12)
        self.assertEqual(svg.count('class="edge"'), 13)
        self.assertEqual(svg.count("anisotropic · d ="), 7)
        self.assertIn("split · d = −1", svg)
        self.assertIn("seven pairwise-incomparable anisotropic classes", svg)

    def test_manual_three_adic_diagram(self):
        path = (
            PROJECT_ROOT
            / "paper_finite_fields_maps"
            / "images"
            / "quadratic-map-poset-q3-adic.svg"
        )
        ElementTree.parse(path)
        svg = path.read_text(encoding="utf-8")
        self.assertEqual(svg.count('class="node '), 8)
        self.assertEqual(svg.count('class="edge"'), 9)
        self.assertEqual(svg.count("anisotropic · d ="), 3)
        self.assertIn("split · d = −1", svg)
        self.assertIn("three pairwise-incomparable anisotropic classes", svg)


class CubicMapPosetTests(unittest.TestCase):
    def test_affine_input_classes_partition_q3_cubic_functions(self):
        classes = cubic_q3_classes("linear")
        self.assertEqual(cubic_q3_map_count(), 6561)
        self.assertEqual(len(classes), 14)
        self.assertEqual(sum(item.size for item in classes), 6561)
        self.assertEqual(
            sorted(item.size for item in classes),
            [
                3,
                24,
                72,
                144,
                162,
                216,
                324,
                432,
                432,
                648,
                648,
                864,
                1296,
                1296,
            ],
        )
        self.assertEqual(len(cubic_q3_covers("linear")), 22)

    def test_quadratic_input_generated_preorder_is_a_three_class_chain(self):
        classes = cubic_q3_classes("quadratic")
        self.assertEqual(
            {item.key: item.size for item in classes},
            {"constant": 3, "two-valued": 504, "surjective": 6054},
        )
        self.assertEqual(
            cubic_q3_covers("quadratic"),
            (("surjective", "two-valued"), ("two-valued", "constant")),
        )

    def test_cubic_cli_generates_both_svg_diagrams(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            result = subprocess.run(
                [str(CUBIC_CLI), "--output-dir", temporary_directory],
                cwd=PROJECT_ROOT,
                text=True,
                capture_output=True,
                check=True,
            )
            self.assertIn("linear\t6561\t14", result.stdout)
            self.assertIn("quadratic\t6561\t3", result.stdout)
            cases = (
                ("linear", 14, 22, "|class| = 1,296"),
                ("quadratic", 3, 2, "|class| = 6,054"),
            )
            for processor_case, node_count, edge_count, size_label in cases:
                path = (
                    Path(temporary_directory)
                    / f"cubic-map-poset-q3-{processor_case}-input.svg"
                )
                ElementTree.parse(path)
                svg = path.read_text(encoding="utf-8")
                self.assertEqual(svg.count('class="node '), node_count)
                self.assertEqual(svg.count('class="edge'), edge_count)
                self.assertIn(size_label, svg)

    def test_q8_quadratic_input_class_sizes_partition_all_maps(self):
        self.assertEqual(cubic_q8_map_count(), 1_073_741_824)
        self.assertEqual(cubic_q8_generated_class_count(), 110)
        self.assertEqual(
            sum(item.total_size for item in CUBIC_Q8_QUADRATIC_CLASSES),
            cubic_q8_map_count(),
        )
        self.assertEqual(len(CUBIC_Q8_QUADRATIC_CLASSES), 23)
        self.assertEqual(len(CUBIC_Q8_QUADRATIC_COVERS), 57)

    def test_q8_compressed_hasse_diagram_is_acyclic(self):
        adjacency = {
            item.key: set() for item in CUBIC_Q8_QUADRATIC_CLASSES
        }
        for source, target in CUBIC_Q8_QUADRATIC_COVERS:
            adjacency[source].add(target)

        def descendants(start):
            seen = set()
            stack = list(adjacency[start])
            while stack:
                node = stack.pop()
                self.assertNotEqual(node, start)
                if node not in seen:
                    seen.add(node)
                    stack.extend(adjacency[node])
            return seen

        for item in CUBIC_Q8_QUADRATIC_CLASSES:
            descendants(item.key)

    def test_cubic_cli_generates_q8_quadratic_svg(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            result = subprocess.run(
                [
                    str(CUBIC_CLI),
                    "--q",
                    "8",
                    "--case",
                    "quadratic",
                    "--output-dir",
                    temporary_directory,
                ],
                cwd=PROJECT_ROOT,
                text=True,
                capture_output=True,
                check=True,
            )
            self.assertIn("quadratic\t1073741824\t110", result.stdout)
            path = (
                Path(temporary_directory)
                / "cubic-map-poset-q8-quadratic-input.svg"
            )
            ElementTree.parse(path)
            svg = path.read_text(encoding="utf-8")
            self.assertEqual(svg.count('class="node"'), 23)
            self.assertEqual(svg.count('class="edge"'), 57)
            self.assertIn("24 classes × 12,644,352 each", svg)


if __name__ == "__main__":
    unittest.main()
