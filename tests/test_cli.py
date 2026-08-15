import subprocess
import tempfile
import unittest
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
CLI = PROJECT_ROOT / "cli" / "kmax_cli"
CLUSTER_CLI = PROJECT_ROOT / "cli" / "cluster_cli"


class KMaxCliTests(unittest.TestCase):
    def test_generates_svg_and_table(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            output = Path(temporary_directory) / "plot.svg"
            result = subprocess.run(
                [
                    str(CLI),
                    "{2, 2}",
                    "3,1",
                    "--n-max",
                    "5",
                    "--output",
                    str(output),
                ],
                cwd=PROJECT_ROOT,
                text=True,
                capture_output=True,
                check=True,
            )
            self.assertIn("5\t3\t0.600000000", result.stdout)
            self.assertIn("C(g -> f)", result.stdout)
            svg = output.read_text(encoding="utf-8")
            self.assertIn("<svg", svg)
            self.assertIn("exact cₙ = kₘₐₓ(n)/n", svg)
            self.assertIn("One maximum-coverage curve per convergent", svg)
            self.assertIn("exact polyline contains every point", svg)
            self.assertIn("log 2 / log 3", svg)

    def test_rejects_unbounded_case(self):
        result = subprocess.run(
            [str(CLI), "2,2", "1"],
            cwd=PROJECT_ROOT,
            text=True,
            capture_output=True,
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("unbounded", result.stderr)

    def test_selection_range_is_independent_of_plot_range(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            output = Path(temporary_directory) / "plot.svg"
            result = subprocess.run(
                [
                    str(CLI),
                    "2,2",
                    "3,1",
                    "--n-max",
                    "100",
                    "--hyperbola-n-max",
                    "99999",
                    "--max-convergents",
                    "7",
                    "--output",
                    str(output),
                    "--quiet",
                ],
                cwd=PROJECT_ROOT,
                text=True,
                capture_output=True,
                check=True,
            )
            self.assertIn("selected using n=1..99999", result.stdout)
            self.assertIn("c_n = 5/8 + 19/(8n)  (22 points", result.stdout)
            svg = output.read_text(encoding="utf-8")
            self.assertIn("cₙ=53/84 − 43/(42n) · 527 pts", svg)

    def test_half_branches_show_formulas(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            output = Path(temporary_directory) / "branches.svg"
            subprocess.run(
                [
                    str(CLI),
                    "2,2",
                    "3,1",
                    "--n-max",
                    "100",
                    "--half-branches",
                    "14",
                    "--output",
                    str(output),
                    "--quiet",
                ],
                cwd=PROJECT_ROOT,
                text=True,
                capture_output=True,
                check=True,
            )
            svg = output.read_text(encoding="utf-8")
            self.assertIn("First 14 upward e=1/2 branches", svg)
            self.assertIn("cₙ=1/2 · 3 pts", svg)
            self.assertIn("cₙ=1/2 + 1/(2n) · 4 pts", svg)
            self.assertIn("cₙ=1/2 + 1/n · 4 pts", svg)
            self.assertIn("cₙ=1/2 + 3/(2n) · 4 pts", svg)
            self.assertIn("cₙ=1/2 + 13/(2n)", svg)

    def test_unit_cover_partitions_all_plotted_points(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            output = Path(temporary_directory) / "unit-cover.svg"
            subprocess.run(
                [
                    str(CLI),
                    "3,1",
                    "2,2",
                    "--n-max",
                    "99",
                    "--unit-cover",
                    "--output",
                    str(output),
                    "--quiet",
                ],
                cwd=PROJECT_ROOT,
                text=True,
                capture_output=True,
                check=True,
            )
            svg = output.read_text(encoding="utf-8")
            self.assertIn("e=1 partition (all plotted points)", svg)
            for offset in range(1, 7):
                self.assertIn(f"cₙ=1 − {offset}/n", svg)
            self.assertNotIn("cₙ=1 − 7/n", svg)

    def test_inverse_limit_gap_plot(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            output = Path(temporary_directory) / "inverse-gap.svg"
            subprocess.run(
                [
                    str(CLI),
                    "3,1",
                    "2,2",
                    "--n-max",
                    "99",
                    "--unit-cover",
                    "--inverse-limit-gap",
                    "--output",
                    str(output),
                    "--quiet",
                ],
                cwd=PROJECT_ROOT,
                text=True,
                capture_output=True,
                check=True,
            )
            svg = output.read_text(encoding="utf-8")
            self.assertIn("Yₙ = −1/(cₙ − C)", svg)
            self.assertIn("exact Yₙ = −1/(cₙ−C)", svg)
            self.assertIn("inverse limit gap Yₙ", svg)
            self.assertIn("Yₙ=−1/(1 − 1/n−C)", svg)
            self.assertNotIn('stroke="#f59e0b" stroke-width="3"', svg)


class ClusterCliTests(unittest.TestCase):
    def test_stops_after_requested_number_of_members(self):
        result = subprocess.run(
            [
                str(CLUSTER_CLI),
                "{3, 1, 1}",
                "--n-max",
                "3",
                "--max-b",
                "10",
                "--grid-size",
                "64",
                "--quiet",
            ],
            cwd=PROJECT_ROOT,
            text=True,
            capture_output=True,
            check=True,
        )
        self.assertIn("target\t{3,1,1}", result.stdout)
        self.assertIn("reported\t3", result.stdout)
        self.assertIn("n_max\t3 (reached)", result.stdout)

    def test_rejects_shell_that_cannot_contain_target(self):
        result = subprocess.run(
            [str(CLUSTER_CLI), "3,1,1", "--max-b", "8", "--quiet"],
            cwd=PROJECT_ROOT,
            text=True,
            capture_output=True,
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("beyond --max-b 8", result.stderr)

    def test_strict_relation_is_reported(self):
        result = subprocess.run(
            [
                str(CLUSTER_CLI),
                "3,1",
                "--n-max",
                "2",
                "--max-b",
                "8",
                "--grid-size",
                "64",
                "--strict",
                "--quiet",
            ],
            cwd=PROJECT_ROOT,
            text=True,
            capture_output=True,
            check=True,
        )
        self.assertIn("relation\tstrict >", result.stdout)


if __name__ == "__main__":
    unittest.main()
