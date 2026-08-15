#!/usr/bin/env python3
"""Generate Appendix B's stabilized 69-signature order and exceptions.

The candidate exhaustion is

    S_B = {a : len(a) >= 2, a_1 > 1, a_1 + 2 len(a) <= B}.

It contains every non-special signature for all sufficiently large ``B``.
The deterministic condensation-DAG order has the same first 69 entries for
``B=18`` and ``B=19`` (verified with the power-sum screening calculation).
Displayed rates and exception relations are recomputed here with the
high-accuracy persistent cache.
"""

from __future__ import annotations

from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from exchange_matrix_extended import (  # noqa: E402
    comparison_graph,
)
from fn_complexity import ExchangeRateCache  # noqa: E402

OUTPUT_PATH = PROJECT_ROOT / "paper" / "appendix_b_signatures.tex"
TABLE_SIZE = 69
STABILITY_BUDGETS = (18, 19)
STABILIZED_ORDER_PREFIX = (
    (2, 1),
    (2, 2),
    (3, 1),
    (2, 1, 1),
    (3, 2),
    (2, 2, 1),
    (3, 3),
    (2, 1, 1, 1),
    (2, 2, 2),
    (4, 1),
    (4, 2),
    (4, 3),
    (4, 4),
    (5, 1),
    (5, 2),
    (5, 3),
    (5, 4),
    (5, 5),
    (6, 1),
    (6, 2),
    (6, 3),
    (6, 4),
    (6, 5),
    (6, 6),
    (7, 1),
    (7, 2),
    (7, 3),
    (7, 4),
    (7, 5),
    (7, 6),
    (7, 7),
    (8, 1),
    (8, 2),
    (8, 3),
    (8, 4),
    (8, 5),
    (8, 6),
    (8, 7),
    (8, 8),
    (9, 1),
    (9, 2),
    (9, 3),
    (9, 4),
    (9, 5),
    (9, 6),
    (9, 7),
    (9, 8),
    (9, 9),
    (10, 1),
    (10, 2),
    (10, 3),
    (10, 4),
    (10, 5),
    (10, 6),
    (10, 7),
    (10, 8),
    (10, 9),
    (10, 10),
    (11, 1),
    (11, 2),
    (11, 3),
    (11, 4),
    (11, 5),
    (11, 6),
    (11, 7),
    (11, 8),
    (11, 9),
    (11, 10),
    (11, 11),
    (12, 1),
    (12, 2),
    (12, 3),
    (12, 4),
    (12, 5),
    (12, 6),
    (12, 7),
    (12, 8),
    (12, 9),
    (12, 10),
    (12, 11),
    (12, 12),
    (13, 1),
    (13, 2),
    (13, 3),
    (13, 4),
    (13, 5),
    (13, 6),
    (13, 7),
    (13, 8),
    (13, 9),
    (13, 10),
    (13, 11),
    (13, 12),
    (13, 13),
    (14, 1),
    (14, 2),
    (14, 3),
    (14, 4),
    (14, 5),
)


def tex_signature(signature: tuple[int, ...]) -> str:
    return r"\(\sig{" + ",".join(map(str, signature)) + r"}\)"


def render_order_table(ordered: list[tuple[int, ...]]) -> list[str]:
    if len(ordered) % 3:
        raise ValueError("the three-block table requires a multiple of three")
    rows_per_block = len(ordered) // 3
    lines = [
        r"\begingroup",
        r"\footnotesize",
        r"\begin{longtable}{@{}r l@{\qquad}r l@{\qquad}r l@{}}",
        rf"  \caption{{The stabilized first {len(ordered)} non-special "
        r"signatures in the",
        r"  deterministic condensation-DAG order.}",
        rf"  \label{{tab:first-{len(ordered)}-signatures}}\\",
        r"    \toprule",
        r"    \(n\) & signature & \(n\) & signature & \(n\) & signature \\",
        r"    \midrule",
    ]
    for row in range(rows_per_block):
        cells: list[str] = []
        for offset in (0, rows_per_block, 2 * rows_per_block):
            index = row + offset
            cells.extend([str(index + 1), tex_signature(ordered[index])])
        lines.append("    " + " & ".join(cells) + r" \\")
    lines.extend(
        [
            r"    \bottomrule",
            r"\end{longtable}",
            r"\endgroup",
        ]
    )
    return lines


def render_exception_table(
    ordered: list[tuple[int, ...]],
    outgoing: dict[tuple[int, ...], set[tuple[int, ...]]],
    rates: dict[tuple[tuple[int, ...], tuple[int, ...]], float],
) -> tuple[list[str], int]:
    number = {signature: index for index, signature in enumerate(ordered, 1)}
    exception_rows: list[tuple[int, tuple[int, ...], list[tuple[int, ...]]]] = []
    for index, signature in enumerate(ordered, 1):
        previous = [
            earlier
            for earlier in ordered[: index - 1]
            if earlier in outgoing[signature]
        ]
        if previous:
            exception_rows.append((index, signature, previous))

    lines = [
        r"\begin{longtable}{@{}r l >{\raggedright\arraybackslash}p{0.67\textwidth}@{}}",
        r"  \caption{Backward comparisons that obstruct a linear order.}",
        r"  \label{tab:signature-order-exceptions}\\",
        r"  \toprule",
        r"  \(n\) & \(x\) & exceptions with previous \(n_a\) values \\",
        r"  \midrule",
        r"  \endfirsthead",
        r"  \multicolumn{3}{@{}l}{\small\itshape Table "
        r"\thetable\ continued from the previous page}\\",
        r"  \toprule",
        r"  \(n\) & \(x\) & exceptions with previous \(n_a\) values \\",
        r"  \midrule",
        r"  \endhead",
        r"  \midrule",
        r"  \multicolumn{3}{r@{}}{\small continued on the next page}\\",
        r"  \endfoot",
        r"  \bottomrule",
        r"  \endlastfoot",
    ]
    for row_index, (index, signature, previous) in enumerate(exception_rows):
        grouped: dict[tuple[str, str], list[int]] = {}
        for earlier in previous:
            rate_pair = (
                f"{rates[signature, earlier]:.6f}",
                f"{rates[earlier, signature]:.6f}",
            )
            grouped.setdefault(rate_pair, []).append(number[earlier])

        rendered_exceptions = []
        for (forward_rate, reverse_rate), earlier_numbers in grouped.items():
            if len(earlier_numbers) == 1:
                earlier_number = earlier_numbers[0]
                earlier_signature = ",".join(
                    map(str, ordered[earlier_number - 1])
                )
                reference = (
                    rf"\(n_a={earlier_number},\ "
                    rf"a=\sig{{{earlier_signature}}}\)"
                )
            else:
                numbers = ",".join(str(value) for value in earlier_numbers)
                reference = rf"\(n_a\in\{{{numbers}\}}\)"
            rendered_exceptions.append(
                reference
                + r":\newline "
                rf"\(\quad C(x\!\to\!a)={forward_rate},\ "
                rf"C(a\!\to\!x)={reverse_rate}\)"
            )
        lines.append(
            f"  {index} & {tex_signature(signature)} & "
            + r"\newline ".join(rendered_exceptions)
            + r" \\"
        )
        if row_index < len(exception_rows) - 1:
            lines.append(r"  \cmidrule(lr){1-3}")
    lines.append(r"\end{longtable}")
    return lines, len(exception_rows)


def main() -> int:
    ordered = list(STABILIZED_ORDER_PREFIX[:TABLE_SIZE])
    if len(ordered) != TABLE_SIZE or len(set(ordered)) != TABLE_SIZE:
        raise AssertionError(
            f"the stabilized table must contain {TABLE_SIZE} signatures"
        )

    cache = ExchangeRateCache()
    rates = {
        (implementer, implemented): cache.get(implementer, implemented)
        for implementer in ordered
        for implemented in ordered
    }
    cache.save()
    outgoing, _ = comparison_graph(ordered, rates)

    exception_lines, exception_count = render_exception_table(
        ordered,
        outgoing,
        rates,
    )
    content = [
        "% Generated by analysis/appendix_b_signatures.py.",
        *render_order_table(ordered),
        "",
        *exception_lines,
        "",
    ]
    OUTPUT_PATH.write_text("\n".join(content), encoding="utf-8")
    print(f"stable budgets: {STABILITY_BUDGETS}")
    print(f"ordered signatures: {len(ordered)}")
    print(f"exception rows: {exception_count}")
    print(
        f"cache: {cache.hits} hits, {cache.misses} misses, "
        f"{len(cache)} total"
    )
    print(OUTPUT_PATH)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
