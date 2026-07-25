#!/usr/bin/env python3
"""Generate Appendix B's 99-signature order and cycle exceptions."""

from __future__ import annotations

from itertools import combinations_with_replacement
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from exchange_matrix_extended import (  # noqa: E402
    comparison_graph,
    condensation_layers,
    signatures as original_signatures,
    strongly_connected_components,
)
from fn_complexity import ExchangeRateCache  # noqa: E402

OUTPUT_PATH = PROJECT_ROOT / "paper" / "appendix_b_signatures.tex"
TABLE_SIZE = 99


def extended_signatures() -> list[tuple[int, ...]]:
    """Return the current catalog plus length-two signatures up to eight."""

    result = original_signatures()
    result.extend(
        tuple(reversed(values))
        for values in combinations_with_replacement(range(1, 9), 2)
        if max(values) >= 7
    )
    # Singletons and all-ones signatures are exceptional infinite families and
    # are deliberately excluded from the finite table.
    return sorted(
        {
            signature
            for signature in result
            if len(signature) > 1 and signature[0] > 1
        },
        key=lambda signature: (len(signature), signature),
    )


def tex_signature(signature: tuple[int, ...]) -> str:
    return r"\(\sig{" + ",".join(map(str, signature)) + r"}\)"


def ordered_signatures(
    candidates: list[tuple[int, ...]],
    rates: dict[tuple[tuple[int, ...], tuple[int, ...]], float],
) -> tuple[
    list[tuple[int, ...]],
    dict[tuple[int, ...], set[tuple[int, ...]]],
]:
    outgoing, _ = comparison_graph(candidates, rates)
    components = strongly_connected_components(candidates, outgoing)
    layers = condensation_layers(components, outgoing)
    ordered = [
        signature
        for layer in layers
        for component_index in layer
        for signature in components[component_index]
    ]
    return ordered[:TABLE_SIZE], outgoing


def render_order_table(ordered: list[tuple[int, ...]]) -> list[str]:
    lines = [
        r"\begingroup",
        r"\footnotesize",
        r"\begin{longtable}{@{}r l@{\qquad}r l@{\qquad}r l@{}}",
        r"  \caption{The first 99 non-special signatures in the deterministic",
        r"  condensation-DAG order.}",
        r"  \label{tab:first-99-signatures}\\",
        r"    \toprule",
        r"    \(n\) & signature & \(n\) & signature & \(n\) & signature \\",
        r"    \midrule",
    ]
    for row in range(33):
        cells: list[str] = []
        for offset in (0, 33, 66):
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
            numbers = ",".join(str(value) for value in earlier_numbers)
            rendered_exceptions.append(
                rf"\(n_a\in\{{{numbers}\}}\):\newline "
                rf"\(\quad C(x\!\mid\!a)={forward_rate},\ "
                rf"C(a\!\mid\!x)={reverse_rate}\)"
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
    candidates = extended_signatures()
    cache = ExchangeRateCache()
    rates = {
        (implementer, implemented): cache.get(implementer, implemented)
        for implementer in candidates
        for implemented in candidates
    }
    cache.save()
    ordered, outgoing = ordered_signatures(candidates, rates)
    if len(ordered) != TABLE_SIZE:
        raise AssertionError(f"expected {TABLE_SIZE} signatures")

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
    print(f"candidates: {len(candidates)}")
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
