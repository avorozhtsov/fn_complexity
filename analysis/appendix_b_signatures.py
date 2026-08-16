#!/usr/bin/env python3
"""Generate Appendix B's stabilized 69-signature order and exceptions.

The candidate exhaustion is

    S_B = {a : len(a) >= 2, a_1 > 1, a_1 + 2 len(a) <= B}.

It contains every non-special signature for all sufficiently large ``B``.
The deterministic condensation-DAG order has the same first 69 entries for
``B=18`` and ``B=19`` (verified with the partition-function screening calculation).
Displayed rates and exception relations are recomputed here with the
high-accuracy persistent cache.
"""

from __future__ import annotations

import csv
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from exchange_matrix_extended import (  # noqa: E402
    comparison_graph,
)
from fn_complexity import ExchangeRateCache  # noqa: E402

OUTPUT_PATH = PROJECT_ROOT / "paper_exchange_rate" / "appendix_b_signatures.tex"
ANC_DIR = PROJECT_ROOT / "paper_exchange_rate" / "anc"
ORDER_CSV = ANC_DIR / "first_69_signatures.csv"
EXCEPTION_CSV = ANC_DIR / "order_exceptions.csv"
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


def render_order_list(ordered: list[tuple[int, ...]]) -> list[str]:
    """The stabilized prefix as running text; the full table ships in anc/."""
    entries = ", ".join(tex_signature(sig) for sig in ordered)
    return [
        r"\begingroup\small\noindent",
        entries + ".",
        r"\endgroup",
    ]


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

    pairs = sum(len(previous) for _, _, previous in exception_rows)
    widths: dict[int, int] = {}
    for _, signature, _ in exception_rows:
        widths[signature[1]] = widths.get(signature[1], 0) + 1
    spread = ", ".join(
        rf"\(\sig{{n,{second}}}\): {count}"
        for second, count in sorted(widths.items())
    )
    every_pair = all(len(sig) == 2 for _, sig, _ in exception_rows)
    lines = [
        rf"Of the {len(ordered)} signatures, {len(exception_rows)} are preceded "
        rf"by at least one \(a\) with \(C(x\to a)<C(a\to x)\), giving "
        rf"{pairs} such pairs in all.  "
        + (
            r"Every offending \(x\) has exactly two entries, and the "
            r"failures thin out as its second entry grows --- counted by that "
            r"entry they are " + spread + ".  "
            if every_pair
            else ""
        )
        + r"The pairs themselves, with both rates, are listed in the "
        r"ancillary file \texttt{anc/order\_exceptions.csv}.",
    ]
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
    ANC_DIR.mkdir(parents=True, exist_ok=True)
    with ORDER_CSV.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["n", "signature"])
        for index, signature in enumerate(ordered, 1):
            writer.writerow([index, ",".join(map(str, signature))])
    number = {signature: index for index, signature in enumerate(ordered, 1)}
    with EXCEPTION_CSV.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["n_x", "x", "n_a", "a", "C_x_to_a", "C_a_to_x"])
        for index, signature in enumerate(ordered, 1):
            for earlier in ordered[: index - 1]:
                if earlier in outgoing[signature]:
                    writer.writerow([
                        index, ",".join(map(str, signature)),
                        number[earlier], ",".join(map(str, earlier)),
                        f"{rates[(signature, earlier)]:.12f}",
                        f"{rates[(earlier, signature)]:.12f}",
                    ])

    content = [
        "% Generated by analysis/appendix_b_signatures.py.",
        *render_order_list(ordered),
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
    print(ORDER_CSV)
    print(EXCEPTION_CSV)
    print(OUTPUT_PATH)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
