#!/usr/bin/env python3
"""Hand-authored SVG figures for the finite-field exchange section.

``quadratic-map-exchange-order-q{q}.svg``
    Two panels in the house style of the poset gallery.  Left: the degeneration
    poset, with its antichains.  Right: the total order the exchange rate puts
    on the fiber signatures.  Dotted connectors carry each class across, and
    two of them land on the same node, which is the picture of the fact that
    the signature merges the linear and parabolic classes.

``quadratic-map-gibbs-regions-q{q}.svg``
    The Gibbs energy--entropy regions of the four non-degenerate classes, with
    the homothety that realises one exchange rate.

No plotting dependency: both files are emitted directly, matching the palette,
typography and card geometry of ``render_quadratic_map_poset_svg``.
"""

from __future__ import annotations

import math
from pathlib import Path
import sys
from xml.sax.saxutils import escape

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from fn_complexity import (  # noqa: E402
    exchange_rate,
    exchange_rate_result,
    gibbs_point,
)
from fn_complexity.finite_field_maps import (  # noqa: E402
    quadratic_map_classes,
    quadratic_map_covers,
)

OUTPUT_DIRECTORY = PROJECT_ROOT / "paper_finite_fields_maps" / "images"
FIELD_SIZES = (2, 4, 8, 9, 16, 25)

SIGNATURE_COLOR = {"S": "#2563eb", "L": "#0f766e", "A": "#be123c", "X": "#7c3aed"}


def signature_of_class(q: int) -> dict[str, str]:
    """Which fiber signature each poset class carries.

    Several classes share one, and that is the point of the badges.  For odd
    ``q`` the linear and parabolic classes both give ``{q,...,q}``; for even
    ``q`` the pure square joins them, because Frobenius is a bijection, and the
    role of the non-flat small class passes to the Artin--Schreier map
    ``x^2 + x``.
    """

    if q == 2:
        return {"linear": "L", "singleton": "X"}
    if q % 2:
        return {
            "linear": "L",
            "parabolic": "L",
            "rank1": "S",
            "split": "X",
            "anisotropic": "A",
        }
    return {
        "linear": "L",
        "parabolic": "L",
        "rank1": "L",
        "separable": "S",
        "split": "X",
        "anisotropic": "A",
    }


def signature_names(q: int) -> dict[str, str]:
    if q == 2:
        return {"L": "balanced", "X": "weight 1 or 3"}
    if q % 2:
        return {
            "S": "pure square x²",
            "L": "linear · parabolic",
            "A": "anisotropic",
            "X": "split",
        }
    return {
        "S": "square + aligned linear",
        "L": "linear · square · parabolic",
        "A": "anisotropic",
        "X": "split",
    }
EXCHANGE_COLOR = "#b45309"
MUTED = "#64748b"
INK = "#172033"

def poset_positions(q: int) -> dict[str, tuple[int, int]]:
    if q == 2:
        return {
            "singleton": (450, 330),
            "linear": (450, 590),
            "constant": (450, 830),
        }
    if q % 2:
        return {
            "parabolic": (170, 330),
            "split": (450, 330),
            "anisotropic": (730, 330),
            "linear": (310, 590),
            "rank1": (590, 590),
            "constant": (450, 830),
        }
    return {
        "parabolic": (170, 330),
        "split": (450, 330),
        "anisotropic": (730, 330),
        "linear": (170, 590),
        "rank1": (450, 590),
        "separable": (730, 590),
        "constant": (450, 830),
    }
LADDER_X = 1155


def ladder_order(q: int) -> list[str]:
    """Signature classes, most valuable first.

    Over ``F_2`` the identities ``x^2 = x`` and ``y^2 = y`` hold as functions,
    so only three classes survive and only two of them are non-degenerate.
    """

    return ["X", "L"] if q == 2 else ["X", "A", "L", "S"]


def ladder_y(q: int) -> dict[str, int]:
    if q == 2:
        return {"X": 400, "L": 700}
    return {"X": 320, "A": 490, "L": 660, "S": 830}


def comparison_pair(q: int) -> tuple[str, str]:
    """The (source, target) pair whose containment the Gibbs figure draws."""

    return ("L", "X") if q == 2 else ("S", "X")

CARD_WIDTH, CARD_HEIGHT = 244, 116
LADDER_WIDTH, LADDER_HEIGHT = 300, 96
WIDTH, HEIGHT = 1420, 1010


def signatures(q: int) -> dict[str, tuple[int, ...]]:
    if q == 2:
        # x has fibers {2,2}; xy has {3,1}.  The formulas below would collapse
        # A onto X and S onto the constant class, so q = 2 is handled apart.
        return {"L": (2, 2), "X": (3, 1)}
    if q % 2:
        pure_square = tuple(sorted([2 * q] * ((q - 1) // 2) + [q], reverse=True))
    else:
        pure_square = tuple([2 * q] * (q // 2))
    return {
        "S": pure_square,
        "L": tuple([q] * q),
        "X": tuple(sorted([2 * q - 1] + [q - 1] * (q - 1), reverse=True)),
        "A": tuple(sorted([q + 1] * (q - 1) + [1], reverse=True)),
    }


def signature_label(signature: tuple[int, ...]) -> str:
    """A short {a,b,…} rendering with repetition counts."""

    runs: list[tuple[int, int]] = []
    for value in signature:
        if runs and runs[-1][0] == value:
            runs[-1] = (value, runs[-1][1] + 1)
        else:
            runs.append((value, 1))
    parts = [f"{value}" if count == 1 else f"{value}^{count}" for value, count in runs]
    return "{" + ", ".join(parts) + "}"


def field_label(q: int) -> str:
    subscript = str(q).translate(str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉"))
    return f"F{subscript}"


def poset_card(
    item,
    position: tuple[int, int],
    signature_key: str | None,
    signature: tuple[int, ...] | None,
) -> str:
    x, y = position
    half_width, half_height = CARD_WIDTH / 2, CARD_HEIGHT / 2
    signature_line = (
        f'<text class="sig" x="0" y="40">{escape(signature_label(signature))}</text>'
        if signature is not None
        else '<text class="sig" x="0" y="40">degenerate resource</text>'
    )
    if signature_key is None:
        badge = ""
    else:
        color = SIGNATURE_COLOR[signature_key]
        badge = (
            f'<circle cx="{half_width - 26:g}" cy="{-half_height + 26:g}" r="17" '
            f'fill="{color}" fill-opacity="0.9"/>'
            f'<text class="badge" x="{half_width - 26:g}" y="{-half_height + 32:g}">'
            f"{signature_key}</text>"
        )
    return (
        f'<g class="node" transform="translate({x} {y})">'
        f'<rect x="{-half_width:g}" y="{-half_height:g}" width="{CARD_WIDTH}" '
        f'height="{CARD_HEIGHT}" rx="18" fill="{item.color}" fill-opacity="0.10" '
        f'stroke="{item.color}" stroke-width="2.5"/>'
        f'<text class="form" x="0" y="-22">{escape(item.representative)}</text>'
        f'<text class="kind" x="0" y="4">{escape(item.name)}</text>'
        f"{signature_line}{badge}"
        "</g>"
    )


def ladder_card(q: int, key: str, signature: tuple[int, ...]) -> str:
    x, y = LADDER_X, ladder_y(q)[key]
    half_width, half_height = LADDER_WIDTH / 2, LADDER_HEIGHT / 2
    color = SIGNATURE_COLOR[key]
    return (
        f'<g class="node" transform="translate({x} {y})">'
        f'<rect x="{-half_width:g}" y="{-half_height:g}" width="{LADDER_WIDTH}" '
        f'height="{LADDER_HEIGHT}" rx="16" fill="{color}" fill-opacity="0.12" '
        f'stroke="{color}" stroke-width="2.5"/>'
        f'<text class="rank" x="{-half_width + 26:g}" y="6" '
        f'fill="{color}">{key}</text>'
        f'<text class="ladder-kind" x="18" y="-14">{escape(signature_names(q)[key])}</text>'
        f'<text class="sig" x="18" y="14">{escape(signature_label(signature))}</text>'
        "</g>"
    )


def bezier_edge(start: tuple[float, float], end: tuple[float, float], css: str) -> str:
    (start_x, start_y), (end_x, end_y) = start, end
    middle_y = (start_y + end_y) / 2
    return (
        f'<path class="{css}" d="M {start_x:g} {start_y:g} '
        f"C {start_x:g} {middle_y:g}, {end_x:g} {middle_y:g}, "
        f'{end_x:g} {end_y:g}" marker-end="url(#arrow)"/>'
    )


def render_exchange_order(q: int) -> Path:
    classes = {item.key: item for item in quadratic_map_classes(q)}
    covers = quadratic_map_covers(q)
    positions = poset_positions(q)
    class_signature = signature_of_class(q)
    sig = signatures(q)

    pieces: list[str] = []

    # Panel separator.
    pieces.append(
        f'<line class="divider" x1="900" y1="180" x2="900" y2="{HEIGHT - 120}"/>'
    )
    pieces.append(
        '<text class="panel" x="450" y="212" text-anchor="middle">'
        "classical degeneration order</text>"
    )
    pieces.append(
        f'<text class="panel" x="{LADDER_X}" y="212" text-anchor="middle">'
        "exchange comparison</text>"
    )
    pieces.append(
        '<text class="panel-note" x="450" y="240" text-anchor="middle">'
        + (
            f"{len(classes)} classes forming a chain — already a total order"
            if q == 2
            else f"{len(classes)} classes, three of them maximal and pairwise "
            "incomparable"
        )
        + "</text>"
    )
    pieces.append(
        f'<text class="panel-note" x="{LADDER_X}" y="240" text-anchor="middle">'
        + (
            "the same order, now with a price on every step"
            if q == 2
            else "a strict total order on the fiber signatures"
        )
        + "</text>"
    )

    # Hasse covers.
    for source, target in covers:
        source_x, source_y = positions[source]
        target_x, target_y = positions[target]
        pieces.append(
            bezier_edge(
                (source_x, source_y + CARD_HEIGHT / 2),
                (target_x, target_y - CARD_HEIGHT / 2 - 10),
                "edge",
            )
        )

    # The ladder itself, with both directed rates on every step.
    order = ladder_order(q)
    heights = ladder_y(q)
    for index in range(len(order) - 1):
        upper, lower = order[index], order[index + 1]
        forward = exchange_rate(sig[lower], sig[upper])
        backward = exchange_rate(sig[upper], sig[lower])
        y_start = heights[lower] - LADDER_HEIGHT / 2
        y_end = heights[upper] + LADDER_HEIGHT / 2 + 10
        pieces.append(
            f'<path class="exchange" d="M {LADDER_X:g} {y_start:g} '
            f'L {LADDER_X:g} {y_end:g}" marker-end="url(#arrow-exchange)"/>'
        )
        middle_y = (y_start + y_end) / 2
        pieces.append(
            f'<text class="rate" x="{LADDER_X + 22:g}" y="{middle_y - 6:g}">'
            f"C({lower}→{upper}) = {forward:.4f}</text>"
        )
        pieces.append(
            f'<text class="rate" x="{LADDER_X + 22:g}" y="{middle_y + 14:g}">'
            f"C({upper}→{lower}) = {backward:.4f}</text>"
        )
        pieces.append(
            f'<text class="verdict" x="{LADDER_X - 26:g}" y="{middle_y + 5:g}" '
            f'text-anchor="end">{lower} ≺ {upper}</text>'
        )

    for key in positions:
        signature_key = class_signature.get(key)
        pieces.append(
            poset_card(
                classes[key],
                positions[key],
                signature_key,
                sig[signature_key] if signature_key else None,
            )
        )
    for key in order:
        pieces.append(ladder_card(q, key, sig[key]))

    field = field_label(q)
    if q == 2:
        headline = "The smallest case: the exchange rate agrees with degeneration"
        strapline = (
            f"quadratic maps {field}² → {field} · the two non-constant classes carry "
            "the signatures {2,2} and {3,1}, the opening example of the companion paper"
        )
    else:
        headline = "What the exchange rate decides that degeneration does not"
        strapline = (
            f"quadratic maps {field}² → {field} · fiber signatures are point counts, "
            "so every rate below is a ratio of logarithms of Weil numbers"
        )
    merged_count = {1: "one", 2: "two", 3: "three", 4: "four"}[
        sum(1 for value in class_signature.values() if value == "L")
    ]
    legend_y = HEIGHT - 46
    legend = (
        f'<line class="edge" x1="60" y1="{legend_y - 5}" x2="104" y2="{legend_y - 5}" '
        'marker-end="url(#arrow)"/>'
        f'<text class="legend" x="116" y="{legend_y}">degeneration cover '
        "(classical)</text>"
        f'<circle cx="392" cy="{legend_y - 5}" r="14" fill="{SIGNATURE_COLOR["L"]}" '
        'fill-opacity="0.9"/>'
        f'<text class="badge" x="392" y="{legend_y}">L</text>'
        f'<text class="legend" x="416" y="{legend_y}">fiber signature of the class'
        + ("" if q == 2 else f" · {merged_count} cards carry L")
        + "</text>"
        f'<line class="exchange" x1="836" y1="{legend_y - 5}" x2="880" y2="{legend_y - 5}" '
        'marker-end="url(#arrow-exchange)"/>'
        f'<text class="legend" x="892" y="{legend_y}">strictly the better resource, '
        "decided by the rate and not by the poset</text>"
    )

    body = "\n  ".join(pieces)
    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH}" height="{HEIGHT}"
     viewBox="0 0 {WIDTH} {HEIGHT}" role="img"
     aria-labelledby="title description">
  <title id="title">Degeneration order and exchange comparison over {field}</title>
  <desc id="description">The degeneration poset of quadratic maps over {field} has
  three pairwise incomparable maximal classes. The asymptotic exchange rate merges
  the linear and parabolic classes, which share a fiber signature, and totally
  orders the four remaining signatures.</desc>
  <defs>
    <marker id="arrow" viewBox="0 0 10 10" refX="8" refY="5"
            markerWidth="8" markerHeight="8" orient="auto-start-reverse">
      <path d="M 0 0 L 10 5 L 0 10 z" fill="{MUTED}"/>
    </marker>
    <marker id="arrow-exchange" viewBox="0 0 10 10" refX="8" refY="5"
            markerWidth="9" markerHeight="9" orient="auto-start-reverse">
      <path d="M 0 0 L 10 5 L 0 10 z" fill="{EXCHANGE_COLOR}"/>
    </marker>
    <style>
      .edge {{ fill: none; stroke: {MUTED}; stroke-width: 2.2; opacity: .74; }}
      .carry {{ fill: none; stroke-width: 1.6; stroke-dasharray: 3 5; opacity: .55; }}
      .exchange {{ fill: none; stroke: {EXCHANGE_COLOR}; stroke-width: 2.6; }}
      .divider {{ stroke: #e2e8f0; stroke-width: 2; }}
      .node text {{ fill: {INK}; text-anchor: middle;
                    font-family: Inter, ui-sans-serif, system-ui, sans-serif; }}
      .form {{ font-size: 23px; font-weight: 700; }}
      .kind {{ font-size: 15px; fill: #475569; }}
      .sig {{ font-size: 14px; fill: #475569;
              font-family: "SF Mono", ui-monospace, Menlo, monospace; }}
      .badge {{ fill: #ffffff; font: 700 17px Inter, ui-sans-serif, system-ui, sans-serif;
                text-anchor: middle; }}
      .rank {{ font-size: 26px; font-weight: 700; text-anchor: middle; }}
      .ladder-kind {{ font-size: 15px; fill: #475569; }}
      .rate {{ fill: {EXCHANGE_COLOR}; font: 13px "SF Mono", ui-monospace, Menlo, monospace; }}
      .verdict {{ fill: {EXCHANGE_COLOR}; font: 700 15px Inter, ui-sans-serif, system-ui, sans-serif; }}
      .title {{ fill: {INK}; font: 700 32px Inter, ui-sans-serif, system-ui, sans-serif; }}
      .subtitle {{ fill: #475569; font: 17px Inter, ui-sans-serif, system-ui, sans-serif; }}
      .panel {{ fill: {INK}; font: 650 19px Inter, ui-sans-serif, system-ui, sans-serif; }}
      .panel-note {{ fill: #64748b; font: 14px Inter, ui-sans-serif, system-ui, sans-serif; }}
      .legend {{ fill: #475569; font: 14px Inter, ui-sans-serif, system-ui, sans-serif; }}
    </style>
  </defs>
  <rect width="{WIDTH}" height="{HEIGHT}" fill="#ffffff"/>
  <text class="title" x="{WIDTH / 2:g}" y="62" text-anchor="middle">
    {headline}
  </text>
  <text class="subtitle" x="{WIDTH / 2:g}" y="94" text-anchor="middle">
    {strapline}
  </text>
  {body}
  {legend}
</svg>
'''
    output = OUTPUT_DIRECTORY / f"quadratic-map-exchange-order-q{q}.svg"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(svg, encoding="utf-8")
    return output


def gibbs_polyline(signature: tuple[int, ...]) -> list[tuple[float, float]]:
    temperatures = (
        [0.0]
        + [10.0 ** (-3.0 + 6.0 * index / 399.0) for index in range(400)]
        + [math.inf]
    )
    points = [gibbs_point(signature, temperature) for temperature in temperatures]
    return [(point.energy, point.entropy) for point in points]


def render_gibbs_regions(q: int) -> Path:
    sig = signatures(q)
    curves = {key: gibbs_polyline(sig[key]) for key in ladder_order(q)}
    source, target = comparison_pair(q)
    rate = exchange_rate(sig[source], sig[target])
    scaled = [(x / rate, y / rate) for x, y in curves[source]]

    width, height = 1200, 820
    left, right, top, bottom = 130, width - 330, 150, height - 90
    energies = [x for curve in curves.values() for x, _ in curve] + [x for x, _ in scaled]
    entropies = [y for curve in curves.values() for _, y in curve] + [y for _, y in scaled]
    min_energy, max_energy = min(energies) * 1.06, 0.0
    max_entropy = max(entropies) * 1.08
    span_energy = max_energy - min_energy

    def project(energy: float, entropy: float) -> tuple[float, float]:
        px = left + (energy - min_energy) / span_energy * (right - left)
        py = bottom - entropy / max_entropy * (bottom - top)
        return px, py

    pieces: list[str] = []
    for index in range(6):
        value = min_energy + span_energy * index / 5
        px, _ = project(value, 0)
        pieces.append(f'<line class="grid" x1="{px:g}" y1="{top}" x2="{px:g}" y2="{bottom}"/>')
        pieces.append(
            f'<text class="tick" x="{px:g}" y="{bottom + 26}" text-anchor="middle">'
            f"{value:.1f}</text>"
        )
    for index in range(5):
        value = max_entropy * index / 4
        _, py = project(min_energy, value)
        pieces.append(f'<line class="grid" x1="{left}" y1="{py:g}" x2="{right}" y2="{py:g}"/>')
        pieces.append(
            f'<text class="tick" x="{left - 14}" y="{py + 5:g}" text-anchor="end">'
            f"{value:.1f}</text>"
        )

    def region_points(curve: list[tuple[float, float]]) -> str:
        """Boundary of the closed region: the Gibbs curve, then E = 0, then H = 0."""

        projected = [project(x, y) for x, y in curve]
        right_x, _ = project(0.0, 0.0)
        _, zero_y = project(0.0, 0.0)
        first_px, _ = projected[0]
        _, last_py = projected[-1]
        closure = [(right_x, last_py), (right_x, zero_y), (first_px, zero_y)]
        return " ".join(f"{px:g},{py:g}" for px, py in projected + closure)

    for key in (target, source):
        color = SIGNATURE_COLOR[key]
        pieces.append(
            f'<polygon class="fill" fill="{color}" points="{region_points(curves[key])}"/>'
        )
        pieces.append(
            f'<polygon class="outline" stroke="{color}" fill="none" '
            f'points="{region_points(curves[key])}"/>'
        )
        pieces.append(
            f'<polyline class="curve" stroke="{color}" points="'
            + " ".join(f"{px:g},{py:g}" for px, py in (project(x, y) for x, y in curves[key]))
            + '"/>'
        )

    pieces.append(
        f'<polygon class="curve scaled" fill="none" stroke="{SIGNATURE_COLOR[source]}" '
        f'points="{region_points(scaled)}"/>'
    )
    for key in (target, source):
        for energy, entropy in (curves[key][0], curves[key][-1]):
            px, py = project(energy, entropy)
            pieces.append(
                f'<circle cx="{px:g}" cy="{py:g}" r="4.5" fill="#ffffff" '
                f'stroke="{SIGNATURE_COLOR[key]}" stroke-width="2.2"/>'
            )

    legend_x = right + 48
    for index, key in enumerate((target, source)):
        y = top + 18 + index * 34
        pieces.append(
            f'<rect x="{legend_x}" y="{y - 12}" width="26" height="14" rx="4" '
            f'fill="{SIGNATURE_COLOR[key]}" fill-opacity="0.25" '
            f'stroke="{SIGNATURE_COLOR[key]}" stroke-width="2"/>'
        )
        pieces.append(
            f'<text class="legend" x="{legend_x + 36}" y="{y}">'
            f"{key} · {escape(signature_names(q)[key])}</text>"
        )
    y = top + 18 + 2 * 34 + 10
    pieces.append(
        f'<line class="curve scaled" stroke="{SIGNATURE_COLOR[source]}" '
        f'x1="{legend_x}" y1="{y - 5}" x2="{legend_x + 26}" y2="{y - 5}"/>'
    )
    pieces.append(
        f'<text class="legend" x="{legend_x + 36}" y="{y}">{source} dilated by</text>'
    )
    pieces.append(
        f'<text class="legend" x="{legend_x + 36}" y="{y + 20}">'
        f"1/C({source}→{target}) = {1 / rate:.4f},</text>"
    )
    pieces.append(
        f'<text class="legend" x="{legend_x + 36}" y="{y + 40}">'
        "the smallest copy</text>"
    )
    pieces.append(
        f'<text class="legend" x="{legend_x + 36}" y="{y + 60}">'
        f"containing {target}</text>"
    )

    # Where the two boundaries touch depends on which endpoint attains the
    # infimum: beta = 0 compares fiber counts (equal heights), beta = infinity
    # compares largest fibers (equal widths).
    contact = exchange_rate_result(implemented=sig[target], implementer=sig[source])
    if contact.beta == 0.0:
        contact_entropy = max(entropy for _, entropy in curves[target])
        contact_x, contact_y = project(0.0, contact_entropy)
        pieces.append(
            f'<line class="contact" x1="{contact_x - 250:g}" y1="{contact_y:g}" '
            f'x2="{contact_x:g}" y2="{contact_y:g}"/>'
        )
        pieces.append(
            f'<text class="note" x="{contact_x - 12:g}" y="{contact_y - 14:g}" '
            'text-anchor="end">contact at T = ∞: equal heights, '
            f'log {len(sig[target])} = {1 / rate:.4f} · log {len(sig[source])}</text>'
        )
    elif math.isinf(contact.beta):
        contact_energy = min(energy for energy, _ in curves[target])
        contact_x, _ = project(contact_energy, 0.0)
        _, base_y = project(contact_energy, 0.0)
        _, top_y = project(contact_energy, max_entropy * 0.62)
        pieces.append(
            f'<line class="contact" x1="{contact_x:g}" y1="{base_y:g}" '
            f'x2="{contact_x:g}" y2="{top_y:g}"/>'
        )
        pieces.append(
            f'<text class="note" x="{contact_x + 12:g}" y="{top_y - 12:g}">'
            "contact at T = 0: equal widths, "
            f"log {max(sig[target])} = {1 / rate:.4f} · log {max(sig[source])}</text>"
        )
    else:
        temperature = 1.0 / contact.beta
        point = gibbs_point(sig[target], temperature)
        contact_x, contact_y = project(point.energy, point.entropy)
        pieces.append(
            f'<circle cx="{contact_x:g}" cy="{contact_y:g}" r="7" fill="none" '
            'stroke="#b45309" stroke-width="3"/>'
        )
        pieces.append(
            f'<text class="note" x="{contact_x + 14:g}" y="{contact_y - 12:g}">'
            f"tangency at T = {temperature:.4f}</text>"
        )

    field = field_label(q)
    names = signature_names(q)
    pair_text = f"{names[source]} ({source}) against {names[target]} ({target})"
    body = "\n  ".join(pieces)
    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}"
     viewBox="0 0 {width} {height}" role="img"
     aria-labelledby="title description">
  <title id="title">Gibbs regions of the quadratic-map classes over {field}</title>
  <desc id="description">Energy-entropy regions of the four non-degenerate fiber
  signatures. The dotted curve is the smallest dilation of the anisotropic region
  about the origin that contains the split region; its scale factor is the
  reciprocal of the exchange rate.</desc>
  <defs>
    <style>
      .grid {{ stroke: #e2e8f0; stroke-width: 1; }}
      .axis {{ stroke: #94a3b8; stroke-width: 1.8; }}
      .curve {{ fill: none; stroke-width: 2.6; stroke-linejoin: round; }}
      .scaled {{ stroke-width: 2; stroke-dasharray: 6 5; }}
      .fill {{ fill-opacity: 0.11; }}
      .outline {{ stroke-width: 1.4; opacity: .45; }}
      .tick {{ fill: #64748b; font: 13px Inter, ui-sans-serif, system-ui, sans-serif; }}
      .axis-label {{ fill: #475569; font: 15px Inter, ui-sans-serif, system-ui, sans-serif; }}
      .legend {{ fill: #475569; font: 14px Inter, ui-sans-serif, system-ui, sans-serif; }}
      .contact {{ stroke: #b45309; stroke-width: 3; opacity: .8; }}
      .note {{ fill: #b45309; font: 600 13px Inter, ui-sans-serif, system-ui, sans-serif; }}
      .title {{ fill: {INK}; font: 700 30px Inter, ui-sans-serif, system-ui, sans-serif; }}
      .subtitle {{ fill: #475569; font: 16px Inter, ui-sans-serif, system-ui, sans-serif; }}
    </style>
  </defs>
  <rect width="{width}" height="{height}" fill="#ffffff"/>
  <text class="title" x="{width / 2:g}" y="58" text-anchor="middle">
    An exchange rate as a dilation of energy–entropy regions
  </text>
  <text class="subtitle" x="{width / 2:g}" y="90" text-anchor="middle">
    {pair_text} over {field} · the rate is the largest factor by which one region fits inside the other
  </text>
  {body}
  <line class="axis" x1="{left}" y1="{bottom}" x2="{right}" y2="{bottom}"/>
  <line class="axis" x1="{left}" y1="{top}" x2="{left}" y2="{bottom}"/>
  <text class="axis-label" x="{(left + right) / 2:g}" y="{bottom + 58}"
        text-anchor="middle">mean energy E</text>
  <text class="axis-label" x="{left - 76}" y="{(top + bottom) / 2:g}"
        text-anchor="middle" transform="rotate(-90 {left - 76} {(top + bottom) / 2:g})">
    entropy H</text>
</svg>
'''
    output = OUTPUT_DIRECTORY / f"quadratic-map-gibbs-regions-q{q}.svg"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(svg, encoding="utf-8")
    return output


def main(argv: list[str]) -> int:
    orders = [int(value) for value in argv] or list(FIELD_SIZES)
    for q in orders:
        for path in (render_exchange_order(q), render_gibbs_regions(q)):
            print("written", path.relative_to(PROJECT_ROOT))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
