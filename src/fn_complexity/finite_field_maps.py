"""Affine degeneration posets for quadratic maps over finite fields."""

from __future__ import annotations

from dataclasses import dataclass
from html import escape
from pathlib import Path


@dataclass(frozen=True)
class QuadraticMapClass:
    """One affine-equivalence class of maps ``F_q^2 -> F_q``."""

    key: str
    name: str
    representative: str
    size: int
    color: str


# An arrow points from a resource class to a class obtainable from it by a
# possibly singular affine input processor and an affine output processor.
ODD_QUADRATIC_MAP_COVERS: tuple[tuple[str, str], ...] = (
    ("parabolic", "linear"),
    ("parabolic", "rank1"),
    ("split", "linear"),
    ("split", "rank1"),
    ("anisotropic", "rank1"),
    ("linear", "constant"),
    ("rank1", "constant"),
)

EVEN_QUADRATIC_MAP_COVERS: tuple[tuple[str, str], ...] = (
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
)

BINARY_QUADRATIC_MAP_COVERS: tuple[tuple[str, str], ...] = (
    ("singleton", "linear"),
    ("linear", "constant"),
)


def _is_prime_power(value: int) -> bool:
    if value < 2:
        return False
    prime = None
    candidate = 2
    while candidate * candidate <= value:
        if value % candidate == 0:
            prime = candidate
            break
        candidate += 1
    if prime is None:
        return True
    remainder = value
    while remainder % prime == 0:
        remainder //= prime
    return remainder == 1


def validate_field_order(q: int) -> None:
    """Require a prime power, hence the order of a finite field."""

    if not _is_prime_power(q):
        raise ValueError(f"q must be a prime power at least 2, got {q}")


def _anisotropic_representative(q: int) -> str:
    if q == 3:
        return "x² + y²"
    if q == 5:
        return "x² − 2y²"
    return "x² − νy²"


def _field_label(q: int) -> str:
    subscript_digits = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
    return f"F{str(q).translate(subscript_digits)}"


def quadratic_map_count(q: int) -> int:
    """Return the number of degree-at-most-two polynomial functions."""

    validate_field_order(q)
    # Over F_2, x^2=x and y^2=y as functions, leaving four Boolean
    # coefficients: xy, x, y, and 1.  For q>2 the six coefficients are unique.
    return 16 if q == 2 else q**6


def quadratic_map_classes(q: int) -> tuple[QuadraticMapClass, ...]:
    """Return affine classes of quadratic polynomial maps ``F_q^2 -> F_q``.

    Input and output equivalence uses invertible affine processors.  For q>2,
    sizes count coefficient vectors ``(a,b,c,d,e,f)``.  For q=2 they count
    distinct polynomial functions, after using ``x^2=x`` and ``y^2=y``.
    """

    validate_field_order(q)
    if q == 2:
        classes = (
            QuadraticMapClass("constant", "constant", "constant", 2, "#64748b"),
            QuadraticMapClass("linear", "balanced", "x", 6, "#0f766e"),
            QuadraticMapClass(
                "singleton", "weight 1 or 3", "xy", 8, "#7c3aed"
            ),
        )
    elif q % 2 == 0:
        classes = (
            QuadraticMapClass("constant", "constant", "constant", q, "#64748b"),
            QuadraticMapClass(
                "linear", "linear", "x", q * (q * q - 1), "#0f766e"
            ),
            QuadraticMapClass(
                "rank1",
                "pure square",
                "x²",
                q * (q * q - 1),
                "#2563eb",
            ),
            QuadraticMapClass(
                "separable",
                "square + aligned linear",
                "x² + x",
                q * (q - 1) * (q * q - 1),
                "#0891b2",
            ),
            QuadraticMapClass(
                "parabolic",
                "parabolic",
                "x² + y",
                q * q * (q - 1) * (q * q - 1),
                "#c2410c",
            ),
            QuadraticMapClass(
                "split",
                "split quadratic",
                "xy",
                q**4 * (q * q - 1) // 2,
                "#7c3aed",
            ),
            QuadraticMapClass(
                "anisotropic",
                "anisotropic · Tr(δ)=1",
                "x² + xy + δy²",
                q**4 * (q - 1) ** 2 // 2,
                "#be123c",
            ),
        )
    else:
        classes = (
            QuadraticMapClass("constant", "constant", "constant", q, "#64748b"),
            QuadraticMapClass(
                "linear", "linear", "x", q * (q * q - 1), "#0f766e"
            ),
            QuadraticMapClass(
                "rank1",
                "rank-1 quadratic",
                "x²",
                q * q * (q * q - 1),
                "#2563eb",
            ),
            QuadraticMapClass(
                "parabolic",
                "parabolic",
                "x² + y",
                q * q * (q - 1) * (q * q - 1),
                "#c2410c",
            ),
            QuadraticMapClass(
                "split",
                "split quadratic",
                "xy",
                q**4 * (q * q - 1) // 2,
                "#7c3aed",
            ),
            QuadraticMapClass(
                "anisotropic",
                "anisotropic",
                _anisotropic_representative(q),
                q**4 * (q - 1) ** 2 // 2,
                "#be123c",
            ),
        )
    if sum(item.size for item in classes) != quadratic_map_count(q):
        raise AssertionError("quadratic-map class sizes must partition all maps")
    return classes


def class_sizes(q: int) -> dict[str, int]:
    """Return ``class key -> orbit size`` for the quadratic-map poset."""

    return {item.key: item.size for item in quadratic_map_classes(q)}


def quadratic_map_covers(q: int) -> tuple[tuple[str, str], ...]:
    """Return the Hasse covers for the field of order ``q``."""

    validate_field_order(q)
    if q == 2:
        return BINARY_QUADRATIC_MAP_COVERS
    if q % 2 == 0:
        return EVEN_QUADRATIC_MAP_COVERS
    return ODD_QUADRATIC_MAP_COVERS


def _diagram_layout(q: int) -> tuple[dict[str, tuple[int, int]], tuple[str, ...]]:
    if q == 2:
        return (
            {
                "singleton": (600, 230),
                "linear": (600, 485),
                "constant": (600, 735),
            },
            ("singleton", "linear", "constant"),
        )
    if q % 2 == 0:
        return (
            {
                "parabolic": (180, 225),
                "split": (600, 225),
                "anisotropic": (1020, 225),
                "linear": (180, 500),
                "rank1": (600, 500),
                "separable": (1020, 500),
                "constant": (600, 750),
            },
            (
                "parabolic",
                "split",
                "anisotropic",
                "linear",
                "rank1",
                "separable",
                "constant",
            ),
        )
    return (
        {
            "parabolic": (190, 235),
            "split": (600, 235),
            "anisotropic": (1010, 235),
            "linear": (385, 500),
            "rank1": (815, 500),
            "constant": (600, 735),
        },
        ("parabolic", "split", "anisotropic", "linear", "rank1", "constant"),
    )


def render_quadratic_map_poset_svg(q: int, output: Path) -> Path:
    """Write a dependency-free SVG Hasse diagram for one prime-power order."""

    classes = {item.key: item for item in quadratic_map_classes(q)}
    covers = quadratic_map_covers(q)
    positions, node_order = _diagram_layout(q)
    width, height = 1200, 875
    node_width, node_height = 256, 112

    def node_markup(key: str) -> str:
        item = classes[key]
        x, y = positions[key]
        return (
            f'<g class="node {escape(key)}" transform="translate({x} {y})">'
            f'<rect x="{-node_width / 2:g}" y="{-node_height / 2:g}" '
            f'width="{node_width}" height="{node_height}" rx="18" '
            f'fill="{item.color}" fill-opacity="0.10" stroke="{item.color}" '
            'stroke-width="2.5"/>'
            f'<text class="form" x="0" y="-20">{escape(item.representative)}</text>'
            f'<text class="kind" x="0" y="8">{escape(item.name)}</text>'
            f'<text class="size" x="0" y="37">|class| = {item.size:,}</text>'
            "</g>"
        )

    outgoing = {
        source: sorted(
            (target for edge_source, target in covers if edge_source == source),
            key=lambda key: positions[key][0],
        )
        for source in positions
    }
    incoming = {
        target: sorted(
            (source for source, edge_target in covers if edge_target == target),
            key=lambda key: positions[key][0],
        )
        for target in positions
    }

    def port_offset(key: str, other: str, adjacency: dict[str, list[str]]) -> float:
        neighbours = adjacency[key]
        if len(neighbours) == 1:
            return 0
        return -48 + 96 * neighbours.index(other) / (len(neighbours) - 1)

    def edge_markup(source: str, target: str) -> str:
        source_x, source_y = positions[source]
        target_x, target_y = positions[target]
        start_x = source_x + port_offset(source, target, outgoing)
        end_x = target_x + port_offset(target, source, incoming)
        start_y = source_y + node_height / 2
        end_y = target_y - node_height / 2 - 10
        middle_y = (start_y + end_y) / 2
        return (
            f'<path class="edge" d="M {start_x:g} {start_y:g} '
            f'C {start_x:g} {middle_y:g}, {end_x:g} {middle_y:g}, '
            f'{end_x:g} {end_y:g}" marker-end="url(#arrow)"/>'
        )

    edges = "".join(edge_markup(*cover) for cover in covers)
    nodes = "".join(node_markup(key) for key in node_order)
    field = _field_label(q)
    map_count = quadratic_map_count(q)
    count_label = "polynomial functions" if q == 2 else "coefficient vectors"
    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}"
     viewBox="0 0 {width} {height}" role="img"
     aria-labelledby="title description">
  <title id="title">Quadratic-map degeneration poset over {field}</title>
  <desc id="description">{len(classes)} affine-equivalence classes. Arrows point downward
  from a quadratic map to a class obtainable using possibly singular affine
  input and output processors. Every node displays its equivalence-class size.</desc>
  <defs>
    <marker id="arrow" viewBox="0 0 10 10" refX="8" refY="5"
            markerWidth="8" markerHeight="8" orient="auto-start-reverse">
      <path d="M 0 0 L 10 5 L 0 10 z" fill="#64748b"/>
    </marker>
    <style>
      .edge {{ fill: none; stroke: #64748b; stroke-width: 2.2; opacity: .74; }}
      .node text {{ fill: #172033; text-anchor: middle;
                    font-family: Inter, ui-sans-serif, system-ui, sans-serif; }}
      .form {{ font-size: 24px; font-weight: 700; }}
      .kind {{ font-size: 16px; fill: #475569; }}
      .size {{ font-size: 18px; font-weight: 650; font-variant-numeric: tabular-nums; }}
      .title {{ fill: #172033; font: 700 32px Inter, ui-sans-serif, system-ui, sans-serif; }}
      .subtitle {{ fill: #475569; font: 17px Inter, ui-sans-serif, system-ui, sans-serif; }}
      .legend {{ fill: #475569; font: 15px Inter, ui-sans-serif, system-ui, sans-serif; }}
    </style>
  </defs>
  <rect width="{width}" height="{height}" fill="#ffffff"/>
  <text class="title" x="600" y="48" text-anchor="middle">
    Quadratic maps {field}² → {field}
  </text>
  <text class="subtitle" x="600" y="78" text-anchor="middle">
    affine-equivalence classes and singular-processor degenerations
  </text>
  <text class="legend" x="600" y="111" text-anchor="middle">
    arrow: resource → implementable degeneration · {map_count:,} {count_label}
  </text>
  {edges}
  {nodes}
</svg>
'''
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(svg, encoding="utf-8")
    return output
