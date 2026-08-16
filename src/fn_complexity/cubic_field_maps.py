"""Cubic-map degeneration posets over the field with three elements."""

from __future__ import annotations

from dataclasses import dataclass
from html import escape
from pathlib import Path
import subprocess


@dataclass(frozen=True)
class CubicMapClass:
    """One class in a cubic-map degeneration poset over ``F_3``."""

    key: str
    name: str
    representative: str
    coefficients: tuple[int, ...]
    size: int
    fiber_label: str
    color: str


# Coefficients use this basis. Since x^3=x and y^3=y as functions on F_3,
# these eight monomials are a basis for degree-at-most-three polynomial
# functions F_3^2 -> F_3. The missing reduced monomial is x^2 y^2.
CUBIC_Q3_BASIS: tuple[tuple[int, int], ...] = (
    (0, 0),
    (1, 0),
    (0, 1),
    (2, 0),
    (1, 1),
    (0, 2),
    (2, 1),
    (1, 2),
)


AFFINE_INPUT_CLASSES: tuple[CubicMapClass, ...] = (
    CubicMapClass(
        "constant",
        "constant",
        "0",
        (0, 0, 0, 0, 0, 0, 0, 0),
        3,
        "fibers 9",
        "#64748b",
    ),
    CubicMapClass(
        "cubic-711",
        "cubic A",
        "xy + xy²",
        (0, 0, 0, 0, 1, 0, 0, 1),
        216,
        "fibers 7+1+1",
        "#9333ea",
    ),
    CubicMapClass(
        "rank1",
        "rank-1 quadratic",
        "x²",
        (0, 0, 0, 1, 0, 0, 0, 0),
        72,
        "fibers 6+3",
        "#2563eb",
    ),
    CubicMapClass(
        "cubic-63",
        "cubic B",
        "y + y² + x²y",
        (0, 0, 1, 0, 0, 1, 1, 0),
        432,
        "fibers 6+3",
        "#9333ea",
    ),
    CubicMapClass(
        "cubic-522-a",
        "cubic C",
        "y² + x²y",
        (0, 0, 0, 0, 0, 1, 1, 0),
        1296,
        "fibers 5+2+2",
        "#9333ea",
    ),
    CubicMapClass(
        "cubic-522-b",
        "cubic D",
        "xy²",
        (0, 0, 0, 0, 0, 0, 0, 1),
        648,
        "fibers 5+2+2",
        "#9333ea",
    ),
    CubicMapClass(
        "cubic-441-a",
        "cubic E",
        "xy + y² + x²y",
        (0, 0, 0, 0, 1, 1, 1, 0),
        1296,
        "fibers 4+4+1",
        "#9333ea",
    ),
    CubicMapClass(
        "split",
        "split quadratic",
        "xy",
        (0, 0, 0, 0, 1, 0, 0, 0),
        324,
        "fibers 5+2+2",
        "#7c3aed",
    ),
    CubicMapClass(
        "cubic-441-b",
        "cubic F",
        "y + xy + xy²",
        (0, 0, 1, 0, 1, 0, 0, 1),
        432,
        "fibers 4+4+1",
        "#9333ea",
    ),
    CubicMapClass(
        "linear",
        "linear",
        "x",
        (0, 1, 0, 0, 0, 0, 0, 0),
        24,
        "fibers 3+3+3",
        "#0f766e",
    ),
    CubicMapClass(
        "cubic-333-a",
        "cubic G",
        "y + x²y",
        (0, 0, 1, 0, 0, 0, 1, 0),
        648,
        "fibers 3+3+3",
        "#9333ea",
    ),
    CubicMapClass(
        "anisotropic",
        "anisotropic quadratic",
        "x² + y²",
        (0, 0, 0, 1, 0, 1, 0, 0),
        162,
        "fibers 4+4+1",
        "#be123c",
    ),
    CubicMapClass(
        "cubic-333-b",
        "cubic H",
        "x + y + y² + x²y",
        (0, 1, 1, 0, 0, 1, 1, 0),
        864,
        "fibers 3+3+3",
        "#9333ea",
    ),
    CubicMapClass(
        "parabolic",
        "parabolic quadratic",
        "x² + y",
        (0, 0, 1, 1, 0, 0, 0, 0),
        144,
        "fibers 3+3+3",
        "#c2410c",
    ),
)


AFFINE_INPUT_COVERS: tuple[tuple[str, str], ...] = (
    ("cubic-711", "linear"),
    ("cubic-711", "rank1"),
    ("cubic-63", "rank1"),
    ("cubic-522-a", "linear"),
    ("cubic-522-a", "rank1"),
    ("cubic-522-b", "linear"),
    ("cubic-522-b", "rank1"),
    ("cubic-441-a", "linear"),
    ("cubic-441-a", "rank1"),
    ("split", "linear"),
    ("split", "rank1"),
    ("cubic-441-b", "linear"),
    ("cubic-441-b", "rank1"),
    ("cubic-333-a", "linear"),
    ("cubic-333-a", "rank1"),
    ("anisotropic", "rank1"),
    ("cubic-333-b", "linear"),
    ("cubic-333-b", "rank1"),
    ("parabolic", "linear"),
    ("parabolic", "rank1"),
    ("linear", "constant"),
    ("rank1", "constant"),
)


QUADRATIC_INPUT_CLASSES: tuple[CubicMapClass, ...] = (
    CubicMapClass(
        "constant",
        "constant functions",
        "0",
        (0, 0, 0, 0, 0, 0, 0, 0),
        3,
        "image size 1",
        "#64748b",
    ),
    CubicMapClass(
        "two-valued",
        "all two-valued cubic functions",
        "x²",
        (0, 0, 0, 1, 0, 0, 0, 0),
        504,
        "fibers 6+3",
        "#2563eb",
    ),
    CubicMapClass(
        "surjective",
        "all surjective cubic functions",
        "x",
        (0, 1, 0, 0, 0, 0, 0, 0),
        6054,
        "image size 3",
        "#9333ea",
    ),
)


QUADRATIC_INPUT_COVERS: tuple[tuple[str, str], ...] = (
    ("surjective", "two-valued"),
    ("two-valued", "constant"),
)


def cubic_q3_map_count() -> int:
    """Return the number of degree-at-most-three polynomial functions."""

    return 3**8


def cubic_q3_classes(processor_case: str) -> tuple[CubicMapClass, ...]:
    """Return classes for ``linear`` or generated ``quadratic`` inputs."""

    if processor_case == "linear":
        return AFFINE_INPUT_CLASSES
    if processor_case == "quadratic":
        return QUADRATIC_INPUT_CLASSES
    raise ValueError("processor_case must be 'linear' or 'quadratic'")


def cubic_q3_covers(processor_case: str) -> tuple[tuple[str, str], ...]:
    """Return Hasse covers for the requested processor case."""

    if processor_case == "linear":
        return AFFINE_INPUT_COVERS
    if processor_case == "quadratic":
        return QUADRATIC_INPUT_COVERS
    raise ValueError("processor_case must be 'linear' or 'quadratic'")


def _node_markup(
    item: CubicMapClass,
    position: tuple[float, float],
    node_width: float,
    node_height: float,
) -> str:
    x, y = position
    return (
        f'<g class="node {escape(item.key)}" transform="translate({x:g} {y:g})">'
        f'<rect x="{-node_width / 2:g}" y="{-node_height / 2:g}" '
        f'width="{node_width:g}" height="{node_height:g}" rx="17" '
        f'fill="{item.color}" fill-opacity="0.10" stroke="{item.color}" '
        'stroke-width="2.5"/>'
        f'<text class="form" x="0" y="-22">{escape(item.representative)}</text>'
        f'<text class="kind" x="0" y="3">{escape(item.name)} · '
        f'{escape(item.fiber_label)}</text>'
        f'<text class="size" x="0" y="33">|class| = {item.size:,}</text>'
        '</g>'
    )


def _edge_markup(
    source: str,
    target: str,
    positions: dict[str, tuple[float, float]],
    node_height: float,
    target_index: int = 0,
    target_count: int = 1,
) -> str:
    source_x, source_y = positions[source]
    target_x, target_y = positions[target]
    start_x = source_x
    if target_count > 1:
        start_x += -32 + 64 * target_index / (target_count - 1)
    end_x = target_x
    start_y = source_y + node_height / 2
    end_y = target_y - node_height / 2 - 10
    middle_y = start_y + 0.54 * (end_y - start_y)
    marker = "arrow-linear" if target == "linear" else "arrow-default"
    edge_class = "edge linear-edge" if target == "linear" else "edge"
    return (
        f'<path class="{edge_class}" d="M {start_x:g} {start_y:g} '
        f'C {start_x:g} {middle_y:g}, {end_x:g} {middle_y:g}, '
        f'{end_x:g} {end_y:g}" marker-end="url(#{marker})"/>'
    )


def render_cubic_q3_poset_svg(processor_case: str, output: Path) -> Path:
    """Render one of the two cubic-map Hasse diagrams as dependency-free SVG."""

    classes = {item.key: item for item in cubic_q3_classes(processor_case)}
    covers = cubic_q3_covers(processor_case)
    if sum(item.size for item in classes.values()) != cubic_q3_map_count():
        raise AssertionError("cubic-map class sizes must partition all maps")

    if processor_case == "linear":
        width, height = 3200, 940
        node_width, node_height = 252, 108
        top_order = (
            "cubic-711",
            "parabolic",
            "split",
            "cubic-522-a",
            "cubic-522-b",
            "cubic-441-a",
            "cubic-441-b",
            "cubic-333-a",
            "cubic-333-b",
            "cubic-63",
            "anisotropic",
        )
        positions = {
            key: (160 + index * 288, 245)
            for index, key in enumerate(top_order)
        }
        positions.update(
            {
                "linear": (1080, 620),
                "rank1": (2120, 620),
                "constant": (1600, 835),
            }
        )
        node_order = top_order + ("linear", "rank1", "constant")
        title = "Cubic maps F₃² → F₃ · affine input processors"
        subtitle = (
            "14 affine-equivalence classes · singular affine processors "
            "define the degeneration order"
        )
        description = (
            "Fourteen classes of degree-at-most-three polynomial functions. "
            "Arrows point from a resource map to an affine degeneration."
        )
    else:
        width, height = 1200, 875
        node_width, node_height = 390, 116
        positions = {
            "surjective": (600, 235),
            "two-valued": (600, 505),
            "constant": (600, 755),
        }
        node_order = ("surjective", "two-valued", "constant")
        title = "Cubic maps F₃² → F₃ · quadratic input processors"
        subtitle = (
            "poset of the generated preorder · linear output processors"
        )
        description = (
            "Three mutual-reachability classes in the transitive closure of "
            "quadratic-input reductions. The classes form a chain."
        )

    outgoing: dict[str, list[str]] = {key: [] for key in positions}
    for source, target in covers:
        outgoing[source].append(target)
    edges = "".join(
        _edge_markup(
            source,
            target,
            positions,
            node_height,
            outgoing[source].index(target),
            len(outgoing[source]),
        )
        for source, target in covers
    )
    nodes = "".join(
        _node_markup(classes[key], positions[key], node_width, node_height)
        for key in node_order
    )
    legend = (
        "arrow: resource → implementable degeneration · "
        f"{cubic_q3_map_count():,} polynomial functions"
    )
    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}"
     viewBox="0 0 {width} {height}" role="img"
     aria-labelledby="title description">
  <title id="title">{escape(title)}</title>
  <desc id="description">{escape(description)}</desc>
  <defs>
    <marker id="arrow-default" viewBox="0 0 10 10" refX="8" refY="5"
            markerWidth="8" markerHeight="8" orient="auto-start-reverse">
      <path d="M 0 0 L 10 5 L 0 10 z" fill="#64748b"/>
    </marker>
    <marker id="arrow-linear" viewBox="0 0 10 10" refX="8" refY="5"
            markerWidth="8" markerHeight="8" orient="auto-start-reverse">
      <path d="M 0 0 L 10 5 L 0 10 z" fill="#0f766e"/>
    </marker>
    <style>
      .edge {{ fill: none; stroke: #64748b; stroke-width: 2.2; opacity: .58; }}
      .linear-edge {{ stroke: #0f766e; }}
      .node text {{ fill: #172033; text-anchor: middle;
                    font-family: Inter, ui-sans-serif, system-ui, sans-serif; }}
      .form {{ font-size: 21px; font-weight: 700; }}
      .kind {{ font-size: 13px; fill: #475569; }}
      .size {{ font-size: 17px; font-weight: 650; font-variant-numeric: tabular-nums; }}
      .title {{ fill: #172033; font: 700 32px Inter, ui-sans-serif, system-ui, sans-serif; }}
      .subtitle {{ fill: #475569; font: 17px Inter, ui-sans-serif, system-ui, sans-serif; }}
      .legend {{ fill: #475569; font: 15px Inter, ui-sans-serif, system-ui, sans-serif; }}
    </style>
  </defs>
  <rect width="{width}" height="{height}" fill="#ffffff"/>
  <text class="title" x="{width / 2:g}" y="48" text-anchor="middle">
    {escape(title)}
  </text>
  <text class="subtitle" x="{width / 2:g}" y="78" text-anchor="middle">
    {escape(subtitle)}
  </text>
  <text class="legend" x="{width / 2:g}" y="108" text-anchor="middle">
    {escape(legend)}
  </text>
  {edges}
  {nodes}
</svg>
'''
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(svg, encoding="utf-8")
    return output


@dataclass(frozen=True)
class CubicQ8GeneratedClass:
    """One node, or one set of order-twin nodes, in the generated poset."""

    key: str
    representative: str
    name: str
    image_size: int
    class_count: int
    class_size: int
    color: str

    @property
    def total_size(self) -> int:
        return self.class_count * self.class_size


# These are the 110 mutual-reachability classes obtained from the 126 affine
# orbits of degree-at-most-three maps F_8^2 -> F_8.  Six sets of pairwise
# incomparable classes have identical strict upper and lower sets; each such
# antichain is drawn as one multiplicity box to keep the Hasse diagram legible.
CUBIC_Q8_QUADRATIC_CLASSES: tuple[CubicQ8GeneratedClass, ...] = (
    CubicQ8GeneratedClass("g0", "0", "constants", 1, 1, 8, "#64748b"),
    CubicQ8GeneratedClass(
        "g1", "x", "surjective SCC · 17 affine orbits", 8, 1, 62_734_392,
        "#0f766e",
    ),
    CubicQ8GeneratedClass("g2", "x²", "pure-square quadratic", 8, 1, 504, "#2563eb"),
    CubicQ8GeneratedClass("g3", "x²+x", "additive quadratic", 4, 1, 3_528, "#0891b2"),
    CubicQ8GeneratedClass(
        "g4", "x²+xy+y²", "anisotropic quadratic", 8, 1, 100_352,
        "#be123c",
    ),
    CubicQ8GeneratedClass("g5", "y³", "homogeneous cube", 8, 1, 4_032, "#c2410c"),
    CubicQ8GeneratedClass("g6", "x+y³", "triple-root cubic", 8, 1, 225_792, "#c2410c"),
    CubicQ8GeneratedClass("g7", "y³+y", "five-valued cubic", 5, 1, 28_224, "#0284c7"),
    CubicQ8GeneratedClass("g8", "x²+y³", "triple-root cubic", 8, 1, 225_792, "#c2410c"),
    CubicQ8GeneratedClass("g9", "x+x²+y³", "triple-root cubic", 8, 1, 790_272, "#c2410c"),
    CubicQ8GeneratedClass("g10", "x+y+x²+y³", "triple-root cubic", 8, 1, 790_272, "#c2410c"),
    CubicQ8GeneratedClass(
        "g11", "x²+xy+y³", "2 order-twin triple-root classes", 8, 2,
        6_322_176, "#c2410c",
    ),
    CubicQ8GeneratedClass("g12", "x²+xy²", "double+simple leading roots", 8, 1, 1_806_336, "#15803d"),
    CubicQ8GeneratedClass(
        "g13", "y+x²+xy²", "41 order twins · leading types 1+8+32", 8,
        41, 12_644_352, "#7c3aed",
    ),
    CubicQ8GeneratedClass(
        "g14", "x+x²+xy²", "15 order-twin double+simple classes", 8, 15,
        6_322_176, "#15803d",
    ),
    CubicQ8GeneratedClass(
        "g15", "x²+xy+xy²", "9 order twins · leading types 1+8", 8, 9,
        6_322_176, "#7c3aed",
    ),
    CubicQ8GeneratedClass(
        "g16", "x+x²+xy+xy²", "2 order twins · leading types 1+1", 8, 2,
        6_322_176, "#7c3aed",
    ),
    CubicQ8GeneratedClass("g17", "x²y+xy²", "three split leading roots", 8, 1, 301_056, "#7c3aed"),
    CubicQ8GeneratedClass("g18", "xy+x²y+xy²", "three split leading roots", 8, 1, 2_107_392, "#7c3aed"),
    CubicQ8GeneratedClass(
        "g19", "x²y+xy²+y³", "linear × irreducible quadratic", 8, 1,
        903_168, "#a21caf",
    ),
    CubicQ8GeneratedClass(
        "g20", "x³+xy²+αy³", "homogeneous irreducible cubic", 8, 1,
        602_112, "#6d28d9",
    ),
    CubicQ8GeneratedClass(
        "g21", "x+x³+xy²+αy³", "24 order-twin irreducible-leading classes",
        8, 24, 12_644_352, "#6d28d9",
    ),
    CubicQ8GeneratedClass(
        "g22", "x³+xy²+αy³+xy+(α²+α+1)x+(α²+1)y",
        "seven-valued irreducible-leading cubic", 7, 1, 4_214_784,
        "#9333ea",
    ),
)


CUBIC_Q8_QUADRATIC_COVERS: tuple[tuple[str, str], ...] = tuple(
    (f"g{source}", f"g{target}")
    for source, target in (
        (1,4), (1,12), (1,15), (1,16), (1,17), (1,18), (1,19),
        (2,0), (3,0), (4,2), (4,3), (5,0),
        (6,1), (6,8), (6,9), (6,10), (6,11), (7,0),
        (8,2), (8,5), (8,7), (9,3), (9,5), (9,7),
        (10,3), (10,7), (11,2), (11,3), (11,5), (11,7),
        (12,2), (12,3), (12,5), (12,7), (13,1), (14,1),
        (15,2), (15,3), (15,5), (15,7),
        (16,3), (16,5), (16,7), (17,3), (17,5), (17,7),
        (18,2), (18,3), (18,7), (19,3), (19,5), (19,7),
        (20,5), (20,7), (21,5), (21,7), (22,7),
    )
)


def cubic_q8_map_count() -> int:
    """Return the number of degree-at-most-three maps over ``F_8``."""

    return 8**10


def cubic_q8_generated_class_count() -> int:
    """Return the number of SCCs in the quadratic-input generated preorder."""

    return sum(item.class_count for item in CUBIC_Q8_QUADRATIC_CLASSES)


def render_cubic_q8_quadratic_poset_svg(output: Path) -> Path:
    """Render the compressed 110-element Hasse diagram using Graphviz."""

    if sum(item.total_size for item in CUBIC_Q8_QUADRATIC_CLASSES) != cubic_q8_map_count():
        raise AssertionError("generated-class sizes must partition all cubic maps")
    if cubic_q8_generated_class_count() != 110:
        raise AssertionError("expected 110 generated-preorder classes")

    levels = (
        ("g6", "g13", "g14"),
        ("g1",),
        (
            "g4", "g8", "g9", "g10", "g11", "g12", "g15", "g16",
            "g17", "g18", "g19", "g20", "g21", "g22",
        ),
        ("g2", "g3", "g5", "g7"),
        ("g0",),
    )
    rows = []
    for item in CUBIC_Q8_QUADRATIC_CLASSES:
        if item.class_count == 1:
            size_label = f"|class| = {item.class_size:,}"
        else:
            size_label = f"{item.class_count} classes × {item.class_size:,} each"
        rows.append(
            f'''  {item.key} [color="{item.color}", fillcolor="{item.color}18", label=<
    <TABLE BORDER="0" CELLBORDER="0" CELLPADDING="2">
      <TR><TD><FONT POINT-SIZE="17"><B>{escape(item.representative)}</B></FONT></TD></TR>
      <TR><TD><FONT POINT-SIZE="11">{escape(item.name)} · image {item.image_size}</FONT></TD></TR>
      <TR><TD><FONT POINT-SIZE="13"><B>{escape(size_label)}</B></FONT></TD></TR>
    </TABLE>>];'''
        )
    edges = "\n".join(
        f"  {source} -> {target};"
        for source, target in CUBIC_Q8_QUADRATIC_COVERS
    )
    ranks = "\n".join(
        "  { rank=same; " + "; ".join(level) + "; }" for level in levels
    )
    dot = f'''digraph cubic_q8_quadratic {{
  graph [rankdir=TB, bgcolor="white", pad="0.24", nodesep="0.20",
         ranksep="0.72 equally", splines="polyline",
         label="Cubic maps F₈² → F₈ · quadratic input processors\n110 generated-preorder classes · 23 boxes after grouping order twins",
         labelloc="t", labeljust="c", fontsize="24", fontname="Arial"];
  node [shape=box, style="rounded,filled", penwidth="1.7",
        fontname="Arial", margin="0.08,0.05"];
  edge [color="#64748b99", penwidth="1.25", arrowsize="0.68"];
{chr(10).join(rows)}
{ranks}
{edges}
}}
'''
    output.parent.mkdir(parents=True, exist_ok=True)
    try:
        subprocess.run(
            ["dot", "-Tsvg", "-o", str(output)],
            input=dot,
            text=True,
            check=True,
        )
    except FileNotFoundError as error:
        raise RuntimeError("Graphviz 'dot' is required for the q=8 diagram") from error
    return output
