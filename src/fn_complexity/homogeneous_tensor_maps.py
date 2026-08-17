"""Homogeneous tensor-map degeneration posets over the field F_3.

The objects are tensors ``Sym^degree(V)^* -> W``.  Equivalence uses
invertible linear changes of coordinates in ``V`` and ``W``; degeneration
allows arbitrary, possibly singular, linear processors on both sides.
"""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass
from functools import lru_cache
from html import escape
from itertools import combinations, product
from pathlib import Path
import subprocess

import numpy as np


P = 3


@dataclass(frozen=True)
class TensorCase:
    number: int
    degree: int
    input_dimension: int
    output_dimension: int
    short_name: str
    title: str

    @property
    def monomial_count(self) -> int:
        return len(monomial_exponents(self.input_dimension, self.degree))

    @property
    def tensor_dimension(self) -> int:
        return self.output_dimension * self.monomial_count

    @property
    def tensor_count(self) -> int:
        return P**self.tensor_dimension


TENSOR_CASES: dict[int, TensorCase] = {
    1: TensorCase(
        1, 2, 3, 1, "quadratic-ternary-form",
        "Quadratic ternary forms over F₃",
    ),
    2: TensorCase(
        2, 2, 2, 2, "quadratic-p1-map",
        "Quadratic homogeneous maps F₃² → F₃²",
    ),
    3: TensorCase(
        3, 2, 3, 3, "quadratic-p2-map",
        "Quadratic homogeneous maps F₃³ → F₃³",
    ),
    4: TensorCase(
        4, 3, 3, 1, "cubic-ternary-form",
        "Cubic ternary forms over F₃",
    ),
    5: TensorCase(
        5, 3, 2, 2, "cubic-p1-map",
        "Cubic homogeneous maps F₃² → F₃²",
    ),
    6: TensorCase(
        6, 3, 3, 3, "cubic-p2-map",
        "Cubic homogeneous maps F₃³ → F₃³",
    ),
}


@dataclass(frozen=True)
class TensorOrbit:
    key: str
    representative: tuple[tuple[int, ...], ...]
    size: int
    output_rank: int
    projective_basepoints: int


@dataclass(frozen=True)
class TensorPoset:
    case: TensorCase
    orbits: tuple[TensorOrbit, ...]
    covers: tuple[tuple[str, str], ...]


@lru_cache(maxsize=None)
def monomial_exponents(dimension: int, degree: int) -> tuple[tuple[int, ...], ...]:
    """Return degree-``degree`` exponent vectors in readable lexicographic order."""

    def visit(prefix: tuple[int, ...], remaining_slots: int, remaining: int):
        if remaining_slots == 1:
            yield prefix + (remaining,)
            return
        for value in range(remaining, -1, -1):
            yield from visit(prefix + (value,), remaining_slots - 1, remaining - value)

    return tuple(visit((), dimension, degree))


def _rref(rows: np.ndarray, width: int | None = None) -> tuple[tuple[int, ...], ...]:
    """Canonical reduced row-echelon basis over F_3."""

    if rows.size == 0:
        return ()
    matrix = np.asarray(rows, dtype=np.int16).copy() % P
    if matrix.ndim == 1:
        matrix = matrix.reshape(1, -1)
    columns = matrix.shape[1] if width is None else width
    pivot_row = 0
    for column in range(columns):
        pivot = next(
            (row for row in range(pivot_row, matrix.shape[0]) if matrix[row, column]),
            None,
        )
        if pivot is None:
            continue
        if pivot != pivot_row:
            matrix[[pivot_row, pivot]] = matrix[[pivot, pivot_row]]
        if matrix[pivot_row, column] == 2:
            matrix[pivot_row] = (2 * matrix[pivot_row]) % P
        for row in range(matrix.shape[0]):
            if row != pivot_row and matrix[row, column]:
                matrix[row] = (
                    matrix[row] - matrix[row, column] * matrix[pivot_row]
                ) % P
        pivot_row += 1
        if pivot_row == matrix.shape[0]:
            break
    return tuple(tuple(int(value) for value in row) for row in matrix[:pivot_row])


def matrix_rank(matrix: tuple[int, ...] | np.ndarray, dimension: int) -> int:
    array = np.asarray(matrix, dtype=np.int16).reshape(dimension, dimension)
    return len(_rref(array))


@lru_cache(maxsize=None)
def all_matrices(dimension: int) -> tuple[tuple[int, ...], ...]:
    return tuple(product(range(P), repeat=dimension * dimension))


@lru_cache(maxsize=None)
def invertible_matrices(dimension: int) -> tuple[tuple[int, ...], ...]:
    return tuple(
        matrix
        for matrix in all_matrices(dimension)
        if matrix_rank(matrix, dimension) == dimension
    )


def _matrix_product(left: tuple[int, ...], right: tuple[int, ...], dimension: int):
    a = np.asarray(left, dtype=np.int16).reshape(dimension, dimension)
    b = np.asarray(right, dtype=np.int16).reshape(dimension, dimension)
    return tuple(int(value) for value in ((a @ b) % P).flat)


@lru_cache(maxsize=None)
def gl_generators(dimension: int) -> tuple[tuple[int, ...], ...]:
    """A checked elementary generating set for GL(dimension, 3)."""

    if dimension == 1:
        return ((2,),)
    generators: list[tuple[int, ...]] = []
    for index in range(dimension - 1):
        matrix = np.eye(dimension, dtype=np.int16)
        matrix[[index, index + 1]] = matrix[[index + 1, index]]
        generators.append(tuple(int(value) for value in matrix.flat))
    scale = np.eye(dimension, dtype=np.int16)
    scale[0, 0] = 2
    generators.append(tuple(int(value) for value in scale.flat))
    shear = np.eye(dimension, dtype=np.int16)
    shear[0, 1] = 1
    generators.append(tuple(int(value) for value in shear.flat))

    identity = tuple(int(value) for value in np.eye(dimension, dtype=np.int16).flat)
    generated = {identity}
    queue = deque([identity])
    while queue:
        current = queue.popleft()
        for generator in generators:
            moved = _matrix_product(current, generator, dimension)
            if moved not in generated:
                generated.add(moved)
                queue.append(moved)
    if generated != set(invertible_matrices(dimension)):
        raise AssertionError("elementary matrices did not generate the full GL group")
    return tuple(generators)


def _poly_multiply(
    left: dict[tuple[int, ...], int],
    right: dict[tuple[int, ...], int],
) -> dict[tuple[int, ...], int]:
    output: dict[tuple[int, ...], int] = {}
    for first, a in left.items():
        for second, b in right.items():
            exponent = tuple(x + y for x, y in zip(first, second))
            output[exponent] = (output.get(exponent, 0) + a * b) % P
    return {key: value for key, value in output.items() if value}


@lru_cache(maxsize=None)
def input_action(
    dimension: int,
    degree: int,
    matrix: tuple[int, ...],
) -> np.ndarray:
    """Matrix for pulling homogeneous forms back along one linear map."""

    exponents = monomial_exponents(dimension, degree)
    exponent_index = {exponent: index for index, exponent in enumerate(exponents)}
    linear = np.asarray(matrix, dtype=np.int16).reshape(dimension, dimension)
    action = np.zeros((len(exponents), len(exponents)), dtype=np.int8)
    zero = (0,) * dimension
    for source_index, source_exponent in enumerate(exponents):
        polynomial: dict[tuple[int, ...], int] = {zero: 1}
        for source_variable, power in enumerate(source_exponent):
            factor = {
                tuple(1 if target == column else 0 for target in range(dimension)):
                int(linear[source_variable, column])
                for column in range(dimension)
                if linear[source_variable, column]
            }
            for _ in range(power):
                polynomial = _poly_multiply(polynomial, factor)
        for exponent, coefficient in polynomial.items():
            action[source_index, exponent_index[exponent]] = coefficient
    return action


def _encode_tensor(tensor: np.ndarray) -> int:
    code = 0
    multiplier = 1
    for value in tensor.flat:
        code += int(value) * multiplier
        multiplier *= P
    return code


def _decode_tensor(code: int, rows: int, columns: int) -> np.ndarray:
    values = np.zeros(rows * columns, dtype=np.int8)
    for index in range(len(values)):
        values[index] = code % P
        code //= P
    return values.reshape(rows, columns)


def _apply_input(tensor: np.ndarray, action: np.ndarray) -> np.ndarray:
    return ((tensor.astype(np.int16) @ action.astype(np.int16)) % P).astype(np.int8)


def _apply_output(tensor: np.ndarray, matrix: tuple[int, ...]) -> np.ndarray:
    rows = tensor.shape[0]
    output = np.asarray(matrix, dtype=np.int16).reshape(rows, rows)
    return ((output @ tensor.astype(np.int16)) % P).astype(np.int8)


def _projective_points(dimension: int) -> tuple[tuple[int, ...], ...]:
    points = []
    for vector in product(range(P), repeat=dimension):
        if not any(vector):
            continue
        first = next(value for value in vector if value)
        normalized = tuple((value * (1 if first == 1 else 2)) % P for value in vector)
        if normalized == vector:
            points.append(vector)
    return tuple(points)


def _evaluate_form(coefficients: tuple[int, ...], point: tuple[int, ...], degree: int) -> int:
    exponents = monomial_exponents(len(point), degree)
    return sum(
        coefficient
        * np.prod([pow(value, power, P) for value, power in zip(point, exponent)])
        for coefficient, exponent in zip(coefficients, exponents)
    ) % P


def _basepoint_count(tensor: tuple[tuple[int, ...], ...], case: TensorCase) -> int:
    return sum(
        all(_evaluate_form(row, point, case.degree) == 0 for row in tensor)
        for point in _projective_points(case.input_dimension)
    )


@lru_cache(maxsize=None)
def all_subspaces(width: int, maximum_dimension: int) -> tuple[tuple[tuple[int, ...], ...], ...]:
    """Enumerate every RREF subspace of F_3^width up to the requested dimension."""

    subspaces: list[tuple[tuple[int, ...], ...]] = [()]
    for dimension in range(1, maximum_dimension + 1):
        for pivots in combinations(range(width), dimension):
            allowed = [
                (row, column)
                for row, pivot in enumerate(pivots)
                for column in range(pivot + 1, width)
                if column not in pivots
            ]
            for values in product(range(P), repeat=len(allowed)):
                matrix = np.zeros((dimension, width), dtype=np.int8)
                for row, pivot in enumerate(pivots):
                    matrix[row, pivot] = 1
                for (row, column), value in zip(allowed, values):
                    matrix[row, column] = value
                subspaces.append(tuple(tuple(int(x) for x in row) for row in matrix))
    return tuple(subspaces)


def gaussian_binomial(n: int, r: int) -> int:
    numerator = 1
    denominator = 1
    for index in range(r):
        numerator *= P**n - P**index
        denominator *= P**r - P**index
    return numerator // denominator


def spanning_tuple_count(output_dimension: int, rank: int) -> int:
    if rank == 0:
        return 1
    result = 1
    for index in range(rank):
        result *= P**output_dimension - P**index
    return result


def _subspaces_of_span(basis: tuple[tuple[int, ...], ...]):
    rank = len(basis)
    if rank == 0:
        return ((),)
    basis_array = np.asarray(basis, dtype=np.int16)
    output = []
    for coefficient_basis in all_subspaces(rank, rank):
        if not coefficient_basis:
            output.append(())
            continue
        coefficients = np.asarray(coefficient_basis, dtype=np.int16)
        output.append(_rref((coefficients @ basis_array) % P))
    return tuple(output)


def _hasse_covers(adjacency: list[set[int]]) -> tuple[tuple[int, int], ...]:
    closure: list[set[int]] = []
    for start in range(len(adjacency)):
        seen = {start}
        stack = [start]
        while stack:
            source = stack.pop()
            for target in adjacency[source]:
                if target not in seen:
                    seen.add(target)
                    stack.append(target)
        closure.append(seen)
    for source, targets in enumerate(closure):
        for target in targets:
            if target != source and source in closure[target]:
                raise AssertionError("mutual reachability exceeded invertible equivalence")
    strict = {
        (source, target)
        for source, targets in enumerate(closure)
        for target in targets
        if source != target
    }
    return tuple(sorted(
        (source, target)
        for source, target in strict
        if not any(
            (source, middle) in strict and (middle, target) in strict
            for middle in range(len(adjacency))
            if middle not in (source, target)
        )
    ))


def _small_tensor_poset(case: TensorCase) -> TensorPoset:
    rows = case.output_dimension
    columns = case.monomial_count
    total = case.tensor_count
    orbit_index = np.full(total, -1, dtype=np.int32)
    representatives: list[int] = []
    orbit_sizes: list[int] = []
    input_actions = tuple(
        input_action(case.input_dimension, case.degree, matrix)
        for matrix in gl_generators(case.input_dimension)
    )
    output_generators = gl_generators(rows)

    for start in range(total):
        if orbit_index[start] >= 0:
            continue
        orbit_id = len(representatives)
        representatives.append(start)
        orbit_index[start] = orbit_id
        queue = deque([start])
        size = 0
        while queue:
            code = queue.popleft()
            size += 1
            tensor = _decode_tensor(code, rows, columns)
            neighbours = [
                _apply_input(tensor, action) for action in input_actions
            ] + [
                _apply_output(tensor, matrix) for matrix in output_generators
            ]
            for neighbour in neighbours:
                target = _encode_tensor(neighbour)
                if orbit_index[target] < 0:
                    orbit_index[target] = orbit_id
                    queue.append(target)
        orbit_sizes.append(size)

    all_actions = tuple(
        input_action(case.input_dimension, case.degree, matrix)
        for matrix in all_matrices(case.input_dimension)
    )
    adjacency = [{orbit_id} for orbit_id in range(len(representatives))]
    tensors = [
        _decode_tensor(code, rows, columns) for code in representatives
    ]
    for orbit_id, tensor in enumerate(tensors):
        transformed_spans = {
            _rref(_apply_input(tensor, action)) for action in all_actions
        }
        for span in transformed_spans:
            for subspace in _subspaces_of_span(span):
                padded = np.zeros((rows, columns), dtype=np.int8)
                if subspace:
                    padded[:len(subspace)] = np.asarray(subspace, dtype=np.int8)
                adjacency[orbit_id].add(int(orbit_index[_encode_tensor(padded)]))

    covers = _hasse_covers(adjacency)
    orbits = []
    for orbit_id, tensor in enumerate(tensors):
        representative = tuple(tuple(int(value) for value in row) for row in tensor)
        rank = len(_rref(tensor))
        orbits.append(TensorOrbit(
            key=f"c{orbit_id}",
            representative=representative,
            size=orbit_sizes[orbit_id],
            output_rank=rank,
            projective_basepoints=_basepoint_count(representative, case),
        ))
    if sum(orbit.size for orbit in orbits) != total:
        raise AssertionError("orbit sizes do not partition the tensor universe")
    return TensorPoset(
        case=case,
        orbits=tuple(orbits),
        covers=tuple((f"c{source}", f"c{target}") for source, target in covers),
    )


def _quadratic_ternary_vector_poset(case: TensorCase) -> TensorPoset:
    width = case.monomial_count
    subspaces = all_subspaces(width, case.output_dimension)
    subspace_set = set(subspaces)
    orbit_for: dict[tuple[tuple[int, ...], ...], int] = {}
    representatives: list[tuple[tuple[int, ...], ...]] = []
    orbit_subspace_sizes: list[int] = []
    generator_actions = tuple(
        input_action(case.input_dimension, case.degree, matrix)
        for matrix in gl_generators(case.input_dimension)
    )
    for start in subspaces:
        if start in orbit_for:
            continue
        orbit_id = len(representatives)
        representatives.append(start)
        orbit_for[start] = orbit_id
        queue = deque([start])
        size = 0
        while queue:
            subspace = queue.popleft()
            size += 1
            basis = np.asarray(subspace, dtype=np.int16)
            for action in generator_actions:
                moved = _rref((basis @ action) % P) if subspace else ()
                if moved not in orbit_for:
                    if moved not in subspace_set:
                        raise AssertionError("input action left the Grassmannian")
                    orbit_for[moved] = orbit_id
                    queue.append(moved)
        orbit_subspace_sizes.append(size)

    all_actions = tuple(
        input_action(case.input_dimension, case.degree, matrix)
        for matrix in all_matrices(case.input_dimension)
    )
    adjacency = [{orbit_id} for orbit_id in range(len(representatives))]
    for orbit_id, representative in enumerate(representatives):
        basis = np.asarray(representative, dtype=np.int16)
        transformed_spans = {
            _rref((basis @ action) % P) if representative else ()
            for action in all_actions
        }
        for span in transformed_spans:
            for subspace in _subspaces_of_span(span):
                adjacency[orbit_id].add(orbit_for[subspace])

    covers = _hasse_covers(adjacency)
    orbits = []
    for orbit_id, representative in enumerate(representatives):
        rank = len(representative)
        padded = representative + ((0,) * width,) * (case.output_dimension - rank)
        size = orbit_subspace_sizes[orbit_id] * spanning_tuple_count(
            case.output_dimension, rank
        )
        orbits.append(TensorOrbit(
            key=f"c{orbit_id}",
            representative=padded,
            size=size,
            output_rank=rank,
            projective_basepoints=_basepoint_count(padded, case),
        ))
    if sum(orbit.size for orbit in orbits) != case.tensor_count:
        raise AssertionError("Grassmannian orbit sizes do not partition all tensors")
    return TensorPoset(
        case=case,
        orbits=tuple(orbits),
        covers=tuple((f"c{source}", f"c{target}") for source, target in covers),
    )


@lru_cache(maxsize=None)
def compute_tensor_poset(case_number: int) -> TensorPoset:
    case = TENSOR_CASES[case_number]
    if case_number == 6:
        raise ValueError("case 6 has too many equivalence classes for a full diagram")
    if case_number == 3:
        return _quadratic_ternary_vector_poset(case)
    return _small_tensor_poset(case)


def _coefficient_label(coefficient: int, monomial: str) -> str | None:
    if coefficient == 0:
        return None
    if monomial == "1":
        return str(coefficient)
    return monomial if coefficient == 1 else f"−{monomial}"


def _monomial_label(exponent: tuple[int, ...]) -> str:
    variables = "xyzuvw"
    superscripts = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")
    pieces = []
    for variable, power in zip(variables, exponent):
        if power == 1:
            pieces.append(variable)
        elif power:
            pieces.append(variable + str(power).translate(superscripts))
    return "".join(pieces) or "1"


def form_label(coefficients: tuple[int, ...], case: TensorCase) -> str:
    terms = [
        _coefficient_label(coefficient, _monomial_label(exponent))
        for coefficient, exponent in zip(
            coefficients,
            monomial_exponents(case.input_dimension, case.degree),
        )
    ]
    return " + ".join(term for term in terms if term).replace("+ −", "− ") or "0"


def tensor_label(orbit: TensorOrbit, case: TensorCase) -> str:
    nonzero_rows = orbit.representative[:orbit.output_rank]
    labels = [form_label(row, case) for row in nonzero_rows]
    if case.output_dimension == 1:
        return labels[0] if labels else "0"
    return "⟨" + "; ".join(labels) + "⟩" if labels else "0"


def rank_counts_for_case6() -> dict[int, int]:
    case = TENSOR_CASES[6]
    width = case.monomial_count
    return {
        rank: gaussian_binomial(width, rank)
        * spanning_tuple_count(case.output_dimension, rank)
        for rank in range(case.output_dimension + 1)
    }


def minimum_orbit_counts_for_case6() -> dict[int, int]:
    group_order = len(invertible_matrices(3)) ** 2
    return {
        rank: (count + group_order - 1) // group_order
        for rank, count in rank_counts_for_case6().items()
    }


def _graphviz_svg(dot: str, output: Path, *, balance: bool = False) -> Path:
    output.parent.mkdir(parents=True, exist_ok=True)
    try:
        if balance:
            balanced = subprocess.run(
                ["unflatten", "-l", "6", "-c", "8"],
                input=dot,
                text=True,
                capture_output=True,
                check=True,
            ).stdout
        else:
            balanced = dot
        subprocess.run(
            ["dot", "-Tsvg", "-o", str(output)],
            input=balanced,
            text=True,
            check=True,
        )
    except FileNotFoundError as error:
        raise RuntimeError("Graphviz 'dot' is required to render these diagrams") from error
    return output


def render_tensor_poset_svg(
    poset: TensorPoset,
    output: Path,
    *,
    show_titles: bool = True,
) -> Path:
    """Render one exact case-1-through-case-5 Hasse diagram."""

    nodes = []
    for orbit in poset.orbits:
        color = ("#64748b", "#2563eb", "#7c3aed", "#be123c")[orbit.output_rank]
        label = tensor_label(orbit, poset.case)
        if poset.case.output_dimension == 1:
            detail = f"projective zeros {orbit.projective_basepoints}"
        else:
            detail = (
                f"coordinate rank {orbit.output_rank} · "
                f"base points {orbit.projective_basepoints}"
            )
        nodes.append(f'''  {orbit.key} [color="{color}", fillcolor="{color}18", label=<
    <TABLE BORDER="0" CELLBORDER="0" CELLPADDING="2">
      <TR><TD><FONT POINT-SIZE="16"><B>{escape(label)}</B></FONT></TD></TR>
      <TR><TD><FONT POINT-SIZE="11">{detail}</FONT></TD></TR>
      <TR><TD><FONT POINT-SIZE="13"><B>|class| = {orbit.size:,}</B></FONT></TD></TR>
    </TABLE>>];''')
    edges = "\n".join(f"  {source} -> {target};" for source, target in poset.covers)
    graph_title = ""
    if show_titles:
        graph_title = (
            f'         label="{escape(poset.case.title)}\\n'
            f'{len(poset.orbits)} classes · {len(poset.covers)} Hasse covers · '
            'singular linear processors allowed",\n'
            '         labelloc="t", labeljust="c", fontsize="23", '
            'fontname="Arial",\n'
        )
    dot = f'''digraph tensor_poset {{
  graph [rankdir=TB, bgcolor="white", pad="0.24", nodesep="0.20",
         ranksep="0.66 equally", splines="polyline",
{graph_title}         fontname="Arial"];
  node [shape=box, style="rounded,filled", penwidth="1.7",
        fontname="Arial", margin="0.08,0.05"];
  edge [color="#64748b99", penwidth="1.25", arrowsize="0.68"];
{chr(10).join(nodes)}
{edges}
}}
'''
    return _graphviz_svg(dot, output, balance=len(poset.orbits) > 25)


def render_case6_scale_svg(
    output: Path,
    exact_rank1_orbits: int,
    *,
    show_titles: bool = True,
) -> Path:
    """Render the honest size obstruction in place of an impossible full graph."""

    counts = rank_counts_for_case6()
    minima = minimum_orbit_counts_for_case6()
    labels = {
        3: ("output rank 3", counts[3], minima[3]),
        2: ("output rank 2", counts[2], minima[2]),
        1: ("output rank 1", counts[1], exact_rank1_orbits),
        0: ("zero tensor", counts[0], 1),
    }
    nodes = []
    for rank in (3, 2, 1, 0):
        name, tensors, classes = labels[rank]
        if rank == 0:
            tensor_label_text = "1 tensor"
            class_label_text = "1 exact class"
        else:
            tensor_label_text = f"{tensors:,} tensors"
            qualifier = "exact classes" if rank == 1 else "classes at minimum"
            class_label_text = f"{classes:,} {qualifier}"
        color = ("#64748b", "#2563eb", "#7c3aed", "#be123c")[rank]
        nodes.append(f'''  r{rank} [color="{color}", fillcolor="{color}18", label=<
    <TABLE BORDER="0" CELLBORDER="0" CELLPADDING="3">
      <TR><TD><FONT POINT-SIZE="18"><B>{name}</B></FONT></TD></TR>
      <TR><TD><FONT POINT-SIZE="13">{tensor_label_text}</FONT></TD></TR>
      <TR><TD><FONT POINT-SIZE="13"><B>{class_label_text}</B></FONT></TD></TR>
    </TABLE>>];''')
    total_minimum = minima[3] + minima[2] + exact_rank1_orbits + 1
    graph_title = ""
    if show_titles:
        graph_title = (
            '         label="Cubic homogeneous maps F₃³ → F₃³\\n'
            'Coarse output-rank stratification, not a Hasse diagram · at least '
            f'{total_minimum:,} classes",\n'
            '         labelloc="t", labeljust="c", fontsize="24", '
            'fontname="Arial",\n'
        )
    dot = f'''digraph case6_scale {{
  graph [rankdir=TB, bgcolor="white", pad="0.28", nodesep="0.35",
         ranksep="0.82 equally", splines="polyline",
{graph_title}         fontname="Arial"];
  node [shape=box, style="rounded,filled", penwidth="1.8",
        fontname="Arial", margin="0.12,0.08"];
  edge [color="#64748b99", penwidth="1.4", arrowsize="0.72"];
{chr(10).join(nodes)}
  r3 -> r2;
  r2 -> r1;
  r1 -> r0;
}}
'''
    return _graphviz_svg(dot, output)
