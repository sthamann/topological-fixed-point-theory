"""Exact finite-rotation controls for the TFPT origin exploration."""

from collections import Counter
from itertools import permutations, product

from sympy import I, Matrix, Rational, ZZ, diag, expand, eye, sign, sqrt, zeros
from sympy.matrices.normalforms import smith_normal_form


def exact(matrix):
    return matrix.applyfunc(expand)


def key(matrix):
    return tuple(matrix)


def square_action(rotation):
    columns = []
    for form in FORMS:
        image = exact(rotation * form * rotation.T)
        columns.append(Matrix([image[0, 0], -image[2, 2], image[0, 1],
                               image[0, 2], image[1, 2]]))
    return Matrix.hstack(*columns)


PHI = (1 + sqrt(5)) / 2
FORMS = [diag(1, -1, 0), diag(0, 1, -1)]
for first, second in [(0, 1), (0, 2), (1, 2)]:
    symmetric = zeros(3)
    symmetric[first, second] = symmetric[second, first] = 1
    FORMS.append(symmetric)


def signed_permutations():
    tetrahedral, octahedral = [], []
    for indices in permutations(range(3)):
        permutation = eye(3)[:, indices]
        for signs in product([-1, 1], repeat=3):
            rotation = permutation * diag(*signs)
            if rotation.det() == 1:
                octahedral.append(rotation)
                if permutation.det() == 1:
                    tetrahedral.append(rotation)
    return tetrahedral, octahedral


def icosahedral_rotations():
    vertices = []
    for first, second in product([-1, 1], repeat=2):
        coordinates = (0, first, expand(second * PHI))
        vertices.extend(Matrix(coordinates[k:] + coordinates[:k]) for k in range(3))

    initial = vertices[0]
    initial_neighbor = next(v for v in vertices if expand(initial.dot(v) - PHI) == 0)
    initial_frame = Matrix.hstack(initial, initial_neighbor, initial.cross(initial_neighbor))
    inverse = initial_frame.inv().applyfunc(lambda entry: expand(entry.ratsimp()).simplify())

    rotations = []
    for vertex in vertices:
        for neighbor in vertices:
            if expand(vertex.dot(neighbor) - PHI) != 0:
                continue
            frame = Matrix.hstack(vertex, neighbor, exact(vertex.cross(neighbor)))
            rotation = exact(frame * inverse)
            rotations.append(rotation)
    vertex_keys = {key(vertex) for vertex in vertices}
    assert len({key(rotation) for rotation in rotations}) == 60
    for rotation in rotations:
        assert exact(rotation * rotation.T) == eye(3)
        assert expand(rotation.det()) == 1
        assert {key(exact(rotation * vertex)) for vertex in vertices} == vertex_keys
    return rotations


def character_data(group):
    actions = [square_action(rotation) for rotation in group]
    norm = sum(expand(action.trace() ** 2) for action in actions) / len(group)
    invariant_rank = sum(action.trace() for action in actions) / len(group)
    return expand(norm), expand(invariant_rank)


def marked_planes(group):
    subgroup = [eye(3), diag(1, -1, -1), diag(-1, 1, -1), diag(-1, -1, 1)]
    base_plane = sum([square_action(rotation) for rotation in subgroup], zeros(5)) / 4
    assert base_plane == diag(1, 1, 0, 0, 0)
    cycle = eye(3)[:, [1, 2, 0]]
    omega = (-1 + I * sqrt(3)) / 2
    cycle_action = square_action(cycle)
    base_line = exact((eye(5) + omega.conjugate() * cycle_action
                       + omega * cycle_action ** 2) * base_plane / 3)
    assert exact(base_line ** 2) == base_line and base_line.trace() == 1
    base_axis_lines = [eye(5)[:, index] * eye(5)[index, :] for index in range(2, 5)]

    frames = {}
    for rotation in group:
        subgroup_key = frozenset(key(exact(rotation * item * rotation.T)) for item in subgroup)
        if subgroup_key in frames:
            assert exact(square_action(rotation) * base_line * square_action(rotation.T)) == frames[subgroup_key][1]
            continue
        action = square_action(rotation)
        inverse_action = square_action(rotation.T)
        plane = exact(action * base_plane * inverse_action)
        line = exact(action * base_line * inverse_action)
        axis_lines = [exact(action * axis_line * inverse_action) for axis_line in base_axis_lines]
        axis_directions = [rotation[:, index] for index in (2, 1, 0)]
        frames[subgroup_key] = plane, line, axis_lines, axis_directions

    assert len(frames) == 5
    planes = [frame[0] for frame in frames.values()]
    lines = [frame[1] for frame in frames.values()]
    assert sum(planes, zeros(5)) == 2 * eye(5)
    assert sum(lines, zeros(5)) == eye(5)
    for first, plane in enumerate(planes):
        assert exact(plane ** 2) == plane and plane.trace() == 2
        opposite_line = exact(plane - lines[first])
        assert lines[first].conjugate() == opposite_line
        for second, other_plane in enumerate(planes):
            if first == second:
                assert expand((lines[second] * opposite_line).trace()) == 0
                continue
            assert exact(plane * other_plane * plane) == plane / 4
            assert exact(lines[first] * lines[second]) == zeros(5)
            assert expand((opposite_line * lines[second]).trace()) == Rational(1, 4)

    effect_transition = Matrix(5, 5, lambda i, j: expand((planes[i] * planes[j]).trace() / 4))
    assert effect_transition.eigenvals() == {Rational(3, 8): 4, 1: 1}
    print('A5 frames:', len(frames), 'ranks: 2+3', 'sum P_i = 2 I', 'P_i P_j P_i = P_i/4')
    print('Conjugate line bases: same orbit orthonormal, cross orbit squared overlaps 0 / 1/4')
    print('POVM transition:', effect_transition[0, :], 'spectrum:', effect_transition.eigenvals())

    # A single real rank-one measurement selects the frame and its axis together.
    axis_lines = [line for frame in frames.values() for line in frame[2]]
    assert len({key(line) for line in axis_lines}) == 15
    assert sum(axis_lines, zeros(5)) == 3 * eye(5)
    joint_transition = Matrix(15, 15, lambda i, j: expand((axis_lines[i] * axis_lines[j]).trace() / 3))
    row_multisets = [Counter(joint_transition[i, :]) for i in range(15)]
    assert all(row == row_multisets[0] for row in row_multisets)
    assert joint_transition == joint_transition.T and joint_transition * Matrix.ones(15, 1) == Matrix.ones(15, 1)
    print('Joint rank-one POVM: 15 uniform labels; transition value multiplicities:', dict(row_multisets[0]))
    print('Joint transition characteristic polynomial:', joint_transition.charpoly().as_expr().factor())

    # The flags form one A5 orbit, so use the first frame and its first axis.
    # A second observed axis can orient the first by a quartic invariant of
    # the unordered perpendicular axis pair. The sign of the second axis cancels.
    axis, first_perpendicular, second_perpendicular = next(iter(frames.values()))[3]
    masses = Counter()
    row = joint_transition[0, :]
    axes = [direction for frame in frames.values() for direction in frame[3]]
    for probability, next_axis in zip(row, axes):
        projected = exact(next_axis - axis.dot(next_axis) * axis)
        gradient = (projected.dot(first_perpendicular) ** 3 * first_perpendicular
                    + projected.dot(second_perpendicular) ** 3 * second_perpendicular)
        orientation = sign(expand(axis.dot(projected.cross(gradient))))
        masses[orientation] += probability
    assert masses == {-1: Rational(1, 3), 0: Rational(1, 3), 1: Rational(1, 3)}
    print('Quarter-turn orientation from ordered flags: masses (-,0,+) = (1/3,1/3,1/3)')


def axis_responses():
    axes = [Matrix([1, y, z]) for y, z in product([-1, 1], repeat=2)]
    values = Matrix([[axis.dot(form * axis) for form in FORMS] for axis in axes])
    projectors = [axis * axis.T / 3 - eye(3) / 3 for axis in axes]
    gram = Matrix(4, 4, lambda i, j: (projectors[i] * projectors[j]).trace())
    assert gram == (8 * eye(4) - 2 * Matrix.ones(4)) / 9

    def jet(axis):
        return Matrix.hstack(*(form * axis for form in FORMS))

    ranks = values.rank(), jet(axes[0]).rank(), Matrix.vstack(jet(axes[0]), jet(axes[1])).rank()
    assert ranks == (3, 3, 5)

    # An order-four lift of each half-turn belongs to Q8.  The order-three
    # generator cycles the three half-turn axes, so it also cycles their lines.
    sign_characters = [square_action(item).diagonal() for item in
                       [diag(1, -1, -1), diag(-1, 1, -1), diag(-1, -1, 1)]]
    assert len({tuple(character) for character in sign_characters}) == 3
    print('T frame: 4 projective marks, family Gram diagonal 2/3 and off diagonal -2/9')
    print('Response ranks (values on 4 marks, jet at 1, jets at 2):', ranks)
    print('Three half-turns lift to conjugate C4 subgroups; 2T/Q8=C3')


def six_axis_matchings(group):
    five_groups = {}
    for rotation in group:
        if expand(rotation.trace()) not in (expand(PHI), expand(1 - PHI)):
            continue
        elements = [eye(3), rotation]
        for _ in range(3):
            elements.append(exact(rotation * elements[-1]))
        assert exact(rotation * elements[-1]) == eye(3)
        five_groups[frozenset(key(item) for item in elements)] = elements
    assert len(five_groups) == 6

    five_indices = {subgroup: index for index, subgroup in enumerate(five_groups)}
    action = {}
    for rotation in group:
        permutation = tuple(five_indices[frozenset(key(exact(rotation * item * rotation.T))
                                                   for item in subgroup)]
                            for subgroup in five_groups.values())
        action[key(rotation)] = permutation
        assert sum(i == j for i, j in enumerate(permutation)) - 1 == square_action(rotation).trace()

    subgroup = [eye(3), diag(1, -1, -1), diag(-1, 1, -1), diag(-1, -1, 1)]
    two_groups = {frozenset(key(exact(rotation * item * rotation.T)) for item in subgroup)
                  for rotation in group}
    assert len(two_groups) == 5
    all_edges = set()
    for two_group in two_groups:
        permutations_of_six = [action[item] for item in two_group]
        orbits = {tuple(sorted({permutation[index] for permutation in permutations_of_six}))
                  for index in range(6)}
        assert len(orbits) == 3 and all(len(orbit) == 2 for orbit in orbits)
        signs = {tuple(1 if permutation[first] == first else -1
                       for permutation in permutations_of_six) for first, _ in orbits}
        assert len(signs) == 3 and all(sum(character) == 0 for character in signs)
        assert all_edges.isdisjoint(orbits)
        all_edges.update(orbits)
    assert all_edges == {(first, second) for first in range(6) for second in range(first + 1, 6)}
    print('Six Sylow-5 labels: exact permutation character 1+W; five V4 matchings partition 15 K6 edges')


def equatorial_marks(group):
    quarter_turns = [rotation for rotation in group
                     if rotation ** 4 == eye(3) and rotation ** 2 != eye(3)]
    axes = {tuple((rotation - eye(3)).nullspace()[0]) for rotation in quarter_turns}
    assert len(quarter_turns) == 6 and len(axes) == 3
    assert all(sum(value != 0 for value in axis) == 1 for axis in axes)

    # Once a quarter-turn axis is selected, the two other axes determine a
    # square on its oriented equator. A stereographic chart sets its radius to 1.
    cross_ratios = set()
    for a, b, c, d in permutations([1, I, -1, -I]):
        cross_ratios.add(((a - c) * (b - d) / ((a - d) * (b - c))).simplify())
    assert cross_ratios == {-1, Rational(1, 2), 2}
    j_values = {256 * (1 - value + value ** 2) ** 3 / (value ** 2 * (1 - value) ** 2)
                for value in cross_ratios}
    assert j_values == {1728}
    print('Octahedral stabilizer: 6 quarter turns on 3 axes, equator cross ratios:', sorted(cross_ratios), 'j:', 1728)

    axis = Matrix([0, 0, 1])
    tetrahedral_directions = [Matrix(signs) for signs in product([-1, 1], repeat=3)]
    square = {key(exact((direction - axis.dot(direction) * axis) / sqrt(2)))
              for direction in tetrahedral_directions}
    assert len(square) == 4
    square = [Matrix(point) for point in sorted(square)]

    def square_values(forms, points):
        return Matrix([[expand(point.dot(axis.cross(form * axis)) + point.dot(form * point))
                        for form in forms] for point in points])

    stabilizer = [rotation for rotation in group if rotation * axis in (axis, -axis)]
    model = square_values(FORMS[2:], square)
    assert len(stabilizer) == 8 and model.rank() == 3
    assert Matrix.ones(1, 4) * model == zeros(1, 3)
    for rotation in stabilizer:
        assert square_values([rotation * form * rotation.T for form in FORMS[2:]], square) == square_values(
            FORMS[2:], [rotation.T * point for point in square])
    print('Projected tetrahedral marks: 4 square points, D4-equivariant zero-sum response rank:', model.rank())


def carrier_and_spin_controls():
    """Compare actual actions and character ranks after obtaining the square."""
    quarter_turn = Matrix([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
    half_turn = diag(1, -1, -1)
    rotation = square_action(quarter_turn)
    reflection = square_action(half_turn)

    # On D={z^4=-1}, ordinary sections are spanned by z^j/(1+z^4).
    section_rotation = diag(*(I ** j for j in range(5)))
    section_reflection = Matrix(5, 5, lambda i, j: int(i + j == 4))
    group_words = list(product(range(4), range(2)))
    matches = []
    for rotation_sign, reflection_sign in product((1, -1), repeat=2):
        if all((rotation ** k * reflection ** l).trace()
               == ((rotation_sign * section_rotation) ** k
                   * (reflection_sign * section_reflection) ** l).trace()
               for k, l in group_words):
            matches.append((rotation_sign, reflection_sign))
    assert rotation.trace() == -1 and section_rotation.trace() == 1
    assert matches == [(-1, 1)]
    print('D4 quarter trace: quadrupole -1, ordinary sections +1; unique character twist (r,s)=(-1,+1)')

    # The positive round-S^2 Dirac level k has spin k+1/2. A scalar phase
    # changes the order-eight lift to order four by shifting m to m+1/2.
    rank_rows = []
    for level in range(4):
        ranks = Counter(weight % 4 for weight in range(-level, level + 2))
        rank_rows.append(tuple(ranks[charge] for charge in range(4)))
    assert rank_rows == [(1, 1, 0, 0), (1, 1, 1, 1), (1, 1, 2, 2), (2, 2, 2, 2)]
    print('Auxiliary round S^2 positive Dirac character ranks (k=0..3):', rank_rows)


def mckay_icosahedral(group):
    angle_counts = Counter(expand(rotation.trace()) for rotation in group)
    assert angle_counts == Counter({3: 1, -1: 15, 0: 20, expand(PHI): 12, expand(1-PHI): 12})
    lift_trace = [2, 0, 1, expand(PHI), expand(PHI - 1)]
    perm_character = [4, 0, 1, -1, -1]  # A5 action on 5 Sylow 2-subgroups, less the trivial line.
    trace_counts = [angle_counts[3], angle_counts[-1], angle_counts[0],
                    angle_counts[expand(PHI)], angle_counts[expand(1 - PHI)]]
    weights = trace_counts + trace_counts
    trace = Matrix(lift_trace + [-item for item in lift_trace])

    def inner(left, right):
        return expand(sum(w * a * b for w, a, b in zip(weights, left, right)) / 120)

    def multiply(left, right):
        return exact(left.multiply_elementwise(right))

    # These are symmetric powers of the defining SU(2) lift; the other
    # characters come from Galois conjugation and the action on five Sylow groups.
    chars = [Matrix.ones(10, 1), trace]
    for _ in range(4):
        chars.append(multiply(trace, chars[-1]) - chars[-2])
    chars += [Matrix(perm_character * 2), chars[2].subs(sqrt(5), -sqrt(5)),
              trace.subs(sqrt(5), -sqrt(5))]
    dimensions = [int(character[0]) for character in chars]
    assert sum(dimension ** 2 for dimension in dimensions) == 120
    assert Matrix(9, 9, lambda i, j: inner(chars[i], chars[j])) == eye(9)
    adjacency = Matrix(9, 9, lambda i, j: inner(multiply(trace, chars[i]), chars[j]))
    assert adjacency == adjacency.T and set(adjacency) == {0, 1}
    gram = 2 * eye(9) - adjacency
    assert gram * Matrix(dimensions) == zeros(9, 1)
    root_gram = gram[1:, 1:]
    assert root_gram.det() == 1 and root_gram.is_positive_definite
    edges = [(i, j) for i in range(9) for j in range(i+1, 9) if adjacency[i, j]]
    assert edges == [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (5, 7), (6, 8)]
    a3_gram = gram.extract([0, 1, 2], [0, 1, 2])
    d5_gram = gram.extract([4, 5, 6, 7, 8], [4, 5, 6, 7, 8])
    assert a3_gram.det() == d5_gram.det() == 4
    root_vectors = [-Matrix(dimensions[1:]), *(eye(8)[:, index] for index in range(8))]
    retained = [0, 1, 2, 4, 5, 6, 7, 8]
    sublattice = Matrix.hstack(*(root_vectors[index] for index in retained))
    assert abs(sublattice.det()) == 4
    assert [abs(item) for item in smith_normal_form(sublattice, domain=ZZ).diagonal()] == [1] * 7 + [4]
    a3_projection = a3_gram.inv() * gram.extract([0, 1, 2], [3])
    d5_projection = d5_gram.inv() * gram.extract([4, 5, 6, 7, 8], [3])
    glue_norms = (a3_projection.dot(a3_gram * a3_projection), d5_projection.dot(d5_gram * d5_projection))
    assert glue_norms == (Rational(3, 4), Rational(5, 4))
    print('Binary I McKay dimensions:', dimensions, 'edges:', edges)
    print('E8 rank/det:', root_gram.rows, root_gram.det(), 'canonical chain deletion: A3 + D5, quotient C4')
    print('Glue projection norms:', glue_norms, 'sum:', sum(glue_norms))


def main():
    tetrahedral, octahedral = signed_permutations()
    icosahedral = icosahedral_rotations()
    data = [character_data(group) for group in (tetrahedral, octahedral, icosahedral)]
    assert data == [(3, 0), (2, 0), (1, 0)]
    print('T, O, I orders:', [len(tetrahedral), len(octahedral), len(icosahedral)])
    print('Quadrupole character norm and trivial multiplicity:', data)
    print('T and O preserve the coordinate line system; I is primitive')
    marked_planes(icosahedral)
    axis_responses()
    six_axis_matchings(icosahedral)
    equatorial_marks(octahedral)
    carrier_and_spin_controls()
    mckay_icosahedral(icosahedral)


if __name__ == '__main__':
    main()
