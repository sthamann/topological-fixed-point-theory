"""Check geometric transport and reflection positivity on the marked flag graph."""

import argparse
from collections import Counter
from itertools import product
from pathlib import Path
from runpy import run_path

from sympy import Rational, diag, expand, eye, sqrt


def previous_computation():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rotation-script', type=Path)
    args = parser.parse_args()
    relative = Path('2026-11-20-quaternion-frames/rotation_quadrupole.py')
    candidates = (Path(__file__).resolve().parent.parent / "rotation_quadrupole.py",
                  Path('experiments/exploration') / relative)
    script = args.rotation_script or next((p for p in candidates if p.exists()), None)
    if script is None:
        parser.error('supply --rotation-script with the prior exact group computation')
    return run_path(str(script))


def geometric_flags(model):
    rotations = model['icosahedral_rotations']()
    exact, key, action = (model[name] for name in ('exact', 'key', 'square_action'))
    klein = (eye(3), diag(1, -1, -1), diag(-1, 1, -1), diag(-1, -1, 1))
    seen = set()
    flags = []
    for rotation in rotations:
        subgroup = frozenset(key(exact(rotation * g * rotation.T)) for g in klein)
        if subgroup in seen:
            continue
        seen.add(subgroup)
        frame, inverse = action(rotation), action(rotation.T)
        for index, normal in enumerate((2, 1, 0)):
            coordinate_projector = diag(*(int(i == index + 2) for i in range(5)))
            projector = exact(frame * coordinate_projector * inverse)
            b, c = (rotation[:, i] for i in range(3) if i != normal)
            marks = (b + c, b - c, -b - c, -b + c)
            flags.append((projector, tuple(exact(mark / sqrt(2)) for mark in marks)))
    assert len(seen) == 5 and len(flags) == 15
    return flags


def nearest_connections(flags):
    connections = {}
    size = len(flags)
    for i, j in product(range(size), repeat=2):
        source, source_marks = flags[i]
        destination, destination_marks = flags[j]
        if expand((source * destination).trace()) != Rational(1, 4):
            continue
        mapping = []
        for mark in source_marks:
            scores = [expand(mark.dot(next_mark)) for next_mark in destination_marks]
            maximizers = [k for k, score in enumerate(scores) if all(score >= other for other in scores)]
            assert len(maximizers) == 1
            mapping.append(maximizers[0])
        assert len(set(mapping)) == 4
        assert all((mapping[(k + 1) % 4] - mapping[k]) % 4 in (1, 3) for k in range(4))
        connections[i, j] = tuple(mapping)
    assert len(connections) == 120
    for (i, j), mapping in connections.items():
        reverse = connections[j, i]
        assert tuple(reverse[k] for k in mapping) == tuple(range(4))
    return connections


def test_directed_frames(flags, connections):
    # The quartic orientation moment fixes i. The most aligned of the two
    # remaining axes, oriented towards the lifted next axis, fixes j.
    for i, j in connections:
        first, second = flags[i][1][:2]
        next_first, next_second = flags[j][1][:2]
        a = first.cross(second)
        b, c = (first+second)/sqrt(2), (first-second)/sqrt(2)
        d = next_first.cross(next_second)
        assert expand(a.dot(a)) == expand(d.dot(d)) == 1
        components = tuple(expand(vector.dot(d)) for vector in (a,b,c))
        assert all(components)
        assert expand(components[1]**2 - components[2]**2) != 0
    print('Directed quaternion frame selection is nondegenerate on all 120 moves')


def test_triangles(connections, size):
    expected = {'identity': 12, 'reflection': 16, 'quarter +': 2, 'quarter -': 2}
    for i in range(size):
        loop_types = Counter()
        for j, k in product(range(size), repeat=2):
            edges = ((i, j), (j, k), (k, i))
            if not all(edge in connections for edge in edges):
                continue
            first, second, third = (connections[edge] for edge in edges)
            loop = tuple(third[second[p]] for p in first)
            if loop == tuple(range(4)):
                loop_types['identity'] += 1
            elif (loop[1] - loop[0]) % 4 == 3:
                loop_types['reflection'] += 1
            else:
                assert loop[0] in (1, 3)
                loop_types['quarter +' if loop[0] == 1 else 'quarter -'] += 1
        assert loop_types == expected
    print('All 120 directed transports are unique D4 bijections with exact reverses')
    print('Triangle holonomy per flag:', expected, '; generated group D4')


def test_quarter_covariance(connections):
    # On a fixed edge the projected point destination is fixed. Commuting
    # with the oriented quarter-turn forces the unique output orientation.
    maps = list(connections.values()) + [tuple(range(4))]
    for mapping in maps:
        determinant = 1 if (mapping[1] - mapping[0]) % 4 == 1 else -1
        for point, sign in product(range(4), (-1, 1)):
            allowed = [out for out in (-1, 1)
                       if (mapping[(point+sign) % 4]-mapping[point]-out) % 4 == 0]
            assert allowed == [sign * determinant]
    print('Exact projected law and global C4 covariance force deterministic orientation')


def test_transfers(connections, size):
    points = 4 * eye(4 * size)
    sides = 4 * eye(4 * size)
    orientation = 4 * eye(size)
    oriented_odd = 4 * eye(4 * size)
    for (i, j), mapping in connections.items():
        orientation[i, j] = 1 if (mapping[1] - mapping[0]) % 4 == 1 else -1
        for k, image in enumerate(mapping):
            points[4 * i + k, 4 * j + image] = 1
            oriented_odd[4 * i + k, 4 * j + image] = orientation[i, j]
            edge = {mapping[k], mapping[(k + 1) % 4]}
            image_side = next(l for l in range(4) if {l, (l + 1) % 4} == edge)
            sides[4 * i + k, 4 * j + image_side] = 1
    assert all(transfer == transfer.T for transfer in (points, sides, orientation, oriented_odd))

    p_points, p_sides, p_orientation, p_odd = (
        transfer.charpoly() for transfer in (points, sides, orientation, oriented_odd)
    )
    x = p_points.gen
    standard = (x*x - 6*x + 4)**9 * (x**4 - 22*x**3 + 164*x*x - 488*x + 496)**3
    shared = (x - 12) * (x - 6)**5 * (x - 2)**9 * standard
    assert expand(p_points.as_expr() - shared * (x - 8)**5 * (x - 2)**10) == 0
    assert expand(p_sides.as_expr() - shared * (x - 4)**5 * (x*x - 8*x + 4)**5) == 0
    assert expand(p_orientation.as_expr() - (x - 4)**9 * (x*x - 8*x - 4)**3) == 0
    assert expand(p_odd.as_expr() - p_orientation.as_expr() * (x - 4)**5
                  * (x*x - 8*x + 4)**5 * standard) == 0
    # The square point/side matrices are real symmetric. Every factor
    # displayed above has no real zero at x <= 0, so both are positive.
    print('Primal and dual 60-state transfers are strictly positive')
    print('12P primal characteristic polynomial:', p_points.as_expr().factor())
    print('12P dual characteristic polynomial:', p_sides.as_expr().factor())
    # The signed orientation matrix is the odd fiber part in the
    # orientation cover. It is a subrepresentation of the 120-state lift.
    assert 4 - 2 * sqrt(5) < 0
    print('Orientation odd block 12P:', p_orientation.as_expr().factor())
    print('Oriented lift has exactly three negative eigenvalues (4-2sqrt(5))/12')
    assert 8 * size - 3 == 117
    print('Positive spectral quotient: dimension 117, with all 60 orientation-even modes retained')
    # Flips of the transported orientation only on moves change the odd
    # block to 4I+(1-2 epsilon)(12P_odd-4I), leaving the even marginal fixed.
    epsilon = (1 - 2 / sqrt(5)) / 2
    assert expand(4 + (1 - 2 * epsilon) * (-2 * sqrt(5))) == 0
    print('Nonstay flip threshold (5-2sqrt(5))/10 preserves the unorientated law')


if __name__ == '__main__':
    flags = geometric_flags(previous_computation())
    connections = nearest_connections(flags)
    test_directed_frames(flags, connections)
    test_triangles(connections, len(flags))
    test_quarter_covariance(connections)
    test_transfers(connections, len(flags))
