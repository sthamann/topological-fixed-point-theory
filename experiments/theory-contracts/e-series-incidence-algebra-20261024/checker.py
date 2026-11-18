#!/usr/bin/env python3
"""Exact E8 cut scheme and its two incidence transforms.

The root definition is the one supplied by the user. All geometry uses
integers in doubled coordinates. SymPy is used only for small exact Cartan
determinants and the five-dimensional adjacency algebra.
"""
import argparse
from collections import Counter, deque
from itertools import combinations, product
import json
from pathlib import Path

import sympy as sp


def need(value, label):
    if not value:
        raise RuntimeError(label)


def indices(mask):
    while mask:
        bit = mask & -mask
        yield bit.bit_length() - 1
        mask ^= bit


def bits(items):
    return sum(1 << i for i in items)


def orbit(seed, generators, action):
    seen, queue = {seed}, deque([seed])
    while queue:
        point = queue.popleft()
        for generator in generators:
            image = action(point, generator)
            if image not in seen:
                seen.add(image)
                queue.append(image)
    return sorted(seen)


def calculate():
    root_set = set()
    for i, j in combinations(range(8), 2):
        for a, b in product((-2, 2), repeat=2):
            v = [0] * 8
            v[i], v[j] = a, b
            root_set.add(tuple(v))
    root_set.update(s for s in product((-1, 1), repeat=8) if s.count(-1) % 2 == 0)
    roots = sorted(root_set)
    lookup = {v: i for i, v in enumerate(roots)}
    n = len(roots)
    need(n == 240, "root count")
    dots = [[sum(a*b for a, b in zip(x, y)) for y in roots] for x in roots]
    neg = [lookup[tuple(-a for a in v)] for v in roots]
    positive = {i for i, v in enumerate(roots) if next(a for a in v if a) > 0}

    def simple_ids(mask):
        pos = positive.intersection(indices(mask))
        decomposable = set()
        for i, j in combinations(sorted(pos), 2):
            k = lookup.get(tuple(a+b for a, b in zip(roots[i], roots[j])))
            if k in pos:
                decomposable.add(k)
        return sorted(pos - decomposable)

    def cartan(mask):
        simple = simple_ids(mask)
        matrix = sp.Matrix([[sp.Rational(dots[i][j], 4) for j in simple] for i in simple])
        need(all(x.q == 1 for x in matrix), "integral Cartan")
        return simple, matrix

    simple, gram = cartan((1 << n) - 1)
    need(len(simple) == 8 and gram.det() == 1, "actual E8 simple roots")
    permutations = []
    for a in range(n):
        images = []
        for i, v in enumerate(roots):
            need(dots[i][a] % 4 == 0, "root reflection integral coefficient")
            k = dots[i][a] // 4
            images.append(lookup[tuple(v[j] - k*roots[a][j] for j in range(8))])
        permutations.append(images)

    def image(mask, permutation):
        return bits(permutation[i] for i in indices(mask))

    weyl_generators = [permutations[i] for i in simple]
    a2 = set()
    for i, j in combinations(range(n), 2):
        if dots[i][j] == -4:
            k = lookup[tuple(a+b for a, b in zip(roots[i], roots[j]))]
            a2.add(bits((i, neg[i], j, neg[j], k, neg[k])))
    normals = sorted(a2)
    need(len(normals) == 1120, "all A2 normals")
    need(orbit(normals[0], weyl_generators, image) == normals, "transitive E8 action")
    normal_ids = {mask: i for i, mask in enumerate(normals)}
    cuts = [bits(k for k in range(n) if all(dots[k][j] == 0 for j in indices(a))) for a in normals]
    need(len(set(cuts)) == 1120 and {a.bit_count() for a in cuts} == {72}, "E6 complements")
    rel_names = ["identity", "D5", "D4", "A4", "A2+A2"]
    rel_sizes = [72, 40, 24, 20, 12]
    size_to_rel = {v: i for i, v in enumerate(rel_sizes)}
    adjacency = [[] for _ in rel_names]
    for x in cuts:
        row = [0] * 5
        for j, y in enumerate(cuts):
            rel = size_to_rel.get((x & y).bit_count())
            need(rel is not None, "complete pair intersection sizes")
            row[rel] |= 1 << j
        for i in range(5):
            adjacency[i].append(row[i])
    valencies = [a[0].bit_count() for a in adjacency]
    need(valencies == [1, 81, 270, 648, 120], "pair valencies")
    need(all([a[j].bit_count() for a in adjacency] == valencies for j in range(1120)), "all row valencies")

    # This explicitly generated subgroup fixes cut0 and its normal A2.
    stabilizer_roots = simple_ids(cuts[0]) + simple_ids(normals[0])
    cut_generators = []
    for a in stabilizer_roots:
        perm = [normal_ids[image(mask, permutations[a])] for mask in normals]
        need(perm[0] == 0, "stabilizer fixes cut0")
        cut_generators.append(perm)
    representatives, cartan_types = [], []
    for relation, row in enumerate(adjacency):
        members = set(indices(row[0]))
        representative = min(members)
        representatives.append(representative)
        need(set(orbit(representative, cut_generators, lambda i, p: p[i])) == members,
             "one stabilizer orbit for " + rel_names[relation])
        simple_cut, matrix = cartan(cuts[0] & cuts[representative])
        cartan_types.append({"relation": rel_names[relation], "simple_root_ids": simple_cut,
                             "cartan": [list(map(int, matrix.row(i))) for i in range(matrix.rows)],
                             "determinant": int(matrix.det())})
    need([(len(q["simple_root_ids"]), q["determinant"]) for q in cartan_types]
         == [(6, 3), (5, 4), (4, 4), (4, 5), (4, 9)], "actual Cartan intersection types")

    # Transitivity on cuts and each pair relation makes these ALL structure
    # constants, not counts for a suggestive selection of pairs.
    p = [[[0]*5 for _ in range(5)] for _ in range(5)]
    for k, representative in enumerate(representatives):
        for i in range(5):
            for j in range(5):
                p[k][i][j] = (adjacency[i][0] & adjacency[j][representative]).bit_count()
    left = [sp.Matrix(5, 5, lambda k, j: p[k][i][j]) for i in range(5)]
    need(all(a*b == b*a for a in left for b in left), "commutative adjacency algebra")
    for i in range(5):
        for j in range(5):
            need(left[i]*left[j] == sum((p[k][i][j]*left[k] for k in range(5)), sp.zeros(5)),
                 "complete multiplication law")
    separating = sum((10**i * left[i] for i in range(5)), sp.zeros(5))
    eigenspaces = separating.eigenvects()
    need(len(eigenspaces) == 5 and all(m == 1 for _, m, _ in eigenspaces), "five primitive summands")
    characters = []
    for _, _, vectors in eigenspaces:
        v = vectors[0]
        pivot = next(i for i in range(5) if v[i])
        values = [(a*v)[pivot]/v[pivot] for a in left]
        need(all(a*v == value*v for a, value in zip(left, values)), "joint eigencharacter")
        scale = sum(v[i]*values[i] for i in range(5))
        rank = 1120 * v[0] / scale
        need(rank.q == 1 and all(x.q == 1 for x in values), "integral multiplicities and eigenvalues")
        characters.append({"dimension": int(rank), "eigenvalues": list(map(int, values))})
    characters.sort(key=lambda q: q["dimension"])
    need([c["dimension"] for c in characters] == [1, 35, 84, 300, 700], "multiplicity-free dimensions")
    for character in characters:
        a1, a2 = character["eigenvalues"][1:3]
        kernel = int(a1 == -3)
        need(kernel*sp.Rational(a2+6, 24) == int(character["dimension"] == 300), "explicit P300")
        need(kernel*sp.Rational(18-a2, 24) == int(character["dimension"] == 700), "explicit P700")
    for i in range(5):
        for j in range(5):
            need(sum(c["dimension"]*c["eigenvalues"][i]*c["eigenvalues"][j] for c in characters)
                 == (1120*valencies[i] if i == j else 0), "exact character orthogonality")

    integer_mask = bits(i for i, v in enumerate(roots) if sum(x != 0 for x in v) == 2)
    d8_frames = orbit(integer_mask, weyl_generators, image)
    need(len(d8_frames) == 135, "actual D8 frame orbit")
    seed_a3 = bits(i for i in indices(integer_mask) if all(roots[i][j] == 0 for j in range(5)))
    tetra_normals = orbit(seed_a3, weyl_generators, image)
    need(len(tetra_normals) == 7560, "actual A3 frame orbit")
    tetra = [tuple(i for i, a in enumerate(normals) if a & normal == a) for normal in tetra_normals]
    need({len(t) for t in tetra} == {4}, "four E6 cuts per tetrahedron")
    tetra_rows = [0] * 1120
    pair_tetra = Counter()
    cores = []
    for j, t in enumerate(tetra):
        core = (1 << n) - 1
        for i in t:
            tetra_rows[i] |= 1 << j
            core &= cuts[i]
        need(core.bit_count() == 40, "D5 core of every tetrahedron")
        cores.append(core)
        pair_tetra.update(combinations(t, 2))
    need(len(set(cores)) == 7560, "D5 core uniquely determines tetrahedron")
    need({r.bit_count() for r in tetra_rows} == {27}, "all 27-fold tetra incidences")
    need(len(pair_tetra) == 45360 and set(pair_tetra.values()) == {1}, "unique tetrahedron per D5 edge")
    frame_tetra = []
    tetra_frame_count = Counter()
    raw_a3_frame_counts = []
    for d8 in d8_frames:
        # Mere A3 containment is insufficient: D8 also contains 560 A3's
        # of the four-coordinate chain type. A compiler frame must contain
        # BOTH the A3 normal and the common D5 core of its four E6 cuts.
        raw_a3_frame_counts.append(sum(a & d8 == a for a in tetra_normals))
        members = [i for i, a in enumerate(tetra_normals)
                   if (a | cores[i]) & d8 == (a | cores[i])]
        need(len(members) == 56, "56 tetrahedra in every D8 frame")
        frame_tetra.append(members)
        tetra_frame_count.update(members)
    need(len(tetra_frame_count) == 7560 and set(tetra_frame_count.values()) == {1},
         "unique D8 frame for every tetrahedron")
    need(set(raw_a3_frame_counts) == {616}, "56 compiler and 560 other A3 subsystems inside D8")

    lines = [bits((i, neg[i])) for i in range(n) if i < neg[i]]
    e7_cuts = [bits(k for k in range(n) if dots[k][next(indices(line))] == 0) for line in lines]
    need(len(lines) == 120 and {c.bit_count() for c in e7_cuts} == {126}, "120 actual E7 cuts")
    n_rows = [bits(j for j, line in enumerate(lines) if line & a == line) for a in normals]
    n_columns = [bits(i for i, row in enumerate(n_rows) if row >> j & 1) for j in range(120)]
    need({row.bit_count() for row in n_rows} == {3}, "three E7 per E6")
    need({col.bit_count() for col in n_columns} == {28}, "28 E6 per E7")
    for i in range(1120):
        for j in range(i, 1120):
            d5 = (adjacency[1][i] >> j) & 1
            need((n_rows[i] & n_rows[j]).bit_count() == (3 if i == j else d5), "NNt identity")
            need((tetra_rows[i] & tetra_rows[j]).bit_count() == (27 if i == j else d5), "MMt identity")
            if d5:
                need(pair_tetra[(i, j)] == 1, "same D5 edge in both incidence geometries")
    b_rows = []
    for i in range(120):
        row = 0
        for j in range(120):
            count = (n_columns[i] & n_columns[j]).bit_count()
            need(count == 28 if i == j else count in (0, 1), "NtN exact entries")
            if i != j and count:
                row |= 1 << j
        b_rows.append(row)
    need({row.bit_count() for row in b_rows} == {56}, "E7 graph degree")
    for i, j in combinations(range(120), 2):
        need((b_rows[i] & b_rows[j]).bit_count() == (28 if b_rows[i] >> j & 1 else 24),
             "strongly regular graph common neighbours")
    d5_spectrum = Counter()
    for c in characters:
        d5_spectrum[c["eigenvalues"][1]] += c["dimension"]
    need(d5_spectrum == {81: 1, 33: 35, 21: 84, -3: 1000}, "D5 graph spectrum")
    need(sum(m for a, m in d5_spectrum.items() if a+3) == 120, "rank N")
    need(all(a+27 > 0 for a in d5_spectrum), "tetra transform full row rank")
    x = sp.Symbol("x")
    inverse_poly = sp.interpolate([(a, sp.Rational(1, a+27)) for a in d5_spectrum], x).expand()
    need(all(inverse_poly.subs(x, a)*(a+27) == 1 for a in d5_spectrum), "exact reconstruction polynomial")
    coefficients = [str(inverse_poly.coeff(x, i)) for i in range(4)]
    results = {
        "research_id": "MATH.SEAM.ESERIES.INCIDENCE.07", "verdict": "PASS",
        "counts": {"roots": 240, "E7": 120, "E6": 1120, "D8_frames": 135, "tetrahedra": 7560,
                   "tetra_E6_flags": 30240, "E6_E7_flags": 3360, "D5_edges": 45360},
        "relation_order": rel_names, "valencies": valencies,
        "cartan_certificates": cartan_types,
        "structure_constants_p_k_i_j": p, "characters": characters,
        "left_multiplication_matrices": [[list(map(int, a.row(i))) for i in range(5)] for a in left],
        "E7_graph_parameters": [120, 56, 28, 24],
        "A3_in_D8": {"all_A3_subsystems": 616, "compiler_D5_A3_frames": 56,
                      "raw_D8_containments_per_A3": 11, "full_flag_D8_containments": 1},
        "D5_spectrum": {str(a): m for a, m in sorted(d5_spectrum.items())},
        "E6_E7_incidence": {"gram": "N Nt=A_D5+3I", "rank": 120, "kernel_dimension": 1000},
        "E6_tetra_incidence": {"gram": "M Mt=A_D5+27I", "rank": 1120,
                              "gram_eigenvalues": {str(a+27): m for a, m in sorted(d5_spectrum.items())},
                              "inverse_gram_polynomial_coefficients_ascending": coefficients},
        "scope": "Exact finite geometry and actual Weyl transport on standard E8 roots. No source selection or physical inference.",
    }
    data = {"coordinate_scale": 2, "roots": roots, "E6_normal_root_ids": [list(indices(a)) for a in normals],
            "E6_root_ids": [list(indices(a)) for a in cuts], "E7_normal_line_ids": [list(indices(a)) for a in lines],
            "E7_root_ids": [list(indices(a)) for a in e7_cuts], "D8_root_ids": [list(indices(a)) for a in d8_frames],
            "tetra_normal_A3_root_ids": [list(indices(a)) for a in tetra_normals], "tetra_E6_ids": tetra,
            "D8_tetra_ids": frame_tetra, "E6_E7_ids": [list(indices(a)) for a in n_rows],
            "relation_bit_rows_hex": [[hex(a) for a in row] for row in adjacency]}
    return results, data


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    result, data = calculate()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for name, value in (("results.json", result), ("geometry.json", data)):
        (args.output_dir/name).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    print(json.dumps({k: result[k] for k in ("research_id", "verdict", "counts", "valencies", "characters",
                                             "E6_E7_incidence", "E6_tetra_incidence")}, indent=2))
