#!/usr/bin/env python3
"""Exact central-lift census for the nine source-compatible Pauli Q actions.

The standard two-qubit Clifford generators and v783's frozen Pauli chart are
used.  By v783 P3.3/P3.7, the Q(i)-rational half of this matrix group is G31,
its kernel over Sp(4,2) is mu4 times Pauli16, and the other half is zeta8 G31.
"""

from collections import deque
from itertools import combinations, product

import sympy as sp


def check(label, condition):
    if not condition:
        raise AssertionError(label)
    print("PASS", label)


def f2_mul(a, b):
    return tuple(tuple(sum(a[i][k] * b[k][j] for k in range(4)) % 2
                       for j in range(4)) for i in range(4))


def f2_apply(a, v):
    return tuple(sum(a[i][j] * v[j] for j in range(4)) % 2
                 for i in range(4))


def f2_inverse(a):
    aug = [list(a[i]) + [int(i == j) for j in range(4)] for i in range(4)]
    for col in range(4):
        pivot = next(row for row in range(col, 4) if aug[row][col])
        aug[col], aug[pivot] = aug[pivot], aug[col]
        for row in range(4):
            if row != col and aug[row][col]:
                aug[row] = [x ^ y for x, y in zip(aug[row], aug[col])]
    return tuple(tuple(row[4:]) for row in aug)


def from_columns(*columns):
    return tuple(tuple(column[i] for column in columns) for i in range(4))


def symplectic(v, w):
    return (v[0] * w[1] + v[1] * w[0] +
            v[2] * w[3] + v[3] * w[2]) % 2


def canonical(expr):
    return sp.expand_complex(sp.radsimp(sp.simplify(expr)))


def is_qi_rational(matrix):
    return all(not canonical(x).has(sp.sqrt(2)) for x in matrix)


def scalar_square(matrix):
    square = matrix * matrix
    value = canonical(square[0, 0])
    if all(canonical(square[r, c] - (value if r == c else 0)) == 0
           for r in range(4) for c in range(4)):
        return value
    return None


def scalar_label(value):
    if value is None:
        return "nonscalar"
    for label, target in (("+1", 1), ("-1", -1), ("+i", sp.I), ("-i", -sp.I)):
        if canonical(value - target) == 0:
            return label
    return str(canonical(value))


def matrix_equal(a, b):
    return all(canonical(a[r, c] - b[r, c]) == 0
               for r in range(a.rows) for c in range(a.cols))


def main():
    identity = tuple(tuple(int(i == j) for j in range(4)) for i in range(4))
    omega = ((0, 1, 0, 0), (1, 0, 0, 0),
             (0, 0, 0, 1), (0, 0, 1, 0))
    sigma_pauli = ((0, 0, 1, 0), (0, 1, 0, 1),
                   (1, 0, 1, 0), (0, 1, 0, 0))
    sigma_pauli_2 = f2_mul(sigma_pauli, sigma_pauli)

    bits = tuple(product((0, 1), repeat=4))
    nonzero = tuple(v for v in bits if any(v))
    planes = set()
    for a, b in combinations(nonzero, 2):
        planes.add(frozenset(((0, 0, 0, 0), a, b,
                              tuple(x ^ y for x, y in zip(a, b)))))
    stable_symplectic_planes = []
    for plane in planes:
        nz = tuple(v for v in plane if any(v))
        if ({f2_apply(sigma_pauli, v) for v in plane} == set(plane)
                and any(symplectic(a, b) for a, b in combinations(nz, 2))):
            stable_symplectic_planes.append(plane)
    check("two sigma-stable nondegenerate planes",
          len(stable_symplectic_planes) == 2)

    sp4 = []
    for entries in product((0, 1), repeat=16):
        matrix = tuple(tuple(entries[4 * r + c] for c in range(4))
                       for r in range(4))
        try:
            f2_inverse(matrix)
        except StopIteration:
            continue
        if f2_mul(f2_mul(tuple(zip(*matrix)), omega), matrix) == omega:
            sp4.append(matrix)
    check("full Sp(4,2) enumeration", len(sp4) == 720)

    q_candidates = [
        matrix for matrix in sp4
        if f2_mul(matrix, matrix) == identity
        and f2_mul(f2_mul(matrix, sigma_pauli), matrix) == sigma_pauli_2
        and all({f2_apply(matrix, v) for v in plane} == set(plane)
                for plane in stable_symplectic_planes)
    ]
    check("nine plane-preserving reversing Q actions", len(q_candidates) == 9)
    sigma_centralizer_sp = [matrix for matrix in sp4
                            if f2_mul(matrix, sigma_pauli) ==
                            f2_mul(sigma_pauli, matrix)]
    check("symplectic centralizer of sigma has order 18",
          len(sigma_centralizer_sp) == 18)
    q0 = q_candidates[0]
    q_orbit = {
        f2_mul(f2_mul(g, q0), f2_inverse(g))
        for g in sigma_centralizer_sp
    }
    check("the sigma centralizer is transitive on all nine Q actions",
          q_orbit == set(q_candidates))

    ii = sp.I
    sqrt2 = sp.sqrt(2)
    identity2 = sp.eye(2)
    identity4 = sp.eye(4)
    x = sp.Matrix(((0, 1), (1, 0)))
    z = sp.diag(1, -1)
    y = sp.Matrix(((0, -ii), (ii, 0)))
    h = sp.Matrix(((1, 1), (1, -1))) / sqrt2
    phase = sp.diag(1, ii)
    cnot = sp.Matrix(((1, 0, 0, 0), (0, 1, 0, 0),
                      (0, 0, 0, 1), (0, 0, 1, 0)))
    gates = (sp.kronecker_product(h, identity2),
             sp.kronecker_product(identity2, h),
             sp.kronecker_product(phase, identity2),
             sp.kronecker_product(identity2, phase), cnot)

    # Exact symplectic images in v783 bit order (x1,z1,x2,z2).
    gate_actions = (
        from_columns((0, 1, 0, 0), (1, 0, 0, 0),
                     (0, 0, 1, 0), (0, 0, 0, 1)),
        from_columns((1, 0, 0, 0), (0, 1, 0, 0),
                     (0, 0, 0, 1), (0, 0, 1, 0)),
        from_columns((1, 1, 0, 0), (0, 1, 0, 0),
                     (0, 0, 1, 0), (0, 0, 0, 1)),
        from_columns((1, 0, 0, 0), (0, 1, 0, 0),
                     (0, 0, 1, 1), (0, 0, 0, 1)),
        from_columns((1, 0, 1, 0), (0, 1, 0, 0),
                     (0, 0, 1, 0), (0, 1, 0, 1)),
    )
    words = {identity: ()}
    frontier = deque((identity,))
    while frontier:
        current = frontier.popleft()
        for index, generator in enumerate(gate_actions):
            image = f2_mul(current, generator)
            if image not in words:
                words[image] = words[current] + (index,)
                frontier.append(image)
    check("standard Clifford generators cover Sp(4,2)", len(words) == 720)

    one_qubit = {(0, 0): identity2, (1, 0): x,
                 (0, 1): z, (1, 1): y}
    paulis = {b: sp.kronecker_product(one_qubit[b[:2]], one_qubit[b[2:]])
              for b in bits}
    mu4 = (1, ii, -1, -ii)
    zeta8 = (1 + ii) / sqrt2

    expected_g31 = {"+1": 8, "-1": 8, "nonscalar": 48}
    expected_coset = {"+i": 8, "-i": 8, "nonscalar": 48}
    all_censuses = []
    joint_survivor_counts = []
    joint_survivors_by_q = {}
    representatives = {}
    frozen_sigma = sp.Matrix(((0, 0, 1, 0), (1, 0, 0, 0),
                              (0, 1, 0, 0), (0, 0, 0, 1)))
    frozen_sigma_inverse = frozen_sigma ** 2
    for q_action in q_candidates:
        representative = identity4
        for index in words[q_action]:
            representative = representative * gates[index]
        representative = representative.applyfunc(canonical)
        # For these nine actions the BFS representatives themselves lie in
        # Q(i), hence in G31 by v783's exact index-two rationality theorem.
        check("Q representative lies in the Q(i)-rational G31 half",
              is_qi_rational(representative))
        representatives[q_action] = representative

        g31_census = {}
        coset_census = {}
        joint_survivors = []
        for scalar in mu4:
            for pauli_bits, pauli in paulis.items():
                lift = scalar * pauli * representative
                label = scalar_label(scalar_square(lift))
                g31_census[label] = g31_census.get(label, 0) + 1
                coset_lift = zeta8 * lift
                label = scalar_label(scalar_square(coset_lift))
                coset_census[label] = coset_census.get(label, 0) + 1
                if (label == "+i" and
                        matrix_equal(coset_lift * frozen_sigma *
                                     coset_lift.conjugate().T,
                                     frozen_sigma_inverse)):
                    joint_survivors.append((scalar, pauli_bits))
        all_censuses.append((g31_census, coset_census))
        joint_survivor_counts.append((len(joint_survivors), joint_survivors))
        joint_survivors_by_q[q_action] = joint_survivors
    check("all nine G31 fibers have only central squares plus or minus one",
          all(g31 == expected_g31 for g31, _ in all_censuses))
    check("all nine zeta8 fibers have central squares plus or minus i",
          all(coset == expected_coset for _, coset in all_censuses))
    check("full alpha sigma alpha^-1=sigma^-1 relation has two lifts per Q",
          all(count == 2 for count, _ in joint_survivor_counts))

    q_displayed = ((1, 0, 1, 0), (0, 1, 1, 0),
                   (0, 0, 1, 0), (1, 1, 1, 1))
    u_displayed = representatives[q_displayed]
    expected_u = sp.Matrix(((1, 0, 0, 0), (0, 0, 0, -1),
                            (0, 0, 1, 0), (0, 1, 0, 0)))
    check("displayed Q lift U is the stated rational signed permutation",
          u_displayed == expected_u)
    check("displayed Q has only plus or minus the explicit XI lift",
          joint_survivors_by_q[q_displayed] ==
          [(1, (1, 0, 0, 0)), (-1, (1, 0, 0, 0))])
    alpha = zeta8 * sp.kronecker_product(x, identity2) * u_displayed
    check("explicit alpha has the imposed compatible native clock relations",
          scalar_label(scalar_square(alpha)) == "+i" and
          (alpha ** 4).applyfunc(canonical) == -identity4 and
          matrix_equal(alpha * frozen_sigma * alpha.conjugate().T,
                       frozen_sigma_inverse))

    print("uniform G31 census:", expected_g31)
    print("uniform zeta8 G31 census:", expected_coset)
    print("joint survivors per Q, before alpha~-alpha:",
          [count for count, _ in joint_survivor_counts])
    print("displayed U =", u_displayed.tolist())
    print("explicit alpha = zeta8 * XI * U; alpha^2=iI, alpha^4=-I, "
          "alpha sigma alpha^-1=sigma^-1")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
