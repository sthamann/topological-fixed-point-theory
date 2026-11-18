"""Exact certificates for the D4 history-channel family and a word test."""

from itertools import permutations

from sympy import Integer, Matrix, Rational, eye, kronecker_product, symbols, zeros


def channel(cycles):
    result = zeros(16)
    for cycle in cycles:
        permutation = eye(4)[:, list(cycle)]
        result += kronecker_product(permutation, permutation)
    return result / len(cycles)


def spectral_projector(operator, eigenvalue, spectrum):
    identity = eye(operator.rows)
    result = identity
    for other in spectrum:
        if other != eigenvalue:
            result = result * (operator - other * identity) / (eigenvalue - other)
    return result


def choi(operator):
    return Matrix(
        16, 16, lambda row, col:
        operator[4 * (row % 4) + col % 4, 4 * (row // 4) + col // 4]
    )


def audit():
    rotation, inverse = (1, 2, 3, 0), (3, 0, 1, 2)
    cycles = [
        p for p in permutations(range(4))
        if all(len({i, p[i], p[p[i]], p[p[p[i]]]}) == 4 for i in range(4))
    ]
    assert len(cycles) == 6
    planar = channel([rotation, inverse])
    crossed = channel([p for p in cycles if p not in (rotation, inverse)])
    assert planar == planar.T and crossed == crossed.T and planar * crossed == crossed * planar

    half = Rational(1, 2)
    a = {v: spectral_projector(planar, v, (1, 0, -1)) for v in (1, 0, -1)}
    b = {
        v: spectral_projector(crossed, v, (1, 0, half, -half))
        for v in (1, 0, half, -half)
    }
    fixed = a[1] * b[1]
    planar_modes = (a[1] + a[-1]) * b[0]
    difference_modes = a[1] * b[-half] + a[-1] * b[half]
    crossed_modes = a[0]
    projectors = (fixed, planar_modes, difference_modes, crossed_modes)
    assert sum(projectors, zeros(16)) == eye(16)
    assert all(p == p.T and p * p == p for p in projectors)
    assert [p.trace() for p in projectors] == [2, 4, 2, 8]

    x = symbols("x", real=True)
    mixture = x * planar + (1 - x) * crossed
    factors = (Integer(1), x**2, ((1 - 3*x)/2)**2, ((1 - x)/2)**2)
    assert (mixture**2 - sum((f*p for f, p in zip(factors, projectors)), zeros(16))).expand() == zeros(16)

    identity_vector = eye(4).reshape(16, 1)
    witnesses = (
        Matrix([0, 1, 0, -1, 1, 0, -1, 0, 0, -1, 0, 1, -1, 0, 1, 0]),
        Matrix([0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0]),
    )
    coeffs = []
    for witness in witnesses:
        assert witness.dot(identity_vector) == 0
        coeffs.append(tuple(
            witness.dot(choi(p) * witness) / witness.dot(witness) for p in projectors[1:]
        ))
    assert coeffs == [(-1, Rational(2, 3), 0), (1, half, -2)]

    # For x > 1/4, the first Rayleigh quotient is negative iff
    # (1-3x)^2 < 4x^3. The exceptional singular x are checked separately.
    assert ((1 - 3*x)**2 - 4*x**3 + (1 - x)**2 * (4*x - 1)).expand() == 0
    # For 0 < x <= 1/4, the other quotient is strictly negative:
    # x^2(1-3x) <= 4/243 by differentiation, while 1-x >= 3/4.
    assert (x*x*(1 - 3*x)).diff(x).expand() == -9*x*x + 2*x
    assert 8 * Rational(4, 243) < Rational(3, 4)**4
    assert all(any(f.subs(x, weight) == 0 for f in factors)
               for weight in (0, Rational(1, 3), 1))

    # A genuine CP semigroup supplies a nonnegative conditional Choi control.
    depolarizing_generator = identity_vector * identity_vector.T / 4 - eye(16)
    complement = eye(16) - identity_vector * identity_vector.T / 4
    assert complement * choi(depolarizing_generator) * complement == complement / 4

    # A mixed compression defect can vanish even when both operators leak.
    v = eye(3)[:, 0]
    u, y = eye(3)[:, [1, 0, 2]], eye(3)[:, [2, 1, 0]]
    outside = eye(3) - v*v.T
    assert v.T * u * outside * y * v == zeros(1)
    assert (outside * u * v).norm() == (outside * y * v).norm() == 1

    print("D4 quantum paired channel excluded for all x in [0,1] under tracial symmetric CP time.")
    print("Separate controls passed: depolarizing CP generator and vanishing mixed defect with leakage.")


if __name__ == "__main__":
    audit()
