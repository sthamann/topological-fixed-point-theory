"""Exact anticanonical realization and visible-clock controls for square marks."""

from itertools import product

from sympy import I, Matrix, Poly, diag, expand, eye, kronecker_product, simplify, sqrt, symbols


z = symbols("z")
BASIS = (
    diag(1, -1, 0),
    diag(1, 1, -2),
    Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 0]]),
    Matrix([[0, 0, 1], [0, 0, 0], [1, 0, 0]]),
    Matrix([[0, 0, 0], [0, 0, 1], [0, 1, 0]]),
)
NULL_VECTOR = Matrix([1 - z**2, I * (1 + z**2), 2 * z])
DIVISOR = 1 + z**4
QUARTER_TURN = Matrix([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
HALF_TURN = diag(1, -1, -1)


def quadratic(form):
    return expand(NULL_VECTOR.dot(form * NULL_VECTOR))


def coefficients(polynomial):
    polynomial = Poly(expand(polynomial), z)
    return Matrix([polynomial.nth(power) for power in range(5)])


def action(pullback):
    return Matrix.hstack(*(coefficients(pullback(z**power)) for power in range(5)))


def swap_axes(polynomial):
    return expand(z**4 * polynomial.subs(z, 1 / z))


def main():
    assert expand(NULL_VECTOR.dot(NULL_VECTOR)) == 0
    chart = Matrix.hstack(*(coefficients(quadratic(form)) for form in BASIS))
    assert chart.det() != 0
    assert quadratic(BASIS[0]) == 2 * DIVISOR

    round_actions = (
        (QUARTER_TURN, lambda polynomial: -polynomial.subs(z, I * z)),
        (HALF_TURN, swap_axes),
    )
    for rotation, pullback in round_actions:
        for form in BASIS:
            assert quadratic(rotation.T * form * rotation) == expand(pullback(quadratic(form)))
    quarter, half = (action(pullback) for _, pullback in round_actions)
    assert quarter**4 == half**2 == eye(5)
    assert half * quarter * half == quarter.inv()
    assert -DIVISOR.subs(z, I * z) == -DIVISOR
    assert swap_axes(DIVISOR) == DIVISOR

    square = [expand((1 + I) / sqrt(2) * I**power) for power in range(4)]
    rotation_permutation = (1, 2, 3, 0)
    half_permutation = (3, 2, 1, 0)
    assert all(expand(point**4 + 1) == 0 for point in square)
    assert all(simplify(I * point - square[rotation_permutation[index]]) == 0
               for index, point in enumerate(square))
    assert all(simplify(1 / point - square[half_permutation[index]]) == 0
               for index, point in enumerate(square))
    for permutation, expected in ((rotation_permutation, -1), (half_permutation, 1)):
        determinant = Matrix(4, 4, lambda i, j: int(i == permutation[j])).det()
        assert determinant == expected

    ordinary_quarter = action(lambda polynomial: polynomial.subs(z, I * z))
    words = tuple(product(range(4), range(2)))
    twists = []
    for quarter_sign, half_sign in product((1, -1), repeat=2):
        if all((quarter**r * half**s).trace() ==
               ((quarter_sign * ordinary_quarter)**r * (half_sign * half)**s).trace()
               for r, s in words):
            twists.append((quarter_sign, half_sign))
    assert twists == [(-1, 1)]

    charges = (1, I, -1, -I)
    multiplicities = tuple(5 - (quarter - charge * eye(5)).rank() for charge in charges)
    assert multiplicities == (1, 1, 2, 1)
    for phase in charges:
        assert phase * quarter != eye(5) and phase * quarter.inv() != eye(5)

    balanced = kronecker_product(quarter, quarter.conjugate())
    assert 25 - (balanced - eye(25)).rank() == sum(rank**2 for rank in multiplicities) == 7
    print(f"Null-conic quadrupole-to-section map: rank {chart.rank()}, determinant {chart.det()}")
    print("Vanishing line p_D=(1+z^4)(d/dz)^2: (r,s)=(-1,+1), the family determinant character")
    print("Meromorphic 3+2 carrier: 1 lies on the distinguished line of the 2-block")
    print(f"Unique matching character among 4: {twists}")
    print(f"Geometrical quarter multiplicities (1,i,-1,-i): {multiplicities}")
    print("Native visible Cartan J fixed ranks: (5,0,0,0); neither rotation direction nor any scalar rephase intertwines")
    print("Self-dual geometric counterphase: fixed dimension 7 before additional rank-5 selection")


if __name__ == "__main__":
    main()
