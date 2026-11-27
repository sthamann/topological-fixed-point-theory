"""Wave propagation, reflection positivity and transversal wave probes."""

import sympy as sp


def causal_coefficient():
    # At (t,0) integrate two delta shells G=delta(t-r)/(4 pi r).
    # First integrate intermediate time; the remaining delta(t-2r) gives 1/2.
    solid_angle = 4 * sp.pi
    jacobian = sp.S.Half
    coefficient = solid_angle * jacobian / (4 * sp.pi) ** 2
    assert coefficient == 1 / (8 * sp.pi)
    # Independent massive tail: -m J_1(m tau)/(4 pi tau).
    m, tau = sp.symbols("m tau", positive=True)
    expansion = sp.series(-m * sp.besselj(1, m * tau) / (4 * sp.pi * tau), m, 0, 4)
    assert expansion.removeO().expand().coeff(m, 2) == -coefficient
    # Ignore separable positive factors e^{-omega(t_i+t_j)}/(4 omega^3).
    # The remaining matrix is 1+omega(t_i+t_j) for the double pole.
    a, b, w = sp.symbols("a b w", positive=True)
    hankel = sp.Matrix([[1 + w * (x + y) for y in (a, b)] for x in (a, b)])
    assert sp.factor(hankel.det()) == -(w**2) * (a - b) ** 2
    simple = sp.Matrix(
        [[sp.exp(-w * (x + y)) / (2 * w) for y in (a, b)] for x in (a, b)]
    )
    assert sp.simplify(simple.det()) == 0
    print(
        "Retarded double propagation 1/(8π), also via the mass derivative. Double pole fails reflection positivity; simple pole has a positive rank-one Gram matrix."
    )


def null_directions():
    # Mobius-normalized quadruple (0, infinity, 1, x+iy).
    x, y = sp.symbols("x y", real=True)
    rows = [
        [1, 0, 0, 1],
        [1, 0, 0, -1],
        [2, 2, 0, 0],
        [1 + x * x + y * y, 2 * x, 2 * y, 1 - x * x - y * y],
    ]
    rays = sp.Matrix(rows)
    metric = sp.diag(1, -1, -1, -1)
    assert all(
        sp.simplify((rays[i, :] * metric * rays[i, :].T)[0]) == 0 for i in range(4)
    )
    assert sp.expand(rays.det()) == -8 * y
    assert rays.subs({x: -1, y: 0}).rank() == 3 and rays.subs({x: 0, y: 1}).rank() == 4
    print(
        "Four null probes: determinant -8 Im(z). Any harmonic quartet has rank at most three; a nonreal cross ratio gives four transverse directions."
    )


if __name__ == "__main__":
    causal_coefficient()
    null_directions()
