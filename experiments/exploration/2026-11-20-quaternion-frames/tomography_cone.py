"""Exact controls for a local cone selected by tomographic composition.

Assume homogeneous Hermitian matrix families over R, C, and H with
multiplied matrix size in composites. These are the associative real
normed division algebras. This is a declared model class, not a census
of all operational theories or all possible composition rules.

An ordered state cone is available within the model. Interpreting
positive differences as physical causal intervals still requires an event
realization. The cone supplies neither gluing of local charts nor a clock,
absolute units, a physical state, curvature, or gravity.
"""

from sympy import I, Matrix, Rational, eye, expand, kronecker_product, trace
from sympy.abc import a, m, n, t, x, y, z


def dimension(capacity, scalar_dimension):
    """Count one real diagonal component and a scalar per off-diagonal pair."""
    return capacity + scalar_dimension * capacity * (capacity - 1) / 2


def test_tomographic_dimension():
    defect = expand(dimension(n * m, a) - dimension(n, a) * dimension(m, a))
    formula = a * (2 - a) * n * m * (n - 1) * (m - 1) / 4
    assert expand(defect - formula) == 0
    assert dimension(1, a) == 1
    survivors = []
    for name, scalar_dimension in (("R", 1), ("C", 2), ("H", 4)):
        local = dimension(2, Rational(scalar_dimension))
        joint = dimension(4, Rational(scalar_dimension))
        product = local**2
        if joint == product:
            survivors.append(name)
        print(
            f"{name}: local={local}, joint={joint}, product={product}, "
            f"defect={joint-product}, cone=1+{local-1}"
        )
    assert survivors == ["C"]
    print("Exact tomography: 1 of 3 matrix families")


def test_rebit_control():
    identity = eye(2)
    sx = Matrix([[0, 1], [1, 0]])
    sy = Matrix([[0, -I], [I, 0]])
    sz = Matrix([[1, 0], [0, -1]])
    witness = kronecker_product(sy, sy)
    rho_plus = (eye(4) + witness) / 4
    rho_minus = (eye(4) - witness) / 4
    for state in (rho_plus, rho_minus):
        assert state.is_symmetric() and state.conjugate() == state
        assert trace(state) == 1
        assert state.eigenvals() == {0: 2, Rational(1, 2): 2}
    for left in (identity, sx, sz):
        for right in (identity, sx, sz):
            observable = kronecker_product(left, right)
            assert trace((rho_plus - rho_minus) * observable) == 0
    assert trace(rho_plus * witness) == 1
    assert trace(rho_minus * witness) == -1
    print("Rebits: all 9 product probes agree, joint witness = +1 and -1")


def test_cone():
    increment = Matrix([[t + z, x - I * y], [x + I * y, t - z]])
    determinant = t**2 - x**2 - y**2 - z**2
    assert expand(increment.det() - determinant) == 0
    assert trace(increment) == 2 * t
    assert (increment - t * eye(2))**2 == (x**2 + y**2 + z**2) * eye(2)
    print(f"Complex determinant: {determinant}")
    print("PSD iff t >= sqrt(x^2 + y^2 + z^2), by trace and determinant")

    boost = Matrix.diag(2, Rational(1, 2))
    boosted = boost * increment * boost.H
    assert boost.det() == 1 and boosted.det() == increment.det()
    new_t = expand(trace(boosted) / 2)
    new_z = expand((boosted[0, 0] - boosted[1, 1]) / 2)
    assert new_t == (17 * t + 15 * z).expand() / 8
    assert new_z == (15 * t + 17 * z).expand() / 8
    assert boosted[1, 0] == increment[1, 0]
    print(f"Rational boost: t'={new_t}, z'={new_z}, x'=x, y'=y")
    print("Invertible congruence preserves the PSD cone")

    assert (a * increment).det().expand() == (a**2 * determinant).expand()
    assert increment.T == increment.subs(y, -y)
    print("Scale counterexample: a>0 preserves PSD order, det scales by a^2")
    print("Orientation counterexample: transpose preserves PSD and flips y")


if __name__ == "__main__":
    test_tomographic_dimension()
    test_rebit_control()
    test_cone()
