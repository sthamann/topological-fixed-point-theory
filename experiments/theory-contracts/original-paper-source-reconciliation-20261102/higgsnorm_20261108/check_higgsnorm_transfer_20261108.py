"""Check the radial contact and the actual v330 generator without a search."""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import sympy as sp


def main(repo):
    sys.path.insert(0, str(repo / "verification"))
    from v330_qgamb_admissible_measure import correlator, symmetric_transfer

    x = sp.symbols("x1:5", real=True)
    a, s = sp.symbols("a s", positive=True)
    u = sum(coordinate**2 for coordinate in x)
    psi = sp.Function("psi")(*x)
    mark = 1 - sp.exp(-s * u / 2)

    def kinetic(f):
        return -a * sum(sp.diff(f, coordinate, 2) for coordinate in x)

    def double(f, g):
        return f * (kinetic(g * psi) - g * kinetic(psi)) - (
            kinetic(g * f * psi) - g * kinetic(f * psi)
        )

    assert sp.simplify(double(mark, mark) - 2 * a * s**2 * u * sp.exp(-s * u) * psi) == 0
    assert sp.simplify(double(1 - sp.exp(-s * u), x[0]) - 4 * a * s * sp.exp(-s * u) * x[0] * psi) == 0

    transfer, stationary, lam2, lam3 = symmetric_transfer()
    p0 = sp.ones(3, 3) / 3
    v2, v3 = sp.Matrix([1, -1, 0]), sp.Matrix([1, 1, -2])
    p2, p3 = v2 * v2.T / 2, v3 * v3.T / 6
    assert p0 + p2 + p3 == sp.eye(3)
    assert all(p * p == p for p in (p0, p2, p3))
    generator = 6 * sp.log(sp.Rational(3, 2)) * p2 + 6 * sp.log(3) * p3
    expected = p0 + sp.Rational(64, 729) * p2 + sp.Rational(1, 729) * p3
    assert np.allclose(transfer, np.array(expected, dtype=float), rtol=0, atol=1e-14)
    assert np.allclose([lam2, lam3], [64 / 729, 1 / 729], rtol=0, atol=1e-14)
    assert generator * sp.ones(3, 1) == sp.zeros(3, 1)
    assert sp.simplify(generator[0, 1] - sp.log(sp.Rational(8, 9))) == 0
    assert sp.simplify(generator[0, 2] + 2 * sp.log(3)) == 0
    diagonal = sp.diag(1, -1, 0)
    centered = diagonal * sp.ones(3, 1)
    first_contact = sp.simplify((centered.T * generator * centered)[0] / 3)
    assert sp.simplify(first_contact - 4 * sp.log(sp.Rational(3, 2))) == 0
    for n in (1, 2):
        actual = correlator(transfer, stationary, np.array(diagonal, dtype=float), n)
        assert abs(actual - (2 / 3) * (64 / 729) ** n) < 1e-14

    print(json.dumps({
        "scope": "Symbolic local kinetic contacts and existing v330, without a primitive marked source.",
        "laplace_test_contact": "2*a*s**2*u*exp(-s*u)",
        "yukawa_coordinate_contact": "4*a*s*exp(-s*u)*x1",
        "v330_step_matrix": str(expected),
        "v330_generator_12": "log(8/9)",
        "v330_generator_13": "-2*log(3)",
        "v330_actual_D_negative_time_derivative": "4*log(3/2)",
        "matrix_max_abs_error": float(np.max(abs(transfer - np.array(expected, dtype=float)))),
        "v330_actual_D_steps_1_2": "pass",
    }, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("repo", type=Path)
    main(parser.parse_args().repo.resolve())
