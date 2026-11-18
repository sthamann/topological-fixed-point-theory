#!/usr/bin/env python3
"""Exact checks of the Hodge bracket under an explicitly supplied 5-slot map.

This checks a conditional tensor identity, not a physical source selection.
"""
from itertools import combinations
from pathlib import Path
import sys

import sympy as sp

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "source_recursion_20261028"))
from check_native_spin5_source_record import hodge_amplitudes

PAIRS = list(combinations(range(5), 2))
C = hodge_amplitudes()


def wedge2(a):
    return sp.Matrix([[a.extract(i, j).det() for j in PAIRS] for i in PAIRS])


def transported(a):
    d = wedge2(a)
    return [d.T * c * d for c in C]


def check(a):
    adj = a.adjugate()
    k = transported(a)
    for m in range(5):
        assert k[m] == sum((adj[n, m] * C[n] for n in range(5)), sp.zeros(10))
    return k


def main():
    general = sp.Matrix([[1, 2, 0, 0, 0], [0, 2, 1, 0, 0],
                         [0, 0, 3, 1, 0], [0, 0, 0, 4, 1],
                         [1, 0, 0, 0, 5]])
    singular = sp.diag(1, 2, 3, 4, 0)
    cycle = sp.zeros(5)
    for j in range(5):
        cycle[(j + 1) % 5, j] = 1
    for a in (general, singular, cycle):
        check(a)
    assert general.det() == 122 and singular.det() == 0
    assert cycle**5 == sp.eye(5)
    cycle_effect = sum((k.T.conjugate() * k for k in transported(cycle)), sp.zeros(10))
    assert cycle_effect == 3 * sp.eye(10)
    anisotropic = sp.diag(1, 1, 1, 1, sp.Rational(1, 2))
    effect = sum((k.T.conjugate() * k for k in check(anisotropic)), sp.zeros(10)) / 3
    assert effect[PAIRS.index((0, 1)), PAIRS.index((0, 1))] == sp.Rational(1, 2)
    assert effect[PAIRS.index((0, 4)), PAIRS.index((0, 4))] == sp.Rational(1, 4)
    assert effect != sp.eye(10) * effect[0, 0]
    # Existing R087 interacting two-mode control, not a new TFPT vertex.
    t = sp.symbols("t", real=True, nonnegative=True)
    h = sp.diag(0, sp.log(4), sp.log(64), sp.log(128))
    v = sp.diag(0, 0, 0, -sp.log(2))
    creator = sp.zeros(4)
    creator[1, 0] = 1
    creator[3, 2] = 1
    u = (-h).exp()
    defect = u * creator - creator * u / 4
    integrand = (-(1 - t) * h).exp() * (v * creator - creator * v) * sp.exp(-t * sp.log(4)) * (-t * h).exp()
    integral = integrand.applyfunc(lambda x: sp.integrate(x, (t, 0, 1)))
    assert sp.simplify(defect + integral) == sp.zeros(4)
    assert defect[3, 2] == sp.Rational(1, 256)
    print("PASS: 15 exact cofactor identities, including singular and 5-cycle cases")
    print("PASS: 5-cycle order 5 and full original Kraus effect 3I")
    print("PASS: anisotropic effect 1/2 versus 1/4 (common scalar cannot fix)")
    print("PASS: full interaction identity; existing two-mode defect 1/256")


if __name__ == "__main__":
    main()
