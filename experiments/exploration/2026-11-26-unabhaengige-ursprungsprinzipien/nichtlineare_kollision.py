#!/usr/bin/env python3
"""Exakter J-Test in der Khan-Penrose-Wechselwirkungszone.

Die Metrikfunktionen sind Gl. (51)--(54) aus M. Halil,
J. Math. Phys. 20 (1979) 120--125, DOI 10.1063/1.523951.
Halil identifiziert diesen Grenzfall ausdrücklich mit Khan--Penrose.

Getestet wird nur die glatte Wechselwirkungszone u>0, v>0 vor der
Fokussierung. Die distributionswertigen Wellenfronten sind nicht Teil der
symbolischen Rechnung.
"""

from __future__ import annotations

import sympy as sp


def khan_penrose_functions():
    """Liefert U,V,M der exakten Vakuummetrik als Funktionen von p,q."""
    p, q = sp.symbols("p q", real=True)
    r = sp.sqrt(1 - p**2)
    w = sp.sqrt(1 - q**2)
    t = sp.sqrt(1 - p**2 - q**2)

    exp_minus_m = t**3 / (r * w * (p * q + r * w) ** 2)
    U = -sp.log(t**2)
    V = sp.log((r + q) / (r - q) * (w + p) / (w - p))
    M = -sp.log(exp_minus_m)
    return p, q, U, V, M, exp_minus_m


def weyl_scalars():
    """Berechnet Psi_0,Psi_2,Psi_4 in der Rosen-Szekeres-NP-Tetrade.

    p=a*u und q=b*v führen die zwei unabhängig schaltbaren Impulsstärken
    a,b ein. Daher ist d/du=a*d/dp und d/dv=b*d/dq.
    """
    a, b = sp.symbols("a b", real=True)
    p, q, U, V, M, exp_minus_m = khan_penrose_functions()
    exp_m = 1 / exp_minus_m

    psi0 = -sp.Rational(1, 2) * exp_m * b**2 * (
        sp.diff(V, q, 2)
        - sp.diff(U, q) * sp.diff(V, q)
        + sp.diff(M, q) * sp.diff(V, q)
    )
    psi2 = sp.Rational(1, 2) * exp_m * a * b * sp.diff(M, p, q)
    psi4 = -sp.Rational(1, 2) * exp_m * a**2 * (
        sp.diff(V, p, 2)
        - sp.diff(U, p) * sp.diff(V, p)
        + sp.diff(M, p) * sp.diff(V, p)
    )
    return p, q, a, b, psi0, psi2, psi4


def exact_interaction_test():
    """Ein algebraisch exakter Punkt tief in der regulären Interaktionszone."""
    p, q, a, b, psi0, psi2, psi4 = weyl_scalars()
    # a=2,b=3 und u=1/4,v=1/6 geben p=q=1/2; t^2=1/2>0.
    point = {p: sp.Rational(1, 2), q: sp.Rational(1, 2), a: 2, b: 3}
    values = tuple(
        sp.radsimp(sp.simplify(x.subs(point))) for x in (psi0, psi2, psi4)
    )
    p0, p2, p4 = values

    # Fuer psi1=psi3=0 lauten die tetradeninvarianten Weyl-Polynome:
    invariant_i = sp.simplify(p0 * p4 + 3 * p2**2)
    invariant_j = sp.radsimp(sp.simplify(p2 * (p0 * p4 - p2**2)))
    discriminant = sp.simplify(invariant_i**3 - 27 * invariant_j**2)

    assert values == (-54 * sp.sqrt(6), 39 * sp.sqrt(2), -24 * sp.sqrt(6))
    assert invariant_i == 16902
    assert invariant_j == 184626 * sp.sqrt(2)
    assert discriminant == 2987837829504
    return values, invariant_i, invariant_j, discriminant


def controls():
    """Ausschaltkontrollen und der bewusst lineare Zweifeld-Gegenfall."""
    p, q, a, b, psi0, psi2, psi4 = weyl_scalars()

    # Abseits der distributionswertigen Front bleibt ein einzelner KP-Impuls
    # in dieser Zone flach. Insbesondere verschwinden I und J.
    only_left = {p: sp.Rational(1, 2), q: 0, a: 2, b: 0}
    only_right = {p: 0, q: sp.Rational(1, 2), a: 0, b: 3}
    single_left = tuple(
        sp.simplify(x.subs(only_left)) for x in (psi0, psi2, psi4)
    )
    single_right = tuple(
        sp.simplify(x.subs(only_right)) for x in (psi0, psi2, psi4)
    )
    assert single_left == (0, 0, 0)
    assert single_right == (0, 0, 0)

    # Lineare Addition zweier Typ-N-Weylanteile: beide transversalen Anteile
    # sind da, aber ohne die nichtlinear erzeugte longitudinale Komponente.
    linear_psi0 = sp.Integer(9)
    linear_psi2 = sp.Integer(0)
    linear_psi4 = sp.Integer(4)
    linear_i = linear_psi0 * linear_psi4 + 3 * linear_psi2**2
    linear_j = linear_psi2 * (linear_psi0 * linear_psi4 - linear_psi2**2)
    assert linear_i == 36 and linear_j == 0
    return single_left, single_right, linear_i, linear_j


def ricci_check_at_sample():
    """Unabhängiger Vakuumcheck R_ab=0 am selben exakten Punkt."""
    u, v, x, y = sp.symbols("u v x y", real=True)
    a, b = sp.Integer(2), sp.Integer(3)
    p, q = a * u, b * v
    r, w = sp.sqrt(1 - p**2), sp.sqrt(1 - q**2)
    t = sp.sqrt(1 - p**2 - q**2)
    exp_minus_m = t**3 / (r * w * (p * q + r * w) ** 2)
    exp_v = (r + q) / (r - q) * (w + p) / (w - p)

    metric = sp.zeros(4)
    metric[0, 1] = metric[1, 0] = exp_minus_m
    metric[2, 2] = -t**2 * exp_v
    metric[3, 3] = -t**2 / exp_v
    inverse = sp.simplify(metric.inv())
    coordinates = (u, v, x, y)

    gamma = [
        [[sp.Integer(0) for _ in range(4)] for _ in range(4)]
        for _ in range(4)
    ]
    for k in range(4):
        for i in range(4):
            for j in range(4):
                gamma[k][i][j] = sp.Rational(1, 2) * sum(
                    inverse[k, ell]
                    * (
                        sp.diff(metric[ell, j], coordinates[i])
                        + sp.diff(metric[ell, i], coordinates[j])
                        - sp.diff(metric[i, j], coordinates[ell])
                    )
                    for ell in range(4)
                )

    sample = {u: sp.Rational(1, 4), v: sp.Rational(1, 6)}
    ricci_at_sample = sp.zeros(4)
    for i in range(4):
        for j in range(4):
            component = sum(
                sp.diff(gamma[k][i][j], coordinates[k])
                - sp.diff(gamma[k][i][k], coordinates[j])
                + sum(
                    gamma[k][k][ell] * gamma[ell][i][j]
                    - gamma[k][j][ell] * gamma[ell][i][k]
                    for ell in range(4)
                )
                for k in range(4)
            )
            ricci_at_sample[i, j] = sp.simplify(component.subs(sample))
    assert ricci_at_sample == sp.zeros(4)
    return ricci_at_sample


def main():
    values, invariant_i, invariant_j, discriminant = exact_interaction_test()
    single_left, single_right, linear_i, linear_j = controls()
    ricci = ricci_check_at_sample()

    print("Khan-Penrose, p=q=1/2, a=2, b=3")
    print(f"Psi0, Psi2, Psi4 = {values}")
    print(f"I = {invariant_i}")
    print(f"J = {invariant_j}")
    print(f"I^3 - 27 J^2 = {discriminant}")
    print(f"Ricci tensor at the same point = {ricci.tolist()}")
    print(f"single-pulse controls = {single_left}, {single_right}")
    print(f"linear two-type-N control: I={linear_i}, J={linear_j}")
    print("VERDICT: exact nonlinear collision does not preserve J=0")


if __name__ == "__main__":
    main()
