#!/usr/bin/env python3
"""Kleinster exakter Test der praemetrischen Kausalrekonstruktion.

Der Test isoliert drei Aussagen:

1. Eine interne, punktweise Reziprozitaet zwischen Feldstaerke F und
   Anregung H kann ohne zusaetzlichen Gradwechsler nur in n = 4 wirken.
2. Die isotrope lokale Antwort erfuellt die Closure Relation und ihre
   physische Fresnel Flaeche ist das Quadrat genau eines Lorentzkegels.
3. Eine uniaxiale Antwort ist die Negativkontrolle. Sie verletzt Closure
   und erzeugt zwei verschiedene quadratische Ausbreitungskegel.

Das ist kein allgemeiner Beweis des Tamm Rubilar Satzes. Es ist der kleinste
exakte Gegenfall, der die im Bericht benannte Voraussetzung unterscheidet.
"""

from __future__ import annotations

import sympy as sp


def constitutive_operator(epsilon: sp.Matrix) -> sp.Matrix:
    """6x6 Antwort auf der Basis (01, 02, 03, 23, 31, 12)."""
    zero = sp.zeros(3)
    identity = sp.eye(3)
    return sp.Matrix.vstack(
        sp.Matrix.hstack(zero, identity),
        sp.Matrix.hstack(-epsilon, zero),
    )


def closure_defect(kappa: sp.Matrix) -> tuple[sp.Expr, sp.Matrix]:
    """Liefert lambda^2 und Kappa^2 + lambda^2 I."""
    kappa_squared = sp.simplify(kappa * kappa)
    lambda_squared = sp.simplify(-sp.trace(kappa_squared) / 6)
    defect = sp.simplify(kappa_squared + lambda_squared * sp.eye(6))
    return lambda_squared, defect


def physical_fresnel(
    epsilon: sp.Matrix,
    omega: sp.Symbol,
    wavevector: sp.Matrix,
) -> sp.Expr:
    """Entfernt den longitudinalen Maxwell Faktor aus det M."""
    spatial_norm = (wavevector.T * wavevector)[0]
    principal_symbol = (
        wavevector * wavevector.T
        - spatial_norm * sp.eye(3)
        + omega**2 * epsilon
    )
    longitudinal_factor = epsilon[0, 0] * omega**2
    return sp.factor(principal_symbol.det() / longitudinal_factor)


def is_polynomial_square(expression: sp.Expr) -> bool:
    """Prueft, ob alle irreduziblen Faktoren gerade Vielfachheit haben."""
    _, factors = sp.factor_list(expression)
    return all(multiplicity % 2 == 0 for _, multiplicity in factors)


def main() -> None:
    n = sp.symbols("n", integer=True)
    dimensions = sp.solve(sp.Eq(2, n - 2), n)
    assert dimensions == [4]

    omega, kx, ky, kz = sp.symbols("omega k_x k_y k_z", real=True)
    a = sp.symbols("a", positive=True)
    wavevector = sp.Matrix([kx, ky, kz])

    isotropic_epsilon = sp.diag(a, a, a)
    isotropic_kappa = constitutive_operator(isotropic_epsilon)
    isotropic_scale, isotropic_defect = closure_defect(isotropic_kappa)
    isotropic_fresnel = physical_fresnel(
        isotropic_epsilon, omega, wavevector
    )
    isotropic_cone = -a * omega**2 + kx**2 + ky**2 + kz**2

    assert isotropic_scale == a
    assert isotropic_defect == sp.zeros(6)
    assert sp.expand(isotropic_fresnel - isotropic_cone**2) == 0
    assert is_polynomial_square(isotropic_fresnel)

    lorentz_form = isotropic_cone.subs(a, 1)
    variables = (omega, kx, ky, kz)
    quadratic_matrix = sp.hessian(lorentz_form, variables) / 2
    assert quadratic_matrix == sp.diag(-1, 1, 1, 1)

    control_epsilon = sp.diag(1, 1, 2)
    control_kappa = constitutive_operator(control_epsilon)
    control_scale, control_defect = closure_defect(control_kappa)
    control_fresnel = physical_fresnel(control_epsilon, omega, wavevector)
    ordinary_cone = -omega**2 + kx**2 + ky**2 + kz**2
    extraordinary_cone = -2 * omega**2 + kx**2 + ky**2 + 2 * kz**2

    assert control_scale == sp.Rational(4, 3)
    assert control_defect != sp.zeros(6)
    assert sp.expand(control_fresnel - ordinary_cone * extraordinary_cone) == 0
    assert not is_polynomial_square(control_fresnel)
    assert sp.expand(ordinary_cone - extraordinary_cone) != 0

    print("PASS: Reziprozitaet ohne Gradwechsler waehlt n = 4.")
    print(f"PASS: isotrope Closure, lambda^2 = {isotropic_scale}.")
    print(f"PASS: ein Lorentzkegel, Fresnel = ({isotropic_cone})^2.")
    print(
        "PASS Gegenfall: uniaxiale Antwort verletzt Closure mit "
        f"lambda^2 = {control_scale}."
    )
    print(f"PASS Gegenfall: Fresnel = ({ordinary_cone}) * ({extraordinary_cone}).")


if __name__ == "__main__":
    main()
