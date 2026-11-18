"""Exact counterexamples for the H0 interpretation and spherical source routes."""
from itertools import product

import sympy as sp


def controls():
    checks = []

    def check(name, condition):
        if not condition:
            raise ArithmeticError(name)
        checks.append(name)

    # Central extensions of C2 with a fixed sign kernel. The normalized
    # two-cocycle table has only one free entry, at (1, 1).
    for sign in (1, -1):
        def multiplier(a, b):
            return sign ** (a * b)

        check(f"C2 cocycle law with sign {sign}", all(
            multiplier(a, b) * multiplier(a ^ b, c)
            == multiplier(b, c) * multiplier(a, b ^ c)
            for a, b, c in product(range(2), repeat=3)
        ))
    check("sign-valued section changes preserve the obstruction", all(
        lift**2 == 1 for lift in (1, -1)
    ))
    check("complex phase redefinition splits the C2 extension", sp.I**2 * (-1) == 1)

    # A commutator obstruction that survives complex rephasing. Its entire
    # twisted group algebra is already fixed by the Pauli representation.
    identity = sp.eye(2)
    x = sp.Matrix([[0, 1], [1, 0]])
    z = sp.diag(1, -1)
    klein = tuple(product(range(2), repeat=2))
    basis = {(a, b): x**a * z**b for a, b in klein}
    check("all sixteen projective Klein products", all(
        basis[a, b] * basis[c, d] == (-1)**(b*c) * basis[a ^ c, b ^ d]
        for (a, b), (c, d) in product(klein, repeat=2)
    ))
    check("commutator phase survives any scalar rephasing", x*z*x*z == -identity)
    check("no algebra directions are lost in reconstruction", sp.Matrix.hstack(*(
        operator.reshape(4, 1) for operator in basis.values()
    )).rank() == 4)

    # Even a fixed faithful trace and the exact same algebra admit two
    # positive CP transfers with different ratios between decay rates.
    pauli = (x, sp.I*x*z, z)
    spectra = []
    for rates in ((1, 1, 1), (1, 1, 2)):
        def transfer(operator):
            for axis, rate in zip(pauli, rates):
                decay = sp.Rational(1, 2)**rate
                operator = ((1 + decay)*operator + (1 - decay)*axis*operator*axis)/2
            return sp.simplify(operator)

        spectrum = [sp.Rational(1, 2)**(sum(rates) - rate) for rate in rates]
        check(f"unital tracial transfer {rates}", transfer(identity) == identity and all(
            sp.trace(transfer(axis)) == 0 for axis in pauli
        ))
        check(f"positive exact decay spectrum {rates}", all(
            0 < eigenvalue < 1 and transfer(axis) == eigenvalue*axis
            for axis, eigenvalue in zip(pauli, spectrum)
        ))
        spectra.append(spectrum)
    check("shared Z decay cannot calibrate the X decay", spectra[0][2] == spectra[1][2]
          and spectra[0][0] != spectra[1][0])

    # Add u/r² to the original local Gaussian action, leaving its kinetic
    # coefficient, reference condition and radial coordinate fixed.
    ell = sp.symbols("ell", integer=True, nonnegative=True)
    u = sp.symbols("u", positive=True)
    r, q = sp.symbols("r q", positive=True)
    half = sp.Rational(1, 2)
    a = ell + half
    frequency = sp.sqrt(a**2 + u)
    for exponent in (frequency - half, -frequency - half):
        check(f"radial Euler equation for exponent {exponent}", sp.simplify(
            exponent*(exponent + 1) - ell*(ell + 1) - u
        ) == 0)
    check("regular and decaying modes sew to precision 2*frequency/R", sp.simplify((
        (frequency - half) - (-frequency - half)
    )/r - 2*frequency/r) == 0)
    check("positive OU covariance coefficient", (1/(2*frequency)).is_positive)
    readout = 1/(8*sp.pi)
    check("undeformed Abel antipodal response", sp.limit(1/(4*sp.pi*(1 + q)), q, 1, dir="-") == readout)
    a_cont = sp.symbols("a", positive=True)
    correction = 1 - a_cont/sp.sqrt(a_cont**2 + u)
    check("alternating remainder has positive terms", sp.simplify(
        u / (sp.sqrt(a_cont**2 + u)*(sp.sqrt(a_cont**2 + u) + a_cont)) - correction
    ) == 0)
    check("alternating terms strictly decrease", sp.simplify(
        sp.diff(correction, a_cont) + u/(a_cont**2 + u)**sp.Rational(3, 2)
    ) == 0)
    check("absolute convergence by inverse square comparison", sp.limit(a_cont**2*correction, a_cont, sp.oo) == u/2)
    b0 = correction.subs({a_cont: half, u: sp.Rational(1, 4)})
    b1 = correction.subs({a_cont: 3*half, u: sp.Rational(1, 4)})
    interval = [sp.simplify((half - term)/(4*sp.pi)) for term in (b0, b0 - b1)]
    check("strict deformed bounds below the native target scalar", 0 < interval[0] < interval[1] < readout)
    ratio_squared = sp.simplify((frequency.subs(ell, 1)/frequency.subs(ell, 0))**2)
    check("ratio three forces u=0 within this action class", sp.simplify(
        ratio_squared - 9 + 32*u/(1 + 4*u)
    ) == 0)
    check("translation Ward test detects the local deformation", sp.diff(u/r**2, r) == -2*u/r**3)
    return {
        "status": "EXACT_CONTROL",
        "checks": checks,
        "scope": "comparison algebras and inverse-square Gaussian action; no native source identification",
        "pauli_transfer_spectra_without_constant": [[str(value) for value in s] for s in spectra],
        "gaussian_antipodal_response_at_zero": str(readout),
        "gaussian_antipodal_strict_bounds_at_u_1_4": [str(value) for value in interval],
        "gaussian_antipodal_bounds_decimal": [float(value) for value in interval],
        "radial_clock_ratio_squared": str(ratio_squared),
        "physical_solution": False,
    }
