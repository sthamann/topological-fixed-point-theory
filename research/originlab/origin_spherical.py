"""Spherical scalar sewing: exact comparison controls, never a TFPT source.

The assumed action is kappa/2 integral_R3 |grad(phi)|^2.  A regular interior
and an exterior decaying to a fixed zero yield a full spherical response
and radial transfer.  The comparison exposes, rather than removes, the
stiffness, absolute reference, source dictionary and clock assumptions.

No existing ledger, paper, native source or index claim is promoted.
Import has no I/O.  ``controls()`` returns JSON-compatible results and
raises ArithmeticError on a failed check, including under Python -O.
"""
from __future__ import annotations

import json

import sympy as sp


def _zero(value):
    if isinstance(value, sp.MatrixBase):
        return all(sp.simplify(entry) == 0 for entry in value)
    return sp.simplify(value) == 0


def _require(checks, name, condition):
    if not bool(condition):
        raise ArithmeticError("spherical comparison failed: " + name)
    checks.append({"name": name, "passed": True})


def controls():
    """Evaluate the fixed analytic comparison and its counterexamples."""
    checks = []
    require = lambda name, condition: _require(checks, name, condition)
    ell = sp.symbols("ell", integer=True, nonnegative=True)
    kappa, radius, r = sp.symbols("kappa R r", positive=True)
    pi = sp.pi
    omega = ell + sp.Rational(1, 2)

    # Harmonic solutions and the two regions' outward normal derivatives.
    u_in = (r / radius) ** ell
    u_out = (radius / r) ** (ell + 1)

    def radial_laplacian(u):
        return sp.diff(u, r, 2) + 2 * sp.diff(u, r) / r - ell * (ell + 1) * u / r**2

    require("regular interior harmonic equation", _zero(radial_laplacian(u_in)))
    require("decaying exterior harmonic equation", _zero(radial_laplacian(u_out)))
    n_in = sp.simplify(sp.diff(u_in, r).subs(r, radius))
    n_out = sp.simplify(-sp.diff(u_out, r).subs(r, radius))
    require("interior outward DtN ell/R", _zero(n_in - ell / radius))
    require("exterior outward DtN (ell+1)/R", _zero(n_out - (ell + 1) / radius))
    sewn = sp.simplify(kappa * (n_in + n_out))
    require("physical sewn precision 2*kappa*Omega/R", _zero(sewn - 2 * kappa * omega / radius))
    require("angular measure precision has extra R squared", _zero(radius**2 * sewn - 2 * kappa * radius * omega))
    require("chi=sqrt(R)*f removes remaining R", _zero(radius**2 * sewn / radius - 2 * kappa * omega))

    # Actual radial change of variables, retaining the boundary derivative.
    tau = sp.symbols("tau", real=True)
    chi = sp.Function("chi")(tau)
    radial_coordinate = radius * sp.exp(tau)
    phi = chi / sp.sqrt(radial_coordinate)
    radial_energy = sp.simplify(radial_coordinate * sp.diff(phi, tau)**2)
    require("radial derivative chain rule", _zero(radial_energy - (sp.diff(chi, tau) - chi / 2)**2))
    cylinder_radial = sp.diff(chi, tau)**2 + chi**2 / 4
    boundary_derivative = -sp.diff(chi**2 / 2, tau)
    require("radial action includes total boundary derivative", _zero(radial_energy - cylinder_radial - boundary_derivative))
    require("angular action transforms with chi squared", _zero(radial_coordinate * phi**2 - chi**2))
    require("cylinder frequency squared", _zero(omega**2 - ell * (ell + 1) - sp.Rational(1, 4)))

    # Full unequal-time covariance from both the Euclidean distance and
    # the Legendre generating function; q=exp(-|tau-tau'|).
    q = sp.symbols("q", positive=True)
    cos_gamma = sp.symbols("cos_gamma", real=True)
    require("full cylinder covariance denominator", _zero((1 - 2 * q * cos_gamma + q**2) / q - (q + 1 / q - 2 * cos_gamma)))
    require("addition theorem cancels 2ell+1", _zero((2 * ell + 1) / (4 * pi) / (2 * kappa * omega) - 1 / (4 * pi * kappa)))
    require("every radial RP spectral coefficient is positive", (1 / (2 * kappa * omega)).is_positive is True)
    f1, f2 = sp.symbols("f1 f2")
    v1, v2 = sp.symbols("v1 v2", positive=True)
    vectors = sp.Matrix([f1, f2])
    rank_one = sp.Matrix([[v1**2, v1*v2], [v1*v2, v2**2]])
    quadratic = (sp.conjugate(vectors).T * rank_one * vectors)[0]
    factorized = sp.conjugate(v1*f1 + v2*f2) * (v1*f1 + v2*f2)
    require("radial RP summands factor into absolute squares", _zero(sp.expand(quadratic - factorized)))
    require("frequency generator is not conservative on constants", omega.subs(ell, 0) == sp.Rational(1, 2))

    # Stiffness is independent of the frequencies, and fixes amplitude.
    stiffness_controls = [
        {"kappa": str(k), "same_frequency": str(omega),
         "dimensionless_antipodal_response": str(sp.simplify(1 / (8*pi*k)))}
        for k in (sp.Rational(1, 2), sp.Integer(1), sp.Integer(2))
    ]
    require("three stiffness choices give distinct response amplitudes", len({item["dimensionless_antipodal_response"] for item in stiffness_controls}) == 3)
    require("frequency is independent of stiffness", sp.diff(omega, kappa) == 0)

    # The local equation alone does not select the exterior empty reference.
    constant, boundary_value = sp.symbols("constant boundary_value", real=True)
    bounded_exterior = constant + radius * (boundary_value - constant) / r
    require("bounded exterior has a free constant reference", _zero(sp.diff(bounded_exterior, r, 2) + 2 * sp.diff(bounded_exterior, r) / r))
    zero_normal = sp.simplify(-sp.diff(bounded_exterior, r).subs(r, radius))
    require("zero mode depends on absolute exterior reference", _zero(zero_normal - (boundary_value - constant) / radius))
    require("constant exterior makes constant DtN vanish", zero_normal.subs(constant, boundary_value) == 0)

    # Point readout and the distinction between a deck and an RP reflection.
    c = 1 / (8 * pi)
    require("dimensionless antipodal free Green response", _zero(radius / (4*pi*kappa*2*radius) - c/kappa))
    off = 1 / (4*pi*sp.sqrt(3))
    anti_gram = sp.Matrix([[c, off], [off, c]])
    anti_minus = sp.simplify(c - off)
    require("antipodal Gram has exact negative eigenvalue", _zero(anti_minus - (3 - 2*sp.sqrt(3))/(24*pi)) and anti_minus < 0)
    require("antipodal difference vector witnesses failure", _zero((sp.Matrix([[1, -1]]) * anti_gram * sp.Matrix([1, -1]))[0] - 2*anti_minus))
    plane_gram = sp.Matrix([[c, off], [off, 1/(4*pi)]])
    require("plane reflection control has positive determinant", sp.simplify(plane_gram.det()) == 1/(96*pi**2))

    abel = 1 / (4*pi*(1 + q))
    require("Abel angular regulator gives antipodal response", sp.limit(abel, q, 1, dir="-") == c)
    sharp = [sp.simplify(sum((-1)**l for l in range(L+1))/(4*pi)) for L in range(6)]
    require("sharp angular truncations do not converge at antipode", sharp == [1/(4*pi), 0, 1/(4*pi), 0, 1/(4*pi), 0])

    dimension_values = {}
    for dimension in range(3, 7):
        area = 2*pi**sp.Rational(dimension, 2)/sp.gamma(sp.Rational(dimension, 2))
        dimension_values[str(dimension)] = str(sp.simplify(1 / ((dimension-2)*area*2**(dimension-2))))
    require("ambient dimension signature is not dimension independent", dimension_values == {
        "3": "1/(8*pi)", "4": "1/(16*pi**2)",
        "5": "1/(64*pi**2)", "6": "1/(64*pi**3)",
    })

    # Local half-plane image formula, relative to full-plane heat diagonal.
    height, heat_time = sp.symbols("height heat_time", positive=True)
    image_term = sp.exp(-height**2/heat_time)
    require("one-sided local heat diagonal is not automatically one half",
            (1-image_term).subs(height, 0) == 0
            and (1+image_term).subs(height, 0) == 2
            and sp.limit(image_term, heat_time, 0, dir="+") == 0)

    # Separate abstract four-sector comparison: these are not native mu4
    # deck sectors and this is not an identification of a Jones index.
    a, b, overall, rate = sp.symbols("a b overall rate", positive=True)
    ps = sp.ones(4) / 4
    pt = sp.eye(4) - ps
    require("four-sector singlet and triplet are orthogonal projectors",
            _zero(ps**2-ps) and _zero(pt**2-pt) and _zero(ps*pt))

    def permutation(images):
        matrix = sp.zeros(4)
        for column, row in enumerate(images):
            matrix[row, column] = 1
        return matrix

    cycle3 = permutation([0, 2, 3, 1])
    double_transposition = permutation([1, 0, 3, 2])
    gram = a*ps + b*pt
    require("four-sector Gram is invariant under A4 generators",
            _zero(cycle3*gram*cycle3.T-gram)
            and _zero(double_transposition*gram*double_transposition.T-gram))
    state = gram/sp.trace(gram)
    require("all normalized sector weights stay one quarter",
            all(_zero(state[i, i]-sp.Rational(1, 4)) for i in range(4)))
    mixed = sp.factor(2*state[0, 1])
    require("mixed response is not fixed by quarter weights", _zero(mixed-(a-b)/(2*(a+3*b))))
    mixed_witnesses = [sp.simplify(mixed.subs({a:aa, b:bb}))
                       for aa, bb in [(1, 1), (2, 1), (sp.Rational(1, 2), 2)]]
    require("three symmetric positive states have distinct mixed responses",
            mixed_witnesses == [0, sp.Rational(1, 10), -sp.Rational(3, 26)])
    require("normalized state does not determine absolute source scale",
            _zero(overall*gram/sp.trace(overall*gram)-state))
    require("four-sector positive dynamics has an unconstrained rate",
            _zero(rate*pt*sp.ones(4, 1))
            and _zero(rate*pt-(rate*pt).T)
            and _zero(cycle3*(rate*pt)*cycle3.T-rate*pt))

    return {
        "status": "EXACT_COMPARISON_MODEL_ONLY",
        "F1_closed": False,
        "native_source_matched": False,
        "physical_3plus1_derived": False,
        "check_count": len(checks),
        "checks": checks,
        "assumptions": [
            "Free scalar Dirichlet action on Euclidean R3",
            "Positive stiffness kappa supplied",
            "Regular interior, exterior fixed to zero at infinity",
            "Spherical interface and radial conformal clock supplied",
            "Antipodal readout pair, not antipodal OS reflection",
            "Compatible smooth angular regulator",
        ],
        "DtN": {"interior": str(n_in), "exterior": str(n_out),
                "physical_precision_dS_R": str(sewn),
                "angular_precision_unrescaled_field": "2*kappa*R*Omega",
                "rescaled_field": "chi=sqrt(r)*phi, r=R*exp(tau)",
                "angular_precision_rescaled_field": "2*kappa*Omega"},
        "radial_action_boundary_term": "-kappa/4 * [integral_S2 chi^2 dOmega]_ends",
        "cylinder_frequency": str(omega),
        "cylinder_covariance": "1/(4*pi*kappa*sqrt(2*(cosh(s)-cos(gamma))))",
        "covariance_modes": "exp(-Omega*abs(s))/(2*kappa*Omega)",
        "radial_RP_proof": "Sum_lm |Sum_i f_i exp(-omega_l*tau_i) Y_lm(n_i)|^2/(2*kappa*omega_l) >= 0",
        "covariance_is_not_transfer": True,
        "frequency_is_not_conservative_Markov_generator": True,
        "dimensionless_antipodal_response": "1/(8*pi*kappa)",
        "stiffness_controls": stiffness_controls,
        "absolute_reference_obstruction": "bounded exterior A+R(f-A)/r gives DtN_0=(f-A)/R; decay chooses A=0",
        "antipodal_RP": {"passed": False, "gram": str(anti_gram),
                         "smallest_eigenvalue": str(anti_minus),
                         "smallest_eigenvalue_numeric": float(anti_minus)},
        "sharp_cutoff_antipodal_values": [str(value) for value in sharp],
        "ambient_dimension_antipodal_responses": dimension_values,
        "half_plane_heat_diagonal_relative_to_full_plane": {"Dirichlet_boundary": 0, "Neumann_boundary": 2, "interior_leading": 1},
        "separate_four_sector_control": {
            "native_index_identified": False,
            "gram": "Q=a*J4/4+b*(I4-J4/4), a,b>0",
            "normalized_diagonals": ["1/4"]*4,
            "Haar_product_probability_density": "1/(8*pi); not a physical response identity",
            "mixed_response": str(mixed),
            "mixed_response_witnesses": [str(value) for value in mixed_witnesses],
            "absolute_source_scale": "arbitrary positive overall factor",
            "RP_time_family": "a*Ps+b*exp(-rate*(t+s))*Pt, a,b,rate>0",
            "conservative_transfer": "Ps+exp(-rate*t)*Pt, arbitrary rate>0",
        },
        "native_matching_gates": [
            "Match actual field, measure and full boundary operator before c3",
            "Derive stiffness and absolute empty reference from original source",
            "Check native winding, Dirac/fermion sector and reflection action",
            "Identify physical clock and interacting 3+1D continuation",
        ],
        "primary_context_urls": [
            "https://arxiv.org/abs/2512.03587",
            "https://arxiv.org/abs/0705.0712",
            "https://arxiv.org/abs/2412.16652",
        ],
    }


if __name__ == "__main__":
    print(json.dumps(controls(), indent=2))
