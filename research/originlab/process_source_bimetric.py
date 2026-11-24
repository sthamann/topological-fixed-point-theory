"""Exact H6h26 source, spectral, and quartic-amplitude checks.

The source Ward identity is not inferred from averaging a scalar Fourier
density.  Two distinct EP tests are kept separate: a species-dependent
MASSLESS pole (the actual Z48 toy ansatz), and a finite-range MASSIVE force.

Conventions: signature (-+++), g=eta+2*h_g/M_g, f=eta+2*h_f/M_f,
Einstein terms M_i^2 R_i/2, and HR potential action -k_edge*V_tau with
V_tau=(-4+e1+e3-4*e4)/2+(tau/2)*det(I-X).  The FP combination is one.

The exact quartic result is the difference of full tree four-point
amplitudes at fixed quadratic/cubic vertices, not the full baseline
amplitude or a gravitational forward-positivity verdict.
"""
from fractions import Fraction

import sympy as sp


Q_ORIGINAL_EXACT = sp.exp(-2*sp.pi)
PRIMARY_SOURCES = {
    "microscope_final": "https://arxiv.org/abs/2209.15487",
    "microscope_short_range": "https://arxiv.org/abs/2102.00022",
    "soft_graviton_universality": "https://doi.org/10.1103/PhysRev.135.B1049",
    "bimetric_eikonal": "https://arxiv.org/abs/1712.10020",
    "bimetric_screening": "https://arxiv.org/abs/2101.08795",
    "massless_pole_positivity": "https://arxiv.org/abs/2007.12667",
    "regge_gravitational_bounds": "https://arxiv.org/abs/2007.15009",
}


def _exact(value, name):
    if isinstance(value, (bool, float)):
        raise TypeError(f"{name} must be exact, not a float or boolean")
    if isinstance(value, Fraction):
        value = sp.Rational(value.numerator, value.denominator)
    elif isinstance(value, str):
        try:
            parsed = Fraction(value)
        except (ValueError, ZeroDivisionError) as error:
            raise ValueError(f"use a rational string or an exact SymPy expression for {name}") from error
        value = sp.Rational(parsed.numerator, parsed.denominator)
    else:
        value = sp.sympify(value)
    if not isinstance(value, sp.Expr) or value.has(sp.Float):
        raise TypeError(f"{name} must be an exact scalar expression")
    return sp.simplify(value)


def _positive(value, name):
    value = _exact(value, name)
    if value.is_positive is not True:
        raise ValueError(f"{name} must be provably positive")
    return value


def _q(value):
    value = _positive(value, "q")
    if (1-value).is_positive is not True:
        raise ValueError("the nondegenerate transfer comparison requires 0<q<1")
    return value


def bimetric_source_data(M_g, M_f, m_squared):
    """Physical TT residue matrices for two minimally coupled metric sources.

    For conserved sources the full tensor contractions use -T*T'/2 for
    the massless numerator and -T*T'/3 for the massive numerator.
    """
    Mg, Mf, m2 = (_positive(v, n) for v, n in (
        (M_g, "M_g"), (M_f, "M_f"), (m_squared, "m_squared")
    ))
    Meff2, alpha = Mg*Mg+Mf*Mf, Mf/Mg
    Meff = sp.sqrt(Meff2)
    u0, um = sp.Matrix([Mg, Mf])/Meff, sp.Matrix([-Mf, Mg])/Meff
    Dinv = sp.diag(1/Mg, 1/Mf)
    c0, cm = sp.simplify(Dinv*u0), sp.simplify(Dinv*um)
    x = sp.Symbol("x", positive=True)
    R0, Rm = sp.simplify(c0*c0.T), sp.simplify(cm*cm.T)
    return {
        "alpha": alpha, "M_eff_squared": Meff2,
        "M_rel_squared": sp.simplify(1/(1/Mg**2+1/Mf**2)),
        "m_squared": m2, "zero_mode": u0, "massive_mode": um,
        "mass_matrix": sp.simplify(m2*um*um.T),
        "massless_source_vector": c0, "massive_source_vector": cm,
        "massless_residue_matrix": R0, "massive_residue_matrix": Rm,
        "physical_TT_covariance": R0/x+Rm/(x+m2),
        "g_source_massive_weight_IR_normalized": alpha**2,
        "g_source_massive_fraction_of_total": sp.simplify(alpha**2/(1+alpha**2)),
        "tensor_scope": "Conserved-source numerators: P0=T.T'-(1/2)*tr(T)*tr(T'); Pm uses 1/3. The residue matrices alone describe the common TT scalar factor.",
    }


def massless_harmonic_eotvos(q, h_a, h_b):
    """Exact eta in the conditional Z48 MASSLESS-residue toy model.

    Assumptions: rho0/inertial_mass is common, g_A is proportional to
    1+q*h_A, and only this long-range force is being compared.  The Earth's
    common source factor and sqrt(1-q^2) cancel.  Physical positive-source
    use requires 1+q*h_a and 1+q*h_b positive.  No massive range enters.
    This identity does not derive the source map or a Ward identity.
    """
    q, ha, hb = _q(q), _exact(h_a, "h_a"), _exact(h_b, "h_b")
    denominator = 1+q*(ha+hb)/2
    if sp.simplify(denominator) == 0:
        raise ValueError("the average source coupling vanishes")
    return sp.simplify(q*(ha-hb)/denominator)


def massive_yukawa_eotvos(chi_e, chi_a, chi_b, z, form_factor=1):
    """Linear unscreened massive spin-2 eta for z=m*distance >=0.

    Distinct from massless_harmonic_eotvos.  For an extended Earth,
    form_factor encodes the source density/charge convolution.  Nonlinear
    screening is not included or claimed.  Equal test charges give zero.
    """
    ce, ca, cb = (_exact(v, n) for v, n in (
        (chi_e, "chi_e"), (chi_a, "chi_a"), (chi_b, "chi_b")
    ))
    z = _exact(z, "z")
    if z.is_nonnegative is not True:
        raise ValueError("z must be provably nonnegative")
    F = _positive(form_factor, "form_factor")
    attenuation = (1+z)*sp.exp(-z)*F
    return sp.simplify(
        sp.Rational(4, 3)*ce*(ca-cb)*attenuation /
        (1+sp.Rational(2, 3)*ce*(ca+cb)*attenuation)
    )


def transfer_band_vs_bimetric(q, mu_squared=1):
    """Compare the MP band with a single massive pole matching w and a.

    For every 0<q<1 the next moment differs.  Jensen gives a strict
    all-positive-x comparison; no numerical closeness test is used.
    """
    q, mu2 = _q(q), _positive(mu_squared, "mu_squared")
    r, Z = sp.simplify(q*q), sp.simplify(1-q*q)
    w = sp.simplify(r/Z)
    a = sp.simplify(r/(mu2*Z**2))
    bmp = sp.simplify(r*(1+r)/(mu2**2*Z**4))
    cmp = sp.simplify(r*(1+3*r+r*r)/(mu2**3*Z**6))
    pole = sp.simplify(mu2*Z)
    bbi = sp.simplify(w*(1+w)/pole**2)
    cbi = sp.simplify(w*(1+w)**2/pole**3)
    return {
        "q": q, "r": r, "Z": Z, "mu_squared": mu2,
        "massive_weight": w, "alpha_squared_matching_weight": w,
        "matching_pole_mass_squared": pole, "a": a,
        "b_MP": bmp, "b_bimetric": bbi,
        "c_MP": cmp, "c_bimetric": cbi,
        "b_difference": sp.simplify(bmp-bbi),
        "c_difference": sp.simplify(cmp-cbi),
        "KL_hankel_MP": sp.simplify(r**3/(mu2**2*Z**5)),
        "KL_hankel_one_massive_atom": sp.Integer(0),
        "band": (sp.simplify(mu2*(1-q)**2), sp.simplify(mu2*(1+q)**2)),
        "weight_match": True, "first_IR_match": True,
        "full_response_match": False,
        "jensen_strict_order": "C_MP,massive(x)<C_bimetric,massive(x) for every x>0, with w and a matched",
        "jensen_proof": "For u=1/s and normalized massive measure, phi(u)=u/(1+x*u) is strictly concave. The MP measure has nonzero inverse-mass variance.",
    }


def quartic_forward_data(tau, m_squared, M_rel_squared):
    """Exact change in the physical massive helicity-0 forward amplitude.

    It is the tau-dependent difference between full tree amplitudes:
    all exchange diagrams are unchanged because V_tau-V_0 starts at
    fourth order.  The common massless t pole is not assigned a finite
    positivity bound by this operation.
    """
    tau = _exact(tau, "tau")
    m2, Mr2 = _positive(m_squared, "m_squared"), _positive(M_rel_squared, "M_rel_squared")
    s = sp.Symbol("s", real=True)
    delta = sp.factor(tau*s*(s-4*m2)/(9*Mr2*m2))
    return {
        "tau": tau, "m_squared": m2, "M_rel_squared": Mr2,
        "s": s, "delta_amplitude_0000_forward": delta,
        "half_second_s_derivative": sp.simplify(sp.diff(delta, s, 2)/2),
        "L4_coefficient_of_det_H_cov": tau*m2/(2*Mr2),
        "scope": "fixed-field-convention difference of complete tree amplitudes; not the baseline amplitude, not a loop result, and not a full gravitational positivity verdict",
    }


def eikonal_family_data(tau):
    tau = _exact(tau, "tau")
    return {
        "beta": ((-4+tau)/2, (1-tau)/2, tau/2, (1-tau)/2, (-4+tau)/2),
        "c3": sp.Rational(1, 4), "d5": (tau-4)/96,
        "source_interval": "-alpha <= chi <= 1/alpha; chi_g=-alpha, chi_f=1/alpha",
        "scope": "necessary leading eikonal time-delay conditions under the paper's EFT and scale-separation assumptions; tau remains free",
        "source": PRIMARY_SOURCES["bimetric_eikonal"],
        "source_equations": "(4.4)-(4.9), (4.53)-(4.55); our massive mode has the opposite sign to that paper's v",
    }


def _analytic_checks():
    checks = {}
    zero = lambda expr: sp.simplify(expr) == 0
    Mg, Mf, m2 = sp.symbols("Mg Mf m2", positive=True)
    data = bimetric_source_data(Mg, Mf, m2)
    Meff = sp.sqrt(data["M_eff_squared"])
    alpha = Mf/Mg
    u0, um = data["zero_mode"], data["massive_mode"]
    checks["orthonormal_eigenmodes"] = zero((u0.T*um)[0]) and zero((u0.T*u0)[0]-1) and zero((um.T*um)[0]-1)
    checks["universal_massless_source_vector"] = sp.simplify(data["massless_source_vector"]-sp.ones(2, 1)/Meff) == sp.zeros(2, 1)
    checks["full_massive_residue_matrix"] = sp.simplify(data["massive_residue_matrix"]-sp.Matrix([[alpha**2, -1], [-1, alpha**-2]])/Meff**2) == sp.zeros(2)

    r, mu2 = sp.symbols("r mu2", positive=True)
    Z, w = 1-r, r/(1-r)
    a, bmp = r/(mu2*Z**2), r*(1+r)/(mu2**2*Z**4)
    pole = mu2*Z
    bbi = w*(1+w)/pole**2
    checks["single_pole_matches_w_and_a"] = zero(w/pole-a)
    checks["b_discriminant_equals_a_squared"] = zero(bmp-bbi-a*a)
    checks["MP_strict_KL_variance"] = zero(w*(bmp-a*a)-a*a-r**3/(mu2**2*Z**5))
    checks["one_pole_saturates_KL_variance"] = zero(w*(bbi-a*a)-a*a)
    u, mean_u, x = sp.symbols("u mean_u x", positive=True)
    phi = lambda v: v/(1+x*v)
    remainder = x*(u-mean_u)**2/((1+x*mean_u)**2*(1+x*u))
    checks["Jensen_exact_remainder"] = zero(phi(mean_u)+(u-mean_u)/(1+x*mean_u)**2-phi(u)-remainder)

    # Exact finite counterexample to replacing the actual quarter symmetry
    # by an unproved continuous rotation symmetry.
    q, harmonic = sp.symbols("q harmonic", positive=True)
    A = sp.zeros(8, 9)
    for j in range(8):
        A[j, j+1] = 1
        if j >= 3:
            A[j, j-3] = -q
    L = A.T*A
    R4 = sp.diag(*(sp.I**j for j in range(9)))
    R8 = sp.diag(*(sp.exp(sp.I*sp.pi*j/4) for j in range(9)))
    v = sp.eye(9)[:, 0]+harmonic*sp.eye(9)[:, 4]
    checks["quarter_symmetry_preserves_operator"] = sp.simplify(R4*L*R4.conjugate().T-L) == sp.zeros(9)
    checks["quarter_symmetry_allows_nonzero_fourth_harmonic"] = R4*v == v
    checks["continuous_rotation_counterexample"] = zero((R8*L*R8.conjugate().T-L)[0, 4]-2*q)
    checks["common_nonzero_harmonic_keeps_EP"] = massless_harmonic_eotvos(sp.Rational(1, 2), 3, 3) == 0
    checks["common_minimal_metric_keeps_EP"] = massive_yukawa_eotvos(-alpha, -alpha, -alpha, x) == 0

    # Independent determinant computation for physical on-shell tensors.
    E, k, m = sp.symbols("E k m", positive=True)
    eta = sp.diag(-1, 1, 1, 1)
    ex, ey = sp.Matrix([0, 1, 0, 0]), sp.Matrix([0, 0, 1, 0])
    lp, lm = sp.Matrix([-k/m, 0, 0, E/m]), sp.Matrix([-k/m, 0, 0, -E/m])
    ep = (ex*ex.T+ey*ey.T-2*lp*lp.T)/sp.sqrt(6)
    em = (ex*ex.T+ey*ey.T-2*lm*lm.T)/sp.sqrt(6)
    shell = lambda expr: sp.factor(sp.expand(expr).subs(E**2, k**2+m**2))
    pp, pm = sp.Matrix([E, 0, 0, k]), sp.Matrix([E, 0, 0, -k])
    checks["physical_helicity0_transversality"] = pp.T*ep == sp.zeros(1, 4) and pm.T*em == sp.zeros(1, 4)
    checks["physical_helicity0_trace_and_norm"] = all(shell(expr) == 0 for expr in (sp.trace(eta*ep), sp.trace(eta*em), sp.trace(eta*ep*eta*ep)-1, sp.trace(eta*em*eta*em)-1))
    a1, a2 = sp.symbols("a1 a2")
    det = sp.factor((a1*ep+a2*em).det())
    checks["physical_polarization_determinant"] = zero(det-4*E**2*k**2*a1*a2*(a1+a2)**2/(9*m**4))
    z1, z2, z3, z4 = sp.symbols("z1 z2 z3 z4")
    coefficient = sp.Poly(sp.expand(det.subs({a1: z1+z3, a2: z2+z4})), z1, z2, z3, z4).coeff_monomial(z1*z2*z3*z4)
    checks["four_leg_combinatorics"] = zero(coefficient-32*E**2*k**2/(9*m**4))
    tau, Mr2, s = sp.symbols("tau Mr2 s", positive=True)
    computed = sp.factor((tau*m**2*coefficient/(2*Mr2)).subs(k**2, E**2-m**2).subs(E**2, s/4))
    checks["physical_quartic_amplitude"] = zero(computed-tau*s*(s-4*m**2)/(9*Mr2*m**2))
    checks["physical_quartic_second_derivative"] = zero(sp.diff(computed, s, 2)/2-tau/(9*Mr2*m**2))
    c3, d5 = sp.Rational(1, 4), (tau-4)/96
    literature_beta = (48*d5+24*c3-6, -48*d5-18*c3+3, 48*d5+12*c3-1, -48*d5-6*c3, 48*d5)
    checks["eikonal_parameter_map"] = all(zero(a-b) for a, b in zip(literature_beta, eikonal_family_data(tau)["beta"]))
    return checks


def _json_exact(value):
    if isinstance(value, dict):
        return {key: _json_exact(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json_exact(item) for item in value]
    if isinstance(value, sp.MatrixBase):
        return [[str(item) for item in row] for row in value.tolist()]
    if isinstance(value, sp.Basic):
        return str(value)
    return value


def source_bimetric_audit():
    """Directly JSON-serializable audit; no floating-point selection gate."""
    q = Q_ORIGINAL_EXACT
    checks = _analytic_checks()
    tau, m2, Mr2 = sp.symbols("tau m2 Mr2", positive=True)
    return _json_exact({
        "status": "EXACT_SOURCE_AND_AMPLITUDE_AUDIT_WITH_CONDITIONAL_PHYSICAL_IDENTIFICATIONS",
        "checks": checks, "all_checks_pass": all(checks.values()),
        "physical_solution": False,
        "native_source": "H6h26 z48_zuordnung_quelle.py, Parts 3 and 4; main base 0172a74bd",
        "q": q,
        "native_massless_EP_claim": {
            "exact_conditional_formula": "eta=q*(h_A-h_B)/(1+q*(h_A+h_B)/2)",
            "leading_order_formula_valid_under_source_ansatz": True,
            "massive_range_relevant_to_this_massless_ansatz": False,
            "common_nonzero_harmonics_also_universal": True,
            "Ward_identity_derived_by_Fourier_average": False,
            "exact_source_e0_forced_by_experiment": False,
            "inertial_energy_normalization_derived": False,
            "illustrative_eta_limit": "1e-14 (coarse order-of-magnitude choice, not a new confidence interval)",
            "illustrative_delta_h_limit_at_zero_mean": str(sp.N(sp.exp(2*sp.pi)/10**14, 16)),
        },
        "rotation_scope": {
            "exact_transfer_covariance": "R_theta A_q R_theta^dagger=exp(-i theta) A_(q exp(4i theta))",
            "fixed_real_q_symmetry": "mu4; full U(1) requires an additional spurion or gauge construction",
            "counterexample": "e0+h*e4 is invariant under mu4 for every h",
            "full_U1_quadratic_densities": "J_(-n) J_n remain invariant for every n",
        },
        "massive_EP_scope": "A separate composition-dependent massive force requires its Earth charge, mass/range, source geometry and screening. Universal minimal g matter gives zero differential acceleration in the classical test-body limit.",
        "transfer_bimetric_comparison": transfer_band_vs_bimetric(q),
        "quartic": quartic_forward_data(tau, m2, Mr2),
        "eikonal": eikonal_family_data(tau),
        "gravitational_positivity_scope": "The common massless t-channel pole prevents an unqualified standard forward bound. The tau amplitude difference is finite, but its positivity does not follow from positivity of either full theory.",
        "next_decisive_tests": [
            "derive a common energy-normalized stress-tensor source and its gravitational Ward identity, including photons and binding energy",
            "decide exact MP-band dynamics versus a stated single-pole approximation by w,a,b or the KL inverse-mass variance",
            "match a physical four-source helicity-0 amplitude; the exact tau-dependent difference is known",
        ],
        "sources": PRIMARY_SOURCES,
    })
