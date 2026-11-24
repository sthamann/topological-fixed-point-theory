"""Scoped exact obstructions and remaining semidirect flavour freedom.

No full-Klein requirement is imposed on semidirect Z2+CP models.  The
thermal-weight theorem applies only to the specified basis and Q(t)
coefficient field, evaluated at a transcendental t.  It is not an exclusion
of general S4 models, theta functions, or nonlinear vacuum selection.
"""
from fractions import Fraction
import mpmath as mp
import sympy as sp


Z49_N_VALUES = (2, 3, 4, 5, 6, 7, 8, 10, 12, 16)
Z49_REACTOR_RANGE = ("0.020", "0.024")
NUFIT61_REACTOR_RANGE = ("0.02064", "0.02418")
NUFIT61_SOLAR_MAX = sp.Rational("0.3295")
NUFIT61_SOURCE = "https://www.nu-fit.org/sites/default/files/v61.tbl-parameters.pdf"


def tm1_mu_tau_matrix(theta):
    """Exact representative; theta is an explicit input, not a prediction."""
    theta = sp.sympify(theta)
    if theta.has(sp.Float):
        raise ValueError("theta must be exact")
    tbm = sp.Matrix([[sp.sqrt(sp.Rational(2, 3)), 1/sp.sqrt(3), 0],
                     [-1/sp.sqrt(6), 1/sp.sqrt(3), -1/sp.sqrt(2)],
                     [-1/sp.sqrt(6), 1/sp.sqrt(3), 1/sp.sqrt(2)]])
    c, s = sp.cos(theta), sp.sin(theta)
    return tbm*sp.Matrix([[1, 0, 0], [0, c, sp.I*s], [0, sp.I*s, c]])


def algebraic_integer_certificate(value):
    """Use the monic minimal polynomial, never a numerical integer fit."""
    value = sp.sympify(value)
    if value.has(sp.Float):
        raise ValueError("value must be exact")
    z = sp.Symbol("z")
    polynomial = sp.Poly(sp.minimal_polynomial(value, z), z, domain=sp.QQ)
    monic = polynomial.monic()
    return {
        "value": str(sp.simplify(value)),
        "primitive_minimal_polynomial": str(polynomial.as_expr()),
        "monic_minimal_polynomial": str(monic.as_expr()),
        "is_algebraic_integer": all(c.q == 1 for c in monic.all_coeffs()),
    }


def full_klein_certificate():
    """The exact 15-degree TM1 point cannot come from a finite full Klein."""
    U = tm1_mu_tau_matrix(sp.pi/12)
    omega = (-1+sp.I*sp.sqrt(3))/2
    T = sp.diag(1, omega, sp.conjugate(omega))
    S = [sp.simplify(2*U[:, i]*U[:, i].conjugate().T-sp.eye(3)) for i in range(3)]
    traces = [sp.simplify(sp.trace(T*Si)) for Si in S]
    certs = [algebraic_integer_certificate(z) for z in traces]
    norm = sp.simplify(traces[2]*traces[2].xreplace({sp.sqrt(3): -sp.sqrt(3)}))
    return {
        "status": "EXACT_FINITE_FULL_KLEIN_OBSTRUCTION",
        "theta_input": "pi/12",
        "traces": certs,
        "quadratic_field_norm_of_nonintegral_traces": str(norm),
        "checks": {
            "unitary_mixing": sp.simplify(U.conjugate().T*U) == sp.eye(3),
            "involutions": all(sp.simplify(Si*Si) == sp.eye(3) for Si in S),
            "commuting_klein": sp.simplify(S[0]*S[1]-S[2]) == sp.zeros(3),
            "third_trace_exact": sp.simplify(traces[2]+(2+sp.sqrt(3))/4) == 0,
            "second_trace_exact": sp.simplify(traces[1]+(2-sp.sqrt(3))/4) == 0,
            "extra_traces_not_integral": not certs[1]["is_algebraic_integer"] and not certs[2]["is_algebraic_integer"],
            "field_norm_is_one_sixteenth": norm == sp.Rational(1, 16),
        },
        "assumptions": [
            "charged Z3 with three distinct characters in an ordinary unitary representation",
            "the full nondegenerate Majorana Klein group is contained in the same finite flavour group",
            "the exact TM1 angle pi/12 and mu-tau reflection",
        ],
        "proof": "Finite-order eigenvalues are roots of unity, so their trace is an algebraic integer. The displayed traces are not.",
        "semidirect_Z2_CP_excluded": False,
        "semidirect_reason": "Only S1 need belong to the specified flavour group; S2 and S3 may be accidental mass symmetries outside it.",
    }


def _in_Qt(value, t):
    value = sp.sympify(value)
    if value.has(sp.Float) or value.free_symbols-{t}:
        raise ValueError("coefficient must lie in Q(t), with exact rational constants")
    numerator, denominator = sp.cancel(value).as_numer_denom()
    try:
        sp.Poly(numerator, t, domain=sp.QQ)
        sp.Poly(denominator, t, domain=sp.QQ)
    except (sp.PolynomialError, sp.polys.polyerrors.CoercionFailed) as exc:
        raise ValueError("coefficient is outside the scoped field Q(t)") from exc
    return sp.cancel(value)


def delta96_majorana_certificate():
    """The old modulus pattern can have a different, valid Klein origin."""
    omega = (-1+sp.I*sp.sqrt(3))/2
    Ue = sp.Matrix([[1, 1, 1], [1, omega, sp.conjugate(omega)],
                    [1, sp.conjugate(omega), omega]])/sp.sqrt(3)
    Un = sp.Matrix([[sp.I/sp.sqrt(2), 0, -sp.I/sp.sqrt(2)],
                    [1/sp.sqrt(2), 0, 1/sp.sqrt(2)], [0, 1, 0]])
    K1 = sp.diag(-1, -1, 1)
    K2 = sp.Matrix([[0, sp.I, 0], [-sp.I, 0, 0], [0, 0, -1]])
    a_gen = sp.Matrix([[0, 1, 0], [0, 0, 1], [1, 0, 0]])
    b_gen = -sp.Matrix([[0, 0, 1], [0, 1, 0], [1, 0, 0]])
    c_gen = sp.diag(sp.I, -sp.I, 1)
    M = sp.simplify(Un.conjugate()*sp.diag(1, 2, 3)*Un.conjugate().T)
    U = sp.simplify(Ue.conjugate().T*Un)
    P = U.applyfunc(lambda z: sp.simplify(z*sp.conjugate(z)))
    a, c = (2-sp.sqrt(3))/6, (2+sp.sqrt(3))/6
    expected = sp.Matrix([[sp.Rational(1, 3)]*3, [a, sp.Rational(1, 3), c],
                          [c, sp.Rational(1, 3), a]])
    clock = Un*sp.diag(-1, sp.I, -sp.I)*Un.conjugate().T
    checks = {
        "Klein_generators_are_Delta96_words": K1 == c_gen**2 and K2 == a_gen*b_gen*c_gen,
        "same_exact_modulus_pattern": sp.simplify(P-expected) == sp.zeros(3),
        "both_Klein_generators_preserve_Majorana_mass": all(sp.simplify(K.T*M*K-M) == sp.zeros(3) for K in (K1, K2)),
        "distinct_Takagi_masses": sp.simplify(Un.T*M*Un) == sp.diag(1, 2, 3),
        "full_Z4_is_not_a_mass_symmetry": sp.simplify(clock.T*M*clock-M) != sp.zeros(3),
    }
    return {"status": "SAME_PATTERN_DIFFERENT_NEUTRAL_RESIDUAL", "checks": checks,
            "generator_dictionary": "n=4, eta=i, K1=c^2, K2=a*b*c in the standard Delta generators",
            "mass_matrix": [[str(v) for v in M.row(i)] for i in range(3)],
            "scope": "The neutral symmetry is Klein, not the original full Z4 clock; this does not revoke the Z4 Majorana degeneracy."}


def rational_weight_angle_certificate(A, B, C, t):
    """For M=[[A,sqrt(2)B],[sqrt(2)B,C]], certify at transcendental t.

    This assumes the real plane basis u=(1,-1,-1,1)/2,
    w=(0,1,-1,0)/sqrt(2).  New algebraic coefficient fields are not admitted.
    """
    if not isinstance(t, sp.Symbol):
        raise ValueError("t must be an indeterminate; evaluation is a separate theorem")
    A, B, C = (_in_Qt(z, t) for z in (A, B, C))
    equation = sp.cancel((A-C)**2-24*B**2)
    numerator = sp.Poly(equation.as_numer_denom()[0], t, domain=sp.QQ)
    degenerate = sp.cancel(A-C) == 0 and B == 0
    if numerator.is_zero and not degenerate:
        raise RuntimeError("24 cannot be a square in Q(t)")
    return {
        "field": "Q(t)",
        "angle_equation": str(equation),
        "polynomial_numerator": str(numerator.as_expr()),
        "degenerate_angle_unselected": degenerate,
        "pi_over_12_excluded_at_transcendental_t": not degenerate,
        "scope": "Only this real basis and finite rational coefficient functions; not arbitrary S4 models.",
    }


def direct_reactor_fractions(n):
    """Exact angle fractions for the standard charged-Z3/full-Klein class.

    They include row/column reassignments of the analytic Delta(6n^2)
    formula.  Odd n has no Klein subgroup, by the order of its 2-Sylow group.
    """
    if isinstance(n, bool) or not isinstance(n, int) or n < 2:
        raise ValueError("n must be an integer >=2")
    if n % 2:
        return ()
    out = set()
    for k in range(1, n//2+1):
        a = Fraction(k, n)
        for f in (a, a+Fraction(1, 2), a+Fraction(1, 6), a+Fraction(2, 3),
                  Fraction(1, 6)-a, Fraction(2, 3)-a):
            f %= 1
            out.add(min(f, 1-f))
    return tuple(sorted(out))


def direct_reactor_scan(max_n=48, bounds=NUFIT61_REACTOR_RANGE):
    """Interval-certified finite diagnostic; this is not a full group scan."""
    if isinstance(max_n, bool) or not isinstance(max_n, int) or max_n < 2:
        raise ValueError("max_n must be an integer >=2")
    previous = mp.iv.dps
    mp.iv.dps = 60
    hits = []
    try:
        lower, upper = (mp.iv.mpf(b) for b in bounds)
        for n in range(2, max_n+1, 2):
            for f in direct_reactor_fractions(n):
                x = mp.iv.mpf(2)/3*mp.iv.sin(mp.iv.pi*f.numerator/f.denominator)**2
                inside = x.a >= lower.b and x.b <= upper.a
                outside = x.b < lower.a or x.a > upper.b
                if not inside and not outside:
                    raise RuntimeError("interval overlaps a window endpoint; no classification made")
                if inside:
                    expression = 2*sp.sin(sp.pi*sp.Rational(f.numerator, f.denominator))**2/3
                    hits.append({"n": n, "group_order": 6*n*n, "angle_over_pi": str(f),
                                 "s13_exact": str(expression), "s13_interval": str(x),
                                 "s12_TM2": str(sp.N(1/(3*(1-expression)), 18))})
    finally:
        mp.iv.dps = previous
    return {"max_n": max_n, "bounds": list(bounds), "hits": hits,
            "scope": "Analytic charged-Z3/full-Klein Delta family; interval test is reactor-only, not a likelihood."}


def semidirect_potential_data():
    """A stable counterfamily with only the named residual Z2+CP imposed."""
    x, y, a, b = sp.symbols("x y a b", real=True)
    v, lam, kappa = sp.symbols("v lambda kappa", positive=True)
    V = lam*(x*x+y*y-v*v)**2+kappa*(b*x-a*y)**2
    vacuum = {x: v*a/sp.sqrt(a*a+b*b), y: v*b/sp.sqrt(a*a+b*b)}
    H = sp.hessian(V, (x, y)).subs(vacuum).applyfunc(sp.simplify)
    radial, angular = sp.Matrix([a, b]), sp.Matrix([-b, a])
    checks = {
        "stationary_vacuum": all(sp.simplify(sp.diff(V, z).subs(vacuum)) == 0 for z in (x, y)),
        "positive_radial_eigenvalue": sp.simplify(H*radial-8*lam*v*v*radial) == sp.zeros(2, 1),
        "positive_angular_eigenvalue": sp.simplify(H*angular-2*kappa*(a*a+b*b)*angular) == sp.zeros(2, 1),
    }
    return {
        "potential": str(V), "conditions": "lambda,kappa,v>0 and (a,b)!=(0,0)",
        "vacuum": "(x,y)=+-v*(a,b)/sqrt(a^2+b^2)",
        "mass_block_after_CP_rephasing": "m*I+g*[[x,y],[y,-x]]",
        "tan_2theta": "y/x=b/a",
        "pi_over_12_requires": "b/a=1/sqrt(3)",
        "checks": checks,
        "first_free_choice": "direction ratio b/a in the allowed effective scalar potential",
        "scope": "Real CP-even mass spurions and only residual Z2+CP; no claim of a full S4-invariant UV scalar action.",
        "angle_dynamically_derived_from_TFPT": False,
    }


def minimal_doublet_alignment_data():
    """A possible angular invariant is not yet a physical mixing operator."""
    x, y = sp.symbols("x y", real=True)
    I3 = x*(x*x-3*y*y)
    rotation = {x: -x/2-sp.sqrt(3)*y/2, y: sp.sqrt(3)*x/2-y/2}
    checks = {
        "cubic_is_D3_invariant": sp.expand(I3.subs(rotation, simultaneous=True)-I3) == 0,
        "cubic_is_reflection_even": sp.expand(I3.subs(y, -y)-I3) == 0,
        "cubic_forbidden_only_after_extra_sign": sp.expand(I3.subs({x: -x, y: -y}, simultaneous=True)+I3) == 0,
        "sextic_has_pi_over_six_direction": sp.simplify(I3.subs({x: sp.sqrt(3)/2, y: sp.Rational(1, 2)})) == 0,
        "sextic_also_has_other_direction": I3.subs({x: 0, y: 1}) == 0,
    }
    return {"checks": checks, "cubic_invariant": str(I3),
            "extra_inputs": ["a genuine scalar doublet", "an extra sign symmetry forbidding I3", "a mass-operator coupling and relative vacuum choice"],
            "minimal_S4_limit": "The S4 doublet factors through S4/V4=S3. Doublet and singlet vevs alone preserve V4, so they cannot evade the full-Klein obstruction.",
            "PMNS_angle_derived": False}


def flavour_obstruction_audit():
    t = sp.Symbol("t", real=True)
    klein = full_klein_certificate()
    delta96 = delta96_majorana_certificate()
    scalar = semidirect_potential_data()
    doublet = minimal_doublet_alignment_data()
    old = direct_reactor_scan(24, Z49_REACTOR_RANGE)
    current = direct_reactor_scan(48, NUFIT61_REACTOR_RANGE)
    return {
        "status": "EXACT_SCOPED_OBSTRUCTIONS_AND_OPEN_SEMIDIRECT_DYNAMICS",
        "physical_solution": False,
        "full_klein": klein,
        "delta96_Majorana_pattern": delta96,
        "thermal_rational_weights": {
            "basis": ["(1,-1,-1,1)/2", "(0,1,-1,0)/sqrt(2)"],
            "matrix": "[[A(t),sqrt(2)*B(t)],[sqrt(2)*B(t),C(t)]], A,B,C in Q(t)",
            "q": "exp(-2*pi)",
            "transcendence": "Gelfond-Schneider applied to (-1)^(2*i)=exp(-2*pi), with log(-1)=i*pi.",
            "embedding": "Evaluation Q(t)->R at a transcendental q is injective; denominators cannot vanish there.",
            "proof": "The necessary equation (A-C)^2=24*B^2 cannot hold nondegenerately in Q(t), where 24 is not a square.",
            "three_corner_example": rational_weight_angle_certificate(sp.Rational(7, 3), -sp.Rational(1, 3), sp.Rational(2, 3), t),
            "degenerate_control": rational_weight_angle_certificate(1/(1-t), 0, 1/(1-t), t),
            "excludes_all_S4_models": False,
            "excluded_extensions": [],
            "outside_theorem": ["new algebraic coefficient fields", "infinite theta-function weights", "nonlinear vacuum selection", "different physical basis or field dictionary"],
        },
        "z49_scope": {"n_values": list(Z49_N_VALUES), "all_integers_through_16_scanned": False,
                      "original_reactor_window": list(Z49_REACTOR_RANGE)},
        "larger_direct_family": {"old_window_first_hit": old["hits"][0],
                                  "NuFIT61_first_hit": current["hits"][0],
                                  "NuFIT61_reactor_scan": current,
                                  "all_n_TM2_solar_exclusion": "s12=1/[3(1-s13)]>1/3>0.3295",
                                  "intervals_are_not_joint_likelihood": True},
        "semidirect_potential": scalar,
        "minimal_doublet_alignment": doublet,
        "all_exact_checks_pass": all(klein["checks"].values()) and all(delta96["checks"].values()) and all(scalar["checks"].values()) and all(doublet["checks"].values()),
        "sources": ["https://arxiv.org/abs/1305.3200", "https://arxiv.org/abs/1409.8005",
                    "https://arxiv.org/abs/1405.3678", NUFIT61_SOURCE,
                    "https://www.mathnet.ru/php/archive.phtml?jrnid=im&option_lang=eng&paperid=4924&wshow=paper"],
    }
