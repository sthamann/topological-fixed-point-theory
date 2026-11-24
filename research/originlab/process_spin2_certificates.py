"""Analytic certificates for the scalar TT form factors used by process_rules.

These helpers replace finite momentum scans for the declared candidate classes.
They certify a scalar quadratic response, not its identification with gravity,
nonlinear Ward identities, or a complete interacting quantum theory.

Coefficient convention: after dividing H by its positive infrared slope,
    H_normalized(x) = x - a*x**2 + b*x**3 - c*x**4 + O(x**5).
``mirror_weight`` means that the *total* positive massive spectral weight is
one when the massless residue is one. It does not assert a mirror-source map.

Inputs are exact integers, Fraction objects, or rational/decimal strings.
Binary floating point inputs are rejected instead of silently rounded.
Outputs retain Fraction values and a set of clauses for the existing lane;
callers may stringify Fraction values and sort sets for JSON serialization.
"""
from fractions import Fraction as Fr


CLAUSE_NAMES = (
    "no_euclidean_zero", "kallen_lehmann", "dirac_low_energy",
    "mirror_weight", "rational",
)
SCOPE = (
    "analytic scalar quadratic TT-form-factor certificate in the declared "
    "class; no physical source identification, nonlinear gravity, or UV claim"
)
DIRAC_A = Fr(3, 10)

# Both examples obey every one of the five scalar lane clauses, but have
# different higher coefficients from the selected single-auxiliary kernel.
EXACT_COUNTEREXAMPLES = (
    {
        "name": "Schur two channels w=1/8,3/8 m=5/6,5/2",
        "weights": (Fr(1, 8), Fr(3, 8)),
        "masses": (Fr(5, 6), Fr(5, 2)),
    },
    {
        "name": "Schur two channels w=1/4,1/4 m=1,5",
        "weights": (Fr(1, 4), Fr(1, 4)),
        "masses": (Fr(1), Fr(5)),
    },
)


def _fraction(value):
    if isinstance(value, bool) or not isinstance(value, (int, str, Fr)):
        raise TypeError("Use an integer, Fraction, or exact rational string; not a float.")
    return Fr(value)


def _certificate(kind, *, a, b, c, no_zero, kl, rational,
                 total_massive_weight=None, slope=Fr(1), proof=(), witnesses=None):
    checks = {
        "no_euclidean_zero": bool(no_zero),
        "kallen_lehmann": bool(kl),
        "dirac_low_energy": a == DIRAC_A,
        "mirror_weight": bool(kl and total_massive_weight == 1),
        "rational": bool(rational),
    }
    return {
        "kind": kind,
        "certificate": "analytic",
        "clauses": {name for name, holds in checks.items() if holds},
        "checks": checks,
        "a": a, "b": b, "c": c,
        "infrared_slope": slope,
        "massless_residue": 1 / slope,
        "total_massive_weight": total_massive_weight,
        "coefficient_convention": "H/infrared_slope = x-a*x^2+b*x^3-c*x^4+O(x^5)",
        "proof": list(proof),
        "witnesses": {} if witnesses is None else witnesses,
        "scope": SCOPE,
    }


def rational_certificate(weights, masses):
    """Certify H=x-sum(w_i*x^2/(x+m_i)) for w_i,m_i>0, sum(w_i)<1.

    The parent is A+x*B, with A=diag(0,m_i) and
    B=[[1,sqrt(w)^T],[sqrt(w),I]]. Its Schur criterion is 1-sum(w)>0.
    Hence A>=0 and B>0. By conjugation with B**(-1/2), the parent inverse
    is the resolvent of a positive operator, and its (0,0) entry is Stieltjes.
    Its infrared residue is one; lim(x*C(x),x->infinity)=1/(1-sum(w)).

    Invalid inputs are outside this sufficient certificate class: they raise
    ValueError rather than receiving an unsupported negative classification.
    Empty lists give Einstein's H=x as the zero-auxiliary limiting case.
    """
    ws = tuple(_fraction(value) for value in weights)
    ms = tuple(_fraction(value) for value in masses)
    if len(ws) != len(ms):
        raise ValueError("weights and masses must have the same length")
    if any(w <= 0 for w in ws) or any(m <= 0 for m in ms):
        raise ValueError("this certificate requires positive weights and positive masses")
    W = sum(ws, Fr(0))
    if W >= 1:
        raise ValueError("this certificate requires sum(weights) < 1")
    a = sum((w / m for w, m in zip(ws, ms)), Fr(0))
    b = sum((w / m**2 for w, m in zip(ws, ms)), Fr(0))
    c = sum((w / m**3 for w, m in zip(ws, ms)), Fr(0))
    massive_weight = W / (1 - W)
    result = _certificate(
        "positive_schur", a=a, b=b, c=c,
        no_zero=True, kl=True, rational=True,
        total_massive_weight=massive_weight,
        proof=(
            "A=diag(0,m_i) is nonnegative and B=[[1,sqrt(w)^T],[sqrt(w),I]] is positive by 1-sum(w)>0.",
            "The exact Schur complement is H; a positive affine-pencil resolvent supplies a positive Stieltjes measure.",
            "H/x -> 1 at zero and H/x -> 1-sum(w) at infinity; the massive measure has weight sum(w)/(1-sum(w)).",
        ),
        witnesses={
            "no_euclidean_zero": {"kind": "positive_affine_pencil", "kinetic_schur_margin": 1-W},
            "kallen_lehmann": {"kind": "positive_operator_resolvent", "massless_residue": Fr(1)},
        },
    )
    result.update({"weights": ws, "masses": ms, "weight_sum": W})
    return result


STANDARD_KINDS = (
    "einstein", "stelle", "dirac_local", "ledger_entire",
    "exp", "step", "laplace3", "exp_plus_einstein",
)


def candidate_certificate(kind, *, kappa=None):
    """Exact standard-candidate certificates; see STANDARD_KINDS.

    exp_plus_einstein defaults to kappa=2 and certifies every exact kappa>=2.
    This restricted range has a simple all-x positivity proof. Unsupported
    kappa values are not guessed from a numerical mesh.

    All spectral-action profiles have their full cosmological volume term
    tuned for a stationary flat background. laplace3 is normalized by
    integral_0^infinity (1+u)^(-3) du = 1/2, as in process_rules.
    """
    if kind not in STANDARD_KINDS:
        raise ValueError(f"unknown candidate kind: {kind!r}")
    if kappa is not None and kind != "exp_plus_einstein":
        raise ValueError("kappa is only supported for exp_plus_einstein")
    if kind == "einstein":
        result = rational_certificate((), ())
        result["kind"] = kind
        result["proof"] = ["H=x and C=1/x: one positive massless atom and no massive measure."]
        return result
    if kind == "stelle":
        return _certificate(
            kind, a=Fr(-1), b=Fr(0), c=Fr(0), no_zero=True, kl=False, rational=True,
            proof=("H=x(1+x)>0 for every x>0, but C=1/x-1/(x+1) has a negative massive residue.",),
            witnesses={"kallen_lehmann": {"kind": "negative_spectral_residue", "mass_squared": Fr(1), "residue": Fr(-1)}},
        )
    if kind == "dirac_local":
        return _certificate(
            kind, a=DIRAC_A, b=Fr(0), c=Fr(0), no_zero=False, kl=False, rational=True,
            proof=("H=x(1-3*x/10) has an exact positive Euclidean zero at x=10/3.",),
            witnesses={"no_euclidean_zero": {"kind": "exact_positive_zero", "x": Fr(10, 3)}},
        )
    if kind == "ledger_entire":
        return _certificate(
            kind, a=Fr(-1), b=Fr(1, 2), c=Fr(-1, 6), no_zero=True, kl=False, rational=False,
            proof=(
                "H=x*exp(x)>0 for x>0 and is nonpolynomial entire, hence not rational.",
                "C=exp(-x)/x has unit massless residue but x*C<1 for x>0; a positive massive measure would imply x*C>=1.",
            ),
            witnesses={"kallen_lehmann": {"kind": "massless_lower_bound_violation", "identity": "x*C=exp(-x)<1"}},
        )
    if kind == "step":
        return _certificate(
            kind, a=DIRAC_A, b=Fr(0), c=Fr(0), no_zero=False, kl=False, rational=False,
            proof=(
                "For f=1_[0,1] and 0<=x<=4, direct integration gives H=x-3*x^2/10.",
                "Thus x=10/3 is an exact positive zero. The full profile has a square-root threshold at x=4 and is not rational.",
            ),
            witnesses={"no_euclidean_zero": {"kind": "exact_positive_zero", "x": Fr(10, 3), "identity_interval": (Fr(0), Fr(4))}},
        )
    if kind == "laplace3":
        a, b, c = Fr(3, 5), Fr(3, 14), Fr(1, 15)
        return _certificate(
            kind, a=a, b=b, c=c, no_zero=False, kl=False, rational=False,
            proof=(
                "For f=(1+u)^(-3), F2(v)=1/[2(1+v)], integral f=1/2, and integral u*f=1/2.",
                "The normalized kernel has H/x->1 at zero and H->-6 at infinity, so continuity gives a positive Euclidean zero.",
                "Its J integral contains an inverse hyperbolic tangent and is not rational; b-a^2<0 also excludes a positive spectral measure.",
            ),
            witnesses={
                "no_euclidean_zero": {"kind": "opposite_endpoint_signs", "H_over_x_at_zero": Fr(1), "H_at_infinity": Fr(-6)},
                "kallen_lehmann": {"kind": "negative_inverse_spectral_moment", "b_minus_a_squared": b-a*a},
            },
        )
    if kind == "exp":
        a, b, c = DIRAC_A, Fr(1, 28), Fr(1, 360)
        return _certificate(
            kind, a=a, b=b, c=c, no_zero=False, kl=False, rational=False,
            proof=(
                "F(x)=integral_0^1 exp(-x*t*(1-t))dt and H=(3*x+12)*F-12.",
                "H/x->1 at zero and H->-6 at infinity, so continuity gives a positive Euclidean zero.",
                "H is nonpolynomial entire and b-a^2=-19/350<0 excludes a positive spectral measure.",
            ),
            witnesses={
                "no_euclidean_zero": {"kind": "opposite_endpoint_signs", "H_over_x_at_zero": Fr(1), "H_at_infinity": Fr(-6)},
                "kallen_lehmann": {"kind": "negative_inverse_spectral_moment", "b_minus_a_squared": b-a*a},
            },
        )
    # Only exp_plus_einstein remains. Preserve normalization explicitly.
    kap = Fr(2) if kappa is None else _fraction(kappa)
    if kap < 2:
        raise ValueError("this all-x positivity certificate covers kappa >= 2 only")
    slope = 1 + kap
    a, b, c = DIRAC_A/slope, Fr(1, 28)/slope, Fr(1, 360)/slope
    return _certificate(
        kind, a=a, b=b, c=c, no_zero=True, kl=False, rational=False, slope=slope,
        proof=(
            "1-exp(-y)<=y implies 1-F(x)<=x/6, hence H_exp+kappa*x>=(kappa-2)*x+3*x*F(x)>0 for kappa>=2 and x>0.",
            "After dividing H by 1+kappa, its coefficients violate the necessary Stieltjes Hankel condition a*c>=b^2.",
            "The exact defect is -13/[29400*(1+kappa)^2]; Euclidean positivity therefore does not imply positive spectral weights.",
        ),
        witnesses={
            "no_euclidean_zero": {"kind": "positive_lower_bound", "kappa": kap, "bound": "(kappa-2)*x+3*x*F(x)>0"},
            "kallen_lehmann": {"kind": "negative_inverse_moment_hankel", "a_c_minus_b_squared": a*c-b*b},
        },
    )
