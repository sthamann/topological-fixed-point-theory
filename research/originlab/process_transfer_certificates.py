"""Exact scalar response certificates for the quarter-step Hardy transfer.

The positive operator is L_q=A_q^* A_q, A_q=S^*(1-q*S^4), on the
infinite unilateral-shift space.  Its e_0 resolvent is divided by its
massless residue Z=1-q**2.  A positive mass-squared rescaling is explicit:

    C_normalized,scale(x) = C_raw(x/scale)/(scale*Z),
    H_scale(x) = scale*H_1(x/scale).

These are scalar operator certificates.  Identifying this positive response
with a connected graviton, choosing the physical sources, and fixing the
mass scale are additional physical tasks.  In particular, an operator's
one-particle massless eigenvector is not its many-body vacuum.

Inputs use exact integers, Fraction, rational strings, or supported SymPy
expressions.  Binary floats, free symbols, and unproved parameter classes
are rejected.  Supported q are nonnegative square roots of rational r<1,
the exact Dirac-matching algebraic q, and exp(-n*pi) for positive rational n.
Supported scales are positive rationals or positive rational multiples of
the unscaled coefficient a(q).  Outputs retain SymPy objects and clause
sets; JSON callers should stringify expressions and sort sets.
"""
from fractions import Fraction

import sympy as sp


DIRAC_A = sp.Rational(3, 10)
EXACT_R_DIRAC = (8-sp.sqrt(55))/3
EXACT_Q_DIRAC = sp.sqrt(EXACT_R_DIRAC)
EXACT_TRANSFER_Q_VALUES = (
    ("1/4", sp.Rational(1, 4)),
    ("1/2", sp.Rational(1, 2)),
    ("3/4", sp.Rational(3, 4)),
    ("q_dirac", EXACT_Q_DIRAC),
    ("exp(-pi)", sp.exp(-sp.pi)),
)
CLAUSE_NAMES = (
    "no_euclidean_zero", "kallen_lehmann", "dirac_low_energy",
    "mirror_weight", "rational", "transfer_memory",
)
SCOPE = (
    "exact positive scalar quadratic boundary resolvent of the declared "
    "infinite Hardy transfer; no derived 4D embedding, helicity, physical "
    "source dictionary, nonlinear gravitational dynamics, G, or UV completion"
)


def _exact(value, name):
    if isinstance(value, bool) or isinstance(value, float):
        raise TypeError(f"{name} must be exact; binary floats and booleans are unsupported")
    if isinstance(value, Fraction):
        result = sp.Rational(value.numerator, value.denominator)
    elif isinstance(value, int):
        result = sp.Integer(value)
    elif isinstance(value, str):
        aliases = {
            "q_dirac": EXACT_Q_DIRAC,
            "exp(-pi)": sp.exp(-sp.pi),
            "e^-pi": sp.exp(-sp.pi),
            "exp(-2*pi)": sp.exp(-2*sp.pi),
            "e^-2pi": sp.exp(-2*sp.pi),
        }
        if value in aliases:
            result = aliases[value]
        else:
            try:
                number = Fraction(value)
            except (ValueError, ZeroDivisionError) as error:
                raise ValueError(f"unsupported exact {name} string: {value!r}") from error
            result = sp.Rational(number.numerator, number.denominator)
    elif isinstance(value, sp.Expr):
        result = value
    else:
        raise TypeError(f"unsupported {name} type; use Fraction or an exact SymPy expression")
    if result.has(sp.Float) or result.free_symbols:
        raise ValueError(f"{name} must contain no floats or free symbols")
    return sp.simplify(result)


def _parameter(q):
    q = _exact(q, "q")
    r = sp.simplify(q*q)
    if r.is_Rational:
        if not 0 <= r < 1 or sp.simplify(q-sp.sqrt(r)) != 0:
            raise ValueError("the certificate requires 0 <= q < 1")
        return q, r, "sqrt_rational", "0<=r<1 is an exact rational inequality."
    if sp.simplify(q-EXACT_Q_DIRAC) == 0:
        return EXACT_Q_DIRAC, EXACT_R_DIRAC, "dirac_algebraic", (
            "7<sqrt(55)<8 implies 0<r=(8-sqrt(55))/3<1/3."
        )
    if q.func == sp.exp:
        n = sp.simplify(-q.args[0]/sp.pi)
        if n.is_Rational and n > 0:
            return q, sp.exp(-2*n*sp.pi), "nome", (
                "n>0 and pi>0 imply 0<exp(-n*pi)<1."
            )
    raise ValueError("q is outside the declared exact certificate classes")


def _scale(scale, base_a):
    scale = _exact(scale, "scale")
    if scale.is_Rational:
        if scale <= 0:
            raise ValueError("scale must be strictly positive")
        return scale, "rational", None
    if base_a != 0:
        ratio = sp.simplify(scale/base_a)
        if ratio.is_Rational and ratio > 0:
            return scale, "proportional_to_a", ratio
    raise ValueError("scale must be a positive rational or positive rational times a(q)")


def dirac_matching_scale(q):
    """Return the exact scale giving a=3/10 for every supported 0<q<1.

    This explicit continuum of matches is why Dirac matching determines q
    only after the relative transfer mass scale has independently been fixed.
    """
    _, r, _, _ = _parameter(q)
    if r == 0:
        raise ValueError("q=0 has a=0 at every positive scale")
    return sp.simplify(10*r/(3*(1-r)**2))


def transfer_certificate(q, *, scale=1):
    """Certify the normalized positive Hardy boundary response analytically.

    Coefficients refer to H=x-a*x**2+b*x**3-c*x**4+O(x**5).
    ``mirror_weight`` is only the equality of the total massive and
    massless spectral weights.  It does not prove a mirror-source rule.
    ``transfer_memory`` records membership in this declared operator class;
    it does not establish a common bundle/4D/gravity identification.

    At q=0 the response is Einstein's scalar 1/x, with zero continuum
    weight.  At q=1 the normalizable zero mode and its residue disappear,
    so q=1 and larger values are outside this normalization and rejected.
    """
    q, r, parameter_class, parameter_proof = _parameter(q)
    Z = sp.simplify(1-r)
    base_a = sp.simplify(r/Z**2)
    scale, scale_class, scale_ratio = _scale(scale, base_a)
    a = sp.simplify(base_a/scale)
    b = sp.simplify(r*(1+r)/(scale**2*Z**4))
    c = sp.simplify(r*(1+3*r+r*r)/(scale**3*Z**6))
    massive_weight = sp.simplify(r/Z)
    nontrivial = r != 0

    # For a nome and rational scale, equality with 3/10 would make r a
    # root of a rational quadratic.  But exp(-2*n*pi) is transcendental
    # for nonzero rational n by Gelfond--Schneider, applied to (-1)^(2*n*i).
    # This avoids using numerical symbolic-equality guesses for a negative.
    if parameter_class == "nome" and scale_class == "rational":
        dirac_match = False
        equality_proof = (
            "For rational scale, a=3/10 would imply a rational quadratic "
            "equation for r=exp(-2*n*pi), contrary to Gelfond-Schneider."
        )
    else:
        dirac_match = sp.simplify(a-DIRAC_A) == 0
        equality_proof = "Dirac matching is decided by exact algebra in the declared parameter class."
    if parameter_class == "nome":
        mirror_match = False
    else:
        mirror_match = sp.simplify(r-sp.Rational(1, 2)) == 0

    checks = {
        "no_euclidean_zero": True,
        "kallen_lehmann": True,
        "dirac_low_energy": bool(dirac_match),
        "mirror_weight": bool(mirror_match),
        "rational": not nontrivial,
        "transfer_memory": True,
    }
    band = None if not nontrivial else (
        sp.simplify(scale*(1-q)**2),
        sp.simplify(scale*(1+q)**2),
    )
    proof = [
        parameter_proof,
        "A_q=S^*(1-q*S^4), with 1-q*S^4 boundedly invertible, has Fredholm index +1; L_q=A_q^*A_q>=0.",
        "The e_0 sector has atom Z=1-q^2 at zero and density sqrt((b-s)*(s-a))/(2*pi*s) on a=(1-q)^2,b=(1+q)^2.",
        "The raw measure has total weight one; its continuum weight is q^2. Dividing by Z makes the massless residue one.",
        "Positive mass scaling pushes the spectral measure to s -> scale*s and preserves all spectral weights.",
        "For x>0, the resolvent is finite and strictly positive. Its reciprocal has no positive Euclidean zero.",
        equality_proof,
    ]
    if nontrivial:
        proof.append("The strictly positive continuum density on a nondegenerate interval gives a branch cut, so the response is not rational.")
    else:
        proof.append("At q=0 the entire e_0 measure is a unit atom at zero, so C=1/x exactly.")
    return {
        "kind": "hardy_quarter_transfer",
        "certificate": "analytic",
        "clauses": {name for name, holds in checks.items() if holds},
        "checks": checks,
        "q": q, "r": r, "Z": Z, "scale": scale,
        "parameter_class": parameter_class,
        "scale_class": scale_class,
        "a": a, "b": b, "c": c,
        "infrared_slope": sp.Integer(1),
        "ultraviolet_slope": Z,
        "massless_residue": sp.Integer(1),
        "raw_massless_residue": Z,
        "raw_total_massive_weight": r,
        "total_massive_weight": massive_weight,
        "total_spectral_weight": sp.simplify(1/Z),
        "spectral_band": band,
        "branch_cut": None if band is None else (-band[1], -band[0]),
        "normalized_spectral_density": (
            "sqrt((band_upper-s)*(s-band_lower))/(2*pi*scale*Z*s) on the band"
            if nontrivial else "zero"
        ),
        "invariant_b_over_a2": sp.simplify((1+r)/r) if nontrivial else None,
        "dirac_matching_scale": dirac_matching_scale(q) if nontrivial else None,
        "coefficient_convention": "H = x-a*x^2+b*x^3-c*x^4+O(x^5), unit massless residue",
        "nontrivial_memory": nontrivial,
        "proof": proof,
        "witnesses": {
            "no_euclidean_zero": {"kind": "positive_bounded_operator_resolvent"},
            "kallen_lehmann": {"kind": "explicit_positive_atom_and_density", "massless_residue": sp.Integer(1)},
            "moment_hankel": {
                "b_minus_a_squared": sp.simplify(b-a*a),
                "a_c_minus_b_squared": sp.simplify(a*c-b*b),
            },
        },
        "source_scope": {
            "one_particle_mass_squared": "A centered free 4D field retains its massless connected two-point pole; centering removes its one-point function only.",
            "vacuum_hamiltonian": "If L is instead a microscopic Hamiltonian with eta its unique vacuum, (O-omega(O))*eta is orthogonal to ker(L), so that connected spectral vector has no vacuum pole.",
            "basis_source": "For raw sources e_(4*j), the massless residue matrix is Z*q^(j+k); universal gravitational coupling needs a separate physical source normalization.",
        },
        "scope": SCOPE,
    }
