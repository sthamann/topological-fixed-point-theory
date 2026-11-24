"""Exact certificates for K = det(I-W*Z)^(-nu) + eps det(Z) conjugate(det(W)).

The domain is the complex 2 by 2 matrix ball; nu and eps are real rationals.
In the Schur expansion eps changes only the (1,1) coefficient. For nu >= 1
all other coefficients are nonnegative. At nu = 0 only the constant survives.
For 0 < nu < 1, the unchanged (2,1) coefficient nu*(nu**2-1)/3 is negative;
for nu < 0, the (1,0) coefficient is negative. This proves the positivity
classification at every degree, including the allowed negative eps bound.
On a normalized rank-one disc K = (1-z*conjugate(w))**(-nu), so Hardy iff nu=1.
Writing D=det(I-W*Z), Cayley differentiation gives
  det(partial_Z) K = conjugate(det(W)) * [nu*(nu-1)*D**(-nu-1) + 2*eps].
The U(2,2) block identity for I-(gW)*(gZ) proves determinant covariance at
eps=0. The exact boost A=D=2I/sqrt(3), B=C=I/sqrt(3) sends 0 to I/2 and
has covariance defect eps/16 at Z=W=0, excluding every nonzero eps.
Noninteger weights use a consistent local/covering-group multiplier, not
a global principal-branch group representation. Weight nu=1 is unambiguous.
These are identities and negative witnesses, not tests at sampled points.
They certify this mathematical family, not a native physical identification.
"""

from fractions import Fraction
from numbers import Integral, Rational


def _rational(value):
    """Accept exact Python or SymPy rationals; never silently round a float."""
    if isinstance(value, Fraction):
        return value
    if isinstance(value, Integral):
        return Fraction(int(value))
    if isinstance(value, Rational):
        return Fraction(int(value.numerator), int(value.denominator))
    if getattr(value, "is_Rational", None) is True:
        return Fraction(int(value.p), int(value.q))
    if isinstance(value, str):
        return Fraction(value)
    raise TypeError("nu and eps must be exact real rationals, not floats")


def matrix_certificate(nu, eps):
    """Return JSON-serializable clauses, proof statements and negative witnesses."""
    nu, eps = _rational(nu), _rational(eps)
    epsilon_lower_bound = -nu * (nu - 1) / 2
    positive = (nu == 0 or nu >= 1) and eps >= epsilon_lower_bound
    hardy = nu == 1
    wave = (nu in (0, 1) and eps == 0) or (nu == -1 and eps == -1)
    covariant = eps == 0
    predicates = {
        "positive": positive,
        "hardy_slice": hardy,
        "wave_equation": wave,
        "conformally_covariant": covariant,
    }
    evidence = {
        "positive": (
            "Schur coefficients at all degrees: nu in {0} union [1,infinity), "
            "and eps >= -nu*(nu-1)/2; bound=" + str(epsilon_lower_bound)
        ),
        "hardy_slice": (
            "K(zE11,wE11)=(1-z*conjugate(w))^(-nu); first coefficient=" + str(nu)
        ),
        "wave_equation": (
            "det(partial_Z)K=conjugate(det(W))*"
            "[nu*(nu-1)*det(I-W*Z)^(-nu-1)+2*eps]"
        ),
        "conformally_covariant": (
            "U(2,2) determinant identity for eps=0; "
            "exact boost at Z=W=0 has defect eps/16=" + str(eps / 16)
        ),
        "multiplier_scope": (
            "Integer determinant weight; the multiplier is single-valued."
            if nu.denominator == 1 else
            "Noninteger determinant weight: use a consistent local multiplier "
            "or a covering group; no global principal-branch group law is claimed."
        ),
    }
    witnesses = {}
    if not positive:
        if nu < 0:
            partition, coefficient = [1, 0], nu
        elif 0 < nu < 1:
            partition, coefficient = [2, 1], nu * (nu * nu - 1) / 3
        else:
            partition, coefficient = [1, 1], nu * (nu - 1) / 2 + eps
        witnesses["positive"] = {
            "kind": "negative_schur_coefficient",
            "partition": partition,
            "coefficient": str(coefficient),
        }
    if not hardy:
        witnesses["hardy_slice"] = {
            "kind": "first_coefficient_mismatch",
            "actual": str(nu), "required": "1",
        }
    if not wave:
        a = nu * (nu - 1)
        value = (a + 2 * eps) / 4
        witnesses["wave_equation"] = {
            "kind": "nonzero_cayley_derivative",
            "W": [["1/2", "0"], ["0", "1/2"]],
        }
        if value:
            witnesses["wave_equation"].update({
                "Z": [["0", "0"], ["0", "0"]],
                "value": str(value),
            })
        else:
            # The constant cancellation at D=1 is not an identity: a!=0,
            # nu!=-1, and (9/16)^(-nu-1)!=1 for real nu!=-1.
            witnesses["wave_equation"].update({
                "Z": [["1/2", "0"], ["0", "1/2"]],
                "value": "(" + str(a) + "*(9/16)^(" + str(-nu-1)
                         + ")+" + str(2*eps) + ")/4",
                "nonzero_reason": "a!=0, exponent!=0, 0<9/16<1",
            })
    if not covariant:
        witnesses["conformally_covariant"] = {
            "kind": "exact_boost_defect",
            "Z": [["0", "0"], ["0", "0"]],
            "W": [["0", "0"], ["0", "0"]],
            "A_and_D": "2I/sqrt(3)", "B_and_C": "I/sqrt(3)",
            "g_zero": "I/2", "defect": str(eps / 16),
        }
    return {
        "nu": str(nu), "eps": str(eps),
        "clauses": sorted(key for key, holds in predicates.items() if holds),
        "evidence": evidence,
        "witnesses": witnesses,
        "integer_determinant_multiplier": nu.denominator == 1,
        "method": "analytic_identity_and_exact_schur_certificate",
    }
