"""Executable scope checks for the R113 scalar and the original regional action."""

from hashlib import sha256
import json
from pathlib import Path

import sympy as sp

from source import NATIVE

ORDERED = Path("experiments/theory-contracts/chat-triangulation-20261108/"
               "R113-ORIGINAL-UND-REFERENZ-20261109.txt")
REGIONAL = Path("docs/research/fundamental-plan/TFPT-REGIONALE-KOMPOSITION-20261112.txt")
TRANSPORT = Path("docs/research/fundamental-plan/TFPT-TRAEGERTRANSPORT-PICARD-KUBIK-20261111.txt")
INPUTS = (ORDERED, REGIONAL, TRANSPORT)


def charged_fibre(root):
    """Enclose the original fixed-history fermion scalar and trace separately."""
    from flint import arb, ctx
    from r113 import squared_singular_values

    raw = (root/NATIVE).read_bytes()
    fixture = json.loads(raw)
    if fixture["domain_id"] != "original-variable-H-quarter-full96-diagnostic-v1":
        raise ValueError("unexpected original conditional domain")
    for ref in fixture["source_originals"]:
        if sha256((root/ref["path"]).read_bytes()).hexdigest() != ref["sha256"]:
            raise ValueError(f"changed original input: {ref['path']}")

    with ctx.workprec(192):
        p, t, x = arb("1/4"), arb("4/5"), arb("1/2")
        vacuum_s, polar, trace_s = arb(0), arb(0), arb(0)
        pairs = 0
        for name, multiplicity in (("Yu", 6), ("Yd", 6), ("Ye", 2)):
            for squared in squared_singular_values(fixture["couplings"][name]):
                pairs += multiplicity
                energy = (p**2+x**2*squared).sqrt()
                gap = 2*((p/energy)*(t*energy/2).sinh()).asinh()
                sea_s = -t*(energy-p)
                vacuum_s += multiplicity*sea_s
                polar += multiplicity*(-t*energy+gap)
                trace_s += multiplicity*(sea_s + ((2+2*gap.cosh())/(2+2*(t*p).cosh())).log())
        if pairs != 42 or not (vacuum_s < trace_s and trace_s < 0 and polar < 0):
            raise ArithmeticError("the original charged block did not separate the conditional factors")
        return {
            "status": "CERTIFIED_CONDITIONAL",
            "source": str(ORDERED),
            "domain_sha256": sha256(raw).hexdigest(),
            "p": "1/4", "t": "4/5", "x": "1/2", "pairs": pairs,
            "log_unmarked_scalar_over_S": str(vacuum_s),
            "log_charged_trace_over_S": str(trace_s),
            "log_sewn_over_polar_scalar": str(polar),
            "log_unmarked_scalar_over_T": "0",
            "log_charged_trace_over_T": "0",
            "scope": "fixed even history, unprojected charged 84-coordinate factor only",
            "missing": ["absolute bosonic endpoints and history weights", "full physical selector",
                        "absolute rest contributions and bound reference normalization",
                        "ordered U and Y insertions on the same original source"],
            "physical_R113": dict.fromkeys(("Z0", "F0", "CUY")),
        }


def regional_transport():
    """Test the full linear span of the four native parabolic vector fields."""
    x = sp.Symbol("x")
    fields = (x**2, (x-1)**2/2, sp.S.One, (x+1)**2/2)

    def coefficients(expression):
        polynomial = sp.Poly(expression, x)
        return sp.Matrix([polynomial.nth(n) for n in range(3)])

    def bracket(left, right):
        return left*sp.diff(right, x) - right*sp.diff(left, x)

    basis = (sp.S.One, x, x**2)
    wedge_map = sp.Matrix.hstack(*(coefficients(bracket(basis[i], basis[j]))
                                   for i, j in ((0, 1), (0, 2), (1, 2))))
    rank = sp.Matrix.hstack(*(coefficients(field) for field in fields)).rank()
    determinant = wedge_map.det()
    if rank != 3 or determinant == 0:
        raise ArithmeticError("the original parabolic wedge test did not replay")
    native_bracket = sp.expand(bracket((1-x)**2, -x**2))
    return {
        "status": "EXCLUDED_SCOPED",
        "sources": [str(TRANSPORT), str(REGIONAL)],
        "original_generator_coefficients_in_1_x_x2": [list(map(str, coefficients(f))) for f in fields],
        "span_rank": rank, "bracket_wedge_determinant": str(determinant),
        "max_independent_commuting_rank": 1,
        "witness": {"bracket": str(native_bracket), "value_at_half": str(native_bracket.subs(x, sp.Rational(1, 2)))},
        "scope": "original Möbius action of the chiral circle and its linear Lie span",
    }


def neutral_return_budget(positive_rank, neutral_rank):
    """Necessary graded Ward test, after quotienting null vectors in each sector."""
    if any(type(rank) is not int or rank < 0 for rank in (positive_rank, neutral_rank)):
        raise ValueError("sector ranks must be nonnegative integers")
    shortfall = max(0, positive_rank-neutral_rank)
    return {"charged_rank": positive_rank, "neutral_rank": neutral_rank,
            "min_additional_neutral_rank": shortfall,
            "status": "EXCLUDED_SCOPED" if shortfall else "NECESSARY_RANK_OK"}


def replay(root):
    return {
        "status": "SCOPED_CHECKS",
        "original_source_sha256": {str(path): sha256((root/path).read_bytes()).hexdigest() for path in INPUTS},
        "charged_fixed_history": charged_fibre(root),
        "regional_transport": regional_transport(),
        "neutral_return": {
            "original_grade_one": neutral_return_budget(1, 1),
            "hypothetical_second_grade_one_at_same_vacuum_rank": neutral_return_budget(2, 1),
            "scope": "full Cartan charge gamma, positive grade-one form and metric adjoints "
                     "on a common invariant domain; F*F=I+E*E on that sector",
        },
        "physical_solution": False,
    }
