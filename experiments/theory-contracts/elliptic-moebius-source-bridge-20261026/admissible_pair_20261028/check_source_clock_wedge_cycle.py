"""Replay the signed E8 Gram and the literal carrier five-cycle.

This uses the existing actual Chevalley source, the existing A4 sector chart,
and the existing ordered E8 lattice cocycle.  It checks a source-clock
compatibility question; it does not introduce a physical P1 clock or reader.
"""

from collections import Counter
from fractions import Fraction
import hashlib
import importlib.util
import json
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
WORKSPACE = HERE.parent
REPO = Path("/Users/stefanhamann/Projekte/tfpt-theoryv4")
A4_CHECKER = WORKSPACE / "research-20261028-hodge-complement-channel" / "check_e8_a4_response.py"
if not A4_CHECKER.is_file():
    A4_CHECKER = REPO / "experiments/theory-contracts/elliptic-moebius-source-bridge-20261026/hodge_complement_20261028/check_e8_a4_response.py"
E8_RESPONSE = REPO / "experiments/theory-contracts/seam-collar-multifield-20261016/e8_boundary_response.py"
OUTPUT = HERE / "source-clock-wedge-cycle.json"


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    need(spec is not None and spec.loader is not None, f"load spec for {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def physical(root):
    return tuple(Fraction(x, 2) for x in root)


RHO_A7 = tuple(Fraction(x, 2) for x in (7, 5, 3, 1, -1, -3, -5, -7))


def h8_affine_grade(root):
    """Return (b,Q) for h_beta(a)=b+aQ."""
    beta = physical(root)
    charge = sum(beta)
    base = Fraction(8) - sum(x*y for x, y in zip(RHO_A7, beta)) - charge/2
    return base, charge


def carrier_cycle(root):
    """Literal e0->e1->...->e4->e0, fixing e5,e6,e7."""
    return (root[4], root[0], root[1], root[2], root[3],
            root[5], root[6], root[7])


def fraction_text(value):
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def counter_rows(counter):
    return [
        {
            "base": fraction_text(base),
            "charge": fraction_text(charge),
            "multiplicity": counter[(base, charge)],
        }
        for base, charge in sorted(counter, key=lambda item: (item[1], item[0]))
    ]


def main():
    sys.path.insert(0, str(REPO))
    a4 = load_module("existing_e8_a4_response", A4_CHECKER)
    response = load_module("existing_e8_boundary_response", E8_RESPONSE)

    source = a4._actual_current_source()
    roots = set(source["chevalley"].roots)
    need(len(roots) == 240, "actual Chevalley source has 240 E8 roots")

    insertion = sorted(root for root in roots if a4.carrier_type(root) == "5")
    need(len(insertion) == 50, "actual (5,10) insertion sector has 50 roots")
    insertion_set = set(insertion)
    need(all(tuple(-x for x in root) not in insertion_set for root in insertion),
         "insertion sector contains no opposite root")

    cycle_images = [carrier_cycle(root) for root in insertion]
    need(set(cycle_images) == insertion_set,
         "literal carrier five-cycle preserves the actual insertion sector")
    need(all(carrier_cycle(carrier_cycle(carrier_cycle(carrier_cycle(carrier_cycle(root))))) == root
             for root in roots), "literal carrier cycle has order five on all E8 roots")

    # It acts on the carrier fundamental label and fixes the centralizer A4
    # root coordinates literally, so this is the actual A4 x A4 chart rather
    # than an arbitrary permutation of fifty labels.
    need(all(a4.carrier_label(carrier_cycle(root), "5")
             == (a4.carrier_label(root, "5") + 1) % 5
             for root in insertion), "cycle advances every carrier 5 label")
    family_roots = [root for root in roots if len(set(root[:5])) == 1]
    need(len(family_roots) == 20, "actual family centralizer has 20 roots")
    need(all(carrier_cycle(root) == root for root in family_roots),
         "carrier cycle fixes the family A4 roots pointwise")

    beta1 = (2, 0, 0, 0, 0, 0, 0, -2)
    beta2 = (0, 2, 0, 0, 0, 0, 0, -2)
    need(beta1 in insertion_set and beta2 in insertion_set,
         "e1-e8 and e2-e8 are actual insertion roots")
    need(carrier_cycle(beta1) == beta2, "literal cycle maps e1-e8 to e2-e8")
    need(a4.carrier_label(beta1, "5") == 0 and a4.carrier_label(beta2, "5") == 1,
         "witnesses have consecutive carrier labels")

    grade_counts = Counter(h8_affine_grade(root) for root in insertion)
    need(sum(grade_counts.values()) == 50 and len(grade_counts) == 26,
         "complete insertion-sector H8 affine grade census")
    need(h8_affine_grade(beta1) == (1, 0), "e1-e8 has H8 grade 1")
    need(h8_affine_grade(beta2) == (2, 0), "e2-e8 has H8 grade 2")

    delta_counts = Counter()
    for root in insertion:
        base, charge = h8_affine_grade(root)
        next_base, next_charge = h8_affine_grade(carrier_cycle(root))
        delta_counts[(next_base-base, next_charge-charge)] += 1
    need(delta_counts == Counter({(Fraction(1), Fraction(0)): 40,
                                  (Fraction(-4), Fraction(0)): 10}),
         "exact [H8,W5] support is +1 on40 and -4 on10 roots")
    need(all(delta != 0 for delta, _ in delta_counts),
         "no insertion root lies in the H8 commutant of the cycle")

    raw_signs = Counter()
    compact_coefficients = Counter()
    for root in insertion:
        beta = physical(root)
        minus_beta = tuple(-x for x in beta)
        raw = response.vertex_response((minus_beta, beta), (2, 1))
        sign = response.adjoint_sign(beta)
        need(raw in (-1, 1), "raw opposite-root two-point coefficient is signed unit")
        need(raw == sign, "raw two-point phase equals the actual adjoint sign")
        raw_signs[int(raw)] += 1
        compact_coefficients[int(sign*raw)] += 1
    need(raw_signs == Counter({1: 26, -1: 24}),
         "signed Gram replay has exact 26/24 split")
    need(compact_coefficients == Counter({1: 50}),
         "compact adjoint cancels every raw cocycle phase")

    # The already read original clock is central on the carrier modules.  The
    # congruence check records the typed distinction from the slot cycle.
    x5 = (-2, -2, -2, 3, 3)
    need({x % 5 for x in x5} == {3}, "U5 is scalar zeta5^3 on the carrier 5")
    wedge2_exponents = {(x5[i]+x5[j]) % 5 for i in range(5) for j in range(i+1, 5)}
    need(wedge2_exponents == {1}, "U5 is scalar zeta5 on Lambda2 carrier")
    need(carrier_cycle(beta1) != beta1, "slot cycle is not the scalar U5 action")

    result = {
        "verdict": "EXACT_NONCOMMUTING_LITERAL_CARRIER_CYCLE; ORIGINAL_U5_IS_A_DISTINCT_SCALAR_CLOCK",
        "sources": {
            "actual_chevalley_source": {
                "path": str(REPO / "tfpt_explorer/current_block_geometry.py"),
                "function": "_actual_current_source",
                "sha256": sha256(REPO / "tfpt_explorer/current_block_geometry.py"),
            },
            "existing_a4_sector_checker": {
                "path": str(A4_CHECKER),
                "sha256": sha256(A4_CHECKER),
            },
            "existing_ordered_cocycle_response": {
                "path": str(E8_RESPONSE),
                "sha256": sha256(E8_RESPONSE),
            },
            "original_u5_audit": {
                "path": str(HERE / "P5-OBSERVABLE-IMAGE-AND-INSTRUMENT.txt"),
                "lines": "377-421",
                "sha256": sha256(HERE / "P5-OBSERVABLE-IMAGE-AND-INSTRUMENT.txt"),
            },
        },
        "actual_root_data": {
            "e8_roots": len(roots),
            "insertion_sector": "(5,10)",
            "insertion_roots": len(insertion),
            "contains_opposite_roots": False,
            "family_centralizer_roots_fixed_pointwise_by_cycle": len(family_roots),
        },
        "signed_gram": {
            "raw_opposite_root_coefficients": {"+1": raw_signs[1], "-1": raw_signs[-1]},
            "compact_adjoint_coefficients": {"+1": compact_coefficients[1]},
            "meaning": "raw V_-beta,V_beta phases; compact V_beta^dagger,V_beta Gram is positive",
        },
        "h8_grades": {
            "formula": "h_beta(a)=base+a*Q(beta)",
            "affine_classes": len(grade_counts),
            "distribution": counter_rows(grade_counts),
        },
        "literal_carrier_cycle": {
            "cycle": "e0->e1->e2->e3->e4->e0; e5,e6,e7 fixed",
            "order": 5,
            "preserves_all_50_insertion_roots": True,
            "advances_carrier_5_label": True,
            "witness": {
                "beta1_doubled": list(beta1),
                "beta1": "e1-e8",
                "beta1_h8_grade": "1",
                "beta2_doubled": list(beta2),
                "beta2": "e2-e8",
                "beta2_h8_grade": "2",
                "cycle_beta1_equals_beta2": True,
                "commutator": "[H8,W5]|beta1> = + W5|beta1>, up to the nonzero lift phase",
            },
            "all_50_commutator_grade_differences": {
                "+1": delta_counts[(Fraction(1), Fraction(0))],
                "-4": delta_counts[(Fraction(-4), Fraction(0))],
            },
            "commutes_with_h8_on_insertion_sector": False,
            "phase_independent": "any nonzero signed Chevalley lift has the same nonzero grade differences",
        },
        "original_u5_clock": {
            "definition": "U5=exp(2*pi*i*X/5)",
            "X": list(x5),
            "on_carrier_5": "zeta5^3 I5",
            "on_carrier_Lambda2_10": "zeta5 I10",
            "permutes_slots": False,
            "equals_literal_carrier_cycle": False,
        },
        "scope": {
            "proved": "the literal carrier slot cycle preserves the actual (5,10) sector but fails the H8 commutator on every one of its 50 roots",
            "distinction": "the original U5 degree clock is scalar and cannot be substituted for that slot cycle",
            "not_claimed": "no physical P1 execution, source-clock replacement, or new E8 reader is selected",
        },
    }
    OUTPUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "verdict": result["verdict"],
        "insertion_roots": len(insertion),
        "signed_gram": result["signed_gram"],
        "h8_affine_classes": len(grade_counts),
        "cycle_commutator_grade_differences": result["literal_carrier_cycle"]["all_50_commutator_grade_differences"],
        "witness": result["literal_carrier_cycle"]["witness"],
        "output": str(OUTPUT),
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
