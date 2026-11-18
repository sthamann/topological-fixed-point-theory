#!/usr/bin/env python3
"""Small exact checks for the bounded quotient/locality questions.

All arithmetic is integral or Fraction based.  The computation checks finite
counts and root energies; it does not supply a physical P1/P2 derivation.
"""

from fractions import Fraction
from itertools import product, combinations
import json
from pathlib import Path

import sympy as sp


OUT = Path(__file__).with_name("results.json")


def fstr(x):
    x = Fraction(x)
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def q(r):
    return Fraction(sum(x * x for x in r), 8)


def root_local_to_d8(r):
    # With q(r)=sum(r_i^2)/8, the induced pairing against a D8 root
    # e_i +/- e_j is (r_i +/- r_j)/2.  All such pairings are integral iff
    # every coordinate has the same parity.
    return all((r[i] - r[j]) % 2 == 0 and (r[i] + r[j]) % 2 == 0
               for i, j in combinations(range(8), 2))


def e8_roots():
    roots = []
    # D8 integer roots: +/- e_i +/- e_j.
    for i, j in combinations(range(8), 2):
        for si, sj in product((-1, 1), repeat=2):
            p = [0] * 8
            p[i], p[j] = si, sj
            roots.append(("D8", tuple(p)))
    # E8 spinor roots: half-integral coordinates, even number of minus signs.
    for signs in product((-1, 1), repeat=8):
        if sum(s < 0 for s in signs) % 2 == 0:
            roots.append(("spinor", tuple(Fraction(s, 2) for s in signs)))
    return roots


def energy(p, delta):
    norm2 = sum(x * x for x in p)
    dot = sum(Fraction(x) * d for x, d in zip(p, delta))
    return 8 * (Fraction(norm2, 2) + dot)


def quotient_locality_checks():
    total_sum_zero = 0
    local = []
    q_integral_local = 0
    for r in product(range(4), repeat=8):
        if sum(r) % 4 != 0:
            continue
        total_sum_zero += 1
        if root_local_to_d8(r):
            local.append(r)
            q_integral_local += q(r).denominator == 1

    parity_counts = {
        "all_even": sum(all(x % 2 == 0 for x in r) for r in local),
        "all_odd": sum(all(x % 2 == 1 for x in r) for r in local),
    }
    full_dual = (1, 3, 0, 0, 0, 0, 0, 0)
    full_dual_root_pairings = sorted({
        Fraction(full_dual[i] + sign * full_dual[j], 2)
        for i, j in combinations(range(8), 2)
        for sign in (-1, 1)
    })
    return {
        "total_Z4_8_labels_with_sum_0_mod_4": total_sum_zero,
        "expected_total_4_pow_7": 4 ** 7,
        "local_labels_root_test": len(local),
        "local_labels_parity_counts": parity_counts,
        "all_local_q_integral": bool(q_integral_local == len(local)),
        "local_q_integral_count": q_integral_local,
        "full_dual_label": list(full_dual),
        "full_dual_sum_mod_4": sum(full_dual) % 4,
        "full_dual_q": fstr(q(full_dual)),
        "full_dual_q_mod_1": fstr(q(full_dual) % 1),
        "full_dual_is_local": root_local_to_d8(full_dual),
        "full_dual_root_pairing_examples": [fstr(x) for x in full_dual_root_pairings[:6]],
        "checks_pass": (total_sum_zero == 4 ** 7 and len(local) == 256
                        and parity_counts == {"all_even": 128, "all_odd": 128}
                        and q_integral_local == 256 and q(full_dual) % 1 == Fraction(1, 4)
                        and not root_local_to_d8(full_dual)),
    }


def alpha_energy_checks():
    roots = e8_roots()
    assert len(roots) == 240
    root_type_counts = {kind: sum(t == kind for t, _ in roots) for kind in ("D8", "spinor")}
    p_star = tuple([Fraction(1, 2)] * 4 + [Fraction(-1, 2)] * 4)
    cases = {}
    for alpha in (Fraction(1, 4), Fraction(3, 4)):
        mu = [Fraction(k) + alpha for k in range(8)]
        mu = [x / 8 for x in mu]
        delta = [x - Fraction(1, 2) for x in mu]
        root_energies = [(kind, p, energy(p, delta)) for kind, p in roots]
        zero_roots = [(kind, p) for kind, p, e in root_energies if e == 0]
        d8_min = min(e for kind, _, e in root_energies if kind == "D8")
        spin_min = min(e for kind, _, e in root_energies if kind == "spinor")

        sum_abs_delta = sum(abs(x) for x in delta)
        # Coordinatewise convex bounds over the two E8 cosets. For |delta_i|<1/2,
        # n^2/2+delta_i*n has unique integer minimum 0 at n=0 and is strictly
        # positive for every nonzero integer n. The half-integer minimum occurs
        # at +/-1/2 with sign opposite delta_i.
        integer_lattice_lower_bound = Fraction(0)
        integer_zero_vector_energy = energy((0,) * 8, delta)
        nonzero_integer_coordinate_energies_strictly_positive = all(
            abs(d) < Fraction(1, 2) for d in delta
        )
        spinor_lattice_lower_bound = 8 * (Fraction(1) - sum_abs_delta / 2)
        coordinate_convex_global_lower_bound = min(integer_lattice_lower_bound,
                                                   spinor_lattice_lower_bound)
        fermion_gap = 8 * min(min(mu), 1 - max(mu))
        cases[str(alpha)] = {
            "mu": [fstr(x) for x in mu],
            "delta": [fstr(x) for x in delta],
            "all_abs_delta_lt_half": all(abs(x) < Fraction(1, 2) for x in delta),
            "sum_abs_delta": fstr(sum_abs_delta),
            "sum_delta_squared": fstr(sum(x * x for x in delta)),
            "root_counts": root_type_counts,
            "p_star": [fstr(x) for x in p_star],
            "p_star_energy": fstr(energy(p_star, delta)),
            "p_star_source_charge_sum": fstr(sum(p_star)),
            "spinor_zero_root_count": sum(t == "spinor" for t, _ in zero_roots),
            "d8_zero_root_count": sum(t == "D8" for t, _ in zero_roots),
            "d8_min_root_energy": fstr(d8_min),
            "other_spinor_min_root_energy": fstr(min(e for t, _, e in root_energies
                                                       if t == "spinor" and e > 0)),
            "coordinate_convex_integer_coset_lower_bound": fstr(integer_lattice_lower_bound),
            "integer_coset_zero_vector": [0] * 8,
            "integer_coset_zero_vector_energy": fstr(integer_zero_vector_energy),
            "analytic_condition_for_strictly_positive_nonzero_integer_coordinates": nonzero_integer_coordinate_energies_strictly_positive,
            "integer_coordinate_argument": "For n>=1, n^2/2+delta*n=n*(n/2+delta)>0; for n<=-1, n/2+delta<0 and n<0, so the product is positive. Thus the unique integer-coordinate minimum is n=0 with value 0.",
            "integer_coordinate_minimizer": 0,
            "integer_coordinate_minimizer_unique": True,
            "coordinate_convex_spinor_coset_lower_bound": fstr(spinor_lattice_lower_bound),
            "coordinate_convex_global_lower_bound": fstr(coordinate_convex_global_lower_bound),
            "global_minimizers_in_the_two_E8_cosets": ["integer vacuum p=0", "unique spinor minimizer p_star"],
            "original_one_fermion_gap": fstr(fermion_gap),
            "checks_pass": (all(abs(x) < Fraction(1, 2) for x in delta)
                            and sum_abs_delta == 2 and energy(p_star, delta) == 0
                            and sum(p_star) == 0
                            and sum(t == "spinor" for t, _ in zero_roots) == 1
                            and sum(t == "D8" for t, _ in zero_roots) == 0
                            and d8_min > 0
                            and min(e for t, _, e in root_energies if t == "spinor" and e > 0) > 0
                            and integer_zero_vector_energy == 0
                            and nonzero_integer_coordinate_energies_strictly_positive
                            and integer_lattice_lower_bound == 0
                            and coordinate_convex_global_lower_bound == 0
                            and fermion_gap == Fraction(1, 4)),
        }

    # Reversing both q and delta, with p -> -p, leaves the energy unchanged.
    delta_for_sign_check = [Fraction(x) for x in cases["1/4"]["delta"]]
    reverse_sign_invariance = all(
        energy(p, delta_for_sign_check)
        == energy(tuple(-x for x in p), [-x for x in delta_for_sign_check])
        for _, p in roots
    )
    return {
        "e8_root_counts": root_type_counts,
        "total_e8_roots": len(roots),
        "cases": cases,
        "simultaneous_sign_reversal_energy_invariant": reverse_sign_invariance,
        "checks_pass": (len(roots) == 240 and root_type_counts == {"D8": 112, "spinor": 128}
                        and all(c["checks_pass"] for c in cases.values())
                        and reverse_sign_invariance),
        "scope_note": "Finite root enumeration checks root counts/zeros. The global E8-lattice minimum is certified by the recorded coordinatewise convex bounds, not by root enumeration.",
    }


def geometry_checks():
    # Compute the stated source geometry and deck identity symbolically.
    theta = sp.symbols("theta", real=True)
    f = sp.cos(2 * theta) + 2 * sp.sin(2 * theta)
    deck_difference = sp.simplify(f.subs(theta, theta + sp.pi) - f)
    points = {"0": sp.Integer(0), "pi/4": sp.pi / 4, "pi": sp.pi}
    symbolic_values = {name: sp.simplify(f.subs(theta, point))
                       for name, point in points.items()}
    values = {name: Fraction(int(value.p), int(value.q))
              for name, value in symbolic_values.items()}
    deck_invariant = deck_difference == 0
    z8_step_changes_value = values["0"] != values["pi/4"]
    return {
        "f(theta)": "cos(2 theta) + 2 sin(2 theta)",
        "exact_values": {k: fstr(v) for k, v in values.items()},
        "deck_difference_simplified": str(deck_difference),
        "deck_Z2_invariance_f(theta_plus_pi)=f(theta)": deck_invariant,
        "not_invariant_under_Z8_step_pi_over_4": z8_step_changes_value,
        "P2_rank5_neutral_geometry": "f(theta) * I_5",
        "P2_trace_fI5": {k: fstr(5 * v) for k, v in values.items()},
        "P2_per_component_eigenvalue": {k: fstr(v) for k, v in values.items()},
        "charged_wilson_phases_before": ["+i", "-i"],
        "charged_wilson_phases_after_neutral_scalar_geometry": ["+i", "-i"],
        "wilson_phases_untouched": True,
        "checks_pass": (deck_invariant and z8_step_changes_value
                        and values["0"] == values["pi"]
                        and {k: fstr(5 * v) for k, v in values.items()}
                        == {"0": "5", "pi/4": "10", "pi": "5"}
                        and {k: fstr(v) for k, v in values.items()}
                        == {"0": "1", "pi/4": "2", "pi": "1"}),
        "scope_note": "This distinguishes the stated Z2 deck-invariant source function from a Z8 mark quotient; it does not select a physical quotient or derive a gauge action.",
    }


def main():
    checks = {
        "1_Z4_labels_and_D8_locality": quotient_locality_checks(),
        "2_alpha_shifted_E8_root_energies": alpha_energy_checks(),
        "3_Z2_deck_vs_Z8_mark_quotient": geometry_checks(),
    }
    all_pass = all(c["checks_pass"] for c in checks.values())
    result = {
        "verdict": "PASS_BOUNDED_EXACT_CHECKS_ONLY" if all_pass else "FAIL",
        "categories": {
            "finite_integer_rational_checks": "exact",
            "trigonometric_values": "computed_exactly_with_sympy",
            "Z4_label_count_interpretation": "KLMVACUUM_REFERENCEV2Z_only; not established as the dual-index label space of the actual quarter-twisted source",
            "actual_quarter_twisted_source_global_sector_or_complement_transport": "not_checked",
            "global_E8_lattice_minimum": "analytic_coordinate_convexity_certificate_with_rational_bounds; vacuum and unique spinor minimizer both attain zero",
            "physical_quotient_selection": "not_established",
            "P1_P2_source_derivation": "not_established",
        },
        "assumptions": [
            "The 16384 count is for the stipulated KLMVACUUM_REFERENCEV2Z label set Z4^8 with total coordinate sum 0 modulo 4; it is not asserted to be the dual-index label space of the actual quarter-twisted source.",
            "The local-label count and parity split are exact mathematical counts within that stipulated KLM reference label set; global-sector/complement transport for the actual quarter-twisted source was not checked.",
            "Quadratic form is q(r)=sum(r_i^2)/8; D8 root locality is tested by integrality of all induced pairings (r_i +/- r_j)/2.",
            "For the E8 energy check, the stated lattice consists of the D8 integer-even-sum coset and the even-parity half-integer spinor coset.",
            "Alpha cases are exactly 1/4 and 3/4, with k=0,...,7 and mu=(k+alpha)/8.",
            "The source-geometry check uses the given f(theta), the deck shift theta->theta+pi, and scalar neutral action f I_5; it assumes no gauge transformation on the charged Wilson phases.",
        ],
        "source_references": [
            {"path": "tfpt_1_architecture_e8.tex", "lines": "167-188", "ids": ["AX.P1.01", "AX.P2.01"], "use": "P1/P2 are declared seam and carrier inputs."},
            {"path": "tfpt_1_architecture_e8.tex", "lines": "1461-1496", "ids": ["QGEO.REALIZE.01"], "use": "Four-mark reduction summary; carrier realization remains open/conditional."},
            {"path": "_archive/tfpt-45/source_extracts/01_boundary_kernel_source.tex", "lines": "42-55;94-149", "ids": ["QGEO.REALIZE.01"], "use": "Ordinary geometric two-sheet quotient from boundary polarization."},
            {"path": "tfpt_research_contracts.tex", "lines": "10634-10656;13492-13515", "ids": ["SEAM.EQUIV.01", "CELEST.WP5D.01"], "use": "Conditional local E8 extension and finite versus continuum Haag-duality scope."},
            {"path": "tfpt_research_contracts.tex", "lines": "12633-12660", "ids": ["SEAM.BULK4D.RECON.01"], "use": "Four-dimensional local-algebra reconstruction remains an open contract."},
            {"path": "verification/status_ledger.csv", "lines": "2;256", "ids": ["AX.P1.01", "QGEO.MARKS.01"], "use": "P1 remains declared; raw-seam four-mark construction remains open."},
            {"path": "experiments/theory-contracts/v16-operator-pullback-20261014/contract_index.json", "lines": "1-119", "ids": ["contract:v16-operator-pullback-20261014"], "use": "The stipulated eightfold physical quotient/localization and its Pin/gauge action remain missing derivations."},
        ],
        "checks": checks,
        "complete_TFPT_solution": False,
        "raw_P1P2_derivation": False,
    }
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "verdict": result["verdict"],
        "check_pass_counts": {
            "passed": sum(bool(c["checks_pass"]) for c in checks.values()),
            "total": len(checks),
        },
        "data": {
            "labels_total": checks["1_Z4_labels_and_D8_locality"]["total_Z4_8_labels_with_sum_0_mod_4"],
            "labels_local": checks["1_Z4_labels_and_D8_locality"]["local_labels_root_test"],
            "alpha_1_4_sum_abs_delta": checks["2_alpha_shifted_E8_root_energies"]["cases"]["1/4"]["sum_abs_delta"],
            "alpha_3_4_sum_abs_delta": checks["2_alpha_shifted_E8_root_energies"]["cases"]["3/4"]["sum_abs_delta"],
            "Z2_vs_Z8_values": checks["3_Z2_deck_vs_Z8_mark_quotient"]["exact_values"],
        },
        "results_file": str(OUT),
        "complete_TFPT_solution": False,
        "raw_P1P2_derivation": False,
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
