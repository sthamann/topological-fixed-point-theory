#!/usr/bin/env python3
"""Exact joint mark / sewing-time gate for the conditional inverse-Casimir filter.

Requires the already available SymPy package for exact symbolic log and matrix
identities. This is a source-scope check, not a physical-origin derivation.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import sympy as sp


def s_matrix(matrix: sp.Matrix) -> list[list[str]]:
    return [[str(sp.simplify(x)) for x in row] for row in matrix.tolist()]


def circulant_mark_gate() -> dict[str, Any]:
    L = sp.log(2)
    c3 = sp.symbols("c3", positive=True)
    angles = [0, sp.pi/2, sp.pi, 3*sp.pi/2]
    C = sp.Matrix([[0, 1, 2, 1],
                   [1, 0, 1, 2],
                   [2, 1, 0, 1],
                   [1, 2, 1, 0]])
    T = sp.Matrix(4, 4, lambda i, j: int(i == (j+1) % 4))
    E = (sp.eye(4) + T**2)/2
    P0 = sp.ones(4, 4)/4
    Ccos = sp.Matrix([[sp.cos(angles[i]-angles[j]) for j in range(4)]
                      for i in range(4)]).applyfunc(sp.simplify)

    # Full scalar Green kernel uses the specified common H_N subtraction,
    # giving diagonal 0. The filtered R sum shares that same counter.
    Kfull = -L*C/2
    KR_shared = Kfull - Ccos  # includes diagonal -1, since n=1 was removed
    KR_diag0 = KR_shared + sp.eye(4)  # separately subtract H_N-1
    G_orbit = -c3*L*C

    checks = {
        "full_kernel_is_minus_log2_over2_C": sp.simplify(Kfull + L*C/2) == sp.zeros(4),
        "shared_filtered_equals_full_minus_cos": sp.simplify(KR_shared - Kfull + Ccos) == sp.zeros(4),
        "even_projection_annihilates_cos_mode": sp.simplify(E*Ccos*E) == sp.zeros(4),
        "even_shared_filtered_matches_full": sp.simplify(E*KR_shared*E - E*Kfull*E) == sp.zeros(4),
        "2c3_even_match_to_original_orbit_G": sp.simplify(2*c3*E*KR_shared*E - E*G_orbit*E) == sp.zeros(4),
        "own_diag0_reference_shifts_by_E": sp.simplify(E*KR_diag0*E - E*Kfull*E - E) == sp.zeros(4),
        "E_is_projector": sp.simplify(E*E-E) == sp.zeros(4),
        "Ccos_left_halfturn_odd": sp.simplify(T**2*Ccos+Ccos) == sp.zeros(4),
        "Ccos_right_halfturn_odd": sp.simplify(Ccos*T**2+Ccos) == sp.zeros(4),
        "Ccos_simultaneous_shift_invariant": sp.simplify(T**2*Ccos*T**2-Ccos) == sp.zeros(4),
    }
    assert all(checks.values()), checks

    # Shared-counter filtered matrix spectrum: it is not a positive covariance.
    eig = KR_shared.eigenvals()
    eig_serialized = {str(sp.simplify(k)): int(v) for k, v in eig.items()}
    # Independent own-diagonal-zero convention adds I and changes the even block.
    own_eig = KR_diag0.eigenvals()
    own_eig_serialized = {str(sp.simplify(k)): int(v) for k, v in own_eig.items()}
    return {
        "angles": ["0", "pi/2", "pi", "3pi/2"],
        "abel_formula_off_diagonal": "sum_{n>=2} cos(n theta)/n = -log(2*abs(sin(theta/2))) - cos(theta)",
        "reference_convention": {
            "common_counter": "Entrywise: K_full,N[i,j]=sum_{n=1}^N cos(n*(theta_i-theta_j))/n - H_N*delta_ij; K_R,N[i,j]=sum_{n=2}^N cos(n*(theta_i-theta_j))/n - H_N*delta_ij.",
            "counter_status": "This fixed diagonal harmonic subtraction is a comparison reference for the conditional test; it is not derived from P1.",
            "diagonals_in_limit": {"K_full": "0", "K_R_shared_counter": "-1"},
            "same_counter_matrix_identity": "K_R_shared = K_full - C_cos",
            "separate_filtered_diag0": "Subtracting H_N-1 instead adds I: K_R_diag0=K_R_shared+I; its even projection shifts by E.",
        },
        "matrices": {"C_original": s_matrix(C), "C_cos_n1": s_matrix(Ccos),
                     "E_plus": s_matrix(E), "K_full": s_matrix(Kfull),
                     "K_R_shared_counter": s_matrix(KR_shared),
                     "K_R_separate_diag0": s_matrix(KR_diag0),
                     "E_K_R_shared_E": s_matrix(E*KR_shared*E),
                     "E_G_orbit_E": s_matrix(E*G_orbit*E)},
        "identities_verified": checks,
        "even_sector_verdict": "With the shared original H_N counter, 2*c3*E*K_R*E = E*G_orbit*E exactly because E*C_cos*E=0.",
        "separate_diag0_verdict": "If K_R is independently renormalized to diagonal zero, the even block gains E and no longer matches the original projected mark Green matrix.",
        "eigenvalues": {"K_R_shared_counter": eig_serialized,
                        "K_R_separate_diag0": own_eig_serialized},
        "positivity_scope": "The shared-counter finite matrix has negative eigenvalues (and the separate-diagonal-zero matrix is indefinite); neither is presented as a positive finite covariance.",
    }


def sewing_gate() -> dict[str, Any]:
    alpha = sp.symbols("alpha", positive=True)
    x = sp.exp(-alpha)
    plane_infinite = sp.atanh(sp.exp(-2*alpha))/2
    native_cutoff8 = sp.exp(-4*alpha) + sp.Rational(1, 2)*sp.exp(-8*alpha)
    plane_first = sp.Rational(1, 2)*sp.exp(-2*alpha) + sp.Rational(1, 6)*sp.exp(-6*alpha)
    return {
        "native_scalar_mark_port": {
            "state": "p0=(1,1,1,1)/2; unweighted constant four-mark average",
            "mode_energy": "native cylinder Fourier mode n has L0=n",
            "selection_rule": "n=0 mod 4",
            "positive_heated_norm": "4*sum_{m>=1} exp(-4*m*alpha)/(4*m) = -log(1-exp(-4*alpha))",
            "abel_sum_closed_form": "-log(1-exp(-4*alpha))",
            "first_terms": "exp(-4*alpha) + exp(-8*alpha)/2 + exp(-12*alpha)/3 + ...",
            "cutoff_N8": str(native_cutoff8),
            "single_grade2_verdict": "No grade-2 contribution at all; grades 4 and 8 already give two positive terms for cutoff N>=8.",
        },
        "plane_R_field_port": {
            "field": "R(z)=sum_n R_n*z^(-n-2); creation coefficient R_-n appears as z^(n-2)",
            "constant_four_mark_average": "(1/4)*sum_{k=0}^3 i^(k*(n-2)) = 1 iff n=2 mod 4",
            "selection_rule": "n=2 mod 4, forced by the weight-2 phase z^2 relative to the scalar Green rule",
            "filtered_positive_heat_tower": "sum_{m>=0} exp(-(4*m+2)*alpha)/(4*m+2)",
            "closed_form_for_alpha_positive": "atanh(exp(-2*alpha))/2 = log((1+exp(-2*alpha))/(1-exp(-2*alpha)))/4",
            "first_terms": "exp(-2*alpha)/2 + exp(-6*alpha)/6 + exp(-10*alpha)/10 + ...",
            "cutoff_N8": str(plane_first),
            "single_grade2_verdict": "The grade-6 coefficient is strictly positive, so the full tower is not a fixed constant times exp(-2*alpha).",
            "unfiltered_R_norm_note": "Before applying s_n, the positive R-mode norm tower has weights kappa_(4m+2)=8*(4m+2)*((4m+2)^2-1), also with multiple grades.",
        },
        "comparison": {
            "scalar_mark_character": "0 mod 4; starts at grade 4",
            "weight_2_field_character": "2 mod 4; starts at grade 2 and continues at 6,10,...",
            "twist_boundary": "The shift is fixed by the weight-2 coefficient z^(n-2); applying the same unshifted P0 rule to scalar Green and R(z) would silently change the field-to-scalar map.",
            "no_single_grade2_for_either_full_tower": True,
        },
        "positive_spectral_measure_scope": "These heated norms are sums of positive per-grade terms. This does not turn the separately renormalized finite Green matrix into a positive covariance.",
    }


def main() -> dict[str, Any]:
    circle = circulant_mark_gate()
    sewing = sewing_gate()
    return {
        "status": "CONDITIONAL_CIRCLE_MATCH; FULL_SEWING_TOWER_NOT_PURE_GRADE2",
        "hypothesis": "Identify four marks with theta_k=k*pi/2, set the sewing time equal to positive alpha, and apply the inverse-Casimir-filtered scalar Green tower.",
        "scope": "This tests only a literal filtered-field evaluation with the full positive tower, the fixed mark projection, and alpha*L0 sewing. It does not identify physical alpha, choose a source vertex, or rule out a separately defined contact port/state.",
        "circle_and_reference_gate": circle,
        "joint_time_and_character_gate": sewing,
        "positive_controls": {
            "rank_one_R_primary_state": {
                "state": "R_-2*Omega in the neutral source",
                "L0_grade": 2,
                "norm_squared": 48,
                "heated_norm": "48*exp(-2*alpha)",
                "identity": "The exact grade-2 source state gives a pure e^(-2*alpha) heat norm when selected by its own rank-one grade projector.",
            },
            "separate_filtered_scalar_green_port": "The full filtered circle/mark tower passes the conditional even-sector Green match with the shared H_N*delta_ij reference; it is a separate scalar port, not the rank-one R-state evaluation.",
            "scope": "These controls show that the rejection concerns using the same literal fixed projection to evaluate the entire positive tower as one pure grade-2 kernel. They leave a common-source multiport model open, but do not select its P_L0=2 projector, filter, coupling, or reference from P1.",
        },
        "operator_homomorphism_gate": "No algebra-homomorphism condition is imposed as a universal requirement on a boundary-to-field port. The heat obstruction here follows from the fixed positive mode weights and mark character under the stated literal-evaluation hypothesis.",
        "selection_firewall": "The circle-sector agreement uses an explicit geometric identification and shared harmonic diagonal reference subtraction. Neither it nor the rank-one positive control derives the mark/field map, P_L0=2 selection, counter, or physical alpha-sewing rule from P1.",
    }


def cli() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outputpath", type=Path, help="write JSON results here; stdout if omitted")
    args = parser.parse_args()
    payload = json.dumps(main(), indent=2, ensure_ascii=False, sort_keys=True) + "\n"
    if args.outputpath:
        args.outputpath.parent.mkdir(parents=True, exist_ok=True)
        args.outputpath.write_text(payload, encoding="utf-8")
    else:
        print(payload, end="")


if __name__ == "__main__":
    cli()
