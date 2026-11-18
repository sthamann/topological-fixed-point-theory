#!/usr/bin/env python3
"""Exact/numerical audit of Claude46's V_adm Gram response."""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
OUT = HERE / "claude46-gram-response-results.json"


def require(ok: object, label: str) -> None:
    if not bool(ok):
        raise RuntimeError(label)


def invsqrt(a: np.ndarray) -> np.ndarray:
    values, vectors = np.linalg.eigh(a)
    require(np.min(values) > 0, "positive input for inverse square root")
    return vectors @ np.diag(values ** -0.5) @ vectors.conj().T


def main() -> None:
    m0 = np.array([
        [0, -(1 + 1j) / 2, (1 - 1j) / 2],
        [-(1 + 1j) / 2, -1j / 2, -1 / 2],
        [(1 - 1j) / 2, -1 / 2, 1j / 2],
    ], dtype=complex)
    p0 = np.diag([1, 0, 0]).astype(complex)
    one = np.eye(3, dtype=complex)
    require(np.allclose(m0.conj().T @ m0, one), "M0 unitary")
    require(np.allclose(np.linalg.matrix_power(m0, 3), one), "M0 order three")

    y = 2 / 3
    delta = 0.593277280133224001394066545420926568161598543351481820436114
    gap = y - delta
    d0 = y * one - delta * m0
    h0 = d0.conj().T @ d0
    values, vectors = np.linalg.eigh(h0)
    expected_high = y * y + delta * delta + y * delta
    require(np.allclose(values, [gap * gap, expected_high, expected_high]),
            "original positive-branch Gram spectrum")

    f = np.array([1, -1, 1], dtype=complex) / math.sqrt(3)
    pf = np.outer(f, f.conj())
    require(np.allclose(m0 @ f, f), "M0 fixed vector")
    require(abs(abs(f[0]) ** 2 - 1 / 3) < 1e-14, "fixed-vector marked weight")

    def frechet_invsqrt(e: np.ndarray) -> np.ndarray:
        ee = vectors.conj().T @ e @ vectors
        divided = np.empty((3, 3), dtype=float)
        for a in range(3):
            for b in range(3):
                if abs(values[a] - values[b]) < 1e-12:
                    divided[a, b] = -0.5 * values[a] ** -1.5
                else:
                    divided[a, b] = ((values[a] ** -0.5 - values[b] ** -0.5)
                                     / (values[a] - values[b]))
        return vectors @ (divided * ee) @ vectors.conj().T

    def derivative(direction: complex) -> np.ndarray:
        v = -direction * m0 @ p0
        dh = v.conj().T @ d0 + d0.conj().T @ v
        return frechet_invsqrt(dh)

    d_re = derivative(1)
    d_im = derivative(1j)
    leading = 1 / (3 * gap * gap)
    re_ff = np.vdot(f, d_re @ f)
    im_ff = np.vdot(f, d_im @ f)
    require(abs(re_ff - leading) < 2e-11, "real resonant divided difference")
    require(abs(im_ff) < 1e-11, "imaginary direction has zero resonant line")

    phase = np.exp(-5j * math.pi / 12)
    d_phase = derivative(phase)
    phase_ff = np.vdot(f, d_phase @ f)
    require(abs(phase_ff - math.cos(5 * math.pi / 12) * leading) < 2e-11,
            "actual Hardy-phase resonant coefficient")

    # As g=y-delta tends to zero, the divided difference between the fixed
    # line and its complement is O(1/g), not O(1).  Check the exact analytic
    # coefficient with a small-gap replay in the same M0,P0 direction.
    q = pf
    r = one - q
    row_norm = float(np.linalg.norm(q @ p0 @ (m0.conj().T - one) @ r, 2))
    require(abs(row_norm - math.sqrt(2 / 3)) < 2e-14,
            "exact cross-block row norm")

    def derivative_at_gap(test_gap: float, direction: complex) -> np.ndarray:
        test_delta = y - test_gap
        test_d = y * one - test_delta * m0
        test_h = test_d.conj().T @ test_d
        test_values, test_vectors = np.linalg.eigh(test_h)
        v = -direction * m0 @ p0
        test_e = v.conj().T @ test_d + test_d.conj().T @ v
        test_ee = test_vectors.conj().T @ test_e @ test_vectors
        divided = np.empty((3, 3), dtype=float)
        for a in range(3):
            for b in range(3):
                if abs(test_values[a] - test_values[b]) < 1e-12:
                    divided[a, b] = -0.5 * test_values[a] ** -1.5
                else:
                    divided[a, b] = (
                        (test_values[a] ** -0.5 - test_values[b] ** -0.5)
                        / (test_values[a] - test_values[b])
                    )
        return test_vectors @ (divided * test_ee) @ test_vectors.conj().T

    asymptotic_gap = 1e-4
    asymptotic_im = derivative_at_gap(asymptotic_gap, 1j)
    asymptotic_cross = q @ asymptotic_im @ r + r @ asymptotic_im @ q
    expected_cross_op_scaled = math.sqrt(2 / 3) / (3 * y)
    expected_cross_fro_scaled = math.sqrt(2) * expected_cross_op_scaled
    observed_cross_op_scaled = asymptotic_gap * float(np.linalg.norm(asymptotic_cross, 2))
    observed_cross_fro_scaled = asymptotic_gap * float(np.linalg.norm(asymptotic_cross, "fro"))
    require(abs(observed_cross_op_scaled - expected_cross_op_scaled) < 3e-5,
            "cross operator 1/g asymptotic")
    require(abs(observed_cross_fro_scaled - expected_cross_fro_scaled) < 5e-5,
            "cross Frobenius 1/g asymptotic")

    def data(a: np.ndarray, coeff: complex) -> dict[str, object]:
        fro = float(np.linalg.norm(a, "fro"))
        op = float(np.linalg.norm(a, 2))
        line = float(abs(coeff))
        return {
            "rank": int(np.linalg.matrix_rank(a, tol=1e-10)),
            "singular_values": [float(x) for x in np.linalg.svd(a, compute_uv=False)],
            "operator_norm": op,
            "frobenius_norm": fro,
            "fixed_line_coefficient_abs": line,
            "fixed_line_fraction_of_frobenius_amplitude": line / fro,
            "fixed_line_fraction_of_frobenius_square": (line / fro) ** 2,
        }

    # The source coefficient is kappa=(4 X_*/pi)I, X_*=1/(2 cosh(ell/2)).
    hardy_i = complex(
        0.1258739334853158906973132600522979652656,
        -0.4697679151157443151336120497675343216716,
    )
    kappa_numerator = 2 * abs(hardy_i) / math.pi
    require(abs(np.angle(hardy_i) + 5 * math.pi / 12) < 2e-15,
            "actual Hardy phase")
    threshold_ratio = kappa_numerator / (gap / 2)
    ell_threshold = 2 * math.acosh(threshold_ratio)

    # Linearization is distinct from the finite perturbation.  Record a
    # dimensionless in-class control at |kappa|/gap=0.49, just inside the
    # admissibility boundary 1/2, without selecting an ell_NS.
    ratio = 0.49
    kappa_control = ratio * gap * phase
    v_control = -kappa_control * m0 @ p0
    finite = invsqrt((d0 + v_control).conj().T @ (d0 + v_control)) - invsqrt(h0)
    linear = ratio * gap * d_phase
    finite_relative_linear_error = float(
        np.linalg.norm(finite - linear, "fro") / np.linalg.norm(finite, "fro")
    )

    # P0 is the marked family endpoint, not an original structural projector.
    # Even the forbidden substitution Pstr=P0 would not satisfy the contract.
    p0_commutator = p0 @ (m0 @ p0) - (m0 @ p0) @ p0
    require(abs(np.linalg.norm(p0_commutator, 2) - 1) < 1e-14,
            "marked endpoint does not commute with the direction")

    result = {
        "scope": (
            "Unregulated plus-branch Gram derivative for the declared finite "
            "direction -kappa M0P0; not a derivation of V_adm or Pstr typing"
        ),
        "original_gap": gap,
        "gram_eigenvalues": [float(x) for x in values],
        "leading_fixed_line_coefficient": leading,
        "derivatives_per_unit_direction": {
            "real": data(d_re, re_ff),
            "imaginary": data(d_im, im_ff),
            "actual_Hardy_phase_minus_5pi_over_12": data(d_phase, phase_ff),
        },
        "gap_asymptotics": {
            "fixed_line": "Re(kappa)/(3 g^2) Q exactly",
            "cross_operator_norm": "|kappa|/(sqrt(6) g)+O(1) at y=2/3",
            "cross_frobenius_norm": "|kappa|/(sqrt(3) g)+O(1) at y=2/3",
            "cross_row_norm_exact": row_norm,
            "numeric_test_gap": asymptotic_gap,
            "numeric_g_times_cross_operator_norm": observed_cross_op_scaled,
            "numeric_g_times_cross_frobenius_norm": observed_cross_fro_scaled,
        },
        "hardy_norm_condition_if_direction_is_Vadm": {
            "kappa_magnitude": (
                f"{kappa_numerator:.16g}/cosh(ell_NS/2)"
            ),
            "required_strict_ell_NS_lower_bound": ell_threshold,
            "uses_original_epsilon_f": "2/3-delta_ph",
        },
        "finite_vs_linear_control": {
            "kappa_over_gap": ratio,
            "finite_response_frobenius_norm": float(np.linalg.norm(finite, "fro")),
            "linear_response_frobenius_norm": float(np.linalg.norm(linear, "fro")),
            "relative_linearization_error": finite_relative_linear_error,
        },
        "projector_typing": {
            "P0_is_not_Pstr": True,
            "norm_commutator_P0_M0P0": float(np.linalg.norm(p0_commutator, 2)),
            "plus_C6_branch_is_not_original_P_Seam_plus": True,
        },
        "verdict": (
            "The Frechet resonance and the conditional ell_NS threshold are "
            "correct.  Claude46's 62.33 and 8.25 are Frobenius norms; the "
            "operator norms are different.  The actual complex Hardy direction "
            "has only 89.03% fixed-line Frobenius amplitude, not 99.3%.  The "
            "remainder after the g^-2 fixed-line term is O(1/g), not O(1). "
            "The calculation remains conditional because the source does not set "
            "V_adm=-kappa M0P0 and no commutation with the actual Pstr is shown."
        ),
    }
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "checks": "PASS",
        "leading": leading,
        "real_frobenius": result["derivatives_per_unit_direction"]["real"]["frobenius_norm"],
        "imag_frobenius": result["derivatives_per_unit_direction"]["imaginary"]["frobenius_norm"],
        "actual_phase_fixed_line_fraction": result["derivatives_per_unit_direction"]["actual_Hardy_phase_minus_5pi_over_12"]["fixed_line_fraction_of_frobenius_amplitude"],
        "ell_threshold": ell_threshold,
    }, sort_keys=True))


if __name__ == "__main__":
    main()
