#!/usr/bin/env python3
"""Exact general-family KO6 and weak-doublet test for the Claude42 proposal."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
OUT = HERE / "claude42-general-family-ko-results.json"
sys.path.insert(0, str(ROOT / "verification"))
import v252_full_finite_triple as v252  # noqa: E402

BASE_PATH = HERE / "check_hardy_bmu_weak_clifford_lift.py"
SPEC = importlib.util.spec_from_file_location("hardy_weak", BASE_PATH)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"cannot load {BASE_PATH}")
BASE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(BASE)


def require(ok: object, label: str) -> None:
    if not bool(ok):
        raise RuntimeError(label)


def family_lift(m: np.ndarray) -> np.ndarray:
    """Lift one 3x3 G-column coefficient identically to all actual 96 slots."""
    require(m.shape == (3, 3), "family matrix shape")
    out = np.zeros((v252.NH, v252.NH), complex)
    for sector in (0, 1):
        offset = sector * v252.NP
        for name, chirality, color, source_generation in v252.P:
            source = offset + v252.IDX[(name, chirality, color, source_generation)]
            for target_generation in range(3):
                target = offset + v252.IDX[(name, chirality, color, target_generation)]
                out[target, source] = m[target_generation, source_generation]
    return out


def family_block(d: np.ndarray, sector: int, row: tuple[str, str],
                 col: tuple[str, str], color: int = 0) -> np.ndarray:
    """Return the 3x3 family block for one fixed carrier matrix entry."""
    offset = sector * v252.NP
    rows = [offset + v252.IDX[(row[0], row[1], color, g)] for g in range(3)]
    cols = [offset + v252.IDX[(col[0], col[1], color, g)] for g in range(3)]
    return d[np.ix_(rows, cols)]


def comm(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return a @ b - b @ a


def main() -> None:
    c = BASE.creators()
    w = BASE.whole96_dictionary()
    gamma = v252.grading().astype(int)
    kswap = v252.Kswap().astype(int)
    pplus = (np.eye(96, dtype=int) + gamma) // 2
    eye3 = np.eye(3, dtype=int)
    x1 = w.T @ np.kron(eye3, c[3]) @ w
    x2 = w.T @ np.kron(eye3, c[4]) @ w

    # First prove the formula for a generic non-symmetric complex matrix.
    m_generic = np.array([
        [1 + 2j, 3 - 1j, 0.5j],
        [-2 + 0.25j, 0.75 - 3j, 4],
        [1.5, -0.5 + 2j, -1j],
    ], complex)
    m96 = family_lift(m_generic)
    b = m96 @ pplus
    a0 = comm(x1, b)
    require(np.allclose(a0, m96 @ x1 @ gamma), "generic commutator")
    d_adj = a0 + a0.conj().T
    d_ko = d_adj + kswap @ d_adj.conj() @ kswap
    s_generic = m_generic + m_generic.T
    s96 = family_lift(s_generic)
    predicted = s96 @ x1 @ gamma - s96.conj().T @ x1.T @ gamma
    require(np.allclose(d_ko, predicted), "generic KO symmetric-part formula")
    require(np.allclose(d_ko, d_ko.conj().T), "generic KO selfadjoint")
    require(np.allclose(kswap @ d_ko.conj() @ kswap, d_ko), "generic KO J-real")
    require(np.allclose(gamma @ d_ko + d_ko @ gamma, 0), "generic KO odd")

    # Claude42's quoted two-entry matrix.  The second value is available only
    # at the displayed decimal precision and its overall sign is uncalibrated.
    a = 0.1258739334853158906973132600522979652656 - 0.4697679151157443151336120497675343216716j
    b21 = -0.421375 + 1.572594j
    m = np.zeros((3, 3), complex)
    m[0, 2] = a
    m[2, 1] = b21
    s = m + m.T
    r = float(np.sqrt(abs(a) ** 2 + abs(b21) ** 2))
    require(np.allclose(np.linalg.svd(m, compute_uv=False),
                        sorted((abs(a), abs(b21), 0), reverse=True)),
            "directed M singular values")
    require(np.allclose(np.linalg.svd(s, compute_uv=False), (r, r, 0)),
            "KO family singular values")
    require(np.allclose(m @ m.conj().T,
                        np.diag((abs(a) ** 2, 0, abs(b21) ** 2))),
            "naive up Gram is diagonal")
    require(np.allclose(m.conj().T @ m,
                        np.diag((0, abs(b21) ** 2, abs(a) ** 2))),
            "naive down Gram is diagonal")

    d_two = family_lift(s) @ x1 @ gamma - family_lift(s).conj().T @ x1.T @ gamma
    singular_d = np.linalg.svd(d_two, compute_uv=False)
    require(np.count_nonzero(singular_d > 1e-10) == 64, "full KO rank 64")
    require(np.allclose(singular_d[:64], r), "full KO repeated singular value")
    require(np.allclose(singular_d[64:], 0), "full KO zero singular values")

    # Exact four-sector family blocks.  LR means rows L, columns R.
    sd = s.conj().T
    expected = {
        "w1": {
            (0, ("nu", "L"), ("e", "R")): -s,
            (0, ("e", "L"), ("nu", "R")): -sd,
            (0, ("u", "L"), ("d", "R")): s,
            (0, ("d", "L"), ("u", "R")): sd,
            (1, ("nu", "L"), ("e", "R")): -sd,
            (1, ("e", "L"), ("nu", "R")): -s,
            (1, ("u", "L"), ("d", "R")): sd,
            (1, ("d", "L"), ("u", "R")): s,
        },
        "w2": {
            (0, ("nu", "L"), ("nu", "R")): sd,
            (0, ("e", "L"), ("e", "R")): -s,
            (0, ("u", "L"), ("u", "R")): -sd,
            (0, ("d", "L"), ("d", "R")): s,
            (1, ("nu", "L"), ("nu", "R")): s,
            (1, ("e", "L"), ("e", "R")): -sd,
            (1, ("u", "L"), ("u", "R")): -s,
            (1, ("d", "L"), ("d", "R")): sd,
        },
    }
    for name, x in (("w1", x1), ("w2", x2)):
        d = family_lift(s) @ x @ gamma - family_lift(s).conj().T @ x.T @ gamma
        for (sector, row, col), value in expected[name].items():
            require(np.allclose(family_block(d, sector, row, col), value),
                    f"{name} {sector} {row} {col} family block")

    # The actual weak raising operator maps the w2 creator to w1.  A family
    # coefficient is a spectator of this doublet action, hence two arbitrary
    # coefficients on its two components violate the fixed action.
    tplus = w.T @ np.kron(eye3, c[3] @ c[4].T) @ w
    require(np.array_equal(comm(tplus, x2), x1), "weak raising x2 to x1")
    require(np.count_nonzero(comm(tplus, x1)) == 0, "weak highest component")
    m1 = family_lift(np.diag((1, 2, 3)))
    m2 = family_lift(np.diag((1, 4, 3)))
    h1 = m1 @ x1 @ gamma
    h2 = m2 @ x2 @ gamma
    require(np.allclose(comm(tplus, h2) - h1, (m2 - m1) @ x1 @ gamma),
            "component-dependent coefficient covariance residual")
    require(np.linalg.norm(comm(tplus, h2) - h1) > 0,
            "different component coefficients break fixed doublet")
    require(np.allclose(comm(tplus, m1 @ x2 @ gamma), m1 @ x1 @ gamma),
            "common coefficient preserves fixed doublet")

    ratio = b21 / a
    output = {
        "scope": "Actual H_F=96 KO6 completion and fixed weak-doublet action; no family-model fit",
        "general_formula": {
            "input": "B=M P_gamma,+",
            "holomorphic": "A=M X_w gamma_F",
            "D_KO": "(M+M^T)X_w gamma_F-(M^dagger+bar(M))X_w^dagger gamma_F",
            "surviving_family_matrix": "S=M+M^T",
        },
        "quoted_two_entry_candidate": {
            "M": "a E_02+b E_21",
            "a": [a.real, a.imag],
            "b_rounded": [b21.real, b21.imag],
            "b_over_a": [ratio.real, ratio.imag],
            "common_phase_residual_from_rounded_b": abs(ratio.imag),
            "directed_M_singular_values": [abs(b21), abs(a), 0.0],
            "KO_S_singular_values": [r, r, 0.0],
            "full_H96_D_KO_singular_values": {"value": r, "multiplicity": 64, "zero_multiplicity": 32},
            "naive_M_Mdagger_mix": "both Gram matrices are diagonal; at most a permutation, J_CP=0",
            "KO_mix": "the two nonzero singular values are degenerate, so no unique CKM matrix is defined",
        },
        "four_sector_pattern": {
            "common_family_matrix": "S=M+M^T",
            "particle_antiparticle_relation": "fixed J_F conjugation",
            "quark_lepton_relation": "fixed relative minus signs",
            "majorana_block": False,
        },
        "weak_doublet_test": {
            "raising_action": "[T_plus,X_w2]=X_w1, [T_plus,X_w1]=0",
            "same_family_coefficient": "covariant",
            "different_family_coefficients_on_bare_components": "nonzero covariance residual",
            "scope_note": "Separate sector intertwiners or source readers may carry distinct matrices; this test only excludes relabeling the two bare components of one fixed doublet.",
        },
        "source_scope": "The original primitive wedge consumer is present. This checker does not equate the Hardy moment with the whole original Hol(nabla_F*) K_f^str coefficient.",
    }
    OUT.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "checks": "PASS",
        "general_KO_uses_symmetric_part": True,
        "KO_family_singular_values": [round(r, 12), round(r, 12), 0.0],
        "full_KO_rank": 64,
        "unique_CKM_from_candidate": False,
        "different_bare_doublet_component_matrices_covariant": False,
    }, sort_keys=True))


if __name__ == "__main__":
    main()
