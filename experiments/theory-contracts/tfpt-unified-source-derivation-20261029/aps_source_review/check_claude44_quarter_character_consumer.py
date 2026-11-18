#!/usr/bin/env python3
"""Exact type check for Claude44's quarter-character source proposal."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
import sympy as sp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
OUT = HERE / "claude44-quarter-character-consumer-results.json"
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


def smat(a: sp.Matrix) -> list[list[str]]:
    return [[str(sp.simplify(x)) for x in a.row(i)] for i in range(a.rows)]


def sseq(a: list[sp.Expr]) -> list[str]:
    return [str(sp.simplify(x)) for x in a]


def dft(sequence: list[sp.Expr]) -> list[sp.Expr]:
    """Claude44 convention: sum_k i^(m k) sequence[k]."""
    return [sp.simplify(sum(sp.I ** (m * k) * sequence[k] for k in range(4)))
            for m in range(4)]


def main() -> None:
    I = sp.I
    kq = sp.Matrix([[0, 1, 0], [-1, 0, 0], [0, 0, 1]])
    cs = sp.Matrix([[0, 0, 1], [-1, 0, 0], [0, -1, 0]])
    e02 = sp.zeros(3)
    e02[0, 2] = 1
    e20 = e02.T
    e22 = sp.diag(0, 0, 1)
    require(kq ** 2 == sp.diag(-1, -1, 1), "KQ square")
    require(kq ** 4 == sp.eye(3), "KQ order four")

    cycles = [sp.simplify(kq ** (-k) * cs * kq ** k) for k in range(4)]
    require(all(c ** 3 == sp.eye(3) for c in cycles), "four order-three cycles")

    # Exact scalar trace sequences claimed by Claude44.
    even = [sp.trace(c ** 2 * e02) for c in cycles]
    kernel = [sp.trace((c + c ** 2) * e02) for c in cycles]
    require(even == [1, 0, -1, 0], "even-word trace sequence")
    require(kernel == [1, 1, -1, -1], "kernel-word trace sequence")
    require(dft(even) == [0, 2, 0, 2], "even-word Fourier characters")
    require(dft(kernel) == [0, 2 + 2 * I, 0, 2 - 2 * I],
            "kernel-word Fourier characters")

    # The actual endomorphism orbit has only odd quarter characters.
    orbit = [sp.simplify(kq ** (-k) * e02 * kq ** k) for k in range(4)]
    projections = [sp.simplify(sum(
        (I ** (m * k) * orbit[k] for k in range(4)), sp.zeros(3)) / 4)
                   for m in range(4)]
    require(projections[0] == sp.zeros(3) and projections[2] == sp.zeros(3),
            "no even character in E02 orbit")
    require(projections[1] != sp.zeros(3) and projections[3] != sp.zeros(3),
            "both odd characters present")

    # Check the already fixed common marked Bstar map at the same time.
    bstar = sp.Matrix([[0, 0, -1 + I], [1, I, 0], [I, 1, 0]]) / sp.sqrt(2)
    u = sp.diag(1, I, -I)
    p0 = (sp.eye(3) + u ** 2) / 2
    m0 = sp.simplify(bstar * cs * bstar.H)
    require(sp.simplify(bstar * bstar.H) == sp.eye(3), "Bstar unitary")
    require(sp.simplify(bstar * kq.inv() * bstar.H) == u, "marked quarter map")
    require(sp.simplify(bstar * e02 * bstar.H) == m0 * p0,
            "Hardy open word M0P0")
    require(sp.simplify(bstar * e22 * bstar.H) == p0, "marked endpoint P0")

    # Actual whole-96 operator test.  Native kappa=-i gamma_F is used by
    # conjugation on End(H_F); its scalar phase cannot be used as a left-only
    # state-vector weight for this operator-valued coefficient.
    gamma = v252.grading().astype(complex)
    one96 = np.eye(96, dtype=complex)
    pplus = (one96 + gamma) / 2
    e02_96 = BASE.family_02().astype(complex)
    b02 = e02_96 @ pplus
    knat = -1j * gamma
    knat_inv = 1j * gamma
    require(np.allclose(knat @ b02 @ knat_inv, b02),
            "native clock fixes B02 under adjoint action")
    require(np.allclose(knat @ b02, -1j * b02),
            "left action alone has the state-vector phase")
    require(np.allclose(b02 @ knat_inv, 1j * b02),
            "right action cancels the state-vector phase")

    creators = BASE.creators()
    w = BASE.whole96_dictionary()
    weak = {}
    for index, name in ((3, "w1"), (4, "w2")):
        x = w.T @ np.kron(np.eye(3), creators[index]) @ w
        a = x @ b02 - b02 @ x
        predicted = e02_96 @ x @ gamma
        require(np.allclose(a, predicted), f"{name} Hardy weak commutator")
        require(np.linalg.matrix_rank(a) == 16, f"{name} rank")
        require(np.allclose(knat @ a @ knat_inv, -a),
                f"{name} native adjoint character minus one")
        weak[name] = {
            "rank": 16,
            "Ad_kappa_nat_on_A": "-A",
            "quarter_character": 2,
        }

    # Multiplying the family orbit by the gamma-odd adjoint sign merely swaps
    # odd characters 1 and 3.  It creates no invariant character 0.
    combined_average = sp.simplify(sum(
        ((-1) ** k * orbit[k] for k in range(4)), sp.zeros(3)) / 4)
    require(combined_average == sp.zeros(3), "no invariant combined adjoint average")

    result = {
        "scope": (
            "Exact family character and actual H_F=96 operator-adjoint test; "
            "no claim about a separately supplied charged Hom consumer"
        ),
        "KQ_squared": smat(kq ** 2),
        "trace_sequences": {
            "Tr(Ck_squared_E02)": sseq(even),
            "Tr((Ck+Ck_squared)_E02)": sseq(kernel),
        },
        "fourier_sum_i_to_mk": {
            "even_word_m0_to_m3": sseq(dft(even)),
            "kernel_word_m0_to_m3": sseq(dft(kernel)),
        },
        "E02_AdKQ_character_projections": {
            f"m{m}": smat(projections[m]) for m in range(4)
        },
        "actual_96": {
            "Ad_kappa_nat_on_B02": "+B02",
            "left_kappa_nat_on_B02": "-i B02",
            "right_kappa_nat_inverse_on_B02": "+i B02",
            "weak_blocks": weak,
            "combined_adjoint_average": "0",
        },
        "common_marked_frame": {
            "Bstar_KQ_inverse_Bstar_dagger": "diag(1,i,-i)",
            "Bstar_E02_Bstar_dagger": "M0 P0",
            "Bstar_E22_Bstar_dagger": "P0",
        },
        "verdict": (
            "Claude44's scalar trace/Fourier table is exact.  Native kappa "
            "does not provide its claimed compensator on the actual operator: "
            "adjoint action fixes B02 and sends the weak lifted block to -A. "
            "A left-only -i phase would require a separately typed charged Hom "
            "consumer, not the present End(H_F) coefficient."
        ),
    }
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "checks": "PASS",
        "even_word_fourier": result["fourier_sum_i_to_mk"]["even_word_m0_to_m3"],
        "Ad_kappa_nat_B02": "+B02",
        "Ad_kappa_nat_weak": "-A",
        "combined_adjoint_average": "0",
    }, sort_keys=True))


if __name__ == "__main__":
    main()
