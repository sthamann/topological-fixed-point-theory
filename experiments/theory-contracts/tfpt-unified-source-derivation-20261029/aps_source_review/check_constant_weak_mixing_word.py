"""Exact targeted check of the first Claude54 constant-background claim.

This imports the actual v252 finite triple.  It compares two down-type Yukawa
matrices with identical singular values and evaluates a neutral closed word
built from the actual charged weak generators and D_F^2.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
sys.path.insert(0, str(REPO / "verification"))

import v252_full_finite_triple as v252  # noqa: E402


def clean(z: complex) -> int | float | list[float]:
    if abs(z.imag) < 1e-12:
        x = float(z.real)
        return int(round(x)) if abs(x - round(x)) < 1e-12 else x
    return [float(z.real), float(z.imag)]


def main() -> int:
    z2 = np.zeros((2, 2), complex)
    z3 = np.zeros((3, 3), complex)
    sx = np.array([[0, 1], [1, 0]], complex)
    sy = np.array([[0, -1j], [1j, 0]], complex)

    # These are the complex charged combinations of the actual quaternionic
    # weak generators in v252.rep.  Wplus maps d_L to u_L on every colour and
    # generation.  They are not newly fitted event matrices.
    tx = v252.rep(0.0, sx, z3)
    ty = v252.rep(0.0, sy, z3)
    wplus = (tx + 1j * ty) / 2
    wminus = wplus.conj().T

    yu = np.diag([1, 2, 3]).astype(complex)
    yd_aligned = np.diag([4, 5, 6]).astype(complex)
    reverse = np.eye(3, dtype=complex)[[2, 1, 0]]
    yd_rotated = reverse @ yd_aligned
    zero = np.zeros((3, 3), complex)

    K = v252.Kswap()
    jwplus = K @ wplus.conj() @ K
    jwminus = jwplus.conj().T
    wplus_ko = wplus + jwplus
    wminus_ko = wminus + jwminus

    def evaluate(yd: np.ndarray) -> dict[str, object]:
        D = v252.full_D(yu, yd, zero, zero, zero)
        H = D @ D
        hu = yu @ yu.conj().T
        hd = yd @ yd.conj().T
        return {
            "yu_singular_values_squared": sorted(clean(x) for x in np.linalg.eigvalsh(hu)),
            "yd_singular_values_squared": sorted(clean(x) for x in np.linalg.eigvalsh(hd)),
            "family_mixed_trace": clean(np.trace(hd @ hu)),
            "particle_charged_word": clean(np.trace(wplus @ H @ wminus @ H)),
            "ko_completed_charged_word": clean(np.trace(wplus_ko @ H @ wminus_ko @ H)),
        }

    aligned = evaluate(yd_aligned)
    rotated = evaluate(yd_rotated)

    # Directly inspect the three colour-zero d_L -> u_L entries.  The other
    # colours are identical, which is the exact factor three in the word.
    charged_entries = []
    for g in range(3):
        i = v252.IDX[("u", "L", 0, g)]
        j = v252.IDX[("d", "L", 0, g)]
        charged_entries.append(clean(wplus[i, j]))

    checks = {
        "charged_generator_is_actual_v252_quaternion_combination": np.linalg.norm(
            wplus - (tx + 1j * ty) / 2
        ) == 0,
        "d_to_u_entries_are_one": charged_entries == [1, 1, 1],
        "same_yu_spectrum": aligned["yu_singular_values_squared"]
        == rotated["yu_singular_values_squared"],
        "same_yd_spectrum": aligned["yd_singular_values_squared"]
        == rotated["yd_singular_values_squared"],
        "mixed_family_trace_changes": aligned["family_mixed_trace"] == 440
        and rotated["family_mixed_trace"] == 280,
        "particle_word_is_three_colour_mixed_trace": aligned["particle_charged_word"]
        == 3 * aligned["family_mixed_trace"]
        and rotated["particle_charged_word"] == 3 * rotated["family_mixed_trace"],
        "particle_word_changes": aligned["particle_charged_word"] == 1320
        and rotated["particle_charged_word"] == 840,
        "ko_completion_doubles_disjoint_particle_word": aligned[
            "ko_completed_charged_word"
        ]
        == 2 * aligned["particle_charged_word"]
        and rotated["ko_completed_charged_word"]
        == 2 * rotated["particle_charged_word"],
        "particle_and_J_particle_weak_actions_are_disjoint": np.linalg.norm(
            wplus @ jwplus
        )
        == 0
        and np.linalg.norm(jwplus @ wplus) == 0,
    }
    checks = {name: bool(value) for name, value in checks.items()}

    result = {
        "verdict": "CONSTANT_FINITE_SPECTRUM_IS_BLIND_BUT_ACTUAL_CHARGED_WEAK_WORD_SEES_RELATIVE_ORIENTATION",
        "source": "verification/v252_full_finite_triple.py",
        "identity": "Tr(W+ D_F^2 W- D_F^2)=3 Tr((Yd Yd*) (Yu Yu*)) on the quark particle block",
        "aligned": aligned,
        "left_rotated_down_sector": rotated,
        "difference": {
            "family_mixed_trace": aligned["family_mixed_trace"]
            - rotated["family_mixed_trace"],
            "particle_charged_word": aligned["particle_charged_word"]
            - rotated["particle_charged_word"],
            "ko_completed_charged_word": aligned["ko_completed_charged_word"]
            - rotated["ko_completed_charged_word"],
        },
        "charged_entries_colour_zero": charged_entries,
        "checks": checks,
        "all_checks_pass": all(checks.values()),
        "scope": [
            "Exact finite v252 operator word at constant Yukawa background.",
            "Shows that the interacting charged weak operator is not determined by the two separate singular spectra.",
            "Does not assert a nonzero coefficient for this word in the unfinished full Gamma_ren caller.",
            "Does not vary or select Yukawa matrices as original dynamical fields.",
        ],
    }

    raw = json.dumps(result, indent=2, sort_keys=True) + "\n"
    out = HERE / "constant-weak-mixing-word-results.json"
    out.write_text(raw, encoding="utf-8")
    print(raw, end="")
    print("sha256", hashlib.sha256(raw.encode()).hexdigest())
    return 0 if result["all_checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
