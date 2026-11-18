#!/usr/bin/env python3
"""Exact KO6 completion of the Hardy Bmu weak-Clifford block on H_F=96."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
OUT = HERE / "hardy-bmu-ko-completion-results.json"
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


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def pair_matrix(p: np.ndarray, q: np.ndarray, rows: list[int], cols: list[int]) -> list[list[str]]:
    """Render exact p*kappa+q*bar(kappa) coefficients."""
    out = []
    for i in rows:
        row = []
        for j in cols:
            a, b = int(p[i, j]), int(q[i, j])
            terms = []
            if a:
                terms.append("kappa" if a == 1 else "-kappa")
            if b:
                terms.append("bar(kappa)" if b == 1 else "-bar(kappa)")
            row.append(" + ".join(terms).replace("+ -", "- ") if terms else "0")
        out.append(row)
    return out


def main() -> None:
    creators = BASE.creators()
    w = BASE.whole96_dictionary()
    e02 = BASE.family_02()
    e20 = e02.T
    family_real = e02 + e20
    family_support = e02 @ e20 + e20 @ e02
    gamma = v252.grading().astype(int)
    kswap = v252.Kswap().astype(int)

    # The positive real Hardy prefactor does not affect any identity.  Use its
    # printed complex I to test the actual non-real phase.
    kappa = complex(
        0.1258739334853158906973132600522979652656,
        -0.4697679151157443151336120497675343216716,
    )
    require(kappa != 0, "actual Hardy phase is nonzero")

    weak = {}
    for k, name in ((3, "w1"), (4, "w2")):
        x = w.T @ np.kron(np.eye(3, dtype=int), creators[k]) @ w
        a0 = e02 @ x @ gamma
        a = kappa * a0
        adj = a.conj().T
        d_adj = a + adj
        j_a = kswap @ a.conj() @ kswap
        j_d_adj = kswap @ d_adj.conj() @ kswap
        d_ko = d_adj + j_d_adj  # epsilon'=+1 in the actual KO6 convention

        predicted_j_a = -np.conjugate(kappa) * e02 @ x.T @ gamma
        predicted_ko = (kappa * family_real @ x @ gamma
                        - np.conjugate(kappa) * family_real @ x.T @ gamma)
        require(np.allclose(j_a, predicted_j_a), f"{name} exact J action")
        require(np.allclose(d_ko, predicted_ko), f"{name} compact KO formula")
        require(np.allclose(d_adj, d_adj.conj().T), f"{name} adjoint selfadjoint")
        require(np.allclose(gamma @ d_adj + d_adj @ gamma, 0), f"{name} adjoint odd")
        require(not np.allclose(j_d_adj, d_adj), f"{name} adjoint is not J-real")
        require(np.allclose(d_ko, d_ko.conj().T), f"{name} KO selfadjoint")
        require(np.allclose(gamma @ d_ko + d_ko @ gamma, 0), f"{name} KO odd")
        require(np.allclose(kswap @ d_ko.conj() @ kswap, d_ko), f"{name} KO J-real")
        require(np.linalg.matrix_rank(d_adj) == 32, f"{name} adjoint rank")
        require(np.linalg.matrix_rank(d_ko) == 64, f"{name} KO rank")
        require(np.allclose(d_ko @ d_ko, abs(kappa) ** 2 * family_support),
                f"{name} KO square")

        # Old-complex-linear geometric C0 acts only on the normal two-space.
        normal_n = np.array([[0, -1], [1, 0]], complex)
        c0 = 1j * np.kron(normal_n, np.eye(96))
        embedded = np.kron(np.eye(2), d_ko)
        require(np.allclose(c0 @ embedded, embedded @ c0), f"{name} C0 finite commutation")

        p_coeff = x @ gamma
        q_coeff = -x.T @ gamma
        sectors = {}
        for sector, offset in (("particle", 0), ("antiparticle", v252.NP)):
            for kind, species in (("lepton", ("nu", "e")), ("quark", ("u", "d"))):
                color = 0
                left = [offset + v252.IDX[(s, "L", color, 0)] for s in species]
                right = [offset + v252.IDX[(s, "R", color, 0)] for s in species]
                sectors[f"{sector}_{kind}_LR"] = pair_matrix(p_coeff, q_coeff, left, right)
                sectors[f"{sector}_{kind}_RL"] = pair_matrix(p_coeff, q_coeff, right, left)

        expected_sectors = {
            "w1": {
                "particle_lepton_LR": [["0", "-kappa"], ["-bar(kappa)", "0"]],
                "particle_lepton_RL": [["0", "-kappa"], ["-bar(kappa)", "0"]],
                "particle_quark_LR": [["0", "kappa"], ["bar(kappa)", "0"]],
                "particle_quark_RL": [["0", "kappa"], ["bar(kappa)", "0"]],
                "antiparticle_lepton_LR": [["0", "-bar(kappa)"], ["-kappa", "0"]],
                "antiparticle_lepton_RL": [["0", "-bar(kappa)"], ["-kappa", "0"]],
                "antiparticle_quark_LR": [["0", "bar(kappa)"], ["kappa", "0"]],
                "antiparticle_quark_RL": [["0", "bar(kappa)"], ["kappa", "0"]],
            },
            "w2": {
                "particle_lepton_LR": [["bar(kappa)", "0"], ["0", "-kappa"]],
                "particle_lepton_RL": [["kappa", "0"], ["0", "-bar(kappa)"]],
                "particle_quark_LR": [["-bar(kappa)", "0"], ["0", "kappa"]],
                "particle_quark_RL": [["-kappa", "0"], ["0", "bar(kappa)"]],
                "antiparticle_lepton_LR": [["kappa", "0"], ["0", "-bar(kappa)"]],
                "antiparticle_lepton_RL": [["bar(kappa)", "0"], ["0", "-kappa"]],
                "antiparticle_quark_LR": [["-kappa", "0"], ["0", "bar(kappa)"]],
                "antiparticle_quark_RL": [["-bar(kappa)", "0"], ["0", "kappa"]],
            },
        }
        require(sectors == expected_sectors[name], f"{name} exact four-sector blocks")

        weak[name] = {
            "J_A_J_inverse": "-bar(kappa) E_02 X_wa^dagger gamma_F",
            "A_adjoint": "-bar(kappa) E_20 X_wa^dagger gamma_F",
            "D_adjoint_rank": 32,
            "D_adjoint_J_real": False,
            "D_KO": "kappa(E_02+E_20)X_wa gamma_F - bar(kappa)(E_02+E_20)X_wa^dagger gamma_F",
            "D_KO_rank": 64,
            "D_KO_square": "|kappa|^2(E_00+E_22) tensor I_32",
            "D_KO_selfadjoint": True,
            "D_KO_gamma_F_odd": True,
            "D_KO_J_F_real": True,
            "four_sector_carrier_blocks": sectors,
        }

    sources = [
        ROOT / "verification/v252_full_finite_triple.py",
        ROOT / "verification/v254_inner_fluctuations.py",
        HERE / "HARDY-BMU-FIRST-NONZERO-ENTRY.txt",
        HERE / "HARDY-BMU-WEAK-CLIFFORD-ODD-LIFT-REVIEW.txt",
        BASE_PATH,
        ROOT / "experiments/theory-contracts/tfpt-unified-source-derivation-20261029/graded_boundary_functor/FINITE-COUPLED.txt",
        ROOT / "experiments/theory-contracts/tfpt-unified-source-derivation-20261029/graded_boundary_functor/charge_transport/DERIVED-PAIRING.txt",
    ]
    result = {
        "scope": "Exact KO6/adjoint completion on actual H_F=96; no source execution",
        "epsilon_prime": 1,
        "kappa_phase_input": [kappa.real, kappa.imag],
        "kappa_note": "The omitted factor 4 X_*/pi is positive real and changes only scale.",
        "C0_note": "Geometric C0=iN tensor I_96 commutes with the finite coefficient in the old complex type; it is distinct from J_F and supplies no finite KO completion.",
        "CAR_parity_note": "gamma_F oddness is finite internal grading; the quadratic CAR background insertion remains CAR-even.",
        "weak_components": weak,
        "source_sha256": {str(path.relative_to(ROOT)): sha256(path) for path in sources},
    }
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "checks": "PASS",
        "D_adjoint_ranks": {name: data["D_adjoint_rank"] for name, data in weak.items()},
        "D_KO_ranks": {name: data["D_KO_rank"] for name, data in weak.items()},
        "J_completion_required": True,
    }, sort_keys=True))


if __name__ == "__main__":
    main()
