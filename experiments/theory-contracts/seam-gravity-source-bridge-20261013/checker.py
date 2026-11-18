#!/usr/bin/env python3
"""Bounded source-native diagnostic for the Dorau–Much seam/gravity bridge.

This is an experiment/theory contract. It audits the finite v323 construction and
symbolically checks two coefficient-selection gaps. It does not establish a local
conformal net, a microscopic replica trace, or a complete TFPT gravity derivation.
"""
from __future__ import annotations

import hashlib
import ast
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import sympy as sp

ROOT = Path("/Users/stefanhamann/Projekte/tfpt-theoryv4")
VERIFY = ROOT / "verification"
HERE = Path(__file__).resolve().parent
SOURCES = {
    "v323": VERIFY / "v323_bw_geometric_modular.py",
    "v73": VERIFY / "v73_k_c3_half.py",
    "v358": VERIFY / "v358_grav_entropy_equilibrium.py",
}


def load_v323():
    sys.path.insert(0, str(VERIFY))
    spec = importlib.util.spec_from_file_location("tfpt_v323_native", SOURCES["v323"])
    if spec is None or spec.loader is None:
        raise RuntimeError("could not load the source-native v323 helper module")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run_v323_diagnostic():
    v323 = load_v323()
    modes = np.arange(-12, 13)
    L = np.diag(modes.astype(float))
    rho = np.diag(np.power(1j, modes))
    C = v323.rotation_covariant_C(modes)
    K = v323.modular_hamiltonian(C)
    errs = {
        "[K,L]": float(v323.comm_norm(K, L)),
        "[rho,K]": float(v323.comm_norm(rho, K)),
        "[rho,C]": float(v323.comm_norm(rho, C)),
        "diagonal_in_L_basis": bool(np.allclose(K, np.diag(np.diag(K)))),
        "K_equals_abs_L_max_error": float(np.max(np.abs(np.diag(K).real - np.abs(modes)))),
    }
    # Bilinear transport across the two adjacent unit mode gaps around n=0.
    # The signs differ, so these modes do not share one scalar quarter-turn phase.
    kdiag = np.diag(K).real
    steps = [(-1, 0), (0, 1)]
    transport = []
    for m, n in steps:
        dm = float(kdiag[n + 12] - kdiag[m + 12])
        transport.append({
            "m": m, "n": n, "mode_gap": n - m, "K_gap": dm,
            "modular_phase_t_pi_over_2": {
                "real": float(np.exp(1j * np.pi * dm / 2).real),
                "imag": float(np.exp(1j * np.pi * dm / 2).imag),
            },
        })
    fixed_rotation_phase = {"real": 0.0, "imag": 1.0}
    rigid_residual = float(max(np.hypot(row["modular_phase_t_pi_over_2"]["real"],
                                         row["modular_phase_t_pi_over_2"]["imag"] - 1.0) for row in transport))

    # The best scalar affine relation K=aL+b still has nonzero residual.
    a, b = np.polyfit(modes.astype(float), kdiag, 1)
    affine_residual = float(np.linalg.norm(kdiag - (a * modes + b)))
    opposite_phases = bool(abs(transport[0]["K_gap"] + 1.0) < 1e-9 and abs(transport[1]["K_gap"] - 1.0) < 1e-9
                           and transport[0]["modular_phase_t_pi_over_2"]["imag"] < -0.999
                           and transport[1]["modular_phase_t_pi_over_2"]["imag"] > 0.999)
    checks_pass = bool(errs["[K,L]"] < 1e-9 and errs["[rho,K]"] < 1e-9
                       and errs["[rho,C]"] < 1e-9 and errs["K_equals_abs_L_max_error"] < 1e-9
                       and errs["diagonal_in_L_basis"] and opposite_phases
                       and affine_residual > 1.0)
    return {
        "pass": checks_pass,
        "criteria": errs,
        "adjacent_bilinear_transport": transport,
        "single_rigid_rotation_phase": fixed_rotation_phase,
        "rigid_rotation_max_phase_residual": rigid_residual,
        "least_squares_K_equals_aL_plus_b": {"a": float(a), "b": float(b), "residual_l2": affine_residual},
        "scope": "finite mode diagnostic only; does not exclude doubled/chiral geometric actions or prove local interval Bisognano-Wichmann",
    }


def run_symbolic_coefficient_diagnostics():
    g, J, lam, A, n, k, z, m2, c3 = sp.symbols(
        "g J lambda A n k z m2 c3", nonzero=True, real=True
    )
    z0 = sp.Function("Z0")
    q = sp.Function("Q")
    before = (z0(g, J) * sp.exp(k * q(g))) / (z0(g, 0) * sp.exp(k * q(g)))
    after = (z0(g, J) * sp.exp((k + lam) * q(g))) / (z0(g, 0) * sp.exp((k + lam) * q(g)))
    normalized_correlation_shift = sp.simplify(after - before)

    replica_delta_q = 4 * sp.pi * (1 - n) * A
    replica_before = k * replica_delta_q
    replica_after = (k + lam) * replica_delta_q
    entropy_unshifted = sp.simplify(-sp.diff(replica_before, n).subs(n, 1))
    entropy_shifted = sp.simplify(-sp.diff(replica_after, n).subs(n, 1))
    entropy_coefficient_shift = sp.simplify(entropy_shifted - entropy_unshifted)

    G_times_m2 = m2 / (16 * sp.pi * (k + z * m2))
    d_Gm2_dz = sp.simplify(sp.diff(G_times_m2, z))
    einstein_coeff_general = sp.Rational(1, 2) / k
    einstein_coeff_v73_substitution = sp.simplify(einstein_coeff_general.subs(k, c3 / 2))
    symbolic_pass = bool(
        normalized_correlation_shift == 0
        and sp.simplify(entropy_shifted - 4 * sp.pi * (k + lam) * A) == 0
        and entropy_coefficient_shift == 4 * sp.pi * lam * A
        and sp.simplify(einstein_coeff_v73_substitution - 1 / c3) == 0
        and sp.simplify(d_Gm2_dz + m2**2 / (16 * sp.pi * (k + z * m2)**2)) == 0
        and d_Gm2_dz != 0
    )
    return {
        "fixed_metric_normalized_correlation_shift": str(normalized_correlation_shift),
        "conical_convention": "Q_n - n Q_1 = 4*pi*(1-n)*A",
        "entropy_from_shifted_replica_functional": str(entropy_shifted),
        "entropy_before_shift": str(entropy_unshifted),
        "coefficient_shift": str(entropy_coefficient_shift),
        "microscopic_trace_caveat": "multiplying the replica functional by exp(lambda Q[g]) changes the replica completion; this does not alter the entropy of the same fixed density matrix trace",
        "four_dimensional_scaling": {"lambda": "z*m2", "G_times_m2": str(G_times_m2), "d_dz": str(d_Gm2_dz)},
        "EH_variation_convention": "delta S_m = -(1/2) int sqrt(g) T_ab delta g^ab",
        "general_einstein_coefficient": str(einstein_coeff_general),
        "after_v73_k_equals_c3_over_2": str(einstein_coeff_v73_substitution),
        "interpretation": "the variation supplies the factor 1/2 conditional on the action normalization k; substituting k=c3/2 yields 1/c3 but does not independently select k",
        "symbolic_pass": symbolic_pass,
    }


def run_original_verifiers():
    env = os.environ.copy()
    env["PYTHONPATH"] = str(VERIFY) + os.pathsep + env.get("PYTHONPATH", "")
    statuses = {}
    for key in ("v73", "v358", "v323"):
        proc = subprocess.run([sys.executable, str(SOURCES[key])], cwd=ROOT, env=env,
                              text=True, capture_output=True, timeout=120)
        statuses[key] = {
            "exit_status": proc.returncode,
            "python_executable": sys.executable,
            "stdout": proc.stdout.strip(),
            "stderr": proc.stderr.strip(),
        }
    return statuses


def source_provenance():
    citations = {
        "v323": "47-50,74-78,92-114",
        "v73": "34-46",
        "v358": "61-86",
    }
    anchors = {
        "v73": ("k",),
        "v358": ("eta",),
    }
    result = {
        key: {
            "path": str(path),
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "cited_lines": citations[key],
        }
        for key, path in SOURCES.items()
    }
    for key, names in anchors.items():
        tree = ast.parse(SOURCES[key].read_text())
        matches = []
        for node in ast.walk(tree):
            if isinstance(node, ast.Assign) and any(
                isinstance(target, ast.Name) and target.id in names for target in node.targets
            ):
                matches.append({"line": node.lineno, "source": ast.unparse(node)})
        result[key]["assignment_anchors"] = matches
    return result


def main():
    result = {
        "contract": HERE.name,
        "checker_runtime": {"python": sys.executable, "optimization_level": sys.flags.optimize},
        "scope": "source-native research diagnostic; no physical-model addition",
        "v323_diagnostic": run_v323_diagnostic(),
        "coefficient_diagnostics": run_symbolic_coefficient_diagnostics(),
        "original_verifier_replays": run_original_verifiers(),
        "source_provenance": source_provenance(),
        "verdict": "PARTIAL",
        "scope_summary_pass": False,
        "theory_gates_moved": [],
        "t1_t8_closed": [],
        "full_tfpt_solution_claimed": False,
    }
    result["scope_summary_pass"] = bool(
        result["v323_diagnostic"]["pass"]
        and result["coefficient_diagnostics"]["symbolic_pass"]
        and all(v["exit_status"] == 0 for v in result["original_verifier_replays"].values())
    )
    (HERE / "results.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "verdict": result["verdict"],
        "scope_summary_pass": result["scope_summary_pass"],
        "v323": result["v323_diagnostic"],
        "symbolic_pass": result["coefficient_diagnostics"]["symbolic_pass"],
        "verifier_exit_statuses": {k: v["exit_status"] for k, v in result["original_verifier_replays"].items()},
    }, indent=2))
    return 0 if result["scope_summary_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
