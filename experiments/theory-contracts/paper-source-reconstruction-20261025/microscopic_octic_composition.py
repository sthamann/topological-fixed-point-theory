#!/usr/bin/env python3
"""Check the actual QWZ modes and their neutral folded Ward dictionary.

The all-mode statements are symbolic identities.  Floating strip calculations
are diagnostics of the existing QWZ source, not a proof of a scaling limit.
No full old suite is declared green: its transitive source pin has drifted.
"""

from __future__ import annotations

import ast
import hashlib
import importlib.util
import itertools
import json
from pathlib import Path
import sys

import numpy as np
import sympy as sp


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "verification"))


def need(condition, message):
    if not condition:
        raise ArithmeticError(message)


def load(name, relative):
    spec = importlib.util.spec_from_file_location(name, REPO / relative)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def sha(relative):
    return hashlib.sha256((REPO / relative).read_bytes()).hexdigest()


def source_nodes(relative, names):
    tree = ast.parse((REPO / relative).read_text())
    found = {}
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and node.name in names:
            found[node.name] = ast.dump(node, include_attributes=False)
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id in names:
                    found[target.id] = ast.dump(node.value, include_attributes=False)
    need(set(found) == set(names), "source symbol absent")
    return found


def main():
    actual_path = "verification/v1033_charged_disorder.py"
    historical_path = "verification/historical_replay/v1033_original.py"
    names = {"strip_at_momentum", "SX", "SY", "SZ", "TX", "TY"}
    need(source_nodes(actual_path, names) == source_nodes(historical_path, names),
         "the actual QWZ strip or its hopping constants changed")
    source = load("paper_actual_qwz", actual_path)
    linear = load("paper_existing_linearization",
                  "experiments/theory-contracts/microscopic-energy-linearization/checker.py")
    charged = load("paper_existing_car",
                   "experiments/theory-contracts/microscopic-charged-car-limit/checker.py")
    native = load("paper_existing_octic",
                  "experiments/theory-contracts/native-octic-seam-origin-20261025/native_folded_octic_origin.py")

    # The first original-CAR identification, before a bosonic dictionary:
    # K = Lcal_-8/8 = sum_r c(r) :a*_(r+8) a_r:.
    r, n = sp.symbols("r n", integer=True)
    alpha = sp.Rational(1, 4)
    coefficient = lambda index: (index + alpha + 4) / 8
    need(sp.expand(coefficient(r + 8 * n) - coefficient(r)) == n,
         "CAR commutator does not give the full current-mode derivation")
    delta = sp.Matrix([(k + alpha) / 8 - sp.Rational(1, 2) for k in range(8)])
    need(sp.Matrix([coefficient(k - 8) for k in range(8)]) == delta,
         "actual sea-crossing pairs do not carry the original delta")
    need(all(-((k + sp.Rational(3, 4)) / 8 - sp.Rational(1, 2))
             == delta[7 - k] for k in range(8)),
         "orientation/counterterm convention mismatch")
    j = sp.symbols("j", integer=True)
    need(sp.expand((-j + alpha) - (sp.Rational(1, 4) - j)) == 0,
         "orientation changes the physical one-particle generator")

    _, orthogonal, positive, fold = native.original_transport()
    a_map = native.clean(positive.T * orthogonal.T * fold)
    b = native.clean(a_map * delta)
    i8_value = sp.expand(sum(z ** 8 for z in b)
        + 14 * sum(b[u] ** 4 * b[v] ** 4
                   for u, v in itertools.combinations(range(4), 2))
        + 168 * sp.prod(z ** 2 for z in b))
    amplitude = sp.simplify(sp.factorial(8) * sp.conjugate(i8_value)
                            / (96 * sp.sqrt(140)))
    expected = sp.sqrt(35) * (16250979 + 56071008 * sp.I) / 128
    need(sp.simplify(8 ** 8 * amplitude - expected) == 0,
         "microscopic source normalization changes the Ward amplitude")

    # Reuse the actual source strip and the existing analytic mode adapter.
    # The old inherited() pin is not bypassed or reported as passed.  These
    # two functions are called explicitly with this identified source, whose
    # relevant AST is checked against the preserved historical source above.
    rows = []
    row = charged.row_spinor()
    for lattice_n in (4096, 8192):
        maxima = dict(raw_generator_residual=0., projected_generator_residual=0.,
                      row_to_polarized_mode=0., wrong_sea_occupation=0.)
        for label in range(-64, 65):
            data = linear.source_strip(lattice_n, label, source)
            p, rho = data["p"], data["rho"]
            k = label - .25
            scale = lattice_n / (2 * np.pi)
            raw = np.linalg.norm((scale * data["h"] + k * np.eye(16)) @ row)
            raw_cap = np.pi * k ** 2 / lattice_n + 2 * np.pi ** 2 * abs(k) ** 3 / (3 * lattice_n ** 2)
            projected = np.linalg.norm((scale * data["h"] + k * np.eye(16)) @ data["j"])
            projected_cap = scale * (2 * rho ** 8 + abs(p) ** 3 / 6)
            distance = np.linalg.norm(data["j"] - row)
            occupation = np.vdot(row, data["negative"] @ row).real
            wrong = abs(occupation - int(label >= 1))
            need(raw <= raw_cap + 1e-10, "actual top-row generator cap failed")
            need(projected <= projected_cap + 1e-10, "polarized generator cap failed")
            need(distance <= 4 * rho + 1e-12, "same-sea mode estimate failed")
            need(wrong <= 9 * rho ** 2 + 1e-12, "actual filled-sea covariance estimate failed")
            for key, value in zip(maxima, (raw, projected, distance, wrong)):
                maxima[key] = max(maxima[key], float(value))
        rows.append({"N": lattice_n, "labels": [-64, 64], "modes_checked": 129,
                     "maxima": maxima, "interval_arithmetic": False})

    report = {
        "status": "PASS_ACTUAL_QWZ_NEUTRAL_OCTIC_COMPOSITION",
        "exact": {
            "source_hopping_and_strip_equal_preserved_original_ast": True,
            "K_commutator_all_modes": "[K,h_(k,-n)]=n*h_(k,-n-1)",
            "K_vacuum_coefficients": [str(z) for z in delta],
            "QWZ_orientation_and_counterterm_identity": True,
            "charge": 0, "original_source_excitation_energy": 64,
            "ward_amplitude_folded": str(amplitude),
            "ward_amplitude_original_Lcal": str(expected),
        },
        "actual_source_diagnostics": rows,
        "source_hashes": {p: sha(p) for p in [actual_path, historical_path,
            "experiments/theory-contracts/microscopic-energy-linearization/checker.py",
            "experiments/theory-contracts/microscopic-charged-car-limit/README.md",
            "experiments/theory-contracts/native-octic-seam-origin-20261025/native_folded_octic_origin.py"]},
        "scope": {
            "fixed_original_QWZ_and_fold": True,
            "new_independent_sources": 0,
            "simultaneous_lattice_UV_limit_proved": False,
            "microscopic_full_stress_unitary_limit_proved": False,
            "P1_source_selection_proved": False,
            "all_old_source_pin_tests_passed": False,
            "complete_TFPT_proved": False,
        },
    }
    (HERE / "composition-results.json").write_text(json.dumps(report, indent=2) + "\n")
    print(report["status"])
    print("Exact original-CAR Ward identification and 258 actual source modes checked.")


if __name__ == "__main__":
    main()
