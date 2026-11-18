#!/usr/bin/env python3
"""Exact charged particle/hole collar control for the retained alpha=1/4 branch.

The raw modes are the retained Fourier modes ``n=0,-1`` of the same chiral
``A_alpha`` component, with signed energies ``(1/4,-3/4)``.  Filling the
negative mode gives the relative particle/hole energies ``(1/4,3/4)``.
With the explicitly declared control
length ``ell = 4 log(4)``, their collar responses are ``(1/4, 1/64)``.
The graph is interpreted as a finite CAR covariance and is unvectorized with
the oriented face pairing implemented in :mod:`boundary_amplitude`.

This is a conditional source reconstruction.  It neither selects alpha or ell
from P1/P2 nor recovers the absolute partition-function scalar.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import importlib
import json
from pathlib import Path
import sys
from typing import Any

import sympy as sp

from boundary_amplitude import graph_projection, state_from_projection, verify_graph_route
from source_sewing import creator, source_lift


T_TOKEN = -1
ALPHA = sp.Rational(1, 4)
ENERGIES = (ALPHA, 1 - ALPHA)
RAW_SIGNED_ENERGIES = (ALPHA, -(1 - ALPHA))
Q_BASE = sp.Rational(1, 4)
RAW_RESPONSE = sp.diag(Q_BASE, Q_BASE ** -3)
RESPONSE = sp.diag(Q_BASE, Q_BASE ** 3)


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def _search_root() -> Path:
    return _repo_root() / "experiments/theory-contracts/seam-positive-kernel-search-20261015"


def _scalar(value: sp.Expr) -> dict[str, str]:
    value = sp.expand(value)
    return {"re": str(sp.re(value)), "im": str(sp.im(value))}


def _matrix_json(matrix: sp.MatrixBase) -> list[list[dict[str, str]]]:
    return [
        [_scalar(matrix[row, column]) for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ]


def majoranas() -> tuple[sp.Matrix, ...]:
    result: list[sp.Matrix] = []
    for mode in range(2):
        creation = creator(2, mode)
        annihilation = creation.adjoint()
        result.extend((creation + annihilation, sp.I * (creation - annihilation)))
    return tuple(result)


def vacuum_density() -> sp.Matrix:
    result = sp.zeros(4)
    result[0, 0] = 1
    return sp.Matrix(result)


def moment(word: tuple[int, ...] | list[int], transfer: sp.MatrixBase) -> sp.Expr:
    gammas = majoranas()
    operator = sp.eye(4)
    for token in word:
        operator *= transfer if token == T_TOKEN else gammas[token]
    return sp.expand(sp.trace(vacuum_density() * operator))


def charged_fixed_moments() -> list[dict[str, Any]]:
    """Vacuum covariance plus exact charged one- and two-field responses."""
    transfer = source_lift(RESPONSE)
    words: list[tuple[int, ...]] = [
        (left, right) for left in range(4) for right in range(left + 1, 4)
    ]
    words.extend((generator, T_TOKEN, generator) for generator in range(4))
    words.append((0, 2, T_TOKEN, 2, 0))
    return [{"word": list(word), **_scalar(moment(word, transfer))} for word in words]


def run_moment_search(transfer_floor: float, *, solver: str = "SCS") -> dict[str, Any]:
    """Use the existing word-moment search without supplying G, Q, H, or rho."""
    search_root = str(_search_root())
    if search_root not in sys.path:
        sys.path.insert(0, search_root)
    kernel_search = importlib.import_module("kernel_search")
    search = kernel_search.MomentSearch(
        m=2,
        degree=3,
        transfer_floor=transfer_floor,
        transfer_ceiling=1,
        transfer_kind="vacuum_effect",
        fixed_moments=charged_fixed_moments(),
        source_observables={
            "label": "actual charged alpha=1/4 particle-hole collar pair",
            "energies": ["1/4", "3/4"],
            "ell": "4*log(4)",
            "scope": "Only exact vacuum covariance and charged field responses are supplied; no graph, transfer, Hamiltonian, or density matrix is optimizer input.",
        },
    )
    return search.solve(
        objective="feasibility", seed=0, solver=solver, tolerance=1e-7
    )


def verify_charged_pair(*, run_search: bool = False, solver: str = "SCS") -> dict[str, Any]:
    # Start with the signed collar, where the negative-energy transfer grows.
    raw_transfer, raw_route_receipt = verify_graph_route(RAW_RESPONSE)
    expected_raw_transfer = sp.diag(1, Q_BASE, 64, 16)
    if raw_transfer != expected_raw_transfer:
        raise ArithmeticError("raw signed collar lift is wrong")
    if raw_transfer != source_lift(RAW_RESPONSE):
        raise ArithmeticError("raw boundary-state and creator-covariance lifts disagree")
    negative_creator = creator(2, 1)
    particle_hole = negative_creator + negative_creator.adjoint()
    if particle_hole.adjoint() * particle_hole != sp.eye(4):
        raise ArithmeticError("negative-mode particle-hole map is not unitary")
    # The filled-sea eigenvalue is 64.  Conjugation followed by division by
    # that reference value gives the bounded relative particle-hole operator.
    relative_from_raw = (
        particle_hole * raw_transfer * particle_hole.adjoint() / 64
    ).applyfunc(sp.simplify)

    boundary_graph = graph_projection(RESPONSE)
    boundary_omega, boundary_norm = state_from_projection(boundary_graph)
    boundary_density = (
        boundary_omega * boundary_omega.adjoint() / boundary_norm
    ).applyfunc(sp.simplify)
    transfer, route_receipt = verify_graph_route(RESPONSE)
    expected_transfer = sp.diag(
        1, Q_BASE, Q_BASE ** 3, Q_BASE ** 4
    )
    if transfer != expected_transfer:
        raise ArithmeticError("charged graph route did not recover Lambda(R)")
    if relative_from_raw != transfer:
        raise ArithmeticError("particle-hole sea normalization did not recover the relative Q")

    positive_creator = creator(2, 0)
    number_plus = positive_creator * positive_creator.adjoint()
    number_minus = negative_creator * negative_creator.adjoint()
    raw_number = number_plus + number_minus
    relative_charge = (
        particle_hole * raw_number * particle_hole.adjoint() - sp.eye(4)
    ).applyfunc(sp.simplify)
    if relative_charge != number_plus - number_minus:
        raise ArithmeticError("particle-hole relative charge identity failed")
    raw_hamiltonian = ALPHA * number_plus - (1 - ALPHA) * number_minus
    relative_hamiltonian = (
        particle_hole * raw_hamiltonian * particle_hole.adjoint()
        + (1 - ALPHA) * sp.eye(4)
    ).applyfunc(sp.simplify)
    expected_relative_hamiltonian = sp.diag(0, ALPHA, 1 - ALPHA, 1)
    if relative_hamiltonian != expected_relative_hamiltonian:
        raise ArithmeticError("particle-hole relative energy identity failed")

    thermal_density = transfer.adjoint() * transfer
    thermal_density = (thermal_density / sp.trace(thermal_density)).applyfunc(sp.simplify)
    expected_density = sp.diag(65536, 4096, 16, 1) / 69649
    if thermal_density != expected_density:
        raise ArithmeticError("Q*Q did not produce the charged thermal density")

    # Extending the same physical collar twice is chronological exterior sewing.
    doubled_response = RESPONSE * RESPONSE
    doubled_transfer, _ = verify_graph_route(doubled_response)
    if transfer * transfer != doubled_transfer:
        raise ArithmeticError("chronological collar sewing failed")
    if doubled_transfer != source_lift(doubled_response):
        raise ArithmeticError("boundary-state and creator-covariance sewing disagree")

    # The actual charged quarter has opposite phases on the two source fields.
    one_particle_quarter = sp.diag(sp.I, -sp.I)
    charged_quarter = source_lift(one_particle_quarter)
    raw_quarter = source_lift(sp.diag(sp.I, sp.I))
    relative_quarter_from_raw = (
        particle_hole * raw_quarter * particle_hole.adjoint() / sp.I
    ).applyfunc(sp.simplify)
    expected_quarter = sp.diag(1, sp.I, -sp.I, 1)
    if charged_quarter != expected_quarter:
        raise ArithmeticError("charged quarter lift is wrong")
    if relative_quarter_from_raw != charged_quarter:
        raise ArithmeticError("raw signed quarter did not become the charged hole quarter")
    creators = (positive_creator, negative_creator)
    expected_phases = (sp.I, -sp.I)
    for field, phase in zip(creators, expected_phases):
        if charged_quarter * field * charged_quarter.adjoint() != phase * field:
            raise ArithmeticError("charged field covariance failed")
    # exp(2 pi i H) on Fock words, evaluated exactly from source quarter energies.
    fock_energies = tuple(relative_hamiltonian[index, index] for index in range(4))
    phase_lookup = {
        sp.Integer(0): sp.Integer(1),
        sp.Rational(1, 4): sp.I,
        sp.Rational(3, 4): -sp.I,
        sp.Integer(1): sp.Integer(1),
    }
    exponential_from_energies = sp.diag(*(phase_lookup[energy] for energy in fock_energies))
    if exponential_from_energies != charged_quarter:
        raise ArithmeticError("W differs from exp(2 pi i H) on the charged Fock words")
    if charged_quarter * transfer != transfer * charged_quarter:
        raise ArithmeticError("charged quarter does not preserve the collar transfer")

    if boundary_density ** 2 != boundary_density or sp.trace(boundary_density) != 1:
        raise ArithmeticError("boundary CAR projection did not yield a pure state")
    expected_first_face = sp.diag(sp.Rational(16, 17), sp.Rational(4096, 4097))
    expected_second_face = sp.diag(sp.Rational(1, 17), sp.Rational(1, 4097))
    if boundary_graph[:2, :2] != expected_first_face:
        raise ArithmeticError("input-dual face covariance is wrong")
    if boundary_graph[2:, 2:] != expected_second_face:
        raise ArithmeticError("output face covariance is wrong")

    proof_path = _repo_root() / "experiments/theory-contracts/seam-full-collar-reconstruction-20261015/PROOF.txt"
    report: dict[str, Any] = {
        "schema": "tfpt.charged-pair-collar-control.v1",
        "ok": True,
        "verdict": "EXACT_CHARGED_PAIR_CHOIRAY_UNDER_DECLARED_CAR_ORIENTATION",
        "declared_source_inputs": {
            "alpha": "1/4",
            "charged_energies": ["alpha=1/4", "1-alpha=3/4"],
            "raw_signed_D_pair": ["1/4", "-3/4"],
            "ell": "4*log(4)",
            "raw_response_R": ["1/4", "64"],
            "response_R": ["1/4", "1/64"],
            "CAR_covariance_interpretation": True,
            "oriented_particle_hole_Choi_pairing": True,
            "absolute_partition_scalar_fixed": False,
            "source_choice_from_P1_P2": False,
        },
        "exact_results": {
            "raw_signed_Q": ["1", "1/4", "64", "16"],
            "filled_sea_reference_eigenvalue": "64",
            "particle_hole_relative_Q": ["1", "1/4", "1/64", "1/256"],
            "Q_ray_spectrum": ["1", "1/4", "1/64", "1/256"],
            "thermal_rho_diag": ["65536/69649", "4096/69649", "16/69649", "1/69649"],
            "input_dual_face_covariance": ["16/17", "4096/4097"],
            "output_face_covariance": ["1/17", "1/4097"],
            "chronological_two_collar_Q": ["1", "1/16", "1/4096", "1/65536"],
            "charged_quarter_W": ["1", "i", "-i", "1"],
            "field_covariance": ["c_plus^* -> i c_plus^*", "c_minus^* -> -i c_minus^*"],
            "W_equals_exp_2pi_i_H_on_source_Fock_words": True,
            "raw_signed_quarter_to_relative_hole_quarter": True,
            "raw_filled_sea": {"charge": "1", "energy": "-3/4", "quarter_phase": "i"},
            "relative_charge_diag": ["0", "1", "-1", "0"],
            "relative_H_diag": ["0", "1/4", "3/4", "1"],
            "common_charge_time_quarter_path": True,
        },
        "exact_matrices": {
            "raw_signed_D_pair": _matrix_json(sp.diag(*RAW_SIGNED_ENERGIES)),
            "raw_response_R": _matrix_json(RAW_RESPONSE),
            "raw_signed_Q": _matrix_json(raw_transfer),
            "negative_mode_particle_hole": _matrix_json(particle_hole),
            "raw_number_Nplus_Nminus": _matrix_json(raw_number),
            "raw_signed_H": _matrix_json(raw_hamiltonian),
            "relative_charge_Nplus_minus_Nminus": _matrix_json(relative_charge),
            "relative_H_after_sea_shift": _matrix_json(relative_hamiltonian),
            "response_R": _matrix_json(RESPONSE),
            "boundary_graph_G": _matrix_json(boundary_graph),
            "Q_ray": _matrix_json(transfer),
            "thermal_rho": _matrix_json(thermal_density),
            "charged_quarter_W": _matrix_json(charged_quarter),
        },
        "route_receipt": route_receipt,
        "raw_signed_route_receipt": raw_route_receipt,
        "fixed_moments": charged_fixed_moments(),
        "quarter_scope": {
            "physical_charged_quarter_proved_exactly": True,
            "moment_search_internal_quarter": "independent CAR gauge invariance; it is not used to prove the opposite charged phases",
        },
        "source_reference": {
            "path": str(proof_path.relative_to(_repo_root())),
            "sha256": sha256(proof_path.read_bytes()).hexdigest(),
            "atlas_source_ids": [
                "premise:actual-chiral-opposite-quarter-domains",
                "premise:charged-car-negative-sea-source",
            ],
            "scope": "The raw pair is the two retained Fourier modes n=0,-1 of the alpha=1/4 chiral A_alpha component, not the two graded components of the doubled collar. The product collar is supplied; alpha and ell are not selected by this control.",
        },
        "proof_boundary": {
            "proved": "raw signed collar lift, the shared charge/time/quarter particle-hole identities with the filled sea, explicit sea-relative normalization, finite exact CAR projection state, oriented Choi operator ray, Q*Q thermal density, chronological exterior sewing, and charged quarter field covariance",
            "conditional": "the full graph is the physical charged CAR boundary covariance and the labeled faces carry the declared oriented particle-hole pairing",
            "not_proved": ["absolute partition normalization", "selection of the sea-relative representation", "selection of alpha or ell from P1/P2", "physical event law", "infinite-dimensional local-net completion"],
        },
    }
    if run_search:
        accepted = run_moment_search(0, solver=solver)
        rejected = run_moment_search(sp.Rational(1, 4), solver=solver)
        report["moment_search_control"] = {
            "floor0_status": accepted.get("solver_status"),
            "floor_quarter_status": rejected.get("solver_status"),
            "floor0_ok": accepted.get("solver_status") in ("optimal", "optimal_inaccurate"),
            "floor_quarter_rejected": rejected.get("solver_status") in ("infeasible", "infeasible_inaccurate"),
            "fixed_moment_count": len(charged_fixed_moments()),
            "no_target_G_Q_H_or_rho_supplied": True,
        }
        if not report["moment_search_control"]["floor0_ok"]:
            raise ArithmeticError("charged source control failed at floor zero")
        if not report["moment_search_control"]["floor_quarter_rejected"]:
            raise ArithmeticError("charged 1/64 response did not reject floor one quarter")
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out",
        type=Path,
        default=Path(__file__).resolve().parent / "results/charged-pair-control.json",
    )
    parser.add_argument("--search-controls", action="store_true")
    parser.add_argument("--solver", choices=("CLARABEL", "SCS"), default="SCS")
    args = parser.parse_args()
    report = verify_charged_pair(run_search=args.search_controls, solver=args.solver)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "ok": report["ok"],
        "verdict": report["verdict"],
        "Q_ray_spectrum": report["exact_results"]["Q_ray_spectrum"],
        "moment_search_control": report.get("moment_search_control"),
        "out": str(args.out),
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
