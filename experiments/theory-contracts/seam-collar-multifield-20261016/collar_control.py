#!/usr/bin/env python3
"""Exact two-face collar control and conditional multifield sewing test.

The collar graph is source data.  When it is declared to be a finite CAR
covariance, positivity forces its unique Slater boundary ray.  A separately
declared oriented particle-hole pairing unvectorizes that ray to the operator
``Q = Lambda(R)``.  No Hamiltonian or density matrix is supplied to the moment
optimizer, and the operator's overall partition-function scalar is not
recovered from the normalized boundary state.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from functools import lru_cache
from hashlib import sha256
import importlib
import json
from pathlib import Path
import sys
from typing import Any

import sympy as sp


T_TOKEN = -1
MODE_COUNT = 2
B = sp.Rational(1, 4)
T_SCALAR = sp.Rational(1, 4)  # exp(-ell B), with ell = 4 log(4)

SOURCE_REFERENCES = [
    {
        "path": "experiments/theory-contracts/seam-full-collar-reconstruction-20261015/PROOF.txt",
        "lines": "63-143",
        "supports": "two-face graph projector, recovery of B, and the explicitly additional CAR/quasifree/free-dynamics premise",
    },
    {
        "path": "_archive/tfpt-45/source_extracts/01_boundary_kernel_source.tex",
        "lines": "480-520",
        "supports": "operational seed and collar generator; it does not state a Grassmann Gaussian or Wick sewing law",
    },
    {
        "path": "website/public/verification/v160_seam_gaussianity_from_pf.py",
        "lines": "1-24",
        "supports": "the full-correlation-cone Perron premise needed to force Gaussianity is explicitly unproved",
    },
]


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def _search_root() -> Path:
    return _repo_root() / "experiments/theory-contracts/seam-positive-kernel-search-20261015"


def _scalar(value: sp.Expr) -> dict[str, str]:
    value = sp.expand(value)
    return {"re": str(sp.re(value)), "im": str(sp.im(value))}


def _matrix_json(matrix: sp.MatrixBase) -> list[list[dict[str, str]]]:
    return [[_scalar(matrix[row, column]) for column in range(matrix.cols)]
            for row in range(matrix.rows)]


def _gamma_action(generator: int, state: int) -> tuple[int, sp.Expr]:
    mode = generator // 2
    phase: sp.Expr = sp.Integer((-1) ** ((state & ((1 << mode) - 1)).bit_count()))
    if generator % 2:
        phase *= -sp.I if state & (1 << mode) else sp.I
    return state ^ (1 << mode), phase


def gamma_matrices() -> tuple[sp.SparseMatrix, ...]:
    result = []
    for generator in range(2 * MODE_COUNT):
        entries = {}
        for state in range(1 << MODE_COUNT):
            destination, phase = _gamma_action(generator, state)
            entries[destination, state] = phase
        result.append(sp.SparseMatrix(4, 4, entries))
    return tuple(result)


def annihilation_matrices() -> tuple[sp.SparseMatrix, ...]:
    gammas = gamma_matrices()
    return tuple(sp.SparseMatrix((gammas[2 * mode] + sp.I * gammas[2 * mode + 1]) / 2)
                 for mode in range(MODE_COUNT))


def exterior_transfer(graph: sp.MatrixBase, scalar: sp.Rational) -> sp.SparseMatrix:
    """Generate Q=scalar*Lambda(K), K=scalar*(1-G)+scalar^-1*G."""
    graph = sp.Matrix(graph)
    identity = sp.eye(2)
    one_particle = (scalar * (identity - graph) + graph / scalar).applyfunc(sp.expand)
    result = sp.zeros(4)
    result[0, 0] = scalar
    result[1:3, 1:3] = scalar * one_particle
    result[3, 3] = sp.expand(scalar * one_particle.det())
    return sp.SparseMatrix(result)


@dataclass(frozen=True)
class CollarControl:
    graph: sp.Matrix
    cross: sp.Rational
    one_particle_h: sp.Matrix
    many_body_h: sp.Matrix
    transfer: sp.SparseMatrix
    omega: sp.Matrix
    state_density: sp.Matrix
    gammas: tuple[sp.SparseMatrix, ...]
    quasiparticle_creators: tuple[sp.Matrix, ...]

    def moment(self, word: tuple[int, ...] | list[int], *, transfer: sp.MatrixBase | None = None) -> sp.Expr:
        operator = sp.eye(4)
        chosen_transfer = self.transfer if transfer is None else transfer
        for token in word:
            operator *= chosen_transfer if token == T_TOKEN else self.gammas[token]
        return sp.expand(sp.trace(self.state_density * operator))


@dataclass(frozen=True)
class SourceChoiControl:
    """Two identical source modes reconstructed from the boundary graph ray."""

    response: sp.Matrix
    boundary_graph: sp.Matrix
    boundary_state_density: sp.Matrix
    transfer: sp.Matrix
    thermal_density: sp.Matrix
    vacuum_density: sp.Matrix
    gammas: tuple[sp.SparseMatrix, ...]
    route_receipt: dict[str, Any]

    def moment(self, word: tuple[int, ...] | list[int]) -> sp.Expr:
        operator = sp.eye(4)
        for token in word:
            operator *= self.transfer if token == T_TOKEN else self.gammas[token]
        return sp.expand(sp.trace(self.vacuum_density * operator))


@lru_cache(None)
def build_control() -> CollarControl:
    graph = sp.Matrix([[16, 4], [4, 1]]) / 17
    cross = sp.Rational(4, 17)
    one_particle_h = (B * (sp.eye(2) - 2 * graph)).applyfunc(sp.expand)
    many_body_h = sp.zeros(4)
    many_body_h[0, 0] = B
    many_body_h[1:3, 1:3] = B * sp.eye(2) + one_particle_h
    many_body_h[3, 3] = B + sp.trace(one_particle_h)
    transfer = exterior_transfer(graph, T_SCALAR)
    sqrt17 = sp.sqrt(17)
    omega = sp.Matrix([0, 4 / sqrt17, 1 / sqrt17, 0])
    state_density = (omega * omega.adjoint()).applyfunc(sp.expand)
    gammas = gamma_matrices()
    annihilators = annihilation_matrices()
    creators = tuple(matrix.adjoint() for matrix in annihilators)
    # Relative to the occupied graph orbital, these create a hole and a particle.
    hole_creator = (4 * annihilators[0] + annihilators[1]) / sqrt17
    particle_creator = (-creators[0] + 4 * creators[1]) / sqrt17
    return CollarControl(
        graph=graph,
        cross=cross,
        one_particle_h=one_particle_h,
        many_body_h=sp.Matrix(many_body_h),
        transfer=transfer,
        omega=omega,
        state_density=sp.Matrix(state_density),
        gammas=gammas,
        quasiparticle_creators=(sp.Matrix(hole_creator), sp.Matrix(particle_creator)),
    )


@lru_cache(None)
def build_source_choi_control() -> SourceChoiControl:
    """Recover the two-field transfer from the projection state and orientation.

    ``verify_graph_route`` solves the exact CAR annihilation conditions of the
    four-mode boundary projection before unvectorizing.  It does not assume a
    Gaussian Hamiltonian or Wick rule.
    """
    from boundary_amplitude import graph_projection, state_from_projection, verify_graph_route

    response = sp.eye(2) / 4
    boundary_graph = graph_projection(response)
    boundary_omega, boundary_norm = state_from_projection(boundary_graph)
    boundary_state_density = (
        boundary_omega * boundary_omega.adjoint() / boundary_norm
    ).applyfunc(sp.simplify)
    transfer, route_receipt = verify_graph_route(response)
    expected = sp.diag(1, sp.Rational(1, 4), sp.Rational(1, 4), sp.Rational(1, 16))
    if transfer != expected:
        raise ArithmeticError("oriented boundary amplitude did not recover Lambda(R)")
    thermal_density = (transfer.adjoint() * transfer)
    thermal_density = (thermal_density / sp.trace(thermal_density)).applyfunc(sp.simplify)
    vacuum_density = sp.zeros(4)
    vacuum_density[0, 0] = 1
    return SourceChoiControl(
        response=sp.Matrix(response),
        boundary_graph=sp.Matrix(boundary_graph),
        boundary_state_density=sp.Matrix(boundary_state_density),
        transfer=sp.Matrix(transfer),
        thermal_density=sp.Matrix(thermal_density),
        vacuum_density=sp.Matrix(vacuum_density),
        gammas=gamma_matrices(),
        route_receipt=dict(route_receipt),
    )


def adversarial_interacting_transfer(control: CollarControl | None = None) -> sp.SparseMatrix:
    """Change only the two-quasiparticle response from 1/16 to 1/8."""
    control = build_control() if control is None else control
    # The two-quasiparticle vector is the unoccupied graph orbital in the
    # absolute one-particle sector.  Its projector is rational.
    two_particle_projector = sp.zeros(4)
    two_particle_projector[1:3, 1:3] = sp.eye(2) - control.graph
    return sp.SparseMatrix(control.transfer + sp.Rational(1, 16) * two_particle_projector)


def thermal_gns_transfer() -> tuple[sp.Matrix, sp.Matrix, sp.Matrix, sp.Matrix]:
    """Return one-face H, Liouville L, exp(-ell L), and the thermal GNS vector.

    The exponential is evaluated by the exact relation exp(-ell B)=1/4,
    avoiding any floating logarithms.
    """
    one_face_h = sp.diag(0, B)
    liouville = sp.kronecker_product(one_face_h, sp.eye(2)) - sp.kronecker_product(
        sp.eye(2), one_face_h.T
    )
    gns_transfer = sp.diag(1, 1 / T_SCALAR, T_SCALAR, 1)
    thermal_omega = sp.Matrix([4 / sp.sqrt(17), 0, 0, 1 / sp.sqrt(17)])
    return (
        sp.Matrix(one_face_h),
        sp.Matrix(liouville),
        sp.Matrix(gns_transfer),
        thermal_omega,
    )


def source_fixed_moments() -> list[dict[str, Any]]:
    """Vacuum covariance and field responses of the reconstructed Choi ray.

    These are observable moments.  The optimizer never receives ``Q``, ``rho``,
    the graph projection, or a Hamiltonian matrix.
    """
    control = build_source_choi_control()
    words: list[tuple[int, ...]] = []
    # The complete off-diagonal Majorana covariance fixes the source vacuum.
    words.extend((left, right) for left in range(4) for right in range(left + 1, 4))
    # All one-field responses and one independent two-field response.
    words.extend((generator, T_TOKEN, generator) for generator in range(4))
    words.append((0, 2, T_TOKEN, 2, 0))
    return [{"word": list(word), **_scalar(control.moment(word))} for word in words]


def thermal_fixed_moments() -> list[dict[str, Any]]:
    """Exact one-mode thermal-GNS moments distinguishing ceiling four from one."""
    return [
        {"word": [0, 1], "re": "0", "im": "15/17"},
        {"word": [0, T_TOKEN, 0], "re": "8/17", "im": "0"},
        {"word": [1, T_TOKEN, 1], "re": "8/17", "im": "0"},
        {"word": [0, T_TOKEN, T_TOKEN, 0], "re": "1", "im": "0"},
    ]


def _word_algebra():
    search_root = str(_search_root())
    if search_root not in sys.path:
        sys.path.insert(0, search_root)
    return importlib.import_module("word_algebra").WordAlgebra(MODE_COUNT)


def exact_source_moments(max_degree: int = 6) -> list[dict[str, Any]]:
    """Exact rational oracle data for independent replay; not optimizer input."""
    control = build_source_choi_control()
    algebra = _word_algebra()
    result = []
    for word in algebra.words(max_degree):
        value = control.moment(word)
        if sp.re(value).is_Rational is not True or sp.im(value).is_Rational is not True:
            raise ValueError(f"nonrational source moment at {word}: {value}")
        result.append({"word": list(word), **_scalar(value)})
    return result


def verify_control() -> dict[str, Any]:
    control = build_control()
    source_control = build_source_choi_control()
    identity2, identity4 = sp.eye(2), sp.eye(4)
    zero4 = sp.zeros(4)
    if control.graph ** 2 != control.graph or control.graph != control.graph.adjoint():
        raise ValueError("collar graph is not an orthogonal projection")
    if control.graph != sp.Matrix([[16, 4], [4, 1]]) / 17:
        raise ValueError("unexpected collar graph")
    if control.one_particle_h * control.graph != -B * control.graph:
        raise ValueError("occupied graph orbital does not have energy -B")
    if control.one_particle_h * (identity2 - control.graph) != B * (identity2 - control.graph):
        raise ValueError("complementary orbital does not have energy +B")
    if (control.omega.adjoint() * control.omega)[0] != 1:
        raise ValueError("Omega is not normalized")
    if control.transfer * control.omega != control.omega:
        raise ValueError("Q Omega != Omega")
    expected_spectrum = {sp.Integer(1): 1, sp.Rational(1, 4): 2, sp.Rational(1, 16): 1}
    if control.transfer.eigenvals() != expected_spectrum:
        raise ValueError("exterior transfer has the wrong exact spectrum")
    expected_h_spectrum = {sp.Integer(0): 1, sp.Rational(1, 4): 2, sp.Rational(1, 2): 1}
    if control.many_body_h.eigenvals() != expected_h_spectrum:
        raise ValueError("normal-ordered H has the wrong exact spectrum")

    for index, gamma in enumerate(control.gammas):
        if gamma != gamma.adjoint():
            raise ValueError(f"gamma_{index} is not selfadjoint")
        for previous in range(index + 1):
            expected = 2 * identity4 if previous == index else zero4
            if gamma * control.gammas[previous] + control.gammas[previous] * gamma != expected:
                raise ValueError(f"CAR failed at {index},{previous}")

    for creator in control.quasiparticle_creators:
        annihilator = creator.adjoint()
        if annihilator * control.omega != sp.zeros(4, 1):
            raise ValueError("quasiparticle annihilator does not kill Omega")
        if control.transfer * creator != T_SCALAR * creator * control.transfer:
            raise ValueError("free push-through covariance failed")
        if annihilator * control.transfer != T_SCALAR * control.transfer * annihilator:
            raise ValueError("Euclidean adjoint push-through covariance failed")
    first, second = control.quasiparticle_creators
    if control.transfer * second * first != T_SCALAR ** 2 * second * first * control.transfer:
        raise ValueError("two-field exterior composition failed")

    half_transfer = exterior_transfer(control.graph, sp.Rational(1, 2))
    if half_transfer * half_transfer != control.transfer:
        raise ValueError("exterior semigroup composition failed")
    if control.moment((0, T_TOKEN, 0)) != sp.Rational(1, 4):
        raise ValueError("one-field response is not 1/4")
    if control.moment((0, 2, T_TOKEN, 2, 0)) != sp.Rational(1, 16):
        raise ValueError("two-field response is not 1/16")

    adversarial = adversarial_interacting_transfer(control)
    grade_vectors = [control.omega]
    grade_vectors.extend(creator * control.omega for creator in control.quasiparticle_creators)
    if any(adversarial * vector != control.transfer * vector for vector in grade_vectors):
        raise ValueError("adversarial transfer changed grade 0/1 data")
    two_vector = second * first * control.omega
    if adversarial * two_vector != sp.Rational(1, 8) * two_vector:
        raise ValueError("adversarial transfer did not change the grade-2 response")
    if all(adversarial * creator == T_SCALAR * creator * adversarial
           for creator in control.quasiparticle_creators):
        raise ValueError("adversarial transfer unexpectedly obeys full field covariance")

    right_face_occupation = control.graph[1, 1]
    if right_face_occupation != sp.Rational(1, 17):
        raise ValueError("right-face occupation is not 1/17")
    if (1 - right_face_occupation) / right_face_occupation != 16:
        raise ValueError("occupation logit ratio is not 16")

    one_face_h, liouville, gns_transfer, thermal_omega = thermal_gns_transfer()
    if liouville.eigenvals() != {
        sp.Integer(0): 2, B: 1, -B: 1
    }:
        raise ValueError("thermal GNS Liouville spectrum is wrong")
    if gns_transfer.eigenvals() != {
        sp.Integer(1): 2, T_SCALAR: 1, 1 / T_SCALAR: 1
    }:
        raise ValueError("thermal GNS transfer spectrum is wrong")
    if liouville * thermal_omega != sp.zeros(4, 1):
        raise ValueError("thermal GNS Liouvillean does not fix its cyclic vector")
    if gns_transfer * thermal_omega != thermal_omega:
        raise ValueError("thermal GNS transfer does not fix its cyclic vector")
    if all((sp.eye(4) - gns_transfer)[index, index] >= 0 for index in range(4)):
        raise ValueError("thermal GNS transfer was incorrectly classified as a contraction")

    expected_source_q = sp.diag(1, T_SCALAR, T_SCALAR, T_SCALAR ** 2)
    if source_control.transfer != expected_source_q:
        raise ValueError("source graph/Choi route has the wrong exterior transfer")
    expected_boundary_graph = sp.Matrix.vstack(sp.eye(2), sp.eye(2) / 4)
    expected_boundary_graph = (
        expected_boundary_graph
        * (expected_boundary_graph.adjoint() * expected_boundary_graph).inv()
        * expected_boundary_graph.adjoint()
    ).applyfunc(sp.simplify)
    if source_control.boundary_graph != expected_boundary_graph:
        raise ValueError("four-mode boundary graph is wrong")
    if source_control.boundary_state_density ** 2 != source_control.boundary_state_density:
        raise ValueError("boundary CAR state is not pure")
    if sp.trace(source_control.boundary_state_density) != 1:
        raise ValueError("boundary CAR state is not normalized")
    if source_control.thermal_density != sp.diag(256, 16, 16, 1) / 289:
        raise ValueError("Q*Q did not induce the expected two-mode thermal state")
    if expected_source_q ** 2 != sp.diag(
        1, sp.Rational(1, 16), sp.Rational(1, 16), sp.Rational(1, 256)
    ):
        raise ValueError("two-step exterior semigroup identity failed")
    if source_control.boundary_graph[:2, :2] != sp.eye(2) * sp.Rational(16, 17):
        raise ValueError("input-dual boundary CAR covariance is wrong")
    if source_control.boundary_graph[2:, 2:] != sp.eye(2) * sp.Rational(1, 17):
        raise ValueError("output boundary CAR covariance is wrong")

    source_hashes = {}
    for reference in SOURCE_REFERENCES:
        path = _repo_root() / reference["path"]
        source_hashes[reference["path"]] = sha256(path.read_bytes()).hexdigest()
    return {
        "schema": "tfpt.collar-multifield-control.v2",
        "ok": True,
        "verdict": "EXACT_GRAPH_STATE_TO_CHOIRAY_CONTROL_UNDER_DECLARED_CAR_ORIENTATION",
        "declared_inputs": {
            "physical_boundary_branch_alpha": "1/4",
            "lowest_mode_B": "1/4",
            "control_mode_selection": "two declared identical lowest-mode copies; not the physical particle/hole pair",
            "ell": "4*log(4)",
            "exp_minus_ell_B": "1/4",
            "field_premise": "the source graph is a finite CAR covariance and the two labeled faces use the declared oriented particle-hole Choi pairing",
        },
        "exact_results": {
            "A": "16/17",
            "Cross": "4/17",
            "right_face_occupation": "1/17",
            "occupation_logit_ratio": "16",
            "occupation_logit": "log(16)=2*ell*B",
            "transfer_spectrum": ["1", "1/4", "1/4", "1/16"],
            "one_field_response": "1/4",
            "two_field_response": "1/16",
            "two_mode_response_R": "diag(1/4,1/4)",
            "source_choi_ray_Q": ["1", "1/4", "1/4", "1/16"],
            "source_thermal_density_diag": ["256/289", "16/289", "16/289", "1/289"],
            "source_boundary_face_covariance": {
                "input_dual_boundary_CAR": "(16/17) I_2",
                "output_boundary_CAR": "(1/17) I_2",
            },
            "source_semigroup_two_steps": "Lambda(R)^2=Lambda(R^2)=diag(1,1/16,1/16,1/256)",
        },
        "time_type_audit": {
            "restricted_one_face_heat": {
                "generator_spectrum": ["0", "1/4"],
                "transfer_spectrum": ["1", "1/4"],
                "role": "source-derived one-face thermal dynamics under the declared CAR/KMS premise",
            },
            "thermal_GNS_Liouville": {
                "generator_spectrum": ["0", "0", "-1/4", "1/4"],
                "transfer_spectrum": ["1", "1", "1/4", "4"],
                "fixes_thermal_cyclic_vector": True,
                "is_positive_contraction": False,
                "conclusion": "the effect search constraint T<=I is not a universal KMS/GNS source class",
            },
            "optional_pure_graph_parent": {
                "generator_spectrum": ["0", "1/4", "1/4", "1/2"],
                "transfer_spectrum": ["1", "1/4", "1/4", "1/16"],
                "fixes_graph_cyclic_vector": True,
                "is_positive_contraction": True,
                "conclusion": "valid positive control, but an extra full parent/time choice",
            },
            "source_graph_choi_ray": {
                "transfer_spectrum": ["1", "1/4", "1/4", "1/16"],
                "is_positive_contraction": True,
                "role": "operator ray reconstructed from the four-mode boundary CAR projection and oriented face pairing; no full-time Hamiltonian is inferred",
                "operator_scalar_recovered": False,
            },
        },
        "exact_matrices": {
            "source_response_R": _matrix_json(source_control.response),
            "source_boundary_graph_G4": _matrix_json(source_control.boundary_graph),
            "source_boundary_state_density_Fock16": _matrix_json(source_control.boundary_state_density),
            "source_choi_ray_Q_Fock4": _matrix_json(source_control.transfer),
            "source_thermal_density_QstarQ": _matrix_json(source_control.thermal_density),
            "optional_pure_graph_parent": {
                "source_derived": False,
                "choice": "h=B(I-2G), filled graph orbital and sea subtraction",
                "graph_G2": _matrix_json(control.graph),
                "one_particle_h": _matrix_json(control.one_particle_h),
                "normal_ordered_H": _matrix_json(control.many_body_h),
                "transfer_Q": _matrix_json(control.transfer),
                "state_density_generated_from_G": _matrix_json(control.state_density),
            },
            "one_face_thermal_H": _matrix_json(one_face_h),
            "thermal_GNS_Liouville": _matrix_json(liouville),
            "thermal_GNS_transfer": _matrix_json(gns_transfer),
        },
        "search_fixed_moments": source_fixed_moments(),
        "source_graph_route_receipt": source_control.route_receipt,
        "adversarial_control": {
            "agrees_on_quasiparticle_grades": [0, 1],
            "free_grade2_response": "1/16",
            "interacting_grade2_response": "1/8",
            "full_push_through_covariance": False,
        },
        "quarter_scope": {
            "physical_alpha": "1/4 is a declared boundary-domain label",
            "internal_CAR_quarter": "an independent gauge rotation used by the moment search",
            "identified": False,
        },
        "source_references": SOURCE_REFERENCES,
        "source_sha256": source_hashes,
        "proof_boundary": {
            "source_derived": "the two-face graph and restricted occupation 1/17 under the declared product-collar premise",
            "conditional": "declaring that graph to be the physical finite CAR covariance forces its unique Slater boundary ray; the separately declared oriented face pairing unvectorizes the ray to Lambda(R)",
            "not_source_derived": "the operator's overall partition-function scalar, the pure graph parent h=B(I-2G), a physical full-time Hamiltonian, and the CAR/Choi interpretation from bare P1",
            "not_claimed": ["P1/P2 selection", "physical event law", "identification of physical alpha with internal CAR gauge quarter"],
        },
    }


def run_moment_search(transfer_floor: float, *, solver: str = "SCS") -> dict[str, Any]:
    """Run the real SDP API with only source observables, never matrices/H/rho."""
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
        fixed_moments=source_fixed_moments(),
        source_observables={
            "label": "exact two-field graph-state/Choi-ray source control",
            "source": SOURCE_REFERENCES[0],
            "scope": "vacuum covariance and responses derived from the four-mode boundary CAR projection plus the declared oriented face pairing; no H, rho, graph, or target T matrix supplied; operator scalar and bare-P1 origin remain open",
        },
    )
    return search.solve(
        objective="feasibility", seed=0, solver=solver, tolerance=1e-7
    )


def run_thermal_moment_search(transfer_ceiling: float, *, solver: str = "SCS") -> dict[str, Any]:
    """Check the separately typed thermal-GNS control against its true bound."""
    search_root = str(_search_root())
    if search_root not in sys.path:
        sys.path.insert(0, search_root)
    kernel_search = importlib.import_module("kernel_search")
    search = kernel_search.MomentSearch(
        m=1,
        degree=2,
        transfer_floor=0,
        transfer_ceiling=transfer_ceiling,
        transfer_kind="thermal_gns",
        fixed_moments=thermal_fixed_moments(),
        source_observables={
            "label": "exact one-mode thermal GNS Liouville control",
            "source": SOURCE_REFERENCES[0],
            "scope": "T=exp(-ell(H tensor 1 - 1 tensor H^T)); ceiling4 is the exact finite control bound, while ceiling1 is deliberately the wrong effect class",
        },
    )
    return search.solve(
        objective="feasibility", seed=0, solver=solver, tolerance=1e-7
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--oracle-moments", action="store_true")
    parser.add_argument("--search-controls", action="store_true")
    parser.add_argument("--solver", choices=("CLARABEL", "SCS"), default="SCS")
    args = parser.parse_args()
    report = verify_control()
    if args.oracle_moments:
        report["oracle_moments_through_degree6"] = exact_source_moments(6)
    if args.search_controls:
        accepted = run_moment_search(0.0, solver=args.solver)
        rejected = run_moment_search(0.25, solver=args.solver)
        report["moment_search_controls"] = {
            "floor0_status": accepted.get("solver_status"),
            "floor_quarter_status": rejected.get("solver_status"),
            "floor0_ok": accepted.get("solver_status") in ("optimal", "optimal_inaccurate"),
            "floor_quarter_rejected": rejected.get("solver_status") in ("infeasible", "infeasible_inaccurate"),
            "fixed_moment_count": len(source_fixed_moments()),
            "no_target_matrix_H_or_rho_supplied": True,
            "scope": "source graph-state/Choi-ray effect control under declared CAR/orientation premises; not the thermal GNS transfer",
        }
        if not report["moment_search_controls"]["floor0_ok"]:
            raise ValueError("floor-0 positive control was not accepted")
        if not report["moment_search_controls"]["floor_quarter_rejected"]:
            raise ValueError("global floor-1/4 negative control was not rejected")
        thermal_accepted = run_thermal_moment_search(4, solver=args.solver)
        thermal_wrong_class = run_thermal_moment_search(1, solver=args.solver)
        report["thermal_GNS_search_controls"] = {
            "ceiling4_status": thermal_accepted.get("solver_status"),
            "ceiling1_status": thermal_wrong_class.get("solver_status"),
            "ceiling4_ok": thermal_accepted.get("solver_status") in ("optimal", "optimal_inaccurate"),
            "ceiling1_rejected": thermal_wrong_class.get("solver_status") in ("infeasible", "infeasible_inaccurate"),
            "fixed_moment_count": len(thermal_fixed_moments()),
            "interpretation": "ceiling1 rejection is a declared transfer-class mismatch, not a physical source exclusion",
        }
        if not report["thermal_GNS_search_controls"]["ceiling4_ok"]:
            raise ValueError("thermal GNS control was not accepted at its exact ceiling 4")
        if not report["thermal_GNS_search_controls"]["ceiling1_rejected"]:
            raise ValueError("thermal GNS control was not rejected by the wrong ceiling-1 class")
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(rendered)
    print(json.dumps({
        "ok": report["ok"],
        "verdict": report["verdict"],
        "one_field_response": report["exact_results"]["one_field_response"],
        "two_field_response": report["exact_results"]["two_field_response"],
        "search_controls": report.get("moment_search_controls"),
        "thermal_GNS_search_controls": report.get("thermal_GNS_search_controls"),
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
