#!/usr/bin/env python3
"""Exact add-on check for clock-pair charge projections and packet norms.

This scratch checker does not register or modify a TFPT contract.  If --repo is
supplied, it also compares its pinned exact Weyl matrix with the actual
v16-operator-pullback fold certificate in that repository.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from collections import Counter
from fractions import Fraction as F
from pathlib import Path

import sympy as sp


HERE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser()
parser.add_argument("--repo", type=Path, help="TFPT repository root for source-certificate comparison")
parser.add_argument("--out", type=Path, default=HERE / "projection-results.json")
args = parser.parse_args()

checks: dict[str, bool] = {}


def require(name: str, condition) -> None:
    if not bool(condition):
        raise AssertionError(name)
    checks[name] = True


def dot(x, y):
    return sum(a * b for a, b in zip(x, y))


def norm2(x):
    return dot(x, x)


def matrix_from_fraction_rows(rows):
    return sp.Matrix([[sp.Rational(x.numerator, x.denominator) if isinstance(x, F) else x for x in row] for row in rows])


def perm(order):
    return sp.Matrix([[int(j == order[i]) for j in range(8)] for i in range(8)])


# Pinned exact matrix from the previously reviewed fold checker.
hadamard4 = sp.Matrix(
    [[1, 1, 1, 1], [1, 1, -1, -1], [1, -1, 1, -1], [1, -1, -1, 1]]
) / 2
perm_p = perm([0, 4, 2, 6, 1, 5, 3, 7])
perm_q = perm([2, 3, 4, 6, 7, 0, 1, 5])
weyl_m = perm_q * sp.diag(hadamard4, hadamard4) * perm_p
require("matrix_orthogonal", weyl_m.T * weyl_m == sp.eye(8))
require("matrix_orientation_positive", weyl_m.det() == 1)

matrix_source = {"compared": False}
if args.repo is not None:
    source_path = args.repo.resolve() / "experiments/theory-contracts/v16-operator-pullback-20261014/fold-results.json"
    require("source_fold_results_exists", source_path.is_file())
    source_payload = json.loads(source_path.read_text())
    source_matrix = sp.Matrix([[sp.Rational(str(x)) for x in row] for row in source_payload["explicit_Weyl_matrix"]])
    require("pinned_matrix_matches_actual_fold_certificate", source_matrix == weyl_m)
    matrix_source = {
        "compared": True,
        "path": source_path.relative_to(args.repo.resolve()).as_posix(),
        "sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
    }


def roots(chirality: str):
    result = []
    for i, j in itertools.combinations(range(8), 2):
        for a, b in itertools.product((-1, 1), repeat=2):
            vector = [F(0)] * 8
            vector[i], vector[j] = F(a), F(b)
            result.append(tuple(vector))
    wanted_minus_parity = 0 if chirality == "s" else 1
    result.extend(
        tuple(F(sign, 2) for sign in signs)
        for signs in itertools.product((-1, 1), repeat=8)
        if signs.count(-1) % 2 == wanted_minus_parity
    )
    return result


roots_by_chirality = {chirality: roots(chirality) for chirality in ("s", "c")}
for chirality, system in roots_by_chirality.items():
    require(f"{chirality}_root_count_240", len(system) == len(set(system)) == 240)
    require(f"{chirality}_all_roots_norm_two", all(norm2(root) == 2 for root in system))

s_rootset = set(roots_by_chirality["s"])
transported_s = {tuple(F(x) for x in weyl_m * sp.Matrix(root)) for root in s_rootset}
require("pinned_matrix_preserves_Ls_roots", transported_s == s_rootset)
q = tuple(F(1) for _ in range(8))
require("pinned_matrix_transports_q_as_recorded", weyl_m * sp.ones(8, 1) == sp.Matrix([0, 0, 2, 0, 0, 2, 0, 0]))
primitive_q = tuple(value/2 for value in q)
require("source_q_half_is_s_root_not_c_root", primitive_q in s_rootset and primitive_q not in set(roots_by_chirality["c"]))
mixed_weight_roots = {root for root in s_rootset if norm2(root[:5]) == 1 and norm2(root[5:]) == 1}
mixed_cartesian_weights = set()
for carrier, family in itertools.product(range(5), range(5, 8)):
    for carrier_sign, family_sign in itertools.product((-1, 1), repeat=2):
        weight = [F(0)]*8
        weight[carrier], weight[family] = F(carrier_sign), F(family_sign)
        mixed_cartesian_weights.add(tuple(weight))
require("mixed_q_half_root_class_is_exact_10_tensor_6", mixed_weight_roots == mixed_cartesian_weights and len(mixed_weight_roots) == 60)
primitive_q_image = tuple(F(value) for value in weyl_m*sp.Matrix(primitive_q))
require("actual_source_q_half_is_mixed_10_6_root", primitive_q_image in mixed_weight_roots)
require("all_primitive_carrier_directions_are_unit_integer_axes", all(sum(value != 0 for value in root[:5]) == 1 and all(value.denominator == 1 for value in root) for root in mixed_weight_roots))
require("primitive_root_source_C_grade2", dot(tuple(F(k, 4) for k in range(8)), primitive_q) % 1 == F(1, 2))
require("primitive_root_source_Cartan_charge_four", dot(q, primitive_q) == 4)
primitive_energies = {}
for alpha, expected in ((F(1, 4), F(7, 8)), (F(3, 4), F(9, 8))):
    delta = tuple((F(k)+alpha)/8-F(1, 2) for k in range(8))
    energy = norm2(primitive_q)/2+dot(delta, primitive_q)
    require(f"primitive_root_quarter_source_energy_{alpha}", energy == expected)
    primitive_energies[str(alpha)] = {"fold_K": str(energy), "source_H_8K": str(8*energy)}


def fixed_components(system, clock):
    remaining = {root for root in system if dot(clock, root).denominator == 1}
    components = []
    while remaining:
        seed = remaining.pop()
        component = {seed}
        frontier = [seed]
        while frontier:
            root = frontier.pop()
            neighbors = {other for other in remaining if dot(root, other) != 0}
            remaining -= neighbors
            component |= neighbors
            frontier.extend(neighbors)
        components.append(component)
    return sorted(components, key=len)


def projection_norm(vector, roots_in_component):
    root_matrix = sp.Matrix.hstack(*(sp.Matrix(root) for root in roots_in_component))
    basis = root_matrix.columnspace()
    basis_matrix = sp.Matrix.hstack(*basis)
    gram = basis_matrix.T * basis_matrix
    projector = basis_matrix * gram.inv() * basis_matrix.T
    v = sp.Matrix(vector)
    return sp.simplify((v.T * projector * v)[0]), basis_matrix.rank()


projection_rows = {}
expected_pairs = {"s": (sp.Integer(4), sp.Integer(4)), "c": (sp.Integer(6), sp.Integer(2))}
for chirality, system in roots_by_chirality.items():
    rows = []
    for k in range(8):
        clock = tuple(F(k + j, 4) for j in range(8))
        fixed = [root for root in system if dot(clock, root).denominator == 1]
        components = fixed_components(system, clock)
        component_by_size = {len(component): component for component in components}
        require(f"{chirality}_k{k}_fixed_root_count", len(fixed) == 52)
        require(f"{chirality}_k{k}_actual_component_sizes", sorted(component_by_size) == [12, 40])

        d5_norm, d5_rank = projection_norm(q, component_by_size[40])
        a3_norm, a3_rank = projection_norm(q, component_by_size[12])
        require(f"{chirality}_k{k}_actual_component_ranks", (d5_rank, a3_rank) == (5, 3))
        require(f"{chirality}_k{k}_projection_norm_pair", (d5_norm, a3_norm) == expected_pairs[chirality])
        rows.append(
            {
                "k": k,
                "clock": [str(x) for x in clock],
                "fixed_root_count": len(fixed),
                "components": [
                    {"root_count": 40, "rank": d5_rank, "q_projection_norm2": str(d5_norm)},
                    {"root_count": 12, "rank": a3_rank, "q_projection_norm2": str(a3_norm)},
                ],
            }
        )
    require(
        f"{chirality}_projection_pair_invariant_over_all_offsets",
        all(
            (row["components"][0]["q_projection_norm2"], row["components"][1]["q_projection_norm2"])
            == tuple(map(str, expected_pairs[chirality]))
            for row in rows
        ),
    )
    projection_rows[chirality] = rows


# Exact A3-to-D3 weight-coordinate isometry on the zero-sum four-plane.
t_map = sp.Matrix([[1, 1, -1, -1], [1, -1, 1, -1], [1, -1, -1, 1]]) / 2
zero_sum_projector = sp.eye(4) - sp.ones(4, 4) / 4
require("T_transpose_T_is_zero_sum_projector", t_map.T * t_map == zero_sum_projector)
require("T_Ttranspose_is_identity", t_map * t_map.T == sp.eye(3))

m_a = sp.Matrix([0, 0, -1, -1])
beta_a = m_a + sp.ones(4, 1) / 2
require("packet_example_sum_beta_zero", sum(beta_a) == 0)
require("packet_example_beta_norm_one", (beta_a.T * beta_a)[0] == 1)
require("packet_example_Tbeta", t_map * beta_a == sp.Matrix([1, 0, 0]))


def tuples_with_sum_and_norm2(target_norm2):
    # For target norm <= 8, |m_j| <= sqrt(8) < 3; [-2,2]^4 is exhaustive.
    return [
        values
        for values in itertools.product(range(-2, 3), repeat=4)
        if sum(values) == -2 and sum(x * x for x in values) == target_norm2
    ]


norm2_class = tuples_with_sum_and_norm2(2)
norm8_class = tuples_with_sum_and_norm2(8)
require(
    "packet_norm2_class_is_two_minus_ones",
    len(norm2_class) == 6 and all(Counter(values) == Counter({-1: 2, 0: 2}) for values in norm2_class),
)
require("packet_norm8_class_empty", norm8_class == [])
minimal_family_weights = {tuple(F(value) for value in t_map*sp.Matrix([F(entry)+F(1, 2) for entry in values])) for values in norm2_class}
vector6_weights = {tuple(F(sign) if axis == j else F(0) for axis in range(3)) for j in range(3) for sign in (-1, 1)}
require("six_minimal_packet_marks_are_actual_family6_weights", minimal_family_weights == vector6_weights)
require(
    "packet_beta_norm_identity_on_exhaustive_box",
    all(
        sum((F(x) + F(1, 2)) ** 2 for x in values) == sum(x * x for x in values) - 1
        for values in itertools.product(range(-2, 3), repeat=4)
        if sum(values) == -2
    ),
)
require(
    "packet_sumsq_even_from_sum_constraint",
    all(sum(x * x for x in values) % 2 == 0 for values in itertools.product(range(-2, 3), repeat=4) if sum(values) == -2),
)


def packet_solutions(target_family_norm2):
    found = []
    # Nonzero integer t obeys t^2 <= target when sum(m^2)-1 >= 1.
    for t in range(-target_family_norm2, target_family_norm2 + 1):
        if t == 0:
            continue
        for values in itertools.product(range(-2, 3), repeat=4):
            if sum(values) != -2:
                continue
            if t * t * (sum(x * x for x in values) - 1) == target_family_norm2:
                found.append((t, values))
    return found


family4_solutions = packet_solutions(4)
family2_solutions = packet_solutions(2)
require(
    "family_norm4_selects_t_plusminus2_and_minimal_m",
    {t for t, _ in family4_solutions} == {-2, 2}
    and {values for _, values in family4_solutions} == set(norm2_class),
)
require("family_norm2_has_no_integer_packet_solution", family2_solutions == [])
require("family_norm2_only_possible_integer_t_control", [t for t in range(-2, 3) if t and 2 % (t * t) == 0] == [-1, 1])

# Trace one does not select the primitive norm-one carrier direction.
x_weak = (1, 0, 0, 0, 0)
x_trace_counterexample = (2, -1, 0, 0, 0)
require("carrier_trace_one_both_examples", sum(x_weak) == sum(x_trace_counterexample) == 1)
require("carrier_norm_hidden_input_witness", norm2(x_weak) == 1 and norm2(x_trace_counterexample) == 5)

# Compare the actual, separately obtained flavor local system only under an
# EXTRA parallel identification of its triplet+singlet with the SU4 lift.
# Its generator B_F is not the eight-dimensional memory root p*.
family_m = sp.Matrix([
    [0, -(1 + sp.I)/2, (1 - sp.I)/2],
    [-(1 + sp.I)/2, -sp.I/2, -sp.Rational(1, 2)],
    [(1 - sp.I)/2, -sp.Rational(1, 2), sp.I/2],
])
family_u = sp.diag(1, sp.I, -sp.I)
family_loops = [sp.simplify(family_u**j * family_m * family_u**(-j)) for j in range(4)]
family_source = {"compared": False}
if args.repo is not None:
    family_path = args.repo.resolve() / "experiments/theory-contracts/source-three-route-closure-20260922/family_certificate.json"
    family_cert = json.loads(family_path.read_text())
    require("flavor_M_matches_actual_certificate", str(family_m) == family_cert["matrices"]["M"])
    require("flavor_deck_matches_actual_certificate", str(family_u) == family_cert["matrices"]["U_deck"])
    v117_path = args.repo.resolve() / "verification/v117_monodromy_weyl_a3.py"
    v117_sha = hashlib.sha256(v117_path.read_bytes()).hexdigest()
    require("actual_flavor_original_v117_pin_current", family_cert["source_hashes"][str(v117_path)] == v117_sha)
    family_source = {"compared": True, "path": family_path.relative_to(args.repo.resolve()).as_posix(),
                     "sha256": hashlib.sha256(family_path.read_bytes()).hexdigest(), "v117_sha256": v117_sha}
require("actual_flavor_loop_product_identity", sp.simplify(sp.prod(family_loops)) == sp.eye(3))
for size, matrices, expected in ((3, family_loops, 1), (4, [sp.diag(h, 1) for h in family_loops], 2)):
    entries = sp.symbols(f"z0:{size*size}")
    unknown = sp.Matrix(size, size, entries)
    equations = [entry for h in matrices for entry in unknown*h-h*unknown]
    coeff = sp.linear_eq_to_matrix(equations, entries)[0]
    require(f"actual_flavor_commutant_dimension_{size}", size*size - coeff.rank() == expected)
family_bridge_rows = []
for values in norm2_class:
    b_f = sp.diag(*[sp.Rational(value) + sp.Rational(1, 2) for value in values])
    comm = b_f*sp.diag(family_m, 1)-sp.diag(family_m, 1)*b_f
    comm_norm = sp.simplify(sp.trace(comm.H*comm))
    require(f"minimal_lift_not_parallel_actual_flavor_{values}", comm_norm > 0)
    family_bridge_rows.append({"m": list(values), "first_loop_commutator_norm2": str(comm_norm)})
for sign in (1, -1):
    b_f = sp.diag(sign*sp.Rational(1, 2), sign*sp.Rational(1, 2),
                 sign*sp.Rational(1, 2), -sign*sp.Rational(3, 2))
    require(f"parallel_triplet_scalar_control_{sign}", all(b_f*sp.diag(h, 1) == sp.diag(h, 1)*b_f for h in family_loops))
    require(f"parallel_triplet_control_norm3_{sign}", sp.trace(b_f*b_f) == 3)
r_sym, t_sym = sp.symbols("r t", integer=True)
b_sym = r_sym + sp.Rational(1, 2)
require("parallel_family_norm_symbolic_identity", sp.expand(t_sym**2*12*b_sym**2 - 3*t_sym**2*(2*r_sym+1)**2) == 0)
require("parallel_family_norm_targets_not_divisible_by3", 4 % 3 != 0 and 2 % 3 != 0)
memory_root = sp.Matrix([sp.Rational(1, 2)]*4 + [-sp.Rational(1, 2)]*4)
memory_image = weyl_m*memory_root
require("memory_root_is_distinct_D5_object", tuple(memory_image) == (1, 0, 0, 1, 0, 0, 0, 0))

# A constructive alternative uses the four-mark A3 Cartan as a residue
# VECTOR space, not a parallel endomorphism of a fundamental 3+1 module.
native_n = sp.Matrix([[1, 0, 0], [0, 5, -2], [0, -2, 1]])
residue_h = sp.Matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1], [-1, -1, -1]])
residue_x = sp.Matrix([[-1, -3, 1], [1, 1, -1], [1, -1, 1]])
b_native = (1 + sp.I)/2
native_s = sp.Matrix([[0, -2, 1], [-b_native, sp.I*b_native, 0], [-sp.I*b_native, b_native, 0]])
residue_l = residue_h*residue_x/2
residue_source = {"compared": False}
if args.repo is not None:
    native_path = args.repo.resolve() / "experiments/theory-contracts/source-native-clock-flavor-bridge-20260922/certificate.json"
    native_cert = json.loads(native_path.read_text())
    for key, matrix in {"N": native_n, "H": residue_h, "X": residue_x, "S": native_s, "M": family_m, "U": family_u}.items():
        require(f"residue_dictionary_{key}_matches_original_certificate", str(matrix) == native_cert["matrices"][key])
    residue_source = {"compared": True, "path": native_path.relative_to(args.repo.resolve()).as_posix(),
                      "sha256": hashlib.sha256(native_path.read_bytes()).hexdigest()}
require("residue_L_is_N_isometry", residue_l.H*residue_l == native_n)
require("residue_S_is_N_isometry", sp.simplify(native_s.H*native_s) == native_n)
require("residue_L_in_traceless_A3_hyperplane", sp.ones(1, 4)*residue_l == sp.zeros(1, 3))
axes = [sp.Matrix([1, 0, 0]), sp.Matrix([0, 0, 1]), sp.Matrix([0, 1, 2])]
axis_matrix = sp.Matrix.hstack(*axes)
h_axes = [sp.simplify(native_s*axis) for axis in axes]
require("three_residue_axes_N_orthonormal", axis_matrix.T*native_n*axis_matrix == sp.eye(3))
require("three_monodromy_axes_orthonormal", sp.simplify(sp.Matrix.hstack(*h_axes).H*sp.Matrix.hstack(*h_axes)) == sp.eye(3))
residue_rows = []
for values in norm2_class:
    beta_f = sp.Matrix([sp.Rational(value) + sp.Rational(1, 2) for value in values])
    axis = native_n.inv()*residue_l.T*beta_f
    matches = [(index, sign) for index, direction in enumerate(axes) for sign in (1, -1) if axis == sign*direction]
    require(f"minimal_lift_residue_axis_{values}", residue_l*axis == beta_f and len(matches) == 1)
    residue_rows.append({"m": list(values), "native_axis": list(map(str, axis)), "axis_index_sign": list(matches[0])})
for index, target in enumerate((2, 0, 1)):
    require(f"monodromy_cycles_minimal_residue_axis_{index}", sp.simplify(family_m*h_axes[index] - h_axes[target]) == sp.zeros(3, 1))
for index, (target, sign) in enumerate(((2, -1), (1, 1), (0, 1))):
    require(f"deck_transports_minimal_residue_axis_{index}", sp.simplify(family_u*h_axes[index] - sign*h_axes[target]) == sp.zeros(3, 1))
require("two_minimal_axes_are_not_U_eigenlines", all(sp.Matrix.hstack(h_axes[index], family_u*h_axes[index]).rank() == 2 for index in (0, 2)))
native_b = sp.Matrix([[1, 0, 0], [0, 0, -1], [0, -1, 0]])
for index, direction in enumerate(h_axes):
    require(f"minimal_residue_axis_J_real_{index}", sp.simplify(native_b*sp.conjugate(direction)-direction) == sp.zeros(3, 1))

def alternating_source(direction):
    h1, h2, h3 = direction
    return sp.Matrix([[0, h3, -h2], [-h3, 0, h1], [h2, -h1, 0]])

generic_c = sp.Matrix(3, 3, sp.symbols("c0:9"))
seed_h = h_axes[1]  # e1, already exhaustively tested in the old one-clock class.
seed_y = generic_c.T*alternating_source(seed_h)-alternating_source(seed_h)*generic_c
for power in (1, 2):
    rotation = family_m**power
    rotated_c = rotation*generic_c*rotation.inv()
    rotated_h = rotation*seed_h
    rotated_y = rotated_c.T*alternating_source(rotated_h)-alternating_source(rotated_h)*rotated_c
    difference = rotated_y-rotation.inv().T*seed_y*rotation.inv()
    require(f"minimal_mark_vertex_congruence_symbolic_{power}", all(sp.simplify(entry) == 0 for entry in difference))
if args.repo is not None:
    seed_scan = native_cert["three_marked_line_scan"]["e1"]
    require("old_seed_full_rank_spectrum_is_twofold_degenerate", seed_scan["full_rank_characteristic_histogram"] == {"(x - 4)*(x - 1)**2": 8})
    require("old_seed_shared_axis_CP_pairs_vanish", seed_scan["nonzero_full_rank_cp_pairs"] == 0)


output = {
    "research_id": "UR.SEAM.JOINT_WARD_STATE.01",
    "addon": "clock-pair-projection-and-packet-norm",
    "verdict": "SCOPED_EXACT_CHECKS_PASS",
    "checks": checks,
    "n_checks": len(checks),
    "matrix_source": matrix_source,
    "clock_projection_results": projection_rows,
    "projection_norm_pairs_D5_A3": {"Ls": [4, 4], "Lc": [6, 2]},
    "T_matrix": [[str(x) for x in row] for row in t_map.tolist()],
    "T_transpose_T": [[str(x) for x in row] for row in (t_map.T * t_map).tolist()],
    "packet_example": {
        "m": list(map(int, m_a)),
        "beta": [str(x) for x in beta_a],
        "beta_norm2": str((beta_a.T * beta_a)[0]),
        "T_beta": [str(x) for x in t_map * beta_a],
    },
    "packet_classification": {
        "sum_m_minus2_norm2_solutions": [list(values) for values in norm2_class],
        "sum_m_minus2_norm8_solutions": [list(values) for values in norm8_class],
        "family_norm4_integer_t": sorted({t for t, _ in family4_solutions}),
        "family_norm2_integer_solutions": [],
    },
    "carrier_trace_counterexample": {
        "primitive_assumed": {"x": list(x_weak), "trace": 1, "norm2": 1},
        "same_trace": {"x": list(x_trace_counterexample), "trace": 1, "norm2": 5},
        "conclusion": "Trace one alone does not derive the primitive norm-one carrier direction.",
    },
    "parallel_flavor_comparison": {
        "extra_premise": "Identify the SU4 lift with the separately given irreducible SU3_F triplet plus singlet, and require the determinant generator to preserve that connection parallelly.",
        "actual_flavor_source": family_source,
        "commutant_dimensions": {"triplet": 1, "triplet_plus_singlet": 2},
        "minimal_lift_witnesses": family_bridge_rows,
        "all_integer_lifts_norm_identity": "3*t^2*(2*r+1)^2",
        "required_source_family_norm2": [4, 2],
        "scope": "Excluded only under the EXTRA parallel identification; not a general TFPT or independent family-connection exclusion. B_F is distinct from the D5 memory root p*.",
    },
    "constructive_cartan_residue_bridge": {
        "original_dictionary_source": residue_source,
        "L": [[str(entry) for entry in row] for row in residue_l.tolist()],
        "minimal_marks_to_three_axes": residue_rows,
        "monodromy_cycle_axis_indices": [2, 0, 1],
        "deck_axis_index_sign": [[2, -1], [1, 1], [0, 1]],
        "scope": "Exact four-mark Cartan residue-vector/isometry bridge using the prior finite metric N and orientation-twisted M/U dictionary. It is not a parallel End(F3) identification. Native monodromy U is not automatically the folded source g, C or physical time.",
        "not_selected": "No choice of axis, physical Higgs field, coupling, actual collar transport, raw E8 extension, source state or event law is derived.",
        "shared_axis_one_clock_corollary": "All minimal directions are signed M-orbit images of e1. Unitary congruence transfers the old complete e1 scan: every full-rank single-C vertex has squared spectrum 1,1,4 up to scale, and a common minimal axis with real clock coefficients gives no CP-odd invariant. This does not exclude multi-vertex or source-derived complex interactions.",
    },
    "intrinsic_primitive_charge_channel": {
        "source_q_half_is_root": True,
        "mixed_root_class": "(10,6)",
        "weight_count": 60,
        "canonical_q_half": list(map(str, primitive_q_image)),
        "primitive_carrier_family_norm2": [1, 1],
        "derived_q_relation": "q=2*lambda, lambda a mixed (10,6) E8 root",
        "coarse_C_grade": 2,
        "root_current_conformal_weight": 1,
        "root_current_source_Cartan_charge": 4,
        "root_state_actual_quarter_time_energies": primitive_energies,
        "scope": "Derived from the retained Ls source q/2 lattice class plus the inherited joint clock projections, not from projection norms alone. This is a charge/weight dictionary, not a physical Chern-clutch, CAR field/statistics, native interaction, coupling, state or time derivation.",
    },
    "scope": (
        "Exact simultaneous inherited g^2 clock and source-charge projections, with the actual D8 currents, positive even-unimodular lattice extensions and determinant-only packet formula retained. "
        "The stronger packet selector uses intrinsic integral carrier trace, not an assumed unit-clutch direction or q=k lambda_m. The dressed clock and raw physical dictionary are still premises. "
        "The flavor comparison adds an explicitly separate parallel-identification premise."
    ),
    "not_claimed": [
        "No new source, E8 field map, physical time, state, or instrument selection.",
        "No universal exclusion beyond the declared D8-containing fold and determinant-only packet comparison.",
        "No contract metadata, Atlas, ledger, paper, or verification update by this checker.",
    ],
}

args.out.write_text(json.dumps(output, indent=2) + "\n")
print(f"{len(checks)}/{len(checks)} checks passed")
print(hashlib.sha256(args.out.read_bytes()).hexdigest())
