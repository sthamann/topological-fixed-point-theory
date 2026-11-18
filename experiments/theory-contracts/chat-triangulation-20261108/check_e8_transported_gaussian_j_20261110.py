"""Exact transport of the source-defined Gaussian J under native E8 maps."""
import ast
from collections import Counter, deque
from itertools import combinations, permutations, product
import json
from pathlib import Path
import sys

import sympy as sp
from sympy.combinatorics import Permutation, PermutationGroup

repo = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(repo / "verification"))
from v55_coxeter_cycle import _coxeter, _e8_cartan
from v973_seam_route_narrowing import _e8_simple_roots

source = repo / "experiments/tfpt-discovery/clock_fiber_product_probe.py"
functions = [
    n for n in ast.parse(source.read_text()).body
    if isinstance(n, ast.FunctionDef) and n.name in ("J_all", "sigma")
]
native = {}
exec(compile(ast.Module(body=functions, type_ignores=[]), str(source), "exec"), native)
identity = sp.eye(8)
basis = sp.Matrix.hstack(*map(sp.Matrix, _e8_simple_roots()))
assert basis.T * basis == _e8_cartan()
coxeter = basis * _coxeter(_e8_cartan()) * basis.inv()
gaussian, family = (
    sp.Matrix.hstack(*(sp.Matrix(native[name](tuple(identity[:, i]))) for i in range(8)))
    for name in ("J_all", "sigma")
)
assert gaussian**2 == -identity and coxeter**15 == -identity and family**3 == identity
assert gaussian * family == family * gaussian
seam = gaussian * family
assert seam**4 == family and seam**9 == gaussian

# Two times each root, with coordinates integral in the native frame.
roots = set()
for a, b in combinations(range(8), 2):
    for x, y in product((-2, 2), repeat=2):
        root = [0] * 8
        root[a], root[b] = x, y
        roots.add(tuple(root))
roots.update(signs for signs in product((-1, 1), repeat=8) if signs.count(-1) % 2 == 0)
assert len(roots) == 240
roots = sorted(roots)
index = {root: i for i, root in enumerate(roots)}


def root_permutation(matrix):
    return bytes(index[tuple(matrix * sp.Matrix(root))] for root in roots)


def compose(left, right):
    return bytes(left[i] for i in right)


def inverse(permutation):
    result = bytearray(len(permutation))
    for i, j in enumerate(permutation):
        result[j] = i
    return bytes(result)


def conjugacy_orbit(origin, generators):
    actions = [(g, inverse(g)) for g in generators]
    distance = {origin: 0}
    pending = deque([origin])
    while pending:
        current = pending.popleft()
        for forward, back in actions:
            result = compose(forward, compose(current, back))
            if result not in distance:
                distance[result] = distance[current] + 1
                pending.append(result)
    return distance


def order(*generators):
    return PermutationGroup(*(Permutation(list(g)) for g in generators)).order()


def commutator(left, right):
    return left * right * left.T * right.T


def structure_readout(transport):
    transported = transport * gaussian * transport.T
    # First 5 coordinates carry D5 and last 3 carry D3=A3 in v634.
    mixing = transported[:5, 5:]
    return {
        "equals_J": transported == gaussian,
        "equals_minus_J": transported == -gaussian,
        "difference_rank": (transported - gaussian).rank(),
        "family_commutator_rank": (transported * family - family * transported).rank(),
        "D5_A3_mixing_rank": mixing.rank(),
        "D5_A3_mixing_squared_frobenius": str(sp.trace(mixing.T * mixing)),
        "first_glue_block_invariant": transport[:5, 5:] == sp.zeros(5, 3),
        "transport_order": order(root_permutation(transport)),
    }


C, F, J, S = map(root_permutation, (coxeter, family, gaussian, seam))
results = {}
for label, gens in {"C": (C,), "C_family": (C, F), "C_seam": (C, S)}.items():
    distances = conjugacy_orbit(J, gens)
    group_order = order(*gens)
    assert group_order % len(distances) == 0
    results[label] = {
        "generated_group_order": group_order,
        "orbit_size": len(distances),
        "stabilizer_order": group_order // len(distances),
        "positive_generator_distance_histogram": dict(sorted(Counter(distances.values()).items())),
        "contains_minus_J": inverse(J) in distances,
    }
results["coxeter_returns_J"] = [
    k for k in range(30) if coxeter**k * gaussian == gaussian * coxeter**k
]
results["coxeter_returns_up_to_sign"] = [
    k for k in range(30) if coxeter**k * gaussian == -gaussian * coxeter**k
]
results["loops"] = {
    label: structure_readout(matrix) for label, matrix in {
        "identity": identity,
        "one_C": coxeter,
        "C6": coxeter**6,
        "C15": coxeter**15,
        "family": family,
        "seam": seam,
        "C_family_commutator": commutator(coxeter, family),
        "C_seam_commutator": commutator(coxeter, seam),
    }.items()
}
# In v634, the D5 and D3 spans are the first five and last three coordinates.
# A full isometry of this split must be signed-permutation on each summand:
# the D5 and D3 root automorphisms both have that form. Preservation of the
# half roots then requires an even number of minus signs over all 8 slots.
def split_family_stabilizer():
    family_map = [list(family[:, i]).index(1) for i in range(8)]
    for d5 in permutations(range(5)):
        for d3 in permutations(range(5, 8)):
            axes = d5 + d3
            if any(axes[family_map[i]] != family_map[axes[i]] for i in range(8)):
                continue
            parity = sum(a > b for i, a in enumerate(axes) for b in axes[i + 1 :])
            if parity % 2:
                continue
            for signs in product((-1, 1), repeat=8):
                if signs.count(-1) % 2:
                    continue
                if any(signs[family_map[i]] != signs[i] for i in range(8)):
                    continue
                matrix = sp.Matrix.hstack(*(signs[i] * identity[:, axes[i]] for i in range(8)))
                assert matrix.det() == 1 and matrix * family == family * matrix
                p = root_permutation(matrix)
                assert family_group.contains(Permutation(list(p)))
                yield matrix, signs[:5].count(-1) % 2 == 0


family_group = PermutationGroup(Permutation(list(C)), Permutation(list(F)))
assert family_group.order() == 348364800
assert family_group.derived_subgroup().order() == family_group.order()
results["joint_group_perfect"] = True
marked_stabilizer = list(split_family_stabilizer())
results["return_to_both_fixed_marks"] = {}
for label, strict in [("preserve_split_and_family", False), ("preserve_glue_labels_and_family", True)]:
    group = [g for g, preserves_glue_labels in marked_stabilizer if preserves_glue_labels or not strict]
    compatible = [g for g in group if g * gaussian == gaussian * g]
    orbit = {tuple(g * gaussian * g.T) for g in group}
    assert all(a * b == b * a for a, b in combinations(group, 2))
    assert len(group) // len(compatible) == len(orbit) == 6
    results["return_to_both_fixed_marks"][label] = {
        "stabilizer_of_returned_marks": len(group),
        "also_restores_J": len(compatible),
        "J_fiber_size": len(orbit),
        "effective_action": "cyclic order six",
    }
d5_cycle = identity[:, [2, 1, 4, 3, 0, 5, 6, 7]]
assert (d5_cycle, True) in marked_stabilizer
assert structure_readout(d5_cycle)["difference_rank"] > 0
results["loops"]["D5_cycle_after_both_marks_return"] = structure_readout(d5_cycle)
weyl_generators = tuple(root_permutation(
    identity - basis[:, i] * basis[:, i].T
) for i in range(8))
weyl_order = order(*weyl_generators)
assert weyl_order == 2 * family_group.order()
assert all(m.det() == 1 for m in (coxeter, family, gaussian, seam))
results["full_Weyl_group_order"] = weyl_order

# The same equations have rational coefficients over R, so nullity one
# excludes even an arbitrary real, fixed complex structure for both clocks.
commutation_rows = []
for matrix in (coxeter, family):
    commutation_rows.append(sp.Matrix.hstack(*(
        (unit * matrix - matrix * unit).reshape(64, 1)
        for i in range(8) for j in range(8)
        for unit in [identity[:, i] * identity[j, :]]
    )))
commutation_system = sp.Matrix.vstack(*commutation_rows)
commutation_system *= sp.ilcm(*(v.q for v in commutation_system))
assert all(v.is_Integer for v in commutation_system)
assert commutation_system * identity.reshape(64, 1) == sp.zeros(128, 1)
# A rank-63 minor nonzero modulo 101 is also nonzero over Q and R.
assert commutation_system.to_DM().convert_to(sp.GF(101)).rank() == 63
common_dimension = 1
results["real_joint_commutant_dimension"] = common_dimension

x = sp.Symbol("x")
results["gaussian_defect"] = {}
for label, matrix in [("Coxeter_C", coxeter), ("seam_c", seam), ("family_sigma", family)]:
    defect = commutator(matrix, gaussian)
    p = Permutation(list(root_permutation(defect)))
    results["gaussian_defect"][label] = {
        "order": int(p.order()),
        "characteristic_polynomial": str(sp.factor(defect.charpoly(x).as_expr())),
        "rank_defect_minus_identity": (defect - identity).rank(),
        "root_cycle_census": dict(sorted(Counter(map(len, p.full_cyclic_form)).items())),
    }
coxeter_defect = commutator(coxeter, gaussian)
assert coxeter_defect.charpoly(x).as_expr() == sp.expand(sp.cyclotomic_poly(10, x)**2)
assert coxeter_defect**5 == -identity

omega = coxeter**10
eisenstein_rotation = 2 * omega + identity
assert omega**2 + omega + identity == sp.zeros(8)
assert eisenstein_rotation**2 == -3 * identity
assert eisenstein_rotation.T == -eisenstein_rotation
assert omega.charpoly(x).as_expr() == sp.expand(sp.cyclotomic_poly(3, x)**4)
results["Eisenstein"] = {
    "C10_order": 3,
    "C10_characteristic_polynomial": str(sp.factor(omega.charpoly(x).as_expr())),
    "C10_fixed_dimension": 8 - (omega - identity).rank(),
    "family_fixed_dimension": 8 - (family - identity).rank(),
    "C10_family_commutator_rank": (omega * family - family * omega).rank(),
    "C10_gaussian_commutator_rank": (omega * gaussian - gaussian * omega).rank(),
    "real_Coxeter_complex_structure": "(I+2*C^10)/sqrt(3)",
}
# Exact common-coordinate diagnostic. v775 gives the role charge but
# warns that its placement in the root model is not a physical embedding.
hypercharge = sp.Matrix([-2, -2, -2, 3, 3, 0, 0, 0]) / 6
marked_charge_stabilizer = [g for g, _ in marked_stabilizer if g * hypercharge == hypercharge]
assert len(marked_charge_stabilizer) == 2
assert all(g * gaussian == gaussian * g for g in marked_charge_stabilizer)
project_A3 = sp.diag(0, 0, 0, 0, 0, 1, 1, 1)
charge_coupling = lambda j: (project_A3 * j * hypercharge).norm() ** 2
shifted = d5_cycle * gaussian * d5_cycle.T
assert gaussian * identity[:, 0] == identity[:, 1]
assert shifted * identity[:, 0] == identity[:, 5]
assert charge_coupling(gaussian) == sp.Rational(1, 4)
assert charge_coupling(shifted) == sp.Rational(1, 9)
results["common_coordinate_hypercharge"] = {
    "preserving_split_family_and_Y": len(marked_charge_stabilizer),
    "all_restore_J": True,
    "initial_J_e1": "e2",
    "D5_cycle_transported_J_e1": "e6",
    "initial_squared_A3_component_of_J_Y": str(charge_coupling(gaussian)),
    "D5_cycle_squared_A3_component_of_J_Y": str(charge_coupling(shifted)),
}
# Vary role assignments on the same v634 split to keep the dependency visible.
role_census = Counter()
for color in combinations(range(5), 3):
    y = sp.Matrix([-2 if i in color else 3 if i < 5 else 0 for i in range(8)])
    fiber = {tuple(g * gaussian * g.T) for g, _ in marked_stabilizer if g * y == y}
    role_census[len(fiber)] += 1
assert role_census == {1: 9, 3: 1}
results["same_split_varying_color_slots_J_fiber_size_census"] = dict(sorted(role_census.items()))
def reconstruct_gaussian(split_projector, family_step, y):
    b3 = 2 * (identity - split_projector) * family_step * y
    b = [family_step * b3, family_step**2 * b3, b3]
    ya = y - sum((v * v.dot(y) for v in b), sp.zeros(8, 1))
    a3 = sp.Rational(6, 5) * (ya - 2 * (identity + family_step + family_step**2) * ya)
    a = [family_step * a3, family_step**2 * a3, a3]
    frame = sp.Matrix.hstack(*(v for pair in zip(a, b) for v in pair))
    assert frame.T * frame == sp.eye(6)
    u, v = sp.GramSchmidt(frame.T.nullspace(), orthonormal=True)
    if sp.Matrix.hstack(frame, u, v).det() < 0:
        v = -v
    return sum((w * z.T - z * w.T for z, w in list(zip(a, b)) + [(u, v)]), sp.zeros(8))


for transport in (identity, coxeter, family, seam, d5_cycle):
    P = transport * (identity - project_A3) * transport.T
    F_transported = transport * family * transport.T
    Y_transported = transport * hypercharge
    assert reconstruct_gaussian(P, F_transported, Y_transported) == transport * gaussian * transport.T
results["equivariant_reconstruction_tests"] = ["identity", "C", "sigma", "seam", "D5_cycle"]
print(json.dumps(results, indent=2))
