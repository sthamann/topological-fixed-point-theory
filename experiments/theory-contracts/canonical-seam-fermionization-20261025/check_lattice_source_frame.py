#!/usr/bin/env python3
"""Exact frame test for the canonical D5+A3 halfway lattice and folded CAR D8."""

from __future__ import annotations

from collections import Counter
from fractions import Fraction as Q
from itertools import combinations, product
import hashlib
import json
from pathlib import Path
import sys

import sympy as sp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT))

from tfpt_explorer import cartan_source as cs  # noqa: E402


def need(value: bool, label: str) -> None:
    if not value:
        raise RuntimeError(label)


def mat_vec(matrix, vector):
    return tuple(sum(a * b for a, b in zip(row, vector)) for row in matrix)


def dot(left, right):
    return sum(a * b for a, b in zip(left, right))


def fstr(value: Q) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def vecstr(vector) -> list[str]:
    return [fstr(Q(value)) for value in vector]


def counterstr(counter: Counter) -> dict[str, int]:
    return {fstr(Q(key)): counter[key] for key in sorted(counter)}


def standard_d8_roots() -> set[tuple[Q, ...]]:
    roots = set()
    for i, j in combinations(range(8), 2):
        for a, b in product((-1, 1), repeat=2):
            vector = [Q(0)] * 8
            vector[i], vector[j] = Q(a), Q(b)
            roots.add(tuple(vector))
    need(len(roots) == 112, "standard D8 root count")
    return roots


def e8_roots() -> set[tuple[Q, ...]]:
    integer = standard_d8_roots()
    spinor = {
        tuple(Q(sign, 2) for sign in signs)
        for signs in product((-1, 1), repeat=8)
        if sum(sign < 0 for sign in signs) % 2 == 0
    }
    need(len(spinor) == 128 and len(integer | spinor) == 240, "standard E8 roots")
    return integer | spinor


def phase_sign(h, root) -> int:
    twice = 2 * dot(h, root)
    need(twice.denominator == 1, "defect character is a sign")
    return 1 if int(twice) % 2 == 0 else -1


def calculate() -> dict:
    lattice = cs._lattice_data()
    source_basis = lattice["basis"]
    target_basis = sp.Matrix.hstack(*(sp.Matrix(root) for root in cs._TARGET_SIMPLE))
    O = target_basis * source_basis.inv() / sp.sqrt(2)
    need(sp.simplify(O.T * O) == sp.eye(8), "original O is orthogonal")

    mapped_doubled = {
        tuple(int(value) for value in target_basis * sp.Matrix(coordinates))
        for coordinates in lattice["root_coordinates"]
    }
    need(mapped_doubled == cs._standard_doubled_roots(), "O maps all original E8 roots")

    # In the marked D5+A3 frame, the glue degree is the first-five sum mod 4.
    grade = lambda root: sum(root[:5]) % 4
    grade_counts = Counter(grade(root) for root in mapped_doubled)
    need([grade_counts[k] for k in range(4)] == [52, 64, 60, 64],
         "marked C4 root grading")
    canonical_roots = {
        tuple(Q(value, 2) for value in root)
        for root in mapped_doubled if grade(root) in (0, 2)
    }
    d8 = standard_d8_roots()
    need(canonical_roots == d8, "order-two glue intermediate is standard target D8")

    # Lattice-level parity identity K union (K+2g)=D8, not only a root count.
    def in_k(bits):
        return sum(bits[:5]) % 2 == 0 and sum(bits[5:]) % 2 == 0

    def in_k_plus_2g(bits):
        shifted = tuple(bit ^ 1 for bit in bits)
        return in_k(shifted)

    for bits in product((0, 1), repeat=8):
        need((in_k(bits) or in_k_plus_2g(bits)) == (sum(bits) % 2 == 0),
             "K plus the unique order-two glue coset equals D8")
    # Parentheses make the intended equivalence explicit after the direct check.
    parity_identity = all(
        (in_k(bits) or in_k_plus_2g(bits)) == (sum(bits) % 2 == 0)
        for bits in product((0, 1), repeat=8)
    )
    need(parity_identity, "full parity identity")

    fold_data = json.loads((ROOT / "experiments/theory-contracts/"
        "v16-operator-pullback-20261014/fold-results.json").read_text())
    M = tuple(tuple(Q(value) for value in row)
              for row in fold_data["explicit_Weyl_matrix"])
    MT = tuple(zip(*M))
    need(all(dot(M[i], M[j]) == (1 if i == j else 0)
             for i in range(8) for j in range(8)), "original M is orthogonal")
    e8 = e8_roots()
    need({mat_vec(M, root) for root in e8} == e8, "original M preserves E8 roots")

    folded_roots = {mat_vec(M, root) for root in d8}
    intersection = canonical_roots & folded_roots
    need(len(folded_roots) == 112 and len(intersection) == 48,
         "two exact D8 frames meet in 48 roots")
    need(canonical_roots != folded_roots, "decisive D8 equality must fail")

    canonical_only = min(canonical_roots - folded_roots)
    folded_only = min(folded_roots - canonical_roots)
    need(mat_vec(MT, canonical_only) not in d8, "canonical-only witness")
    need(mat_vec(MT, folded_only) in d8, "fold-only witness")

    # The two D8s are the +1 root sectors of two distinct torus characters.
    h_canonical = tuple([Q(1)] + [Q(0)] * 7)
    h_folded = tuple(row[0] for row in M)  # M e_0
    joint = Counter((phase_sign(h_canonical, root), phase_sign(h_folded, root))
                    for root in e8)
    need(joint == {(1, 1): 48, (1, -1): 64, (-1, 1): 64, (-1, -1): 64},
         "distinct charged-root actions")
    need({root for root in e8 if phase_sign(h_canonical, root) == 1} == canonical_roots,
         "canonical defect fixes canonical D8")
    need({root for root in e8 if phase_sign(h_folded, root) == 1} == folded_roots,
         "fold defect fixes transported CAR D8")

    # The source-frame spinor singled out by the quarter shift. It belongs
    # to the standard E8 spinor neighbour D8+s, but not to the opposite
    # odd-parity convention called c in ODD-SECTOR-REFIELDING-AND-NAHT.
    p_spinor = tuple([Q(1, 2)] * 4 + [Q(-1, 2)] * 4)
    s_spinor = tuple([Q(1, 2)] * 8)
    p_minus_s = tuple(a - b for a, b in zip(p_spinor, s_spinor))
    need(p_spinor in e8 and p_spinor not in d8, "source spinor is an E8 root outside CAR D8")
    need(all(value.denominator == 1 for value in p_minus_s), "spinor coset representative")
    need(sum(p_minus_s) % 2 == 0, "source spinor lies in the even-parity E8 spinor coset")
    q_target = mat_vec(M, p_spinor)
    q_expected = tuple([Q(1), Q(0), Q(0), Q(1)] + [Q(0)] * 4)
    need(q_target == q_expected, "fixed M sends p to the already known e1+e4 line")
    need(q_target in canonical_roots and q_target not in folded_roots,
         "the transported spinor is a canonical-only D8 root")

    time_cases = {}
    for alpha in (Q(1, 4), Q(3, 4)):
        delta = tuple(Q(k, 8) + alpha / 8 - Q(1, 2) for k in range(8))
        d = mat_vec(M, delta)
        canonical_phases = Counter(dot(d, root) % 1 for root in canonical_roots)
        folded_phases = Counter(dot(d, root) % 1 for root in folded_roots)
        source_phases = Counter(dot(delta, root) % 1 for root in d8)
        need(folded_phases == source_phases, "M transports the actual quarter time")
        need(canonical_phases != folded_phases, "quarter time distinguishes the two D8 frames")
        source_spinor_energy = dot(p_spinor, p_spinor) / 2 + dot(delta, p_spinor)
        target_root_energy = dot(q_target, q_target) / 2 + dot(d, q_target)
        need(dot(delta, p_spinor) == -1, "quarter shift pairing with p")
        need(source_spinor_energy == 0 and target_root_energy == 0,
             "M transports the zero-energy spinor line")
        time_cases[fstr(alpha)] = {
            "delta_in_folded_CAR_frame": vecstr(delta),
            "d_equals_M_delta_in_target_frame": vecstr(d),
            "canonical_D8_phase_counts": counterstr(canonical_phases),
            "transported_fold_D8_phase_counts": counterstr(folded_phases),
            "canonical_only_witness_pairing": fstr(dot(d, canonical_only)),
            "fold_only_witness_pairing": fstr(dot(d, folded_only)),
            "source_spinor_p_pairing": fstr(dot(delta, p_spinor)),
            "source_spinor_shifted_energy": fstr(source_spinor_energy),
            "target_root_q_pairing": fstr(dot(d, q_target)),
            "target_root_shifted_energy": fstr(target_root_energy),
        }

    return {
        "status": "SCOPED_FALSE_CANONICAL_D8_IS_NOT_THE_M_TRANSPORTED_CAR_D8",
        "marked_C4_glue": {
            "generator": "g=(1/2)^8 modulo K=D5+A3",
            "unique_order_two_subgroup": "<2g>",
            "lattice_identity": "K union (K+2g)=D8",
            "root_grade_counts_0_1_2_3": [grade_counts[k] for k in range(4)],
            "intermediate_root_count": len(canonical_roots),
            "O_maps_original_E8_to_target_frame": True,
        },
        "frame_comparison": {
            "canonical_D8_roots": len(canonical_roots),
            "M_transported_CAR_D8_roots": len(folded_roots),
            "intersection": len(intersection),
            "canonical_only": len(canonical_roots - folded_roots),
            "fold_only": len(folded_roots - canonical_roots),
            "equal": False,
            "canonical_only_root": vecstr(canonical_only),
            "canonical_only_fold_preimage": vecstr(mat_vec(MT, canonical_only)),
            "fold_only_root": vecstr(folded_only),
            "fold_only_fold_preimage": vecstr(mat_vec(MT, folded_only)),
        },
        "legitimacy_boundary": {
            "fixed_simultaneous_marking": "The identity comparison with both supplied O and M frozen is false.",
            "common_Weyl_transport": "A single Weyl map applied to the D5+A3 flag, O frame, time vector, source embedding, field words and cocycles is a coordinate-equivalent presentation.",
            "partial_alignment": "Rotating only one D8 while retaining the other marked data is a new identification and is not licensed by abstract D8 conjugacy.",
            "TFPT_no_go": False,
        },
        "Cartan_character_comparison": {
            "canonical_character_vector": vecstr(h_canonical),
            "fold_character_vector_Me0": vecstr(h_folded),
            "joint_root_sign_counts": {f"{a},{b}": n for (a, b), n in sorted(joint.items())},
            "both_act_trivially_on_the_eight_Cartan_currents": True,
            "charged_root_actions_equal": False,
        },
        "quarter_time": time_cases,
        "spinor_zero_ground": {
            "source_spinor_p": vecstr(p_spinor),
            "norm_squared": fstr(dot(p_spinor, p_spinor)),
            "p_minus_standard_spinor_representative": vecstr(p_minus_s),
            "in_standard_E8_even_spinor_coset": True,
            "in_ODD_SECTOR_odd_parity_c_convention": False,
            "in_bare_integer_CAR_charge_lattice": False,
            "M_p_target": vecstr(q_target),
            "target_line_one_based_name": "C e^(e1+e4)",
            "in_canonical_D8_but_not_M_transported_CAR_D8": True,
            "already_known_equal_energy_target_root": True,
            "consequence": "The full E8 charge extension has an additional zero-energy line absent from the bare integer CAR charge lattice; state/time compatibility requires an explicit sector rule.",
        },
        "gated_followups": {
            "nonisotropic_vector_coset_and_odd_neighbor_Z8": "NOT_RUN_AFTER_FAILED_EQUALITY_GATE",
            "E8_Z2_defect_sector_v_plus_c": "NOT_IDENTIFIED_ACROSS_THE_TWO_FRAMES",
            "reason": "Those discriminant labels belong to a chosen D8 embedding; the two actual embeddings are unequal.",
        },
        "sources": {
            "O_and_marked_glue": {
                "path": "tfpt_explorer/cartan_source.py",
                "sha256": hashlib.sha256((ROOT / "tfpt_explorer/cartan_source.py").read_bytes()).hexdigest(),
            },
            "M": {
                "path": "experiments/theory-contracts/v16-operator-pullback-20261014/fold-results.json",
                "sha256": hashlib.sha256((ROOT / "experiments/theory-contracts/"
                    "v16-operator-pullback-20261014/fold-results.json").read_bytes()).hexdigest(),
            },
            "known_locality_exclusion": {
                "path": "experiments/theory-contracts/mixed-current-origin-20261014/ODD-SECTOR-REFIELDING-AND-NAHT.txt",
                "sha256": hashlib.sha256((ROOT / "experiments/theory-contracts/"
                    "mixed-current-origin-20261014/ODD-SECTOR-REFIELDING-AND-NAHT.txt").read_bytes()).hexdigest(),
            },
            "known_equal_energy_target_root": {
                "path": "experiments/theory-contracts/joint-native-algebra-selection-20261016/FOLD-COVARIANCE-TO-NATIVE-FIELD-RESPONSE.txt",
                "atlas_source_id": "research-source:1a8c7c3f539e779301a7377d",
                "sha256": hashlib.sha256((ROOT / "experiments/theory-contracts/"
                    "joint-native-algebra-selection-20261016/FOLD-COVARIANCE-TO-NATIVE-FIELD-RESPONSE.txt").read_bytes()).hexdigest(),
            },
        },
    }


if __name__ == "__main__":
    print(json.dumps(calculate(), indent=2, sort_keys=True))
