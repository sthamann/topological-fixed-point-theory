#!/usr/bin/env python3
"""Exact two leg Choi witness for the Hodge and old twenty event channels.

The old Kraus matrices are reconstructed by the pinned original source.  The
Hodge matrices and the six duad charts come from the pinned preceding local
checker.  This file tests an operator statement in that common chart.  It
does not add an instrument, a preparation rule, or a physical identification
of the old Bell module with an SU(5) module.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import itertools as it
import json
from fractions import Fraction
from pathlib import Path

import numpy as np
import sympy as sp


HERE = Path(__file__).resolve().parent
WORK = HERE.parent
REPO = Path("/Users/stefanhamann/Projekte/tfpt-theoryv4")
OLD_SOURCE = REPO / "experiments/theory-contracts/compiler-correlated-event-clock-20260919/source_channel.py"
OLD_PROOF = REPO / "experiments/theory-contracts/compiler-correlated-event-clock-20260919/PROOF.txt"
HODGE_DIR = WORK / "research-20261028-hodge-complement-channel"
if not HODGE_DIR.is_dir():
    HODGE_DIR = REPO / "experiments/theory-contracts/elliptic-moebius-source-bridge-20261026/hodge_complement_20261028"
HODGE_CHECKER = HODGE_DIR / "hodge_complement_channel.py"
HODGE_RESULT = HODGE_DIR / "result.json"
A4_REVIEW = HODGE_DIR / "E8-A4-PARENT-REVIEW.txt"
PINS = {
    str(OLD_SOURCE): "3852fabaf10c143fe20238b67406df6396f85c931dfd4834fe91a2d7ef900593",
    str(OLD_PROOF): "bf1f85c99f5f89a0a15da1af5e3ff321c2756fda132a9f8d644f9ce0e08c2f5c",
    str(HODGE_CHECKER): "72b04a5bf8b3bd805568af49dec76bd97926fcff4009a830137e2d318807522b",
    str(HODGE_RESULT): "811d82cbc63efee9551887ba3ddecb70c7330934d322069ae19f10de25025062",
    str(A4_REVIEW): "7199b91ec2f24c4d17a5519758c49e19a11a4bac7fd59b34f134a38de15eb221",
}

checks: list[str] = []


def require(ok: bool, name: str) -> None:
    if not bool(ok):
        raise RuntimeError(name)
    checks.append(name)


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, f"import specification {name}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def exterior_two_differential(X: np.ndarray, pairs, pair_index) -> np.ndarray:
    """Differential of Lambda squared on a Gaussian integer matrix."""
    out = np.zeros((10, 10), dtype=np.complex128)
    for col, (i, j) in enumerate(pairs):
        for a in range(5):
            if X[a, i] and a != j:
                pair = tuple(sorted((a, j)))
                out[pair_index[pair], col] += (1 if a < j else -1) * X[a, i]
            if X[a, j] and a != i:
                pair = tuple(sorted((i, a)))
                out[pair_index[pair], col] += (1 if i < a else -1) * X[a, j]
    return out


def su5_generators() -> list[tuple[str, np.ndarray]]:
    """A Gaussian integer anti Hermitian basis of su(5)."""
    generators = []
    for a, b in it.combinations(range(5), 2):
        X = np.zeros((5, 5), dtype=np.complex128)
        X[a, b], X[b, a] = 1, -1
        generators.append((f"A{a}{b}", X))
        X = np.zeros((5, 5), dtype=np.complex128)
        X[a, b], X[b, a] = 1j, 1j
        generators.append((f"S{a}{b}", X))
    for a in range(4):
        X = np.zeros((5, 5), dtype=np.complex128)
        X[a, a], X[4, 4] = 1j, -1j
        generators.append((f"D{a}", X))
    require(len(generators) == 24, "24 anti Hermitian su5 generators")
    return generators


def sm_generators() -> list[tuple[str, np.ndarray]]:
    """One fixed su(3)+su(2)+u(1) inside su(5)."""
    generators = []
    for a, b in it.combinations(range(3), 2):
        X = np.zeros((5, 5), dtype=np.complex128)
        X[a, b], X[b, a] = 1, -1
        generators.append((f"su3-A{a}{b}", X))
        X = np.zeros((5, 5), dtype=np.complex128)
        X[a, b], X[b, a] = 1j, 1j
        generators.append((f"su3-S{a}{b}", X))
    X = np.zeros((5, 5), dtype=np.complex128)
    X[3, 4], X[4, 3] = 1, -1
    generators.append(("su2-A34", X))
    X = np.zeros((5, 5), dtype=np.complex128)
    X[3, 4], X[4, 3] = 1j, 1j
    generators.append(("su2-S34", X))
    generators.extend([
        ("su3-D0", np.diag([1j, 0, -1j, 0, 0]).astype(np.complex128)),
        ("su3-D1", np.diag([0, 1j, -1j, 0, 0]).astype(np.complex128)),
        ("su2-D", np.diag([0, 0, 0, 1j, -1j]).astype(np.complex128)),
        ("Y", np.diag([-2j, -2j, -2j, 3j, 3j]).astype(np.complex128)),
    ])
    require(len(generators) == 12, "12 generators of one fixed SM subgroup")
    return generators


def two_leg_generator(X: np.ndarray, pairs, pair_index) -> np.ndarray:
    r = exterior_two_differential(X, pairs, pair_index)
    eye = np.eye(10, dtype=np.complex128)
    return np.kron(r, eye) + np.kron(eye, r)


def exact_complex_integer(a: np.ndarray) -> np.ndarray:
    require(np.array_equal(a.real, np.rint(a.real)), "Gaussian integer real part")
    require(np.array_equal(a.imag, np.rint(a.imag)), "Gaussian integer imaginary part")
    return np.rint(a.real).astype(np.int64) + 1j * np.rint(a.imag).astype(np.int64)


def old_kraus_for_mark(old, hodge, mark, bell, six, mapping):
    idx, *_ = old.marked_data(mark, bell, six, mapping)
    chart, rename = hodge.duad_chart(mark, mapping)
    groups = {m: [] for m in range(5)}
    all_kraus = []
    for k in idx:
        U = bell[k]
        permutation = old.bell_permutation(U)
        jump_support = np.eye(10) - np.diag([int(permutation[a] == a) for a in range(10)])
        K = exact_complex_integer(chart.T @ (U @ jump_support) @ chart)
        m = rename[six[k][mark]]
        groups[m].append(K)
        all_kraus.append(K)
    require(len(all_kraus) == 20, f"mark {mark}: twenty actual old Kraus matrices")
    require(all(len(groups[m]) == 4 for m in range(5)), f"mark {mark}: four lifts per complement")
    return groups, all_kraus


def partial_gains(kraus: list[np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    left = sum((K @ K.conj().T for K in kraus), start=np.zeros((10, 10), dtype=np.complex128))
    right = sum((K.conj().T @ K for K in kraus), start=np.zeros((10, 10), dtype=np.complex128))
    return exact_complex_integer(left), exact_complex_integer(right)


def main(output: Path | None) -> dict:
    checks.clear()
    for path, digest in PINS.items():
        require(hashlib.sha256(Path(path).read_bytes()).hexdigest() == digest, "source pin " + path)

    old = load_module("pinned_old_source", OLD_SOURCE)
    hodge = load_module("pinned_hodge_checker", HODGE_CHECKER)
    Cs = hodge.hodge_amplitudes()
    V = np.stack([C.reshape(-1, order="F") for C in Cs], axis=1)
    require(np.array_equal(V.T @ V, 6 * np.eye(5, dtype=np.int64)), "Hodge Choi vectors have Gram 6I5")
    P5_numerator = V @ V.T
    require(np.array_equal(P5_numerator @ P5_numerator, 6 * P5_numerator), "P5 is a projector after division by six")
    require(int(sp.Matrix(P5_numerator).rank()) == 5, "P5 has rank five")
    identity100 = np.eye(100, dtype=np.int64)
    swap = np.zeros((100, 100), dtype=np.int64)
    for a in range(10):
        for b in range(10):
            swap[b + 10 * a, a + 10 * b] = 1
    # Q_d=6P_d keeps every check integral.
    Q5 = P5_numerator
    Q45 = 3 * (identity100 - swap)
    Q50 = 3 * (identity100 + swap) - Q5
    sector_data = [(5, Q5), (45, Q45), (50, Q50)]
    for dimension, Q in sector_data:
        require(np.array_equal(Q @ Q, 6 * Q), f"P{dimension} exact projector")
        require(int(np.trace(Q)) == 6 * dimension, f"P{dimension} exact rank")
    require(np.array_equal(Q5 + Q45 + Q50, 6 * identity100), "P5 plus P45 plus P50 is identity")
    require(not np.any(Q5 @ Q45) and not np.any(Q5 @ Q50) and not np.any(Q45 @ Q50),
            "three SU5 projectors are pairwise orthogonal")

    full_generators = su5_generators()
    for name, X in full_generators:
        L = two_leg_generator(X, hodge.PAIRS, hodge.PAIR_INDEX)
        require(np.array_equal(L @ P5_numerator, P5_numerator @ L), f"P5 commutes with su5 generator {name}")
    sm = sm_generators()
    for name, X in sm:
        L = two_leg_generator(X, hodge.PAIRS, hodge.PAIR_INDEX)
        require(np.array_equal(L @ P5_numerator, P5_numerator @ L), f"P5 commutes with fixed SM generator {name}")

    # J_Phi=(1/3)VV^dagger=2P5 and rho_Phi=J_Phi/10.
    # Therefore Tr(P5 rho_Phi)=Tr(P5)/5=1.
    require(np.trace(P5_numerator) == 30, "trace P5 numerator is thirty")
    p_phi = Fraction(int(np.trace(P5_numerator @ P5_numerator)), 6 * 3 * 10)
    require(p_phi == 1, "normalized Hodge Choi witness is one")
    phi_sector_weights = [Fraction(int(np.trace(Q @ P5_numerator)), 6 * 3 * 10)
                          for _dimension, Q in sector_data]
    require(phi_sector_weights == [1, 0, 0], "Hodge Choi lies entirely in P5")
    phi_left, phi_right = partial_gains(Cs)
    require(np.array_equal(phi_left, 3 * np.eye(10)), "Phi unital gain 3I")
    require(np.array_equal(phi_right, 3 * np.eye(10)), "Phi trace preserving gain 3I")

    rays = old.source_rays()
    bell, raw_six, _ = old.reflection_actions(rays)
    six, _ = old.outer_mark_action(raw_six)
    mapping = old.find_bell_partition_intertwiner(bell, six)
    mark_results = []
    for mark in range(6):
        groups, kraus = old_kraus_for_mark(old, hodge, mark, bell, six, mapping)
        left, right = partial_gains(kraus)
        require(np.array_equal(left, 12 * np.eye(10)), f"mark {mark}: Psi unital gain 12I")
        require(np.array_equal(right, 12 * np.eye(10)), f"mark {mark}: Psi trace preserving gain 12I")

        event_overlaps = []
        for K in kraus:
            amplitudes = V.conj().T @ K.reshape(-1, order="F")
            event_overlaps.append(int(np.vdot(amplitudes, amplitudes).real))
        require(sorted(event_overlaps) == [4] * 18 + [36] * 2,
                f"mark {mark}: exact twenty event P5 overlap distribution")
        overlap_sum = sum(event_overlaps)
        p_psi = Fraction(overlap_sum, 6 * 12 * 10)
        require(overlap_sum == 144, f"mark {mark}: total P5 numerator 144")
        require(p_psi == Fraction(1, 5), f"mark {mark}: normalized old Choi witness is one fifth")

        group_sums = [sum(int(np.vdot(V.conj().T @ K.reshape(-1, order="F"),
                                         V.conj().T @ K.reshape(-1, order="F")).real)
                          for K in groups[m]) for m in range(5)]
        require(sorted(group_sums) == [16, 16, 16, 48, 48],
                f"mark {mark}: complement resolved overlap sums")

        unscaled_choi = sum((np.outer(K.reshape(-1, order="F"),
                                              K.reshape(-1, order="F").conj()) for K in kraus),
                            start=np.zeros((100, 100), dtype=np.complex128))
        unscaled_choi = exact_complex_integer(unscaled_choi)
        sector_weights = []
        projected_ranks = []
        kraus_vectors = np.stack([K.reshape(-1, order="F") for K in kraus], axis=1)
        for dimension, Q in sector_data:
            numerator = np.trace(Q @ unscaled_choi)
            require(numerator.imag == 0, f"mark {mark}: P{dimension} weight is real")
            sector_weights.append(Fraction(int(numerator.real), 6 * 12 * 10))
            projected_ranks.append(int(sp.Matrix(Q @ kraus_vectors).rank()))
        require(sector_weights == [Fraction(1, 5), Fraction(2, 5), Fraction(2, 5)],
                f"mark {mark}: actual SU5 twirl weights one fifth two fifths two fifths")
        symmetric_rank = int(sp.Matrix((identity100 + swap) @ kraus_vectors).rank())
        antisymmetric_rank = int(sp.Matrix((identity100 - swap) @ kraus_vectors).rank())
        require([symmetric_rank, antisymmetric_rank] == [9, 6],
                f"mark {mark}: old Choi support has nine symmetric and six antisymmetric directions")

        # The SU5 twirl is sum_d p_d (10/d)P_d.  Its population
        # entries depend only on whether two duads are equal, intersecting,
        # or disjoint.  Evaluate every one of the hundred entries.
        population_values = {"stay": set(), "intersect": set(), "disjoint": set()}
        for a, e in enumerate(hodge.PAIRS):
            for b, f in enumerate(hodge.PAIRS):
                kind = "stay" if a == b else "disjoint" if set(e).isdisjoint(f) else "intersect"
                coordinate = b + 10 * a
                value = sum((weight * Fraction(10, dimension) * Fraction(int(Q[coordinate, coordinate]), 6)
                             for weight, (dimension, Q) in zip(sector_weights, sector_data)),
                            start=Fraction(0))
                population_values[kind].add(value)
        require(population_values == {
                    "stay": {Fraction(2, 25)},
                    "intersect": {Fraction(19, 225)},
                    "disjoint": {Fraction(31, 225)},
                }, f"mark {mark}: exact SU5 twirl population leakage")
        require(Fraction(2, 25) + 6 * Fraction(19, 225) + 3 * Fraction(31, 225) == 1,
                f"mark {mark}: twirled population column normalization")
        full_commuting = 0
        for _name, X in full_generators:
            L = two_leg_generator(X, hodge.PAIRS, hodge.PAIR_INDEX)
            full_commuting += int(np.array_equal(L @ unscaled_choi, unscaled_choi @ L))
        sm_commuting = 0
        for _name, X in sm:
            L = two_leg_generator(X, hodge.PAIRS, hodge.PAIR_INDEX)
            sm_commuting += int(np.array_equal(L @ unscaled_choi, unscaled_choi @ L))
        require(full_commuting < 24, f"mark {mark}: old Choi matrix is not fully SU5 invariant")
        require(sm_commuting < 12, f"mark {mark}: old Choi matrix is not invariant under fixed SM subgroup")

        mark_results.append({
            "mark": mark,
            "p_Psi": "1/5",
            "su5_twirl_weights": {"P5": "1/5", "P45": "2/5", "P50": "2/5"},
            "projected_kraus_span_ranks": {"P5": projected_ranks[0], "P45": projected_ranks[1], "P50": projected_ranks[2]},
            "symmetric_antisymmetric_choi_support_ranks": [symmetric_rank, antisymmetric_rank],
            "su5_twirl_population_per_target": {kind: str(next(iter(values)))
                                                  for kind, values in population_values.items()},
            "su5_twirl_population_total": {
                "stay": "2/25",
                "intersect": "38/75",
                "disjoint": "31/75",
            },
            "event_overlap_squared_distribution": {"4": 18, "36": 2},
            "total_overlap_numerator": overlap_sum,
            "complement_group_overlap_sums": group_sums,
            "normalized_choi_marginals": "I_10/10 on both legs",
            "old_choi_su5_generators_commuting": f"{full_commuting}/24",
            "old_choi_fixed_sm_generators_commuting": f"{sm_commuting}/12",
        })

    require(all(item["p_Psi"] == "1/5" for item in mark_results), "all six marks give p_Psi one fifth")
    result = {
        "verdict": [
            "P5_IS_A_RANK_FIVE_FULL_SU5_INVARIANT_TWO_LEG_PROJECTOR",
            "NORMALIZED_HODGE_CHOI_HAS_P5_PROBABILITY_ONE",
            "ACTUAL_OLD_TWENTY_EVENT_CHOI_HAS_P5_PROBABILITY_ONE_FIFTH_FOR_ALL_SIX_MARKS",
            "ACTUAL_OLD_CHOI_SU5_TWIRL_WEIGHTS_ARE_ONE_FIFTH_TWO_FIFTHS_TWO_FIFTHS",
            "SU5_TWIRL_INTRODUCES_EXACT_STAY_AND_INTERSECTION_POPULATION_LEAKAGE",
            "TWO_LEG_SU5_AND_FIXED_SM_NEUTRAL_OBSERVABLES_CAN_DISTINGUISH_PHASE_ERASURE",
            "ALL_SINGLE_LEG_OBSERVABLES_HAVE_EQUAL_MAXIMALLY_MIXED_CHOI_MARGINALS",
            "OLD_BELL10_TO_DUAD_CHART_IS_NOT_A_PHYSICAL_SU5_MODULE_IDENTIFICATION",
            "NO_PREPARATION_OR_EVENT_SELECTION_RULE_IS_DERIVED",
        ],
        "source_sha256": PINS,
        "witness": {
            "definition": "P5=(1/6) sum_m |C_m>><<C_m|",
            "space_in_fixed_chart": "Lambda^2 C^5 tensor Lambda^2 C^5, equivalently the stated dual bar10 tensor bar10 convention",
            "rank": 5,
            "projector": True,
            "su5_invariant": True,
            "su5_generators_checked": 24,
            "fixed_sm_invariant": True,
            "fixed_sm_generators_checked": 12,
        },
        "normalized_choi_witness": {
            "p_Phi": "1",
            "p_Psi_all_six_marks": "1/5",
            "derivation_Phi": "J_Phi=(1/3)VV^dagger=2P5 and rho_Phi=J_Phi/10",
            "derivation_Psi": "sum_l ||V^dagger vec(K_l)||^2=144, divided by 6*12*10",
        },
        "actual_old_su5_twirl": {
            "weights": {"P5": "1/5", "P45": "2/5", "P50": "2/5"},
            "choi_support_ranks_under_swap": {"symmetric": 9, "antisymmetric": 6},
            "projected_span_ranks": {"P5": 5, "P45": 6, "P50": 7},
            "population_per_target": {"stay": "2/25", "intersect": "19/225", "disjoint": "31/225"},
            "population_total_per_input": {"stay": "2/25", "intersect": "38/75", "disjoint": "31/75"},
            "interpretation": "twirling the noncovariant old channel preserves the P5 weight but creates non-Petersen stay and intersecting transitions",
        },
        "commutant_scope": {
            "two_leg_full_su5_commutant_restrictions_equal": False,
            "two_leg_fixed_sm_commutant_restrictions_equal": False,
            "reason": "P5 belongs to both commutants and its two normalized Choi expectations are 1 and 1/5",
            "single_leg_restrictions_equal": True,
            "single_leg_reason": "both normalized Choi marginals are I_10/10 because both channels are trace preserving and unital",
            "old_choi_itself_su5_invariant": False,
            "old_choi_itself_fixed_sm_invariant": False,
        },
        "marks": mark_results,
        "scope": {
            "operator_claim": "exact in the fixed duad chart used by the old source comparison",
            "physical_module_claim": False,
            "caveat": "the old Bell10 is Sym^2 C^4; the event label duad chart alone is not a linear SU5 intertwiner",
            "preparation_claim": False,
        },
        "checks": len(checks),
    }
    if output is not None:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    print(json.dumps(main(args.out), indent=2, sort_keys=True))
