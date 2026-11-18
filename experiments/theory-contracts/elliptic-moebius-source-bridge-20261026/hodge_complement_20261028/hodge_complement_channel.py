#!/usr/bin/env python3
"""Exact finite checks for the five separated Hodge complement amplitudes.

The calculation compares the resulting CPTP map with the *actual* marked
twenty-event jump map reconstructed by the pinned source contract.  It proves
finite exterior-algebra and channel identities only; it does not select a
physical TFPT instrument.
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
OLD_DIR = REPO / "experiments/theory-contracts/compiler-correlated-event-clock-20260919"
OLD_SOURCE = OLD_DIR / "source_channel.py"
OLD_PROOF = OLD_DIR / "PROOF.txt"
PREVIOUS = WORK / "research-20261027-dynamic-reconstruction/HODGE-BORN-SOURCE-TEST.txt"
V810 = REPO / "verification/v810_cutcode_moment.py"
CURRENT_DESCENDANT_PROOF = REPO / "experiments/theory-contracts/compiler-current-descendant-20260919/PROOF.txt"
CURRENT_PRODUCT_PROOF = REPO / "experiments/theory-contracts/compiler-current-product-20260919/PROOF.txt"
PINS = {
    str(OLD_SOURCE): "3852fabaf10c143fe20238b67406df6396f85c931dfd4834fe91a2d7ef900593",
    str(OLD_PROOF): "bf1f85c99f5f89a0a15da1af5e3ff321c2756fda132a9f8d644f9ce0e08c2f5c",
    str(PREVIOUS): "49d0c4fe47baa2284e95d5600a8591d614b0ae24266443ea824a8bf98e9d3344",
    str(V810): "8dcf729eb9be10b4d51ecdc61f41624f8cf87a991cd5a2c72eaeb2269212afae",
    str(CURRENT_DESCENDANT_PROOF): "ad45198c44b84365a862b709e8b7e49af708f2e1a58adff67dc55d7cf4fa459b",
    str(CURRENT_PRODUCT_PROOF): "940131acfceba1937d36b44447ecce39b25571ef7c95e870667b838f4e9ae0db",
}

PAIRS = list(it.combinations(range(5), 2))
PAIR_INDEX = {p: i for i, p in enumerate(PAIRS)}
checks: list[str] = []


def require(ok: bool, name: str) -> None:
    if not bool(ok):
        raise RuntimeError(name)
    checks.append(name)


def parity(word: tuple[int, ...]) -> int:
    inversions = sum(word[i] > word[j] for i in range(len(word)) for j in range(i + 1, len(word)))
    return -1 if inversions % 2 else 1


def hodge_amplitudes() -> list[np.ndarray]:
    """C_m(e_i wedge e_j)=epsilon(m,i,j,k,l)e_k wedge e_l."""
    out = []
    universe = set(range(5))
    for m in range(5):
        C = np.zeros((10, 10), dtype=np.int64)
        for e in PAIRS:
            if m in e:
                continue
            f = tuple(sorted(universe - {m, *e}))
            C[PAIR_INDEX[f], PAIR_INDEX[e]] = parity((m,) + e + f)
        out.append(C)
    return out


def exterior_two_differential(X: np.ndarray) -> np.ndarray:
    """Differential of Lambda^2 on an integer 5 by 5 matrix."""
    r = np.zeros((10, 10), dtype=np.int64)
    for col, (i, j) in enumerate(PAIRS):
        for a in range(5):
            if X[a, i] and a != j:
                pair = tuple(sorted((a, j)))
                sign = 1 if a < j else -1
                r[PAIR_INDEX[pair], col] += sign * int(X[a, i])
            if X[a, j] and a != i:
                pair = tuple(sorted((i, a)))
                sign = 1 if i < a else -1
                r[PAIR_INDEX[pair], col] += sign * int(X[a, j])
    return r


def verify_sl5_covariance(Cs: list[np.ndarray]) -> int:
    """Check the infinitesimal epsilon-congruence identity on a basis of sl5."""
    generators = []
    for a in range(5):
        for b in range(5):
            if a != b:
                X = np.zeros((5, 5), dtype=np.int64)
                X[a, b] = 1
                generators.append(X)
    for a in range(4):
        X = np.zeros((5, 5), dtype=np.int64)
        X[a, a] = 1
        X[4, 4] = -1
        generators.append(X)
    require(len(generators) == 24, "24 element basis of sl5")
    for gi, X in enumerate(generators):
        r = exterior_two_differential(X)
        for m in range(5):
            lhs = r @ Cs[m] + Cs[m] @ r.T
            rhs = -sum((int(X[m, n]) * Cs[n] for n in range(5)), start=np.zeros((10, 10), dtype=np.int64))
            require(np.array_equal(lhs, rhs), f"sl5 epsilon congruence generator {gi} index {m}")
    return len(generators)


def verify_su5_sector_spectrum(Cs: list[np.ndarray], Snew: np.ndarray) -> dict:
    """Identify the 1+24+75 decomposition, not merely its dimensions."""
    generators = []
    for a in range(5):
        for b in range(5):
            if a != b:
                X = np.zeros((5, 5), dtype=np.int64)
                X[a, b] = 1
                generators.append(X)
    for a in range(4):
        X = np.zeros((5, 5), dtype=np.int64)
        X[a, a] = 1
        X[4, 4] = -1
        generators.append(X)
    adjoint = [exterior_two_differential(X) for X in generators]
    adjoint_columns = np.stack([A.reshape(-1, order="F") for A in adjoint], axis=1)
    require(int(sp.Matrix(adjoint_columns).rank()) == 24,
            "dGamma(sl5) is the actual 24-dimensional adjoint subspace")
    identity_vector = np.eye(10, dtype=np.int64).reshape(-1, order="F")
    require(all(np.dot(identity_vector, column) == 0 for column in adjoint_columns.T),
            "adjoint subspace is traceless and orthogonal to identity")
    for i, A in enumerate(adjoint):
        image = (Snew @ A.reshape(-1, order="F")).reshape((10, 10), order="F")
        require(np.array_equal(image, -2 * A.T),
                f"Phi sends dGamma generator {i} to minus two thirds transpose")

    square = Snew @ Snew
    require(np.array_equal(square @ identity_vector, 9 * identity_vector),
            "scaled Phi squared eigenvalue nine on identity")
    require(all(np.array_equal(square @ column, 4 * column) for column in adjoint_columns.T),
            "scaled Phi squared eigenvalue four on actual adjoint")
    identity_adjoint = np.column_stack((identity_vector, adjoint_columns))
    require(int(sp.Matrix(identity_adjoint).rank()) == 25,
            "identity plus actual adjoint has dimension 25")
    require(int(sp.Matrix(square - np.eye(100, dtype=np.int64)).rank()) == 25,
            "orthogonal complement is the 75-dimensional Phi-squared one-eigenspace")

    # The actual population restriction to distance shells about one Petersen
    # vertex has shell sizes 1,3,6 and the signed three-rate quotient.
    adjacency = np.zeros((10, 10), dtype=np.int64)
    for a, e in enumerate(PAIRS):
        for b, f in enumerate(PAIRS):
            adjacency[b, a] = int(a != b and set(e).isdisjoint(f))
    shell0 = [0]
    shell1 = list(np.flatnonzero(adjacency[:, 0]))
    shell2 = [i for i in range(10) if i not in shell0 + shell1]
    shells = [shell0, shell1, shell2]
    require([len(shell) for shell in shells] == [1, 3, 6], "Petersen radial shell sizes 1,3,6")
    quotient = np.zeros((3, 3), dtype=np.int64)
    for source_shell, source_vertices in enumerate(shells):
        reference_counts = None
        for vertex in source_vertices:
            counts = tuple(int(adjacency[target_shell, vertex].sum()) for target_shell in shells)
            if reference_counts is None:
                reference_counts = counts
            require(counts == reference_counts, "Petersen radial partition is equitable")
        quotient[:, source_shell] = reference_counts
    require(np.array_equal(quotient, np.array([[0, 1, 0], [3, 0, 1], [0, 2, 2]], dtype=np.int64)),
            "Petersen radial quotient is the exact three-shell transfer")
    require(multiplicities_from_moments(quotient, (3, 1, -2)) == [1, 1, 1],
            "Petersen radial quotient has scaled rates 3,1,-2")
    return {
        "endomorphism_decomposition": "End(Lambda^2 C^5)=1+24+75",
        "adjoint_identification": "24=dGamma(sl5), with Phi(A)=-(2/3) A^T",
        "seventy_five_identification": "orthogonal complement of identity plus dGamma(sl5)",
        "Phi_squared": "Pi_1+(4/9)Pi_24+(1/9)Pi_75",
        "Phi_sixth_spectrum": {"1": 1, "64/729": 24, "1/729": 75},
        "radial_population_shells": [1, 3, 6],
        "radial_population_transfer": "(1/3)*[[0,1,0],[3,0,1],[0,2,2]] in shell-mass column convention",
        "radial_spectrum": ["1", "1/3", "-2/3"],
    }


def load_old_source():
    spec = importlib.util.spec_from_file_location("pinned_correlated_source", OLD_SOURCE)
    require(spec is not None and spec.loader is not None, "old source import specification")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def gaussian_integer(a: np.ndarray) -> np.ndarray:
    require(np.array_equal(a.real, np.rint(a.real)), "Gaussian integer real part")
    require(np.array_equal(a.imag, np.rint(a.imag)), "Gaussian integer imaginary part")
    return np.rint(a.real).astype(np.int64) + 1j * np.rint(a.imag).astype(np.int64)


def exact_gaussian_matrix(a: np.ndarray) -> sp.Matrix:
    return sp.Matrix([[sp.Integer(int(round(z.real))) + sp.I * sp.Integer(int(round(z.imag))) for z in row]
                      for row in a])


def mu4_exponent(z: complex) -> int:
    values = (1, 1j, -1, -1j)
    hits = [k for k, value in enumerate(values) if abs(z - value) < 1e-9]
    require(len(hits) == 1, "phase is in mu4")
    return hits[0]


def projective_phase_code(values: np.ndarray) -> tuple[int, ...]:
    exponents = [mu4_exponent(z) for z in values]
    origin = exponents[0]
    return tuple((x - origin) % 4 for x in exponents)


def multiplicities_from_moments(matrix: np.ndarray, roots: tuple[int, ...]) -> list[int]:
    moments = [int(np.trace(np.linalg.matrix_power(matrix, k))) for k in range(len(roots))]
    vandermonde = sp.Matrix([[r**k for r in roots] for k in range(len(roots))])
    return [int(x) for x in vandermonde.inv() * sp.Matrix(moments)]


def duad_chart(mark: int, mapping: dict[int, tuple[int, ...]]) -> tuple[np.ndarray, dict[int, int]]:
    others = [i for i in range(6) if i != mark]
    rename = {old: new for new, old in enumerate(others)}
    P = np.zeros((10, 10), dtype=np.int64)
    for bell_label, part in mapping.items():
        side = set(part)
        if mark not in side:
            side = set(range(6)) - side
        duad = tuple(sorted(rename[x] for x in side - {mark}))
        P[bell_label, PAIR_INDEX[duad]] = 1
    require(np.array_equal(P.T @ P, np.eye(10, dtype=np.int64)), f"mark {mark}: duad chart permutation")
    return P, rename


def old_mark_comparison(old, mark, bell, six, mapping, Cs, Snew):
    idx, counts, adjacency_bell, stay, jump, _c20, *_ = old.marked_data(mark, bell, six, mapping)
    P, rename = duad_chart(mark, mapping)
    Q = np.kron(P, P)
    J = gaussian_integer(Q.T @ jump @ Q)
    require(not np.any(J.imag), f"mark {mark}: old jump is real in duad chart")
    J = J.real.astype(np.int64)
    require(np.array_equal(J, J.T), f"mark {mark}: old jump Hilbert-Schmidt self-adjoint")
    vec_identity = np.eye(10, dtype=np.int64).reshape(-1, order="F")
    require(np.array_equal(J @ vec_identity, 12 * vec_identity), f"mark {mark}: old normalized jump unital")

    adjacency = np.zeros((10, 10), dtype=np.int64)
    for a, e in enumerate(PAIRS):
        for b, f in enumerate(PAIRS):
            adjacency[b, a] = int(a != b and set(e).isdisjoint(f))
    require(np.array_equal(P.T @ adjacency_bell @ P, adjacency), f"mark {mark}: Petersen chart")

    # New channel is Snew/3; normalized old jump channel is J/12.
    for a in range(10):
        col = a + 10 * a
        require(np.array_equal(J[:, col], 4 * Snew[:, col]), f"mark {mark}: equal populations column {a}")
    require(not np.array_equal(J, 4 * Snew), f"mark {mark}: full channels differ")
    require(int(sp.Matrix(J).rank()) == 40, f"mark {mark}: old normalized jump superoperator rank 40")
    old_roots = (12, 4, -8, -4, 0)
    poly = np.eye(100, dtype=np.int64)
    for root in old_roots:
        poly = poly @ (J - root * np.eye(100, dtype=np.int64))
    require(not np.any(poly), f"mark {mark}: old jump spectral polynomial")
    old_mult = multiplicities_from_moments(J, old_roots)
    require(old_mult == [1, 20, 4, 15, 60], f"mark {mark}: old jump spectrum multiplicities")

    groups: dict[int, list[np.ndarray]] = {m: [] for m in range(5)}
    for k in idx:
        U = bell[k]
        perm = old.bell_permutation(U)
        D_jump = np.eye(10) - np.diag([int(perm[a] == a) for a in range(10)])
        K = gaussian_integer(P.T @ (U @ D_jump) @ P)
        target = six[k][mark]
        groups[rename[target]].append(K)

    all_kraus = []
    group_ranks = []
    phase_options: dict[int, list[np.ndarray]] = {}
    phase_edge_options: dict[int, list[list[tuple[int, int, int, int]]]] = {}
    phase_group_codes = []
    for m in range(5):
        require(len(groups[m]) == 4, f"mark {mark}: four old lifts in complement group {m}")
        effect = sum((K.conj().T @ K for K in groups[m]), start=np.zeros((10, 10), dtype=complex))
        require(np.array_equal(effect, 4 * (Cs[m].T @ Cs[m])),
                f"mark {mark}: grouped old effect equals four P_m")
        vectors = np.stack([K.reshape(-1, order="F") for K in groups[m]], axis=1)
        gram = gaussian_integer(vectors.conj().T @ vectors)
        require(np.array_equal(gram @ gram, 8 * gram), f"mark {mark}: grouped old Choi Gram polynomial {m}")
        rank = int(exact_gaussian_matrix(gram).rank())
        require(rank == 3, f"mark {mark}: grouped old Choi rank three {m}")
        group_ranks.append(rank)
        support = [i for i, pair in enumerate(PAIRS) if m not in pair]
        factor_data = []
        for K in groups[m]:
            D = K @ Cs[m]
            require(np.array_equal(D, np.diag(np.diag(D))),
                    f"mark {mark}: old Kraus is diagonal phase times C_{m}")
            require(all(abs(D[i, i]) == 1 for i in support) and
                    all(D[i, i] == 0 for i in range(10) if i not in support),
                    f"mark {mark}: phase factor has exact P_{m} support")
            require(np.array_equal(K, D @ Cs[m]), f"mark {mark}: exact phase factorization {m}")
            factor_data.append((projective_phase_code(np.diag(D)[support]), D))
        factor_data.sort(key=lambda item: item[0])
        phase_options[m] = [item[1] for item in factor_data]
        phase_edge_options[m] = []
        for D in phase_options[m]:
            edges = []
            for e in range(10):
                hits = np.flatnonzero(Cs[m][:, e])
                if len(hits) and e < int(hits[0]):
                    f = int(hits[0])
                    edges.append((e, f, m, mu4_exponent(D[f, f])))
            phase_edge_options[m].append(edges)
        D0 = factor_data[0][1]
        relative_codes = []
        relative_diagonals = []
        for _code, D in factor_data:
            relative = D @ D0.conj().T
            code = projective_phase_code(np.diag(relative)[support])
            relative_codes.append(code)
            relative_diagonals.append(relative)
        canonical_klein = [
            (0, 0, 0, 0, 0, 0),
            (0, 0, 2, 2, 0, 0),
            (0, 2, 0, 0, 2, 0),
            (0, 2, 2, 2, 2, 0),
        ]
        require(relative_codes == canonical_klein,
                f"mark {mark}: relative phase records form canonical Klein four on support {m}")
        require(all(tuple((a + b) % 4 for a, b in zip(x, y)) in relative_codes
                    for x in relative_codes for y in relative_codes),
                f"mark {mark}: relative phase codes close projectively {m}")
        characters = [tuple(code[j] for code in relative_codes) for j in range(6)]
        character_multiplicities = sorted(characters.count(c) for c in set(characters))
        require(character_multiplicities == [2, 2, 2],
                f"mark {mark}: Klein phase representation has three characters twice {m}")
        phase_group_codes.append(relative_codes)
        G = sum((np.kron(K.conj(), K) for K in groups[m]), start=np.zeros((100, 100), dtype=complex))
        N = np.kron(Cs[m], Cs[m])
        for a in range(10):
            col = a + 10 * a
            require(np.array_equal(G[:, col], 4 * N[:, col]),
                    f"mark {mark}: grouped population equals Hodge outcome {m},{a}")
        require(not np.array_equal(G, 4 * np.kron(Cs[m], Cs[m])),
                f"mark {mark}: grouped old operation differs from one Hodge Kraus {m}")
        phase_twirl = sum((np.kron((H @ D0).conj(), H @ D0) @ N
                           for H in relative_diagonals),
                          start=np.zeros((100, 100), dtype=complex))
        require(np.array_equal(G, phase_twirl),
                f"mark {mark}: old block is conditional Klein phase twirl after Hodge outcome {m}")
        for D in phase_options[m]:
            # Each individual three-transposition lift has no cycle flux on
            # its own support.  The obstruction below only appears when all
            # five complement classes must share one ten-label gauge.
            for e in support:
                f = int(np.flatnonzero(Cs[m][:, e])[0])
                require(D[e, e] * D[f, f] == 1,
                        f"mark {mark}: individual phase lift has trivial two-cycle flux {m},{e}")
        all_kraus.extend(groups[m])

    vectors = np.stack([K.reshape(-1, order="F") for K in all_kraus], axis=1)
    gram = gaussian_integer(vectors.conj().T @ vectors)
    require(np.array_equal(gram @ gram, 8 * gram), f"mark {mark}: old Choi Gram polynomial")
    require(int(exact_gaussian_matrix(gram).rank()) == 15, f"mark {mark}: old Choi rank fifteen")

    # Test whether one diagonal gauge on all ten duad labels could remove a
    # choice of one baseline lift in every complement class.  Global phases
    # of Hermitian monomial lifts reduce to signs.  Exhaust all 4^5 choices
    # and 2^5 class signs; propagation on the connected Petersen graph is an
    # exact Z4 linear consistency test.
    def common_gauge_exists(choice: tuple[int, ...], class_signs: tuple[int, ...]) -> bool:
        equations = [equation for m, option in enumerate(choice)
                     for equation in phase_edge_options[m][option]]
        labels: list[int | None] = [None] * 10
        labels[0] = 0
        changed = True
        while changed:
            changed = False
            for e, f, m, phase in equations:
                sign = class_signs[m]
                if labels[e] is not None and labels[f] is None:
                    labels[f] = (phase - sign + labels[e]) % 4
                    changed = True
                elif labels[f] is not None and labels[e] is None:
                    labels[e] = (sign + labels[f] - phase) % 4
                    changed = True
        return all(labels[e] is not None and
                   (class_signs[m] + labels[f] - labels[e] - phase) % 4 == 0
                   for e, f, m, phase in equations)

    gauge_solutions = 0
    for choice in it.product(range(4), repeat=5):
        for signs in it.product((0, 2), repeat=5):
            gauge_solutions += int(common_gauge_exists(choice, signs))
    require(gauge_solutions == 0,
            f"mark {mark}: no common diagonal gauge removes the five baseline phase lifts")

    # A concrete coherence witness: the new map preserves two terms for every
    # distinct intersecting pair; the old jump map kills all of them.
    witness_count = 0
    for a, e in enumerate(PAIRS):
        for b, f in enumerate(PAIRS):
            if a != b and len(set(e) & set(f)) == 1:
                col = a + 10 * b
                require(np.count_nonzero(Snew[:, col]) == 2, f"mark {mark}: new intersecting coherence has two terms")
                require(not np.any(J[:, col]), f"mark {mark}: old intersecting coherence vanishes")
                witness_count += 1
    require(witness_count == 60, f"mark {mark}: sixty ordered intersecting coherence witnesses")
    return {
        "mark": mark,
        "population_equal": True,
        "old_superoperator_rank": 40,
        "old_trace_preserving": True,
        "old_unital": True,
        "old_hilbert_schmidt_self_adjoint": True,
        "old_superoperator_spectrum_scaled_by_12": {"12": 1, "4": 20, "-8": 4, "-4": 15, "0": 60},
        "old_choi_rank": 15,
        "old_nonzero_choi_eigenvalue": "2/3 (multiplicity 15)",
        "group_sizes": [4] * 5,
        "group_choi_ranks": group_ranks,
        "group_effects": "P_m/3 after normalization",
        "phase_factorization": "R_l=D_l C_m on P_m",
        "relative_phase_group": "Klein four projectively; three characters each with multiplicity two",
        "conditional_channel_factorization": "old block = outcome-conditioned Klein phase twirl after one Hodge outcome",
        "individual_baseline_two_cycle_flux": "trivial and removable with an m-conditioned diagonal gauge",
        "common_ten_label_gauge_solutions": gauge_solutions,
        "common_edge_flux": "nontrivial: none of 4^5 baseline choices and 2^5 class signs admits one common diagonal gauge",
        "full_channel_equal": False,
        "intersecting_coherence_witnesses": witness_count,
    }


def fock_annihilators() -> list[np.ndarray]:
    """Jordan-Wigner signed annihilators on the 32 subset basis."""
    out = []
    for m in range(5):
        a = np.zeros((32, 32), dtype=np.int64)
        for mask in range(32):
            if mask & (1 << m):
                target = mask ^ (1 << m)
                sign = -1 if (mask & ((1 << m) - 1)).bit_count() % 2 else 1
                a[target, mask] = sign
        out.append(a)
    return out


def verify_full_clifford_path(Cs: list[np.ndarray]):
    annihilation = fock_annihilators()
    identity32 = np.eye(32, dtype=np.int64)
    total = sum((a.T @ a + a @ a.T for a in annihilation), start=np.zeros((32, 32), dtype=np.int64))
    require(np.array_equal(total, 5 * identity32), "ten Clifford Kraus effects sum to 5I on Fock space")

    # Keep the two-particle basis in the same lexicographic pair order used by
    # the Hodge matrices.  Integer mask order is a different ordering.
    grade2 = [sum(1 << i for i in pair) for pair in PAIRS]
    grade3 = [mask for mask in range(32) if mask.bit_count() == 3]
    grade1 = [mask for mask in range(32) if mask.bit_count() == 1]
    E2 = np.eye(32, dtype=np.int64)[:, grade2]
    E3 = np.eye(32, dtype=np.int64)[:, grade3]
    E1 = np.eye(32, dtype=np.int64)[:, grade1]
    create_effect = sum((E2.T @ a @ a.T @ E2 for a in annihilation), start=np.zeros((10, 10), dtype=np.int64))
    annihilate_effect = sum((E2.T @ a.T @ a @ E2 for a in annihilation), start=np.zeros((10, 10), dtype=np.int64))
    require(np.array_equal(create_effect, 3 * np.eye(10, dtype=np.int64)), "creation branch effect 3I")
    require(np.array_equal(annihilate_effect, 2 * np.eye(10, dtype=np.int64)), "annihilation branch effect 2I")

    # Complement triples to pairs.  This is the same epsilon transport as C_m,
    # up to one harmless fixed basis-sign convention, checked at CP-map level.
    D = np.zeros((10, 10), dtype=np.int64)
    for col, mask in enumerate(grade3):
        triple = tuple(i for i in range(5) if mask & (1 << i))
        pair = tuple(i for i in range(5) if not (mask & (1 << i)))
        D[PAIR_INDEX[pair], col] = parity(triple + pair)
    require(np.array_equal(D.T @ D, np.eye(10, dtype=np.int64)), "triple-complement isometry")
    creation_kraus = [D @ E3.T @ a.T @ E2 for a in annihilation]
    creation_super = sum((np.kron(K, K) for K in creation_kraus), start=np.zeros((100, 100), dtype=np.int64))
    hodge_super = sum((np.kron(C, C) for C in Cs), start=np.zeros((100, 100), dtype=np.int64))
    require(np.array_equal(creation_super, hodge_super), "conditional creation channel equals Hodge complement channel")

    # The ten Kraus population walk toggles one bit uniformly: Q_5/5.
    population = np.zeros((32, 32), dtype=np.int64)
    for a in annihilation:
        for K in (a, a.T):
            population += K * K
    expected = np.zeros((32, 32), dtype=np.int64)
    for mask in range(32):
        for m in range(5):
            expected[mask ^ (1 << m), mask] = 1
    require(np.array_equal(population, expected), "ten Clifford Kraus populations are the five-cube adjacency")

    # Quotient the five-cube by antipodes.  Since complement reverses parity
    # in odd dimension, each antipodal class has one even representative.
    even = [mask for mask in range(32) if mask.bit_count() % 2 == 0]
    even_index = {mask: i for i, mask in enumerate(even)}
    quotient = np.zeros((16, 16), dtype=np.int64)
    for mask in even:
        for m in range(5):
            odd_neighbor = mask ^ (1 << m)
            representative = odd_neighbor ^ 31
            quotient[even_index[representative], even_index[mask]] = 1
    require(np.all(quotient.sum(axis=0) == 5), "antipodal five-cube quotient is five-regular")
    require(np.array_equal(quotient @ quotient,
                           5 * np.eye(16, dtype=np.int64) + 2 * (np.ones((16, 16), dtype=np.int64)
                                                                  - np.eye(16, dtype=np.int64) - quotient)),
            "antipodal five-cube quotient is Clebsch SRG(16,5,0,2)")
    return {
        "full_fock_channel": "(1/5) sum_m [Ad(a_m)+Ad(a_m^dagger)]",
        "population": "Q5/5 on 32 subset states",
        "grade2_creation_probability": "3/5",
        "grade2_annihilation_probability": "2/5",
        "conditioned_creation_after_complement": "five Hodge channel",
        "annihilation_output_dimension": 5,
        "old_stay_identification": False,
        "antipodal_quotient": "Clebsch SRG(16,5,0,2), the already checked v810 object",
    }


def main(output: Path | None) -> dict:
    checks.clear()
    for path, digest in PINS.items():
        require(hashlib.sha256(Path(path).read_bytes()).hexdigest() == digest, "source pin " + path)

    Cs = hodge_amplitudes()
    effects = []
    for m, C in enumerate(Cs):
        Pm = np.diag([int(m not in pair) for pair in PAIRS]).astype(np.int64)
        require(np.array_equal(C, C.T), f"C_{m} symmetric")
        require(int(sp.Matrix(C).rank()) == 6, f"C_{m} rank six")
        require(np.array_equal(C @ C, Pm), f"C_{m} squared is pair-exclusion projector")
        effects.append(C.T @ C)
    require(np.array_equal(sum(effects), 3 * np.eye(10, dtype=np.int64)), "sum C_m dagger C_m is 3I")
    require(np.array_equal(sum(C @ C.T for C in Cs), 3 * np.eye(10, dtype=np.int64)), "five-Hodge channel unital")
    gram = np.array([[np.trace(A.T @ B) for B in Cs] for A in Cs], dtype=np.int64)
    require(np.array_equal(gram, 6 * np.eye(5, dtype=np.int64)), "five Hodge Kraus vectors orthogonal")

    Snew = sum((np.kron(C, C) for C in Cs), start=np.zeros((100, 100), dtype=np.int64))
    require(np.array_equal(Snew, Snew.T), "five-Hodge superoperator Hilbert-Schmidt self-adjoint")
    require(int(sp.Matrix(Snew).rank()) == 100, "five-Hodge superoperator invertible")
    roots = (3, 2, 1, -1, -2)
    poly = np.eye(100, dtype=np.int64)
    for root in roots:
        poly = poly @ (Snew - root * np.eye(100, dtype=np.int64))
    require(not np.any(poly), "five-Hodge spectral polynomial")
    multiplicities = multiplicities_from_moments(Snew, roots)
    require(multiplicities == [1, 10, 40, 35, 14], "five-Hodge spectrum multiplicities")
    sl5_generators = verify_sl5_covariance(Cs)
    su5_sectors = verify_su5_sector_spectrum(Cs, Snew)

    old = load_old_source()
    rays = old.source_rays()
    bell, raw_six, _ = old.reflection_actions(rays)
    six, _ = old.outer_mark_action(raw_six)
    mapping = old.find_bell_partition_intertwiner(bell, six)
    marks = [old_mark_comparison(old, mark, bell, six, mapping, Cs, Snew) for mark in range(6)]
    fock = verify_full_clifford_path(Cs)

    result = {
        "verdict": [
            "FIVE_SEPARATED_HODGE_AMPLITUDES_DEFINE_A_CANONICAL_CPTP_CHANNEL",
            "POPULATION_IS_EXACTLY_PETERSEN_ADJACENCY_OVER_THREE",
            "OLD_NORMALIZED_J20_HAS_THE_SAME_FIVE_COARSE_EFFECTS_AND_POPULATIONS",
            "FULL_QUANTUM_CHANNELS_AND_OUTCOME_OPERATIONS_ARE_NOT_EQUAL",
            "OLD_J20_BLOCKS_ARE_HODGE_OUTCOMES_FOLLOWED_BY_CONDITIONAL_KLEIN4_PHASE_TWIRLS",
            "PHI_SQUARED_HAS_THE_EXACT_SU5_PROJECTOR_SPECTRUM_1_PLUS_24_PLUS_75",
            "PETERSEN_POPULATIONS_RESTRICT_TO_THE_ACTUAL_RADIAL_SHELLS_1_3_6",
            "CHANNEL_IS_THE_NORMALIZED_ONE_PARTICLE_CREATION_OR_DUAL_FERMION_LOSS_CHANNEL",
            "PHYSICAL_TFPT_SELECTION_OR_RECORD_OF_THE_COMPLEMENT_INDEX_IS_NOT_DERIVED",
        ],
        "source_sha256": PINS,
        "hodge_channel": {
            "space": "Lambda^2 C^5",
            "dimension": 10,
            "kraus_operators": "K_m=C_m/sqrt(3), m=0,...,4",
            "kraus_ranks": [6] * 5,
            "effect_sum": "sum C_m^dagger C_m=3 I_10",
            "C_m_square": "projector onto pairs not containing m",
            "trace_preserving": True,
            "unital": True,
            "hilbert_schmidt_self_adjoint": True,
            "choi_rank": 5,
            "choi_spectrum": {"2": 5, "0": 95},
            "superoperator_rank": 100,
            "superoperator_spectrum": {"1": 1, "2/3": 10, "1/3": 40, "-1/3": 35, "-2/3": 14},
            "population": "A_Petersen/3",
            "sl5_infinitesimal_generators_checked": sl5_generators,
            "covariance": "Phi(conj(Lambda2 U) rho conj(Lambda2 U)^dagger)=Lambda2 U Phi(rho) Lambda2 U^dagger for U in SU(5)",
            "su5_sectors": su5_sectors,
        },
        "old_J20_comparisons": marks,
        "full_clifford_fock_comparison": fock,
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
    payload = main(args.out)
    print(json.dumps(payload, indent=2, sort_keys=True))
