#!/usr/bin/env python3
"""Exact whole-96 type check for the Hardy Bmu weak-Clifford odd lift."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
OUT = HERE / "hardy-bmu-weak-clifford-odd-lift-results.json"
sys.path.insert(0, str(ROOT / "verification"))
import v252_full_finite_triple as v252  # noqa: E402


def require(ok: object, label: str) -> None:
    if not bool(ok):
        raise RuntimeError(label)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def creators() -> list[np.ndarray]:
    out = []
    for k in range(5):
        c = np.zeros((32, 32), dtype=int)
        for mask in range(32):
            if not ((mask >> k) & 1):
                sign = (-1) ** ((mask & ((1 << k) - 1)).bit_count())
                c[mask | (1 << k), mask] = sign
        out.append(c)
    return out


def compact_conjugation() -> np.ndarray:
    compact = np.zeros((32, 32), dtype=int)
    for mask in range(32):
        complement = 31 ^ mask
        source = [i for i in range(5) if (mask >> i) & 1]
        target = [i for i in range(5) if (complement >> i) & 1]
        inversions = sum(i > j for i in source for j in target)
        sign = (-1) ** (inversions + len(target) * (len(target) - 1) // 2)
        compact[complement, mask] = sign
    return compact


def whole96_dictionary() -> np.ndarray:
    compact = compact_conjugation()
    w = np.zeros((v252.NH, v252.NH), dtype=int)
    for i, (name, chirality, color, generation) in enumerate(v252.P):
        if chirality == "L":
            if name == "nu":
                mask, sign = 15, 1
            elif name == "e":
                mask, sign = 23, 1
            else:
                mask, sign = (1 << color) | (8 if name == "u" else 16), 1
            column = i
        else:
            if name == "nu":
                mask, sign = 0, 1
            elif name == "e":
                mask, sign = 24, 1
            else:
                mask, sign = 7 ^ (1 << color), (-1) ** color
                if name == "d":
                    mask |= 24
            column = v252.NP + i
        w[32 * generation + mask, column] = sign
        paired = (column + v252.NP) % v252.NH
        w[32 * generation:32 * (generation + 1), paired] = compact[:, mask] * sign
    require(np.array_equal(w.T @ w, np.eye(96, dtype=int)), "whole96 signed unitary")
    return w


def family_02() -> np.ndarray:
    """G-column row 0, column 2, with all actual finite labels retained."""
    e02 = np.zeros((96, 96), dtype=int)
    for sector in (0, 1):
        offset = sector * v252.NP
        for local, (name, chirality, color, generation) in enumerate(v252.P):
            if generation == 2:
                source = offset + local
                target = offset + v252.IDX[(name, chirality, color, 0)]
                e02[target, source] = 1
    return e02


def label(index: int) -> str:
    sector = "particle" if index < v252.NP else "antiparticle"
    name, chirality, color, generation = v252.P[index % v252.NP]
    return f"{sector}:{name}_{chirality}:c{color}:G{generation}"


def main() -> None:
    c = creators()
    w = whole96_dictionary()
    gamma = v252.grading().astype(int)
    pplus = (np.eye(96, dtype=int) + gamma) // 2
    e02 = family_02()
    require(np.array_equal(e02 @ gamma, gamma @ e02), "family coefficient is neutral/even")

    parity32 = np.diag([(-1) ** mask.bit_count() for mask in range(32)])
    gammas = []
    for k in range(5):
        gammas.extend((c[k] + c[k].T, 1j * (c[k] - c[k].T)))
    d5 = [(gammas[i] @ gammas[j] - gammas[j] @ gammas[i]) / 4
          for i in range(10) for j in range(i + 1, 10)]
    require(len(d5) == 45, "forty-five D5 bilinears")
    require(all(np.array_equal(m @ parity32, parity32 @ m) for m in d5),
            "all D5 bilinears commute with parity")

    six_y_ext = np.diag([
        sum((-2, -2, -2, 3, 3)[i] for i in range(5) if (mask >> i) & 1)
        for mask in range(32)
    ] * 3)
    two_t3_ext = np.diag([
        ((mask >> 3) & 1) - ((mask >> 4) & 1) for mask in range(32)
    ] * 3)
    color_ext = [
        np.diag([((mask >> i) & 1) - ((mask >> j) & 1) for mask in range(32)] * 3)
        for i, j in ((0, 1), (1, 2))
    ]
    six_y = w.T @ six_y_ext @ w
    two_t3 = w.T @ two_t3_ext @ w
    color = [w.T @ h @ w for h in color_ext]

    results = {}
    for k, weak_name, t3_weight in ((3, "w1", 1), (4, "w2", -1)):
        x = w.T @ np.kron(np.eye(3, dtype=int), c[k]) @ w
        block = x @ (e02 @ pplus) - (e02 @ pplus) @ x
        expected = e02 @ x @ gamma
        require(np.array_equal(block, expected), f"{weak_name} commutator identity")
        require(np.count_nonzero(block) == 16, f"{weak_name} sixteen entries")
        require(np.linalg.matrix_rank(block) == 16, f"{weak_name} rank sixteen")
        require(np.array_equal(gamma @ block + block @ gamma, np.zeros((96, 96), dtype=int)),
                f"{weak_name} gamma odd")
        require(np.array_equal(six_y @ block - block @ six_y, 3 * block),
                f"{weak_name} hypercharge plus one half")
        require(np.array_equal(two_t3 @ block - block @ two_t3, t3_weight * block),
                f"{weak_name} weak weight")
        require(all(np.array_equal(h @ block, block @ h) for h in color),
                f"{weak_name} color neutral")
        mappings = [
            {"input": label(int(j)), "output": label(int(i)), "sign": int(block[i, j])}
            for i, j in zip(*np.nonzero(block))
        ]
        results[weak_name] = {
            "rank": 16,
            "nonzero_entries": 16,
            "six_times_hypercharge": 3,
            "two_times_T3": t3_weight,
            "gamma_F_odd": True,
            "color_neutral": True,
            "mappings": mappings,
        }

    zero = np.zeros((3, 3), complex)
    d_zero = v252.full_D(zero, zero, zero, zero, zero)
    require(np.count_nonzero(d_zero) == 0, "zero Yukawa gives zero finite Dirac")

    source_files = [
        ROOT / "verification/v252_full_finite_triple.py",
        ROOT / "verification/v474_entropic_hodge_carrier.py",
        HERE / "HARDY-BMU-FIRST-NONZERO-ENTRY.txt",
        HERE / "check_actual_chiral_dictionary.py",
    ]
    output = {
        "scope": "Exact actual-v252/P2 operator type; no source execution or Higgs identification",
        "hardy_coefficient": "kappa E_02 P_gamma,+ with kappa=(4 X_*/pi) I != 0",
        "d5_zero_mode_result": "all 45 bilinears commute with P_gamma,+",
        "weak_clifford_identity": "[X_wa,E_02 P_gamma,+]=E_02 X_wa gamma_F",
        "zero_yukawa_one_forms": "zero because D_F=0, while both weak blocks retain rank 16",
        "weak_blocks": results,
        "source_sha256": {str(path.relative_to(ROOT)): sha256(path) for path in source_files},
    }
    OUT.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "d5_bilinears": 45,
        "weak_ranks": {key: value["rank"] for key, value in results.items()},
        "checks": "PASS",
    }, sort_keys=True))


if __name__ == "__main__":
    main()
