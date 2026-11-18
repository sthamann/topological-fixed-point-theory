"""Replay the actual v252 signed dictionary; no original run()/ledger call."""
from pathlib import Path
import argparse
import json
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "verification"))
import v252_full_finite_triple as v252


def replay():
    dim = 32
    creators = []
    for i in range(5):
        creator = np.zeros((dim, dim), dtype=int)
        for mask in range(dim):
            if not (mask >> i) & 1:
                sign = (-1) ** ((mask & ((1 << i) - 1)).bit_count())
                creator[mask | (1 << i), mask] = sign
        creators.append(creator)
    compact = np.zeros((dim, dim), dtype=int)
    for mask in range(dim):
        complement = 31 ^ mask
        source = [i for i in range(5) if (mask >> i) & 1]
        target = [i for i in range(5) if (complement >> i) & 1]
        inversions = sum(i > j for i in source for j in target)
        sign = (-1) ** (inversions + len(target) * (len(target) - 1) // 2)
        compact[complement, mask] = sign

    W = np.zeros((v252.NH, v252.NH), dtype=int)
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
        W[32 * generation + mask, column] = sign
        paired = (column + v252.NP) % v252.NH
        W[32 * generation : 32 * (generation + 1), paired] = compact[:, mask] * sign

    identity = np.eye(dim, dtype=int)
    cartan2 = [2 * c @ c.T - identity for c in creators]
    roots = [creators[i] @ creators[j].T for i in range(5) for j in range(5) if i != j]
    upper = [creators[i] @ creators[j] for i in range(5) for j in range(i + 1, 5)]
    roots += upper + [root.T for root in upper]
    checks = {
        "W_real_signed_orthogonal": np.array_equal(W.T @ W, np.eye(v252.NH, dtype=int)),
        "whole96_grading": np.array_equal(W @ v252.grading() @ W.T,
                                         np.diag([(-1) ** m.bit_count() for m in range(32)] * 3)),
        "whole96_JF_compact": np.array_equal(W @ v252.Kswap() @ W.T,
                                            np.kron(np.eye(3, dtype=int), compact)),
        "compact_square_plus_one": np.array_equal(compact @ compact, identity),
        "five_Cartan_compact_flip": all(np.array_equal(compact @ h @ compact, -h) for h in cartan2),
        "forty_D5_roots_compact_conjugation": all(np.array_equal(compact @ e @ compact, -e.T) for e in roots),
        "forty_D5_roots_split_reality": all(np.array_equal(e.conj(), e) for e in roots),
    }
    basis = []
    for column in range(v252.NH):
        name, chirality, color, generation = v252.P[column % v252.NP]
        row, = np.flatnonzero(W[:, column])
        mask = int(row % dim)
        subset = [i for i in range(5) if (mask >> i) & 1]
        basis.append({
            "original_index": column,
            "sector": "particle" if column < v252.NP else "antiparticle",
            "field": name, "chirality": chirality, "color_zero_based": color,
            "generation_zero_based": generation,
            "G_line": ("L0", "L2", "L1")[generation],
            "internal_gamma": (-1) ** len(subset),
            "wedge_indices_one_based": [i + 1 for i in subset],
            "sign": int(W[row, column]),
            "six_times_hypercharge": sum((-2, -2, -2, 3, 3)[i] for i in subset),
        })
    assert all(checks.values()), checks
    return {
        "scope": "Exact actual-v252 fibre/SM/grade/real-structure comparison; no P1 source promotion",
        "basis_order": ["c1", "c2", "c3", "w1", "w2"],
        "original_functions_used": ["P", "grading", "Kswap"],
        "original_run_called": False,
        "checks": {name: bool(value) for name, value in checks.items()},
        "basis": basis,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = replay()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"checks": result["checks"], "basis_rows": len(result["basis"])}))
