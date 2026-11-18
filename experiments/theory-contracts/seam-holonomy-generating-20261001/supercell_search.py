#!/usr/bin/env python3
"""Part SC: BCC.SUPERCELL.SEARCH.01 -- periodic bcc supercells that contain the light node Gamma AND the mirror P.

Superlattices L = n B (n integer 3x3 in Hermite normal form, index N = det n = number of bcc sites) of the primitive bcc
basis b1 = (1,1,-1), b2 = (1,-1,1), b3 = (-1,1,1).  P = (pi/2)(1,1,1) is an allowed momentum iff P.R in 2 pi Z for all
superlattice vectors, i.e. every row sum of n is = 0 mod 4.
  SC1 the index-4 debug cell R1 = b1 - b2, R2 = b2 - b3, R3 = b1 + b2 + 2 b3 contains Gamma and P; its four momenta are
      exactly {Gamma, P, -P, H} (all nodes of the Weyl block) -- the Hamiltonian walk vanishes identically on it, so it
      tests the interaction alone (strong geometric folding).
  SC2 census for N <= 32: the number of P-compatible cells per index, the best ones by residual point group (subgroup
      of O_h, 48 elements), shortest vector and anisotropy (longest / shortest reduced basis vector); the cubic box of
      side 4 (N = 16) and side 8 (N = 128) appear with the full O_h.
  SC3 (statement) full-O_h cells at N = 4, 16, 32 (and 128); 4 and 8 only for debugging, 16 and 32 for the phase study.
    python3 -B experiments/theory-contracts/seam-holonomy-generating-20261001/supercell_search.py
"""
from __future__ import annotations

import json
import sys
from itertools import permutations, product
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
LOG: list[dict] = []


def step(name: str, ok: bool, detail: str = "") -> None:
    LOG.append({"step": name, "pass": bool(ok), "detail": detail})
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))
    sys.stdout.flush()


B = np.array([[1, 1, -1], [1, -1, 1], [-1, 1, 1]])                  # rows b1, b2, b3
P = np.pi / 2 * np.ones(3)
fcc = lambda k: np.allclose(np.cos(np.array([k @ np.array(s) for s in product((-1, 1), repeat=3)])), 1, atol=1e-9)
PG = []
for p in permutations(range(3)):
    for sg in product((-1, 1), repeat=3):
        R = np.zeros((3, 3), int)
        for i in range(3):
            R[i, p[i]] = sg[i]
        PG.append(R)


def momenta(Rv):
    """allowed momenta of the supercell modulo the fcc reciprocal lattice."""
    G = 2 * np.pi * np.linalg.inv(Rv).T                              # rows: reciprocal superlattice basis (cubic coords)
    N = int(round(abs(np.linalg.det(Rv)) / 4))
    ks = []
    for m in product(range(-N, N + 1), repeat=3):
        k = np.array(m) @ G
        if not any(fcc(k - q) for q in ks):
            ks.append(k)
        if len(ks) == N:
            break
    return ks


def residual_group(Rv):
    inv = np.linalg.inv(Rv)
    cnt = 0
    for g in PG:
        c = (g @ Rv.T).T @ inv                                       # images of the rows in the R basis
        if np.allclose(c, np.round(c), atol=1e-9):
            cnt += 1
    return cnt


def reduced_lengths(Rv):
    vecs = [np.array(m) @ Rv for m in product(range(-2, 3), repeat=3) if any(m)]
    L = sorted(set(round(float(np.linalg.norm(v)), 9) for v in vecs))
    return L[0], L


# SC1 ----------------------------------------------------------------------------------------------------
R_dbg = np.array([B[0] - B[1], B[1] - B[2], B[0] + B[1] + 2 * B[2]])
N_dbg = int(round(abs(np.linalg.det(R_dbg)) / 4))
ks = momenta(R_dbg)
names = {"Gamma": np.zeros(3), "P": P, "-P": -P, "H": np.array([np.pi, 0, 0])}
found = {n: any(fcc(k - q) for k in ks) for n, q in names.items()}
sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.diag([1.0 + 0j, -1]); I2 = np.eye(2)
Wk = lambda k: (np.cos(k[0]) * I2 - 1j * np.sin(k[0]) * sx) @ (np.cos(k[1]) * I2 - 1j * np.sin(k[1]) * sy) @ (np.cos(k[2]) * I2 - 1j * np.sin(k[2]) * sz)
Hk = lambda k: 0.5j * (Wk(k) - Wk(k).conj().T)
walk_zero = all(np.allclose(Hk(k), 0) for k in ks)
step("SC1.index-4 debug cell R1 = b1 - b2, R2 = b2 - b3, R3 = b1 + b2 + 2 b3 has 4 sites and contains Gamma and P; its momenta are exactly {Gamma, P, -P, H}; the Hamiltonian walk vanishes identically on it (interaction-only test, strong folding)",
     N_dbg == 4 and len(ks) == 4 and all(found.values()) and walk_zero, f"momenta found {found}")

# SC2 ----------------------------------------------------------------------------------------------------
census = {}
best = {}
for N in range(1, 33):
    cells = []
    for a in range(1, N + 1):
        if N % a:
            continue
        for d in range(1, N // a + 1):
            if (N // a) % d:
                continue
            f_ = N // (a * d)
            for b_ in range(a):                                    # row-HNF: n[1][0], n[2][0] mod a, n[2][1] mod d
                for c in range(a):
                    for e in range(d):
                        n = np.array([[a, 0, 0], [b_, d, 0], [c, e, f_]])   # lower-triangular HNF (rows = superlattice vectors)
                        if all(int(r.sum()) % 4 == 0 for r in n):
                            cells.append(n)
    census[N] = len(cells)
    if cells:
        scored = []
        for n in cells:
            Rv = n @ B
            grp = residual_group(Rv)
            lmin, Ls = reduced_lengths(Rv)
            vecs = sorted(np.linalg.norm(np.array(m) @ Rv) for m in product(range(-1, 2), repeat=3) if any(m))
            aniso = vecs[2] / vecs[0] if len(vecs) > 2 else np.inf
            scored.append((-grp, -lmin, aniso, n.tolist()))
        scored.sort()
        best[N] = scored[0]
cubic16 = np.array([[4, 0, 0], [0, 4, 0], [0, 0, 4]])                  # cubic box side 4 in cubic coordinates
grp16 = residual_group(cubic16.astype(float))
n16 = int(round(abs(np.linalg.det(cubic16)) / 4))
cubic128 = 2 * cubic16
grp128 = residual_group(cubic128.astype(float))
step("SC2.census N <= 32: P-compatible cells exist only for N = 0 mod 4; best cell per index by residual point group, shortest vector and anisotropy; the cubic boxes of side 4 (16 sites) and 8 (128 sites) keep the full O_h (48)",
     all(census[N] == 0 for N in census if N % 4) and all(census[N] > 0 for N in (4, 8, 12, 16, 20, 24, 28, 32)) and grp16 == 48 and n16 == 16 and grp128 == 48,
     "; ".join(f"N={N}: {census[N]} cells, best |G|={-best[N][0]}, l_min={-best[N][1]:.2f}" for N in (4, 8, 12, 16, 24, 32)))
rec = [N for N in (16, 24, 32) if -best[N][0] >= 8]
full_Oh = [N for N in best if -best[N][0] == 48]
step("SC3.(statement) full-O_h P-compatible cells exist for N = 4, 16, 32 (and the cubic 128); recommended: 4 and 8 for debugging only (their folding makes the phase vanish, Part RW), 16 and 32 for the phase study, then 128",
     full_Oh == [4, 16, 32], f"full-O_h indices up to 32: {full_Oh}; best cells: { {N: best[N][3] for N in (8, 16, 24, 32)} }")

out = {"steps": LOG, "passed": sum(x["pass"] for x in LOG), "total": len(LOG), "census": census,
       "best": {str(N): {"group": -v[0], "l_min": -v[1], "anisotropy": v[2], "n": v[3]} for N, v in best.items()}, "debug_cell": R_dbg.tolist()}
(HERE / "supercell_search_results.json").write_text(json.dumps(out, indent=2, default=float))
print(f"{out['passed']}/{out['total']} steps passed")
sys.exit(0 if out["passed"] == out["total"] else 1)
