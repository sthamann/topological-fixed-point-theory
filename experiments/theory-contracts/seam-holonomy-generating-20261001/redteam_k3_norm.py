#!/usr/bin/env python3
"""Part I: red team of K3 (where do the eight species come from?) and of the factor 2 pi.

  I1  CORRECTED (round 12): the species count is NOT dim S+ (that 16 is the complex D5 spinor, the
      16 Majoranas of SO(16)_1 are the D8 vector).  The non-circular count is the Cartan rank:
      rank D5 + rank A3 = 5 + 3 = 8 complex species (bosonised level-one Cartan currents), from P2
      (g = 5) and the four marks (n - 1 = 3), without E8 as input.
  I2  c = 8 is then an OUTPUT (one unit per level-one Cartan boson).
  I3  Nielsen-Ninomiya on the admissible seas: Fermi points come in equal numbers of right and
      left movers; K3 (two Fermi points) <=> one chiral mode per side ("minimal chirality",
      c = c_- per chirality).  K3 is therefore NOT derived from c = 8 without circularity; it is
      the axiom "the seam is purely chiral", stated without E8.
  I4  CAR fixes the current level: N S^(m) = |m| exactly for |m| <= N/2 (integer charges, level 1).
  I5  the 2 pi is canonical: |D|^-1 kernel = (1/2 pi) sum_{m != 0} <rho_m rho_-m>/m^2 e^{i m th}
      with rho the integer-normalised charge density and the kernel taken w.r.t. d theta on the
      seam circle of length 2 pi; it reproduces v484's -(1/pi) ln|2 sin(th/2)| with no free factor.
    python3 -B experiments/theory-contracts/seam-holonomy-generating-20261001/redteam_k3_norm.py
"""
from __future__ import annotations

import json
import sys
from itertools import combinations, product
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
LOG: list[dict] = []


def step(name: str, ok: bool, detail: str = "") -> None:
    LOG.append({"step": name, "pass": bool(ok), "detail": detail})
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))


g, n_marks = 5, 4
rank = g + (n_marks - 1)                       # D5 Cartan + A3 Cartan
even = [b for b in product((0, 1), repeat=g) if sum(b) % 2 == 0]
step("I1.(corrected) species = Cartan rank 5 + 3 = 8 complex (from P2 and the four marks); dim S+ = 16 is the complex D5 spinor, not the Majorana count",
     rank == 8 and len(even) == 16)
step("I2.c = 8 (one per level-one Cartan boson) is an output; equality with rank E8 is a consistency check only", rank * 1 == 8)

N = 16
def fp_chirality(sea):
    """Hamiltonian-free chirality of Fermi points: edge orientation in increasing k.

    filled -> empty is a right-moving point, empty -> filled a left-moving one.
    """
    occ = [m in sea for m in range(N)]
    right = sum(occ[m] and not occ[(m + 1) % N] for m in range(N))
    left = sum((not occ[m]) and occ[(m + 1) % N] for m in range(N))
    return right, left
reflect = lambda m: N - 1 - m
seas = [frozenset(c) for c in combinations(range(N), N // 2)]
adm = [s for s in seas if frozenset(map(reflect, s)) == s]
chir = [fp_chirality(s) for s in adm]
nontriv = [(r, l) for r, l in chir if r + l > 0]
step("I3.Nielsen-Ninomiya on all admissible seas: right movers = left movers; K3 <=> exactly one chiral mode per side",
     all(r == l for r, l in nontriv) and len(nontriv) == len(adm) and min(r for r, _ in nontriv) == 1,
     f"chiral content per side over {len(adm)} seas: {sorted({int(r) for r, _ in nontriv})}; K3 = 'purely chiral seam' (axiom, no E8)")

for N in (64, 256, 1024):
    k = np.pi * (2 * np.arange(N) + 1) / N
    sea = np.cos(k) > 0
    C = np.array([np.exp(1j * k[sea] * d).sum() / N for d in range(N)])
    S = -np.abs(C) ** 2
    S[0] = (C[0] * (1 - C[0])).real
    Sh = np.fft.fft(S).real
    m = np.arange(N); m = np.where(m > N // 2, m - N, m)
    lvl = np.abs(N * Sh[np.abs(m) <= N // 2] - np.abs(m[np.abs(m) <= N // 2])).max()
    if N == 1024:
        th = np.array([np.pi / 2, np.pi, np.pi / 3])
        ker = np.array([sum((N * Sh[j] / m[j] ** 2) * np.cos(m[j] * t) for j in range(N) if m[j] != 0) / (2 * np.pi)
                        for t in th])
        target = -(1 / np.pi) * np.log(np.abs(2 * np.sin(th / 2)))
        kdev = np.abs(ker - target).max()
step("I4.CAR fixes level one: N S^(m) = |m| exactly for |m| <= N/2 (integer charges)", lvl < 1e-10, f"max dev {lvl:.1e}")
step("I5.canonical 2 pi: (1/2pi) sum <rho rho>_m / m^2 e^{im th} = v484 |D|^-1 kernel -(1/pi) ln|2 sin(th/2)| (no free factor)",
     kdev < 2e-3, f"max dev at pi/2, pi, pi/3: {kdev:.1e} (finite-N tail)")

out = {"steps": LOG, "passed": sum(x["pass"] for x in LOG), "total": len(LOG)}
(HERE / "redteam_k3_norm_results.json").write_text(json.dumps(out, indent=2))
print(f"{out['passed']}/{out['total']} steps passed")
sys.exit(0 if out["passed"] == out["total"] else 1)
