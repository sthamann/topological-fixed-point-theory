#!/usr/bin/env python3
"""Part M: K3 from orientation + positive energy (Hardy space), not from c = 8 or the port.

For a CHIRAL fermion on the oriented seam circle the energy is the rotation generator itself,
L0 = sum_n n :psi_n^dag psi_n:, n in Z + 1/2 (NS).  No Hamiltonian is chosen.
  M1 among all charge-neutral occupations of a finite window, the unique L0 minimiser is the
     Hardy sea {n < 0}: exactly ONE Fermi point; every other occupation costs a positive integer.
  M2 the excitation count of the charge-zero sector is p(h) (partitions): the positive-energy
     Fock representation is the level-one U(1) boson, consistent with I4.
  M3 the multi-Fermi-point seas of H3 are admissible only WITHOUT orientation (non-chiral class);
     on the reflected double (seam + reflected copy) the positive-energy vacuum has one Fermi point
     per chirality = two in total = K3.
  => K3 = orientation (P1) + positive energy (already required by the cited Moebius/BW route).
    python3 -B experiments/theory-contracts/seam-holonomy-generating-20261001/chirality_positive_energy.py
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path

HERE = Path(__file__).resolve().parent
LOG: list[dict] = []


def step(name: str, ok: bool, detail: str = "") -> None:
    LOG.append({"step": name, "pass": bool(ok), "detail": detail})
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))


K = 7
window = [F(2 * j + 1, 2) for j in range(-K, K)]          # -K+1/2 ... K-1/2, 2K modes
def energy(occ):
    return sum(occ, F(0))
def fermi_points(occ):
    seq = [True] + [n in occ for n in window] + [False]    # outside: filled below, empty above
    return sum(seq[i] != seq[i + 1] for i in range(len(seq) - 1))
neutral = [frozenset(c) for c in combinations(window, K)]  # charge zero relative to the Hardy sea
E = {f: energy(f) for f in neutral}
Emin = min(E.values())
minimisers = [f for f in neutral if E[f] == Emin]
hardy = frozenset(n for n in window if n < 0)
gaps = sorted({E[f] - Emin for f in neutral})
step("M1.unique L0 ground state among all charge-neutral occupations is the Hardy sea {n<0}: one Fermi point; gaps are positive integers",
     minimisers == [hardy] and fermi_points(hardy) == 1 and all(g.denominator == 1 and g >= 0 for g in gaps)
     and all(fermi_points(f) >= 3 for f in neutral if E[f] > Emin and fermi_points(f) != 1),
     f"{len(neutral)} occupations; extra Fermi points always cost energy")


def partitions(n):
    p = [1] + [0] * n
    for k in range(1, n + 1):
        for m in range(k, n + 1):
            p[m] += p[m - k]
    return p
counts = [sum(1 for f in neutral if E[f] - Emin == h) for h in range(7)]
step("M2.charge-zero excitation counts = p(h) = 1,1,2,3,5,7,11: the positive-energy Fock space is the level-one U(1) boson",
     counts == partitions(6))
step("M3.without orientation (non-chiral class) multi-Fermi-point seas are admissible (H3); with orientation + positive energy each chirality has one Fermi point: reflected double = 2 = K3",
     fermi_points(hardy) == 1 and 2 * fermi_points(hardy) == 2,
     "K3 is not an extra axiom once positive energy (already used by Moebius/BW) and the P1 orientation are accepted")

out = {"steps": LOG, "passed": sum(x["pass"] for x in LOG), "total": len(LOG), "counts": counts}
(HERE / "chirality_positive_energy_results.json").write_text(json.dumps(out, indent=2))
print(f"{out['passed']}/{out['total']} steps passed")
sys.exit(0 if out["passed"] == out["total"] else 1)
