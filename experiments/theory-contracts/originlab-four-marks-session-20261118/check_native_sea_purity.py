#!/usr/bin/env python3
"""The purity of round twenty-five evaluated on the native tangential operator (18 November 2026).

Uses, unchanged, the assembler of the original Q contract
(original-paper-source-reconciliation-20261102/contribution_20261106Q/runtime/assemble.py) with
the domain native-diagnostic-v1.json: 96 fibres, four circle modes, the recorded quarter shift.
Four modes are four places on the circle, the four marks.  The state is the filled negative
spectrum of the assembled operator h_N.  This choice of state is an assumption of the hypothesis;
the assembler itself labels h_N as the tangential one-particle operator, not the joint generator.

Sections
    N1  the native operator: modes, quarter shift, spectrum
    N2  purity of the four places and of pairs of places; the 4 x 4 table of ratios
    N3  the same without any potential (free comparison) and for the baseline coefficients
"""
import sys
from pathlib import Path

import numpy as np

FAIL = []


def check(name, ok, info=None):
    print(("  ok   " if ok else "  FAIL ") + name + (("   " + str(info)) if info is not None else ""))
    if not ok:
        FAIL.append(name)


root = Path(sys.argv[1]).resolve()
base = root/"experiments/theory-contracts/original-paper-source-reconciliation-20261102/contribution_20261106Q"
sys.path.insert(0, str(base/"runtime"))
from assemble import Assembler  # noqa: E402

A = Assembler(base/"domain/native-diagnostic-v1.json")


def candidate(genes):
    return {"schema": "tfpt.field-configuration.v1", "domain_sha256": A.domain["_sha256"], "genes": genes}


def table(h, modes, alpha):
    """Ratios pur({a,b}) / (pur{a} pur{b}) for the four places, and the single purities."""
    m = len(modes)
    w, v = np.linalg.eigh(h)
    gap = float(np.min(np.abs(w)))
    P = v[:, w < 0] @ v[:, w < 0].conj().T
    theta = 2*np.pi*np.arange(m)/m
    F = np.array([[np.exp(1j*(n+alpha)*t) for n in modes] for t in theta])/np.sqrt(m)   # place x mode
    U = np.kron(F, np.eye(96))
    C = U @ P @ U.conj().T

    def logpur(places):
        idx = np.concatenate([np.arange(96*p, 96*(p+1)) for p in places])
        nu = np.clip(np.linalg.eigvalsh(C[np.ix_(idx, idx)]), 0, 1)
        return float(np.sum(np.log(nu**2+(1-nu)**2)))
    single = [logpur([a]) for a in range(m)]
    T = np.zeros((m, m))
    for a in range(m):
        for b in range(a+1, m):
            T[a, b] = T[b, a] = logpur([a, b])-single[a]-single[b]
    return {"negative": int(np.sum(w < 0)), "gap": gap, "single": single, "log_ratio": T,
            "whole": logpur(list(range(m)))}


# ---------------------------------------------------------------- N1
print("N1  the native operator")
full = A.assemble(candidate([2, 2, 2]))
h, modes, alpha = full["h_N"], full["modes"], full["alpha"]
print("      modes:", modes, " quarter shift:", alpha, " dimension:", h.shape[0])
check("four modes, 96 fibres, quarter shift 1/4", len(modes) == 4 and h.shape == (384, 384) and abs(alpha-0.25) < 1e-15)
check("assembled operator is self-adjoint", max(full["kinematic_residuals"].values()) < 1e-10)
res = table(h, modes, alpha)
check("half of the 384 levels are negative, with a gap", res["negative"] == 192 and res["gap"] > 1e-6,
      (res["negative"], round(res["gap"], 6)))
check("the whole circle is in a pure state (purity 1)", abs(res["whole"]) < 1e-9)

# ---------------------------------------------------------------- N2
print("N2  places and pairs of places, coefficients [2,2,2]")
np.set_printoptions(precision=6, suppress=True)
print("      log purity of one place:", np.round(res["single"], 6))
print("      log of the ratio for pairs of places:")
print(res["log_ratio"])
T = res["log_ratio"]
adj = [T[0, 1], T[1, 2], T[2, 3], T[3, 0]]
opp = [T[0, 2], T[1, 3]]
check("ratios are larger than one for every pair (the places are correlated)", min(adj+opp) > 0)
print("      neighbours:", np.round(adj, 6), " opposite:", np.round(opp, 6))
check("complementary pairs have equal purity (pure state): pur{0,2} = pur{1,3}, pur{0,1} = pur{2,3}",
      abs((T[0, 2]+res["single"][0]+res["single"][2])-(T[1, 3]+res["single"][1]+res["single"][3])) < 1e-8
      and abs((T[0, 1]+res["single"][0]+res["single"][1])-(T[2, 3]+res["single"][2]+res["single"][3])) < 1e-8)

# ---------------------------------------------------------------- N3
print("N3  comparisons")
z = np.zeros((96, 96), complex)
free = np.block([[(n+alpha)*A.gamma if n == m else z for m in modes] for n in modes])
rf = table(free, modes, alpha)
print("      no potential: one place", round(rf["single"][0], 6), " neighbours", round(rf["log_ratio"][0, 1], 6),
      " opposite", round(rf["log_ratio"][0, 2], 6))
per = np.array([rf["single"][0], rf["log_ratio"][0, 1], rf["log_ratio"][0, 2]])/96
print("      per fibre:", np.round(per, 8), " exp:", np.round(np.exp(per), 8))
check("without potential the 96 fibres contribute equally and the numbers are those of one shifted mode set",
      abs(rf["single"][0]-rf["single"][1]) < 1e-9)
rb = table(A.assemble(candidate([1, 1, 1]))["h_N"], modes, alpha)
print("      baseline [1,1,1]: one place", np.round(rb["single"], 6))
print("      baseline [1,1,1]: neighbours", round(rb["log_ratio"][0, 1], 6), " opposite", round(rb["log_ratio"][0, 2], 6))
print("      [2,2,2] against no potential, neighbours:", round(float(np.mean(adj))-rf["log_ratio"][0, 1], 6),
      " opposite:", round(float(np.mean(opp))-rf["log_ratio"][0, 2], 6))
check("the native potential changes the table (the functional reads the native couplings)",
      abs(float(np.mean(opp))-rf["log_ratio"][0, 2]) > 1e-6)

print()
if FAIL:
    print("FAILED:", FAIL)
    sys.exit(1)
print("all checks passed")
