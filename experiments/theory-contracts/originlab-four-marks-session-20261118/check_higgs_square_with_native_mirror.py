#!/usr/bin/env python3
"""Does the Higgs square give the negative mixed answer with the native mirror? (18 November 2026)

argv[1] = repository root.  Native objects through the unchanged assembler; pair currents carried
to the native fibres by the bijection of round thirty-eight (imported from that script's logic,
rebuilt here).  Forms on the flat glued state of the 96 fibres:
    plain mirror   Q[F,G]   = Tr(F^+ G)/96
    K-mirror       Q_K[F,G] = Tr((K conj(F) K)^+ G)/96
Connected mixed answer of a marker U through a hermitian cut X:  Q[XU, UX] - Q[XU, XU].

Sections
    H1  parity of the cuts under the real structure
    H2  Higgs square through the direction of the native potential
    H3  Higgs square through the 60 carried pair currents, both mirrors
    H4  hypercharge through the same cuts, both mirrors
"""
import itertools
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

Asm = Assembler(base/"domain/native-diagnostic-v1.json")
orig, g, K = Asm.original, Asm.gamma, Asm.K
I96 = np.eye(96)

# ---------------------------------------------------------------- chart and bijection (round thirty-eight)
bits = np.array(list(itertools.product((0, 1), repeat=8)))
ann = []
for i in range(8):
    a = np.zeros((256, 256))
    for k, b in enumerate(bits):
        if b[i]:
            c = b.copy()
            c[i] = 0
            a[int("".join(map(str, c)), 2), k] = (-1)**int(b[:i].sum())
    ann.append(a)
spin2 = 2*bits-1
matter = np.array([list(v).count(-1) % 2 == 0 for v in spin2])
even_fam = sorted({tuple(v[5:]) for v in spin2 if list(v[5:]).count(-1) % 2 == 0})
kept = even_fam[1:]
three = matter & np.array([tuple(v[5:]) not in (even_fam[0], tuple(-x for x in even_fam[0])) for v in spin2])
idx = np.where(three)[0]
pos = {int(k): n for n, k in enumerate(idx)}


def fibre(s, perm):
    particle = int(np.sum(s[:3] > 0)) in (1, 3)
    t = s if particle else -s
    quark = int(np.sum(t[:3] > 0)) == 1
    colour = int(np.argmax(t[:3] > 0)) if quark else 0
    name = ("u" if t[3] > 0 else "d") if quark else ("nu" if t[3] > 0 else "e")
    ch = "L" if t[3] != t[4] else "R"
    fam = tuple(int(x) for x in t[5:])
    fam = fam if fam in kept else tuple(-x for x in fam)
    return (0 if particle else 48)+orig.IDX[(name, ch, colour, perm[kept.index(fam)])]


def carried(perm):
    Phi = np.zeros((96, 96))
    for k in idx:
        Phi[fibre(spin2[k], perm), pos[int(k)]] = 1
    out = {}
    for i in range(5):
        for j in range(5, 8):
            for si in (1, -1):
                for sj in (1, -1):
                    X = ((ann[i].T if si > 0 else ann[i])@(ann[j].T if sj > 0 else ann[j]))[np.ix_(idx, idx)]
                    out[(i, j, si, sj)] = Phi@X@Phi.T
    return out


def hermitian_cuts(Xn):
    cuts = []
    for (i, j, si, sj), X in Xn.items():
        if si > 0:                                               # one representative of each pair of opposite roots
            cuts.append(("weak" if i >= 3 else "colour", (X+X.T)/np.sqrt(2)))
            cuts.append(("weak" if i >= 3 else "colour", 1j*(X-X.T)/np.sqrt(2)))
    return cuts


Q = lambda F, G: (np.trace(F.conj().T@G)/96).real
QK = lambda F, G: (np.trace((K@F.conj()@K).conj().T@G)/96).real
conn = lambda form, U, X: form(X@U, U@X)-form(X@U, X@U)
kpar = lambda X: 1 if np.allclose(K@X.conj()@K, X) else (-1 if np.allclose(K@X.conj()@K, -X) else 0)

DF = Asm.assemble({"schema": "tfpt.field-configuration.v1", "domain_sha256": Asm.domain["_sha256"], "genes": [1, 1, 1]})["phi_fourier"][0]
wv, vv = np.linalg.eigh(-1j*g@DF)
Ypot = (vv*np.sign(wv))@vv.conj().T
H2 = sum(v@v for v in Asm.V)
U = H2/np.sqrt(np.trace(H2@H2).real/96)                          # Higgs square with <U^2> = 1
YT = {("u", "L"): 1/6, ("d", "L"): 1/6, ("u", "R"): 2/3, ("d", "R"): -1/3, ("nu", "L"): -1/2, ("e", "L"): -1/2, ("nu", "R"): 0.0, ("e", "R"): -1.0}
yv = [YT[(n, c)] for (n, c, _, _) in orig.P]
Yh = np.diag(yv+[-v for v in yv])

# ---------------------------------------------------------------- H1
print("H1  parity of the cuts")
cuts = hermitian_cuts(carried((0, 1, 2)))
pars = [kpar(X) for _, X in cuts]
check("the 60 hermitian combinations of the carried currents have a definite parity under K: 32 even and 28 odd "
      "(the first run expected 30 and 30)", len(cuts) == 60 and pars.count(1) == 32 and pars.count(-1) == 28,
      (pars.count(1), pars.count(-1), pars.count(0)))
check("the Higgs square and the direction of the native potential are even, the hypercharge is odd",
      kpar(U) == 1 and kpar(Ypot) == 1 and kpar(Yh) == -1)

# ---------------------------------------------------------------- H2
print("H2  Higgs square through the direction of the native potential")
c_plain, c_k = conn(Q, U, Ypot), conn(QK, U, Ypot)
com = U@Ypot-Ypot@U
print("      connected answer: plain mirror %.6f, K-mirror %.6f; -(1/2)<|[U,Y]|^2> = %.6f" % (c_plain, c_k, -0.5*Q(com, com)))
check("the Higgs square does not commute with the direction of the potential; the connected answer is negative and "
      "the same for both mirrors, = -(1/2)<|[U,Y]|^2>", c_plain < -1e-6 and abs(c_plain-c_k) < 1e-10 and abs(c_plain+0.5*Q(com, com)) < 1e-10)

# ---------------------------------------------------------------- H3, H4
print("H3  Higgs square through the carried pair currents")
rows = []
for perm in itertools.permutations(range(3)):
    for kind, X in hermitian_cuts(carried(perm)):
        rows.append((perm, kind, kpar(X), conn(Q, U, X), conn(QK, U, X), conn(Q, Yh, X), conn(QK, Yh, X), Q(X, X)))
arr = np.array([[r[2], r[3], r[4], r[5], r[6], r[7]] for r in rows], float)
kinds = np.array([r[1] for r in rows])
print("      plain mirror, Higgs square: weak %.6f to %.6f, colour %.6f to %.6f"
      % (arr[kinds == "weak", 1].min(), arr[kinds == "weak", 1].max(), arr[kinds == "colour", 1].min(), arr[kinds == "colour", 1].max()))
check("plain mirror: the Higgs square gives a negative connected answer through every carried current, in all six "
      "orders of the generations", arr[:, 1].max() < -1e-9)
check("K-mirror: the answer keeps its value for the cuts that are even under K and changes sign for the odd ones",
      np.allclose(arr[:, 2], arr[:, 0]*arr[:, 1]))
check("so with the native mirror the Higgs square gives the negative answer through the even cuts only (32 of 60): the sign "
      "is decided by the parity of the cut, not by the marker being non-negative",
      (arr[arr[:, 0] > 0, 2].max() < 0) and (arr[arr[:, 0] < 0, 2].min() > 0))
vals = sorted({round(float(v), 6) for v in arr[:, 1]})
check("the values depend on the current and on the order of the generations; they are not the recorded -1/48 and "
      "-1/108", len(vals) > 2, (vals[0], vals[-1], len(vals)))
print("H4  hypercharge through the same cuts")
check("plain mirror: -1/48 for the weak and -1/108 for the colour currents, for every order of the generations",
      np.allclose(arr[kinds == "weak", 3], -1/48) and np.allclose(arr[kinds == "colour", 3], -1/108))
check("K-mirror: the sign of the hypercharge answer is the opposite of the parity of the cut (the hypercharge is odd): "
      "negative through the odd cuts, positive through the even ones", np.allclose(arr[:, 4], -arr[:, 0]*arr[:, 3]))

print()
if FAIL:
    print("FAILED:", FAIL)
    sys.exit(1)
print("all checks passed")
