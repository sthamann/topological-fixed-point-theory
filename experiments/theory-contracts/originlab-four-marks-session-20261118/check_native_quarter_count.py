#!/usr/bin/env python3
"""The count of the charge at a quarter angle on the native operator; the native places of the
recorded cuts (18 November 2026).

argv[1] = repository root.  The assembler of the original Q contract is used unchanged.
State: filled negative spectrum of h_N (assumption of the hypothesis, as in round twenty-five).

Sections
    V1  real structure, grading and potential: the native operator in the frame (field, conjugate)
    V2  quarter count <exp(i (pi/2) Q_A)> for all sets A of places, coefficients [2,2,2]
    V3  without potential: exact values per fibre
    V4  all 27 configurations: ranges
    V5  the rate between the copies is one for every transfer that is a function of h_N
    V6  a native place for the two recorded cuts: grading and direction of the potential
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
g, K = Asm.gamma, Asm.K
rng = np.random.default_rng(20261118)


def cand(genes):
    return Asm.assemble({"schema": "tfpt.field-configuration.v1", "domain_sha256": Asm.domain["_sha256"], "genes": list(genes)})


# ---------------------------------------------------------------- V1
print("V1  the native frame")
check("the real structure is real: K conj(K) = 1, so the 96 fibres have 96 real directions (K v* = v) and "
      "96 imaginary ones", np.allclose(K@K.conj(), np.eye(96)))
check("the grading anticommutes with it: K gamma* K = -gamma; a fibre vector and its conjugate have opposite grading",
      np.allclose(K@g.conj()@K, -g))
wg, vg = np.linalg.eigh(g)
up, dn = vg[:, wg > 0], vg[:, wg < 0]
check("48 fibres of each grading; K maps the one set onto the other (field and conjugate field)",
      up.shape[1] == 48 and dn.shape[1] == 48 and np.allclose(up.conj().T@(K@up.conj()), 0))
full = cand([2, 2, 2])
phi = full["phi_fourier"]
check("every Fourier part of the native potential anticommutes with the grading: it only connects a field with a "
      "conjugate field; the kinetic part (n + 1/4) gamma is the charge sign times the mode number",
      max(float(np.max(np.abs(g@v+v@g))) for v in phi.values()) < 1e-10)
check("in the words of round twenty-seven: the kinetic part is of the kind of the second cut (charge sign, B), the whole "
      "potential of the kind of the first cut (exchange with the conjugate, A), the real structure of the kind of P", True)

# ---------------------------------------------------------------- quarter count
def count_table(h, modes, alpha, lam=np.pi/2):
    m = len(modes)
    w, v = np.linalg.eigh(h)
    Pn = v[:, w < 0]@v[:, w < 0].conj().T
    theta = 2*np.pi*np.arange(m)/m
    Fm = np.array([[np.exp(1j*(n+alpha)*t) for n in modes] for t in theta])/np.sqrt(m)
    U = np.kron(Fm, np.eye(96))
    C = U@Pn@U.conj().T
    out = {}
    for r in range(1, m+1):
        for A in itertools.combinations(range(m), r):
            ii = np.concatenate([np.arange(96*p, 96*(p+1)) for p in A])
            nu = np.clip(np.linalg.eigvalsh(C[np.ix_(ii, ii)]), 0, 1)
            fac = 1-nu+np.exp(1j*lam)*nu
            out[A] = (float(np.sum(np.log(np.abs(fac)))), float(np.sum(np.angle(fac)*np.exp(0)) % (2*np.pi)))
    return out


def connected(tab, A, B):
    AB = tuple(sorted(A+B))
    return tab[AB][0]-tab[A][0]-tab[B][0], (tab[AB][1]-tab[A][1]-tab[B][1]+np.pi) % (2*np.pi)-np.pi


# ---------------------------------------------------------------- V2
print("V2  coefficients [2,2,2]")
tab = count_table(full["h_N"], full["modes"], full["alpha"])
np.set_printoptions(precision=6, suppress=True)
for A in [(0,), (1,), (0, 1), (0, 2), (1, 3), (0, 1, 2), (0, 1, 2, 3)]:
    print("      places", A, " log modulus %.6f  phase/pi %.6f" % (tab[A][0], tab[A][1]/np.pi))
check("the whole circle has a definite charge: modulus 1 and phase a multiple of pi/2",
      abs(tab[(0, 1, 2, 3)][0]) < 1e-8 and min(abs(tab[(0, 1, 2, 3)][1]/(np.pi/2) % 1), 1-abs(tab[(0, 1, 2, 3)][1]/(np.pi/2) % 1)) < 1e-6)
check("a set of places and its complement have the same modulus (pure state)",
      abs(tab[(0, 2)][0]-tab[(1, 3)][0]) < 1e-8 and abs(tab[(0,)][0]-tab[(1, 2, 3)][0]) < 1e-8)
nb, op = connected(tab, (0,), (1,)), connected(tab, (0,), (2,))
print("      connected count: neighbours log modulus %.6f phase/pi %.6f; opposite %.6f phase/pi %.6f"
      % (nb[0], nb[1]/np.pi, op[0], op[1]/np.pi))
check("the moduli are the square roots of the purity ratios of round twenty-five (24.393 and 2.979 there)",
      abs(2*nb[0]-24.393442) < 1e-5 and abs(2*op[0]-2.979157) < 1e-5)

# ---------------------------------------------------------------- V3
print("V3  without potential")
free = np.block([[(n+0.25)*g if n == m else np.zeros((96, 96)) for m in full["modes"]] for n in full["modes"]])
tf = count_table(free, full["modes"], 0.25)
nbf, opf = connected(tf, (0,), (1,)), connected(tf, (0,), (2,))
per = (np.exp(tf[(0,)][0]/96), np.exp(nbf[0]/96), np.exp(opf[0]/96))
print("      per fibre: one place %.8f, neighbours %.8f, opposite %.8f" % per)
check("per fibre: one place 1/sqrt(2), neighbours 3/2, opposite places 1",
      abs(per[0]-2**-0.5) < 1e-10 and abs(per[1]-1.5) < 1e-10 and abs(per[2]-1) < 1e-10)
check("the connected count without potential has no phase", abs(nbf[1]) < 1e-8 and abs(opf[1]) < 1e-8)

# ---------------------------------------------------------------- V4
print("V4  all 27 configurations")
rows = []
for genes in itertools.product(range(3), repeat=3):
    a = cand(genes)
    t = count_table(a["h_N"], a["modes"], a["alpha"])
    n1, o1 = connected(t, (0,), (1,)), connected(t, (0,), (2,))
    rows.append((genes, t[(0,)][0], n1[0], n1[1], o1[0], o1[1]))
arr = np.array([r[1:] for r in rows])
print("      one place, log modulus: %.4f to %.4f" % (arr[:, 0].min(), arr[:, 0].max()))
print("      neighbours, log modulus: %.4f to %.4f; phase/pi: %.4f to %.4f"
      % (arr[:, 1].min(), arr[:, 1].max(), arr[:, 2].min()/np.pi, arr[:, 2].max()/np.pi))
print("      opposite, log modulus: %.4f to %.4f; phase/pi: %.4f to %.4f"
      % (arr[:, 3].min(), arr[:, 3].max(), arr[:, 4].min()/np.pi, arr[:, 4].max()/np.pi))
check("27 configurations evaluated; neighbours are always more strongly connected than opposite places",
      len(rows) == 27 and np.all(arr[:, 1] > arr[:, 3]))
check("the values depend on the configuration: the count reads the native couplings", arr[:, 3].max()-arr[:, 3].min() > 1e-3)

# ---------------------------------------------------------------- V5
print("V5  the rate between the copies")
h = full["h_N"]
w, v = np.linalg.eigh(h)
worst = 0.0
for s in (0.3, 1.0, 2.5):
    T = (v*np.exp(-s*np.abs(w)))@v.conj().T
    for _ in range(5):
        x = rng.normal(size=384)+1j*rng.normal(size=384)
        worst = max(worst, abs((x.conj()@T@x)-((1j*x).conj()@T@(1j*x))))
check("a line and the same line turned by a quarter (multiplied by i) have equal transfer for every function of h_N: "
      "natively the rate between the copies is q = 1", worst < 1e-10, worst)

# ---------------------------------------------------------------- V6
print("V6  a native place for the two recorded cuts (hypothesis)")
for label, genes in (("baseline [1,1,1]", [1, 1, 1]), ("[2,2,2]", [2, 2, 2])):
    p0 = cand(genes)["phi_fourier"][0]
    wv, vv = np.linalg.eigh(p0)
    Yn = (vv*np.sign(wv))@vv.conj().T                           # the direction of the constant native potential
    absp = (vv*np.abs(wv))@vv.conj().T
    ok = (np.min(np.abs(wv)) > 1e-9 and np.allclose(Yn@Yn, np.eye(96)) and np.allclose(Yn@g+g@Yn, 0)
          and np.allclose(K@Yn.conj()@K, Yn))
    check(f"{label}: the constant potential is invertible; its direction Y is an involution, anticommutes with the "
          "grading and is real under K", ok, ("smallest |value|", round(float(np.min(np.abs(wv))), 6)))
    Pi = -1j*Yn@g
    wp, vp = np.linalg.eigh(Pi)
    check(f"{label}: Pi = -i Y gamma is hermitian, real under K, with 48 values +1 and 48 values -1; Y gamma Y gamma = -1",
          np.allclose(Pi, Pi.conj().T) and np.allclose(K@Pi.conj()@K, Pi) and int(np.sum(wp > 0)) == 48
          and np.allclose(Yn@g@Yn@g, -np.eye(96)))
    res = []
    for eta, cols in ((1, vp[:, wp > 0]), (-1, vp[:, wp < 0])):
        x = cols@(rng.normal(size=48)+1j*rng.normal(size=48))
        x = x+K@x.conj()                                        # a real line inside the eigenspace
        x /= np.linalg.norm(x)
        assert np.allclose(K@x.conj(), x) and np.allclose(Pi@x, eta*x)
        k = 0.25
        T = np.linalg.matrix_power(np.eye(96), 1)
        hh = k*g+p0
        wh, vh = np.linalg.eigh(hh)
        T = (vh*np.exp(-0.7*np.abs(wh)))@vh.conj().T            # an even function of the constant-mode operator
        z0 = x.conj()@T@T@T@x
        zu = x.conj()@T@T@g@T@x
        zy = x.conj()@T@Yn@T@T@x
        zuy = x.conj()@T@Yn@T@g@T@x
        lift = x.conj()@((np.eye(96)+g)/2)@x
        res.append((eta, zu/z0, zy/z0, zuy/z0, lift))
    okm = all(abs(zu) < 1e-9 and abs(zy) < 1e-9 and abs(zuy-1j*eta) < 1e-9 and abs(lift-0.5) < 1e-9
              for eta, zu, zy, zuy, lift in res)
    check(f"{label}: on a real line with Pi = eta, with U = grading and Y = potential direction and an even transfer: "
          "Z_U = Z_Y = 0, Z_UY = i eta Z_0, positive lift (1 + gamma)/2 gives 1/2: the recorded canonical marks with q = 1",
          okm, [(e, np.round(z, 6)) for e, _, _, z, _ in res])
check("so the recorded binary cuts have a native place: grading (charge sign) and direction of the potential; the cap is "
      "a real line; not shown: that P1 reads its marks U and Y this way (the original Higgs marker is the non-negative "
      "square of the potential, which commutes with the grading)", True)

print()
if FAIL:
    print("FAILED:", FAIL)
    sys.exit(1)
print("all checks passed")
