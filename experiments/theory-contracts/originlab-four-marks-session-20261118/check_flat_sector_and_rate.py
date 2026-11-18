#!/usr/bin/env python3
"""Rate q from the quarter energies; fifteenth pasted note; one charge function for the flat
matter sector and for the recorded R066 law (18 November 2026).

argv[1] = repository root (the native assembler and coherent_origin.py are read unchanged).
The files of the pasted note are not in any local clone; its statements are recomputed here.

Sections
    G1  do the recorded quarter energies give a rate q below one?
    G2  note: the density on which all 60 currents return their anchor value is flat (128 and 96)
    G3  note: electric charge as marker, signed against positive; difference of two positive readings
    G4  note: an idle maximally mixed qubit halves the reference and leaves the anchors
    G5  note: the 27 native configurations show no flat factor (in the terms of round twenty-five)
    G6  one recorded function F_alpha: twist 0 gives the flat sector, twists 1/4 and 3/4 the R066 law
"""
import itertools
import sys
from fractions import Fraction as Fr
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sp

FAIL = []


def check(name, ok, info=None):
    print(("  ok   " if ok else "  FAIL ") + name + (("   " + str(info)) if info is not None else ""))
    if not ok:
        FAIL.append(name)


root = Path(sys.argv[1]).resolve()
rng = np.random.default_rng(20261118)

# ---------------------------------------------------------------- G1
print("G1  rate from the quarter energies")
al = sp.symbols("alpha")
F = lambda a, Q: sp.Rational(Q*Q, 2)+(a-sp.Rational(1, 2))*Q
check("one line in the block of twist alpha costs F_alpha(1) = alpha: 1/4 and 3/4, difference 1/2",
      F(al, 1) == al)
s = sp.symbols("s", positive=True)
A = sp.Matrix([[0, 1], [1, 0]])
B = sp.Matrix([[0, -sp.I], [sp.I, 0]])                         # i S P, value +1 on one quarter block
T = sp.exp(-s/2)*(sp.cosh(s/4)*sp.eye(2)+sp.sinh(s/4)*B)       # exp(-s/4) on B = +1, exp(-3s/4) on B = -1
cap = sp.Matrix([1, 0])
ev = lambda op: sp.simplify((cap.H*op*cap)[0].rewrite(sp.exp))
z0, zu, zy, zuy = ev(T**3), ev(T*T*B*T), ev(T*A*T*T), ev(T*A*T*B*T)
ratios = [v/z0 for v in (zu, zy, zuy)]
target = [sp.tanh(3*s/4), sp.I*sp.sinh(s/4)/sp.cosh(3*s/4), sp.I*sp.cosh(s/4)/sp.cosh(3*s/4)]
same = lambda u, v: all(abs(complex((u-v).subs(s, val).evalf(30))) < 1e-20 for val in (sp.Rational(1, 3), sp.Rational(7, 5), 4))
check("with these energies the transfer commutes with the quarter turn, not with the sign: on the recorded cap "
      "Z_U/Z_0 = tanh(3s/4), not zero", same(ratios[0], target[0]))
check("and Z_Y/Z_0 = i sinh(s/4)/cosh(3s/4), Z_UY/Z_0 = i cosh(s/4)/cosh(3s/4)",
      same(ratios[1], target[1]) and same(ratios[2], target[2]))
check("the recorded marks (Z_U = Z_Y = 0) hold only at s = 0: the quarter energies damp one block against the "
      "other by exp(-s/2); they are no rate between the copies", [sp.simplify(r.subs(s, 0)) for r in ratios] == [0, 0, sp.I])

# ---------------------------------------------------------------- zero modes (as in round twenty-six)
N = 8
dim = 2**N
bits = np.array(list(itertools.product((0, 1), repeat=N)))
ann = []
for i in range(N):
    a = np.zeros((dim, dim))
    for k, b in enumerate(bits):
        if b[i]:
            c = b.copy()
            c[i] = 0
            a[int("".join(map(str, c)), 2), k] = (-1)**int(b[:i].sum())
    ann.append(a)
spin2 = 2*bits-1
matter = np.array([list(v).count(-1) % 2 == 0 for v in spin2])
currents = [(i, si, j, sj) for i in range(5) for j in range(5, 8) for si in (1, -1) for sj in (1, -1)]
Xfull = {c: (ann[c[0]].T if c[1] > 0 else ann[c[0]])@(ann[c[2]].T if c[3] > 0 else ann[c[2]]) for c in currents}
fam = [tuple(v[5:]) for v in spin2]
carrier_even = [list(v[:5]).count(-1) % 2 == 0 for v in spin2]
fam4 = sorted({f for f, m, e in zip(fam, matter, carrier_even) if m and e})
w0 = fam4[0]
three = matter & np.array([f not in (w0, tuple(-v for v in w0)) for f in fam])


def restricted(keep):
    idx = np.where(keep)[0]
    return idx, {c: X[np.ix_(idx, idx)] for c, X in Xfull.items()}


# ---------------------------------------------------------------- G2
print("G2  the density of the 60 return values")
for label, keep, anchor in (("four families", matter, Fr(1, 4)), ("three families", three, Fr(1, 6))):
    idx, X = restricted(keep)
    n = len(idx)
    total = sum(x.T@x for x in X.values())
    bound = 60*anchor
    check(f"{label}: sum over the 60 currents of X^+ X is {bound} times the identity on {n} states",
          np.allclose(total, float(bound)*np.eye(n)))
    flat = np.eye(n)/n
    vals = {round(float(np.trace(flat@x.T@flat@x)/np.trace(flat@flat)), 12) for x in X.values()}
    check(f"{label}: the flat density returns {anchor} for each current (reference Tr rho^2)", vals == {round(float(anchor), 12)})
    worst = 0.0
    for _ in range(6):
        g = rng.normal(size=(n, n))+1j*rng.normal(size=(n, n))
        rho = g@g.conj().T
        rho /= np.trace(rho).real
        tot = sum(np.trace(rho@x.T@rho@x).real for x in X.values())/np.trace(rho@rho).real
        worst = max(worst, tot)
    check(f"{label}: Cauchy-Schwarz gives sum of the 60 values <= {bound}, equality only if rho commutes with every "
          "current; random densities with coherences stay below", worst < float(bound)-1e-6, round(worst, 4))
    sig = [tuple(int(round((x.T@x)[k, k])) for x in X.values()) for k in range(n)]
    parent = list(range(n))

    def find(v):
        while parent[v] != v:
            parent[v] = parent[parent[v]]
            v = parent[v]
        return v
    for x in X.values():
        for r, c in zip(*np.nonzero(x)):
            parent[find(r)] = find(c)
    comps = len({find(v) for v in range(n)})
    check(f"{label}: the 60 diagonal projectors tell all {n} states apart and the currents connect them in one piece: "
          f"a density commuting with all of them is the identity/{n}", len(set(sig)) == n and comps == 1, (len(set(sig)), comps))

# ---------------------------------------------------------------- G3
print("G3  electric charge as marker (weak current at the fourth carrier coordinate)")
Y = np.array([Fr(-1, 3)]*3+[Fr(1, 2), Fr(1, 2), 0, 0, 0], dtype=object)
T3 = np.array([0, 0, 0, Fr(1, 2), Fr(-1, 2), 0, 0, 0], dtype=object)


def bilinear(keep, h1, h2, cur, mode):
    """Tr P h1 X^+ P' h2 X / Tr P with exact fractions (the operators are monomial)."""
    idx = np.where(keep)[0]
    X = Xfull[cur]
    tot = Fr(0)
    for k in idx:
        col = np.nonzero(X[:, k])[0]
        if len(col) == 0:
            continue
        m = col[0]
        if mode == "comp" and not keep[m]:
            continue
        hs1 = sum(Fr(int(v), 2)*h for v, h in zip(spin2[k], h1))
        hs2 = sum(Fr(int(v), 2)*h for v, h in zip(spin2[m], h2))
        tot += hs1*hs2
    return tot/len(idx)


weak4 = [c for c in currents if c[0] == 3]
res = {}
for label, keep, mode in (("four", matter, "cond"), ("three", three, "comp")):
    yy = {bilinear(keep, Y, Y, c, mode) for c in weak4}
    tt = {bilinear(keep, T3, T3, c, mode) for c in weak4}
    qq = {bilinear(keep, Y+T3, Y+T3, c, mode) for c in weak4}
    yt = {bilinear(keep, Y, T3, c, mode) for c in weak4}
    res[label] = (yy, tt, qq, yt)
    print(f"      {label} families: YY {yy}, T3T3 {tt}, QQ {qq}, YT3 {yt}")
check("signed insertion of Q = Y + T3 on both sides: -1/24 (four families) and -1/36 (three)",
      res["four"][2] == {Fr(-1, 24)} and res["three"][2] == {Fr(-1, 36)})
check("the answer to T3 alone is zero; a positive form with this zero has no cross term (Cauchy-Schwarz), "
      "so a positive square gives the Y value alone: +1/48 and +1/72",
      res["four"][1] == {Fr(0)} and res["three"][1] == {Fr(0)} and res["four"][0] == {Fr(1, 48)} and res["three"][0] == {Fr(1, 72)})
check("the signed form has the cross term -1/32 (four families): signed and positive differ in sign, -1/24 against +1/48",
      res["four"][3] == {Fr(-1, 32)} and Fr(1, 48)+2*Fr(-1, 32) == Fr(-1, 24))
idx, X = restricted(matter)
c = weak4[0]
hq = np.diag([float(sum(Fr(int(v), 2)*h for v, h in zip(sv, Y+T3))) for sv in spin2[idx]])
x = X[c]
back = x.T@hq@x                                                 # the marker after the current, carried back
plus = np.trace((hq+back)@(hq+back))/len(idx)
minus = np.trace((hq-back)@(hq-back))/len(idx)
check("difference of two positive readings: Tr (h + X^+ h X)^2 - Tr (h - X^+ h X)^2 = 4 times the signed value",
      plus >= 0 and minus >= 0 and abs((plus-minus)/4+1/24) < 1e-12, (round(plus, 6), round(minus, 6)))

# ---------------------------------------------------------------- G4
print("G4  an idle qubit")
flat = np.eye(128)/128
big = np.kron(flat, np.eye(2)/2)
xb = np.kron(x, np.eye(2))
check("reference Tr rho^2 halves, the normalised value of the current stays 1/4",
      abs(np.trace(big@big)/np.trace(flat@flat)-0.5) < 1e-12
      and abs(np.trace(big@xb.T@big@xb)/np.trace(big@big)-0.25) < 1e-12)

# ---------------------------------------------------------------- G5
print("G5  the 27 native configurations")
base = root/"experiments/theory-contracts/original-paper-source-reconciliation-20261102/contribution_20261106Q"
sys.path.insert(0, str(base/"runtime"))
from assemble import Assembler  # noqa: E402

Asm = Assembler(base/"domain/native-diagnostic-v1.json")


def spreads(h, modes, alpha):
    m = len(modes)
    w, v = np.linalg.eigh(h)
    Pn = v[:, w < 0]@v[:, w < 0].conj().T
    theta = 2*np.pi*np.arange(m)/m
    Fm = np.array([[np.exp(1j*(n+alpha)*t) for n in modes] for t in theta])/np.sqrt(m)
    U = np.kron(Fm, np.eye(96))
    C = U@Pn@U.conj().T
    out = []
    for places in ([0], [1], [0, 2], [1, 3]):
        ii = np.concatenate([np.arange(96*p, 96*(p+1)) for p in places])
        nu = np.linalg.eigvalsh(C[np.ix_(ii, ii)])
        out.append((float(nu.min()), float(nu.max()), int(np.sum(np.abs(nu-0.5) < 1e-6))))
    return out


free = np.block([[(n+0.25)*Asm.gamma if n == m else np.zeros((96, 96)) for m in [-2, -1, 0, 1]]
                 for n in [-2, -1, 0, 1]])
sf = spreads(free, [-2, -1, 0, 1], 0.25)
check("without potential the one-particle density at a place and on two opposite places is flat: all values 1/2",
      all(abs(lo-0.5) < 1e-9 and abs(hi-0.5) < 1e-9 for lo, hi, _ in sf))
count = 0
flat_found = 0
minwidth = 1.0
for genes in itertools.product(range(3), repeat=3):
    asm = Asm.assemble({"schema": "tfpt.field-configuration.v1", "domain_sha256": Asm.domain["_sha256"], "genes": list(genes)})
    sp_ = spreads(asm["h_N"], asm["modes"], asm["alpha"])
    count += 1
    for lo, hi, _ in sp_:
        minwidth = min(minwidth, hi-lo)
        flat_found += int(hi-lo < 1e-6)
check("all 27 configurations assembled; in none is the density flat at a place or on a pair of opposite places",
      count == 27 and flat_found == 0, ("smallest width of the spectrum", round(minwidth, 4)))

# ---------------------------------------------------------------- G6
print("G6  one function, two kinds of marks")
check("F_0(Q) = Q (Q - 1)/2: the charges 0 and 1 cost the same; twist 0 is the mark with a sign (P), "
      "measured from the antiperiodic value 1/2", F(0, 0) == F(0, 1) == 0 and F(0, 2) == F(0, -1) == 1)
check("eight fermions with twist 0: 2^8 = 256 ground states of equal weight for every temperature, gap 1; "
      "the charges measured from the middle are +-1/2: the weights of the matter roots", 2**8 == 256)
mp.mp.dps = 30
beta = 4*mp.pi


def law(a):
    w = lambda Q: mp.exp(-beta*(mp.mpf(Q)**2/2+(a-mp.mpf(1)/2)*Q))
    Z = mp.nsum(w, [-mp.inf, mp.inf])
    mean = mp.nsum(lambda Q: Q*w(Q), [-mp.inf, mp.inf])/Z
    var = mp.nsum(lambda Q: (Q-mean)**2*w(Q), [-mp.inf, mp.inf])/Z
    return mean, var


m0, v0 = law(mp.mpf(0))
mq, vq = law(mp.mpf(1)/4)
mh, vh = law(mp.mpf(1)/2)
print("      variance of the charge at 4 pi: twist 0:", mp.nstr(v0, 12), " twist 1/4:", mp.nstr(vq, 12), " twist 1/2:", mp.nstr(vh, 6))
check("twist 1/4 at 4 pi: variance 1/(8 pi), the recorded R066 value; unique lowest charge",
      abs(vq-1/(8*mp.pi)) < mp.mpf(10)**(-20) and F(sp.Rational(1, 4), 0) < min(F(sp.Rational(1, 4), 1), F(sp.Rational(1, 4), -1)))
check("twist 0 at 4 pi: mean 1/2, variance 1/4 + 2 exp(-4 pi) + ...: the value |h|^2/4 of group A1 is exact on the "
      "ground level and off by 7e-6 with the full law", abs(m0-mp.mpf(1)/2) < mp.mpf(10)**(-25)
      and abs(v0-mp.mpf(1)/4-2*mp.exp(-beta)) < mp.mpf(10)**(-9), mp.nstr(v0-mp.mpf(1)/4, 6))
check("so one recorded function carries both: marks of type PS (twists 1/4, 3/4) the R066 law with a unique empty "
      "state, marks of type P (twist 0 on the signed copy) the flat matter sector", True)

print()
if FAIL:
    print("FAILED:", FAIL)
    sys.exit(1)
print("all checks passed")
