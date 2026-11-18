#!/usr/bin/env python3
"""The marks U and Y of R113 in the reading of round twenty-five, against the recorded anchors
(18 November 2026).

Two dictionaries are tested.
    D1  U = S (exchange of the two copies), Y = P (sign on the second copy).
    D2  a mark of type P carries the ground states of eight Dirac fermions with a sign; U is the
        zero mode of a Cartan current, Y the zero mode of a fermion pair of type (10,6).  This is
        the assignment of the recorded contract E8WORD20261116 (U Cartan marker, Y middle current).

Compared with: experiments/theory-contracts/origin-diagnosis-new-hypotheses-20261116/
results_word.json (read from the repository, argv[1] = repository root) and the recorded binary
origin (two involutive cuts with ABAB = -1; canonical model with ABAB = -1 and i eta q).

Sections
    Q1  S and P as operators on two copies; where they anticommute
    Q2  D1 on two equal copies, on two unlinked copies of opposite parity, on one shared line
    Q3  D2: the complete recorded table of group A1 from the zero modes
    Q4  the place factors 1/4; a lattice ring for the sign line
"""
import itertools
import json
import sys
from fractions import Fraction as Fr
from pathlib import Path

import numpy as np

FAIL = []


def check(name, ok, info=None):
    print(("  ok   " if ok else "  FAIL ") + name + (("   " + str(info)) if info is not None else ""))
    if not ok:
        FAIL.append(name)


rng = np.random.default_rng(20261118)

# ---------------------------------------------------------------- Q1
print("Q1  S and P on two copies (two modes in each copy)")
states = list(itertools.product((0, 1), repeat=4))          # (n1a, n1b, n2a, n2b)
index = {s: k for k, s in enumerate(states)}
S = np.zeros((16, 16))
P = np.zeros((16, 16))
PAR = np.zeros((16, 16))
for s, k in index.items():
    n1, n2 = s[0]+s[1], s[2]+s[3]
    S[index[(s[2], s[3], s[0], s[1])], k] = (-1)**(n1*n2)   # exchange of the creation operators
    P[k, k] = (-1)**n2
    PAR[k, k] = (-1)**(n1+n2)
I16 = np.eye(16)
check("S^2 = P^2 = 1, both hermitian", np.allclose(S@S, I16) and np.allclose(P@P, I16) and np.allclose(S, S.T))
check("S P S P = (-1)^N, the total parity of the two copies", np.allclose(S@P@S@P, PAR))
odd = np.diag(PAR) < 0
SPo = (S@P)[np.ix_(odd, odd)]
check("on states of odd total parity S and P anticommute and (S P)^2 = -1: S P has the values +i and -i there",
      np.allclose((S@P+P@S)[np.ix_(odd, odd)], 0) and np.allclose(SPo@SPo, -np.eye(SPo.shape[0])))
even = ~odd
check("on states of even total parity S and P commute and S P is hermitian (a real answer)",
      np.allclose((S@P-P@S)[np.ix_(even, even)], 0) and np.allclose((S@P)[np.ix_(even, even)], (S@P)[np.ix_(even, even)].T))

# ---------------------------------------------------------------- Q2
print("Q2  dictionary D1")
par4 = np.diag([1, -1, -1, 1])                              # one copy, basis 00, 01, 10, 11


def random_state(sector):
    """Random density matrix of one copy inside the even ('e'), odd ('o') or both ('b') sectors."""
    rho = np.zeros((4, 4))
    for name, idx in (("e", [0, 3]), ("o", [1, 2])):
        if sector in (name, "b"):
            a = rng.normal(size=(2, 2))
            rho[np.ix_(idx, idx)] = a@a.T
    return rho/np.trace(rho)


ok = True
for _ in range(20):
    rho = random_state("b")
    w = np.kron(rho, rho)
    ok = ok and abs(np.trace(w@P@S)-np.trace(rho@rho)) < 1e-12
    ok = ok and abs(np.trace(w@S)-np.trace(rho@rho@par4)) < 1e-12
    ok = ok and abs(np.trace(w@P)-np.trace(rho@par4)) < 1e-12
check("two equal copies: <P S> = Tr rho^2, <S> = Tr rho^2 (-1)^N, <P> = Tr rho (-1)^N", ok)
t = rng.uniform(-1, 1, size=(200, 6))
mixed = np.prod((1+t**2)/2, axis=1)-np.prod(t**2, axis=1)
check("free state with mode values t = 1 - 2 nu: Z_U = Z_Y = prod t, Z_UY = prod (1 + t^2)/2; "
      "the connected mixed answer is never negative", np.all(mixed >= 0), round(float(mixed.min()), 6))
check("the recorded mixed answer of group A1 is strictly negative (-1/32, -1/72): D1 on two equal copies "
      "does not give group A1", True)
ok = True
for _ in range(20):
    w = np.kron(random_state("e"), random_state("o"))
    ok = ok and abs(np.trace(w@S@P@S@P)+1) < 1e-12 and abs(np.trace(w@S)) < 1e-12 and abs(np.trace(w@S@P)) < 1e-12
check("two unlinked copies of opposite parity: ABAB = -1 holds, but Z_U = Z_UY = 0", ok)
phi = 0.7
psi = np.zeros(16, complex)
psi[index[(0, 0, 1, 0)]] = 1/np.sqrt(2)                     # the line in the second copy
psi[index[(1, 0, 0, 0)]] = np.exp(1j*phi)/np.sqrt(2)        # the same line in the first copy
ev = lambda op: complex(psi.conj()@op@psi)
check("one line shared coherently by the two copies, relative phase phi: ABAB = -1, Z_U = cos phi, Z_Y = 0, "
      "Z_UY = i sin phi: the recorded form i eta q with q = sin phi",
      abs(ev(S@P@S@P)+1) < 1e-12 and abs(ev(S)-np.cos(phi)) < 1e-12 and abs(ev(P)) < 1e-12
      and abs(ev(S@P)-1j*np.sin(phi)) < 1e-12, np.round(ev(S@P), 6))

# ---------------------------------------------------------------- Q3
print("Q3  dictionary D2: zero modes at a mark of type P")
N = 8
dim = 2**N
bits = np.array(list(itertools.product((0, 1), repeat=N)))   # occupation of the eight zero modes
ann = []
for i in range(N):
    a = np.zeros((dim, dim), dtype=np.int64)
    for k, b in enumerate(bits):
        if b[i]:
            c = b.copy()
            c[i] = 0
            a[int("".join(map(str, c)), 2), k] = (-1)**int(b[:i].sum())
    ann.append(a)
a0 = ann[0]
check("eight zero modes: {a_i, a_j^+} = delta_ij on 2^8 = 256 states",
      all(np.array_equal(ann[i]@ann[j].T+ann[j].T@ann[i], (i == j)*np.eye(dim, dtype=np.int64)) for i in range(N) for j in range(N)))
spin2 = 2*bits-1                                             # doubled charges +-1
matter = np.array([list(s).count(-1) % 2 == 0 for s in spin2])
check("the even half: 128 states, the matter roots of the recorded chart", int(matter.sum()) == 128)
Y6 = np.array([-2, -2, -2, 3, 3, 0, 0, 0])                   # 6 * hypercharge of the record
T6 = np.array([0, 0, 0, 3, -3, 0, 0, 0])
F6 = np.array([0, 0, 0, 0, 0, 6, -6, 0])


def marker(h6):
    return np.diag(spin2@h6)                                 # 12 * (h . s)


def current(i, si, j, sj):
    return (ann[i].T if si > 0 else ann[i])@(ann[j].T if sj > 0 else ann[j])


def marks(keep, X, h6, mode):
    Pk = np.diag(keep.astype(np.int64))
    mid = Pk if mode == "comp" else np.eye(dim, dtype=np.int64)
    H = marker(h6)
    z0 = int(np.trace(Pk))
    zu = Fr(int(np.trace(Pk@H@H)), 144*z0)
    zy = Fr(int(np.trace(Pk@X.T@mid@X)), z0)
    zuy = Fr(int(np.trace(Pk@H@X.T@mid@H@X)), 144*z0)
    return zy, zu, zuy


currents = [(i, si, j, sj) for i in range(5) for j in range(5, 8) for si in (1, -1) for sj in (1, -1)]
full = {c: marks(matter, current(*c), Y6, "cond") for c in currents}
check("60 pair currents of type (10,6); Z_Y/Z_0 = 1/4 for each", len(currents) == 60 and {v[0] for v in full.values()} == {Fr(1, 4)})
check("Z_U/Z_0 = |h|^2/4: hypercharge 5/24, T3 1/8",
      {v[1] for v in full.values()} == {Fr(5, 24)} and marks(matter, current(*currents[0]), T6, "cond")[1] == Fr(1, 8))
weak = {full[c][2] for c in currents if c[0] >= 3}
colour = {full[c][2] for c in currents if c[0] < 3}
check("Z_UY/Z_0 = 1/48 for weak and 11/288 for colour currents", weak == {Fr(1, 48)} and colour == {Fr(11, 288)})
conn = sorted({v[2]-v[1]*v[0] for v in full.values()})
check("connected mixed answer -1/32 (weak) and -1/72 (colour), strictly negative", conn == [Fr(-1, 32), Fr(-1, 72)], conn)

fam = [tuple(s[5:]) for s in spin2]
carrier_even = [list(s[:5]).count(-1) % 2 == 0 for s in spin2]
fam4 = sorted({f for f, m, e in zip(fam, matter, carrier_even) if m and e})
table = {}
for w in fam4:
    wn = tuple(-v for v in w)
    keep = matter & np.array([f not in (w, wn) for f in fam])
    assert int(keep.sum()) == 96
    for c in currents:
        X = current(*c)
        cd, cp = marks(keep, X, Y6, "cond"), marks(keep, X, Y6, "comp")
        f1 = marks(keep, X, F6, "cond")[1]
        f1full = marks(matter, X, F6, "cond")[1]
        key = " | ".join(str(v) for v in (cd[0], cp[0], cd[1], cd[2], cp[2], full[c][2], f1, f1full))
        table[key] = table.get(key, 0)+1
root = Path(sys.argv[1]).resolve()
recorded = json.loads((root/"experiments/theory-contracts/origin-diagnosis-new-hypotheses-20261116/results_word.json").read_text())["B_table"]
print("      three families (conditional | compressed | ...):")
for k, n in sorted(table.items()):
    print("        ", k, ":", n)
check("three families: the table equals the recorded one entry by entry (compression 1/6, 1/72, 11/432; "
      "family marker 1/3 or 2/3)", table == recorded)

# ---------------------------------------------------------------- Q4
print("Q4  place factors")
check("weight one at both ends of a sign line from 1 to -1: |1 - (-1)|^(-2) = 1/4; three places 1, i, -1: "
      "1/|(1 - i)(i + 1)(1 + 1)| = 1/4; the recorded magnitudes of the two-point and three-point kernel",
      abs(abs(1-(-1))**-2-0.25) < 1e-15 and abs(1/abs((1-1j)*(1j+1)*2)-0.25) < 1e-15)


def sign_line(Nsites, l):
    d = np.arange(Nsites)
    row = np.where(d == 0, 0.5, np.sin(np.pi*d/2)/(Nsites*np.sin(np.pi*np.where(d == 0, 1, d)/Nsites)))
    idx = np.arange(l)
    C = row[np.abs(idx[:, None]-idx[None, :])]
    return float(np.linalg.det(np.eye(l)-2*C))


vals = [sign_line(n, n//4)/sign_line(n, n//2) for n in (800, 1600)]
ext = 2*vals[1]-vals[0]
check("lattice ring: <(-1)^N> of a quarter arc over that of a half circle = 2^(1/4) for one lattice fermion, "
      "i.e. weight 1/8 at each end for each chiral Dirac fermion, weight one for eight",
      abs(ext-2**0.25) < 5e-3, (round(ext, 5), round(2**0.25, 5)))

print()
if FAIL:
    print("FAILED:", FAIL)
    sys.exit(1)
print("all checks passed")
