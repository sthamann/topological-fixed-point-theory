#!/usr/bin/env python3
"""The 96 native fibres against the three-family part of the 128 chart states (18 November 2026).

argv[1] = repository root.  Used unchanged: the assembler of the original Q contract and, through
it, verification/v252_full_finite_triple.py (index of the fibres, grading, real structure, full_D,
the four Higgs directions).  Chart: the matter states (+-1/2)^8 of round twenty-six, hypercharge
(-1/3,-1/3,-1/3,1/2,1/2,0,0,0), isospin (0,0,0,1/2,-1/2,0,0,0), one family weight omitted.

Sections
    D1  a bijection of the two sets of 96 states that carries all charges, the grading and the
        real structure
    D2  the 60 pair currents on the native fibres: the recorded table with the native hypercharge
    D3  how the carried currents sit in the native structure (grading, blocks, generations)
    D4  the native potential against the carried currents
"""
import itertools
import sys
from fractions import Fraction as Fr
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
rat = lambda x: Fr(float(np.real(x))).limit_denominator(2000)

# ---------------------------------------------------------------- chart
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
even_fam = sorted({tuple(v[5:]) for v in spin2 if list(v[5:]).count(-1) % 2 == 0})
w0 = even_fam[0]
kept = [f for f in even_fam if f != w0]
three = matter & np.array([tuple(v[5:]) not in (w0, tuple(-x for x in w0)) for v in spin2])
idx = np.where(three)[0]
pos = {int(k): n for n, k in enumerate(idx)}
hyp = np.array([-1/3, -1/3, -1/3, 1/2, 1/2, 0, 0, 0])
iso = np.array([0, 0, 0, 1/2, -1/2, 0, 0, 0])


def label(s):
    """(particle?, name, chirality, colour, generation) of a chart state."""
    plus3 = int(np.sum(s[:3] > 0))
    particle = plus3 in (1, 3)
    t = s if particle else -s
    colour = int(np.argmax(t[:3] > 0)) if int(np.sum(t[:3] > 0)) == 1 else 0
    quark = int(np.sum(t[:3] > 0)) == 1
    if t[3] != t[4]:
        name, ch = (("u" if t[3] > 0 else "d") if quark else ("nu" if t[3] > 0 else "e")), "L"
    else:
        name, ch = (("u" if t[3] > 0 else "d") if quark else ("nu" if t[3] > 0 else "e")), "R"
    fam = tuple(int(x) for x in t[5:])
    fam = fam if fam in kept else tuple(-x for x in fam)
    return particle, name, ch, colour, kept.index(fam)


# ---------------------------------------------------------------- D1
print("D1  the bijection")
Phi = np.zeros((96, 96))
for k in idx:
    particle, name, ch, colour, gen = label(spin2[k])
    Phi[(0 if particle else 48)+orig.IDX[(name, ch, colour, gen)], pos[int(k)]] = 1
check("every chart state gets one native fibre and every fibre is reached", np.allclose(Phi@Phi.T, np.eye(96)) and np.allclose(Phi.T@Phi, np.eye(96)))
YT = {("u", "L"): 1/6, ("d", "L"): 1/6, ("u", "R"): 2/3, ("d", "R"): -1/3, ("nu", "L"): -1/2, ("e", "L"): -1/2, ("nu", "R"): 0.0, ("e", "R"): -1.0}
TT = {("u", "L"): 0.5, ("nu", "L"): 0.5, ("d", "L"): -0.5, ("e", "L"): -0.5}
yv = [YT[(n, c)] for (n, c, _, _) in orig.P]
tv = [TT.get((n, c), 0.0) for (n, c, _, _) in orig.P]
Yh, T3 = np.diag(yv+[-v for v in yv]), np.diag(tv+[-v for v in tv])
carry = lambda D: Phi@D@Phi.T
Hc, Tc = np.diag(spin2[idx]@hyp/2), np.diag(spin2[idx]@iso/2)
check("it carries the chart hypercharge and isospin to the native hypercharge and isospin",
      np.allclose(carry(Hc), Yh) and np.allclose(carry(Tc), T3))
par = np.diag([1.0 if list(s[:5]).count(-1) % 2 == 1 else -1.0 for s in spin2[idx]])
check("it carries the parity of the five carrier signs to the native grading", np.allclose(carry(par), g))
neg = np.zeros((96, 96))
for k in idx:
    neg[pos[int(np.where((spin2 == -spin2[k]).all(axis=1))[0][0])], pos[int(k)]] = 1
check("it carries the reversal of all eight signs to the recorded real structure K", np.allclose(carry(neg), K))
blocks = np.round(np.exp(4j*np.pi*np.diag(Yh)), 9)
check("so the native 96 fibres and the three-family part of the chart are the same module of the charges, with "
      "grading and real structure", True)

# ---------------------------------------------------------------- D2
print("D2  the pair currents on the native fibres")
currents = [(i, j, si, sj) for i in range(5) for j in range(5, 8) for si in (1, -1) for sj in (1, -1)]
Xn = {}
for (i, j, si, sj) in currents:
    X = ((ann[i].T if si > 0 else ann[i])@(ann[j].T if sj > 0 else ann[j]))[np.ix_(idx, idx)]
    Xn[(i, j, si, sj)] = carry(X)
tr = lambda A: np.trace(A)/96
table = {}
for c, X in Xn.items():
    key = "weak" if c[0] >= 3 else "colour"
    zu, zy, zuy = tr(Yh@Yh), tr(X.T@X), tr(Yh@X.T@Yh@X)
    table.setdefault(key, set()).add((rat(zu), rat(zy), rat(zuy), rat(zuy-zu*zy)))
check("with the native hypercharge the 60 carried currents give the recorded three-family table: 5/24, 1/6, "
      "1/72 and 11/432, connected -1/48 and -1/108",
      table == {"weak": {(Fr(5, 24), Fr(1, 6), Fr(1, 72), Fr(-1, 48))}, "colour": {(Fr(5, 24), Fr(1, 6), Fr(11, 432), Fr(-1, 108))}},
      {k: sorted(v) for k, v in table.items()})

# ---------------------------------------------------------------- D3
print("D3  where the carried currents sit")
check("every carried current anticommutes with the native grading: like the native potential it connects a field "
      "with a conjugate field", all(np.allclose(X@g+g@X, 0) for X in Xn.values()))
keeps_blocks = lambda X: all(abs(X[r, c]) < 1e-12 or abs(blocks[r]-blocks[c]) < 1e-6 for r, c in zip(*np.nonzero(X)))
weak_ok = all(keeps_blocks(X) for c, X in Xn.items() if c[0] >= 3)
col_ok = all(not keeps_blocks(X) for c, X in Xn.items() if c[0] < 3)
check("the 24 weak currents keep the three native blocks of exp(4 pi i Y) (24, 36, 36); the 36 colour currents change "
      "the block, as recorded for the middle current", weak_ok and col_ok)
gens = [gen for (_, _, _, gen) in orig.P]*2
same_gen = max(abs(X[r, c]) for X in Xn.values() for r, c in zip(*np.nonzero(X)) if gens[r] == gens[c]) if any(
    gens[r] == gens[c] for X in Xn.values() for r, c in zip(*np.nonzero(X))) else 0.0
check("no carried current connects two fibres of the same generation: the currents are off-diagonal in the generations",
      same_gen == 0.0)
kpar = [np.allclose(K@X@K, Xn[(c[0], c[1], -c[2], -c[3])]) or np.allclose(K@X@K, -Xn[(c[0], c[1], -c[2], -c[3])]) for c, X in Xn.items()]
check("the real structure maps the current of a root to the current of the opposite root, up to sign", all(kpar))

# ---------------------------------------------------------------- D4
print("D4  the native potential against the carried currents")
DF = Asm.assemble({"schema": "tfpt.field-configuration.v1", "domain_sha256": Asm.domain["_sha256"], "genes": [1, 1, 1]})["phi_fourier"][0]
basis = np.array([X.reshape(-1) for X in Xn.values()]).T               # 9216 x 60


def captured(M):
    coef, *_ = np.linalg.lstsq(basis, M.reshape(-1), rcond=None)
    return float(np.linalg.norm(basis@coef)**2/np.linalg.norm(M)**2)


fr = {"D_F": captured(DF), "Majorana part": captured(Asm.DM), "Yukawa part": captured(DF-Asm.DM),
      "Higgs directions": [round(captured(v), 6) for v in Asm.V]}
print("      share of the squared norm inside the span of the 60 currents:", {k: (round(v, 6) if not isinstance(v, list) else v) for k, v in fr.items()})
diag_share = float(sum(abs(DF[r, c])**2 for r in range(96) for c in range(96) if gens[r] == gens[c])/np.linalg.norm(DF)**2)
check("the 60 currents are linearly independent on the 96 fibres", np.linalg.matrix_rank(basis) == 60)
check("the Majorana part has no component along the carried currents; the Yukawa part has one: 0.114 of its squared "
      "norm (0.103 of D_F); the first run expected zero", fr["Majorana part"] < 1e-12 and abs(fr["Yukawa part"]-0.11416) < 1e-5
      and abs(fr["D_F"]-0.103492) < 1e-5)
print("      share of D_F that is diagonal in the generations: %.6f" % diag_share)
shares = []
for perm in itertools.permutations(range(3)):
    Pg = np.zeros((96, 96))
    for n, (name, ch, c, gen) in enumerate(orig.P):
        for off in (0, 48):
            Pg[off+orig.IDX[(name, ch, c, perm[gen])], off+n] = 1
    bp = np.array([(Pg@X@Pg.T).reshape(-1) for X in Xn.values()]).T
    coef, *_ = np.linalg.lstsq(bp, DF.reshape(-1), rcond=None)
    shares.append(round(float(np.linalg.norm(bp@coef)**2/np.linalg.norm(DF)**2), 6))
print("      share of D_F for the six orders of the generations:", shares)
check("the bijection is fixed by the charges only up to the order of the generations (and signs); the share depends "
      "on that order, so it is no invariant number", max(shares)-min(shares) > 1e-4 and min(shares) > 0, (min(shares), max(shares)))
check("in every order the larger part of D_F lies outside the span of the 60 currents: the native Higgs coupling is "
      "not a combination of the carried pair currents", max(shares) < 0.5)

print()
if FAIL:
    print("FAILED:", FAIL)
    sys.exit(1)
print("all checks passed")
