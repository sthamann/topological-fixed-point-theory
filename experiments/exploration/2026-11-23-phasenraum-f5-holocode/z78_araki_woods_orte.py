"""Z78: Vakuum gegen thermisch und Orte aus Korrelationen -- Prüfung der sechsten externen Analyse (24.11.).

Teil 1  Araki-Woods: Aus den nativen Stromoszillatoren a_i (zwei Quartersektoren, CCR aus z58) und ihren Tomita-Partnern b_i
        sind d_i = (a_i - q b_i^+)/sqrt(1-q^2), e_i = (b_i - q a_i^+)/sqrt(1-q^2) kanonisch (CCR) und vernichten die Purifikation
        des thermischen Originalzustands (q = e^-2pi). Das 'Vakuum' der Zeit- und Raumstrukturen ist damit die Purifikation des
        Originals selbst: Die Temperatur muss nicht auf null gesetzt werden.
Teil 2  c_i = d_i + e_i^+ = kappa (a_i + b_i^+), kappa^2 = (1-q)/(1+q), vertauscht mit allen c_j^+; P_mu = (1/2) c^+ sigma_mu c
        vertauschen, P_0 >= 0, P^2 = 0 (Rang eins).
Teil 3  Orte: Für |X> = exp(i P(X)) |Omega> gilt exakt (M_X - M_X^+)/i = X mit (M_X)_ij = <X| a_i b_j |X>, unabhängig von q;
        die linke Algebra allein sieht nur X^2 (unterscheidet X und -X nicht).
"""
import json
from pathlib import Path

import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spla
import sympy as sp

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val):
    BEFUNDE[name] = bool(val)
    say(f"   [Befund: {'ja' if val else 'nein'}] {name}")


# ------------------------------------------------------------------ Teil 1/2 (exakt)
say("Teil 1/2: Araki-Woods-Oszillatoren und Translationen (exakt)")
q = sp.symbols("q", positive=True)
s = sp.sqrt(1 - q**2)
nbar = q**2/(1 - q**2)                     # <a^+ a> in der Purifikation sqrt(1-q^2) sum q^n |n,n>
ab = q/(1 - q**2)                          # <a b>
kap2 = (1 - q)/(1 + q)
# Kommutatoren linearer Kombinationen x = u1 a + u2 a^+ + v1 b + v2 b^+: [x, y] = u1 y2 - u2 y1 + v1 w2 - v2 w1
def comm(x, y):
    return sp.simplify(x[0]*y[1] - x[1]*y[0] + x[2]*y[3] - x[3]*y[2])
d = (1/s, 0, 0, -q/s)                      # (a, a^+, b, b^+)
e = (0, -q/s, 1/s, 0)
dag = lambda x: (x[1], x[0], x[3], x[2])
ccr = comm(d, dag(d)) == 1 and comm(e, dag(e)) == 1 and comm(d, e) == 0 and comm(d, dag(e)) == 0
# Vernichtung der Purifikation: a|Omega> = q b^+|Omega> und b|Omega> = q a^+|Omega> (TFD), also d|Omega> = e|Omega> = 0
nmax = 40
qq = 0.3
tfd = np.zeros((nmax, nmax))
for n in range(nmax):
    tfd[n, n] = np.sqrt(1 - qq**2)*qq**n
a1 = np.diag(np.sqrt(np.arange(1, nmax)), 1)
I = np.eye(nmax)
A, Bm = np.kron(a1, I), np.kron(I, a1)
v = tfd.reshape(-1)
dv = (A - qq*Bm.T) @ v
ev = (Bm - qq*A.T) @ v
keep = [i for i in range(nmax*nmax) if i//nmax < nmax - 3 and i % nmax < nmax - 3]
annih = np.max(np.abs(dv[keep])) < 1e-12 and np.max(np.abs(ev[keep])) < 1e-12
check("d_i, e_i erfüllen die CCR exakt (symbolisch in q) und vernichten die Purifikation des thermischen Zustands (numerisch, "
      "q = 0,3): das Oszillatorvakuum ist die Purifikation des Originals", ccr and annih)
c = tuple(sp.simplify(d[k] + dag(e)[k]) for k in range(4))
c_prop = all(abs(float((c[k] - sp.sqrt(kap2)).subs(q, qv))) < 1e-14 for k in (0, 3) for qv in (0.3, float(sp.exp(-2*sp.pi)))) and c[1] == 0 and c[2] == 0
check("c_i = d_i + e_i^+ = kappa (a_i + b_i^+) mit kappa^2 = (1-q)/(1+q), und [c_i, c_i^+] = 0: c und c^+ vertauschen, also "
      "vertauschen alle P_mu = (1/2) c^+ sigma_mu c, P_0 = (1/2) c^+ c >= 0 und P^2 = det(c c^+) = 0 (Rang eins)",
      c_prop and comm(c, dag(c)) == 0)

# ------------------------------------------------------------------ Teil 3 (exakt)
say("\nTeil 3: Orte aus der gemischten Korrelation")
# Heisenberg: e^{-iP} a e^{iP} = a + i(kappa^2/2)(X c~), e^{-iP} b e^{iP} = b + i(kappa^2/2)(c~^+ X), c~ = a + b^+ (vertauscht mit P)
# (M - M^+)/i = kappa^2 X (1 + 2 nbar + 2 <ab>)
fac = sp.simplify(kap2*(1 + 2*nbar + 2*ab))
check("(M_X - M_X^+)/i = kappa^2 (1 + 2 nbar + 2<ab>) X und der Faktor ist exakt 1 für jedes q: die vier Koordinaten X in Herm(2) "
      "sind aus der gemischten Korrelation <a_i b_j> exakt und temperaturunabhängig auslesbar", sp.simplify(fac - 1) == 0)


def fock4(n):
    a = sps.diags(np.sqrt(np.arange(1, n)), 1, format="csr")
    Id = sps.identity(n, format="csr")
    ops = []
    for k in range(4):
        mats = [Id]*4
        mats[k] = a
        op = mats[0]
        for M in mats[1:]:
            op = sps.kron(op, M, format="csr")
        ops.append(op)
    return ops                              # a1, a2, b1, b2


sig = [np.eye(2), np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.diag([1, -1])]
rng = np.random.default_rng(78)
ok_num, ok_left = True, True
n = 9
a1_, a2_, b1_, b2_ = fock4(n)
aa, bb = [a1_, a2_], [b1_, b2_]
for qv in (0.2, float(np.exp(-2*np.pi))):
    om = np.zeros([n]*4)
    for i in range(n):
        for j in range(n):
            om[i, j, i, j] = (1 - qv**2)*qv**(i + j)
    om = om.reshape(-1).astype(complex)
    kap = np.sqrt((1 - qv)/(1 + qv))
    ct = [kap*(aa[i] + bb[i].T) for i in range(2)]
    for _ in range(3):
        xv = rng.normal(size=4)*0.15
        X = sum(xv[k]*sig[k] for k in range(4))
        P = 0.5*sum(ct[i].conj().T @ ct[j]*X[i, j] for i in range(2) for j in range(2))
        st = spla.expm_multiply(1j*P, om)
        stm = spla.expm_multiply(-1j*P, om)
        M = np.array([[np.vdot(st, aa[i] @ (bb[j] @ st)) for j in range(2)] for i in range(2)])
        Xrec = (M - M.conj().T)/1j
        ok_num &= np.allclose(Xrec, X, atol=3e-3)
        L = np.array([[np.vdot(st, aa[i].T @ (aa[j] @ st)) for j in range(2)] for i in range(2)])
        Lm = np.array([[np.vdot(stm, aa[i].T @ (aa[j] @ stm)) for j in range(2)] for i in range(2)])
        Mm = np.array([[np.vdot(stm, aa[i] @ (bb[j] @ stm)) for j in range(2)] for i in range(2)])
        ok_left &= np.allclose(L, Lm, atol=1e-6) and not np.allclose(M, Mm, atol=1e-4)
check("numerisch auf dem Fockraum (vier Moden, je 9 Stufen; q = 0,2 und e^-2pi; drei zufällige X): (M_X - M_X^+)/i = X",
      ok_num)
check("die linke Algebra allein sieht |X> und |-X> gleich (<a_i^+ a_j> gleich), die gemischte Korrelation unterscheidet sie: "
      "die Ortsinformation steckt in der Beziehung zum Tomita-Partner", ok_left)
befund("Vakuum gegen thermisch (H6h35 (c)) ist strukturell beantwortet: Die positiven Translationen und ihr Vakuum leben in der "
       "Araki-Woods-Darstellung, deren Vakuum die Purifikation des thermischen Originals ist; die linke thermische Algebra allein "
       "hat keine halbseitigen Inklusionen (z72). I_pos hat damit einen konkreten Kandidaten: Lagen = relative Zustände "
       "exp(i P(X)) Omega, ausgelesen über die gemischte Korrelation mit dem Partner. Offen bleibt, dass die Oszillatoren a_1, a_2 "
       "der beiden Quarter das Dublett der Naht sind und die Lorentzwirkung auf diesem Index nativ ist (Lesart)", True)

urteil = ("Die sechste externe Analyse hält in den geprüften Punkten. Die Araki-Woods-Oszillatoren d_i, e_i aus den nativen "
          "Stromoszillatoren und ihren Tomita-Partnern sind kanonisch und haben die Purifikation des thermischen Originals als "
          "Vakuum. c_i = d_i + e_i^+ = kappa (a_i + b_i^+) vertauscht mit c^+, daher vertauschen die P_mu = (1/2) c^+ sigma_mu c, "
          "P_0 >= 0 und P^2 = 0. Für |X> = exp(i P(X)) Omega gilt exakt (M_X - M_X^+)/i = X für jedes q; die linke Algebra "
          "allein sieht nur X^2. Damit ist 'Vakuum gegen thermisch' strukturell beantwortet (das Vakuum ist die Purifikation) und "
          "I_pos hat einen konkreten Kandidaten; offen ist, dass der Quarterindex das Dublett der Naht trägt.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z78_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
