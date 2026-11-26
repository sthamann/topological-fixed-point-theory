"""Z160: Kopplung zweier nativer lemniskatischer Ketten: Entsteht das wechselwirkende Gatter i^(Q_A Q_B) mit Viertelwinkel?

Leitlinie 'Rekonstruktion statt Auswahl'; H6 und die gelösten Blöcke (Flavour, Massen, alpha, 16, E8) bleiben unberührt,
nichts wird befördert. Statuskürzel: E exakt, C bedingt (Prämisse genannt), L Lesart bzw. Literatur, N numerisch, O offen.

Anlass (z157, VERDICT Block B4): (i) Die Nahtkette aus z152 (Hüpfen t, Wechselpotential m = 2t) liegt wörtlich auf der
Markenkurve E_i: w^2 = 1 - s^4; ihre R-Matrix ist der Grenzfall reiner Paarung der Bazhanov-Stroganov-Familie mit
k^2 = -1, frei in den Fermionen, mit U(1)-Ladung; am Uhrpunkt ist sie die Paaruhr diag(-i, i). (ii) In der symmetrischen
elliptischen Klasse (Acht-Vertex) gilt auf E_i für jeden Kreuzungsparameter lam am Punkt u = i lam: unitär, delta_y = -i.
Die Familie hat die Enden lam -> 0 (z120-Gatter (1 + iP)/(1 + i)) und lam = K (freie Nahtkette); Prämissen P_aut, P_lam.
(iii) Das Gatter des Kandidaten ist in zirkularen Moden die gesteuerte Uhr, i^(Q_A Q_B) = 'Ladung um Fluss' (z120, z145).

Frage: Entsteht das wechselwirkende Gatter (gesteuerte Uhr zwischen zwei Ladungen, Viertelwinkel) aus der Kopplung zweier
nativer lemniskatischer Ketten, mit einem durch E_i festgelegten Kopplungswert, und ist der Limes lam -> 0 nativ
ausgezeichnet?

Literatur [L] (gelesen über arXiv-Zusammenfassung bzw. ar5iv; Formeln daraus sind unten nachgerechnet, nicht übernommen):
  Shastry, Phys. Rev. Lett. 56 (1986) 1529 und 2453; J. Stat. Phys. 50 (1988) 57: Hubbard-Modell als zwei gekoppelte
    XX-Ketten, L = I(h) L_auf L_ab I(h), Bedingung sinh(2h)/sin(2 theta) = U/4 (aus dem Gedächtnis zitiert; unten in B2
    numerisch reproduziert). Shiroishi/Wadati, J. Phys. Soc. Jpn. 64 (1995) 57 (Yang-Baxter-Gleichung für Shastrys
    R-Matrix) und 64 (1995) 2795 ('Decorated star-triangle relations for the free-fermion model and a new solvable
    bilayer vertex model'; nur als Titel bekannt, nicht geprüft). Umeno/Shiroishi/Wadati, J. Phys. Soc. Jpn. 67 (1998)
    1930 und 2242 (fermionischer R-Operator).
  Melikyan, arXiv:1912.10786 (über ar5iv gelesen): fermionischer Bazhanov-Stroganov-R-Operator; die dekorierte
    Yang-Baxter-Gleichung gilt für die ganze elliptische Familie (Gl. 3.1: (2n-1) R (2n-1) = R mit um 2K verschobenen
    Farben und k -> -k); die tetraedrische Zamolodchikov-Algebra ist dort nur für k = 0 (XX-Modell) gelöst, die
    gekoppelte R-Matrix für den elliptischen Fall ausdrücklich offen ('future work').
  De Marco/Tolle/Halati/Sheikhan/Läuchli/Kollath, arXiv:2203.05379 (Zusammenfassung): Das ionische Hubbard-Modell
    (Hubbard-Kette mit Wechselpotential) hat Niveaustatistik einer Gaußschen Zufallsmatrix, das Hubbard-Modell Poisson.
  Schulz/Shastry, Phys. Rev. Lett. 80 (1998) 1924, cond-mat/9707135 (aus dem Gedächtnis): Fermionen mit korreliertem
    Hüpfen, lösbar durch eine nichtlokale unitäre Transformation; hier in B6 als abelscher Twist nachgerechnet.
  Reshetikhin, Lett. Math. Phys. 20 (1990) 331 (abelscher Twist einer R-Matrix; aus dem Gedächtnis).
  Baxter, Exactly Solved Models (1982): Kap. 10.4, 10.14 (Acht-Vertex, XYZ), Kap. 13 (Eckentransfermatrix: Dichtematrix der
    Halbkette = Produkt von diag(1, x^(2j))), Abschnitt 10.12 und Johnson/Krinsky/McCoy, Phys. Rev. A 8 (1973) 2526
    (Korrelationslänge; Normierung unten am freien Punkt geeicht). Kadanoff/Wegner, Phys. Rev. B 4 (1971) 3989 und
    Wu, Phys. Rev. B 4 (1971) 2312: Acht-Vertex-Modell = zwei Ising-Modelle mit Vierspinkopplung. Kohmoto/den Nijs/
    Kadanoff, Phys. Rev. B 24 (1981) 5229 (Ashkin-Teller-Kette). Destri/de Vega, Nucl. Phys. B 290 (1987) 363
    (Lichtkegelgitter); Vanicat/Zadnik/Prosen, arXiv:1712.00431 (integrable Trotterisierung).
  Thetafunktionen: Jacobi-Transformation theta_1(v/tau | -1/tau) = -i sqrt(-i tau) e^(i pi v^2/tau) theta_1(v | tau);
    Stufe-2-Thetafunktionen und Weil-Darstellung (Kac/Peterson; unten numerisch geprüft).

Teil A  Die native R-Matrix in Normalform auf der Lemniskatenkurve x^2 + y^2 + x^2 y^2 = 1; Symmetrien, Torsionspunkte;
        neu: Die Verschiebung um den reellen 2-Torsionspunkt ist die Uhr-Eichung.
Teil B  Kopplung zweier Ketten: Shastry-Form (Kontrolle XX, nativ), lokale Erhaltungsgröße der Reichweite 3, Halbketten-
        leiter, abelscher Twist = Torsionstranslation, Ising-Form = Acht-Vertex-Familie.
Teil C  Der Limes lam -> 0: geschlossene Form der Familie u = i lam, Uhrwerte, stetige Symmetrie, Masse, Leiter,
        Lichtkegelgitter.
Teil D  Die zwei Lesarten des Viertels (z155 Villain/Theta, z157 komplexe Multiplikation) sind eine.
Teil E  Bilanz und Urteil.
"""
import itertools
import json
import os
import time
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")                     # deterministische Rundung (einfädig)

import mpmath as mp
import numpy as np
import sympy as sp
from scipy.linalg import eigvalsh as eigvalsh_teil
from scipy.optimize import minimize
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import LinearOperator, eigsh

mp.mp.dps = 50
HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
T0 = time.time()
TOL = mp.mpf(10) ** -40


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val, text):
    BEFUNDE[name] = val
    say(f"   [Befund: {val}] {text}")


def ns(x, n=10):
    return mp.nstr(x, n)


def fl(x, n=4):
    return f"{round(float(x), n) + 0.0:.{n}f}"


# --------------------------------------------------------------------------------------------------------------------
# Werkzeuge (mpmath)
PI = mp.pi
I = mp.mpc(0, 1)


def kron(A, B):
    n, m = A.rows, B.rows
    M = mp.zeros(n * m, n * m)
    for i in range(n):
        for j in range(n):
            if A[i, j] == 0:
                continue
            for k in range(m):
                for l in range(m):
                    M[i * m + k, j * m + l] = A[i, j] * B[k, l]
    return M


def nrm(A):
    return mp.sqrt(mp.fsum(abs(A[i, j]) ** 2 for i in range(A.rows) for j in range(A.cols)))


E2 = mp.eye(2)
E4 = mp.eye(4)
SX = mp.matrix([[0, 1], [1, 0]])
SY = mp.matrix([[0, -I], [I, 0]])
SZ = mp.matrix([[1, 0], [0, -1]])
XX, YY, ZZ = kron(SX, SX), kron(SY, SY), kron(SZ, SZ)
PERM = mp.matrix([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]])        # Tausch P
P23 = kron(E2, PERM)


def e12(R):
    return kron(R, E2)


def e23(R):
    return kron(E2, R)


def e13(R):
    return P23 * e12(R) * P23


def sn(x, k):
    return mp.ellipfun("sn", x, k * k)


def cn(x, k):
    return mp.ellipfun("cn", x, k * k)


def dn(x, k):
    return mp.ellipfun("dn", x, k * k)


def sl(x):
    """Lemniskatischer Sinus: Umkehrung von Integral dz/sqrt(1 - z^4) (Jacobi-sn zum Modul k = i)."""
    return mp.ellipfun("sn", x, -1)


def cl(x):
    """Lemniskatischer Kosinus cl = cn/dn zum Modul k = i; es gilt sl^2 + cl^2 + sl^2 cl^2 = 1."""
    return mp.ellipfun("cn", x, -1) / mp.ellipfun("dn", x, -1)


VARPI = mp.gamma(mp.mpf(1) / 4) ** 2 / (2 * mp.sqrt(2 * PI))            # Lemniskatenkonstante
KB = (mp.sqrt(2) - 1) ** 2                                              # Baxters Modul der Ising-Form (z157)
K_ = mp.ellipk(KB ** 2)


def wQ(u, k=I):
    """Gewichte der Ladungskette aus z157 A15 (Grenzfall reiner Paarung der Bazhanov-Stroganov-Familie)."""
    eu = cn(u, k) + I * sn(u, k)
    return eu, (I / 2) * (1 - eu) / sn(u / 2, k), (k / 2) * (1 + eu) * sn(u / 2, k)


def RcQ(u, k=I):
    eu, c, d = wQ(u, k)
    return mp.matrix([[1, 0, 0, d], [0, c, 0, 0], [0, 0, c, 0], [d, 0, 0, eu]])


def RQ(u, k=I):
    return PERM * RcQ(u, k)


def Bn(u):
    """Normalform des Ein-Teilchen-Blocks auf der Lemniskatenkurve: s = sl(u/2), c = cl(u/2)."""
    s, c = sl(u / 2), cl(u / 2)
    return mp.matrix([[c - I * s, I * s * c], [I * s * c, c + I * s]])


def Gn(u):
    """Gatter der Ladungskette in Normalform (Paarungsbild): gerader Block Bn(u), ungerader Block Eins."""
    b = Bn(u)
    return mp.matrix([[b[0, 0], 0, 0, b[0, 1]], [0, 1, 0, 0], [0, 0, 1, 0], [b[1, 0], 0, 0, b[1, 1]]])


def npm(M):
    return np.array([[complex(M[i, j]) for j in range(M.cols)] for i in range(M.rows)])


NX = np.array([[0, 1], [1, 0]], complex)
NY = np.array([[0, -1j], [1j, 0]])
NZ = np.diag([1.0, -1.0]).astype(complex)
NI = np.eye(2, dtype=complex)
NP4 = np.array([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]], complex)


def op(o, j, L):
    m = np.array([[1]], complex)
    for i in range(L):
        m = np.kron(m, o if i == j else NI)
    return m


def emb(G, pos, n):
    """Operator G auf den Qubits pos (Liste, in dieser Reihenfolge) in n Qubits; Qubit 0 ist das höchstwertige."""
    k = len(pos)
    T = np.asarray(G, complex).reshape([2] * (2 * k))
    full = np.eye(2 ** n, dtype=complex).reshape([2] * (2 * n))
    res = np.tensordot(T, full, axes=(list(range(k, 2 * k)), list(pos)))
    res = np.moveaxis(res, list(range(k)), list(pos))
    return res.reshape(2 ** n, 2 ** n)


def rel(A, B):
    return float(np.linalg.norm(A - B) / max(np.linalg.norm(A), 1e-300))


# --------------------------------------------------------------------------------------------------------------------
say("Z160: Kopplung zweier lemniskatischer Ketten, das Gatter i^(Q_A Q_B) und der Limes lam -> 0")
say("=" * 100)
say("\nTeil A  Die native R-Matrix der Ladungskette, exakt")
pts = [(mp.mpc("0.37", "0.2"), mp.mpc("-0.21", "0.05")), (mp.mpf("1.3"), mp.mpf("0.45")),
       (mp.mpc("-0.8", "0.6"), mp.mpc("0.3", "-1.1"))]
ok = True
for uq, vq in pts:
    Lh = e12(RQ(uq - vq)) * e13(RQ(uq)) * e23(RQ(vq))
    Rh = e23(RQ(vq)) * e13(RQ(uq)) * e12(RQ(uq - vq))
    ok &= nrm(Lh - Rh) < TOL * nrm(Lh)
check("A1 [E] Yang-Baxter-Gleichung der R-Matrix der Ladungskette (z157 A15, k = i) in Differenzform an drei Punkten "
      "(zwei komplex), Rest < 1e-40", ok)
ok = True
for u in (mp.mpf("0.37"), mp.mpf("2.1"), mp.mpc("0.6", "0.8"), mp.mpc("-1.2", "0.3")):
    eu, c, d = wQ(u)
    ok &= nrm(RcQ(u) / c - Gn(u)) < TOL
    s, c_ = sl(u / 2), cl(u / 2)
    ok &= abs(s * s + c_ * c_ + s * s * c_ * c_ - 1) < TOL and abs(mp.det(Bn(u)) - 1) < TOL
xs, ys = sp.symbols("x y")
kurve = xs ** 2 + ys ** 2 + xs ** 2 * ys ** 2 - 1
aut = sp.simplify(kurve.subs({xs: sp.I * xs, ys: 1 / ys}, simultaneous=True) * ys ** 2 + kurve)
check("A2 [E] Normalform: Bis auf den Skalar c ist das Gatter im ungeraden Block die Eins und im geraden Block "
      "B(u) = [[cl - i sl, i sl cl], [i sl cl, cl + i sl]] mit sl = sl(u/2), cl = cl(u/2). det B = cl^2 + sl^2 + sl^2 cl^2 "
      "= 1 ist die Gleichung der Lemniskatenkurve x^2 + y^2 + x^2 y^2 = 1 (Edwards-Form von E_i); (x, y) -> (i x, 1/y) "
      "ist ihr Automorphismus der Ordnung 4", ok and aut == 0)
u0 = mp.mpf("0.83")
G0 = Gn(u0)
Zm, Zp = kron(SZ, E2) - kron(E2, SZ), kron(SZ, E2) + kron(E2, SZ)
ok = (nrm(G0 * Zm - Zm * G0) < TOL and nrm(G0 * ZZ - ZZ * G0) < TOL and nrm(PERM * G0 * PERM - G0) < TOL
      and nrm(G0.T - G0) < TOL and nrm(Gn(-u0) * G0 - E4) < TOL and nrm(G0 * G0.H - E4) < TOL
      and nrm(G0 * Zp - Zp * G0) > mp.mpf("0.1"))
check("A3 [E] Symmetrien: U(1) (im Paarungsbild Z_1 - Z_2, die gestaffelte Ladung), Parität Z_1 Z_2, Platztausch, "
      "Transposition, G(-u) = G(u)^(-1); für reelles u unitär. Leeres und doppelt besetztes Paar tragen die Eins und "
      "det B = 1: Das Gatter ist einkörperig (gaußsch)", ok)
ur = mp.mpf("0.9")
Bi = Bn(I * ur)
s, c_ = sl(ur / 2), cl(ur / 2)
ok = (abs(sl(I * ur / 2) - I * s) < TOL and abs(cl(I * ur / 2) - 1 / c_) < TOL
      and nrm(Bi - mp.matrix([[1 / c_ + s, -s / c_], [-s / c_, 1 / c_ - s]])) < TOL
      and nrm(Bi - Bi.H) < TOL and mp.re(Bi[0, 0]) > 0 and abs(mp.det(Bi) - 1) < TOL)
check("A4 [E] Z_4-Automorphie: u -> i u ist (sl, cl) -> (i sl, 1/cl). Sie führt das unitäre Gatter (reelle Zeit) in ein "
      "hermitesches positives mit Determinante 1 über (euklidische Richtung): B(i u) = [[1/cl + sl, -sl/cl], "
      "[-sl/cl, 1/cl - sl]]", ok)
B1, B2, Bh = Bn(VARPI), Bn(2 * VARPI), Bn(VARPI / 2)
Rkreuz = RcQ(I * VARPI)
s2 = sp.symbols("s2", positive=True)
pmax = sp.solve(sp.diff(s2 * (1 - s2) / (1 + s2), s2), s2)
pval = sp.simplify((s2 * (1 - s2) / (1 + s2)).subs(s2, pmax[0]) - (sp.sqrt(2) - 1) ** 2)
ok = (nrm(B1 - mp.diag([-I, I])) < TOL and nrm(B2 + E2) < TOL
      and abs(Bh[0, 1] - I * (mp.sqrt(2) - 1)) < TOL and abs(mp.arg(Bh[0, 0]) + PI / 4) < TOL
      and nrm(Rkreuz - mp.matrix([[1, 0, 0, -1], [0, 0, 0, 0], [0, 0, 0, 0], [-1, 0, 0, 1]])) < TOL
      and len(pmax) == 1 and pval == 0)
check("A5 [E] Torsionspunkte (Gitter 2 varpi Z[i]): u = varpi (reelle 2-Torsion): B = diag(-i, i), die Paaruhr; "
      "u = 2 varpi: B = -1 (Parität); u = i varpi (Kreuzung): Rang 1; u = varpi/2 (4-Torsion): Transferamplitude "
      "i (sqrt 2 - 1), Phase e^(-i pi/4). Die Transferwahrscheinlichkeit sl^2 cl^2 ist auf der reellen Achse höchstens "
      "(sqrt 2 - 1)^2 (exakt): Der symmetrische Strahlteiler P_sym ist kein Wert (z157 C15/C16 reproduziert)", ok)
S1 = mp.matrix([[1, 0], [0, I]])                                        # Uhr i^n eines Platzes
S2 = kron(E2, S1)
ok = True
for u in (mp.mpf("0.41"), mp.mpc("0.3", "0.7"), mp.mpc("-1.9", "0.2")):
    ok &= nrm(Gn(u + VARPI) - Gn(VARPI) * S2 * Gn(u) * S2 ** -1) < TOL
    ok &= nrm(Bn(u + 2 * I * VARPI) + SZ * Bn(u) * SZ) < TOL
    ok &= nrm(S2 ** -1 * Gn(u) * S2 - Gn(VARPI) * Gn(u - VARPI)) < TOL
ok &= nrm(S1 ** 4 - E2) < TOL and nrm(S1 ** 2 - SZ) < TOL
check("A6 [E, neu] Die Verschiebung um den reellen 2-Torsionspunkt ist die Uhr: G(u + varpi) = G(varpi) S G(u) S^(-1) "
      "mit S = i^n auf dem zweiten Platz (Ordnung 4, S^2 = sigma_z) und der Paaruhr G(varpi). Gleichwertig "
      "S^(-1) G(u) S = G(varpi) G(u - varpi). Die imaginäre Richtung gibt mit 2 i varpi nur S^2", ok)


def detBp(u):
    d = mp.diff(lambda x: Bn(x)[0, 0], u), mp.diff(lambda x: Bn(x)[0, 1], u), mp.diff(lambda x: Bn(x)[1, 1], u)
    return d[0] * d[2] - d[1] * d[1]


ok = True
for u in (mp.mpf("0.3"), mp.mpc("0.7", "0.4"), mp.mpc("-0.2", "1.1")):
    p = (sl(u / 2) * cl(u / 2)) ** 2
    ok &= abs(detBp(u) - (1 - p) / 2) < mp.mpf(10) ** -30
ok &= abs(detBp(mp.mpf(0)) - mp.mpf(1) / 2) < mp.mpf(10) ** -30
nullst = [abs(sl(a / 2) * cl(a / 2)) < TOL for a in (VARPI, 2 * VARPI, 3 * VARPI, 2 * I * VARPI, VARPI + 2 * I * VARPI)]
nicht = [abs(sl(a / 2) * cl(a / 2)) > mp.mpf("0.1") for a in (VARPI / 2, I * VARPI * (1 - mp.mpf(10) ** -20),
                                                               (1 + I) * VARPI / 2, mp.mpf("0.7"))]
check("A7 [E] Eindeutigkeit: Eine Translation u -> u + a wirkt nur dann als Eichung mal festem Faktor, wenn "
      "det B'(a) = det B'(0). Es ist det B'(u) = (1 - sl^2 cl^2)/2, also muss sl(a/2) cl(a/2) = 0 sein: "
      "a in varpi Z + 2 i varpi Z. Die zugehörigen Eichungen sind genau die Potenzen der Uhr S. Eine Eichung "
      "diag(1, e^(i phi)) ist also nur für phi in (pi/2) Z eine Translation auf der Kurve", ok and all(nullst) and all(nicht))
befund("A", "E",
       "Die native R-Matrix ist ein Punkt der Lemniskatenkurve x^2 + y^2 + x^2 y^2 = 1 in SU(2): "
       "B(u) = cl - i sl (sigma_z - cl sigma_x). Sie ist einkörperig, trägt U(1), Parität und die Z_4-Automorphie "
       "(reelle gegen euklidische Zeit). Neu: Die Uhr S = i^n ist die Eichung der Halbperiode varpi; andere "
       "Eichwinkel als Vielfache von pi/2 sind keine Translationen der Kurve.")

# --------------------------------------------------------------------------------------------------------------------
say("\nTeil B  Kopplung zweier Ketten")
say("  (a) Hubbard-Kopplung der Ladungsform = ionisches Hubbard-Modell")
L4 = 4
nq = 2 * L4                                                             # Qubits: A1..A4, B1..B4
sm_ = np.array([[0, 1], [0, 0]], complex)
cs = []
for j in range(nq):
    m_ = np.array([[1]], complex)
    for i in range(nq):
        m_ = np.kron(m_, NZ if i < j else (sm_ if i == j else NI))
    cs.append(m_)
Uh = 0.7
Hf = np.zeros((2 ** nq, 2 ** nq), complex)
for spz in (0, 1):
    o = spz * L4
    for j in range(L4 - 1):
        Hf += -(cs[o + j].conj().T @ cs[o + j + 1] + cs[o + j + 1].conj().T @ cs[o + j])
    for j in range(L4):
        Hf += -2.0 * (-1) ** j * cs[o + j].conj().T @ cs[o + j]
E256 = np.eye(2 ** nq)
for j in range(L4):
    Hf += Uh * (cs[j].conj().T @ cs[j] - 0.5 * E256) @ (cs[L4 + j].conj().T @ cs[L4 + j] - 0.5 * E256)
Uph = np.eye(2 ** nq, dtype=complex)
for j in range(1, L4, 2):
    Uph = Uph @ op(NX, j, nq) @ op(NX, L4 + j, nq)
Hp = np.zeros((2 ** nq, 2 ** nq), complex)
for spz in (0, 1):
    o = spz * L4
    for j in range(L4 - 1):
        Hp += -0.5 * (op(NX, o + j, nq) @ op(NX, o + j + 1, nq) - op(NY, o + j, nq) @ op(NY, o + j + 1, nq))
    for j in range(L4):
        Hp += op(NZ, o + j, nq)
for j in range(L4):
    Hp += Uh / 4 * op(NZ, j, nq) @ op(NZ, L4 + j, nq)
check("B1 [E] Zwei Ladungsketten (t = 1, m = -2) mit Dichtekopplung U (n_A - 1/2)(n_B - 1/2) sind nach Jordan-Wigner "
      "und Teilchen-Loch-Tausch auf den ungeraden Plätzen exakt zwei Paarungsketten aus z157 A13 mit der Kopplung "
      "(U/4) Z_A Z_B je Platz (vier Plätze, offen). Das ist das ionische Hubbard-Modell bei Wechselpotential = 2t. "
      "Z_A Z_B ist die einzige Platzkopplung, die beide Ladungen erhält",
      np.linalg.norm(Uph @ Hf @ Uph.conj().T - Hp) < 1e-10)

say("  (b) Shastry-Form L = I(h) L_A L_B I(h), I(h) = exp(h Z_A Z_B/2) im Hilfsraum")
ZZn = np.kron(NZ, NZ)


def r_xx(th):
    a, b = np.cos(th), np.sin(th)
    return np.array([[a, 0, 0, 0], [0, b, 1, 0], [0, 1, b, 0], [0, 0, 0, a]], complex)


def r_nat(u):
    return npm(RQ(mp.mpmathify(u)))


def L_shastry(R, u, h):
    """Indizes (Hilfsraum A, B; Quantenraum A, B), je 4 x 4."""
    Ra = R(u).reshape(2, 2, 2, 2)
    Lm = np.einsum("iqjs,krlt->ikqrjlst", Ra, Ra).reshape(16, 16)
    Ia = np.kron(np.diag(np.exp(h * np.diag(ZZn) / 2)), np.eye(4))
    return (Ia @ Lm @ Ia).reshape(4, 4, 4, 4)


def verflechter(R, x1, x2):
    """Kleinste Eigenwerte von M^+ M für X (L(x1) x L(x2)) = (L(x2) x L(x1)) X; Verhältnis sqrt(w0/w1)."""
    La, Lb = L_shastry(R, *x1), L_shastry(R, *x2)
    A = np.einsum("aqcs,bsdt->abqcdt", La, Lb).reshape(16, 4, 16, 4)
    B = np.einsum("aqcs,bsdt->abqcdt", Lb, La).reshape(16, 4, 16, 4)
    M = np.zeros((16, 4, 16, 4, 16, 16), complex)
    for i in range(16):
        M[i, :, :, :, i, :] += A.transpose(1, 2, 3, 0)
    for k in range(16):
        M[:, :, k, :, :, k] -= B.transpose(0, 1, 3, 2)
    M = M.reshape(4096, 256)
    w = np.linalg.eigvalsh(M.conj().T @ M)
    return float(np.sqrt(abs(w[0]) / w[1]))


U4 = 0.7
th1, th2 = 0.4, 0.9
h1s, h2s = np.arcsinh(U4 * np.sin(2 * th1)) / 2, np.arcsinh(U4 * np.sin(2 * th2)) / 2
r_auf = verflechter(r_xx, (th1, h1s), (th2, h2s))
r_neben = verflechter(r_xx, (th1, h1s), (th2, h2s + 0.1))
fx = minimize(lambda v: verflechter(r_xx, (th1, h1s), (th2, v[0] + 1j * v[1])), (0.1, 0.05), method="Nelder-Mead",
              options={"xatol": 1e-9, "fatol": 1e-12, "maxiter": 90})
ZAHLEN["Shastry-Kontrolle XX: Verhältnis sqrt(w0/w1) neben der Kurve (auf der Kurve < 1e-6)"] = fl(r_neben, 3)
ZAHLEN["Shastry-Kontrolle XX: gefundenes h2 gegen arsinh((U/4) sin 2 theta)/2"] = [fl(fx.x[0], 5), fl(h2s, 5)]
check("B2 [N] Kontrolle (Shastry): Für zwei XX-Ketten gibt es zu L(theta_1, h_1), L(theta_2, h_2) genau dann einen "
      "Verflechter, wenn sinh(2h)/sin(2 theta) für beide gleich ist (= U/4). Auf der Kurve fällt der kleinste "
      "Singulärwert auf Rundungsniveau, daneben nicht; die freie Suche in komplexem h_2 findet Shastrys Wert",
      r_auf < 1e-6 and r_neben > 1e-2 and abs(fx.x[0] - h2s) < 1e-5 and abs(fx.x[1]) < 1e-5 and fx.fun < 1e-6)
x1n = (0.4, 0.2)
r_frei = verflechter(r_nat, (0.4, 0.0), (0.9, 0.0))
r_triv = verflechter(r_nat, x1n, x1n)
suchen = {}
for u2, start in ((0.9, (0.5, 0.5)), (0.9, (-0.7, 0.3)), (1.7, (1.0, -1.0))):
    f_ = minimize(lambda v: verflechter(r_nat, x1n, (u2, v[0] + 1j * v[1])), start, method="Nelder-Mead",
                  options={"xatol": 1e-6, "fatol": 1e-9, "maxiter": 150})
    suchen[(u2, start)] = (float(f_.fun), float(f_.x[0]), float(f_.x[1]))
fj = minimize(lambda v: verflechter(r_nat, x1n, (v[0] + 1j * v[1], v[2] + 1j * v[3])), (0.9, 0.1, 0.5, 0.2),
              method="Nelder-Mead", options={"xatol": 1e-8, "fatol": 1e-12, "maxiter": 500})
ZAHLEN["Shastry-Form nativ, h_1 = 0,2: kleinstes Verhältnis und h_2 dort für u_2 = 0,9 (zwei Starts) und u_2 = 1,7"] = [
    [fl(v[0], 3), fl(v[1], 3)] for _, v in sorted(suchen.items())]
say("   nativ, (u_1, h_1) = (0,4; 0,2): kleinstes Verhältnis sqrt(w0/w1) über komplexes h_2: "
    + "; ".join(f"u_2 = {k[0]}: {fl(v[0], 3)} bei h_2 = {fl(v[1], 3)}" for k, v in sorted(suchen.items())))
check("B3 [N negativ] Dieselbe Form mit der nativen R-Matrix (Paarungsbild): Ohne Kopplung (h = 0) gibt es den "
      "Verflechter. Mit h_1 = 0,2 hat das Verhältnis über komplexes h_2 für u_2 = 0,9 ein einziges Minimum 0,092 "
      "(bei reellem h_2 = 0,254; zwei Starts), für u_2 = 1,7 das Minimum 0,128: keine Nullstelle. Die gemeinsame Suche "
      "in komplexem (u_2, h_2) läuft in den trivialen Punkt x_2 = x_1. Es gibt keine Kurve h(u) durch diesen Punkt: "
      "Shastrys Form trägt die native Kette nicht (gescopt: diese Form, diese Startwerte)",
      r_frei < 1e-6 and r_triv < 1e-6 and all(0.05 < v[0] < 0.2 and abs(v[2]) < 1e-3 for v in suchen.values())
      and abs(suchen[(0.9, (0.5, 0.5))][0] - suchen[(0.9, (-0.7, 0.3))][0]) < 1e-4
      and fj.fun < 1e-5 and max(abs(fj.x[0] - 0.4), abs(fj.x[1]), abs(fj.x[2] - 0.2), abs(fj.x[3])) < 1e-3)

say("  (c) Lokale Erhaltungsgröße der Reichweite 3 (unabhängig vom Ansatz für die R-Matrix)")
# Pauli-Strings: Platz = zwei Qubits (A, B), lokaler Index 4 pA + pB mit p in (1, X, Y, Z)
_pm = [np.eye(2), NX, NY, NZ]
PMUL = {}
for p_ in range(4):
    for q_ in range(4):
        Mm = _pm[p_] @ _pm[q_]
        for r_ in range(4):
            c_ = np.trace(Mm @ _pm[r_].conj().T) / 2
            if abs(c_) > 0.5:
                PMUL[(p_, q_)] = (r_, complex(c_))
SMUL = {}
for a_ in range(16):
    for b_ in range(16):
        ra, pa = PMUL[(a_ // 4, b_ // 4)]
        rb, pb = PMUL[(a_ % 4, b_ % 4)]
        SMUL[(a_, b_)] = (4 * ra + rb, pa * pb)


def pcomm(s1, s2):
    out, p, q = [], 1, 1
    for a, b in zip(s1, s2):
        r, ph = SMUL[(a, b)]
        out.append(r)
        p *= ph
        q *= SMUL[(b, a)][1]
    if abs(p - q) < 1e-12:
        return None
    return tuple(out), p - q


def pcanon(s):
    i = 0
    while i < len(s) and s[i] == 0:
        i += 1
    if i == len(s):
        return (), 0
    j = len(s)
    while s[j - 1] == 0:
        j -= 1
    return s[i:j], i


def ppar(s):
    return (sum(1 for x in s if x // 4 in (1, 2)) % 2, sum(1 for x in s if x % 4 in (1, 2)) % 2)


def pbasis(r, sektor):
    out = []
    for rr in range(1, r + 1):
        for s in itertools.product(range(16), repeat=rr):
            if s[0] != 0 and s[-1] != 0 and ppar(s) == sektor:
                out.append(s)
    return out


def kommutant(hs, r, sektor, mom):
    """Matrizen der Abbildung q -> [H, Q] für Q = Summe_j (+-1)^j tau^j(q), je Komponente von H, gemeinsame Zeilen."""
    Bs = pbasis(r, sektor)
    rows = {}
    datas = [([], [], []) for _ in hs]
    for col, q in enumerate(Bs):
        lq = len(q)
        for ih, h in enumerate(hs):
            acc = {}
            for hstr, hc in h:
                lh = len(hstr)
                for d in range(-lh + 1, lq):
                    lo, hi = min(0, d), max(lq, d + lh)
                    a, b = [0] * (hi - lo), [0] * (hi - lo)
                    for i, x in enumerate(q):
                        a[i - lo] = x
                    for i, x in enumerate(hstr):
                        b[d + i - lo] = x
                    c = pcomm(tuple(b), tuple(a))
                    if c is None:
                        continue
                    cs_, sh = pcanon(c[0])
                    sgn = (-1) ** ((lo + sh) % 2) if mom else 1
                    acc[cs_] = acc.get(cs_, 0) + hc * c[1] * sgn
            for cs_, v in acc.items():
                if abs(v) > 1e-13:
                    if cs_ not in rows:
                        rows[cs_] = len(rows)
                    datas[ih][0].append(rows[cs_])
                    datas[ih][1].append(col)
                    datas[ih][2].append(v)
    return [coo_matrix((dt[2], (dt[0], dt[1])), shape=(len(rows), len(Bs))).toarray() for dt in datas]


def psite(pA, pB):
    return 4 * pA + pB


def komponenten(art="paar", jxy=(1.0, 1.0)):
    hop, pot, cou = [], [], [((psite(3, 3),), 0.25)]
    for spz in (0, 1):
        s_ = (lambda p: psite(p, 0)) if spz == 0 else (lambda p: psite(0, p))
        if art == "paar":                                             # -(t/2)(XX - YY) - (m/2) Z
            hop += [((s_(1), s_(1)), -0.5), ((s_(2), s_(2)), 0.5)]
        else:                                                         # Jx XX + Jy YY
            hop += [((s_(1), s_(1)), jxy[0]), ((s_(2), s_(2)), jxy[1])]
        pot += [((s_(3),), -0.5)]
    return hop, pot, cou


def kern(M, skala=1.0, gram=None):
    """Zahl der Nullvektoren (unter den zwölf kleinsten Eigenwerten von M^+ M) und der kleinste weitere Eigenwert."""
    w = eigvalsh_teil(M.conj().T @ M if gram is None else gram, subset_by_index=[0, 11])
    n = int(np.sum(w < 1e-9 * skala))
    return n, float(w[n])


KOMM = {(mom, sek): kommutant(komponenten(), 3, sek, mom) for mom in (0, 1) for sek in ((0, 0), (1, 1), (0, 1))}
Ugitter = (0.05, 0.5, 1.0, 4.0, 12.0, 40.0, -2.0)
res_h, res_n = {}, {}
for mom in (0, 1):
    Mh, Mp_, Mu = KOMM[(mom, (0, 0))]
    res_h[mom] = [kern(Mh + U * Mu, max(1.0, U * U)) for U in Ugitter]
    res_n[mom] = [kern(Mh + 2.0 * Mp_ + U * Mu, max(1.0, U * U)) for U in Ugitter]
Mh, Mp_, Mu = KOMM[(0, (0, 0))]
M0 = Mh + 2.0 * Mp_
GA_, GB_, GC_ = M0.conj().T @ M0, M0.conj().T @ Mu, Mu.conj().T @ Mu
GB_ = GB_ + GB_.conj().T
Ufein = [float(x) for x in np.logspace(-1.3, 1.7, 31)]
Ufein = Ufein + [-x for x in Ufein]
fein = [(U, kern(None, max(1.0, U * U), GA_ + U * GB_ + U * U * GC_)) for U in Ufein]
norm_f = [x[1] * (1 + U ** 4) / U ** 2 for U, x in fein]
n_frei = [kern(KOMM[(mom, (0, 0))][0] + 2.0 * KOMM[(mom, (0, 0))][1])[0] for mom in (0, 1)]
n_andere = [kern(KOMM[(mom, sek)][0] + 2.0 * KOMM[(mom, sek)][1] + 1.0 * KOMM[(mom, sek)][2])[0]
            for mom in (0, 1) for sek in ((1, 1), (0, 1))]
mscan = [(m, kern(Mh + m * Mp_ + Mu)) for m in (0.0, 0.01, 0.1, 0.5, 1.0, 2.0, 3.0)]
ZAHLEN["Reichweite 3, Hubbard (m = 0): Zahl der Ladungen (Impuls 0, pi) für sieben Werte von U"] = [
    sorted(set(x[0] for x in res_h[0])), sorted(set(x[0] for x in res_h[1]))]
ZAHLEN["Reichweite 3, nativ (m = 2t): Zahl der Ladungen (Impuls 0, pi) für sieben Werte von U"] = [
    sorted(set(x[0] for x in res_n[0])), sorted(set(x[0] for x in res_n[1]))]
ZAHLEN["nativ, Impuls 0, 62 Werte 0,05 <= |U| <= 50: Zahl der Ladungen; kleinster und größter Wert von "
       "lambda_1 (1 + U^4)/U^2"] = [sorted(set(x[0] for _, x in fein)), fl(min(norm_f), 3), fl(max(norm_f), 1)]
ZAHLEN["nativ: kleinster Eigenwert jenseits von H bei U = 0,05; 1; 4; 12; 40 (Impuls 0)"] = [
    f"{res_n[0][Ugitter.index(U)][1]:.2e}" for U in (0.05, 1.0, 4.0, 12.0, 40.0)]
ZAHLEN["U = 1: kleinster Eigenwert jenseits von H gegen m (Impuls 0)"] = [(m, f"{x[1]:.2e}") for m, x in mscan[1:]]
say("   Zahl der Ladungen (Impuls 0, pi): frei (U = 0, m = 2t) %s; Hubbard (m = 0) mit U: %s; nativ mit U: %s"
    % (n_frei, [sorted(set(x[0] for x in res_h[mom])) for mom in (0, 1)],
       [sorted(set(x[0] for x in res_n[mom])) for mom in (0, 1)]))
say("   nativ, 62 Werte von U: lambda_1 (1 + U^4)/U^2 zwischen %s und %s" % (fl(min(norm_f), 3), fl(max(norm_f), 1)))
check("B4 [N] Kontrolle: Das Hubbard-Modell (m = 0) hat im geraden Sektor bei Impuls 0 für jedes U im Raster (sieben "
      "Werte) genau zwei translationsinvariante Ladungen der Reichweite <= 3 (H und Shastrys Dreiplatzladung), bei "
      "Impuls pi die zwei gestaffelten Ladungen", all(x[0] == 2 for x in res_h[0]) and all(x[0] == 2 for x in res_h[1]))
check("B5 [N negativ] Nativ (m = 2t) bleibt für jedes U ungleich 0 bei Impuls 0 nur H: 62 Werte, 0,05 <= |U| <= 50, "
      "beide Vorzeichen. Der kleinste weitere Eigenwert lambda_1 von M^+ M ist glatt in U, 0,044 U^2 für kleine und "
      "rund 2/U^2 für große U, ohne Nullstelle (lambda_1 (1 + U^4)/U^2 liegt zwischen 0,044 und 17). Bei Impuls pi "
      "bleiben nur die zwei Ladungen, in den Sektoren ungerader Parität nichts. Der Rest des Kommutators wächst "
      "linear in m (U = 1: lambda_1 = 7e-5 bei m = 0,01, 6,7e-3 bei m = 0,1). Die Hubbard-Kopplung zweier nativer "
      "Ketten ist in diesem Sinn nicht integrabel; integrabel ist sie nur bei m = 0, wo die Kurve entartet",
      all(x[0] == 1 and x[1] > 1e-5 for x in res_n[0]) and all(x[0] == 2 and x[1] > 1e-5 for x in res_n[1])
      and all(x[0] == 1 for _, x in fein) and 0.04 < min(norm_f) < 0.05 and max(norm_f) < 20
      and n_frei == [8, 8] and n_andere == [0, 0, 0, 0] and mscan[0][1][0] == 2
      and all(x[0] == 1 for _, x in mscan[1:]) and 1e-5 < mscan[1][1][1] < 2e-4
      and abs(mscan[2][1][1] / mscan[1][1][1] - 100) < 10)
KXY = kommutant(komponenten("xy", (1.0, 2 ** -0.5)), 3, (0, 0), 0)
n_xy0, n_xy1 = kern(KXY[0])[0], kern(KXY[0] + KXY[2])[0]
check("B6 [N negativ] Auch zwei XY-Ketten ohne Feld mit J_y/J_x = 1/sqrt 2 (je eine Nahtkette und ihre duale, z157 A6) "
      "haben mit der Kopplung U Z_A Z_B keine Dreiplatzladung (U = 1: nur H; U = 0: acht)", n_xy0 == 8 and n_xy1 == 1)


def ionisch(L, t, m, U):
    """Zwei Spezies auf L Plätzen, Moden nach Platz geordnet (A1, B1, A2, ...); Sektor halber Füllung je Spezies.
    Rückgabe: ladungsaufgelöste Niveaus -ln p der linken Hälfte in Einheiten von pi, vom tiefsten aus gezählt."""
    confs = []
    for a in itertools.combinations(range(L), L // 2):
        ma = sum(1 << (2 * j) for j in a)
        for b in itertools.combinations(range(L), L // 2):
            confs.append(ma | sum(1 << (2 * j + 1) for j in b))
    confs.sort()
    index = {c: i for i, c in enumerate(confs)}
    N = len(confs)
    ri, ci, vv = [], [], []
    for i, c in enumerate(confs):
        e = 0.0
        for j in range(L):
            na, nb = (c >> (2 * j)) & 1, (c >> (2 * j + 1)) & 1
            e += m * (-1) ** j * (na + nb) + U * (na - 0.5) * (nb - 0.5)
        ri.append(i)
        ci.append(i)
        vv.append(e)
        for j in range(L - 1):
            for s_ in (0, 1):
                p, q = 2 * j + s_, 2 * (j + 1) + s_
                if ((c >> p) & 1) != ((c >> q) & 1):
                    ri.append(index[c ^ ((1 << p) | (1 << q))])
                    ci.append(i)
                    vv.append(-t * (-1) ** ((c >> (p + 1)) & 1))
    H = coo_matrix((vv, (ri, ci)), shape=(N, N)).toarray()
    g = np.linalg.eigh(H)[1][:, 0]
    Mg = np.zeros((2 ** L, 2 ** L))
    for i, c in enumerate(confs):
        Mg[c & (2 ** L - 1), c >> L] = g[i]
    lv = []
    for na in range(L // 2 + 1):
        for nb in range(L // 2 + 1):
            rows = [x for x in range(2 ** L) if bin(x & 0x5555).count("1") == na and bin(x & 0xAAAA).count("1") == nb]
            cols = [x for x in range(2 ** L) if bin(x & 0x5555).count("1") == L // 2 - na
                    and bin(x & 0xAAAA).count("1") == L // 2 - nb]
            for x in np.linalg.svd(Mg[np.ix_(rows, cols)], compute_uv=False)[:4]:
                if x > 1e-9:
                    lv.append((-2 * np.log(x), na, nb))
    lv.sort()
    return [((e - lv[0][0]) / np.pi, na, nb) for e, na, nb in lv[:6]], (lv[0][1], lv[0][2])


leiter = {}
for U in (0.0, 0.5, 1.0, -1.0, 2.0):
    lv, (a0, b0) = ionisch(6, 1.0, 2.0, U)
    eins = [e for e, na, nb in lv if abs(na - a0) + abs(nb - b0) == 1][0]
    beide = [e for e, na, nb in lv if abs(na - a0) == 1 and abs(nb - b0) == 1 and (na - a0) == (nb - b0)][0]
    leiter[U] = (eins, beide - 2 * eins)
ZAHLEN["Halbkette zweier gekoppelter Ladungsketten (6 Plätze): erstes Niveau/pi und Additivitätsfehler/pi gegen U"] = {
    str(U): [fl(v[0], 3), fl(v[1], 3)] for U, v in sorted(leiter.items())}
say("   Halbkettenleiter (erstes Niveau/pi, Fehler der Additivität/pi): "
    + ", ".join(f"U = {U}: {fl(v[0], 3)}, {fl(v[1], 3)}" for U, v in sorted(leiter.items())))
check("B7 [N] Die Halbkette der Hubbard-gekoppelten Ketten verliert die native Leiter schon in erster Ordnung in U: "
      "Bei U = 0 liegt das erste Niveau bei pi und die Niveaus addieren sich (sechs Plätze, Fehler 2e-3); bei U = +-1 "
      "verschiebt es sich um -+0,1 pi, und das Paarniveau weicht von der Summe ab. Das wiederholt z155 (Eichkopplung) "
      "für die Dichtekopplung",
      abs(leiter[0.0][0] - 1) < 5e-3 and abs(leiter[0.0][1]) < 5e-3 and leiter[1.0][0] < 0.93 and leiter[-1.0][0] > 1.07
      and abs(leiter[1.0][1]) > 0.05 and abs((leiter[0.5][0] - 1) / (leiter[1.0][0] - 1) - 0.5) < 0.1)

say("  (d) Abelscher Twist: die gesteuerte Uhr als Kopplung, die die Yang-Baxter-Gleichung erhält")
# drei Plätze, je zwei Qubits (A_j, B_j): Qubit 2j = A_j, 2j + 1 = B_j
GA = lambda u, j: emb(npm(Gn(mp.mpmathify(u))), [2 * j, 2 * j + 2], 6)
GB = lambda u, j: emb(npm(Gn(mp.mpmathify(u))), [2 * j + 1, 2 * j + 3], 6)


def Dtw(phi, j):
    return emb(np.diag(np.exp(-1j * phi * np.diag(ZZn) / 4)), [2 * j, 2 * j + 3], 6)


def Gtw(u, j, phi):
    D = Dtw(phi, j)
    return D @ GA(u, j) @ GB(u, j) @ D.conj().T


ut, vt = 0.37 + 0.2j, -0.21 + 0.05j
tw_ok, tw_rest = True, []
for phi in (0.0, np.pi / 2, 0.7, np.pi):
    lhs = Gtw(ut - vt, 0, phi) @ Gtw(ut, 1, phi) @ Gtw(vt, 0, phi)
    rhs = Gtw(vt, 1, phi) @ Gtw(ut, 0, phi) @ Gtw(ut - vt, 1, phi)
    tw_rest.append(rel(lhs, rhs))
    tw_ok &= tw_rest[-1] < 1e-12
Om = np.eye(64, dtype=complex)
for (j, k_) in ((0, 1), (0, 2), (1, 2)):
    Om = Om @ emb(np.diag(np.exp(1j * (-1) ** (j + k_) * (np.pi / 2) * np.diag(ZZn) / 4)), [2 * j, 2 * k_ + 1], 6)
om_ok = all(rel(Om @ GA(0.6, j) @ GB(0.6, j) @ Om.conj().T, Gtw(0.6, j, np.pi / 2)) < 1e-12 for j in (0, 1))


def osr_AB(G4):
    """Operator-Schmidt-Rang eines Operators auf (A_1, B_1, A_2, B_2) bezüglich der Teilung A gegen B."""
    T = G4.reshape([2] * 8).transpose(0, 2, 4, 6, 1, 3, 5, 7).reshape(16, 16)
    return int(np.sum(np.linalg.svd(T, compute_uv=False) > 1e-10))


def Gtw4(u, phi):
    D = emb(np.diag(np.exp(-1j * phi * np.diag(ZZn) / 4)), [0, 3], 4)
    g = npm(Gn(mp.mpmathify(u)))
    return D @ emb(g, [0, 2], 4) @ emb(g, [1, 3], 4) @ D.conj().T


osr_frei, osr_tw, osr_uhr = osr_AB(Gtw4(0.6, 0.0)), osr_AB(Gtw4(0.6, np.pi / 2)), osr_AB(Gtw4(float(VARPI), np.pi / 2))
uhr_gleich = rel(Gtw4(float(VARPI), np.pi / 2), Gtw4(float(VARPI), 0.0))
QQ = np.diag([0.0, 0.0, 0.0, 1.0])                                      # Q_A Q_B mit Q = (1 - Z)/2
cuhr = np.diag(np.exp(1j * np.pi / 2 * np.diag(QQ)))                    # i^(Q_A Q_B) = diag(1, 1, 1, i)
lokal = np.kron(np.diag(np.exp(-1j * np.pi / 8 * np.diag(NZ))), np.diag(np.exp(-1j * np.pi / 8 * np.diag(NZ))))
uhr_ist_D = rel(cuhr, np.exp(1j * np.pi / 8) * lokal @ np.diag(np.exp(-1j * (-np.pi / 2) * np.diag(ZZn) / 4)))
ZAHLEN["Twist: Operator-Schmidt-Rang A gegen B bei u = 0,6 (frei, phi = pi/2) und am Uhrpunkt"] = [osr_frei, osr_tw, osr_uhr]
check("B8 [E/N] Für jeden Winkel phi erfüllt das verdrillte Gatter D G_A(u) G_B(u) D^(-1) mit "
      "D = exp(-i phi Z^A_j Z^B_(j+1)/4) die Yang-Baxter-Gleichung auf drei Plätzen (phi = 0, pi/2, 0,7, pi; Rest < 1e-12). "
      "Grund: Es ist das freie Gatter, konjugiert mit Omega = Produkt über j < k von exp(i phi q^A_j q^B_k), "
      "q_j = (-1)^j Z_j/2 die Ladung im Paarungsbild (geprüft für phi = pi/2 auf beiden Bindungen). Das ist ein "
      "Reshetikhin-Twist; die Ketten sind nichtlokal frei, lokal korreliertes Hüpfen wie bei Schulz und Shastry",
      tw_ok and om_ok)
check("B9 [E] Der Twistfaktor D ist bei phi = -+pi/2 bis auf Einplatzphasen die gesteuerte Uhr i^(+-Q_A Q_B) "
      "zwischen dem A-Platz j und dem B-Platz j + 1. Das verdrillte Gatter wechselwirkt lokal (Operator-Schmidt-Rang "
      "zwischen den Ketten 4 statt 1), aber am Uhrpunkt u = varpi ist es gleich dem unverdrillten Produkt der zwei Paaruhren (Rang 1): Die Paaruhr ist "
      "diagonal und sieht den Twist nicht. Die gesteuerte Uhr sitzt im Rahmen, nicht im Wert der R-Matrix",
      osr_frei == 1 and osr_tw == 4 and osr_uhr == 1 and uhr_gleich < 1e-12 and uhr_ist_D < 1e-12)
g_p = lambda phi: mp.diag([mp.exp(-I * phi / 4), mp.exp(I * phi / 4)])   # exp(-i phi Z/4) auf dem zweiten B-Platz
ok = True
for u in (mp.mpf("0.41"), mp.mpc("0.3", "0.7")):
    for phi in (PI / 2, mp.mpf("0.7")):
        gp, gm = kron(E2, g_p(phi)), kron(E2, g_p(-phi))
        zweig_p, zweig_m = gp * Gn(u) * gp ** -1, gm * Gn(u) * gm ** -1         # Z^A = +1 und Z^A = -1
        verschoben = gp * Gn(VARPI) * Gn(u - VARPI) * gp ** -1
        d_ = nrm(zweig_m - verschoben)
        ok &= (d_ < TOL) if phi == PI / 2 else (d_ > mp.mpf("0.01"))
check("B10 [E, neu] Der Twist ist eine ladungsgesteuerte Verschiebung des Spektralparameters: Die zwei Ladungszustände "
      "des A-Platzes sehen B-Gatter, die sich um die Eichung diag(1, e^(-i phi)) unterscheiden. Nach A6/A7 ist das "
      "genau für phi in (pi/2) Z eine Translation auf E_i; bei phi = pi/2 ist es die Verschiebung um den reellen "
      "2-Torsionspunkt varpi mal Paaruhr. phi = pi ist die Periode 2 varpi (Eichung S^2 = sigma_z, "
      "die Parität), phi = +-pi/2 sind die zwei nichtentarteten Bicharaktere i^(+-ab) der Uhrgruppe (Bestand z144)", ok)

say("  (e) Ising-Form: Die Acht-Vertex-Familie ist die integrable Kopplung der Nahtkette an ihre duale Kette")
Xs = [op(NX, j, 4) for j in range(4)]
Ys = [op(NY, j, 4) for j in range(4)]
Zs = [op(NZ, j, 4) for j in range(4)]
z1a, z1b, b1 = Xs[0] @ Xs[1], Xs[2] @ Xs[3], -Ys[1] @ Ys[2]              # Kette 1 (Nahtkette): Plätze, Bindung
z2a, z2b, b2 = Ys[0] @ Ys[1], Ys[2] @ Ys[3], Xs[1] @ Xs[2]               # Kette 2 (dual)
check("B11 [E] In der XYZ-Kette ist der Term J_z ZZ die Kopplung der zwei Ising-Ketten aus z157 A6: auf ungeraden "
      "Bindungen Z Z = -(Platzterm Kette 1)(Platzterm Kette 2), auf geraden Z Z = (Bindung 1)(Bindung 2). Das ist die "
      "Vierspinkopplung von Kadanoff/Wegner bzw. die Ashkin-Teller-Kopplung. Die Familie aus z157 (lam von K nach 0) "
      "ist also die integrable Kopplung der Nahtkette (Feld sqrt 2) an ihre Kramers-Wannier-duale Kette (Feld "
      "1/sqrt 2), mit J_z/J_x = cn lam dn lam/(1 + k sn^2 lam) als freiem Kopplungswert",
      np.linalg.norm(Zs[0] @ Zs[1] + z1a @ z2a) < 1e-12 and np.linalg.norm(Zs[2] @ Zs[3] + z1b @ z2b) < 1e-12
      and np.linalg.norm(Zs[1] @ Zs[2] - b1 @ b2) < 1e-12)
befund("B", "nein (Hubbard-Form), bedingt (Twist: P_schub; Ising-Form: P_lam frei)",
       "Drei Kopplungen sind gerechnet. (1) Dichtekopplung der Ladungsform (ionisches Hubbard-Modell): keine "
       "Shastry-Lösung und keine Dreiplatzladung für U ungleich 0; integrabel nur bei m = 0. Die Naht bleibt nicht die "
       "Naht. (2) Abelscher Twist: für jeden Winkel integrabel, aber nichtlokal frei; die Kurve zeichnet die Winkel "
       "k pi/2 aus, weil nur sie Translationen um Torsionspunkte sind. Am Uhrpunkt ist der Twist unsichtbar. "
       "(3) Ising-Form: Die Acht-Vertex-Familie koppelt die Nahtkette integrabel an ihre duale Kette; der Wert ist lam.")

# --------------------------------------------------------------------------------------------------------------------
say("\nTeil C  Der Limes lam -> 0 in der symmetrischen Klasse")


def w8(u, lam, k=KB):
    s, su, sul = sn(lam, k), sn(u, k), sn(u + lam, k)
    return sul, su, s, k * s * su * sul


def R8(u, lam, k=KB):
    a, b, c, d = w8(u, lam, k)
    return mp.matrix([[a, 0, 0, d], [0, b, c, 0], [0, c, b, 0], [d, 0, 0, a]])


def Rc8(u, lam, k=KB):
    a, b, c, d = w8(u, lam, k)
    return mp.matrix([[a, 0, 0, d], [0, c, b, 0], [0, b, c, 0], [d, 0, 0, a]]) / c


def deltas(u, lam, k=KB):
    a, b, c, d = w8(u, lam, k)
    pp, pm, sp_, sm = a + d, a - d, c + b, c - b
    return pm * sm / (pp * sp_), pp * sm / (pm * sp_), sp_ * sm / (pp * pm)


def Ffun(x):
    """F(x) = e^(-i pi/4) sl(e^(i pi/4) (1 + k) x/2) = slh((1 + k) x/2), der hyperbolische lemniskatische Sinus."""
    z = mp.exp(I * PI / 4)
    return sl(z * (1 + KB) * x / 2) / z


def theta_fam(lam):
    return mp.atan(mp.re(Ffun(lam)) ** 2)


ok = True
fracs = [mp.mpf(x) for x in ("0.001", "0.1", "0.25", "0.5", "0.75", "0.9", "1", "1.3", "1.7", "1.95")]
ths = []
for fr in fracs:
    lam = fr * K_
    dx, dy, dz = deltas(I * lam, lam)
    th = theta_fam(lam)
    ths.append(th)
    ok &= (abs(mp.im(Ffun(lam))) < TOL and abs(dy + I) < TOL and abs(dx + I * mp.exp(-2 * I * th)) < mp.mpf(10) ** -35
           and abs(dz + I * mp.exp(2 * I * th)) < mp.mpf(10) ** -35)
w_ = (1 + I) * mp.mpf("0.37")
ok &= abs(Ffun(-I * w_) + I * Ffun(w_)) < TOL
ok &= abs(Ffun(K_) - 1) < TOL and abs(1 / Ffun(2 * K_ - mp.mpf(10) ** -20)) < mp.mpf(10) ** -15
ok &= all(ths[i] < ths[i + 1] for i in range(len(ths) - 1)) and abs(ths[6] - PI / 4) < TOL
check("C1 [E] Geschlossene Form der Familie u = i lam (Baxters Modul k = tan^2(pi/8)): (delta_x, delta_y, delta_z) = "
      "(-i e^(-2 i theta), -i, -i e^(2 i theta)), also (c1, c2, c3) = (pi/8 + theta/2, pi/8, pi/8 - theta/2), mit "
      "tan theta = slh^2((1 + k) lam/2), slh der hyperbolische lemniskatische Sinus. theta wächst streng von 0 (lam = 0) "
      "über pi/4 (lam = K) auf pi/2 (lam = 2K). delta_y = F(-i w)/F(w) = -i folgt allein aus sl(i z) = i sl(z) "
      "(zehn Werte von lam, 35 Stellen)", ok)
uhr_th = sp.symbols("th", real=True)
G120 = (E4 + I * PERM) / (1 + I)
lsg = sp.solveset(sp.Eq(sp.cos(8 * uhr_th), 1), uhr_th, sp.Interval(0, sp.pi / 2))
ok = set(lsg) == {sp.Integer(0), sp.pi / 4, sp.pi / 2}
dz8 = deltas(I * K_ / 2, K_ / 2)[2]


def prop_mp(A, B):
    i, j = max(((i, j) for i in range(B.rows) for j in range(B.cols)), key=lambda t: abs(B[t[0], t[1]]))
    return nrm(A - A[i, j] / B[i, j] * B) / nrm(A)


eps_l = mp.mpf(10) ** -12 * K_
Y1 = kron(SY, E2)
ende_ok = (prop_mp(Rc8(I * eps_l, eps_l), G120) < mp.mpf(10) ** -10
           and prop_mp(Rc8(I * (2 * K_ - eps_l), 2 * K_ - eps_l), XX * Y1 * G120 * Y1) < mp.mpf(10) ** -10
           and prop_mp(Rc8(I * K_, K_), G120) > mp.mpf("0.1"))
check("C2 [E] Uhrwerte um alle drei Achsen (delta^4 = 1) verlangen e^(8 i theta) = 1, also theta in {0, pi/4, pi/2}: "
      "lam = 0 (z120-Gatter), lam = K (freie Nahtkette), lam = 2K (das z120-Gatter in der Eichung sigma_y auf einem "
      "Platz, mal dem lokalen Produkt XX). z157 C14 (dort numerisch, nur an den Enden exakt) ist damit exakt; "
      "dazwischen ist delta_z keine vierte Einheitswurzel (bei lam = K/2 ist Re delta_z = sqrt 2 - 1)",
      ok and abs(mp.re(dz8) - (mp.sqrt(2) - 1)) < TOL and ende_ok)


def expw(c1, c2, c3):
    def ex(c, M):
        return mp.cos(c) * E4 + I * mp.sin(c) * M
    return ex(c1, XX) * ex(c2, YY) * ex(c3, ZZ)


def uhrkomm(G, achse, vz):
    """Kommutator mit der Produktuhr i^(Q_1 + vz Q_2) um eine Achse, Q = (1 - sigma)/2."""
    S_ = {"x": SX, "y": SY, "z": SZ}[achse]
    c1 = mp.cos(PI / 4) * E2 - I * mp.sin(PI / 4) * S_
    c2 = c1 if vz > 0 else c1.H
    C = kron(c1, c2)
    return nrm(G * C - C * G) / nrm(G)


def gatter_fam(lam):
    G = Rc8(I * lam, lam)
    return G / mp.sqrt(mp.sqrt(mp.det(G)))


ok = True
kom_tab = {}
G0_ = gatter_fam(mp.mpf(10) ** -12 * K_)
G2_ = gatter_fam(2 * K_ - mp.mpf(10) ** -12 * K_)
for fr in (mp.mpf("0.5"), mp.mpf("1"), mp.mpf("1.5")):
    G = gatter_fam(fr * K_)
    th = theta_fam(fr * K_)
    W = expw(PI / 8 + th / 2, PI / 8, PI / 8 - th / 2)
    ph = G[1, 1] / W[1, 1]
    ok &= nrm(G - ph * W) < mp.mpf(10) ** -35 * nrm(G)
    vals = [uhrkomm(G, a, v) for a in "xyz" for v in (1, -1)]
    kom_tab[ns(fr, 3)] = fl(min(vals), 3)
    ok &= min(vals) > mp.mpf("0.05")
ok &= all(uhrkomm(G0_, a, 1) < mp.mpf(10) ** -10 for a in "xyz") and uhrkomm(G2_, "x", -1) < mp.mpf(10) ** -10
ok &= sum(1 for a in "xyz" for v in (1, -1) if uhrkomm(G2_, a, v) < mp.mpf(10) ** -10) == 1
ok &= all(0 < x < PI for t in ths for x in (PI / 4 + t / 2, PI / 4, PI / 4 - t / 2, t / 2, t))


def schiefe_uhr(G, n1, n2):
    c1 = (np.eye(2) - 1j * (n1[0] * NX + n1[1] * NY + n1[2] * NZ)) / np.sqrt(2)
    c2 = (np.eye(2) - 1j * (n2[0] * NX + n2[1] * NY + n2[2] * NZ)) / np.sqrt(2)
    C = np.kron(c1, c2)
    return float(np.linalg.norm(G @ C - C @ G) / np.linalg.norm(G))


kugel = [(np.sin(a) * np.cos(b), np.sin(a) * np.sin(b), np.cos(a))
         for a in np.linspace(0, np.pi, 13) for b in np.linspace(0, 2 * np.pi, 24, endpoint=False)]
schief = {}
for fr in (mp.mpf("0.5"), mp.mpf("1"), mp.mpf("1.5")):
    Gs = npm(gatter_fam(fr * K_))
    schief[ns(fr, 3)] = fl(min(schiefe_uhr(Gs, n, tuple(v * x for x in n)) for n in kugel for v in (1, -1)), 3)
    ok &= float(schief[ns(fr, 3)]) > 0.05
ZAHLEN["dasselbe für 312 Achsen auf der Kugel (gleicher oder entgegengesetzter Drehsinn)"] = schief
ZAHLEN["kleinster Kommutator des Gatters mit einer Produktuhr bei lam/K = 0,5; 1; 1,5"] = kom_tab
check("C3 [E für die Koordinatenachsen, N für schiefe Achsen; neu] Das Gatter bei u = i lam vertauscht mit einer "
      "Produktuhr i^(Q_1 +- Q_2) um eine Koordinatenachse genau dann, wenn die zwei anderen Koordinaten gleich oder "
      "entgegengesetzt sind: c_a -+ c_b in pi Z. Auf der Familie sind die Differenzen theta/2, theta/2, theta und die "
      "Summen pi/4 + theta/2, pi/4, pi/4 - theta/2: Das geht nur bei theta = 0 (lam = 0: alle drei Achsen, also die "
      "volle SU(2)) und theta = pi/2 (lam = 2K: eine gestaffelte U(1)). Für schiefe Achsen (312 Richtungen) bleibt "
      "der Kommutator im Inneren über 0,05. Ladungserhaltung des Gatters legt lam also fest; am freien Ende lam = K "
      "gibt es keine erhaltene Ladung",
      ok)
jz_rate = []
for eps_ in (mp.mpf("0.1"), mp.mpf("0.01")):
    s_ = sn(eps_ * K_, KB)
    jx, jy, jz = 1 + KB * s_ * s_, 1 - KB * s_ * s_, cn(eps_ * K_, KB) * dn(eps_ * K_, KB)
    jz_rate.append(((jx - jy) / (eps_ * K_) ** 2, (jy - jz) / (eps_ * K_) ** 2))
ok = (abs(jz_rate[1][0] - 2 * KB) < mp.mpf("1e-4") and abs(jz_rate[1][1] - (1 - KB) ** 2 / 2) < mp.mpf("1e-4")
      and abs(sn(K_, KB) - 1) < TOL)
for fr in (mp.mpf("0.3"), mp.mpf("1"), mp.mpf("1.6")):
    s_ = sn(fr * K_, KB)
    ok &= abs((1 - KB * s_ * s_) ** 2 - (cn(fr * K_, KB) * dn(fr * K_, KB)) ** 2 - (1 - KB) ** 2 * s_ * s_) < TOL
dG = [nrm(gatter_fam(mp.mpf(10) ** -n * K_) / gatter_fam(mp.mpf(10) ** -n * K_)[0, 0]
          - (E4 + I * PERM) / (1 + I)) for n in (2, 3)]
ok &= abs(dG[0] / dG[1] - 100) < 1 and nrm((E4 + I * PERM) / (1 + I) - (1 * E4 + I * PERM) / (1 + I)) == 0
check("C4 [E; Kritikalität L] Am Ende lam -> 0 ist die Kette die XXX-Kette: J_x - J_y = 2 k lam^2, J_y - J_z = "
      "(1 - k)^2 lam^2/2. Für 0 < lam < 2K ist J_x > J_y > |J_z|, denn J_x - J_y = 2 k sn^2 lam und "
      "J_y^2 - J_z^2 = (1 - k)^2 sn^2 lam; zwei Kopplungen sind nur an den Enden betragsgleich. Masselos ist die "
      "Familie [L] nur dort: lam -> 0 antiferromagnetisch kritisch (su(2)_1, c = 1), lam -> 2K dasselbe mit "
      "umgekehrtem Vorzeichen von H (ferromagnetisch). Das Gatter geht quadratisch gegen (1 + iP)/(1 + i) = Yangs R(u) = u + iP bei u = 1, das Destri-de-Vega-Gitter Theta = 1/2 aus "
      "z121; P_aut (u/lam = i) ist dort dieselbe Aussage wie Theta = 1/2", ok)


def k1_von(fr):
    q = mp.exp(-PI * fr)
    return (mp.jtheta(2, 0, q) / mp.jtheta(3, 0, q)) ** 2


masse = {ns(fr, 3): ns(-mp.log(k1_von(fr)), 6) for fr in (mp.mpf(1), mp.mpf("0.75"), mp.mpf("0.5"), mp.mpf("0.25"),
                                                           mp.mpf("0.1"))}
ZAHLEN["Masse m(lam) = -ln k_1 mit Nome e^(-pi lam/K) [L, am freien Punkt geeicht], lam/K = 1; 0,75; 0,5; 0,25; 0,1"] = masse
jx1, jy1 = 1 + KB, 1 - KB
ok = (abs(-mp.log(k1_von(1)) - mp.log(mp.sqrt(2))) < TOL and abs(jy1 / jx1 - 1 / mp.sqrt(2)) < TOL
      and abs(-mp.log(k1_von(mp.mpf("0.1"))) / (8 * mp.exp(-10 * PI)) - 1) < mp.mpf("1e-10"))
say("   Masse je Platz [L]: " + ", ".join(f"lam/K = {a}: {b}" for a, b in masse.items()))
check("C5 [E am freien Punkt, sonst L] Korrelationslänge entlang der Familie (antiferromagnetisches Vorzeichen, "
      "Johnson/Krinsky/McCoy, Baxter 10.12): Masse m(lam) = -ln k_1, k_1 der Modul zum Nome e^(-pi lam/K). Bei lam = K "
      "ist k_1 = 1/sqrt 2 und m = ln sqrt 2, die exakte inverse Korrelationslänge der XY-Kette mit "
      "J_y/J_x = 1/sqrt 2. Für lam -> 0 ist m = 8 e^(-pi K/lam): wesentliche Singularität, kein Potenzgesetz. Die "
      "Formel für lam < K ist Literatur und hier nicht unabhängig geprüft (Korrelationen bis 18 Plätze lösen sie in "
      "der Erkundung nicht auf)", ok)


def xyz_grund(L, jx, jy, jz):
    N = 2 ** L
    idx = np.arange(N, dtype=np.int64)
    terms, diag = [], np.zeros(N)
    for j in range(L - 1):
        s1, s2 = 1 - 2 * ((idx >> j) & 1), 1 - 2 * ((idx >> (j + 1)) & 1)
        terms.append((idx ^ ((1 << j) | (1 << (j + 1))), jx - jy * s1 * s2))
        diag += jz * s1 * s2

    def mv(v):
        out = diag * v
        for flip, co in terms:
            out += co * v[flip]
        return out

    v0 = np.cos(0.37 * np.arange(N) + 0.11) + 0.5                      # fester Startvektor
    w, v = eigsh(LinearOperator((N, N), matvec=mv, dtype=float), k=2, which="SA", v0=v0, tol=1e-11)
    g = v[:, np.argmin(w)]
    p = np.linalg.svd(g.reshape(2 ** (L // 2), 2 ** (L - L // 2)), compute_uv=False) ** 2
    return -np.log(p[:8] / p[0])


ed = {}
for fr, Led in ((mp.mpf("1.5"), 18), (mp.mpf("1.25"), 14)):
    lam = fr * K_
    s_ = sn(lam, KB)
    lv = xyz_grund(Led, float(1 + KB * s_ * s_), float(1 - KB * s_ * s_), float(cn(lam, KB) * dn(lam, KB)))
    ed[ns(fr, 3)] = [float(x) / float(PI * fr) for x in (lv[2], lv[4], lv[6])]
ZAHLEN["XYZ-Halbkette: Niveaus in Einheiten pi lam/K bei lam/K = 1,5 (18 Plätze) und 1,25 (14 Plätze)"] = {
    a: [fl(x, 2) for x in b] for a, b in ed.items()}
say("   Halbkettenniveaus der XYZ-Kette in Einheiten pi lam/K: "
    + "; ".join(f"lam/K = {a}: {', '.join(fl(x, 2) for x in b)}" for a, b in ed.items()))
check("C6 [L, N] Leiter der Halbkette entlang der Familie (Baxter Kap. 13: Dichtematrix = Produkt von diag(1, x^(2j)); "
      "hier mit x^2 = e^(-pi lam/K)): Niveaus j pi lam/K. Geprüft durch exakte Diagonalisierung bei lam = 1,5 K "
      "(18 Plätze: 1,00; 2,00; 2,98) und 1,25 K (14 Plätze, auf 5 Prozent). Bei lam = K sind es die Vielfachen von pi: "
      "die Niveaus pi (2l + 1) der Nahtkette und 2 pi l ihrer dualen Kette (z152 B1, B2). Die Einheit ist streng "
      "monoton in lam: Die native Leiter gibt es nur bei lam = K",
      all(abs(ed["1.5"][i] - (i + 1)) < 0.02 * (i + 1) for i in range(3))
      and all(abs(ed["1.25"][i] - (i + 1)) < 0.05 * (i + 1) for i in range(3)))


def t_inh(v, xis, lam):
    Mt = np.eye(2, dtype=complex).reshape(2, 2, 1, 1)
    for xi in xis:
        R4 = npm(R8(mp.mpmathify(v - xi), lam)).reshape(2, 2, 2, 2)
        Mt = np.einsum("xsbt,baPQ->xaPsQt", R4, Mt).reshape(2, 2, Mt.shape[2] * 2, Mt.shape[3] * 2)
    return np.einsum("aaPQ->PQ", Mt)


def prop(A, B):
    i = np.unravel_index(np.argmax(abs(B)), B.shape)
    return float(np.linalg.norm(A - A[i] / B[i] * B) / np.linalg.norm(A))


Lk = 6
lk_ok, lk_rest = True, []
for fr in (mp.mpf(1), mp.mpf("0.5"), mp.mpf("0.001")):
    lam = fr * K_
    ug = complex(I * lam)
    Gg = npm(gatter_fam(lam))
    xis = [(-1) ** j * ug / 2 for j in range(Lk)]
    F = t_inh(ug / 2, xis, lam) @ np.linalg.inv(t_inh(-ug / 2, xis, lam))
    ev_, od_ = np.eye(2 ** Lk, dtype=complex), np.eye(2 ** Lk, dtype=complex)
    for j in range(0, Lk, 2):
        ev_ = ev_ @ emb(Gg, [j, (j + 1) % Lk], Lk)
    for j in range(1, Lk, 2):
        od_ = od_ @ emb(Gg, [j, (j + 1) % Lk], Lk)
    tv = t_inh(0.3 + 0.2j, xis, lam)
    r1, r2 = prop(F, od_ @ ev_), float(np.linalg.norm(tv @ F - F @ tv) / np.linalg.norm(tv) / np.linalg.norm(F))
    r3 = float(np.linalg.norm(Gg @ Gg.conj().T - np.eye(4)))
    lk_rest.append(max(r1, r2, r3))
    lk_ok &= max(r1, r2, r3) < 1e-10
check("C7 [N] Für jedes lam ist die Ziegelschaltung aus dem unitären Gatter bei u = i lam integrabel: Sie ist "
      "t(u/2) t(-u/2)^(-1) der inhomogenen Acht-Vertex-Transfermatrix mit Inhomogenitäten +-u/2 und vertauscht mit "
      "t(v) (Ring mit sechs Plätzen; lam = K, K/2, K/1000). Das ist das elliptische Lichtkegelgitter; bei lam -> 0 ist "
      "es das Gitter aus z121. Das Verhältnis u/lam = i bleibt fest: lam -> 0 ist der rationale Keim der Kurve, nicht "
      "der Doppellimes Theta -> unendlich von Destri und de Vega", lk_ok)
befund("C", "bedingt: P_lam folgt aus der Ladungserhaltung des Gatters (in der symmetrischen Klasse, unter P_aut)",
       "Die Familie ist geschlossen lösbar: tan theta = slh^2((1 + k) lam/2). Drei unabhängige Bedingungen wählen "
       "Punkte: Uhrwerte um alle Achsen (lam = 0, K, 2K), erhaltene Ladung bzw. Vertauschen mit der Uhr (lam = 0: "
       "SU(2); lam = 2K: gestaffelte U(1)), native Leiter (nur lam = K). Masselos ist die Familie nur an den "
       "rationalen Enden; die Naht als Zustand liegt in der Mitte. P_lam ist damit keine freie Zahl mehr, aber "
       "Zustand (lam = K) und Gatter (lam -> 0) bleiben zwei verschiedene Punkte der Familie.")

# --------------------------------------------------------------------------------------------------------------------
say("\nTeil D  Die zwei Lesarten des Viertels: Villain/Theta (z155) und komplexe Multiplikation (z157)")


def th2(n, tau, z, M=14):
    """Thetafunktion der Stufe 2: Summe über r = n mod 4 von q^(r^2/8) zeta^r."""
    return mp.fsum(mp.exp(2 * PI * I * tau * mp.mpf(r * r) / 8 + 2 * PI * I * z * r) for r in range(n - 4 * M, n + 4 * M + 1, 4))


def villain(n, a, M=14):
    return mp.fsum(mp.exp(-a * mp.mpf((n + 4 * m) ** 2) / 2) for m in range(-M, M + 1))


Smat = mp.matrix(4, 4)
for n in range(4):
    for m in range(4):
        Smat[n, m] = (-I) ** (n * m) / 2
wv = mp.matrix([villain(n, PI / 2) for n in range(4)])
ok = all(abs(wv[n] - th2(n, I, 0)) < TOL for n in range(4)) and nrm(Smat * wv - wv) < TOL
a_ = mp.mpf("0.9")
dual = mp.matrix([villain(n, PI ** 2 / (4 * a_)) for n in range(4)]) * mp.sqrt(2 * PI / a_) / 2
ok &= nrm(Smat * mp.matrix([villain(n, a_) for n in range(4)]) - dual) < TOL
check("D1 [E] z155 nachgerechnet: Das Villain-Gewicht auf Z_4 ist die Thetakonstante der Stufe 2, "
      "w_a(n) = theta_n(tau = 2 i a/pi, 0). Die Z_4-Fourier-Transformation S = (1/2) i^(-n n') schickt a nach "
      "pi^2/(4a); selbstdual ist genau a = pi/2, also tau = i", ok)
zt = mp.mpc("0.13", "0.07")
lhs = mp.matrix([th2(n, I, -I * zt) for n in range(4)])
rhs = mp.exp(4 * PI * zt * zt) * Smat * mp.matrix([th2(n, I, zt) for n in range(4)])
check("D2 [E] Dieselbe Matrix S wirkt bei tau = i auf die Thetafunktionen mit elliptischem Argument: "
      "theta_n(i, -i z) = e^(4 pi z^2) Summe_n' S_(n n') theta_n'(i, z). Die Modultransformation S hat tau = i als "
      "Fixpunkt und wirkt dort auf z als die komplexe Multiplikation mit -i: Die Automorphie von E_i und die "
      "Z_4-Fourier-Transformation sind dieselbe Abbildung", nrm(lhs - rhs) < TOL * nrm(lhs))
Cmat = mp.matrix(4, 4)
for n in range(4):
    Cmat[n, (-n) % 4] = 1
v_odd = mp.matrix([0, 1, 0, -1])
evs = sorted([complex(x) for x in mp.eig(Smat)[0]], key=lambda z: (round(z.real, 6), round(z.imag, 6)))
ok = (nrm(Smat * Smat - Cmat) < TOL and nrm(Smat ** 4 - mp.eye(4)) < TOL and nrm(Smat * v_odd + I * v_odd) < TOL
      and max(abs(a - b) for a, b in zip(evs, [-1, -1j, 1, 1])) < 1e-12)
check("D3 [E] S hat die Ordnung 4, S^2 ist die Ladungskonjugation. Eigenwerte: 1, 1, -1 auf dem geraden Sektor und -i "
      "auf dem ungeraden Vektor (0, 1, 0, -1). Die Villain-Gewichte sind gerade (Eigenwert 1: Selbstdualität, z155); "
      "die ungerade Kombination theta_1 - theta_3 trägt den Eigenwert -i", ok)
vt_ = mp.mpc("0.21", "0.16")
q_i = mp.exp(-PI)
t1 = lambda v: mp.jtheta(1, PI * v, q_i)
t3 = lambda v: mp.jtheta(3, PI * v, q_i)
ok = (abs(th2(1, I, vt_ / 2) - th2(3, I, vt_ / 2) - I * t1(vt_)) < TOL
      and abs(t1(-I * vt_) + I * mp.exp(PI * vt_ ** 2) * t1(vt_)) < TOL
      and abs(t3(-I * vt_) - mp.exp(PI * vt_ ** 2) * t3(vt_)) < TOL
      and abs(sl(VARPI * vt_) - t1(vt_) / t3(vt_)) < TOL
      and abs(sl(-I * VARPI * vt_) + I * sl(VARPI * vt_)) < TOL)
check("D4 [E] Die ungerade Kombination ist Jacobis theta_1: theta_1 - theta_3 (Stufe 2, Argument z) = i theta_1(2z | i), "
      "und theta_1(-i v | i) = -i e^(pi v^2) theta_1(v | i), theta_3(-i v | i) = e^(pi v^2) theta_3(v | i). Der "
      "lemniskatische Sinus ist ihr Quotient, sl(varpi v) = theta_1(v | i)/theta_3(v | i) (ohne Vorfaktor), also "
      "sl(-i x) = -i sl(x): Der Multiplikator -i aus z157 ist der Eigenwert von S auf dem ungeraden Sektor", ok)
Tm = mp.diag([mp.exp(-I * PI * mp.mpf(n * n) / 4) for n in range(4)])   # Gatter e^(-i (pi/2) Q^2/2) = inverse T-Matrix
pol_ok = all(abs(Tm[(n + m) % 4, (n + m) % 4] / (Tm[n, n] * Tm[m, m]) - (-I) ** (n * m)) < TOL
             for n in range(4) for m in range(4))
phase = mp.exp(-I * PI / 4)
gauss4 = mp.fsum(mp.exp(-I * PI * mp.mpf(n * n) / 4) for n in range(4))
ok = (pol_ok and nrm((Smat * Tm) ** 3 - phase * mp.eye(4)) < TOL
      and nrm((Smat ** -1 * Tm) ** 3 - phase * Smat * Smat) < TOL and abs(gauss4 - 2 * phase) < TOL)
check("D5 [E] Die gesteuerte Uhr ist beides zugleich: i^(-n n') ist (bis auf 1/2) der Kern von S und die Polarisation "
      "der Gatterphase T(n) = e^(-i pi n^2/4): T(n + n')/(T(n) T(n')) = i^(-n n'). Es gilt (S T)^3 = e^(-i pi/4) und "
      "(S^(-1) T)^3 = e^(-i pi/4) S^2 mit der Gauß-Summe Summe_n T(n) = 2 e^(-i pi/4) (Bestand z145, z155). S (Ordnung 4, fixiert tau = i) und T (tau -> tau - 1) sind "
      "verschiedene Elemente derselben Weil-Darstellung", ok)
Had = mp.matrix([[1, 1], [1, -1]]) / mp.sqrt(2)
g2 = lambda x: 1 + I * x                                                  # Summe über n in Z_2 von i^(n^2) x^n
G120 = (E4 + I * PERM) / (1 + I)
ok = (nrm(Had * mp.diag([1, -I]) * Had - (E2 + I * SX) / (1 + I)) < TOL
      and abs(g2(-1) / g2(1) + I) < TOL and abs(g2(1) / abs(g2(1)) - mp.exp(2 * PI * I / 8)) < TOL
      and nrm(G120 * (E4 + PERM) / 2 - (E4 + PERM) / 2) < TOL and nrm(G120 * (E4 - PERM) / 2 + I * (E4 - PERM) / 2) < TOL)
check("D6 [E; Lesart L] Arithmetik des z120-Gatters: (1 + iP)/(1 + i) = g(P)/g(1) mit der Gauß-Summe "
      "g(x) = Summe über n in Z_2 von i^(n^2) x^n; auf dem antisymmetrischen Sektor ist es g(-1)/g(1) = (1 - i)/(1 + i) "
      "= -i, und g(1)/|g(1)| = e^(2 pi i/8). Es ist die mit der Z_2-Fourier-Matrix (Hadamard) konjugierte inverse "
      "Uhr: H diag(1, -i) H = (1 + iX)/(1 + i). Lesart: Der Nenner 1 + i ist die Gauß-Summe einer Ladung mit c = 1", ok)
befund("D", "E: dieselbe Aussage in komplementären Sektoren; der reelle Schritt T ist nicht mit erfasst",
       "Beide Lesarten sind die Wirkung des Stabilisators von tau = i (erzeugt von S, Ordnung 4) in der "
       "Weil-Darstellung der Uhrgruppe Z_4. z155 sieht den geraden Sektor bei z = 0 (Eigenwert 1: Selbstdualität des "
       "Villain-Gewichts), z157 den ungeraden Sektor bei z ungleich 0 (Eigenwert -i: sl(-i x) = -i sl(x), delta_y = "
       "-i). Kein Zufall, sondern ein Satz. Nicht dasselbe sind die Prämissen: P_aut benutzt S; P_viertel als Laufzeit "
       "pi/2 unter Q^2/2 (z145, Form iii) ist T. Beide liefern denselben Bicharakter i^(-n n'), weil er Kern von S "
       "und Polarisation von T ist.")

# --------------------------------------------------------------------------------------------------------------------
say("\nTeil E  Bilanz")
praemissen = [
    ("P_aut", "benutzt (z157)", "Der Abstand der zwei Linienscharen ist das Aut(E_i)-Bild der Kreuzung, u = i lam. "
     "Gleichwertig Theta = 1/2 im Lichtkegelgitter (C4, C7). In Teil D als Wirkung von S erkannt"),
    ("P_lam", "reduziert", "In der symmetrischen Klasse unter P_aut keine freie Zahl mehr: Uhrwerte um alle Achsen "
     "lassen lam = 0, K, 2K (C2); Vertauschen des Gatters mit der Produktuhr (P_ladung) lässt nur lam = 0 (und die "
     "Eichkopie 2K) (C3)"),
    ("P_ladung", "neu, Folge von P_uhr", "Das Gatter erhält die Ladung der zwei Plätze, d. h. es vertauscht mit der "
     "Produktuhr i^(Q_1 + Q_2). Das folgt aus P_uhr ('die Phase ist eine Funktion der Uhrwerte') und ist schwächer "
     "als P_null. Um alle drei "
     "Achsen verlangt, ist es die SU(2)-Invarianz des Gatters; als Lesart das Gittergegenstück zu P_sl2 aus z133 "
     "(dort die native konforme sl(2), hier die SU(2) der zwei Zustände je Platz; nicht dieselbe Gruppe)"),
    ("P_schub", "neu, Alternative zu P_viertel", "Die Kopplung zweier Ketten ist eine ladungsgesteuerte Translation auf "
     "E_i (abelscher Twist). Dann ist der Kopplungsfaktor eine Funktion der Uhrwerte (P_null) und der Winkel in "
     "(pi/2) Z. Dass die Verschiebung genau eine Halbperiode varpi ist (nichtentartet, +-pi/2) und nicht 0 oder "
     "2 varpi, gehört zur Prämisse (B10)"),
    ("P_platz, P_buendel, P_wurzel, P_spek", "benutzt (z157)", "Rahmen der symmetrischen elliptischen Klasse; unverändert"),
    ("P_CTM", "benutzt (z152)", "Die Naht ist die Halbkette; für B7 und C6"),
    ("P_sym", "unverändert offen", "Der symmetrische Strahlteiler ist kein Wert der nativen R-Matrix (A5); der Twist "
     "liefert ihn nicht"),
]
for nme, st, txt in praemissen:
    say(f"   {nme} [{st}]: {txt}")
bilanz = [
    ("Teil 1", "E", "Native R-Matrix in Normalform auf x^2 + y^2 + x^2 y^2 = 1: B(u) = cl - i sl (sigma_z - cl sigma_x), "
     "sl, cl bei u/2. Yang-Baxter 40 Stellen. U(1), Parität, Z_4: (sl, cl) -> (i sl, 1/cl). 2-Torsion: Paaruhr, Parität, "
     "Rang 1; 4-Torsion: Transfer i tan(pi/8). Neu: G(u + varpi) = G(varpi) S G(u) S^(-1), die Uhr ist die Eichung der "
     "Halbperiode, und nur Winkel k pi/2 sind Translationen."),
    ("Teil 2, Shastry", "N negativ", "Dichtekopplung = ionisches Hubbard-Modell (exakte Operatoridentität). Shastrys "
     "Form reproduziert für XX sinh(2h)/sin(2 theta) = U/4; nativ kein Verflechter. Keine Dreiplatzladung für U "
     "ungleich 0 (62 Werte), Bruch in erster Ordnung in m und in U. U ist nicht frei, nicht diskret, sondern 0. Gilt "
     "für Platzkopplung, "
     "Reichweite 3, Shastry-Form im Hilfsraum."),
    ("Teil 2, Twist", "E/C", "Abelscher Twist mit der gesteuerten Phase: Yang-Baxter für jeden Winkel, nichtlokal frei. "
     "Auf E_i sind genau phi = k pi/2 Translationen (2-Torsion = Uhr); +-pi/2 nichtentartet. Am Uhrpunkt unsichtbar: "
     "Das Gatter ist dort die Paaruhr, nicht diag(1,1,1,-i)."),
    ("Teil 2, Ising-Form", "E/L", "Die Acht-Vertex-Familie ist die Kopplung der Nahtkette an ihre duale Kette "
     "(Kadanoff/Wegner). Der Wert lam ist in der Yang-Baxter-Gleichung frei."),
    ("Teil 3", "E/L/N", "tan theta = slh^2((1 + k) lam/2); Uhrwerte nur bei lam = 0, K, 2K; erhaltene Ladung nur bei "
     "lam = 0 (SU(2)) und 2K (gestaffelte U(1)); masselos nur dort; Masse -ln k_1(e^(-pi lam/K)), bei lam -> 0 wie "
     "8 e^(-pi K/lam) [L]; Leitereinheit pi lam/K, nativ nur bei lam = K [L, N]. Jedes lam gibt ein integrables Lichtkegelgitter; lam -> 0 ist z121."),
    ("Teil 4", "E", "Eine Aussage: S der Weil-Darstellung bei tau = i; gerader Sektor = Selbstdualität (z155), ungerader "
     "= Multiplikator -i (z157). Der Kern von S ist die gesteuerte Uhr und zugleich die Polarisation von T."),
    ("P_uhr", "C", "P_null: für den Kopplungsfaktor Folge von P_schub (Twist); in der symmetrischen Klasse genügt die "
     "schwächere P_ladung. P_viertel: folgt aus P_aut (z157) oder aus P_schub (Torsion). Beide Wege ersetzen P_uhr "
     "durch andere Prämissen (P_ladung + P_aut bzw. P_schub); die Zahl der Prämissen sinkt nicht."),
    ("P_lam, P_aut, P_sym", "C", "P_lam folgt aus P_ladung. P_aut bleibt, jetzt gleich 'Theta = 1/2' und gleich 'S'. "
     "P_sym bleibt offen: Die native R-Matrix erreicht höchstens Transfer (sqrt 2 - 1)^2."),
    ("K13, K15, K16", "O", "K13: nichts (1+1 Dimensionen). K15: Die Dichtekopplung ist nicht integrabel und zerstört "
     "die Naht; der Twist ist nichtlokal frei, also keine Dynamik; die symmetrische Familie hat SU(2) nur am "
     "masselosen Ende. K16: Masselos ist die Familie nur an den rationalen Enden, dort aber nicht die Naht."),
    ("Neun Prämissen", "C", "A*, I_step, P_lift, P_IF, P_univ, P_ord, P_curv, P_G: unberührt. P_uhr bleibt, mit zwei "
     "bedingten Ersatzformen (P_ladung + P_aut in der symmetrischen Klasse; P_schub für zwei Ladungsketten). Die Zahl "
     "bleibt neun; OriginLab ist nicht zu ändern."),
    ("Nächster Test", "O", "Die zwei Enden verbinden: Halbkette der kritischen XXX-Schaltung (z121, lam -> 0) gegen die "
     "Naht. Konkret: Ist die Naht (lam = K) der stationäre Zustand bzw. die Randbedingung des Lichtkegelgitters bei "
     "lam -> 0, d. h. gibt es eine integrable Grenzfläche zwischen lam = K und lam -> 0 auf derselben Kurve "
     "(Baxters Z-invariantes Gitter mit zwei Linienscharen verschiedener Kreuzungsparameter)? Ja: Zustand und Gatter "
     "liegen in einem Modell. Nein: P_uhr bleibt endgültig ein eigener Schritt."),
]
for nme, st, txt in bilanz:
    say(f"   [{st}] {nme}: {txt}")

urteil = ("Bedingt, in der wörtlichen Form nein. (1) Die Shastry-artige Kopplung zweier nativer Ladungsketten ist das "
          "ionische Hubbard-Modell bei Wechselpotential 2t. Sie ist nicht integrabel: kein Verflechter in Shastrys "
          "Form, keine Dreiplatzladung für U ungleich 0, und die Halbkette verliert die native Leiter in erster "
          "Ordnung. Integrabel ist die Kopplung nur bei m = 0, wo die Kurve E_i entartet. Die Kurve legt also kein U "
          "fest; sie lässt nur U = 0 zu (gescopt: Platzkopplung, die beide Ladungen erhält; Reichweite 3; "
          "Shastry-Form). (2) Die gesteuerte Uhr koppelt zwei native Ketten trotzdem exakt, aber als abelscher Twist: "
          "für jeden Winkel gilt die Yang-Baxter-Gleichung, und das Modell ist nichtlokal frei. Neu ist, dass die "
          "Kurve hier den Winkel auszeichnet: Die Uhr S = i^n ist die Eichung der Verschiebung um den reellen "
          "2-Torsionspunkt, G(u + varpi) = G(varpi) S G(u) S^(-1), und nur Eichwinkel k pi/2 sind Translationen auf "
          "E_i. Unter der Prämisse P_schub (Kopplung = ladungsgesteuerte Torsionstranslation) folgen P_null für den "
          "Kopplungsfaktor und der Viertelwinkel bis auf das Vorzeichen. Am Uhrpunkt selbst ist das Gatter dann aber die Paaruhr, nicht "
          "diag(1,1,1,-i); die gesteuerte Uhr sitzt im Rahmen. (3) In der Ising-Form ist die Acht-Vertex-Familie aus "
          "z157 die integrable Kopplung der Nahtkette an ihre duale Kette. Die Familie u = i lam ist geschlossen "
          "lösbar, tan theta = slh^2((1 + k) lam/2). Uhrwerte um alle Achsen gibt es exakt nur bei lam = 0, K, 2K; "
          "eine erhaltene Ladung (Vertauschen mit der Produktuhr) nur an den rationalen Enden, und dort ist die Kette "
          "masselos; bei lam -> 0 ist das Gatter SU(2)-symmetrisch. P_lam ist damit keine freie Zahl: Es folgt aus der "
          "Ladungserhaltung des "
          "Gatters (P_ladung, eine Folge von P_uhr). Der Limes lam -> 0 ist in diesem Sinn nativ ausgezeichnet, aber als "
          "rationaler Keim bei festem u/lam = i (Theta = 1/2), nicht als Kontinuumslimes der Naht: Die native Leiter "
          "gibt es nur bei lam = K. (4) Die zwei Lesarten des Viertels sind eine: Selbstdualität des Villain-Gewichts "
          "(z155) und Multiplikator -i der komplexen Multiplikation (z157) sind der gerade und der ungerade Sektor "
          "derselben Matrix S bei tau = i; ihr Kern ist die gesteuerte Uhr. (5) Für den Kandidaten: P_uhr bleibt "
          "Prämisse, ist aber enger gefasst (P_ladung + P_aut oder P_schub); P_lam entfällt als eigene Zahl; P_sym, "
          "P_ord und die neun Prämissen bleiben. K13, K15, K16 sind nicht gelöst. Nichts ist befördert.")
say("\nUrteil\n   " + urteil)
print(f"\n   (Laufzeit {time.time() - T0:.1f} s, nicht Teil der Ausgabe)", flush=True)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z160_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "praemissen": [list(x) for x in praemissen],
     "bilanz": [list(x) for x in bilanz], "urteil": urteil}, indent=1, ensure_ascii=False, default=str) + "\n")
(HERE / "ausgabe_z160.txt").write_text("\n".join(OUT) + "\n")
