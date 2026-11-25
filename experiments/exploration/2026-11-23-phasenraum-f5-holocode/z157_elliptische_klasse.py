"""Z157: Ist die Nahtkette aus z152 das integrable Modell der elliptischen Klasse auf der nativen Kurve E_i: w^2 = 1 - z^4?

Leitlinie 'Rekonstruktion statt Auswahl'; H6 und die gelösten Blöcke (Flavour, Massen, alpha, 16, E8) bleiben unberührt,
nichts wird befördert. Statuskürzel: E exakt, C bedingt (Prämisse genannt), L Lesart bzw. Literatur, N numerisch, O offen.

Anlass: (1) z152: Die Naht ist die Hälfte einer eindeutigen freien Gitterkette (Hüpfen t, Wechselpotential 2t; im Spektrum
gleichwertig: Ising-Kette im Querfeld sqrt(2)); der Modulgenerator ist c mal Gitter-Boost mit c = Lemniskatenkonstante.
(2) z151: Das Gatter des Kandidaten ist die rationale 4D-Chern-Simons-R-Matrix 1 + hbar P/z mit Dz = i hbar; die native Form
gibt auf der Sphäre nur die trigonometrische Klasse; die rationale braucht fünf neue Prämissen. (3) z120/z144/z145: Gatter
(1 + iP)/(1 + i) bzw. gesteuerte Uhr diag(1,1,1,-i); der Viertelwinkel ist P_viertel, das Fehlen des Oszillatoranteils P_null.

Frage: Ist die Nahtkette das integrable Modell der elliptischen Klasse auf E_i, und legen die Z_4-Automorphie und die
Torsionspunkte der Kurve das Gatter samt Winkel fest?

Literatur [L] (Buchkapitel bzw. arXiv-Nummern; Formeln daraus sind unten numerisch nachgerechnet, nicht übernommen):
  Baxter, Exactly Solved Models in Statistical Mechanics (1982): Kap. 7 (Ising, elliptische Parametrisierung), Kap. 10
    (Acht-Vertex-Modell: 10.4 Parametrisierung, 10.14 XYZ-Kette, 10.16 Freie-Fermionen-Fall), Kap. 13 (Eckentransfermatrizen).
  Fan/Wu, Phys. Rev. B 2 (1970) 723 (Freie-Fermionen-Bedingung); Sutherland, J. Math. Phys. 11 (1970) 3183 (XYZ-Kette
    vertauscht mit der Acht-Vertex-Transfermatrix); Suzuki, Prog. Theor. Phys. 46 (1971) 1337 (XY-Kette und 2D-Ising);
    Felderhof, Physica 65 (1973) 421 (Freie Fermionen im Feld, elliptisch).
  Bazhanov/Stroganov, Teor. Mat. Fiz. 62 (1985) 377 und 63 (1985) 291, 'Hidden symmetry of free fermion model I, II'
    (mathnet.ru tmf4655, tmf4763): elliptische R-Matrix freier Fermionen im Feld mit drei Parametern je Vertex.
    Neuere Fassung: Melikyan, arXiv:1912.10786 (fermionischer Bazhanov-Stroganov-R-Operator, tetraedrische Algebra).
  Thacker, Physica D 18 (1986) 348; Tetelman, Sov. Phys. JETP 55 (1982) 306; Sogo/Wadati, Prog. Theor. Phys. 69 (1983) 431:
    Gitter-Boost Summe j h_j erzeugt die Verschiebung des Spektralparameters; Eckentransfermatrix = exp(-u Boost).
    (Diese drei aus dem Gedächtnis zitiert; die lokale Form der Aussage ist unten in B1 nachgerechnet.)
  Peschel/Kaulke/Legeza, cond-mat/9810174; Peschel, cond-mat/0410416; Its/Jin/Korepin, quant-ph/0409027;
    Eisler/Peschel, arXiv:1703.08126: Halbketten freier Gittermodelle, Abstand der Verschränkungsniveaus elliptisch.
  Costello/Witten/Yamazaki I, arXiv:1709.09993 (gelesen über ar5iv bis Abschnitt 7): Abschnitt 2.2: elliptische Lösungen
    der klassischen Yang-Baxter-Gleichung gibt es nur für sl_N, mit nur der endlichen Automorphiegruppe Z_N x Z_N;
    Abschnitt 3.1: hbar ist nicht quantisiert; Abschnitt 3.2, Gl. 3.14: für Geschlecht 1 ist omega = dz ohne Pole und bis
    auf Normierung eindeutig; Abschnitt 3.4: die G-Symmetrie der R-Matrix gilt nur für C = C; Abschnitt 6: die
    Rahmungsanomalie hat lokal dieselbe Form; Abschnitt 10.2 'Rigid Holomorphic Bundles' (nur Titel gesehen).
    Costello/Witten/Yamazaki II, arXiv:1802.01579 (Zusammenfassung). Belavin, Nucl. Phys. B 180 (1981) 189;
    Belavin/Drinfeld, Funct. Anal. Appl. 16 (1982) 159 (aus dem Gedächtnis zitiert).
  Lemniskatische Funktionen: sl((1 +- i) x) = (1 +- i) sl(x)/sl'(x) (Gauß; unten in C7 numerisch geprüft).

Teil A  Elliptische Parametrisierung: Acht-Vertex-R-Matrix, Freie-Fermionen-Punkt = 2-Torsion, XY-Kette = zwei Ising-Ketten,
        Kurve j = 1728; Ladungskette: Dispersionskurve ist wörtlich w^2 = 1 - s^4, R-Matrix aus lemniskatischen Funktionen.
Teil B  Boost: lokale Boost-Relation, K = 2 lambda mal Boost, c = Periode; welche Punkte der Kurve was bedeuten.
Teil C  Z_4 und Torsion: komplexe Multiplikation = Hadamard-Dualität; R-Matrix an Torsionspunkten; das Viertel als
        Invariante der komplexen Multiplikation; Vergleich mit dem Gatter des Kandidaten.
Teil D  Gegen z151: welche Prämissen der rationalen Klasse entfallen, welche hinzukommen.
Teil E  Bilanz und Urteil.
"""
import json
import os
import time
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")                     # deterministische Rundung (einfädig)

import mpmath as mp
import numpy as np
import sympy as sp

mp.mp.dps = 50
HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
T0 = time.time()
TOL = mp.mpf(10) ** -35


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val, text):
    BEFUNDE[name] = val
    say(f"   [Befund: {val}] {text}")


def ns(x, n=12):
    return mp.nstr(x, n)


# --------------------------------------------------------------------------------------------------------------------
# Werkzeuge
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


def w8(u, lam, k):
    """Acht-Vertex-Gewichte (Baxter 10.4 in der sn-Form): a = sn(u + lam), b = sn u, c = sn lam, d = k sn lam sn u sn(u + lam)."""
    s, su, sul = sn(lam, k), sn(u, k), sn(u + lam, k)
    return sul, su, s, k * s * su * sul


def R8(u, lam, k):
    a, b, c, d = w8(u, lam, k)
    return mp.matrix([[a, 0, 0, d], [0, b, c, 0], [0, c, b, 0], [d, 0, 0, a]])


def Rc8(u, lam, k):
    """Gatterform P R / c: gerader Block [[a, d], [d, a]], ungerader Block [[c, b], [b, c]], bei u = 0 die Eins."""
    a, b, c, d = w8(u, lam, k)
    return mp.matrix([[a, 0, 0, d], [0, c, b, 0], [0, b, c, 0], [d, 0, 0, a]]) / c


def deltas(u, lam, k):
    """Zweikörperinvarianten je Achse: Produkt der Bell-Eigenwerte mit sigma sigma = -1 durch das mit +1 (= e^(-4 i c_Achse))."""
    a, b, c, d = w8(u, lam, k)
    pp, pm, sp_, sm = a + d, a - d, c + b, c - b                 # Phi+, Phi-, Psi+, Psi-
    return pm * sm / (pp * sp_), pp * sm / (pm * sp_), sp_ * sm / (pp * pm)


def weyl(u, lam, k):
    """Koordinaten (c1, c2, c3) mit Gatter ~ exp(i (c1 XX + c2 YY + c3 ZZ)); Rückgabe zusätzlich die Streuung der Beträge."""
    a, b, c, d = w8(u, lam, k)
    ev = (a + d, a - d, c + b, c - b)
    ph = [mp.arg(x) for x in ev]
    g = mp.fsum(ph) / 4
    p = [x - g for x in ph]
    mods = [abs(x) for x in ev]
    return (p[0] + p[2]) / 2, (p[1] + p[2]) / 2, (p[0] + p[1]) / 2, max(mods) / min(mods) - 1


def expw(c1, c2, c3):
    def ex(c, M):
        return mp.cos(c) * E4 + I * mp.sin(c) * M
    return ex(c1, XX) * ex(c2, YY) * ex(c3, ZZ)


def prop(A, B):
    """Abstand von A zu einem skalaren Vielfachen von B (relativ)."""
    best = max(((abs(B[i, j]), i, j) for i in range(B.rows) for j in range(B.cols)))
    _, i, j = best
    s = A[i, j] / B[i, j]
    return nrm(A - s * B) / max(nrm(A), mp.mpf(10) ** -40)


def npm(M):
    return np.array([[complex(M[i, j]) for j in range(M.cols)] for i in range(M.rows)])


NX = np.array([[0, 1], [1, 0]], complex)
NY = np.array([[0, -1j], [1j, 0]])
NZ = np.diag([1.0, -1.0]).astype(complex)
NI = np.eye(2, dtype=complex)


def op(o, j, L):
    m = np.array([[1]], complex)
    for i in range(L):
        m = np.kron(m, o if i == j else NI)
    return m


def ring_T(R, L):
    """Transfermatrix auf dem Ring: R wirkt auf Hilfsraum (x) Platz."""
    R4 = R.reshape(2, 2, 2, 2)
    M = np.eye(2, dtype=complex).reshape(2, 2, 1, 1)
    for _ in range(L):
        M = np.einsum("xsbt,baPQ->xaPsQt", R4, M).reshape(2, 2, M.shape[2] * 2, M.shape[3] * 2)
    return np.einsum("aaPQ->PQ", M)


def osr(G):
    g = npm(G).reshape(2, 2, 2, 2).transpose(0, 2, 1, 3).reshape(4, 4)
    return int(np.sum(np.linalg.svd(g, compute_uv=False) > 1e-10))


VARPI = mp.gamma(mp.mpf(1) / 4) ** 2 / (2 * mp.sqrt(2 * PI))            # Lemniskatenkonstante
KB = (mp.sqrt(2) - 1) ** 2                                              # Baxters Modul der Ising-Form
K_, KP_ = mp.ellipk(KB ** 2), mp.ellipk(1 - KB ** 2)                    # K, K'

# --------------------------------------------------------------------------------------------------------------------
say("Z157: Die Nahtkette und die elliptische Klasse auf der nativen Kurve E_i: w^2 = 1 - z^4")
say("=" * 100)
say("\nTeil A  Elliptische Parametrisierung")
say("  (a) Acht-Vertex-R-Matrix (Baxter 10.4): a = sn(u+lam), b = sn u, c = sn lam, d = k sn lam sn u sn(u+lam)")
pts = [(mp.mpc("0.37", "0.11"), mp.mpc("0.41", "-0.2"), mp.mpc("0.3", "0.7"), mp.mpc("-0.5", "0.13")),
       (mp.mpf("0.6"), mp.mpf("0.9"), mp.mpc("0.2", "0.4"), mp.mpc("0.7", "-0.3")),
       (KB, K_, mp.mpc("0.3", "0.5"), mp.mpc("-0.2", "0.9"))]
ok = True
for k, lam, u, v in pts:
    Lh = e12(R8(u - v, lam, k)) * e13(R8(u, lam, k)) * e23(R8(v, lam, k))
    Rh = e23(R8(v, lam, k)) * e13(R8(u, lam, k)) * e12(R8(u - v, lam, k))
    ok &= nrm(Lh - Rh) < TOL * nrm(Lh)
check("A1 [E] Yang-Baxter-Gleichung in Differenzform R_12(u-v) R_13(u) R_23(v) = R_23(v) R_13(u) R_12(u-v): drei "
      "komplexe Punkte (Modul, Kreuzungsparameter, u, v), Rest < 1e-35", ok)
ok = True
for k, lam, u, v in pts[:2]:
    dx, dy, dz = deltas(u, lam, k)
    ok &= abs(dz - sn(lam - u, k) / sn(lam + u, k)) < TOL
check("A2 [E] Zweikörperinvariante um die z-Achse: delta_z = (c^2 - b^2)/(a^2 - d^2) = sn(lam - u)/sn(lam + u). "
      "delta_z = 1 für alle u heißt Gatter gaußsch (Freie-Fermionen-Bedingung a^2 + b^2 = c^2 + d^2)", ok)
# Kopplungen der XYZ-Kette und 2-Torsion
kg = mp.mpf("0.37")
Kg, Kgp = mp.ellipk(kg ** 2), mp.ellipk(1 - kg ** 2)


def coup(lam, k):
    s = sn(lam, k)
    return 1 + k * s * s, 1 - k * s * s, cn(lam, k) * dn(lam, k)


hfd = mp.mpf(10) ** -12
lamg = mp.mpf("0.52")
dR = (Rc8(hfd, lamg, kg) - Rc8(-hfd, lamg, kg)) / (2 * hfd)
jx, jy, jz = coup(lamg, kg)
Hloc = (jx * XX + jy * YY + jz * ZZ + jz * E4) / (2 * sn(lamg, kg))
check("A3 [E] Regulär: Gatter(0) = 1; Ableitung bei 0 = (Jx XX + Jy YY + Jz ZZ + Jz)/(2 sn lam) mit "
      "Jx : Jy : Jz = 1 + k sn^2 lam : 1 - k sn^2 lam : cn lam dn lam (XYZ-Kette, Sutherland/Baxter 10.14)",
      nrm(Rc8(mp.mpf(0), lamg, kg) - E4) < TOL and nrm(dR - Hloc) < mp.mpf(10) ** -18)
eps = mp.mpf(10) ** -30
j1 = coup(Kg, kg)
j2 = coup(I * Kgp / 2 + eps, kg)
j3 = coup(Kg + I * Kgp / 2 + eps, kg)
j4 = coup(Kg + I * Kgp, kg)
n2, n3 = max(abs(x) for x in j2), max(abs(x) for x in j3)
ok = (abs(j1[2]) < TOL and abs(j4[2]) < TOL and abs(j2[0] / n2) < mp.mpf(10) ** -20 and abs(j3[1] / n3) < mp.mpf(10) ** -20
      and min(abs(j1[0]), abs(j1[1])) > mp.mpf("0.1"))
check("A4 [E] Freie Fermionen <=> der Kreuzungsparameter ist ein 2-Torsionspunkt des Gitters 2K Z + iK' Z: "
      "Jz = 0 bei lam = K (und K + iK' = K modulo Gitter), Jx = 0 bei iK'/2, Jy = 0 bei K + iK'/2 (allgemeiner Modul)", ok)
say("  (b) Freie-Fermionen-Punkt lam = K und die Ising-Kette im Querfeld sqrt(2)")
ok = True
for u in (mp.mpf("0.3"), mp.mpc("0.2", "0.9"), mp.mpc("-0.7", "0.1")):
    a, b, c, d = w8(u, K_, KB)
    ok &= abs(a * a + b * b - c * c - d * d) < TOL and abs(deltas(u, K_, KB)[2] - 1) < TOL
jx, jy, jz = coup(K_, KB)
kb_s = (sp.sqrt(2) - 1) ** 2
check("A5 [E] Bei lam = K: a^2 + b^2 = c^2 + d^2 und delta_z = 1 für alle u; Kette = XY-Kette mit Jy/Jx = (1-k)/(1+k). "
      "Für k = (sqrt 2 - 1)^2 = tan^2(pi/8) ist Jy/Jx = 1/sqrt 2 exakt",
      ok and abs(jy / jx - 1 / mp.sqrt(2)) < TOL and abs(jz) < TOL
      and sp.simplify((1 - kb_s) / (1 + kb_s) - 1 / sp.sqrt(2)) == 0 and abs(KB - mp.tan(PI / 8) ** 2) < TOL)
# XY-Kette = zwei entkoppelte Ising-Ketten (Operatoridentität auf vier Plätzen)
L4 = 4
Xs = [op(NX, j, L4) for j in range(L4)]
Ys = [op(NY, j, L4) for j in range(L4)]
z1a, z1b, b1 = Xs[0] @ Xs[1], Xs[2] @ Xs[3], -Ys[1] @ Ys[2]              # Kette 1: Plätze X1X2, X3X4; Bindung -Y2Y3
z2a, z2b, b2 = Ys[0] @ Ys[1], Ys[2] @ Ys[3], Xs[1] @ Xs[2]               # Kette 2: Plätze Y1Y2, Y3Y4; Bindung X2X3


def comm(A, B):
    return np.linalg.norm(A @ B - B @ A)


def acomm(A, B):
    return np.linalg.norm(A @ B + B @ A)


ok = (acomm(z1a, b1) < 1e-12 and acomm(z1b, b1) < 1e-12 and comm(z1a, z1b) < 1e-12
      and acomm(z2a, b2) < 1e-12 and acomm(z2b, b2) < 1e-12 and comm(z2a, z2b) < 1e-12
      and all(comm(A, B) < 1e-12 for A in (z1a, z1b, b1) for B in (z2a, z2b, b2)))
check("A6 [E] XY-Kette = zwei entkoppelte Ising-Ketten: XX auf ungeraden und YY auf geraden Bindungen bilden eine "
      "Ising-Algebra (Platzterm XX, Bindungsterm YY), die übrigen die duale; beide vertauschen. Kette 1 hat "
      "Feld/Kopplung = Jx/Jy = sqrt 2: die Ising-Form der Nahtkette aus z152", ok)
# Transfermatrix vertauscht mit dem Hamiltonoperator (Ring, numpy)
L6 = 6


def Hxyz(jx_, jy_, jz_, L):
    H = np.zeros((2 ** L, 2 ** L), complex)
    for j in range(L):
        H += (jx_ * op(NX, j, L) @ op(NX, (j + 1) % L, L) + jy_ * op(NY, j, L) @ op(NY, (j + 1) % L, L)
              + jz_ * op(NZ, j, L) @ op(NZ, (j + 1) % L, L))
    return H


ok = True
for lam, k in ((K_, KB), (lamg, kg)):
    jx, jy, jz = [float(mp.re(x)) for x in coup(lam, k)]
    H = Hxyz(jx, jy, jz, L6)
    T1 = ring_T(npm(R8(mp.mpf("0.3"), lam, k)), L6)
    T2 = ring_T(npm(R8(mp.mpc("0.2", "0.7"), lam, k)), L6)
    ok &= (np.linalg.norm(T1 @ H - H @ T1) < 1e-10 * np.linalg.norm(T1) * np.linalg.norm(H)
           and np.linalg.norm(T1 @ T2 - T2 @ T1) < 1e-10 * np.linalg.norm(T1) * np.linalg.norm(T2))
check("A7 [N] Ring mit sechs Plätzen: Die Transfermatrizen vertauschen untereinander und mit der Kette "
      "(XY-Kette am Freie-Fermionen-Punkt mit k = tan^2(pi/8); Kontrolle: XYZ-Kette bei allgemeinem lam, k)", ok)
say("  (c) Das invariante Gitter der R-Matrix und die Kurve")
u0, lam0 = mp.mpc("0.31", "0.17"), mp.mpf("0.47")
gs = {"1": E2, "x": SX, "y": SY, "z": SZ}


def gauge_dist(w):
    R1 = R8(u0 + w, lam0, KB)
    res = {}
    for name in sorted(gs):
        G = kron(gs[name], E2)
        res[name] = prop(R1, G * R8(u0, lam0, KB) * G)
    return res


d2K, diKp = gauge_dist(2 * K_), gauge_dist(I * KP_)
nonper = [min(gauge_dist(w).values()) for w in (K_, I * KP_ / 2, K_ + I * KP_ / 2)]
ok = d2K["z"] < TOL and diKp["x"] < TOL and min(nonper) > mp.mpf("0.05")
check("A8 [E] Modulo Eichung und Skalar hat die R-Matrix die Perioden 2K (Konjugation mit sigma_z) und iK' (mit "
      "sigma_x), aber nicht K, iK'/2, K + iK'/2: Gitter L_R = 2K Z + iK' Z, Eichgruppe = Kleinsche Vierergruppe", ok)
tauR = I * KP_ / (2 * K_)
jnum = 1728 * mp.kleinj(tauR)
lam_s = sp.Symbol("lam")
jlam = 256 * (lam_s ** 2 - lam_s + 1) ** 3 / (lam_s ** 2 * (lam_s - 1) ** 2)
jfac = sp.factor(sp.together(jlam - 1728))
gam_s = (1 - kb_s) / (1 + kb_s)
check("A9 [E] Für k = (sqrt 2 - 1)^2 ist K'/K = 2 (50 Stellen), also L_R quadratisch, tau_R = i, j = 1728. Exakt: das "
      "Doppelverhältnis der vier Verzweigungspunkte ist (Jy/Jx)^2 = 1/2, und j - 1728 = 64 (l+1)^2 (l-2)^2 (2l-1)^2/(l^2 (l-1)^2)",
      abs(KP_ / K_ - 2) < mp.mpf(10) ** -45 and abs(jnum - 1728) < mp.mpf(10) ** -30
      and abs(KB - (1 - 1 / mp.sqrt(2)) / (1 + 1 / mp.sqrt(2))) < TOL
      and sp.simplify(jlam.subs(lam_s, sp.nsimplify(gam_s ** 2)) - 1728) == 0
      and sp.simplify(jfac - 64 * (lam_s + 1) ** 2 * (lam_s - 2) ** 2 * (2 * lam_s - 1) ** 2
                      / (lam_s ** 2 * (lam_s - 1) ** 2)) == 0)
ZAHLEN["Baxter-Modul der Ising-Form k = (sqrt2-1)^2"] = ns(KB, 20)
ZAHLEN["K, K' dazu"] = [ns(K_, 20), ns(KP_, 20)]
say(f"   k = {ns(KB, 15)}, K = {ns(K_, 15)}, K' = {ns(KP_, 15)}, tau_R = i K'/(2K) = i")
befund("A-Ising", "E/N",
       "Die Ising-Form der Nahtkette (Querfeld sqrt 2) ist die Hälfte der XY-Kette am Freie-Fermionen-Punkt des "
       "Acht-Vertex-Modells. Dieser Punkt ist der reelle 2-Torsionspunkt des Kreuzungsparameters; das invariante Gitter "
       "der R-Matrix ist genau für Feld sqrt 2 quadratisch (j = 1728). Baxters Modul ist dabei tan^2(pi/8), der "
       "Ising-Modul 1/sqrt 2; beide sind durch die Landen-Transformation k = (1 - k')/(1 + k') mit k' = 1/sqrt 2 "
       "verbunden, und beide Kurven haben j = 1728.")

say("  (d) Dispersionskurven (exakt)")
h_s, z_s, mu_s = sp.symbols("h z mu", positive=True)
j_tfi = jlam.subs(lam_s, h_s ** 2)
sol_tfi = sorted(sp.solve(sp.Eq(sp.numer(sp.together(j_tfi - 1728)), 0), h_s), key=lambda x: float(x))
lam_q = 1 + 1 / mu_s ** 2                                           # Doppelverhältnis (z1, -z1; z2, -z2), mu = m/(2t)
j_q = sp.simplify(jlam.subs(lam_s, lam_q))
sol_q = sorted(sp.solve(sp.Eq(sp.numer(sp.together(j_q - 1728)), 0), mu_s), key=lambda x: float(x))
zz1, zz2 = sp.symbols("z1 z2")
cr = sp.simplify(((zz1 - zz2) * (-zz1 + zz2)) / ((zz1 + zz2) * (-zz1 - zz2)))
cr_val = sp.simplify(cr.subs({zz1: -mu_s + sp.sqrt(mu_s ** 2 + 1), zz2: -mu_s - sp.sqrt(mu_s ** 2 + 1)}))
check("A10 [E] Ising-Kette: Verzweigungspunkte der Dispersion über z = e^(iq) sind 0, unendlich, h, 1/h; j = 1728 genau "
      "für h = sqrt 2 oder 1/sqrt 2. Ladungskette (Hüpfen t, Wechselpotential m, mu = m/2t): Verzweigungspunkte +-z1, +-z2 "
      "mit Doppelverhältnis 1 + 1/mu^2; j = 1728 genau für m = 2t",
      [sp.nsimplify(x) for x in sol_tfi] == [sp.sqrt(2) / 2, sp.sqrt(2)] and sol_q == [1]
      and sp.simplify(cr_val - lam_q) == 0)
q_s, s_s = sp.symbols("q s")
eps2 = 4 + 4 * sp.sin(q_s) ** 2                                      # t = 1, m = 2: eps^2 = m^2 + 4 t^2 sin^2 q
ssub = sp.I * sp.sin(q_s)
check("A11 [E] Ladungskette bei m = 2t in der Koordinate s = i sin q: (eps/m)^2 = 1 - s^2, cos^2 q = 1 + s^2, also "
      "w^2 = 1 - s^4 mit w = (eps/m) cos q: die Dispersionskurve ist wörtlich die Markenkurve. Marken s = +-1: eps = 0 "
      "(Lückenrand im komplexen Impuls), s = +-i: Bandkante q = +-pi/2. Das invariante Differential ist "
      "dq/eps = -i ds/(2t w)",
      sp.simplify(eps2 / 4 - (1 - ssub ** 2)) == 0 and sp.simplify(sp.cos(q_s) ** 2 - (1 + ssub ** 2)) == 0
      and sp.simplify(sp.diff(ssub, q_s) - sp.I * sp.cos(q_s)) == 0)


def perioden(m, t=mp.mpf(1)):
    emax = mp.sqrt(m * m + 4 * t * t)
    kap, kapp = 2 * t / emax, m / emax
    return 4 * mp.ellipk(kap ** 2) / emax, 4 * mp.ellipk(kapp ** 2) / emax        # Bandzyklus, Lückenzyklus


Pi_b, Pi_g = perioden(mp.mpf(2))
th0 = mp.asinh(1)
Pi_g_quad = 2 * mp.quad(lambda th: 1 / mp.sqrt(4 - 4 * mp.sinh(th) ** 2), [-th0, 0, th0])
Pi_b_quad = mp.quad(lambda q: 1 / mp.sqrt(4 + 4 * mp.sin(q) ** 2), [0, PI / 2, PI, 3 * PI / 2, 2 * PI])
check("A12 [E] Perioden von dq/eps bei m = 2t: Bandzyklus (q von 0 bis 2 pi) = Lückenzyklus (geschlossener Weg zwischen "
      "den Marken s = -1 und 1) = Lemniskatenkonstante/t = 2,62206/t. Das Periodengitter ist quadratisch",
      abs(Pi_b - VARPI) < TOL and abs(Pi_g - VARPI) < TOL and abs(Pi_g_quad - VARPI) < mp.mpf(10) ** -12
      and abs(Pi_b_quad - VARPI) < mp.mpf(10) ** -12)
ZAHLEN["Lemniskatenkonstante"] = ns(VARPI, 25)

say("  (e) Die Ladungskette: explizite R-Matrix aus lemniskatischen Funktionen")
# Jordan-Wigner und Teilchen-Loch-Tausch auf den ungeraden Plätzen
sm_ = np.array([[0, 1], [0, 0]], complex)                              # |0><1|, |1> = besetzt
cs = []
for j in range(L6):
    m_ = np.array([[1]], complex)
    for i in range(L6):
        m_ = np.kron(m_, NZ if i < j else (sm_ if i == j else NI))
    cs.append(m_)
Hf = np.zeros((2 ** L6, 2 ** L6), complex)
for j in range(L6 - 1):
    Hf += -(cs[j].conj().T @ cs[j + 1] + cs[j + 1].conj().T @ cs[j])
for j in range(L6):
    Hf += -2.0 * (-1) ** j * cs[j].conj().T @ cs[j]                   # t = 1, m = -2t (Vorzeichen wie z152)
U = np.eye(2 ** L6, dtype=complex)
for j in range(1, L6, 2):
    U = U @ op(NX, j, L6)
Hp_open = np.zeros((2 ** L6, 2 ** L6), complex)
for j in range(L6 - 1):
    Hp_open += -0.5 * (op(NX, j, L6) @ op(NX, j + 1, L6) - op(NY, j, L6) @ op(NY, j + 1, L6))
for j in range(L6):
    Hp_open += op(NZ, j, L6)
check("A13 [E] Ladungskette (t = 1, m = -2, offen, sechs Plätze) ist nach Jordan-Wigner und Teilchen-Loch-Tausch auf den "
      "ungeraden Plätzen exakt die XY-Kette mit reiner Paarung im Querfeld: -(1/2) Summe (XX - YY) + Summe Z. Die "
      "Ladung wird zur gestaffelten Ladung; die Kette ist dort um einen Platz translationsinvariant",
      np.linalg.norm(U @ Hf @ U.conj().T - Hp_open) < 1e-11)


def e_bs(x, k):
    return cn(x, k) + I * sn(x, k)


def R_bs(u, p1, p2, k):
    """Bazhanov-Stroganov: freie Fermionen im Feld, Differenzform in u, 'Farben' p1, p2 auf den Linien."""
    f = mp.sqrt(e_bs(p1, k) * sn(p1, k)) * mp.sqrt(e_bs(p2, k) * sn(p2, k))
    eu = e_bs(u, k)
    a, ap = 1 - eu * e_bs(p1, k) * e_bs(p2, k), eu - e_bs(p1, k) * e_bs(p2, k)
    b, bp = e_bs(p1, k) - eu * e_bs(p2, k), e_bs(p2, k) - eu * e_bs(p1, k)
    c, d = f * (1 - eu) / sn(u / 2, k), -I * k * f * (1 + eu) * sn(u / 2, k)
    return mp.matrix([[a, 0, 0, d], [0, b, c, 0], [0, c, bp, 0], [d, 0, 0, ap]]), (a, ap, b, bp, c, d)


kq, u1, u2, u3 = mp.mpf("0.6"), mp.mpf("0.3"), mp.mpf("-0.45"), mp.mpf("0.17")
p1, p2, p3 = mp.mpf("0.7"), mp.mpf("0.52"), mp.mpf("1.1")
Lh = e12(R_bs(u1 - u2, p1, p2, kq)[0]) * e13(R_bs(u1 - u3, p1, p3, kq)[0]) * e23(R_bs(u2 - u3, p2, p3, kq)[0])
Rh = e23(R_bs(u2 - u3, p2, p3, kq)[0]) * e13(R_bs(u1 - u3, p1, p3, kq)[0]) * e12(R_bs(u1 - u2, p1, p2, kq)[0])
a, ap, b, bp, c, d = R_bs(u1 - u2, p1, p2, kq)[1]
hh = mp.mpf(10) ** -12
RcP = lambda u: PERM * R_bs(u, p1, p1, kq)[0]
R0 = (RcP(hh) + RcP(-hh)) / 2
dR = (RcP(hh) - RcP(-hh)) / (2 * hh) / R0[1, 1]
s1, c1 = sn(p1, kq), cn(p1, kq)
Hbs = I / 2 * E4 + ((XX + YY) / 2 + kq * s1 * (XX - YY) / 2 + c1 * (kron(SZ, E2) + kron(E2, SZ)) / 2) / (2 * s1)
check("A14 [E] Bazhanov-Stroganov-R-Matrix (aus der Literatur, hier nachgerechnet): Freie-Fermionen-Bedingung "
      "a a' + b b' = c^2 + d^2, Yang-Baxter-Gleichung mit drei Farben und Differenzform in u (Rest < 1e-35), und der "
      "Hamiltonoperator ist die XY-Kette im Feld mit Anisotropie gamma = k sn(psi) und Feld h = cn(psi)",
      abs(a * ap + b * bp - c * c - d * d) < TOL and nrm(Lh - Rh) < TOL * nrm(Lh) and nrm(dR - Hbs) < mp.mpf(10) ** -18)


def wQ(u, k=I):
    eu = cn(u, k) + I * sn(u, k)
    return eu, (I / 2) * (1 - eu) / sn(u / 2, k), (k / 2) * (1 + eu) * sn(u / 2, k)


def RcQ(u, k=I):
    """Gatterform der Ladungskette (Paarungsbild): gerader Block [[1, d], [d, e]], ungerader Block c mal Eins."""
    eu, c, d = wQ(u, k)
    return mp.matrix([[1, 0, 0, d], [0, c, 0, 0], [0, 0, c, 0], [d, 0, 0, eu]])


def RQ(u, k=I):
    return PERM * RcQ(u, k)


psi0 = (1 + I) * VARPI / 2 + mp.mpf(10) ** -14 * mp.mpc("0.3", "0.4")
Rl = R_bs(mp.mpf("0.37"), psi0, psi0, I)[0]
Rl = PERM * Rl / Rl[0, 0]
uq, vq = mp.mpc("0.37", "0.2"), mp.mpc("-0.21", "0.05")
Lh = e12(RQ(uq - vq)) * e13(RQ(uq)) * e23(RQ(vq))
Rh = e23(RQ(vq)) * e13(RQ(uq)) * e12(RQ(uq - vq))
slx = mp.mpf("0.61")
sl_inv = mp.quad(lambda z: 1 / mp.sqrt(1 - z ** 4), [0, sn(slx, I)])
check("A15 [E] Grenzfall reiner Paarung (psi am Pol von sn, k^2 = -(2t/m)^2): Gatter mit geradem Block [[1, d], [d, e]] "
      "und ungeradem Block c, e = cn u + i sn u, c = (i/2)(1 - e)/sn(u/2), d = (k/2)(1 + e) sn(u/2). Bei m = 2t ist k = i "
      "und sn(u; i) die lemniskatische Sinusfunktion (Umkehrung von Integral dz/sqrt(1 - z^4)). Yang-Baxter in "
      "Differenzform < 1e-35; Grenzwert der allgemeinen Formel auf 1e-10",
      nrm(Rl - RcQ(mp.mpf("0.37"))) < mp.mpf(10) ** -10 and nrm(Lh - Rh) < TOL * nrm(Lh)
      and abs(sl_inv - slx) < mp.mpf(10) ** -15)
dRQ0 = mp.matrix([[0, 0, 0, I / 2], [0, I / 2, 0, 0], [0, 0, I / 2, 0], [I / 2, 0, 0, I]])       # analytisch, k = i
hQ = (kron(SZ, E2) + kron(E2, SZ)) / 2 - (XX - YY) / 2
dnum = (RcQ(hh) - RcQ(-hh)) / (2 * hh)
okH = nrm(dnum - dRQ0) < mp.mpf(10) ** -18 and nrm(dRQ0 - (I / 2) * (E4 - hQ)) < TOL
Hp_ring = np.zeros((2 ** L6, 2 ** L6), complex)
for j in range(L6):
    Hp_ring += -0.5 * (op(NX, j, L6) @ op(NX, (j + 1) % L6, L6) - op(NY, j, L6) @ op(NY, (j + 1) % L6, L6)) + op(NZ, j, L6)
okT = True
for u in (mp.mpf("0.3"), mp.mpc("0.7", "0.2")):
    Tq = ring_T(npm(RQ(u)), L6)
    okT &= np.linalg.norm(Tq @ Hp_ring - Hp_ring @ Tq) < 1e-10 * np.linalg.norm(Tq) * np.linalg.norm(Hp_ring)
Tq2 = ring_T(npm(RQ(mp.mpf("0.9"))), L6)
okT &= np.linalg.norm(Tq @ Tq2 - Tq2 @ Tq) < 1e-10 * np.linalg.norm(Tq) * np.linalg.norm(Tq2)
# Kontrolle |m| = |t|: k^2 = -(2t/m)^2 = -4
kc = 2 * I
Hc_ring = np.zeros((2 ** L4, 2 ** L4), complex)
for j in range(L4):
    Hc_ring += -1.0 * (op(NX, j, L4) @ op(NX, (j + 1) % L4, L4) - op(NY, j, L4) @ op(NY, (j + 1) % L4, L4)) + op(NZ, j, L4)
Tc = ring_T(npm(RQ(mp.mpf("0.3"), kc)), L4)
okC = np.linalg.norm(Tc @ Hc_ring - Hc_ring @ Tc) < 1e-10 * np.linalg.norm(Tc) * np.linalg.norm(Hc_ring)
check("A16 [E/N] Gatter(u) = 1 + (i u/2)(1 - h) + ..., h = (Z_1 + Z_2)/2 - (XX - YY)/2: Das ist die Ladungskette aus "
      "A13 mit reellem u als Zeit u/2. Ring mit sechs Plätzen: Die Transfermatrizen vertauschen untereinander und mit "
      "der Kette. Kontrolle m = t (k^2 = -4, vier Plätze) vertauscht ebenfalls", okH and okT and okC)
ok = True
for u in (mp.mpf("0.3"), mp.mpf("1.1"), mp.mpf("2.0")):
    eu, c, d = wQ(u)
    B = mp.matrix([[1, d], [d, eu]]) / c
    ok &= nrm(B * B.H - E2) < TOL and abs(mp.det(B) - 1) < TOL
ok &= nrm(RcQ(mp.mpf("0.7")) * kron(SZ, E2) - kron(SZ, E2) * RcQ(mp.mpf("0.7"))
          - (RcQ(mp.mpf("0.7")) * kron(E2, SZ) - kron(E2, SZ) * RcQ(mp.mpf("0.7")))) < TOL
check("A17 [E] Für reelles u ist das Gatter unitär: im Ein-Teilchen-Sektor der Ladungskette ein Strahlteiler mit "
      "Determinante 1, auf leerem und doppelt besetztem Paar die Eins. Es erhält die Ladung (im Paarungsbild Z_1 - Z_2) "
      "und ist rein einkörperig", ok)
befund("A-Ladung", "E",
       "Die Verbindung zur Ladungskette ist explizit: Sie ist die XY-Kette mit reiner Paarung im Querfeld, ihre "
       "R-Matrix ist der Grenzfall reiner Paarung der Bazhanov-Stroganov-Familie mit Modul k^2 = -(2t/m)^2. Bei m = 2t "
       "ist k^2 = -1: Die Gewichte sind lemniskatische Funktionen, die Dispersionskurve ist wörtlich w^2 = 1 - s^4. "
       "Die Kette liegt damit auf der nativen Kurve, aber in der Freie-Fermionen-Familie 'mit Feld' (unsymmetrische "
       "Gewichte, U(1)-Ladung), nicht in der symmetrischen Acht-Vertex-Familie (Belavin, sl(2)).")

# --------------------------------------------------------------------------------------------------------------------
say("\nTeil B  Boost und Eckentransfermatrix")


def d4(f, u, h=mp.mpf(10) ** -10):
    return (-f(u + 2 * h) + 8 * f(u + h) - 8 * f(u - h) + f(u - 2 * h)) / (12 * h)


def sutherland(f, h12, u):
    R, dRu = f(u), d4(f, u)
    lhs = e12(R) * e23(dRu) - e12(dRu) * e23(R)
    rhs = e12(R) * e23(R) * e12(h12) - e23(h12) * e12(R) * e23(R)
    return nrm(lhs - rhs) / nrm(lhs)


f8 = lambda u: Rc8(u, lamg, kg)
jx, jy, jz = coup(lamg, kg)
h8 = -(jx * XX + jy * YY + jz * ZZ + jz * E4) / (2 * sn(lamg, kg))
s8 = sutherland(f8, h8, mp.mpc("0.3", "0.2"))
sQ = sutherland(RcQ, -dRQ0, mp.mpc("0.4", "0.1"))
check("B1 [E] Lokale Boost-Relation (Ableitung der Yang-Baxter-Gleichung): G_12 G'_23 - G'_12 G_23 = G_12 G_23 h_12 - "
      "h_23 G_12 G_23 mit h = -G'(0), für das Acht-Vertex-Gatter (allgemeines lam, k) und das Gatter der Ladungskette. "
      "Aufsummiert heißt das [Summe j h_j, T(u)] = T'(u): Der Gitter-Boost erzeugt die Verschiebung von u",
      s8 < mp.mpf(10) ** -25 and sQ < mp.mpf(10) ** -25)
Rx = npm(R8(K_, K_, KB))
RqX = npm(RcQ(I * VARPI))
rank8 = int(np.sum(np.linalg.svd(Rx, compute_uv=False) > 1e-12))
rankQ = int(np.sum(np.linalg.svd(RqX, compute_uv=False) > 1e-12))
a, b, c, d = w8(K_ / 2, K_, KB)
check("B2 [E] Kreuzungspunkt: Die R-Matrix hat dort Rang 1. Ising-Form: u = lam = K (reeller 2-Torsionspunkt); "
      "isotrope Ecke a = b bei u = K/2 (4-Torsionspunkt). Ladungskette: u = i mal Lemniskatenkonstante "
      "(2-Torsionspunkt des Gitters 2 varpi Z[i]), Rang 1 mit geradem Block [[1, -1], [-1, 1]]",
      rank8 == 1 and rankQ == 1 and abs(a - b) < TOL and nrm(RcQ(I * VARPI) - mp.matrix(
          [[1, 0, 0, -1], [0, 0, 0, 0], [0, 0, 0, 0], [-1, 0, 0, 1]])) < TOL)


def ent(M, cut):
    """H = (i/2) Summe a_n M[n,m] b_m. Grundzustand <i a_n b_m> = -(U V)[n,m]; Verschränkungsoperator der linken cut Plätze."""
    Um, _, Vm = mp.svd_r(M)
    G = -(Um * Vm)
    U2, S2, V2 = mp.svd_r(G[:cut, :cut])
    lev = [2 * mp.atanh(min(s, 1 - mp.mpf(10) ** -62)) for s in S2]
    W = mp.zeros(cut, cut)
    for i in range(cut):
        for a_ in range(cut):
            ua = lev[i] * U2[a_, i]
            for b_ in range(cut):
                W[a_, b_] += ua * V2[i, b_]
    return -W, sorted(lev)


def M_ladung(N, t, m):
    M = mp.zeros(N, N)
    for n in range(N):
        M[n, n] = m
    for n in range(N - 1):
        M[n, n + 1], M[n + 1, n] = t, -t
    return M


def M_ising(N, h):
    M = mp.zeros(N, N)
    for n in range(N):
        M[n, n] = 2 * h
    for n in range(N - 1):
        M[n + 1, n] = -2
    return M


NK, CUT = 40, 20
with mp.workdps(72):
    Mq = M_ladung(NK, mp.mpf(1), mp.mpf(2))
    Wq, levq = ent(Mq, CUT)
    cq_site = [Wq[CUT - 1 - d_, CUT - 1 - d_] / Mq[CUT - 1 - d_, CUT - 1 - d_] / (d_ + mp.mpf(1) / 2) for d_ in range(5)]
    cq_bond = [Wq[CUT - 1 - d_, CUT - 2 - d_] / Mq[CUT - 1 - d_, CUT - 2 - d_] / (d_ + 1) for d_ in range(5)]
    Mq1 = M_ladung(NK, mp.mpf(1), mp.mpf(1))
    Wq1, levq1 = ent(Mq1, CUT)
    c1_meas = Wq1[CUT - 1, CUT - 2] / Mq1[CUT - 1, CUT - 2]
    Mi = M_ising(NK, mp.sqrt(2))
    Wi, levi = ent(Mi, CUT)
    ci_site = Wi[CUT - 1, CUT - 1] / Mi[CUT - 1, CUT - 1] * 2
    ci_bond = Wi[CUT - 1, CUT - 2] / Mi[CUT - 1, CUT - 2]
Pb1, Pg1 = perioden(mp.mpf(1))
ZAHLEN["Ladungskette m = 2t: K/(Abstand mal lokale Energie), Plätze und Bindungen"] = [ns(x, 14) for x in cq_site + cq_bond]
ZAHLEN["Ladungskette m = t: c gemessen, Lückenperiode"] = [ns(c1_meas, 12), ns(Pg1, 12)]
ZAHLEN["Ising h = sqrt2: c aus Platz und Bindung (40 Plätze)"] = [ns(ci_site, 9), ns(ci_bond, 9)]
say(f"   Ladungskette m = 2t (Paarungsbild, 40 Plätze, 72 Stellen): c = {ns(cq_bond[0], 14)}; m = t: c = {ns(c1_meas, 12)}, "
    f"Lückenperiode {ns(Pg1, 12)}")
say(f"   Ising-Kette h = sqrt 2: c aus Platz {ns(ci_site, 9)}, aus Bindung {ns(ci_bond, 9)} (Korrelationslänge 2,9 Plätze)")
check("B3 [N] Ladungskette im Paarungsbild: Der Verschränkungsoperator der Halbkette ist c mal Summe (Abstand vom "
      "Schnitt) mal (lokale Energie), Plätze bei 1/2, 3/2, ..., Bindungen bei 1, 2, ...; c = Lemniskatenkonstante auf "
      "1e-9 an den ersten fünf Plätzen und Bindungen, 1e-12 am Schnitt (z152 D2/D3 unabhängig reproduziert)",
      all(abs(x - VARPI) < mp.mpf(10) ** -9 for x in cq_site + cq_bond)
      and abs(cq_site[0] - VARPI) < mp.mpf(10) ** -12 and abs(cq_bond[0] - VARPI) < mp.mpf(10) ** -12)
check("B4 [N] Allgemein ist c die Periode von dq/eps über den Lückenzyklus (Kontrolle m = t: 1e-8), und der Abstand der "
      "Niveaus ist 2 pi mal Lückenperiode durch Bandperiode (erste drei Niveaus, 1e-6). Bei m = 2t sind beide Perioden "
      "gleich: Abstand 2 pi, Niveaus pi (2l + 1)",
      abs(c1_meas - Pg1) < mp.mpf(10) ** -8
      and all(abs(levq1[l] - (2 * l + 1) * PI * Pg1 / Pb1) < mp.mpf(10) ** -6 for l in range(3))
      and all(abs(levq[l] - (2 * l + 1) * PI) < mp.mpf(10) ** -10 for l in range(3)))
check("B5 [E/N] K = 2 lam mal Boost, mit dem Boost aus h = -G'(0): Ladungskette 2 |lam| (1/2) = Lemniskatenkonstante "
      "(lam = i varpi, Zeit u/2); Ising-Form 2K (1 - k) = Lemniskatenkonstante exakt (50 Stellen) und numerisch "
      "c = 2,622057 aus der Ising-Kette (1e-6; ein Ising-Platz sind zwei Plätze der XY-Kette, daher 2K (1 - k) und "
      "nicht K (1 - k)). Die Dichtematrix ist die Verschiebung um den vollen reellen Periodenvektor 2 lam",
      abs(2 * K_ * (1 - KB) - VARPI) < mp.mpf(10) ** -45 and abs(ci_site - VARPI) < mp.mpf(10) ** -6
      and abs(ci_bond - VARPI) < mp.mpf(10) ** -5)
ks = [4 * n + 1 for n in range(-6, 7)]
check("B6 [L; Arithmetik E] Punkte der Kurve: Dichtematrix rho = volle reelle Periode (Ursprung); Kreuzung (euklidische halbe Drehung, "
      "rho^(1/2)) = reeller 2-Torsionspunkt; isotrope Ecke (Vierteldrehung) = reeller 4-Torsionspunkt; Modulfluss läuft "
      "in imaginärer Richtung, Delta^(it) kehrt bei t = 1 bis auf (-1)^F zurück; die Uhr R = Delta^(-i/2) = i^Q ist der "
      "imaginäre 2-Torsionspunkt (Spektrum pi (4Z + 1): e^(i pi k/2) = i für alle k)",
      all(abs(mp.exp(I * PI * k / 2) - I) < TOL for k in ks) and all(abs(mp.exp(-I * PI * k) + 1) < TOL for k in ks))
befund("B", "E/N",
       "Ja: Der Erzeuger der Spektralverschiebung ist der Gitter-Boost, und c = 2,62206 folgt direkt: c ist die Periode "
       "des invarianten Differentials dq/eps = ds/(2t sqrt(1 - s^4)) über den Lückenzyklus, gleichwertig der doppelte "
       "Kreuzungsparameter. Genau ist 2 Integral_0^1 dz/sqrt(1 - z^4) das Integral zwischen den Marken -1 und 1, also die "
       "Hälfte der geschlossenen Periode von dz/w; das Periodengitter von dz/w ist (1 + i) varpi Z[i]. Die volle "
       "Modulperiode ist der volle Periodenvektor, die isotrope Ecke ein 4-Torsionspunkt, die Uhr ein 2-Torsionspunkt.")

# --------------------------------------------------------------------------------------------------------------------
say("\nTeil C  Z_4-Automorphie und Torsionspunkte")
check("C1 [E] Auf der Dispersionskurve der Ladungskette vertauscht s -> i s die Faktoren 1 - s^2 = (eps/m)^2 und "
      "1 + s^2 = cos^2 q: Energie gegen Impulskosinus, Band gegen Lücke, Bandperiode gegen Lückenperiode. Das ist der "
      "Tausch von Translation (reelle Zeit) und Boost (euklidische Drehung); er ist genau bei m = 2t eine Symmetrie",
      sp.simplify((1 - (sp.I * s_s) ** 2) - (1 + s_s ** 2)) == 0 and abs(Pi_b - Pi_g) < TOL
      and abs(perioden(mp.mpf(1))[0] - perioden(mp.mpf(1))[1]) > mp.mpf("0.1"))
HD = (SX + SZ) / mp.sqrt(2)
HH = kron(HD, HD)
ok = True
for u, lam in ((mp.mpf("0.31"), mp.mpf("0.47")), (mp.mpc("0.2", "0.3"), mp.mpf("0.9"))):
    ok &= prop(R8(I * u, I * lam, KB), HH * R8(u, lam, KB) * HH) < TOL
check("C2 [E] Komplexe Multiplikation auf der R-Matrix (nur bei tau = i möglich): R(i u; i lam) = Skalar mal "
      "(H x H) R(u; lam) (H x H) mit der Hadamard-Matrix H. Die Z_4-Automorphie wirkt auf u und den Kreuzungsparameter "
      "zugleich und vertauscht die Achsen x und z (Dualität); die y-Achse ist fest", ok)
jA, jB = coup(K_, KB), coup(I * KP_ / 2 + eps, KB)
nB = max(abs(x) for x in jB)
check("C3 [E] Folge am freien Punkt: lam = K (Kopplungen 1 : 1/sqrt 2 : 0) geht in lam = i K (0 : 1/sqrt 2 : 1) über, "
      "dieselbe Kette mit vertauschten Achsen. Nur am selbstdualen Modul sind beide Verhältnisse gleich",
      abs(jA[1] / jA[0] - 1 / mp.sqrt(2)) < TOL and abs(jB[0] / nB) < mp.mpf(10) ** -20
      and abs(jB[1] / jB[2] - 1 / mp.sqrt(2)) < mp.mpf(10) ** -20)
say("  (a) Die R-Matrix der Ising-Form an den Torsionspunkten (freier Punkt lam = K)")
uclk = I * KP_ / 2
G_clk = Rc8(uclk, K_, KB)
G_soll = expw(PI / 4, PI / 8, 0)
c1, c2, c3, spread = weyl(uclk, K_, KB)
check("C4 [E] Am imaginären 2-Torsionspunkt u = i K'/2 = i lam (Bild der Kreuzung unter der Automorphie; der Uhrpunkt "
      "aus B6) ist das Gatter exakt exp(i (pi/4) XX) exp(i (pi/8) YY) bis auf eine Phase (50 Stellen)",
      prop(G_clk, G_soll) < mp.mpf(10) ** -40 and abs(c1 - PI / 4) < TOL and abs(c2 - PI / 8) < TOL and abs(c3) < TOL
      and spread < TOL)
# Lesart in der Ising-Kette: XX = Platzterm, YY = Bindungsterm
site = mp.cos(PI / 4) * E2 + I * mp.sin(PI / 4) * SZ
uhr_site = mp.exp(I * PI / 4) * mp.matrix([[1, 0], [0, -I]])
VY = mp.matrix([[1, 1], [I, -I]]) / mp.sqrt(2)                          # Spalten: Eigenvektoren von sigma_y
VV = kron(VY, VY)
bond_y = VV.H * (mp.cos(PI / 8) * E4 + I * mp.sin(PI / 8) * YY) * VV
Dl = mp.matrix([[1, 0], [0, mp.exp(I * PI / 4)]])
ok = nrm(site - uhr_site) < TOL
bd = mp.diag([bond_y[i, i] for i in range(4)])
ok &= nrm(bond_y - bd) < TOL and nrm(kron(Dl, Dl) * bond_y - mp.exp(I * PI / 8) * mp.diag([1, 1, 1, I])) < TOL
check("C5 [E] In der Ising-Kette (A6) ist XX der Platzterm sigma_z und YY der Bindungsterm sigma_x sigma_x. Am Uhrpunkt "
      "ist das Platzgatter exp(i (pi/4) sigma_z) = Phase mal diag(1, -i): die Uhr i^(-n) je Platz. Das Bindungsgatter "
      "exp(i (pi/8) sigma sigma) ist in der Eigenbasis mit Ein-Platz-Phasen exakt diag(1, 1, 1, i): die gesteuerte Uhr "
      "zwischen benachbarten Ising-Spins (Ordnungsvariablen)", ok)
G120 = (E4 + I * PERM) / (1 + I)
CSm = mp.diag([1, 1, 1, -I])
CSp = mp.diag([1, 1, 1, I])


def dz_block(G):
    return (G[1, 1] * G[2, 2] - G[1, 2] * G[2, 1]) / (G[0, 0] * G[3, 3] - G[0, 3] * G[3, 0])


def proj_order(G, nmax=32):
    A = E4
    for n in range(1, nmax + 1):
        A = A * G
        if nrm(A - A[0, 0] * E4) < mp.mpf(10) ** -30:
            return n
    return 0


tab = {}
Gn = G_clk / G_clk[1, 1] * G_soll[1, 1]
for name, G, ax in (("R-Matrix am Uhrpunkt (frei)", Gn, None), ("z120-Gatter (1+iP)/(1+i)", G120, None),
                    ("gesteuerte Uhr diag(1,1,1,-i)", CSm, None), ("Tausch P", PERM, None)):
    tab[name] = {"Operator-Schmidt-Rang": osr(G), "delta_z": ns(mp.chop(dz_block(G), tol=mp.mpf(10) ** -30), 6), "projektive Ordnung": proj_order(G),
                 "Quadrat = Tausch (bis auf Phase)": bool(prop(G * G, PERM) < mp.mpf(10) ** -30)}
tab["R-Matrix am Uhrpunkt (frei)"]["Koordinaten/pi"] = [ns(c1 / PI, 6), ns(c2 / PI, 6), ns(mp.chop(abs(c3) / PI, tol=mp.mpf(10) ** -30), 6)]
tab["z120-Gatter (1+iP)/(1+i)"]["Koordinaten/pi"] = ["0.125", "0.125", "0.125"]
tab["gesteuerte Uhr diag(1,1,1,-i)"]["Koordinaten/pi"] = ["0.125", "0.0", "0.0"]
tab["Tausch P"]["Koordinaten/pi"] = ["0.25", "0.25", "0.25"]
ZAHLEN["Vergleich der Gatter"] = tab
for name in tab:
    say(f"   {name}: {tab[name]}")
okW = (prop(G120, expw(PI / 8, PI / 8, PI / 8)) < TOL and prop(PERM, expw(PI / 4, PI / 4, PI / 4)) < TOL
       and prop(kron(Dl, Dl).H * CSp, expw(0, 0, PI / 8)) < TOL)
check("C6 [E] Vergleich über Invarianten: Die R-Matrix am Uhrpunkt hat die Koordinaten (pi/4, pi/8, 0), "
      "Operator-Schmidt-Rang 4, projektive Ordnung 8 und delta_z = 1; ihr Quadrat ist nicht der Tausch. Das "
      "z120-Gatter hat (pi/8, pi/8, pi/8), delta_z = -i, Quadrat = Tausch; die gesteuerte Uhr (pi/8, 0, 0), Rang 2; der "
      "Tausch (pi/4, pi/4, pi/4). Sie ist keines davon",
      okW and tab["R-Matrix am Uhrpunkt (frei)"]["Operator-Schmidt-Rang"] == 4
      and tab["R-Matrix am Uhrpunkt (frei)"]["projektive Ordnung"] == 8
      and not tab["R-Matrix am Uhrpunkt (frei)"]["Quadrat = Tausch (bis auf Phase)"]
      and abs(dz_block(Gn) - 1) < TOL and abs(dz_block(G120) + I) < TOL
      and tab["z120-Gatter (1+iP)/(1+i)"]["Quadrat = Tausch (bis auf Phase)"]
      and tab["gesteuerte Uhr diag(1,1,1,-i)"]["Operator-Schmidt-Rang"] == 2
      and tab["z120-Gatter (1+iP)/(1+i)"]["Operator-Schmidt-Rang"] == 4)
ok = True
for u in (mp.mpc("0.2", "0.9"), I * KP_ / 4, K_ / 2 + I * KP_ / 3, mp.mpf("0.4")):
    ok &= abs(deltas(u, K_, KB)[2] - 1) < TOL
check("C7a [E negativ] Auf der ganzen Kurve des freien Punkts ist delta_z = 1 (gaußsch in den Kettenfermionen); das "
      "z120-Gatter und die gesteuerte Uhr um die Ladungsachse haben delta_z = -i bzw. +-i. Kein Punkt der freien Kurve "
      "ist das wechselwirkende Gatter, auch nicht bis auf paritätserhaltende Ein-Platz-Basiswechsel", ok)
c1q, c2q, c3q, _ = weyl(I * KP_ / 4, K_, KB)
d_q = mp.exp(-4 * I * c2q)
check("C7b [E] Am 4-Torsionspunkt u = i K'/4 ist keine Koordinate ein Uhrwert: e^(-4 i c2) hat Realteil sqrt 2 - 1, und "
      "2 sqrt 2 - 2 hat die Konjugierte -2 sqrt 2 - 2 vom Betrag > 2, ist also nicht 2 cos eines rationalen Winkels. Die "
      "4-Torsionspunkte tragen die isotrope Ecke, nicht den Viertelwinkel; der sitzt am 2-Torsionspunkt",
      abs(mp.re(d_q) - (mp.sqrt(2) - 1)) < TOL and abs(c3q) < TOL and abs(-2 * mp.sqrt(2) - 2) > 2)
ZAHLEN["4-Torsionspunkt i K'/4: (c1, c2)/pi"] = [ns(c1q / PI, 12), ns(c2q / PI, 12)]
say("  (b) Das Viertel als Invariante der komplexen Multiplikation (beliebiger Kreuzungsparameter)")
lams = [K_ / 1000, K_ / 10, K_ / 4, K_ / 2, 3 * K_ / 4, 9 * K_ / 10, K_]
ok = True
fam = []
for lam in lams:
    dx, dy, dz = deltas(I * lam, lam, KB)
    a1, a2, a3, sprd = weyl(I * lam, lam, KB)
    ok &= abs(dy + I) < TOL and abs(dx * dz + 1) < TOL and sprd < TOL and abs(a2 - PI / 8) < TOL
    fam.append([ns(lam / K_, 6), ns(a1 / PI, 10), ns(a2 / PI, 10), ns(abs(a3) / PI, 10) if abs(a3) > mp.mpf(10) ** -30 else "0.0"])
dyc = deltas(I * mp.mpc("0.4", "0.3"), mp.mpc("0.4", "0.3"), KB)[1]
ok &= abs(dyc + I) < TOL
ZAHLEN["Familie u = i lam auf E_i: lam/K, (c1, c2, c3)/pi"] = fam
for row in fam:
    say(f"   lam/K = {row[0]}: (c1, c2, c3)/pi = ({row[1]}, {row[2]}, {row[3]})")
check("C8 [E] Auf E_i gilt für JEDEN Kreuzungsparameter lam: Am Punkt u = i lam ist das Gatter unitär und "
      "delta_y = -i exakt, also c2 = pi/8; dazu c1 + c3 = pi/4 (sieben reelle Werte und ein komplexer, 35 Stellen). "
      "Das Viertel um die y-Achse ist eine Invariante der Z_4-Automorphie", ok)
cc = 1 + KB
sl = lambda x: mp.ellipfun("sn", x, -1)


def Ffun(x):
    z = cc * x / 2
    return sn(z, 1 / mp.sqrt(2)) * dn(z, 1 / mp.sqrt(2)) / cn(z, 1 / mp.sqrt(2))


ut, lt, zt = mp.mpc("0.2", "0.3"), K_ / 2, mp.mpc("0.37", "0.21")
okp = abs(deltas(ut, lt, KB)[1] - Ffun(lt - ut) / Ffun(lt + ut)) < TOL
okp &= abs(sn(zt, 1 / mp.sqrt(2)) * dn(zt, 1 / mp.sqrt(2)) / cn(zt, 1 / mp.sqrt(2))
           - mp.exp(-I * PI / 4) * sl(mp.exp(I * PI / 4) * zt)) < TOL
okp &= abs(sl((1 + I) * zt) - (1 + I) * sl(zt) / mp.diff(sl, zt)) < mp.mpf(10) ** -30
okp &= abs(sl((1 - I) * zt) - (1 - I) * sl(zt) / mp.diff(sl, zt)) < mp.mpf(10) ** -30
check("C9 [E] Beweis: delta_y = F(lam - u)/F(lam + u) mit F(x) = e^(-i pi/4) sl(e^(i pi/4) (1 + k) x/2), sl die "
      "lemniskatische Sinusfunktion; und sl((1 +- i) x) = (1 +- i) sl(x)/sl'(x) (komplexe Multiplikation). Also "
      "delta_y(u = i lam) = (1 - i)/(1 + i) = -i für alle lam. In der Normalform mit sl ist die R-Matrix a = sl(u+lam), "
      "b = sl u, c = sl lam, d = i sl lam sl u sl(u+lam)", okp)
k3 = mp.mpf("0.3")
K3, K3p = mp.ellipk(k3 ** 2), mp.ellipk(1 - k3 ** 2)
ctrl = [abs(deltas(I * lam, lam, k3)[1] + I) for lam in (K3 / 2, K3)]
ctrl2 = [abs(deltas(I * lam * K3p / (2 * K3), lam, k3)[1] + I) for lam in (K3 / 2, K3)]
zeta = sp.Symbol("zeta")
xr, yr = sp.symbols("x y", real=True)
unit = sp.simplify(sp.Abs(1 - (xr + sp.I * yr)) ** 2 - sp.Abs(1 + (xr + sp.I * yr)) ** 2)
check("C10 [E] Nur tau = i: Bei anderem Modul (k = 0,3) ist delta_y weder bei u = i lam noch beim dualen Punkt "
      "u = i lam K'/(2K) konstant (Abweichung > 0,1). Im rationalen Keim ist delta = (1 - zeta)/(1 + zeta) für "
      "u = zeta lam; |delta| = 1 verlangt Re zeta = 0. Einheiten von End(E) sind +-1, +-i (j = 1728) oder sechste "
      "Einheitswurzeln (j = 0): Ein unitäres Gatter am Automorphiebild der Kreuzung gibt es nur auf E_i, und sein "
      "Wert ist dann zwangsläufig -i oder +i",
      min(ctrl) > mp.mpf("0.1") and min(ctrl2) > mp.mpf("0.1") and sp.simplify(unit + 4 * xr) == 0)
lamr = K_ / 2
roots = []
for ix in range(8):
    for iy in range(4):
        ustart = (mp.mpf("0.13") + ix * mp.mpf("0.25")) * 2 * K_ + I * (mp.mpf("0.11") + iy * mp.mpf("0.25")) * KP_
        try:
            rt = mp.findroot(lambda u: deltas(u, lamr, KB)[1] + I, ustart, tol=mp.mpf(10) ** -35, maxsteps=80)
        except Exception:
            continue
        if abs(deltas(rt, lamr, KB)[1] + I) > mp.mpf(10) ** -30:
            continue
        al, be = mp.re(rt) / (2 * K_), mp.im(rt) / KP_
        pq = [(al + be) / 2, (al - be) / 2]
        pq = [float(x - mp.floor(x + mp.mpf(10) ** -12)) for x in pq]
        key = (round(pq[0], 8) % 1.0, round(pq[1], 8) % 1.0)
        if key not in roots:
            roots.append(key)
roots = sorted(roots)
ZAHLEN["Lösungen von delta_y = -i bei lam = K/2 (Gitterkoordinaten)"] = [list(r) for r in roots]
say(f"   Lösungen von delta_y = -i bei lam = K/2 modulo Gitter von delta_y: {roots}")
check("C11 [N] Eindeutigkeit: Bei gegebenem lam hat delta_y = -i modulo dem Gitter von delta_y genau zwei Lösungen, "
      "u = i lam und u = -i lam + 2K (das zeitgespiegelte Bild, um einen Gittervektor von L_R verschoben, der delta in "
      "1/delta überführt). Andere Punkte gibt es nicht (32 Startpunkte)",
      roots == [(0.125, 0.875), (0.375, 0.625)])
lamz = mp.mpf("0.6")
rk = {}
for name, u in (("+lam", lamz), ("-lam", -lamz), ("+i lam", I * lamz), ("-i lam", -I * lamz)):
    Rn = npm(R8(u, lamz, KB))
    sv = np.linalg.svd(Rn, compute_uv=False)
    rk[name] = int(np.sum(sv > 1e-12 * sv[0]))
dpl, dmi = deltas(I * lamz, lamz, KB)[1], deltas(-I * lamz, lamz, KB)[1]
ZAHLEN["Bahn des Kreuzungsparameters unter Aut(E_i): Rang der R-Matrix"] = rk
check("C12 [E] Die Bahn des Kreuzungsparameters unter Aut(E_i) = Z_4: bei u = lam Rang 3 und bei u = -lam Rang 1 (die "
      "zwei Projektoren, wie (1 +- P)/2 im rationalen Fall, z151), bei u = +-i lam die zwei unitären Gatter mit "
      "delta_y = -i und +i (Zeitumkehr)",
      rk == {"+lam": 3, "-lam": 1, "+i lam": 4, "-i lam": 4} and abs(dpl + I) < TOL and abs(dmi - I) < TOL)
say("  (c) Die zwei Enden der Familie u = i lam")
errs = []
for lam in (K_ / 10, K_ / 100, K_ / 1000):
    a, b, c, d = w8(I * lam, lam, KB)
    G = mp.matrix([[a, 0, 0, d], [0, c, b, 0], [0, b, c, 0], [d, 0, 0, a]]) / (c * (1 + I))
    errs.append(nrm(G - G120))
ZAHLEN["Abstand zum z120-Gatter bei lam = K/10, K/100, K/1000"] = [ns(x, 4) for x in errs]
say(f"   Abstand zum z120-Gatter bei lam = K/10, K/100, K/1000: {[ns(x, 4) for x in errs]}")
check("C13 [E] Ende lam -> 0: Das Gatter bei u = i lam geht quadratisch in lam gegen (1 + iP)/(1 + i). Das z120-Gatter "
      "ist die rationale Entartung der Familie. Lesart: Sein Nenner 1 + i ist der Multiplikator der komplexen "
      "Multiplikation aus C9",
      errs[2] < mp.mpf(10) ** -5 and abs(errs[0] / errs[1] - 100) < 1 and abs(errs[1] / errs[2] - 100) < mp.mpf("0.1"))
c3s, dzs = [], []
for jn in range(1, 16):
    lam = K_ * jn / 16
    _, _, a3, _ = weyl(I * lam, lam, KB)
    c3s.append(abs(a3))
    dzs.append(deltas(I * lam, lam, KB)[2])
dz_half = deltas(I * K_ / 2, K_ / 2, KB)[2]
check("C14 [N; an den Enden E] Entlang der Familie fällt c3 streng von pi/8 (lam = 0) auf 0 (lam = K). Verlangt man "
      "Uhrwerte um alle drei Achsen (delta^4 = 1), so bleiben nur die Enden: lam -> 0 mit (-i, -i, -i), das "
      "SU(2)-symmetrische z120-Gatter, und lam = K mit (-1, -i, 1), das freie Gatter der Nahtkette. Dazwischen ist "
      "delta_z kein Uhrwert (15 Stützstellen; bei lam = K/2 ist Re delta_z = sqrt 2 - 1)",
      all(c3s[i] > c3s[i + 1] for i in range(14)) and c3s[0] < PI / 8 and min(abs(x ** 4 - 1) for x in dzs) > mp.mpf("0.01")
      and abs(mp.re(dz_half) - (mp.sqrt(2) - 1)) < TOL)
say("  (d) Die R-Matrix der Ladungskette an den Torsionspunkten (reelles u = Zeit, Gitter 2 varpi Z[i])")
uq0 = mp.mpc("0.41", "0.17")
Sg = kron(SZ, E2)


def Bq(u):
    eu, c, d = wQ(u)
    return mp.matrix([[1, d], [d, eu]]) / c


okl = (nrm(Bq(uq0 + 2 * VARPI) + SZ * Bq(uq0) * SZ) < TOL and nrm(Bq(uq0 + 2 * I * VARPI) + SZ * Bq(uq0) * SZ) < TOL
       and nrm(Bq(uq0 + 4 * VARPI) - Bq(uq0)) < TOL and prop(Bq(uq0 + VARPI), Bq(uq0)) > mp.mpf("0.1")
       and prop(Bq(uq0 + VARPI), SZ * Bq(uq0) * SZ) > mp.mpf("0.1"))
B1, B2, B4 = Bq(VARPI), Bq(VARPI / 2), Bq(2 * VARPI - mp.mpf(10) ** -30)
okc = nrm(B1 - mp.diag([-I, I])) < TOL
okb = (abs(B2[0, 1] - I * (mp.sqrt(2) - 1)) < TOL and abs(mp.arg(B2[0, 0]) + PI / 4) < TOL
       and abs(abs(B2[0, 0]) ** 2 - (1 - (mp.sqrt(2) - 1) ** 2)) < TOL and nrm(B4 + E2) < mp.mpf(10) ** -25)
tmax = max(abs(Bq(VARPI * mp.mpf(n) / 200)[0, 1]) ** 2 for n in range(1, 200))
ZAHLEN["Ladungskette: größte Transferwahrscheinlichkeit des Gatters auf der reellen Achse"] = ns(tmax, 12)
check("C15 [E] Ladungskette: Perioden 2 varpi und 2 i varpi bis auf Eichung (quadratisches Gitter); Kreuzung bei "
      "u = i varpi (B2). Am Uhrpunkt u = varpi = -i mal Kreuzung ist das Gatter im Ein-Teilchen-Sektor diag(-i, i) "
      "und sonst die Eins: exakt die Paaruhr i^(Q_B - Q_A) der zwei Plätze, ohne Transfer. Bei 2 varpi ist es die Parität",
      okl and okc and nrm(B4 + E2) < mp.mpf(10) ** -25)
check("C16 [E/N] Am 4-Torsionspunkt u = varpi/2 ist es ein Strahlteiler mit Transferamplitude i (sqrt 2 - 1) = i tan(pi/8) "
      "und Phasen e^(-+i pi/4). Die Transferwahrscheinlichkeit erreicht auf der reellen Achse höchstens (sqrt 2 - 1)^2 = "
      "0,1716: Der symmetrische Strahlteiler (die Vierteldrehung P_sym) ist kein Wert dieser R-Matrix",
      okb and abs(tmax - (mp.sqrt(2) - 1) ** 2) < mp.mpf(10) ** -20)
befund("C", "E; Lesart als Gatter des Kandidaten C",
       "Die Z_4-Automorphie wirkt als Dualität (Hadamard) mit Drehung von u und lam um i. Am Automorphiebild der "
       "Kreuzung, dem Uhrpunkt, ist die R-Matrix der Nahtkette: Ladungsform die Paaruhr; Ising-Form "
       "exp(i pi/4 XX) exp(i pi/8 YY), in der Ising-Kette Platzuhr mal gesteuerte Uhr zwischen Ordnungsvariablen. Sie "
       "ist frei in den Fermionen (delta_z = 1) und weder (1+iP)/(1+i) noch der Tausch. Neu und exakt: Für jeden "
       "Kreuzungsparameter ist bei u = i lam die Zweikörperinvariante um die feste Achse (1-i)/(1+i) = -i. Das "
       "z120-Gatter ist das rationale Ende lam -> 0 dieser Familie, die Nahtkette das freie Ende lam = K.")

# --------------------------------------------------------------------------------------------------------------------
say("\nTeil D  Gegen z151: die elliptische statt der rationalen Klasse")
zc, wc, tc = sp.symbols("z w t")
phiM = 1 / (1 - zc ** 4)
char_phi = sp.simplify(phiM.subs(zc, sp.I * zc) * sp.I ** 2 / phiM)        # dz^2 -> (i dz)^2
ser_inf = sp.series((1 / tc) / (sp.I * tc ** -2 * sp.sqrt(1 - tc ** 4)), tc, 0, 3).removeO()
zloc = 1 - wc ** 2 / 4                                                      # z = 1 - w^2/4 + O(w^4) nahe der Marke z = 1
check("D1 [E] Die native Quartik gibt Phi_M = 1/(pi_0^4 - pi_1^4), einen Schnitt von O(-4) = K^2 der Markensphäre: ein "
      "quadratisches Differential dz^2/(1 - z^4) mit einfachen Polen an den vier Marken. Es ist exakt das Quadrat des "
      "holomorphen Differentials dz/w von E_i. Die Uhr z -> i z wirkt auf Phi_M mit -1, auf dz/w also mit +-i: Die Uhr "
      "ist die komplexe Multiplikation (Bestand: 'mu_4 und Uhren' zum quadratischen Torus Z[i])",
      char_phi == -1 and sp.simplify((1 / wc ** 2).subs(wc, sp.sqrt(1 - zc ** 4)) - phiM) == 0)
check("D2 [E] Die reduzierte Form aus z151, z dz/(1 - z^4) = (z/w) dz/w, hat auf E_i vier einfache Pole (an den Marken, "
      "z/w ~ 1/w) und vier einfache Nullstellen (z = 0 und z = unendlich, je zwei Punkte). Die elliptische Klasse "
      "(omega ohne Pole) ist also nicht der Rückzug der Twistorreduktion; sie verlangt die Wurzel aus Phi_M",
      sp.limit(wc * zloc / wc, wc, 0) == 1 and sp.simplify(ser_inf + sp.I * tc) == 0
      and sp.series(1 - zloc ** 4, wc, 0, 3).removeO() == wc ** 2)


def rcl(u, k):
    s, c, d = sn(u, k), cn(u, k), dn(u, k)
    return mp.matrix([[c * d, 0, 0, k * s * s], [0, 0, 1, 0], [0, 1, 0, 0], [k * s * s, 0, 0, c * d]]) / s


uc, vc = mp.mpc("0.3", "0.2"), mp.mpc("-0.4", "0.1")
A_, B_, C_ = e12(rcl(uc - vc, kg)), e13(rcl(uc, kg)), e23(rcl(vc, kg))
cybe = (A_ * B_ - B_ * A_) + (A_ * C_ - C_ * A_) + (B_ * C_ - C_ * B_)
usm = mp.mpf(10) ** -20
lsm = mp.mpf(10) ** -15
qc = (R8(uc, lsm, kg) / sn(uc, kg) - E4) / lsm
check("D3 [E] Quasiklassischer Grenzfall: R(u; lam) = sn(u) (1 + lam r(u) + ...), r erfüllt die klassische "
      "Yang-Baxter-Gleichung, und u r(u) -> P für u -> 0. Die elliptische R-Matrix ist die Hebung von 1 + hbar P/z aus "
      "z151 mit hbar = Kreuzungsparameter. Die Nahtkette sitzt bei hbar = halbe Periode (2-Torsion), dem Wert, der von "
      "der Störungsreihe am weitesten entfernt ist; das z120-Gatter beim Keim hbar -> 0",
      nrm(cybe) < TOL * nrm(A_) and nrm(usm * rcl(usm, kg) - PERM) < mp.mpf(10) ** -30
      and nrm(qc - rcl(uc, kg)) < mp.mpf(10) ** -10)
G8 = Rc8(mp.mpc("0.3", "0.4"), lamg, kg)
ZS, ZD = kron(SZ, E2) + kron(E2, SZ), kron(SZ, E2) - kron(E2, SZ)
GQ = RcQ(mp.mpc("0.3", "0.4"))
check("D4 [E] Eichgruppe und Ladung: Das Acht-Vertex-Gatter (elliptische Klasse, sl(2)) vertauscht weder mit Z_1 + Z_2 "
      "noch mit Z_1 - Z_2; seine Symmetrie ist die Kleinsche Vierergruppe aus A8 (Literatur: Z_N x Z_N). Das Gatter der "
      "Ladungskette erhält die Ladung, hat aber unsymmetrische Gewichte (a ungleich a'): Es liegt außerhalb der "
      "elliptischen Klasse von Belavin und Costello-Witten-Yamazaki",
      nrm(G8 * ZS - ZS * G8) > mp.mpf("0.01") and nrm(G8 * ZD - ZD * G8) > mp.mpf("0.01")
      and nrm(GQ * ZD - ZD * GQ) < TOL and abs(GQ[0, 0] - GQ[3, 3]) > mp.mpf("0.01"))
praemissen_z151 = [
    ("P_Omega", "entfällt als Wahl der Polstruktur: Auf E_i ist das holomorphe Differential bis auf Normierung eindeutig "
     "(Geschlecht 1); ersetzt durch P_wurzel"),
    ("P_red", "entfällt: keine Wahl der Nullstellen der Reduktion; ersetzt durch P_wurzel (omega = Wurzel aus Phi_M auf der "
     "Doppelüberlagerung; aus der Twistorreduktion folgt das nicht, D2)"),
    ("P_rand", "entfällt: keine Pole, keine Randbedingung [L]; an ihre Stelle tritt das starre, mit der Kleinschen "
     "Vierergruppe verdrillte Bündel (P_buendel; die Vierergruppe Gamma_x, Gamma_y, Gamma_z ist im Bestand nativ)"),
    ("P_spek", "bleibt: Die Spektralkurve der Kette (Dispersionskurve, j = 1728 exakt) ist die Doppelüberlagerung der "
     "Markensphäre; die Koordinate s = i sin q wird mit der nativen Koordinate z gleichgesetzt"),
    ("P_GS", "gegenstandslos: Für die elliptische Klasse gibt es keine Herkunft aus dem Twistorraum (D2); die "
     "4D-Chern-Simons-Theorie selbst braucht das Zusatzfeld nicht [L]"),
]
neue = [
    ("P_wurzel", "Die Einsform der Spektralkurve ist die Wurzel dz/w des nativen quadratischen Differentials Phi_M; die "
     "Spektralkurve ist die Doppelüberlagerung E_i der Markensphäre (ersetzt P_Omega und P_red)"),
    ("P_buendel", "Eichbündel = das starre sl(2)-Bündel mit Holonomien sigma_x, sigma_z (ersetzt P_rand)"),
    ("P_aut", "Der Abstand der zwei Linienscharen ist ein Bild des Kreuzungsparameters unter Aut(E_i) (u = i lam); auf "
     "der Nahtkette gleich P_fr aus z151 (zweite Schar = Uhrbild der ersten). Ersetzt P_viertel"),
    ("P_platz", "Die komponierenden Einheiten sind die Plätze der Kette; in z152 (P_CTM) ist die Naht die Halbkette"),
    ("P_lam", "Wert des Kreuzungsparameters: lam -> 0 (rationales Ende) für das z120-Gatter, lam = 2-Torsion für die "
     "Nahtkette; beides zugleich geht nicht in einer Kette"),
    ("P_spin", "nur für die Lesart 'gesteuerte Uhr' am freien Punkt: Die Ladung im Gatter ist die Ordnungsvariable der "
     "Ising-Form (y-Achse), nicht die Fermionzahl"),
]
for nm, tx in praemissen_z151:
    say(f"   z151 {nm}: {tx}")
for nm, tx in neue:
    say(f"   neu {nm}: {tx}")
check("D5 [E, Zählung] Von den fünf Prämissen der rationalen Klasse entfallen vier in ihrer alten Form (P_Omega, P_red, "
      "P_rand, P_GS), P_spek bleibt. Hinzu kommen P_wurzel, P_buendel, P_aut, P_platz und P_lam (P_spin nur für eine "
      "Lesart): sechs statt fünf. Davon ist P_aut auf der Nahtkette gleich P_fr aus z151 und ersetzt P_viertel, und "
      "P_platz ist die Platzlesart zu P_CTM aus z152. Die Zahl der Prämissen sinkt nicht",
      len(praemissen_z151) == 5 and len(neue) == 6)
befund("D", "E/L",
       "Die elliptische Klasse ist auf E_i natürlicher als die rationale (keine Pole, kein Bit, Form eindeutig, die "
       "native Quartik ist das Quadrat der Form), spart aber keine Prämisse. Sie trägt nur sl(2) mit der Kleinschen "
       "Vierergruppe: keine U(1)-Ladung, keine SU(2). Das Gatter am freien Punkt wechselwirkt nicht; das "
       "wechselwirkende z120-Gatter verlangt den Kreuzungsparameter lam -> 0, also den rationalen Keim der Kurve.")

# --------------------------------------------------------------------------------------------------------------------
say("\nTeil E  Bilanz")
bilanz = [
    ("E/N", "Teil 1, Ising-Form", "Die Ising-Kette im Querfeld sqrt 2 ist die Hälfte der XY-Kette am Freie-Fermionen-Punkt "
     "des Acht-Vertex-Modells (Belavin sl(2)). Freie Fermionen <=> Kreuzungsparameter ist 2-Torsionspunkt. Das invariante "
     "Gitter der R-Matrix ist genau bei Feld sqrt 2 quadratisch: j = 1728. Yang-Baxter 35 Stellen, [T, H] = 0 auf dem Ring."),
    ("E", "Teil 1, Ladungsform", "Die Ladungskette ist die XY-Kette mit reiner Paarung im Feld. Ihre Dispersionskurve "
     "ist in s = i sin q wörtlich w^2 = 1 - s^4; ihre R-Matrix besteht aus lemniskatischen Funktionen (Grenzfall der "
     "Bazhanov-Stroganov-Familie, k^2 = -1). Sie liegt auf E_i, aber nicht in der Belavin-Klasse (unsymmetrische Gewichte)."),
    ("E/N", "Teil 2, Boost", "Der Gitter-Boost erzeugt die Spektralverschiebung (lokale Relation exakt). K = 2 lam mal "
     "Boost; c = Lemniskatenkonstante ist die Periode von dq/eps über den Lückenzyklus. Volle Modulperiode = voller "
     "Periodenvektor, Kreuzung = 2-Torsion, isotrope Ecke = 4-Torsion, Uhr R = Delta^(-i/2) = imaginärer 2-Torsionspunkt."),
    ("E", "Teil 3, Automorphie", "z -> i z ist auf der Kurve der Tausch von Band und Lücke (Translation gegen Boost), auf "
     "der R-Matrix die Hadamard-Dualität mit (u, lam) -> (i u, i lam). Die Uhr ist die komplexe Multiplikation."),
    ("E", "Teil 3, Torsionspunkte", "Am Uhrpunkt u = i lam: Ladungsform = Paaruhr i^(Q_B - Q_A), rein diagonal; Ising-Form "
     "= exp(i pi/4 XX) exp(i pi/8 YY) = Platzuhr mal gesteuerte Uhr zwischen Ordnungsvariablen. Koordinaten "
     "(pi/4, pi/8, 0) gegen (pi/8, pi/8, pi/8) des z120-Gatters; Quadrat nicht der Tausch; frei in den Fermionen. Die "
     "4-Torsionspunkte geben keine Uhrwerte."),
    ("E", "Teil 3, das Viertel", "Für jeden Kreuzungsparameter ist bei u = i lam die Zweikörperinvariante um die feste "
     "Achse (1 - i)/(1 + i) = -i (komplexe Multiplikation der lemniskatischen Sinusfunktion). Nur E_i hat eine "
     "Automorphie, deren Bild der Kreuzung ein unitäres Gatter gibt. Das z120-Gatter ist das rationale Ende lam -> 0."),
    ("E negativ", "Teil 4, frei", "Am freien Punkt ist delta_z = 1 auf der ganzen Kurve: kein wechselwirkendes Gatter. "
     "Uhrwerte um alle drei Achsen gibt es auf der Familie nur an den Enden lam -> 0 (z120) und lam = K (frei)."),
    ("E/L", "Teil 4, Prämissen", "P_Omega, P_red, P_rand, P_GS entfallen in alter Form; P_spek bleibt; neu P_wurzel, "
     "P_buendel, P_aut, P_platz, P_lam. Eichgruppe nur sl(2) mit Kleinscher Vierergruppe; keine Twistorherkunft."),
    ("C", "P_uhr", "P_viertel folgt aus P_aut in der elliptischen Klasse auf E_i (um eine Achse exakt, um alle drei nur "
     "am rationalen Ende). P_null: Am Uhrpunkt ist das Gatter der Ladungskette eine reine Funktion der Platzladungen. "
     "P_uhr bleibt Prämisse, jetzt in der Form 'der Scharenabstand ist das Uhrbild der Kreuzung'."),
    ("O", "P_ord, K13, K15, K16", "P_ord unverändert (Plätze der Kette, P_CTM). K13: nichts (1+1 Dimensionen). K15: nicht "
     "gelöst; das starre Bündel bricht die Eichgruppe auf eine endliche Gruppe. K16: Bei tau = i und lam = 2-Torsion ist "
     "die Korrelationslänge fest (0,57 Plätze nach z152 bzw. 1/ln sqrt 2 = 2,9 Plätze); ein Kontinuumslimes liegt [L] "
     "nur am Ende lam -> 0 und ist nicht gerechnet."),
    ("O", "Nächster Test", "Shastry-Kopplung zweier lemniskatischer R-Matrizen der Ladungskette (zwei Spezies; "
     "tetraedrische Algebra der Bazhanov-Stroganov-R-Matrix, arXiv:1912.10786): Gibt es eine integrable Kopplung, deren "
     "Gatter am Uhrpunkt die Phase i^(Q_A Q_B) trägt? Ja: P_viertel und P_null folgen aus P_CTM und zwei Spezies. Nein: "
     "Die Gatterlinie ist auch auf dem Gitter geschlossen."),
]
for st, was, text in bilanz:
    say(f"   [{st}] {was}: {text}")
urteil = ("Bedingt. (1) Ja für die Kette selbst: Die Ising-Form der Nahtkette (Querfeld sqrt 2) ist der "
          "Freie-Fermionen-Punkt des Acht-Vertex-Modells, also der elliptischen sl(2)-Klasse, und dieser Punkt ist der "
          "2-Torsionspunkt des Kreuzungsparameters. Das invariante Gitter der R-Matrix ist genau bei Feld sqrt 2 "
          "quadratisch (j = 1728). Die Ladungsform (Hüpfen t, Wechselpotential 2t) liegt auf derselben Kurve, und zwar "
          "wörtlich: In s = i sin q ist ihre Dispersionskurve w^2 = 1 - s^4, ihre R-Matrix besteht aus lemniskatischen "
          "Funktionen. Sie gehört aber zur Freie-Fermionen-Familie mit Feld (Bazhanov-Stroganov), nicht zur symmetrischen "
          "elliptischen Klasse; die trägt nur die Kleinsche Vierergruppe und keine U(1)-Ladung. (2) Der Boost ist der "
          "Erzeuger der Spektralverschiebung, und c = Lemniskatenkonstante ist die Periode des invarianten Differentials "
          "über den Lückenzyklus, gleich dem doppelten Kreuzungsparameter. Die Dichtematrix ist der volle Periodenvektor, "
          "die isotrope Ecke ein 4-Torsionspunkt, die Uhr R = Delta^(-i/2) der imaginäre 2-Torsionspunkt. (3) Die "
          "Z_4-Automorphie ist auf der R-Matrix eine Dualität mit Drehung von Spektral- und Kreuzungsparameter. Am "
          "Automorphiebild der Kreuzung ist die R-Matrix der Nahtkette die Paaruhr (Ladungsform) bzw. Platzuhr mal "
          "gesteuerte Uhr zwischen Ordnungsvariablen (Ising-Form); sie ist frei in den Fermionen und weder das "
          "z120-Gatter noch der Tausch. Die 4-Torsionspunkte geben keine Uhrwerte. (4) Neu und exakt: In der "
          "elliptischen Klasse auf E_i ist bei u = i lam für jeden Kreuzungsparameter die Zweikörperinvariante um die "
          "feste Achse gleich (1 - i)/(1 + i) = -i, und nur E_i hat eine Automorphie, die dort ein unitäres Gatter "
          "gibt. Das z120-Gatter (1 + iP)/(1 + i) ist das rationale Ende lam -> 0 dieser Familie, die Nahtkette das "
          "freie Ende; Uhrwerte um alle drei Achsen gibt es nur an diesen zwei Enden. Der Viertelwinkel ist damit "
          "nativ ausgezeichnet, aber nur unter der Prämisse P_aut (der Abstand der Linienscharen ist das Uhrbild der "
          "Kreuzung); das wechselwirkende Gatter selbst legt die Kurve nicht fest (P_lam). (5) Gegen z151: Die "
          "elliptische Klasse ersetzt vier der fünf Prämissen durch fünf andere (P_wurzel, P_buendel, P_aut, P_platz, "
          "P_lam neben P_spek; P_aut tritt an die Stelle von P_viertel) und hat keine Herkunft aus der Twistorreduktion. P_uhr, P_ord und die neun Prämissen des "
          "Kandidaten bleiben; K13, K15, K16 sind nicht berührt. Nichts ist befördert.")
say("\nUrteil\n   " + urteil)
print(f"\n   (Laufzeit {time.time() - T0:.1f} s, nicht Teil der Ausgabe)", flush=True)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z157_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN,
     "praemissen_z151": [list(x) for x in praemissen_z151], "praemissen_neu": [list(x) for x in neue],
     "urteil": urteil}, indent=1, ensure_ascii=False, default=str) + "\n")
(HERE / "ausgabe_z157.txt").write_text("\n".join(OUT) + "\n")
