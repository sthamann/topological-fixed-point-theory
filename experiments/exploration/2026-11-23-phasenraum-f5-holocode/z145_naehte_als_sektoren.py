"""Z145: Nähte als Sektoren eines Netzes, das Gatter als Statistikoperator? Test des Ein-Prinzip-Kandidaten für P_uhr.

Leitlinie 'Rekonstruktion statt Auswahl'; H6 und die gelösten Blöcke (Flavour, Massen, alpha, 16, E8) bleiben unberührt,
nichts wird befördert. Statuskürzel: E exakt, C bedingt (Prämisse genannt), L Lesart, N numerisch, O offen.

Frage (Schritt 1 nach H6h50): Folgt P_uhr ('die Kompositionsphase ist eine Funktion der Uhrwerte', z144) aus einem
Prinzip, nämlich daraus, dass Nähte lokalisierte Sektoren eines Netzes sind und die Komposition ihre Statistik ist?

Bestand: z113 (Ladungssektor als Spektralfluss, L_0 = (Q -+ 1/4)^2/2 - 1/32 + Oszillatoren, e^(2 pi i L_0) = R), z120
(Gatter U = B^+ e^(-i phi H_Q) B, B = e^(-i pi J_x/2)), z121 (nativ ist das Gatter trivial), z133 (U_in U_out einkörperig),
z144 (P_uhr <=> phi in (pi/2) Z und kein Oszillatoranteil; Wechselphase i^(Q_A Q_B); konstante Zopfrelation nur
einkörperig; Diskriminantenformen 5a^2/8 von D5 und 3a^2/8 von A3).

Teil A  Das native Netz: die Uhr als Eichgruppe des Nahtfermions
  A1  Viertelverdrilltes Fermion im Fenster: e^(2 pi i L_0) = i^Q, L_0 - (Q^2/2 - Q/4) ganz und >= 0, Zustandszahlen.
  A2  Uhrinvarianter Teil: Gitter 4Z in Z, Sektorgruppe Z_16 mit q(n) = n^2/32; die Naht ist der Sektor n = 4Q - 1.
  A3  Statistik je Einheit: Fermion (n = 4), Viertelpaar (n = 2), Viertel (n = 1); Ladung um Fluss.
  A4  Kriterium: Die Wechselphase hat Ordnung 4 genau für Einheiten mit Gewicht ungerade/8.
Teil B  Das z120-Gatter im Fockmodell zweier Nähte, umgeschrieben
  B1  Zirkulare Moden = Eigenmoden der Vierteldrehung P_sym; B ist der Basiswechsel; Q_+- = N/2 +- J_y.
  B2  U(phi) = e^(-i phi (Q_+^2 + Q_-^2)/2); bei pi/2: gesteuerte Uhr i^(Q_+ Q_-) mal Gesamtladungsphase;
      gleichwertig eine Paarphase: -i je Paar von Anregungen in derselben zirkularen Mode.
  B3  Heisenberg-Form von P_uhr: U X U^+ = X mal Uhrwert genau für phi in (pi/2) Z.
  B4  Die acht Uhrgatter sind die achten Wurzeln der KMS-Periode beta = 4 pi des Ladungsgitters.
Teil C  Quadratische Formen auf der Uhrgruppe
  C1  Die vier nichtentarteten Formen c a^2/8 (c = 1, 3, 5, 7) aus D1, A3, D5, D7; Zuordnung zu den Gattern.
  C2  Weil-Darstellung: S unitär genau für ungerade k; (ST)^3 = S^2.
  C3  Fünfeck und Sechseck: Assoziator (-1)^(a Übertrag(b,c)), Klasse 2 in H^3(Z_4, U(1)); Paare sind anomaliefrei.
  C4  Unter P_lift: Form 5a^2/8, Gatter = (-1)^N mal z120-Gatter in umgekehrter Zeitrichtung; D5 mal A3 = 1.
Teil D  E8 als Probe: Drehungen des D5-Faktors auf Ketten von Klebewurzeln (2 pi: halbes Gatter; 4 pi: das Gatter)
Teil E  Bilanz und Urteil.
"""
import cmath
import itertools
import json
import os
import time
from fractions import Fraction as Fr
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")                     # deterministische Rundung (einfädig)

import numpy as np
import sympy as sp
from scipy.linalg import expm

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
T0 = time.time()
TOL = 1e-10


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val, text):
    BEFUNDE[name] = val
    say(f"   [Befund: {val}] {text}")


def frac1(x):
    """Bruchteil in [0, 1)."""
    x = Fr(x)
    return x - (x.numerator // x.denominator)


def ph(x):
    return cmath.exp(2j * cmath.pi * float(x))


def dag(a):
    return a.conj().T


def nrm(a):
    return float(np.linalg.norm(a))


I2 = np.eye(2, dtype=complex)
Z2 = np.diag([1, -1]).astype(complex)


def kron(*ms):
    out = np.eye(1, dtype=complex)
    for m in ms:
        out = np.kron(out, m)
    return out


def fock(n):
    low = np.array([[0, 1], [0, 0]], complex)
    return [kron(*([Z2] * k + [low] + [I2] * (n - k - 1))) for k in range(n)]


# ------------------------------------------------------------------------------------------------------------------
say("Z145: Nähte als Sektoren eines Netzes, das Gatter als Statistikoperator?")
say("=" * 100)
say("\nTeil A: Das native Netz. Die Uhr als Eichgruppe des Nahtfermions")

# A1: Teilchen (Ladung +1) mit Energien 1/4 + j, Löcher (Ladung -1) mit Energien 3/4 + j, j = 0..M-1
M = 5
teil = [Fr(1, 4) + j for j in range(M)]
loch = [Fr(3, 4) + j for j in range(M)]
ok_uhr, ok_boson, zahl = True, True, {}
for bits in itertools.product((0, 1), repeat=2 * M):
    Q = sum(bits[:M]) - sum(bits[M:])
    E = sum(e for e, b in zip(teil + loch, bits) if b)
    ok_uhr &= frac1(E) == frac1(Fr(Q, 4))
    lev = E - (Fr(Q * Q, 2) - Fr(Q, 4))
    ok_boson &= lev.denominator == 1 and lev >= 0
    zahl[(Q, int(lev))] = zahl.get((Q, int(lev)), 0) + 1
check("[E] Viertelverdrilltes Fermion (1024 Zustände): e^(2 pi i L_0) = i^Q für jeden Zustand; die Uhr ist linear in der "
      "Ladung", ok_uhr)
check("[E] L_0 = Q^2/2 - Q/4 + Niveau mit ganzzahligem Niveau >= 0 (Spektralfluss aus z113); im Sektor Q = 0, 1, -1 "
      "sind die Zustandszahlen der Niveaus 0 bis 3 die Partitionszahlen 1, 1, 2, 3",
      ok_boson and all([zahl.get((Q, l), 0) for l in range(4)] == [1, 1, 2, 3] for Q in (0, 1, -1)))


# A2: Z_16 mit q(n) = n^2/32
def q16(n):
    return frac1(Fr(n * n, 32))


def b16(n, m):
    return frac1(Fr(n * m, 16))


ok = all(q16(n + 16) == q16(n) and frac1(q16(n + m) - q16(n) - q16(m)) == b16(n, m)
         for n in range(-20, 20) for m in range(-20, 20))
check("[E] Uhrinvarianter Teil des Fermions: Gitter 4Z (Norm 16), Dual (1/4)Z, Sektorgruppe Z_16 mit q(n) = n^2/32 und "
      "Bicharakter e^(2 pi i n m/16); ganzzahlige Ladung Q ist n = 4Q", ok)
ok = all(frac1(q16(4 * Q - 1) - Fr(1, 32)) == frac1(Fr(Q * Q, 2) - Fr(Q, 4)) for Q in range(-12, 13))
check("[E] Die Naht mit Ladung Q ist der Sektor n = 4Q - 1 (Viertelverdrillung): Sein Twist ist e^(2 pi i/32) i^Q, also "
      "bis auf eine Konstante die Uhr aus A1", ok)
ok = all(frac1(b16(4 * a - 1, 4 * b - 1) - Fr(1, 16) + Fr(a + b, 4)) == 0 for a in range(-8, 9) for b in range(-8, 9))
check("[E negativ] Monodromie zweier Nähte (4 Q_A - 1, 4 Q_B - 1) = e^(2 pi i/16) (-i)^(Q_A + Q_B): ein Produkt von "
      "Einnahtfunktionen ohne Term in Q_A Q_B. Die native Statistik der Nahtsektoren ist einkörperig", ok)

# A3: Einheiten
tab = {}
for s, name in ((4, "Fermion"), (2, "Viertelpaar"), (1, "Viertel")):
    tab[name] = sorted({b16(s * a, s * b) for a in range(16) for b in range(16)})
check("[E] Wechselphase von a und b Einheiten der Größe s (n = s a, s b): Fermion s = 4 immer 1; Viertelpaar s = 2 "
      "genau i^(a b); Viertel s = 1 sechzehnte Wurzeln",
      tab["Fermion"] == [0] and tab["Viertelpaar"] == [0, Fr(1, 4), Fr(1, 2), Fr(3, 4)]
      and all(b16(2 * a, 2 * b) == frac1(Fr(a * b, 4)) for a in range(8) for b in range(8))
      and len(tab["Viertel"]) == 16)
check("[E] Ladung um Fluss: Q_A Fermionen (n = 4 Q_A) um Q_B Viertel (n = Q_B) geben i^(Q_A Q_B) = R_B^(Q_A). Die "
      "Wechselphase des z120-Gatters ist die Uhr der einen Naht, gesteuert durch die Ladung der anderen",
      all(b16(4 * a, b) == frac1(Fr(a * b, 4)) for a in range(-8, 9) for b in range(-8, 9)))

# A4: Kriterium über das Gewicht der Einheit
gew = {"Fermion": Fr(1, 2), "Boson (Strom)": Fr(1), "Viertel": Fr(1, 32), "Viertelpaar": Fr(1, 8),
       "4 von SU(4) (A3)": Fr(3, 8), "16 von Spin(10) (D5)": Fr(5, 8), "Spinor von D7": Fr(7, 8)}
ordn = {}
for name, h in gew.items():
    w = frac1(2 * h)                                   # Wechselphase zweier Einheiten: e^(2 pi i 2h)
    ordn[name] = next(k for k in range(1, 65) if frac1(k * w) == 0)
ZAHLEN["Ordnung der Wechselphase je Einheit"] = ordn
check("[E] Kriterium: Die Wechselphase e^(2 pi i 2 h a b) hat Ordnung 4 genau für Gewicht h = ungerade/8. Fermion und "
      "Strom geben 1 (freie Komposition); Viertelpaar, die 4 von SU(4), die 16 von Spin(10) geben i^(+-ab)",
      ordn["Fermion"] == 1 and ordn["Boson (Strom)"] == 1 and ordn["Viertel"] == 16
      and all(ordn[k] == 4 for k in ("Viertelpaar", "4 von SU(4) (A3)", "16 von Spin(10) (D5)", "Spinor von D7")))
befund("native Statistik trägt das wechselwirkende Gatter", "nein",
       "Teil A: Im uhrgeeichten Netz des Nahtfermions ist die Statistik ganzzahliger Ladungen linear (Uhr, Parität). "
       "Die Phase i^(Q_A Q_B) braucht eine Einheit mit Gewicht ungerade/8 oder die Lesart Ladung um Fluss.")

# ------------------------------------------------------------------------------------------------------------------
say("\nTeil B: Das z120-Gatter im Fockmodell zweier Nähte, umgeschrieben")
res_B = {}
for K in (3, 4):
    c = fock(2 * K)
    cA, cB = c[:K], c[K:]
    NA = sum(dag(x) @ x for x in cA)
    NB = sum(dag(x) @ x for x in cB)
    nA, nB = np.rint(np.diag(NA).real).astype(int), np.rint(np.diag(NB).real).astype(int)
    Jp = sum(dag(a) @ b for a, b in zip(cA, cB))
    Jx, Jy = (Jp + dag(Jp)) / 2, (Jp - dag(Jp)) / 2j
    Bq = expm(-1j * np.pi / 2 * Jx)
    Ntot = NA + NB

    def fn(f, _B=Bq, _a=nA, _b=nB):
        """Funktion der Ladungen der zirkularen Moden: B^+ f(N_A, N_B) B."""
        return dag(_B) @ np.diag(f(_a, _b).astype(complex)) @ _B

    def U(phi, _fn=fn):
        return _fn(lambda a, b: np.exp(-1j * phi * (a * a + b * b) / 2))

    r = {}
    dplus = [dag(Bq) @ x @ Bq for x in cA]
    dminus = [dag(Bq) @ x @ Bq for x in cB]
    r["zirkular"] = max(nrm(dplus[y] - (cA[y] - 1j * cB[y]) / np.sqrt(2)) for y in range(K)) \
        + max(nrm(dminus[y] + 1j * (cA[y] + 1j * cB[y]) / np.sqrt(2)) for y in range(K))
    Qp, Qm = dag(Bq) @ NA @ Bq, dag(Bq) @ NB @ Bq
    r["Q+-"] = nrm(Qp - (Ntot / 2 + Jy)) + nrm(Qm - (Ntot / 2 - Jy))
    Gy = expm(1j * np.pi * Jy)
    r["Vierteldrehung"] = max(nrm(Gy @ cA[y] @ dag(Gy) + cB[y]) + nrm(Gy @ cB[y] @ dag(Gy) - cA[y]) for y in range(K))
    r["Eigenmoden"] = max(nrm(Gy @ dplus[y] @ dag(Gy) + 1j * dplus[y]) + nrm(Gy @ dminus[y] @ dag(Gy) - 1j * dminus[y])
                          for y in range(K))
    r["Paaruhr"] = nrm(Gy - fn(lambda a, b: 1j ** ((a - b) % 4)))
    # B2
    Uq = U(np.pi / 2)
    r["Gatter"] = nrm(Uq - dag(Bq) @ expm(-1j * np.pi / 2 * (NA @ NA + NB @ NB) / 2) @ Bq)
    CR = fn(lambda a, b: 1j ** ((a * b) % 4))
    r["gesteuerte Uhr"] = nrm(Uq - expm(-1j * np.pi / 4 * Ntot @ Ntot) @ CR)
    r["Uhr hoch Ladung"] = nrm(CR - sum(
        dag(Bq) @ np.diag((nA == a).astype(complex)) @ Bq @ np.linalg.matrix_power(fn(lambda x, y: 1j ** (y % 4)), a)
        for a in range(K + 1)))
    paare = fn(lambda a, b: (-1j) ** ((a * (a - 1) // 2 + b * (b - 1) // 2) % 4))
    r["Paarphase"] = nrm(Uq - expm(-1j * np.pi / 4 * Ntot) @ paare)
    # B3
    X = dag(dplus[0]) @ dminus[0]
    r["Heisenberg"] = nrm(Uq @ X @ dag(Uq) + 1j * X @ dag(Gy))
    cols = np.array([(X @ np.linalg.matrix_power(Gy, k)).ravel() for k in range(4)]).T

    def rest(phi, _U=U, _X=X, _cols=cols):
        v = (_U(phi) @ _X @ dag(_U(phi))).ravel()
        return float(np.linalg.norm(v - _cols @ np.linalg.lstsq(_cols, v, rcond=None)[0]))

    r["Normalisator"] = max(rest(k * np.pi / 2) for k in range(8))
    r["kein Normalisator"] = min(rest(p) for p in (0.3, 0.7, np.pi / 4, 1.0, 2.0))
    # B4
    r["Periode"] = nrm(U(4 * np.pi) - np.eye(len(NA))) + (0 if nrm(U(2 * np.pi) - np.eye(len(NA))) > 1 else 1)
    r["Uhrfunktion"] = max(nrm(U(k * np.pi / 2) - fn(lambda a, b, k=k: np.exp(-1j * np.pi * k * ((a % 4) ** 2 + (b % 4) ** 2) / 4)))
                           for k in range(8))
    # C4 (Vorgriff): Formen 5a^2/8 und 3a^2/8 auf den zirkularen Ladungen
    par = fn(lambda a, b: (-1.0) ** (a + b))
    T5 = fn(lambda a, b: np.exp(2j * np.pi * 5 * ((a % 4) ** 2 + (b % 4) ** 2) / 8))
    T3 = fn(lambda a, b: np.exp(2j * np.pi * 3 * ((a % 4) ** 2 + (b % 4) ** 2) / 8))
    r["D5"] = nrm(T5 - par @ U(-np.pi / 2)) + nrm(par - expm(1j * np.pi * Ntot))
    r["A3"] = nrm(T3 - par @ U(np.pi / 2))
    r["D5 mal A3"] = nrm(T5 @ T3 - np.eye(len(NA)))
    res_B[K] = r

mx = lambda key: max(res_B[K][key] for K in res_B)
check("[E] Zirkulare Moden d_+- = (c_A -+ i c_B)/sqrt 2 sind B^+ c_A B und (bis auf -i) B^+ c_B B; ihre Ladungen sind "
      "Q_+- = N/2 +- J_y (K = 3, 4)", mx("zirkular") < TOL and mx("Q+-") < TOL)
check("[E] Die Vierteldrehung A -> B -> -A (P_sym) ist e^(i pi J_y); die zirkularen Moden sind ihre Eigenmoden zu -+i, "
      "und sie ist die Paaruhr i^(Q_+ - Q_-) der zirkularen Moden", mx("Vierteldrehung") < TOL
      and mx("Eigenmoden") < TOL and mx("Paaruhr") < TOL)
check("[E] U(phi) = B^+ e^(-i phi H_Q) B = e^(-i phi (Q_+^2 + Q_-^2)/2); bei phi = pi/2 ist U = e^(-i pi N^2/4) i^(Q_+ Q_-): "
      "bis auf die Gesamtladungsphase die gesteuerte Uhr R_-^(Q_+) zwischen den zirkularen Moden",
      mx("Gatter") < TOL and mx("gesteuerte Uhr") < TOL and mx("Uhr hoch Ladung") < TOL)
check("[E] Paarphase: U(pi/2) = e^(-i pi N/4) (-i)^P mit P = Zahl der Paare von Anregungen in derselben zirkularen Mode "
      "(Binomial(Q_+, 2) + Binomial(Q_-, 2)). Das Gatter ist eine Hubbard-artige Paarwechselwirkung der zirkularen Moden "
      "über die Zeit pi/2", mx("Paarphase") < TOL)
check("[E] Heisenberg-Form: Für den Einheitstransfer X = d_+^+ d_- gilt U X U^+ = X e^(-i phi (2 J_y + 1)); bei pi/2 ist das "
      "-i X mal inverse Paaruhr. U X U^+ liegt genau für phi in (pi/2) Z in X mal Uhrgruppe (Kontrollen 0,3; 0,7; pi/4; 1; 2)",
      mx("Heisenberg") < TOL and mx("Normalisator") < 1e-8 and min(res_B[K]["kein Normalisator"] for K in res_B) > 0.1)
check("[E] e^(-i t H_Q) hat auf dem Ladungsgitter genau die Periode 4 pi = beta (bei 2 pi bleibt die Parität). Die acht "
      "Gatter U(k pi/2) = U(k beta/8) sind genau die Funktionen der Ladungen modulo 4",
      mx("Periode") < TOL and mx("Uhrfunktion") < 1e-8)
ZAHLEN["Rest von U X U^+ außerhalb der Uhrgruppe (Minimum der Kontrollen)"] = {
    K: round(res_B[K]["kein Normalisator"], 6) for K in res_B}
befund("P_uhr umformuliert", "ja",
       "Teil B: P_uhr sagt gleichwertig (i) die Phase ist eine Funktion der Uhrwerte, (ii) das Gatter bildet einen "
       "Einheitstransfer auf Transfer mal Uhrwert ab (Clifford-Bedingung), (iii) die Laufzeit des Ladungsgitters ist "
       "ein Vielfaches von beta/8.")

# ------------------------------------------------------------------------------------------------------------------
say("\nTeil C: Quadratische Formen auf der Uhrgruppe Z_4")


def cartan_A(n):
    return sp.Matrix(n, n, lambda i, j: 2 if i == j else (-1 if abs(i - j) == 1 else 0))


def cartan_D(n):
    m = sp.zeros(n, n)
    for i in range(n):
        m[i, i] = 2
    for i in range(n - 2):
        m[i, i + 1] = m[i + 1, i] = -1
    m[n - 3, n - 1] = m[n - 1, n - 3] = -1
    return m


formen = {"D1": Fr(1, 8)}                                # D1 = 2Z, Dual (1/2)Z, Erzeuger 1/2 mit Norm 1/4
for name, C, node in (("A3", cartan_A(3), 0), ("D5", cartan_D(5), 4), ("D7", cartan_D(7), 6)):
    inv = C.inv()
    assert C.det() == 4
    formen[name] = frac1(Fr(int(inv[node, node].p), int(inv[node, node].q)) / 2)
check("[E] Die vier nichtentarteten quadratischen Formen c a^2/8 auf Z_4 sind die Diskriminantenformen von D1 (c = 1), "
      "A3 (c = 3), D5 (c = 5) und D7 (c = 7), aus den Cartan-Matrizen",
      formen == {"D1": Fr(1, 8), "A3": Fr(3, 8), "D5": Fr(5, 8), "D7": Fr(7, 8)})


def gatterform(k):
    """Exponent/(2 pi) der Phase von U(k pi/2) je zirkularer Mode: -k a^2/8."""
    return [frac1(Fr(-k * a * a, 8)) for a in range(4)]


zu = {k: next((n for n, c in formen.items() if gatterform(k) == [frac1(c * a * a) for a in range(4)]), None)
      for k in range(8)}
ZAHLEN["Gatter k -> Form"] = {k: (v or "entartet") for k, v in zu.items()}
check("[E] Zuordnung: U(k pi/2) trägt für k = 1, 3, 5, 7 die Form von D7, D5, A3, D1; für gerades k ist die Form "
      "entartet (k = 0 trivial, k = 4 Parität, k = 2, 6 Rang 2)",
      zu == {0: None, 1: "D7", 2: None, 3: "D5", 4: None, 5: "A3", 6: None, 7: "D1"})

ok_weil, raenge = True, {}
for k in range(8):
    th = np.array([ph(-Fr(k * a * a, 8)) for a in range(4)])
    Mm = np.array([[th[(a + b) % 4] / (th[a] * th[b]) for b in range(4)] for a in range(4)])
    Sm = Mm.conj() / 2                                 # S_ab = e^(-2 pi i b(a,b))/2
    raenge[k] = int(np.linalg.matrix_rank(Mm, tol=1e-9))
    if k % 2:
        gauss = th.sum() / 2
        Tm = np.diag(th) / gauss ** (1 / 3)
        Cc = np.array([[1.0 if (a + b) % 4 == 0 else 0.0 for b in range(4)] for a in range(4)])
        ok_weil &= nrm(Sm @ dag(Sm) - np.eye(4)) < TOL and abs(abs(gauss) - 1) < TOL
        ok_weil &= nrm(np.linalg.matrix_power(Sm @ Tm, 3) - Sm @ Sm) < 1e-9 and nrm(Sm @ Sm - Cc) < TOL
        ok_weil &= abs(gauss - ph(Fr((-k) % 8, 8))) < TOL
ZAHLEN["Rang der Monodromiematrix von U(k pi/2)"] = raenge
check("[E] Modularität: Die Monodromiematrix hat Rang 1, 4, 2, 4, 1, 4, 2, 4 für k = 0..7. Für ungerades k ist S unitär, "
      "S^2 die Ladungskonjugation, (ST)^3 = S^2 und die Gauß-Summe e^(2 pi i c/8) mit c = -k mod 8",
      ok_weil and [raenge[k] for k in range(8)] == [1, 4, 2, 4, 1, 4, 2, 4])


def carry(b, c):
    return (b + c) // 4


ok_f, inv_klasse = True, {}
for cc in (1, 3, 5, 7):
    Rm = lambda a, b, cc=cc: ph(Fr(cc * a * b, 8))
    Fm = lambda a, b, c: (-1.0) ** (a * carry(b, c))
    for a, b, c in itertools.product(range(4), repeat=3):
        bc, ab = (b + c) % 4, (a + b) % 4
        ok_f &= abs(Rm(a, bc) * Fm(b, c, a) * Fm(a, b, c) - Rm(a, b) * Rm(a, c) * Fm(b, a, c)) < TOL
        ok_f &= abs(Rm(ab, c) * Fm(c, a, b) * Fm(a, b, c) - Rm(a, c) * Rm(b, c) * Fm(a, c, b)) < TOL
        for d in range(4):
            ok_f &= abs(Fm(a, b, c) * Fm(a, bc, d) * Fm(b, c, d) - Fm(ab, c, d) * Fm(a, b, (c + d) % 4)) < TOL
    triv = all(abs(Rm(a, (b + c) % 4) - Rm(a, b) * Rm(a, c)) < TOL for a, b, c in itertools.product(range(4), repeat=3))
    inv_klasse[cc] = (int(round(float(np.prod([Fm(1, j, 1) for j in range(4)])))), triv)
check("[E] Fünfeck und beide Sechsecke gelten für R(a,b) = e^(2 pi i c a b/8) mit dem Assoziator F = (-1)^(a Übertrag(b,c)); "
      "mit trivialem Assoziator scheitert das Sechseck. Die Invariante prod_j F(1, j, 1) = -1: Klasse 2 in "
      "H^3(Z_4, U(1)) = Z_4 für alle vier Formen", ok_f and all(v == (-1, False) for v in inv_klasse.values()))
check("[E] Zwei Faktoren zusammen (zwei Nähte, oder D5 mit A3) haben auf der diagonalen Uhr die Klasse 2 + 2 = 0: Die "
      "gemeinsame Uhr ist anomaliefrei, die einer einzelnen Naht nicht", (2 + 2) % 4 == 0)
check("[E] Unter P_lift (Uhr = Zentrum von Spin(10)) ist die Form 5 a^2/8. Die Funktion theta_D5(Q_+) theta_D5(Q_-) ist "
      "(-1)^N U(-pi/2), theta_A3(Q_+) theta_A3(Q_-) ist (-1)^N U(+pi/2), und ihr Produkt ist 1 (K = 3, 4). Das z120-Gatter "
      "ist bis auf die zentrale Parität der Twist des A3-Faktors und der inverse Twist des D5-Faktors",
      mx("D5") < 1e-8 and mx("A3") < 1e-8 and mx("D5 mal A3") < 1e-8)
befund("Hebung legt die Phase fest", "bedingt",
       "Teil C: Die acht Uhrgatter sind eine Z_8; die vier wechselwirkenden entsprechen D1, A3, D5, D7. P_lift wählt "
       "D5. Die zwei Zeitrichtungen des Gatters sind die zwei Faktoren D5 und A3 der E8-Verklebung.")

# ------------------------------------------------------------------------------------------------------------------
say("\nTeil D: E8 als Probe. Die 2 pi-Drehung des D5-Faktors auf Ketten von Klebewurzeln")
H = Fr(1, 2)
wurz = []
for i, j in itertools.combinations(range(8), 2):
    for si, sj in itertools.product((1, -1), repeat=2):
        v = [Fr(0)] * 8
        v[i], v[j] = Fr(si), Fr(sj)
        wurz.append(tuple(v))
for sg in itertools.product((1, -1), repeat=8):
    if sum(1 for s in sg if s < 0) % 2 == 0:
        wurz.append(tuple(H * s for s in sg))
n5 = lambda v: sum(x * x for x in v[:5])
n3 = lambda v: sum(x * x for x in v[5:])
typ = {}
for v in wurz:
    typ[(n5(v), n3(v))] = typ.get((n5(v), n3(v)), 0) + 1
kleb = [v for v in wurz if n5(v) == Fr(5, 4)]
check("[E] E8 = D8 mit Spinorklasse, D5 auf den ersten fünf, A3 = D3 auf den letzten drei Koordinaten: 240 Wurzeln = 40 "
      "(D5) + 12 (A3) + 60 ((10,6)) + 128 Klebewurzeln (16,4) + (16-quer,4-quer) mit Norm 5/4 + 3/4",
      len(wurz) == 240 and typ == {(Fr(2), Fr(0)): 40, (Fr(0), Fr(2)): 12, (Fr(1), Fr(1)): 60, (Fr(5, 4), Fr(3, 4)): 128})

rng = np.random.default_rng(145)
punkte = []
for _ in range(60):
    ks, ws = rng.integers(-2, 3, 6), [wurz[t] for t in rng.integers(0, 240, 6)]
    punkte.append(tuple(sum(int(k) * w[i] for k, w in zip(ks, ws)) for i in range(8)))
ok_quad, ok_lin, ok_ganz, ok_rel, werte, werte_rel = True, True, True, True, set(), set()
for p in punkte:
    ok_ganz &= frac1(sum(x * x for x in p) / 2) == 0
    for g in kleb[::8]:
        f = lambda m: frac1(-sum((p[i] + m * g[i]) ** 2 for i in range(5)) / 2)           # 2 pi-Drehung von D5
        f3 = lambda m: frac1(-sum((p[i] + m * g[i]) ** 2 for i in range(5, 8)) / 2)       # 2 pi-Drehung von A3
        ok_quad &= all(frac1(f(m + 1) + f(m - 1) - 2 * f(m)) == Fr(3, 4) for m in range(-3, 4))
        lin = frac1(sum(p[i] * g[i] for i in range(5)))
        ok_lin &= frac1(4 * lin) == 0
        ok_ganz &= all(frac1(f(m) + f3(m)) == 0 for m in range(-3, 4))
        rel = lambda m: frac1(f(m) - f3(m))                                               # D5 um 2 pi gegen A3 um -2 pi
        ok_rel &= all(rel(m) == frac1(2 * f(m)) and frac1(rel(m + 1) + rel(m - 1) - 2 * rel(m)) == Fr(1, 2)
                      and rel(m + 2) == rel(m) for m in range(-3, 4))
        werte.add(lin)
        werte_rel.add(frac1(rel(1) - rel(0)))
check("[E] Auf jeder Kette p + m gamma (60 Gitterpunkte, 16 Klebewurzeln) ist die 2 pi-Drehung des D5-Faktors allein "
      "e^(-2 pi i L_0^D5) eine quadratische Phase in m mit zweiter Differenz -i (achte Wurzeln, Periode 4 in m); der "
      "lineare Anteil ist ein Uhrwert. Das ist die Einmodenfunktion, entlang des Transfers also das halbe Gatter",
      ok_quad and ok_lin and werte == {Fr(0), Fr(1, 4), Fr(1, 2), Fr(3, 4)})
check("[E] E8 ist gerade: e^(-2 pi i L_0^D5) e^(-2 pi i L_0^A3) = 1 auf dem ganzen Gitter. Die Drehungen von D5 und A3 "
      "sind zueinander invers, wie U_in und U_out in z133", ok_ganz)
check("[E] Die relative Drehung, D5 um 2 pi gegen A3 um -2 pi, ist die 4 pi-Drehung von D5. Auf jeder Kette hat sie die "
      "zweite Differenz -1 und Periode 2 in m, wie das z120-Gatter entlang des Transfers (e^(-i pi m^2/2), also a + b mal "
      "Tausch); der Schritt ist +-i",
      ok_rel and werte_rel <= {Fr(1, 4), Fr(3, 4)} and len(werte_rel) == 2)
ok_vier = all(frac1(2 * Fr(5 * a * a, 8)) == frac1(Fr(a * a, 4))
              and frac1(Fr(a * a, 4)) == frac1(Fr(a, 4) + Fr(a * (a - 1) // 2 % 2, 2)) for a in range(8))
ok_ab = all(frac1(Fr(5 * a * a, 8) + Fr(5 * b * b, 8)) == frac1(Fr(a * a, 4)) for a in range(8) for b in range(8)
            if (a + b) % 4 == 0)
check("[E] Auf den Klassen ist die 4 pi-Drehung von D5 gleich i^(a^2) = i^a (-1)^(a(a-1)/2): Uhr mal Vorzeichen der "
      "Umkehrung von a Fermionen. Nativ ist e^(4 pi i L_0) = (-1)^F; die gehobene 4 pi-Drehung ist eine quadratische "
      "Wurzel davon. Im niveauabgeglichenen Sektor N = 0 mod 4 ist theta_D5(Q_+) theta_D5(Q_-) = i^(Q_+^2)",
      ok_vier and ok_ab)
klass = [frac1(Fr(5 * a * a, 8)) for a in range(4)]
check("[E] Die Drehung ist keine Symmetrie von E8: Die Phasen 1, e^(2 pi i 5/8), -1, e^(2 pi i 5/8) der vier Klassen sind "
      "kein Charakter der Uhrgruppe (Klasse 1 + 1 = 2 verlangte e^(2 pi i 10/8) = i statt -1). Sie wirkt also echt auf den "
      "Klebeströmen", frac1(2 * klass[1]) != klass[2])
befund("Phase ohne Winkelprämisse im E8-Bild", "ja (Lesart)",
       "Teil D: Liest man die zwei zirkularen Moden als die Faktoren D5 und A3 und den Transfer als Klebestrom (16,4), "
       "dann ist das z120-Gatter die 4 pi-Drehung des D5-Faktors (D5 um 2 pi gegen A3): kein freier Winkel, P_null von "
       "selbst, Niveauabgleich N = 0 mod 4 eingebaut. Die Lesart ist nicht hergeleitet.")

# ------------------------------------------------------------------------------------------------------------------
say("\nTeil E: Bilanz und Urteil")
bilanz = [
    ("E negativ", "Prinzip 'Nähte sind Sektoren, das Gatter ist ihre Statistik' im nativen Netz",
     "Die Statistik der Nahtsektoren (n = 4Q - 1 in Z_16) ist ein Produkt von Uhren. Fermionische Einheiten komponieren "
     "frei; das bestätigt z121 (Gatter nativ trivial) und z144 (Zopfrelation nur einkörperig) von der Netzseite."),
    ("E", "Was das wechselwirkende Gatter verlangt",
     "Eine Einheit mit Gewicht ungerade/8 (Z_4-Anyon) oder die Lesart Ladung um Fluss: i^(Q_A Q_B) = R_B^(Q_A). Im "
     "Fockmodell ist das Gatter genau die gesteuerte Uhr zwischen den Eigenmoden der Vierteldrehung P_sym, "
     "gleichwertig die Paarphase -i je Paar in derselben zirkularen Mode."),
    ("E", "P_uhr in drei gleichwertigen Formen",
     "Funktion der Uhrwerte; Clifford-Bedingung U X U^+ = X mal Uhrwert; Laufzeit in (beta/8) Z mit beta = 4 pi als "
     "Periode des Ladungsgitters."),
    ("E/C", "Hebung",
     "Die vier wechselwirkenden Uhrgatter sind die Formen von D1, A3, D5, D7. Unter P_lift ist es D5; das z120-Gatter "
     "ist bis auf (-1)^N der Twist von A3 und der inverse Twist von D5. Zeitrichtung = Wahl des Faktors."),
    ("E", "Anomalie",
     "Jede der vier Formen hat den Assoziator der Klasse 2 in H^3(Z_4, U(1)); Paare sind anomaliefrei."),
    ("L", "E8-Bild",
     "Die 4 pi-Drehung des D5-Faktors (D5 um 2 pi gegen A3) wirkt auf Ketten von Klebewurzeln wie das z120-Gatter "
     "entlang des Transfers, ohne Winkel und ohne Oszillatoranteil; auf den Klassen ist sie Uhr mal Umkehrvorzeichen. "
     "Als Mechanismus für P_uhr ist das eine Lesart (Faktoren D5 und A3 statt zweier gleicher Nähte)."),
    ("E negativ", "Dynamik",
     "Eine punktierte Kategorie hat nur eindimensionale Fusionsräume und vertauschende Statistikphasen. Der Transfer "
     "B und jede Kettendynamik liegen außerhalb des Prinzips; K15 wird nicht berührt."),
]
for st, was, text in bilanz:
    say(f"   [{st}] {was}: {text}")
urteil = ("Der Ein-Prinzip-Kandidat ist in seiner strengen Form widerlegt und in einer schwächeren Form geschärft. (1) "
          "Streng: Im nativen Netz (Nahtfermion mit geeichter Uhr, Sektorgruppe Z_16, Naht = Sektor n = 4Q - 1) ist die "
          "Statistik der Nähte ein Produkt von Uhren. Aus 'Nähte sind Sektoren' folgt also das triviale Gatter, nicht "
          "das wechselwirkende. Fermionische Transfereinheiten komponieren frei. (2) Was das wechselwirkende Gatter "
          "verlangt, ist jetzt genau benannt: eine transferierte Einheit mit Gewicht ungerade/8, oder gleichwertig die "
          "gesteuerte Uhr, bei der die Ladung der einen zirkularen Mode die Uhr der anderen dreht. Die zirkularen "
          "Moden sind die Eigenmoden der Vierteldrehung P_sym. (3) P_uhr hat drei gleichwertige Formen: Funktion der "
          "Uhrwerte, Clifford-Bedingung (ein Transfer wird zu Transfer mal Uhrwert), Laufzeit des Ladungsgitters in "
          "Vielfachen von beta/8, wobei beta = 4 pi genau die Periode des Ladungsgitters ist. (4) Die acht Uhrgatter "
          "bilden eine Z_8; die vier wechselwirkenden tragen die Diskriminantenformen von D1, A3, D5 und D7. Unter "
          "P_lift ist die Form die von D5, und das z120-Gatter ist bis auf die zentrale Parität der Twist des "
          "A3-Faktors und der inverse Twist des D5-Faktors. Die zwei Zeitrichtungen des Gatters entsprechen den zwei "
          "Faktoren der E8-Verklebung. (5) Im E8-Gitter ist die 2 pi-Drehung des D5-Faktors allein entlang einer Kette "
          "von Klebewurzeln (16,4) nur das halbe Gatter (achte Wurzeln). Das z120-Gatter entspricht der relativen "
          "Drehung, D5 um 2 pi gegen A3 um -2 pi, also der 4 pi-Drehung des D5-Faktors: zweite Differenz -1 entlang "
          "des Transfers, kein freier Winkel, kein Oszillatoranteil, Niveauabgleich N = 0 mod 4 von selbst. Auf den "
          "Klassen ist sie i^(a^2), die Uhr mal das Vorzeichen der Umkehrung von a Fermionen; nativ ist die 4 pi-Drehung "
          "dagegen nur (-1)^F. Das ist ein Mechanismus für P_uhr nur in der Lesart, dass die zwei zirkularen Moden die "
          "Faktoren D5 und A3 sind und der Transfer ein Klebestrom ist; hergeleitet ist diese Lesart nicht. Im "
          "Fockmodell ist dasselbe Gatter eine Paarphase: -i je Paar von Anregungen in derselben zirkularen Mode. (6) Dynamik folgt aus dem Prinzip in "
          "keiner Form: Statistikphasen vertauschen, der Transfer liegt außerhalb. Für den Kandidaten: P_uhr bleibt "
          "Prämisse. Sie ist entweder durch P_lift plus die Lesart 'Transfereinheit = Klebeklasse' ersetzbar (keine "
          "Ersparnis) oder als Clifford-Abschluss zu lesen. P_ord, K15, P_curv und G bleiben offen. Nichts ist "
          "befördert.")
say("\nUrteil\n   " + urteil)
print(f"\n   (Laufzeit {time.time() - T0:.1f} s, nicht Teil der Ausgabe)", flush=True)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z145_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
(HERE / "ausgabe_z145.txt").write_text("\n".join(OUT) + "\n")
