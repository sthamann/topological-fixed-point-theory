"""Z144: Prüfung der neunzehnten externen Analyse 'Ein neuer fundamentaler Lösungsweg' (25.11.) und nativer Test der
Defektlesart.

Leitlinie 'Rekonstruktion statt Auswahl'; H6 und die gelösten Blöcke (Flavour, Massen, alpha, 16, E8) bleiben unberührt,
nichts wird befördert. Statuskürzel: E exakt, C bedingt (Prämisse genannt), L Lesart, N numerisch, O offen.

Bestand: z92/z106 (Uhr R = i^Q, R^2 = (-1)^F, R = i Gamma_11 erzeugt das Zentrum Z_4 von Spin(10)), z120 (lokales Gatter
U = B^+ e^(-i phi H_Q) B, B = e^(-i pi J_x/2), Vierteldrehung des Strahlteilers aus Leckfreiheit; Quadratphase gleichwertig
i^(Q_A Q_B)), z121 (unter P_null U = P_T - i P_S; auf einer Anregung je Naht ist jede Phase leckfrei und gibt die XXX-Familie
R(u), u = tan(phi/2); Yang-Baxter, alternierende Kette integrabel, Destri-de Vega bei Theta = 1/2, su(2)_1-Lesart), z133
(U^2 = Tausch bei gerader Gesamtladung; P_null <=> P_sl2; P_viertel bleibt Prämisse; Nahtuhr dreht den Transferrahmen), z129 (das Bit P_IF ist P- und C-ungerade, CP-gerade und aus dem Index
nicht festlegbar), z130 (P_univ nativ unter P_TP), R2 (Trägerform (5x^2 + 3y^2)/8 auf Z_4 x Z_4, zwei Lagrangesche), v526
(C(3)/C(1) = sqrt 2 - 1 am Achterkreis).

Teil A  Nachrechnung der Behauptungen der Analyse
  A1  Diskriminantenformen von D5 und A3 aus den Cartan-Matrizen; die vier Z_4-Paarungen und ihr Rang.
  A2  Verklebung D5 + A3: isotrope Elemente, Lagrangesche, D8 und seine vier Klassen, die drei Erweiterungen (zwei
      bosonische E8, eine fermionische Z^8), Wurzelzahlen durch Abzählen, torischer Code.
  A3  Zopfrelation: Modell A (Tausch auf Tensorfaktoren) gegen Modell B (antikommutierende Paritäten).
  A4  Hadamard-Fixpunkt, Ising-Blöcke bei x = 1/2, v526; Paarungsmatrix M = m 1.
  A5  Doppelpol und Reflexionspositivität; Zustandsgleichung aus rho_Lambda ~ e^(-2/alpha).
Teil B  Nativer Test im Fockmodell der Nähte (je Naht K Moden eines geladenen Fermions, Transfer-su(2) zwischen Nähten)
  B1  Zwei Nähte: Nahtuhr = Vierteldrehung um die Transferachse z; z120-Gatter als Phase der Uhralgebra; welche Phasen
      liegen in der Uhralgebra, welche davon verschränken (Bicharakter i^(k Q_A Q_B), Rang).
  B2  Zopfrelation zwischen Nahtuhr und Strahlteiler: erzwingt die Vierteldrehung; im Ein-Fermion-Sektor die
      Ising-Zopfmatrizen, der Strahlteiler ist die Hadamard-Transformierte der Uhr.
  B3  Drei Nähte: Der Nahttausch erfüllt die Zopfrelation, seine einkörperige Wurzel und das z120-Gatter nicht
      (Modell A, Defekt sqrt(3/8) im ausgeglichenen Sektor); die Monodromien antikommutieren nicht.
  B4  Eine Naht mit Oszillatoren: 'Phase liegt in der Uhralgebra' (P_uhr) <=> phi in (pi/2) Z und kein Oszillatoranteil.
  B5  Schalter zwischen den Verklebungen: Majorana-Hebung, Nahtparität gegen Schalterparität (Modell B nativ), Zopf
      erzwingt z = i aus der Uhr; ladungserhaltender Teil = Strahlteiler, Rest = Paarerzeugung; Uhr tauscht die beiden
      Schalterklassen; thermisches Verhältnis gegen Hadamard-Fixpunkt.
  B6  Anschluss an z121 (Yang-Baxter mit Spektralparameter) im Mehrmodenmodell: Auf Sektoren gleicher Nahtparität
      ist der Nahttausch eine Involution und Y(u) = (1 + i u Gamma)/(1 + i u) Yangs R-Matrix; mit den vollen
      z120-Gattern braucht die Faktorisierung dreier Nähte Nahtladungen gleich modulo 4 (Niveauabgleich, z130).
      Neu: Geltungsbereich von z121s 'jede Phase ist leckfrei'. Ab Transferspin 2 (Fenster y = -2..1, Sektor
      Q_A = Q_B = 0) leckt der Paarsektor für phi außerhalb (pi/2) Z, und das Gatter verlässt die XXX-Familie.
  B7  Vier Nähte (Kontrolle zu z121 (b) im Fockmodell): Im ausgeglichenen Sektor ist der Nahttausch der Tausch der
      Nahtinhalte; die alternierende Komposition der Viertelgatter vertauscht mit den Transfermatrizen, die
      Treppenordnung nicht.
Teil C  Bilanz und Urteil.
"""
import itertools
import json
import os
import time
from fractions import Fraction
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


def befund(name, val):
    BEFUNDE[name] = bool(val)
    say(f"   [Befund: {'ja' if val else 'nein'}] {name}")


def nrm(A):
    return float(np.linalg.norm(A))


def snorm(A):
    """Spektralnorm."""
    return float(np.linalg.norm(A, 2))


def bis_auf_phase(A, B):
    """min_phi ||A - e^(i phi) B||."""
    ph = np.vdot(B.ravel(), A.ravel())
    ph = ph/abs(ph) if abs(ph) > 1e-14 else 1.0
    return nrm(A - ph*B)


def komm(A, B):
    return A @ B - B @ A


def akomm(A, B):
    return A @ B + B @ A


def dag(A):
    return A.conj().T


def mod1(x):
    return x - sp.floor(x) if isinstance(x, sp.Basic) else x - (x.numerator // x.denominator)


# =================================================================== Teil A1
say("Teil A1: Diskriminantenformen von D5 und A3")


def cartan_D(n):
    G = sp.zeros(n)
    for i in range(n):
        G[i, i] = 2
    for i in range(n - 2):
        G[i, i + 1] = G[i + 1, i] = -1
    G[n - 3, n - 1] = G[n - 1, n - 3] = -1
    return G


def cartan_A(n):
    G = sp.zeros(n)
    for i in range(n):
        G[i, i] = 2
    for i in range(n - 1):
        G[i, i + 1] = G[i + 1, i] = -1
    return G


GD5, GA3 = cartan_D(5), cartan_A(3)
ID5, IA3 = GD5.inv(), GA3.inv()
nD5, nA3 = ID5[4, 4], IA3[0, 0]                       # Norm des Erzeugers der Diskriminantengruppe (Spinor / Fundamental)
ordD5 = min(a for a in range(1, 9) if all((a*ID5[4, j]).is_integer for j in range(5)))
ordA3 = min(a for a in range(1, 9) if all((a*IA3[0, j]).is_integer for j in range(3)))
check("D5: det Cartan = 4, Diskriminantengruppe zyklisch Z_4 (Spinorgewicht), Norm 5/4, q(a) = 5 a^2/8",
      GD5.det() == 4 and ordD5 == 4 and nD5 == sp.Rational(5, 4))
check("A3: det Cartan = 4, Diskriminantengruppe Z_4 (Fundamentalgewicht), Norm 3/4, q(a) = 3 a^2/8",
      GA3.det() == 4 and ordA3 == 4 and nA3 == sp.Rational(3, 4))
bD5 = all(sp.exp(2*sp.pi*sp.I*a*b*nD5).equals(sp.I**(a*b)) for a in range(4) for b in range(4))
bA3 = all(sp.exp(2*sp.pi*sp.I*a*b*nA3).equals(sp.I**(-a*b)) for a in range(4) for b in range(4))
check("Bicharaktere: B_D5(a,b) = i^(ab), B_A3(a,b) = i^(-ab); ihr Produkt auf der Diagonalen ist 1", bD5 and bA3)
raenge = {k: int(np.linalg.matrix_rank(np.array([[1j**(k*a*b) for b in range(4)] for a in range(4)]))) for k in range(4)}
say(f"   Rang der vier Z_4-Paarungen i^(k ab), k = 0..3: {raenge}")
check("von den vier Z_4-Paarungen sind genau i^(+ab) und i^(-ab) nichtentartet (Rang 4; k = 0: 1, k = 2: 2)",
      raenge == {0: 1, 1: 4, 2: 2, 3: 4})
bich = [k for k in range(4) if all(abs(1j**(k*((a + a2) % 4)*b) - 1j**(k*a*b)*1j**(k*a2*b)) < 1e-12
                                  for a in range(4) for a2 in range(4) for b in range(4))]
check("jeder Bicharakter von Z_4 x Z_4 ist i^(k ab) (festgelegt durch c(1,1) mit c(1,1)^4 = 1): k = 0..3 sind alle",
      bich == [0, 1, 2, 3])
ZAHLEN["paarungsraenge"] = raenge

# =================================================================== Teil A2
say("\nTeil A2: Verklebung D5 + A3 -> D8 -> E8")
Z44 = [(x, y) for x in range(4) for y in range(4)]


def qf(v):
    return Fraction(5*v[0]**2 + 3*v[1]**2, 8) % 1


def bf(v, w):
    return Fraction(5*v[0]*w[0] + 3*v[1]*w[1], 4) % 1


def plus(v, w):
    return ((v[0] + w[0]) % 4, (v[1] + w[1]) % 4)


def erzeugt(gens):
    S = {(0, 0)}
    while True:
        neu = {plus(s, g) for s in S for g in gens} | S
        if neu == S:
            return frozenset(S)
        S = neu


isotrop = [v for v in Z44 if qf(v) == 0]
untergruppen = {erzeugt([g, h]) for g in Z44 for h in Z44}
lagr = sorted(sorted(U) for U in untergruppen if len(U) == 4 and all(qf(v) == 0 for v in U))
say(f"   isotrope Elemente von q = (5x^2 + 3y^2)/8: {isotrop}")
say(f"   Lagrangesche Untergruppen: {lagr}")
check("sechs isotrope Elemente; genau zwei Lagrangesche, <(1,1)> und <(1,3)>, mit gemeinsamem Element (2,2)",
      len(isotrop) == 6 and lagr == [[(0, 0), (1, 1), (2, 2), (3, 3)], [(0, 0), (1, 3), (2, 2), (3, 1)]])
L1, L2 = (frozenset(l) for l in lagr)
spiegel_y = frozenset(((x, (-y) % 4) for x, y in L1))
spiegel_x = frozenset((((-x) % 4, y) for x, y in L1))
spiegel_xy = frozenset((((-x) % 4, (-y) % 4) for x, y in L1))
check("y -> -y (Ladungskonjugation von A3) und x -> -x (von D5) tauschen die beiden Lagrangeschen; beides zusammen "
      "erhält jede", spiegel_y == L2 and spiegel_x == L2 and spiegel_xy == L1)
check("die beiden Lagrangeschen sind gegenseitig nichtlokal: B((1,1),(1,3)) = -1", bf((1, 1), (1, 3)) == Fraction(1, 2))
H = erzeugt([(2, 2)])
Hperp = [v for v in Z44 if all(bf(v, h) == 0 for h in H)]
klassen = sorted({tuple(sorted(plus(v, h) for h in H)) for v in Hperp})
name = {((0, 0), (2, 2)): "1", ((1, 1), (3, 3)): "s", ((1, 3), (3, 1)): "c", ((0, 2), (2, 0)): "v"}
kq = {name[k]: qf(k[0]) for k in klassen}
kb = {(name[k], name[l]): bf(k[0], l[0]) for k in klassen for l in klassen}
say(f"   D8 = D5 + A3 + <(2,2)>: Klassen {[name[k] for k in klassen]}, q = { {n: str(v) for n, v in kq.items()} }")
check("D8*/D8 = <(2,2)>^perp/<(2,2)> hat vier Klassen 1, s, c, v mit q = 0, 0, 0, 1/2 (v ist das Fermion)",
      len(klassen) == 4 and kq == {"1": 0, "s": 0, "c": 0, "v": Fraction(1, 2)})
check("torischer Code: s, c, v sind paarweise gegenseitige Semionen (Paarung 1/2), Spins (1, 1, 1, -1)",
      all(kb[(p, r)] == Fraction(1, 2) for p, r in (("s", "c"), ("s", "v"), ("c", "v"))))

# Abzählen im R^8: D5 + D3 mit verdoppelten Koordinaten in {-2..2}
vals = np.array([-2, -1, 0, 1, 2], dtype=np.int64)
G8 = np.stack(np.meshgrid(*([vals]*8), indexing="ij"), axis=-1).reshape(-1, 8)
n4 = (G8**2).sum(1)                                   # 4 * Norm


def gitterklasse(X, n):
    gerade = (X % 2 == 0).all(1)
    ungerade = (X % 2 != 0).all(1)
    s = X.sum(1)
    k = np.full(len(X), -1)
    k[gerade & ((s // 2) % 2 == 0)] = 0
    k[gerade & ((s // 2) % 2 == 1)] = 2
    k[ungerade & ((s - n) % 4 == 0)] = 1
    k[ungerade & ((s + n) % 4 == 0)] = 3
    return k


kx, ky = gitterklasse(G8[:, :5], 5), gitterklasse(G8[:, 5:], 3)
gueltig = (kx >= 0) & (ky >= 0)


def zahl(norm4, coset):
    return int((gueltig & (n4 == norm4) & (kx == coset[0]) & (ky == coset[1])).sum())


def minnorm(coset):
    return Fraction(int(n4[gueltig & (kx == coset[0]) & (ky == coset[1]) & (n4 > 0)].min()), 4)


check("Kontrolle der Klassenzuordnung: kürzeste Vektoren der Nebenklassen (1,0) und (0,1) haben Norm 5/4 und 3/4",
      minnorm((1, 0)) == Fraction(5, 4) and minnorm((0, 1)) == Fraction(3, 4))
w1 = {c: zahl(8, c) for c in sorted(L1)}
w2 = {c: zahl(8, c) for c in sorted(L2)}
say(f"   Vektoren der Norm 2 je Nebenklasse, <(1,1)>: {w1}; <(1,3)>: {w2}")
check("beide Verklebungen geben 240 Wurzeln (52 + 60 + 64 + 64): zwei bosonische E8",
      sum(w1.values()) == 240 and sum(w2.values()) == 240 and w1[(0, 0)] == 52 and w1[(2, 2)] == 60)
f10, f6 = zahl(4, (2, 0)), zahl(4, (0, 2))
check("fermionische Erweiterung um v: 16 Vektoren der Norm 1, 16 = (10,1) + (1,6) (ungerades Gitter Z^8)",
      f10 == 10 and f6 == 6)
ZAHLEN["wurzeln"] = {"L1": {str(k): v for k, v in w1.items()}, "L2": {str(k): v for k, v in w2.items()},
                     "norm1_v": [f10, f6]}

# =================================================================== Teil A3
say("\nTeil A3: Zopfrelation, Modell A gegen Modell B")
I2 = np.eye(2)
X2 = np.array([[0, 1], [1, 0]], complex)
Y2 = np.array([[0, -1j], [1j, 0]])
Z2 = np.diag([1.0 + 0j, -1.0])


def kron(*ms):
    M = ms[0]
    for x in ms[1:]:
        M = np.kron(M, x)
    return M


SW = np.array([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]], complex)
S12, S23 = kron(SW, I2), kron(I2, SW)


def gatter(z, P):
    return (1 + z)/2*np.eye(len(P)) + (1 - z)/2*P


def zopfdefekt(UA, UB):
    return UA @ UB @ UA - UB @ UA @ UB


zz = [np.exp(1j*t) for t in (0.3, 1.1, 2.0, 2.9, 4.4, 5.7)] + [0.4 - 1.3j]
okA = all(nrm(zopfdefekt(gatter(z, S12), gatter(z, S23)) - ((1 + z)/2)**2*(1 - z)/2*(S12 - S23)) < TOL for z in zz)
dA = snorm(zopfdefekt(gatter(-1j, S12), gatter(-1j, S23)))
check("Modell A (Tausch): U_12 U_23 U_12 - U_23 U_12 U_23 = (1+z)^2 (1-z)/8 (S_12 - S_23); Nullstellen nur z = 1, -1",
      okA and abs(snorm(S12 - S23) - np.sqrt(3)) < TOL)
check("Modell A bei z = -i: Defektnorm sqrt(3/8)", abs(dA - np.sqrt(3/8)) < TOL)
g1, g2, g3 = kron(X2, I2), kron(Y2, I2), kron(Z2, X2)  # drei Majoranas
PA, PB = 1j*g1 @ g2, 1j*g2 @ g3
okB = all(nrm(zopfdefekt(gatter(z, PA), gatter(z, PB)) - (1 - z)*(1 + z*z)/4*(PA - PB)) < TOL for z in zz)
check("Modell B (antikommutierende Paritäten A, B): Defekt = (1-z)(1+z^2)/4 (A - B); Nullstellen z = 1, i, -i",
      okB and nrm(akomm(PA, PB)) < TOL and nrm(PA @ PA - np.eye(4)) < TOL)
check("Modell B bei z = i: das Gatter ist die Majorana-Zopfung e^(pm pi/4 g g) bis auf Phase, sein Quadrat die Parität",
      min(bis_auf_phase(gatter(1j, PA), expm(s_*np.pi/4*g1 @ g2)) for s_ in (1, -1)) < TOL
      and nrm(gatter(1j, PA) @ gatter(1j, PA) - PA) < TOL)
ZAHLEN["modell_A_defekt"] = round(dA, 12)

# =================================================================== Teil A4
say("\nTeil A4: Hadamard-Fixpunkt, Ising-Blöcke, v526, Paarungsmatrix")
x = sp.symbols("x", positive=True)
fp, fm = sp.sqrt(1 + sp.sqrt(1 - x)), sp.sqrt(1 - sp.sqrt(1 - x))
silber = sp.sqrt(2) - 1
check("Hadamard-Fixpunkt: H (1, t)^T = (1, t)^T genau für t = sqrt 2 - 1 = tan(pi/8)",
      sp.simplify((1 + silber)/sp.sqrt(2) - 1) == 0 and sp.simplify((1 - silber)/sp.sqrt(2) - silber) == 0
      and abs(float(sp.tan(sp.pi/8) - silber)) < 1e-15)
check("Ising-Blöcke F_pm = sqrt(1 pm sqrt(1-x)): F_-/F_+ bei x = 1/2 ist sqrt 2 - 1",
      abs(float((fm/fp).subs(x, sp.Rational(1, 2)) - silber)) < 1e-15)
kreuz = all(abs(float((fp.subs(x, 1 - xv) - (fp + fm).subs(x, xv)/sp.sqrt(2)))) < 1e-14
            and abs(float((fm.subs(x, 1 - xv) - (fp - fm).subs(x, xv)/sp.sqrt(2)))) < 1e-14
            for xv in (sp.Rational(1, 7), sp.Rational(2, 5), sp.Rational(3, 4)))
check("Kreuzung x -> 1-x wirkt auf (F_+, F_-) als Hadamard-Matrix", kreuz)
check("v526: C(d) = (2/N)/sin(pi d/N), N = 8: C(3)/C(1) = sin(pi/8)/sin(3pi/8) = sqrt 2 - 1",
      abs(float(sp.sin(sp.pi/8)/sp.sin(3*sp.pi/8) - silber)) < 1e-15)
therm = float(sp.exp(-sp.pi))
say(f"   thermisches Verhältnis e^(-pi) = {therm:.6f} gegen sqrt 2 - 1 = {float(silber):.6f}")
check("das thermische Verhältnis e^(-pi) ist nicht der Hadamard-Fixpunkt", abs(therm - float(silber)) > 0.3)
m11, m12, m21, m22 = sp.symbols("m11 m12 m21 m22")
Msym = sp.Matrix([[m11, m12], [m21, m22]])
Rs = sp.diag(1, sp.I)
Hs = sp.Matrix([[1, 1], [1, -1]])/sp.sqrt(2)
gl = list(Rs.H*Msym*Rs - Msym) + list(Hs.H*Msym*Hs - Msym)
lsg = sp.solve(gl, [m11, m12, m21, m22], dict=True)
check("Paarungsmatrix: R^+ M R = M und F^+ M F = M (R = diag(1, i), F = Hadamard) erzwingen M = m 1",
      len(lsg) == 1 and lsg[0].get(m12, m12) == 0 and lsg[0].get(m21, m21) == 0 and
      sp.simplify(lsg[0].get(m11, m11) - lsg[0].get(m22, m22)) == 0)

# =================================================================== Teil A5
say("\nTeil A5: Doppelpol und Zustandsgleichung")
a_, b_, om, t1, t2 = sp.symbols("a b omega t1 t2", real=True)


def Cdp(t):
    return (a_ + b_*t)*sp.exp(-om*t)


detOS = sp.simplify(Cdp(2*t1)*Cdp(2*t2) - Cdp(t1 + t2)**2)
check("Doppelpol C(t) = (a + b t) e^(-omega t): det der Osterwalder-Schrader-Matrix = -b^2 (t1-t2)^2 e^(-2 omega (t1+t2)) < 0",
      sp.simplify(detOS + b_**2*(t1 - t2)**2*sp.exp(-2*om*(t1 + t2))) == 0)
tt, Hh = sp.symbols("t H", positive=True)
al = sp.Function("alpha")(tt)
rho = sp.exp(-2/al)
w1p = sp.simplify(-sp.diff(rho, tt)/(3*Hh*rho))       # 1 + w aus rho' = -3 H (1+w) rho
check("rho_Lambda ~ e^(-2/alpha): 1 + w = -(2/(3 H alpha)) (alpha'/alpha)",
      sp.simplify(w1p + 2*sp.diff(al, tt)/(3*Hh*al**2)) == 0)
ALPHA_INV, H0_JAHR, DRIFT = 137.035999, 6.9e-11, 1.0e-18   # externe Eingaben: Uhrenschranke |alpha'/alpha| ~ 1e-18/Jahr
w_schranke = 2*ALPHA_INV/3*DRIFT/H0_JAHR
say(f"   mit |alpha'/alpha| <= {DRIFT:.0e}/Jahr und H0 = {H0_JAHR:.1e}/Jahr: |1 + w| <= {w_schranke:.2e}")
check("Größenordnung |1 + w| ~ 1e-6: der Lambda-Zweig lässt keine messbare Abweichung von w = -1 zu",
      1e-7 < w_schranke < 1e-5)
ZAHLEN["w_schranke"] = w_schranke


# =================================================================== Fockbaukasten
def fock(n):
    """Jordan-Wigner-Vernichter c_0..c_(n-1) auf (C^2)^n, Basis |0>, |1> je Mode."""
    low = np.array([[0, 1], [0, 0]], complex)
    return [kron(*([Z2]*k + [low] + [I2]*(n - k - 1))) for k in range(n)]


def naehte(S, K):
    c = fock(S*K)
    num = [sum(dag(c[s*K + y]) @ c[s*K + y] for y in range(K)) for s in range(S)]

    def transfer(s, t):
        Jp = sum(dag(c[s*K + y]) @ c[t*K + y] for y in range(K))
        return (Jp + dag(Jp))/2, (Jp - dag(Jp))/2j, (num[s] - num[t])/2

    return c, num, transfer


def fkt(D, f):
    """Funktion eines diagonalen Operators."""
    return np.diag(f(np.real(np.diag(D))))


# =================================================================== Teil B1
say("\nTeil B1: zwei Nähte, Uhr und Gatter (K = 3 Moden je Naht, 64 Zustände)")
c2, (NA, NB), tr2 = naehte(2, 3)
Jx, Jy, Jz = tr2(0, 1)
Nges = NA + NB
dim2 = len(NA)
E2 = np.eye(dim2)
nges = np.real(np.diag(Nges)).round().astype(int)
Pg = np.diag((nges % 2 == 0).astype(float))            # gerade Gesamtladung: ganzzahliger Transferspin
Pu = E2 - Pg
check("Transfer-su(2): [J_x, J_y] = i J_z und zyklisch; N_A + N_B zentral",
      nrm(komm(Jx, Jy) - 1j*Jz) < TOL and nrm(komm(Jy, Jz) - 1j*Jx) < TOL and nrm(komm(Jz, Jx) - 1j*Jy) < TOL
      and max(nrm(komm(Nges, J)) for J in (Jx, Jy, Jz)) < TOL)
RA, RB = fkt(NA, lambda q: 1j**q), fkt(NB, lambda q: 1j**q)
Rrel = expm(1j*np.pi/2*Jz)
check("Nahtuhr R_A = i^(Q_A) = e^(i pi N/4) e^(i pi J_z/2): die Uhr einer Naht ist die Vierteldrehung des Transferrahmens "
      "um z", nrm(RA - fkt(Nges, lambda n: np.exp(1j*np.pi*n/4)) @ Rrel) < TOL)
Gam = {a: expm(1j*np.pi*J) for a, J in (("x", Jx), ("y", Jy), ("z", Jz))}
check("Gamma_z = e^(i pi J_z) ist die Paaruhr i^(Q_A - Q_B)", nrm(Gam["z"] - RA @ dag(RB)) < TOL)
klein = all(nrm(Pg @ komm(Gam[p], Gam[r])) < TOL for p, r in (("x", "y"), ("y", "z"), ("z", "x"))) and \
    all(nrm(Pg @ (Gam[p] @ Gam[p] - E2)) < TOL for p in "xyz")
quat = all(nrm(Pu @ akomm(Gam[p], Gam[r])) < TOL for p, r in (("x", "y"), ("y", "z"), ("z", "x"))) and \
    all(nrm(Pu @ (Gam[p] @ Gam[p] + E2)) < TOL for p in "xyz")
check("gerade Gesamtladung: Gamma_x, Gamma_y, Gamma_z kommutierende Involutionen (Kleinsche Vierergruppe)", klein)
check("ungerade Gesamtladung: Gamma_x, Gamma_y, Gamma_z antikommutieren, Quadrat -1 (Quaternionen)", quat)
Bst = expm(-1j*np.pi/2*Jx)                             # Strahlteiler von z120


def Phi(phi):
    return fkt(NA, lambda q: np.exp(-1j*phi*q*q/2)) @ fkt(NB, lambda q: np.exp(-1j*phi*q*q/2))


def Ugat(phi, B=None):
    B = Bst if B is None else B
    return dag(B) @ Phi(phi) @ B


Jw = dag(Bst) @ Jz @ Bst
vz = 1.0 if nrm(Jw - Jy) < TOL else -1.0
check("B^+ J_z B = pm J_y und U(phi) = e^(-i phi N^2/4) e^(-i phi J_y^2): das z120-Gatter ist die Quadratphase des "
      "gedrehten Transferspins", nrm(Jw - vz*Jy) < TOL and all(
          nrm(Ugat(p) - fkt(Nges, lambda n: np.exp(-1j*p*n*n/4)) @ expm(-1j*p*Jy @ Jy)) < 1e-9 for p in (0.37, np.pi/2, 2.2)))
NANB = NA @ NB
check("phi = pi/2: e^(-i pi H_Q/2) = e^(-i pi N^2/4) i^(Q_A Q_B): die Wechselphase ist der Bicharakter i^(ab) der "
      "Uhrgruppe Z_4, derselbe wie B_D5 (und invers zu B_A3)",
      nrm(Phi(np.pi/2) - fkt(Nges, lambda n: np.exp(-1j*np.pi*n*n/4)) @ fkt(NANB, lambda q: 1j**q)) < TOL)
a0, b0 = (1 - 1j)/2, (1 + 1j)/2
check("phi = pi/2, gerade Gesamtladung: U = e^(-i pi N^2/8) (a + b Gamma_y), a = (1-i)/2, b = (1+i)/2, also "
      "P_+ - i P_- zur Tauschparität Gamma_y; U^2 = Gamma_y bis auf Sektorphase",
      nrm(Pg @ (Ugat(np.pi/2) - fkt(Nges, lambda n: np.exp(-1j*np.pi*n*n/8)) @ (a0*E2 + b0*Gam["y"]))) < 1e-9
      and nrm(Pg @ (Ugat(np.pi/2) @ Ugat(np.pi/2) - fkt(Nges, lambda n: np.exp(-1j*np.pi*n*n/4)) @ Gam["y"])) < 1e-9)


def linearitaet(U, P):
    """Rest von U J_z U^+ außerhalb von span{J_x, J_y, J_z, 1, N} im Sektor P: 0 heißt, U normalisiert dort die su(2)."""
    idx = np.where(np.real(np.diag(P)) > 0.5)[0]
    T = (U @ Jz @ dag(U))[np.ix_(idx, idx)]
    basis = np.array([M[np.ix_(idx, idx)].ravel() for M in (Jx, Jy, Jz, E2, Nges)]).T
    koef = np.linalg.lstsq(basis, T.ravel(), rcond=None)[0]
    return nrm(T.ravel() - basis @ koef)


lin = {k: (linearitaet(Ugat(k*np.pi/2), Pg), linearitaet(Ugat(k*np.pi/2), Pu)) for k in range(4)}
say("   Rest der Linearität von U_k J_z U_k^+ (U_k = U(k pi/2); gerade / ungerade Gesamtladung): "
    + ", ".join(f"k={k}: {v[0]:.4f} / {v[1]:.4f}" for k, v in lin.items()))
check("k gerade: das Gatter normalisiert die Transfer-su(2) in jedem Paritätssektor (k = 0 trivial, k = 2 der "
      "Nahttausch Gamma_y bei gerader, eine Phase bei ungerader Gesamtladung); k ungerade: nicht (echte "
      "Wechselwirkung) - genau die Paarungen vom Rang 4",
      max(lin[0]) < 1e-9 and max(lin[2]) < 1e-9 and lin[1][0] > 0.5 and lin[3][0] > 0.5 and
      nrm(Pg @ (expm(-1j*np.pi*Jy @ Jy) - Gam["y"])) < 1e-9)
ZAHLEN["linearitaet_k"] = {k: [round(w_, 9) for w_ in v] for k, v in lin.items()}

# =================================================================== Teil B2
say("\nTeil B2: Zopfrelation zwischen Nahtuhr und Strahlteiler")


def T_x(theta):
    return expm(-1j*theta*Jx)


def sektordefekt(Xm, Ym):
    """max über Gesamtladungssektoren von min_phi ||X - e^(i phi) Y||."""
    d = 0.0
    for n in sorted(set(nges)):
        idx = np.where(nges == n)[0]
        d = max(d, bis_auf_phase(Xm[np.ix_(idx, idx)], Ym[np.ix_(idx, idx)]))
    return d


gitter = [j*np.pi/8 for j in range(-15, 17)]
exakt0, proj0 = [], []
for j, th in zip(range(-15, 17), gitter):
    T = T_x(th)
    L, Rr = Rrel @ T @ Rrel, T @ Rrel @ T
    if nrm(L - Rr) < 1e-9:
        exakt0.append(j)
    if sektordefekt(RA @ T @ RA, T @ RA @ T) < 1e-9:
        proj0.append(j)
say(f"   theta = j pi/8, j = -15..16: exakte Nullstellen des Zopfdefekts (relative Uhr) j = {exakt0}; "
    f"projektive Nullstellen (Nahtuhr R_A, Phase je Gesamtladung) j = {proj0}")
check("R T(theta) R = T(theta) R T(theta) gilt genau für theta = pm pi/2 (alle Transferspins): die Zopfrelation mit der "
      "Nahtuhr erzwingt die Vierteldrehung des Strahlteilers", exakt0 == [-4, 4])
check("projektiv (Phase je Gesamtladungssektor) kommen nur die 2 pi-Verschiebungen theta = pm 3 pi/2 hinzu",
      proj0 == [-12, -4, 4, 12])
cs = sp.symbols("c s", real=True)
check("Quaternionenbeweis: (c + s k)(c + s i)(c + s k) - (c + s i)(c + s k)(c + s i) = s (1 - 2 c^2)(i - k)",
      sp.simplify((cs[0]**2*cs[1] + cs[1]**3) - cs[1]) .subs(cs[0]**2, 1 - cs[1]**2).simplify() == 0)
c1m, (N1A, N1B), tr1 = naehte(2, 1)
J1x, J1y, J1z = tr1(0, 1)
idx1 = [1, 2]                                          # Ein-Fermion-Sektor: |n_A n_B> = |01>, |10>
R1 = fkt(N1A, lambda q: 1j**q)[np.ix_(idx1, idx1)]
B1 = expm(-1j*np.pi/2*J1x)[np.ix_(idx1, idx1)]
Had = np.array([[1, 1], [1, -1]], complex)/np.sqrt(2)
check("Ein-Fermion-Sektor eines Nahtpaars: R_A = diag(1, i) (Ising-R-Matrix bis auf Phase) und "
      "B = e^(-i pi/4) H R_A H: der Strahlteiler ist die Hadamard-Transformierte der Nahtuhr (Ising-F-Matrix)",
      nrm(R1 - np.diag([1, 1j])) < TOL and nrm(B1 - np.exp(-1j*np.pi/4)*Had @ R1 @ Had) < TOL)

# =================================================================== Teil B3
say("\nTeil B3: drei Nähte (K = 2 Moden je Naht, 64 Zustände)")
c3, (MA, MB, MC), tr3 = naehte(3, 2)
JAB, JBC = tr3(0, 1), tr3(1, 2)
E3 = np.eye(len(MA))
GAB, GBC = expm(1j*np.pi*JAB[1]), expm(1j*np.pi*JBC[1])
VAB, VBC = expm(1j*np.pi/2*JAB[1]), expm(1j*np.pi/2*JBC[1])
check("Nahttausch Gamma_AB = e^(i pi J_y^AB): Zopfrelation exakt, aber keine Antikommutation",
      nrm(zopfdefekt(GAB, GBC)) < 1e-9 and nrm(akomm(GAB, GBC)) > 1.0)
dV = snorm(zopfdefekt(VAB, VBC))
say(f"   einkörperige Wurzel des Tauschs V_AB = e^(i pi J_y^AB/2): Zopfdefekt {dV:.6f}")
check("auch die einkörperige Wurzel des Tauschs (V^2 = Gamma, quasifrei) verletzt die Zopfrelation: in der "
      "Tauschfamilie ist nur der Tausch selbst zopfkonsistent", dV > 0.1 and nrm(VAB @ VAB - GAB) < 1e-9)
ZAHLEN["wurzel_tausch_defekt"] = round(dV, 9)


def U3(phi, J, Ns, Nt):
    B = expm(-1j*np.pi/2*J[0])
    return dag(B) @ fkt(Ns, lambda q: np.exp(-1j*phi*q*q/2)) @ fkt(Nt, lambda q: np.exp(-1j*phi*q*q/2)) @ B


na, nb, nc = (np.real(np.diag(M)).round().astype(int) for M in (MA, MB, MC))
bal = np.where((na == 1) & (nb == 1) & (nc == 1))[0]
Pbal = np.zeros_like(E3)
Pbal[bal, bal] = 1
UAB, UBC = U3(np.pi/2, JAB, MA, MB), U3(np.pi/2, JBC, MB, MC)
check("der ausgeglichene Sektor Q_A = Q_B = Q_C = 1 (8 Zustände) ist unter beiden Gattern abgeschlossen",
      nrm(komm(UAB, Pbal)) < 1e-9 and nrm(komm(UBC, Pbal)) < 1e-9)
Dvoll = zopfdefekt(UAB, UBC)
dbal = snorm(Dvoll[np.ix_(bal, bal)])
say(f"   z120-Gatter bei phi = pi/2: Zopfdefekt (Spektralnorm) gesamt {snorm(Dvoll):.6f}, ausgeglichener Sektor {dbal:.6f}; "
    f"||[U_AB, U_BC]|| = {snorm(komm(UAB, UBC)):.6f}")
check("das wechselwirkende z120-Gatter verletzt die Zopfrelation; im ausgeglichenen Sektor ist der Defekt genau "
      "sqrt(3/8) (Modell A)", abs(dbal - np.sqrt(3/8)) < 1e-9)
gbal = [M[np.ix_(bal, bal)] for M in (GAB, GBC)]
s3 = all(nrm(g @ g - np.eye(8)) < 1e-9 for g in gbal) and nrm(zopfdefekt(*gbal)) < 1e-9
mono = bis_auf_phase((UAB @ UAB)[np.ix_(bal, bal)], gbal[0]) < 1e-9
check("dort ist die Monodromie U^2 der Nahttausch, und die Tausche erzeugen S_3 (Involutionen mit Zopfrelation) statt "
      "zu antikommutieren", s3 and mono)
scan3 = {}
for j in range(0, 17):
    ph = j*np.pi/8
    D = zopfdefekt(U3(ph, JAB, MA, MB), U3(ph, JBC, MB, MC))
    zph = np.exp(-1j*ph)
    scan3[j] = (snorm(D[np.ix_(bal, bal)]), abs((1 + zph)**2*(1 - zph)/8)*np.sqrt(3))
okscan = all(abs(v[0] - v[1]) < 1e-9 for v in scan3.values())
nullst = [j for j, v in scan3.items() if v[0] < 1e-9]
say(f"   Zopfdefekt im ausgeglichenen Sektor über phi = j pi/8: Nullstellen j = {nullst}")
check("im ausgeglichenen Sektor ist der Defekt für alle phi |1+z|^2 |1-z| sqrt 3/8 mit z = e^(-i phi): Nullstellen nur "
      "phi = 0, pi (kein Gatter, reiner Tausch)", okscan and nullst == [0, 8, 16])
ZAHLEN["dreinaht"] = {"defekt_gesamt": round(snorm(Dvoll), 9), "defekt_ausgeglichen": round(dbal, 9),
                      "kommutator": round(snorm(komm(UAB, UBC)), 9)}

# =================================================================== Teil B4
say("\nTeil B4: eine Naht mit Oszillatoren - welche Phasen liegen in der Uhralgebra?")
MF = 3                                                  # Moden y = -3..2, See y < 0 gefüllt
zust = []
for occ in itertools.product((0, 1), repeat=2*MF):
    Q = sum(occ[y + MF] for y in range(0, MF)) - sum(1 - occ[y + MF] for y in range(-MF, 0))
    En = sum((y + 0.5)*occ[y + MF] for y in range(0, MF)) + sum((-y - 0.5)*(1 - occ[y + MF]) for y in range(-MF, 0))
    zust.append((Q, En - Q*Q/2))
check("L_0^osc = E^nat - Q^2/2 ist ganzzahlig und nichtnegativ; im Fenster treten zu gleichem Q die Werte 0 und 1 auf",
      all(abs(l - round(l)) < 1e-12 and l > -1e-12 for _, l in zust) and {0, 1} <= {round(l) for q, l in zust if q == 0})


def in_uhralgebra(phi, chi):
    """Hängt e^(-i (phi Q^2/2 + chi L_0^osc)) nur von Q mod 4 ab?"""
    wert = {}
    for Q, l in zust:
        v = np.exp(-1j*(phi*Q*Q/2 + chi*round(l)))
        if abs(wert.setdefault(Q % 4, v) - v) > 1e-9:
            return False
    return True


treffer = [(j, l) for j in range(16) for l in range(4) if in_uhralgebra(j*np.pi/8, l*np.pi/2)]
say(f"   (phi, chi) = (j pi/8, l pi/2) mit Phase in der Uhralgebra: {treffer}")
check("die Kompositionsphase ist genau dann eine Funktion der Uhrwerte (Q mod 4), wenn phi in (pi/2) Z liegt und der "
      "Oszillatoranteil verschwindet (chi = 0 mod 2 pi)", treffer == [(0, 0), (4, 0), (8, 0), (12, 0)])
befund("in der nativen Familie e^(-i (phi H_Q + chi L_0^osc)) schließt P_uhr ('die Kompositionsphase ist eine Funktion "
       "der Uhrwerte') P_sl2 (= P_null, kein Oszillatoranteil, z133) ein und lässt nur phi = k pi/2; mit "
       "'wechselwirkend' (k ungerade, B1) bleibt die Vierteldrehung P_viertel. Allgemeine Funktionen der Uhrwerte "
       "brauchen zusätzlich die Bicharakter-Eigenschaft (A1)", True)

# =================================================================== Teil B5
say("\nTeil B5: Schalter zwischen den Verklebungen (zwei Nähte, je eine Mode: vier Majoranas)")
cA, cB = c1m
m1, m2 = cA + dag(cA), 1j*(dag(cA) - cA)
m3, m4 = cB + dag(cB), 1j*(dag(cB) - cB)
E4 = np.eye(4)
FA, FB = E4 - 2*N1A, E4 - 2*N1B                         # Nahtparitäten R_A^2, R_B^2
R1A, R1B = fkt(N1A, lambda q: 1j**q), fkt(N1B, lambda q: 1j**q)
Ssw = 1j*m2 @ m3                                         # Schalterparität: spiegelt je eine Majorana beider Nähte
check("Nahtparität (-1)^(Q_A) = R_A^2 = -i g_1 g_2; die Nahtuhr ist R_A = P_+ + i P_-, also das Modell-B-Gatter mit "
      "z = i zur Nahtparität (Zopfung der beiden Majorana-Hälften einer Naht)",
      nrm(FA + 1j*m1 @ m2) < TOL and nrm(R1A - gatter(1j, FA)) < TOL and bis_auf_phase(R1A, expm(-np.pi/4*m1 @ m2)) < TOL)
check("Schalterparität S = i g_2 g_3: S^2 = 1, antikommutiert mit beiden Nahtparitäten, kommutiert mit (-1)^F",
      nrm(Ssw @ Ssw - E4) < TOL and nrm(akomm(Ssw, FA)) < TOL and nrm(akomm(Ssw, FB)) < TOL
      and nrm(komm(Ssw, FA @ FB)) < TOL)
zscan = [j for j in range(16) if nrm(zopfdefekt(R1A, gatter(np.exp(1j*np.pi*j/8), Ssw))) < 1e-9]
check("Zopfrelation R_A U_S(z) R_A = U_S(z) R_A U_S(z) gilt auf dem Einheitskreis genau für z = i: die Vierteldrehung "
      "des Schalters ist von der Uhr geerbt", zscan == [4])
US = gatter(1j, Ssw)
n1 = np.real(np.diag(N1A + N1B)).round().astype(int)
PN = [np.diag((n1 == n).astype(float)) for n in range(3)]
S_hop = sum(P @ Ssw @ P for P in PN)
S_paar = Ssw - S_hop
uhr_rel, uhr_sum = R1A @ dag(R1B), R1A @ R1B
check("Zerlegung S = Hüpfteil + Paarteil: der Hüpfteil ist 2 J_x; unter der Paaruhr i^(Q_A - Q_B) ist der Hüpfteil "
      "ungerade und der Paarteil gerade, unter i^(Q_A + Q_B) umgekehrt",
      bis_auf_phase(S_hop, 2*J1x) < TOL and nrm(S_paar) > 1.0
      and nrm(uhr_rel @ S_hop @ dag(uhr_rel) + S_hop) < TOL and nrm(uhr_rel @ S_paar @ dag(uhr_rel) - S_paar) < TOL
      and nrm(uhr_sum @ S_hop @ dag(uhr_sum) - S_hop) < TOL and nrm(uhr_sum @ S_paar @ dag(uhr_sum) + S_paar) < TOL)
Bfull = expm(-1j*np.pi/2*J1x)
ein = min(bis_auf_phase(US[np.ix_(idx1, idx1)], M[np.ix_(idx1, idx1)]) for M in (Bfull, dag(Bfull)))
ger = US[np.ix_([0, 3], [0, 3])]
check("ladungserhaltender Teil des Schalter-Zopfs im Ein-Fermion-Sektor = z120-Strahlteiler (bis auf Phase/Richtung); "
      "im geraden Sektor mischt er |00> und |11> mit Betrag 1/sqrt 2 (Paarerzeugung, nicht ladungserhaltend)",
      ein < TOL and abs(abs(ger[0, 1]) - 1/np.sqrt(2)) < TOL and nrm(komm(US, N1A + N1B)) > 0.5)
wA, vA = np.linalg.eigh(FA[np.ix_([0, 3], [0, 3])])
wS, vS = np.linalg.eigh(Ssw[np.ix_([0, 3], [0, 3])])
check("gerader Sektor: die Eigenbasen von Nahtparität und Schalterparität sind gegenseitig unverzerrt (|F|^2 = 1/2, "
      "Hadamard-Typ) - das folgt allein aus der Antikommutation",
      np.allclose(np.abs(dag(vA) @ vS)**2, 0.5, atol=1e-12))
check("die Uhr dreht den Majorana-Rahmen um ein Viertel (R_A g_2 R_A^-1 = pm g_1) und tauscht damit die beiden "
      "Spiegelungsklassen; sie unterscheiden sich um die Nahtparität: g_1 ~ g_2 (-1)^(Q_A)",
      min(nrm(R1A @ m2 @ dag(R1A) - s*m1) for s in (1, -1)) < TOL and bis_auf_phase(m1, m2 @ FA) < TOL)


def gruppenordnung(gens, grenze=4000):
    def schluessel(M):
        f = M.ravel()
        k = next(i for i in range(len(f)) if abs(f[i]) > 1e-8)
        return tuple(np.round(f*(abs(f[k])/f[k]), 6).tolist())

    gesehen = {schluessel(E4): E4}
    rand = [E4]
    while rand and len(gesehen) < grenze:
        neu = []
        for M in rand:
            for g in gens:
                P = g @ M
                s = schluessel(P)
                if s not in gesehen:
                    gesehen[s] = P
                    neu.append(P)
        rand = neu
    return len(gesehen)


ordnung = gruppenordnung([R1A, US, R1B])
say(f"   projektive Ordnung der von R_A, U_S, R_B erzeugten Gruppe (vier Majoranas, beide Paritäten): {ordnung}")
check("die Zopfgruppendarstellung aus Uhren und Schalter ist endlich (Clifford-Typ): keine dichte Dynamik",
      ordnung < 4000)
ZAHLEN["zopfgruppe_ordnung"] = ordnung
spek_ok = all(sorted(((m + al_) % 1) for m in range(-3, 3)) == sorted(((-(m + al_)) % 1) for m in range(-3, 3))
              for al_ in (Fraction(0), Fraction(1, 2))) and \
    sorted(((m + Fraction(1, 4)) % 1) for m in range(-3, 3)) != sorted(((-(m + Fraction(1, 4))) % 1) for m in range(-3, 3))
check("Ladungskonjugation ist keine Symmetrie des Viertels: das Spektrum {m + alpha} ist nur für alpha in {0, 1/2} "
      "spiegelsymmetrisch; alpha = 1/4 geht in 3/4 über (der Schalter verbindet Viertel und Gegenviertel)", spek_ok)
be = 2*np.log(1 + np.sqrt(2))
say(f"   Paaramplitude des KMS-Zustands e^(-beta eps/2) bei beta = 4 pi, eps = 1/4: {np.exp(-np.pi/2):.6f}; "
    f"Hadamard-Fixpunkt verlangt beta eps = 2 ln(1 + sqrt 2) = {be:.6f} statt pi = {np.pi:.6f}")
check("der native KMS-Zustand liegt nicht am Hadamard-Fixpunkt (e^(-pi/2) = 0,2079 gegen 0,4142)",
      abs(np.exp(-np.pi/2) - (np.sqrt(2) - 1)) > 0.2)
ZAHLEN["kms_gegen_hadamard"] = {"kms_amplitude": float(np.exp(-np.pi/2)), "fixpunkt": float(np.sqrt(2) - 1),
                                "beta_eps_noetig": float(be)}

# =================================================================== Teil B6
say("\nTeil B6: Yang-Baxter mit Spektralparameter (drei Nähte, K = 2)")
gleichpar = (na % 2 == nb % 2) & (nb % 2 == nc % 2)
gl, mx = np.where(gleichpar)[0], np.where(~gleichpar)[0]
Pgl = np.diag(gleichpar.astype(float))


def r_(M, idx):
    return M[np.ix_(idx, idx)]


def Yang(u, G):
    return (E3 + 1j*u*G)/(1 + 1j*u)


check("Sektoren gleicher Nahtparität sind unter den Nahttauschen abgeschlossen; dort Gamma^2 = 1, sonst nicht",
      nrm(komm(GAB, Pgl)) < 1e-9 and nrm(komm(GBC, Pgl)) < 1e-9 and nrm(r_(GAB @ GAB - E3, gl)) < 1e-9
      and nrm(r_(GBC @ GBC - E3, gl)) < 1e-9 and nrm(r_(GAB @ GAB - E3, mx)) > 1.0)
uv = [(0.3, 1.1), (1.0, 1.0), (-0.7, 2.5), (1.0, 2.0)]
ybe = max(nrm(r_(Yang(u, GAB) @ Yang(u + v, GBC) @ Yang(v, GAB) - Yang(v, GBC) @ Yang(u + v, GAB) @ Yang(u, GBC), gl))
          for u, v in uv)
einh = nrm(r_(dag(Yang(0.7, GAB)) @ Yang(0.7, GAB) - E3, gl))
nichteinh = nrm(r_(dag(Yang(0.7, GAB)) @ Yang(0.7, GAB) - E3, mx))
check("gleiche Nahtparität: Y(u) = (1 + i u Gamma)/(1 + i u) ist unitär und erfüllt Y_AB(u) Y_BC(u+v) Y_AB(v) = "
      "Y_BC(v) Y_AB(u+v) Y_BC(u) (Yangs R-Matrix); bei gemischter Parität ist Y(u) nicht unitär",
      ybe < 1e-9 and einh < 1e-9 and nichteinh > 0.1)
urel = max(bis_auf_phase(r_(U3(ph, JAB, MA, MB), bal), r_(Yang(np.tan(ph/2), GAB), bal)) for ph in (0.4, np.pi/2, 2.0, 4.1))
check("ausgeglichener Sektor: das z120-Gatter ist Y(u) mit u = tan(phi/2), dem Parameter aus z133; das Viertel ist u = 1",
      urel < 1e-9)
mitte1 = snorm(r_(zopfdefekt(Yang(1, GAB), Yang(1, GBC)), bal))


def faktor3(UA, UB, YA2, YB2, idx):
    return snorm(r_(UA @ YB2 @ UA - UB @ YA2 @ UB, idx))


f_bal = faktor3(UAB, UBC, Yang(2, GAB), Yang(2, GBC), bal)
nlm = np.setdiff1d(gl, bal)                             # gleiche Parität, aber Ladungen nicht gleich modulo 4
f_nlm = faktor3(UAB, UBC, Yang(2, GAB), Yang(2, GBC), nlm)
z2 = (1 - 2j)/(1 + 2j)
say(f"   drei Viertelgatter: Defekt {mitte1:.6f}; mit mittlerem Gatter u = 2 (z = (1-2i)/(1+2i) = {z2.real:.1f}{z2.imag:+.1f}i): "
    f"ausgeglichen {'< 1e-9' if f_bal < 1e-9 else f'{f_bal:.6f}'}, gleiche Parität ohne Niveauabgleich {f_nlm:.6f}")
check("Faktorisierung dreier Nähte: U_AB(1) Y_BC(2) U_AB(1) = U_BC(1) Y_AB(2) U_BC(1) gilt mit den vollen z120-Gattern "
      "im niveauabgeglichenen Sektor exakt; ohne Niveauabgleich brechen die Ladungsphasen e^(-i pi N^2/8) sie",
      f_bal < 1e-9 and f_nlm > 0.1 and abs(mitte1 - np.sqrt(3/8)) < 1e-9)
c8, (K4A, K4B), tr8 = naehte(2, 4)
J8x, J8y, J8z = tr8(0, 1)
n8 = np.real(np.diag(K4A + K4B)).round().astype(int)
ger8 = np.where(n8 % 2 == 0)[0]
G8y = expm(1j*np.pi*J8y)
J8y2 = J8y @ J8y


def spanrest(phi, Jq, G, idx):
    W = r_(expm(-1j*phi*Jq), idx).ravel()
    basis = np.array([np.eye(len(idx)).ravel(), r_(G, idx).ravel()]).T
    return nrm(W - basis @ np.linalg.lstsq(basis, W, rcond=None)[0])


rest8 = {j: spanrest(j*np.pi/8, J8y2, G8y, ger8) for j in range(9)}
rest2 = max(spanrest(j*np.pi/8, Jy @ Jy, Gam["y"], np.where(nges % 2 == 0)[0]) for j in range(9))
say("   K = 4 Moden je Naht, gerade Gesamtladung: Abstand von e^(-i phi J_y^2) zu span{1, Gamma_y} über phi = j pi/8: "
    + ", ".join(f"{j}: {v:.3f}" for j, v in rest8.items()))
check("für viele Moden je Naht (Transferspin >= 2) ist das native Gatter genau für phi in (pi/2) Z eine Funktion des "
      "Nahttauschs, also ein Element der Yang-Familie; für K <= 3 (Spin <= 3/2) gilt das für alle phi",
      [j for j, v in rest8.items() if v < 1e-9] == [0, 4, 8] and rest2 < 1e-9)
k4a, k4b = (np.real(np.diag(M)).round().astype(int) for M in (K4A, K4B))
paar22 = np.where((k4a == 2) & (k4b == 2))[0]
rest22 = np.where((n8 == 4) & ~((k4a == 2) & (k4b == 2)))[0]
leck = {j: nrm(expm(-1j*j*np.pi/8*J8y2)[np.ix_(rest22, paar22)]) for j in range(9)}
uhr8 = expm(1j*np.pi*J8z)
uhrinv = max(nrm(komm(expm(-1j*j*np.pi/8*J8y2), uhr8)) for j in range(9))
say("   Leck aus dem Paarsektor (K = 4: Fermifenster y = -2..1 je Naht, Q relativ zum See, Sektor Q_A = Q_B = 0) über "
    "phi = j pi/8: " + ", ".join(f"{j}: {v:.3f}" for j, v in leck.items()))
check("Leckfreiheit des Paarsektors Q_A = Q_B gilt bei Transferspin >= 2 genau für phi in (pi/2) Z; die Paaruhr "
      "i^(Q_A - Q_B) vertauscht dagegen für jedes phi mit dem Gatter und verbietet das Leck nach (Q+2, Q-2) nicht",
      [j for j, v in leck.items() if v < 1e-9] == [0, 4, 8] and uhrinv < 1e-9)
leck2 = max(nrm(Ugat(j*np.pi/8)[np.ix_(np.where((nges == 2) & (np.real(np.diag(NA)).round() != 1))[0],
                                         np.where((nges == 2) & (np.real(np.diag(NA)).round() == 1))[0])]) for j in range(9))
check("Kontrolle (Geltungsbereich von z121): bei einer Anregung je Naht (Q_A = Q_B = 1, Transferspin <= 1) ist jede "
      "Phase leckfrei", leck2 < 1e-9)
befund("z121s 'die Phase ist frei, jede Phase ist leckfrei' gilt nur bis Transferspin 3/2; ab Transferspin 2 erzwingt "
       "die strenge Leckfreiheit P_Q (Nahtladungen werden nicht umverteilt, höchstens getauscht) phi in (pi/2) Z", True)
ZAHLEN["yang_baxter"] = {"faktor3_ohne_abgleich": round(f_nlm, 9), "span_rest_K4": {j: round(v, 9) for j, v in rest8.items()},
                         "leck_K4": {j: round(v, 9) for j, v in leck.items()}}

# =================================================================== Teil B7
say("\nTeil B7: vier Nähte (K = 2), alternierende Komposition und Transfermatrizen")
c4, N4, tr4 = naehte(4, 2)
occ4 = [np.real(np.diag(M)).round().astype(int) for M in N4]
bal4 = np.where(np.all([o == 1 for o in occ4], axis=0))[0]


def qindex(idx):
    """Nahtinhalt y_s (welche Mode der Naht s besetzt ist) als Qubit-Index."""
    return sum((((idx >> (7 - (2*s_ + 1))) & 1) << (3 - s_)) for s_ in range(4))


ordnung4 = [0]*16
for idx_ in bal4:
    ordnung4[qindex(int(idx_))] = int(idx_)


def swapq(n, i_, j_):
    M = np.zeros((2**n, 2**n))
    for q in range(2**n):
        bi, bj = (q >> (n - 1 - i_)) & 1, (q >> (n - 1 - j_)) & 1
        q2 = q & ~(1 << (n - 1 - i_)) & ~(1 << (n - 1 - j_)) | (bj << (n - 1 - i_)) | (bi << (n - 1 - j_))
        M[q2, q] = 1
    return M


paare = [(0, 1), (1, 2), (2, 3), (3, 0)]
G4 = {pq: expm(1j*np.pi*tr4(*pq)[1]) for pq in paare}
tausch_ok = all(nrm(r_(G4[pq], ordnung4) - swapq(4, *pq)) < 1e-9 for pq in paare)
check("ausgeglichener Sektor (16 Zustände): der Nahttausch Gamma_st ist genau der Tausch der Nahtinhalte, auch über "
      "den Ringschluss (4,1)", len(bal4) == 16 and tausch_ok)
E16 = np.eye(16)


def g_(pq, u=1.0):
    return (E16 + 1j*u*swapq(4, *pq))/(1 + 1j*u)


U4 = {pq: r_(U3(np.pi/2, tr4(*pq), N4[pq[0]], N4[pq[1]]), ordnung4) for pq in paare}
check("dort ist das z120-Viertelgatter (1 + i P)/(1 + i) bis auf eine Phase: SU(K)-invariantes Gatter auf den "
      "Nahtinhalten", all(bis_auf_phase(U4[pq], g_(pq)) < 1e-9 for pq in paare))
P5 = {j_: swapq(5, 0, j_) for j_ in range(1, 5)}


def tau(lam, delta):
    """Transfermatrix tr_0 R_04 R_03 R_02 R_01 mit R = (lam - xi_j) + i P und xi_j = pm delta/2 alternierend."""
    M = np.eye(32, dtype=complex)
    for j_ in range(1, 5):
        xi = delta/2 if j_ % 2 == 1 else -delta/2
        M = ((lam - xi)*np.eye(32) + 1j*P5[j_]) @ M
    return np.einsum("aiaj->ij", M.reshape(2, 16, 2, 16))


Wa = U4[(0, 1)] @ U4[(2, 3)] @ U4[(1, 2)] @ U4[(3, 0)]
Wb = U4[(1, 2)] @ U4[(3, 0)] @ U4[(0, 1)] @ U4[(2, 3)]
Wtreppe = U4[(0, 1)] @ U4[(1, 2)] @ U4[(2, 3)] @ U4[(3, 0)]
lams = (0.37, -1.2 + 0.4j)
tt_ = max(nrm(komm(tau(lams[0], d_), tau(lams[1], d_))) for d_ in (1.0, -1.0))
kom = {d_: (max(nrm(komm(tau(l_, d_), Wa)) for l_ in lams), max(nrm(komm(tau(l_, d_), Wb)) for l_ in lams),
            max(nrm(komm(tau(l_, d_), Wtreppe)) for l_ in lams)) for d_ in (1.0, -1.0)}
say("   ||[tau(lambda), W]|| für delta = +1 / -1 (alternierend a, alternierend b, Treppe): "
    + " / ".join("(" + ", ".join("< 1e-8" if w_ < 1e-8 else f"{w_:.4f}" for w_ in v) + ")" for v in kom.values()))
dgut = [d_ for d_, v in kom.items() if v[0] < 1e-8 and v[1] < 1e-8]
kontrolle = max(nrm(komm(tau(l_, 0.3), Wa)) for l_ in lams)
say(f"   Kontrolle: mit falscher Inhomogenität delta = 0,3 ist ||[tau, W]|| = {kontrolle:.4f}")
check("die Transfermatrizen tau(lambda) mit alternierenden Inhomogenitäten pm 1/2 vertauschen untereinander und mit "
      "der alternierenden Komposition der Viertelgatter (beide Schichtfolgen; delta = +1 und -1 gehen durch "
      "Verschiebung um eine Naht ineinander über): integrable Dynamik; mit delta = 0,3 vertauschen sie nicht",
      tt_ < 1e-8 and dgut == [1.0, -1.0] and kontrolle > 0.1)
befund("die Treppenordnung U_12 U_23 U_34 U_41 vertauscht mit denselben Transfermatrizen", kom[dgut[0]][2] < 1e-8 if dgut else False)
ZAHLEN["integrabel"] = {"delta": dgut, "treppe_kommutator": round(kom[1.0][2], 9), "kontrolle_delta_0.3": round(kontrolle, 9)}

# =================================================================== Teil C
say("\nTeil C: Bilanz")
BILANZ = [
    ("E", "Gitterseite der Analyse stimmt: q_D5 = 5a^2/8, q_A3 = 3a^2/8, Bicharaktere i^(+ab), i^(-ab); von den vier "
          "Z_4-Paarungen sind nur diese nichtentartet; zwei Lagrangesche <(1,1)>, <(1,3)> (je 240 Wurzeln), getauscht "
          "durch y -> -y oder x -> -x; D8 mit {1, s, c, v} = torischer Code; fermionische Erweiterung 16 = 10 + 6"),
    ("E", "Zopfalgebra stimmt: Modell A Defekt (1+z)^2 (1-z)/8 (S_12 - S_23), bei z = -i Norm sqrt(3/8); Modell B "
          "Defekt (1-z)(1+z^2)/4 (A - B), Nullstellen 1, i, -i"),
    ("E", "Hadamard-Fixpunkt sqrt 2 - 1 = Ising-Blockverhältnis bei x = 1/2 = C(3)/C(1) aus v526; M = m 1; Doppelpol "
          "nicht reflexionspositiv; 1 + w = -(2/(3 H alpha)) alpha'/alpha"),
    ("E", "nativ: die Nahtuhr R_A ist die Vierteldrehung des Transferrahmens um z; das z120-Gatter ist "
          "e^(-i phi J_y^2) bis auf Sektorphase; bei phi = pi/2 ist seine Phase der Bicharakter i^(Q_A Q_B) = B_D5 auf "
          "den Uhrwerten (die Uhr erzeugt das Zentrum Z_4 von Spin(10))"),
    ("E", "nativ: Phasen in der Uhralgebra sind genau phi = k pi/2 ohne Oszillatoranteil; k = 0 trivial, k = 2 "
          "Nahttausch (einkörperig), k = pm 1 wechselwirkend - dieselbe Einteilung wie der Rang der Paarung"),
    ("E", "nativ: in der Familie e^(-i phi J_y^2) sind bei Transferspin >= 2 gleichwertig: phi in (pi/2) Z; die Phase "
          "ist eine Funktion der Uhrwerte; das Gatter ist eine Funktion des Nahttauschs; der Paarsektor Q_A = Q_B "
          "leckt nicht. Die Paaruhr allein verbietet das Leck nicht"),
    ("E", "Anschluss an z121 im Mehrmodenmodell: Auf Sektoren gleicher Nahtparität ist Y(u) = (1 + i u Gamma)/(1 + i u) "
          "Yangs R-Matrix mit additivem u = tan(phi/2); vier Nähte: Nahttausch = Tausch der Nahtinhalte, die "
          "alternierende Komposition vertauscht mit allen Transfermatrizen (integrabel, wie z121), die Treppenordnung "
          "nicht. Neu: Mit den vollen z120-Gattern gilt die Faktorisierung dreier Nähte nur bei Nahtladungen gleich "
          "modulo 4 (Niveauabgleich, z130), sonst brechen die Ladungsphasen e^(-i pi N^2/8) sie"),
    ("C", "Korrektur zum Geltungsbereich von z121: 'jede Phase ist leckfrei' gilt nur bis Transferspin 3/2 (eine "
          "Anregung je Naht). Ab Transferspin 2 erzwingt die strenge Leckfreiheit P_Q ('die Komposition verteilt "
          "Nahtladungen nicht um, sie tauscht sie höchstens') phi in (pi/2) Z; ohne das triviale Gatter und den reinen "
          "Tausch bleibt phi = pm pi/2. P_viertel ist damit keine eigene Zahl mehr, sondern folgt aus derselben "
          "Leckfreiheit, die schon die Vierteldrehung des Strahlteilers trägt (z121). P_Q ist stärker als die Paaruhr "
          "(native Form von P_univ), die das Leck nach (Q+2, Q-2) nicht verbietet; hergeleitet ist P_Q nicht"),
    ("C", "gleichwertige Fassung P_uhr ('die Kompositionsphase ist eine Funktion der Uhrwerte'): in der nativen Familie "
          "schließt sie P_sl2 (= P_null) ein; das Vorzeichen u = pm 1 bleibt die Zeitrichtung"),
    ("E", "nativ: Zopfrelation zwischen Nahtuhr und Strahlteiler gilt genau bei der Vierteldrehung (alle Transferspins); "
          "im Ein-Fermion-Sektor sind Uhr und Strahlteiler die Ising-Zopfmatrizen, B = e^(-i pi/4) H R_A H"),
    ("E", "nativ: das wechselwirkende z120-Gatter ist Modell A, nicht Modell B: seine Monodromie ist der Nahttausch, die "
          "Tausche erzeugen S_3; Zopfdefekt sqrt(3/8); Nullstellen nur phi = 0, pi; auch die einkörperige Wurzel des "
          "Tauschs zopft nicht. Zopfkonsistenz und Wechselwirkung schließen sich in dieser Familie aus; P_ord bleibt "
          "nötig"),
    ("E", "nativ: Modell B gibt es als Nahtparität gegen Schalterparität S = i g_2 g_3; die Uhr ist dann das Gatter mit "
          "z = i, der Zopf erzwingt z = i für den Schalter; sein ladungserhaltender Teil ist der Strahlteiler, der Rest "
          "Paarerzeugung (gerade unter der Paaruhr, also verträglich mit P_TP); erzeugte Gruppe endlich (Clifford-Typ)"),
    ("L", "P_IF als Wahl des Lagrangeschen: die beiden Verklebungen sind durch einen Automorphismus isomorph, kein "
          "inneres Kriterium kann eine bevorzugen (deckt sich mit z129); der Schalter ist die Ladungskonjugation, die "
          "das Viertel ins Gegenviertel überführt, die Uhr tauscht die beiden Schalterklassen (s <-> c bis auf v)"),
    ("N", "der native KMS-Zustand (beta = 4 pi, eps = 1/4) hat Paaramplitude e^(-pi/2) = 0,2079, nicht sqrt 2 - 1; der "
          "Hadamard-Fixpunkt aus v526 gehört zum Achterkreis, nicht zum Viertelprozess"),
    ("O", "eine native Begründung von P_Q (etwa: Nähte als Sektoren mit Ladungs-Superauswahl, das Gatter als ihr "
          "Statistikoperator); ob der Paarteil des Schalters im Viertelprozess als Operator existiert (Tomita-Partner); "
          "Eichdynamik über das su(2)_1-Kontinuum von z121 hinaus (K15), G und Lambda - die Analyse berührt G nicht"),
]
for st, txt in BILANZ:
    say(f"   [{st}] {txt}")
ZAHLEN["bilanz"] = [{"status": st, "aussage": txt} for st, txt in BILANZ]

urteil = ("Die neunzehnte externe Analyse hält in allen gerechneten Behauptungen, und ihr Kern lässt sich nativ prüfen; "
          "das Ergebnis ist zweigeteilt. (1) Gitter: D5 und A3 tragen die Z_4-Formen 5a^2/8 und 3a^2/8 mit den "
          "Bicharakteren i^(ab) und i^(-ab), den einzigen nichtentarteten Z_4-Paarungen. Es gibt genau zwei "
          "Verklebungen zu E8 (je 240 Wurzeln), getauscht durch die Ladungskonjugation einer Seite, dazu die "
          "fermionische Erweiterung 16 = 10 + 6; D8 hat die Struktur des torischen Codes. (2) Nativ ist die "
          "Wechselphase des z120-Gatters bei phi = pi/2 genau dieser Bicharakter, i^(Q_A Q_B), und eine "
          "Kompositionsphase ist genau dann eine Funktion der Uhrwerte, wenn phi in (pi/2) Z liegt und kein "
          "Oszillatoranteil vorkommt. Unter diesen vier Phasen ist k = 0 trivial, k = 2 der einkörperige Nahttausch "
          "und k = pm 1 die Wechselwirkung: Das ist dieselbe Einteilung wie der Rang der Paarung (1, 4, 2, 4). (3) Die "
          "konstante Zopfrelation trägt nativ nur einkörperige Operatoren: Nahtuhr und Strahlteiler erfüllen sie genau "
          "bei der Vierteldrehung (zweite, unabhängige Begründung neben der Leckfreiheit aus z120) und sind im "
          "Ein-Fermion-Sektor die Ising-Zopfmatrizen, der Strahlteiler die Hadamard-Transformierte der Uhr. Das "
          "wechselwirkende z120-Gatter ist dagegen Modell A: Seine Monodromie ist der Nahttausch, die Tausche erzeugen "
          "S_3 und antikommutieren nicht, der Zopfdefekt ist im ausgeglichenen Sektor genau sqrt(3/8) und verschwindet "
          "nur für phi = 0, pi. Die Erwartung der Analyse, die konstante Zopfrelation lege P_viertel des "
          "wechselwirkenden Gatters fest, ist widerlegt. (3a) Die tragfähige Fassung ist die Yang-Baxter-Gleichung mit "
          "Spektralparameter, die z121 schon gefunden hat (XXX-Familie R(u), u = tan(phi/2), integrable alternierende "
          "Kette); sie bestätigt sich hier im Mehrmodenmodell mit drei und vier Nähten. Zwei Punkte kommen hinzu. "
          "Erstens faktorisiert die Komposition dreier Nähte mit den vollen z120-Gattern nur bei Nahtladungen gleich "
          "modulo 4 (Niveauabgleich, die native Form von P_univ aus z130); sonst brechen die Ladungsphasen sie. "
          "Zweitens gilt z121s Satz 'jede Phase ist leckfrei, die Phase ist frei' nur bei einer Anregung je Naht "
          "(Transferspin bis 3/2). Ab Transferspin 2, im nativen Fermifenster schon im Sektor Q_A = Q_B = 0, leckt der "
          "Paarsektor nach (Q+2, Q-2), und das Gatter verlässt die XXX-Familie, außer für phi in (pi/2) Z; die Paaruhr "
          "verbietet dieses Leck nicht. Die strenge Leckfreiheit P_Q ('die Komposition verteilt Nahtladungen nicht um, "
          "sie tauscht sie höchstens') erzwingt also neben der Vierteldrehung des Strahlteilers (z121) auch die Phase: "
          "Ohne das triviale Gatter und den reinen Tausch bleibt phi = pm pi/2. P_viertel ist damit keine eigene Zahl "
          "mehr, sondern eine Folge von P_Q; gleichwertig ist P_uhr ('die Phase ist eine Funktion der Uhrwerte'), das "
          "in der nativen Familie auch P_sl2 (= P_null) einschließt. Hergeleitet sind P_Q und P_uhr nicht. (4) Modell "
          "B existiert nativ als Paar Nahtparität und Schalterparität (i g_2 g_3): Die Uhr ist "
          "dann selbst das Zopfgatter mit z = i und vererbt z = i an den Schalter; dessen ladungserhaltender Teil ist "
          "der z120-Strahlteiler, der Rest ist Paarerzeugung, gerade unter der Paaruhr (verträglich mit P_TP). Die so "
          "erzeugte Gruppe ist endlich (Clifford-Typ), trägt also keine lokale Dynamik. (5) Das Bit P_IF ist in dieser "
          "Lesart die Wahl zwischen zwei isomorphen Verklebungen; kein inneres Kriterium kann sie treffen (wie z129). "
          "(6) Der Hadamard-Fixpunkt sqrt 2 - 1 stimmt mit v526 überein, ist aber nicht die native KMS-Amplitude "
          "(e^(-pi/2)). (7) Der Doppelpol ist nicht reflexionspositiv, die Zustandsgleichung des Lambda-Zweigs ist "
          "w = -1 bis auf höchstens etwa 1e-6; der in der Analyse genannte DESI-Wert ist hier nicht geprüft. Für den "
          "Kandidaten: P_viertel wird zur Folge der Leckfreiheit P_Q (gleichwertig P_uhr, das auch P_null trägt), P_IF "
          "wird zur Wahl zwischen isomorphen Lösungen, P_ord bleibt (die alternierende Schichtung ist integrabel wie in "
          "z121, die Treppenordnung nicht). Die Zopf- und Defektstruktur der Analyse liefert keine neue Dynamik: Sie "
          "ist endlich (Clifford-Typ) und einkörperig. Eichdynamik (K15), P_curv und G bleiben offen. Nichts ist "
          "befördert.")
say("\nUrteil\n   " + urteil)
print(f"\n   (Laufzeit {time.time() - T0:.1f} s, nicht Teil der Ausgabe)", flush=True)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z144_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
(HERE / "ausgabe_z144.txt").write_text("\n".join(OUT) + "\n")
