"""Z165: Einundzwanzigste externe Analyse nachgerechnet; Folgerechnungen zum Kandidaten 'relative graduierte Determinante'.

Leitlinie 'Rekonstruktion statt Auswahl'; H6 und die gelösten Blöcke (Flavour, Massen, alpha, 16, E8) bleiben unberührt,
nichts wird befördert. Statuskürzel: E exakt, C bedingt (Prämisse genannt), L Lesart, N numerisch, O offen.

Anlass: Die externe Analyse (zum Stand z163) prüft die Wege W1 bis W5 und schlägt als ein Objekt eine relative
graduierte Determinante über dem Raum der Twistorgeraden vor. Sie rechnet: (1) Doppelpol verletzt
Reflexionspositivität, (2) Vertexkern (1 - z-quer w)^2 ist nicht positiv, (3) q = dz^2/(z^4 - 1) hat Charakter -1,
(4) C = R^2 U B U für die beiden Markenmatrizen, (5) Lambda G = 3/(32 pi^3) e^(-2/alpha), (6) magnetische Sektoren in
der SU(5)-Einbettung, (7) f_0 f_4 >= f_2^2 für ein positives Spektralprofil, (8) ALE-Auflösung und AdS sind nicht
konform verträglich.

Bestand: z164 (eine freie Zahl; nativer Unendlichkeitstwistor vom Einstein-Typ), z163, z150 (Yang-Mills-artige
Gravitation), PRUEFUNG-BESTAND-H6h34 (W1: q hat Charakter -1), Ledger ALPHA.QUILLEN.EXACT.01 und SEAM.CONTACT.UNIT.01
(R_0 = 16 c3 circ(2,-1,0,-1), G_marks = -4 c3 ln2 circ(0,1,2,1), 'swapped mu4 channels'),
SEAM.DETLINE.UNIFICATION.01, GRAV.SPIN2.EMERGENCE.01.

Teil A  Nachrechnung der acht Aussagen
Teil B  Folgerechnungen: (B1) Gravitation ohne Doppelpol: Vertex B.B mit so(3,2)-Krümmung; (B2) Lambda-Stufe als
        relative Nahtantwort; (B3) Ungleichung gegen die Ledgerzahlen; (B4) native Graduierung: Uhr-Spuren auf E8
Teil C  Bilanz und Urteil.
"""
import itertools
import json
import math
import os
import time
from fractions import Fraction as Fr
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import mpmath as mp
import sympy as sp

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
T0 = time.time()
mp.mp.dps = 40


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val, text):
    BEFUNDE[name] = val
    say(f"   -> {name}: {text}")


say("Z165: externe Analyse 21 nachgerechnet; relative graduierte Determinante")
say("=" * 100)
say("\nTeil A  Nachrechnung")
# A1 Doppelpol
A, B, om, t1, t2 = sp.symbols("A B omega t1 t2", real=True)
Ck = lambda t: sp.exp(-om * t) * (A + B * t)
detM = sp.simplify(Ck(2 * t1) * Ck(2 * t2) - Ck(t1 + t2) ** 2)
check("A1 Kern e^(-omega t)(A + B t): det der gespiegelten 2x2-Matrix = -B^2 (t1 - t2)^2 e^(-2 omega (t1 + t2)); "
      "für B = 0 ist sie null (Rang 1, positiv semidefinit) [E]",
      sp.simplify(detM + B ** 2 * (t1 - t2) ** 2 * sp.exp(-2 * om * (t1 + t2))) == 0)
# A2 Vertexkern
K1 = lambda z, w: 1 / (1 - z * w) ** 2
K2 = lambda z, w: (1 - z * w) ** 2
pts = [sp.Integer(0), sp.Rational(1, 2)]
d1 = sp.Matrix(2, 2, lambda i, j: K1(pts[i], pts[j])).det()
d2 = sp.Matrix(2, 2, lambda i, j: K2(pts[i], pts[j])).det()
xx = sp.Symbol("x")
c_minus = [sp.series((1 - xx) ** (-2), xx, 0, 5).removeO().coeff(xx, n) for n in range(5)]
c_plus = [sp.expand((1 - xx) ** 2).coeff(xx, n) for n in range(3)]
check("A2 Punkte 0 und 1/2: Stromkern 1/(1 - z-quer w)^2 hat det 7/9; mit Zähler, (1 - z-quer w)^2, det -7/16. "
      "Grund: Taylor-Koeffizienten (1, 2, 3, ...) positiv gegen (1, -2, 1) [E]",
      d1 == sp.Rational(7, 9) and d2 == sp.Rational(-7, 16) and all(c > 0 for c in c_minus) and c_plus == [1, -2, 1])
# A3 Charakter von q
z, u = sp.symbols("z u")
qcoef = lambda zz: 1 / (zz ** 4 - 1)                    # q = qcoef(z) dz^2
rot = sp.simplify(qcoef(sp.I * u) * sp.diff(sp.I * u, u) ** 2 / qcoef(u))
inv = sp.simplify(qcoef(1 / u) * sp.diff(1 / u, u) ** 2 / qcoef(u))
check("A3 q = dz^2/(z^4 - 1): unter z -> i z und unter z -> 1/z geht q in -q über [E; Bestand: Prüfung H6h34, W1]",
      rot == -1 and inv == -1)
# A4 Markenmatrizen
circ = lambda c: sp.Matrix(4, 4, lambda i, j: c[(j - i) % 4])
Bm, Cm = circ([2, -1, 0, -1]), circ([0, 1, 2, 1])
Rm = circ([0, 1, 0, 0])
Um = sp.diag(1, -1, 1, -1)
specB = sorted(Bm.eigenvals(multiple=True))
specC = sorted(Cm.eigenvals(multiple=True))
check("A4 C = R^2 U B U und C^2 = U B^2 U für B = circ(2,-1,0,-1), C = circ(0,1,2,1); Spektren {0,2,2,4} und "
      "{-2,-2,0,4}: die explizite Form des 'Kanaltauschs' aus ALPHA.QUILLEN.EXACT.01 [E]",
      Cm == Rm ** 2 * Um * Bm * Um and Cm ** 2 == Um * Bm ** 2 * Um and specB == [0, 2, 2, 4]
      and specC == [-2, -2, 0, 4])
check("A4b U ist die Multiplikation mit dem Charakter -1 der Vierteldrehung (U R U = -R): derselbe Charakter wie "
      "der von q [E]", Um * Rm * Um == -Rm)
# A5 Lambda G
AINV = mp.mpf("137.0359992168407")
LG = 3 / (32 * mp.pi ** 3) * mp.e ** (-2 * AINV)
ZAHLEN["Lambda_G"] = mp.nstr(LG, 6)
check("A5 Lambda G = 3/(32 pi^3) e^(-2/alpha) = 2,83508e-122 (aus rho_Lambda/M^4 = (3/4 pi^2) e^(-2/alpha), "
      "M^2 = 1/(8 pi G)) [E gegeben die Ledgerformel]", abs(LG / mp.mpf("2.83508e-122") - 1) < 1e-5)
# A6 magnetische Sektoren: Felder (Q in Einheiten e, Trialität)
fields = {"u": (Fr(2, 3), 1), "d": (Fr(-1, 3), 1), "e": (Fr(-1), 0), "nu": (Fr(0), 0), "W": (Fr(1), 0)}
ok_sector = lambda m, c: all((m * Q + Fr(c * t, 3)).denominator == 1 for Q, t in fields.values())
table = {m: [c for c in range(3) if ok_sector(m, c)] for m in range(1, 7)}
ZAHLEN["magnetische_sektoren"] = {str(m): v for m, v in table.items()}
say(f"   Zulässige Farbflüsse c (mod 3) je elektromagnetischem Fluss m (e g = 2 pi m): {table}")
check("A6 Reiner elektromagnetischer Fluss (c = 0) ist nur für m = 0 mod 3 zulässig; m = 1, 2, 4 nur mit Farbfluss. "
      "Die Eichgruppe zeichnet m = 2 (alpha_m = 1/alpha) nicht aus [E]",
      table[1] == [1] and table[2] == [2] and table[3] == [0] and table[4] == [1])
# A7 Cauchy-Schwarz
w1, w2, s1, s2 = sp.symbols("w1 w2 s1 s2", positive=True)
f0, f2, f4 = w1 + w2, w1 / s1 + w2 / s2, w1 / s1 ** 2 + w2 / s2 ** 2
check("A7 Positives Laplace-Profil f(x) = Summe w_i e^(-s_i x): f_0 f_4 - f_2^2 = w1 w2 (1/s1 - 1/s2)^2 >= 0; "
      "mit einem negativen Gewicht (graduiert) kann die Differenz null oder negativ sein [E]",
      sp.simplify(f0 * f4 - f2 ** 2 - w1 * w2 * (1 / s1 - 1 / s2) ** 2) == 0
      and (f0 * f4 - f2 ** 2).subs({w1: 1, w2: -1, s1: 1, s2: 2}) < 0)
# A8 konforme Verträglichkeit: Eguchi-Hanson (A1) als Modell der Auflösung
r, th, ph, ps, a = sp.symbols("r theta phi psi a", positive=True)
W = 1 - a ** 4 / r ** 4
X = [r, th, ph, ps]
g = sp.zeros(4)
g[0, 0] = 1 / W
g[1, 1] = r ** 2 / 4
g[2, 2] = r ** 2 / 4 * sp.sin(th) ** 2 + r ** 2 / 4 * W * sp.cos(th) ** 2
g[3, 3] = r ** 2 / 4 * W
g[2, 3] = g[3, 2] = r ** 2 / 4 * W * sp.cos(th)
gi = g.inv()
Gam = [[[sum(gi[i, d] * (sp.diff(g[d, k], X[j]) + sp.diff(g[d, j], X[k]) - sp.diff(g[j, k], X[d])) for d in range(4)) / 2
         for k in range(4)] for j in range(4)] for i in range(4)]
pt = {r: sp.Rational(3, 2), th: sp.pi / 3, a: 1}
Gn = [[[Gam[i][j][k] for k in range(4)] for j in range(4)] for i in range(4)]


def riem(i, j, k, l):
    e = sp.diff(Gn[i][l][j], X[k]) - sp.diff(Gn[i][k][j], X[l]) + sum(
        Gn[i][k][m] * Gn[m][l][j] - Gn[i][l][m] * Gn[m][k][j] for m in range(4))
    return mp.mpf(sp.N(e.subs(pt), 30))


Rup = [[[[riem(i, j, k, l) for l in range(4)] for k in range(4)] for j in range(4)] for i in range(4)]
gN = [[mp.mpf(sp.N(g[i, j].subs(pt), 30)) for j in range(4)] for i in range(4)]
giN = [[mp.mpf(sp.N(gi[i, j].subs(pt), 30)) for j in range(4)] for i in range(4)]
ricci = max(abs(sum(Rup[k][i][k][j] for k in range(4))) for i in range(4) for j in range(4))
Rdn = [[[[sum(gN[i][m] * Rup[m][j][k][l] for m in range(4)) for l in range(4)] for k in range(4)] for j in range(4)]
       for i in range(4)]
Rupall = [[[[sum(giN[j][n] * giN[k][p] * giN[l][q] * Rup[i][n][p][q] for n in range(4) for p in range(4)
                 for q in range(4)) for l in range(4)] for k in range(4)] for j in range(4)] for i in range(4)]
kret = sum(Rdn[i][j][k][l] * Rupall[i][j][k][l]
           for i in range(4) for j in range(4) for k in range(4) for l in range(4))
kret_ref = 384 / mp.mpf(3) ** 12 * 2 ** 12
ZAHLEN.update(EH_ricci=mp.nstr(ricci, 3), EH_kretschmann=mp.nstr(kret, 12))
say(f"   Eguchi-Hanson bei r = 3/2, a = 1: max |Ricci| = {mp.nstr(ricci, 3)}, Riem^2 = {mp.nstr(kret, 12)} "
    f"(384 a^8/r^12 = {mp.nstr(kret_ref, 12)})")
check("A8 Die aufgelöste Geometrie ist Ricci-flach und nicht flach (Riem^2 = 384 a^8/r^12), also Weyl ungleich null "
      "und nicht konform flach. S^4 und AdS_4 (z164) sind konform flach: konform verträglich ist nur der "
      "Orbifoldpunkt a = 0, nicht die Auflösung [E für A1; A3 analog, L]",
      ricci < mp.mpf("1e-20") and abs(kret / kret_ref - 1) < mp.mpf("1e-20"))
befund("A_nachrechnung", "8 von 8", "Alle acht gerechneten Aussagen der externen Analyse stimmen. Zwei sind Bestand "
       "(Charakter von q; Kanaltausch der Markenmatrizen), die Operatorform C = R^2 U B U ist neu. Präzisierung zu "
       "(8): Am Orbifoldpunkt sind beide Geometrien verträglich, erst die Auflösung trennt sie.")

# ------------------------------------------------------------------------------------------------------------------
say("\nTeil B  Folgerechnungen")
say("   B1 Gravitation ohne Doppelpol")
Lam, G, Rsc, cc, Eul = sp.symbols("Lambda G R c E")
# eps F F mit F = R - (Lambda/3) e e: eps R R = E (Euler), eps R e e = 2 R sqrt(g), eps e e e e = 24 sqrt(g)
dens = Eul - 2 * (Lam / 3) * 2 * Rsc + (Lam / 3) ** 2 * 24
check("B1a eps F F mit F = R - (Lambda/3) e e ist Euler - (4 Lambda/3)(R - 2 Lambda): Einstein-Hilbert mit "
      "kosmologischem Term im richtigen Verhältnis, zweiter Ordnung [E]",
      sp.expand(dens - (Eul - sp.Rational(4, 3) * Lam * (Rsc - 2 * Lam))) == 0)
csol = sp.solve(sp.Eq(cc * (-sp.Rational(4, 3) * Lam), 1 / (16 * sp.pi * G)), cc)[0]
cval = 3 / (64 * mp.pi * LG)
ZAHLEN["log10_kopplung_MM"] = mp.nstr(mp.log10(cval), 6)
check("B1b Koeffizient c = -3/(64 pi G Lambda): Die einzige Kopplung ist die reine Zahl 1/(G Lambda), wie A^2/V "
      "bei Spin 1 (z164); mit der Leiter |c| = 10^119,7 [E]", sp.simplify(csol + 3 / (64 * sp.pi * G * Lam)) == 0
      and abs(mp.log10(cval) - mp.mpf("119.72")) < 0.01)
bb, FF, VV = sp.symbols("b F V")
Sbf = bb * FF - VV * bb ** 2
check("B1c Erste Ordnung plus Vertex b^2: Ausintegrieren gibt F^2/(4V), quadratisch in der Krümmung, also zweiter "
      "Ordnung in der Metrik (F enthält R linear): kein Doppelpol. Voraussetzung ist Lambda ungleich null "
      "(nicht entarteter Unendlichkeitstwistor, z164 D3) [E für die Algebra, C unter P_MM]",
      sp.simplify(Sbf.subs(bb, sp.solve(sp.diff(Sbf, bb), bb)[0]) - FF ** 2 / (4 * VV)) == 0)
say("   B2 Lambda-Stufe als relative Nahtantwort")
c3 = 1 / (8 * mp.pi)
dtop = 48 * c3 ** 4
qf = lambda x: dtop * mp.e ** (-2 * x)
devf = lambda x: qf(x) * (1 - qf(x)) ** (-mp.mpf(5) / 4)      # phi_seam(x) - phi_base, ohne Auslöschung
rel_dual = devf(AINV)
rho = 3 / (4 * mp.pi ** 2) * mp.e ** (-2 * AINV)
rel_el = devf(1 / AINV) * 6 * mp.pi
ZAHLEN.update(rel_antwort_elektrisch=mp.nstr(rel_el, 6), q_elektrisch=mp.nstr(qf(1 / AINV), 6))
check("B2 rho_Lambda c3^2/M^4 = phi_seam(1/alpha) - phi_base bis auf relativ 1e-120: Die Lambda-Stufe ist die "
      "Abweichung der Nahtantwort von ihrem Basiswert, am dualen Punkt. Am elektrischen Punkt ist dieselbe "
      "Abweichung 0,22 % des Basiswerts [E für die Identität, L für die Lesart, eine Wahl]",
      abs(rel_dual / (rho * c3 ** 2) - 1) < mp.mpf("1e-30") and abs(rel_el - mp.mpf("0.00223")) < 1e-4)
say(f"   (phi_seam - phi_base)/phi_base bei alpha: {mp.nstr(rel_el, 6)}; q(alpha) = {mp.nstr(qf(1 / AINV), 6)}")
m_ = sp.Symbol("m", positive=True)
al = sp.Symbol("alpha", positive=True)
S2d = sp.pi * m_ ** 2 / (2 * 4 * sp.pi * al)            # (1/2e^2) int F^2 auf der Einheitssphäre, Fluss 2 pi m
check("B2b Zwei Lesarten desselben Exponenten: 2 alpha_m mit m = 2, oder Maxwell-Wirkung des Flusses m = 4 auf der "
      "Einheitssphäre (m^2/(8 alpha) = 2/alpha). Keine ist als reiner elektromagnetischer Fluss zulässig (A6): die "
      "duale Lesart ist nicht eindeutig [E/L]", sp.simplify(S2d.subs(m_, 4) - 2 / al) == 0
      and 0 not in table[2] and 0 not in table[4])
say("   B3 Ungleichung gegen die Ledgerzahlen")
ratio = rho * 16 * mp.pi                                  # rho_Lambda/M^4 mal 1/g^2 an der Planckmasse (1/g^2 ~ 4, z164)
ZAHLEN["log10_f0f4_durch_f2quadrat"] = mp.nstr(mp.log10(rho * 4), 6)
check("B3 [N] Als Momente eines positiven Profils müssten rho_Lambda (f_4), M^2 (f_2) und 1/g^2 (f_0) "
      "f_0 f_4 >= f_2^2 bis auf Zählfaktoren der Ordnung eins erfüllen. Mit den Ledgerzahlen ist das Verhältnis "
      "10^-119,5: Die kosmologische Konstante kann kein absolutes Moment sein", mp.log10(rho * 4) < -119)
say("   B4 Native Graduierung: Spuren der Uhr auf der adjungierten Darstellung von E8")


def e8_roots():
    out = []
    for i, j in itertools.combinations(range(8), 2):
        for si in (1, -1):
            for sj in (1, -1):
                w = [Fr(0)] * 8
                w[i], w[j] = Fr(si), Fr(sj)
                out.append(tuple(w))
    for eps in itertools.product((1, -1), repeat=8):
        if eps.count(-1) % 2 == 0:
            out.append(tuple(Fr(e, 2) for e in eps))
    return sorted(out)


E8 = e8_roots()
cls = lambda w: int(2 * sum(w[:5])) % 4
dims = [sum(1 for w in E8 if cls(w) == k) + (8 if k == 0 else 0) for k in range(4)]
trR = [sum(sp.I ** (n * k) * dims[k] for k in range(4)) for n in range(4)]
p2c = [sum(w[0] ** 2 for w in E8 if cls(w) == k) for k in range(4)]      # Koeffizient von x1^2 je Klasse
q2c = [sum(w[7] ** 2 for w in E8 if cls(w) == k) for k in range(4)]
trRX = [(sum(sp.I ** (n * k) * p2c[k] for k in range(4)), sum(sp.I ** (n * k) * q2c[k] for k in range(4)))
        for n in range(4)]
ZAHLEN.update(uhr_dims=dims, Sp_R_n=[str(x) for x in trR], Sp_R_n_X2=[[str(sp.nsimplify(a_)), str(sp.nsimplify(b_))]
                                                                    for a_, b_ in trRX])
say(f"   Klassen {dims}; Sp(R^n) = {[str(x) for x in trR]}; Sp(R^n X^2) = (p-Teil, q-Teil): "
    f"{[(str(a_), str(b_)) for a_, b_ in trRX]}")
check("B4a Sp_248(R) = Sp_248(R^3) = 0, Sp_248(R^2) = -8 = -Rang; quadratische Spuren je Klasse (16, 16, 12, 16) "
      "und (8, 16, 20, 16), Summe 60 = 2 h [E]", trR == [248, 0, -8, 0] and p2c == [16, 16, 12, 16]
      and q2c == [8, 16, 20, 16])
check("B4b Im Uhr-verdrehten Sektor (n = 1) verschwindet die Zählspur (Volumenterm), die quadratische Spur nicht: "
      "Sp(R X^2) = 4 p2 - 12 q2. Mit der Uhr als Graduierung hat die Determinante am symmetrischen Punkt keinen "
      "Volumenterm, aber einen Eichterm [E; Bedeutung L]", trRX[1] == (4, -12) and trR[1] == 0)
say("   B5 Die beiden Markenmatrizen des alpha-Blocks als graduiertes Paar")
xs_ = sp.Symbol("x")
UBU = Um * Bm * Um
detC = sp.factor((sp.eye(4) - xs_ * Cm).det())
detB = sp.factor((sp.eye(4) - xs_ * UBU).det())
R2 = Rm ** 2
check("B5a R^2 vertauscht mit U B U, und C = R^2 (U B U): C ist dieselbe Matrix mit dem Vorzeichen R^2 je Kanal. "
      "In H6 ist R^2 = (-1)^F (Spin-Z4-Struktur der Uhr) [E; Identifikation mit (-1)^F ist Bestand von H6, L]",
      R2 * UBU == UBU * R2 and Cm == R2 * UBU and R2 ** 2 == sp.eye(4))
check("B5b det(1 - x C) = (1 - 4x)(1 + 2x)^2 und det(1 - x U B U) = (1 - 4x)(1 - 2x)^2: Die Kanäle mit R^2 = +1 "
      "(Eigenwerte 4 und 0) sind gemeinsam, das Dublett mit R^2 = -1 wechselt das Vorzeichen. Das Verhältnis "
      "((1 + 2x)/(1 - 2x))^2 ist eine relative Determinante, die nur den ungeraden Kanal sieht [E]",
      sp.expand(detC - (1 - 4 * xs_) * (1 + 2 * xs_) ** 2) == 0
      and sp.expand(detB - (1 - 4 * xs_) * (1 - 2 * xs_) ** 2) == 0)
befund("B_folge", {"MM_kopplung": "1/(G Lambda)", "uhr_volumen": 0},
       "(1) Der Doppelpol entsteht nicht, wenn die Gravitation wie Spin 1 als erste Ordnung plus Vertex b^2 mit der "
       "so(3,2)-Krümmung eingeht (MacDowell-Mansouri, nativ im Bestand); die einzige Kopplung ist 1/(G Lambda). "
       "(2) Die Lambda-Stufe ist eine relative Antwort, kein absolutes Moment; die positive Ungleichung wäre um 119 "
       "Größenordnungen verletzt. (3) Die native Graduierung ist die Uhr, nicht Supersymmetrie: Sp(R) = 0 auf E8. "
       "(4) Im alpha-Block des Ledgers steht die graduierte Struktur schon: Die Greensche Matrix der Marken ist "
       "die mit R^2 graduierte Sprungmatrix.")

# ------------------------------------------------------------------------------------------------------------------
say("\nTeil C  Bilanz und Urteil")
bilanz = [
    ("E", "Nachrechnung", "Alle acht Aussagen der externen Analyse stimmen (A1 bis A8)."),
    ("E negativ", "W3", "Der Vertexkern ist kein positiver Zweipunktkern. Positivität gehört an den Propagator "
     "(Serre-Paarung mit Spiegelung, z123), nicht an den Vertex. Der erste Test von W3 war falsch gestellt."),
    ("E/C", "W2", "Kein Doppelpol bei erster Ordnung plus b^2 mit so(3,2)-Krümmung (P_MM); Kopplung 1/(G Lambda)."),
    ("E", "Hintergrund", "ALE-Auflösung und nativer Unendlichkeitstwistor sind nur am Orbifoldpunkt verträglich."),
    ("E/L", "W5", "Lambda-Stufe = relative Nahtantwort am dualen Punkt; die Flusszahl ist nicht ausgezeichnet."),
    ("E/L", "Graduierung", "Sp(R) = 0 auf E8: Uhr als native Graduierung, kein Volumenterm am symmetrischen Punkt."),
    ("E", "alpha-Block", "C = R^2 (U B U): Die beiden Markenmatrizen sind ein mit R^2 graduiertes Paar; ihre "
     "relative Determinante sieht nur das ungerade Dublett."),
    ("O", "Objekt", "Der Operator, sein Zustandsraum, das Maß und die Metrik der Determinantenlinie fehlen."),
]
for st, was, text in bilanz:
    say(f"   [{st}] {was}: {text}")
urteil = ("Die externe Analyse ist in allen gerechneten Punkten richtig und korrigiert zwei Stellen des Stands z163: "
          "Der erste Positivitätstest von W3 war am falschen Objekt gestellt, und der Hintergrund des "
          "Celest-Vertrags (aufgelöste ALE-Geometrie) ist mit dem nativen Unendlichkeitstwistor nur am "
          "Orbifoldpunkt verträglich. Ihr Kandidat, eine relative graduierte Determinante, wird durch drei eigene "
          "Rechnungen gestützt und eingeengt. (1) Relativ: Die Lambda-Stufe des Ledgers ist exakt die Abweichung "
          "der Nahtantwort von ihrem Basiswert am dualen Punkt, und als absolutes Moment eines positiven Profils "
          "wäre sie um 119 Größenordnungen ausgeschlossen. (2) Graduiert: Die native Graduierung ist die Uhr; "
          "ihre Spur auf E8 verschwindet, die quadratische nicht. (3) Einstein statt konform: Der Doppelpol "
          "entsteht nicht, wenn die Gravitation wie die Eichfelder als erste Ordnung plus quadratischem Vertex "
          "mit der so(3,2)-Krümmung eingeht; das verlangt Lambda ungleich null, was der native "
          "Unendlichkeitstwistor liefert, und lässt genau die Zahl G Lambda frei. Nicht gelöst ist das Objekt "
          "selbst: der Operator, sein Maß und die Metrik der Determinantenlinie (Bismut-Freed), aus der die "
          "Kopplung folgen müsste, sowie die Auswahl des zweiten Sektors. Eine volle Lösung liegt nicht vor. "
          "Nichts ist befördert.")
say("\nUrteil\n   " + urteil)
print(f"\n   (Laufzeit {time.time() - T0:.1f} s, nicht Teil der Ausgabe)", flush=True)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z165_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
(HERE / "ausgabe_z165.txt").write_text("\n".join(OUT) + "\n")
