"""Z151: Ist die integrable Nahtkette die Wilson-Linien-Konfiguration einer 4D-Chern-Simons-Theorie aus dem Twistorbild?

Leitlinie 'Rekonstruktion statt Auswahl'; H6 und die gelösten Blöcke (Flavour, Massen, alpha, 16, E8) bleiben unberührt,
nichts wird befördert. Statuskürzel: E exakt, C bedingt (Prämisse genannt), L Lesart oder Literaturaussage (nicht
nachgerechnet), N numerisch, O offen.

Frage: In H6 ist die Komposition von Nähten unter P_uhr und P_ord die XXX-Familie R(u) = (1 + i u P)/(1 + i u) bei u = 1,
die alternierende Kette das Destri-de-Vega-Lichtkegelgitter (z121, z144 B6/B7). In der Literatur entstehen R-Matrizen
mit Spektralparameter aus der 4D-Chern-Simons-Theorie (Costello; Costello-Witten-Yamazaki), und diese aus holomorpher
Chern-Simons-Theorie auf dem Twistorraum mit meromorpher Volumenform (Costello; Bittleston-Skinner). Ist die Nahtkette
eine solche Wilson-Linien-Konfiguration, folgt sie aus dem nativen Twistorbild (Naht = Twistorgerade, Marken z^4 = 1,
I_pi, I_omega), und welche 4D-Theorie wäre das?

Teil L  Literaturbefunde mit Quelle (Status L).
Teil A  R-Matrix: Yang-Baxter und klassische r-Matrix der rationalen 4D-CS-Klasse; das H6-Gatter in dieser Sprache;
        was u = 1 heißt; besondere Punkte; Kreuzung und Rahmung; rational gegen trigonometrisch.
Teil B  Die alternierende Kette als Gitter zweier Scharen von Wilson-Linien; die Treppe als eine Linie durch alle.
Teil C  Meromorphe Formen auf der Naht: Klassifikation, Uhrcharaktere, O(-4), Residuen, Reduktion, Randbedingungen.
Teil D  Anomalie-Bedingung (quartische Casimir-Identität) für die Gruppen des Bestands, mit und ohne Materie.
Teil E  Bewertung der Strategie, erste entscheidende Rechnung, Urteil.
"""
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

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
LITERATUR = []
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


def lit(kennung, quelle, gelesen, aussage):
    LITERATUR.append({"kennung": kennung, "quelle": quelle, "gelesen": gelesen, "aussage": aussage})
    say(f"   [L] {kennung} ({quelle}; gelesen: {gelesen}): {aussage}")


# ------------------------------------------------------------------------------------------------------------------
say("Z151: Nahtkette als Wilson-Linien einer 4D-Chern-Simons-Theorie aus dem Twistorbild?")
say("=" * 100)
say("\nTeil L: Literaturbefunde (Status L; nur Gelesenes; arXiv-Nummern)")

say("  (a) R-Matrix aus der 4D-Chern-Simons-Theorie")
lit("L1", "Costello, arXiv:1303.2632", "Zusammenfassung",
    "Eine verdrillte, deformierte N = 1-Eichtheorie wird vom Yangian regiert; ein Netz aus n + m Wilson-Operatoren "
    "gibt die Zustandssumme eines integrablen Spinkettenmodells auf dem n x m-Gitter.")
lit("L2", "Costello-Witten-Yamazaki I, arXiv:1709.09993", "Volltext (ar5iv), Abschnitte 3 bis 7",
    "Wirkung S = (1/2 pi) Integral omega ^ CS(A) auf Sigma x C, Gewicht exp(i S/hbar) (Gl. 3.3, 3.5, 3.13). Drei Fälle "
    "(Gl. 3.14): rational C = C, omega = dz, Doppelpol bei unendlich; trigonometrisch C = C*, omega = dz/z, zwei "
    "einfache Pole; elliptisch. omega hat keine Nullstellen; am Doppelpol gilt A = 0. hbar ist nicht quantisiert "
    "(Abschnitt 3.1). Für C = C ist die Theorie invariant unter gemeinsamer Skalierung von z und hbar: R hängt nur "
    "von (z_1 - z_2)/hbar ab (Abschnitt 3.4), nicht von der Lage in Sigma. Führende Ordnung (Gl. 4.14): R = 1 + hbar "
    "c/(z_1 - z_2) + O(hbar^2) mit c = sum t_a (x) t_a (Gl. 4.12). Für V (x) V von GL_N: R = A + P B mit B/A = hbar/z "
    "(Gl. 3.24); für V (x) V*: W(z) = -hbar/(z + hbar N/2) (Gl. 3.29). Rahmungsanomalie (Abschnitt 3.6, 7.4): Dreht "
    "sich eine Linie um den Winkel phi, verschiebt sich ihr Spektralparameter um hbar h^v phi/(2 pi); eine halbe "
    "Drehung gibt hbar N/2 (Kreuzung, Gl. 3.31). Die höheren Ordnungen sind bis auf Umparametrisierung eindeutig "
    "(Drinfeld, Abschnitt 5.5); die Kategorie der Linien sind Darstellungen des Yangians.")
lit("L3", "Costello-Witten-Yamazaki II, arXiv:1802.01579", "Zusammenfassung",
    "Die Quantengruppen der rationalen, trigonometrischen und elliptischen Klasse folgen zu allen Ordnungen in hbar "
    "über die RTT-Darstellung; am direktesten für gl_N, erweiterbar auf alle einfachen Algebren außer e_8.")
lit("L4", "Costello-Yamazaki III, arXiv:1908.02289", "Volltext (ar5iv), erste Hälfte",
    "Ordnungsdefekte: freie chirale und antichirale Fermionen in einer reellen Darstellung an zwei Punkten z_0, z_1 "
    "von C = C (omega = dz). Ausintegrieren des Eichfelds gibt eine Strom-Strom-Kopplung proportional 1/(z_0 - z_1) "
    "(Gl. 4.16): Gross-Neveu (so_N) bzw. masseloses Thirring-Modell (sl_N). Anomaliefrei nur mit gleichem "
    "Dynkin-Index beider Seiten. Unordnungsdefekte: Pole des Eichfelds sind nur an Nullstellen von omega erlaubt.")
lit("L5", "Ashwinkumar-Sakamoto-Yamazaki, arXiv:2309.14412", "Volltext (ar5iv), erste Hälfte",
    "Chirale und antichirale Flächendefekte werden entlang ihrer Lichtkegelrichtung in Wilson-Linien zerlegt: ein "
    "Lichtkegelgitter mit chiralen Linien auf geraden und antichiralen auf ungeraden Plätzen (Abschnitt 4.3), "
    "Spektralparameter z -+ nu (Gl. 4.20 bis 4.24). Für Thirring ist die 2D-Kopplung hbar/(z_+ - z_-) bis auf "
    "Normierung (Gl. 5.16). Der Kontinuumslimes braucht Gitterabstand -> 0 mit Skalierung von hbar (Gl. 4.28, "
    "Fußnoten 10, 11) und Korrelationslänge >> Gitterabstand.")

say("  (b) Meromorphe Volumenformen auf dem Twistorraum")
lit("L6", "Bittleston-Skinner, arXiv:2011.04638", "Volltext (ar5iv)",
    "S = (1/2 pi i) Integral Omega ^ hCS(A-quer), Omega = Phi D^3Z mit Phi meromorpher Schnitt von O(-4), nur von pi "
    "abhängig (Abschnitt 3.1). Zwei Doppelpole (Costellos Wahl, Gl. 2.4): WZW_4, Bewegungsgleichung = Yangs Gleichung "
    "= ASD-Yang-Mills; Randbedingung: A-quer durch (Z.A)(Z.B) teilbar. Ein Pol vierter Ordnung: "
    "Leznov-Mukhtarov-Parkes-Wirkung (Gl. 3.10). Vier einfache Pole (Abschnitt 3.3): 'trigonometrische' Wirkung, "
    "Randwerte in Lagrangeschen Unteralgebren (Manin-Tripel); die Wirkungen entsprechen wieder WZW_4. Ist Phi "
    "nullstellenfrei, sind die Bewegungsgleichungen ASD-Yang-Mills; mit Nullstellen entstehen gekoppelte "
    "Sigma-Modelle mit 4D-Lax-Verbindung (Abschnitt 3.4). Reduktion längs zweier Translationen gibt 4D-CS mit omega "
    "= Phi <pi gamma><pi gamma-Dach><d pi pi> (Abschnitt 2.5): Die zwei Nullstellen kommen aus der Reduktion; in "
    "euklidischer Signatur sind sie antipodal. Für Costellos Phi: omega = (zeta - zeta_0)(zeta - zeta_0-Dach) d zeta/"
    "((zeta_0 - zeta_0-Dach) zeta^2), Hauptchirales Modell mit WZ-Term k = <alpha gamma><alpha-Dach gamma-Dach> + "
    "<alpha gamma-Dach><alpha-Dach gamma>. Realität (Abschnitt 4): euklidisch oder ultrahyperbolisch (2,2); in "
    "lorentzscher Signatur lassen reelle Eichgruppen keine ASD-Verbindungen zu. Die Arbeit ist klassisch.")
lit("L7", "Penna, arXiv:2011.05831", "Zusammenfassung",
    "Eine Chern-Simons-Wirkung auf dem Twistorraum liefert Lagrangefunktionen und Lax-Operatoren integrabler "
    "Sigma-Modelle aus dimensional reduzierter Gravitation und Supergravitation.")

say("  (c) Anomalien der Twistortheorien")
lit("L8", "Costello, arXiv:2111.08879", "Volltext (ar5iv)",
    "WZW_4 = holomorphe Chern-Simons-Theorie mit Doppelpolen der Volumenform bei z = 0 und z = unendlich; selbstduales "
    "Yang-Mills = holomorphe BF-Theorie. Die Einschleifen-Anomalie (Kastendiagramm, Gl. 2.1.1) ist proportional zur "
    "adjungierten Spur von (X^4). Ein Green-Schwarz-Feld hebt sie genau dann, wenn Tr_adj(X^4) = lambda^2 tr(X^2)^2 "
    "(Gl. 1.2.2, 2.2.8); das gilt für sl_2, sl_3, so_8 und alle exzeptionellen Algebren. Für selbstduales Yang-Mills "
    "ist das Zusatzfeld ein Axion rho mit Term (Laplace rho)^2 und Kopplung an CS(A) (Gl. 1.2.4); mit Materie in "
    "einer reellen Darstellung R wird tr_R(X^4) abgezogen (Fußnote 3), etwa sl(N_c) mit N_f = N_c oder so(N_c) mit "
    "N_f = N_c - 8 (Fußnote 4). Für WZW_4 ist die Quantisierung nur für SO(8) bewiesen (Satz 2.2.2); das Zusatzfeld "
    "ist ein Kähler-Potential mit (Laplace rho)^2 (Gl. 1.1.3, 1.1.4). Die Theorien sind nicht unitär, ihre "
    "Korrelatoren rational. Volles Yang-Mills entsteht erst durch die Deformation -(g^2/4) B^2 (Abschnitt 1.4). Die "
    "dort genannten Werte von lambda^2 (8 für sl_2, 9 für sl_3, 3/2 für so_8) stimmen für sl_2 und sl_3 mit Teil D "
    "überein; der Wert für so_8 steht in einer Spurnormierung, die aus dem gelesenen Text nicht zu klären war.")
lit("L9", "Costello-Paquette, arXiv:2204.05301", "Volltext (ar5iv)",
    "Die Einschleifen-Kollinearsingularitäten der selbstdualen Eichtheorie bilden nur dann eine assoziative chirale "
    "Algebra, wenn die Twistoranomalie gehoben ist. Bedingung (Gl. 15, 20): Tr_adj(X^4) = lambda_g^2 tr(X^2)^2 mit "
    "lambda_g^2 = 10 (h^v)^2/(dim g + 2), tr die minimale Spur. Gruppen: SU(2), SU(3), SO(8), exzeptionelle. "
    "Alternative mit Materie: SU(2) mit N_f = 8 (acht Fundamentale und acht Antifundamentale).")
lit("L10", "Costello-Paquette, arXiv:2201.02595", "Zusammenfassung",
    "Formfaktoren 4D-selbstdualer Theorien sind Korrelatoren einer chiralen Algebra auf der Twistorgeraden; "
    "Parke-Taylor und CSW folgen daraus.")

say("  (d) Was volles Yang-Mills oder Gravitation erreicht")
lit("L11", "Costello-Paquette-Sharma, arXiv:2208.14233", "Zusammenfassung",
    "WZW_4 mit SO(8), gekoppelt an skalarflache Kähler-Gravitation auf dem Burns-Raum, ist dual zu einer chiralen "
    "Algebra aus geeichten beta-gamma-Systemen. Das ist keine Einstein-Gravitation.")
lit("L12", "Bittleston-Sharma-Skinner, arXiv:2208.12701", "Zusammenfassung",
    "Die Poisson-BF-Theorie (selbstduale Gravitation) ist auf dem Twistorraum anomal; ein weiteres Twistorfeld (auf der "
    "Raumzeit ein Axion) hebt die Anomalie, danach verschwinden alle Amplituden ab vier Teilchen.")
lit("L13", "Mason, arXiv:hep-th/0507269", "Zusammenfassung",
    "Twistorwirkungen für volles Yang-Mills und für konforme Gravitation, mit nicht integrablen fast komplexen "
    "Strukturen; Grundlage der Twistorstring-Baumamplituden.")
lit("L14", "Berkovits-Witten, arXiv:hep-th/0406051", "Zusammenfassung",
    "Im Twistorstring koppelt konforme Supergravitation mit derselben Stärke wie Yang-Mills und trägt zu den "
    "Schleifen bei; reines Super-Yang-Mills entsteht so nicht.")
lit("L15", "Adamo-Mason, arXiv:1307.5043", "Zusammenfassung",
    "Einstein-Gravitation mit kosmologischer Konstante ist in die konforme Gravitation eingebettet; die "
    "MHV-Amplituden folgen aus der Twistorwirkung der konformen Gravitation.")

# ------------------------------------------------------------------------------------------------------------------
say("\nTeil A: Das H6-Gatter als rationale 4D-CS-R-Matrix")


def perm_sym(n):
    """Tauschoperator P auf C^n (x) C^n als sympy-Matrix."""
    m = sp.zeros(n * n, n * n)
    for i in range(n):
        for j in range(n):
            m[i * n + j, j * n + i] = 1
    return m


def drei(n):
    """P_12, P_23, P_13 auf (C^n)^(x3)."""
    p, e = perm_sym(n), sp.eye(n)
    p12, p23 = sp.kronecker_product(p, e), sp.kronecker_product(e, p)
    return p12, p23, p23 * p12 * p23


def nullmat(m):
    return all(sp.expand(x) == 0 for x in m)


z1, z2, z3, hb, u, v, x, th = sp.symbols("z1 z2 z3 hbar u v x theta")
I = sp.I
P2 = perm_sym(2)
E4, E8m = sp.eye(4), sp.eye(8)
P12, P23, P13 = drei(2)

# A1: Yang-Baxter in Differenzform für R(z) = z + hbar P (gleich 1 + hbar P/z bis auf den Faktor z)
Rt = lambda P, z: z * sp.eye(P.shape[0]) + hb * P
lhs = Rt(P12, z1 - z2) * Rt(P13, z1 - z3) * Rt(P23, z2 - z3)
rhs = Rt(P23, z2 - z3) * Rt(P13, z1 - z3) * Rt(P12, z1 - z2)
ok_n2 = nullmat(lhs - rhs)
Q12, Q23, Q13 = drei(3)
ok_n3 = True
for pt in ((sp.Rational(3, 7), sp.Rational(-5, 4), sp.Rational(2, 9), sp.Rational(7, 3)),
           (sp.Rational(-11, 5), sp.Rational(1, 6), sp.Rational(13, 8), sp.Rational(-2, 7))):
    a1, a2, a3, hh = pt
    r = lambda P, z: z * sp.eye(27) + hh * P
    ok_n3 &= (r(Q12, a1 - a2) * r(Q13, a1 - a3) * r(Q23, a2 - a3)
              - r(Q23, a2 - a3) * r(Q13, a1 - a3) * r(Q12, a1 - a2)) == sp.zeros(27, 27)
check("[E] R(z) = 1 + hbar P/z (CWY I, Gl. 3.24: B/A = hbar/z) erfüllt die Yang-Baxter-Gleichung in Differenzform "
      "R_12(z_1 - z_2) R_13(z_1 - z_3) R_23(z_2 - z_3) = R_23 R_13 R_12: symbolisch für N = 2 (alle z_i und hbar frei), "
      "für N = 3 an zwei rationalen Punkten exakt", ok_n2 and ok_n3)

# A2: klassische r-Matrix r = c/z, c = sum t_a (x) t_a
a, b = sp.symbols("a b")
com = lambda A, B: A * B - B * A
cybe = ((com(P12 / a, P13 / (a + b)) + com(P12 / a, P23 / b) + com(P13 / (a + b), P23 / b)) * a * b * (a + b)).applyfunc(
    sp.cancel)
cas_ok = True
for n in (2, 3):
    c = sp.zeros(n * n, n * n)
    for i in range(n):
        for j in range(n):
            eij, eji = sp.zeros(n, n), sp.zeros(n, n)
            eij[i, j], eji[j, i] = 1, 1
            c += sp.kronecker_product(eij, eji)                 # Basis E_ij, dual E_ji bezüglich der Spur
    cas_ok &= c == perm_sym(n)
    c0 = c - sp.eye(n * n) / n                                    # sl_n-Anteil
    ptr = sp.zeros(n, n)
    for i in range(n):
        for j in range(n):
            ptr[i, j] = sum(c0[i * n + k, j * n + k] for k in range(n))
    cas_ok &= ptr == sp.zeros(n, n)
check("[E] Führende Ordnung (CWY I, Gl. 4.12, 4.14): Die klassische r-Matrix r(z) = c/z mit c = sum t_a (x) t_a erfüllt "
      "die klassische Yang-Baxter-Gleichung. Für gl_N mit der Spur der Fundamentalen ist c = P, für sl_N c = P - 1/N "
      "(N = 2, 3). R = 1 + hbar r ist für die Fundamentale schon eine exakte Lösung; die volle Quanten-R-Matrix weicht "
      "davon nur um einen skalaren Faktor und eine Umparametrisierung von hbar ab (CWY I, 5.5) [L]",
      nullmat(cybe) and cas_ok)

# A3: das H6-Gatter und die zwei Wörterbücher
G = lambda P, uu: (sp.eye(P.shape[0]) + I * uu * P) / (1 + I * uu)            # z121, Zopfform
Rn = lambda P, xx: (sp.eye(P.shape[0]) + xx * P) / (1 + xx)                    # normierte Linien-R-Matrix, x = hbar/Dz
Bz = lambda P, uu: sp.eye(P.shape[0]) + I * uu * P
zopf = Bz(P12, u) * Bz(P23, u + v) * Bz(P12, v) - Bz(P23, v) * Bz(P12, u + v) * Bz(P23, u)
platz = (G(P2, u) - P2 * Rn(P2, hb / (I * u * hb))).applyfunc(sp.simplify)       # Naht = Platz: Dz = i u hbar
linie = (G(P2, u) - Rn(P2, hb / (hb / (I * u)))).applyfunc(sp.simplify)          # Naht = Linie: Dz = hbar/(i u)
umschr = (G(P12, u) * G(P23, u + v) * G(P12, v)
          - P13 * Rn(P23, 1 / (I * u)) * Rn(P13, 1 / (I * (u + v))) * Rn(P12, 1 / (I * v))).applyfunc(sp.simplify)
inst = Bz(P12, 1) * Bz(P23, 2) * Bz(P12, 1) - Bz(P23, 1) * Bz(P12, 2) * Bz(P23, 1)
lin_yb = lambda u12, u13, u23: (Bz(P12, u12) * Bz(P13, u13) * Bz(P23, u23) - Bz(P23, u23) * Bz(P13, u13) * Bz(P12, u12))
harm = nullmat(lin_yb(1, sp.Rational(1, 2), 1)) and not nullmat(lin_yb(1, 2, 1))
check("[E] Das H6-Gatter R(u) = (1 + i u P)/(1 + i u) (z121) ist die rationale 4D-CS-R-Matrix in zwei Lesarten: 'Naht = "
      "Platz' (Gatter = P mal Linien-R-Matrix) mit Dz = i u hbar, und 'Naht = Linie' mit Dz = hbar/(i u). Die "
      "H6-Beziehungen sind in u additiv, und das passt nur zur ersten: Die Zopfgleichung R_12(u) R_23(u+v) R_12(v) = "
      "R_23(v) R_12(u+v) R_23(u) ist exakt P_13 mal der Liniengleichung mit Abständen (u, u+v, v) i hbar; die "
      "Faktorisierung dreier Nähte aus z144 B6 mit (1, 2, 1) ist der Fall z_1 - z_2 = z_2 - z_3 = i hbar. Als Linien "
      "gelesen wäre 1/u additiv: Drei Linien mit u_12 = u_23 = 1 verlangen u_13 = 1/2, mit u_13 = 2 scheitert die "
      "Gleichung",
      nullmat(zopf) and nullmat(platz) and nullmat(linie) and nullmat(umschr) and nullmat(inst) and harm)

# A4: besondere Punkte der Familie R(x) = (1 + x P)/(1 + x), x = hbar/Dz
ar, br = sp.symbols("a_r b_r", real=True)
xc = ar + I * br
unit = ((E4 + xc * P2) * (E4 + sp.conjugate(xc) * P2) - sp.expand((1 + xc) * (1 + sp.conjugate(xc))) * E4
        - 2 * ar * (P2 - E4)).applyfunc(sp.expand)
inv_ok = nullmat((E4 + x * P2) * (E4 - x * P2) - (1 - x ** 2) * E4)
dual_ok = nullmat(P2 * (E4 + x * P2) - x * (E4 + P2 / x))
wurzel_ok = nullmat((Rn(P2, I) ** 2 - P2).applyfunc(sp.simplify)) and nullmat((Rn(P2, -I) ** 2 - P2).applyfunc(sp.simplify))
proj_p, proj_m = (E4 + P2) / 2, (E4 - P2) / 2
proj_ok = proj_p ** 2 == proj_p and proj_m ** 2 == proj_m and proj_p.rank() == 3 and proj_m.rank() == 1
w_lsg = set(sp.solve(sp.Eq(((1 - x) / (1 + x)) ** 4, 1), x))
fix = {"x -> -x": sorted(sp.solve(sp.Eq(x, -x), x), key=str),
       "x -> 1/x": sorted(sp.solve(sp.Eq(x, 1 / x), x), key=lambda t: sp.re(t)),
       "x -> -1/x": sorted(sp.solve(sp.Eq(x, -1 / x), x), key=lambda t: sp.im(t))}
okta = set(fix["x -> 1/x"]) | set(fix["x -> -1/x"])
check("[E] Besondere Punkte von R(x) = (1 + x P)/(1 + x), x = hbar/Dz: R(x) R(-x) = 1; P R(x) = R(1/x); unitär genau für "
      "x rein imaginär (R R^+ - 1 ist proportional zu 2 Re x (P - 1)). Die drei Involutionen x -> -x, 1/x, -1/x haben die "
      "Fixpunkte {0, unendlich}, {1, -1}, {i, -i}: bei 0 und unendlich Identität und Tausch, bei +-1 die Projektoren "
      "(1 + P)/2 (Rang 3) und (1 - P)/2 (Rang 1, das Singulett von P_ident aus z121), bei +-i die unitären Gatter mit "
      "R^2 = P. Das Verhältnis der Eigenwerte (1 - x)/(1 + x) ist ein Uhrwert genau für x in {0, i, -i, unendlich}",
      nullmat(unit) and inv_ok and dual_ok and wurzel_ok and proj_ok and w_lsg == {0, -I, I}
      and fix["x -> -x"] == [0] and fix["x -> 1/x"] == [-1, 1] and fix["x -> -1/x"] == [-I, I]
      and okta == {1, -1, I, -I})
ZAHLEN["u = 1 in 4D-CS-Sprache"] = {"Naht = Platz": "Dz = i hbar, x = hbar/Dz = -i",
                                    "Naht = Linie": "Dz = -i hbar, x = +i",
                                    "Unterschied": "Vorzeichen von Dz = Zeitrichtung (R(-1) = R(1)^-1)"}
befund("u = 1 in der Sprache der 4D-CS", "Dz = +-i hbar",
       "P_uhr (z144: die Phase ist ein Uhrwert) heißt x = hbar/Dz in {0, unendlich, i, -i}; ohne Identität und Tausch "
       "bleibt Dz = +-i hbar. Das ist der Fixpunkt von x -> -1/x: Dort sind die Lesarten 'Naht = Platz' und 'Naht = "
       "Linie' Zeitumkehrungen voneinander. Die rationale Theorie hat die Skalensymmetrie (z, hbar) -> (l z, l hbar) "
       "(CWY I, 3.4) und damit je Paar von Scharen genau eine reine Zahl Dz/hbar; sie zeichnet keinen Wert aus. hbar ist "
       "also keine neue freie Zahl (es ist u), und P_uhr bleibt die Prämisse, die sie auf i setzt.")

# A5: Kreuzung und Rahmung
eps = sp.Matrix([[0, 1], [-1, 0]])


def pt2(m, n):
    """Teiltransposition im zweiten Faktor von C^n (x) C^n."""
    out = sp.zeros(n * n, n * n)
    for i in range(n):
        for j in range(n):
            for k in range(n):
                for l in range(n):
                    out[i * n + l, k * n + j] = m[i * n + j, k * n + l]
    return out


z, cc = sp.symbols("z c")
e2 = sp.kronecker_product(sp.eye(2), eps)
kreuz2 = (e2 * pt2(z * E4 + hb * P2, 2) * e2.inv() + ((-z - hb) * E4 + hb * P2)).applyfunc(sp.expand)
kreuz_ok, w_ok = True, True
for n in (2, 3):
    K = pt2(perm_sym(n), n)
    en = sp.eye(n * n)
    prod = ((-z - cc) * en + hb * K) * ((z - cc) * en + hb * K)                # R_VV*(z) R_VV*(-z), R_VV* = R(-z - c)^t2
    rest = (prod - (cc ** 2 - z ** 2) * en).applyfunc(sp.expand)
    lsg = sp.solve([sp.expand(e) for e in rest if sp.expand(e) != 0], cc)
    werte = [lsg[cc]] if isinstance(lsg, dict) else [t[0] if isinstance(t, tuple) else t for t in lsg]
    kreuz_ok &= K * K == n * K and werte == [n * hb / 2]
    wz = -hb / (z + n * hb / 2)
    w_ok &= nullmat((((-z - n * hb / 2) * en + hb * K) / (-(z + n * hb / 2)) - (en + wz * K)).applyfunc(sp.simplify))
rap_ok = True
for n in (2, 3):
    Pn, en = perm_sym(n), sp.eye(n * n)
    S = (th * en - (2 * sp.pi * I / n) * Pn) / (th - 2 * sp.pi * I / n)
    rap_ok &= nullmat((G(Pn, n * th / (2 * sp.pi)) - Pn * S).applyfunc(sp.simplify))
check("[E] Kreuzung: Für N = 2 ist (1 (x) eps) R(z)^(t_2) (1 (x) eps)^-1 = -R(-z - hbar). Allgemein ist R_VV*(z) = "
      "R(-z - c)^(t_2) genau für c = hbar N/2 kreuzungsunitär (R_VV*(z) R_VV*(-z) proportional zu 1; N = 2, 3) und dann "
      "gleich 1 + W K mit W = -hbar/(z + hbar N/2), CWY I Gl. 3.29. Die Kreuzungsverschiebung d_X = hbar N/2 ist die "
      "Rahmungsverschiebung einer halben Drehung (CWY I, Gl. 3.31) und von der Normierung der Spur unabhängig. In der "
      "Rapidität theta mit Dz = i d_X theta/pi ist das Gatter P (theta - (2 pi i/N) P)/(theta - 2 pi i/N) mit u = N "
      "theta/(2 pi)", nullmat(kreuz2) and kreuz_ok and w_ok and rap_ok)
ZAHLEN["Kreuzungsverschiebung"] = "d_X = hbar N/2; N = 2: d_X = hbar"
ZAHLEN["u = 1 in Rapidität"] = "theta = 2 pi/N; N = 2: theta = pi, also Dz = i d_X"
befund("Ist i hbar = Dz eine ausgezeichnete Stelle?", "nur als Lesart",
       "Exakt gilt für N = 2 (zwei Moden je Naht): Dz = i d_X, der Abstand der zwei Scharen ist i mal die Verschiebung "
       "einer halben Drehung. Die euklidische halbe Drehung selbst (Dz = +-d_X, x = +-1) gibt die Projektoren, das "
       "Singulett von P_ident. Lesart [L]: Eine Schar ist das Bild der anderen unter einer lorentzschen halben "
       "Drehung (Rapidität pi), und das ist im Bestand die Uhr R = Delta^(-i/2) (z127, z130). Das würde P_uhr durch "
       "P_fr ersetzen ('die zweite Schar ist das Uhrbild der ersten bei gleichem nacktem Spektralparameter'). "
       "Hergeleitet ist es nicht: Der native Modulfluss ist elliptisch (z136), seine Gleichsetzung mit einer Drehung "
       "der Linien ist Teil der Prämisse; für K Moden je Naht wäre der Wert 2 pi/K und nur für K = 2 die halbe Drehung.")

# A6: rational gegen trigonometrisch
eta, nu = sp.symbols("eta nu")
ident = sp.simplify((sp.sinh(eta + nu) - sp.sinh(eta) - sp.sinh(nu)
                     - 4 * sp.sinh(eta / 2) * sp.sinh(nu / 2) * sp.sinh((eta + nu) / 2)).rewrite(sp.exp))
wa, wb, wc = sp.symbols("w_a w_b w_c")
R6 = sp.Matrix([[wa, 0, 0, 0], [0, wb, wc, 0], [0, wc, wb, 0], [0, 0, 0, wa]])
sp_plus = sp.Matrix([[0, 1], [0, 0]])
Splus = sp.kronecker_product(sp_plus, sp.eye(2)) + sp.kronecker_product(sp.eye(2), sp_plus)
komm = (R6 * Splus - Splus * R6).applyfunc(sp.factor)
eintr = sorted({str(e) for e in komm if e != 0})
check("[E] Die trigonometrische Klasse trägt das H6-Gatter nicht: Die Sechs-Vertex-R-Matrix mit Gewichten a, b, c "
      "vertauscht mit der SU(2) genau für a = b + c. Mit a = sinh(eta + nu), b = sinh(nu), c = sinh(eta) ist a - b - c = "
      "4 sinh(eta/2) sinh(nu/2) sinh((eta + nu)/2): SU(2)-invariant nur an entarteten Stellen (c = 0, b = 0 oder a = 0). "
      "Das exakt SU(2)-invariante Gatter (z121) verlangt die rationale Klasse, also einen Doppelpol von omega",
      ident == 0 and eintr == ["-w_a + w_b + w_c", "w_a - w_b - w_c"])

# ------------------------------------------------------------------------------------------------------------------
say("\nTeil B: Die alternierende Kette als Gitter zweier Scharen von Wilson-Linien")


def perm_np(L, i, j, d=2):
    """Tausch der Plätze i, j auf (C^d)^(xL), numpy."""
    D = d ** L
    m = np.zeros((D, D), complex)
    for s in range(D):
        dig = [(s // d ** (L - 1 - k)) % d for k in range(L)]
        dig[i], dig[j] = dig[j], dig[i]
        t = sum(dg * d ** (L - 1 - k) for k, dg in enumerate(dig))
        m[t, s] = 1
    return m


def kette(L, uu=1.0):
    """Beide Schichtfolgen der alternierenden Kette aus Gattern R(u) auf einem Ring."""
    D = 2 ** L
    g = lambda i, j: (np.eye(D) + 1j * uu * perm_np(L, i, j)) / (1 + 1j * uu)
    ger, ung = np.eye(D, dtype=complex), np.eye(D, dtype=complex)
    for k in range(0, L, 2):
        ger = ger @ g(k, k + 1)
        ung = ung @ g(k + 1, (k + 2) % L)
    return {"gerade zuerst": ung @ ger, "ungerade zuerst": ger @ ung}


def transfer(L, lam, xi):
    """t(lam) = Spur_a prod_j (lam - xi_j + i P_aj), Hilfsraum a = Platz L."""
    D = 2 ** (L + 1)
    m = np.eye(D, dtype=complex)
    for j in range(L):
        m = m @ ((lam - xi[j]) * np.eye(D) + 1j * perm_np(L + 1, L, j))
    return m.reshape(2 ** L, 2, 2 ** L, 2).trace(axis1=1, axis2=3)


def prop(A, B):
    """Abstand von A zu einem Vielfachen von B (Frobenius, relativ)."""
    cfac = np.vdot(B, A) / np.vdot(B, B)
    return float(np.linalg.norm(A - cfac * B) / np.linalg.norm(A)), cfac


treffer = {}
for L in (4, 6, 8):
    Fs = kette(L)
    best = None
    for muster, sgn in (("xi = (+1/2, -1/2, ...)", 1), ("xi = (-1/2, +1/2, ...)", -1)):
        xi = [sgn * (0.5 if j % 2 == 0 else -0.5) for j in range(L)]
        cand = transfer(L, -0.5, xi) @ np.linalg.inv(transfer(L, 0.5, xi))
        for name in sorted(Fs):
            d, cfac = prop(Fs[name], cand)
            if best is None or d < best[0] - 1e-9:
                best = (d, muster, name, cfac)
    falsch = [0.3 * (1 if j % 2 == 0 else -1) for j in range(L)]
    dk = min(prop(Fs[name], transfer(L, -0.5, falsch) @ np.linalg.inv(transfer(L, 0.5, falsch)))[0] for name in sorted(Fs))
    treffer[L] = (best, dk)
ok_B1 = all(treffer[L][0][0] < TOL and abs(abs(treffer[L][0][3]) - 1) < 1e-9 and treffer[L][1] > 1e-3 for L in treffer)
ZAHLEN["Lichtkegelgitter"] = {f"L = {L}": {"Muster": treffer[L][0][1], "Schichtfolge": treffer[L][0][2],
                                           "Abstand": "< 1e-10", "Kontrolle xi = +-0,3": round(treffer[L][1], 3)}
                              for L in sorted(treffer)}
check("[E] Alternierende Kette der Gatter R(1) auf dem Ring (L = 4, 6, 8) = t(-1/2) t(1/2)^-1 bis auf eine Phase, mit der "
      "inhomogenen Transfermatrix t(lam) = Spur_a prod_j (lam - xi_j + i P_aj), xi_j = +-1/2 alternierend (z121 "
      "reproduziert; Kontrolle mit +-0,3 scheitert). In 4D-CS-Sprache (hbar/z = i/lam): zwei Scharen paralleler Linien "
      "in der topologischen Ebene, Rechts- und Linksläufer des Lichtkegelgitters, an den Stellen z = +-i hbar/2; Theta "
      "= 1/2 heißt z_R - z_L = +-i hbar. Die Nähte sind die Plätze der Kette; an jedem Gatter kreuzen sich zwei Linien "
      "und wechseln den Platz",
      ok_B1)

L = 4
D = 2 ** L
gn = lambda i, j: (np.eye(D) + 1j * perm_np(L, i, j)) / (1 + 1j)
rn = lambda i, j, xx: (np.eye(D) + xx * perm_np(L, i, j)) / (1 + xx)
treppe = gn(0, 1) @ gn(1, 2) @ gn(2, 3)
zykl = perm_np(L, 0, 1) @ perm_np(L, 1, 2) @ perm_np(L, 2, 3)
mono = rn(0, 3, -1j) @ rn(1, 3, -1j) @ rn(2, 3, -1j)
check("[E] Die Treppenordnung U_12 U_23 U_34 (z144 B7) ist exakt eine zyklische Verschiebung mal R_14 R_24 R_34: eine "
      "einzelne Linie, die alle anderen beim selben Abstand i hbar kreuzt (eine ungespurte Monodromie). Alternierende "
      "Schichtung und Treppe sind zwei verschiedene Lagen von Linien in der topologischen Ebene; die 4D-CS lässt beide "
      "zu und wählt keine", float(np.linalg.norm(treppe - zykl @ mono)) < TOL)
befund("P_ord als Ordnung in der topologischen Ebene", "Ort, keine Herleitung",
       "P_ord bekommt einen geometrischen Ort: Die Reihenfolge der Kompositionen ist die Lage der Linien in Sigma, und "
       "die R-Matrix hängt von nichts anderem ab (CWY I, 3.4). Welche Lage vorliegt (zwei Scharen oder eine Linie "
       "durch alle), ist aber die Wahl der Observablen. P_ord bleibt Prämisse.")

# ------------------------------------------------------------------------------------------------------------------
say("\nTeil C: Meromorphe Formen auf der Naht (Marken z^4 = 1, Uhrpole 0 und unendlich)")
zz, tt = sp.symbols("z t")
marken = [sp.Integer(1), I, sp.Integer(-1), -I]


def ordnungen(aa, bb):
    """Ordnungen von z^a (1 - z^4)^b dz bei 0, an jeder Marke, bei unendlich."""
    return aa, bb, -aa - 4 * bb - 2


ok_inf = True
for aa, bb in ((0, 0), (1, -1), (-2, 0), (3, -1), (2, -2)):
    f_inf = sp.simplify((zz ** aa * (1 - zz ** 4) ** bb).subs(zz, 1 / tt) * (-1 / tt ** 2))
    lead = sp.limit(f_inf / tt ** ordnungen(aa, bb)[2], tt, 0)
    ok_inf &= lead not in (0, sp.oo, -sp.oo, sp.zoo)
grad_ok = all(sum((o[0], 4 * o[1], o[2])) == -2 for o in (ordnungen(aa, bb) for aa in range(-6, 7) for bb in range(-3, 4)))
nullfrei = sorted((aa, bb) for aa in range(-6, 7) for bb in range(-3, 4) if max(ordnungen(aa, bb)) <= 0)
mit_marken = {aa: sum(max(o, 0) * mult for o, mult in zip(ordnungen(aa, -1), (1, 4, 1))) for aa in range(-6, 7)}
char = {aa: sp.simplify((I * zz) ** aa * I / zz ** aa) for aa in (0, -1, -2, 1)}
check("[E] Formen z^a (1 - z^4)^b dz haben die Ordnungen a (bei 0), b (an jeder Marke), -a - 4b - 2 (bei unendlich), "
      "Summe -2. Nullstellenfrei sind genau dz, dz/z, dz/z^2: Pole nur an den Uhrpolen. Unter der Uhr z -> i z tragen "
      "sie die Charaktere i, 1, -i. Mit einfachen Polen an den vier Marken (b = -1) hat jede Form mindestens zwei "
      "Nullstellen. In der Klasse ohne Nullstellen (CWY I: Gittermodelle) liegen die Marken also nie als Pole; mit "
      "Polen an den Marken ist man in der Klasse mit Unordnungsdefekten (CY III, Bittleston-Skinner)",
      ok_inf and grad_ok and nullfrei == [(-2, 0), (-1, 0), (0, 0)] and min(mit_marken.values()) == 2
      and char[0] == I and char[-1] == 1 and char[-2] == -I)

p0, p1, lam_s, sa = sp.symbols("pi0 pi1 lambda s")
q_m = p0 ** 4 - p1 ** 4
phi_m, phi_u = 1 / q_m, 1 / (p0 * p1) ** 2
hom_ok = all(sp.simplify(f.subs({p0: lam_s * p0, p1: lam_s * p1}, simultaneous=True) * lam_s ** 4 / f) == 1
             for f in (phi_m, phi_u))
Zv = sp.Matrix([sa * 1, sa * sp.Symbol("a1"), sa * sp.Symbol("a2"), sa * sp.Symbol("a3")])
var = [sp.Symbol("a1"), sp.Symbol("a2"), sp.Symbol("a3")]
d3 = sp.Matrix.hstack(Zv, *[Zv.diff(w) for w in var]).det()
d3s = sp.Matrix.hstack(Zv, Zv.diff(sa), Zv.diff(var[0]), Zv.diff(var[1])).det()
res = {j: [sp.simplify(sp.residue(zz ** j / (1 - zz ** 4), zz, m)) for m in marken] for j in range(4)}
res_ok = all(sp.simplify(res[j][k] + marken[k] ** (j + 1) / 4) == 0 for j in range(4) for k in range(4))
summen = {j: sp.simplify(sum(res[j])) for j in range(4)}
check("[E] Die native Quartik der Marken q = pi_0^4 - pi_1^4 gibt Phi_M = 1/q: homogen vom Grad -4, vier einfache Pole, "
      "keine Nullstelle; die Zahl der Marken ist genau der Grad, den D^3Z (Gewicht +4: D^3Z = s^4 da db dc, ohne "
      "ds-Anteil) braucht. Phi_U = 1/(pi_0 pi_1)^2 hat denselben Grad mit Doppelpolen an den Uhrpolen. Residuen von z^j "
      "dz/(1 - z^4) an der Marke zeta: -zeta^(j+1)/4; ihre Summe ist 0 für j = 0, 1, 2 (für j = 3 ist sie -1, dort liegt "
      "ein Pol bei unendlich). Für j = 0 sind die Residuen -1/4 mal die Uhrwerte 1, i, -1, -i",
      hom_ok and sp.simplify(d3 - sa ** 4) == 0 and d3s == 0 and res_ok
      and [summen[j] for j in range(4)] == [0, 0, 0, -1])
ZAHLEN["Residuen z dz/(1 - z^4) an 1, i, -1, -i"] = [str(r) for r in res[1]]

A_, B_, C_ = sp.symbols("A B C")
quad = [A_, B_, C_]                                    # A pi_0^2 + B pi_0 pi_1 + C pi_1^2, Bild unter pi_1 -> i pi_1
bild = [A_, I * B_, -C_]
inv_div = []
for mu in (1, I, -1, -I):
    lsg = sp.solve([sp.Eq(bi, mu * qi) for bi, qi in zip(bild, quad)], quad, dict=True)
    frei = sorted(str(s) for s in quad if lsg and lsg[0].get(s, s) != 0)
    if frei:
        inv_div.append((str(mu), frei))
dlog = sp.simplify(sp.diff(sp.log((1 + zz ** 2) / (1 - zz ** 2)), zz) / 4 - zz / (1 - zz ** 4))
ww = sp.Symbol("w")
dw = sp.simplify((zz / (1 - zz ** 4)) / sp.diff(zz ** 2, zz) - (sp.Rational(1, 2) / (1 - ww ** 2)).subs(ww, zz ** 2))
check("[E] Reduktion (Bittleston-Skinner 2.5): omega = Phi <pi gamma><pi gamma-Dach><d pi pi> hat zwei Nullstellen aus "
      "der Reduktion. Uhrkovariante Wahlen des Nullstellenpaars gibt es genau drei (pi_0^2, pi_0 pi_1, pi_1^2); "
      "symmetrisch unter dem Tausch der Uhrpole ist nur pi_0 pi_1. Mit Phi_M gibt das omega_1 = z dz/(1 - z^4), "
      "Residuen -1/4, +1/4, -1/4, +1/4, ungerade unter der Uhr. In w = z^2 ist omega_1 = dw/(2 (1 - w^2)) = (1/4) "
      "dlog((1 + w)/(1 - w)): auf dem Quotienten nach z -> -z die trigonometrische Form dx/(4x) mit zwei einfachen "
      "Polen, die Nullstellen sind die Verzweigungspunkte",
      inv_div == [("1", ["A"]), ("I", ["B"]), ("-1", ["C"])] and dlog == 0 and dw == 0
      and [str(r) for r in res[1]] == ["-1/4", "1/4", "-1/4", "1/4"])

g0, g1, h0, h1, phi_s = sp.symbols("g0 g1 h0 h1 varphi", real=True)
kl = lambda p, q: p[0] * q[1] - p[1] * q[0]                               # <p q>
pi_aff = (1, zz)
om_red = lambda phi, ga, gb: sp.simplify(phi.subs({p0: 1, p1: zz}) * kl(pi_aff, ga) * kl(pi_aff, gb))
max_ordn = 0
for f in (1, zz, zz ** 2, (1 - zz) * (1 + zz), (zz - I) * (zz + I), (1 - zz) * (zz - I)):
    for m in marken:
        max_ordn = max(max_ordn, 0 if sp.limit(f / (1 - zz ** 4) * (zz - m) ** 2, zz, m) == 0 else 2)
rat_a = om_red(phi_u, (1, 0), (1, 0))                    # Nullstellen doppelt bei z = 0: dz
rat_b = om_red(phi_u, (0, 1), (0, 1))                    # Nullstellen doppelt bei unendlich: dz/z^2
pcm = sp.simplify(om_red(phi_u, (1, 1), (1, -1)))        # Nullstellen an den Marken +-1
alpha, alphad = (1, 0), (0, 1)
k_wz = lambda ga, gb: sp.simplify(kl(alpha, ga) * kl(alphad, gb) + kl(alpha, gb) * kl(alphad, ga))
ph = sp.exp(I * phi_s)
k_werte = [k_wz((1, 1), (1, -1)), k_wz((1, I), (1, -I)), k_wz((1, ph), (1, -ph))]
check("[E] Aus Phi_M entsteht nie die rationale Klasse: f(z) dz/(1 - z^4) mit quadratischem f hat nur einfache Pole, "
      "auch wenn die Nullstellen auf Marken liegen (dann bleiben zwei einfache Pole, trigonometrisch). Die rationale "
      "Form verlangt einen Doppelpol: Phi_U = 1/(pi_0 pi_1)^2, Costellos WZW_4-Form mit den Uhrpolen als Polen, und "
      "eine Reduktion mit doppelter Nullstelle an einem Uhrpol gibt omega = dz bzw. dz/z^2. Diese Reduktion ist nicht "
      "euklidisch (dort sind die Nullstellen antipodal). Mit Nullstellen an einem antipodalen Markenpaar gibt Phi_U "
      "omega = (1 - z^2) dz/z^2 bis auf Vorzeichen, und der WZ-Koeffizient k = <alpha gamma><alpha^ gamma^> + "
      "<alpha gamma^><alpha^ gamma> verschwindet für jedes Paar auf dem Äquator der Uhr: das Hauptchirale Modell ohne "
      "WZ-Term", max_ordn == 0 and rat_a == 1 and sp.simplify(rat_b - 1 / zz ** 2) == 0
      and sp.simplify(pcm + (1 - zz ** 2) / zz ** 2) == 0 and k_werte == [0, 0, 0])

pb0, pb1 = sp.symbols("pib0 pib1")                      # komplex Konjugierte als unabhängige Symbole
q_anti = q_m.subs({p0: -pb1, p1: pb0}, simultaneous=True)
q_konj = q_m.subs({p0: pb0, p1: pb1}, simultaneous=True)
u_anti = ((p0 * p1) ** 2).subs({p0: -pb1, p1: pb0}, simultaneous=True)
anti_paare = sorted(str(sp.simplify(-1 / sp.conjugate(m))) for m in marken)
res1 = res[1]
passend = [pm for pm in (((0, 1), (2, 3)), ((0, 3), (1, 2)), ((0, 2), (1, 3)))
           if all(sp.simplify(res1[i] + res1[j]) == 0 for i, j in pm)]
uhr = lambda pm: tuple(sorted(tuple(sorted(((i + 1) % 4, (j + 1) % 4))) for i, j in pm))
check("[E] Realität und Randbedingungen: Unter der Antipodenabbildung ist q -> -q-quer und (pi_0 pi_1)^2 -> sein "
      "Konjugiertes (Phi_M ist bis auf den Faktor i, Phi_U ohne Faktor mit der euklidischen Realstruktur verträglich), "
      "und die Marken bilden zwei antipodale Paare. Für omega_1 sind die Residuen +-1/4 "
      "abwechselnd; eine diagonale Lagrangesche Randbedingung (Eichfeld an zwei Polen entgegengesetzten Residuums "
      "gleichgesetzt, Form kappa - kappa) gibt es in genau zwei Paarungen, {1 mit i, -1 mit -i} und {1 mit -i, -1 mit "
      "i}. Die Uhr vertauscht sie: Die Randbedingung ist ein Bit",
      sp.expand(q_anti + q_konj) == 0 and sp.expand(u_anti - (pb0 * pb1) ** 2) == 0
      and anti_paare == sorted(str(m) for m in marken) and len(passend) == 2
      and uhr(passend[0]) == tuple(sorted(tuple(sorted(pp)) for pp in passend[1])))
befund("Kanonische meromorphe Form aus den nativen Daten", "ja, aber die falsche Klasse",
       "I_pi legt die Faserung Z -> pi fest, die Marken die Quartik q: Omega_M = D^3Z/(pi_0^4 - pi_1^4) ist kanonisch, "
       "homogen vom Grad 0, mit vier einfachen Polen und ohne Nullstelle. Die Literatur gibt dafür ASD-Yang-Mills in "
       "der 'trigonometrischen' WZW_4-Form an (Bittleston-Skinner 3.3), mit Randbedingungen aus Lagrangeschen "
       "Unteralgebren. Die Reduktion ist trigonometrisch mit Unordnungsdefekten, nie rational. Die XXX-Kette verlangt "
       "die rationale Klasse, also Omega_U = D^3Z/(pi_0 pi_1)^2 (Doppelpole an den Uhrpolen, Costellos WZW_4) mit "
       "einer nicht euklidischen Nullreduktion. I_omega geht in keine der beiden Formen ein; die Faserung über I_omega "
       "gäbe die gespiegelte Konstruktion des Partners.")

# ------------------------------------------------------------------------------------------------------------------
say("\nTeil D: Anomalie-Bedingung (quartische Casimir-Identität) für die Gruppen des Bestands")


def T4(ws):
    w = np.array(ws, dtype=np.int64)
    return np.einsum("ni,nj,nk,nl->ijkl", w, w, w, w)


def T2(ws):
    w = np.array(ws, dtype=np.int64)
    return np.einsum("ni,nj->ij", w, w)


def T3(ws):
    w = np.array(ws, dtype=np.int64)
    return np.einsum("ni,nj,nk->ijk", w, w, w)


def sym22(A, B):
    """Dreifaches des symmetrisierten Produkts zweier quadratischer Formen (als Tensor vierter Stufe)."""
    s = (np.einsum("ij,kl->ijkl", A, B) + np.einsum("ik,jl->ijkl", A, B) + np.einsum("il,jk->ijkl", A, B)
         + np.einsum("ij,kl->ijkl", B, A) + np.einsum("ik,jl->ijkl", B, A) + np.einsum("il,jk->ijkl", B, A))
    return s // 2 if np.all(s % 2 == 0) else None


def verhaeltnis(Q4, Q2):
    """c mit Q4 = c (Q2)^2 als Polynome, sonst None. Tensorvergleich: 3 Q4 = c sym22(Q2, Q2)."""
    S = sym22(Q2, Q2)
    idx = np.argwhere(S != 0)
    if len(idx) == 0:
        return None
    i0 = tuple(idx[0])
    c = Fr(3 * int(Q4[i0]), int(S[i0]))
    return c if np.array_equal(3 * Q4 * c.denominator, c.numerator * S) else None


def wurzeln_A(n):
    out = []
    for i in range(n + 1):
        for j in range(n + 1):
            if i != j:
                vv = [0] * (n + 1)
                vv[i], vv[j] = 1, -1
                out.append(tuple(vv))
    return out


def wurzeln_D(n):
    out = []
    for i, j in itertools.combinations(range(n), 2):
        for si in (1, -1):
            for sj in (1, -1):
                vv = [0] * n
                vv[i], vv[j] = si, sj
                out.append(tuple(vv))
    return out


def wurzeln_E8():
    """Mit 2 skaliert: ganzzahlig."""
    out = [tuple(2 * c for c in r) for r in wurzeln_D(8)]
    for sg in itertools.product((1, -1), repeat=8):
        if sum(1 for s_ in sg if s_ < 0) % 2 == 0:
            out.append(tuple(sg))
    return out


dot = lambda p, q: sum(a_ * b_ for a_, b_ in zip(p, q))
E8 = wurzeln_E8()
r0 = E8[0]
E7 = [r for r in E8 if dot(r, r0) == 0]
r1 = next(r for r in E8 if dot(r, r0) == -4)              # Skalarprodukt -1 in der Norm-2-Skala
E6 = [r for r in E8 if dot(r, r0) == 0 and dot(r, r1) == 0]
G2 = wurzeln_A(2) + [tuple(s_ * (3 * (k == i) - 1) for k in range(3)) for i in range(3) for s_ in (1, -1)]
F4 = ([tuple(2 * c for c in r) for r in wurzeln_D(4)]
      + [tuple(2 * s_ * (k == i) for k in range(4)) for i in range(4) for s_ in (1, -1)]
      + list(itertools.product((1, -1), repeat=4)))
algebren = [("A1 = su(2)", wurzeln_A(1), 1, 2), ("A2 = su(3)", wurzeln_A(2), 2, 3), ("A3 = su(4)", wurzeln_A(3), 3, 4),
            ("A4 = su(5)", wurzeln_A(4), 4, 5), ("D4 = so(8)", wurzeln_D(4), 4, 6), ("D5 = so(10)", wurzeln_D(5), 5, 8),
            ("D6 = so(12)", wurzeln_D(6), 6, 10), ("G2", G2, 2, 4), ("F4", F4, 4, 9), ("E6", E6, 6, 12),
            ("E7", E7, 7, 18), ("E8", E8, 8, 30)]
tab, ok_formel = {}, True
for name, ws, rang, hv in algebren:
    dim = rang + len(ws)
    c = verhaeltnis(T4(ws), T2(ws))
    if c is None:
        tab[name] = {"dim": dim, "Identität": "nein"}
    else:
        lam2 = 4 * hv * hv * c
        ok_formel &= c == Fr(5, 2 * (dim + 2)) and lam2 == Fr(10 * hv * hv, dim + 2)
        tab[name] = {"dim": dim, "Identität": "ja", "Tr_adj X^4/(Tr_adj X^2)^2": str(c),
                     "lambda^2 (minimale Spur)": str(lam2)}
ZAHLEN["quartische Identität"] = tab
so8 = wurzeln_D(4)
v8 = [tuple(s_ * (k == i) for k in range(4)) for i in range(4) for s_ in (1, -1)]
c_so8_v = verhaeltnis(T4(so8), T2(v8))                   # Tr_adj X^4 = c (tr_v X^2)^2
ZAHLEN["so(8): lambda^2 in der Vektorspur"] = str(c_so8_v)
ja = sorted(k for k in tab if tab[k]["Identität"] == "ja")
check("[E] Tr_adj(X^4) = lambda^2 tr(X^2)^2 auf der Cartan-Algebra (Summen über die Wurzeln, exakter Tensorvergleich): "
      "erfüllt genau für su(2), su(3), so(8), G2, F4, E6, E7, E8 unter den geprüften Algebren, mit Tr_adj X^4 = 5/(2 "
      "(dim "
      "+ 2)) (Tr_adj X^2)^2 und lambda^2 = 10 (h^v)^2/(dim + 2) in der minimalen Spur (Costello-Paquette Gl. 15, 20 "
      "nachgerechnet): 8, 9, 12, 10, 15, 18, 24, 36. In der Vektorspur ist der Wert für so(8) gleich 3. Nicht erfüllt "
      "für su(4), su(5), so(10), so(12)",
      ok_formel and ja == sorted(["A1 = su(2)", "A2 = su(3)", "D4 = so(8)", "G2", "F4", "E6", "E7", "E8"])
      and [tab[k]["lambda^2 (minimale Spur)"]
           for k in ("A1 = su(2)", "A2 = su(3)", "D4 = so(8)", "G2", "F4", "E6", "E7", "E8")]
      == ["8", "9", "12", "10", "15", "18", "24", "36"] and c_so8_v == 3)

# Produktgruppen: Standardmodell-Algebra und D5 x A3
a1s, a2s, bs, ys = sp.symbols("a1 a2 b y")
w3 = [a1s, a2s, -a1s - a2s]
adj_sm = [wi - wj for wi in w3 for wj in w3 if wi is not wj] + [2 * bs, -2 * bs]
P_sm = sp.expand(sum(f ** 4 for f in adj_sm))
q3, q2 = sp.expand(sum(wi ** 2 for wi in w3)), 2 * bs ** 2
sm_ident = sp.expand(P_sm - 9 * q3 ** 2 - 8 * q2 ** 2) == 0


def ist_quadrat(poly):
    """Ist das Polynom (bis auf eine Konstante) das Quadrat eines Polynoms?"""
    if poly == 0:
        return "null"
    konst, fakt = sp.factor_list(poly)
    return "Quadrat" if all(m % 2 == 0 for _, m in fakt) else "kein Quadrat"


d5 = [tuple(list(r) + [0, 0, 0]) for r in wurzeln_D(5)]
a3 = [tuple([0] * 5 + list(r)) for r in wurzeln_D(3)]
c_d5a3 = verhaeltnis(T4(d5 + a3), T2(d5 + a3))
check("[E] Produktgruppen: Für su(3) + su(2) + u(1) ist Tr_adj X^4 = 9 q_3^2 + 8 q_2^2 (q_i = minimale Spuren), kein "
      "Quadrat einer quadratischen Form: Ein einziges Green-Schwarz-Feld hebt die Anomalie nicht, je ein Feld für su(3) "
      "und su(2) schon (das u(1) trägt nichts bei). Für D5 x A3 ist die Identität nicht erfüllt, auch nicht "
      "faktorweise", sm_ident and ist_quadrat(P_sm) == "kein Quadrat" and c_d5a3 is None)

# Materie: Lesart 'ein Weyl-Fermion in r trägt -tr_r(X^4) bei' (Costello 2111.08879, Fn. 3 für reelles R)
bb_ = sp.Symbol("b")
kal = {}
kal["SU(2), N_f = 8"] = sp.expand(2 * (2 * bb_) ** 4 - 16 * 2 * bb_ ** 4)
x3 = sp.symbols("x1:4")
kal["SU(3), N_f = 9"] = sp.expand(sum((xi - xj) ** 4 for xi in x3 for xj in x3 if xi is not xj)
                                  - 18 * sum(xi ** 4 for xi in x3)).subs(x3[2], -x3[0] - x3[1]).expand()
for n in (4, 5):
    xs = sp.symbols(f"x1:{n + 1}")
    sub = {xs[-1]: -sum(xs[:-1])}
    Pn = sp.expand((sum((xi - xj) ** 4 for xi in xs for xj in xs if xi is not xj) - 2 * n * sum(xi ** 4 for xi in xs)
                    - 6 * sum(xi ** 2 for xi in xs) ** 2).subs(sub))
    kal[f"SU({n}), N_f = {n}: Rest minus 6 t_2^2"] = Pn
check("[E] Eichung der Zählregel an der Literatur: Mit 'jedes Weyl-Fermion in r trägt -tr_r(X^4) bei' verschwindet das "
      "Polynom exakt für SU(2) mit N_f = 8 (Costello-Paquette) und SU(3) mit N_f = 9, und für SU(N) mit N_f = N bleibt "
      "genau 6 (tr_F X^2)^2 (N = 4, 5; Costello Fn. 4: mit Axion)", all(val == 0 for val in kal.values()))

v10 = [tuple(2 * s_ * (k == i) for k in range(5)) for i in range(5) for s_ in (1, -1)]       # mit 2 skaliert
s16 = [sg for sg in itertools.product((1, -1), repeat=5) if sum(1 for s_ in sg if s_ < 0) % 2 == 0]
adj45 = [tuple(2 * c for c in r) for r in wurzeln_D(5)]
t4v, t2v = T4(v10), T2(v10)
t22 = sym22(t2v, t2v)                                    # 3 mal (tr_v X^2)^2 als Tensor


def zerlege(Q4):
    """Q4 = a tr_v X^4 + b (tr_v X^2)^2 für Weyl-invariante Quartiken von D5; exakt."""
    i_a, i_b = (0, 0, 0, 0), (0, 0, 1, 1)
    M = sp.Matrix([[int(t4v[i_a]), sp.Rational(int(t22[i_a]), 3)], [int(t4v[i_b]), sp.Rational(int(t22[i_b]), 3)]])
    aa, bb = M.solve(sp.Matrix([int(Q4[i_a]), int(Q4[i_b])]))
    ok = np.array_equal(3 * Q4 * aa.q * bb.q, 3 * aa.p * bb.q * t4v + bb.p * aa.q * t22)
    return (aa, bb) if ok else None


zer = {"45": zerlege(T4(adj45)), "16": zerlege(T4(s16)), "10": zerlege(t4v)}
faelle = {}
for n16, n10 in ((0, 0), (1, 0), (3, 0), (3, 5), (4, 6)):
    aa = zer["45"][0] - n16 * zer["16"][0] - n10 * zer["10"][0]
    bb = zer["45"][1] - n16 * zer["16"][1] - n10 * zer["10"][1]
    faelle[f"{n16} x 16 + {n10} x 10"] = (aa, bb)
ZAHLEN["Spin(10): a tr_v X^4 + b (tr_v X^2)^2"] = {k: [str(val[0]), str(val[1])] for k, val in faelle.items()}
check("[E] Spin(10): tr_45 X^4 = 2 tr_v X^4 + 3 (tr_v X^2)^2, tr_16 X^4 = -tr_v X^4 + (3/4) (tr_v X^2)^2, und die kubische "
      "Spur der 16 verschwindet. Mit n Weyl-Fermionen in der 16 und m in der 10 ist das Kastenpolynom (2 + n - m) tr_v "
      "X^4 + (3 - 3n/4) (tr_v X^2)^2. Der irreduzible Anteil 2 + n verschwindet mit der 16 allein nie: Die 16 "
      "verschlechtert die Bilanz. Für den Bestand (3 x 16) bleibt 5 tr_v X^4 + (3/4) (tr_v X^2)^2, kein Quadrat",
      zer["45"] == (2, 3) and zer["16"] == (-1, sp.Rational(3, 4)) and zer["10"] == (1, 0)
      and not np.any(T3(s16)) and faelle["3 x 16 + 0 x 10"] == (5, sp.Rational(3, 4))
      and faelle["3 x 16 + 5 x 10"] == (0, sp.Rational(3, 4)) and faelle["4 x 16 + 6 x 10"] == (0, 0))

typ = {}
for r in E8:
    key = (sum(c * c for c in r[:5]), sum(c * c for c in r[5:]))
    typ[key] = typ.get(key, 0) + 1
kleb = [r for r in E8 if all(abs(c) == 1 for c in r)]
k16_4 = [r for r in kleb if sum(1 for c in r[:5] if c < 0) % 2 == 0]
k10_6 = [r for r in E8 if sum(c * c for c in r[:5]) == 4 and sum(c * c for c in r[5:]) == 4]
adj_d5a3 = [r for r in E8 if sum(c * c for c in r[5:]) == 0 or sum(c * c for c in r[:5]) == 0]
Q_e8 = T4(adj_d5a3) - T4(k16_4) - T4(k10_6)
auf_d5 = Q_e8[:5, :5, :5, :5]
auf_a3 = Q_e8[5:, 5:, 5:, 5:]
gemischt = Q_e8[:5, :5, 5:, 5:]
check("[E] Beobachtung im E8-Bild: 248 = (45,1) + (1,15) + (10,6) + (16,4) + (16-quer,4-quer) (Wurzeln 40 + 12 + 60 + 64 + "
      "64). Weyl-Fermionen in (16,4) + (10,6), also vier 16 und sechs 10, heben das Kastenpolynom von Spin(10) exakt "
      "(identisch null auf der Cartan-Algebra von D5), ohne Zusatzfeld. Auf dem A3-Faktor und in den gemischten "
      "Termen verschwindet es nicht. Das ist nicht der Inhalt des Bestands (drei Familien, keine Fermionen in der 10)",
      typ == {(8, 0): 40, (0, 8): 12, (4, 4): 60, (5, 3): 128} and len(k16_4) == 64 and len(k10_6) == 60
      and len(adj_d5a3) == 52 and not np.any(auf_d5) and np.any(auf_a3) and np.any(gemischt))

sechs = sp.Rational(1, 6)
gen = ([wi + s_ * bs + sechs * ys for wi in w3 for s_ in (1, -1)]                 # Q  (3, 2, 1/6)
       + [-wi - sp.Rational(2, 3) * ys for wi in w3]                             # u^c (3-quer, 1, -2/3)
       + [-wi + sp.Rational(1, 3) * ys for wi in w3]                             # d^c (3-quer, 1, 1/3)
       + [s_ * bs - ys / 2 for s_ in (1, -1)]                                    # L  (1, 2, -1/2)
       + [ys])                                                                  # e^c (1, 1, 1)
kub = sp.expand(sum(f ** 3 for f in gen))
P_sm3 = sp.expand(P_sm - 3 * sum(f ** 4 for f in gen))
P_sm1 = sp.expand(P_sm - sum(f ** 4 for f in gen))
ZAHLEN["Standardmodell, drei Familien: Kastenpolynom"] = str(sp.Poly(P_sm3, a1s, a2s, bs, ys).as_expr())
check("[E] Standardmodell mit Materie (je Familie Q, u^c, d^c, L, e^c; die kubische Anomalie verschwindet identisch): "
      "Das Kastenpolynom Tr_adj X^4 - n_Fam sum_Materie (Gewicht)^4 ist für eine und für drei Familien weder null "
      "noch das Quadrat einer quadratischen Form", kub == 0 and ist_quadrat(P_sm3) == "kein Quadrat"
      and ist_quadrat(P_sm1) == "kein Quadrat")

rng = np.random.default_rng(151)


def spurstrukturen(n, proben=40):
    rows = []
    for _ in range(proben):
        X = []
        for _k in range(4):
            m = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
            m = m + m.conj().T
            X.append(m - np.trace(m) / n * np.eye(n))
        tr = lambda *ms: np.trace(np.linalg.multi_dot(ms))
        einz = [tr(X[0], X[p[0]], X[p[1]], X[p[2]]) for p in itertools.permutations((1, 2, 3))]
        dopp = [tr(X[0], X[1]) * tr(X[2], X[3]), tr(X[0], X[2]) * tr(X[1], X[3]), tr(X[0], X[3]) * tr(X[1], X[2])]
        rows.append(einz + dopp)
    M = np.array(rows)
    sv = np.linalg.svd(M, compute_uv=False)
    rel = float(np.max(np.abs(M[:, :6].sum(axis=1) - M[:, 6:].sum(axis=1))) / np.max(np.abs(M)))
    return int(np.sum(sv > 1e-9 * sv[0])), rel


raenge = {n: spurstrukturen(n) for n in (2, 3, 4)}
ZAHLEN["Rang der neun n = 4-Spurstrukturen"] = {f"su({n})": raenge[n][0] for n in sorted(raenge)}
check("[N] Anschluss an z124: Der Rang der neun n = 4-Spurstrukturen (sechs Einspur, drei Doppelspur) ist 3, 8, 9 für "
      "su(2), su(3), su(4) (z124 reproduziert). Für su(2) und su(3) gilt die polarisierte Identität 'Summe der sechs "
      "Einspurterme = Summe der drei Doppelspurterme' (tr X^4 = (1/2) (tr X^2)^2), für su(4) nicht. Die Rangdefekte aus "
      "z124 (3 und 8 statt 9) enthalten also die quartische Identität; unter den su(N) treten sie genau bei den "
      "Gruppen der Green-Schwarz-Liste auf",
      [raenge[n][0] for n in (2, 3, 4)] == [3, 8, 9] and raenge[2][1] < 1e-9 and raenge[3][1] < 1e-9
      and raenge[4][1] > 1e-3)
befund("Passt die Anomalie-Bedingung zu den Gruppen des Bestands?", "nur E8 und die Transfer-su(2)",
       "Spin(10) und SU(4) haben eine unabhängige quartische Casimir-Invariante und fallen durch, mit der 16 erst "
       "recht. Die Standardmodell-Algebra besteht nur faktorweise (zwei Zusatzfelder) und mit ihrer Materie gar "
       "nicht. E8 besteht (lambda^2 = 36), ebenso su(2), die Transfergruppe der Kette (lambda^2 = 8). Für WZW_4, die "
       "Theorie, aus der die 4D-CS käme, ist die Quantisierung nur für SO(8) bewiesen (Costello, Satz 2.2.2); SO(8) "
       "ist keine Gruppe des Bestands. In jedem Fall braucht es ein Zusatzfeld mit Term (Laplace rho)^2.")

# ------------------------------------------------------------------------------------------------------------------
say("\nTeil E: Bewertung und Urteil")
D12, D23, D13, t1, t2, t3 = sp.symbols("D12 D23 D13 t1 t2 t3")
cybe_D = ((com(P12, P13) / (D12 * D13) + com(P12, P23) / (D12 * D23) + com(P13, P23) / (D13 * D23)) * D12 * D13 * D23
          - (D13 - D12 - D23) * com(P12, P23)).applyfunc(sp.cancel)
quad_def = sp.factor(((t1 - t3) ** 2 - (t1 - t2) ** 2 - (t2 - t3) ** 2))
h11, h22, hre, him = sp.symbols("h11 h22 h_re h_im", real=True)
herm = sp.Matrix([[h11, hre + I * him], [hre - I * him, h22]])
det_quad = sp.simplify((t1 * herm).det() - t1 ** 2 * herm.det()) == 0
check("[E] Vorprüfung zur ersten entscheidenden Rechnung: Für einen Austausch r_AB = c/D_AB zwischen drei Nähten ist der "
      "Rest der klassischen Yang-Baxter-Gleichung exakt (D_AC - D_AB - D_BC) [c_12, c_23]/(D_AB D_AC D_BC), und [c_12, "
      "c_23] verschwindet nicht. In diesem Ansatz tragen Nähte also genau dann einen Spektralparameter, wenn der "
      "Austausch von einem additiven Abstand abhängt. Der native Cauchykern 1/det(x - y) (z111) ist längs jeder "
      "Geraden quadratisch im "
      "Abstand (det(t n) = t^2 det(n)); sein Defekt ist 2 (t_1 - t_2)(t_2 - t_3). Unreduziert erfüllt er die Bedingung "
      "nicht", nullmat(cybe_D) and not nullmat(com(P12, P23)) and det_quad
      and sp.simplify(quad_def - 2 * (t1 - t2) * (t2 - t3)) == 0)
bilanz = [
    ("E", "Gatter",
     "R(u) aus z121 ist die rationale 4D-CS-R-Matrix. Konsistent mit der Additivität ist 'Naht = Platz': Dz = i u hbar. "
     "u = 1 heißt Dz = i hbar, der Fixpunkt von x -> -1/x; für zwei Moden je Naht ist das i mal die "
     "Kreuzungsverschiebung (Rapidität pi)."),
    ("E", "Kette",
     "Die alternierende Kette ist das Lichtkegelgitter zweier Scharen von Wilson-Linien bei z = +-i hbar/2 (Theta = "
     "1/2); die Treppe ist eine Linie durch alle. Die Nähte sind die Plätze, nicht die Linien."),
    ("L", "Kontinuum (K16)",
     "Nach CY III und Ashwinkumar-Sakamoto-Yamazaki ist der Kontinuumspartner ein Paar chiraler Flächendefekte an "
     "z_+, z_- mit Strom-Strom-Kopplung hbar/(z_+ - z_-): masseloses Thirring- bzw. chirales Gross-Neveu-Modell. Der "
     "Limes braucht eine Skalierung; bei festem u = 1 fehlt der Parameter (wie z121)."),
    ("E", "Native Form",
     "Die Marken geben kanonisch Omega_M = D^3Z/(pi_0^4 - pi_1^4): vier einfache Pole, Grad genau -4, Residuen "
     "proportional zu Uhrwerten, Summe null. Das ist die trigonometrische Klasse (ASD-Yang-Mills in WZW_4-Form), "
     "nicht die rationale. Die XXX-Kette verlangt Doppelpole an den Uhrpolen (Costellos WZW_4-Form) und eine nicht "
     "euklidische Nullreduktion."),
    ("E", "Anomalie",
     "Die quartische Identität erfüllen su(2), su(3), so(8) und die exzeptionellen Algebren. Vom Bestand also E8 und "
     "die Transfer-su(2); Spin(10), SU(4) und das Standardmodell mit Materie nicht. Auf der Seite der 4D-CS ist e_8 "
     "die eine Algebra, für die die Herleitung der Quantengruppe ausgenommen ist (CWY II)."),
    ("L", "Welche 4D-Theorie",
     "Selbstduales Yang-Mills bzw. WZW_4 mit der Transfergruppe (su(2), allgemein gl_K) als Eichgruppe, in "
     "euklidischer oder (2,2)-Signatur, nicht unitär, mit einem Zusatzfeld vierter Ordnung. Nicht Spin(10), nicht das "
     "Standardmodell, keine Materie, keine Einstein-Gravitation."),
    ("L", "Strukturvorbehalt",
     "In der Reduktion von Bittleston-Skinner ist C die pi-Faser, und eine Twistorgerade ist {Punkt von Sigma} x C. "
     "Wilson-Linien sitzen an Punkten der Faser und laufen in Sigma. Eine Naht als ganze Twistorgerade ist dort ein "
     "Punkt von Sigma, keine Linie. Dass die Spektralebene der Kette (u aus der Laufzeit phi) die Faser ist, ist "
     "nicht gezeigt (P_spek)."),
    ("O", "Erste entscheidende Rechnung",
     "Ein-Gluon-Austausch zwischen zwei Nähten mit Omega_U: den a-b-Kern aus z111/z123 mit b = Phi_U a lesen "
     "(holomorphe Chern-Simons- statt BF-Paarung), linear an die Transferströme zweier Nähte koppeln, längs der "
     "Nullebene reduzieren. Ja: Der Term erster Ordnung ist hbar mal Casimir durch D mit einem für drei Nähte "
     "additiven D (D_AC = D_AB + D_BC), und D = i hbar tan(phi/2) gibt den Term erster Ordnung von z121. Nein: Der "
     "Kern bleibt quadratisch im Abstand wie 1/det(x - y) (z111); dann tragen Nähte keinen Spektralparameter, und "
     "die Strategie ist nur eine Analogie."),
    ("O", "Zweite Rechnung",
     "Ist der native Doppelspuranteil aus z124 (Skalar mit 1/p^4, Stärke k = 1) der Green-Schwarz-Austausch, den "
     "Costello für su(2) verlangt (lambda^2 = 8)? Ja: Der native Twistorsektor ist die anomaliefreie selbstduale "
     "Theorie mit Axion. Nein: Er ist als lokale holomorphe Theorie anomal."),
]
for st, was, text in bilanz:
    say(f"   [{st}] {was}: {text}")
praemissen = {
    "benutzt": ["P_uhr", "P_ord", "T/T' (Naht = Twistorgerade)", "P_X (Austauschfeld auf dem Twistorraum, z123)"],
    "neu nötig": ["P_Omega (Volumenform Omega_U mit Doppelpolen an den Uhrpolen; Faserung über I_pi)",
                  "P_red (Nullreduktion mit doppelter Nullstelle an einem Uhrpol, nicht euklidisch)",
                  "P_rand (Randbedingung an den Polen; bei einfachen Polen ein Bit)",
                  "P_spek (die Spektralebene der Kette ist die reduzierte Faser C)",
                  "P_GS (Zusatzfeld rho mit (Laplace rho)^2 für die Anomaliefreiheit)",
                  "P_fr (nur für die Lesart: zweite Schar = Uhrbild der ersten)"],
}
say("   Prämissen: benutzt " + ", ".join(praemissen["benutzt"]) + "; neu nötig " + "; ".join(praemissen["neu nötig"]))
befund("Urteil zur Strategie", "bedingt lohnend",
       "Eine Rechnung lohnt (Ein-Gluon-Austausch zwischen zwei Nähten mit Omega_U). Fällt sie negativ aus, ist der Weg "
       "geschlossen.")
urteil = ("Bedingt lohnend. (1) Die Einordnung ist exakt: Das H6-Gatter ist die rationale R-Matrix der "
          "4D-Chern-Simons-Theorie, und die alternierende Kette ist das Lichtkegelgitter zweier Scharen von "
          "Wilson-Linien an den Stellen z = +-i hbar/2. Mit der Additivität der Spektralparameter verträglich ist nur "
          "die Lesart 'Naht = Platz' (Dz = i u hbar); die Nähte sind nicht die Wilson-Linien. (2) u = 1 heißt Dz = i "
          "hbar. Die rationale Theorie hat eine Skalensymmetrie, also je Paar von Scharen genau eine reine Zahl "
          "Dz/hbar; hbar ist keine neue freie Zahl, aber die Theorie zeichnet den Wert i auch nicht aus. P_uhr bleibt "
          "Prämisse. Exakt ist: Dz = i hbar ist der Fixpunkt von x -> -1/x, und für zwei Moden je Naht ist es i mal "
          "die Kreuzungsverschiebung. Die Lesart 'zweite Schar = Uhrbild der ersten' (P_fr) ist nicht hergeleitet. (3) "
          "P_ord bekommt einen Ort (Lage der Linien in der topologischen Ebene), wird aber nicht hergeleitet: "
          "Alternierende Schichtung und Treppe sind beide zulässige Lagen. (4) Die nativen Marken geben kanonisch die "
          "Volumenform D^3Z/(pi_0^4 - pi_1^4), Grad genau -4, vier einfache Pole, Residuen proportional zu den "
          "Uhrwerten. Das ist die trigonometrische Klasse (ASD-Yang-Mills in WZW_4-Form) und führt nie auf die "
          "rationale; die exakt SU(2)-invariante XXX-Kette verlangt Doppelpole an den Uhrpolen (Costellos WZW_4-Form) "
          "und eine nicht euklidische Reduktion. Das native Twistorbild und die Kette passen also nicht über die "
          "Marken zusammen. (5) Die Anomalie-Bedingung erfüllen vom Bestand nur E8 und die Transfer-su(2); Spin(10), "
          "SU(4) und das Standardmodell mit seiner Materie nicht. Die Eichgruppe einer solchen 4D-Theorie wäre die "
          "Transfergruppe, und die Theorie wäre selbstduales Yang-Mills bzw. WZW_4 mit einem Zusatzfeld vierter "
          "Ordnung: nicht unitär, euklidisch oder (2,2), ohne volles Yang-Mills, ohne chirale Materie, ohne "
          "Einstein-Gravitation. (6) Für den Kandidaten: P_uhr und P_ord bleiben, keine der neun Prämissen entfällt, "
          "und es kämen P_Omega, P_red, P_rand, P_spek und P_GS hinzu. K15 und K13 werden nicht gelöst; für K16 "
          "liefert die Literatur den Kontinuumspartner (Thirring bzw. chirales Gross-Neveu mit Kopplung hbar/(z_+ - "
          "z_-)), aber keinen Skalenparameter bei festem u = 1. (7) Die erste entscheidende Rechnung ist der "
          "Ein-Gluon-Austausch zwischen zwei Nähten mit Omega_U: Gibt er hbar mal Casimir durch einen additiven "
          "Abstand, trägt die Strategie; bleibt der Kern 1/det(x - y), ist sie eine Analogie. Nichts ist befördert.")
say("\nUrteil\n   " + urteil)
print(f"\n   (Laufzeit {time.time() - T0:.1f} s, nicht Teil der Ausgabe)", flush=True)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z151_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "literatur": LITERATUR, "praemissen": praemissen,
     "urteil": urteil}, indent=1, ensure_ascii=False, default=str) + "\n")
(HERE / "ausgabe_z151.txt").write_text("\n".join(OUT) + "\n")
