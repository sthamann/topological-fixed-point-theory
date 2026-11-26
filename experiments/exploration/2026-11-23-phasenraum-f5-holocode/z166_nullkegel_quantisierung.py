"""Z166: Prüfung der Arbeit 'An algebraic quantization of gravity from the null cone' (V. Chandrasekaran,
arXiv:2611.34924, eingereicht 25.11.2026) gegen den Bestand.

Leitlinie 'Rekonstruktion statt Auswahl'; H6 und die gelösten Blöcke (Flavour, Massen, alpha, 16, E8) bleiben unberührt,
nichts wird befördert. Statuskürzel: E exakt, C bedingt (Prämisse genannt), L Lesart, N numerisch, O offen,
Lit Aussage der Arbeit, hier nicht geprüft.

Gelesen (HTML-Fassung, v1): Einleitung, 3.1 (Wirkung), 4.3 und 4.4 (Uhroperator, Endlichkeit, Paaramplitude), Anfang
von 5.1, 5.2 (endliche Dimension, Boostgewicht), 5.3 (Entropie, Logarithmus), 6 (halbklassischer Grenzfall),
Danksagung, Anfang von Anhang E (Axiome). Nicht gelesen: Abschnitt 2 im Detail, 3.2 bis 3.6, 4.1, 4.2, 5.4, Anhänge A
bis D. Die Arbeit ist drei Tage alt, nicht begutachtet, und nach eigener Angabe wurden mehrere Rechnungen und Beweise
zusammen mit einem Sprachmodell entwickelt. Sie zählt hier deshalb als ungeprüfte Literatur; geprüft sind nur die
unten nachgerechneten Formeln.

Was die Arbeit tut (Lit): Vakuumgravitation in vier Dimensionen mit Lambda, Start von der Plebanski-Form erster
Ordnung (3.2), (3.8); reduzierter Phasenraum auf zwei sich schneidenden Nullkegeln eines Kausaldiamanten mit Ecke S;
Boostladung der Ecke K = A/(8 pi G); Gaußsche Zustände, Algebra vom Typ I; Summen über Zwischenzustände konvergieren
absolut (Cauchy-Schwarz); Paaramplitude (4.33); endlichdimensionaler Hilbertraum eines Diamanten beschränkter Größe;
Dichteoperator exp(-A/(4 G hbar))/Z; Logarithmus-Koeffizient 106/45; Einstein-Gleichungen mit Lambda bis eine
Schleife. G_N und Lambda sind Eingaben. Keine Materie, keine Eichfelder, keine Spinoren oder Twistoren.

Teil A  Nachrechnung von fünf Formeln der Arbeit
Teil B  Brücken zum Bestand: (B1) Startwirkung = MacDowell-Mansouri-Form aus z165 minus Euler-Term; (B2) Lambda-Stufe
        als Zustandszahl des kosmischen Diamanten; (B3) Art der Nichtanalytizität; (B4) Nahtfaktor als Modensumme;
        (B5) Boostperiode 2 pi gegen KMS 4 pi; (B6) Scherung als Modul tau, scherungsfrei = tau = i
Teil C  Bilanz und Urteil.
"""
import json
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
AINV = mp.mpf("137.0359992168407")          # Wurzel der alpha-Gleichung des Ledgers (Bestand, Eingabe)


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val, text):
    BEFUNDE[name] = val
    say(f"   -> {name}: {text}")


say("Z166: Nullkegel-Quantisierung der Gravitation (arXiv:2611.34924) gegen den Bestand")
say("=" * 100)
say("\nTeil A  Nachrechnung")
# A1 Paaramplitude (4.33)
say("   A1 Paaramplitude: Varianz (Delta q_k)^2 = 1/(2 s k), Paarphase e^(-2ik/s), s = A/(4 G hbar)")
ok = True
for s_ in (mp.mpf(3), mp.mpf(10), mp.mpf(50)):
    direct = mp.nsum(lambda k: 2 * (1 / (2 * s_ * k)) ** 2 * mp.e ** (-2j * k / s_), [1, mp.inf])
    closed = mp.polylog(2, mp.e ** (-2j / s_)) / (2 * s_ ** 2)
    ok = ok and abs(direct - closed) < mp.mpf("1e-25")
check("A1 Summe über alle Radialmoden = Li_2(e^(-2i/s))/(2 s^2), für s = 3, 10, 50 auf 25 Stellen [N]; die Reihe "
      "konvergiert absolut für jedes s > 0, weil die Gewichte wie 1/k^2 fallen", ok)
# A2 Entwicklung in G vor der Summe (4.34)
g = sp.Symbol("g", positive=True)            # g = 1/s, proportional zu G
N_ = 7
ser = sp.series(sum(g ** 2 / 2 * sp.exp(-2 * sp.I * k * g) / k ** 2 for k in range(1, N_ + 1)), g, 0, 4).removeO()
c3g = sp.simplify(ser.coeff(g, 3))
check("A2 Entwickelt man bei endlicher Modenzahl N zuerst in g = 1/s, ist der Koeffizient von g^3 gleich -i H_N "
      "(harmonische Zahl), wächst also wie log N: die logarithmische Divergenz der Arbeit (4.34), hier für N = 7 [E]",
      sp.simplify(c3g + sp.I * sp.harmonic(N_)) == 0)
# A3 Summe vor der Entwicklung: endlicher Term -i log s / s^3
s_big = mp.mpf(10) ** 6
amp = mp.polylog(2, mp.e ** (-2j / s_big)) / (2 * s_big ** 2)
im_rest = mp.im(amp) * s_big ** 3 + mp.log(s_big)        # soll gegen ln 2 - 1 gehen
re_rest = mp.re(amp) * 2 * s_big ** 2                    # soll pi^2/6 - pi/s + 1/s^2 sein
ZAHLEN.update(A3_im_rest=mp.nstr(im_rest, 12), A3_soll=mp.nstr(mp.log(2) - 1, 12))
check("A3 Summiert man zuerst, ist Im(Amplitude) = (-log s + ln 2 - 1)/s^3 + O(1/s^5) und Re(Amplitude) 2 s^2 = "
      "pi^2/6 - pi/s + 1/s^2 exakt: aus log N wird log s. Die Amplitude ist für jedes G > 0 endlich, ihre dritte "
      "Ableitung nach G divergiert bei G -> 0+ [N bei s = 10^6; die Realteilformel ist das Bernoulli-Polynom, E]",
      abs(im_rest - (mp.log(2) - 1)) < mp.mpf("1e-10")
      and abs(re_rest - (mp.pi ** 2 / 6 - mp.pi / s_big + 1 / s_big ** 2)) < mp.mpf("1e-30"))
# A4 Boostgewicht und Dichteoperator
A_, G_, hb = sp.symbols("A G hbar", positive=True)
K_ = A_ / (8 * sp.pi * G_)
check("A4 Boostladung K = A/(8 pi G): halbe Drehung (i pi) gibt e^(-A/(8 G hbar)), zwei Kopien geben "
      "e^(-2 pi K/hbar) = e^(-A/(4 G hbar)) wie (5.23) [E]",
      sp.simplify(sp.exp(-sp.pi * K_ / hb) - sp.exp(-A_ / (8 * G_ * hb))) == 0
      and sp.simplify(sp.exp(-2 * sp.pi * K_ / hb) - sp.exp(-A_ / (4 * G_ * hb))) == 0)
spec = [Fr(3, 2), Fr(5, 2), Fr(5, 2), Fr(4), Fr(11, 2)]                    # beliebiges endliches Spektrum von s
wts = [mp.e ** (-mp.mpf(x.numerator) / x.denominator) for x in spec]
Z_ = mp.fsum(wts)
p_ = [w / Z_ for w in wts]
S_vn = -mp.fsum(p * mp.log(p) for p in p_)
mean_s = mp.fsum(p * mp.mpf(x.numerator) / x.denominator for p, x in zip(p_, spec))
dim_ = Z_ * mp.fsum(p / w for p, w in zip(p_, wts))
check("A4b Für rho = e^(-s)/Z auf endlich vielen Zuständen gilt S = <s> + log Z (5.50) und dim = Z Sp(rho e^s) "
      "(5.39c): Der Flächenterm A/(4 G hbar) ist der Mittelwert des Boostgewichts, die Korrektur ist log Z [E, "
      "an einem Beispielspektrum N]",
      abs(S_vn - (mean_s + mp.log(Z_))) < mp.mpf("1e-30") and abs(dim_ - len(spec)) < mp.mpf("1e-30"))
# A5 Logarithmus-Koeffizient
qv = Fr(257, 302)
coef = qv / (2 * (1 - qv)) - Fr(1, 2)
check("A5 Mit q = 257/302 (5.68, Lit) ist q/(2(1 - q)) - 1/2 = 106/45 (5.69) [E für die Arithmetik]; das ist "
      "Sens Koeffizient 212/45 je log des Horizontradius, also 106/45 je log der Fläche. Woher 257/302 kommt "
      "(Anhang D), ist hier nicht geprüft", coef == Fr(106, 45) and 2 * coef == Fr(212, 45))
befund("A_nachrechnung", {"geprueft": 5, "paaramplitude": "Li_2(e^(-2i/s))/(2 s^2)", "log_koeff": "106/45"},
       "Die fünf nachgerechneten Formeln stimmen. Nicht geprüft: die Herleitung des reduzierten Phasenraums, die "
       "Typ-I-Aussage, der Beweis der endlichen Dimension (5.2), q = 257/302 und der Schleifenterm.")

# ------------------------------------------------------------------------------------------------------------------
say("\nTeil B  Brücken zum Bestand")
say("   B1 Startwirkung gegen die MacDowell-Mansouri-Form aus z165")
E_, X_, V_, Lam, Gn = sp.symbols("E X V Lambda G", positive=True)     # E = eps R R, X = eps e e R, V = eps e e e e
a_ = Lam / 3
mm = E_ - 2 * a_ * X_ + a_ ** 2 * V_                                  # eps (R - a ee)(R - a ee), 2-Formen vertauschen
paper = (X_ - Lam / 6 * V_) / (32 * sp.pi * Gn)                       # (3.8)
c_mm = -3 / (64 * sp.pi * Gn * Lam)                                   # Kopplung aus z165 B1
Rs = sp.Symbol("R")
check("B1 c eps F F mit F = R - (Lambda/3) e e und c = -3/(64 pi G Lambda) (z165 B1) ist genau die Wirkung (3.8) "
      "der Arbeit plus c mal Euler-Dichte: Die Arbeit quantisiert die Form, die der Kandidat unter P_MM verlangt, "
      "bis auf einen topologischen Term [E]", sp.expand(c_mm * mm - paper - c_mm * E_) == 0)
check("B1b (3.8) = (R - 2 Lambda)/(16 pi G) mal Volumenform, mit eps e e R = 2 R vol und eps e e e e = 24 vol [E]",
      sp.simplify(paper.subs({X_: 2 * Rs, V_: 24}) - (Rs - 2 * Lam) / (16 * sp.pi * Gn)) == 0)
say("   B2 Lambda-Stufe als Zustandszahl")
c3 = 1 / (8 * mp.pi)
dtop = 48 * c3 ** 4
q_dual = dtop * mp.e ** (-2 * AINV)                                   # q(1/alpha)
LG = 3 / (32 * mp.pi ** 3) * mp.e ** (-2 * AINV)                      # Lambda G (z165 A5)
S_dS = 3 * mp.pi / LG                                                 # A/(4G) mit A = 12 pi/Lambda
rho = 3 / (4 * mp.pi ** 2) * mp.e ** (-2 * AINV)                      # rho_Lambda/Mquer^4
ZAHLEN.update(S_dS=mp.nstr(S_dS, 8), log10_S_dS=mp.nstr(mp.log10(S_dS), 8), q_dual=mp.nstr(q_dual, 8),
              eins_durch_s=mp.nstr(1 / S_dS, 8))
check("B2 Mit der Lambda-Stufe des Ledgers ist das Boostgewicht des kosmischen Horizonts s = A/(4G) = 3 pi/(Lambda G) "
      "= 32 pi^4 e^(2/alpha) = (3/8)/q(1/alpha) = 3,3243e122; gleichwertig s rho_Lambda/Mquer^4 = 24 pi^2 = "
      "(3/8)/c3^2. Der Schleifenparameter der Arbeit, 1/s, ist also (8/3) mal das Nahtnomen am dualen Punkt [E für "
      "die Identitäten; die Lambda-Stufe selbst bleibt Lesart des Ledgers]",
      abs(S_dS / (32 * mp.pi ** 4 * mp.e ** (2 * AINV)) - 1) < mp.mpf("1e-30")
      and abs(S_dS * q_dual / (mp.mpf(3) / 8) - 1) < mp.mpf("1e-30")
      and abs(S_dS * rho / (24 * mp.pi ** 2) - 1) < mp.mpf("1e-30")
      and abs(24 * mp.pi ** 2 * c3 ** 2 - mp.mpf(3) / 8) < mp.mpf("1e-30")
      and abs(mp.log10(S_dS) - mp.mpf("122.5217")) < mp.mpf("1e-3"))
say(f"   s = {mp.nstr(S_dS, 8)}, 1/s = {mp.nstr(1 / S_dS, 8)}, q(1/alpha) = {mp.nstr(q_dual, 8)}")
say("   B3 Art der Nichtanalytizität")
h_ = sp.exp(-2 / g)
check("B3 e^(-2/g) hat bei g -> 0+ lauter verschwindende Ableitungen (Taylorreihe null), die Paaramplitude dagegen "
      "eine Entwicklung mit Potenzen und einem log: verschiedene Arten. Der Mechanismus der Arbeit erzeugt die "
      "Hierarchie e^(-2/alpha) nicht [E]",
      all(sp.limit(sp.diff(h_, g, n), g, 0, "+") == 0 for n in range(0, 4)))
say("   B4 Nahtfaktor als Modensumme")
q_el = dtop * mp.e ** (-2 / AINV)
lhs = mp.e ** (mp.mpf(5) / 4 * mp.polylog(1, q_el))
check("B4 (1 - q)^(-5/4) = exp((5/4) Li_1(q)): Der Faktor der Nahtantwort ist die Exponentialfunktion einer Summe "
      "mit Gewichten 1/k, also einer Varianzsumme (Delta q_k)^2 ~ 1/k mit Faktor q je Modenzahl; die "
      "Paaramplitude der Arbeit ist dieselbe Summe mit 1/k^2 (zwei Kontraktionen). Formgleichheit, keine "
      "Herleitung: q ist dort eine Phase e^(-2i/s), bei uns reell 48 c3^4 e^(-2 alpha) [E für die Identität, L]",
      abs(lhs - (1 - q_el) ** (-mp.mpf(5) / 4)) < mp.mpf("1e-30"))
say("   B5 Boostperiode")
check("B5 Das Gewicht der Arbeit gehört zur Boostperiode 2 pi (zwei halbe Drehungen). Der KMS-Wert 4 pi des "
      "Viertelprozesses wäre in derselben Normierung e^(-2s), das Gewicht von Sp(rho^2). Ob 4 pi die doppelte "
      "Überlagerung (Spinoren, R^2 = (-1)^F) oder eine andere Normierung des Generators ist, ist ein offener "
      "Abgleich [E für die Arithmetik, O für die Zuordnung]",
      sp.simplify(sp.exp(-4 * sp.pi * K_ / hb) - sp.exp(-A_ / (4 * G_ * hb)) ** 2) == 0)
say("   B6 Scherung als Modul: die Gravitonkoordinaten der Arbeit und der Punkt tau = i")
z1, z2, dz1, dz2 = sp.symbols("z1 z2 dz1 dz2", real=True)
qbar = sp.Matrix([[sp.exp(z1), sp.exp(z1) * z2], [sp.exp(z1) * z2, sp.exp(-z1) + sp.exp(z1) * z2 ** 2]])   # (4.21a)
t1_, t2_ = z2, sp.exp(-z1)                                            # tau = z2 + i e^(-z1)
torus = sp.Matrix([[1, t1_], [t1_, t1_ ** 2 + t2_ ** 2]]) / t2_       # |dx + tau dy|^2 / Im tau
dq = qbar.diff(z1) * dz1 + qbar.diff(z2) * dz2
qi = qbar.inv()
ds2 = sp.simplify(sp.Rational(1, 2) * (qi * dq * qi * dq).trace())    # (4.21b)
poinc = (sp.diff(t1_, z1) * dz1 + sp.diff(t1_, z2) * dz2) ** 2 / t2_ ** 2 \
    + (sp.diff(t2_, z1) * dz1 + sp.diff(t2_, z2) * dz2) ** 2 / t2_ ** 2
check("B6a Die konforme Metrik (4.21a) der Arbeit ist |dx + tau dy|^2 / Im tau mit tau = z2 + i e^(-z1), det 1; ihre "
      "Feldraummetrik (4.21b) ist dz1^2 + e^(2 z1) dz2^2, die Poincaré-Metrik der oberen Halbebene. Der "
      "scherungsfreie Punkt z = 0 ist tau = i [E]",
      sp.simplify(qbar - torus) == sp.zeros(2, 2) and sp.simplify(qbar.det() - 1) == 0
      and sp.simplify(ds2 - (dz1 ** 2 + sp.exp(2 * z1) * dz2 ** 2)) == 0 and sp.simplify(ds2 - poinc) == 0)
Smat = sp.Matrix([[0, -1], [1, 0]])
check("B6b tau = i ist Fixpunkt von S: tau -> -1/tau; S hat Ordnung 4 mit S^2 = -1, das auf tau trivial wirkt. Das "
      "ist die Struktur der Uhr (R von Ordnung 4, R^2 = (-1)^F) [E; Gleichsetzung mit der Uhr ist Lesart, L]",
      sp.simplify(-1 / sp.I - sp.I) == 0 and Smat ** 2 == -sp.eye(2) and Smat ** 4 == sp.eye(2))
jf = lambda x, y: mp.kleinj(tau=mp.mpc(x, y))                         # J = j/1728
dj = mp.diff(lambda x: jf(x, 1), 0)
chi2 = lambda x, y: mp.mpf(2) / 3 * mp.log(abs(1728 * jf(x, y)))      # log |chi|^2, chi = j^(1/3) (E8, Niveau 1)
mp.mp.dps = 30
fxx = mp.diff(lambda x: chi2(x, 1), 0, 2)
fyy = mp.diff(lambda y: chi2(0, y), 1, 2)
mp.mp.dps = 40
ZAHLEN.update(B6_fxx=mp.nstr(fxx, 8), B6_fyy=mp.nstr(fyy, 8))
check("B6c Der E8-Charakter j^(1/3) ist bei tau = i stationär (j'(i) = 0, J(i) = 1). log |chi|^2 ist dort harmonisch: "
      "die zweiten Ableitungen in Re tau und Im tau sind entgegengesetzt gleich (Sattel, keine positive Form). Der "
      "positive Teil einer Graviton-Normierung kann nur aus der Metrik der Determinantenlinie kommen, also aus "
      "M1 [N für die Zahlen, E für die Harmonizität; Lesart als Graviton L]",
      abs(jf(0, 1) - 1) < mp.mpf("1e-25") and abs(dj) < mp.mpf("1e-15")
      and abs(fxx + fyy) < mp.mpf("1e-10") and abs(fxx) > 1)
say(f"   Hesse-Einträge von log |chi|^2 bei tau = i: {mp.nstr(fxx, 8)} (Re tau), {mp.nstr(fyy, 8)} (Im tau)")
befund("B_bruecken", {"s_dS": mp.nstr(S_dS, 6), "wirkung": "gleich bis auf Euler-Term"},
       "(1) Die Arbeit quantisiert genau die Gravitationsform des Kandidaten (erste Ordnung, Lambda ungleich null). "
       "(2) In ihrer Sprache ist die Lambda-Stufe eine Zustandszahl: log dim = (3/8)/q(1/alpha) für den kosmischen "
       "Diamanten. (3) Ihre Nichtanalytizität ist ein Logarithmus, unsere Hierarchie ist e^(-2/alpha): sie liefert "
       "die Hierarchie nicht. (4) G und Lambda sind dort Eingaben; die zwei freien reinen Zahlen bleiben frei.")

# ------------------------------------------------------------------------------------------------------------------
say("\nTeil C  Bilanz und Urteil")
bilanz = [
    ("E", "Nachrechnung", "Fünf Formeln der Arbeit stimmen: Paaramplitude, log-Divergenz bei Entwicklung, endlicher "
     "log s-Term, Boostgewicht und Entropieformel, Arithmetik des Koeffizienten 106/45."),
    ("E", "Pflichtenheft Punkt 3", "Die Startwirkung der Arbeit ist die MacDowell-Mansouri-Form aus z165 minus "
     "Euler-Term. Für den Gravitationsteil des Kandidaten gibt es damit eine ausgearbeitete Vorlage."),
    ("Lit", "Pflichtenheft Punkt 5, W3", "Positivität steht dort an den Zuständen: Präparation, verklebt mit ihrer "
     "Zeitspiegelung, und Cauchy-Schwarz für die Zwischensummen (Anhang E, Axiome). Das ist die Form, die z165 "
     "nach der Korrektur von W3 verlangt."),
    ("E/L", "W5", "Lambda-Stufe als Zustandszahl des kosmischen Diamanten: log dim = (3/8)/q(1/alpha)."),
    ("E negativ", "Hierarchie", "Der Mechanismus der Arbeit (erst summieren, dann entwickeln) gibt einen log, "
     "nicht e^(-2/alpha)."),
    ("O", "Zwei reine Zahlen", "G und Lambda sind Eingaben der Arbeit. W1 (Kopplung) und die Zahl G Lambda bleiben "
     "offen."),
    ("O", "Materie", "Keine Eichfelder, keine Fermionen, keine Spinoren oder Twistoren. T3 bis T8 (gemeinsamer "
     "Elter) werden nicht berührt."),
    ("O", "Abgleich", "Boostperiode 2 pi der Arbeit gegen KMS 4 pi des Viertelprozesses."),
    ("E/L", "Graviton als Modul", "Die zwei Polarisationen der Arbeit sind je Lichtstrahl ein Punkt tau der oberen "
     "Halbebene; scherungsfrei ist tau = i, der Punkt der Naht, mit der Uhr als Stabilisator. Der E8-Charakter ist "
     "dort stationär, sein Betrag gibt einen Sattel."),
]
for st, was, text in bilanz:
    say(f"   [{st}] {was}: {text}")
urteil = ("Die Arbeit hilft an drei Stellen und an den zwei entscheidenden nicht. Sie hilft, weil sie (1) genau die "
          "Gravitationsform quantisiert, die der Kandidat 'relative graduierte Determinante' verlangt (erste "
          "Ordnung mit Lambda ungleich null; ihre Wirkung ist unsere Form aus z165 minus Euler-Term), (2) die "
          "Positivität dort ansetzt, wo z165 sie nach der Korrektur von W3 hingelegt hat (Zustände und "
          "Zeitspiegelung, nicht der Vertex), und (3) eine Sprache liefert, in der die Lambda-Stufe des Ledgers "
          "eine prüfbare Zustandszahl ist: log dim des kosmischen Diamanten = (3/8)/q(1/alpha). Sie hilft nicht bei "
          "den zwei freien reinen Zahlen (G und Lambda sind dort Eingaben), nicht bei der Hierarchie (ihre "
          "Nichtanalytizität ist ein Logarithmus) und nicht beim fehlenden gemeinsamen Objekt: Sie enthält keine "
          "Materie, keine Eichfelder und keine Twistoren. Sie ist drei Tage alt, nicht begutachtet und teils mit "
          "einem Sprachmodell entwickelt; hier zählt sie als ungeprüfte Literatur, nur fünf Formeln sind "
          "nachgerechnet. Am Stand ändert sich nichts: 13 P / 7 L, keine volle Lösung, nichts befördert. Neu als "
          "Rechnungen: M5 (Naht als Eckdaten: Boostladung und Uhr der Arbeit gegen Uhr und KMS-Zustand der Naht, "
          "dabei der Abgleich 2 pi gegen 4 pi), M6 (Zustandszahl des kosmischen Diamanten aus der Nahtantwort), "
          "M7 (Graviton als Modul der Naht: zweite Ableitung der Nahtdeterminante mit Linienmetrik bei tau = i), "
          "M8 (Wechselwirkung aus der Uhrgleichung statt aus dem Zustand) und M9 (Verklebungsbit P_IF aus der "
          "Klebregel mit Außennormalen).")
say("\nUrteil\n   " + urteil)
print(f"\n   (Laufzeit {time.time() - T0:.1f} s, nicht Teil der Ausgabe)", flush=True)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z166_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
(HERE / "ausgabe_z166.txt").write_text("\n".join(OUT) + "\n")
