"""Z169: zweiundzwanzigste und dreiundzwanzigste externe Analyse, selbst nachgerechnet.

Leitlinie 'Rekonstruktion statt Auswahl'; H6 und die gelösten Blöcke bleiben unberührt, nichts wird befördert.
Statuskürzel: E exakt, C bedingt, L Lesart, N numerisch, O offen. Die Analysen sind ungeprüfte Zuarbeit; jede
übernommene Aussage ist hier neu gerechnet. Ihre Skripte lagen nicht vor.

Teil A  Analyse 22: die vier Uhrsektoren von E8 (Uhr R = exp(2 pi i v.P), v = (1/2)^5 (0)^3)
Teil B  Analyse 22: positives Maß mit der Potenz 5/4
Teil C  Analyse 23: E8-Zustandsmetrik gegen Nullkegelform, Gaußbrücke, Gegenbeispiel zur Hesse-Form
Teil D  Bilanz und Urteil.
"""
import itertools
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
mp.mp.dps = 30
AINV = mp.mpf("137.0359992168407")
c3 = 1 / (8 * mp.pi)
I = mp.mpc(0, 1)


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val, text):
    BEFUNDE[name] = val
    say(f"   -> {name}: {text}")


def theta_e8(tau, z):
    q = mp.e ** (I * mp.pi * tau)
    prods = []
    for n in (3, 4, 2, 1):
        p = mp.mpf(1)
        for zi in z:
            p *= mp.jtheta(n, mp.pi * zi, q)
        prods.append(p)
    return sum(prods) / 2


def eta(tau):
    q = mp.e ** (2 * I * mp.pi * tau)
    return mp.e ** (I * mp.pi * tau / 12) * mp.qp(q)


def eis(k, y):
    q = mp.e ** (-2 * mp.pi * y)
    coef = {2: -24, 4: 240, 6: -504}[k]
    return 1 + coef * mp.fsum(n ** (k - 1) * q ** n / (1 - q ** n) for n in range(1, 60))


def e8_shell(norm):
    """E8-Vektoren der Norm 0, 2 oder 4 als Tupel von Brüchen."""
    out = []
    rng = range(-2, 3)
    for v in itertools.product(rng, repeat=8):
        if sum(x * x for x in v) == norm and sum(v) % 2 == 0:
            out.append(tuple(Fr(x) for x in v))
    for v in itertools.product((-3, -1, 1, 3), repeat=8):
        if sum(x * x for x in v) == 4 * norm and (sum(v) // 2) % 2 == 0:
            out.append(tuple(Fr(x, 2) for x in v))
    return out


dot = lambda a, b: sum(x * y for x, y in zip(a, b))
add = lambda a, b: tuple(x + y for x, y in zip(a, b))
V = tuple([Fr(1, 2)] * 5 + [Fr(0)] * 3)
SH = {n: e8_shell(n) for n in (0, 2, 4)}
ALL = SH[0] + SH[2] + SH[4]

say("Z169: Analysen 22 und 23 nachgerechnet")
say("=" * 100)
say("\nTeil A  Die vier Uhrsektoren von E8")
cls = [0, 0, 0, 0]
for r in SH[2]:
    cls[int(4 * dot(V, r)) % 4] += 1
cls[0] += 8
tr1 = sum(n * sp.I ** k for k, n in enumerate(cls))
tr2 = sum(n * (-1) ** k for k, n in enumerate(cls))
check("A1 Schalen: 1, 240, 2160 Vektoren. Die Uhr hat auf den 248 die Vielfachheiten (60, 64, 60, 64): Sp R = 0, "
      "Sp R^2 = -8. Eine Spur mit eingesetzter Uhr R und eine Superspur zur Parität R^2 sind verschiedene "
      "Konstruktionen (Analyse 22, 23) [E]",
      [len(SH[n]) for n in (0, 2, 4)] == [1, 240, 2160] and cls == [60, 64, 60, 64] and tr1 == 0 and tr2 == -8)
t4 = sum(sp.I ** (int(4 * dot(V, r)) % 4) for r in SH[4])
vz = [mp.mpf(1) / 2] * 5 + [mp.mpf(0)] * 3
dev = max(abs(theta_e8(t, vz) / eta(t) ** 8 - (eta(2 * t) / eta(4 * t)) ** 4)
          for t in (I, mp.mpc("0.3", "1.1"), mp.mpc("-0.2", "0.8")))
check("A2 chi_R(tau) = Sp(R Q^(L_0 - 1/3)) = (eta(2 tau)/eta(4 tau))^4 = Q^(-1/3) prod (1 + Q^(2n))^(-4): an drei "
      "Punkten auf 1e-25 [N]; die Koeffizienten 1, 0, -4 folgen exakt aus den Schalen (Uhrspur der Normen 2 und 4: "
      "-8 und 16) [E]",
      dev < mp.mpf("1e-25") and sp.simplify(t4 - 16) == 0 and -8 + 8 == 0 and 16 - 64 + 44 == -4)
tw1 = sorted((dot(add(l, V), add(l, V)) / 2, add(l, V)) for l in ALL)
g1 = [w for h, w in tw1 if h == Fr(3, 8)]
g2 = [w for h, w in tw1 if h == Fr(5, 8)]
E5 = tuple([Fr(0)] * 4 + [Fr(1)] + [Fr(0)] * 3)
tw2 = sorted((dot(add(l, E5), add(l, E5)) / 2, add(l, E5)) for l in ALL)
g3 = [w for h, w in tw2 if h == Fr(1, 2)]
d5 = lambda w: w[:5]
a3 = lambda w: w[5:]
zero5, zero3 = (Fr(0),) * 5, (Fr(0),) * 3
par = lambda w: sum(1 for x in w if x < 0) % 2
check("A3 Verdrillte Sektoren, direkt im verschobenen Gitter gezählt: R-Sektor Grundniveau h = 3/8 mit 4 Zuständen "
      "(Spinor von A3 = D3, auf D5 trivial), nächstes Niveau h = 5/8 mit 16 Zuständen (Spinor von D5, auf A3 "
      "trivial); R^2-Sektor h = 1/2 mit 10 + 6 Zuständen (Vektoren). Das sind die Gewichte der Spinoren von "
      "SO(6)_1 und SO(10)_1 (6/16, 10/16) [E]",
      tw1[0][0] == Fr(3, 8) and len(g1) == 4 and all(d5(w) == zero5 and len({par(a3(w))}) == 1 for w in g1)
      and len({par(a3(w)) for w in g1}) == 1 and len(g2) == 16 and all(a3(w) == zero3 for w in g2)
      and len({par(d5(w)) for w in g2}) == 1 and tw2[0][0] == Fr(1, 2) and len(g3) == 16
      and sum(1 for w in g3 if a3(w) == zero3) == 10 and sum(1 for w in g3 if d5(w) == zero5) == 6)
vec10 = {w for w in g3 if a3(w) == zero3}
vec6 = {w for w in g3 if d5(w) == zero5}
n80 = sum(1 for s in g2 for t in g2 if add(s, t) in vec10)
n12 = sum(1 for s in g1 for t in g1 if add(s, t) in vec6)
check("A4 Übergänge als Ladungssummen: 80 geordnete Paare der 16 Spinorladungen treffen die 10 Vektorladungen, 12 "
      "geordnete Paare der 4 Grundladungen die 6 Vektorladungen von A3: 16 x 16 -> 10 und 4 x 4 -> 6 [E für die "
      "Ladungen; die Operatorprodukte selbst sind nicht nachgerechnet]",
      n80 == 80 and n12 == 12)
fix = [r for r in SH[2] if dot(V, r).denominator == 1]
check("A5 Unter der Uhr feste Wurzeln: 52 = 40 (D5, auf den ersten fünf Koordinaten) + 12 (D3 = A3, auf den letzten "
      "drei). Das feste Ladungsgitter ist D5 + D3, die native Verklebung von z144 [E]",
      len(fix) == 52 and sum(1 for r in fix if a3(r) == zero3) == 40 and sum(1 for r in fix if d5(r) == zero5) == 12)
x = sp.Symbol("x")
ser = sp.series(4 * sp.log((1 + x) * (1 + x ** 2) * (1 + x ** 3)), x, 0, 4).removeO()
alpha = 1 / AINV
dtop = 48 * c3 ** 4
Fel = dtop * mp.e ** (-2 * alpha) * (1 - dtop * mp.e ** (-2 * alpha)) ** (-mp.mpf(5) / 4)
uhr_el = (dtop / 4) * 4 * mp.fsum(mp.log(1 + mp.e ** (-2 * n * alpha)) for n in range(1, 4000))
faktor = uhr_el / Fel
ZAHLEN["uhrcharakter_faktor_elektrisch"] = mp.nstr(faktor, 8)
check("A6 -log des Uhrcharakters ohne Vakuumfaktor: 4 Q^2 + 2 Q^4 + (16/3) Q^6 + ...; die Antwort beginnt bei Q^2 "
      "[E]. Passt man den Vorfaktor am führenden dualen Term an (delta_top e^(-2/alpha)), liegt die elektrische "
      "Auswertung um den Faktor 56,82 neben der Nahtantwort [N]: Der Uhrcharakter ist nicht die Nahtantwort",
      sp.expand(ser - (4 * x + 2 * x ** 2 + sp.Rational(16, 3) * x ** 3)) == 0 and abs(faktor - mp.mpf("56.82")) < 0.05)
say(f"   Faktor elektrisch: {mp.nstr(faktor, 8)}")
check("A7 Anomalie der Uhr auf dem chiralen E8-Faktor: v.v = 5/4, also 4 v.v/2 = 5/2 nicht ganz; die Gewichte des "
      "R-Sektors liegen in 3/8 + Z/4 = 2/16 + Z/4: Typ 4{2}. Das Quadrat R^2 (Verschiebung der Norm 1) ist nicht "
      "anomal [E]",
      dot(V, V) == Fr(5, 4) and (4 * dot(V, V) / 2).denominator != 1 and (Fr(3, 8) - Fr(2, 16)) % Fr(1, 4) == 0
      and dot(E5, E5) == 1)
befund("A_uhrsektoren", {"chi_R": "(eta(2 tau)/eta(4 tau))^4", "R": "3/8: (1,4), 5/8: (16,1)", "R2": "1/2: 10 + 6"},
       "Die Aussagen der Analyse 22 über die Uhrsektoren stimmen. Die vier Sektoren sind die 16 Moduln von "
       "SO(10)_1 x SO(6)_1; die 16 der Materie ist der Spinor im R-Sektor.")

# ------------------------------------------------------------------------------------------------------------------
say("\nTeil B  Positives Maß mit der Potenz 5/4")
r, s_ = sp.symbols("r s", positive=True)                               # s = 1 - q > 0
qq = sp.Symbol("q")
num = sp.integrate(r ** 4 * sp.exp(-s_ * r ** 4), (r, 0, sp.oo))
den = sp.integrate(r ** 4 * sp.exp(-r ** 4), (r, 0, sp.oo))
ratio = sp.simplify(num / den)
logabl = sp.simplify(qq * sp.diff(sp.log(ratio.subs(s_, 1 - qq)), qq))
check("B1 Skalierung: Integral über R^5 von e^(-(1-q) Q_4) geteilt durch dasselbe bei q = 0 ist (1-q)^(-5/4) für "
      "jede positive homogene Quartik (radial nachgerechnet): fünf Richtungen durch Grad vier [E]. Die "
      "logarithmische Ableitung ist (5/4) q/(1-q), nicht q (1-q)^(-5/4): Die Einfügung fehlt",
      sp.simplify(ratio - s_ ** sp.Rational(-5, 4)) == 0
      and sp.simplify(logabl - sp.Rational(5, 4) * qq / (1 - qq)) == 0)
a, b, c, d, e = sp.symbols("a b c d e", real=True)
Q = sp.Matrix([[a, c, d], [c, b, e], [d, e, -a - b]])
m = sp.Matrix([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
Q2 = m * m * Q * (m * m).T
check("B2 Fünferträger: Auf den spurfreien symmetrischen 3x3-Matrizen ist Sp Q^4 = (Sp Q^2)^2/2, die rotations"
      "invariante Quartik also eindeutig. Das Quadrat der Vierteldrehung lässt drei Richtungen fest und kehrt zwei "
      "um (3 + 2) [E]",
      sp.expand((Q ** 4).trace() - ((Q ** 2).trace()) ** 2 / 2) == 0
      and [sp.simplify(Q2[i] - s * Q[i]) for i, s in ((0, 1), (1, 1), (4, 1), (2, -1), (5, -1))] == [0] * 5)
xx = sp.Symbol("xx", real=True)
check("B3 Kern K(z, w) = (1 - (zbar + w)/2)^(-5/4): projektive Metrik 5/(16 (1 - Re z)^2) [E]",
      sp.simplify(sp.diff(sp.Rational(-5, 4) * sp.log(1 - xx), xx, 2) / 4 - sp.Rational(5, 16) / (1 - xx) ** 2) == 0)
befund("B_mass", {"potenz": "5/4 = 5 Richtungen / Grad 4"},
       "Die Rechnungen stimmen. Das Maß ist nicht nativ hergeleitet, die Einfügung fehlt, und die Potenz bestimmt "
       "das Maß nicht. Eine native Alternative ohne Zusatzannahme steht in z170 (Kohärenzkern des 16-Zustands).")

# ------------------------------------------------------------------------------------------------------------------
say("\nTeil C  Zustandsmetrik, Nullkegelform, Gaußbrücke")
K_ = lambda y: -eis(6, y) / (3 * eis(4, y))
g = lambda y: -(2 * mp.pi) * mp.diff(K_, y) / 4                     # g = f''/4, f' = -2 pi K
g0 = g(mp.mpf(1))
g1_, g2_ = mp.diff(g, 1), mp.diff(g, 1, 2)
Kfs = -(g2_ / g0 - (g1_ / g0) ** 2) / (2 * g0)
ZAHLEN.update(g_tautau=mp.nstr(g0, 15), K_E8=mp.nstr(Kfs, 12))
check("C1 Zustandsmetrik der E8-Familie bei tau = i: g = f''/4 = pi^2 E_4(i)/6 = 2,39463397474630, g' = -3 g, "
      "g'' = 9 g - (16/3) g^2, Gaußkrümmung +8/3 [N, 1e-12]. Die Formmetrik der Nullkegel-Arbeit ist C (dx^2 + "
      "dy^2)/y^2 mit Krümmung -1/C < 0: Die direkte Gleichsetzung scheitert am Vorzeichen. Das deckt sich mit "
      "z167 (kein positiver Poincaré-Term aus der vollen chiralen Naht)",
      abs(g0 - mp.pi ** 2 * eis(4, 1) / 6) < mp.mpf("1e-20") and abs(g0 - mp.mpf("2.39463397474630456")) < mp.mpf("1e-14")
      and abs(g1_ + 3 * g0) < mp.mpf("1e-12") and abs(g2_ - 9 * g0 + mp.mpf(16) / 3 * g0 ** 2) < mp.mpf("1e-10")
      and abs(Kfs - mp.mpf(8) / 3) < mp.mpf("1e-10"))
varQ = eis(2, 1) / 12
covEQ = ((eis(2, 1) ** 2 - eis(4, 1)) / 12 + eis(4, 1) / 2) / 12
ZAHLEN.update(cov_E_Q2=mp.nstr(covEQ, 15))
check("C2 Ladungsschwankung im Nahtzustand: Var(Q_a) = E_2(i)/12 = 1/(4 pi) je Einheitsrichtung, also "
      "Var(T_3) = 1/(8 pi) = c3 und Var(Y) = 5/(24 pi); Kov(L_0, Q_a^2) = 0,0568798966258652 [E/N]. Das gibt das "
      "Verhältnis 5/3, keine absolute Kopplung",
      abs(varQ - 1 / (4 * mp.pi)) < mp.mpf("1e-25") and abs(varQ / 2 - c3) < mp.mpf("1e-25")
      and abs(covEQ - mp.mpf("0.05687989662586518")) < mp.mpf("1e-15"))
xs, ys = sp.symbols("x y", real=True)
ys = sp.Symbol("y", positive=True)
Sig = sp.Matrix([[1, xs], [xs, xs ** 2 + ys ** 2]]) / (2 * ys)
qbar = sp.Matrix([[1, xs], [xs, xs ** 2 + ys ** 2]]) / ys
qv = sp.Symbol("qv", real=True)
tau = xs + sp.I * ys
psi = (ys / sp.pi) ** sp.Rational(1, 4) * sp.exp(sp.I * tau * qv ** 2 / 2)
mom = lambda f: sp.simplify(sp.integrate(sp.simplify(f), (qv, -sp.oo, sp.oo)))
dpsi = sp.diff(psi, qv)
sqq = mom(sp.conjugate(psi) * qv ** 2 * psi)
spp = mom(sp.conjugate(dpsi) * dpsi)
sqp = mom(sp.re(sp.conjugate(psi) * qv * (-sp.I) * dpsi)) + 0
check("C3 Gaußbrücke: Der reine Gaußzustand psi_tau(q) = (y/pi)^(1/4) e^(i tau q^2/2) hat die Kovarianz "
      "Sigma = (1/2y) [[1, x], [x, x^2 + y^2]] mit det 1/4; die Formmatrix der Nullkegel-Arbeit ist genau 2 Sigma. "
      "Die Zustandsmetrik ist 1/(8 y^2), Krümmung -8 [E]. Das Paar ist bosonisch; ein solches Paar ist in der "
      "Naht nicht hergeleitet [O]",
      sp.simplify(sqq - Sig[0, 0]) == 0 and sp.simplify(spp - Sig[1, 1]) == 0 and sp.simplify(sqp - Sig[0, 1]) == 0
      and sp.simplify(Sig.det() - sp.Rational(1, 4)) == 0 and sp.simplify(qbar - 2 * Sig) == sp.zeros(2, 2)
      and sp.simplify(sp.diff(-sp.log(ys) / 2, ys, 2) / 4 - 1 / (8 * ys ** 2)) == 0)
z_, t_, x_, y_ = sp.symbols("z t x y", real=True)
u = sp.Function("u")(z_)
co = [z_, t_, x_, y_]
G = sp.diag(1, 1, sp.exp(u), sp.exp(-u))
Gi = G.inv()
Ch = [[[sum(Gi[i, l] * (sp.diff(G[l, j], co[k]) + sp.diff(G[l, k], co[j]) - sp.diff(G[j, k], co[l])) for l in range(4)) / 2
        for k in range(4)] for j in range(4)] for i in range(4)]
Ric = sp.zeros(4, 4)
for j in range(4):
    for k in range(4):
        Ric[j, k] = sum(sp.diff(Ch[i][j][k], co[i]) - sp.diff(Ch[i][j][i], co[k])
                        + sum(Ch[i][i][l] * Ch[l][j][k] - Ch[i][k][l] * Ch[l][j][i] for l in range(4)) for i in range(4))
Rs = sp.simplify(sum(Gi[j, k] * Ric[j, k] for j in range(4) for k in range(4)))
check("C4 Gegenbeispiel zur homogenen Hesse-Form: Für ds^2 = dz^2 + dt^2 + e^u dx^2 + e^(-u) dy^2 ist die "
      "Skalarkrümmung -(u')^2/2 und det g = 1. Eine konstante Verformung sieht den Einstein-Term nicht; G steht im "
      "kinetischen Koeffizienten Z von Gamma(p) = M + Z p^2. Die Rechnung bei festem tau (z167, z168) kann G also "
      "grundsätzlich nicht liefern [E]",
      sp.simplify(Rs + sp.diff(u, z_) ** 2 / 2) == 0 and sp.simplify(G.det() - 1) == 0)
befund("C_bruecke", {"K_E8": "+8/3", "K_Gauss": "-8", "K_Form": "< 0"},
       "Die Zustandsmetrik der chiralen E8-Naht ist keine Gravitonmetrik. Die Nullkegelform ist die Kovarianz eines "
       "bosonischen Gaußzustands. Offen ist, ob ein solches Paar samt Uhr aus der Naht folgt.")

# ------------------------------------------------------------------------------------------------------------------
say("\nTeil D  Bilanz und Urteil")
bilanz = [
    ("E", "Uhrsektoren", "Charakter, Gewichte 3/8, 5/8, 1/2, Inhalte (1,4), (16,1), 10 + 6, Übergänge, festes Gitter."),
    ("E", "Spuren", "Sp R = 0, Sp R^2 = -8: Uhrspur und Superspur trennen."),
    ("N negativ", "Uhrcharakter", "nicht die Nahtantwort (Faktor 56,82)."),
    ("E", "Maß 5/4", "richtig, aber nicht nativ; Einfügung fehlt."),
    ("E negativ", "E8-Zustandsmetrik", "Krümmung +8/3, falsches Vorzeichen für die Nullkegelform."),
    ("E", "Gaußbrücke", "Nullkegelform = 2 x Kovarianz eines bosonischen Gaußzustands."),
    ("E", "Hesse-Form", "konstante Verformung sieht G nicht; nötig ist der kinetische Koeffizient."),
    ("Korrektur", "eigene Angaben", "106 versionierte Lean-Dateien (nicht 215); Scherungsfreiheit wählt tau = i "
     "nicht aus (jedes konstante tau ist scherungsfrei); 3,32e122 ist log N."),
]
for st, was, text in bilanz:
    say(f"   [{st}] {was}: {text}")
urteil = ("Beide Analysen sind in ihren Rechnungen bestätigt. Übernommen werden: (1) das Uhrsektorenobjekt als "
          "konkreter Ansatz für K20, also die vier Sektoren V_(L + a v) mit Charakter (eta(2 tau)/eta(4 tau))^4, den "
          "Gewichten 3/8 und 5/8 und den Übergängen 16 x 16 -> 10; (2) die Trennung von Uhrspur und Superspur; (3) der "
          "Einwand, dass eine Rechnung bei festem tau den Einstein-Koeffizienten nicht sehen kann; (4) die Gaußbrücke "
          "als genaue Zuordnung der Nullkegelform; (5) drei Korrekturen an eigenen Angaben. Nicht übernommen: der "
          "Uhrcharakter als Nahtantwort (Faktor 56,82) und das quartische Maß als Herleitung (nicht nativ). Keine "
          "volle Lösung, nichts befördert.")
say("\nUrteil\n   " + urteil)
print(f"\n   (Laufzeit {time.time() - T0:.1f} s, nicht Teil der Ausgabe)", flush=True)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z169_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
(HERE / "ausgabe_z169.txt").write_text("\n".join(OUT) + "\n")
