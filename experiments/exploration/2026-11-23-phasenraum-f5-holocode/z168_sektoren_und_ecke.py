"""Z168: Stufe 2 und Stufe 3 des Lösungswegs. Zwei Sektoren, die Naht als Ecke, Uhrgleichung, Verklebung, Materie.

Leitlinie 'Rekonstruktion statt Auswahl'; H6 und die gelösten Blöcke (Flavour, Massen, alpha, 16, E8) bleiben unberührt,
nichts wird befördert. Statuskürzel: E exakt, C bedingt (Prämisse genannt), L Lesart, N numerisch, O offen.

Bestand: z167 (Linienmetrik; t = e^2/c3), z166 (Nullkegel-Arbeit: Boostgewicht e^(-2 pi K), K = A/(8 pi G); Klebregel
mit Außennormalen), z165 (Lambda-Stufe als relative Nahtantwort; Markenmatrizen), z144 (zwei Lagrangesche Untergruppen),
VERDICT (Gesetzestorus = Nahttorus, tau_law = i beta/(4 pi); 4 pi L_0 = 2 pi L_0^NS als Erzeugerwahl).

Teil A  M5: Naht als Ecke. KMS-Gewicht, mittlere Boostladung, Entropie, Spiegelung
Teil B  M8: die Uhrgleichung der Naht (Boostladung als Funktion von Verdrillung und Temperatur)
Teil C  Stufe 2: Hierarchie. Torsionssektoren, Markenmatrizen (M2), Zustandszahl (M6)
Teil D  Stufe 3: Verklebungsbit (M9), Helizität der Materie (W4)
Teil E  Bilanz und Urteil.
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
mp.mp.dps = 40
AINV = mp.mpf("137.0359992168407")          # Wurzel der alpha-Gleichung des Ledgers (Bestand, Eingabe)
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


def eis(k, tau):
    q = mp.e ** (2 * I * mp.pi * tau)
    coef = {2: -24, 4: 240, 6: -504}[k]
    return 1 + coef * mp.nsum(lambda n: n ** (k - 1) * q ** n / (1 - q ** n), [1, mp.inf])


zero = [mp.mpf(0)] * 8
Zf = lambda tau, z=zero: theta_e8(tau, z) / eta(tau) ** 8
Kmean = lambda tau: -eis(6, tau) / (3 * eis(4, tau))          # <L_0 - 1/3> = q d/dq log Z
varpi = mp.gamma(mp.mpf(1) / 4) ** 2 / (2 * mp.sqrt(2 * mp.pi))

say("Z168: Stufe 2 und 3. Zwei Sektoren, die Naht als Ecke, Uhrgleichung, Verklebung")
say("=" * 100)
say("\nTeil A  M5: die Naht als Ecke")
# Charakter-Koeffizienten von j^(1/3) = q^(-1/3)(1 + 248 q + 4124 q^2 + ...)
NT = 14


def pmul(a, b):
    c = [0] * NT
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if i + j < NT:
                    c[i + j] += x * y
    return c


dn = [1] + [240 * sum(d ** 3 for d in range(1, n + 1) if n % d == 0) for n in range(1, NT)]
for n in range(1, NT):
    inv = [1 if k % n == 0 else 0 for k in range(NT)]
    for _ in range(8):
        dn = pmul(dn, inv)
mp.mp.dps = 30
Zsum = mp.fsum(d * mp.e ** (-2 * mp.pi * (n - mp.mpf(1) / 3)) for n, d in enumerate(dn))
mp.mp.dps = 40
ZAHLEN.update(grade=dn[:5], Z_i=mp.nstr(mp.re(Zf(I)), 15))
check("A1 Bei tau = i ist die Zustandssumme der Naht Sp e^(-2 pi K) mit K = L_0 - 1/3 (c/24 = 1/3): "
      "Summe d_n e^(-2 pi (n - 1/3)) mit d_n = 1, 248, 4124, 34752, ... gibt 12 = j(i)^(1/3) [N, 14 Stufen auf 1e-25; "
      "exakt 12]. Das Gewicht hat die Periode 2 pi im Erzeuger L_0",
      dn[:4] == [1, 248, 4124, 34752] and abs(Zsum - 12) < mp.mpf("1e-25") and abs(Zf(I) - 12) < mp.mpf("1e-30"))
be, L0 = sp.symbols("beta L0", positive=True)
check("A2 Abgleich 2 pi gegen 4 pi: Mit tau_law = i beta/(4 pi) (Bestand) ist q^(L_0) = e^(-beta L_0/2). beta = 4 pi "
      "ist die Periode im Erzeuger H_Q = L_0/2 und gibt dasselbe Gewicht e^(-2 pi L_0) wie die Boostperiode 2 pi "
      "der Nullkegel-Arbeit. Kein Faktor zwei, kein Widerspruch [E unter der Bestandsdefinition]",
      sp.simplify(sp.exp(2 * sp.pi * sp.I * (sp.I * be / (4 * sp.pi)) * L0).subs(be, 4 * sp.pi)
                  - sp.exp(-2 * sp.pi * L0)) == 0)
S_naht = 2 * mp.pi * Kmean(I) + mp.log(Zf(I))
check("A3 Mittlere Boostladung <K> = -E_6/(3 E_4) = 0 bei tau = i (das ist <T> = 0 aus H6). Entropie "
      "S = 2 pi <K> + log Z = log 12 = 2,4849, ganz aus log Z, ohne Flächenterm [E]",
      abs(Kmean(I)) < mp.mpf("1e-30") and abs(S_naht - mp.log(12)) < mp.mpf("1e-30"))
ys = [mp.mpf("0.7"), mp.mpf("1.3"), mp.mpf("2.1")]
check("A4 Spiegelung: <K>(i/y) = -y^2 <K>(i y) für y = 0,7; 1,3; 2,1 [N, 25 Stellen; folgt aus den Gewichten 6 "
      "und 4]. Unter tau -> -1/tau wechselt die Boostladung das Vorzeichen, am Fixpunkt ist sie null. Das ist "
      "die Form der Eckbedingung der Arbeit (Boostladungen mit Außennormalen heben sich auf) [E für die "
      "Identität, L für die Zuordnung K <-> A/(8 pi G)]",
      all(abs(Kmean(I / y) + y ** 2 * Kmean(I * y)) < mp.mpf("1e-25") for y in ys))
befund("A_ecke", {"Z": 12, "K_mittel": 0, "S": "log 12"},
       "Die Naht am symmetrischen Punkt ist ein Zustand e^(-2 pi K)/12 mit mittlerer Boostladung null. In der "
       "Sprache der Nullkegel-Arbeit ist das eine Ecke mit mittlerer Fläche null: der Scheitel eines Lichtkegels, "
       "passend zur Twistorgeraden als Himmel eines Punkts [L].")

# ------------------------------------------------------------------------------------------------------------------
say("\nTeil B  M8: die Uhrgleichung der Naht")
E4i = eis(4, I)
mp.mp.dps = 30
k_y1 = mp.diff(lambda y: mp.re(Kmean(I * y)), 1)
k_x2 = mp.diff(lambda x: mp.re(Kmean(mp.mpc(x, 1))), 0, 2) / 2
k_x1 = mp.diff(lambda x: mp.re(Kmean(mp.mpc(x, 1))), 0)
mp.mp.dps = 40
ZAHLEN.update(K_lin_temperatur=mp.nstr(k_y1, 12), K_quad_verdrillung=mp.nstr(k_x2, 12))
check("B1 K(tau) = -E_6/(3 E_4) ist eine exakte nichtlineare Funktion von Verdrillung x = Re tau und "
      "Temperaturmodul y = Im tau. Um tau = i: Re K = -(pi E_4(i)/3) (y - 1) - (pi E_4(i)/2) x^2 + ..., mit "
      "pi E_4(i)/3 = varpi^4/pi^3 = 1,5245 [N, 12 Stellen]",
      abs(k_y1 + mp.pi * E4i / 3) < mp.mpf("1e-12") and abs(k_x2 + mp.pi * E4i / 2) < mp.mpf("1e-12")
      and abs(k_x1) < mp.mpf("1e-12") and abs(mp.pi * E4i / 3 - varpi ** 4 / mp.pi ** 3) < mp.mpf("1e-30"))
check("B2 Vorzeichen: Verdrillung senkt die Boostladung quadratisch, ohne linearen Term. Das ist das Muster der "
      "Fokussierung (bei fester Laufstrecke senkt Scherung die Fläche). Die Temperaturrichtung geht linear ein, "
      "wie ein Boost. Die Nichtlinearität steckt in der Funktion K(tau), nicht im Zustand, der Gaußsch bleibt: "
      "S-II ist für die eine Naht bestätigt [E für die Funktion, L für die Lesart als Fokussierung]. Eine "
      "Wechselwirkung zwischen verschiedenen Geraden ist das noch nicht",
      k_x2 < 0 and abs(k_x1) < mp.mpf("1e-12") and k_y1 < 0)
befund("B_uhr", {"K": "-E6/(3 E4)", "quadratisch": mp.nstr(k_x2, 8)},
       "Die Naht hat eine eigene Uhrgleichung: K(tau) = -E_6/(3 E_4). Eine Polarisation (Verdrillung) verhält sich "
       "wie Scherung, die andere (Temperatur) wie ein Boost. Die Gleichsetzung von tau mit beiden "
       "Gravitonpolarisationen (z166, S-I) passt deshalb nur zur Hälfte.")

# ------------------------------------------------------------------------------------------------------------------
say("\nTeil C  Stufe 2: Hierarchie und zwei Sektoren")
Yv = [mp.mpf(-1) / 3] * 3 + [mp.mpf(1) / 2] * 2 + [mp.mpf(0)] * 3


def n2(tau, a, b, vec=None):
    """||Z||^2 am Punkt z = (a + b tau) vec."""
    vec = Yv if vec is None else vec
    z = [(a + b * tau) * v for v in vec]
    y2 = mp.fsum(mp.im(x) ** 2 for x in z)
    return abs(Zf(tau, z)) ** 2 * mp.e ** (-2 * mp.pi * y2 / mp.im(tau))


Yf = [Fr(-1, 3)] * 3 + [Fr(1, 2)] * 2 + [Fr(0)] * 3
T3 = [Fr(0)] * 3 + [Fr(1, 2), Fr(-1, 2)] + [Fr(0)] * 3
diff = [3 * a - b for a, b in zip(Yf, T3)]
roots = []
for i, j in itertools.combinations(range(8), 2):
    for sa in (1, -1):
        for sb in (1, -1):
            v = [Fr(0)] * 8
            v[i], v[j] = Fr(sa), Fr(sb)
            roots.append(v)
for sg in itertools.product((1, -1), repeat=8):
    if sum(1 for x in sg if x < 0) % 2 == 0:
        roots.append([Fr(x, 2) for x in sg])
gerade = sum(1 for r in roots if (3 * sum(a * b for a, b in zip(Yf, r))).denominator == 1)
check("C0 Die Zwei-Torsion der Hyperladung, e^(2 pi i 3Y), ist auf der Naht das Zentrum der schwachen SU(2): "
      "3Y - T_3 = (-1,-1,-1,1,2,0,0,0) liegt im E8-Gitter. Der Zentralisator hat 128 Wurzeln (E7 x SU(2), "
      "136 = 133 + 3), 112 = (56,2) sind ungerade. Die Verschiebung hat Norm 1/2: das ist die anomale Involution "
      "von E8 (Spin 1/4 im verdrillten Sektor) [E]",
      all(x.denominator == 1 for x in diff) and sum(diff) % 2 == 0 and len(roots) == 240 and gerade == 128
      and sum(x * x for x in T3) == Fr(1, 2))
mp.mp.dps = 30
sek = {k: n2(I, *ab) for k, ab in {"0": (0, 0), "3": (3, 0), "3i": (0, 3), "3+3i": (3, 3)}.items()}
weg = n2(mp.mpc(0, "1.2"), 3, 3)
e1 = [mp.mpf(1)] + [mp.mpf(0)] * 7
d8 = n2(I, 1, 1, e1)
per = abs(n2(I, mp.mpf("0.37") + 6, mp.mpf("0.21")) - n2(I, mp.mpf("0.37"), mp.mpf("0.21"))) \
    + abs(n2(I, mp.mpf("0.37"), mp.mpf("0.21") + 6) - n2(I, mp.mpf("0.37"), mp.mpf("0.21")))
eps = mp.mpf("1e-3")
steig = n2(I * (1 + eps), 3, 3) / eps ** 2
mp.mp.dps = 40
ZAHLEN["sektoren_norm2"] = {k: mp.nstr(v, 10) for k, v in sek.items()}
ZAHLEN.update(doppelt_bei_1_2i=mp.nstr(weg, 10), d8_doppelt=mp.nstr(d8, 10), nullstelle_steigung=mp.nstr(steig, 6))
check("C1 Die vier Sektoren der Zwei-Torsion bei tau = i: ||Z||^2 = 144 (unverdrillt), 72 und 72 (einfach "
      "verdrillt, exakt die Hälfte), 0 (doppelt verdrillt) [N, 25 Stellen]. Der doppelt verdrillte Sektor "
      "verschwindet nur am symmetrischen Punkt (bei tau = 1,2 i ist er 2,47) und nur für die anomale Involution "
      "(die nicht anomale D8-Involution gibt dort einen endlichen Wert). Perioden 6 und 6 tau geprüft",
      abs(sek["0"] - 144) < mp.mpf("1e-25") and abs(sek["3"] - 72) < mp.mpf("1e-25")
      and abs(sek["3i"] - 72) < mp.mpf("1e-25") and sek["3+3i"] < mp.mpf("1e-50") and weg > 1 and d8 > 1
      and per < mp.mpf("1e-20"))
say("   ||Z||^2: " + ", ".join(f"{k}: {mp.nstr(v, 8)}" for k, v in sek.items())
    + f"; doppelt verdrillt bei 1,2 i: {mp.nstr(weg, 8)}; D8-Involution doppelt: {mp.nstr(d8, 8)}")
check("C1b Folge für die Hierarchie: Die Verhältnisse der Sektoren sind 1, 1/2, 1/2 und 0. Eine kleine endliche "
      "Zahl wie 1e-122 kommt als Verhältnis von Torsionssektoren bei festem tau = i nicht vor; die Null ist eine "
      "doppelte Nullstelle in tau - i (||Z||^2 = 76,6 |tau - i|^2), keine Unterdrückung. Stufe 2 in der "
      "angekündigten Form ist negativ [E negativ]",
      abs(steig - mp.mpf("76.6")) < 1 and sek["3+3i"] < mp.mpf("1e-50"))
alpha = 1 / AINV
e2 = 4 * mp.pi * alpha
check("C2 Die Hierarchie braucht die Kopplung im Exponenten: e^(-2/alpha) = exp(-(1/c3)/e^2), mit "
      "1/c3 = 8 pi = |Z_2| 2 pi chi(S^2) = Spur der Hesse-Form (z167). Das ist das Gewicht einer Konfiguration mit "
      "Wirkung (Gesamtkrümmung der zweiblättrigen Naht)/e^2, linear in der Krümmung und skalenfrei [E für die "
      "Identität; die Lesart als Instantongewicht der Naht selbst ist L und ersetzt die skalenabhängige Lesart "
      "'Fluss 4 auf der Einheitssphäre' aus z165]",
      abs(2 * AINV - (1 / c3) / e2) < mp.mpf("1e-25") and abs(1 / c3 - 2 * 2 * mp.pi * 2) < mp.mpf("1e-30"))
xm = c3 * mp.log(2)
Bm = sp.Matrix(4, 4, lambda i, j: [2, -1, 0, -1][(j - i) % 4])
Cm = sp.Matrix(4, 4, lambda i, j: [0, 1, 2, 1][(j - i) % 4])
Um = sp.diag(1, -1, 1, -1)
xsym = sp.Symbol("x")
detC = sp.factor((sp.eye(4) - xsym * Cm).det())
detB = sp.factor((sp.eye(4) - xsym * Um * Bm * Um).det())
rel = mp.log(((1 + 2 * xm) / (1 - 2 * xm)) ** 2)
ldC = -mp.log((1 - 4 * xm) * (1 + 2 * xm) ** 2)
dtop = 48 * c3 ** 4
q_el = dtop * mp.e ** (-2 * alpha)
resp = q_el * (1 - q_el) ** (-mp.mpf(5) / 4)
ZAHLEN.update(M2_relativ=mp.nstr(rel, 8), M2_logdetC=mp.nstr(ldC, 8), nahtantwort=mp.nstr(resp, 8))
check("C3 M2: Die relative Determinante der Markenmatrizen bei x = c3 ln 2 ist log((1+2x)/(1-2x))^2 = 0,2209, "
      "-log det(1 - x C) = 0,009508. Die Nahtantwort q (1-q)^(-5/4) ist 1,1858e-4. Keine der Größen stimmt, auch "
      "nicht bis auf einfache Faktoren (Verhältnisse 1863 und 80,2): M2 ist negativ [N negativ]",
      sp.expand(detC - (1 - 4 * xsym) * (1 + 2 * xsym) ** 2) == 0
      and sp.expand(detB - (1 - 4 * xsym) * (1 - 2 * xsym) ** 2) == 0
      and abs(rel - mp.mpf("0.2209")) < mp.mpf("1e-3") and abs(rel / resp - 1863) < 5 and abs(ldC / resp - 80.2) < 0.5)
S_dS = 3 * mp.pi / (3 / (32 * mp.pi ** 3) * mp.e ** (-2 * AINV))
N12 = S_dS / mp.log(12)
ZAHLEN.update(naehte_je_horizont_log12=mp.nstr(N12, 6))
check("C4 M6: Trägt jede Naht die Entropie log 12 (A3), braucht der kosmische Horizont 1,338e122 Nähte. Der "
      "Bestand nennt 9,6e122 (das entspricht (ln 2)/2 je Naht). Die Zustandszahl folgt nicht aus der Naht; sie "
      "bleibt die Zahl der Nähte je Fläche, dieselbe freie Zahl wie G/R_vac^2 [N; E negativ für eine Herleitung]",
      abs(N12 / mp.mpf("1.338e122") - 1) < mp.mpf("1e-3")
      and abs(S_dS / (mp.log(2) / 2) / mp.mpf("9.59e122") - 1) < mp.mpf("1e-2"))
befund("C_hierarchie", {"torsion": "1, 1/2, 1/2, 0", "M2": "negativ", "exponent": "(1/c3)/e^2"},
       "Die Hierarchie ist keine relative Determinante der Naht bei fester Kopplung. Ihr Exponent ist die "
       "Gesamtkrümmung geteilt durch e^2: das Gewicht der Naht selbst als klassische Konfiguration eines "
       "Eichfelds mit Kopplung e. Dafür braucht es die Dynamik der Kopplung, nicht nur die Naht.")

# ------------------------------------------------------------------------------------------------------------------
say("\nTeil D  Stufe 3: Verklebung und Materie")
qD = lambda x: Fr(5 * x * x, 8) % 1
qA = lambda y: Fr(3 * y * y, 8) % 1
anti = all((qD(x) + qA(x)) % 1 == 0 for x in range(4))
elems = list(itertools.product(range(4), range(4)))
qtot = lambda v: (qD(v[0]) + qA(v[1])) % 1


def span(g):
    return {((k * g[0]) % 4, (k * g[1]) % 4) for k in range(4)}


lagr = []
for g in elems:
    s = span(g)
    if len(s) == 4 and all(qtot(v) == 0 for v in s) and s not in lagr:
        lagr.append(s)
graphs = [{(x, (sgn * x) % 4) for x in range(4)} for sgn in (1, -1)]
check("D1 M9: q(A3) = -q(D5) mod 1 auf Z_4 (3/8 = -5/8 mod 1). Die zwei Lagrangeschen Untergruppen sind genau "
      "die Graphen der zwei Anti-Isometrien y = x und y = -x. Sie unterscheiden sich um die "
      "Ladungskonjugation auf einer Seite [E]",
      anti and len(lagr) == 2 and all(g in lagr for g in graphs))
check("D2 Klebregel der Nullkegel-Arbeit: Die zweite Seite ist die Zeitspiegelung der ersten (Orientierung "
      "umgekehrt, Quellen komplex konjugiert, keine zusätzliche Ladungskonjugation), und die Summe der Ladungen "
      "mit Außennormalen ist null. Bezeichnet man die Klassen der zweiten Seite durch die Spiegelung selbst, ist "
      "die Verklebung der Graph der Identität, y = x. Die andere Verklebung gehört zur Spiegelung mit "
      "Ladungskonjugation. P_IF ist damit die Wahl der Spiegelung (J oder CJ), und unter P_TP (zweite Naht = "
      "Tomita-Partner, Bestand) ist sie getroffen [E für die Algebra, C unter P_TP, L für die Zuordnung]",
      graphs[0] in lagr and all((v[0] - v[1]) % 4 == 0 for v in graphs[0])
      and all((v[0] + v[1]) % 4 == 0 for v in graphs[1]))
hel = lambda n: Fr(-(n + 2), 2)                                 # Penrose: O(n) <-> Helizität -(n+2)/2
check("D3 W4: Die 16 Fermionen der Naht sind Schnitte von K^(1/2) = O(-1) auf der Geraden. Nach der "
      "Penrose-Zuordnung O(n) <-> Helizität -(n+2)/2 ist das einheitlich Helizität -1/2; der Serre-Partner "
      "O(-3) hat +1/2. Eichfeld: O(0) und O(-4) (Helizität -1, +1); Gravitation: O(2) und O(-6) (z164). Die "
      "Chiralität der 16 ist das Gewicht des Spinbündels der Naht [E für die Gewichte, L für die Lesart der "
      "Nahtfermionen als Raumzeitfelder]",
      hel(-1) == Fr(-1, 2) and hel(-3) == Fr(1, 2) and hel(0) == -1 and hel(-4) == 1 and hel(2) == -2
      and hel(-6) == 2 and all(n + m == -4 for n, m in ((-1, -3), (0, -4), (2, -6))))
befund("D_stufe3", {"P_IF": "Wahl der Spiegelung", "helizitaet_16": "-1/2"},
       "Die Verklebung ist auf die Wahl der Spiegelung zurückgeführt; unter P_TP entfällt P_IF als eigene "
       "Prämisse. Die 16 tragen eine einheitliche Helizität. Nicht gerechnet: Wechselwirkung zwischen Geraden, "
       "Einstein-Form in Twistorvariablen (M3), Positivität über mehrere Geraden (W3).")

# ------------------------------------------------------------------------------------------------------------------
say("\nTeil E  Bilanz und Urteil")
bilanz = [
    ("E", "M5", "Naht = e^(-2 pi K)/12 mit K = L_0 - 1/3; kein Faktor zwei gegen die Boostperiode 2 pi."),
    ("E/L", "Ecke", "<K> = 0 und K wechselt unter Spiegelung das Vorzeichen: Form der Eckbedingung; Lesart als "
     "Kegelscheitel."),
    ("E/L", "M8", "Uhrgleichung K = -E_6/(3 E_4): Verdrillung quadratisch mit Fokussierungsvorzeichen, Temperatur "
     "linear. S-II für eine Naht bestätigt; S-I nur zur Hälfte."),
    ("E negativ", "Stufe 2", "Torsionssektoren geben 1, 1/2, 1/2, 0; M2 negativ; Zustandszahl nicht aus der Naht."),
    ("E", "Schwaches Zentrum", "e^(2 pi i 3Y) = Zentrum der schwachen SU(2) = anomale Involution von E8; ihr doppelt "
     "verdrillter Sektor verschwindet genau bei tau = i."),
    ("E/L", "Hierarchie", "Exponent = Gesamtkrümmung/e^2: Gewicht der Naht als klassische Konfiguration."),
    ("E/C", "M9", "P_IF = Wahl der Spiegelung; unter P_TP getroffen."),
    ("E/L", "W4", "Helizität der 16 aus dem Spinbündel der Naht."),
    ("O", "Dynamik", "Wechselwirkung zwischen Geraden, M3, W3 über mehrere Geraden, Kopplungswert, G."),
]
for st, was, text in bilanz:
    say(f"   [{st}] {was}: {text}")
urteil = ("Stufe 2 ist in der angekündigten Form negativ: Die Hierarchie e^(-2/alpha) ist kein Verhältnis von "
          "Determinanten der Naht bei fester Kopplung, die Markenmatrizen geben die Nahtantwort nicht, und die "
          "Zustandszahl des Horizonts folgt nicht aus der einzelnen Naht. Positiv und exakt sind fünf Dinge. "
          "(1) Die Naht am symmetrischen Punkt ist ein Boostzustand e^(-2 pi K)/12 mit K = L_0 - 1/3 und "
          "mittlerer Boostladung null; der vermeintliche Widerspruch 2 pi gegen 4 pi ist eine Erzeugerwahl. "
          "(2) Die Naht hat eine eigene nichtlineare Uhrgleichung K(tau) = -E_6/(3 E_4), mit dem Vorzeichen der "
          "Fokussierung in der Verdrillung: Nichtlinearität ohne nicht-Gaußschen Zustand. (3) Der Exponent der "
          "Hierarchie ist die Gesamtkrümmung der Naht geteilt durch e^2. (4) Das Verklebungsbit ist die Wahl der "
          "Spiegelung und unter P_TP keine eigene Prämisse mehr. (5) Die Zwei-Torsion der Hyperladung ist das Zentrum "
          "der schwachen SU(2) und die anomale Involution von E8; ihr doppelt verdrillter Sektor verschwindet genau "
          "am symmetrischen Punkt, die einfach verdrillten tragen exakt das halbe Gewicht. Zusammen mit z167 heißt das: Alles, was an einer "
          "einzelnen Naht rechenbar ist, ist jetzt gerechnet und stimmt in der Form; es liefert Einheiten, "
          "Verhältnisse und Vorzeichen. Die drei offenen Zahlen (gemeinsame Kopplung, G, Hierarchie) und die "
          "Wechselwirkung liegen nachweislich nicht in der einzelnen Naht, sondern in der Kopplung zwischen "
          "Geraden. Keine volle Lösung, nichts befördert.")
say("\nUrteil\n   " + urteil)
print(f"\n   (Laufzeit {time.time() - T0:.1f} s, nicht Teil der Ausgabe)", flush=True)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z168_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
(HERE / "ausgabe_z168.txt").write_text("\n".join(OUT) + "\n")
