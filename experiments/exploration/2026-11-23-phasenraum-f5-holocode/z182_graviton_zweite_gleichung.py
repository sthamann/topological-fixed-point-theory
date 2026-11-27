"""Z182: F4 und F5, Graviton als kollektive Anregung und die zweite Gravitationsgleichung. Ziele: P_curv, P_G, Gravitonklausel.

Leitlinie 'Rekonstruktion statt Auswahl'; H6 und die gelösten Bausteine bleiben unberührt, nichts wird befördert.
Statuskürzel: E exakt, N numerisch, C bedingt, L Lesart oder Literatur, O offen.
Bestand: z83 (unter I_step trägt die Standardform einer Naht jede masselose Helizität h = -2..2 genau einmal, als
Tiefstgewichtsdarstellung von su(2,2) bei E = |h| + 1); die erste Gravitationsrelation eps^2 G/a^2 = 4 c_3^3 e^(-2/alpha)
mit c_3 = 1/(8 pi); P_G ist der eine freie Zahlenwert der Bilanz.

Fragen. F4: Entsteht ein masseloser Helizität-2-Pol als kollektive Anregung mehrerer Nähte, mit positiven Freiheitsgraden
und universeller Kopplung, im Einklang mit Weinberg-Witten? F5: Gibt dasselbe normierte Spektrum eine zweite Gleichung, die
G/a^2 und eps^2 trennt, ohne Zahlenjagd?

Teil A  Charaktere von su(2,2): Produkte masseloser Darstellungen (Zweinahtzustände) exakt zerlegt.
Teil B  Weinberg-Witten in Charakterform und die Folgen für den Pol.
Teil C  Spektrum der Naht: keine induzierte Einstein-Wirkung bei festem tau; Entropie des Zustands von A*.
Teil D  Zweite Gleichung: Buchführung, bedingte Form, Mehrdeutigkeit.
Teil E  Bilanz.
"""
import json
import time
from fractions import Fraction as Fr
from pathlib import Path

import mpmath as mp

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
T_START = time.time()
mp.mp.dps = 40


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val):
    BEFUNDE[name] = bool(val)
    say(f"   [Befund: {'ja' if val else 'nein'}] {name}")


# =================================================================== Teil A
say("Teil A: Produkte masseloser Darstellungen von su(2,2)")
ORD = 28                                                                   # Potenzreihen in t = s^(1/2) bis t^ORD
U, V = Fr(2), Fr(3)                                                        # x = U^2, y = V^2 (exakte rationale Probe)
U2, V2 = Fr(5, 3), Fr(7, 2)                                                # zweite Probe


def reihe(c0, grad=0):
    out = [Fr(0)] * (ORD + 1)
    if grad <= ORD:
        out[grad] = Fr(c0)
    return out


def mul(a, b):
    out = [Fr(0)] * (ORD + 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b[:ORD + 1 - i]):
                if y:
                    out[i + j] += x * y
    return out


def add(*rs):
    return [sum(col) for col in zip(*rs)]


def skal(c, a):
    return [c * x for x in a]


def chi(j2, u):
    """SU(2)-Charakter zum Spin j = j2/2 bei x = u^2: Summe u^(2m)."""
    if j2 < 0:
        return Fr(0)
    return sum(u ** (j2 - 2 * k) for k in range(j2 + 1))


def geo(w):
    return [Fr(w) ** (n // 2) if n % 2 == 0 else Fr(0) for n in range(ORD + 1)]   # 1/(1 - s w), s = t^2


def impuls(u, v):
    out = reihe(1)
    for w in (u * v, u / v, v / u, 1 / (u * v)):
        out = mul(out, geo(w))
    return out


def masselos(j2, links, u, v, P):
    """Masselose Darstellung der Helizität j = j2/2: (j + 1; j, 0) für links, (j + 1; 0, j) für rechts."""
    a, b = (u, v) if links else (v, u)
    kopf = add(reihe(chi(j2, a), j2 + 2), reihe(-chi(j2 - 1, a) * chi(1, b), j2 + 4), reihe(chi(j2 - 2, a), j2 + 6))
    if j2 == 0:
        kopf = add(reihe(1, 2), reihe(-1, 6))
    return mul(kopf, P)


def strom(a2, b2, u, v, P):
    """Erhaltener Strom (j1 + j2 + 2; j1, j2), j1 = a2/2, j2 = b2/2 > 0."""
    kopf = add(reihe(chi(a2, u) * chi(b2, v), a2 + b2 + 4), reihe(-chi(a2 - 1, u) * chi(b2 - 1, v), a2 + b2 + 6))
    return mul(kopf, P)


def lang(delta2, a2, b2, u, v, P):
    return mul(reihe(chi(a2, u) * chi(b2, v), delta2), P)


zerlegung = {}
for name, (u, v) in (("Probe 1", (U, V)), ("Probe 2", (U2, V2))):
    P = impuls(u, v)
    # Skalar (x) Skalar = (2; 0, 0) + Summe der Ströme (l + 2; l/2, l/2), l >= 1
    sk = mul(masselos(0, True, u, v, P), masselos(0, True, u, v, P))
    rhs = lang(4, 0, 0, u, v, P)
    for l in range(1, ORD):
        rhs = add(rhs, strom(l, l, u, v, P))
    zerlegung[(name, "Skalar")] = sk == rhs
    # (j + 1; j, 0) (x) (j + 1; 0, j) = Summe über k der Ströme (2j + 2 + k; j + k/2, j + k/2)
    for j2 in (1, 2, 3, 4):
        lhs = mul(masselos(j2, True, u, v, P), masselos(j2, False, u, v, P))
        rhs = reihe(0)
        for k in range(0, ORD):
            rhs = add(rhs, strom(j2 + k, j2 + k, u, v, P))
        zerlegung[(name, f"Helizität +{Fr(j2, 2)} mal -{Fr(j2, 2)}")] = lhs == rhs
zaehl = masselos(0, True, Fr(1), Fr(1), impuls(Fr(1), Fr(1)))
ref = [0, 1, 4, 9, 16, 25, 36]                                             # m^2 Zustände auf der Schale E = m (z83)
ZAHLEN["A1"] = {"Zerlegungen bis s^14 (zwei rationale Proben)": {f"{k[0]}: {k[1]}": v for k, v in zerlegung.items()},
                "Skalar: Zustände je Schale": [int(zaehl[2 * m]) for m in range(7)]}
check("A1 [E, Potenzreihen bis s^14, zwei rationale Proben für (x, y)] Charaktere: Die masselose Darstellung der "
      "Helizität j ist (j + 1; j, 0); der Skalar zählt m^2 Zustände auf der Schale E = m wie in z83. Die Produkte "
      "zerfallen exakt: Skalar mal Skalar = (2; 0, 0) + Summe der erhaltenen Ströme (l + 2; l/2, l/2); Helizität +j mal "
      "-j = Summe über k >= 0 der erhaltenen Ströme (2j + 2 + k; j + k/2, j + k/2), geprüft für j = 1/2, 1, 3/2, 2. "
      "Zweinahtzustände sind also erhaltene Ströme, keine masselosen Teilchen",
      all(zerlegung.values()) and [int(zaehl[2 * m]) for m in range(7)] == ref)

# gleiche Händigkeit: kleinste Energie
tief = {}
for a2 in range(0, 5):
    for b2 in range(a2, 5):
        P = impuls(U, V)
        pr = mul(masselos(a2, True, U, V, P), masselos(b2, True, U, V, P))
        erste = next(n for n, x in enumerate(pr) if x)
        tief[(a2, b2)] = (Fr(erste, 2), Fr(a2 + b2, 2) + 1)                  # (kleinste Energie, Energie des masselosen Spins j + j')
ZAHLEN["A2"] = {f"j = {Fr(a, 2)}, j' = {Fr(b, 2)}": [str(v[0]), str(v[1])] for (a, b), v in tief.items()}
check("A2 [E] Kein masseloser Pol aus zwei Nähten: Das Produkt zweier masseloser Darstellungen der Helizitäten j und j' "
      "beginnt bei der Energie j + j' + 2. Eine masselose Darstellung des Spins J <= j + j' läge bei J + 1 <= j + j' + 1, "
      "also darunter. Im Produkt kommt daher keine masselose Darstellung vor, für kein Paar von Helizitäten 0..2 (alle 15 "
      "Paare geprüft). Insbesondere enthält Helizität +1 mal +1 bei der Energie 4 die lange Darstellung (4; 2, 0) und "
      "nicht das Graviton (3; 2, 0): Der Spin-2-Kanal des Viererkerns ist ein Kontinuum mit Spektraldichte ~ mu^4 und hat "
      "keinen Pol 1/p^2", all(v[0] == v[1] + 1 for v in tief.values()) and len(tief) == 15)

# =================================================================== Teil B
say("\nTeil B: Weinberg-Witten in Charakterform")
ww = {}
for j2 in range(0, 5):
    spins = [j2 + k for k in range(0, 6)] if j2 else [k for k in range(1, 6)]   # Spin des Stroms = 2 j1 = j2 + k
    ww[str(Fr(j2, 2))] = {"kleinster Strom (Spin)": min(spins), "Ladungsstrom (Spin 1)": 1 in spins,
                          "Energie-Impuls-Tensor (Spin 2)": 2 in spins}
ZAHLEN["B1"] = ww
check("B1 [E] Weinberg-Witten als Darstellungssatz: Die Bilinearen eines masselosen Feldes der Helizität j sind die Ströme "
      "vom Spin >= 2j (A1). Ein erhaltener Ladungsstrom (3; 1/2, 1/2) kommt nur für j <= 1/2 vor, ein erhaltener "
      "Energie-Impuls-Tensor (4; 1, 1) nur für j <= 1. Für j = 2 beginnt die Reihe beim Spin 4, Energie 6. Das Graviton der "
      "Standardform (h = +-2, z83) trägt also keinen bilinearen kovarianten Energie-Impuls-Tensor, und sein Energie-Impuls "
      "ist kein Integral eines lokalen Stroms der Nähte",
      ww["2"] == {"kleinster Strom (Spin)": 4, "Ladungsstrom (Spin 1)": False, "Energie-Impuls-Tensor (Spin 2)": False}
      and ww["1"]["Energie-Impuls-Tensor (Spin 2)"] and not ww["1"]["Ladungsstrom (Spin 1)"]
      and ww["1/2"]["Ladungsstrom (Spin 1)"] and not ww["3/2"]["Energie-Impuls-Tensor (Spin 2)"])
befund("Kriterium 'echter masseloser Pol': als kollektive Anregung mehrerer Nähte (Pol im Viererkern)", False)
befund("Kriterium 'echter masseloser Pol': als elementarer Zustand h = +-2 der Standardform einer Naht (bedingt unter "
       "I_step, z83)", True)
befund("Kriterium 'positive Freiheitsgrade': die Darstellung (3; 2, 0) ist unitär", True)
befund("Kriterium 'universelle Kopplung': aus einem lokalen Energie-Impuls-Tensor der Nähte hergeleitet", False)

# =================================================================== Teil C
say("\nTeil C: Spektrum der Naht")


def chi_e8(y):
    q = mp.exp(-mp.pi * y)
    eta = mp.exp(-2 * mp.pi * y / 24) * mp.qp(mp.exp(-2 * mp.pi * y))
    return sum(mp.jtheta(k, 0, q) ** 8 for k in (2, 3, 4)) / (2 * eta ** 8)


def chi_diag(y):
    q = mp.exp(-mp.pi * y)
    eta = mp.exp(-2 * mp.pi * y / 24) * mp.qp(mp.exp(-2 * mp.pi * y))
    t2, t3, t4 = (mp.jtheta(k, 0, q) for k in (2, 3, 4))
    g = lambda n: (t2 ** (2 * n) + t3 ** (2 * n) + t4 ** (2 * n)) / (2 * eta ** (2 * n))
    return g(5) * g(3)


rest = {}
for y in (mp.mpf('0.15'), mp.mpf('0.25'), mp.mpf('0.4')):
    w = mp.exp(-2 * mp.pi / y)
    rest[mp.nstr(y, 3)] = (abs(mp.log(chi_e8(y)) - 2 * mp.pi / (3 * y) - mp.log(1 + 248 * w + 4124 * w ** 2 + 34752 * w ** 3)),
                           abs((mp.log(chi_diag(y)) - 4 * mp.pi / (3 * y)) / (32 * mp.exp(-3 * mp.pi / (2 * y))) - 1), w ** 4,
                           mp.exp(-mp.pi / (2 * y)))
ZAHLEN["C1"] = {k: [mp.nstr(x, 4) for x in v] for k, v in rest.items()}
check("C1 [N, 40 Stellen] Keine induzierte Einstein-Wirkung aus dem Spektrum einer Naht: Für kleine y (hohe Temperatur, "
      "tau = i y) ist log Z(i y) = 2 pi/(3 y) + log(1 + 248 e^(-2 pi/y) + ...), für das geeichte Paar 4 pi/(3 y) + "
      "32 e^(-3 pi/(2 y)) + ... (die 32 Zustände der so(6)-Spinorsektoren mit h + h-quer = 3/4). Hinter dem führenden Term (Volumen, Cardy) steht kein Potenzterm: kein y^0, kein "
      "log y, kein y^1. Die modulare Invarianz streicht alle Wärmekernkoeffizienten, aus denen ein Term vom Einstein-Typ "
      "entstehen könnte. Aus dem Spektrum einer Naht bei festem tau ist G nicht induzierbar",
      all(v[0] < max(3e5 * v[2], mp.mpf(10) ** -36) and v[1] < 10 * v[3] for v in rest.values()))

# Entropie des Zustands von A*: S = log Z + 2 pi <L_0 - c/24>, und der zweite Term ist null
def koeff_e8(nmax):
    """q-Entwicklung von E_4/eta^8 = q^(-1/3) (1 + 248 q + ...) mit ganzzahliger Arithmetik."""
    sig3 = lambda n: sum(d ** 3 for d in range(1, n + 1) if n % d == 0)
    e4 = [1] + [240 * sig3(n) for n in range(1, nmax + 1)]
    inv = [1] + [0] * nmax                                                  # 1/prod (1 - q^n)^8
    for n in range(1, nmax + 1):
        for _ in range(8):
            for k in range(n, nmax + 1):
                inv[k] += inv[k - n]
    return [sum(e4[i] * inv[n - i] for i in range(n + 1)) for n in range(nmax + 1)]


a_n = koeff_e8(60)
Z = sum(a * mp.exp(-2 * mp.pi * (n - mp.mpf(1) / 3)) for n, a in enumerate(a_n))
p_n = [mp.exp(-2 * mp.pi * (n - mp.mpf(1) / 3)) / Z for n in range(len(a_n))]
S_naht = -sum(a * p * mp.log(p) for a, p in zip(a_n, p_n))
E_naht = sum(a * p * (n - mp.mpf(1) / 3) for n, (a, p) in enumerate(zip(a_n, p_n)))
S_paar = mp.log(80 + 48 * mp.sqrt(2))
entropien = {"eine Naht, log 12": mp.log(12), "Paar E8 (x) E8-quer je Naht": mp.log(144) / 2,
             "geeichtes Paar je Naht, (1/2) log(80 + 48 sqrt(2))": S_paar / 2,
             "eine Naht minus topologischer Term log 4": mp.log(12) - mp.log(4)}
ZAHLEN["C2"] = {"Koeffizienten": a_n[:5], "Z(i)": mp.nstr(Z, 25), "S der Naht": mp.nstr(S_naht, 25),
                "<L_0 - 1/3>": mp.nstr(E_naht, 5), "Kandidaten je Naht": {k: mp.nstr(v, 15) for k, v in entropien.items()}}
check("C2 [E/N] Entropie des Zustands von A*: Weil die Bilanz <L_0 - c/24> bei tau = i verschwindet, ist die "
      "von-Neumann-Entropie des Zustands genau log Z. Für die Naht: S = log 12 (aus den ersten 60 Niveaus der E8-Reihe "
      f"1 + 248 q + 4124 q^2 + ...: S = {mp.nstr(S_naht, 15)}, log 12 = {mp.nstr(mp.log(12), 15)}). Für das geeichte Paar "
      "ist sie log(80 + 48 sqrt(2)), für E8 (x) E8-quer 2 log 12. Das ist die lokale Antwort desselben Spektrums, dessen "
      "globale Antwort Z ist, ohne freien Parameter",
      a_n[:4] == [1, 248, 4124, 34752] and abs(Z - 12) < mp.mpf(10) ** -30 and abs(S_naht - mp.log(12)) < mp.mpf(10) ** -30
      and abs(E_naht) < mp.mpf(10) ** -30)

# =================================================================== Teil D
say("\nTeil D: Zweite Gleichung")
alpha_inv = mp.mpf('137.035999')
c3 = 1 / (8 * mp.pi)
X = 4 * c3 ** 3 * mp.exp(-2 * alpha_inv)
G_Lambda = 12 * X
S_dS = 3 * mp.pi / G_Lambda
ZAHLEN["D1"] = {"X = eps^2 G/a^2": mp.nstr(X, 6), "G Lambda = 12 X": mp.nstr(G_Lambda, 6), "S_dS = pi/(4 X)": mp.nstr(S_dS, 6)}
check("D1 [E] Buchführung: Die erste Relation X = eps^2 G/a^2 = 4 c_3^3 e^(-2/alpha) legt das dimensionslose Produkt "
      f"G Lambda = 12 X = {mp.nstr(G_Lambda, 4)} und die de-Sitter-Entropie pi/(4 X) = {mp.nstr(S_dS, 4)} bereits fest "
      "(mit 1/alpha = 137.035999 als Messwert eingesetzt). Offen ist allein, wie sich X auf den lokalen Faktor G/a^2 und "
      "den globalen Faktor eps^2 verteilt. P_G ist genau diese eine Zahl",
      abs(G_Lambda / mp.mpf('2.83e-122') - 1) < mp.mpf('0.02') and abs(S_dS - mp.pi / (4 * X)) / S_dS < mp.mpf(10) ** -30)

tabelle = {}
for k, s1 in entropien.items():
    tabelle[k] = {"s_1": mp.nstr(s1, 10), "G/a^2 = 1/(4 s_1)": mp.nstr(1 / (4 * s1), 10),
                  "eps^2 = 4 s_1 X": mp.nstr(4 * s1 * X, 6), "Nähte auf dem Horizont": mp.nstr(S_dS / s1, 6)}
werte = sorted(mp.nstr(1 / (4 * s1), 8) for s1 in entropien.values())
ZAHLEN["D2"] = tabelle
check("D2 [C, Mehrdeutigkeit] Bedingte zweite Gleichung: Unter der Flächenregel (Entropie = Fläche/4G, eine Naht je Zelle "
      "a^2; Prämisse D nach Jacobson, nicht hergeleitet) ist G/a^2 = 1/(4 s_1) mit der Entropie s_1 je Naht, und "
      "eps^2 = 4 s_1 X. Das Spektrum liefert s_1 ohne freien Parameter, aber nicht eindeutig: Je nachdem, ob die Zelle "
      "eine Naht (log 12), ein geeichtes Paar oder eine Naht mit abgezogenem topologischem Term trägt, ist G/a^2 = "
      f"{werte}. E8 (x) E8-quer gibt je Naht denselben Wert wie die einzelne Naht. Drei physikalisch verschiedene Modelle "
      "erfüllen dieselben Bedingungen; das ist als Mehrdeutigkeit zu führen, nicht als Wert von G",
      len(set(werte)) == 3 and abs(entropien["eine Naht, log 12"] - entropien["Paar E8 (x) E8-quer je Naht"]) < mp.mpf(10) ** -30)

# =================================================================== Teil E
say("\nTeil E: Bilanz")
befund("Gravitonklausel geschlossen (masseloser Pol, positive Freiheitsgrade und universelle Kopplung zusammen)", False)
befund("P_curv hergeleitet", False)
befund("P_G: freier Zahlenwert beseitigt", False)
befund("P_G: unter der Flächenregel auf drei diskrete Werte eingeengt", True)
befund("Komplettlösung", False)

urteil = (
    "F4 und F5 sind gerechnet. F4, negativ mit scharfem Geltungsbereich: Das Produkt zweier masseloser Darstellungen von "
    "su(2,2) enthält nie eine masselose Darstellung; Zweinahtzustände sind erhaltene Ströme, und der Spin-2-Kanal des "
    "Viererkerns ist ein Kontinuum ohne Pol 1/p^2. Ein Graviton als kollektive Anregung mehrerer Nähte gibt es im "
    "konformen Rahmen von I_step also nicht. Es bleibt der elementare Zustand h = +-2 der Standardform einer Naht (z83): "
    "unitär und masselos, aber ohne bilinearen Energie-Impuls-Tensor (die Ströme beginnen beim Spin 4). Das ist "
    "Weinberg-Witten als Darstellungssatz und zugleich der Grund, warum die universelle Kopplung nicht aus einem lokalen "
    "Strom der Nähte folgen kann; sie braucht ein Eichprinzip, also P_curv. F5: Die erste Relation legt G Lambda und die "
    "de-Sitter-Entropie schon fest; P_G ist nur die Verteilung auf G/a^2 und eps^2. Aus dem Spektrum einer Naht bei "
    "festem tau entsteht kein Einstein-Term (die modulare Invarianz streicht alle Potenzkorrekturen). Dasselbe Spektrum "
    "gibt aber die Entropie des Zustands von A* exakt: log 12 je Naht, weil <T> = 0 ist. Unter der Flächenregel wird "
    "daraus G/a^2 = 1/(4 s_1), allerdings mit drei physikalisch verschiedenen Kandidaten für s_1. P_G ist damit von "
    "einer freien reellen Zahl auf eine diskrete Mehrdeutigkeit unter einer benannten Prämisse eingeengt, nicht "
    "beseitigt. Gravitonklausel, P_curv und P_G bleiben offen. Keine Komplettlösung, nichts befördert.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z182_ergebnis.json").write_text(json.dumps({"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil},
                                                    indent=1, ensure_ascii=False, default=str))
(HERE / "ausgabe_z182.txt").write_text("\n".join(OUT) + "\n")
say(f"(Laufzeit {time.time() - T_START:.0f} s)")
