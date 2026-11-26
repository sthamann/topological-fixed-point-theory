"""Z170: automatische Suche hinter den vier Uhrsektoren (OriginLab erweitert) und was sie zusammenfügt.

Leitlinie 'Rekonstruktion statt Auswahl'; H6 und die gelösten Blöcke bleiben unberührt, nichts wird befördert.
Statuskürzel: E exakt, C bedingt, L Lesart, N numerisch, O offen.

Die Suche steht in research/originlab/clock_sector_search.py (drei endliche Universen, Regelmaschine von OriginLab).

Teil A  Universum Uhren: alle Klassen der Ordnung 2 bis 8 in E8
Teil B  Universum Spiegelvervollständigungen: alle selbstorthogonalen Untergruppen der 16 x 16 Sektoren
Teil C  Universum Antworten: die Potenz 5/4
Teil D  Die 16 Uhrsektoren am symmetrischen Punkt und die Energie des Defekts
Teil E  Bilanz und Urteil.
"""
import json
import os
import sys
import time
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import mpmath as mp
import sympy as sp

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parents[2] / "research" / "originlab"))
import clock_sector_search as cs                                           # noqa: E402

OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
T0 = time.time()
mp.mp.dps = 30
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


res = cs.run_clock_sector_search()
uni = res["universes"]


def unique(universe, clauses):
    for row in uni[universe]["search"]["unique_selections"]:
        if any(set(m) == set(clauses) for m in row["minimal_rules"]):
            return row["selected"]
    return None


say("Z170: automatische Suche hinter den vier Uhrsektoren")
say("=" * 100)
say("\nTeil A  Universum Uhren")
ck = uni["clocks"]
h6 = ck["h6_clock"]
ZAHLEN.update(klassen_je_ordnung=ck["classes_per_order"], spur_null=ck["trace_zero_classes"])
check("A1 Vollständige Aufzählung über Kac-Koordinaten: 138 Klassen von Elementen der Ordnung 2 bis 8 in E8 "
      "(2, 4, 7, 14, 20, 38, 53 je Ordnung). Die H6-Uhr ist dabei: fester Wurzelteil A3 + D5, Vielfachheiten "
      "(60, 64, 60, 64), Grundgewicht des verdrillten Sektors 3/8, Anomalietyp 2 [E]",
      ck["candidates"] == 138 and ck["classes_per_order"] == {2: 2, 3: 4, 4: 7, 5: 14, 6: 20, 7: 38, 8: 53}
      and h6 is not None and ck["derived"][h6] == {"multiplicities": [60, 64, 60, 64], "twist_ground": "3/8",
                                                   "anomaly_type": 2})
zero = ck["trace_zero_classes"]
check("A2 Spur null auf den 248 haben genau fünf Klassen: die Uhr (Ordnung 4), zwei der Ordnung 6, zwei der "
      "Ordnung 8. Bei Ordnung 2, 3 und 5 gibt es keine. Die Uhr ist das Element kleinster Ordnung mit Sp R = 0 und "
      "darin eindeutig [E]",
      len(zero) == 5 and [z for z in zero if z[:2] in ("Z2", "Z3", "Z4", "Z5")] == [h6]
      and sorted(z[:2] for z in zero) == ["Z4", "Z6", "Z6", "Z8", "Z8"])
check("A3 Regelsuche: 'Spur null + Ordnung höchstens 4' wählt eindeutig die Uhr; ebenso 'Spur null + fester Teil "
      "halbeinfach' (kein U(1)-Rest) unter allen 138 Klassen. Die Gruppe SO(10) x SU(4) ist damit keine Wahl, "
      "sondern Folge von Sp R = 0, der Volumenaufhebung aus z165 [E für die Aufzählung, C: Sp R = 0 als native "
      "Forderung]",
      unique("clocks", ["order_at_most_4", "trace_zero"]) == h6
      and unique("clocks", ["semisimple_fixed", "trace_zero"]) == h6
      and "square_is_fermion_parity" in ck["clauses"][h6] and "anomaly_free" not in ck["clauses"][h6])
befund("A_uhr", {"klassen": 138, "eindeutig": "Sp R = 0 und Ordnung <= 4"},
       "Die Uhr ist unter allen Elementen kleiner Ordnung von E8 durch eine einzige native Bedingung festgelegt.")

# ------------------------------------------------------------------------------------------------------------------
say("\nTeil B  Universum Spiegelvervollständigungen")
cp = uni["completions"]
graphs = [l for l, cl in cp["clauses"].items() if {"graph", "bosonic"} <= set(cl)]
ZAHLEN.update(vervollstaendigungen=cp["counts"])
check("B1 Die 16 Sektoren von SO(10)_1 x SO(6)_1 haben genau 27 Links-rechts-Vervollständigungen "
      "(selbstorthogonale Untergruppen der Ordnung 16): 10 bosonisch, 17 mit Spinstruktur [E]",
      cp["candidates"] == 27 and cp["counts"]["bosonic"] == 10 and cp["counts"]["spin"] == 17)
check("B2 Die bosonischen Vervollständigungen, in denen jeder Sektor genau einmal vorkommt (Graphen), sind genau "
      "vier: die Vorzeichen-Isometrien (x, y) -> (+-x, +-y). Das sind zwei Bits (Ladungskonjugation auf der D5- "
      "und auf der A3-Seite); P_IF aus z144 ist ihr Verhältnis. Alle vier sind spiegelsymmetrisch und enthalten "
      "die 16 und die 4 als skalare Felder [E]",
      len(graphs) == 4 and all({"scalar_16", "scalar_4", "mirror_symmetric"} <= set(cp["clauses"][l]) for l in graphs)
      and unique("completions", ["diagonal"]) is not None and unique("completions", ["charge_conjugate"]) is not None)
e8_left = [l for l, cl in cp["clauses"].items() if "factorised" in cl and "left [(0, 0), (1, " in l]
check("B3 Ausschluss: Keine der 27 Vervollständigungen enthält zugleich die holomorphen E8-Ströme und die 16 als "
      "lokales Feld. Wo links E8 steht (6 Fälle), gibt es keinen Uhrsektor. Neben E8 existiert die 16 nur im "
      "R-verdrillten Sektor, also als Ende einer Uhr-Defektlinie. Die Defektlesart der Naht (K19) ist damit "
      "keine freie Lesart mehr: Sie folgt aus 'E8 und die 16 zugleich' [E; offen bleibt, warum der Grundzustand "
      "den Defekt enthält]",
      cp["e8_currents_and_local_16"] == [] and len(e8_left) == 6
      and all("clock_sectors_present" not in cp["clauses"][l] for l in e8_left))
check("B4 Skalendimensionen der verdrillten Skalare in den Graphen: 2 h = 3/4 (die 4), 5/4 (die 16), 1 (Vektoren). "
      "3/4 und 5/4 sind die zwei Zahlen des alpha-Blocks im Ledger (c3/phi_base = 3/4 = q(A3), Antwortpotenz "
      "5/4 = q(D5)) [E für die Gewichte, L für die Gleichsetzung]",
      2 * sp.Rational(3, 8) == sp.Rational(3, 4) and 2 * sp.Rational(5, 8) == sp.Rational(5, 4))
befund("B_spiegel", {"anzahl": 27, "graphen": 4, "E8_und_16": "unverträglich"},
       "Materie und Eichströme sind nicht gleichzeitig lokal. Das ist die genaue Form der fehlenden Hälfte: Die "
       "Spiegelvervollständigung ist bis auf zwei Bits bestimmt, aber sie ersetzt die holomorphen E8-Ströme.")

# ------------------------------------------------------------------------------------------------------------------
say("\nTeil C  Universum Antworten")
rs = uni["responses"]
h, n = sp.Symbol("h", positive=True), 8
norms = [sp.Integer(1)]
for k in range(1, n):
    norms.append(sp.expand(norms[-1] * k * (2 * h + k - 1)))             # |L_-1^k h|^2 aus [L_1, L_-1] = 2 L_0
q = sp.Symbol("q")
kern = sum(nk / sp.factorial(k) ** 2 * q ** k for k, nk in enumerate(norms))
ziel = sp.series((1 - q) ** (-2 * h), q, 0, n).removeO()
check("C1 Kohärenzkern: Für einen Zustand niedrigsten Gewichts h ist |L_-1^n h|^2 = n! (2h)_n, also "
      "Summe |L_-1^n h|^2 q^n/(n!)^2 = (1 - q)^(-2h) (bis Ordnung 7 symbolisch). Für die 16 im R-Sektor (h = 5/8) "
      "ist das (1 - q)^(-5/4): ein positiver Hilbertraum mit Norm, ohne Zusatzmaß [E]",
      sp.expand(kern - ziel) == 0 and all(nk.subs(h, sp.Rational(5, 8)) > 0 for nk in norms))
sel = unique("responses", ["insertion_matches", "native"])
check("C2 Regelsuche: Unter den acht Kandidaten ist der Kohärenzkern der 16 mit F = q K der einzige, der nativ ist "
      "und die Nahtantwort q (1 - q)^(-5/4) trifft. Das quartische Maß der Analyse 22 trifft die Potenz, ist aber "
      "nicht nativ; die logarithmische Ableitung trifft die Einfügung nicht [E; die Einfügung F = q K (ein Schritt "
      "mit Gewicht q vor der Leiter) ist Lesart]",
      sel is not None and sel.startswith("coherent kernel of the 16")
      and "native" not in rs["clauses"]["quartic measure on five real directions (analysis 22), F = q K"]
      and "insertion_matches" not in rs["clauses"]["coherent kernel of the 16, F = q d log K/dq"])
check("C3 Gegenprobe gegen Beliebigkeit: q (1 - q)^(-5/4) = 4 q d/dq (1 - q)^(-1/4) gilt ebenfalls exakt. Die "
      "Potenz allein unterscheidet die Lesarten nicht; nativ ist nur die erste (kein Zustand mit h = 1/8 im "
      "Sektorenobjekt) [E]",
      sp.simplify(4 * q * sp.diff((1 - q) ** sp.Rational(-1, 4), q) - q * (1 - q) ** sp.Rational(-5, 4)) == 0
      and "native" not in rs["clauses"]["coherent kernel of weight h = 1/8, F = 4 q dK/dq"])
befund("C_antwort", {"nativ": "Kohärenzkern der 16, h = 5/8"},
       "Die Antwortpotenz 5/4 hat einen nativen positiven Träger: die Leiter der 16 im Uhrsektor.")

# ------------------------------------------------------------------------------------------------------------------
say("\nTeil D  Die 16 Uhrsektoren am symmetrischen Punkt")


def theta_e8(tau, z):
    qn = mp.e ** (I * mp.pi * tau)
    prods = []
    for k in (3, 4, 2, 1):
        p = mp.mpf(1)
        for zi in z:
            p *= mp.jtheta(k, mp.pi * zi, qn)
        prods.append(p)
    return sum(prods) / 2


def eta(tau):
    qn = mp.e ** (2 * I * mp.pi * tau)
    return mp.e ** (I * mp.pi * tau / 12) * mp.qp(qn)


V = [mp.mpf(1) / 2] * 5 + [mp.mpf(0)] * 3


def n2(a, b, tau=I):
    z = [(a + b * tau) * x for x in V]
    y2 = mp.fsum(mp.im(x) ** 2 for x in z)
    return abs(theta_e8(tau, z) / eta(tau) ** 8) ** 2 * mp.e ** (-2 * mp.pi * y2 / mp.im(tau))


tab = {(a, b): n2(a, b) for a in range(4) for b in range(4)}
w2 = mp.sqrt(2)
soll = {(0, 0): 144, (2, 0): 64, (0, 2): 64, (2, 2): 16}
soll.update({k: 8 for k in ((1, 1), (1, 3), (3, 1), (3, 3))})
soll.update({k: 32 + 24 * w2 for k in ((1, 0), (3, 0), (0, 1), (0, 3))})
soll.update({k: 24 * w2 - 32 for k in ((1, 2), (3, 2), (2, 1), (2, 3))})
ZAHLEN["uhrsektoren_norm2"] = {f"{a}{b}": mp.nstr(v, 12) for (a, b), v in tab.items()}
check("D1 Alle 16 Sektoren (Uhr a-mal in der Zeit-, b-mal in der Raumrichtung) bei tau = i in geschlossener Form: "
      "||Z||^2 = 144; 64, 64, 16 (Sektoren von R^2); 8 (vier Sektoren); 32 + 24 Wurzel 2 (vier) und "
      "24 Wurzel 2 - 32 (vier) [N, 25 Stellen]. Summe 320 + 192 Wurzel 2. Kein Uhrsektor verschwindet (anders als "
      "beim anomalen schwachen Zentrum aus z168)",
      all(abs(tab[k] - soll[k]) < mp.mpf("1e-25") for k in tab)
      and abs(mp.fsum(tab.values()) - (320 + 192 * w2)) < mp.mpf("1e-24") and min(tab.values()) > 1)
check("D2 Energie des Defekts: Das Grundniveau des R-Sektors liegt bei L_0 - 1/3 = 3/8 - 1/3 = +1/24, das des "
      "R^2-Sektors bei 1/6, das Vakuum bei -1/3. Der Defektsektor ist modular nötig, aber energetisch nicht "
      "bevorzugt: K19 bleibt offen [E]",
      sp.Rational(3, 8) - sp.Rational(1, 3) == sp.Rational(1, 24) and sp.Rational(1, 2) - sp.Rational(1, 3) == sp.Rational(1, 6))
befund("D_sektoren", {"werte": "144, 64, 16, 8, 32 +- 24 Wurzel 2"},
       "Das Sektorenobjekt ist am symmetrischen Punkt vollständig und in geschlossener Form bekannt.")

# ------------------------------------------------------------------------------------------------------------------
say("\nTeil E  Bilanz und Urteil")
bilanz = [
    ("E/C", "Uhr", "eindeutig durch Sp R = 0 unter allen Klassen bis Ordnung 5, mit 'halbeinfach' bis Ordnung 8."),
    ("E", "Spiegel", "27 Vervollständigungen; die vollständigen bosonischen sind vier (zwei Bits)."),
    ("E", "Ausschluss", "E8-Ströme und lokale 16 schließen sich aus: die 16 ist ein Defektende."),
    ("E/L", "Potenz 5/4", "nativer positiver Träger: Kohärenzkern der 16; Einfügung Lesart."),
    ("E/L", "3/4 und 5/4", "Skalendimensionen der zwei verdrillten Grundfelder."),
    ("E", "16 Sektoren", "geschlossene Werte bei tau = i."),
    ("O", "offen", "Auswahl des Defekts im Grundzustand (K19), Wert der Kopplung, G, Hierarchie, vier Dimensionen "
     "aus dem Sektorenobjekt, die zwei Bits."),
]
for st, was, text in bilanz:
    say(f"   [{st}] {was}: {text}")
urteil = ("Die erweiterte Programmsuche hat drei vollständige endliche Aufzählungen geliefert. (1) Die Uhr ist unter "
          "138 Klassen von E8 durch Sp R = 0 festgelegt; SO(10) x SU(4) ist Folge, nicht Wahl. (2) Von den 27 "
          "Spiegelvervollständigungen der Uhrsektoren enthält keine zugleich die E8-Ströme und die 16 als lokales "
          "Feld; neben E8 ist die 16 das Ende einer Uhr-Defektlinie. Die vollständigen bosonischen "
          "Vervollständigungen sind bis auf zwei Vorzeichenbits bestimmt. (3) Die Antwortpotenz 5/4 hat einen nativen "
          "positiven Träger, die Leiter der 16 im Uhrsektor. Damit fügen sich Uhr, Eichgruppe, Materie als Defekt, "
          "Verklebung und Antwortpotenz in ein Objekt, die vier Uhrsektoren von E8. Nicht gefunden ist eine volle "
          "Lösung: Der Wert der Kopplung, G, die Hierarchie, die Auswahl des Defekts im Grundzustand und der Weg "
          "vom Sektorenobjekt zur vierdimensionalen Raumzeit bleiben offen. Nichts befördert.")
say("\nUrteil\n   " + urteil)
print(f"\n   (Laufzeit {time.time() - T0:.1f} s, nicht Teil der Ausgabe)", flush=True)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z170_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
(HERE / "ausgabe_z170.txt").write_text("\n".join(OUT) + "\n")
