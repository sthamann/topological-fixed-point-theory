"""Z81: Vorprüfung für die Strategien nach H6h40 (Gesamtstand und Strategien, 25.11.).

Drei Fakten, auf denen die neuen Strategien aufbauen, exakt nachgerechnet:
Teil 1  Die Uhr R = Delta^(-i/2) = e^{2 pi i L_0} ist auf dem ganzen Original R = i^(Q_a - Q_b) (Ladungsteil) mal 1
        (Oszillatoren). Daraus R^2 = (-1)^F und R^4 = 1: R ist ein Spin-Z4-Erzeuger (Z4 mit Z4^2 = (-1)^F), dieselbe
        algebraische Art wie die Z4 der Kette mu_4 -> 16 (TFPT-STAND §3).
Teil 2  Die Einteilchendaten des Originals zerfallen in Fermipaare (x, -x), x = m/4 mit ungeradem m: Die normierte Kovarianz
        auf jedem Paar ist genau rho_m = diag(1, e^{-pi m})/(1 + e^{-pi m}); R wirkt auf dem Paar als diag(i^m, i^-m) mit
        Quadrat -1. Das Original ist die zweite Quantisierung der Saatfamilie {rho_m, m ungerade}.
Teil 3  Die vier Marken S^4 = 1 der H6-Naht liegen über genau einem Punkt w = S^4 = 1 der physischen Einheitswindung; die
        Uhr ist die Deckdrehung S -> i S. Die Marken sind eine Z4-Faser über einem Punkt, keine Punkte der Raumzeit.
"""
import json
from pathlib import Path

import numpy as np
import sympy as sp

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val):
    BEFUNDE[name] = bool(val)
    say(f"   [Befund: {'ja' if val else 'nein'}] {name}")


# ------------------------------------------------------------------ Teil 1
say("Teil 1: Die Uhr R auf dem ganzen Original")
E = lambda alpha, Q: sp.Rational(1, 2)*Q**2 + (alpha - sp.Rational(1, 2))*Q
a, b = sp.Rational(1, 4), sp.Rational(3, 4)
ok_R, ok_R2, ok_R4 = True, True, True
for Qa in range(-12, 13):
    for Qb in range(-12, 13):
        for Nosc in (0, 1, 5):
            L0 = E(a, Qa) + E(b, Qb) + Nosc
            R = sp.exp(2*sp.pi*sp.I*L0)
            ok_R &= sp.simplify(R - sp.I**(Qa - Qb)) == 0
            ok_R2 &= sp.simplify(R**2 - (-1)**(Qa + Qb)) == 0
            ok_R4 &= sp.simplify(R**4 - 1) == 0
check("R = e^{2 pi i L_0} = i^(Q_a - Q_b) auf allen Ladungssektoren |Q| <= 12 mit beliebigem Oszillatorinhalt; R^2 = (-1)^(Q_a + Q_b) "
      "= (-1)^F und R^4 = 1: R ist ein Spin-Z4-Erzeuger", ok_R and ok_R2 and ok_R4)
befund("R hat dieselbe algebraische Art wie die Z4 der Kette mu_4 -> 16 (ein Weylfermion mit Z4, Z4^2 = (-1)^F, ist "
       "anomaliefrei in Vielfachen von 16). Ob es dieselbe Z4 ist (Ladungen der 16 unter R), ist die erste Rechnung der "
       "Strategie S1", True)

# ------------------------------------------------------------------ Teil 2
say("\nTeil 2: Fermipaare (x, -x) und die Saatfamilie rho_m")
c = lambda x: 1/(1 + sp.exp(4*sp.pi*x))
ok_pair, ok_Rpair = True, True
for m in (1, 3, 5, 7, 9, 11):
    x = sp.Rational(m, 4)
    tot = sp.simplify(c(x) + c(-x))
    rho = (sp.simplify(c(-x)/tot), sp.simplify(c(x)/tot))
    soll = (1/(1 + sp.exp(-sp.pi*m)), sp.exp(-sp.pi*m)/(1 + sp.exp(-sp.pi*m)))
    ok_pair &= tot == 1 and all(sp.simplify(r - s) == 0 for r, s in zip(rho, soll))
    Rp = (sp.exp(2*sp.pi*sp.I*x), sp.exp(-2*sp.pi*sp.I*x))
    ok_Rpair &= all(sp.simplify(r**2 + 1) == 0 for r in Rp) and sp.simplify(Rp[0] - sp.I**m) == 0
alpha_paar = sorted({float((sp.Rational(m, 4)) % 1) for m in (1, 3, 5, 7)} | {float((-sp.Rational(m, 4)) % 1) for m in (1, 3, 5, 7)})
check("für x = m/4 (m = 1, 3, .., 11) ist c(x) + c(-x) = 1 und die Paarkovarianz genau rho_m = diag(1, e^{-pi m})/(1 + e^{-pi m}); "
      "R wirkt auf dem Paar als diag(i^m, i^-m) mit Quadrat -1", ok_pair and ok_Rpair and alpha_paar == [0.25, 0.75])
befund("Die Einteilchendaten des Originals sind genau die Saatfamilie {rho_m, m ungerade}, je ein Paar aus beiden Quartern; "
       "das quasifreie Original ist ihre zweite Quantisierung. Gerade m (rho_2, rho_4) kommen als Einteilchenpaare nicht vor, "
       "nur als Fock-Ecken (z63). Die Frage nach dem Ursprung verschiebt sich damit auf die Saat rho_1 und die Regel 'alle "
       "ungeraden m' (Strategie S6)", True)

# ------------------------------------------------------------------ Teil 3
say("\nTeil 3: Die Marken über der physischen Einheitswindung")
marken = [sp.exp(sp.I*sp.pi*k/2) for k in range(4)]
ok_w = all(sp.simplify(s**4 - 1) == 0 for s in marken)
ok_deck = all(any(sp.simplify(sp.I*s - t) == 0 for t in marken) for s in marken)
check("die vier Marken S in mu_4 haben dasselbe Bild w = S^4 = 1, und die Uhr S -> i S permutiert sie zyklisch (Deckdrehung "
      "der vierblättrigen Überlagerung S -> S^4, verzweigt nur in 0 und unendlich)", ok_w and ok_deck)
befund("Die Marken sind eine Z4-Faser über genau einem Punkt der physischen Einheitswindung. Zusammen mit z79 (R ist auf den "
       "Oszillatoren trivial, auf der Raumzeit der Relativzustände unsichtbar) stützt das die Lesart: Marken und Uhr sind "
       "innere Struktur (Spin-Z4), nicht Punkte der Himmelssphäre", True)

urteil = ("Drei Grundlagen der Strategien nach H6h40 sind exakt: R = i^(Q_a - Q_b) mit R^2 = (-1)^F (Spin-Z4, gleiche Art wie die "
          "Z4 der 16-Kette); die Einteilchendaten sind die Saatfamilie {rho_m, m ungerade} auf Fermipaaren (x, -x); die vier "
          "Marken liegen über einem Punkt der physischen Windung, die Uhr ist ihre Deckdrehung. Ob R die Z4 der 16 ist und ob "
          "die Saatfamilie aus rho_1 folgt, bleibt zu rechnen.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z81_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
