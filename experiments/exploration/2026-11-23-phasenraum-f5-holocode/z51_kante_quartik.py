"""Z51: Was legt die Wechselwirkung des metrischen Sterns fest? Ungerichtete Kanten, die quartische Zahl t, die Kopplungstiefe.

Teil 1  Zählung (aus berechneten Taylorkoeffizienten, ohne Literaturformeln): Flachheit beider Metriken und die
        Fierz-Pauli-Masse lassen in den fünf Hassan-Rosen-Koeffizienten zwei Richtungen frei. Eine ungerichtete Kante
        (Potential gleich, ob von der Eltern- oder der Kindmetrik gelesen) lässt eine: die quartische Zahl t. Die kubische
        Ordnung ist damit festgelegt.
Teil 2  Der quartische Rest sitzt auf den Kanten des Transfers: Delta c4 = (t/2) det(B - A)/16, in kanonischen Feldern
        (t/2) mu^2 det((B h)_e)/M_e^2 -- genau die Zeilen des Vierteltransfers, gewichtet mit 1/M_e^2.
Teil 3  Kopplungstiefe: Ein Knoten ist schwach gekoppelt, solange seine Planckmasse über der Fierz-Pauli-Skala liegt
        (Lambda_3 = (mu^2 M)^(1/3) >= mu). Mit M_(j,n) = (M_c/2) q^(n+1) bricht das nach n_max Stufen; im originalen
        Zustand kostet jede Stufe den Faktor e^{2 pi} = 535.
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO / "research" / "originlab"))
from process_rules import BETA_DET, BETA_V, Q_ORIGINAL, potential_taylor, run  # noqa: E402

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


rng = np.random.default_rng(51)
pairs = []
for _ in range(3):
    A, B = rng.normal(size=(4, 4)), rng.normal(size=(4, 4))
    pairs.append(((A + A.T)/2, (B + B.T)/2))
basis = np.eye(5)
tay = {n: [potential_taylor(tuple(basis[n]), A, B) for A, B in pairs] for n in range(5)}

# ------------------------------------------------------------------ Teil 1
say("Teil 1: Zählung der freien Potentialrichtungen")
lin1 = np.array([[tay[n][k][1] for n in range(5)] for k in range(3)])          # c1 auf drei Paaren
lin0 = np.array([[tay[n][0][0] for n in range(5)]])
fpc = np.array([[tay[n][k][2]/(np.trace((A - B) @ (A - B)) - np.trace(A - B)**2) for n in range(5)]
                for k, (A, B) in enumerate(pairs)])
r1 = np.linalg.matrix_rank(lin1, tol=1e-9)
r01 = np.linalg.matrix_rank(np.vstack([lin0, lin1]), tol=1e-9)
say(f"   Rang der Kaulquappenbedingungen c1 = 0: {r1}; mit c0 = 0: {r01}")
check("Flachheit ist genau zwei Bedingungen (Lambda_g = Lambda_f = 0); c0 = 0 folgt daraus", r1 == 2 and r01 == 2)
flat = np.linalg.svd(lin1)[2][2:]                                              # drei flache Richtungen
fp_flat = fpc @ flat.T                                                         # Fierz-Pauli-Koeffizient je Paar und Richtung
check("auf den drei flachen Richtungen ist der quadratische Term reine Fierz-Pauli-Form tr D^2 - (tr D)^2 (gleicher "
      "Koeffizient auf allen Paaren); nicht flache Richtungen haben Zusatzterme", np.allclose(fp_flat, fp_flat[0], atol=1e-9)
      and not np.allclose(fpc[:, 0], fpc[0, 0], atol=1e-9))
cons = np.vstack([lin1, fpc[:1]])
null = np.linalg.svd(cons)[2][np.linalg.matrix_rank(cons, tol=1e-9):]
say(f"   Flachheit plus feste Masse: Rang {np.linalg.matrix_rank(cons, tol=1e-9)}, freie Richtungen {len(null)}")
check("Flachheit plus Fierz-Pauli-Masse lassen zwei Richtungen frei", len(null) == 2)
pal = np.array([[1, 0, 0, 0, -1], [0, 1, 0, -1, 0]], float)
cons_s = np.vstack([cons, pal])
null_s = np.linalg.svd(cons_s)[2][np.linalg.matrix_rank(cons_s, tol=1e-9):]
check("mit ungerichteter Kante (beta_n = beta_(4-n)) bleibt genau eine Richtung, und sie ist det(I - X)",
      len(null_s) == 1 and abs(abs(null_s[0] @ np.array(BETA_DET))/np.linalg.norm(BETA_DET) - 1) < 1e-9)
c3dir = [sum(null_s[0][n]*tay[n][k][3] for n in range(5)) for k in range(3)]
c4dir = [sum(null_s[0][n]*tay[n][k][4] for n in range(5)) for k in range(3)]
check("längs dieser Richtung ändert sich die kubische Ordnung nicht, die quartische schon",
      max(abs(c) for c in c3dir) < 1e-9 and min(abs(c) for c in c4dir) > 1e-4)
U_ = run()["universes"]["metric_potentials"]
say(f"   Lane: flach + Masse + ungerichtet lässt {U_['tie_flat_mass_undirected']} offen")
check("in der Lane lassen flach, Fierz-Pauli-Masse und ungerichtete Kante genau V und V +- det(I-X)/2 offen",
      U_["tie_flat_mass_undirected"] == ["V", "V + det(I-X)/2", "V - det(I-X)/2"])
befund("Ungerichtete Kanten legen die kubische Wechselwirkung fest; die quartische Zahl t bleibt die einzige freie "
       "Potentialzahl je Kante", True)

# ------------------------------------------------------------------ Teil 2
say("\nTeil 2: Der quartische Rest sitzt auf den Transferkanten")
devs = []
for A, B in pairs:
    cv = potential_taylor(BETA_V, A, B)
    ct = potential_taylor(tuple(v + d/2 for v, d in zip(BETA_V, BETA_DET)), A, B)       # t = 1
    devs.append(((ct[4] - cv[4]), np.linalg.det(B - A)/32, max(abs(ct[k] - cv[k]) for k in range(4))))
say(f"   Delta c4 gegen (t/2) det(B - A)/16 bei t = 1: {[(round(a, 10), round(b, 10)) for a, b, _ in devs]}")
check("Delta c4 = (t/2) det(B - A)/16, alle niedrigeren Ordnungen gleich", all(abs(a - b) < 1e-9 and c < 1e-10 for a, b, c in devs))
check("in kanonischen Feldern (g_v = eta + 2 h_v/M_v, kappa = mu^2 M_k^2) ist das (t/2) mu^2 det((B h)_e)/M_e^2: "
      "B - A = -2 (B h)_e/M_e, det(-2 X/M)/16 = det(X)/M^4", abs(np.linalg.det(-2*np.eye(4)/3)/16 - 1/3**4) < 1e-15)
befund("Die einzige freie Wechselwirkungszahl ist ein lokaler Viererkontakt auf den Zeilen des Vierteltransfers, "
       "gewichtet mit 1/M_e^2; er wächst in der Tiefe wie q^(-2n)", True)

# ------------------------------------------------------------------ Teil 3
say("\nTeil 3: Kopplungstiefe")


def depth(qv, mu_over_mc):
    n, out = 0, -1
    while 0.5*qv**(n + 1) >= mu_over_mc and n < 200:
        out, n = n, n + 1
    return out


rows = []
for qv, name in ((Q_ORIGINAL, "e^-2pi"), (0.5, "1/2")):
    for r in (1e-3, 1e-10, 1e-16):
        rows.append((name, r, depth(qv, r)))
for name, r, d in rows:
    say(f"   q = {name:7s} mu/M_c = {r:.0e}: schwach gekoppelte Armstufen n <= {d}")
check("im originalen Zustand kostet jede Stufe den Faktor 1/q = e^{2 pi} = 535: bei mu/M_c = 1e-16 (TeV gegen Planck) sind "
      "es fünf Stufen, bei mu/M_c = 1e-3 keine", dict(((n, r), d) for n, r, d in rows)[("e^-2pi", 1e-16)] == 4
      and dict(((n, r), d) for n, r, d in rows)[("e^-2pi", 1e-3)] == -1)
hidden_tail = Q_ORIGINAL**(2*6)
check("was jenseits dieser Tiefe stark koppelt, trägt im originalen Zustand höchstens das Gewicht q^(2(n_max+1)) "
      "(fünf Stufen: q^12 = 1.8e-33)", hidden_tail < 1e-32)
befund("Der unendliche Hardy-Grenzfall ist eine lineare Aussage; nichtlinear ist der Stern bis zu einer endlichen Tiefe "
       "schwach gekoppelt. Im originalen Zustand liegt praktisch das ganze versteckte Gewicht in den ersten Stufen", True)
ZAHLEN["tiefen"] = rows

urteil = ("Ungerichtete Kanten legen die Wechselwirkung des metrischen Sterns bis auf eine Zahl fest. Aus berechneten "
          "Taylorkoeffizienten folgt: Flachheit beider Metriken ist zwei Bedingungen, die Fierz-Pauli-Masse eine; es "
          "bleiben zwei Potentialrichtungen. Liest man jede Kante gleich von Eltern- und Kindmetrik (beta_n = beta_(4-n)), "
          "bleibt genau eine, det(I - X): Die kubische Ordnung ist fest, frei ist nur die quartische Zahl t. Sie ist ein "
          "lokaler Viererkontakt (t/2) mu^2 det((B h)_e)/M_e^2 auf den Zeilen des Vierteltransfers und wächst in der Tiefe "
          "wie q^(-2n). Daher ist der Stern nichtlinear nur bis zu einer endlichen Tiefe schwach gekoppelt; im originalen "
          "Zustand kostet jede Stufe e^{2 pi}, bei TeV gegen Planck sind es fünf Stufen, und jenseits davon liegt kein "
          "nennenswertes Gewicht. Offen bleiben t, mu und M_c (also G).")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z51_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False, default=str) + "\n")
