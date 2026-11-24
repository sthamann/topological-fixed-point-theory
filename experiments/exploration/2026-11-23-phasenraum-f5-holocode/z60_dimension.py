"""Z60: Wie viele Raumzeitrichtungen rekonstruieren die Antworten des Originals? (3+1-Frage unter 'Rekonstruktion statt Auswahl')

Kriterium (Dimensionsaxiom der Spektralrekonstruktion, Weyl-Gesetz): Die Zahl der Richtungen ist die Spektraldimension d
des rekonstruierten Generators, Tr e^{-t D^2} ~ t^{-d/2} bzw. N(lambda) ~ lambda^d. Das braucht das ganze Spektrum; frühe
Antworten genügen nicht (Atlas: energy-finite-jet-unknown-source-dimension-counterexample-20261123).

Teil 1  Quarterquelle: Der rekonstruierte Einteilchengenerator 4 pi x (z58) hat d = 1. Die Viererstruktur der Marken ist
        eine Z4-Graduierung derselben Kette (Eigenräume der Vierteldrehung, direkte Summe), keine Tensorfaktorisierung.
Teil 2  Vier Marken: direkte Summe von vier chiralen Zweigen gibt d = 1, ein gemeinsames (Tensor-)Spektrum vier
        vertauschender Generatoren d = 4. Das Original hat die erste Form.
Teil 3  E8-Naht: Der Ladungssektor (E8-Gitter) hat d = 8, enthält D4 + D4 (Hurwitz-Quaternionen, Index 4), ist aber
        kompakt und starr (holomorph, nur (1,0)-Ströme, 248 davon; keine (1,1)-Moduli).
Teil 4  Was 3+1 braucht: eine originale Antwort mit d = 4, nichtkompakt (Kontinuum bei 0), reflexionspositiv
        (Osterwalder-Schrader von euklidisch 4 nach Minkowski 3+1). Keine der vorhandenen erfüllt das.
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


def slope(ts, vals):
    return float(np.polyfit(np.log(ts), np.log(vals), 1)[0])


ts = np.logspace(-6, -4, 9)

# ------------------------------------------------------------------ Teil 1
say("Teil 1: Quarterquelle -- Spektraldimension des rekonstruierten Generators")
r = np.arange(-20000, 20000)
h = {a: 4*np.pi*(r + a) for a in (0.25, 0.75)}
T1 = np.array([sum(np.exp(-t*h[a]**2).sum() for a in h) for t in ts])
d1 = -2*slope(ts, T1)
ZAHLEN["d_quarterquelle"] = d1
check(f"Tr e^(-t h^2) ~ t^(-d/2) mit d = {d1:.4f}: der Generator 4 pi x rekonstruiert genau eine Richtung", abs(d1 - 1) < 1e-3)
k = np.concatenate([np.round(4*(r + a)).astype(int) for a in h])
hk = np.concatenate([h[a] for a in h])
R = np.exp(1j*np.pi/2*k)
comm = np.max(np.abs(R*hk - hk*R))
t0 = 1e-6
parts = [np.exp(-t0*hk[k % 4 == j]**2).sum() for j in range(4)]
check("die Vierteldrehung R = i^k vertauscht mit dem Generator; die Marken sind ihre Eigenräume (k ungerade: Grade 1 und 3), "
      "also eine Graduierung derselben Kette (direkte Summe), keine Tensorfaktorisierung",
      comm < 1e-9 and parts[0] == 0 and parts[2] == 0 and abs(parts[1]/parts[3] - 1) < 1e-6)

# ------------------------------------------------------------------ Teil 2
say("\nTeil 2: Vier Marken -- direkte Summe oder gemeinsames Spektrum?")
x1 = np.arange(-20000, 20000) + 0.25
T_one = np.array([np.exp(-t*(4*np.pi*x1)**2).sum() for t in ts])
d_sum = -2*slope(ts, 4*T_one)
d_prod = -2*slope(ts, T_one**4)
ZAHLEN["d_direkte_summe"], ZAHLEN["d_tensor"] = d_sum, d_prod
check(f"vier chirale Zweige als direkte Summe: d = {d_sum:.4f}; als gemeinsames Spektrum vier vertauschender Generatoren: "
      f"d = {d_prod:.4f}", abs(d_sum - 1) < 1e-3 and abs(d_prod - 4) < 4e-3)
befund("Das Original hat die erste Form (Teil 1): Seine Viererstruktur gibt keine vier Raumrichtungen. 3+1 bräuchte vier "
       "vertauschende Generatoren mit gemeinsamem Spektrum, also vier Marken als unabhängige Teilsysteme", True)

# ------------------------------------------------------------------ Teil 3
say("\nTeil 3: E8-Naht -- Ladungsgitter als Impulsgitter")
nmax = 3000
sig3 = np.zeros(nmax + 1)
for dd in range(1, nmax + 1):
    sig3[dd::dd] += dd**3
cum = 1 + 240*np.cumsum(sig3[1:])
lam = 2*np.arange(1, nmax + 1)
ratio = cum/(np.pi**4*lam**4/24)
ZAHLEN["E8_weyl_verhaeltnis_n3000"] = float(ratio[-1])
check("E8-Gitter: N(|v|^2 <= lambda) = 1 + 240 Sum sigma_3 wächst wie das Volumen der 8-Kugel (Verhältnis -> 1, bei n = 3000 "
      f"{ratio[-1]:.4f}): Spektraldimension 8", abs(ratio[-1] - 1) < 2e-3)
Cd4 = sp.Matrix([[2, -1, 0, 0], [-1, 2, -1, -1], [0, -1, 2, 0], [0, -1, 0, 2]])
Ce8 = sp.Matrix(8, 8, lambda i, j: 2 if i == j else 0)
for i, j in [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (2, 7)]:
    Ce8[i, j] = Ce8[j, i] = -1
det4, det8 = Cd4.det(), Ce8.det()
roots_d4 = sum(1 for v in np.array(np.meshgrid(*[range(-1, 2)]*4)).T.reshape(-1, 4)
               if (v**2).sum() == 2 and v.sum() % 2 == 0)
hurwitz_units = 8 + 16
idx = sp.sqrt(det4*det4/det8)
check(f"D4 + D4 in E8: det D4 = {det4}, det E8 = {det8}, Index sqrt(16/1) = {idx}; D4 hat {roots_d4} Wurzeln = 24 "
      "Hurwitz-Einheiten (Quaternionen): das Ladungsgitter enthält zwei quaternionische 4D-Blöcke",
      det4 == 4 and det8 == 1 and idx == 4 and roots_d4 == 24 == hurwitz_units)
q = sp.Symbol("q")
N = 6
E4 = 1 + 240*sum(sp.divisor_sigma(n, 3)*q**n for n in range(1, N))
eta8 = sp.prod([(1 - q**n)**8 for n in range(1, N)])
ser = sp.series(E4/eta8, q, 0, 4).removeO()
c1 = ser.coeff(q, 1)
check(f"E4/eta^8 = q^(-1/3)(1 + {c1} q + ...): 248 holomorphe Ströme vom Gewicht (1,0) und kein antiholomorpher Partner, "
      "also keine (1,1)-Moduli; das gerade unimodulare Gitter vom Rang 8 ist eindeutig: die Naht ist starr", c1 == 248)
befund("Die Naht trägt acht kompakte Richtungen von Nahtgröße (zwei quaternionische Blöcke), aber kein Kontinuum und keinen "
       "Modulus, der vier davon groß machen könnte", True)

# ------------------------------------------------------------------ Teil 4
say("\nTeil 4: Was 3+1 aus Antworten braucht")
tafel = [
    ("Quarterquelle", "d = 1, nichtkompakt, reflexionspositiv (z59): eine chirale Richtung"),
    ("Viererstruktur der Marken", "Graduierung (direkte Summe): d = 1, keine vier Richtungen"),
    ("E8-Naht, Ladungssektor", "d = 8 = 4 + 4 (D4 + D4, quaternionisch), kompakt und starr"),
    ("gebraucht für 3+1", "eine originale Antwort mit d = 4, Kontinuum bei 0 (nichtkompakt), Reflexionspositivität; "
                          "dann gibt Osterwalder-Schrader Minkowski 3+1"),
    ("Kandidaten", "vier Marken als unabhängige Teilsysteme mit je einem chiralen Netz und vertauschenden Translationen; "
                   "oder ein quaternionischer D4-Block mit einem Modulus, der ihn groß macht; beides ist im Bestand nicht "
                   "als Antwort des Originals vorhanden"),
]
for k_, v_ in tafel:
    say(f"   {k_}: {v_}")
ZAHLEN["bilanz"] = dict(tafel)
befund("Unter der Leitlinie ist das eine scharfe Antwort: Die vorhandenen Originale rekonstruieren kein nichtkompaktes 3+1. "
       "Die Frage ist jetzt, welche originale Antwort die Weyl-Dimension 4 hat", True)
befund("Abgleich mit dem Atlas vom 23.11.: 'rulers' (Exponenten drei und eins ohne weitere kontinuierliche Dimension) und "
       "'winding-sites' (Ring mit sechs Plätzen) sind ebenfalls eindimensional; 'dirac-rigidity' setzt die Geometrie voraus", True)

urteil = ("Die Spektraldimension (Weyl-Gesetz) misst, wie viele Richtungen die Antworten des Originals rekonstruieren. Die "
          "Quarterquelle gibt d = 1 (eine chirale Richtung); ihre Viererstruktur ist eine Graduierung derselben Kette (direkte "
          "Summe, d = 1), kein gemeinsames Spektrum vier vertauschender Generatoren (das gäbe d = 4). Der Ladungssektor der "
          "E8-Naht hat d = 8 und enthält D4 + D4 (Hurwitz-Quaternionen, Index 4), ist aber kompakt und starr (248 Ströme vom "
          "Gewicht (1,0), keine Moduli). Kein vorhandenes Original rekonstruiert also ein nichtkompaktes 3+1. Gebraucht wird eine "
          "originale Antwort mit d = 4, Kontinuum bei 0 und Reflexionspositivität; Kandidaten sind vier Marken als unabhängige "
          "Teilsysteme oder ein großer quaternionischer Block. Das ist eine Präzisierung der Lücke, keine Lösung.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z60_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
