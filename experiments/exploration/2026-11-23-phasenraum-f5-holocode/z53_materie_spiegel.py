"""Z53: Materie -- Gaußzustand aus der freien Naht, Z16-Bilanz der Familie, Spiegellücke, Flavour mit Dirac-Neutrinos.

Teil 1  Gaußzustand: Der KMS-Zustand eines quadratischen Majorana-Hamiltonoperators erfüllt Wick exakt, mit Wechselwirkung
        nicht. Die Prämisse 'freie Gaußquantisierung' der Quelle ist damit die Prämisse 'freie Naht' des Bestands.
Teil 2  Z16: Eine Familie mit nu_R hat 16 Weyl-Fermionen; ihre Ladungen X = 5(B-L) - 4Y sind alle 1 mod 4, also wirkt
        exp(i pi X/2) als i und ist eine Spin-Z4-Struktur. Alle störungstheoretischen Anomalien verschwinden, die globale
        Zahl ist 16 = 0 mod 16 je Familie (ohne nu_R 15). Nach der Klassifikation (Garcia-Etxebarria, Montero) ist das die
        Bedingung, unter der Spiegel ohne Symmetriebrechung eine Lücke bekommen können.
Teil 3  Spielzeug in 0+1 Dimensionen (Fidkowski-Kitaev): Mit erlaubter, nur quartischer Wechselwirkung bleiben 4 Majoranas
        stets entartet; 8 Majoranas bekommen einen eindeutigen Grundzustand mit Lücke.
Teil 4  Flavour: Auch mit Dirac-Neutrinos (zyklische Reste auf beiden Seiten) gibt Delta(6 n^2) kein TM1 mit Mischung.
"""
import json
import sys
from fractions import Fraction as Fr
from itertools import combinations
from pathlib import Path

import numpy as np
from scipy.linalg import expm

HERE = Path(__file__).parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO / "research" / "originlab"))
from process_rules import flavour_groups  # noqa: E402

OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
rng = np.random.default_rng(53)


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val):
    BEFUNDE[name] = bool(val)
    say(f"   [Befund: {'ja' if val else 'nein'}] {name}")


def majoranas(n_modes):
    X = np.array([[0, 1], [1, 0]], complex)
    Y = np.array([[0, -1j], [1j, 0]])
    Z = np.diag([1.0, -1.0]).astype(complex)
    I = np.eye(2, dtype=complex)
    out = []
    for k in range(n_modes):
        for P in (X, Y):
            ops = [Z]*k + [P] + [I]*(n_modes - k - 1)
            M = ops[0]
            for o in ops[1:]:
                M = np.kron(M, o)
            out.append(M)
    return out


# ------------------------------------------------------------------ Teil 1
say("Teil 1: Gaußzustand aus der freien Naht")
g = majoranas(4)
N = len(g)
A = rng.normal(size=(N, N))
A = A - A.T
H0 = sum(0.25j*A[a, b]*g[a] @ g[b] for a in range(N) for b in range(N))


def wick_error(H, beta=1.3):
    rho = expm(-beta*H)
    rho /= np.trace(rho)
    G = np.array([[np.trace(rho @ g[a] @ g[b]) for b in range(N)] for a in range(N)])
    worst = 0.0
    for a, b, c, d in combinations(range(N), 4):
        four = np.trace(rho @ g[a] @ g[b] @ g[c] @ g[d])
        wick = G[a, b]*G[c, d] - G[a, c]*G[b, d] + G[a, d]*G[b, c]
        worst = max(worst, abs(four - wick))
    return worst


e0 = wick_error(H0)
H1 = H0 + 0.7*g[0] @ g[1] @ g[2] @ g[3] + 0.5*g[2] @ g[5] @ g[6] @ g[7]
e1 = wick_error(H1)
say(f"   8 Majoranas: Wick-Abweichung der Vierpunktfunktion im KMS-Zustand, frei {e0:.1e}, mit Quartikterm {e1:.1e}")
check("der KMS-Zustand eines quadratischen Majorana-Hamiltonoperators ist exakt gaußsch (Wick), mit Wechselwirkung nicht",
      e0 < 1e-12 and e1 > 1e-3)
befund("Die Prämisse 'freie Gaußquantisierung' der elliptischen Quelle ist die Prämisse 'freie Naht auf zwei Kopien' des "
       "Bestands (premise:free-seam-on-two-copies-20261118) plus KMS; b = 4 pi legt die elliptische Geometrie fest", True)

# ------------------------------------------------------------------ Teil 2
say("\nTeil 2: Z16-Bilanz einer Familie")
fam = [  # (Name, Anzahl Weyl, B-L, Y, SU(3)-Index T3, SU(2)-Index T2, Farben, Dublett)
    ("Q", 6, Fr(1, 3), Fr(1, 6), 3, 2), ("u^c", 3, Fr(-1, 3), Fr(-2, 3), 3, 1), ("d^c", 3, Fr(-1, 3), Fr(1, 3), 3, 1),
    ("L", 2, Fr(-1), Fr(-1, 2), 1, 2), ("e^c", 1, Fr(1), Fr(1), 1, 1), ("nu^c", 1, Fr(1), Fr(0), 1, 1)]
X = {f[0]: 5*f[2] - 4*f[3] for f in fam}
say(f"   X = 5(B-L) - 4Y: {dict((k, str(v)) for k, v in X.items())}")
check("alle X-Ladungen sind 1 mod 4 (10: 1, 5bar: -3, 1: 5): exp(i pi X/2) wirkt als i, sein Quadrat ist die "
      "Fermionparität (Spin-Z4)", all(v.denominator == 1 and int(v) % 4 == 1 for v in X.values()))
anom = {
    "X": sum(f[1]*X[f[0]] for f in fam),
    "X^3": sum(f[1]*X[f[0]]**3 for f in fam),
    "X SU(3)^2": sum(Fr(f[1], f[4])*X[f[0]]*Fr(1, 2) for f in fam if f[4] == 3),
    "X SU(2)^2": sum(Fr(f[1], f[5])*X[f[0]]*Fr(1, 2) for f in fam if f[5] == 2),
    "X Y^2": sum(f[1]*X[f[0]]*f[3]**2 for f in fam),
    "X^2 Y": sum(f[1]*X[f[0]]**2*f[3] for f in fam),
    "Y": sum(f[1]*f[3] for f in fam), "Y^3": sum(f[1]*f[3]**3 for f in fam),
}
say(f"   Anomaliekoeffizienten je Familie: {dict((k, str(v)) for k, v in anom.items())}")
check("alle störungstheoretischen Anomalien (X, X^3, X-SU(3)^2, X-SU(2)^2, X-Y^2, X^2-Y, Y, Y^3) verschwinden",
      all(v == 0 for v in anom.values()))
weyl = sum(f[1] for f in fam)
say(f"   globale Zahl nu = Anzahl Weyl mit Ladung 1 mod 4: je Familie {weyl} = {weyl % 16} mod 16; ohne nu^c {weyl - 1} = "
    f"{(weyl - 1) % 16} mod 16; drei Familien {3*weyl % 16} bzw. {3*(weyl - 1) % 16} mod 16")
check("mit nu_R ist nu = 16 = 0 mod 16 je Familie, ohne nu_R 15 (drei Familien: 0 gegen 13): nur mit nu_R ist die "
      "Spin-Z4-Anomalie frei", weyl % 16 == 0 and (weyl - 1) % 16 == 15 and (3*(weyl - 1)) % 16 == 13)
befund("Die 16er-Naht des Bestands (16 Majorana, Kitaev nu = 16 = 0 mod 16, ohne Anyonen) und die 16 Weyl einer Familie "
       "sind dieselbe mod-16-Bilanz; die notwendige Bedingung für Spiegel ohne Symmetriebrechung ist erfüllt, die "
       "Dynamik der Lücke in 4D ist damit nicht bewiesen", True)

# ------------------------------------------------------------------ Teil 3
say("\nTeil 3: Spielzeug 0+1 D -- quartische Lücke")
g4 = majoranas(2)
H4 = g4[0] @ g4[1] @ g4[2] @ g4[3]
ev4 = np.linalg.eigvalsh((H4 + H4.conj().T)/2)
say(f"   4 Majoranas, einziger Quartikterm: Spektrum {np.round(ev4, 6)}")
check("4 Majoranas: der einzige erlaubte Term ist die Parität, jedes Niveau ist zweifach entartet", np.allclose(ev4, [-1, -1, 1, 1]))
g8 = majoranas(4)
quads = list(combinations(range(8), 4))
uniq, trials = 0, 30
gaps = []
for _ in range(trials):
    w = rng.normal(size=len(quads))
    H = sum(wi*g8[a] @ g8[b] @ g8[c] @ g8[d] for wi, (a, b, c, d) in zip(w, quads))
    H = (H + H.conj().T)/2
    ev = np.linalg.eigvalsh(H)
    gaps.append(ev[1] - ev[0])
    uniq += ev[1] - ev[0] > 1e-8
say(f"   8 Majoranas, zufällige reelle Quartikterme: eindeutiger Grundzustand in {uniq} von {trials} Fällen, kleinste Lücke "
    f"{min(gaps):.3f}")
check("8 Majoranas: reelle Quartikterme (zeitumkehrinvariant, quadratische Terme verboten) geben generisch einen "
      "eindeutigen Grundzustand mit Lücke (Fidkowski-Kitaev, Z8)", uniq == trials)
befund("Die Spiegelhälfte kann durch Wechselwirkung allein eine Lücke bekommen, wenn die Zahl der Moden zur Anomalie passt "
       "(0+1 D: 8, 2+1 D: 16, 3+1 D: 16 Weyl je Familie mit Spin-Z4); TFPTs Zählung passt in allen drei Stufen", True)

# ------------------------------------------------------------------ Teil 4
say("\nTeil 4: Flavour mit Dirac-Neutrinos")
fg = flavour_groups(ns=(2, 4, 8), neutrinos="cyclic")
mixed = [o for o in fg if "mixes" in o["clauses"]]
say(f"   Delta(24), Delta(96), Delta(384): {len(fg)} Muster, davon {len(mixed)} gemischt; TM1 mit Mischung: "
    f"{sum('tm1' in o['clauses'] for o in mixed)}, im Datenbereich: {sum('fits_data' in o['clauses'] for o in fg)}")
check("auch mit zyklischen Resten auf beiden Seiten (Dirac-Neutrinos) gibt es kein gemischtes TM1-Muster und keinen "
      "Datentreffer", not any("tm1" in o["clauses"] for o in mixed) and not any("fits_data" in o["clauses"] for o in fg))
befund("TM1 bei 15 Grad ist in Delta(6 n^2) weder mit Majorana- noch mit Dirac-Neutrinos ein Restsymmetriemuster; es braucht "
       "eine Brechung oder eine CP-Regel, die den Winkel festlegt", True)

urteil = ("Drei Lücken auf der Materieseite sind enger. Erstens ist die Gaußquantisierung keine eigene Prämisse: Der "
          "KMS-Zustand eines freien Majorana-Hamiltonoperators erfüllt Wick exakt, also folgt sie aus der freien Naht des "
          "Bestands. Zweitens hat eine Familie mit nu_R 16 Weyl-Fermionen, deren Ladungen X = 5(B-L) - 4Y alle 1 mod 4 sind: "
          "exp(i pi X/2) ist eine Spin-Z4-Struktur, alle störungstheoretischen Anomalien verschwinden, und die globale Zahl "
          "ist 16 = 0 mod 16 (ohne nu_R 15). Das ist die notwendige Bedingung, unter der die Spiegelhälfte ohne "
          "Symmetriebrechung eine Lücke bekommen kann; dieselbe mod-16-Bilanz trägt die 16er-Naht. Im 0+1-D-Spielzeug bleiben "
          "4 Majoranas entartet, 8 bekommen eine Lücke. Die Dynamik der Lücke in 4D ist damit nicht gezeigt. Drittens gibt "
          "Delta(6 n^2) auch mit Dirac-Neutrinos kein TM1; die Familienregel braucht eine Brechung oder CP-Regel.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z53_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False, default=str) + "\n")
