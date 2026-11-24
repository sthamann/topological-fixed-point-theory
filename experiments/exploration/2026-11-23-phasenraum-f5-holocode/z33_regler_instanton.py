"""Z33: Halbgruppenregler und Positivität, Instantonen als kontinuierliche Wandquelle, axiale Anomalien, xi.

Prüft die sechste externe Analyse (24.11.):
Teil 1  Exponentiell regulierter Propagator C(p) = e^{-eps (p^2 + m^2)}/(p^2 + m^2): Hankel-Matrizen C(t_i + t_j) bei
        positiven Zeiten (Reflexionspositivität) und Stieltjes-Bedingung (x C(x) nichtfallend).
Teil 2  Nullmoden je Einheitsinstanton einer Wurzel-SU(2) von Spin(10): 10 -> 2, 16 -> 4, Eltern 4, Familien 12;
        Sättigung der Eltern durch Phi^2, Unterdrückung e^{-2 pi / alpha} bei den Kopplungen aus Z25.
Teil 3  Axiale Anomalien: Spin(10)^2 U(1)_A = 2 + 2 sum r, Gravitation 32 + 16 sum r: keine gemeinsame Lösung.
Teil 4  Nichtminimale Kopplung: 1/(4 G_ind) = (1 - 6 xi)/(48 pi eps^2), flache Verschränkung xi-unabhängig.
"""
import itertools
import json
from fractions import Fraction as Fr
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.integrate import quad

OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
HERE = Path(__file__).parent


def say(s=""):
    print(s)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val):
    BEFUNDE[name] = bool(val)
    say(f"   [Befund: {'ja' if val else 'nein'}] {name}")


# ------------------------------------------------------------------ Teil 1
say("Teil 1: exponentiell regulierter Propagator und Reflexionspositivität")
m = 1.0
def C(t, eps):
    f = lambda p: np.exp(-eps * (p * p + m * m)) / (p * p + m * m) / np.pi
    return quad(f, 0, np.inf, weight="cos", wvar=t, limit=400)[0]
ts = [0.05, 0.1, 0.2, 0.4, 0.8, 1.6, 3.2]
mins = {}
for eps in (0.0, 0.01, 0.1, 1.0):
    worst = 0.0
    for t1, t2 in itertools.combinations(ts, 2):
        M = np.array([[C(2 * t1, eps), C(t1 + t2, eps)], [C(t1 + t2, eps), C(2 * t2, eps)]])
        ev = np.linalg.eigvalsh(M / np.trace(M))
        worst = min(worst, ev[0])
    mins[eps] = worst
    say(f"   eps = {eps}: kleinster normierter Eigenwert der 2x2-Hankel-Matrizen {worst:.2e}")
check("eps = 0: Hankel-Matrizen positiv (bis Rundung)", mins[0.0] > -1e-9)
check("eps > 0: negative Eigenwerte, also Reflexionspositivität verletzt", all(mins[e] < -1e-5 for e in (0.01, 0.1, 1.0)))
xs = np.array([1.0, 2.0, 5.0, 10.0, 50.0])
xC = xs * np.exp(-0.1 * xs) / xs
check("Stieltjes verlangt x C(x) nichtfallend; e^{-eps x} fällt, also keine positive Källén-Lehmann-Darstellung", np.all(np.diff(xC) < 0))
say("   Lesart: Der Wärmeleitungsregler ist als Eigenzeitregularisierung der effektiven Wirkung (log det, Anomalien,")
say("   Z31) zulässig, nicht als endlicher physischer Schnitt des Propagators.")
befund("Halbgruppenregler als physischer Propagatorschnitt verletzt Reflexionspositivität", True)

# ------------------------------------------------------------------ Teil 2
say("\nTeil 2: Instantonnullmoden unter einer Wurzel-SU(2) von Spin(10)")
w10 = [tuple(s * (1 if i == k else 0) for i in range(5)) for k in range(5) for s in (1, -1)]
w16 = [w for w in itertools.product((Fr(1, 2), Fr(-1, 2)), repeat=5) if sum(1 for x in w if x < 0) % 2 == 0]
dbl = lambda ws: sum(1 for w in ws if w[0] - w[1] == 1)      # Wurzel alpha = e1 - e2
n10, n16 = dbl(w10), dbl(w16)
zm = {"Vektor 10": n10, "Spinor 16": n16, "Eltern (psi, chi je (10,1)+(1,6))": 2 * n10, "drei Familien": 3 * n16}
zm["gesamt"] = zm["Eltern (psi, chi je (10,1)+(1,6))"] + zm["drei Familien"]
say(f"   Dubletts (= Nullmoden je Einheitsinstanton, Weyl): {zm}")
check("10 -> 2, 16 -> 4, Eltern 4, Familien 12, gesamt 16", (n10, n16, zm["gesamt"]) == (2, 4, 16))
check("Eltern-Nullmoden: zwei Einsetzungen der Paarung Phi psi chi sättigen 4, also Phi^2 im 't-Hooft-Vertex", zm["Eltern (psi, chi je (10,1)+(1,6))"] == 2 * 2)
supp = {a_inv: float(np.exp(-2 * np.pi * a_inv)) for a_inv in (4.0, 23.0, 60.0)}
say(f"   Unterdrückung e^(-2 pi alpha^-1): {dict((k, f'{v:.1e}') for k, v in supp.items())} (Z25: induziert 4 bzw. <= 23, nackt SO(10) 60)")
check("bei nackter SO(10) (alpha^-1 = 60) ist der Instantonbeitrag ~1e-164 gegen den Gitteroperator O(1) (Z24)", supp[60.0] < 1e-160)
say("   Ungesättigte Familiennullmoden (masselose Familien) geben keinen Vakuumbeitrag; gesättigt (Yukawas, EW-Brechung) ein")
say("   periodisches Potential ~ Re[e^(i theta) Phi^2 C_fam]: die Nahtphase wäre axionartig, die Wand auch im Kontinuum möglich.")
befund("Wand im Kontinuum möglich über Spin(10)-Instantonen, aber nur bei gesättigten Familiennullmoden und stark unterdrückt", True)

# ------------------------------------------------------------------ Teil 3
say("\nTeil 3: axiale Anomalien (Eltern psi, chi axial +1, Familien r_1, r_2, r_3)")
r = sp.symbols("r1:4")
T10, T16 = 1, 2
mixed = 2 * T10 + T16 * sum(r)
grav = 2 * 16 + 16 * sum(r)
sol_m = sp.solve(sp.Eq(mixed, 0), r[0])
sol_g = sp.solve(sp.Eq(grav, 0), r[0])
say(f"   Spin(10)^2 U(1)_A: {mixed} = 0 => sum r = {sp.simplify(sol_m[0] + r[1] + r[2])}; Gravitation: {grav} = 0 => sum r = "
    f"{sp.simplify(sol_g[0] + r[1] + r[2])}")
check("keine gemeinsame Lösung (sum r = -1 gegen -2)", sp.simplify(sol_m[0] - sol_g[0]) != 0)
say(f"   sum r = -1: Spin(10)-Anteil 0, Gravitation {grav.subs({r[0]: -1, r[1]: 0, r[2]: 0})}; sum r = 0 (H6: Familien axial neutral): "
    f"Spin(10)-Anteil {mixed.subs({r[0]: 0, r[1]: 0, r[2]: 0})}, Gravitation {grav.subs({r[0]: 0, r[1]: 0, r[2]: 0})}")
check("mit axial neutralen Familien: Spin(10)-Anomalie 2 (Instantonen brechen U(1)_A), Gravitation 32 (c = 8)",
      mixed.subs({x: 0 for x in r}) == 2 and grav.subs({x: 0 for x in r}) == 32)
befund("U(1)_A ist keine anomaliefreie interne Symmetrie; der Zufluss ist notwendiger Bestandteil", True)

# ------------------------------------------------------------------ Teil 4
say("\nTeil 4: nichtminimale Kopplung")
xi, eps = sp.symbols("xi epsilon", positive=True)
s_ = sp.symbols("s", positive=True)
invG4 = sp.simplify(16 * sp.pi * sp.Rational(1, 2) * (sp.Rational(1, 6) - xi) * sp.integrate(1 / (16 * sp.pi**2 * s_**2), (s_, eps**2, sp.oo)) / 4)
Sent = sp.integrate(1 / (12 * 4 * sp.pi * s_**2), (s_, eps**2, sp.oo))
say(f"   1/(4 G_ind) = {sp.factor(invG4)}, S_ent/A (flach) = {Sent}, R = S_ent / (A/4G) = {sp.simplify(Sent / invG4)}")
check("1/(4 G_ind) = (1 - 6 xi)/(48 pi eps^2); R = 1/(1 - 6 xi); xi = 0 nötig, xi = 1/6 gibt keinen Beitrag",
      sp.simplify(invG4 - (1 - 6 * xi) / (48 * sp.pi * eps**2)) == 0 and sp.simplify(Sent / invG4 - 1 / (1 - 6 * xi)) == 0)

urteil = ("Ein exponentiell regulierter Propagator verletzt die Reflexionspositivität: Die 2x2-Hankel-Matrizen C(t_i + t_j) "
          "werden für eps > 0 indefinit (für eps = 0 positiv), und x C(x) fällt, eine positive Källén-Lehmann-Darstellung "
          "verlangt das Gegenteil. Der Wärmeleitungsregler aus Z31 ist daher nur als Eigenzeitregularisierung der "
          "effektiven Wirkung zulässig, nicht als physischer Schnitt. Unter einer Wurzel-SU(2) von Spin(10) hat ein "
          "Einheitsinstanton 4 Eltern- und 12 Familiennullmoden; die Eltern sättigt Phi^2, also kann ein periodisches "
          "Potential Re[e^(i theta) Phi^2 C_fam] entstehen, sobald die Familiennullmoden gesättigt sind: Die Wand ist auch "
          "ohne Gitter möglich, die Nahtphase axionartig. Bei nackter SO(10) (alpha^-1 = 60, Z25) ist der Beitrag ~1e-164, "
          "gegen den Gitteroperator aus Z24 vernachlässigbar. Die axiale Ladung kann mit einheitlich geladenen Familien "
          "nicht zugleich Spin(10)- und gravitativ anomaliefrei sein (sum r = -1 gegen -2); mit neutralen Familien bleibt "
          "die Spin(10)-Anomalie 2 und die gravitative 32 (c = 8). Nichtminimale Kopplung gibt 1/(4 G_ind) = (1 - 6 xi)/"
          "(48 pi eps^2), also R = 1/(1 - 6 xi); die Halbgruppe allein legt xi nicht fest.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z33_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": {"hankel_min": mins, "nullmoden": zm, "unterdrueckung": supp},
     "urteil": urteil}, indent=1, ensure_ascii=False, default=str) + "\n")
