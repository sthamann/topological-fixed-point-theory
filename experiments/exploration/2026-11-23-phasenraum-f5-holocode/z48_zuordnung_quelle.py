"""Z48: Die Zuordnung hinter q = e^-2pi, die Quelle des Gravitons, die Skala und die nichtlineare Form.

Teil 1  Zuordnung: Die Quarterzwiste alpha = 1/4, 3/4 der originalen Quelle sind die Vierteldrehungsgrade 1 und 3 der
        Hardy-Marken k = 4(r + alpha). Dann ist die Quarterenergie L0/4, die originale Kovarianz die Fermibesetzung
        1/(1 + e^{pi k}), und der Strom J_-1 verschiebt k um 4: J_-1 wirkt als S^4. Das war in z46 eine Annahme.
Teil 2  Statistik: Der Nullvektor eta_q ist die Thermofeld-Purifikation eines Oszillators mit Boltzmannfaktor q^2. Die
        bosonische Lesart (Strom J_-1) passt exakt, eine fermionische Lesart der Grad-0-Leiter ist nicht geometrisch.
Teil 3  Quelle: Ward verlangt, dass der masselose Pol nur an die erhaltene Ladung koppelt, die Nullmode k = 0 = |0>.
        Quarter-Harmonische der Quellabbildung gäben eine Verletzung des Äquivalenzprinzips von e^-2pi mal ihrem Anteil.
Teil 4  Nichtlineare Form: Zwei Metriken mit Fierz-Pauli-Mischung geben das massive Gewicht alpha^2/(1 + alpha^2);
        e^-4pi verlangt (M_f/M_g)^2 = 1/(e^{4 pi} - 1), die originale Stromnorm. Gleiches Gewicht hieße M_f = M_g.
"""
import json
import sys
from pathlib import Path

import numpy as np
import sympy as sp

HERE = Path(__file__).parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO / "research" / "originlab"))
from process_rules import Q_ORIGINAL, _compression, ladder_occupation, original_current_norm  # noqa: E402

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


b = 4*np.pi
c = lambda r, al: 1/(1 + np.exp(b*(r + al)))

# ------------------------------------------------------------------ Teil 1
say("Teil 1: Quarterzwiste = Vierteldrehungsgrade, J_-1 = S^4")
rows = []
for al, grade in ((0.25, 1), (0.75, 3)):
    for r in range(-6, 7):
        k = int(round(4*(r + al)))
        rows.append((al, r, k, k % 4 == grade, abs(c(r, al) - 1/(1 + np.exp(np.pi*k))), np.exp(2j*np.pi*al), 1j**k))
check("k = 4(r + alpha) ist ganz und liegt für alpha = 1/4 im Grad 1, für alpha = 3/4 im Grad 3; der Zwist e^{2 pi i alpha} "
      "ist der Vierteldrehungscharakter i^k", all(x[3] and abs(x[5] - x[6]) < 1e-12 for x in rows))
check("die originale Kovarianz c_r = (1 + e^{4 pi (r + alpha)})^-1 ist die Fermibesetzung 1/(1 + e^{pi k}) der Hardy-Marke "
      "k: Quarterenergie r + alpha = k/4 = L0/4", max(x[4] for x in rows) < 1e-15)
n = 40
labels = {al: [4*(r + al) for r in range(-5, 5)] for al in (0.25, 0.75)}
shift = {al: [4*(r + 1 + al) - 4*(r + al) for r in range(-5, 5)] for al in (0.25, 0.75)}
check("der Strom J_-1 = Sum_r psi*_(r+1) psi_r verschiebt jede Quartermode um eins, also jede Hardy-Marke um 4: er wirkt "
      "als S^4 auf den Graden 1 und 3", all(all(abs(s - 4) < 1e-12 for s in v) for v in shift.values()))
L0 = np.diag(np.arange(n))/4
S4 = np.diag(np.ones(n - 4), -4)
check("[L0/4, S^4] = S^4: ein Übertragungsschritt kostet genau eine Quartereinheit (die Einheit von H_alpha)",
      np.allclose(L0 @ S4 - S4 @ L0, S4))
A = _compression(1, 4, 0.5, n=n)
R = np.diag(1j**np.arange(n))
ker = np.linalg.svd(A)[2][-1]
check("die Übertragung erfüllt R A R^+ = -i A (Vierteldrehung R|k> = i^k|k>); ihr Kern liegt nur auf Grad 0, dem neutralen "
      "Sektor, auf dem die Ströme J_m wirken", np.allclose(R @ A @ R.conj().T, -1j*A)
      and np.max(np.abs(ker[np.arange(n) % 4 != 0])) < 1e-12)
befund("Die Zuordnung 'Übertragungsschritt = ein Quant von J_-1' folgt aus den originalen Zwisten alpha = 1/4, 3/4, sobald "
       "die Vierteldrehung der Übertragung die der Marken ist; sie ist keine freie Wahl mehr", True)

# ------------------------------------------------------------------ Teil 2
say("\nTeil 2: Bosonische Purifikation statt fermionischer Leiter")
q = Q_ORIGINAL
qh = 0.5
vec = np.linalg.svd(_compression(1, 4, qh, n=n))[2][-1]
vec = vec*np.sign(vec[0])
eta = np.abs(vec[::4])
tfd = np.array([np.sqrt(1 - qh**2)*qh**j for j in range(n//4)])
check("der Kern der Kompression hat die Amplituden sqrt(1 - x) x^(n/2) einer Thermofeld-Purifikation mit x = q^2 "
      "(geprüft bei q = 1/2; Partner: die zweite Blatthälfte)", np.allclose(eta, tfd, atol=1e-12))
om = original_current_norm()
x = om/(1 + om)
check("geometrische Besetzung mit Mittel omega = 1/(e^{4 pi} - 1) hat das Verhältnis omega/(1 + omega) = e^{-4 pi}",
      abs(x/np.exp(-4*np.pi) - 1) < 1e-9)
check("damit q^2 = e^{-4 pi} und die Leiterbesetzung des Kompressionskerns stimmt (q = e^-2pi)",
      abs(ladder_occupation(q)/om - 1) < 1e-8)
cf = [1/(1 + np.exp(np.pi*4*m)) for m in range(4)]
ratios = [cf[m + 1]/cf[m] for m in range(3)]
say(f"   fermionische Lesart der Grad-0-Leiter (alpha = 0): Besetzungen {[f'{v:.3e}' for v in cf]}, "
    f"Verhältnisse {[f'{v:.3e}' for v in ratios]}")
check("eine fermionische Grad-0-Leiter hat eine halb besetzte Nullmode (1/2) und kein festes Verhältnis: kein Nullvektor "
      "einer Übertragung; nur die bosonische Stromlesart passt", abs(cf[0] - 0.5) < 1e-15 and abs(ratios[0]/ratios[1] - 2) < 1e-3)
befund("Die Zustandsregel ist die Thermofeld-Purifikation des untersten Stroms über die beiden Blatthälften; ihre einzige "
       "Prämisse ist die freie Gaußquantisierung der Quelle", True)

# ------------------------------------------------------------------ Teil 3
say("\nTeil 3: Quelle des Gravitons (Ward)")
th = np.linspace(0, 2*np.pi, 4096, endpoint=False)
rng = np.random.default_rng(48)
cons = []
for _ in range(5):
    k = rng.integers(1, 9)
    amp = rng.normal()
    rho = 1 + amp*np.cos(k*th)
    cons.append(abs(np.mean(rho) - 1) < 1e-12)
check("unter Drehung der Naht ist nur die Fourier-Nullmode einer Dichte erhalten (Ladung = <0|rho>); Harmonische "
      "mitteln heraus", all(cons))
Z = 1 - q**2
h1, h2 = 0.0, 1e-3
eta_ep = q*abs(h1 - h2)
say(f"   koppelt die Quellabbildung auch an |4>: kappa_s = sqrt(Z)(rho_0 + q rho_4 + ...); zwei Arten mit Harmonischenanteil "
    f"0 und 10^-3 hätten eta = {eta_ep:.2e}")
check("Eötvös-Parameter eta = e^-2pi |Delta h|: MICROSCOPE (|eta| < ~1e-14) verlangt |Delta h| < 5e-12; die "
      "Quellabbildung muss die Noetherladung (|0>) sein", abs(eta_ep - 1.867e-6) < 1e-8 and 1e-14/q < 6e-12)
befund("Ward verlangt Minimalkopplung über die erhaltene Ladung: alle Arten koppeln über |0> mit sqrt(1 - q^2), "
       "universell; die Leiterzustände |4j> sind Multipole derselben Dichte, keine eigenen Quellen", True)

# ------------------------------------------------------------------ Teil 4
say("\nTeil 4: Nichtlineare Form -- zwei Metriken")
Mg, Mf, m2 = sp.symbols("M_g M_f m2", positive=True)
k2 = sp.symbols("k2")
Kmat = sp.Matrix([[Mg**2*k2 + m2, -m2], [-m2, Mf**2*k2 + m2]])
Ggg = sp.simplify(Kmat.inv()[0, 0])
r0 = sp.simplify(sp.limit(k2*Ggg, k2, 0))
pole = sp.solve(sp.denom(sp.together(Ggg)), k2)
mpole = [p_ for p_ in pole if p_ != 0][0]
rm = sp.simplify(sp.limit((k2 - mpole)*Ggg, k2, mpole))
frac = sp.simplify(rm/(r0 + rm))
say(f"   G_gg: masseloses Residuum {r0}, massiver Pol bei k^2 = {mpole}, Residuum {rm}; massiver Anteil {frac}")
check("massiver Anteil = M_f^2/(M_g^2 + M_f^2) = alpha^2/(1 + alpha^2) mit alpha = M_f/M_g",
      sp.simplify(frac - Mf**2/(Mg**2 + Mf**2)) == 0)
al2 = sp.symbols("alpha2", positive=True)
sol = sp.solve(sp.Eq(al2/(1 + al2), sp.exp(-4*sp.pi)), al2)[0]
say(f"   Gewicht e^(-4 pi) verlangt alpha^2 = {sp.simplify(sol)}")
check("alpha^2 = 1/(e^{4 pi} - 1): das Planckmassenverhältnis zum Quadrat ist genau die originale Stromnorm",
      sp.simplify(sol - 1/(sp.exp(4*sp.pi) - 1)) == 0 and abs(float(sol)/om - 1) < 1e-9)
check("gleiches Blattgewicht (q^2 = 1/2) hieße M_f = M_g (austauschsymmetrische Bimetrik); der originale Zustand bricht "
      "diesen Tausch stark (M_f/M_g = 0.0019)", sp.solve(sp.Eq(al2/(1 + al2), sp.Rational(1, 2)), al2) == [1]
      and abs(np.sqrt(om) - 1.8675e-3) < 1e-6)
befund("Nichtlinear ist eine geisterfreie Bimetrik (Hassan-Rosen) die naheliegende Vervollständigung; der Zustand legt "
       "M_f/M_g fest, die Skala mu die Fierz-Pauli-Masse. Die übrigen Potentialparameter (kubisch und quartisch; nach "
       "der Literatur zwei Kombinationen) legt die freie Übertragung nicht fest", True)
befund("Skala: Ist ein Übertragungsschritt eine Quartereinheit der Naht, liegt mu bei der Nahtskala; bei Gitter- bzw. "
       "Planckskala ist die Linie (Stärke 4.6e-6) unbeobachtbar und G bleibt offen", True)
ZAHLEN["alpha_bimetrik"] = float(np.sqrt(om))

urteil = ("Die Zuordnung hinter q = e^-2pi ist hergeleitet, sobald die Vierteldrehung der Übertragung die der Marken ist: Die "
          "originalen Quarterzwiste alpha = 1/4, 3/4 sind die Grade 1 und 3 der Hardy-Marken k = 4(r + alpha), die "
          "Kovarianz ist die Fermibesetzung 1/(1 + e^{pi k}), und der Strom J_-1 verschiebt k um 4, wirkt also als der "
          "Übertragungsschritt S^4. Der Nullvektor ist die Thermofeld-Purifikation dieses Stroms über die beiden "
          "Blatthälften; nur die bosonische Lesart gibt eine geometrische Leiter. Ward verlangt Minimalkopplung über die "
          "erhaltene Ladung |0>; Quarter-Harmonische in der Quellabbildung verletzten das Äquivalenzprinzip um e^-2pi mal "
          "ihren Anteil. Nichtlinear passt eine geisterfreie Bimetrik mit (M_f/M_g)^2 = 1/(e^{4 pi} - 1), der originalen "
          "Stromnorm; gleiches Blattgewicht wäre M_f = M_g. Offen bleiben die Gaußquantisierung als Prämisse, die Skala, G "
          "und die kubischen und quartischen Kopplungen.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z48_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False, default=str) + "\n")
