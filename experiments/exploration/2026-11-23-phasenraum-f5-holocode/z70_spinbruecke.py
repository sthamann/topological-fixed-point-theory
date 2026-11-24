"""Z70: Spinbrücke -- Prüfung der externen Analyse 'Der nächste Schritt zur vollständigen Theorie' (24.11.).

Behauptungen der Analyse, unabhängig nachgerechnet (ihr Paket 'TFPT Spin Bridge' liegt nicht im Repo):
Teil 1  Für die Doppelüberlagerung pi: E -> CP^1 (w^2 = 1 - z^4) und eine nichttriviale 2-Torsions-Linie eta gilt
        pi_* eta = O(-1) + O(-1), also pi_*(eta (x) pi^* O(2)) = O(1) + O(1); ohne eta: pi_* O_E = O + O(-2).
Teil 2  Von den drei nichttrivialen eta ist genau eines unter der Vierteldrehung invariant: das der Antipodenpaarung der
        Marken {1, -1 | i, -i}; als Spinstrukturdifferenz ist es (A,A) - (P,P), der 2-Torsionspunkt (1 + i)/2.
Teil 3  Das 'Quarterdublett' 4,1424 % / 95,8576 % ist der Zustand rho_1 = diag(1, e^-pi)/(1 + e^-pi) (z63, m = 1). Als
        Einteilchen-Einschränkung des Fermipaars +-1/4 ist es nicht KMS für 4 pi x; KMS ist es auf {Omega, Twistzustand}.
Teil 4  Ausschluss: spec(K) in pi Z gibt einen periodischen Modulfluss; ein Lorentz-kovariantes Spektrum positiver Energie
        mit diskreter Energie ist trivial. Der Modulfluss ist also keine physikalische Zeittranslation.
"""
import itertools
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
say("Teil 1: Spaltungstyp der Bildbündel (Riemann-Roch auf E und CP^1)")


def h0_E(deg, trivial):
    """h0 eines Geradenbündels auf einer elliptischen Kurve."""
    if deg > 0:
        return deg
    if deg == 0:
        return 1 if trivial else 0
    return 0


def h0_P1(n):
    return max(n + 1, 0)


def splitting(deg_L, trivial, ks=range(-5, 6)):
    """Spaltungstyp (a, b) von pi_* L aus h0(P1, pi_*L (x) O(k)) = h0(E, L (x) pi^*O(k)), deg(pi^*O(k)) = 2k."""
    seq = {k: h0_E(deg_L + 2*k, trivial and k == 0 or (trivial and deg_L + 2*k == 0)) for k in ks}
    sols = [(a, b) for a in range(-8, 8) for b in range(a, 8)
            if a + b == deg_L - 2 and all(h0_P1(a + k) + h0_P1(b + k) == seq[k] for k in ks)]
    return sols


s_O = splitting(0, True)
s_eta = splitting(0, False)
s_eta2 = splitting(4, False)                         # eta (x) pi^*O(2): Grad 4, nichttrivial
ZAHLEN["spaltung"] = {"pi_* O_E": s_O, "pi_* eta": s_eta, "pi_*(eta (x) pi^*O(2))": s_eta2}
check(f"pi_* O_E = O + O(-2) {s_O}, also O(1) (x) pi_* O_E = O(1) + O(-1) (H6h33); pi_* eta = O(-1) + O(-1) {s_eta}; "
      f"pi_*(eta (x) pi^* O(2)) = O(1) + O(1) {s_eta2}: h0 = 4, h1 = 0",
      s_O == [(-2, 0)] and s_eta == [(-1, -1)] and s_eta2 == [(1, 1)])

# ------------------------------------------------------------------ Teil 2
say("\nTeil 2: Welche 2-Torsion ist unter der Vierteldrehung invariant?")
marks = [1, 1j, -1, -1j]
pairings = [((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))]


def rot(i):            # z -> i z permutiert die Marken: 1 -> i -> -1 -> -i
    return (i + 1) % 4


def norm_pairing(p):
    return tuple(sorted(tuple(sorted(x)) for x in p))


inv = [p for p in pairings if norm_pairing([(rot(a), rot(b)) for a, b in p]) == norm_pairing(p)]
names = [" | ".join(",".join(str(marks[i]).replace("(", "").replace(")", "").replace("j", "i") for i in x) for x in p)
         for p in inv]
z = sp.Symbol("z")
w2 = 1 - z**4
zeros = sp.solve(w2, z)
deg_inf = sp.degree(w2, z)/2
check("die drei nichttrivialen eta entsprechen den drei Paarungen der vier Verzweigungspunkte (P_a - P_b ~ P_c - P_d, "
      "weil div(w) = Summe P - 2 inf_+ - 2 inf_-: w hat vier einfache Nullstellen und Ordnung 2 an beiden Punkten über "
      "unendlich); genau eine Paarung ist unter z -> i z invariant: " + str(names),
      len(inv) == 1 and inv[0] == ((0, 2), (1, 3)) and len(zeros) == 4 and deg_inf == 2)
tors = {"1/2": sp.Rational(1, 2), "i/2": sp.I/2, "(1+i)/2": (1 + sp.I)/2}


def mod_lattice(c):
    re, im = sp.re(c), sp.im(c)
    return sp.nsimplify(re - sp.floor(re)) + sp.I*sp.nsimplify(im - sp.floor(im))


fixed = [k for k, v in tors.items() if sp.simplify(mod_lattice(sp.I*v) - mod_lattice(v)) == 0]
check("auf E = C/(Z + iZ) (tau = i) wirkt die Vierteldrehung als Multiplikation mit i: von den 2-Torsionspunkten 1/2, i/2, "
      "(1+i)/2 ist genau (1+i)/2 fest; Charakteristikdifferenz (1/2, 1/2) = (A,A) - (P,P)", fixed == ["(1+i)/2"])
befund("Die Analyse hält: eta = (A,A) (x) (P,P)^-1 ist wahlfrei (vierteldrehungsinvariant) und beide Spinstrukturen sind "
       "TFPT-Bestand ((A,A) für den Zustand, (P,P) für den flachen E8-Sektor). Damit gibt die Doppelüberlagerung "
       "wahlfrei O(1) + O(1); Hecke-Modifikationen an zwei Marken (OriginLab-Suche) brauchen dagegen eine Wahl", True)

# ------------------------------------------------------------------ Teil 3
say("\nTeil 3: Das 'Quarterdublett' der Analyse")
c = lambda x: 1/(1 + np.exp(4*np.pi*x))
wp, wn = c(0.25), c(-0.25)
rho1 = np.array([1, np.exp(-np.pi)])/(1 + np.exp(-np.pi))
ZAHLEN["quarterdublett"] = [wp, wn]
check(f"die Gewichte {wp:.6f} / {wn:.6f} sind genau rho_1 = diag(1, e^-pi)/(1 + e^-pi) aus z63 (m = 1)",
      abs(wp - rho1[1]) < 1e-12 and abs(wn - rho1[0]) < 1e-12)
pair_N1 = np.array([c(-0.25)*(1 - c(0.25)), c(0.25)*(1 - c(-0.25))])
check("als Einteilchen-Einschränkung des Fermipaars x = -1/4, +1/4 ist das Verhältnis e^-pi, KMS für 4 pi x verlangte "
      "e^-(4 pi * 1/2) = e^-2 pi (so der N=1-Sektor des Paars); KMS mit e^-pi ist es auf {Omega, Twistzustand} (Lücke pi)",
      abs(wp/wn - np.exp(-np.pi)) < 1e-12 and abs(pair_N1[1]/pair_N1[0] - np.exp(-2*np.pi)) < 1e-12)

# ------------------------------------------------------------------ Teil 4
say("\nTeil 4: Ausschluss der periodischen Zeit")
spec = [np.pi*n for n in range(0, 40)]
per = all(abs(np.exp(2j*s) - 1) < 1e-12 for s in spec)
etas = np.linspace(0, 0.2, 2001)
E0 = np.pi                                          # Ruheenergie pi (kleinste Modulanregung)
energies = E0*np.cosh(etas)
off_grid = np.any(np.abs(energies/np.pi - np.round(energies/np.pi)) > 1e-3)
check("spec(K) in pi Z gibt exp(2 i K) = 1 (Periode 2); eine Lorentz-kovariante Energie wird unter Boosts kontinuierlich "
      "(E cosh eta), verlässt also jedes diskrete Gitter: diskretes, boostinvariantes Spektrum positiver Energie ist {0}",
      per and off_grid)
befund("Der periodische Modulfluss ist keine physikalische Zeittranslation (stärker als der Ausschluss für beschränkte "
       "Generatoren, gilt auch unbeschränkt). Mit z62 (Modulfluss auf dem Dublett = Drehung) und GRAV.NONCIRCULAR.01 "
       "heißt das: Zeit der Raumzeit muss aus relativen Modulflüssen (Kozyklen zwischen Zuständen) oder aus T kommen", True)

urteil = ("Die Spinbrücke der externen Analyse hält. Für die Doppelüberlagerung gilt pi_* O_E = O + O(-2), aber pi_* eta = "
          "O(-1) + O(-1) für jede nichttriviale 2-Torsion eta, also pi_*(eta (x) pi^* O(2)) = O(1) + O(1). Genau ein eta ist "
          "unter der Vierteldrehung invariant: das der Antipodenpaarung, als Spinstrukturdifferenz (A,A) - (P,P) = (1+i)/2; "
          "beide Spinstrukturen sind TFPT-Bestand, die Konstruktion ist wahlfrei. Das 'Quarterdublett' 4,14 %/95,86 % ist "
          "rho_1 aus z63 (KMS auf {Omega, Twistzustand}); als Einteilchen-Einschränkung des Fermipaars wäre es nicht KMS. "
          "Der periodische Modulfluss (spec in pi Z) ist keine physikalische Zeittranslation. Offen: die Lorentz-Realform "
          "(die eigenen Realstrukturen der Naht geben nur (2,2) oder (4,0)) und die physische Gleichsetzung der Dubletts.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z70_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
