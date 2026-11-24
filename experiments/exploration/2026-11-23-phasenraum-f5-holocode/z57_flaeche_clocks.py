"""Z57: Plebański-Flächenregel für die quartische Zahl, Familienmischung aus den originalen Clocks, G.

Teil 1  Flächenregel: Ist die Kantenwirkung quadratisch in der Differenz der Plebański-Zweiformen B = e ^ e,
        eps_abcd (B_p - B_k)^ab ^ (B_p - B_k)^cd, so ist das Potential (-3/2) + e_2/2 - (3/2) e_4 = V + (1/2) det(I - X), also
        t = 1 -- dieselbe Form wie das arithmetische Mittel Q = (e_p + e_k)^2, denn (a^2 - b^2)^2 = (a - b)^2 (a + b)^2.
        Numerisch mit der vollen Epsilon-Kontraktion für zufällige Vielbeine geprüft; Lane-Klausel quadratic_in_area.
Teil 2  Clocks: F_C, F_J (Kompression der Coxeter- und Vierteldrehungsprojektoren laut Analyse, Zertifikat
        compiler-integral-triality-20260918) haben die gemeinsame Nulllinie u = (0,1,1)/sqrt 2; in der zyklischen Dreierbasis
        trägt sie die TM1-Gewichte (2/3, 1/6, 1/6). Im Rest ist tan 2 theta = 2 sqrt 2/(3 + 16 gamma); gamma bleibt frei.
        i[F_C, F_J] vernichtet u ebenfalls: Er gibt eine Phase, bricht TM1 aber nicht.
Teil 3  Eingefrorene Werte: sin^2 theta12 = 1/3 - phi0/2 = 0.306747 und sin^2 theta13 = phi0 e^(-5/6) = 0.023108 sind kein TM1
        (TM1 verlangte 0.317563); delta = 240 Grad ist ein zugewiesener Texturwert. Die Lösung muss die TM1-Brechung tragen.
Teil 4  G: 8 pi G = 1/(c3 v_geo^2) mit c3 = 1/(8 pi) gibt G = v_geo^-2.
"""
import json
import sys
from itertools import permutations
from pathlib import Path

import numpy as np
import sympy as sp

HERE = Path(__file__).parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO / "research" / "originlab"))
from process_rules import BETA_DET, BETA_V, potential_taylor, run  # noqa: E402

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
say("Teil 1: Plebański-Flächenregel")
c = (1, 0, -2, 0, 1)                                    # (a^2 - b^2)^2 in eps(e_p^(4-n) e_k^n)
beta = np.array([c[n]/[1, 4, 6, 4, 1][n] for n in range(5)])
beta = beta/(beta[1] + 2*beta[2] + beta[3])
t = 2*(beta[0] - BETA_V[0])
say(f"   beta der Flächenregel (Fierz-Pauli-normiert): {np.round(beta, 6).tolist()}, t = {t}")
check("Flächenregel = (-3/2) + e_2/2 - (3/2) e_4 = V + (1/2) det(I - X): t = 1",
      np.allclose(beta, [-1.5, 0, 0.5, 0, -1.5]) and np.allclose(beta, np.array(BETA_V) + 0.5*np.array(BETA_DET)))
perm4 = list(permutations(range(4)))
eps = np.zeros((4, 4, 4, 4))
for p in perm4:
    eps[p] = np.linalg.det(np.eye(4)[list(p)])


def wedge(Bp, Bk):
    D = Bp - Bk
    return np.einsum("abcd,mnrs,abmn,cdrs->", eps, eps, D, D)


def area(e):
    return np.einsum("am,bn->abmn", e, e) - np.einsum("an,bm->abmn", e, e)


rng = np.random.default_rng(57)
ratios = []
for _ in range(24):
    ep = rng.normal(size=(4, 4)) + 3*np.eye(4)
    S = rng.normal(size=(4, 4))
    S = S @ S.T/4 + np.eye(4)                              # symmetrische positive Relativmatrix (Deser-van-Nieuwenhuizen)
    ek = ep @ S
    lam = np.linalg.eigvalsh(S)
    e = [1, sum(lam), sum(lam[i]*lam[j] for i in range(4) for j in range(i + 1, 4)),
         sum(lam[i]*lam[j]*lam[k] for i in range(4) for j in range(i + 1, 4) for k in range(j + 1, 4)), np.prod(lam)]
    hr = np.linalg.det(ep)*sum(b*en for b, en in zip(beta, e))
    ratios.append(wedge(area(ep), area(ek))/hr)
say(f"   volle Epsilon-Kontraktion / HR-Potential für 24 zufällige Vielbeinpaare: Mittel {np.mean(ratios):.6f}, "
    f"relative Streuung {np.ptp(ratios)/abs(np.mean(ratios)):.1e}")
check("die volle Kontraktion eps eps (B_p - B_k)(B_p - B_k) ist für alle 24 Paare ein fester Vielfacher des HR-Potentials mit t = 1",
      np.ptp(ratios)/abs(np.mean(ratios)) < 1e-10)
a_, b_ = sp.symbols("a b", positive=True)
conf = sp.factor(sum(sp.Rational(bv).limit_denominator(8)*sp.binomial(4, n)*a_**(4 - n)*b_**n for n, bv in enumerate(beta)))
say(f"   konformer Teilraum e_p = a I, e_k = b I: V = {conf}")
check("auf dem konformen Teilraum ist V_1 = -(3/2)(a^2 - b^2)^2", sp.simplify(conf + sp.Rational(3, 2)*(a_**2 - b_**2)**2) == 0)
res = run()
ef = res["universes"]["edge_forms"]
sel = {r["selected"]: r["minimal_rules"] for r in ef["unique_selections"]}
say(f"   Lane edge_forms: {ef['derived']['Q = (e_p + e_k)^2 (arithmetic mean)']['clauses']}")
check("in der Lane wählt 'quadratisch in den Flächen' allein die Mittelform (t = 1); sie ist identisch mit dem "
      "arithmetischen Mittel", ["quadratic_in_area"] in sel.get("Q = (e_p + e_k)^2 (arithmetic mean)", []))
A0, B0 = rng.normal(size=(4, 4)), rng.normal(size=(4, 4))
A0, B0 = (A0 + A0.T)/2, (B0 + B0.T)/2
tay, sw = potential_taylor(tuple(beta), A0, B0), potential_taylor(tuple(beta), B0, A0)
check("die Flächenform ist flach und ungerichtet (Taylor bis Ordnung vier)", abs(tay[0]) < 1e-12 and abs(tay[1]) < 1e-12
      and np.allclose(tay, sw, atol=1e-10))
befund("Mit Flächen (Plebański) als Grundvariablen ist t = 1 die natürliche Regel; mit Vielbeinen und Minimalgrad wäre es t = 4. "
       "Die Wahl der Grundvariablen ist damit die eine offene Entscheidung; die Simplizität der gemittelten Fläche bleibt "
       "außerhalb des konformen Teilraums eine Zusatzbedingung", True)

# ------------------------------------------------------------------ Teil 2
say("\nTeil 2: Familienmischung aus den originalen Clocks")
FC = sp.Matrix([[13, 1, -1], [1, 5, -5], [-1, -5, 5]])/16
FJ = sp.diag(1, 0, 0)
u = sp.Matrix([0, 1, 1])/sp.sqrt(2)
comm = FC*FJ - FJ*FC
hs = sum(abs(v)**2 for v in comm)
check("F_C und F_J vernichten beide u = (0,1,1)/sqrt 2, kommutieren auf dem Rest aber nicht (|[F_C, F_J]|_HS^2 = 1/64)",
      FC*u == sp.zeros(3, 1) and FJ*u == sp.zeros(3, 1) and sp.nsimplify(hs) == sp.Rational(1, 64))
w = sp.exp(2*sp.pi*sp.I/3)
cyc = [sp.Matrix([1, w**k, w**(2*k)])/sp.sqrt(3) for k in range(3)]
wts = [sp.nsimplify(sp.simplify(abs((f.H*u)[0])**2)) for f in cyc]
say(f"   Gewichte von u in der zyklischen Dreierbasis: {wts}")
check("in der zyklischen Dreierbasis trägt u genau die TM1-Gewichte (2/3, 1/6, 1/6)", wts == [sp.Rational(2, 3), sp.Rational(1, 6),
                                                                                         sp.Rational(1, 6)])
e1 = sp.Matrix([1, 0, 0])
wv = sp.Matrix([0, 1, -1])/sp.sqrt(2)
gam = sp.symbols("gamma")
HF = sp.Matrix([[(e.H*(FC + gam*FJ)*f)[0] for f in (e1, wv)] for e in (e1, wv)])
check("im Rest {e_1, (0,1,-1)/sqrt 2} ist H_F = [[13/16 + gamma, sqrt 2/16], [sqrt 2/16, 10/16]]",
      sp.simplify(HF - sp.Matrix([[sp.Rational(13, 16) + gam, sp.sqrt(2)/16], [sp.sqrt(2)/16, sp.Rational(10, 16)]])) == sp.zeros(2, 2))
tan2 = sp.simplify(2*HF[0, 1]/(HF[0, 0] - HF[1, 1]))
check("tan 2 theta = 2 sqrt 2/(3 + 16 gamma)", sp.simplify(tan2 - 2*sp.sqrt(2)/(3 + 16*gam)) == 0)
th = 0.5*np.arctan(2*np.sqrt(2)/5)
s13 = np.sin(th)**2/3
g15 = sp.nsimplify((2*sp.sqrt(6) - 3)/16)
say(f"   gamma = 1/8: theta = {np.degrees(th):.3f} Grad, sin^2 theta13 = {s13:.6f}; theta = 15 Grad verlangte gamma = {g15} = "
    f"{float(g15):.4f}")
check("gamma = 1/8 gibt theta = 14.748 Grad und sin^2 theta13 = 0.021602 (Angabe der Analyse); gamma bleibt eine freie Zahl",
      abs(np.degrees(th) - 14.748) < 1e-3 and abs(s13 - 0.021602) < 1e-6)
icomm = sp.I*comm
check("i[F_C, F_J] ist hermitesch und vernichtet u: er liefert eine Phase im Rest, bricht die TM1-Spalte aber nicht",
      sp.simplify(icomm - icomm.H) == sp.zeros(3, 3) and icomm*u == sp.zeros(3, 1) and icomm != sp.zeros(3, 3))
befund("Die TM1-Struktur steckt in den originalen Clocks (gemeinsame Nulllinie mit Gewichten 2/3, 1/6, 1/6); der Winkel hängt "
       "an einem relativen Gewicht gamma, das nicht hergeleitet ist. Die Zuordnung des Stromindexraums zu den Leptonfamilien "
       "und die Kompression aus dem Zertifikat sind aus der Analyse übernommen", True)

# ------------------------------------------------------------------ Teil 3
say("\nTeil 3: Eingefrorene Vorhersagen")
frozen = json.loads((REPO / "verification/predictions_frozen.json").read_text())
pr = {p["id"]: p for p in frozen["predictions"]}
tex = {p["id"]: p for p in frozen["assigned_texture_values"]}
s12f, s13f = float(pr["SIN2_THETA12_SEED"]["frozen_value"]), float(pr["SIN2_THETA13"]["frozen_value"])
tm1_s12 = (1 - 3*s13f)/(3*(1 - s13f))
say(f"   eingefroren: sin^2 theta12 = {s12f:.6f} (1/3 - phi0/2), sin^2 theta13 = {s13f:.6f} (phi0 e^-5/6); TM1 verlangte "
    f"{tm1_s12:.6f}; delta = {tex['DELTA_CP_NU_DEG']['frozen_value']} Grad ({tex['DELTA_CP_NU_DEG']['layer']})")
check("die eingefrorenen Werte sind kein TM1 (Abstand 0.0108 in sin^2 theta12); delta = 240 Grad ist ein zugewiesener "
      "Texturwert, nicht die mu-tau-Reflexion", abs(tm1_s12 - 0.317563) < 1e-6 and abs(tm1_s12 - s12f) > 0.01
      and tex["DELTA_CP_NU_DEG"]["frozen_value"] == "240")
nf61 = (0.3088, 0.0067)
say(f"   NuFIT 6.1 sin^2 theta12 = {nf61[0]} +- {nf61[1]}: TFPT {abs(s12f - nf61[0])/nf61[1]:.2f} sigma, TM1(15 Grad) "
    f"{abs(0.318107 - nf61[0])/nf61[1]:.2f} sigma")
check("der eingefrorene TFPT-Wert 0.3067 liegt näher an NuFIT 6.1 (0.3 sigma) als TM1 bei 15 Grad (1.4 sigma): exaktes TM1 "
      "nicht erzwingen", abs(s12f - nf61[0])/nf61[1] < 0.5 and abs(0.318107 - nf61[0])/nf61[1] > 1.3)
befund("Die vollständige Familienregel muss TM1 aus den Clocks leicht brechen (der Analyse zufolge über die inverse "
       "Coxeter-Kompression, die die Nulllinie nicht erhält) und dabei die eingefrorenen Werte tragen; ein Gewicht dafür ist "
       "nicht hergeleitet", True)

# ------------------------------------------------------------------ Teil 4
say("\nTeil 4: G")
v, c3s = sp.symbols("v_geo c3", positive=True)
G = sp.solve(sp.Eq(8*sp.pi*sp.Symbol("G"), 1/(c3s*v**2)), sp.Symbol("G"))[0].subs(c3s, 1/(8*sp.pi))
check("8 pi G = 1/(c3 v_geo^2) mit c3 = 1/(8 pi) gibt G = v_geo^-2", sp.simplify(G - 1/v**2) == 0)
befund("G ist damit im Bestand relativ zur Einheit v_geo festgelegt; ein dimensionsloses Verhältnis wie M_Pl/m_e braucht die "
       "gemeinsame Zuordnung von Geometrie- und Materieskala", True)

urteil = ("Drei Punkte der Analyse halten und präzisieren den Stand. Die Plebański-Flächenregel (Kantenwirkung quadratisch in "
          "B_p - B_k) gibt exakt das Potential (-3/2) + e_2/2 - (3/2) e_4, also t = 1; die volle Epsilon-Kontraktion ist für 24 "
          "zufällige Vielbeinpaare ein fester Vielfacher davon, und auf dem konformen Teilraum ist V = -(3/2)(a^2 - b^2)^2. Das "
          "ist dieselbe Form wie das arithmetische Mittel; die Lane wählt sie mit 'quadratisch in den Flächen'. Die originalen "
          "Clocks F_C, F_J haben eine gemeinsame Nulllinie mit TM1-Gewichten (2/3, 1/6, 1/6), der Winkel hängt an einem freien "
          "Gewicht gamma (gamma = 1/8 gibt 14.748 Grad), und i[F_C, F_J] gibt eine Phase ohne TM1-Brechung. Die eingefrorenen "
          "Werte 0.306747 und 0.023108 sind kein TM1; sie liegen näher an NuFIT 6.1, also muss die Familienregel TM1 leicht "
          "brechen. G = v_geo^-2 folgt aus 8 pi = 1/c3.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z57_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False, default=str) + "\n")
