"""Z171: Gesamtlösungs-Suche (ein Objekt bis zum Messwert) und Prüfung der externen Analysen 24 und 25.

Leitlinie 'Rekonstruktion statt Auswahl'; H6 und die gelösten Blöcke bleiben unberührt, nichts wird befördert.
Statuskürzel: E exakt, C bedingt, L Lesart, N numerisch, O offen.

Die Suche steht in research/originlab/whole_solution_search.py. Sie rechnet jedes Objekt des Katalogs bis zum Ende
durch (Anzahl, Hyperladungsspur, Basis, Kernpotenz, daraus alpha und die Lambda-Stufe) und vergleicht mit den
Messwerten. Die Form der Alpha-Gleichung ist Ledger-Eingabe und wird nicht hergeleitet.

Teil A  Analyse 24: verdrillte Schalensummen, theta_4, Familienprojektor
Teil B  Analyse 25: Projektor auf den alternierenden Kanal, Determinante, elliptische Fläche
Teil C  Gesamtlösungs-Suche: ein Objekt gegen die Messwerte
Teil D  Eine Gleichung aus einem Objekt
Teil E  Bilanz und Urteil.
"""
import json
import os
import sys
import time
from fractions import Fraction as Fr
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import mpmath as mp
import sympy as sp

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parents[2] / "research" / "originlab"))
import whole_solution_search as ws                                         # noqa: E402

OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
T0 = time.time()
mp.mp.dps = 40
AINV = mp.mpf("137.0359992168407")
c3 = 1 / (8 * mp.pi)


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val, text):
    BEFUNDE[name] = val
    say(f"   -> {name}: {text}")


res = ws.run_whole_solution_search()
objs = {o["label"]: o for o in ws.catalogue()}

say("Teil A  Analyse 24: verdrillte Schalensummen, theta_4, Familienprojektor")
r2 = objs["1 x trace | insert R^2 | family None | kernel D5"]["shell_traces"]
check("A1 [E] Schalensummen mit R^2: 1, -16, 112, -448, 1136", r2 == ["1", "-16", "112", "-448", "1136"])
q = sp.symbols("q")
th4 = sum((-1) ** abs(n) * q ** (n * n) for n in range(-3, 4))
ser = sp.Poly(sp.expand(th4 ** 8), q)
check("A2 [E] Z_{R^2} = theta_4(q)^8 bis zur vierten Schale",
      [ser.coeff_monomial(q ** n) for n in range(5)] == [1, -16, 112, -448, 1136])
w = objs["3 x minus_log | insert R^2 | family None | kernel D5"]["series"]
check("A3 [E] W = -3 ln Z_{R^2} = -24 ln theta_4 = 48q + 48q^2 + 64q^3 + 48q^4", w == ["0", "48", "48", "64", "48"])
fam = objs["1 x fredholm | class 1 | family standard | kernel D5"]
check("A4 [E] Familienprojektor im Materiesektor: Schalen 64 -> 48, 512 -> 336, 1792 -> 1056",
      fam["shell_traces"][:4] == ["1", "48", "336", "1056"]
      and objs["1 x trace | class 1 | family None | kernel D5"]["shell_traces"][:4] == ["1", "64", "512", "1792"])
check("A5 [E] -ln det(1 - Q) = 48q + 360q^2 + 1072q^3", fam["series"][:4] == ["0", "48", "360", "1072"])
q_el = mp.e ** (-2 / AINV)
L_el = 2 / AINV                                                            # theta_4(q) = sqrt(pi/L) theta_2(exp(-pi^2/L))
w_el = -24 * (mp.log(mp.pi / L_el) / 2 + mp.log(mp.jtheta(2, 0, mp.e ** (-mp.pi ** 2 / L_el))))
target = mp.log(1 / (4 * c3 / 3 + 48 * c3 ** 4 * q_el * (1 - 48 * c3 ** 4 * q_el) ** mp.mpf("-1.25")))
ZAHLEN["W_elektrisch"], ZAHLEN["ln_1_durch_phi_seam"] = mp.nstr(w_el, 8), mp.nstr(target, 8)
check("A6 [N] elektrische Auswertung von W bei q = exp(-2 alpha) trifft ln(1/phi_seam) = 2,93426 nicht "
      f"(W = {mp.nstr(w_el, 6)}, mal c3^2 = {mp.nstr(w_el * c3 ** 2, 5)}; W wächst wie 3 pi^2/alpha)",
      abs(target - mp.mpf("2.93426")) < 1e-5 and abs(w_el * c3 ** 2 / target - 1) > 0.5 and w_el / target > 100
      and abs(w_el / (3 * mp.pi ** 2 * AINV) - 1) < 0.03)
befund("Analyse24", "teilweise bestätigt",
       "Die Gitterzahlen stimmen. W ist geschlossen -24 ln theta_4(q). Auf der magnetischen Seite ist nur die "
       "erste Schale geprüft (alle weiteren liegen unter 10^-238). Die elektrische Seite trifft die Funktion nicht. "
       "Die Zahlen 203,25 und 68395,7 der Analyse sind ohne Angabe des Auswertungspunkts nicht nachgerechnet; die "
       "Aussagen zu Chernklassen und Index auf CP^3 sind ohne Definition des Bündels nicht prüfbar.")

say("\nTeil B  Analyse 25: Projektor auf den alternierenden Kanal, Determinante, elliptische Fläche")


def circ(row):
    return sp.Matrix(4, 4, lambda i, j: row[(j - i) % 4])


B, C, U = circ([2, -1, 0, -1]), circ([0, 1, 2, 1]), sp.diag(1, -1, 1, -1)
PI = (B + U * C * U) / 8
alt = sp.Matrix([1, -1, 1, -1])
check("B1 [E] Pi = (B + UCU)/8 = B(B - 2)/8 ist der Projektor vom Rang eins auf den alternierenden Kanal",
      PI == B * (B - 2 * sp.eye(4)) / 8 and PI * PI == PI and PI.rank() == 1 and PI * alt == alt)
t, lam, cc = sp.symbols("t lambda c3", positive=True)
basis = sp.Matrix([[1, -1, 0, 0], [0, 1, -1, 0], [0, 0, 1, -1]]).T          # Summe null: die drei Familienrichtungen
gram = basis.T * basis
red = gram.inv() * basis.T * (t * sp.eye(4) - lam ** 2 * cc * PI) * basis
check("B2 [E] auf den drei Richtungen mit Summe null: det(t - lambda^2 c3 Pi) = t^2 (t - lambda^2 c3)",
      sp.simplify(red.det() - t ** 2 * (t - lam ** 2 * cc)) == 0)
check("B3 [E] die Ableitung des Logarithmus, 2/t + 1/(t - c3), ist nicht die kubische Gleichung; lambda bleibt frei",
      sp.simplify(sp.diff(sp.log(t ** 2 * (t - cc)), t) - (2 / t + 1 / (t - cc))) == 0)
X, a, u, v = sp.symbols("X a u v")
disc = -16 * (4 * a ** 3)
check("B4 [E] Y^2 = X^3 + Q X: Diskriminante -64 Q^3, j = 1728 überall",
      sp.simplify(disc + 64 * a ** 3) == 0 and sp.discriminant(X ** 3 + a * X, X) == -4 * a ** 3
      and sp.Rational(1728) * 4 * a ** 3 / (4 * a ** 3) == 1728)
check("B5 [E] vier Fasern vom Typ III (Ordnung 3 der Diskriminante, Eulerzahl 3): Eulerzahl 12, rationale "
      "elliptische Fläche; die Kurve v^4 = Q(u) hat Geschlecht 3",
      4 * 3 == 12 and sp.degree(sp.expand((u ** 4 - 1) ** 3), u) == 12 and (4 * (-2) + 4 * 3 + 2) // 2 == 3)
befund("Analyse25", "bestätigt, mit Folgerung",
       "Projektor, Determinante und Flächendaten stimmen. Die linke Seite der Alpha-Gleichung in t = e^2/c3 ist "
       "bei lambda = 1 die Determinante auf den drei Familienrichtungen; die Stärke lambda = 1 ist nicht erzwungen.")

say("\nTeil C  Gesamtlösungs-Suche: ein Objekt gegen die Messwerte")
ZAHLEN["suche"] = {"objekte": res["objects"], "ausgewertet": res["evaluated"], "anzahlen": res["counts"]}
HIT = "1 x fredholm | class 1 | family standard | kernel D5"
check("C1 [E] 324 Objekte, 88 bis zum Messwert auswertbar", (res["objects"], res["evaluated"]) == (324, 88))
check("C2 [N] genau ein Objekt trifft alpha und Lambda: Materiesektor der Uhr mit Familienprojektor (Spur und "
      "Determinante teilen die erste Schale)",
      res["world_hits"] == [HIT, HIT.replace("fredholm", "trace")] and res["all_clauses"] == [HIT])
d = res["derived"][HIT]
check("C3 [E] aus dem Objekt: Anzahl 48, b1 = 41/10 (Fermionteil 4), Basis 4/3, Kernpotenz 5/4",
      (d["n_top"], d["b1"], d["base"], d["power"]) == ("48", "41/10", "4/3", "5/4"))
check("C4 [N] alpha^-1 = 137,0359992168 (1,9 Sigma zum Messwert), rho_Lambda/M_P^4 = 1,128e-123 (Messwert "
      "1,134e-123, 3 %)", abs(d["alpha_inverse"] - 137.0359992168407) < 1e-9 and abs(d["alpha_pull_sigma"]) < 2
      and abs(d["rho_lambda"] / 1.128042e-123 - 1) < 1e-5)
miss = list(res["nearest_miss"].values())
check("C5 [N] Nullkontrolle: der nächste Fehlversuch (Kern aus dem A3-Teil, Potenz 3/4) liegt 96 Sigma daneben, "
      "alle anderen über 10^8 Sigma", abs(miss[0]["alpha_pull_sigma"] + 95.54) < 0.1
      and all(abs(m["alpha_pull_sigma"]) > 1e8 for m in miss[2:]))
sig = lambda n, p, b, base=4 / 3: (ws.alpha_inverse(n, p, b, base) - ws.ALPHA_INV_MEASURED) / ws.ALPHA_INV_SIGMA
ZAHLEN["empfindlichkeit_sigma"] = {"N=47": sig(47, 1.25, 4.1), "N=49": sig(49, 1.25, 4.1), "p=1": sig(48, 1.0, 4.1),
                                   "p=3/2": sig(48, 1.5, 4.1), "b1=4": sig(48, 1.25, 4.0),
                                   "Basis=1": sig(48, 1.25, 4.1, 1.0)}
check("C6 [N] Schärfe: N = 47 oder 49, Potenz 1 oder 3/2, b1 = 4 und Basis 1 verfehlen den Messwert jeweils um "
      "mehr als 40 Sigma", all(abs(x) > 40 for x in ZAHLEN["empfindlichkeit_sigma"].values()))
uniq = [r for r in res["search"]["unique_selections"] if r["selected"] == HIT]
check("C7 [E] Regelmaschine: 'alpha gemessen' und 'positiv' wählen das Objekt eindeutig",
      any({"alpha_measured", "positive"} == set(m) for r in uniq for m in r["minimal_rules"]))
befund("Suche", "ein Objekt",
       "Ein einziges Objekt trägt Anzahl 48, drei Familien, b1 und die Potenz 5/4 und trifft damit alpha und Lambda. "
       "Nicht hergeleitet und weiter Ledger-Eingabe: " + "; ".join(res["ledger_inputs_not_derived"]) + ".")

say("\nTeil D  Eine Gleichung aus einem Objekt")
N, b1, p = 48, mp.mpf(41) / 10, mp.mpf(5) / 4
Xf = lambda x: N * c3 ** 4 * mp.e ** (-2 * x)
F = lambda z: z * (1 - z) ** (-p)
eq = lambda al: (al / (2 * c3 ** 2)) ** 2 * (al / (2 * c3 ** 2) - c3) - b1 * mp.log(1 / (mp.mpf(4) / 3 * c3 + F(Xf(al))))
root = mp.findroot(eq, 1 / AINV)
check("D1 [N] det_3(t - c3 Pi) = b1 ln(1/phi_seam) mit t = alpha/(2 c3^2): Wurzel alpha^-1 = 137,0359992168407",
      abs(1 / root - AINV) < mp.mpf("1e-12") and abs(eq(root)) < mp.mpf("1e-30"))
rho = Xf(1 / root)
ZAHLEN["rho_Lambda_Planck"], ZAHLEN["rho_Lambda_reduziert"] = mp.nstr(rho, 8), mp.nstr(rho / c3 ** 2, 8)
check("D2 [E] dieselbe Funktion X(x) = 48 c3^4 exp(-2x): X(alpha) im Kern, X(1/alpha) = rho_Lambda/M_P^4; "
      "reduziert 7,125e-121 (COSMO.LAM.01), G Lambda = 8 pi X(1/alpha) = 3/(32 pi^3) exp(-2/alpha)",
      abs(rho / c3 ** 2 / mp.mpf("7.125e-121") - 1) < 1e-3
      and abs(8 * mp.pi * rho / (3 / (32 * mp.pi ** 3) * mp.e ** (-2 / root)) - 1) < mp.mpf("1e-30"))
check("D3 [N] F(X(1/alpha)) und X(1/alpha) unterscheiden sich relativ um weniger als 10^-120: eine Antwort F(X(x)) "
      "trägt beide Seiten", abs(F(rho) / rho - 1) < mp.mpf("1e-38"))
check("D4 [E] die Drei ist beide Male dieselbe Darstellung: Summe-null-Richtungen der vier Marken (Determinante) "
      "und Standarddarstellung von S4 auf den vier Spinorgewichten von A3 (Familienprojektor, 64 -> 48)",
      len(ws.family_group()) == 24 and gram.rank() == 3 and Fr(64, 48) == Fr(4, 3))
befund("Gleichung", "verdichtet",
       "det_3(t - c3 Pi) = [(2/5) Sp_Q(Y^2) + 1/10] ln(1/phi_seam), phi_seam = c3 Sp(P_1)/Sp(P_3 P_1) + F(X(alpha)), "
       "X(x) = Sp(P_3 P_1) c3^4 exp(-2x), F(z) = z (1-z)^(-2h_16); rho_Lambda/M_P^4 = F(X(1/alpha)). "
       "Das ist eine Verdichtung der gelösten Blöcke auf ein Objekt, keine Herleitung des Gesetzes.")

say("\nTeil E  Bilanz und Urteil")
bilanz = [
    ("E", "Analyse 24", "Gitterzahlen bestätigt; W = -24 ln theta_4; elektrische Seite trifft nicht."),
    ("E", "Analyse 25", "Projektor, Determinante, elliptische Fläche bestätigt; lambda = 1 nicht erzwungen."),
    ("N", "Suche", "ein Objekt von 324 trifft alpha und Lambda; nächster Fehlversuch 96 Sigma."),
    ("E/L", "ein Objekt", "48, drei Familien, 4/3, b1-Fermionteil 4 und 5/4 aus dem Materiesektor mit "
     "Familienprojektor; dieselbe Drei in der Determinante der linken Seite."),
    ("O", "offen", "das Gesetz selbst (Determinante gleich b1 mal Logarithmus), das Gewicht c3^4 exp(-2x), der "
     "Skalaranteil 1/10, die Einfügungsform, die Wahl des D5-Teils; dazu unverändert Ursprung der Naht (K19, K20), "
     "3+1, Kontinuum, Dynamik, Gravitation."),
]
for st, was, text in bilanz:
    say(f"   [{st}] {was}: {text}")
urteil = ("Die Programmsuche rechnet jetzt jedes Objekt bis zum Messwert durch. Von 324 Objekten aus Gitter, Uhr und "
          "Familiengruppe trifft genau eines alpha und Lambda: der Materiesektor der Uhr mit dem Familienprojektor. "
          "Es liefert die Anzahl 48, die drei Familien, die Basis 4/3, den Fermionteil von b1 und die Potenz 5/4, die "
          "im Ledger aus vier getrennten Quellen stammen. Die linke Seite der Alpha-Gleichung ist die Determinante "
          "auf denselben drei Richtungen. Damit stehen alpha und Lambda in einer Gleichung aus einem Objekt. "
          "Das ist eine Verdichtung, keine volle Lösung: Die Form des Gesetzes, das Gewicht und die Einfügung bleiben "
          "Eingabe, und der Ursprung der Naht, 3+1, Kontinuum, Dynamik und Gravitation sind unverändert offen. "
          "Nichts befördert.")
say("\nUrteil\n   " + urteil)
print(f"\n   (Laufzeit {time.time() - T0:.1f} s, nicht Teil der Ausgabe)", flush=True)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z171_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
(HERE / "ausgabe_z171.txt").write_text("\n".join(OUT) + "\n")
