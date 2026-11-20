"""H9: Auswahl durch scharfe Schranken. Geprüft an den zwei Stellen, an denen der Korpus eine exakte Schranke hat."""
import numpy as np, sympy as sp
from scipy.linalg import expm, logm
from common import Checks
c = Checks("h09_scharf")
# --- A. gemeinsame Registerhierarchie (Gesamtsynthese 28.09.2026): S3 Gewichte w(q), Phi_q(X) = 2/3 X, Phi_q(Y) = q Y, Phi_q(Z) = 1/3 Z
X = np.array([[0, 1], [1, 0]], complex); Y = np.array([[0, -1j], [1j, 0]]); Z = np.diag([1, -1]).astype(complex); I2 = np.eye(2)
rotY = lambda t: np.cos(t / 2) * I2 - 1j * np.sin(t / 2) * Y
axis = lambda a: np.cos(a) * Z + np.sin(a) * X                 # Halbdrehung um eine Achse in der X Z Ebene
U = [I2, rotY(2 * np.pi / 3), rotY(-2 * np.pi / 3), axis(0), axis(2 * np.pi / 3), axis(-2 * np.pi / 3)]
w = lambda q: np.array([1 / 2 + q / 6, q / 6, q / 6, 1 / 18 - q / 6, 2 / 9 - q / 6, 2 / 9 - q / 6])
def chan(q, P): return sum(wi * u @ P @ u.conj().T for wi, u in zip(w(q), U))
okc = all(np.allclose(chan(q, X), 2 / 3 * X) and np.allclose(chan(q, Y), q * Y) and np.allclose(chan(q, Z), 1 / 3 * Z) for q in (0.1, 2 / 9, 0.24, 1 / 3))
c.ok("Rekonstruktion der sechs Wörter: der Kanal wirkt als X -> 2/3 X, Y -> q Y, Z -> 1/3 Z (wie in der Quelle)", okc)
q = sp.symbols("q", positive=True); r, s = sp.Rational(2, 3), sp.Rational(1, 3)
rates = {"e": (sp.log(q) + 2 * sp.log(r) + 2 * sp.log(s)) / 6, "Drehung": (sp.log(q) - sp.log(r) - sp.log(s)) / 6,
         "Spiegelung r59": (-sp.log(q) - 2 * sp.log(r) + 2 * sp.log(s)) / 6, "Spiegelung r14, r13": (-sp.log(q) + sp.log(r) - sp.log(s)) / 6}
lo = sp.solve(sp.Eq(rates["Drehung"], 0), q)[0]; hi = sp.solve(sp.Eq(rates["Spiegelung r59"], 0), q)[0]
c.ok("der Generator in der Gruppenalgebra hat genau für 2/9 <= q <= 1/4 nur nichtnegative Sprungraten (Reproduktion des Intervalls)", (lo, hi) == (sp.Rational(2, 9), sp.Rational(1, 4)), (lo, hi))
zero = lambda e: abs(float(sp.N(e, 30))) < 1e-20
pos = lambda e: float(sp.N(e, 30)) > 1e-12
def numeric_rates(qv):
    V = [np.kron(u, u) for u in U]; S = lambda M: np.kron(M, M.conj())                 # Superoperatoren auf zwei Registern
    Psi = sum(wi * S(v) for wi, v in zip(w(qv), V)); Lg = logm(Psi)
    A = np.array([S(v).flatten() for v in V]).T
    co = np.linalg.lstsq(A, Lg.flatten(), rcond=None)[0]
    return co.real, np.linalg.norm(A @ co - Lg.flatten())
ok_num = True
for qv in (2 / 9, 0.235, 0.25):
    co, res = numeric_rates(qv); ex = [float(rates[k].subs(q, qv)) for k in ("e", "Drehung", "Drehung", "Spiegelung r59", "Spiegelung r14, r13", "Spiegelung r14, r13")]
    ok_num &= res < 1e-8 and np.allclose(co, ex, atol=1e-8)
c.ok("Gegenprobe auf zwei Registern: der Matrixlogarithmus des Kanals hat genau diese sechs Koeffizienten", ok_num)
c.ok("untere Schranke scharf (q = 2/9): die Drehungen springen nicht mehr, es bleiben nur die drei Spiegelungen. Dann ist q = (2/3)(1/3), die dritte Rate ist das Produkt der beiden anderen",
     zero(rates["Drehung"].subs(q, lo)) and lo == r * s and all(pos(rates[k].subs(q, lo)) for k in ("Spiegelung r59", "Spiegelung r14, r13")))
c.ok("obere Schranke scharf (q = 1/4): die Spiegelung r59 springt nicht mehr, Drehungen und zwei Spiegelungen bleiben", zero(rates["Spiegelung r59"].subs(q, hi)) and pos(rates["Drehung"].subs(q, hi)) and pos(rates["Spiegelung r14, r13"].subs(q, hi)))
c.ok("der historische Punkt q = 1/3 liegt außerhalb (Rate von r59 negativ)", float(sp.N(rates["Spiegelung r59"].subs(q, sp.Rational(1, 3)))) < 0)
c.ok("EINDEUTIGKEIT FÄLLT: das Schärfeprinzip lässt zwei Punkte übrig, 2/9 und 1/4, nicht einen", lo != hi)
c.ok("ZUSATZ, der zwischen beiden wählt: die nativen Ereignisse des Originals sind Spiegelungen (60 Reflexionen erzeugen G31). Nur q = 2/9 hat ausschließlich Spiegelungssprünge", zero(rates["Drehung"].subs(q, lo)) and pos(rates["Drehung"].subs(q, hi)))
# --- B. Viererring (OPHCP20261119): gamma >= m, erreicht vom Generator mit Sprüngen J, J+, B, B+
E = lambda i, j: np.outer(np.eye(4)[i], np.eye(4)[j]).astype(complex)
Jm = E(1, 0) + E(2, 1) + E(3, 2); B = E(0, 3); m = 1.0
def Lsup(ops, g):
    L = np.zeros((16, 16), complex)
    for A in ops:
        AdA = A.conj().T @ A
        L += g * (np.kron(A, A.conj()) - 0.5 * np.kron(AdA, np.eye(4)) - 0.5 * np.kron(np.eye(4), AdA.T))
    return L
L4 = Lsup([Jm, Jm.conj().T, B, B.conj().T], m / 2)
E03 = E(0, 3).flatten(); out = (L4 @ E03).reshape(4, 4)
gam = -(out[0, 3] / 1).real
c.ok("Viererring: der Übertrag E03 zerfällt genau mit gamma = m; die Schranke gamma >= m wird erreicht (Reproduktion)", abs(gam - m) < 1e-12 and np.allclose(out, -m * E(0, 3)), gam)
pops = expm(L4 * 0.7) @ E(0, 0).flatten(); p00 = pops.reshape(4, 4)[0, 0].real
c.ok("dieselbe Entwicklung mischt die Besetzungen wie verlangt: p00(t) = (1 + 2 e^(-mt) + e^(-2mt))/4", abs(p00 - (1 + 2 * np.exp(-0.7) + np.exp(-1.4)) / 4) < 1e-10, p00)
c.ok("NICHT VERBUNDEN: A und B sind verschiedene Modelle (zwei Zustände mit S3, vier Stationen mit Zeitkovarianz). Es gibt im Korpus keine Abbildung, die gamma = m in einen Wert von q übersetzt", True)
c.ok("NICHT PRÜFBAR: für kappa der Kugelrandantwort und für das Seitenverhältnis des Torus nennt der Korpus keine Schranke, die scharf werden könnte", True)
c.done({"intervall": [str(lo), str(hi)], "raten_bei_2/9": {k: round(float(sp.N(v.subs(q, lo))), 12) for k, v in rates.items()}, "raten_bei_1/4": {k: round(float(sp.N(v.subs(q, hi))), 12) for k, v in rates.items()},
        "urteil": "Teilweise: scharfe Schranken geben Randpunkte (q = 2/9 oder 1/4; gamma = m wird erreicht), aber keinen eindeutigen Punkt, und zwei der drei freien Größen haben keine Schranke. Neu: q = 2/9 ist der einzige Punkt mit reinen Spiegelungssprüngen, und dort gilt q = (2/3)(1/3)."})
