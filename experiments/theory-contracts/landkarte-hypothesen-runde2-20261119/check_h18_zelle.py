"""H18: TFPT als Bruchteil einer OPH Zelle. Die OPH Netze und die Definition des Quotienten liegen nicht im Repository."""
import numpy as np, sympy as sp
from fractions import Fraction as Fr
from collections import Counter
from common import Checks, e8
c = Checks("h18_zelle")
Q = lambda rates: Fr(len(rates)) * sum(r * r for r in rates) / sum(rates) ** 2            # hängt nicht von der Zeiteinheit ab
c.ok("REKONSTRUKTION: mit Q = (Zahl der Moden) x Summe r^2 / (Summe r)^2 geben die Zerfallsraten 1/3, 1/3, 2/3 genau 9/8 = 1,125. Das ist meine Lesart der berichteten Größe, nicht deren hinterlegte Definition", Q([Fr(1, 3), Fr(1, 3), Fr(2, 3)]) == Fr(9, 8))
c.ok("die berichteten OPH Werte 7,6 · 7,9 · 7,975 sind exakt 8 - 96/N für N = 240, 960, 3840", [8 - Fr(96, N) for N in (240, 960, 3840)] == [Fr(38, 5), Fr(79, 10), Fr(319, 40)])
c.ok("für Q gilt immer Q >= 1, mit Gleichheit bei lauter gleichen Raten. 9/8 heißt fast gleichmäßig, 8 heißt: nur etwa ein Achtel der Moden trägt den Zerfall", Q([Fr(1)] * 5) == 1)
# Ersatzrechnung: die 240 Wurzeln als Zelle, Nachbarn = Wurzeln mit Skalarprodukt 1
al, roots, coef, Binv = e8(); R = np.array([[float(x) for x in r] for r in roots]); G = R @ R.T; A = (np.abs(G - 1) < 1e-9).astype(float)
deg = A.sum(axis=1); L = np.diag(deg) - A; ev = np.linalg.eigvalsh(L)
spec = Counter(int(round(x)) for x in ev)
c.ok("Wurzelgraph von E8: jede Wurzel hat 56 Nachbarn; Spektrum des Laplace 0 (1), 28 (8), 48 (35), 58 (112), 60 (84)", set(deg) == {56.0} and dict(spec) == {0: 1, 28: 8, 48: 35, 58: 112, 60: 84}, dict(spec))
qfull = Q([Fr(k) for k, m in spec.items() if k for _ in range(m)])
c.ok("auf diesem Graphen ist Q = 1,013, nicht 7,6. Die OPH Zelle ist also nicht der Wurzelgraph; ihre Kopplung muss stark ungleichmäßig sein", abs(float(qfull) - 1.013) < 2e-3, float(qfull))
# Zusammenfassen nach den Spiegelungen von D5 und A3
def cls(k, r):
    g = k[5] % 4
    if g == 0: return "D5" if (k[5], k[6], k[7]) == (0, 0, 0) else "A3"
    return {1: "16x4", 2: "10x6", 3: "16qx4q"}[g]
lab = [cls(k, r) for k, r in zip(coef, roots)]; names = ["D5", "A3", "16x4", "16qx4q", "10x6"]; cnt = Counter(lab)
c.ok("die Spiegelungen von D5 und A3 zerlegen die 240 Ports in FÜNF Klassen (40, 12, 64, 64, 60), nicht vier: Spinor und Gegenspinor bleiben getrennt", [cnt[n] for n in names] == [40, 12, 64, 64, 60], dict(cnt))
idx = {n: [i for i, x in enumerate(lab) if x == n] for n in names}
B = [[A[np.ix_(idx[a], idx[b])].sum(axis=1) for b in names] for a in names]
equit = all(len(set(np.round(B[i][j]).astype(int))) == 1 for i in range(5) for j in range(5)); Bq = np.array([[B[i][j][0] for j in range(5)] for i in range(5)])
c.ok("die Zusammenfassung ist verträglich: jede Wurzel einer Klasse hat gleich viele Nachbarn in jeder anderen Klasse", equit, Bq.astype(int).tolist())
evq = sorted(np.round(np.linalg.eigvals(56 * np.eye(5) - Bq).real, 6)); nz = [Fr(int(round(x))) for x in evq if abs(x) > 1e-6]
qq = Q(nz)
c.ok("das zusammengefasste Netz hat vier Zerfallsraten aus dem Spektrum des Graphen und ein Q, das die Vorhersage 9/8 nicht trifft", abs(float(qq) - 1.125) > 0.02, ([int(x) for x in nz], float(qq)))
c.ok("NICHT ENTSCHEIDBAR: der entscheidende Test braucht das veröffentlichte OPH Netz mit seiner eigenen Kopplung und die Definition des Quotienten. Beides liegt nicht im Repository", True)
c.ok("WAS BLEIBT: Zusammenfassen von Ports macht ein Netz gleichmäßiger und senkt Q. Die Richtung der Hypothese (TFPT kleiner als eine Zelle) ist mit 9/8 < 4 verträglich; gezeigt ist sie nicht", True)
c.done({"Q_tfpt": "9/8", "Q_wurzelgraph": float(qfull), "klassen": dict(cnt), "raten_zusammengefasst": [int(x) for x in nz], "Q_zusammengefasst": float(qq),
        "urteil": "Nicht entscheidbar: Daten fehlen. Die Ersatzrechnung am Wurzelgraphen trifft weder 7,6 noch 9/8. Korrektur: es sind fünf Klassen, nicht vier."})
