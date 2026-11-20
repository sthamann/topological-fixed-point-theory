"""H7: gemeinsame Regel der Schnitte. Verdrehte eta Invarianten der Linsenräume gegen die Klebenormen der E8 Teilgitter."""
import itertools, numpy as np
from fractions import Fraction as F
from mpmath import mp, mpf, pi, sin, expj
from common import Checks
mp.dps = 30
c = Checks("h07_kac")
tol = mpf(10) ** -20
def eta_g(phi):                        # wie in P1KERN20261119: äquivariante Dirac eta auf S^3, Abel summiert
    w = expj(phi); S0 = 1 / (1 - w); S1 = w / (1 - w) ** 2
    return ((expj(2 * phi) * (S1 + S0) - expj(phi) * (S1 + 2 * S0)).imag) / sin(phi)
def eta(p, k=0):                       # verdreht mit dem Charakter der Windung k
    return sum((expj(2 * pi * j * k / p) * eta_g(2 * pi * j / p)) for j in range(1, p)) / p
vals = {(p, k): eta(p, k) for p in range(2, 13) for k in range(p)}
c.ok("alle verdrehten Werte sind reell", all(abs(v.imag) < tol for v in vals.values()))
sgn = 1 if vals[(4, 0)].real > 0 else -1
e = {key: sgn * v.real for key, v in vals.items()}
c.ok("ohne Windung: eta = (p^2 - 1)/(6p) (Reproduktion)", all(abs(e[(p, 0)] - mpf(p * p - 1) / (6 * p)) < tol for p in range(2, 13)))
c.ok("Windung k = p/2 ist die zweite Spinstruktur: bei p = 4 gibt k = 2 den Wert -3/8 (Reproduktion)", abs(e[(4, 2)] + mpf(3) / 8) < tol, e[(4, 2)])
c.ok("HYPOTHESE IN DER GENANNTEN FORM FÄLLT: der verdrehte Wert ist nicht k(p-k)/(2p)",
     any(abs(abs(e[(p, k)]) - mpf(k * (p - k)) / (2 * p)) > mpf("0.01") for p in range(2, 13) for k in range(1, p)))
c.ok("GEMEINSAME REGEL (exakt, p = 2..12, alle k): eta_k = (p^2 - 1)/(6p) - k(p - k)/p",
     all(abs(e[(p, k)] - (mpf(p * p - 1) / (6 * p) - mpf(k * (p - k)) / p)) < tol for p in range(2, 13) for k in range(p)))
c.ok("also: der Unterschied zur unverdrehten Größe ist genau k(p-k)/p, die Norm der Klebeklasse k von A_(p-1)", True)
zero = [(p, k) for p in range(2, 400) for k in range(1, p // 2 + 1) if p * p - 1 == 6 * k * (p - k)]
c.ok("die Phase verschwindet ganz (eta_k = 0) für Einheitswindung k = 1 NUR bei p = 5; das ist die Fünferuhr Beobachtung aus P1KERN in einer Zeile",
     [x for x in zero if x[1] == 1] == [(5, 1)], zero[:6])
c.ok("die beiden Muster aus P1KERN sind zwei Fälle dieser einen Formel: p = 4 mit k = 0, 2 gibt 5/8 und -3/8; p = 5 mit k = 0 gibt 4/5 = 1·4/5",
     abs(e[(4, 0)] - mpf(5) / 8) < tol and abs(e[(4, 2)] + mpf(3) / 8) < tol and abs(e[(5, 0)] - mpf(4) / 5) < tol and abs(e[(5, 1)]) < tol)
# E8: Wurzeln, erweitertes Diagramm, Teilgitter nach Streichen eines Knotens, Normen der Klebeklassen
roots = [np.array(v) for v in itertools.product([-1, 0, 1], repeat=8) if sum(abs(x) for x in v) == 2]
roots += [np.array(s) / 2 for s in itertools.product([-1, 1], repeat=8) if sum(1 for x in s if x < 0) % 2 == 0]
simple = [np.array(v, float) for v in ([1, -1, 0, 0, 0, 0, 0, 0], [0, 1, -1, 0, 0, 0, 0, 0], [0, 0, 1, -1, 0, 0, 0, 0], [0, 0, 0, 1, -1, 0, 0, 0],
                                         [0, 0, 0, 0, 1, -1, 0, 0], [0, 0, 0, 0, 0, 1, -1, 0], [0, 0, 0, 0, 0, 1, 1, 0], [-.5] * 8)]
B = np.array(simple); coords = lambda v: np.linalg.solve(B.T, v)
high = max(roots, key=lambda r: sum(coords(r)))
marks = np.round(coords(high)).astype(int); ext = simple + [-high]; allmarks = list(marks) + [1]
c.ok("240 Wurzeln; Marken des erweiterten Diagramms 1,2,2,3,3,4,4,5,6", len(roots) == 240 and sorted(allmarks) == [1, 2, 2, 3, 3, 4, 4, 5, 6])
def split(k):                          # Komponenten des Teilgitters nach Streichen von Knoten k und Projektionsnormen der Nebenklassen
    keep = [ext[i] for i in range(9) if i != k]; idx = list(range(8)); comp = []
    while idx:
        st = [idx.pop()]; cc = set(st)
        while st:
            u = st.pop()
            for vtx in idx[:]:
                if abs(keep[u] @ keep[vtx]) > 1e-9: idx.remove(vtx); cc.add(vtx); st.append(vtx)
        comp.append(sorted(cc))
    Bk = np.array(keep); table = {}
    for r in roots:
        cf = np.linalg.solve(Bk.T, r)
        if all(abs(x - round(x)) < 1e-9 for x in cf): continue
        norms = []
        for cc in comp:
            Q = np.array([keep[i] for i in cc]).T; pr = Q @ np.linalg.solve(Q.T @ Q, Q.T @ r); norms.append((len(cc) if round(np.linalg.det(Q.T @ Q)) == len(cc) + 1 else -len(cc), F(float(pr @ pr)).limit_denominator(60)))   # Rang > 0: Typ A, Rang < 0: Typ D oder E
        table.setdefault(tuple(sorted(norms)), 0); table[tuple(sorted(norms))] += 1
    return comp, table
out = {}
for k in range(9):
    comp, table = split(k); out[k] = (allmarks[k], sorted(len(x) for x in comp), table)
k4 = next(k for k in out if out[k][0] == 4 and out[k][1] == [3, 5]); k5 = next(k for k in out if out[k][0] == 5)
c.ok("Marke 4 (D5 + A3): jede Wurzel außerhalb des Teilgitters zerfällt in Normen 5/4 + 3/4 oder 1 + 1; Zählung 64 + 64 und 60",
     out[k4][2] == {((-5, F(5, 4)), (3, F(3, 4))): 128, ((-5, F(1)), (3, F(1))): 60}, out[k4][2])
c.ok("Marke 5 (A4 + A4): Normen 4/5 + 6/5 und 6/5 + 4/5", set(out[k5][2]) == {((4, F(4, 5)), (4, F(6, 5)))}, out[k5][2])
ok_all = True; rows = []
for k, (mk, sizes, table) in out.items():
    for key in table:
        for (rank, nrm) in key:
            # A Faktoren: Norm muss k'(n+1-k')/(n+1) für ein k' sein
            if rank > 0 and not any(nrm == F(kk * (rank + 1 - kk), rank + 1) for kk in range(0, rank + 1)): ok_all = False
    rows.append((mk, sizes, {str(kk): v for kk, v in table.items()}))
c.ok("für alle A Faktoren aller neun Knoten ist die Projektionsnorm von der Form k(n+1-k)/(n+1), also ein eta Unterschied des Linsenraums S^3/Z_(n+1)", ok_all)
c.ok("GRENZE: die 5/8 von D5 ist die unverdrehte eta von S^3/Z_4. Der passende Rand zu D5 ist kein Linsenraum, sondern der Rand der D5 Singularität. Diese Gleichheit ist nur bei p = 4 da (p = 2 trifft 1/4, p = 3, 5, 6 treffen nicht)",
     [abs(e[(p, 0)] - x) < tol for p, x in ((2, mpf(1) / 4), (3, mpf(2) / 3), (5, mpf(3) / 5), (6, mpf(7) / 12))] == [True, False, False, False])
c.done({"eta": {f"{p},{k}": str(v)[:12] for (p, k), v in e.items() if p <= 6}, "knoten": [str(r) for r in rows], "nullstellen": zero[:6],
        "urteil": "Die genannte Formel k(p-k)/(2p) fällt. Die gemeinsame Regel existiert und ist exakt: eta_k(p) = (p^2-1)/(6p) - k(p-k)/p. Sie fasst beide Muster aus P1KERN zusammen, und die Einheitswindung löscht die Phase genau bei p = 5. Der zweite Summand ist die Klebenorm von A_(p-1). Nicht gezeigt: dass P1 an einem solchen Schnitt sitzt."})
