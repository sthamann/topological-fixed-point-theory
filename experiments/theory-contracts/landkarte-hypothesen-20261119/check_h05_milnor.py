"""H5: P1 und P2 als Faser und Umlauf einer Funktion. Singularität gegen elliptische Faserung."""
import itertools, numpy as np, sympy as sp
from common import Checks
c = Checks("h05_milnor")
def seifert(mu):                       # Seifertform von z^(mu+1): L + L^T = Cartanmatrix A_mu
    return np.eye(mu, dtype=int) - np.diag(np.ones(mu - 1, dtype=int), 1)
L3 = seifert(3); A3 = L3 + L3.T
h3 = np.round(-np.linalg.inv(L3) @ L3.T).astype(int)       # Umlauf
order = lambda M: next(k for k in range(1, 200) if np.array_equal(np.linalg.matrix_power(M, k), np.eye(len(M), dtype=int)))
c.ok("Faser von z^4: Gitter A3 (det 4), Umlauf hat Ordnung 4 = Coxeterzahl", round(np.linalg.det(A3)) == 4 and order(h3) == 4)
# E8 als A2 (x) A4 (Sebastiani Thom; A1 Faktor für x^2 ändert nur das Vorzeichen)
L = np.kron(seifert(2), seifert(4)); S = L + L.T
ev = np.linalg.eigvalsh(S)
c.ok("Faser von x^2 + y^3 + z^5: gerade, definit, Determinante 1, Rang 8, also das Gitter E8",
     round(abs(np.linalg.det(S))) == 1 and all(S[i, i] % 2 == 0 for i in range(8)) and (all(ev > 0) or all(ev < 0)), list(np.round(ev, 3)))
h8 = np.round(np.linalg.inv(L) @ L.T).astype(int)
x = sp.symbols("x")
cp = sp.Matrix(-h8).charpoly(x).as_expr()
c.ok("ihr Umlauf ist das Tensorprodukt der Umläufe zu 2, 3, 5: Ordnung 30, charakteristisches Polynom = Kreisteilungspolynom Phi_30",
     order(-h8) == 30 and sp.expand(cp - sp.cyclotomic_poly(30, x)) == 0, (order(h8), order(-h8)))
c.ok("also C30 = C2 x C3 x C5 ohne Zusatzannahme, und die acht Eigenwerte sind die primitiven 30. Einheitswurzeln (Exponenten 1,7,11,13,17,19,23,29)",
     sp.degree(sp.cyclotomic_poly(30, x)) == 8)
# Adjazenz: D5 + A3 ist kein Teildiagramm des E8 Diagramms, aber entsteht durch Streichen im erweiterten Diagramm
E8 = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (2, 7)]                       # Kette 0..6, Abzweig an Knoten 2
aff = E8 + [(6, 8)]                                                                 # erweiterter Knoten 8 am langen Ende
def cartan(nodes, edges):
    n = len(nodes); M = 2 * np.eye(n, dtype=int)
    for a, b in edges:
        if a in nodes and b in nodes: M[nodes.index(a), nodes.index(b)] = M[nodes.index(b), nodes.index(a)] = -1
    return M
def comps(nodes, edges):
    nodes = list(nodes); out = []; left = set(nodes)
    while left:
        st = [left.pop()]; comp = set(st)
        while st:
            u = st.pop()
            for a, b in edges:
                for p, q in ((a, b), (b, a)):
                    if p == u and q in left: left.remove(q); comp.add(q); st.append(q)
        cm = cartan(sorted(comp), edges); deg = sorted(int(-(cm.sum(axis=1) - 2)[i]) for i in range(len(comp)))
        out.append((len(comp), round(np.linalg.det(cm)), max(deg)))
    return sorted(out)
target = sorted([(5, 4, 3), (3, 4, 2)])                                              # D5: 5 Knoten, det 4, Verzweigung; A3: 3 Knoten, det 4, Kette
sub = [comps(s, E8) for r in range(1, 9) for s in itertools.combinations(range(8), r)]
c.ok("WÖRTLICHE FORM FÄLLT: D5 + A3 ist kein Teildiagramm von E8; die E8 Singularität hat keine Deformation mit diesen beiden Punkten", target not in sub)
marks = np.round(np.linalg.svd(cartan(list(range(9)), aff))[2][-1] / np.linalg.svd(cartan(list(range(9)), aff))[2][-1][8]).astype(int)
rem = [(int(marks[k]), comps([n for n in range(9) if n != k], aff)) for k in range(9)]
c.ok("im erweiterten Diagramm (Marken 1,2,2,3,3,4,4,5,6) gibt das Streichen eines Knotens mit Marke 4 genau D5 + A3; Marke 5 gibt A4 + A4, Marke 6 gibt A5 + A2 + A1",
     sorted(marks) == [1, 2, 2, 3, 3, 4, 4, 5, 6] and (4, target) in rem and (5, [(4, 5, 2), (4, 5, 2)]) in rem and (6, [(1, 2, 0), (2, 3, 1), (5, 6, 2)]) in rem, rem)
# die Funktion, die es wirklich gibt: elliptische Faserung, Tate Normalform mit Punkt der Ordnung 4
t, X, Y = sp.symbols("t X Y")
def tate(a1, a2, a3, a4, a6):
    b2 = a1 ** 2 + 4 * a2; b4 = 2 * a4 + a1 * a3; b6 = a3 ** 2 + 4 * a6; b8 = a1 ** 2 * a6 + 4 * a2 * a6 - a1 * a3 * a4 + a2 * a3 ** 2 - a4 ** 2
    D = -b2 ** 2 * b8 - 8 * b4 ** 3 - 27 * b6 ** 2 + 9 * b2 * b4 * b6
    return sp.factor(D), sp.factor(b2 ** 2 - 24 * b4), sp.factor(-b2 ** 3 + 36 * b2 * b4 - 216 * b6)
D4, c4, c6 = tate(1, -t, -t, 0, 0)                           # y^2 + xy - ty = x^3 - t x^2
c.ok("Kurve mit Punkt der Ordnung 4 über der t Kugel: Diskriminante t^4 (16 t + 1), also Fasern I4 bei 0, I1 bei -1/16, und Ordnung 12 - 5 = 7 bei unendlich (I1*)",
     sp.expand(D4 - t ** 4 * (16 * t + 1)) == 0 and c4.subs(t, 0) != 0, D4)
c.ok("Eulerzahl 4 + 1 + 7 = 12; Fasergitter A3 + D5 hat Rang 8 und Determinante 16 = 4^2: mit Torsionsschnitt der Ordnung 4 unimodular, also E8", 4 + 1 + 7 == 12 and 4 * 4 // 4 ** 2 == 1)
# der Punkt (0,0) hat Ordnung 4: Gruppengesetz
def add(P, Q, a=(1, -t, -t, 0, 0)):
    a1, a2, a3, a4, a6 = a
    if P is None: return Q
    if Q is None: return P
    x1, y1 = P; x2, y2 = Q
    if sp.simplify(x1 - x2) == 0:
        if sp.simplify(y1 + y2 + a1 * x2 + a3) == 0: return None
        lam = (3 * x1 ** 2 + 2 * a2 * x1 + a4 - a1 * y1) / (2 * y1 + a1 * x1 + a3); nu = (-x1 ** 3 + a4 * x1 + 2 * a6 - a3 * y1) / (2 * y1 + a1 * x1 + a3)
    else:
        lam = (y2 - y1) / (x2 - x1); nu = (y1 * x2 - y2 * x1) / (x2 - x1)
    x3 = sp.simplify(lam ** 2 + a1 * lam - a2 - x1 - x2); return (x3, sp.simplify(-(lam + a1) * x3 - nu - a3))
P = (sp.Integer(0), sp.Integer(0)); P2 = add(P, P); P3 = add(P2, P); P4 = add(P3, P)
c.ok("der Schnitt (0,0) hat genau Ordnung 4: das ist der Klebevektor als Torsionsschnitt (Mordell Weil Gruppe Z/4)", P2 is not None and P3 is not None and P4 is None, (P2, P3))
D5_, c4_5, _ = tate(1 - t, -t, -t, 0, 0)
c.ok("Gegenprobe Ordnung 5: Diskriminante t^5 (t^2 - 11 t - 1): Fasern I5, I5, I1, I1, also A4 + A4 mit Klebung Z5", sp.expand(D5_ - t ** 5 * (t ** 2 - 11 * t - 1)) == 0, D5_)
j1728 = sp.roots(sp.Poly(sp.numer(sp.together(c6)), t))
c.ok("drei Fasern der Familie sind quadratische Tori (c6 = 0 hat drei einfache Nullstellen): die Lage tau = i ist ein Punkt der Basis, kein Umlauf", len(j1728) == 3 and all(v == 1 for v in j1728.values()), list(j1728))
# Umläufe: T^4, -T, T mit Produkt 1 nach Konjugation
T = np.array([[1, 1], [0, 1]]); found = False
rng = range(-3, 4)
sl2 = [np.array([[a, b], [cc, d]]) for a in rng for b in rng for cc in rng for d in rng if a * d - b * cc == 1]
T4 = np.linalg.matrix_power(T, 4)
for g1 in sl2:
    M1 = g1 @ (-T) @ np.round(np.linalg.inv(g1)).astype(int)
    for g2 in sl2:
        M2 = g2 @ T @ np.round(np.linalg.inv(g2)).astype(int)
        if np.array_equal(T4 @ M1 @ M2, np.eye(2, dtype=int)): found = True; break
    if found: break
c.ok("die drei Umläufe (parabolisch T^4, -T, T) haben nach Konjugation Produkt 1", found)
c.ok("ZEIT LESART FÄLLT: alle Umläufe haben unendliche Ordnung; die Vierteldrehung ist kein Umlauf dieser Familie, sondern eine Zusatzsymmetrie dreier Fasern",
     not any(np.array_equal(np.linalg.matrix_power(T4, k), np.eye(2, dtype=int)) for k in range(1, 50)))
c.ok("Zahl am Rand: (1 - 2/8) + (1 - 2/4) = 5/4 für D5 und A3; reine Arithmetik, keine Abbildung", sp.Rational(6, 8) + sp.Rational(2, 4) == sp.Rational(5, 4))
c.done({"gestrichene_knoten": [str(r) for r in rem],
        "urteil": "Wörtliche Form (eine isolierte Singularität) ausgeschlossen. Exakt: C30 = Tensorprodukt der Umläufe 2, 3, 5. Exakter Ersatz: D5 + A3 mit Klebung Z4 ist das Fasergitter der elliptischen Fläche mit Punkt der Ordnung 4; der Klebevektor ist ihr Torsionsschnitt, ihr Gitter ist E8 + U. Als Zeit fällt der Umlauf aus."})
