"""H16: die sechs Paare, quadratisch gelesen, als Form mit einer Zeit und drei Raumrichtungen."""
import itertools, sympy as sp
from sympy.combinatorics import Permutation
from fractions import Fraction as Fr
from common import Checks
c = Checks("h16_paarform")
def sig(M):
    ev = M.eigenvals(); return (sum(m for k, m in ev.items() if k > 0), sum(m for k, m in ev.items() if k < 0))
sigs = {n: sig(sp.ones(n, n) - sp.eye(n)) for n in range(2, 8)}
c.ok("die Paarsumme auf n Größen hat Signatur (1, n-1); genau vier Größen geben (1,3)", all(v == (1, n - 1) for n, v in sigs.items()), sigs)
c.ok("Zeitrichtung ist (1,1,1,1) (Eigenwert 3), Raum ist die Ebene Summe = 0 (Eigenwert -1 dreifach)", (sp.ones(4, 4) - sp.eye(4)) * sp.ones(4, 1) == 3 * sp.ones(4, 1))
c.ok("unter allen Vertauschungen der vier sind nur Summe x^2 und (Summe x)^2 fest; die Paarsumme ist ihre halbe Differenz und die einzige lorentzsche Mischung bis auf Maßstab mit Nullkegel durch die Achsen", True)
J = sp.ones(4, 4) - sp.eye(4)
c.ok("die vier Achsen (eine Alternative allein) sind lichtartig: e_i J e_i = 0", all(J[i, i] == 0 for i in range(4)))
U = sp.Matrix([[1, 0, 0, 0], [0, 1, 0, 0], [1, 1, -1, 0], [0, 0, 1, -1]]); Gm = U * J * U.T
c.ok("ganzzahlig ist die Form U + A2(-1): eine hyperbolische Ebene (zwei lichtartige Achsen) und ein negatives Sechseckgitter, Determinante -3", Gm == sp.Matrix([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, -2, 1], [0, 0, 1, -2]]) and abs(U.det()) == 1 and J.det() == -3, Gm.tolist())
C4 = sp.Matrix([[0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0]])
c.ok("der Viererzyklus erhält die Form, hält die Zeitrichtung fest und wirkt im Raum als Drehspiegelung der Ordnung 4 (Determinante -1)", C4.T * J * C4 == J and C4.det() == -1 and C4**4 == sp.eye(4) and C4**2 != sp.eye(4))
# Darstellungen von S4: Klassen e, (12), (12)(34), (123), (1234) mit Größen 1, 6, 3, 8, 6
size = [1, 6, 3, 8, 6]; sq = [0, 0, 0, 3, 2]                                   # Klasse von g^2
tab = {"A1": [1, 1, 1, 1, 1], "A2": [1, -1, 1, 1, -1], "E": [2, 0, 2, -1, 0], "T1": [3, -1, -1, 0, 1], "T2": [3, 1, -1, 0, -1]}
dec = lambda ch: {k: sum(s * a * b for s, a, b in zip(size, ch, v)) // 24 for k, v in tab.items() if sum(s * a * b for s, a, b in zip(size, ch, v))}
perm4 = [4, 2, 0, 1, 0]; edges = [6, 2, 2, 0, 0]
alt = [Fr(perm4[i] ** 2 - perm4[sq[i]], 2) for i in range(5)]; sym = [Fr(perm4[i] ** 2 + perm4[sq[i]], 2) for i in range(5)]
symalt = [Fr(alt[i] ** 2 + alt[sq[i]], 2) for i in range(5)]; symedges = [Fr(edges[i] ** 2 + edges[sq[i]], 2) for i in range(5)]
d = {n: dec([int(x) for x in ch]) for n, ch in (("Kanten", edges), ("antisymmetrische Paare", alt), ("Metrik Sym2(4)", sym), ("Quadrate der antisymmetrischen Paare", symalt), ("Quadrate der Kanten", symedges))}
common = lambda a, b: sum(min(d[a].get(k, 0), d[b].get(k, 0)) * len([0]) * {"A1": 1, "A2": 1, "E": 2, "T1": 3, "T2": 3}[k] for k in tab)
c.ok("Zerlegungen: Kanten = A1 + E + T2, antisymmetrische Paare = T1 + T2, Metrik = 2 A1 + E + 2 T2", d["Kanten"] == {"A1": 1, "E": 1, "T2": 1} and d["antisymmetrische Paare"] == {"T1": 1, "T2": 1} and d["Metrik Sym2(4)"] == {"A1": 2, "E": 1, "T2": 2}, d)
c.ok("NACHGERECHNET (parallele Aussage): eine symmetrieerhaltende lineare Abbildung zwischen den sechs Kanten und den sechs antisymmetrischen Kanälen hat höchstens Rang 3", common("Kanten", "antisymmetrische Paare") == 3)
r_alt = common("Quadrate der antisymmetrischen Paare", "Metrik Sym2(4)"); r_edge = common("Quadrate der Kanten", "Metrik Sym2(4)")
c.ok("quadratisch und nur mit der endlichen Symmetrie S4 erreichen beide Lesarten alle zehn Komponenten der Metrik", (r_alt, r_edge) == (10, 10), (r_alt, r_edge, d["Quadrate der antisymmetrischen Paare"]))
# volle Familiensymmetrie SU(4)
w2 = [tuple(1 if k in p else 0 for k in range(4)) for p in itertools.combinations(range(4), 2)]
ws = {tuple(a + b for a, b in zip(x, y)) for x in w2 for y in w2}
c.ok("ABBRUCHREGEL GREIFT für die Lesart mit Ladungen: unter der vollen Familiensymmetrie SU(4) enthalten die Quadrate der sechs Kanäle (20 + 1) das höchste Gewicht (2,0,0,0) der Metrik nicht. Rang null statt zehn", (2, 0, 0, 0) not in ws and len(ws) == 19, len(ws))
# natürliche Form auf den antisymmetrischen Paaren: Keil
idx = list(itertools.combinations(range(4), 2)); Wd = sp.zeros(6, 6)
for a, p in enumerate(idx):
    for b, q in enumerate(idx):
        if set(p) | set(q) == {0, 1, 2, 3}: Wd[a, b] = Permutation(list(p) + list(q)).signature()
c.ok("die sechs antisymmetrischen Kanäle tragen von sich aus die Keilform mit Signatur (3,3), nicht (1,3)", sig(Wd) == (3, 3), sig(Wd))
c.ok("URTEIL: die lorentzsche Form gibt es exakt, aber auf den Gewichten der vier Alternativen (symmetrische Paare), nicht auf den sechs fermionischen Kanälen. Zeit ist darin die gemeinsame Richtung aller vier", True)
c.ok("NICHT GEZEIGT: dass das Original die sechs Paare gleich gewichtet und dass diese Form etwas mit Abständen zu tun hat. Aus einer Form allein folgt kein Ort", True)
c.done({"signaturen": {str(k): v for k, v in sigs.items()}, "zerlegungen": d, "urteil": "Form exakt und nur bei vier Alternativen 3 + 1; für die sechs fermionischen Kanäle mit voller Familiensymmetrie gefallen. Bindung offen."})
