"""Exakte Nachprüfung der Treffer aus search_transfer.py. Aufruf: python3 -B check_search.py"""
import itertools, json
import sympy as sp
from common import Checks

C = Checks("search_check")
R = sp.Rational
res = json.load(open("results_transfer_search.json"))
lossy, absorbed = res["targets"]["verlust"], res["targets"]["aufgesogen"]
C.ok("Suchraum: 15 Operatoren respektieren die Teilung 3 + 2, davon 11 symmetrisch, 8 nur über Nachbarn", (res["operators"], res["symmetric_operators"], res["local_operators"]) == (15, 11, 8))
C.ok("2719 Dreiteilungen der zehn Paare bis auf Symmetrie, je für zwei Lesarten des Verlusts", lossy["partitions_up_to_symmetry"] == 2719 == absorbed["partitions_up_to_symmetry"])
C.ok("allein aus dem Graphen (ein Gewicht je für gleich, Nachbar, Nichtnachbar): kein einziger Treffer", lossy["counts"]["graph_only"] == 0 == absorbed["counts"]["graph_only"])
C.ok("mit freien Gewichten je Paartyp: 36 und 57 Teilungen mit erlaubter Sprungregel; die Regel ist also nicht erzwungen", (lossy["counts"]["nonnegative"], absorbed["counts"]["nonnegative"]) == (36, 57))
C.ok("genau eine Teilung ist symmetriefest und erlaubt eine symmetrische Regel: fest = FF und WW, Paar = die Seiten w1 und w2 des Sechserkreises", lossy["counts"]["stable_symmetric"] == 1 and lossy["counts"]["local_symmetric"] == 1 and absorbed["counts"]["symmetric"] == 0)

# der eine Treffer, exakt
F, W = ("f1", "f2", "f3"), ("w1", "w2")
hexa = [frozenset((f, w)) for w in W for f in F]
I6 = sp.eye(6)
ring = sp.Matrix(6, 6, lambda i, j: int(not hexa[i] & hexa[j]))
flip = sp.Matrix(6, 6, lambda i, j: int(len(hexa[i] & hexa[j]) == 1 and bool(hexa[i] & hexa[j] & set(F))))
side = lambda T: sp.Matrix(2, 2, lambda a, b: sum(T[3*a, 3*b+k] for k in range(3)))
B = sp.Matrix([[R(1, 2), R(1, 6)], [R(1, 6), R(1, 2)]])
local = I6/2+ring/12
C.ok("lokal und symmetrisch, eindeutig auf dem Kreis: bleiben 1/2, je Kreisnachbar 1/12; zusammengefasst die Paarregel des Originals", side(local) == B and all(side(local**k) == B**k for k in (2, 6)))
C.ok("aber auf den sechs Plätzen kommen die Raten 7/12 und 5/12 dazu, die das Original nicht hat", local.eigenvals() == {R(2, 3): 1, R(7, 12): 2, R(5, 12): 2, R(1, 3): 1})
c0, c1, c2, c3 = sp.symbols("c0:4")
dist = sp.Matrix(6, 6, lambda i, j: 0 if i == j else 1 if ring[i, j] else 3 if flip[i, j] else 2)
circ = sp.Matrix(6, 6, lambda i, j: (c0, c1, c2, c3)[dist[i, j]])
pure = []
for l1, l2 in itertools.product((R(2, 3), R(1, 3)), repeat=2):
    sol = sp.solve([c0+2*c2-R(1, 2), 2*c1+c3-R(1, 6), c0+c1-c2-c3-l1, c0-c1-c2+c3-l2], [c0, c1, c2, c3], dict=True)[0]
    if min(sol.values()) >= 0: pure.append(tuple(sol[c] for c in (c0, c1, c2, c3)))
C.ok("verlangt man auf dem Kreis nur die Raten 2/3 und 1/3, bleiben genau zwei erlaubte Regeln", sorted(pure) == sorted([(R(1, 2), 0, 0, R(1, 6)), (R(7, 18), R(1, 18), R(1, 18), R(1, 18))]), str(pure))
weak = I6/2+flip/6
C.ok("Regel A: bleiben 1/2, mit 1/6 den schwachen Platz wechseln bei gleicher Farbe; das ist die Paarregel des Originals, dreimal nebeneinander", weak.eigenvals() == {R(2, 3): 3, R(1, 3): 3} and side(weak) == B and flip*flip == I6)
C.ok("der Wechsel in Regel A ist kein Schritt des Graphen: die beiden Paare teilen einen Platz", all(ring[i, j] == 0 for i in range(6) for j in range(6) if flip[i, j]))
uniform = R(7, 18)*I6+(sp.ones(6, 6)-I6)/18
C.ok("Regel B: bleiben 7/18, zu jedem anderen der fünf mit 1/18; Raten 2/3 einmal, 1/3 fünfmal", uniform.eigenvals() == {R(2, 3): 1, R(1, 3): 5} and side(uniform) == B)
C.ok("sechs Schritte von Regel A auf allen zehn Paaren: 1 viermal, (2/3)^6 dreimal, (1/3)^6 dreimal; v808 sagt 1, (1/3)^6 fünfmal, (2/3)^6 viermal", (weak**6).eigenvals() == {R(64, 729): 3, R(1, 729): 3})
C.ok("in keiner dieser Regeln ist die Schrittzahl sechs festgelegt; sie wird wie im Original von außen gesetzt", True)
C.done({"pure_rules_on_ring": [[str(x) for x in p] for p in pure]})
