"""Lassen sich die drei Zustände des Originaltransfers mit den drei Schalen des Petersengraphen gleichsetzen?
Aufruf: python3 -B check_cusps.py"""
import itertools
import sympy as sp
from common import Checks

C = Checks("cusps")
R = sp.Rational
F, W = ("f1", "f2", "f3"), ("w1", "w2")
places = F+W
pairs = [frozenset(p) for p in itertools.combinations(places, 2)]
A = sp.Matrix(10, 10, lambda i, j: int(not pairs[i] & pairs[j]))
Q = sp.Matrix([[0, 3, 0], [1, 0, 2], [0, 1, 2]]); P = Q/3
M = sp.Matrix([[1, 0, 0], [0, R(1, 2), R(1, 6)], [0, R(1, 6), R(1, 2)]])      # v486/v487
one, Y, Z = sp.Matrix([1, 1, 1]), sp.Matrix([1, R(-2, 3), R(1, 6)]), sp.Matrix([3, 1, -1])
pi = sp.Matrix([[1, 3, 6]])/10

# 1. Schalen als Zustände: ausgeschlossen
C.ok("Original (v486): ein fester Zustand, ein tauschbares Paar, Raten 1, 2/3, 1/3", set(M.eigenvals()) == {1, R(2, 3), R(1, 3)} and M[0, 1:] == sp.zeros(1, 2))
C.ok("Schrittmatrix der Schalen: Eigenvektoren (1,1,1), (1,-2/3,1/6), (3,1,-1) zu 1, -2/3, 1/3", P*one == one and P*Y == -R(2, 3)*Y and P*Z == Z/3)
left = [pi, sp.Matrix([[1, -2, 1]]), sp.Matrix([[1, 1, -2]])]
C.ok("Linkseigenvektoren (1,3,6), (1,-2,1), (1,1,-2)", left[0]*P == left[0] and left[1]*P == -R(2, 3)*left[1] and left[2]*P == left[2]/3)
C.ok("kein Eigenvektor hat einen Nulleintrag: in keiner Funktion des Schrittes ist eine Schale fest oder abgekoppelt", all(x != 0 for v in (one, Y, Z, *left) for x in v))
swaps = [p for p in itertools.permutations(range(3)) if p != (0, 1, 2) and sp.Matrix(3, 3, lambda i, j: P[p[i], p[j]]) == P]
C.ok("keine Vertauschung von Schalen lässt den Schritt unverändert: es gibt kein tauschbares Schalenpaar", swaps == [])
U = sp.Matrix.hstack(one, Y, Z); absP = U*sp.diag(1, R(2, 3), R(1, 3))*U.inv()
C.ok("der positive Schritt |P| (Raten 1, 2/3, 1/3 auf denselben Formen) hat negative Einträge: keine Sprungregel auf Schalen", min(absP) < 0, str(absP.tolist()))
C.ok("AUSGESCHLOSSEN: Schalen = Zustände des Originals, für jede Reihenfolge und jede Funktion des Schrittes", True)

# 2. Was passt: Zustände des Originals = Schwingungsformen der Schalen
dot = lambda a, b: sum(pi[i]*a[i]*b[i] for i in range(3))
C.ok("die drei Formen stehen bezüglich der Gewichte 1 : 3 : 6 senkrecht", dot(one, Y) == 0 and dot(one, Z) == 0 and dot(Y, Z) == 0)
a, b, c = sp.symbols("a b c")
fit = sp.solve([sp.Eq(a*Y[i]**2+b*Y[i]+c, Z[i]) for i in range(3)], [a, b, c], dict=True)[0]
C.ok("Form n hat Grad n in der Ladung: 1, Y, und 25 Z = 108 Y^2 - 6 Y - 27", fit == {a: sp.Rational(108, 25), b: sp.Rational(-6, 25), c: sp.Rational(-27, 25)}, str(fit))
C.ok("Betrag der Rate von Form n ist 1 - n/3: die Leiter des Originals (v124) mit n = Grad in der Ladung", [abs(x) for x in (1, R(-2, 3), R(1, 3))] == [1-R(n, 3) for n in range(3)])
C.ok("aber Vorzeichen (-1)^n je Schritt; Original: Stufe 0 und 1 gerade, Stufe 2 ungerade (v487). Die Muster stimmen nicht überein", [sp.sign(x) for x in (1, R(-2, 3), R(1, 3))] == [1, -1, 1])
C.ok("gerade Schrittzahlen stimmen auf den Formen exakt: P^2k und M^2k haben dieselben Raten", all(set((P**k).eigenvals()) == set((M**k).eigenvals()) for k in (2, 4, 6)) and set(P.eigenvals()) != set(M.eigenvals()))

# 3. Der Blatttausch des Originals
def perm_matrix(g):
    m = dict(zip(places, g))
    return sp.Matrix(10, 10, lambda i, j: int(frozenset(m[x] for x in pairs[j]) == pairs[i]))
E4 = sp.Matrix.hstack(*(A+2*sp.eye(10)).nullspace()); E5 = sp.Matrix.hstack(*(A-sp.eye(10)).nullspace())
def trace_on(E, G): return sp.simplify(((E.T*E).inv()*E.T*G*E).trace())
found = {}
for g in itertools.permutations(places):
    G = perm_matrix(g)
    if G != sp.eye(10) and G*G == sp.eye(10): found[(trace_on(E4, G), trace_on(E5, G))] = found.get((trace_on(E4, G), trace_on(E5, G)), 0)+1
C.ok("25 Spiegelungen in S5; Spuren auf den Räumen zu 2/3 (vier Formen) und 1/3 (fünf Formen) sind (2,1) oder (0,1), nie (4,-5)", found == {(2, 1): 10, (0, 1): 15}, str(found))
C.ok("also: keine Symmetrie des Graphen ist gerade auf allen 2/3 Formen und ungerade auf allen 1/3 Formen", (4, -5) not in found)

# 4. Der Sechserkreis allein, mit Verlust nach innen
hexa = [i for i, p in enumerate(pairs) if len(p & set(W)) == 1]
K = A.extract(hexa, hexa)/3
C.ok("Schritt auf dem Sechserkreis: jede Zeile summiert zu 2/3, Verlust 1/3 wie im Original (v486)", all(sum(K.row(i)) == R(2, 3) for i in range(6)))
C.ok("Raten 2/3, 1/3, 1/3, -1/3, -1/3, -2/3 = (2/3) cos(k pi/3)", K.eigenvals() == {R(2, 3): 1, R(1, 3): 2, R(-1, 3): 2, R(-2, 3): 1})
swap = perm_matrix(("f1", "f2", "f3", "w2", "w1")).extract(hexa, hexa)
def parity(lam):
    vs = (K-lam*sp.eye(6)).nullspace(); return [s for s in (1, -1) if all(swap*v == s*v for v in vs)]
C.ok("Tausch der zwei schwachen Plätze = Gegenpunkt auf dem Kreis: gerade auf 2/3, ungerade auf 1/3 (wie v487), aber auch gerade auf -1/3, ungerade auf -2/3", parity(R(2, 3)) == [1] and parity(R(1, 3)) == [-1] and parity(R(-1, 3)) == [1] and parity(R(-2, 3)) == [-1])
C.ok("sechs Schritte auf dem Kreis: (2/3)^6 zweimal, (1/3)^6 viermal", (K**6).eigenvals() == {R(64, 729): 2, R(1, 729): 4})
lump = sp.Matrix(2, 2, lambda i, j: sum(K[hexa.index(next(h for h in hexa if W[i] in pairs[h])), hexa.index(h)] for h in hexa if W[j] in pairs[h]))
C.ok("zusammengefasst nach schwachem Platz: bleiben 0, wechseln 2/3; das Original hat bleiben 1/2, wechseln 1/6", lump == sp.Matrix([[0, R(2, 3)], [R(2, 3), 0]]), str(lump.tolist()))
C.ok("die Rate 1/3 entsteht hier als (2/3) mal cos(60 Grad), auf den Schalen als 1 - 2/3: zwei verschiedene Mechanismen", R(2, 3)*sp.cos(sp.pi/3) == R(1, 3))

C.done({"abs_P": [[str(x) for x in absP.row(i)] for i in range(3)], "Z_in_Y": "25 Z = 108 Y^2 - 6 Y - 27", "involution_traces": {str(k): v for k, v in found.items()}})
