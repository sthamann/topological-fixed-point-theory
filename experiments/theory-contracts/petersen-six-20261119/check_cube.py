"""Lesart: drei Familienbits, sechs gleich wahrscheinliche Züge (Bit i auf Wert v setzen). Was folgt daraus exakt,
und wie eindeutig ist das? Aufruf: python3 -B check_cube.py"""
import itertools
from fractions import Fraction as Fr
import sympy as sp
from common import Checks

C = Checks("cube")
R = sp.Rational

def walk(d, q):
    """Wähle eine von d Stellen und einen von q Werten gleich wahrscheinlich, setze die Stelle auf den Wert."""
    states = list(itertools.product(range(q), repeat=d)); idx = {s: i for i, s in enumerate(states)}
    L = sp.zeros(len(states), len(states))
    for s in states:
        for i in range(d):
            for v in range(q): L[idx[s], idx[s[:i]+(v,)+s[i+1:]]] += R(1, d*q)
    return states, idx, L

states, idx, L = walk(3, 2)

# 1. Die Regel und was das Original getrennt festhält
C.ok("sechs Züge = 3 Bits mal 2 Werte; je Zug Wahrscheinlichkeit 1/6", 3*2 == 6)
C.ok("von jedem Zustand: bleiben 1/2, zu jedem der drei Nachbarn 1/6 (v486: s = 1/2, h = 1/6, dort aus dem Spektrum gewählt)", all(L[i, i] == R(1, 2) for i in range(8)) and {x for x in L if x not in (0, R(1, 2))} == {R(1, 6)})
C.ok("Raten 1, 2/3, 1/3, 0 mit Vielfachheit 1, 3, 3, 1: die Leiter 1 - n/3 (v124)", L.eigenvals() == {1: 1, R(2, 3): 3, R(1, 3): 3, 0: 1})
C.ok("die Leiter endet bei n = 3 mit Rate genau 0: die Wand bei der Familienzahl (v124 Punkt 4)", L.eigenvals()[0] == 1)
a, b = idx[(0, 0, 0)], idx[(1, 0, 0)]
pair = L.extract([a, b], [a, b])
C.ok("ein Paar entlang einer Kante: bleiben 1/2, wechseln 1/6, Verlust 1/3 an die zwei anderen Nachbarn (v486)", pair == sp.Matrix([[R(1, 2), R(1, 6)], [R(1, 6), R(1, 2)]]) and 1-sum(pair.row(0)) == R(1, 3))
C.ok("gerade Paarform Rate 2/3, ungerade 1/3 (v487)", pair*sp.Matrix([1, 1]) == R(2, 3)*sp.Matrix([1, 1]) and pair*sp.Matrix([1, -1]) == R(1, 3)*sp.Matrix([1, -1]))
n = sp.symbols("n", positive=True)
C.ok("Verlust (N-1)/(2N) gleich dem 1/N des Originals nur bei N = 3", sp.solve(sp.Eq((n-1)/(2*n), 1/n), n) == [3])
flip = sp.Matrix(8, 8, lambda i, j: int(all(x != y for x, y in zip(states[i], states[j]))))
walsh = lambda S: sp.Matrix([(-1)**sum(s[i] for i in S) for s in states])
C.ok("Umklappen aller Bits vertauscht mit dem Schritt und hat auf Stufe n das Vorzeichen (-1)^n: jede Stufe wechselt die Parität (v487)", flip*L == L*flip and all(L*walsh(S) == (1-R(len(S), 3))*walsh(S) and flip*walsh(S) == (-1)**len(S)*walsh(S) for k in range(4) for S in itertools.combinations(range(3), k)))
C.ok("der feste Zustand ist die Gleichverteilung (v327: der Attraktor ist das demokratische Gesetz)", sp.ones(1, 8)*L == sp.ones(1, 8))

# 2. Der symmetrische Teil: Urne mit drei Kugeln
weight = lambda k: [i for i, s in enumerate(states) if sum(s) == k]
urn = sp.Matrix(4, 4, lambda j, k: sum(L[weight(j)[0], t] for t in weight(k)))
C.ok("zusammengefasst nach Zahl gesetzter Bits: vier Zustände = Familienzahl + 1 (im Original: vier Marken)", urn.shape == (4, 4) and urn == sp.Matrix([[R(1, 2), R(1, 2), 0, 0], [R(1, 6), R(1, 2), R(1, 3), 0], [0, R(1, 3), R(1, 2), R(1, 6)], [0, 0, R(1, 2), R(1, 2)]]))
C.ok("dort jede Rate genau einmal: 1, 2/3, 1/3 und die Wand 0", urn.eigenvals() == {1: 1, R(2, 3): 1, R(1, 3): 1, 0: 1})

# 3. Die Sechs
C.ok("eine Hand = so viele Züge wie es Züge gibt (sechs): Raten 1, (2/3)^6, (1/3)^6, 0", (L**6).eigenvals() == {1: 1, R(64, 729): 3, R(1, 729): 3, 0: 1})
C.ok("L = 1 + G/6 mit G = -2n: die lineare Uhr hat Steigung 2 je Stufe (v124 Punkt 3)", all(6*((1-R(k, 3))-1) == -2*k for k in range(4)))
cyc = lambda m: ((m[0]+1) % 3, 1-m[1])
orbit = [(0, 0)]
while cyc(orbit[-1]) != orbit[0]: orbit.append(cyc(orbit[-1]))
C.ok("Bits im Kreis weiterrücken und Wert umklappen hat Ordnung 6 und durchläuft alle sechs Züge (v118: -M0 hat Ordnung 6)", len(orbit) == 6 and len(set(orbit)) == 6)

# 4. Dasselbe auf Operatoren: sechs Majoranagrößen dreier Moden
X, Yp, Zp, I2 = sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -sp.I], [sp.I, 0]]), sp.Matrix([[1, 0], [0, -1]]), sp.eye(2)
kron = lambda *m: sp.kronecker_product(*m)
gam = [kron(X, I2, I2), kron(Yp, I2, I2), kron(Zp, X, I2), kron(Zp, Yp, I2), kron(Zp, Zp, X), kron(Zp, Zp, Yp)]
C.ok("sechs Majoranagrößen: Quadrat 1, paarweise antivertauschend", all((gam[i]*gam[j]+gam[j]*gam[i]) == 2*int(i == j)*sp.eye(8) for i in range(6) for j in range(6)))
def channel(A): return sum((A/2+g*A*g/2 for g in gam), sp.zeros(8, 8))/6
spec = {}
for k in (0, 2, 4, 6):
    for S in itertools.combinations(range(6), k):
        A = sp.eye(8)
        for i in S: A = A*gam[i]
        img = channel(A); lam = 1-R(k, 6)
        assert sp.simplify(img-lam*A) == sp.zeros(8, 8)
        spec[lam] = spec.get(lam, 0)+1
C.ok("Kanal 'eine der sechs mit Wahrscheinlichkeit 1/2 anwenden': auf geraden Produkten Rate 1 - |S|/6, Vielfachheit 1, 15, 15, 1", spec == {1: 1, R(2, 3): 15, R(1, 3): 15, 0: 1})
C.ok("15 = Dimension von A3 (Drehungen in sechs Richtungen), 6 = Zahl seiner positiven Wurzeln (v124: p_2 = 6 = |R+(A3)|)", sp.binomial(6, 2) == 15 and 4*3//2 == 6)

# 5. Wie eindeutig? Alle Regeln 'd Stellen, q Werte' gegen fünf Angaben des Originals
target = {"stay": Fr(1, 2), "hop": Fr(1, 6), "leak": Fr(1, 3), "second": Fr(2, 3), "hand": 6}
table = {}
for d in range(1, 7):
    for q in range(2, 6):
        got = {"stay": Fr(1, q), "hop": Fr(1, d*q), "leak": Fr(q-1, q)-Fr(1, d*q), "second": 1-Fr(1, d), "hand": d*q}
        table[(d, q)] = sum(got[k] == target[k] for k in target)
C.ok("unter 24 Regeln mit 1 bis 6 Stellen und 2 bis 5 Werten trifft nur (3, 2) alle fünf Angaben; die zweitbeste trifft zwei", [k for k, v in table.items() if v == 5] == [(3, 2)] and sorted(table.values())[-2] <= 2, str(sorted(table.items(), key=lambda kv: -kv[1])[:4]))
C.ok("die Raten 1 - n/3 legen nur die drei Stellen fest; die zwei Werte kommen aus bleiben 1/2 und aus der Hand sechs", {d for (d, q) in table if 1-Fr(1, d) == Fr(2, 3)} == {3} and Fr(1, 2) == target["stay"])
for d, q in ((2, 2), (4, 2), (3, 3), (2, 3)):
    _, _, M = walk(d, q)
    C.ok(f"Gegenprobe {d} Stellen, {q} Werte mit eigener Hand {d*q}: nicht beide Werte 64/729 und 1/729", not {R(64, 729), R(1, 729)} <= set((M**(d*q)).eigenvals()))

C.done({"regeltabelle": {f"{d},{q}": v for (d, q), v in table.items()}})
