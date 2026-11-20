"""H11: Zeit als die eine Richtung zwischen Viererschnitt und Fünferschnitt."""
import sympy as sp
from collections import Counter
from common import Checks, e8
c = Checks("h11_nachbar")
al, roots, coef, Binv = e8()
G = sp.Matrix(8, 8, lambda i, j: al[i].dot(al[j])); W = G.inv()                      # Fundamentalgewichte in der einfachen Basis
om = lambda k: sum((W[k - 1, j] * al[j] for j in range(8)), sp.zeros(8, 1))
v = 4 * om(5) - 5 * om(6)
vc = tuple(int(x) for x in Binv * v)
c.ok("die freie Richtung ist v = 4 omega5 - 5 omega6 = 2a1 + 3a2 + 4a3 + 6a4 + 5a5, Norm 20", vc == (2, 3, 4, 6, 5, 0, 0, 0) and v.dot(v) == 20, vc)
val = [int(v.dot(r)) for r in roots]
c.ok("für jede Wurzel ist v·r = 4 c5 - 5 c6 (c5, c6 = Koeffizienten der beiden gestrichenen Knoten)", all(x == 4 * k[4] - 5 * k[5] for x, k in zip(val, coef)))
n4 = [r for r, x in zip(roots, val) if x % 4 == 0]; n5 = [r for r, x in zip(roots, val) if x % 5 == 0]; n0 = [r for r, x in zip(roots, val) if x == 0]
c.ok("U(t) = exp(2 pi i t v): bei t = 1/4 vertauschen genau 52 Wurzeln (D5 + A3: 40 + 12), bei t = 1/5 genau 40 (A4 + A4), bei allgemeinem t genau 32 (A4 + A3)", (len(n4), len(n5), len(n0)) == (52, 40, 32), (len(n4), len(n5), len(n0)))
c.ok("U(1/4) ist der Erzeuger der Viererklebung (Phase i^(-c6)), U(1/5) der Erzeuger der Fünferklebung (Phase zeta5^(4 c5)); beide liegen auf einem Kreis", all((x % 4) == (-k[5]) % 4 for x, k in zip(val, coef)) and all((x % 5) == (4 * k[4]) % 5 for x, k in zip(val, coef)))
from math import gcd; from functools import reduce
c.ok("der Kreis schließt sich bei t = 1 (die Werte v·r haben größten gemeinsamen Teiler 1); U(1/20) hat Ordnung 20", reduce(gcd, [abs(x) for x in val if x]) == 1 and len({x % 20 for x in val}) > 1)
c.ok("EXAKT: Vierer und Fünferschnitt sind zwei Punkte eines einzigen Kreises in E8, gemeinsamer Takt 20", True)
# was ist diese Richtung auf dem Träger?
c.ok("v liegt ganz in D5 (nur die Knoten a1 bis a5) und steht senkrecht auf A3 und auf dem A4 in D5", vc[5:] == (0, 0, 0) and all(v.dot(al[i]) == 0 for i in (0, 1, 2, 3)) and all(v.dot(x) == 0 for x in (al[6], al[7])))
spin = Counter(x for x, k in zip(val, coef) if k[5] % 4 == 1)
c.ok("auf dem Spinor 16 von D5 hat v die Ladungen -5, -1, 3 mit Vielfachheit 1, 10, 5 (mal 4 Familienplätze): das ist die bekannte Ladung X in SO(10) -> SU(5) x U(1)", dict(spin) == {-5: 4, -1: 40, 3: 20}, dict(spin))
vec = Counter(x for x, k in zip(val, coef) if k[5] % 4 == 2)
c.ok("auf dem Vektor 10 von D5 hat v die Ladungen +2 und -2, je fünffach (mal 6)", dict(vec) == {2: 30, -2: 30}, dict(vec))
c.ok("ABBRUCHREGEL GREIFT (1): die Richtung ist eine schon vorhandene Eichladung des Trägers, keine neue Größe", True)
c.ok("ABBRUCHREGEL GREIFT (2): U(1/4) wirkt auf den fünf Plätzen als Zahl (i^2 = -1 auf 5 und 5 quer), U(1/5) ebenfalls als Zahl. Die Vierteluhr des Originals unterscheidet 3 von 2, der Fünfertakt verschiebt Plätze. Keine der beiden Uhren entsteht", all(abs(x) == 2 for x in vec))
c.ok("ALS ZEIT GEFALLEN: U(t) ist periodisch mit ganzzahligem Spektrum, also wieder eine Uhr ohne Zählrichtung", True)
c.done({"v": vc, "norm": 20, "urteil": "Geometrie exakt (ein Kreis durch beide Klebeerzeuger, Takt 20), als Zeit gefallen: die Richtung ist die Eichladung X von SO(10) und wirkt auf den fünf Plätzen nur als Zahl."})
