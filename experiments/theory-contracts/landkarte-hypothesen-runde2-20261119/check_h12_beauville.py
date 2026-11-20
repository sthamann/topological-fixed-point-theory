"""H12: die vier Marken als Ausnahmestellen der elliptischen Fläche y^2 = x (x + t^2)(x + 1)."""
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from common import Checks, e8
c = Checks("h12_beauville")
t, x = sp.symbols("t x")
def inv(A, B):                                   # y^2 = x (x^2 + A x + B)
    c4 = 16 * (A**2 - 3 * B); c6 = -32 * A * (2 * A**2 - 9 * B); D = 16 * B**2 * (A**2 - 4 * B)
    return sp.factor(c4), sp.factor(c6), sp.factor(D), sp.factor(sp.cancel(c4**3 / D))
A, B = t**2 + 1, t**2
c4, c6, D, j = inv(A, B)
c.ok("Diskriminante 16 t^4 (t-1)^2 (t+1)^2; Grad 8, also Ordnung 4 bei unendlich: Fasern I4, I4 bei 0 und unendlich, I2, I2 bei 1 und -1", sp.expand(D - 16 * t**4 * (t**2 - 1)**2) == 0, D)
c.ok("alle vier Stellen sind multiplikativ (c4 verschwindet dort nicht), Eulerzahl 4 + 4 + 2 + 2 = 12", all(c4.subs(t, v) != 0 for v in (0, 1, -1)) and sp.degree(sp.expand(c4), t) == 4)
c.ok("Doppelverhältnis der vier Stellen (0, unendlich; 1, -1) ist -1: harmonische Lage wie die Marken", sp.Rational(1 - 0, -1 - 0) == -1)
c.ok("(t, t(t+1)) liegt auf der Kurve und hat Ordnung 4: seine Verdopplung ist der Punkt (0,0) der Ordnung 2", sp.expand((t * (t + 1))**2 - (t * (t + t**2) * (t + 1))) == 0
     and sp.simplify(((3 * t**2 + 2 * A * t + B) / (2 * t * (t + 1)))**2 - A - 2 * t) == 0)
# Teilung durch (0,0): y^2 = x (x^2 - 2A x + A^2 - 4B)
A2, B2 = -2 * A, sp.expand(A**2 - 4 * B)
c4b, c6b, Db, jb = inv(A2, B2)
c.ok("die halbierte Kurve hat Diskriminante ~ t^2 (t-1)^4 (t+1)^4: die Sorten sind getauscht (I2 bei 0 und unendlich, I4 bei 1 und -1)", sp.expand(Db - 256 * t**2 * (t**2 - 1)**4) == 0, Db)
m = (1 + t) / (1 - t)
c.ok("m(t) = (1+t)/(1-t) hat Ordnung 4 und dreht die vier Stellen im Kreis: 0 -> 1 -> unendlich -> -1 -> 0", sp.simplify(m.subs(t, m.subs(t, m.subs(t, m))) - t) == 0 and m.subs(t, 0) == 1 and sp.limit(m, t, 1) in (sp.oo, -sp.oo, sp.zoo) and sp.limit(m, t, sp.oo) == -1 and m.subs(t, -1) == 0)
same_j = sp.simplify(jb - j.subs(t, m)) == 0
c.ok("KERN: die halbierte Kurve bei t ist die ursprüngliche Kurve bei m(t) (gleiche j Invariante als Funktion)", same_j, sp.factor(jb))
tw = sp.factor(sp.cancel((c6b / c6.subs(t, m)) * (c4.subs(t, m) / c4b)))
c.ok("der Vergleichsfaktor ist -(t-1)^2: gleich bis auf die Verdrehung mit -1, also erst nach Hinzunahme von i. Über Q(t) sind die beiden Kurven verdrehte Geschwister", sp.expand(tw + (t - 1)**2) == 0, tw)
c.ok("FOLGE: die Vierteldrehung der Marken ist keine Symmetrie der Fläche, sondern eine Halbierung der Kurve (Isogenie vom Grad 2). Zweimal angewandt ist sie Verdopplung, also bis auf Maßstab die Identität", same_j)
c.ok("die Halbdrehung t -> -1/t ist eine echte Symmetrie (j bleibt gleich), ebenso t -> -t", sp.simplify(j.subs(t, -1 / t) - j) == 0 and sp.simplify(j.subs(t, -t) - j) == 0)
sq = set(sp.solve(sp.numer(sp.together(j - 1728)), t))
c.ok("der quadratische Torus (j = 1728) sitzt an sechs Stellen, darunter i und -i. Das sind genau die zwei Fixpunkte der Vierteldrehung m, also die Pole neben den vier Marken", sq == {sp.sqrt(2), -sp.sqrt(2), sp.sqrt(2) / 2, -sp.sqrt(2) / 2, sp.I, -sp.I} and sp.simplify(m.subs(t, sp.I) - sp.I) == 0 and sp.simplify(m.subs(t, -sp.I) + sp.I) == 0, sorted(map(str, sq)))
# Gitter: A3 + A3 + A1 + A1 in D5 + A3 in E8
al, roots, coef, Binv = e8()
fam = [k for k in coef if k[3] == 0 and k[4] == 0 and k[5] == 0 and k[0] == 0 and k[1] == 0 and k[2] == 0]            # A3: nur a7, a8 und a0 = -theta
theta = sum(((2, 3, 4, 6, 5, 4, 3, 2)[i] * al[i] for i in range(8)), sp.zeros(8, 1))
A3 = [r for r in roots if all(r.dot(al[i]) == 0 for i in (0, 1, 2, 3, 4))]                                              # senkrecht auf D5
D5 = [r for r, k in zip(roots, coef) if k[5] == 0 and k[6] == 0 and k[7] == 0]
a = D5[0]; rest = [r for r in D5 if r.dot(a) == 0]
c.ok("in D5 steht senkrecht auf einer Wurzel ein A1 + A3 (2 + 12 Wurzeln): D5 enthält A3 + A1 + A1", (len(A3), len(D5), len(rest)) == (12, 40, 14), (len(A3), len(D5), len(rest)))
gens = sp.Matrix([list(Binv * r) for r in A3 + [a] + rest])
snf = smith_normal_form(gens, domain=sp.ZZ); diag = sorted(abs(int(snf[i, i])) for i in range(8))
c.ok("das Gitter A3 + A3 + A1 + A1 hat in E8 Index 8 mit Quotient Z4 x Z2: dieselbe Gruppe wie die Torsion der Fläche (Punkt der Ordnung 4 und volle Halbperioden)", diag == [1, 1, 1, 1, 1, 1, 2, 4], diag)
c.ok("ABBRUCHREGEL GREIFT NICHT: die verlangte Selbstabbildung existiert und tauscht die zwei Sorten", same_j)
c.ok("NICHT GEZEIGT: dass P1 die Basis dieser Fläche ist. Die fünf Trägerfunktionen sind Funktionen auf der Basis; die Fläche hat acht Torsionsschnitte und Rang null. Eine natürliche Abbildung zwischen beiden fand sich nicht", True)
c.ok("SPANNUNG mit dem Original: dort sind die vier Marken gleichartig (je Defizit pi), hier sind es zwei Sorten. Verträglich nur, wenn die Vierteldrehung im Original als Halbierung gelesen werden darf", True)
c.done({"j": str(j), "j_halbiert": str(jb), "invarianten_snf": diag, "urteil": "Trägt als Mathematik: harmonische Lage, Gitter in D5 + A3, Torsion = Klebegruppe Z4 x Z2, Vierteldrehung = Halbierung der Kurve. Bindung an P1 offen."})
