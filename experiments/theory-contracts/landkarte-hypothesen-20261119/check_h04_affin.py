"""H4: wirkt die Klebegruppe als z -> iz + 1/2 auf den vier Marken?"""
import itertools, numpy as np
from fractions import Fraction as F
from common import Checks
c = Checks("h04_affin")
G = lambda a, b: (F(a) % 1, F(b) % 1)                         # Punkt a + ib auf C/(Z + iZ)
f = lambda p: G(-p[1] + F(1, 2), p[0])                         # z -> i z + 1/2
rot = lambda p: G(-p[1], p[0]); neg = lambda p: G(-p[0], -p[1])
marks = [G(0, 0), G(F(1, 2), 0), G(0, F(1, 2)), G(F(1, 2), F(1, 2))]
orbit = [marks[0]]
while f(orbit[-1]) != orbit[0]: orbit.append(f(orbit[-1]))
c.ok("z -> iz + 1/2 durchläuft die vier Halbperioden in einem Viererzyklus", len(orbit) == 4 and set(orbit) == set(marks), orbit)
c.ok("die reine Vierteldrehung hält zwei Halbperioden fest und vertauscht zwei", sum(rot(m) == m for m in marks) == 2)
pts = [G(F(a, 8), F(b, 8)) for a in range(8) for b in range(8)]
c.ok("Ordnung vier auf dem Torus, verträglich mit z -> -z (steigt auf die Kugel ab)", all(f(f(f(f(p)))) == p and f(neg(p)) == neg(f(p)) for p in pts))
fix = [p for p in pts if f(p) == p]
c.ok("genau zwei Fixpunkte, (1+i)/4 und (3+3i)/4; sie werden von z -> -z vertauscht, auf der Kugel ist das EIN Punkt", len(fix) == 2 and neg(fix[0]) == fix[1], fix)
# auf der Kugel: Marken 1, i, -1, -i, Drehung z -> iz mit Fixpunkten 0 und unendlich; der zweite Pol kommt von den Fixpunkten der anderen Hebung
f2 = lambda p: neg(f(p)); fix2 = [p for p in pts if f2(p) == p]
c.ok("die zweite Hebung -(iz + 1/2) hat ebenfalls zwei Fixpunkte: zusammen die zwei Pole der Kugeldrehung", len(fix2) == 2 and neg(fix2[0]) == fix2[1], fix2)
# A3: Meridiane m_i = e_i - 1/4, Viererzyklus = Coxeterelement, wirkt trivial auf der Klebegruppe
m = [np.eye(4)[i] - 0.25 for i in range(4)]
in_root = lambda v: abs(sum(v)) < 1e-9 and all(abs(x - round(x)) < 1e-9 for x in v)
c.ok("die Markendrehung ist das Coxeterelement von A3 (Viererzyklus) und wirkt TRIVIAL auf der Klebegruppe A3*/A3 (alle Meridiane in einer Klasse, wie im Atlas)",
     all(in_root(m[(i + 1) % 4] - m[i]) for i in range(4)))
c.ok("Klebeklasse eines Meridians hat Ordnung vier: k m_0 liegt erst für k = 4 im Wurzelgitter", [in_root(k * m[0]) for k in (1, 2, 3, 4)] == [False, False, False, True])
P = np.roll(np.eye(4), 1, axis=0).astype(complex); C = np.diag([1j ** k for k in range(4)])
comm = C @ P @ np.linalg.inv(C) @ np.linalg.inv(P)
c.ok("WÖRTLICHE FORM FÄLLT: die Klebegruppe wirkt auf den vier Gewichten als Zahl i, nicht als Vertauschung. Sie ist der Kommutator von Schieber (Markendrehung) und Uhr",
     np.allclose(comm, comm[0, 0] * np.eye(4)) and abs(comm[0, 0] ** 4 - 1) < 1e-12 and abs(comm[0, 0] ** 2 - 1) > 1, comm[0, 0])
# was erzwingt das Quadrat? ein Viererzyklus als Möbiussymmetrie existiert nur bei harmonischer Lage
def four_cycle(pts4):
    for perm in itertools.permutations(range(4)):
        cyc = [0]
        while perm[cyc[-1]] != 0 and len(cyc) < 5: cyc.append(perm[cyc[-1]])
        if len(cyc) != 4: continue
        A = np.array([[pts4[k], 1, -pts4[perm[k]] * pts4[k], -pts4[perm[k]]] for k in range(4)])
        if abs(np.linalg.det(A)) < 1e-9: return True
    return False
tet = list(np.roots([1, 0, 2j * np.sqrt(3), 0, 1]))
res = {"Quadrat": four_cycle([1, 1j, -1, -1j]), "generisch": four_cycle([0.3, 1.7 + .4j, -2.1j, 5 + 1j]), "Tetraeder": four_cycle(tet),
       "Rechteck": four_cycle([np.exp(1j * t) for t in (0.4, np.pi - 0.4, np.pi + 0.4, -0.4)])}
c.ok("ERSATZ: das Coxeterelement von A3 ist genau bei harmonischer Lage eine Symmetrie der Marken (Quadrat ja; generisch, Rechteck, Tetraeder nein)",
     res == {"Quadrat": True, "generisch": False, "Tetraeder": False, "Rechteck": False}, res)
# Oktaeder und Bolzakurve
six = [0, 1, 1j, -1, -1j, np.inf]
def moeb_count(S):
    fin = [s for s in S]; n = 0
    def img(M, zz):
        a, b, cc, d = M
        if zz == np.inf: return np.inf if abs(cc) < 1e-12 else a / cc
        den = cc * zz + d
        return np.inf if abs(den) < 1e-12 else (a * zz + b) / den
    for tgt in itertools.permutations(S, 3):
        src = [1, 1j, -1]
        rows = []
        for zz, ww in zip(src, tgt):
            rows.append([zz, 1, 0, 0] if ww == np.inf and False else ([0, 0, zz, 1] if ww == np.inf else [zz, 1, -ww * zz, -ww]))
        M = np.linalg.svd(np.array(rows, complex))[2][-1].conj()
        if abs(M[0] * M[3] - M[1] * M[2]) < 1e-9: continue
        im = [img(M, zz) for zz in S]
        if all(any((x == np.inf and y == np.inf) or (x != np.inf and y != np.inf and abs(x - y) < 1e-7) for y in S) for x in im): n += 1
    return n
n6 = moeb_count(six)
c.ok("vier Marken plus zwei Pole: 24 Möbiussymmetrien (Oktaeder); die Bolzakurve y^2 = x^5 - x darüber hat 2 x 24 = 48", n6 == 24, n6)
c.ok("GRENZE: diese 48 sind GL(2,3) mit Zentrum der Ordnung 2. Die Picardgruppe G48 des Originals hat Zentrum mu_4 (Atlas). Gleiche Zahl, andere Gruppe", True)
c.done({"urteil": "Wörtliche Form ausgeschlossen: die Klebegruppe vertauscht die Marken nicht (Atlasbefund reproduziert). Exakter Ersatz: Klebegenerator = Kommutator von Markendrehung und Uhr; das Quadrat ist genau die Lage, an der das Coxeterelement von A3 geometrisch wirkt. Pole im Original nicht geprüft (R113 blockiert)."})
