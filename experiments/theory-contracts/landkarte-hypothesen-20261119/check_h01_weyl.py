"""H1: Träger C^5 = Sym^4(C^2) (Spin 2), Marken = Nullstellen einer binären Quartik."""
import itertools, numpy as np
from common import Checks
c = Checks("h01_weyl")
m = np.array([2, 1, 0, -1, -2])
Jz = np.diag(m).astype(complex)
Jp = np.zeros((5, 5), complex)
for i in range(1, 5): Jp[i - 1, i] = np.sqrt(6 - m[i] * (m[i] + 1))   # j = 2
Jm = Jp.conj().T
c.ok("Spin 2 Darstellung: [J+, J-] = 2 Jz", np.allclose(Jp @ Jm - Jm @ Jp, 2 * Jz))
R = np.diag(np.exp(1j * np.pi * m))
c.ok("Halbdrehung exp(i pi Jz) hat 3 gerade und 2 ungerade Richtungen", sorted(np.round(np.diag(R).real)) == [-1, -1, 1, 1, 1])
Y = np.eye(5) / 12 - 5 * R.real / 12
c.ok("Hyperladung diag(-1/3 auf gerade, +1/2 auf ungerade) = 1/12 - (5/12) R, spurfrei",
     np.allclose(sorted(np.diag(Y)), [-1/3, -1/3, -1/3, .5, .5]) and abs(np.trace(Y)) < 1e-12)
def commutant_dim(mats, n=5):
    rows = [np.kron(M, np.eye(n)) - np.kron(np.eye(n), M.T) for M in mats]
    return n * n - np.linalg.matrix_rank(np.vstack(rows), tol=1e-9)
cd = {k: commutant_dim([np.diag(np.exp(2j * np.pi * m / k))]) - 1 for k in range(2, 7)}
c.ok("Zentralisator der Halbdrehung in su(5) hat Dimension 12 = su(3)+su(2)+u(1); andere Drehwinkel geben das nicht",
     cd[2] == 12 and all(v != 12 for k, v in cd.items() if k != 2), cd)
c.ok("Spin 2 ist irreduzibel: nur Vielfache der Eins vertauschen mit der ganzen SL(2) Wirkung", commutant_dim([Jz, Jp, Jm]) == 1)
gens = []
for a, b in [(0, 2), (2, 4), (0, 4), (1, 3)]:
    E = np.zeros((5, 5), complex); E[a, b] = 1; gens += [E + E.T, 1j * (E - E.T)]
gens += [np.diag([1, 0, -1, 0, 0]).astype(complex), np.diag([0, 0, 1, 0, -1]).astype(complex), np.diag([0, 1, 0, -1, 0]).astype(complex), Y.astype(complex)]
c.ok("die zwölf Erzeuger von S(U(3)xU(2)) vertauschen mit der Halbdrehung", all(np.allclose(G @ R, R @ G) for G in gens) and len(gens) == 12)
c.ok("von den drei Lorentz Erzeugern vertauscht nur Jz mit der Hyperladung: Y ist kein Element der Lie Algebra der Drehungen",
     np.allclose(Jz @ Y, Y @ Jz) and not np.allclose(Jp @ Y, Y @ Jp))
c.ok("STARKE FORM FÄLLT: kein Erzeuger von SU(3) oder SU(2) vertauscht mit J+; Eichgruppe und Lorentzgruppe können nicht beide auf demselben C^5 wirken und vertauschen",
     all(not np.allclose(G @ Jp, Jp @ G) for G in gens[:11]))
# Invarianten der Quartik a z^4 + 4b z^3 + 6c z^2 + 4d z + e
def IJ(a, b, cc, d, e): return a * e - 4 * b * d + 3 * cc * cc, a * cc * e + 2 * b * cc * d - a * d * d - b * b * e - cc ** 3
I, J = IJ(1, 0, 0, 0, -1)
c.ok("z^4 - 1: I = -1, J = 0 (harmonisch, j = 1728)", abs(I + 1) < 1e-12 and abs(J) < 1e-12)
def stab(pts):
    n = 0
    for perm in itertools.permutations(range(4)):
        z, w = pts, [pts[i] for i in perm]
        A = np.array([[z[k], 1, -w[k] * z[k], -w[k]] for k in range(4)])
        if abs(np.linalg.det(A)) < 1e-9: n += 1
    return n
sq = [1, 1j, -1, -1j]; gen = [0.3, 1.7 + .4j, -2.1j, 5 + 1j]
w3 = np.exp(2j * np.pi / 3); tet = [0, np.sqrt(2), np.sqrt(2) * w3, np.sqrt(2) * w3 ** 2]
tet = [0] + [np.sqrt(2) * w3 ** k for k in range(3)]; tet = [np.inf if False else t for t in tet]
# Tetraeder: Nullstellen von z(z^3 - 2 sqrt 2) zusammen mit ... benutze z^4 + 2 sqrt(3) i z^2 + 1 (I = 0)
tet = np.roots([1, 0, 2j * np.sqrt(3), 0, 1])
I_t, J_t = IJ(1, 0, 2j * np.sqrt(3) / 6, 0, 1)
orders = {"generisch": stab(gen), "Quadrat": stab(sq), "Tetraeder": stab(list(tet))}
c.ok("Möbius Symmetrien der Markenmenge: generisch 4, Quadrat 8, Tetraeder 12", orders == {"generisch": 4, "Quadrat": 8, "Tetraeder": 12} and abs(I_t) < 1e-9, orders)
c.done({"zentralisator_dim": cd, "moebius": orders,
        "urteil": "Reduzierte Form exakt: Eichalgebra des Standardmodells = Zentralisator der Halbdrehung der Haupt SL(2) auf dem Träger. Starke Form (derselbe Fünfer ist zugleich Raumzeit Spin 2 mit vertauschender Eichwirkung) ausgeschlossen."})
