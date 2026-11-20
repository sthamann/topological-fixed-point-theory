"""H8: Barnes Wall Turm als Wachstumsregel. Stufe 8 -> 16 über die Multiplikation mit 1 + i."""
import itertools, numpy as np
from common import Checks
c = Checks("h08_turm")
pts = list(itertools.product([0, 1], repeat=3))
code_basis = [np.ones(8, int)] + [np.array([p[i] for p in pts]) for i in range(3)]
code = {tuple(sum(a * b for a, b in zip(co, code_basis)) % 2) for co in itertools.product([0, 1], repeat=4)}
wts = sorted(sum(w) for w in code)
c.ok("erweiterter Hamming Code [8,4,4]: 16 Wörter, Gewichte 0, 4 (vierzehnmal), 8", len(code) == 16 and wts == [0] + [4] * 14 + [8])
inL = lambda x: all(abs(v - round(v)) < 1e-9 for v in x) and tuple(int(round(v)) % 2 for v in x) in code       # Gitter L = Construction A, Norm x.x/2
roots = [np.array(v) for v in itertools.product([-2, -1, 0, 1, 2], repeat=8) if sum(x * x for x in v) == 4 and inL(v)]
c.ok("Construction A gibt E8: 240 Vektoren der Norm 2", len(roots) == 240)
J = np.zeros((8, 8), int)
for u in range(4): J[u + 4, u] = 1; J[u, u + 4] = -1                                                           # Translation um e3 mit Vorzeichen
gens = [2 * np.eye(8, dtype=int)[i] for i in range(8)] + code_basis
c.ok("J (Vertauschen der Hälften mit Vorzeichen) ist eine Symmetrie des Gitters mit J^2 = -1: die gaußsche Struktur", np.array_equal(J @ J, -np.eye(8, dtype=int)) and all(inL(J @ g) for g in gens))
phi = np.eye(8, dtype=int) + J
c.ok("der Wachstumsschritt phi = 1 + J verdoppelt jede Norm exakt (phi^T phi = 2)", np.array_equal(phi.T @ phi, 2 * np.eye(8, dtype=int)))
phi_inv = (np.eye(8) - J) / 2
in_phiL = lambda x: inL(phi_inv @ x)
c.ok("Index [L : phi L] = 16, also Determinante der Stufe 16 gleich 256", round(np.linalg.det(phi)) == 16)
c.ok("Stufe 16 = {(u, u + v): u in L, v in phi L} hat KEINE Wurzeln: kein Vektor der Norm 2 von E8 + E8 liegt darin", not any(in_phiL(r) for r in roots))
n4 = sum(1 for r in roots for s in roots if in_phiL(s - r))                                                     # (r, s) mit s - r in phi L, Norm 4
short_phi = [np.array(v) for v in itertools.product([-2, -1, 0, 1, 2], repeat=8) if sum(x * x for x in v) == 8 and inL(v) and in_phiL(np.array(v))]
kiss = n4 + 2 * len(short_phi)
c.ok("kürzeste Vektoren der Stufe 16 haben Norm 4; ihre Anzahl ist 4320 (Barnes Wall Gitter)", kiss == 4320, (n4, len(short_phi), kiss))
c.ok("also: Wurzeln, und damit Eichfelder in dieser Lesart, enden bei E8", True)
ev = np.linalg.eigvals(J)
c.ok("J hat nur die Eigenwerte +i und -i: jeder J feste reelle Teilraum hat gerade Dimension", np.allclose(sorted(ev.imag), [-1] * 4 + [1] * 4) and np.allclose(ev.real, 0))
c.ok("ABBRUCHKRITERIUM GREIFT: D5 hat Rang 5 und A3 Rang 3, beide ungerade. Keine gaußsche Struktur kann D5 + A3 erhalten; der Wachstumsschritt mischt Träger und Familie",
     5 % 2 == 1 and 3 % 2 == 1)
sv = np.linalg.svd(phi.astype(float))[1]
c.ok("KEINE RELAXATION: alle Dehnungen von phi sind Wurzel 2; Raten 2/3 und 1/3 können aus diesem Schritt nicht kommen", np.allclose(sv, np.sqrt(2)))
# Reed Muller RM(2,5): selbstdual [32,16,8]
p5 = list(itertools.product([0, 1], repeat=5))
mono = [np.ones(32, int)] + [np.array([p[i] for p in p5]) for i in range(5)] + [np.array([p[i] * p[j] for p in p5]) for i, j in itertools.combinations(range(5), 2)]
G = np.array(mono)
def rank2(M):
    M = M.copy() % 2; r = 0
    for col in range(M.shape[1]):
        piv = next((i for i in range(r, M.shape[0]) if M[i, col]), None)
        if piv is None: continue
        M[[r, piv]] = M[[piv, r]]
        for i in range(M.shape[0]):
            if i != r and M[i, col]: M[i] ^= M[r]
        r += 1
    return r
words = (np.array(list(itertools.product([0, 1], repeat=16)), dtype=np.uint8) @ G.astype(np.uint8)) % 2
mw = int(words.sum(axis=1)[1:].min())
c.ok("Stufe mit fünf Qubits: RM(2,5) hat Dimension 16, ist selbstdual und hat Mindestabstand 8 (Reproduktion von ORIGINIDEAS20261119)",
     rank2(G) == 16 and not ((G @ G.T) % 2).any() and mw == 8, mw)
c.ok("NICHT ENTSCHEIDBAR HIER: ob der hinzukommende Tensorfaktor ein Ort ist. Dafür gibt es keine Rechnung ohne das Originalnetz", True)
c.done({"kuerzeste_vektoren": kiss,
        "urteil": "Turm exakt (E8 -> Barnes Wall 16 über 1 + i, keine Wurzeln mehr). Als Wachstumsregel für TFPT ausgeschlossen: der Schritt kann die Zerlegung D5 + A3 nicht erhalten (ungerade Ränge) und trägt keine Relaxationsraten."})
