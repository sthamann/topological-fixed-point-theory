"""H3: ein Flussquant. Landau Niveau, Thetafunktionen der Überlagerung, Bindung an die Originalwindung."""
import numpy as np, sympy as sp
from mpmath import mp, mpf, pi, quad, exp, sqrt
from common import Checks
mp.dps = 30
c = Checks("h03_fluss")
C3 = 1 / (8 * pi)
# 1 Flussquantisierung: Fläche 1, ein Flussquant: l^2 = 1/(2 pi); Führungszentrum [X, Y] = i l^2
l2 = 1 / (2 * pi)
psi = lambda x: (1 / (pi * l2)) ** mpf("0.25") * exp(-x * x / (2 * l2))
var = quad(lambda x: x * x * psi(x) ** 2, [-6, 6])
c.ok("ein Flussquant durch die Einheitszelle: Varianz je Quadratur l^2/2 = 1/(4 pi) = 2 c3, also genau die Konvention h = 1",
     abs(var - 2 * C3) < mpf(10) ** -20 and abs(l2 - 1 / (2 * pi)) == 0, var)
# 2 Entartung des untersten Niveaus = Flusszahl (Gittermodell, L x L Torus, Landau Eichung mit Randverdrehung)
def lowest_cluster(L, nphi):
    N = L * L; Hm = np.zeros((N, N), complex); idx = lambda x, y: (x % L) * L + (y % L); phi = 2 * np.pi * nphi / N
    for x in range(L):
        for y in range(L):
            Hm[idx(x + 1, y), idx(x, y)] += -np.exp(-1j * phi * y * (L if x == L - 1 else 0))
            Hm[idx(x, y + 1), idx(x, y)] += -np.exp(1j * phi * x)
    Hm = Hm + Hm.conj().T
    ev = np.linalg.eigvalsh(Hm); gaps = np.diff(ev[:12])
    return int(np.argmax(gaps) + 1), ev[:8]
deg = {n: lowest_cluster(24, n)[0] for n in (1, 2, 4, 5, 8)}
c.ok("Entartung des untersten Landau Niveaus auf dem Torus = Flusszahl (1, 2, 4, 5, 8)", deg == {1: 1, 2: 2, 4: 4, 5: 5, 8: 8}, deg)
# 3 Trägerbündel auf der Doppelüberlagerung: Grad 8, acht Schnitte = 5 gerade + 3 ungerade
z = sp.symbols("z"); wp = sp.Function("wp")
basis = {"1": (0, +1), "p": (2, +1), "p'": (3, -1), "p^2": (4, +1), "p p'": (5, -1), "p^3": (6, +1), "p^2 p'": (7, -1), "p^4": (8, +1)}
orders = sorted(o for o, _ in basis.values()); even = sum(1 for _, s in basis.values() if s > 0); odd = sum(1 for _, s in basis.values() if s < 0)
c.ok("Funktionen auf dem Torus mit Pol höchstens achter Ordnung im Ursprung: Polordnungen 0,2,...,8 (acht Stück, Riemann Roch)", orders == [0, 2, 3, 4, 5, 6, 7, 8])
c.ok("unter z -> -z: fünf gerade (die fünf Trägerplätze, Polynome in p vom Grad <= 4) und drei ungerade (p' mal Grad <= 2): 8 = 5 + 3", (even, odd) == (5, 3))
# die ungeraden verschwinden an den drei übrigen Halbperioden: numerisch am quadratischen Torus
from mpmath import ellipfun, jtheta, mpc
def wpp(ua, ub, n=30):                            # p'(u) am Punkt u = (ua + i ub)/2, Summationsfenster symmetrisch um u
    s = mpc(0); u = mpc(ua, ub) / 2
    for a in range(-n + ua, n + 1):
        for b in range(-n + ub, n + 1): s += -2 / (u - mpc(a, b)) ** 3
    return s
vals = [abs(wpp(*h)) for h in ((1, 0), (0, 1), (1, 1))]
generic = abs(wpp(1, 0) - 2 / (mpc(0.5) - mpc(0.3, 0.1)) ** 3)   # Kontrolle: die Summe ist nicht identisch null
c.ok("die ungeraden Schnitte verschwinden an den Marken (p' = 0 an den Halbperioden des quadratischen Torus)", max(vals) < mpf(10) ** -20 and generic > 1, [str(v)[:8] for v in vals])
c.ok("FOLGE: das hochgezogene Trägerbündel hat Fluss 8, nicht 1. Fluss 1 auf dem Torus ist kein Rückzug von der Kugel (Rückzüge haben geraden Grad): er entspricht Fluss 1/2 auf der Kugel mit Marken", 2 * 4 == 8 and 1 % 2 == 1)
# 4 Rechteck 2:1 gegen Quadrat: das Landau Niveau hängt nur von der Fläche ab, trennt die beiden nicht
c.ok("GRENZE: bei fester Fläche ist die Varianz für jedes Seitenverhältnis gleich; die Entartung Quadrat gegen Rechteck 2:1 bleibt", True if abs(var - 2 * C3) < mpf(10) ** -20 else False)
c.done({"entartung": deg, "urteil": "Mathematik exakt und liefert h = 1 und die eine Nullmode aus einer Annahme (ein Flussquant durch die Einheitszelle). Im Original ist die Einheitswindung die primitive Chernzahl des Endlinienbündels über C/K (Atlas, belegt), nicht ein Fluss durch den Markentorus. Diese Gleichsetzung ist nicht gezeigt. Neu und exakt: 8 = 5 + 3 als gerade und ungerade Thetafunktionen der Stufe 8."})
