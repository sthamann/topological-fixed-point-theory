"""Zweite Erzeugerdatei, zeigt den Einbau: jede Datei gen_*.py mit Funktionen gen_*() wird von suche.py von selbst geladen."""
import math
from fractions import Fraction as Fr
from generators import Z, B, T, obj

def gen_polytope():
    data = {"Simplex 5-Zeller": ([5, 10, 10, 5], 120), "Tesserakt": ([16, 32, 24, 8], 384), "16-Zeller": ([8, 24, 32, 16], 384), "24-Zeller": ([24, 96, 96, 24], 1152), "120-Zeller": ([600, 1200, 720, 120], 14400), "600-Zeller": ([120, 720, 1200, 600], 14400),
            "Ikosaeder": ([12, 30, 20], 120), "Dodekaeder": ([20, 30, 12], 120), "Oktaeder": ([6, 12, 8], 48), "Würfel": ([8, 12, 6], 48), "Tetraeder": ([4, 6, 4], 24),
            "Gosset 4_21": ([240, 6720, 60480], 696729600), "Gosset 3_21": ([56, 756, 4032], 2903040), "Gosset 2_21": ([27, 216, 720], 51840), "Halbwürfel 5": ([16, 80, 160], 1920), "Kreuzpolytop 5": ([10, 40, 80, 80, 32], 3840)}
    for name, (f, sym) in data.items():
        yield obj(f"polytop:{name}", "Polytop", f"Polytop {name}", [Z(f[0], "Ecken"), Z(f[1], "Kanten"), Z(sym, "Ordnung der Symmetriegruppe"), Z(2 * f[1] // f[0], "Nachbarn je Ecke"), Z(len(f), "Dimension")] + [Z(x, "Seitenzahl") for x in f[2:]])

def gen_darstellungen():
    reps = {"A4 = SU(5)": [5, 10, 24, 15, 40, 45, 50], "D5 = SO(10)": [10, 16, 45, 54, 120, 126, 144, 210], "E6": [27, 78, 351, 650], "E7": [56, 133, 912, 1539], "E8": [248, 3875, 30380], "A3 = SU(4)": [4, 6, 10, 15, 20, 36, 45],
            "A2 = SU(3)": [3, 6, 8, 10, 15, 27], "D4 = SO(8)": [8, 28, 35, 56], "A1 = SU(2)": [2, 3, 4, 5], "G2": [7, 14, 27, 64, 77], "F4": [26, 52, 273, 324], "A5 = SU(6)": [6, 15, 20, 35, 21], "B4 = SO(9)": [9, 16, 36, 84, 126], "C3 = Sp(6)": [6, 14, 21]}
    for name, d in reps.items():
        yield obj(f"darst:{name}", "Darstellungen", f"kleine Darstellungen von {name}", [Z(x, "Dimension einer Darstellung") for x in d] + [T(name.split(" ")[0], "Algebra")])

def gen_modulformen():
    for k in range(2, 26, 2):
        dim = k // 12 + (0 if k % 12 == 2 else 1)
        yield obj(f"modform:k={k}", "Modulformen", f"Modulformen vom Gewicht {k} zur vollen Gruppe", [Z(k, "Gewicht"), Z(dim, "Dimension"), B(Fr(k, 12), "Grad des Bündels"), Z(dim - 1, "Spitzenformen") if dim > 1 else None])
    for m in range(3, 10): yield obj(f"cft:M({m+1},{m})", "Randtheorie", f"Minimalmodell M({m+1},{m})", [B(1 - Fr(6, m * (m + 1)), "zentrale Ladung"), Z(m * (m - 1) // 2, "Zahl der Sektoren"), Z(m, "Stufe")])
    for k in range(2, 13): yield obj(f"cft:PF{k}", "Randtheorie", f"Parafermionen Z{k}", [B(Fr(2 * (k - 1), k + 2), "zentrale Ladung"), Z(k * (k + 1) // 2, "Zahl der Sektoren"), Z(k, "Ordnung"), T(f"Z{k}", "Symmetrie")])
    for N in range(3, 6):
        for k in range(1, 7): yield obj(f"cft:SU({N})_{k}", "Randtheorie", f"SU({N}) auf Stufe {k}", [B(Fr(k * (N * N - 1), k + N), "zentrale Ladung"), Z(math.comb(N + k - 1, N - 1), "Zahl der Sektoren"), Z(k, "Stufe"), Z(k + N, "Nenner"), T(f"A{N-1}", "Algebra")])
    for name, t, c, n in (("SO(10) Stufe 1", "D5", 5, 4), ("E6 Stufe 1", "E6", 6, 3), ("E7 Stufe 1", "E7", 7, 2), ("E8 Stufe 1", "E8", 8, 1)):
        yield obj(f"cft:{name}", "Randtheorie", f"Randtheorie {name}", [Z(c, "zentrale Ladung"), Z(n, "Zahl der Sektoren"), T(t, "Algebra")])
