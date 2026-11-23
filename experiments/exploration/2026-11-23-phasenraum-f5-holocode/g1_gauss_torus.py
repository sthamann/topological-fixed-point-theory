"""G: Ein Objekt für P1 und P2 — der Torus C/O_K mit komplexer Multiplikation.

Klasse der Startobjekte: flache 2-Tori, deren Gitter der Ganzheitsring O_K eines imaginär-
quadratischen Körpers ist (Klassenzahl 1: d = 1, 2, 3, 7, 11, 19, 43, 67, 163). Kein TFPT Objekt
wird vorausgesetzt; ausgewählt wird nur durch Regeln.

Regeln
  (A) Die Einheiten U von O_K wirken auf einem Torsionsstück O_K/p als die volle Einheitengruppe
      eines Primkörpers F_p (Reduktion U -> F_p^x bijektiv). Das kleinste solche p bestimmt die
      Plätze: Funktionen auf O_K/p = F_p (Schrödinger-Modell der Weil-Darstellung).
  (B) Der fermionische Teil der Weil-Darstellung von SL(2,F_p) ist ein treues SU(2) Dublett.
  (C) Es gibt eine U-invariante Spinstruktur ohne Dirac-Nullmoden.
Ausgaben
  - die Randnormierung c = covol(O_K) / (2 pi |U|) (Weylterm der Dirac-Wärmespur, rank 2, gezählt
    pro U-Bahn), zu vergleichen mit c3 = 1/(8 pi);
  - für den ausgewählten Torus: Vierteldrehung auf der 5-Torsion = Weyl-Element = Fourier,
    Eigengeraden E[2 +- i], Uhrenrelation, Wärmespur numerisch.
"""
import itertools
from fractions import Fraction

import numpy as np
import sympy as sp

exec(open(__file__.replace("g1_gauss_torus.py", "w1_weil_darstellung.py")).read().split("def main")[0])

D_LIST = [1, 2, 3, 7, 11, 19, 43, 67, 163]


def ring(d):
    """Basis (1, tau) of O_K in C, unit group order, covolume."""
    if d % 4 == 3:
        tau = (1 + 1j * np.sqrt(d)) / 2
        covol = np.sqrt(d) / 2
    else:
        tau = 1j * np.sqrt(d)
        covol = np.sqrt(d)
    units = {1: 4, 3: 6}.get(d, 2)
    return tau, units, covol


def splits_or_ramifies(p, d):
    disc = -d if d % 4 == 3 else -4 * d
    return disc % p == 0 or pow(disc % p, (p - 1) // 2, p) == 1


def first_full_unit_prime(d):
    _, U, _ = ring(d)
    p = U + 1
    return p if sp.isprime(p) and splits_or_ramifies(p, d) else None


def spin_structures(d):
    """Shifts s in (1/2) Lambda* / Lambda*; invariant under the generator of U; invertible if s not in Lambda*."""
    tau, U, _ = ring(d)
    B = np.array([[1, tau.real], [0, tau.imag]])            # columns: lattice basis in R^2
    Bd = np.linalg.inv(B).T                                  # dual basis columns
    rot = 2 * np.pi / U
    Rm = np.array([[np.cos(rot), -np.sin(rot)], [np.sin(rot), np.cos(rot)]])
    out = []
    for e in itertools.product((0, 1), repeat=2):
        s = Bd @ (np.array(e) / 2)
        diff = np.linalg.solve(Bd, Rm @ s - s)               # coordinates in the dual basis
        invariant = np.allclose(diff, np.round(diff), atol=1e-9)
        invertible = e != (0, 0)
        out.append((e, invariant, invertible))
    return out


weil_cache = {}


def weil_doublet(p):
    if p not in weil_cache:
        r = analyse(p)
        weil_cache[p] = bool(r["dim_ferm"] == 2 and r["kernel_ferm"] == 1 and r.get("det1_ferm"))
    return weil_cache[p]


print("1. Auswahl über alle CM-Tori der Klassenzahl 1\n")
print(f"   {'d':>4} {'|U|':>4} {'p (A)':>6} {'Dublett (B)':>12} {'Spinstr. (C)':>13} {'c = covol/(2π|U|)':>20}")
selected = []
for d in D_LIST:
    tau, U, covol = ring(d)
    p = first_full_unit_prime(d)
    B = weil_doublet(p) if p else False
    C = any(inv and invt for _, inv, invt in spin_structures(d))
    c = covol / (2 * np.pi * U)
    tag = "  <- c3 = 1/(8π)" if abs(c - 1 / (8 * np.pi)) < 1e-12 else ""
    print(f"   {d:>4} {U:>4} {str(p):>6} {str(B):>12} {str(C):>13} {c:>20.10f}{tag}")
    if p and B and C:
        selected.append(d)
print(f"\n   (A) und (B) und (C) erfüllt nur: d = {selected}")

print("\n2. Der ausgewählte Torus C/Z[i] und seine 5-Torsion Z[i]/5")
Mi = sp.Matrix([[0, -1], [1, 0]])                          # multiplication by i in basis (1, i)
J = sp.Matrix([[0, 1], [-1, 0]])
print(f"   Multiplikation mit i in der Basis (1, i): {Mi.tolist()}, symplektisch: "
      f"{(Mi.T * J * Mi - J) % 5 == sp.zeros(2)}, Ordnung {next(k for k in range(1, 9) if (Mi ** k) % 5 == sp.eye(2))}")
print("   Das ist das Weyl-Element w = [[0,-1],[1,0]] von SL(2,5); sein Weil-Bild ist die Fourier-"
      "\n   transformation S mit Spektrum {1, -1, -1, i, -i} (w1), also TFPTs R5 bis auf das Vorzeichen.")
for lam in range(5):
    v = [(x, y) for x in range(5) for y in range(5) if (x, y) != (0, 0)
         and ((Mi * sp.Matrix([x, y])) % 5 - lam * sp.Matrix([x, y])) % 5 == sp.zeros(2, 1)]
    if v:
        x, y = v[0]
        gauss = sp.expand(x + y * sp.I)
        # which Gaussian prime kills this line: (2+i) or (2-i)
        killer = [g for g in (2 + sp.I, 2 - sp.I) if all(
            sp.simplify(((g * (a + b * sp.I)).expand().as_real_imag()[0]) % 5) == 0 and
            sp.simplify(((g * (a + b * sp.I)).expand().as_real_imag()[1]) % 5) == 0 for a, b in [v[0]])]
        print(f"   Eigenwert {lam}: Eigengerade erzeugt von {gauss}, Torsionsgerade E[{killer[0] if killer else '?'}]")
print("   Plätze = Funktionen auf der Eigengeraden L = E[2-i] (Eigenwert 2; E[2+i] mit Eigenwert 3 ist"
      "\n   die Gegenorientierung). Translation um L ist C5, die Vierteldrehung wirkt auf L als x -> 2x = G.")
L = [(t, 3 * t % 5) for t in range(5)]                    # line spanned by 1 + 3i
gen = (1, 3)
mult_i = {pt: (int((-pt[1]) % 5), int(pt[0] % 5)) for pt in L}
coord = {pt: t for t, pt in enumerate(L)}
on_line = all(mult_i[pt] in L for pt in L)
dil = all(coord[mult_i[pt]] == (2 * coord[pt]) % 5 for pt in L)
# clock relation on the line: G T G^-1 = T^2 with T = translation by the generator
T = lambda t: (t + 1) % 5
Gm = lambda t: (2 * t) % 5
Ginv = lambda t: (3 * t) % 5
clock = all(Gm(T(Ginv(t))) == (t + 2) % 5 for t in range(5))
print(f"   Kontrolle: i bildet L = {L} auf sich ab: {on_line}; wirkt dort als t -> 2t: {dil}; "
      f"G C5 G^-1 = C5^2: {clock}")

print("\n3. Randnormierung aus der Wärmespur von C/Z[i] (Spinstruktur (-,-), Kantenlänge 1)")
for s in (0.02, 0.01, 0.005):
    ks = np.arange(-60, 60) + 0.5
    total = 2 * sum(np.exp(-4 * np.pi ** 2 * s * (kx ** 2 + ky ** 2)) for kx in ks for ky in [ks])
    total = float(np.sum(total))
    print(f"   s = {s}: s * (1/4) Tr exp(-s D^2) = {s * total / 4:.12f}   (1/(8π) = {1 / (8 * np.pi):.12f})")
assert selected == [1]
