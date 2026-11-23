"""C: Vergleich des zeitneutralen Zustands (N) mit der operativen Saat des Originals.

Original (KERNSIX20261120 Schritte 1, 3, 4; P1BIND20261120): Die Naht Sigma ist der Raum, die Zeit
läuft normal dazu. Kragen D = gamma_n (d_n + B), B = Dirac auf Sigma. Der Randzustand des Kragens
der Breite l ist rein; auf einer Fläche ist er thermisch für B bei beta = 2 l:
C_l = (1 + exp(-2 l B))^-1. Breite null gibt den flachen Kernzustand C = 1/2. l ist im Original
frei ("T chi = 1/(2 l lambda_1) frei"). Nullmoden von B bleiben für jede Breite flach; dort lebt das
aufgezeichnete E8 Wort. Modell für die Marken: B auf dem quadratischen Torus, Spinstruktur (A, A).

N auf die volle euklidische Geometrie gehoben: thermischer Kreis 2l mal Naht T^2(L, L), Fermionen
überall antiperiodisch (KMS). Zeitneutral heißt: die Rekonstruktion entlang des thermischen Kreises
und die entlang eines Nahtzyklus geben dasselbe orientierte thermische System.

Prüfungen
 1. Für welche Breite l sind die beiden Rekonstruktionen dasselbe System (gleiche dimensionslose
    Spektren beta * lambda)? Erwartet: genau 2 l = L.
 2. Was folgt dann für die offene Größe des Originals: T / lambda_1, Verschränkung der Kragenflächen.
 3. Nullmoden: wie viele hat B bei (A, A) und bei den anderen Spinstrukturen?
 4. Lage der Marken: N gegen KERNSIX Schritt 4 (Quadrat nur Minimum in der rechteckigen Familie,
    unter Verdrillung Sattel).
"""
import numpy as np

L = 1.0
M = 40
grid = np.arange(-M, M) + 0.5                                 # antiperiodic momenta (in units 2 pi / length)


def dirac_levels(lx, ly, sx="A", sy="A", count=None):
    ax = np.arange(-M, M) + (0.5 if sx == "A" else 0.0)
    ay = np.arange(-M, M) + (0.5 if sy == "A" else 0.0)
    kx, ky = np.meshgrid(2 * np.pi * ax / lx, 2 * np.pi * ay / ly)
    lev = np.sqrt(kx ** 2 + ky ** 2).ravel()
    lev = np.sort(lev)
    return lev if count is None else lev[:count]


def entropy_face(levels, beta):
    """One face of the pure collar state: occupations 1/(1+exp(beta*lam)) on +lam and -lam."""
    S = 0.0
    for lam in levels:
        for e in (lam, -lam):
            n = 1.0 / (1.0 + np.exp(np.clip(beta * e, -700, 700)))
            if 0 < n < 1:
                S += -(n * np.log(n) + (1 - n) * np.log(1 - n))
    return S


print("1. Zeitneutralität auf thermischer Kreis x Naht: für welche Breite l?")
for ell in (0.25, 0.4, 0.5, 0.6, 0.75, 1.0):
    beta_A = 2 * ell
    A = beta_A * dirac_levels(L, L, count=60)                  # time = thermal circle, space = T^2(L, L)
    B = L * dirac_levels(2 * ell, L, count=60)                 # time = a seam cycle, space = T^2(2l, L)
    d = np.max(np.abs(A - B))
    print(f"   l = {ell:4.2f}: max Abweichung der Spektren beta*lambda (60 tiefste) = {d:.3e}"
          + ("   <- dasselbe System" if d < 1e-9 else ""))

print("\n2. Folgen für die offene Größe des Originals (l = L/2)")
lam1 = dirac_levels(L, L, count=1)[0]
print(f"   lambda_1 = {lam1:.10f} = sqrt(2) pi / L: {abs(lam1 - np.sqrt(2) * np.pi) < 1e-12}")
print(f"   T / lambda_1 = 1/(2 l lambda_1) = {1 / (L * lam1):.10f} = 1/(sqrt(2) pi) = {1 / (np.sqrt(2) * np.pi):.10f}")
lev = dirac_levels(L, L)
for ell in (0.5,):
    print(f"   Verschränkung der beiden Kragenflächen bei l = L/2: S = {entropy_face(lev, 2 * ell):.10f} (endlich, ohne Abschneiden)")
for M_cut in (10, 20, 40):
    ax = np.arange(-M_cut, M_cut) + 0.5
    nmodes = 2 * len(ax) ** 2
    print(f"   Breite null (flacher Kernzustand): S = log 2 je Mode = {nmodes * np.log(2):.1f} bei {nmodes} Moden (wächst mit dem Abschneiden)")

print("\n3. Nullmoden von B auf dem quadratischen Torus")
for sx, sy in (("A", "A"), ("A", "P"), ("P", "A"), ("P", "P")):
    z = int(np.sum(dirac_levels(L, L, sx, sy) < 1e-12))
    print(f"   Spinstruktur ({sx},{sy}): {2 * z} Nullmoden (zwei Spinorkomponenten je Impuls null)")
print("   Bei (A,A) hat das geometrische B keine Nullmoden. Ein flacher Sektor für das E8 Wort bleibt nur,"
      "\n   wenn er nicht an das geometrische B koppelt (direkte Summe statt Tensorprodukt).")

print("\n4. Lage der Marken: Spektralwirkung mit lambda_1 normiert, Torus mit Fläche 1 (KERNSIX Schritt 4 nachgebaut)")


def action(tau, f=lambda x: np.exp(-x ** 2)):
    # lattice (1, tau) scaled to area 1, spin structure antiperiodic on both basis cycles
    a = 1 / np.sqrt(tau.imag)
    w1, w2 = a * 1.0, a * tau
    # dual lattice momenta k = 2 pi (n1 b1 + n2 b2) with b dual to (w1, w2), shifted by 1/2
    Wm = np.array([[w1.real, w2.real], [w1.imag, w2.imag]])
    Bd = 2 * np.pi * np.linalg.inv(Wm).T
    n = np.arange(-25, 25) + 0.5
    N1, N2 = np.meshgrid(n, n)
    k = Bd @ np.vstack([N1.ravel(), N2.ravel()])
    lev = np.sqrt((k ** 2).sum(axis=0))
    lam1 = lev.min()
    return 2 * np.sum(f(lev / lam1)), lam1


base, _ = action(1j)
print(f"   Quadrat tau = i: S = {base:.6f}")
for tau in (1.2j, 0.9j, 0.1 + 1j, -0.1 + 1j):
    s, _ = action(tau)
    print(f"   tau = {tau}: S - S(i) = {s - base:+.6f}")
print("   Rechteckige Familie: Quadrat Minimum; Verdrillung senkt S: Sattel (wie KERNSIX). N schließt die"
      "\n   Verdrillung durch die orientierte Zeitneutralität aus und braucht keine Minimierung.")
