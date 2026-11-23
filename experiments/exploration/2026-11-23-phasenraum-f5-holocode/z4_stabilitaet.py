"""Z4: Ist die Lage tau = i stabil?

Z3 fand: Für die chirale E8 Naht ist tau = i die einzige stationäre Form, aber F = -log|Z| ist dort ein
Sattel (stabil gegen Verdrillung, "instabil" gegen Streckung). Hier wird geprüft, ob das eine
Instabilität ist oder die gewöhnliche Struktur zweier konjugierter Paare.

In Rekonstruktion A (Raum L, Zeit beta) ist tau = x + i y mit y = beta/L. x ist die Verdrillung und
koppelt an den Impuls; y ist die inverse Temperatur und koppelt an die Energie. Z = Tr exp(-beta H + i theta P).

 1. Verdrillung: Bei positivem Spektrum ist |Z(x + i y)| <= Z(i y), Gleichheit nur bei x ganz.
    Zero Verdrillung ist also für jede unitäre Naht ein echtes Maximum von |Z| (stabil).
 2. Temperatur: log Z ist in beta konvex, denn die zweite Ableitung ist die Energievarianz.
    -log Z ist längs y daher für jede unitäre Naht konkav. Ein Minimum über die Temperatur gibt es
    nie; die Temperatur ist eine Legendre-Variable, die durch die Energiebedingung festgelegt wird.
    Stabilität heißt hier positive Wärmekapazität, also eine konkave mikrokanonische Entropie S(E).
 3. Was die Energiebedingung bei c = 8 bedeutet: im Mittel oder als Operatorgleichung.
"""
import json
from pathlib import Path

import numpy as np

OK = []


def check(name, cond):
    OK.append((name, bool(cond)))
    return "ja" if cond else "NEIN"


NMAX = 40


def e8_coeffs(nmax=NMAX):
    """chi_E8 = q^(-1/3) * sum a_n q^n with integer a_n (E4 / prod (1 - q^n)^8)."""
    sigma3 = [0] + [sum(d ** 3 for d in range(1, n + 1) if n % d == 0) for n in range(1, nmax + 1)]
    e4 = [1] + [240 * sigma3[n] for n in range(1, nmax + 1)]
    inv = [1] + [0] * nmax                                   # prod (1 - q^n)^(-8)
    for n in range(1, nmax + 1):
        for _ in range(8):
            for k in range(n, nmax + 1):
                inv[k] += inv[k - n]
    return [sum(e4[i] * inv[n - i] for i in range(n + 1)) for n in range(nmax + 1)]


a = np.array(e8_coeffs(), dtype=float)
lv = np.arange(NMAX + 1)
print(f"Koeffizienten von chi_E8: {[int(v) for v in a[:6]]} (1, 248, 4124, 34752, 213126, 1057504 erwartet)")
check("E8 Koeffizienten", list(a[:6].astype(int)) == [1, 248, 4124, 34752, 213126, 1057504])
C24 = 8 / 24


def stats(y, coeffs=a, shift=C24):
    w = coeffs * np.exp(-2 * np.pi * lv * y)
    m = np.sum(lv * w) / np.sum(w)
    v = np.sum(lv ** 2 * w) / np.sum(w) - m * m
    return m - shift, v, np.log(np.sum(w)) + 2 * np.pi * shift * y


def Z(tau, coeffs=a, shift=C24):
    return np.exp(-2j * np.pi * shift * tau) * np.sum(coeffs * np.exp(2j * np.pi * lv * tau))


print("\n1. Verdrillung (x) bei y = 1\n")
xs = np.linspace(-0.5, 0.5, 201)
ratio = np.array([abs(Z(x + 1j)) / abs(Z(1j)) for x in xs])
print(f"   |Z(x + i)|/Z(i): Maximum {ratio.max():.12f} bei x = {xs[ratio.argmax()]:+.3f}; "
      f"bei x = 0.1: {abs(Z(0.1 + 1j)) / abs(Z(1j)):.6f}; bei x = 0.5: {abs(Z(0.5 + 1j)) / abs(Z(1j)):.6f}")
check("Verdrillung null ist Maximum", xs[ratio.argmax()] == 0 and np.all(ratio[np.abs(xs) > 1e-9] < 1))
h = 1e-4
Fxx = -(np.log(abs(Z(h + 1j))) - 2 * np.log(abs(Z(1j))) + np.log(abs(Z(-h + 1j)))) / h ** 2
print(f"   F_xx = {Fxx:+.5f} > 0: stabil. Grund: positive Koeffizienten, |Tr e^(-beta H) e^(i theta P)| <= Tr e^(-beta H).")

print("\n2. Temperatur (y): Energie, Varianz, Wärmekapazität\n")
for y in (0.8, 1.0, 1.25):
    E, var, lz = stats(y)
    print(f"   beta/L = {y:4.2f}: <L0> - c/24 = {E:+.6f}, Var(L0) = {var:.6f}, d^2 log Z/dy^2 = (2 pi)^2 Var = "
          f"{(2 * np.pi) ** 2 * var:.5f}")
E1, var1, lz1 = stats(1.0)
check("Energie null bei y = 1", abs(E1) < 1e-12)
Fyy = -(2 * np.pi) ** 2 * var1
print(f"   Bei y = 1 ist <L0> = c/24 = 1/3 genau: Casimir- und Anregungsenergie heben sich auf.")
print(f"   F_yy = -(2 pi)^2 Var = {Fyy:+.5f}; Z3 fand -9.57854: gleich {check('F_yy = -Varianz', abs(Fyy + 9.57854) < 1e-4)}")
ys = np.linspace(0.3, 3.0, 271)
logz = np.array([stats(y)[2] for y in ys])
second = np.diff(logz, 2) / (ys[1] - ys[0]) ** 2
print(f"   log Z(y) ist auf [0.3, 3] konvex: kleinste zweite Ableitung {second.min():.4f} > 0 "
      f"{check('log Z konvex', second.min() > 0)}")
print("   -log Z ist längs der Temperatur also für jede unitäre Naht konkav. Das ist positive Wärmekapazität,")
print("   keine Instabilität; ein Minimum über die Temperatur gibt es grundsätzlich nicht.")

print("\n   Mikrokanonische Entropie S(E) = min_beta [log Z(beta) + beta E] (Legendre):")
Es = np.linspace(-0.25, 0.25, 11)
S = []
for E in Es:
    vals = [stats(y)[2] + 2 * np.pi * y * E for y in np.linspace(0.2, 6, 2000)]
    S.append(min(vals))
S = np.array(S)
dS2 = np.diff(S, 2)
print("   E (Einheiten 2 pi/L): " + " ".join(f"{e:+.2f}" for e in Es))
print("   S(E):                 " + " ".join(f"{s:5.3f}" for s in S))
print(f"   S(0) = {S[5]:.6f}, log 12 = {np.log(12):.6f} (chi_E8(i) = 12); S konkav: "
      f"{check('S konkav', np.all(dS2 < 0))}")
check("S(0) = log 12", abs(S[5] - np.log(12)) < 1e-4)
print("   Bei E = 0 ist die Temperatur beta = L, also genau der Punkt aus Z3. Die Lage ist der stabile")
print("   mikrokanonische Gleichgewichtszustand bei Energie null.")

print("\n3. Die Energiebedingung: im Mittel oder als Operatorgleichung?\n")
print(f"   c = 8: Niveaus L0 - c/24 liegen bei n - 1/3; ein Zustand mit Energie genau null existiert nicht.")
print(f"   Also gilt E = 0 bei c = 8 nur im Mittel (thermischer Zustand), nicht als H|psi> = 0.")
j_coeffs = np.convolve(np.convolve(a, a), a)[:NMAX + 1]       # chi^3 = j: shift c/24 = 1
print(f"   c = 24 (E8 hoch 3, Z = j): Niveau L0 - 1 = 0 hat {int(j_coeffs[1])} Zustände (3 x 248 Ströme).")
check("E8^3 Nullniveau 744", int(j_coeffs[1]) == 744)
E24, v24, _ = stats(1.0, j_coeffs, 1.0)
print(f"   Auch bei c = 24 ist im thermischen Zustand bei y = 1: <L0> - 1 = {E24:+.2e}")
check("E8^3 Energie null bei i", abs(E24) < 1e-12)
print("   Die Gabel aus Z1 erscheint hier wieder: Gilt die Energiebedingung als Operatorgleichung")
print("   (Wheeler-DeWitt-artig), braucht es c durch 24 teilbar (F'); gilt sie im Mittel des Zustands,")
print("   genügt c = 8 mit dem 3D Zufluss (F). TFPT nutzt ein E8 mit 248 = 78 + 8 + 2 x 27 x 3")
print("   (SEAM.EQUIV.02), also die Lesart im Mittel.")

print("\nZusammenfassung")
for n, c in OK:
    print(f"   {'bestanden' if c else 'NICHT bestanden'}: {n}")
Path(__file__).with_name("z4_ergebnis.json").write_text(json.dumps({
    "checks": {n: c for n, c in OK}, "Fxx": float(Fxx), "Fyy": float(Fyy), "var_L0_at_i": float(var1),
    "S0": float(S[5]), "e8cubed_zero_level": int(j_coeffs[1])}, indent=1))
assert all(c for _, c in OK)
