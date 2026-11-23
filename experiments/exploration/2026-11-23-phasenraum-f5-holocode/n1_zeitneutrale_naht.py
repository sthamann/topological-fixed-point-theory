"""N: Ursprung von P1 — die zeitneutrale Fermionnaht.

Prinzip: Eine geschlossene, orientierte, flache 2D-Naht mit freien Fermionen zeichnet keine
Zeitrichtung aus. Rekonstruiert man (Osterwalder-Schrader) entlang jedes der beiden Basiszyklen
ein Quantensystem, entsteht dasselbe orientierte thermische System, und jede Rekonstruktion ist
ein echter KMS-Zustand (Fermionen in der Zeit antiperiodisch, keine (-1)^F Einfügung) mit
invertierbarem Hamiltonoperator (keine Nullmoden).

Prüfungen
 1. Für Gitter (1, tau), tau im Fundamentalbereich, und alle vier Spinstrukturen: Rekonstruktion
    A (Raum entlang 1, Zeit entlang tau) gibt (beta/L, Drall/L) = (Im tau, Re tau); Rekonstruktion
    B (orientierte Vertauschung (tau, -1)) gibt (Im tau/|tau|^2, -Re tau/|tau|^2). Gesucht:
    gleiche Daten, gleiche Randbedingungen, KMS und invertierbar.
 2. Die Zeitneutralitätsgruppe und die Normierung c = (Diracrang/(4 pi)) / |Gruppe|.
 3. Der ausgewählte Zustand: Spektrum, Lücke, Modularinvarianz der Zustandssumme am Fixpunkt,
    Viertelzählung als exakte Bahnzählung.
"""
import numpy as np
from mpmath import mp, jtheta, exp, pi, mpf, mpc, sqrt, qp

mp.dps = 30


def recon(tau):
    a = (tau.imag, tau.real)
    t2 = -1 / tau
    b = (t2.imag, t2.real)
    return a, b


def spin_ok(s_space, s_time):
    """KMS needs antiperiodic time; invertible massless Dirac needs antiperiodic space."""
    return s_time == "A" and s_space == "A"


print("1. Suche über Moduli und Spinstrukturen\n")
solutions, rhombic = [], []
xs = np.linspace(-0.5, 0.5, 201)
ys = np.linspace(0.8, 2.5, 341)
for x in xs:
    for y in ys:
        tau = complex(x, y)
        if abs(tau) < 1 - 1e-12:
            continue
        (bA, dA), (bB, dB) = recon(tau)
        same_or = abs(bA - bB) < 1e-9 and abs(dA - dB) < 1e-9
        same_unor = abs(bA - bB) < 1e-9 and abs(dA + dB) < 1e-9
        for s1 in "AP":
            for s2 in "AP":
                # A: space = cycle 1 (s1), time = cycle 2 (s2); B: space = cycle 2, time = cycle 1
                okA, okB = spin_ok(s1, s2), spin_ok(s2, s1)
                if same_or and okA and okB and s1 == s2:
                    solutions.append((round(x, 6), round(y, 6), s1 + s2))
                if same_unor and okA and okB and s1 == s2 and not same_or:
                    rhombic.append((round(x, 4), round(y, 4)))
print(f"   orientiert zeitneutral, KMS und invertierbar: {sorted(set(solutions))}")
print(f"   nur bis auf Parität zeitneutral (rhombische Familie |tau| = 1): {len(set(rhombic))} Gitterpunkte, "
      f"z. B. {sorted(set(rhombic))[:3]}")
hexa = complex(0.5, np.sqrt(3) / 2)
print(f"   hexagonaler Punkt tau = e^(i pi/3): Rekonstruktionen {recon(hexa)} -> Drall wechselt das Vorzeichen,"
      f" orientiert nicht zeitneutral")

print("\n2. Zeitneutralitätsgruppe und Normierung")
S = np.array([[0, -1], [1, 0]])
order = next(k for k in range(1, 9) if np.array_equal(np.linalg.matrix_power(S, k), np.eye(2, dtype=int)))
rank = 2
for name, G in (("quadratisch, orientiert (S)", order), ("rhombisch, nur mit Parität", 2), ("ohne Neutralität", 1)):
    print(f"   {name:30s}: |Gruppe| = {G}, c = 2/(4π·{G}) = {rank / (4 * np.pi * G):.12f}"
          + ("   = 1/(8π) = c3" if G == 4 else ""))

print("\n3. Der ausgewählte Zustand am Fixpunkt tau = i, Spinstruktur (A, A)")
# single-particle spectrum on the antiperiodic unit circle
print("   Einteilchenenergien 2π|n + 1/2| auf dem Kreis der Länge 1: Lücke π, β/L = 1 (selbstduale Temperatur)")
# Dirac partition function |theta_3 / eta|^2 is invariant under tau -> -1/tau at tau = i (check numerically)


def Z(tau):
    q = exp(2j * pi * tau)
    nome = exp(1j * pi * tau)
    eta = exp(1j * pi * tau / 12) * qp(q)
    return abs(jtheta(3, 0, nome) / eta) ** 2


z1, z2 = Z(mpc(0, 1)), Z(-1 / mpc(0, 1))
print(f"   Z(i) = {mp.nstr(z1, 15)},  Z(-1/i) = {mp.nstr(z2, 15)}  (gleich: {abs(z1 - z2) < mpf(10) ** -20})")
zA, zB = Z(mpc(0.1, 1.2)), Z(-1 / mpc(0.1, 1.2))
print(f"   Kontrolle abseits: Z(0.1+1.2i) = {mp.nstr(zA, 12)}, Z(-1/(0.1+1.2i)) = {mp.nstr(zB, 12)}"
      f"  (Z selbst ist S-invariant; zeitneutral ist erst die Gleichheit der Rekonstruktionsdaten)")
# quarter counting = exact orbit counting: no momentum in (Z+1/2)^2 is fixed by the quarter turn
ks = [(a + 0.5, b + 0.5) for a in range(-6, 6) for b in range(-6, 6)]
fixed = [k for k in ks if abs(-k[1] - k[0]) < 1e-12 and abs(k[0] - k[1]) < 1e-12]
print(f"   Impulse (Z+1/2)^2 mit Fixpunkt unter der Vierteldrehung: {len(fixed)} -> alle Bahnen haben Länge 4,"
      f"\n   die Viertelzählung ist exakt die Zählung der zeitneutralen Bahnen")
s = 0.002
total = sum(np.exp(-4 * np.pi ** 2 * s * (a * a + b * b))
            for a in np.arange(-80, 80) + 0.5 for b in np.arange(-80, 80) + 0.5)
print(f"   s · (Bahnzahl der Diracmoden) bei s = {s}: {s * 2 * total / 4:.12f}  (1/(8π) = {1 / (8 * np.pi):.12f})")
assert sorted(set(solutions)) == [(0.0, 1.0, "AA")]
