"""Z2: Gegentest für Prämisse D, die Verschränkungsdichte des Gittervakuums.

Jacobsons Zustandsgleichung braucht die Verschränkungsentropie je Fläche eta, und 8 pi G = 2 pi/eta.
TFPT setzt eta = 1/|mu4| = 1/4 (GRAV.ENTROPY.EQUILIBRIUM.01). In H6 ist die Raumzeit das Gitter Z^4 mit
Gitterabstand eins. Dann ist die Gitterskala genau dann die Planckskala, wenn das Vakuum je
Gitterfläche die Entropie kappa = 1/4 trägt. Diese Zahl ist auf dem Gitter berechenbar.

Rechnung: freies naives Dirac Fermion auf Z^3 (Hamiltonsch), Schnitt entlang einer Gitterebene.
Für jeden Querimpuls k_perp zerfällt das Problem in zwei Kopien einer Kette
h(k) = sigma_x sin k + sigma_z mu mit mu^2 = sin^2 k1 + sin^2 k2 + m^2. Dann ist
kappa = 2 * Integral d^2k_perp/(2 pi)^2 S_Kette(mu).
"""
import json
from pathlib import Path

import numpy as np


def chain_entropy(mu, ring=400):
    """Entanglement entropy per cut of half a ring of the 2-band chain (antiperiodic ring)."""
    k = 2 * np.pi * (np.arange(ring) + 0.5) / ring
    E = np.sqrt(np.sin(k) ** 2 + mu ** 2)
    # projector on the filled band: (1 - h/E)/2 with h = sin k sigma_x + mu sigma_z
    P = np.zeros((ring, 2, 2), complex)
    P[:, 0, 0] = 0.5 * (1 - mu / E)
    P[:, 1, 1] = 0.5 * (1 + mu / E)
    P[:, 0, 1] = -0.5 * np.sin(k) / E
    P[:, 1, 0] = -0.5 * np.sin(k) / E
    half = ring // 2
    x = np.arange(half)
    phase = np.exp(1j * np.outer(x, k))                       # (half, ring)
    C = np.zeros((half, 2, half, 2), complex)
    for a in range(2):
        for b in range(2):
            C[:, a, :, b] = (phase * P[:, a, b]) @ phase.conj().T / ring
    nu = np.linalg.eigvalsh(C.reshape(2 * half, 2 * half))
    nu = nu[(nu > 1e-14) & (nu < 1 - 1e-14)]
    S = -np.sum(nu * np.log(nu) + (1 - nu) * np.log(1 - nu))
    return S / 2                                               # two cuts on the ring


def kappa(m=0.0, ring=200, nk=480):
    """Area coefficient per naive Dirac field; log-spaced table resolves the gapless points."""
    top = np.sqrt(2 + m * m)
    mus = np.unique(np.concatenate([np.geomspace(1e-6, 0.05, 300), np.linspace(0.05, top, 500)]))
    table = np.array([chain_entropy(mu, ring) for mu in mus])
    k = 2 * np.pi * (np.arange(nk) + 0.5) / nk
    K1, K2 = np.meshgrid(k, k)
    mu = np.sqrt(np.sin(K1) ** 2 + np.sin(K2) ** 2 + m * m)
    S = np.interp(np.log(np.maximum(mu, 1e-6)).ravel(), np.log(mus), table)
    return 2 * S.mean()                                        # per lattice area, per naive Dirac field


print("Verschränkungsdichte kappa des naiven Dirac Feldes auf Z^3 (je Gitterfläche, in nat)\n")
print(f"   Kette bei mu = 1: S = {chain_entropy(1.0, 200):.8f} (= ln2/2 = {np.log(2) / 2:.8f})")
print("   Konvergenz (masselos):")
vals = {}
for ring, nk in ((100, 240), (200, 480), (400, 480), (200, 960)):
    vals[(ring, nk)] = kappa(0.0, ring, nk)
    print(f"   Ring {ring:4d}, Querimpulse {nk:4d}^2: kappa = {vals[(ring, nk)]:.6f}")
k0 = vals[(200, 960)]
kt = k0 / 8
print(f"   pi/4 = {np.pi / 4:.6f}: Verhältnis {k0 / (np.pi / 4):.5f}, also keine Sonderzahl")
print(f"\n   Je Geschmack (8 Verdoppler im Hamiltonschen Z^3): kappa_t = {kt:.6f}; je Weylgeschmack {kt / 2:.6f}")
for m in (0.1, 0.3, 1.0):
    print(f"   Masse {m}: kappa = {kappa(m):.5f}")

print("\n   Prämisse D: 8 pi G = 2 pi/eta mit eta = kappa_gesamt/a^2, also l_P^2 = G = a^2/(4 kappa_gesamt)")
print("   und a/l_P = 2 sqrt(kappa_gesamt). Die Gitterskala ist genau dann die Planckskala, wenn kappa_gesamt = 1/4.")
print(f"   1/4 entspricht {0.25 / kt:.2f} Dirac-Geschmäckern. Szenarien (nur Fermionen, Eichfelder und Higgs kämen hinzu):")
for name, n_taste in (("16 Majorana-Weyl der Naht (8 Dirac-Geschmäcker)", 8),
                      ("drei Generationen mit nu_R, Verdoppler entfernt (24 Dirac)", 24),
                      ("drei Generationen als naive Felder (24 x 8 Geschmäcker)", 192)):
    kap = n_taste * kt
    print(f"   {name}: kappa = {kap:.3f}, a/l_P >= {2 * np.sqrt(kap):.2f}")
assert abs(k0 / (np.pi / 4) - 1) > 1e-3
print("\n   Ergebnis: Mit dem Fermioninhalt von H6 trägt das Gittervakuum weit mehr als 1/4 je Gitterfläche."
      "\n   Unter D liegt der Gitterabstand daher bei mindestens etwa drei Plancklängen; der Satz \"Gitter in"
      "\n   Planckeinheiten\" ist falsch. eta = 1/4 ist unter D die Definition der Planckeinheit (tautologisch),"
      "\n   und a/l_P = 2 sqrt(kappa_gesamt) wird eine Vorhersage, sobald der volle Feldinhalt samt Eichfeldern"
      "\n   und deren Kontaktterm feststeht. Das ist Sakharovs Artenproblem in Gitterform.")
Path(__file__).with_name("z2_ergebnis.json").write_text(json.dumps({
    "kappa_naive_dirac": float(k0), "kappa_per_dirac_taste": float(kt),
    "a_over_lP_three_generations_no_doublers": float(2 * np.sqrt(24 * kt)),
    "a_over_lP_three_generations_naive": float(2 * np.sqrt(192 * kt)),
    "tastes_for_quarter": float(0.25 / kt)}, indent=1))
