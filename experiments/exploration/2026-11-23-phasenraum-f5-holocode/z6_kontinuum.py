"""Z6: Kontinuumslimes (K16) des H6 Gitters.

Nach Z2 ist das Lipschitzgitter physisch: Der Abstand ist mindestens etwa drei Plancklängen. Verlangt
ist daher nicht a -> 0, sondern Universalität im Infraroten: Bei Energien weit unter 1/a muss das
Kontinuum mit Lorentzsymmetrie herauskommen, ohne Feinabstimmung. Geprüft wird:
 1. Lichtgeschwindigkeit: Auf dem isotropen Gitter (alle vier Achsen gleichwertig) ist sie für jede Art
    exakt eins; die führende Abweichung ist Dimension 6 und richtungsabhängig.
 2. Ohne Gleichwertigkeit der Achsen (Zeit ausgezeichnet, Arten mit verschiedener Zeitsprungstärke nu)
    bekommt das induzierte Photon eine eigene Geschwindigkeit: Lorentzverletzung bei Dimension 4.
    Gerechnet über die Ein-Schleifen-Vakuumpolarisation naiver Gitterfermionen.
 3. Isotropes Gitter: Die Vakuumpolarisation ist bei Ordnung q^2 richtungsunabhängig (O(4)).
 4. Reflexionspositivität: Transfermatrix naiver Fermionen; T hat Eigenwerte +-e^(-E), T^2 ist positiv.
 5. Größenordnungen: Lorentzverletzung bei a = 3 .. 9 l_P und der Landau-Pol von U(1)_Y.
"""
import json
from pathlib import Path

import numpy as np

OK = []


def check(name, cond):
    OK.append((name, bool(cond)))
    return "ja" if cond else "NEIN"


print("1. Dispersion auf dem isotropen Gitter (Gitterabstand 1)\n")


def e_naive(p):
    return np.arcsinh(np.sqrt(np.sum(np.sin(p) ** 2)))


def e_scalar(p):
    return 2 * np.arcsinh(np.sqrt(np.sum(np.sin(np.asarray(p) / 2) ** 2)))


for name, ef in (("naives Fermion", e_naive), ("Skalar (Laplace)", e_scalar)):
    out = []
    for d in (np.array([1, 0, 0]), np.array([1, 1, 0]) / np.sqrt(2), np.array([1, 1, 1]) / np.sqrt(3)):
        ps = np.array([0.01, 0.02, 0.04])
        v = np.array([ef(p * d) / p for p in ps])
        c6 = np.polyfit(ps ** 2, v - 1, 1)[0]
        out.append((v[0], c6))
    print(f"   {name}: E/|p| bei p -> 0 = " + ", ".join(f"{o[0]:.6f}" for o in out)
          + "; Koeffizient von p^2 (Achse, Flächen-, Raumdiagonale) = " + ", ".join(f"{o[1]:+.4f}" for o in out))
    check(f"c = 1 für {name}", all(abs(o[0] - 1) < 1e-3 for o in out))
print("   Lichtgeschwindigkeit eins für jede Art und Richtung; die Abweichung ist ~ p^2 (Dimension 6)")
print("   und hängt von der Richtung ab (Würfelsymmetrie, nicht O(3)).")

print("\n2./3. Ein-Schleifen-Vakuumpolarisation naiver Fermionen, transversale Komponente Pi_33\n")
N = 36
k = 2 * np.pi * (np.arange(N) + 0.5) / N
P0, P1, P2, P3 = np.meshgrid(k, k, k, k, indexing="ij")


def pi33(q, nu=1.0, M=0.5):
    """Pi_33(q) - Pi_33(0) for a naive Dirac fermion with time hopping nu; q in the (0,1) plane."""
    def s(p0, p1):
        return nu * np.sin(p0), np.sin(p1)
    s3 = np.sin(P3)
    S2 = np.sin(P2)
    def val(qq):
        a0, a1 = s(P0, P1)
        b0, b1 = s(P0 + qq[0], P1 + qq[1])
        ss = a0 ** 2 + a1 ** 2 + S2 ** 2 + s3 ** 2
        tt = b0 ** 2 + b1 ** 2 + S2 ** 2 + s3 ** 2
        dot = a0 * b0 + a1 * b1 + S2 * S2 + s3 * s3
        tr = -4 * (2 * s3 * s3 - dot) + 4 * M * M
        return np.mean(np.cos(P3) ** 2 * tr / ((ss + M * M) * (tt + M * M)))
    return val(q) - val((0.0, 0.0))


qs = (0.15, 0.3)
iso = {}
for label, d in (("Zeitachse", (1, 0)), ("Raumachse", (0, 1)), ("Diagonale", (1 / np.sqrt(2), 1 / np.sqrt(2)))):
    iso[label] = [pi33((q * d[0], q * d[1])) / q ** 2 for q in qs]
    print(f"   isotrop, q entlang {label:9s}: Pi/q^2 = " + ", ".join(f"{v:+.6f}" for v in iso[label]) + f"  (q = {qs})")
spread = max(abs(iso["Diagonale"][0] - iso["Raumachse"][0]), abs(iso["Zeitachse"][0] - iso["Raumachse"][0]))
shrink = abs(iso["Diagonale"][0] - iso["Raumachse"][0]) / abs(iso["Diagonale"][1] - iso["Raumachse"][1])
print(f"   Richtungsunterschied bei q = 0,15: {spread:.2e}; Verhältnis zu q = 0,3: {shrink:.3f} (q^2-Skalierung: 0,25)")
check("Zeit gleich Raum isotrop", abs(iso["Zeitachse"][0] - iso["Raumachse"][0]) < 1e-12)
check("Diagonalabweichung ~ q^2", 0.2 < shrink < 0.3)

print("\n   Gitter mit ausgezeichneter Zeit (Zeitsprung nu):")
res = {}
for nu in (1.0, 0.9):
    t = pi33((0.15, 0.0), nu) / 0.15 ** 2
    x = pi33((0.0, 0.15), nu) / 0.15 ** 2
    res[nu] = (t, x)
    print(f"   nu = {nu}: Pi(Zeit) = {t:+.6f}, Pi(Raum) = {x:+.6f}, c_gamma^2 = Pi(Raum)/Pi(Zeit) = {x / t:.4f}"
          f" (Fermion: 1/nu^2 = {1 / nu ** 2:.4f})")
two = (res[1.0][1] + res[0.9][1]) / (res[1.0][0] + res[0.9][0])
print(f"   zwei Arten (nu = 1 und 0,9): c_gamma^2 = {two:.4f}, die Arten haben 1 und {1 / 0.81:.4f}")
trend = [pi33((0.0, 0.15), 0.9, M) / pi33((0.15, 0.0), 0.9, M) for M in (0.5, 0.3, 0.2)]
print("   eine Art, nu = 0,9, Masse M = 0,5 / 0,3 / 0,2: c_gamma^2 = " + ", ".join(f"{v:.4f}" for v in trend)
      + f" -> 1/nu^2 = {1 / 0.81:.4f}: Der logarithmische (universelle) Anteil folgt der Art, der Rest ist ein")
print("   nicht-universeller Gitteranteil von etwa 2 %, der selbst schon Lorentz bei Dimension 4 verletzt.")
check("eine Art: Photon folgt bis auf Gitteranteil", abs(trend[-1] - 1 / 0.81) < abs(trend[0] - 1 / 0.81) < 0.04)
check("zwei Arten: Fehlanpassung", abs(two - 1) > 0.05 and abs(two - 1 / 0.81) > 0.05)
print("   Mit ausgezeichneter Zeit bekommt jede Art ihre eigene Geschwindigkeit, und das induzierte Photon")
print("   liegt dazwischen: Lorentzverletzung bei Dimension 4, nicht unterdrückt, nur durch Feinabstimmung")
print("   zu beheben. Auf dem H6 Gitter verbietet die Vierteldrehung zwischen Zeit- und Raumachse jedes")
print("   nu != 1 für jede Art und in jeder Ordnung. Die Gleichwertigkeit der Achsen, die die Kritik als")
print("   Hindernis für Lorentz las, ist genau die Bedingung für Lorentz ohne Feinabstimmung.")

print("\n4. Reflexionspositivität: zeitliche Korrelation eines naiven Fermions (Raumimpuls fest, omega = 0,4)\n")
omega = 0.4
p0 = 2 * np.pi * (np.arange(4096) + 0.5) / 4096
Ep = np.arcsinh(omega)
gs, gv = [], []
for t in range(0, 7):
    g_scalar = np.mean(np.cos(p0 * t) * omega / (np.sin(p0) ** 2 + omega ** 2))
    g_vec = np.mean(np.sin(p0 * t) * np.sin(p0) / (np.sin(p0) ** 2 + omega ** 2))
    gs.append(g_scalar)
    gv.append(g_vec)
    print(f"   t = {t}: skalarer Teil {g_scalar:+.6f}, Vektorteil {g_vec:+.6f}, e^(-E t) = {np.exp(-Ep * t):.6f}")
print("   Der skalare Teil lebt nur auf geraden, der Vektorteil nur auf ungeraden Zeiten: Die Pole bei p0 = iE")
print("   und p0 = pi + iE geben die Transfereigenwerte e^(-E) und -e^(-E). T selbst ist nicht positiv,")
print("   T^2 schon; die OS-Rekonstruktion läuft über Doppelschritte, und der Zeitverdoppler ist eine eigene Art.")
ratios = [gs[t + 2] / gs[t] for t in (0, 2, 4)]
print(f"   G(t+2)/G(t) = " + ", ".join(f"{r:.6f}" for r in ratios) + f"; e^(-2E) = {np.exp(-2 * Ep):.6f}")
check("Transfer T^2 positiv", all(abs(gs[t]) < 1e-12 for t in (1, 3, 5)) and all(abs(gv[t]) < 1e-12 for t in (0, 2, 4, 6)) and all(abs(r - np.exp(-2 * Ep)) < 1e-6 for r in ratios))

print("\n5. Größenordnungen\n")
z2 = json.loads(Path(__file__).with_name("z2_ergebnis.json").read_text())
EP = 1.22e19
for ratio in (z2["a_over_lP_three_generations_no_doublers"], z2["a_over_lP_three_generations_naive"]):
    E_LV = EP / (ratio * np.sqrt(1 / 3))
    print(f"   a = {ratio:.2f} l_P: quadratische Lorentzverletzung bei E_LV ~ {E_LV:.1e} GeV; bei 100 TeV ist"
          f" delta v ~ {(1e5 / E_LV) ** 2:.1e}")
print("   Beobachtungsschranken für quadratische Verletzung liegen höchstens um 10^12 GeV, also weit darunter.")
alpha1_inv, b1 = 59.0, 41 / 10
mu_landau = 91.19 * np.exp(alpha1_inv * 2 * np.pi / b1)
print(f"   U(1)_Y (GUT-normiert, alpha1^-1(M_Z) ~ {alpha1_inv}, b1 = 41/10): Landau-Pol bei ~ {mu_landau:.1e} GeV,"
      f" weit über der Gitterskala ~ 1e18 GeV.")
check("Landau-Pol über Gitterskala", mu_landau > 1e20)
print("   Die nichtabelschen Faktoren sind asymptotisch frei. Weil das Gitter physisch ist, braucht U(1)_Y")
print("   keinen Limes a -> 0; die bekannte Trivialität von U(1) ist kein Hindernis.")

print("\nNicht gezeigt: nichtperturbative Universalität (Confinement, Higgs-Sektor) und die chirale")
print("Gitterformulierung. Naive Fermionen haben 16 Geschmäcker in Paaren entgegengesetzter Chiralität")
print("(Nielsen-Ninomiya). Bekannter Weg aus der Literatur, hier nicht gerechnet: Spiegelpartner symmetrisch")
print("massiv machen; das gelingt nur, wenn die Z16-Anomalie einer Generation verschwindet, also mit 16")
print("Weyl-Fermionen je Generation (mit nu_R).")
print("\nZusammenfassung")
for n, c in OK:
    print(f"   {'bestanden' if c else 'NICHT bestanden'}: {n}")
Path(__file__).with_name("z6_ergebnis.json").write_text(json.dumps({
    "checks": {n: c for n, c in OK}, "pi33_isotropic": iso, "anisotropic": {str(k): v for k, v in res.items()},
    "c_gamma2_two_species": two, "landau_pole_GeV": mu_landau}, indent=1))
assert all(c for _, c in OK)
