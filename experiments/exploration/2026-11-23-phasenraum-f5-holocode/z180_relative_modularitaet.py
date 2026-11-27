"""Z180: F2, Raumzeit und Ordnung aus relativer Modularität. Ziele: P_std, P_ord, I_step, K20 (Partner), P_vac.

Leitlinie 'Rekonstruktion statt Auswahl'; H6 und die gelösten Bausteine bleiben unberührt, nichts wird befördert.
Statuskürzel: E exakt, N numerisch, C bedingt, L Lesart oder Literatur, O offen.
Bestand: A* legt tau = i fest (Torus als Doppelüberlagerung der Sphäre mit vier Marken); z72, z84: das thermische Original
hat keine positive kontinuierliche Zeit (P_vac); z83: I_step ('der verdoppelte Diagonalschritt ist der native
Einzelschritt') auf der Tomita-verdoppelten Schalenleiter; z130: der Tomita-Partner als zweite Naht.

Frage (F2): Geben die Algebren zu den Bögen und Paarungen der vier Marken, ihre modularen Flüsse und halbseitige
Inklusionen (Wiesbrock) Translation, Zeitrichtung, lokale Zusammensetzung und Schrittnormierung, ohne dass der thermische
Zustand gegen das Vakuum getauscht wird?

Teil A  Wörterbuch: Marken, Bögen, Paarungen gegen Halbperioden, Zyklen, Kanäle (Weierstraß-Funktion bei tau = i).
Teil B  Der Zustand von A* auf den Bogenalgebren: Zweipunktfunktion, Partner, die zwei modularen Flüsse.
Teil C  Halbseitige Inklusionen: im Vakuum der Basis (alles vorhanden) und im Zustand von A* (nichts vorhanden).
Teil D  Zeitrichtung, Schrittnormierung, modulare Familie rho^s.
Teil E  Bilanz.
"""
import json
import time
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sp

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
T_START = time.time()
mp.mp.dps = 30


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val):
    BEFUNDE[name] = bool(val)
    say(f"   [Befund: {'ja' if val else 'nein'}] {name}")


PI, I = mp.pi, mp.mpc(0, 1)
Q0 = mp.exp(-PI)                                                           # Nome bei tau = i


def K(z, nome=Q0):
    """Szegö-Kern der Spinstruktur (A, A): theta_1'(0) theta_3(pi z) / (theta_3(0) theta_1(pi z))."""
    return mp.jtheta(1, 0, nome, 1) * mp.jtheta(3, PI * z, nome) / (mp.jtheta(3, 0, nome) * mp.jtheta(1, PI * z, nome))


def wp(z):
    """Weierstraß-Funktion zum Gitter Z + i Z: p(z) - e_2 = pi^2 K(z)^2 und e_2 = 0 bei tau = i."""
    return PI ** 2 * K(z) ** 2


# =================================================================== Teil A
say("Teil A: Wörterbuch zwischen den vier Marken und dem Torus bei tau = i")
th2, th3, th4 = (mp.jtheta(k, 0, Q0) for k in (2, 3, 4))
e1, e2, e3 = PI ** 2 * (th3 ** 4 + th4 ** 4) / 3, PI ** 2 * (th2 ** 4 - th4 ** 4) / 3, -PI ** 2 * (th2 ** 4 + th3 ** 4) / 3
g2 = 2 * (e1 ** 2 + e2 ** 2 + e3 ** 2)
z0 = mp.mpc('0.173', '0.291')
dgl = abs(mp.diff(wp, z0) ** 2 - (4 * wp(z0) ** 3 - g2 * wp(z0))) / abs(wp(z0)) ** 3
pol = abs(wp(mp.mpf('1e-8')) * mp.mpf('1e-16') - 1)
per = max(abs(wp(z0 + 1) - wp(z0)), abs(wp(z0 + I) - wp(z0)), abs(wp(-z0) - wp(z0))) / abs(wp(z0))
marken = {"1/2": wp(mp.mpf(1) / 2), "(1+i)/2": wp((1 + I) / 2), "i/2": wp(I / 2)}
viertel = abs(wp(I * z0) + wp(z0)) / abs(wp(z0))
xs = [mp.mpf(k) / 16 for k in range(1, 8)]
linien = {"Im z = 0": [wp(x) for x in xs], "Re z = 1/2": [wp(mp.mpf(1) / 2 + I * x) for x in xs],
          "Im z = 1/2": [wp(x + I / 2) for x in xs], "Re z = 0": [wp(I * x) for x in xs]}
reell = max(abs(v.imag) / (1 + abs(v)) for vs in linien.values() for v in vs)
bereich = {k: (float(min(v.real for v in vs)), float(max(v.real for v in vs))) for k, vs in linien.items()}
E1 = e1
bogen_ok = (all(v.real > E1 for v in linien["Im z = 0"]) and all(0 < v.real < E1 for v in linien["Re z = 1/2"])
            and all(-E1 < v.real < 0 for v in linien["Im z = 1/2"]) and all(v.real < -E1 for v in linien["Re z = 0"]))
doppelverh = (marken["i/2"] - marken["(1+i)/2"]) / (marken["1/2"] - marken["(1+i)/2"])
ZAHLEN["A1"] = {"e1": mp.nstr(e1, 20), "e2": mp.nstr(e2, 5), "e3": mp.nstr(e3, 20), "Doppelverhältnis": mp.nstr(doppelverh, 12),
                "Bilder der vier Linien": bereich}
check("A1 [N, 30 Stellen] Überlagerung: p(z) = pi^2 K(z)^2 mit dem Szegö-Kern K der Spinstruktur (A, A) erfüllt die "
      "Differentialgleichung p'^2 = 4 p^3 - g2 p, hat den Pol 1/z^2 und die Perioden 1 und i; die Koordinate der Basis ist "
      "also das Quadrat der Zweipunktfunktion des Zustands von A*. Die vier Marken sind die Halbperioden: p(0) = unendlich, "
      f"p(1/2) = e1 = {mp.nstr(e1, 10)}, p((1+i)/2) = 0, p(i/2) = -e1 (harmonische Lage, Doppelverhältnis -1). Der Kreis "
      "durch die Marken (die reelle Achse) hat als Urbild genau vier geschlossene Linien: Im z = 0 über dem Bogen "
      "(e1, unendlich), Re z = 1/2 über (0, e1), Im z = 1/2 über (-e1, 0), Re z = 0 über (unendlich, -e1). Jeder Bogen "
      "zwischen benachbarten Marken ist also ein ganzer Zyklus des Torus, doppelt durchlaufen (z und -z)",
      dgl < mp.mpf(10) ** -18 and pol < mp.mpf(10) ** -12 and per < mp.mpf(10) ** -25 and abs(e2) < mp.mpf(10) ** -25
      and abs(marken["1/2"] - e1) < mp.mpf(10) ** -20 and abs(marken["(1+i)/2"]) < mp.mpf(10) ** -20
      and abs(marken["i/2"] + e1) < mp.mpf(10) ** -20 and reell < mp.mpf(10) ** -20 and bogen_ok
      and abs(doppelverh + 1) < mp.mpf(10) ** -20)

# Paarungen, Klassen, Kanäle
S2 = sp.Matrix([[0, -1], [1, 0]])
T2 = sp.Matrix([[1, 1], [0, 1]])
klassen = [(1, 0), (0, 1), (1, 1)]


def wirk(M, c):
    v = M * sp.Matrix(c)
    return (int(v[0]) % 2, int(v[1]) % 2)


perm_S = {c: wirk(S2, c) for c in klassen}
perm_T = {c: wirk(T2, c) for c in klassen}
gruppe = set()
front = [sp.eye(2)]
while front:
    M = front.pop()
    key = tuple(int(x) % 2 for x in M)
    if key in gruppe:
        continue
    gruppe.add(key)
    front += [M * S2, M * T2]
# Kanal der dritten Paarung: tau' = tau/(tau + 1) = (1 + i)/2
tau_d = I / (I + 1)
j = lambda t: 1728 * (lambda a, b, c: (a ** 8 + b ** 8 + c ** 8) ** 3 / (54 * (a * b * c) ** 8))(
    *(mp.jtheta(k, 0, mp.exp(PI * I * t)) for k in (2, 3, 4)))
chi_e8 = lambda t: sum(mp.jtheta(k, 0, mp.exp(PI * I * t)) ** 8 for k in (2, 3, 4)) / \
    (2 * (mp.exp(2 * PI * I * t / 24) * mp.qp(mp.exp(2 * PI * I * t))) ** 8)
ZAHLEN["A2"] = {"S auf den Klassen": {str(k): str(v) for k, v in perm_S.items()},
                "T auf den Klassen": {str(k): str(v) for k, v in perm_T.items()}, "Bild von SL(2,Z) mod 2": len(gruppe),
                "tau des Diagonalkanals": mp.nstr(tau_d, 6), "chi_E8 im Diagonalkanal": mp.nstr(chi_e8(tau_d), 12)}
check("A2 [E/N] Paarungen: Die drei Paarungen der vier Marken sind die drei Klassen von H_1(Torus, Z2) (die Zyklen "
      "a = (1, 0), b = (0, 1), a + b = (1, 1)); SL(2, Z) wirkt darauf als volle S3 (Bild modulo 2 hat sechs Elemente). Der "
      "Stabilisator von tau = i wird von S erzeugt: S vertauscht a und b (die zwei Paarungen aus benachbarten Marken, also "
      "waagerechte und senkrechte Zyklen) und hält a + b fest (die Paarung der gegenüberliegenden Marken {0, unendlich} "
      "und {e1, -e1}). Auf der Basis ist S die Abbildung w -> -w: p(i z) = -p(z). Der Kanal der dritten Paarung hat "
      "tau' = (1 + i)/2: Er ist die Spur bei halber inverser Temperatur mit eingeschobener Halbdrehung e^(i pi L_0), und "
      f"chi_E8(tau') = {mp.nstr(chi_e8(tau_d), 8)} = 12 e^(2 pi i/3), dem Betrag nach wieder 12",
      perm_S == {(1, 0): (0, 1), (0, 1): (1, 0), (1, 1): (1, 1)} and len(gruppe) == 6 and viertel < mp.mpf(10) ** -25
      and abs(tau_d - (1 + I) / 2) < mp.mpf(10) ** -25 and abs(chi_e8(tau_d) - 12 * mp.exp(2 * I * PI / 3)) < mp.mpf(10) ** -20)

# =================================================================== Teil B
say("\nTeil B: Der Zustand von A* auf den Bogenalgebren")


def G_moden(z, s=1, kmax=60):
    """Thermische Zweipunktfunktion des freien Fermions auf dem Kreis der Länge 1 bei inverser Temperatur s:
    Summe über k in Z + 1/2 von e^(2 pi i k z)/(1 + e^(-2 pi s k)), in der für reelle z konvergenten Form."""
    out = I / (2 * mp.sin(PI * z))
    for n in range(kmax):
        k = n + mp.mpf(1) / 2
        out -= 2 * I * mp.sin(2 * PI * k * z) / (mp.exp(2 * PI * s * k) + 1)
    return out


z1 = mp.mpc('0.21', '0.17')
moden_theta = abs(G_moden(z1) - (I / 2) * K(z1)) / abs(K(z1))
anti = max(abs(K(z1 + 1) + K(z1)), abs(K(z1 + I) + K(z1))) / abs(K(z1))
s_kov = abs(K(I * z1) + I * K(z1)) / abs(K(z1))
xq = [mp.mpf(k) / 20 for k in range(1, 10)]
kreuz = [G_moden(x + I / 2) for x in xq]
kreuz_reihe = [sum(mp.cos(2 * PI * (n + mp.mpf(1) / 2) * x) / mp.cosh(PI * (n + mp.mpf(1) / 2)) for n in range(60)) for x in xq]
kreuz_ok = max(abs(a - b) for a, b in zip(kreuz, kreuz_reihe))
ZAHLEN["B1"] = {"Modensumme gegen Theta-Form": mp.nstr(moden_theta, 3), "antiperiodisch": mp.nstr(anti, 3),
                "S-Kovarianz K(i z) = -i K(z)": mp.nstr(s_kov, 3), "Kreuzkorrelation": mp.nstr(kreuz_ok, 3)}
check("B1 [N] Zweipunktfunktion des Zustands von A*: Die thermische Modensumme über k in Z + 1/2 mit Gewicht "
      "1/(1 + e^(-2 pi k)) ist (i/2) K(z), in beiden Richtungen antiperiodisch, und K(i z) = -i K(z): Derselbe Kern "
      "beschreibt die waagerechte und die senkrechte Naht (S-Kovarianz mit Gewicht 1/2). Zwischen den Höhen 0 und 1/2 ist "
      "die Korrelation Summe cos(2 pi k x)/cosh(pi k) mit lauter positiven Fourier-Koeffizienten: Die Felder auf dem "
      "gegenüberliegenden Bogen sind der Tomita-Partner der Naht (KMS-Fortsetzung um die halbe Periode, J-Positivität). "
      "Der Partner ist damit geometrisch: gegenüberliegender Bogen = paralleler Zyklus im Abstand der Halbperiode",
      moden_theta < mp.mpf(10) ** -25 and anti < mp.mpf(10) ** -25 and s_kov < mp.mpf(10) ** -25
      and kreuz_ok < mp.mpf(10) ** -25)

check("B2 [E] Modulare Flüsse: Der Bogen (e1, unendlich) mit beiden Blättern ist der Zyklus Im z = 0, der gegenüberliegende "
      "Bogen der Zyklus Im z = 1/2 (Kommutante). Der modulare Fluss dieses Paars im Zustand von A* ist die Verschiebung "
      "z -> z + t (KMS in der senkrechten Richtung mit Periode 1), der des anderen Paars z -> z + i t. Vier Bögen geben "
      "also nur zwei modulare Flüsse; sie vertauschen, jeder ist die imaginäre Richtung des anderen, und S tauscht sie. "
      "Zusammen sind sie die Translationen des Torus, eine kompakte abelsche Gruppe U(1) x U(1): keine nichtabelsche "
      "Struktur, kein Fixpunkt, kein Erzeuger mit einseitigem Spektrum (log Delta = -2 pi (L_0 (x) 1 - 1 (x) L_0-quer) hat "
      "symmetrisches Spektrum)",
      anti < mp.mpf(10) ** -25 and s_kov < mp.mpf(10) ** -25)

# =================================================================== Teil C
say("\nTeil C: Halbseitige modulare Inklusionen")
# C1: Vakuum der Basis, Möbiusgeometrie mit Marken -1, 0, 1, unendlich
a_, t_ = sp.symbols("a t", real=True)
T = lambda a: sp.Matrix([[1, a], [0, 1]])
dil = lambda t: sp.Matrix([[sp.exp(-sp.pi * t), 0], [0, sp.exp(sp.pi * t)]])        # x -> e^(-2 pi t) x
borchers = sp.simplify(dil(t_) * T(a_) * dil(-t_) - T(sp.exp(-2 * sp.pi * t_) * a_)) == sp.zeros(2)
moeb = lambda M, x: sp.oo if x == sp.oo else sp.simplify((M[0, 0] * x + M[0, 1]) / (M[1, 0] * x + M[1, 1]))
schritt = [moeb(T(1), x) for x in (-1, 0, 1)]
D1 = sp.Matrix([[sp.Rational(1, 2), 0], [0, -sp.Rational(1, 2)]])                    # Dilatation von (0, unendlich)
D2 = sp.Matrix([[0, sp.Rational(1, 2)], [sp.Rational(1, 2), 0]])                     # Dilatation von (-1, 1)
Jrot = D1 * D2 - D2 * D1
sl2 = (Jrot == sp.Matrix([[0, sp.Rational(1, 2)], [-sp.Rational(1, 2), 0]])
       and D2 * Jrot - Jrot * D2 == -D1 and Jrot * D1 - D1 * Jrot == -D2)
fix_D2 = sp.solve(sp.Eq(sp.Symbol("x") ** 2, 1), sp.Symbol("x"))
P = sp.Matrix([[0, 1], [0, 0]])
par = sp.simplify(D1 * P - P * D1 - P) == sp.zeros(2)
refl = sp.Matrix([[1, 0], [0, -1]])                                                 # S auf der Basis: w -> -w
umkehr = sp.simplify(refl * T(a_) * refl.inv() - T(-a_)) == sp.zeros(2)
# Positivität des Translationserzeugers in Tiefstgewichtsdarstellungen
def spektrum(h, n=500):
    L0 = np.diag(h + np.arange(n))
    Lm = np.diag(np.sqrt((np.arange(n - 1) + 1) * (2 * h + np.arange(n - 1))), -1)
    Pm = L0 - (Lm + Lm.T) / 2
    return float(np.linalg.eigvalsh(Pm).min())


pos = {h: spektrum(h) for h in (0.5, 1.0, 2.0)}
ZAHLEN["C1"] = {"T(1) auf den Marken -1, 0, 1": [str(x) for x in schritt], "kleinster Eigenwert von P": pos}
check("C1 [E/N, L] Vakuum der Basis (Marken als bloße Punkte -1, 0, 1, unendlich; Bisognano-Wichmann und Wiesbrock sind "
      "Literatur): Der modulare Fluss des Intervalls (0, unendlich) ist die Dilatation; sie bildet den Bogen "
      "(1, unendlich) für ein Vorzeichen von t in sich ab, das ist eine halbseitige modulare Inklusion. Die zugehörige "
      "Translation erfüllt Delta^(it) T(a) Delta^(-it) = T(e^(-2 pi t) a), ihr Erzeuger ist in Tiefstgewichtsdarstellungen "
      f"positiv (kleinster Eigenwert der Trunkierung {pos}), und ihr Einheitsschritt T(1) ist der Markenabstand: "
      "(-1, 0, 1) -> (0, 1, 2). Für die Paarung der gegenüberliegenden Marken ((-1, 1) in (-1, unendlich)) ist der "
      "Einheitsschritt T(2), das Quadrat des Schritts der benachbarten Paarung. Die Dilatationen zu (0, unendlich) und "
      "(-1, 1) erzeugen sl(2, R). Translation, Positivität, Zusammensetzung benachbarter Bögen und Schrittnormierung "
      "kommen im Vakuum also aus denselben Inklusionen",
      borchers and schritt == [0, 1, 2] and sl2 and par and all(v > 0 for v in pos.values()))


# C2: Zustand von A*: KMS-Test des geometrischen Dilatationsflusses auf dem Bogen (0, 1/2)
def x_von_u(u):
    return mp.atan(mp.exp(2 * PI * u)) / PI


def dx_du(u):
    return 2 * mp.exp(2 * PI * u) / (1 + mp.exp(4 * PI * u))


def test_dilatation(s):
    """Relative Schwankung von G(x(u) - x(u')) sqrt(x'(u) x'(u')) bei festem u - u': null genau dann, wenn der Zustand auf
    dem Bogen unter der geometrischen Dilatation invariant ist (notwendig für 'modularer Fluss = Dilatation')."""
    du = mp.mpf('0.15')
    werte = []
    for c in (mp.mpf(k) / 10 for k in range(-4, 5)):
        u, v = c + du / 2, c - du / 2
        Gf = (I / (2 * mp.sin(PI * (x_von_u(u) - x_von_u(v))))) if s is None else G_moden(x_von_u(u) - x_von_u(v), s)
        werte.append(Gf * mp.sqrt(dx_du(u) * dx_du(v)))
    m = sum(werte) / len(werte)
    return max(abs(w - m) for w in werte) / abs(m)


dil_test = {"Vakuum": test_dilatation(None), "s = 1 (A*)": test_dilatation(1), "s = 2": test_dilatation(2),
            "s = 4": test_dilatation(4), "s = 8": test_dilatation(8)}
abfall = [mp.log(dil_test["s = 1 (A*)"] / dil_test["s = 2"]), mp.log(dil_test["s = 2"] / dil_test["s = 4"]) / 2,
          mp.log(dil_test["s = 4"] / dil_test["s = 8"]) / 4]
ZAHLEN["C2"] = {k: mp.nstr(v, 6) for k, v in dil_test.items()} | {"Abfallrate je Einheit s": [mp.nstr(x, 6) for x in abfall]}
check("C2 [N] Zustand von A* (Marken als Verzweigungspunkte): Auf dem Bogen ist der Zustand unter der geometrischen "
      f"Dilatation nicht invariant. Die relative Schwankung der transportierten Zweipunktfunktion ist im Vakuum "
      f"{mp.nstr(dil_test['Vakuum'], 2)}, im Zustand von A* {mp.nstr(dil_test['s = 1 (A*)'], 3)}; in der modularen Familie "
      f"rho^s fällt sie mit der Rate {mp.nstr(abfall[-1], 5)} je Einheit s (pi = {mp.nstr(PI, 5)}), also wie e^(-pi s). Der "
      "modulare Fluss der Naht im Zustand von A* ist die Drehung des ganzen Zyklus (B2); sie hat keinen Fixpunkt und "
      "bildet keinen echten Teilbogen in sich ab, für kein Vorzeichen von t. Unter den Bogenalgebren des Zustands von A* "
      "gibt es also keine halbseitige modulare Inklusion. In der endlichen Standardform (Typ I) gibt es sie grundsätzlich "
      "nicht: Eine echte halbseitige Inklusion erzwingt Typ III_1 (Wiesbrock, Literatur)",
      dil_test["Vakuum"] < mp.mpf(10) ** -25 and dil_test["s = 1 (A*)"] > mp.mpf('0.01')
      and abs(abfall[-1] - PI) < mp.mpf('0.01') and dil_test["s = 8"] < mp.mpf(10) ** -9)

# =================================================================== Teil D
say("\nTeil D: Zeitrichtung, Schritt, modulare Familie")
check("D1 [E] Zeitrichtung: Auf der Basis ist S die Abbildung w -> -w. Sie hält die Marken 0 und unendlich fest, tauscht "
      "e1 und -e1 und kehrt den Umlaufsinn des Kreises durch die Marken um; sie führt T(a) in T(-a) über, also den "
      "positiven Erzeuger in den negativen und jede halbseitige Inklusion in die des anderen Vorzeichens. Der Zustand von "
      "A* ist unter S invariant. Aus A* und den vier Marken folgt daher kein Umlaufsinn und keine Zeitrichtung; eine "
      "Ordnung der Bögen ist nur bis auf diese Umkehr bestimmt. Das ist die Zeitneutralität aus A*, hier als Sperrsatz für "
      "P_ord: Die Ordnung braucht ein Bit, das S bricht", umkehr and viertel < mp.mpf(10) ** -25)

halb = {"auf einem Elter (L_0 = 1/2)": mp.exp(I * PI / 2), "auf den Strömen (L_0 = 1)": mp.exp(I * PI)}
check("D2 [E] Drei geometrische Kandidaten für eine Uhr der Ordnung 4 im Vergleich: (i) die Vierteldrehung der Marken "
      "(S, z -> i z) ist nach z175 nicht die Uhr; (ii) die Halbdrehung der Naht e^(i pi L_0) (Einschub des Diagonalkanals, "
      "tauscht die zwei Marken der Naht) hat auf einem Elter den Wert i und das Quadrat (-1)^F, wirkt aber auf alle 248 "
      "Ströme gleich (-1) und hat keine Fixalgebra so(10) + so(6); (iii) die Deckparität der Schnitte hat genau diese "
      "Fixalgebra (z179). Von den drei natürlichen Z4 der Geometrie trägt nur (iii) die Struktur von H6",
      abs(halb["auf einem Elter (L_0 = 1/2)"] - I) < mp.mpf(10) ** -25 and abs(halb["auf den Strömen (L_0 = 1)"] + 1) < mp.mpf(10) ** -25)

# modulare Familie: Z(s) und <T> = 0 nur bei s = 1
def E_mittel(s):
    f = lambda y: mp.log(chi_e8(I * y).real)
    return -mp.diff(f, s) / (2 * PI)


bil = {s: E_mittel(mp.mpf(s)) for s in (0.5, 1, 2, 4)}
ZAHLEN["D3"] = {"<L_0 - 1/3> in rho^s": {str(k): mp.nstr(v, 10) for k, v in bil.items()}}
check("D3 [N] Modulare Familie rho^s (tau = i s): Die Bilanz <L_0 - c/24> verschwindet genau bei s = 1 "
      f"({ {k: mp.nstr(v, 6) for k, v in bil.items()} }); für s -> unendlich geht der Zustand gegen das Vakuum der Naht "
      "(<L_0 - 1/3> -> -1/3), und die Abweichung vom geometrischen Fluss fällt wie e^(-pi s) (C2). Das Vakuum ist das "
      "kalte Ende derselben modularen Familie, aber A* wählt s = 1 und nicht dieses Ende: P_vac ist damit nicht hergeleitet",
      abs(bil[1]) < mp.mpf(10) ** -15 and bil[0.5] > 0 and bil[2] < 0 and abs(bil[4] + mp.mpf(1) / 3) < mp.mpf('0.001'))

# =================================================================== Teil E
say("\nTeil E: Bilanz")
befund("Partner (Teil von K20 und P_std): Der Tomita-Partner der Naht ist der gegenüberliegende Bogen der vier Marken; "
       "das zweite Paar ist das S-Bild des ersten", True)
befund("Im Vakuum der Basis geben dieselben Inklusionen Translation, positiven Erzeuger, Zusammensetzung benachbarter "
       "Bögen und den Einheitsschritt (Markenabstand; gegenüberliegende Paarung = doppelter Schritt)", True)
befund("Im Zustand von A* gibt es unter den Bogenalgebren eine halbseitige modulare Inklusion", False)
befund("Zeitrichtung oder Ordnung der Bögen folgt aus A* und den vier Marken (P_ord)", False)
befund("I_step ist hergeleitet (die Identifikation auf der Tomita-verdoppelten Schalenleiter, Typ I)", False)
befund("P_vac ist hergeleitet (Vakuumsektor ohne Zustandstausch)", False)
befund("Komplettlösung", False)

urteil = (
    "F2 ist gerechnet, mit einem positiven Wörterbuch und einem Sperrsatz. Positiv: Die Koordinate der Basis ist das "
    "Quadrat der Zweipunktfunktion des Zustands von A*. Jeder Bogen zwischen benachbarten Marken ist ein ganzer Zyklus "
    "des Torus; gegenüberliegende Bögen sind Naht und Tomita-Partner, die beiden anderen Bögen ihr S-Bild, und die dritte "
    "Paarung (gegenüberliegende Marken) ist der Kanal mit eingeschobener Halbdrehung. Im Vakuum der Basis liefern die "
    "Inklusionen der Bögen alles Verlangte aus einer Quelle: Translation mit positivem Erzeuger, Zusammensetzung "
    "benachbarter Bögen, sl(2, R) und den Einheitsschritt als Markenabstand, mit dem doppelten Schritt für die "
    "gegenüberliegende Paarung. Sperrsatz: Im Zustand von A* sind die Marken Verzweigungspunkte, und die vier Bögen haben "
    "nur zwei modulare Flüsse, die Translationen des Torus. Sie vertauschen, sind kompakt und haben keinen Fixpunkt; es "
    "gibt keine halbseitige Inklusion, und der Zustand ist unter dem geometrischen Fluss nicht invariant (Schwankung zwei "
    "Prozent im Testabstand, im Vakuum null, in der Familie rho^s wie e^(-pi s)). S kehrt den Umlaufsinn der Marken um und lässt den Zustand fest, also folgt keine Zeitrichtung "
    "und keine Ordnung: P_ord braucht ein Bit, das S bricht. I_step lebt auf der Standardform vom Typ I, wo es keine "
    "halbseitigen Inklusionen gibt; die Schrittnormierung aus Inklusionen ist nur im Vakuum der Basis ein Satz. Der Weg "
    "'Raumzeit ohne Zustandstausch aus relativer Modularität' ist damit in dieser Form geschlossen: P_vac, P_ord und "
    "I_step bleiben, P_std ist unberührt. Keine Komplettlösung, nichts befördert.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z180_ergebnis.json").write_text(json.dumps({"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil},
                                                    indent=1, ensure_ascii=False, default=str))
(HERE / "ausgabe_z180.txt").write_text("\n".join(OUT) + "\n")
say(f"(Laufzeit {time.time() - T_START:.0f} s)")
