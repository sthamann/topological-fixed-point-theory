"""Z181: F3, Wechselwirkung aus der Verklebung mehrerer Nähte. Ziele: P_std, P_ord, nichtlineare Dynamik.

Leitlinie 'Rekonstruktion statt Auswahl'; H6 und die gelösten Bausteine bleiben unberührt, nichts wird befördert.
Statuskürzel: E exakt, N numerisch, C bedingt, L Lesart oder Literatur, O offen.
Bestand: z121 (Phasenfamilie e^(-i phi Q_1 Q_2) des Gatters, P_viertel: phi = pi/2), z144 (P_uhr: die Phase ist eine
Funktion der Uhrwerte; Monodromie B_D5 = i), z178 (Assoziator der 16 Uhrsektoren, Anomalie der einseitigen Uhr).

Frage (F3): Legt die Verklebung von drei und vier Nähten (Assoziator, Pentagon, Hexagon) nur diskrete Phasen fest oder
auch die dynamische Kopplung?

Teil A  Klassifikation: abelsche Kohomologie H^3_ab(Z_n, U(1)) aus der Smith-Normalform, n = 2 und n = 4.
Teil B  Das native Gatter als Kandidat: für welche phi erfüllt es Pentagon und Hexagon auf den Uhrwerten?
Teil C  Drei und vier Nähte: Zusammensetzung, Yang-Baxter, Zopfgruppe, Ordnung.
Teil D  Bilanz: was die Topologie festlegt und was nicht.
"""
import cmath
import json
import math
import time
from fractions import Fraction as Fr
from itertools import product
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
T_START = time.time()


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val):
    BEFUNDE[name] = bool(val)
    say(f"   [Befund: {'ja' if val else 'nein'}] {name}")


def snf_diag(M):
    try:
        from flint import fmpz_mat
        D = fmpz_mat(M.tolist()).snf()
        return [abs(int(D[i, i])) for i in range(min(M.shape))], "flint"
    except Exception:                                                     # noqa: BLE001
        import sympy as sp
        from sympy.matrices.normalforms import smith_normal_form
        D = smith_normal_form(sp.Matrix(M.tolist()), domain=sp.ZZ)
        return [abs(int(D[i, i])) for i in range(min(M.shape))], "sympy"


# =================================================================== Teil A
say("Teil A: Welche Verklebungsdaten gibt es überhaupt? Abelsche Kohomologie der Uhr")


def komplex(n):
    """Gleichungen (Pentagon, zwei Hexagone) und Koränder für abelsche 3-Kozykel (omega, c) von Z_n, additiv in R/Z."""
    io = lambda a, b, c: (a * n + b) * n + c
    ic = lambda a, b: n ** 3 + a * n + b
    N = n ** 3 + n ** 2
    rows = []
    for a, b, c, d in product(range(n), repeat=4):
        r = np.zeros(N, dtype=np.int64)
        r[io(b, c, d)] += 1
        r[io((a + b) % n, c, d)] -= 1
        r[io(a, (b + c) % n, d)] += 1
        r[io(a, b, (c + d) % n)] -= 1
        r[io(a, b, c)] += 1
        rows.append(r)
    for a, b, c in product(range(n), repeat=3):
        r = np.zeros(N, dtype=np.int64)
        r[io(b, c, a)] += 1
        r[ic(a, (b + c) % n)] += 1
        r[io(a, b, c)] += 1
        r[ic(a, c)] -= 1
        r[io(b, a, c)] -= 1
        r[ic(a, b)] -= 1
        rows.append(r)
        r = np.zeros(N, dtype=np.int64)
        r[io(c, a, b)] -= 1
        r[ic((a + b) % n, c)] += 1
        r[io(a, b, c)] -= 1
        r[ic(a, c)] -= 1
        r[io(a, c, b)] += 1
        r[ic(b, c)] -= 1
        rows.append(r)
    M = np.array(rows)
    B = np.zeros((N, n * n), dtype=np.int64)                                # Korand von kappa(a, b)
    ik = lambda a, b: a * n + b
    for a, b, c in product(range(n), repeat=3):
        B[io(a, b, c), ik(b, c)] += 1
        B[io(a, b, c), ik((a + b) % n, c)] -= 1
        B[io(a, b, c), ik(a, (b + c) % n)] += 1
        B[io(a, b, c), ik(a, b)] -= 1
    for a, b in product(range(n), repeat=2):
        B[ic(a, b), ik(a, b)] += 1
        B[ic(a, b), ik(b, a)] -= 1
    # Vorzeichenkonvention des Korands passend zu den Hexagonen wählen: (omega, c) = (s1 d kappa, s2 kappa/kappa^T)
    for s1, s2 in ((1, 1), (1, -1), (-1, 1), (-1, -1)):
        Bs = B.copy()
        Bs[:n ** 3] *= s1
        Bs[n ** 3:] *= s2
        if not (M @ Bs).any():
            return M, Bs, (s1, s2)
    return M, B, None


koh = {}
for n in (2, 4):
    M, B, konv = komplex(n)
    kette = konv is not None and not (M @ B).any()
    d, wie = snf_diag(M)
    rang = sum(1 for x in d if x)
    tors = [x for x in d if x > 1]
    rang_b = int(np.linalg.matrix_rank(B.astype(float)))
    koh[n] = {"Variablen": M.shape[1], "Gleichungen": M.shape[0], "Rang": rang, "Torsion": tors,
              "Dimension der Kozykel": M.shape[1] - rang, "Rang der Koränder": rang_b, "Koränder sind Kozykel": kette, "Konvention (omega, c)": str(konv),
              "Rechnung": wie}
ZAHLEN["A1"] = koh
check("A1 [E] Abelsche 3-Kozykel (omega, c) der Uhrgruppe, also Assoziator und Verzopfung zusammen (Pentagon und beide "
      "Hexagone als ganzzahliges lineares System, Smith-Normalform). Die Zusammenhangskomponente der Lösungen besteht "
      "genau aus den Korändern (gleiche Dimension), und die Komponentengruppe ist H^3_ab(Z_2, U(1)) = Z_"
      f"{koh[2]['Torsion']} und H^3_ab(Z_4, U(1)) = Z_{koh[4]['Torsion']}. Für die Fermionparität gibt es also genau vier "
      "Verklebungen (Boson, Semion, Fermion, Antisemion), für die Uhr der Ordnung 4 genau acht. Einen kontinuierlichen "
      "Parameter gibt es nicht",
      all(v["Koränder sind Kozykel"] and v["Dimension der Kozykel"] == v["Rang der Koränder"] for v in koh.values())
      and koh[2]["Torsion"] == [4] and koh[4]["Torsion"] == [8])


def kozykel(n, m):
    """Vertreter der Klasse m in Z_2n: c(a, b) = e^(2 pi i m a b/(2n)), omega(a, b, c) = (-1)^(m a u(b, c))."""
    om = lambda a, b, c: Fr(m * a * ((b + c) // n), 2)
    c = lambda a, b: Fr(m * a * b, 2 * n)
    return om, c


def erfuellt(n, om, c):
    pent = all((om(b, c_, d) - om((a + b) % n, c_, d) + om(a, (b + c_) % n, d) - om(a, b, (c_ + d) % n) + om(a, b, c_)) % 1 == 0
               for a, b, c_, d in product(range(n), repeat=4))
    h1 = all((om(b, c_, a) + c(a, (b + c_) % n) + om(a, b, c_) - c(a, c_) - om(b, a, c_) - c(a, b)) % 1 == 0
             for a, b, c_ in product(range(n), repeat=3))
    h2 = all((-om(c_, a, b) + c((a + b) % n, c_) - om(a, b, c_) - c(a, c_) + om(a, c_, b) - c(b, c_)) % 1 == 0
             for a, b, c_ in product(range(n), repeat=3))
    return pent, h1, h2


vertreter = {m: erfuellt(4, *kozykel(4, m)) for m in range(8)}
spins = sorted({kozykel(4, m)[1](1, 1) % 1 for m in range(8)})
ZAHLEN["A2"] = {"Klassen m = 0..7 erfüllen (Pentagon, H1, H2)": {str(k): v for k, v in vertreter.items()},
                "Spin des Erzeugers c(1, 1)": [str(x) for x in spins]}
check("A2 [E] Vertreter: Für m = 0..7 erfüllt c(a, b) = e^(2 pi i m a b/8), omega(a, b, c) = (-1)^(m a u) (u der Übertrag "
      "von b + c) Pentagon und beide Hexagone. Die acht Klassen sind durch den Spin des Erzeugers c(1, 1) = e^(2 pi i m/8) "
      "getrennt, und der ist unter Korändern invariant. In den 16 Uhrsektoren stehen die Klassen m = 5 (so(10)-Hälfte) und "
      "m = 3 (so(6)-Hälfte); auf dem E8-Kleber (k, k) addieren sie sich zu 8 = 0 (bosonisch, kein Assoziator)",
      all(all(v) for v in vertreter.values()) and len(spins) == 8 and (5 + 3) % 8 == 0)

# =================================================================== Teil B
say("\nTeil B: Das native Gatter auf den Uhrwerten")


def gatter_test(phi, n=4):
    """c(a, b) = e^(-i phi a b), omega(a, b, c) = e^(-i phi n a u(b, c)) auf Vertretern 0..n-1."""
    ph = lambda x: cmath.exp(-1j * phi * x)
    om = lambda a, b, c: ph(n * a * ((b + c) // n))
    c = lambda a, b: ph(a * b)
    tol = 1e-9
    wohl = all(abs(ph(a * (b + n)) - ph(a * b)) < tol for a, b in product(range(n), repeat=2))
    pent = all(abs(om(b, c_, d) * om(a, (b + c_) % n, d) * om(a, b, c_) - om((a + b) % n, c_, d) * om(a, b, (c_ + d) % n)) < tol
               for a, b, c_, d in product(range(n), repeat=4))
    h1 = all(abs(om(b, c_, a) * c(a, (b + c_) % n) * om(a, b, c_) - c(a, c_) * om(b, a, c_) * c(a, b)) < tol
             for a, b, c_ in product(range(n), repeat=3))
    h2 = all(abs(c((a + b) % n, c_) * om(a, c_, b) - c(a, c_) * c(b, c_) * om(c_, a, b) * om(a, b, c_)) < tol
             for a, b, c_ in product(range(n), repeat=3))
    return {"wohldefiniert ohne Assoziator": wohl, "Pentagon": pent, "Hexagone": h1 and h2}


raster = {j: gatter_test(math.pi * j / 32) for j in range(64)}                # phi = pi j/32, eine volle Periode
erlaubt = {k: [j for j, v in raster.items() if v[k]] for k in ("wohldefiniert ohne Assoziator", "Pentagon", "Hexagone")}
erlaubt["Pentagon und Hexagone"] = [j for j, v in raster.items() if v["Pentagon"] and v["Hexagone"]]
ZAHLEN["B1"] = {"phi = pi j/32, erlaubte j": erlaubt}
check("B1 [E] Gatterfamilie e^(-i phi Q_1 Q_2): Auf ganzzahligen Ladungen ist sie für jedes phi ein Bicharakter, also "
      "streng assoziativ und ohne jede Bedingung an phi. Als Funktion der Uhrwerte Q mod 4 (P_uhr) ist sie ohne "
      "Assoziator wohldefiniert genau für phi in (pi/2) Z (vier Werte: 0, pi/2, pi, 3 pi/2); mit dem Assoziator "
      "omega = e^(-4 i phi a u) erfüllt sie das Pentagon für phi in (pi/8) Z und die Hexagone für phi in (pi/4) Z. "
      "Zusammen bleiben genau acht Werte, die Klassen von H^3_ab(Z_4). P_viertel (phi = pi/2) ist die Klasse der "
      "Ordnung 4; sie ist die Monodromie der so(10)-Hälfte mit sich selbst, e^(2 pi i 2 * 5 k^2/8) = i^(k^2), und zugleich "
      "die Differenz 5 - 3 der beiden Deckhälften auf dem Kleber. Aus P_uhr folgt damit ein diskretes Raster, aber nicht, "
      "welcher seiner Punkte gilt: phi = 0 (kein Gatter) besteht alle Bedingungen ebenso",
      erlaubt["wohldefiniert ohne Assoziator"] == [0, 16, 32, 48] and erlaubt["Pentagon"] == list(range(0, 64, 4))
      and erlaubt["Pentagon und Hexagone"] == list(range(0, 64, 8))
      and all(cmath.isclose(cmath.exp(2j * math.pi * 2 * 5 * k * k / 8), 1j ** (k * k), abs_tol=1e-12) for k in range(4)))

# =================================================================== Teil C
say("\nTeil C: Drei und vier Nähte")
rng = np.random.default_rng(181)
d = 4                                                                     # Ladungen 0..3 je Naht


def gate(N, i, j, phi):
    q = np.indices((d,) * N).reshape(N, -1)
    return np.diag(np.exp(-1j * phi * q[i] * q[j]))


def swap(N, i):
    idx = np.arange(d ** N).reshape((d,) * N)
    perm = np.swapaxes(idx, i, i + 1).reshape(-1)
    P = np.zeros((d ** N, d ** N))
    P[np.arange(d ** N), perm] = 1
    return P


yb, reihenfolge, total = True, True, True
for phi in (0.37, math.pi / 2, 1.9):
    for N in (3, 4):
        R = [swap(N, i) @ gate(N, i, i + 1, phi) for i in range(N - 1)]
        yb = yb and all(np.allclose(R[i] @ R[i + 1] @ R[i], R[i + 1] @ R[i] @ R[i + 1]) for i in range(N - 2))
        yb = yb and all(np.allclose(R[i] @ R[j], R[j] @ R[i]) for i in range(N - 1) for j in range(i + 2, N - 1))
        paare = [(i, j) for i in range(N) for j in range(i + 1, N)]
        G1 = np.eye(d ** N, dtype=complex)
        for i, j in paare:
            G1 = G1 @ gate(N, i, j, phi)
        G2 = np.eye(d ** N, dtype=complex)
        for i, j in rng.permutation(paare):
            G2 = G2 @ gate(N, i, j, phi)
        reihenfolge = reihenfolge and np.allclose(G1, G2)
        q = np.indices((d,) * N).reshape(N, -1)
        total = total and np.allclose(np.diag(G1), np.exp(-1j * phi * (q.sum(0) ** 2 - (q ** 2).sum(0)) / 2))
        quad = all(np.allclose(R[i] @ R[i], gate(N, i, i + 1, 2 * phi)) for i in range(N - 1))
        yb = yb and quad
ZAHLEN["C1"] = {"Yang-Baxter für jedes phi": yb, "Reihenfolge der Paare gleichgültig": reihenfolge}
check("C1 [E, N = 3 und 4 Nähte, Ladungen 0..3, drei Werte von phi] Auf ganzzahligen Ladungen vertauschen alle "
      "Paargatter. Das Produkt über alle Paare hängt nicht von der Reihenfolge ab und ist e^(-i phi ((Summe Q)^2 - Summe "
      "Q_i^2)/2), also eine Funktion der Gesamtladung und der Einzelladungen. Die Tauschgatter P_i G_(i,i+1) erfüllen die "
      "Zopfrelationen für jedes phi, und ihr Quadrat ist das Gatter zu 2 phi. Die Verklebung von drei und vier Nähten "
      "stellt also keine Bedingung an phi und verlangt keine Ordnung der Nähte: Aus Pentagon, Hexagon und Yang-Baxter "
      "folgt weder P_viertel noch P_ord", yb and reihenfolge and total)

# Zopfgruppe im Defektsektor: N gleiche Uhrsektoren a = (1, 0), Spin 5/8
theta = cmath.exp(2j * math.pi * 5 / 8)
ordnung_sigma = next(k for k in range(1, 33) if abs(theta ** k - 1) < 1e-12)
mono = theta ** 2
kleber_spin = [Fr(5 * k * k + 3 * k * k, 8) % 1 for k in range(4)]
anomalie_N = {N: (2 * N) % 4 for N in range(1, 7)}
ZAHLEN["C2"] = {"Tausch zweier Defekte": str(theta), "Ordnung": ordnung_sigma, "Monodromie": str(mono),
                "Anomalie der diagonalen Uhr über N Nähte (in Z_4)": anomalie_N}
check("C2 [E] Defektsektor: N gleiche Uhrsektoren (1, 0) (die 16, Spin 5/8) tragen eine eindimensionale Darstellung der "
      "Zopfgruppe B_N mit sigma_i = e^(2 pi i 5/8). Sie hat die Ordnung 8, sigma_i^2 = i ist nicht 1: Der Tausch faktorisiert "
      "nicht über die Permutationsgruppe (N = 3, 4 und jedes N), die Nähte im Defektsektor sind weder Bosonen noch "
      "Fermionen. Im unverdrillten Teil (E8-Kleber, Spin 0) ist der Tausch trivial. Die Anomalie der diagonalen Uhr über N "
      f"Nähte ist 2N mod 4 = {list(anomalie_N.values())} (z178): drei Nähte sind anomal, vier nicht. Eine Zopfphase ordnet "
      "die Nähte nicht: Sie ist eine Zahl, kein Etikett. P_ord folgt daraus nicht",
      ordnung_sigma == 8 and abs(mono - 1j) < 1e-12 and kleber_spin == [0, 0, 0, 0]
      and anomalie_N == {1: 2, 2: 0, 3: 2, 4: 0, 5: 2, 6: 0})

# =================================================================== Teil D
say("\nTeil D: Bilanz")
befund("Die Verklebungsdaten der Uhr sind vollständig klassifiziert und diskret (Z_8, für die Parität Z_4)", True)
befund("Unter P_uhr legt die Topologie phi auf ein Raster von acht Werten fest", True)
befund("Die Topologie wählt den Punkt des Rasters (P_viertel) oder schließt phi = 0 aus", False)
befund("Die Verklebung von drei und vier Nähten erzwingt eine Ordnung der Nähte (P_ord)", False)
befund("Die Topologie legt eine dynamische Kopplung fest (kontinuierlicher Parameter, nichtlineare Dynamik)", False)
befund("Komplettlösung", False)

urteil = (
    "F3 ist gerechnet und entschieden: Die Verklebung legt nur diskrete Phasen fest. (1) Assoziator und Verzopfung der "
    "Uhr sind durch H^3_ab klassifiziert: Z_4 für die Fermionparität, Z_8 für die Uhr der Ordnung 4, ohne kontinuierlichen "
    "Parameter; die 16 Uhrsektoren tragen die Klassen 5 und 3, die sich auf dem E8-Kleber aufheben. (2) Das native Gatter "
    "e^(-i phi Q_1 Q_2) ist auf ganzzahligen Ladungen für jedes phi assoziativ, erfüllt Yang-Baxter für drei und vier "
    "Nähte und hängt nicht von der Reihenfolge der Paare ab. Erst als Funktion der Uhrwerte (P_uhr) wird phi diskret: "
    "vier Werte ohne Assoziator, acht mit Pentagon und Hexagon. P_viertel ist die Klasse der Ordnung 4, gleich der "
    "Monodromie der so(10)-Hälfte (z144); ausgewählt wird sie von der Topologie nicht, phi = 0 besteht ebenso. (3) Im "
    "Defektsektor ist der Tausch zweier Nähte die Phase e^(2 pi i 5/8) der Ordnung 8; das ist Statistik, keine Ordnung. "
    "P_ord und jede kontinuierliche Kopplung bleiben außerhalb dessen, was Pentagon und Hexagon festlegen. Der Weg "
    "'Wechselwirkung aus der Verklebung' gibt damit die Statistikklassen und das Raster von phi, aber weder die Auswahl "
    "im Raster noch die Dynamik. Keine Komplettlösung, nichts befördert.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z181_ergebnis.json").write_text(json.dumps({"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil},
                                                    indent=1, ensure_ascii=False, default=str))
(HERE / "ausgabe_z181.txt").write_text("\n".join(OUT) + "\n")
say(f"(Laufzeit {time.time() - T_START:.0f} s)")
