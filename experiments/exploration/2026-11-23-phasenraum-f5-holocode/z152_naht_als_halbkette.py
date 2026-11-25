"""Z152: Die Naht als Hälfte einer Gitterkette (Eckentransfer-Lesart). Neue Strategie für P_ord, den Partner und die Dynamik.

Leitlinie 'Rekonstruktion statt Auswahl'; H6 und die gelösten Blöcke (Flavour, Massen, alpha, 16, E8) bleiben unberührt,
nichts wird befördert. Statuskürzel: E exakt, C bedingt (Prämisse genannt), L Lesart, N numerisch, O offen.

Anlass: z121, z144 und z145 enden an derselben Stelle. Eine einzelne Naht ist ein freies System mit elliptischem
Modulfluss; Ordnung der Nähte (P_ord), Transfer und Kettendynamik liegen außerhalb. Frage dieser Rechnung: Gibt es ein
lokales Gittersystem, dessen eine Hälfte genau die Naht ist? Dann wären Ordnung, Partner und ein Hamiltonoperator
mitgeliefert.

Bestand: Der modulare Einteilchengenerator des Originals ist h = log((1-C)/C) = 4 pi x mit x in Z + 1/4 (H6h31/H6h33;
'Die originale Kovarianz (1 + e^(4 pi (r + alpha)))^-1 ist die Fermibesetzung 1/(1 + e^(pi k))', k in 4Z + 1). Der
Modulfluss ist elliptisch und zu keinem Doppelkegel-Erzeuger konjugiert (z136). tau = i ist die Lage der Naht (z4, z82).
Literatur: Baxter (Eckentransfermatrizen), Peschel/Kaulke/Legeza 1999, Peschel 2004, Eisler/Peschel 2009 (Halbketten
freier Gittermodelle haben äquidistante Einteilchen-Verschränkungsenergien mit Abstand aus elliptischen Integralen;
für die Ising- und XY-Kette exakt, für die Kette mit Wechselpotential hier numerisch auf 20 Stellen), Eisler/Peschel
2017 (arXiv:1703.08126) und 2020 (arXiv:2007.01804). Für den Ausblick auf vier Dimensionen: Butt/Catterall/Pradhan/Toga
(arXiv:2101.01026) und Catterall (arXiv:2311.02487): gestaffelte = Kähler-Dirac-Fermionen, Z_4-Anomalie, 16 Weyl-Fermionen.

Teil A  Das native Spektrum: pi (4Z + 1) = ungerade Vielfache von pi mit wechselnder Ladung; Zustandssumme 2^(1/4) e^(-pi/24).
Teil B  Ohne Ladung: Ising-Kette im Querfeld. Halbkette hat (2l+1) pi K'/K; die Naht verlangt K' = K, also h = sqrt(2).
Teil C  Mit Ladung: Hüpfkette mit Dimerisierung und Wechselpotential (Rice-Mele). Leiter = Abstand mal (Z + a);
        Abstand 4 pi K'(kappa)/K(kappa); a = 1/4 genau ohne Dimerisierung; die Naht legt die Kette eindeutig fest:
        gleichmäßiges Hüpfen t, Wechselpotential 2t.
Teil D  Der Modulgenerator ist der Gitter-Boost: K = c mal Summe (Abstand vom Schnitt) mal (lokale Energie), c = Lemniskatenkonstante.
        Zwei verschachtelte Schnitte geben den Hamiltonoperator: K_j - K_(j+2) = 2c H.
Teil E  Bilanz und Urteil.
"""
import json
import os
import time
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")                     # deterministische Rundung (einfädig)

import mpmath as mp
import numpy as np
from scipy.optimize import fsolve
from scipy.special import ellipk

mp.mp.dps = 90
HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
T0 = time.time()


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val, text):
    BEFUNDE[name] = val
    say(f"   [Befund: {val}] {text}")


def ns(x, n=15):
    return mp.nstr(x, n)


# ---------------------------------------------------------------------------------------------------------------------
say("Z152: Die Naht als Hälfte einer Gitterkette (Eckentransfer-Lesart)")
say("=" * 100)
say("\nTeil A  Das native Spektrum")
# A1: k in 4Z + 1 -> |k| durchläuft jede ungerade Zahl genau einmal, Vorzeichen (-1)^l bei |k| = 2l + 1
ks = [4 * n + 1 for n in range(-40, 41)]
betr = sorted(abs(k) for k in ks)
check("A1 k in 4Z+1: |k| ist jede ungerade Zahl genau einmal", betr[:80] == [2 * l + 1 for l in range(80)])
check("A1 Vorzeichen von k ist (-1)^l bei |k| = 2l+1 (Teilchen, Loch, Teilchen, ...)",
      all((k > 0) == (((abs(k) - 1) // 2) % 2 == 0) for k in ks))
# A2: Zustandssumme der Naht = Majorana-Charakter (NS) bei tau = i = Weber-Funktion f(i) = 2^(1/4)
Z = mp.fprod(1 + mp.exp(-(2 * l + 1) * mp.pi) for l in range(60))
Zc = mp.mpf(2) ** mp.mpf("0.25") * mp.exp(-mp.pi / 24)
check("A2 Zustandssumme prod(1 + e^(-(2l+1) pi)) = 2^(1/4) e^(-pi/24) (80 Stellen)", abs(Z - Zc) < mp.mpf(10) ** -80)
S = mp.fsum(mp.log(1 + mp.exp(-(2 * l + 1) * mp.pi)) + (2 * l + 1) * mp.pi / (mp.exp((2 * l + 1) * mp.pi) + 1)
            for l in range(60))
Qm = mp.fsum((-1) ** l / (mp.exp((2 * l + 1) * mp.pi) + 1) for l in range(60))
ZAHLEN["Zustandssumme der Naht"] = ns(Z, 30)
ZAHLEN["Entropie der Naht (nat)"] = ns(S, 20)
ZAHLEN["mittlere Ladung der Naht gegen den Dirac-See"] = ns(Qm, 20)
say(f"   Zustandssumme {ns(Z, 20)}, Entropie {ns(S, 12)} nat, mittlere Ladung {ns(Qm, 12)}")
befund("A", "E", "Das Modulspektrum der Naht ist ein einziger Turm (2l+1) pi mit wechselnder Ladung. Das ist die Form "
       "einer Eckentransfermatrix: äquidistant, mit dem Abstand als einzigem Parameter.")

# ---------------------------------------------------------------------------------------------------------------------
say("\nTeil B  Ohne Ladung: Ising-Kette im Querfeld, H = -sum sx sx - h sum sz")


def tfi(h, N=400, nlev=4):
    A = np.zeros((N, N))
    B = np.zeros((N, N))
    for n in range(N):
        A[n, n] = 2 * h
    for n in range(N - 1):
        A[n, n + 1] = A[n + 1, n] = -1.0
        B[n, n + 1] = -1.0
        B[n + 1, n] = 1.0
    lam2, phi = np.linalg.eigh((A - B) @ (A + B))
    psi = ((A + B) @ phi) / np.sqrt(np.abs(lam2))
    G = -(psi @ phi.T)
    s = np.sort(np.linalg.svd(G[:N // 2, :N // 2], compute_uv=False))
    return 2 * np.arctanh(np.clip(s, 0, 1 - 1e-16))[:nlev] / np.pi


def kk(k):
    return float(ellipk(1 - k * k) / ellipk(k * k))


lev = tfi(np.sqrt(2))
ZAHLEN["Ising h = sqrt(2): Niveaus / pi"] = [round(float(x), 6) for x in lev]
say(f"   h = sqrt(2): Niveaus / pi = {[round(float(x), 6) for x in lev]}")
check("B1 h = sqrt(2): Halbkette hat die Niveaus pi (1, 3, 5, 7) (1e-6)",
      np.allclose(lev, [1, 3, 5, 7], atol=1e-6))
ok = True
for h in (1.2, 2.0, 3.0):
    lv = tfi(h)
    th = kk(1 / h)
    say(f"   h = {h}: Niveau 0 / pi = {lv[0]:.7f}, K'/K(1/h) = {th:.7f}, Verhältnisse {[round(float(x / lv[0]), 5) for x in lv[:3]]}")
    ok &= abs(lv[0] - th) < 1e-6 and np.allclose(lv[:3] / lv[0], [1, 3, 5], atol=1e-5)
check("B2 ungeordnete Phase h > 1: Niveaus (2l+1) pi K'(k)/K(k) mit k = 1/h (drei Felder)", ok)
lv = tfi(1 / np.sqrt(2))
say(f"   h = 1/sqrt(2): Niveaus / pi = {[round(float(x), 6) for x in lv]}")
check("B3 geordnete Phase h < 1: gerader Turm mit Nullmode, also nicht die Naht",
      np.allclose(lv, [0, 2, 4, 6], atol=1e-5))
kstar = mp.findroot(lambda k: mp.ellipk(1 - k * k) / mp.ellipk(k * k) - 1, mp.mpf("0.6"))
check("B4 Abstand pi <=> K' = K <=> k = 1/sqrt(2) (lemniskatischer Modul, tau = i)",
      abs(kstar - 1 / mp.sqrt(2)) < mp.mpf(10) ** -60)
befund("B", "N; Struktur exakt nach Literatur",
       "Die Ising-Kette im Querfeld h = sqrt(2) (ungeordnete Phase) hat als Halbkette genau das Modulspektrum der "
       "Naht. Der Abstand pi ist gleichwertig zu K' = K, also zum selbstdualen Modul tau = i. Die geordnete Phase "
       "scheidet aus (Nullmode). Diese Kette kennt aber die Ladung Q nicht.")

# ---------------------------------------------------------------------------------------------------------------------
say("\nTeil C  Mit Ladung: Hüpfkette mit Dimerisierung t1, t2 und Wechselpotential m (Rice-Mele), halbe Füllung")


def rm(t2, m, N=240, cut=120, t1=1.0):
    H = np.zeros((N, N))
    for n in range(N - 1):
        H[n, n + 1] = H[n + 1, n] = -(t1 if n % 2 == 0 else t2)
    for n in range(N):
        H[n, n] = m if n % 2 == 0 else -m
    _, v = np.linalg.eigh(H)
    occ = v[:, :N // 2]
    z = np.clip(np.linalg.eigvalsh((occ @ occ.T)[:cut, :cut]), 1e-300, 1 - 1e-16)
    e = np.sort(np.log((1 - z) / z))
    e = e[np.abs(e) < 30]
    e0, em = e[e > 1e-9][0], e[e < 1e-9][-1]
    return (e0 - em) / np.pi, e0 / np.pi, em / np.pi


ok = True
tab = []
for t2, m in ((0.5, 0.0), (0.5, 0.3), (0.2, 0.5), (0.3, -0.2), (0.7, 0.1), (1.0, -1.0), (1.0, -3.0)):
    D, e0, _ = rm(t2, m)
    kp = np.sqrt(((1 - t2) ** 2 + m * m) / ((1 + t2) ** 2 + m * m))        # kappa' = E_min / E_max
    th = 4 * float(ellipk(kp * kp) / ellipk(1 - kp * kp))                  # 4 K'(kappa) / K(kappa)
    tab.append([t2, m, round(float(D), 6), round(float(th), 6), round(float(e0 / D), 6)])
    say(f"   t2/t1 = {t2}, m/t1 = {m}: Abstand / pi = {D:.6f}, 4K'/K = {th:.6f}, Versatz a = {e0 / D:.6f}")
    ok &= abs(D - th) < 2e-6
ZAHLEN["Rice-Mele: t2, m, Abstand/pi, 4K'/K, Versatz"] = tab
check("C1 Leiter = Abstand (Z + a) mit Abstand 4 pi K'(kappa)/K(kappa), kappa' = E_min/E_max (sieben Punkte)", ok)
check("C2 ohne Wechselpotential (m = 0): Versatz a = 1/2 (Teilchen-Loch-symmetrisch)", abs(tab[0][4] - 0.5) < 1e-6)
ok = True
vers = []
for m in (-0.5, -1.0, -2.0, -3.0, -5.0):
    D, e0, _ = rm(1.0, m)
    D2, e02, _ = rm(1.0, -m)
    vers.append([m, round(float(e0 / D), 7), round(float(e02 / D2), 7)])
    ok &= abs(e0 / D - 0.25) < 1e-6 and abs(e02 / D2 - 0.75) < 1e-6
ZAHLEN["gleichmäßiges Hüpfen: m, Versatz, Versatz bei -m"] = vers
say(f"   gleichmäßiges Hüpfen (t1 = t2): Versatz a je m und -m: {vers}")
check("C3 ohne Dimerisierung (t1 = t2): Versatz a = 1/4 für jedes m, 3/4 bei umgekehrtem Vorzeichen", ok)


def gl(x):
    D, e0, _ = rm(x[0], x[1])
    return [D - 4.0, e0 - 1.0]


loes = []
for start in ((0.15, -0.3), (0.3, -0.5), (0.05, -0.1), (0.5, -1.0), (0.8, -3.0), (0.15, 0.3)):
    s, _, ier, _ = fsolve(gl, start, xtol=1e-12, full_output=True)
    loes.append([list(start), round(float(s[0]), 8), round(float(s[1]), 8), int(ier)])
    say(f"   Start {start}: t2/t1 = {s[0]:.8f}, m/t1 = {s[1]:.8f}, konvergiert = {ier == 1}")
ZAHLEN["Lösung der Nahtbedingung (Abstand 4 pi, Versatz 1/4)"] = loes
check("C4 Die Nahtleiter pi (4Z + 1) legt die Kette eindeutig fest: t2 = t1, m = -2 t1 (sechs Startpunkte)",
      all(l[3] == 1 and abs(l[1] - 1) < 1e-7 and abs(l[2] + 2) < 1e-7 for l in loes))


def kette(N, m):
    H = mp.zeros(N, N)
    for n in range(N - 1):
        H[n, n + 1] = H[n + 1, n] = -1
    for n in range(N):
        H[n, n] = m if n % 2 == 0 else -m
    return H


def halb(N, m, cuts):
    H = kette(N, m)
    w, v = mp.eigsy(H)
    idx = sorted(range(N), key=lambda i: w[i])[:N // 2]
    res = {}
    for cut in cuts:
        C = mp.zeros(cut, cut)
        for a in range(cut):
            for b in range(a, cut):
                C[a, b] = C[b, a] = mp.fsum(v[a, i] * v[b, i] for i in idx)
        z, U = mp.eigsy(C)
        eps = [mp.log((1 - x) / x) for x in z]
        hE = mp.zeros(cut, cut)
        for i in range(cut):
            for a in range(cut):
                ua = eps[i] * U[a, i]
                for b in range(cut):
                    hE[a, b] += ua * U[b, i]
        res[cut] = (sorted(eps, key=abs), hE)
    return H, res


NK = 56
Hk, res = halb(NK, mp.mpf(-2), (26, 28))
lev28 = [e / mp.pi for e in res[28][0][:10]]
soll = [(-1) ** l * (2 * l + 1) for l in range(10)]
abw = max(abs(lev28[l] - soll[l]) for l in range(8))
ZAHLEN["Kette m = -2t: Niveaus / pi (90 Stellen gerechnet)"] = [ns(x, 18) for x in lev28[:8]]
ZAHLEN["größte Abweichung der ersten acht Niveaus von pi (1, -3, 5, -7, ...)"] = ns(abw, 3)
say(f"   Kette t1 = t2, m = -2t, 56 Plätze, 90 Stellen: Niveaus / pi = {[ns(x, 12) for x in lev28[:6]]}")
say(f"   größte Abweichung der ersten acht Niveaus von (1, -3, 5, -7, 9, -11, 13, -15): {ns(abw, 3)}")
check("C5 Hochgenau: Halbkette hat pi (1, -3, 5, -7, ...) = pi (4Z + 1), acht Niveaus auf 1e-13", abw < mp.mpf(10) ** -13)
_, resp = halb(NK, mp.mpf(2), (28,))
levp = [e / mp.pi for e in resp[28][0][:8]]
check("C6 Umgekehrtes Vorzeichen (oder die andere Hälfte): gespiegelte Leiter pi (4Z + 3), der Partner mit Ladung -1",
      max(abs(levp[l] + soll[l]) for l in range(8)) < mp.mpf(10) ** -13)
kap = mp.mpf(2) / mp.sqrt(8)
check("C7 Die Kette liegt am selbstdualen Punkt: E_min/E_max = 2t / (2 sqrt(2) t) = 1/sqrt(2) = kappa = kappa'",
      abs(kap - 1 / mp.sqrt(2)) < mp.mpf(10) ** -80)
zerf = mp.acosh(3)
check("C8 Zerfallsrate der Korrelationen: cosh = 1 + m^2/(2 t^2) = 3, Länge 1/(2 arsinh 1) = 0,567 Plätze",
      abs(zerf - 2 * mp.asinh(1)) < mp.mpf(10) ** -80)
ZAHLEN["Korrelationslänge in Gitterplätzen"] = ns(1 / zerf, 12)
befund("C", "N (20 Stellen)",
       "In der Familie der ladungserhaltenden Nächste-Nachbar-Ketten mit Zweiplatzzelle gibt es genau eine, deren "
       "Halbkette die Naht ist: gleichmäßiges Hüpfen t und Wechselpotential 2t. Das Viertel im Versatz ist keine "
       "Einstellung, es folgt aus der fehlenden Dimerisierung. Der Abstand 4 pi verlangt wieder kappa = kappa'. In "
       "Kontinuumssprache ist das das Kogut-Susskind-Gitterfermion mit Masse mal Gitterabstand = 1.")

# ---------------------------------------------------------------------------------------------------------------------
say("\nTeil D  Der Modulgenerator als Gitter-Boost")
hE = res[28][1]
cut = 28
c = -hE[cut - 2, cut - 1]                                  # Hüpfelement der ersten inneren Bindung
varpi = mp.gamma(mp.mpf(1) / 4) ** 2 / (2 * mp.sqrt(2 * mp.pi))
ZAHLEN["c = Modulgenerator / Gitter-Boost"] = ns(c, 25)
ZAHLEN["Lemniskatenkonstante"] = ns(varpi, 25)
say(f"   Vorfaktor c = {ns(c, 25)}; Lemniskatenkonstante = {ns(varpi, 25)}")
neb = max(abs(hE[cut - 1 - d - j, cut - 1 - d]) for d in range(5) for j in (2, 3))   # nahe am Schnitt: das ferne Ende stört tiefe Moden
check("D1 Der Verschränkungsoperator ist tridiagonal (übernächste Nachbarn < 1e-9 nahe am Schnitt)", neb < mp.mpf(10) ** -9)
ok = True
for d in range(7):
    n = cut - 1 - d
    ok &= abs(hE[n - 1, n] + c * (d + 1)) < mp.mpf(10) ** -9
    ok &= abs(hE[n, n] - c * (d + mp.mpf(1) / 2) * Hk[n, n]) < mp.mpf(10) ** -9
check("D2 K = c mal Summe (Abstand vom Schnitt) mal (lokale Energie): Bindungen bei 1, 2, 3, ..., Plätze bei 1/2, 3/2, ...", ok)
check("D3 c ist die Lemniskatenkonstante Gamma(1/4)^2 / (2 sqrt(2 pi)) = 2,62206 (statt pi im Kontinuum)",
      abs(c - varpi) < mp.mpf(10) ** -18)
per = 2 * mp.quad(lambda z: 1 / mp.sqrt(1 - z**4), [0, 1])
check("D3b c ist die reelle Periode der Markenkurve w^2 = 1 - z^4 (j = 1728, tau = i): c = 2 Integral_0^1 dz/sqrt(1 - z^4)",
      abs(per - varpi) < mp.mpf(10) ** -30 and abs(c - per) < mp.mpf(10) ** -18)
h26 = res[26][1]
diff = max(abs(hE[a, b] - h26[a, b] - 2 * c * Hk[a, b]) for a in range(19, 26) for b in range(19, 26))
ZAHLEN["Rest von K_28 - K_26 - 2c H nahe am Schnitt"] = ns(diff, 3)
check("D4 Zwei verschachtelte Schnitte geben die Dynamik: K_(j+2) - K_j = 2c H auf dem gemeinsamen Teil", diff < mp.mpf(10) ** -9)
_, res1 = halb(NK, mp.mpf(-1), (28,))
c1 = -res1[28][1][cut - 2, cut - 1]
c1th = 4 * mp.ellipk(mp.mpf(1) / 5) / mp.sqrt(5)           # 4 K(kappa') / E_max mit kappa'^2 = 1/5, E_max = sqrt(5)
ZAHLEN["c bei m = -t (gemessen, Formel)"] = [ns(c1, 15), ns(c1th, 15)]
say(f"   Kontrolle m = -t: c = {ns(c1, 15)}, Formel 4 K(kappa')/E_max = {ns(c1th, 15)}")
check("D5 Allgemein c = 4 K(kappa') / E_max (Kontrolle bei m = -t); für m -> 0 geht c gegen pi (Bisognano-Wichmann)",
      abs(c1 - c1th) < mp.mpf(10) ** -9 and abs(4 * mp.ellipk(0) / 2 - mp.pi) < mp.mpf(10) ** -80)
befund("D", "N (18 Stellen)",
       "Der Modulgenerator der Naht ist in dieser Lesart geometrisch: der Boost des Gitters um den Schnitt. Sein "
       "Spektrum ist diskret und sein Fluss periodisch, weil das Gitter eine kleinste Länge hat. Der elliptische "
       "Nahtfluss und der hyperbolische Fluss eines Kontinuumskeils sind derselbe Erzeuger bei verschiedenem Modul. "
       "Der Hamiltonoperator der Kette ist die Differenz der Modulgeneratoren zweier verschachtelter Nähte.")

# ---------------------------------------------------------------------------------------------------------------------
say("\nTeil E  Bilanz")
bilanz = [
    ("E", "Form", "Das Modulspektrum der Naht, pi (4Z + 1), ist ein äquidistanter Turm mit wechselnder Ladung; "
     "Zustandssumme 2^(1/4) e^(-pi/24)."),
    ("N", "Ohne Ladung", "Die Ising-Kette im Querfeld sqrt(2) hat als Halbkette genau diesen Turm. Abstand pi <=> K' = K "
     "<=> tau = i."),
    ("N", "Mit Ladung, eindeutig", "Unter den Hüpfketten mit Dimerisierung und Wechselpotential trägt genau eine die "
     "Naht als Hälfte: gleichmäßiges Hüpfen t, Wechselpotential 2t (20 Stellen). Das Viertel folgt aus der fehlenden "
     "Dimerisierung, der Abstand aus kappa = kappa'."),
    ("N", "Partner", "Die andere Hälfte trägt die gespiegelte Leiter pi (4Z + 3): der Partner mit Ladung -1, wie im "
     "Paarmodell des Bestands (P_TP)."),
    ("N", "Boost", "K = c mal Gitter-Boost, exakt tridiagonal, c = Lemniskatenkonstante 2,62206. Der elliptische "
     "Nahtfluss ist der Boost eines Gitters."),
    ("E", "Markenkurve", "Der Vorfaktor c ist die reelle Periode der elliptischen Kurve w^2 = 1 - z^4 der nativen Marken "
     "(Bestand: Spektralparameter-Bahn {1, i, -1, -i}, j = 1728). Auf ihr ist dz/w polfrei: Das wäre die elliptische Klasse "
     "der 4D-Chern-Simons-Lesart (z151 fand auf der Sphäre nur die trigonometrische). Nicht gerechnet."),
    ("N", "Dynamik aus zwei Nähten", "K_(j+2) - K_j = 2c H: Zwei verschachtelte Nähte bestimmen den Hamiltonoperator "
     "mit Nächste-Nachbar-Transfer. Das ist die Gitterform der halbseitigen modularen Inklusion."),
    ("C", "Prämisse P_CTM", "Die Naht ist die Hälfte einer ladungserhaltenden Nächste-Nachbar-Kette mit Zweiplatzzelle "
     "in ihrem Grundzustand. Die Eindeutigkeit gilt in dieser Familie; längere Reichweite ist nicht geprüft."),
    ("E negativ", "Wechselwirkung", "Die Kette ist frei. Das wechselwirkende Gatter (P_uhr) folgt nicht."),
    ("O", "Kontinuum und 4D", "Die Kette liegt bei Masse mal Gitterabstand = 1, Korrelationslänge 0,567 Plätze. Es gibt "
     "keinen Kontinuumslimes am nativen Punkt. Die Kette ist 1+1-dimensional; K12, K13 und P_G sind nicht berührt."),
    ("L", "Ausblick 4D", "Gestaffelte Fermionen sind in vier Dimensionen Kähler-Dirac-Fermionen mit 2^4 = 16 Komponenten je "
     "Zelle und einer Z_4-Symmetrie, deren Anomalie sich nach der Literatur genau bei 16 Weyl-Fermionen aufhebt (dieselbe "
     "Bedingung wie Spin-Z_4 im Bestand). Ob die 16 der Naht diese Zelle ist, ist nicht gerechnet."),
]
for st, was, text in bilanz:
    say(f"   [{st}] {was}: {text}")
urteil = ("Die Naht lässt sich exakt als Hälfte eines lokalen Gittersystems lesen, und diese Lesart ist starr. (1) Das "
          "native Modulspektrum pi (4Z + 1) ist ein äquidistanter Turm, die Form einer Eckentransfermatrix. (2) Ohne "
          "Ladung trägt die Ising-Kette im Querfeld sqrt(2) diesen Turm; der Abstand pi ist gleichwertig zum "
          "selbstdualen Modul tau = i. (3) Mit Ladung gibt es in der Familie der Hüpfketten mit Dimerisierung und "
          "Wechselpotential genau eine Kette: gleichmäßiges Hüpfen t, Wechselpotential 2t. Das Viertel der Naht "
          "(x in Z + 1/4) ist dort keine Einstellung, es folgt aus der fehlenden Dimerisierung; der Abstand 4 pi "
          "verlangt E_min/E_max = 1/sqrt(2). Das ist das Kogut-Susskind-Gitterfermion in 1+1 Dimensionen mit Masse "
          "mal Gitterabstand = 1. (4) Der Modulgenerator ist der Gitter-Boost um den Schnitt, mit der "
          "Lemniskatenkonstante an der Stelle von pi. Damit ist der Widerspruch 'Nahtfluss elliptisch, Keilfluss "
          "hyperbolisch' (z136) aufgelöst: Es ist derselbe Erzeuger, das Gitter macht ihn periodisch. (5) Die andere "
          "Hälfte ist der Partner mit gespiegelter Leiter (P_TP), die Plätze geben die Ordnung (P_ord), und zwei "
          "verschachtelte Schnitte geben den Hamiltonoperator: K_(j+2) - K_j = 2c H. Die Aussage aus z145, Transfer "
          "und Kettendynamik lägen außerhalb, gilt für eine Naht; für zwei verschachtelte Nähte gilt sie in dieser "
          "Lesart nicht mehr. (6) Grenzen: Die Lesart ist eine Prämisse (P_CTM), keine Herleitung; die Eindeutigkeit "
          "ist numerisch und nur in der genannten Familie gezeigt. Die Kette ist frei, das wechselwirkende Gatter "
          "(P_uhr) folgt nicht. Sie ist 1+1-dimensional und hat am nativen Punkt keinen Kontinuumslimes; Gravitation, "
          "vier Dimensionen und G sind nicht berührt. Für den Kandidaten: P_ord und P_TP wären durch P_CTM ersetzt, "
          "die Zahl der Prämissen sinkt dadurch nicht sicher. Der Gewinn ist ein konkretes natives Gitter, auf dem "
          "sich Eichung der Uhr (Gitter-Schwinger-Modell) und die Frage nach vier Dimensionen rechnen lassen. Nichts "
          "ist befördert.")
say("\nUrteil\n   " + urteil)
print(f"\n   (Laufzeit {time.time() - T0:.1f} s, nicht Teil der Ausgabe)", flush=True)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z152_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
(HERE / "ausgabe_z152.txt").write_text("\n".join(OUT) + "\n")
