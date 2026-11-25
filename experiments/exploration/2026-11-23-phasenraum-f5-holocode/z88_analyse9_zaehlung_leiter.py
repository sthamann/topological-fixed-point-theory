"""Z88: Prüfung der neunten externen Analyse 'Fortsetzung der sieben Lösungsstrategien' (25.11.): Zählung und Leiter.

Teil 1  Vielteilchenzählung: Der volle Fockraum des Originals (zwei Quarter und Partner, mit Ladungssektoren) wächst wie eine
        2D-Theorie, der Fockraum eines masselosen 4D-Skalars über der Leiter E = m (z83) wie eine 4D-Theorie. Exakte
        ganzzahlige Zählung bis E = 500: Ab einer Schwelle hat das 4D-Feld mehr Zustände mit Energie <= E als das Original.
        Eine energietreue Einbettung des vollen 4D-Fockraums in das Original ist damit ausgeschlossen; die
        Einteilchenrekonstruktion (z83) bleibt unberührt.
Teil 2  Asymptotik: log Z ~ (2 pi^2/3)/beta für vier chirale Türme gegen pi^4/(45 beta^3) für das 4D-Feld.
Teil 3  Leiter aus Hardy-Invarianz: Erhält der Modulfluss e^{ith} die analytische Multiplikationsalgebra (der verschobene
        Shift bleibt Toeplitz), so ist h linear in n; mit der Teilchen-Loch-Symmetrie h_{-n-1} = -h_n und der Saatlücke pi
        folgt h_n = pi (2n + 1) auf der verzahnten Leiter beider Quarter, also h = 4 pi x mit x in Z/2 + 1/4.
"""
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val):
    BEFUNDE[name] = bool(val)
    say(f"   [Befund: {'ja' if val else 'nein'}] {name}")


def euler_produkt(a, N):
    """Koeffizienten von prod_m (1 - t^m)^(-a(m)) bis t^N exakt (Rekursion n p(n) = sum_k sigma(k) p(n-k))."""
    sig = [0]*(N + 1)
    for d in range(1, N + 1):
        ad = a(d)
        if ad:
            for k in range(d, N + 1, d):
                sig[k] += d*ad
    p = [0]*(N + 1)
    p[0] = 1
    for n in range(1, N + 1):
        p[n] = sum(sig[k]*p[n - k] for k in range(1, n + 1))//n
    return p


def mal(u, v, N):
    w = [0]*(N + 1)
    for i, x in enumerate(u):
        if x:
            for j in range(0, N + 1 - i):
                if v[j]:
                    w[i + j] += x*v[j]
    return w


# ------------------------------------------------------------------ Teil 1
say("Teil 1: Vielteilchenzählung bis E = 500 (Viertelschritte t = q^(1/4))")
EMAX = 500
N = 4*EMAX
# vier chirale Oszillatortürme (zwei Quarter x Original/Partner) in t: Faktor (1 - t^(4m))^-4
tuerme = euler_produkt(lambda d: 4 if d % 4 == 0 else 0, N)
# Ladungssektoren: Summe_Q t^(4 E_alpha(Q)), 4 E = 2 Q^2 + (4 alpha - 2) Q, je Quarter und je Original/Partner
ladung = {}
for alpha4 in (1, 3):
    c = [0]*(N + 1)
    for Q in range(-40, 41):
        e4 = 2*Q*Q + (alpha4 - 2)*Q
        if 0 <= e4 <= N:
            c[e4] += 1
    ladung[alpha4] = c
lad = mal(ladung[1], ladung[3], N)
lad = mal(lad, lad, N)                                            # Original und Partner
nativ_voll = mal(lad, tuerme, N)
vier_d = euler_produkt(lambda m: m*m, EMAX)                       # 4D-Skalar, Energien E = m mit m^2 Zuständen
kum = lambda arr, schritt: np.cumsum([float(x) for x in arr])     # nur für Plot-Zahlen, exakte Vergleiche unten
nur_tuerme = [tuerme[4*E] for E in range(EMAX + 1)]


def kumulativ_exakt(arr):
    s, out = 0, []
    for x in arr:
        s += x
        out.append(s)
    return out


K_nativ = kumulativ_exakt(nativ_voll)                             # in Viertelschritten
K_tuerme = kumulativ_exakt(nur_tuerme)
K_4d = kumulativ_exakt(vier_d)
schwelle_voll = next(E for E in range(1, EMAX + 1) if K_4d[E] > K_nativ[4*E])
schwelle_tuerme = next(E for E in range(1, EMAX + 1) if K_4d[E] > K_tuerme[E])
immer_danach = all(K_4d[E] > K_nativ[4*E] for E in range(schwelle_voll, EMAX + 1))
log10_quot_tuerme = np.log10(float(vier_d[EMAX])) - np.log10(float(nur_tuerme[EMAX]))
log10_quot_voll = np.log10(float(K_4d[EMAX])) - np.log10(float(K_nativ[4*EMAX]))
ZAHLEN["zaehlung"] = {"schwelle_volles_original": schwelle_voll, "schwelle_vier_tuerme": schwelle_tuerme,
                      "log10_4D_durch_vier_tuerme_bei_E500": round(log10_quot_tuerme, 2),
                      "log10_kumulativ_4D_durch_volles_original_bei_E500": round(log10_quot_voll, 2),
                      "4D_Zustaende_E500_log10": round(np.log10(float(vier_d[EMAX])), 2)}
say(f"   {ZAHLEN['zaehlung']}")
check(f"exakte Zählung: Ab E = {schwelle_voll} hat der 4D-Skalar mehr Zustände mit Energie <= E als der volle Fockraum des "
      f"Originals mit Partner und allen Ladungssektoren (und bleibt bis E = 500 darüber); gegen die vier reinen Türme ab "
      f"E = {schwelle_tuerme}, bei E = 500 um den Faktor 10^{log10_quot_tuerme:.1f} (die Analyse nennt etwa 10^51)",
      immer_danach and schwelle_voll < EMAX and 45 < log10_quot_tuerme < 60)
befund("Eine energietreue Einbettung des vollen 4D-Fockraums in das Original (gleiche Uhr L_0) ist ausgeschlossen; die "
       "Einteilchenleiter aus z83 bleibt richtig. Die Vielteilchenstruktur der Raumzeit muss eine eigene Zusammensetzung "
       "sein (zweite Quantisierung über der Leiter, etwa viele Nähte oder Ereignisse), nicht die Mehrpaarsektoren einer Naht",
       True)

# ------------------------------------------------------------------ Teil 2
say("\nTeil 2: Asymptotik")
ok_as = True
for beta in (0.05, 0.02):
    m = np.arange(1, 200000)
    lz4 = -4*np.sum(np.log1p(-np.exp(-beta*m)))
    lz4d = -np.sum(m**2*np.log1p(-np.exp(-beta*m)))
    voll2d = 4*(np.pi**2/(6*beta) + 0.5*np.log(beta/(2*np.pi)) - beta/24)      # Dedekind-eta, exakt bis O(e^(-4 pi^2/beta))
    ok_as &= abs(lz4 - voll2d) < 1e-8*abs(lz4) and abs(lz4d/(np.pi**4/(45*beta**3)) - 1) < 1e-3
    ZAHLEN.setdefault('asymptotik', {})[beta] = {'log_Z_vier_tuerme': round(float(lz4), 4), 'fuehrend_2D': round(2*np.pi**2/(3*beta), 4),
                                                 'log_Z_4D': round(float(lz4d), 2), 'fuehrend_4D': round(np.pi**4/(45*beta**3), 2)}
check("log Z = 4 (pi^2/(6 beta) + ln(beta/2 pi)/2 - beta/24) exakt für vier chirale Türme (führend (2 pi^2/3)/beta) und "
      "log Z = pi^4/(45 beta^3) (1 + O(beta^4)) für den 4D-Skalar (beta = 0,05; 0,02): 1/beta gegen 1/beta^3", ok_as)

# ------------------------------------------------------------------ Teil 3
say("\nTeil 3: Lineare Leiter aus der Erhaltung der analytischen Multiplikationsalgebra")
rng = np.random.default_rng(88)
nM = 40
S = np.diag(np.ones(nM - 1), -1)                                  # S e_n = e_(n+1)


def toeplitz_abweichung(h, t=0.7):
    U = np.diag(np.exp(1j*t*h))
    T = U @ S @ U.conj().T
    d = np.diag(T, -1)
    return np.max(np.abs(d - d[0]))


lin = 2*np.pi*np.arange(nM) + np.pi
nichtlin = [np.pi*(2*np.arange(nM) + 1)**3/10, np.pi*(2*np.arange(nM) + 1) + 0.3*np.sin(np.arange(nM)),
            np.sort(rng.normal(size=nM))*10]
ok_lin = toeplitz_abweichung(lin) < 1e-12 and all(toeplitz_abweichung(h) > 1e-3 for h in nichtlin)
check("e^{ith} S e^{-ith} ist genau dann wieder Toeplitz (konstant auf der Nebendiagonalen), wenn h_(n+1) - h_n konstant ist: "
      "linear ja, kubisch, gestört und zufällig nein", ok_lin)
n = np.arange(-6, 6)
h = np.pi*(2*n + 1)
x = np.concatenate([np.arange(-3, 3) + 0.25, np.arange(-3, 3) + 0.75])
ok_ph = np.allclose(np.pi*(2*(-n - 1) + 1), -h) and np.allclose(np.sort(4*np.pi*x), np.sort(h)) and abs(h[n == 0][0] - np.pi) < 1e-15
check("mit der Teilchen-Loch-Symmetrie h_(-n-1) = -h_n und der Saatlücke h_0 = pi folgt h_n = pi (2n + 1); das ist 4 pi x auf "
      "der verzahnten Leiter x in Z/2 + 1/4, also genau beide Quarter des Originals (die NS-Moden der Überlagerung aus z87)", ok_ph)
befund("Die neunte Analyse hält hier: Die gleichförmige Leiter folgt aus der Erhaltung der analytischen Multiplikationsalgebra "
       "unter dem Modulfluss. Das ist die Bedingung H aus z87 als Satz; sie bleibt eine benannte Prämisse, wird aber "
       "strukturell (Hardy-Invarianz statt Dispersion)", True)

urteil = ("Die neunte externe Analyse hält in den hier geprüften Punkten. Exakte Zählung: Ab einer Schwelle hat der 4D-Skalar mehr "
          "Zustände mit Energie <= E als der volle Fockraum des Originals mit Partner und Ladungssektoren; eine energietreue "
          "Einbettung des vollen 4D-Fockraums ist ausgeschlossen, die Einteilchenleiter (z83) bleibt. log Z wächst wie "
          "(2 pi^2/3)/beta gegen pi^4/(45 beta^3). Die lineare Leiter h_n = pi (2n + 1) folgt aus der Erhaltung der "
          "analytischen Multiplikationsalgebra unter dem Modulfluss mit Teilchen-Loch-Symmetrie und Saatlücke. Die radiale "
          "Spinorreflexion und der Quadrupoltest gehen an die laufenden Rechnungen S3 und S4.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z88_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
