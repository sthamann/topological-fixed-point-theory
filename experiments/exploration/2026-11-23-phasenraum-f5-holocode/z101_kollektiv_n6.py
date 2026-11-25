"""Z101: Kollektive Nichtgaußizität ab N = 6 als verborgene Wechselwirkung? (Folgeauftrag N2, Hypothese H2).

Leitlinie 'Rekonstruktion statt Auswahl'. Firewalls: DIMENSION.UPLIFT.FIREWALL.01 (Saaten und Nähte sind 2D- bzw.
Fockkinematik, 4D nur als Lesart über die radiale Leiter von z83 und die bosonische Komposition von z91; keine
wechselwirkende 3+1-Schließung), GRAV.NONCIRCULAR.01 (der Modulfluss ist Uhr und Helizitätsphase, nicht die Zeit).
Gelöste Blöcke (Flavour, Massen, alpha, 16, E8) bleiben unberührt; nichts wird befördert. Bestand: z94 (N = 6 gleiche
Saaten: diagonaler Zustand y^(|n| + 4 f(n)), y = e^-pi, Ränder und Uhr exakt; N = 2..5 starr im Fenster |f| <= 8),
z95 Teil 3 (symmetrische f ohne Lösung), z91 (hartkernige Komposition, kappa_4 = -1/(2N); 'lokal frei, kollektiv
korreliert' als Zustand ohne Wechselwirkung), z99 (verbundene 1->2-Amplitude nur mit paarzahländerndem Zweikörperterm;
Singulett-Vierpunkt -2/N rein statistisch), z98 (Tomita-Familie = AdS4, MacDowell-Mansouri-Krümmung null), z83
(radiale Leiter, Schritte 'Teilchen hoch' und 'Loch runter'). Statuskürzel: E exakt, C bedingt, L Lesart, O offen.

Teil 1  Die Familie aller exakten Diagonallösungen. Struktursatz für jedes N (Stufen A_L = {n: e(n) = L} sind
        gleichmäßig überdeckt, gamma_(L+1) = n_L - gamma_L; keine Lücken, also endlich), vollständige Suche ohne
        Fenster für N = 2..5, Profilsuche für N = 2..8. Für N = 6 vollständige Klassifikation: Das Stufenprofil ist
        erzwungen, jede Lösung ist eine Wahl von Modenfolge, zwei geordneten 1-Faktorisierungen von K6 und einer
        geordneten Tripelpartition; z94 ist ein Element. Isospektralität und universelle S_6-invariante Daten,
        k-Körper-Zerlegung von log rho - log rho_Produkt (Möbius, ganzzahlig), Symmetrie (S_6, Komplementspiegelung,
        Uhr), Kumulanten des Kollektivfelds (transversal wie z91, longitudinal).
Teil 2  Radiale Lesart: jede Naht mit der radialen Leiter (Stufen m = 1, 2, hartkernig, z83/z91), die kollektive
        Korrelation auf den Nahtbesetzungen. z99-Kriterium für verbundene 1->2-Amplituden (native Einkörperoperatoren,
        modularer Hamiltonoperator, Kontrolle), Modulfluss der Erzeuger (Kleidung und Uhrzeiten), Drei- und
        Vierpunktantworten gegen die Statistik (-2/N), Transport unter freier radialer Dynamik gegen einen Zweikörperterm.
Teil 3  Tomita-Familie und MM-Krümmung: J ist zustandsunabhängig (Standardform), Familie unter der Einkörpergruppe
        U(6) (Stabilisator, Bahn, Maurer-Cartan, Torsion), Connes-Kozyklus; Bilanz und Urteil H2 (Folgen für K20, T3/T5).
"""
import itertools
import json
import math
import time
from collections import Counter
from pathlib import Path

import mpmath as mp
import numpy as np
import scipy.linalg as sla
import scipy.sparse as sps
from scipy.sparse.linalg import expm_multiply

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
T0 = time.time()
rng = np.random.default_rng(101)
mp.mp.dps = 50
YM = mp.exp(-mp.pi)
YN = float(YM)
P1M = YM/(1 + YM)                                                      # Randbesetzung einer Saat
P1 = float(P1M)
N = 6
NC = 1 << N
VOLL = NC - 1
POP = np.array([bin(m).count("1") for m in range(NC)])
BIT = np.array([[(m >> k) & 1 for k in range(N)] for m in range(NC)])  # Bit k = Saat (Naht) k
EMAX = 22


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")
    return bool(cond)


def befund(name, val):
    BEFUNDE[name] = bool(val)
    say(f"   [Befund: {'ja' if val else 'nein'}] {name}")


def laufzeit():
    say(f"   (Laufzeit bisher {time.time() - T0:.1f} s)")


def maske(s):
    return sum(int(ch) << k for k, ch in enumerate(s))


def wort(m, n=N):
    return "".join(str((m >> k) & 1) for k in range(n))


def multimengen_ok(e, n=N):
    """Randverteilung jeder Saat exakt <=> {e(b): b_k = 1} = {e(b) + 1: b_k = 0} als Multimengen (e^pi transzendent)."""
    nc = 1 << n
    return all(Counter(e[m] for m in range(nc) if (m >> k) & 1) == Counter(e[m] + 1 for m in range(nc) if not (m >> k) & 1)
               for k in range(n))


def uhr_ok(e, n=N):
    pop = [bin(m).count("1") for m in range(1 << n)]
    return all((e[m] - pop[m]) % 4 == 0 for m in range(1 << n))


def gewichte(e):
    w = YN**np.asarray(e, dtype=float)
    return w/w.sum()


def paarkorrelationen(p):
    nb = BIT.T @ p
    return {(i, j): float(p @ (BIT[:, i]*BIT[:, j]) - nb[i]*nb[j]) for i, j in itertools.combinations(range(N), 2)}


def moebius(f):
    """c_S mit f(n) = sum_S c_S prod_(i in S) n_i (Möbiusumkehr auf dem Teilmengenverband, ganzzahlig)."""
    c = np.array(f, dtype=np.int64).copy()
    for k in range(N):
        for m in range(NC):
            if (m >> k) & 1:
                c[m] -= c[m ^ (1 << k)]
    return c


# =================================================================== Teil 1
say("Teil 1: Die Familie der exakten Diagonallösungen für N gleiche Saaten")

# --- 1.1 Reproduktion von z94 (Gegenbeispiel N = 6)
z94 = json.loads((HERE / "z94_ergebnis.json").read_text())
f94roh = z94["zahlen"]["aussenbehauptung_gleiche_saaten"]["f(n) für N = 6"]
S94_alt = z94["zahlen"]["aussenbehauptung_gleiche_saaten"]["Beispiel N = 6 (drei Paare m = 1)"]["S"]
F94 = np.zeros(NC, dtype=np.int64)
for s_, v_ in f94roh.items():
    F94[maske(s_)] = v_
E94 = POP + 4*F94
pm94 = [YM**int(x) for x in E94]
Zm94 = mp.fsum(pm94)
pm94 = [x/Zm94 for x in pm94]
rand94 = max(abs(mp.fsum(pm94[m] for m in range(NC) if (m >> k) & 1) - P1M) for k in range(N))
nb94 = [mp.fsum(pm94[m] for m in range(NC) if (m >> k) & 1) for k in range(N)]
C94m = {(i, j): mp.fsum(pm94[m] for m in range(NC) if (m >> i) & 1 and (m >> j) & 1) - nb94[i]*nb94[j]
        for i, j in itertools.combinations(range(N), 2)}
S94 = -mp.fsum(x*mp.log(x) for x in pm94)
S0 = -N*(P1M*mp.log(P1M) + (1 - P1M)*mp.log(1 - P1M))
n_paar94 = sum(abs(v) > mp.mpf(10)**-45 for v in C94m.values())
ZAHLEN["z94_reproduktion"] = {"Randfehler (50 Stellen)": float(rand94), "Paarkorrelationen != 0": n_paar94,
                              "S": float(S94), "S_Produkt": float(S0), "f-Bereich": [int(F94.min()), int(F94.max())],
                              "verschobene Konfigurationen": int(np.sum(F94 != 0))}
check(f"Reproduktion z94 (N = 6, f aus z94_ergebnis.json): Multimengen je Saat exakt ({multimengen_ok(E94)}), Uhr "
      f"e = |n| mod 4 exakt ({uhr_ok(E94)}), Randfehler {float(rand94):.0e} bei 50 Stellen, {n_paar94} von 15 "
      f"Paarkorrelationen != 0, S = {float(S94):.6f} gegen S_Produkt = {float(S0):.6f}, f in [{F94.min()}, {F94.max()}], "
      f"{int(np.sum(F94 != 0))} Konfigurationen verschoben (wie z94)",
      multimengen_ok(E94) and uhr_ok(E94) and rand94 < mp.mpf(10)**-45 and n_paar94 == 15
      and abs(float(S94) - S94_alt) < 1e-12 and np.sum(F94 != 0) == 51)

# --- 1.2 Struktursatz: Stufen gleichmäßig überdeckt, Rekursion, Endlichkeit; vollständige Suche N = 2..5; Profile
STUFEN94 = {L: [m for m in range(NC) if E94[m] == L] for L in range(EMAX + 1)}
gam, gam_ok = [0], True
for L in range(1, EMAX + 2):
    vor = STUFEN94.get(L - 1, [])
    soll = [sum(1 for m in vor if not (m >> k) & 1) for k in range(N)]
    ist = [sum(1 for m in STUFEN94.get(L, []) if (m >> k) & 1) for k in range(N)]
    gam_ok &= soll == ist and len(set(ist)) == 1 and ist[0] == len(vor) - gam[-1]
    gam.append(ist[0])
luecken_ok = all(len(STUFEN94[L]) > 0 for L in range(EMAX + 1)) and STUFEN94[0] == [0] and STUFEN94[EMAX] == [VOLL]


def vollsuche(n, tl=60):
    """Alle Lösungen ohne Fenster: Stufe für Stufe, A_L mit Spaltensummen c_k = #{b in A_(L-1): b_k = 0}, Gewichte
    = L mod 4, 1^n nur als letzte Stufe. Endlich, weil die Stufen lückenlos sind (e <= 2^n - 1)."""
    nc, voll = 1 << n, (1 << n) - 1
    pc = [bin(b).count("1") for b in range(nc)]
    bits = [[(b >> k) & 1 for k in range(n)] for b in range(nc)]
    loes, t0 = [], time.time()

    def teilmengen(kand, c):
        out, r, gew = [], list(c), []
        suf = [[0]*n for _ in range(len(kand) + 1)]
        for i in range(len(kand) - 1, -1, -1):
            suf[i] = [suf[i + 1][k] + bits[kand[i]][k] for k in range(n)]

        def rek(i):
            if not any(r):
                out.append(list(gew))
                return
            if i == len(kand) or any(r[k] > suf[i][k] for k in range(n)):
                return
            bb = bits[kand[i]]
            if all(r[k] >= bb[k] for k in range(n)):
                for k in range(n):
                    r[k] -= bb[k]
                gew.append(kand[i])
                rek(i + 1)
                gew.pop()
                for k in range(n):
                    r[k] += bb[k]
            rek(i + 1)
        rek(0)
        return out

    def dfs(L, vor, benutzt, e):
        if time.time() - t0 > tl:
            raise TimeoutError
        c = [sum(1 for b in vor if not bits[b][k]) for k in range(n)]
        if not any(c):
            if benutzt == (1 << nc) - 1:
                loes.append(dict(e))
            return
        rest = [b for b in range(nc) if not (benutzt >> b) & 1]
        if rest == [voll] and all(x == 1 for x in c) and pc[voll] % 4 == L % 4:
            dfs(L + 1, [voll], benutzt | (1 << voll), {**e, voll: L})
            return
        for S in teilmengen([b for b in rest if b != voll and pc[b] % 4 == L % 4], c):
            u = benutzt
            for b in S:
                u |= 1 << b
            dfs(L + 1, S, u, {**e, **{b: L for b in S}})
    dfs(1, [0], 1, {0: 0})
    return loes


voll_klein = {n: vollsuche(n) for n in range(2, 6)}
voll_klein_ok = all(len(v) == 1 and all(v[0][b] == bin(b).count("1") for b in range(1 << n)) for n, v in voll_klein.items())


def profile(n, tl=60):
    """Notwendige Bedingungen auf Anzahlen: je Stufe Gewichte w = L mod 4 mit sum w x_w = n gamma_L,
    gamma_(L+1) = sum x_w - gamma_L, 1^n allein als oberste Stufe. Ergebnis: mögliche Stufenprofile (Spektren)."""
    rem0 = {w: math.comb(n, w) for w in range(n + 1)}
    rem0[0] -= 1
    res, t0 = [], time.time()

    def rek(L, g, rem, prof):
        if time.time() - t0 > tl:
            raise TimeoutError
        if g == 0:
            if not any(rem.values()):
                res.append(tuple(prof))
            return
        ws = [w for w in range(1, n + 1) if w % 4 == L % 4 and rem[w] > 0]
        if n in ws and g == 1 and not any(v for w, v in rem.items() if w != n):
            rek(L + 1, 0, {**rem, n: 0}, prof + [1])
        wl = [w for w in ws if w != n]

        def komp(i, rest):
            if i == len(wl):
                if rest == 0:
                    yield ()
                return
            for x in range(min(rem[wl[i]], rest//wl[i]) + 1):
                for c in komp(i + 1, rest - wl[i]*x):
                    yield ((wl[i], x),) + c
        for c in komp(0, n*g):
            nl = sum(x for _, x in c)
            if nl == 0 or nl < g:
                continue
            r2 = dict(rem)
            for w, x in c:
                r2[w] -= x
            rek(L + 1, nl - g, r2, prof + [nl])
    rek(1, 1, rem0, [1])
    return res


PROF = {n: profile(n) for n in range(2, 9)}
prof_tab = {n: {"Profile": len(v), "e_max": sorted({len(p) - 1 for p in v})[:6]} for n, v in PROF.items()}
prof6 = sorted(PROF[6], key=len)
PROFIL_KOLL = list(prof6[-1]) if len(prof6) == 2 else []
ZAHLEN["struktursatz"] = {"gamma_L (z94)": gam[:EMAX + 2], "Stufengrößen n_L (z94)": [len(STUFEN94[L]) for L in range(EMAX + 1)],
                          "vollständige Suche ohne Fenster N = 2..5 (Lösungen)": {n: len(v) for n, v in voll_klein.items()},
                          "Profilsuche (notwendige Anzahlen)": prof_tab, "kollektives Profil N = 6": PROFIL_KOLL}
check(f"Struktursatz (jedes N): Die Randbedingung je Saat lautet stufenweise #{{b in A_L: b_k = 1}} = #{{b in A_(L-1): b_k = "
      f"0}}; aus A_0 = {{0}} folgt induktiv, dass jede Stufe jede Saat gleich oft enthält, gamma_L-mal mit gamma_(L+1) = n_L - "
      f"gamma_L (an z94 geprüft: {gam_ok}, gamma = {gam[:9]}...); unterste Stufe {{0}}, oberste {{1^N}}, keine Lücke "
      f"({luecken_ok}), also e <= 2^N - 1 und die Lösungsmenge endlich. Vollständige Suche ohne Fenster: N = 2..5 haben "
      f"genau eine Lösung, das Produkt ({voll_klein_ok}; z94 hatte |f| <= 8). Profilsuche (notwendige Anzahlbedingungen): "
      f"N = 2..5 nur das Produktprofil, N = 6 genau zwei Profile mit e_max = {prof_tab[6]['e_max']}; N = 7: "
      f"{prof_tab[7]['Profile']}, N = 8: {prof_tab[8]['Profile']} Profile (Realisierbarkeit ab N = 7 nicht geprüft)",
      gam_ok and luecken_ok and voll_klein_ok and all(prof_tab[n]["Profile"] == 1 for n in range(2, 6))
      and prof_tab[6]["Profile"] == 2 and prof_tab[6]["e_max"] == [6, EMAX]
      and PROFIL_KOLL == [len(STUFEN94[L]) for L in range(EMAX + 1)])
befund("E: Die Schwelle N = 6 ist 5 + 1. Auf Stufe 1 muss jede Saat genau einmal überdeckt werden, mit Blöcken der Größe "
       "1 mod 4 (Uhr). Bis N = 5 gibt es nur Einzelquanten (ein Fünferblock wäre 1^5 selbst, also die Spitze), und jede "
       "Restklasse mod 4 enthält außer 0 und 1^N nur ein Gewicht: Jede Stufe ist voll, der Zustand das Produkt. Ab N = 6 "
       "darf ein Fünferkompositum auf derselben Modulstufe sitzen wie ein Einzelquant (5 = 1 mod 4)", True)
laufzeit()

# --- 1.3 Vollständige Klassifikation N = 6
MODI = list(range(N))
PAARE = list(itertools.combinations(MODI, 2))
TRIPEL = list(itertools.combinations(MODI, 3))


def matchings(vs):
    if not vs:
        yield ()
        return
    for b in vs[1:]:
        for m in matchings([v for v in vs if v not in (vs[0], b)]):
            yield ((vs[0], b),) + m


MATCH = [tuple(sorted(m)) for m in matchings(MODI)]
FAKT = []


def _fakt(gew, benutzt):
    rest = [p for p in PAARE if p not in benutzt]
    if not rest:
        FAKT.append(tuple(sorted(gew)))
        return
    for m in MATCH:
        if rest[0] in m and not set(m) & benutzt:
            _fakt(gew + [m], benutzt | set(m))


_fakt([], set())
QUAD = [q for q in itertools.combinations(TRIPEL, 4) if all(sum(v in t for t in q) == 2 for v in MODI)]
TPART = []


def _tpart(gew, benutzt):
    rest = [t for t in TRIPEL if t not in benutzt]
    if not rest:
        TPART.append(tuple(sorted(gew)))
        return
    for q in QUAD:
        if rest[0] in q and not set(q) & benutzt:
            _tpart(gew + [q], benutzt | set(q))


_tpart([], set())
msk = lambda xs: sum(1 << x for x in xs)


def baue(folge, fp, fv, tp):
    """Lösung aus den Bausteinen: Stufe 4j+1 = {e_k, 1 - e_k} (k = folge[j]), 4j+2 = Matching fp[j], 4j+3 = Tripel-
    quadrupel tp[j], 4j+4 = Komplemente des Matchings fv[j] (j = 0..4), Stufe 22 = 1^6."""
    e = np.full(NC, -1, dtype=np.int64)
    e[0], e[VOLL] = 0, EMAX
    for j, k in enumerate(folge):
        e[1 << k] = e[VOLL ^ (1 << k)] = 4*j + 1
    for j in range(5):
        for p in fp[j]:
            e[msk(p)] = 4*j + 2
        for t in tp[j]:
            e[msk(t)] = 4*j + 3
        for p in fv[j]:
            e[VOLL ^ msk(p)] = 4*j + 4
    return e


def zufall_loesung():
    folge = [int(x) for x in rng.permutation(N)]
    a_, b_, c_ = FAKT[rng.integers(len(FAKT))], FAKT[rng.integers(len(FAKT))], TPART[rng.integers(len(TPART))]
    fp = [a_[i] for i in rng.permutation(5)]
    fv = [b_[i] for i in rng.permutation(5)]
    tp = [c_[i] for i in rng.permutation(5)]
    e = baue(folge, fp, fv, tp)
    assert np.all(e >= 0)
    return e


# z94 zerlegen und nachbauen
z_folge = [next(k for k in MODI if (1 << k) in STUFEN94[4*j + 1]) for j in range(6)]
z_fp = [tuple(sorted(tuple(k for k in MODI if (m >> k) & 1) for m in STUFEN94[4*j + 2])) for j in range(5)]
z_tp = [tuple(sorted(tuple(k for k in MODI if (m >> k) & 1) for m in STUFEN94[4*j + 3])) for j in range(5)]
z_fv = [tuple(sorted(tuple(k for k in MODI if not (m >> k) & 1) for m in STUFEN94[4*j + 4])) for j in range(5)]
z94_form = (all(sorted(STUFEN94[4*j + 1]) == sorted([1 << k, VOLL ^ (1 << k)]) for j, k in enumerate(z_folge))
            and tuple(sorted(z_fp)) in FAKT and tuple(sorted(z_fv)) in FAKT and tuple(sorted(z_tp)) in TPART
            and np.array_equal(baue(z_folge, z_fp, z_fv, z_tp), E94))
N_KOLL = math.factorial(6)*(len(FAKT)*math.factorial(5))**2*(len(TPART)*math.factorial(5))
STICH = [zufall_loesung() for _ in range(300)]
stich_ok = all(multimengen_ok(e) and uhr_ok(e) and e[0] == 0 for e in STICH)
korr_fail = 0
for e in STICH[:150]:
    e2 = e.copy()
    a_, b_ = rng.choice([m for m in range(NC) if POP[m] == 2 and e[m] == 2]), rng.choice(
        [m for m in range(NC) if POP[m] == 2 and e[m] == 6])
    e2[a_], e2[b_] = e[b_], e[a_]
    korr_fail += not multimengen_ok(e2) and uhr_ok(e2)
for e in STICH[150:]:
    e2 = e.copy()
    a_, b_ = rng.choice([m for m in range(NC) if POP[m] == 3 and e[m] == 3]), rng.choice(
        [m for m in range(NC) if POP[m] == 3 and e[m] == 11])
    e2[a_], e2[b_] = e[b_], e[a_]
    korr_fail += not multimengen_ok(e2) and uhr_ok(e2)
ZAHLEN["klassifikation_N6"] = {"perfekte Matchings K6": len(MATCH), "1-Faktorisierungen K6": len(FAKT),
                               "2-reguläre Tripelquadrupel": len(QUAD), "Tripelpartitionen in fünf Quadrupel": len(TPART),
                               "kollektive Lösungen (f(0) = 0)": N_KOLL, "z94: Modenfolge der Stufen 4j+1": z_folge,
                               "Stichprobe exakt": f"{len(STICH)}/{len(STICH)}" if stich_ok else "nein",
                               "Störungen (Tausch zwischen Stufen) verletzen die Ränder": f"{korr_fail}/{len(STICH)}"}
check(f"Klassifikation N = 6 (exakt): Das kollektive Profil n_L = 1, 2, 3, 4, 3, 2, 3, 4, 3, ..., 2, 1 (e_max = 22) ist "
      f"erzwungen, gamma = 1 auf allen Stufen 4j+1, 4j+2 und 2 auf 4j+3, 4j+4: Stufe 4j+1 = {{e_k, 1 - e_k}}, 4j+2 = "
      f"perfektes Matching, 4j+3 = vier Tripel mit jeder Saat zweimal, 4j+4 = Komplemente eines Matchings. Die Restklassen "
      f"mod 4 entkoppeln, also ist jede Lösung eine freie Wahl aus Modenfolge (6!), zwei geordneten 1-Faktorisierungen von "
      f"K6 ({len(FAKT)} x 5!) und einer geordneten Partition der 20 Tripel in fünf 2-reguläre Quadrupel ({len(TPART)} x 5!, "
      f"aus {len(QUAD)} Quadrupeln): {N_KOLL:.4e} kollektive Lösungen. z94 hat genau diese Form und wird nachgebaut "
      f"({z94_form}); {len(STICH)} Zufallslösungen sind exakt (Multimengen und Uhr), {korr_fail} von {len(STICH)} "
      f"Störungen (Tausch zweier Konfigurationen gleichen Gewichts zwischen Stufen, Uhr bleibt) verletzen die Ränder",
      len(MATCH) == 15 and len(FAKT) == 6 and len(QUAD) == 75 and len(TPART) == 1350 and z94_form and stich_ok
      and korr_fail == len(STICH) and N_KOLL == 60466176000000)

# --- 1.4 Isospektralität, universelle Daten, k-Körper-Zerlegung
FW94 = [int(sum(F94[m] for m in range(NC) if POP[m] == w)) for w in range(N + 1)]
ORD = lambda c: [int(sum(c[m] for m in range(NC) if POP[m] == k)) for k in range(1, N + 1)]
C94 = moebius(F94)
ord94 = ORD(C94)
gemeinsam = Counter((int(POP[m]), int(E94[m])) for m in range(NC))
iso_ok, ord_ok, eins_ok, sechs_ok, int_ok, fuenf_ok = True, True, True, True, True, True
nz_k = {k: [] for k in range(1, N + 1)}
cmin = {k: 10**9 for k in range(1, N + 1)}
cmax = {k: -10**9 for k in range(1, N + 1)}
for e in STICH + [E94]:
    f = (e - POP)//4
    c = moebius(f)
    iso_ok &= Counter((int(POP[m]), int(e[m])) for m in range(NC)) == gemeinsam
    ord_ok &= ORD(c) == ord94
    eins_ok &= sorted(int(c[1 << k]) for k in range(N)) == list(range(6))
    sechs_ok &= c[VOLL] == 0
    fuenf_ok &= all(c[VOLL ^ (1 << k)] == 4 for k in range(N))
    int_ok &= np.all((e - POP) % 4 == 0)
    for k in range(1, N + 1):
        sel = [c[m] for m in range(NC) if POP[m] == k]
        nz_k[k].append(sum(1 for x in sel if x != 0))
        cmin[k], cmax[k] = min(cmin[k], int(min(sel))), max(cmax[k], int(max(sel)))
prof_zahl = Counter(int(x) for x in E94)
Zpoly = [prof_zahl[L] for L in range(EMAX + 1)]
Sp = -sum(nl*(YM**L/mp.fsum(Zpoly[k]*YM**k for k in range(EMAX + 1)))*mp.log(YM**L/mp.fsum(
    Zpoly[k]*YM**k for k in range(EMAX + 1))) for L, nl in enumerate(Zpoly))
PM = [mp.fsum(YM**int(E94[m]) for m in range(NC) if POP[m] == w) for w in range(N + 1)]
ZM = mp.fsum(PM)
PM = [x/ZM for x in PM]
ZAHLEN["universell"] = {"Stufenprofil n_L": Zpoly, "S (aus dem Profil)": float(Sp),
                        "F_w = Summe f über |n| = w": FW94, "Ordnungssummen sum_(|S|=k) c_S, k = 1..6": ord94,
                        "P(M = w) kollektiv": [float(x) for x in PM],
                        "P(M = 5)/P(M = 1)": float(PM[5]/PM[1]),
                        "Möbius z94: Zahl der c_S != 0 je k": [sum(1 for m in range(NC) if POP[m] == k and C94[m]) for k in
                                                               range(1, N + 1)],
                        "Stichprobe: Zahl der c_S != 0 je k (min, max)": {k: [min(v), max(v)] for k, v in nz_k.items()},
                        "Stichprobe: Koeffizientenbereich je k": {k: [cmin[k], cmax[k]] for k in range(1, N + 1)}}
check(f"Isospektral und universell (exakt, z94 und {len(STICH)} Zufallslösungen): Alle kollektiven Lösungen haben dasselbe "
      f"Stufenprofil, also dasselbe Spektrum (S = {float(Sp):.6f} für alle), und sogar dieselbe gemeinsame Verteilung von "
      f"(|n|, e) ({iso_ok}): P(M = w), F_w = {FW94} und die S_6-invarianten Ordnungssummen der Möbiuskoeffizienten "
      f"sum_(|S|=k) c_S = {ord94} (k = 1..6) sind universell ({ord_ok}). log rho - log rho_Produkt = -4 pi sum_S c_S n_S "
      f"(+ const) mit ganzen c_S ({int_ok}): k = 1..5 treten in jeder Lösung auf, k = 6 nie (c_123456 = 0: {sechs_ok}); der "
      f"Einkörperteil ist stets eine Treppe c_k = Permutation von 0..5 ({eins_ok}), die Fünfkörperkoeffizienten sind "
      f"universell c_S = 4 für alle sechs Fünfermengen ({fuenf_ok}; aus der Regularität: (15 - j_k) - 20 + 20 - 10 + "
      f"(j_k - 1) = 4), die Ordnungen 2..4 hängen von der Lösung ab. Exakt P(M = 5) = P(M = 1) = {float(PM[1]):.6f} "
      f"(Produkt: Verhältnis y^4)",
      iso_ok and ord_ok and eins_ok and sechs_ok and fuenf_ok and int_ok and abs(Sp - S94) < mp.mpf(10)**-40
      and ord94 == [15, -45, 70, -60, 24, 0] and abs(PM[5] - PM[1]) < mp.mpf(10)**-45)
befund("E: Die Nichtgaußizität hat ein festes, quantisiertes Spektrum (eindeutig für N = 6) und ganzzahlige "
       "Vielkörperkoeffizienten in Einheiten 4 pi; frei ist nur ihre Lage, eine diskrete Wahl aus 6.05e13 Lösungen. "
       "Lesart: Einkörperig sind die sechs Saaten zu einer Treppe 4 pi j (j = 0..5) verstimmt, vielkörperig ist das "
       "kompensiert; führend ist das Fünferkompositum 1 - e_k mit dem Modulgewicht eines Einzelquants (e = 1)", True)

# --- 1.5 Symmetrie: S_6, Komplementspiegelung, Uhr; S_6-Mittel
PERMS = list(itertools.permutations(MODI))
PMAP = np.array([[sum(((m >> k) & 1) << p[k] for k in range(N)) for m in range(NC)] for p in PERMS])
stab_groessen = []
spieg_ok, spieg_stab = True, 0
for e in [E94] + STICH[:100]:
    stab_groessen.append(int(np.sum(np.all(e[PMAP] == e[None, :], axis=1))))
    eR = EMAX - e[VOLL ^ np.arange(NC)]
    spieg_ok &= multimengen_ok(eR) and uhr_ok(eR) and eR[0] == 0
    spieg_stab += int(np.any(np.all(eR[PMAP] == e[None, :], axis=1)))
folge_R = [next(k for k in MODI if EMAX - E94[VOLL ^ (1 << k)] == 4*j + 1) for j in range(6)]
# S_6-Mittel: rho-quer(n) = P(M = |n|)/C(6, |n|)
rq = [PM[int(POP[m])]/math.comb(N, int(POP[m])) for m in range(NC)]
rand_q = max(abs(mp.fsum(rq[m] for m in range(NC) if (m >> k) & 1) - P1M) for k in range(N))
phase_q = [mp.exp(-0.5j*mp.log(rq[m]))*(1j)**(-int(POP[m])) for m in range(NC)]
uhr_q = max(abs(ph - phase_q[0]) for ph in phase_q)
phase_k = [mp.exp(-0.5j*mp.log(pm94[m]))*(1j)**(-int(POP[m])) for m in range(NC)]
uhr_k = max(abs(ph - phase_k[0]) for ph in phase_k)
ZAHLEN["symmetrie"] = {"S_6-Stabilisator (z94 und 100 Lösungen)": dict(Counter(stab_groessen)),
                       "Komplementspiegelung e -> 22 - e(1 - n) ist Lösung": spieg_ok,
                       "Lösungen mit Stabilisator in S_6 x Spiegelung": spieg_stab,
                       "z94: Modenfolge vorher / nach Spiegelung": [z_folge, folge_R],
                       "S_6-Mittel: Randfehler": float(rand_q), "S_6-Mittel: Uhrabweichung": float(uhr_q),
                       "kollektiv: Uhrabweichung": float(uhr_k)}
check(f"Symmetrie: Die S_6 wirkt frei (Stabilisatorgrößen {dict(Counter(stab_groessen))} bei z94 und 100 Lösungen; Beweis: "
      f"die Stufen 4j+1 ordnen die Saaten), die Komplementspiegelung e -> 22 - e(1 - n) bildet Lösungen auf Lösungen ab "
      f"({spieg_ok}) und kehrt die Modenfolge um ({z_folge} -> {folge_R}); {spieg_stab} der 101 Lösungen sind unter "
      f"Spiegelung mal Permutation invariant. Die Uhr gilt für jede Lösung exakt (Abweichung {float(uhr_k):.0e}); das "
      f"S_6-Mittel hat exakte Ränder ({float(rand_q):.0e}), verletzt aber die Uhr (Abweichung {float(uhr_q):.3f}): "
      f"S_N-Symmetrie und Uhr schließen sich für kollektive Zustände aus, das Produkt ist der einzige symmetrische",
      set(stab_groessen) == {1} and spieg_ok and folge_R == z_folge[::-1] and uhr_k < 1e-40 and rand_q < 1e-45
      and uhr_q > 0.1)

# --- 1.6 Kumulanten des Kollektivfelds
sx = sps.csr_matrix(np.array([[0., 1.], [1., 0.]]))


def lokal2(op, k):
    mats = [sps.identity(2, format="csr")]*N
    mats[N - 1 - k] = op                                               # Bit k = Saat k
    out = mats[0]
    for M_ in mats[1:]:
        out = sps.kron(out, M_, format="csr")
    return out


XS = sum(lokal2(sx, k) for k in range(N))/np.sqrt(2*N)
XP = [sps.identity(NC, format="csr")]
for _ in range(4):
    XP.append((XP[-1] @ XS).tocsr())


def kappa_x(p):
    mu = [float(p @ XP[r].diagonal()) for r in range(5)]
    return mu[2], mu[3], mu[4] - 3*mu[2]**2


kap_x = {"Produkt": kappa_x(gewichte(POP)), "z94": kappa_x(gewichte(E94))}
for i_, e in enumerate(STICH[:3]):
    kap_x[f"Lösung {i_ + 1}"] = kappa_x(gewichte(e))
kap_x_ok = all(abs(v[0] - 0.5) < 1e-13 and abs(v[1]) < 1e-13 and abs(v[2] + 1/(2*N)) < 1e-13 for v in kap_x.values())


def kumulanten(P):
    mu = [mp.fsum(P[w]*mp.mpf(w)**r for w in range(N + 1)) for r in range(5)]
    return [mu[1], mu[2] - mu[1]**2, mu[3] - 3*mu[2]*mu[1] + 2*mu[1]**3,
            mu[4] - 4*mu[3]*mu[1] - 3*mu[2]**2 + 12*mu[2]*mu[1]**2 - 6*mu[1]**4]


PP = [math.comb(N, w)*P1M**w*(1 - P1M)**(N - w) for w in range(N + 1)]
kP, kK = kumulanten(PP), kumulanten(PM)
sig2 = N*P1M*(1 - P1M)
SUMC = kK[1] - sig2                                                   # sum_(i != j) C_ij
ZAHLEN["kumulanten"] = {"transversal x = sum sigma_x/sqrt(2N) (kappa_2, kappa_3, kappa_4)": {k: [round(x, 14) for x in v] for
                                                                                         k, v in kap_x.items()},
                        "-1/(2N)": -1/(2*N), "longitudinal M: kappa_1..4 Produkt": [float(x) for x in kP],
                        "longitudinal M: kappa_1..4 kollektiv": [float(x) for x in kK],
                        "standardisiert kollektiv kappa_3/kappa_2^1.5, kappa_4/kappa_2^2": [float(kK[2]/kK[1]**1.5),
                                                                                        float(kK[3]/kK[1]**2)],
                        "standardisiert Produkt": [float(kP[2]/kP[1]**1.5), float(kP[3]/kP[1]**2)],
                        "sum_(i != j) C_ij": float(SUMC)}
check(f"Kumulanten: Das transversale Kollektivfeld von z91, x = sum sigma_x/sqrt(2N), hat in jedem Diagonalzustand "
      f"kappa_2 = 1/2, kappa_3 = 0, kappa_4 = -1/(2N) = {-1/(2*N):.6f} (Produkt, z94, drei Lösungen; nur gepaarte "
      f"Umklappungen zählen): Es sieht die kollektive Korrelation nicht. Das longitudinale Feld M = sum n_i sieht sie "
      f"universell: kappa_1 = {float(kK[0]):.6f} (= N p, Ränder), kappa_2 = {float(kK[1]):.5f} gegen {float(kP[1]):.5f}, "
      f"kappa_3 = {float(kK[2]):.4f} gegen {float(kP[2]):.5f}, kappa_4 = {float(kK[3]):.3f} gegen {float(kP[3]):.5f}; "
      f"standardisiert kappa_4/kappa_2^2 = {float(kK[3]/kK[1]**2):.3f} gegen {float(kP[3]/kP[1]**2):.3f} im Produkt "
      f"(unabhängige Saaten, reine Statistik) und -2/N = {-2/N:.4f} für x (z91)",
      kap_x_ok and abs(kK[0] - N*P1M) < mp.mpf(10)**-45 and kK[1] > 3*kP[1] and kK[3] > 10*kP[3]
      and abs(kP[1] - sig2) < mp.mpf(10)**-45)
befund("E: Einkörperig unsichtbar (Ränder, Uhr, Einteilchendaten, transversales Feld von z91), zweikörperig sichtbar: "
       "alle 15 Paarkorrelationen, die Zahlverteilung P(M) und die longitudinalen Kumulanten. 'Lokal unsichtbar' gilt "
       "nur für eine Saat allein", True)
laufzeit()

# =================================================================== Teil 2
say("\nTeil 2: Radiale Lesart (Naht = hartkerniges Quant auf der radialen Leiter, kollektive Besetzung)")
# lokale Basis: 0 Vakuum, 1 |m=1>, 2 |m=2, Loch runter>, 3 |m=2, Teilchen hoch> (z83, nackte Schritte)
DL = 4
DIM = DL**N
E_LOK = np.diag([0., 1., 2., 2.])


def ket(i, j, d=DL):
    M_ = np.zeros((d, d))
    M_[i, j] = 1
    return sps.csr_matrix(M_)


def lokal4(op, i):
    mats = [sps.identity(DL, format="csr")]*N
    mats[i] = sps.csr_matrix(op)                                       # Naht 0 ist die höchste Stelle
    out = mats[0]
    for M_ in mats[1:]:
        out = sps.kron(out, M_, format="csr")
    return out


ZUST = np.array(list(itertools.product(range(DL), repeat=N)))         # Zeile = Nahtzustände (Naht 0 zuerst)
BES = (ZUST > 0).astype(int)
BMASK = BES @ (1 << np.arange(N))
NQ = BES.sum(1)
ZWEI = NQ == 2
OM = np.zeros(DIM)
OM[0] = 1
IDX_SAAT = np.array([sum(((m >> k) & 1)*DL**(N - 1 - k) for k in range(N)) for m in range(NC)])   # Saat m -> m = 1-Stufe
ERAD = np.array([E_LOK.diagonal()[ZUST[:, i]] for i in range(N)]).sum(0)


def K_rad(f):
    return sps.diags(np.pi*ERAD + 4*np.pi*np.asarray(f)[BMASK]).tocsr()


BD = {A: ket(A, 0) for A in (1, 2, 3)}
QUELLEN = [(i, A, lokal4(BD[A], i)) for i in range(N) for A in (1, 2, 3)]
NATIV = [sps.diags(ERAD).tocsr()]
for i in range(N):
    for op in (ket(2, 1), ket(3, 1)):
        NATIV += [lokal4(op, i), lokal4(op.T, i)]
for i, j in itertools.permutations(range(N), 2):
    for A, B in itertools.product((1, 2, 3), repeat=2):
        NATIV.append((lokal4(ket(A, 0), i) @ lokal4(ket(0, B), j)).tocsr())


def verbunden_1_2(X):
    best = 0.0
    xo = X @ OM
    for _, _, Bq in QUELLEN:
        d_ = X @ (Bq @ OM) - Bq @ xo
        best = max(best, float(np.abs(d_[ZWEI]).max()))
    return best


v_nativ = max(verbunden_1_2(X) for X in NATIV)
K94 = K_rad(F94)
v_K = max(verbunden_1_2(K_rad((e - POP)//4)) for e in [E94] + STICH[:5])
V_ctrl = (lokal4(ket(1, 0), 1) @ lokal4(ket(1, 0), 2) @ lokal4(ket(0, 1), 0)).tocsr()
v_ctrl = verbunden_1_2(V_ctrl)
kd = K94.diagonal()
e1 = [IDX_SAAT[1 << i] for i in range(N)]
zz = {(i, j): kd[IDX_SAAT[(1 << i) | (1 << j)]] - kd[e1[i]] - kd[e1[j]] + kd[0] for i, j in itertools.combinations(range(N), 2)}
zz_ok = all(abs(v - 4*np.pi*C94[(1 << i) | (1 << j)]) < 1e-10 for (i, j), v in zz.items())
n_zz = sum(1 for v in zz.values() if abs(v) > 1e-9)
ZAHLEN["z99_kriterium_radial"] = {"Dimension": DIM, "native Einkörperoperatoren": len(NATIV),
                                  "max verbunden 1->2 nativ": v_nativ, "max verbunden 1->2 K_koll (z94 + 5)": v_K,
                                  "Kontrolle b+_1 b+_2 b_0": v_ctrl,
                                  "K_koll 2->2 diagonal verbunden / 4 pi (z94)": {f"{i}{j}": round(v/(4*np.pi), 9) for (i, j), v
                                                                                  in zz.items()}}
check(f"z99-Kriterium radial (sechs Nähte, Stufen m = 1, 2, Dimension {DIM}): Der Zweiquantenanteil von X b+_A Omega - "
      f"b+_A X Omega verschwindet für alle {len(NATIV)} nativen Einkörperoperatoren (Energie, Heber 'Teilchen hoch'/'Loch "
      f"runter', Transfer zwischen Nähten; max {v_nativ:.0e}) und für den modularen Hamiltonoperator K_koll = pi E + 4 pi "
      f"f(n) (z94 und fünf Lösungen; max {v_K:.0e}): K_koll ist je Naht besetzungsdiagonal. Er hat verbundene "
      f"2->2-Diagonalamplituden (Dichte-Dichte) genau 4 pi c_ij ({zz_ok}; {n_zz} von 15 Paaren bei z94), ohne Transfer. "
      f"Kontrolle: b+_1 b+_2 b_0 gibt verbunden {v_ctrl:g}",
      v_nativ < 1e-14 and v_K < 1e-14 and abs(v_ctrl - 1) < 1e-14 and zz_ok and n_zz > 0)

# --- 2.2 Modulfluss der Erzeuger: Kleidung mit Zuschauerbesetzungen, Uhrzeiten
BS = []
for k in range(N):
    M_ = np.zeros((NC, NC))
    for m in range(NC):
        if not (m >> k) & 1:
            M_[m | (1 << k), m] = 1
    BS.append(M_)


def fluss_verbunden(e, t):
    D_ = np.exp(-1j*np.pi*t*np.asarray(e, dtype=float))
    best, formel = 0.0, 0.0
    c = moebius((np.asarray(e) - POP)//4)
    for k, i in itertools.permutations(range(N), 2):
        sig = (D_[:, None]*BS[k])/D_[None, :]
        a_ = sig[(1 << k) | (1 << i), 1 << i] - (BS[i] @ sig)[(1 << k) | (1 << i), 0]
        best = max(best, abs(a_))
        formel = max(formel, abs(abs(a_) - abs(np.exp(-4j*np.pi*t*c[(1 << k) | (1 << i)]) - 1)))
    return best, formel


TS_GEN, TS_UHR = (0.1, 0.25, 0.37), (0.5, 1.0, 1.5)
fl_gen = {t: fluss_verbunden(E94, t) for t in TS_GEN}
fl_uhr = {t: fluss_verbunden(E94, t) for t in TS_UHR}
fl_prod = max(fluss_verbunden(POP, t)[0] for t in TS_GEN + TS_UHR)
frq = sorted(int(E94[1 << k] - E94[0]) for k in range(N))
ZAHLEN["modulfluss"] = {"max |verbunden 1->2| von sigma_t(b+_k), generisch": {t: v[0] for t, v in fl_gen.items()},
                        "Uhrzeiten t in Z/2": {t: v[0] for t, v in fl_uhr.items()},
                        "Abweichung von |e^(-4 pi i t c_ki) - 1|": max(v[1] for v in list(fl_gen.values()) + list(fl_uhr.values())),
                        "Produkt (alle t)": fl_prod, "Modulfrequenzen der Erzeuger auf dem Vakuum (Einheiten pi)": frq}
check(f"Modulfluss (Helizitätsphase, nicht Zeit): sigma_t(b+_k) = rho^(it) b+_k rho^(-it) ist mit den Zuschauerbesetzungen "
      f"gekleidet; seine verbundene 1->2-Amplitude im Sinn von z99 ist exakt |e^(-4 pi i t c_ki) - 1| (Abweichung "
      f"{max(v[1] for v in list(fl_gen.values()) + list(fl_uhr.values())):.0e}), bei t = 0,1, 0,25, 0,37 bis "
      f"{max(v[0] for v in fl_gen.values()):.4f}, bei den Uhrzeiten t = 1/2, 1, 3/2 null "
      f"({max(v[0] for v in fl_uhr.values()):.0e}); im Produkt für alle t null ({fl_prod:.0e}). Auf dem Vakuum hat der "
      f"Erzeuger der Saat k die Modulfrequenz pi (1 + 4 j_k) = pi x {frq}: sechs 'gleiche' Saaten mit sechs Modulenergien",
      max(v[0] for v in fl_gen.values()) > 1 and max(v[0] for v in fl_uhr.values()) < 1e-12 and fl_prod < 1e-12
      and max(v[1] for v in list(fl_gen.values()) + list(fl_uhr.values())) < 1e-12 and frq == [1, 5, 9, 13, 17, 21])
befund("E/L: Die kollektive Struktur wirkt wie eine verborgene Wechselwirkung genau im Modulfluss: Er kleidet jeden "
       "Erzeuger mit den Besetzungen der anderen Nähte (Frequenzen pi (1 + 4j), j ganz) und erzeugt dabei verbundene "
       "1->2-Amplituden mit festen ganzzahligen Ladungen c_ki. Zu den Uhrzeiten t in Z/2 verschwindet die Kleidung exakt, "
       "die Uhr sieht nichts. Der Modulfluss ist aber Uhr und Helizitätsphase, nicht die Zeit (GRAV.NONCIRCULAR.01)", True)

# --- 2.3 Drei- und Vierpunktantworten gegen die Statistik
BSd = [b.T for b in BS]
NS = [np.diag(BIT[:, k].astype(float)) for k in range(N)]
OPS = BS + BSd
p94 = gewichte(E94)
drei_max = 0.0
for a_, b_, c_ in itertools.product(range(2*N), repeat=3):
    drei_max = max(drei_max, abs(float(np.trace(np.diag(p94) @ OPS[a_] @ OPS[b_] @ OPS[c_]))))
C94f = paarkorrelationen(p94)
misch_ok = True
for i, j in itertools.permutations(range(N), 2):
    w_ = float(p94 @ np.diag(BS[i].T @ NS[j] @ BS[i]))                 # <b_i n_j b+_i>
    v_ = w_ - float(p94 @ np.diag(BS[i].T @ BS[i]))*float(p94 @ BIT[:, j])
    misch_ok &= abs(v_ + C94f[tuple(sorted((i, j)))]) < 1e-15
nb = BIT.T @ p94
kap3 = {}
for i, j, k in itertools.combinations(range(N), 3):
    Eijk = p94 @ (BIT[:, i]*BIT[:, j]*BIT[:, k])
    Eij, Eik, Ejk = (p94 @ (BIT[:, a]*BIT[:, b]) for a, b in ((i, j), (i, k), (j, k)))
    kap3[(i, j, k)] = Eijk - Eij*nb[k] - Eik*nb[j] - Ejk*nb[i] + 2*nb[i]*nb[j]*nb[k]
n_kap3 = sum(1 for v in kap3.values() if abs(v) > 1e-14)
PHId = sum(BS)/np.sqrt(N)                                              # Phi^+ (Erzeuger), Phi = PHId^T


def vierpunkt(p):
    rho = np.diag(p)
    zwei = np.trace(rho @ PHId.T @ PHId.T @ PHId @ PHId)               # <Phi Phi Phi+ Phi+>
    ein = np.trace(rho @ PHId.T @ PHId)
    return float(zwei - 2*ein**2)


vak = np.zeros(NC)
vak[0] = 1
G_vak, G_prod, G_94 = vierpunkt(vak), vierpunkt(gewichte(POP)), vierpunkt(p94)
G_st = [vierpunkt(gewichte(e)) for e in STICH[:20]]
XI = G_94 - G_prod
XI_formel = 2/N**2*float(SUMC)
ZAHLEN["drei_vierpunkt"] = {"Paarfeld-Dreipunkt max (1728 Wörter)": drei_max, "<b_i n_j b+_i>_c = -C_ij": misch_ok,
                            "Dichte-Dreipunkt kappa(n_i, n_j, n_k) != 0 (von 20)": n_kap3,
                            "Dichte-Dreipunkt Bereich": [float(min(kap3.values())), float(max(kap3.values()))],
                            "Singulett-Vierpunkt Vakuum": G_vak, "Produkt": G_prod, "kollektiv (z94)": G_94,
                            "Überschuss Xi = (2/N^2) sum_(i != j) C_ij": XI, "Xi Formel": XI_formel,
                            "Streuung über 20 Lösungen": max(G_st) - min(G_st),
                            "Statistik -2(1-p)^2/N": -2*(1 - P1)**2/N}
check(f"Drei- und Vierpunkt (Saaten auf m = 1, Nähte bosonisch komponiert): Paarfeld-Dreipunkte verschwinden exakt "
      f"(1728 Wörter, max {drei_max:.0e}, Nahtzahl). Verbunden und über die Statistik hinaus sind nur Dichteantworten: "
      f"<b_i n_j b+_i>_c = -C_ij ({misch_ok}), {n_kap3} von 20 Dichte-Dreipunkten kappa(n_i, n_j, n_k) != 0. "
      f"Singulett-Vierpunkt (Phi = sum b_i/sqrt N): Vakuum {G_vak:.6f} = -2/N (z99), Produkt {G_prod:.6f} = -2(1-p)^2/N, "
      f"kollektiv {G_94:.6f}: Überschuss Xi = (2/N^2) sum_(i != j) C_ij = {XI:.6f} (Formel {XI_formel:.6f}), universell "
      f"über 20 Lösungen (Streuung {max(G_st) - min(G_st):.0e}), positiv und {XI/abs(G_prod):.3f} der statistischen Größe",
      drei_max < 1e-15 and misch_ok and n_kap3 > 0 and abs(G_vak + 2/N) < 1e-14 and abs(G_prod + 2*(1 - P1)**2/N) < 1e-14
      and abs(XI - XI_formel) < 1e-12 and max(G_st) - min(G_st) < 1e-12 and XI > 0)

# --- 2.4 Transport unter freier radialer Dynamik gegen einen Zweikörperterm
LAM, GKOP = 0.35, 0.6
h_lok = E_LOK + LAM*(ket(2, 1) + ket(3, 1) + ket(1, 2) + ket(1, 3)).toarray()
H_frei = sum(lokal4(h_lok, i) for i in range(N)).tocsr()
V2 = sum((lokal4(ket(2, 1), i) @ lokal4(ket(2, 1), j)) for i, j in itertools.combinations(range(N), 2))
H_ww = (H_frei + GKOP*(V2 + V2.T)).tocsr()
PHIR = sum(lokal4(ket(0, 1), i) for i in range(N)).tocsr()/np.sqrt(N)
PHIRd = PHIR.T.tocsr()
PHI2, PHI2d = (PHIR @ PHIR).tocsr(), (PHIRd @ PHIRd).tocsr()
BASIS = sps.csr_matrix((np.ones(NC), (IDX_SAAT, np.arange(NC))), shape=(DIM, NC))


def g_t(t):
    return complex(sla.expm(-1j*h_lok*t)[1, 1])


def vier_t(H, p, t):
    """G_c(t) = <Phi(t) Phi(t) Phi+ Phi+> - 2 <Phi(t) Phi+>^2 im Diagonalzustand p (Träger m <= 1)."""
    A = expm_multiply(-1j*H*t, BASIS.toarray())
    B2 = expm_multiply(-1j*H*t, (PHI2d @ BASIS).toarray())
    B1 = expm_multiply(-1j*H*t, (PHIRd @ BASIS).toarray())
    zwei = np.sum(p*np.einsum("ij,ij->j", A.conj(), PHI2 @ B2))
    ein = np.sum(p*np.einsum("ij,ij->j", A.conj(), PHIR @ B1))
    return zwei - 2*ein**2


TT = (0.4, 1.3, 2.2)
pP = gewichte(POP)
quot = {}
for name, H in (("frei", H_frei), ("Zweikörper", H_ww)):
    for zn, p in (("kollektiv", p94), ("Produkt", pP)):
        quot[(name, zn)] = [complex(vier_t(H, p, t)/g_t(t)**2) for t in TT]
abw = {k: max(abs(q - (G_94 if k[1] == "kollektiv" else G_prod)) for q in v) for k, v in quot.items()}
ZAHLEN["transport"] = {"lambda (Leitermischung)": LAM, "g (Zweikörperterm)": GKOP, "t": TT,
                       "|G_c(t)/g(t)^2 - G_c(0)|": {f"{k[0]}, {k[1]}": v for k, v in abw.items()},
                       "|g(t)|": [abs(g_t(t)) for t in TT]}
check(f"Transport, nicht Erzeugung: Unter der freien radialen Dynamik h = E + lambda (Heber + Senker) je Naht (lambda = "
      f"{LAM}, |g(t)| = {', '.join(f'{abs(g_t(t)):.3f}' for t in TT)}) ist G_c(t)/g(t)^2 für t = {TT} konstant gleich dem "
      f"statischen Wert, kollektiv ({abw[('frei', 'kollektiv')]:.0e}) wie im Produkt ({abw[('frei', 'Produkt')]:.0e}): Der "
      f"Überschuss Xi wird mit dem freien Propagator transportiert, mit festem Koeffizienten. Mit einem Zweikörperterm "
      f"g (Heber_i Heber_j + h.c.), g = {GKOP}, ändert sich der Quotient ({abw[('Zweikörper', 'kollektiv')]:.3f} bzw. "
      f"{abw[('Zweikörper', 'Produkt')]:.3f}): Erzeugung braucht Dynamik",
      abw[("frei", "kollektiv")] < 1e-9 and abw[("frei", "Produkt")] < 1e-9 and abw[("Zweikörper", "kollektiv")] > 1e-3
      and abw[("Zweikörper", "Produkt")] > 1e-3)
befund("E/C: In der radialen Lesart ist die kollektive Nichtgaußizität eine feste Anfangskorrelation der Nahtbesetzungen: "
       "keine verbundene 1->2-Amplitude im Erzeuger der Zeit (K_koll und alle nativen Operatoren sind je Naht "
       "zahlerhaltend), verbundene Dichteantworten mit festem, quantisiertem Koeffizienten (Singulett-Überschuss Xi "
       "universell), vom freien Propagator transportiert. Wechselwirkung im Sinn von z99 zeigt nur der Modulfluss, und nur "
       "außerhalb der Uhrzeiten", True)
laufzeit()

# =================================================================== Teil 3
say("\nTeil 3: Tomita-Familie, MacDowell-Mansouri-Krümmung und Urteil H2")
# --- 3.1 J ist zustandsunabhängig (Standardform), blockweise polare Zerlegung
JSOLL = np.array([[0, 0, 1, 0], [0, 0, 0, -1], [1, 0, 0, 0], [0, -1, 0, 0]], dtype=float)


def tomita_bloecke(p):
    """S(A x 1) xi = (A+ x 1) xi mit xi = sum sqrt(p_b) |b, b>: S|a, b> = sqrt(p_a/p_b) |b, a> (antilinear). Reelle
    4 x 4-Blöcke auf (Re z_ab, Im z_ab, Re z_ba, Im z_ba), polar zerlegt (SVD)."""
    J_abw, D_max = 0.0, 0.0
    for a, b in itertools.combinations(range(NC), 2):
        r = math.sqrt(p[a]/p[b])
        M_ = np.array([[0, 0, 1/r, 0], [0, 0, 0, -1/r], [r, 0, 0, 0], [0, -r, 0, 0]])
        W, s, Vt = np.linalg.svd(M_)
        J_ = W @ Vt
        J_abw = max(J_abw, float(np.abs(J_ - JSOLL).max()))
        D_max = max(D_max, abs(math.log(r)))
    return J_abw, D_max


tom = {"Produkt": tomita_bloecke(pP), "z94": tomita_bloecke(p94)}
for i_, e in enumerate(STICH[:2]):
    tom[f"Lösung {i_ + 1}"] = tomita_bloecke(gewichte(e))
U_ = sla.expm(1j*(lambda A: (A + A.conj().T)/2)(rng.normal(size=(NC, NC)) + 1j*rng.normal(size=(NC, NC))))
Vz = rng.normal(size=(NC, NC)) + 1j*rng.normal(size=(NC, NC))
kov = float(np.abs((U_ @ (U_.conj().T @ Vz @ U_).conj().T @ U_.conj().T) - Vz.conj().T).max())
z98 = json.loads((HERE / "z98_ergebnis.json").read_text())
FA98 = z98["zahlen"]["berry"]["|F_A| (MM-Krümmung, Maurer-Cartan)"]
FG98 = z98["zahlen"]["berry"]["|F_Gamma| (Berry-Krümmung)"]
ZAHLEN["tomita"] = {"max |J - Tausch-Konjugation| je Zustand": {k: v[0] for k, v in tom.items()},
                    "max |log Delta^(1/2)| je Zustand": {k: v[1] for k, v in tom.items()},
                    "Kovarianz J (U x U-quer) = (U x U-quer) J": kov, "z98: |F_A|, |F_Gamma|": [FA98, FG98]}
check(f"J ist zustandsunabhängig (Araki, Standardform): Die polare Zerlegung S = J Delta^(1/2) der Tomita-Operatoren im "
      f"verdoppelten Saatraum (2016 reelle 4 x 4-Blöcke |a, b>, |b, a>) gibt für Produkt, z94 und zwei Lösungen dieselbe J = Tausch mal "
      f"Konjugation (max Abweichung {max(v[0] for v in tom.values()):.0e}), während Delta verschieden ist (max |log "
      f"Delta^(1/2)| = {tom['Produkt'][1]/np.pi:.1f} pi gegen {tom['z94'][1]/np.pi:.1f} pi); J vertauscht mit jeder Einkörper"
      f"unitären U x U-quer ({kov:.0e}). Die Tomita-Familie J_g = U(g) J U(g)^-1 hängt also nur von der Gruppe ab: für den "
      f"kollektiven Zustand dieselbe wie in z98 (AdS4, Levi-Civita, MM-Krümmung |F_A| = {FA98:.0e})",
      max(v[0] for v in tom.values()) < 1e-12 and abs(tom["Produkt"][1] - 3*np.pi) < 1e-9
      and abs(tom["z94"][1] - 11*np.pi) < 1e-9 and kov < 1e-10 and FA98 < 1e-8)

# --- 3.2 Familie des modularen Hamiltonoperators unter der Einkörpergruppe U(6) (gemeinsame CAR)
Zp, sm = np.diag([1., -1.]), np.array([[0., 1.], [0., 0.]])
CJW = []
for k in range(N):
    mats = [np.eye(2)]*N
    for q in range(k):
        mats[N - 1 - q] = Zp
    mats[N - 1 - k] = sm
    op = mats[0]
    for M_ in mats[1:]:
        op = np.kron(op, M_)
    CJW.append(op)
EINS = {(i, j): CJW[i].T @ CJW[j] for i in range(N) for j in range(N)}
nok = all(np.allclose(np.diag(EINS[(k, k)]), BIT[:, k]) for k in range(N))


def stabilisator_dim(e, ops=None):
    ops = EINS if ops is None else ops
    K = np.diag(np.pi*np.asarray(e, dtype=float))
    A = np.array([(X @ K - K @ X).ravel() for X in ops.values()]).T
    return int(len(ops) - np.linalg.matrix_rank(A, tol=1e-9))


stab_dim = {"Produkt": stabilisator_dim(POP), "z94": stabilisator_dim(E94)}
for i_, e in enumerate(STICH[:3]):
    stab_dim[f"Lösung {i_ + 1}"] = stabilisator_dim(e)
SEKT = [(0, 1, 2), (3, 4, 5)]                                          # Ladungsauswahl je Quarter: drei + drei Moden
EINS33 = {k: v for k, v in EINS.items() if any(k[0] in s_ and k[1] in s_ for s_ in SEKT)}
stab33 = {"Produkt": stabilisator_dim(POP, EINS33), "z94": stabilisator_dim(E94, EINS33)}
xa = rng.normal(size=(N, N)) + 1j*rng.normal(size=(N, N))
xa = (xa - xa.conj().T)/2
Ubig = sla.expm(sum(xa[i, j]*EINS[(i, j)] for i in range(N) for j in range(N)))
rho_s = Ubig @ np.diag(p94) @ Ubig.conj().T
Gam_s = np.array([[np.trace(rho_s @ EINS[(i, j)]) for j in range(N)] for i in range(N)])
w_s = np.linalg.eigvalsh(rho_s)
Rop = np.diag(1j**POP.astype(float))
uhr_s_abw = float(np.abs(Ubig @ Rop - Rop @ Ubig).max())              # [U, R] = 0: rho_s^(-i/2) R^-1 = U (const) U^+
offd = float(np.abs(rho_s - np.diag(np.diag(rho_s))).max())
spek_ok = np.allclose(np.sort(w_s), np.sort(p94), rtol=0, atol=1e-15)
# Maurer-Cartan, Berry-Krümmung, Torsion auf g = e^(sX) e^(tY), h = Cartan (Stabilisator), m = Nichtdiagonale
Ph_ = lambda X: np.diag(np.diag(X))
Pm_ = lambda X: X - np.diag(np.diag(X))


def anti(A):
    A = (A - A.conj().T)/2
    return Pm_(A)


Xg = 0.5*anti(rng.normal(size=(N, N)) + 1j*rng.normal(size=(N, N)))
Yg = 0.5*anti(rng.normal(size=(N, N)) + 1j*rng.normal(size=(N, N)))


def A_s(t):
    return sla.expm(-t*Yg) @ Xg @ sla.expm(t*Yg)


t0_, hh = 0.41, 1e-5
dAs = (A_s(t0_ + hh) - A_s(t0_ - hh))/(2*hh)
As0 = A_s(t0_)
F_A = -dAs + (As0 @ Yg - Yg @ As0)
G_s, e_s = Ph_(As0), Pm_(As0)
F_G = -Ph_(dAs)
riem = -Ph_(e_s @ Yg - Yg @ e_s)
tors = -Pm_(dAs) + (G_s @ Yg - Yg @ G_s)
tors_soll = -Pm_(e_s @ Yg - Yg @ e_s)
ZAHLEN["einkoerperfamilie"] = {"Stabilisator in u(6) (Dimension)": stab_dim, "Bahn U(6)/T^6 (Dimension)": 36 - stab_dim["z94"],
                               "Stabilisator in u(3) + u(3) (Ladungsauswahl je Quarter)": stab33,
                               "Bahn unter U(3) x U(3) (Dimension)": len(EINS33) - stab33["z94"],
                               "Bahnzustand: max |Gamma - p 1|": float(np.abs(Gam_s - P1*np.eye(N)).max()),
                               "Bahnzustand: Uhrabweichung": uhr_s_abw, "Bahnzustand: Nichtdiagonale": offd,
                               "|F_A| (Maurer-Cartan)": float(np.abs(F_A).max()), "|F_Gamma|": float(np.abs(F_G).max()),
                               "|F_Gamma + P_h[e, e]|": float(np.abs(F_G - riem).max()), "|Torsion|": float(np.abs(tors).max()),
                               "|Torsion + P_m[e, e]|": float(np.abs(tors - tors_soll).max())}
check(f"Einkörperfamilie (Fermilesart von z94: sechs Saaten in gemeinsamer CAR, Erzeuger c_i^+ c_j): Der Stabilisator des "
      f"modularen Hamiltonoperators in u(6) hat die Dimension {stab_dim} (Produkt: ganz u(6), kollektiv: nur der "
      f"Cartan-Teil, weil keine Transposition f festhält); die Bahn U(6)/T^6 hat die Dimension {36 - stab_dim['z94']}, "
      f"mit Ladungsauswahl je Quarter (U(3) x U(3), Stabilisator {stab33['z94']}) die Dimension "
      f"{len(EINS33) - stab33['z94']}: kontinuierliche Moduli isospektraler Zustände mit denselben Einteilchendaten "
      f"(|Gamma - p 1| = "
      f"{np.abs(Gam_s - P1*np.eye(N)).max():.0e}) und exakter Uhr ([U, R] = {uhr_s_abw:.0e}), nicht diagonal "
      f"({offd:.3f}). Auf der "
      f"Bahn (z98-Konstruktion, g = e^(sX) e^(tY)): Maurer-Cartan-Krümmung {np.abs(F_A).max():.0e}, Berry-Krümmung der "
      f"Nullmoden F_Gamma = -P_h[e, e] != 0 ({np.abs(F_G).max():.3f}, Rest {np.abs(F_G - riem).max():.0e}), Torsion "
      f"-P_m[e, e] != 0 ({np.abs(tors).max():.3f}, Rest {np.abs(tors - tors_soll).max():.0e}): reduktiv, nicht symmetrisch, "
      f"als Cartan-Geometrie flach",
      nok and stab_dim["Produkt"] == 36 and all(v == 6 for k, v in stab_dim.items() if k != "Produkt")
      and stab33 == {"Produkt": 18, "z94": 6} and np.abs(Gam_s - P1*np.eye(N)).max() < 1e-12 and uhr_s_abw < 1e-12 and offd > 1e-3 and spek_ok
      and np.abs(F_A).max() < 1e-7 and np.abs(F_G).max() > 1e-2 and np.abs(F_G - riem).max() < 1e-7
      and np.abs(tors).max() > 1e-2 and np.abs(tors - tors_soll).max() < 1e-7)

# --- 3.3 Connes-Kozyklus kollektiv gegen Produkt
def kozyklus(t):
    u = np.exp(-1j*np.pi*t*(E94 - POP).astype(float))                   # rho_koll^(it) rho_prod^(-it) bis auf Phase
    return u/u[0]


koz = {t: float(np.abs(kozyklus(t) - 1).max()) for t in (0.1, 0.25, 0.5, 1.0)}
s_, t_ = 0.13, 0.29
koz_char = float(np.abs(kozyklus(s_ + t_) - kozyklus(s_)*kozyklus(t_)).max())
ZAHLEN["connes_kozyklus"] = {"max |u_t - 1| (normiert)": koz, "Charakter u_(s+t) = u_s u_t": koz_char,
                             "Ladungen f (ganzzahlig)": sorted(set(int(x) for x in F94))}
check(f"Connes-Kozyklus (D rho_koll : D rho_prod)_t = e^(-4 pi i t f(n)) bis auf Phase: ein Charakter mit ganzen "
      f"Ladungen f ({koz_char:.0e}), abelsch und flach; trivial genau zu den Uhrzeiten (|u - 1| bei t = 1/2, 1: "
      f"{koz[0.5]:.0e}, {koz[1.0]:.0e}), sonst nicht (t = 0,1: {koz[0.1]:.3f}, t = 1/4: {koz[0.25]:.0f})",
      koz[0.5] < 1e-12 and koz[1.0] < 1e-12 and koz[0.1] > 0.1 and abs(koz[0.25] - 2) < 1e-12 and koz_char < 1e-12)
befund("E: Keine native Tomita-Familie des kollektiven Zustands trägt Krümmung. J ist zustandsunabhängig, also ist die "
       "z98-Familie dieselbe (F = 0); die Nichtgaußizität steckt allein in Delta. Gruppenbahnen von Delta sind homogen "
       "und als Cartan-Geometrie flach (unter U(6) eine reduktive, nicht symmetrische Bahn mit Torsion); die relative "
       "Modulkette ist ein Charakter, der zu den Uhrzeiten trivial ist; die Diagonalfamilie ist diskret. F != 0 bräuchte "
       "eine nicht homogene, stetige Familie mit ortsabhängigem Koeffizienten, also Dynamik", True)
laufzeit()

# --- Bilanz
say("\nBilanz")
STATUS = {
    "Familie N = 6": ("E", f"vollständig klassifiziert: erzwungenes Profil, {N_KOLL:.3e} Lösungen = Modenfolge x zwei "
                           "geordnete 1-Faktorisierungen von K6 x Tripelpartition; N = 2..5 nur das Produkt (ohne Fenster)"),
    "Schwelle N = 6": ("E", "5 + 1: Fünferkompositum auf der Modulstufe eines Einzelquants (5 = 1 mod 4)"),
    "Spektrum und S_6-invariante Daten": ("E", "eindeutig (isospektral), S = 0,3720, P(M = 5) = P(M = 1), Ordnungssummen "
                                               "(15, -45, 70, -60, 24, 0)"),
    "k-Körper-Zerlegung": ("E", "k = 1..5 in jeder Lösung, k = 6 nie; ganzzahlig in Einheiten 4 pi; Einkörperteil = Treppe, "
                               "Fünfkörperteil universell c_S = 4, Ordnungen 2..4 lösungsabhängig"),
    "Symmetrie": ("E", "Uhr exakt; S_6 frei (kein Stabilisator), Komplementspiegelung; S_6-Mittel verletzt die Uhr"),
    "Kollektivfeld": ("E", "transversal (z91) blind: -1/(2N) für jeden Diagonalzustand; longitudinal universell sichtbar"),
    "verbundene 1->2 (Zeit)": ("E, negativ", "K_koll und native Einkörperoperatoren je Naht zahlerhaltend"),
    "Modulfluss": ("E/L", "kleidet Erzeuger, verbundene 1->2-Amplitude |e^(-4 pi i t c) - 1|, null zu den Uhrzeiten"),
    "Drei- und Vierpunkt": ("E", "Paarfelder null; Dichteantworten verbunden; Singulett-Überschuss Xi universell, fest"),
    "Transport": ("E", "freie radiale Dynamik transportiert Xi, ein Zweikörperterm verändert es"),
    "Tomita-Familie und F": ("E, negativ", "J zustandsunabhängig, z98-Familie identisch, F = 0; U(6)-Bahn flach mit Torsion"),
    "Moduli": ("E", "U(6)-Bahn der Dimension 30 (mit Ladungsauswahl je Quarter 12): isospektral, Einteilchendaten und Uhr "
                    "exakt, nicht diagonal"),
    "Auswahl": ("E/O", "Maximalentropie wählt das Produkt (S = 1,0347 gegen 0,3720); ein natives Auswahlprinzip für eine "
                       "kollektive Lösung ist nicht bekannt"),
}
for k_, (s_, t_) in STATUS.items():
    say(f"   [{s_}] {k_}: {t_}")
ZAHLEN["status"] = STATUS
urteil = ("H2 trägt nur als Zustandsstruktur, nicht als Wechselwirkung (bedingt, nicht nativ). (1) Die Familie ist "
          "vollständig bestimmt: Für N gleiche Saaten ist jede Stufe A_L = {n: e(n) = L} gleichmäßig überdeckt "
          "(gamma_(L+1) = n_L - gamma_L), die Stufen sind lückenlos, die Lösungsmenge endlich. N = 2..5 haben ohne Fenster "
          "nur das Produkt. Die Schwelle N = 6 ist 5 + 1: Die Uhr mod 4 erlaubt ein Fünferkompositum auf der Modulstufe "
          "eines Einzelquants. Für N = 6 ist das Stufenprofil erzwungen; jede der 6,05e13 kollektiven Lösungen ist eine "
          "Wahl von Modenfolge, zwei geordneten 1-Faktorisierungen von K6 und einer geordneten Partition der 20 Tripel in "
          "2-reguläre Quadrupel, z94 ist eine davon. Alle sind isospektral (S = 0,3720), und alle S_6-invarianten Daten "
          "sind universell: P(M), die Ordnungssummen der Möbiuskoeffizienten (15, -45, 70, -60, 24, 0), P(M = 5) = "
          "P(M = 1). log rho - log rho_Produkt = -4 pi sum c_S n_S hat ganze Koeffizienten mit k = 1..5, nie k = 6; der "
          "Einkörperteil ist eine Treppe 4 pi j (j = 0..5), die vielkörperig kompensiert wird, und der Fünfkörperteil ist "
          "universell (c_S = 4 für jede Fünfermenge). Die Uhr gilt exakt, die "
          "S_6 wirkt frei, und das S_6-Mittel verletzt die Uhr: Kollektivität bricht die Saatsymmetrie notwendig. Das "
          "transversale Kollektivfeld von z91 ist blind (kappa_4 = -1/(2N) für jeden Diagonalzustand), das longitudinale "
          "sieht die Korrelation universell. 'Lokal unsichtbar' gilt nur für eine einzelne Saat; alle 15 Paarkorrelationen "
          "sind sichtbar. (2) Radial gelesen ist die Struktur eine feste Anfangskorrelation der Nahtbesetzungen: Der "
          "modulare Hamiltonoperator und alle nativen Einkörperoperatoren sind je Naht zahlerhaltend, also keine verbundene "
          "1->2-Amplitude (z99-Kriterium exakt null); K_koll hat nur Dichte-Dichte-Terme 4 pi c_ij ohne Transfer. "
          "Verbundene Drei- und Vierpunktantworten über die Statistik hinaus gibt es nur im Dichtekanal; der "
          "Singulett-Vierpunkt erhält zum statistischen -2(1-p)^2/N einen universellen, festen Überschuss Xi = (2/N^2) "
          "sum C_ij, den die freie radiale Dynamik transportiert und erst ein Zweikörperterm verändert. Wechselwirkungsartig "
          "ist allein der Modulfluss: Er kleidet die Erzeuger mit den Zuschauerbesetzungen (Frequenzen pi (1 + 4j)) und hat "
          "verbundene 1->2-Amplituden |e^(-4 pi i t c_ki) - 1|, die zu den Uhrzeiten exakt verschwinden; nach "
          "GRAV.NONCIRCULAR.01 ist er aber Uhr und Helizitätsphase, nicht Zeit. (3) Die Tomita-Familie gibt kein F != 0: "
          "J ist zustandsunabhängig, also ist die z98-Familie identisch (AdS4, F = 0); die Nichtgaußizität steckt nur in "
          "Delta. Unter der Einkörpergruppe U(6) bildet der kollektive Zustand eine 30-dimensionale Bahn isospektraler, "
          "exakter Zustände (mit Ladungsauswahl je Quarter 12-dimensional; neue stetige Moduli außerhalb der "
          "Diagonalfamilie), deren Cartan-Verbindung flach ist "
          "(reduktiv, nicht symmetrisch, mit Torsion); der Connes-Kozyklus ist ein Charakter mit ganzen Ladungen f, "
          "trivial zu den Uhrzeiten. Urteil H2: verborgene Wechselwirkung nein als Dynamik [E, gescoptes Negativ]; als "
          "Zustand mit festem, quantisiertem Spektrum und universellem Singulett-Überschuss ja, aber nicht nativ gewählt "
          "(Maximalentropie wählt das Produkt, die Lage ist eine diskrete Wahl plus U(6)-Moduli) [C mit einer "
          "Auswahlprämisse]. Folgen: K20 erhält keinen gemeinsamen Elternprozess, weder F != 0 noch eine verbundene "
          "1->2-Amplitude in der Zeit; der Zustand liefert aber einen nativen Mechanismus für feste, ganzzahlige "
          "Vielkörperkoeffizienten (Uhr mod 4), der in eine Dynamik übertragen werden müsste. T3/T5 bleiben offen: Eine "
          "quasilokale wechselwirkende Hamiltonfamilie folgt nicht. Nichts ist befördert.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z101_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
