"""Z177: Weitere Kandidaten, die alles zusammenführen könnten, und die Arbeit arXiv:2611.34924 gegen den Stand nach z174.

Leitlinie 'Rekonstruktion statt Auswahl'; H6 und die gelösten Bausteine bleiben unberührt, nichts wird befördert.
Statuskürzel: E exakt, N numerisch, C bedingt, L Lesart, Lit ungeprüfte Literatur, O offen.
Bestand: z113 (Modultest), z174 (Standardform; die Z4-Verklebung ist auf einem Paar Naht + Partner anomal und über zwei
Paare modulinvariant und positiv; neue Annahmen P_std und P_paar; beide Negativbefunde zeigen auf 'ein Paar von Nähten'),
z166 (erste Durchsicht der Arbeit; offener Abgleich 2 pi gegen 4 pi), z175 (Zuordnung) und z176 (Ursprung), die parallel
zu dieser Rechnung liefen.

Teil A  Wie viele Nähte, und wie verklebt? Modultest der Methode z113/z174 mit robuster Auswertung (Thetafunktionen und
        eta über die exakten Transformationsgesetze in den Fundamentalbereich): Nahtzahl N = 1..12 und die zyklische
        Verklebung zweier Kopien (eine Spur durch beide Kopien, Öffnungswinkel 4 pi).
Teil B  Standardform des Nahtpaars im endlichen Modell: kanonische Implementierung des Gatters zwischen zwei Nähten.
        Trägt das Paar die Wechselwirkung, die eine Naht mit Partner nicht trägt (z174 A7)?
Teil C  arXiv:2611.34924 (Chandrasekaran, 'An algebraic quantization of gravity from the null cone'), jetzt mit den
        Abschnitten 4.3 bis 6 und den Anhängen D, E gelesen: nachgerechnete Formeln, Wörterbuch zur Naht, Abgleich
        2 pi gegen 4 pi, erweiterte Spintabelle der Verklebung.
Teil D  Kandidatentabelle, Eintrag ins Suchprogramm (OriginLab reconstruction_search) und Abstand zur Komplettlösung.

Die Hilfsfunktionen des Modultests und des endlichen Modells werden aus z174 (gleicher Ordner) übernommen.
"""
import json
import re
import sys
import time
from fractions import Fraction as Fr
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.linalg import expm

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parents[2] / "research" / "originlab"))
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


def zeit():
    return f"(Laufzeit bisher {time.time() - T_START:.0f} s)"


SRC174 = (HERE / "z174_standardform_statistik_komposition.py").read_text()
Z174 = json.loads((HERE / "z174_ergebnis.json").read_text())

# =================================================================== Teil A
say("Teil A: Nahtzahl und Art der Verklebung im Modultest (Methode z113/z174, robuste Auswertung)")
ns = {"mp": mp, "Fr": Fr, "say": lambda *a: None}
exec(SRC174.split("# =================================================================== Teil B")[1]
     .split("# --- B1 Bausteine und Normalformen")[0], ns)
I, PI, R_ = ns["I"], ns["PI"], ns["R_"]
theta_dir, eta_dir = ns["theta"], ns["eta"]
HALB = Fr(1, 2)


def theta_r(a, b, w):
    """theta[a;b](w) für beliebiges w der oberen Halbebene: mit theta[a;b](w+1) = e^(-i pi a(a+1)) theta[a;a+b+1/2](w) und
    theta[a;b](-1/w) = sqrt(-i w) e^(2 pi i a b) theta[b;-a](w) in den Fundamentalbereich, dort die Reihe."""
    a, b = Fr(a), Fr(b)
    f = mp.mpc(1)
    for _ in range(400):
        m = int(mp.nint(w.real))
        if m:
            f *= mp.exp(-I*PI*R_(a*(a + 1))*m)
            b = b + m*(a + HALB)
            w = w - m
        a = a - round(a)
        k = round(b)
        if k:
            f *= mp.exp(2*I*PI*R_(a)*k)
            b = b - k
        if abs(w) < mp.mpf('0.999'):
            w = -1/w
            f *= mp.sqrt(-I*w)*mp.exp(2*I*PI*R_(a*b))
            a, b = b, -a
        else:
            break
    return f*theta_dir(a, b, w)


def eta_r(w):
    f = mp.mpc(1)
    for _ in range(400):
        m = int(mp.nint(w.real))
        if m:
            f *= mp.exp(I*PI*m/12)
            w = w - m
        if abs(w) < mp.mpf('0.999'):
            w = -1/w
            f *= mp.sqrt(-I*w)
        else:
            break
    return f*eta_dir(w)


rob = mp.mpf(0)
for (a_, b_) in [(Fr(1, 4), Fr(1, 4)), (Fr(3, 4), Fr(-1, 4)), (Fr(0), HALB), (Fr(1, 4), Fr(0))]:
    for w_ in [mp.mpc('0.137', '0.21'), mp.mpc('-2.31', '0.09'), mp.mpc('5.4', '0.5'), mp.mpc('0.49', '0.13')]:
        mp.mp.dps = 130
        d_, e_ = theta_dir(a_, b_, w_), eta_dir(w_)
        mp.mp.dps = 50
        rob = max(rob, abs(theta_r(a_, b_, w_) - d_)/abs(d_), abs(eta_r(w_) - e_)/abs(e_))
ns["theta"], ns["eta"] = theta_r, eta_r
naht, bahn, charaktere, konsistenz, act = (ns[k] for k in ("naht", "bahn", "charaktere", "konsistenz", "act"))
NAMEN = {"1": "1", "G": "G = U R^2", "X": "X = R i^(-2h)"}
B = {k: bahn(naht(v)) for k, v in NAMEN.items()}
YS = (mp.mpf('0.35'), mp.mpf('0.6'), mp.mpf('1.0'))
KF = 24
GRID = [mp.mpc(mp.mpf(k)/KF, Y) for Y in YS for k in range(KF)]
EPS = -1e-8
NULLCH = {"T": Fr(0), "S": Fr(0)}


def gitter(Bx, dps=50):
    mp.mp.dps = dps
    out = [[Bx["F"](act(g, t)) for t in GRID] for g in Bx["reps"]]
    mp.mp.dps = 50
    return out


def moden(Bx, gw, N=1):
    """Spinmoden |s| <= 5 der Bahnsumme von F^N mit trivialem Charakter. Die Bahn von F^N hat dieselben Elemente wie die
    von F (sind zwei Bahnelemente nicht proportional, dann auch ihre N-ten Potenzen nicht); die Phasen werden N-fach."""
    BN = dict(Bx)
    BN["ph"] = {g: {k: (N*r) % 1 for k, r in d.items()} for g, d in Bx["ph"].items()}
    V, bad = konsistenz(BN, NULLCH)
    if bad:
        return None
    tot = [mp.mpc(0)]*len(GRID)
    for k in range(Bx["n"]):
        ph = mp.exp(2*I*PI*R_(V[k]))
        tot = [t + ph*x**N for t, x in zip(tot, gw[k])]
    return {Y: {s: mp.fsum(x*mp.exp(-2*I*PI*s*k/KF) for k, x in enumerate(tot[iy*KF:(iy + 1)*KF]))/KF
                for s in range(-5, 6)} for iy, Y in enumerate(YS)}


def comb(ms, cs):
    return {Y: {s: sum(c*m_[Y][s] for m_, c in zip(ms, cs)) for s in ms[0][Y]} for Y in YS}


def konj(m_):
    return {Y: {s: mp.conj(m_[Y][s]) for s in m_[Y]} for Y in YS}


def minimum(m_, m0):
    gross = max(abs(m_[Y][s]) for Y in YS for s in m_[Y])
    return (float(min(m_[Y][s].real/abs(m0[Y][0]) for Y in YS for s in m_[Y])),
            float(max(abs(m_[Y][s].imag) for Y in YS for s in m_[Y])/gross))


def tafel(GWx, NS):
    out = {}
    for N in NS:
        M = {k: moden(B[k], GWx[k], N) for k in B}
        m0 = M["1"]
        z = {"G = +1": minimum(comb([m0, M["G"]], [1, 1]), m0), "G = -1": minimum(comb([m0, M["G"]], [1, -1]), m0)}
        if M["X"]:
            for j in range(4):
                z[f"X = i^{j}"] = minimum(comb([m0, M["G"], M["X"], konj(M["X"])], [1, (-1)**j, (-I)**j, I**j]), m0)
        out[N] = z
    return out


GW = {k: gitter(B[k]) for k in B}
say(f"   Gitterwerte (50 Stellen) berechnet {zeit()}")
NS = tuple(range(1, 13))
TAF = tafel(GW, NS)
ref = Z174["zahlen"]["B3_positivitaet (min normierter Spinmode, rel. Imaginärteil; Y = 0.35, 0.6, 1.0, |s| <= 5)"]
paare = [(TAF[N]["G = +1"][0], ref[str(N)]["Z_1 + Z_G"][0]) for N in (1, 2, 3, 4)] + \
        [(TAF[N][f"X = i^{j}"][0], ref[str(N)][f"X = i^{j}"][0]) for N in (2, 4) for j in range(4)]
kontrolle = max(abs(x - y)/max(1.0, abs(y)) for x, y in paare)
bahnen = {k: Bx["n"] for k, Bx in B.items()}
ZAHLEN["A1"] = {"robust_gegen_direkt": float(rob), "Bahnen": bahnen, "Kontrolle_gegen_z174": kontrolle}
check("A1 [E/N] Robuste Auswertung: theta[a;b] und eta werden mit den exakten Transformationsgesetzen in den "
      f"Fundamentalbereich gebracht und stimmen mit den direkten Reihen bei 130 Stellen überein (Fehler {mp.nstr(rob, 2)}). "
      f"Bahnen wie in z174 ({bahnen}); die kleinsten Spinmoden für N = 1..4 reproduzieren z174, das die langen Nähte über "
      f"Normalformen auswertet (größte relative Abweichung {kontrolle:.0e})",
      rob < mp.mpf(10)**-45 and bahnen == {"1": 6, "G": 24, "X": 96} and kontrolle < 1e-6)

GW34 = {k: gitter(B[k], 34) for k in B}
TAF34 = tafel(GW34, NS)
stabil = all((TAF[N][k][0] > EPS) == (TAF34[N][k][0] > EPS) for N in NS for k in TAF[N])
relabw = max(abs(TAF[N][k][0] - TAF34[N][k][0])/max(1.0, abs(TAF[N][k][0])) for N in NS for k in TAF[N])
frei = [N for N in NS if "X = i^0" in TAF[N]]
positiv = {N: [j for j in range(4) if TAF[N][f"X = i^{j}"][0] > EPS] for N in frei}
gtest = {"G = +1": [N for N in NS if TAF[N]["G = +1"][0] > EPS], "G = -1": [N for N in NS if TAF[N]["G = -1"][0] > EPS]}
reell = max(v[1] for z in TAF.values() for v in z.values())
for N in NS:
    say(f"      N = {N}: " + "; ".join(f"{k}: {v[0]:.2e}" for k, v in TAF[N].items()))
ZAHLEN["A2_tafel (min normierter Spinmode, rel. Imaginärteil)"] = {str(N): z for N, z in TAF.items()}
ZAHLEN["A2"] = {"anomaliefrei": frei, "positive_Eigenwerte_i^j": {str(N): v for N, v in positiv.items()}, "G_Test": gtest,
                "Vorzeichen_stabil_34_gegen_50_Stellen": stabil, "rel_Abweichung_34_gegen_50": relabw, "reell": reell}
normal = [N for N in frei if 2 in positiv[N]]
check("A2 [N, entscheidend] Nahtzahl: Für N = 1..12 Paare Naht + Partner ist die diagonal summierte Verklebung X genau für "
      f"gerades N anomaliefrei ({frei}). Der Eigenwert -1 (R = Z_SU^-1, normale Statistik) gibt positive Spinmoden nur für "
      f"N = {normal}. Vertauschte Statistik besteht für N = 2 (-i), N = 4 (+i) und N = 6 (-i); für N = 8, 10, 12 besteht kein "
      "Eigenwert, und X = +1 (Sektor des Skalarvakuums) besteht nie. Die Vorzeichen sind bei 34 und 50 Stellen dieselben "
      f"(relative Abweichung der Minima {relabw:.0e}); alle Summen sind reell ({reell:.0e})",
      frei == [2, 4, 6, 8, 10, 12] and normal == [2] and positiv == {2: [2, 3], 4: [1], 6: [3], 8: [], 10: [], 12: []}
      and stabil and reell < 1e-12)
say(f"   {zeit()}")


def zyklisch(name):
    """Zyklische Verklebung zweier Kopien mit Einschub g = sigma (X (x) 1): Tr g q^(L_0) = Z_X(2 tau)."""
    F = naht(NAMEN[name])
    Bz = bahn(lambda tau: F(2*tau), maxn=3000)
    _, bad = konsistenz(Bz, NULLCH)
    ch = [(str(c["T"]), str(c["S"])) for c, _ in charaktere(Bz)]
    return Bz, sorted(str(x) for x in bad), ch


ZY = {k: zyklisch(k) for k in B}
M2 = {k: moden(B[k], GW[k], 2) for k in B}
m0 = M2["1"]
Z1c = moden(ZY["1"][0], gitter(ZY["1"][0]))
ZGc = moden(ZY["G"][0], gitter(ZY["G"][0]))
sym = {"Sym^2 (sigma = +1)": minimum(comb([m0, Z1c], [1, 1]), m0), "sigma = -1": minimum(comb([m0, Z1c], [1, -1]), m0)}
sg = {"sigma G = +1": minimum(comb([m0, M2["G"], ZGc], [1, 1, 2]), m0),
      "sigma G = -1": minimum(comb([m0, M2["G"], ZGc], [1, 1, -2]), m0),
      "sigma G = +-i": minimum(comb([m0, M2["G"]], [1, -1]), m0)}
ZAHLEN["A3"] = {"Bahnen": {k: v[0]["n"] for k, v in ZY.items()}, "Widersprüche trivialer Charakter": {k: v[1] for k, v in ZY.items()},
                "Charaktere": {k: v[2] for k, v in ZY.items()}, "ohne Uhr": sym, "mit G": sg}
check("A3 [E/N, entscheidend] Zyklische Verklebung zweier Kopien (eine Spur durch beide, Öffnungswinkel 4 pi; in "
      "arXiv:2611.34924 Tr[(C^+ C)^2] gegenüber [Tr C^+ C]^2): Die Spur mit sigma (x) Einschub ist Z(2 tau). Kontrolle ohne "
      f"Uhr: Bahn {ZY['1'][0]['n']}, widerspruchsfrei, positiv (symmetrisches Quadrat). Mit G = U R^2: Bahn "
      f"{ZY['G'][0]['n']}, widerspruchsfrei, positiv nur für sigma G = -1. Mit der Z4-Verklebung X: Bahn {ZY['X'][0]['n']}, "
      "kein widerspruchsfreier Charakter. Die Verklebung, die Chiralität und Spintabelle trägt, verträgt die zyklische "
      "Verklebung nicht; sie braucht zwei getrennt geschlossene Nähte (diagonale Summe)",
      ZY["1"][0]["n"] == 3 and ZY["1"][1] == [] and sym["Sym^2 (sigma = +1)"][0] > EPS
      and ZY["G"][1] == [] and sg["sigma G = -1"][0] > EPS and sg["sigma G = +1"][0] < EPS and sg["sigma G = +-i"][0] < EPS
      and ZY["X"][1] != [] and ZY["X"][2] == [])
befund("Teil A: Normale Statistik wählt im Modultest genau zwei Paare (N = 2 unter N = 1..12). Die Zahl der Nähte ist damit "
       "keine eigene Wahl mehr; in P_paar bleibt die Aussage, dass die Verklebung einseitig ist und über getrennte Nähte "
       "summiert wird. Die zyklische Verklebung (Replik, Öffnungswinkel 4 pi) scheidet als Herkunft des Paars aus", True)
say(f"   {zeit()}")

# =================================================================== Teil B
say("\nTeil B: Standardform des Nahtpaars im endlichen Modell")
nsA = {"np": np, "expm": expm}
for fn in ("jw", "nrm", "modell", "gatter", "linearitaet"):
    exec(re.search(rf"^def {fn}\(.*?(?=^\S)", SRC174, re.S | re.M).group(0), nsA)
modell, gatter, linearitaet, nrm = (nsA[k] for k in ("modell", "gatter", "linearitaet", "nrm"))
THS, PHS = np.linspace(0, np.pi, 13), np.linspace(0, 2*np.pi, 25)


def gewichte(K, art):
    """Diagonale von rho einer Naht (K Moden): entartet (K = pi Q), modenaufgelöst (beta = 4 pi, Energien c + 1/4) oder
    ladungsquadratisch (K = 2 pi Q^2 - pi Q ohne Oszillatoren, die Karikatur aus z174)."""
    d = 2**K
    bits = np.array([[(i >> (K - 1 - c)) & 1 for c in range(K)] for i in range(d)])
    q = bits.sum(axis=1)
    if art == "entartet":
        e = np.pi*q
    elif art == "modenaufgelöst":
        e = 4*np.pi*(bits*(np.arange(K) + 0.25)).sum(axis=1)
    else:
        e = 2*np.pi*q**2 - np.pi*q
    p = np.exp(-(e - e.min()))
    return p/p.sum()


rng = np.random.default_rng(177)
b1, b2 = {}, {}
for K in (2, 3):
    m = modell(K, lambda q: np.pi*q)                       # Seiten A, B sind hier zwei Nähte; Hüpfen Mode c <-> Mode c
    d = m["d"]
    PA = np.diag((-1.0)**(m["QA"] % 2))
    U_all = {(it, ip): gatter(m, th, ph) for it, th in enumerate(THS) for ip, ph in enumerate(PHS)}
    for art in ("entartet", "modenaufgelöst", "ladungsquadratisch"):
        p = gewichte(K, art)
        w = np.sqrt(np.kron(p, p))                         # Omega = vec(rho^(1/2)), rho = rho_A (x) rho_B
        fix = sorted(k for k, U in U_all.items() if nrm(U @ np.diag(w) @ U.conj().T - np.diag(w)) < 1e-10)
        b1[(K, art)] = fix
    Ui = gatter(m, np.pi/2, np.pi/2)                       # ein wechselwirkendes Mitglied (z144: k = 1)
    M_ = rng.normal(size=(d*d, d*d)) + 1j*rng.normal(size=(d*d, d*d))
    P_ = M_ @ M_.conj().T                                  # positives Element des Kegels vec(P)
    kegel = float(np.linalg.eigvalsh(Ui @ P_ @ Ui.conj().T).min())
    jrest = nrm((Ui @ M_ @ Ui.conj().T).conj().T - Ui @ M_.conj().T @ Ui.conj().T)   # W J = J W, J vec(M) = vec(M^+)
    lin = {kk: round(linearitaet(m, gatter(m, np.pi/2, kk*np.pi/4)), 6) for kk in range(8)}
    par = sorted({it for (it, ip), U in U_all.items() if ip == 6 and nrm(U @ PA - PA @ U) < 1e-10})
    b2[K] = {"Kegel_min_Eigenwert": kegel, "J_Rest": jrest, "Linearitätsrest (phi = k pi/4)": lin,
             "theta-Index mit Paritätserhaltung je Naht (phi = pi/2)": par}
alle = len(THS)*len(PHS)
zahl = {f"K={K}, {art}": len(v) for (K, art), v in b1.items()}
quad3 = b1[(3, "ladungsquadratisch")]
quad3_rest = [k for k in quad3 if k[0] not in (0, 12) and k[1] not in (0, 24)]
ZAHLEN["B1_feste_Gatter (von 325)"] = zahl
ZAHLEN["B1_ladungsquadratisch_K3_Rest (theta-Index, phi-Index)"] = quad3_rest
ZAHLEN["B2"] = {str(K): v for K, v in b2.items()}
check("B1 [E] Kanonische Implementierung W = U (x) J U J des Gatters zwischen zwei Nähten in der Standardform des Paars "
      "(Omega = vec(rho^(1/2)), rho = rho_A (x) rho_B; K = 2, 3 Moden je Naht, 13 Winkel, 25 Phasen): W erhält den Kegel "
      "und vertauscht mit J; W Omega = Omega genau, wenn U den Zustand erhält. Im Zustand freier Quarter (entartet und "
      f"modenaufgelöst mit Energien c + 1/4 bei beta = 4 pi) gilt das für alle {alle} Gatter. In der ladungsquadratischen "
      f"Karikatur ohne Oszillatoren bleiben bei K = 3 nur theta in {{0, pi}}, phi in {{0, 2 pi}} und {quad3_rest} "
      "(theta = pi/2, phi = pi), alles einkörperig oder nahtdiagonal",
      all(len(b1[(K, art)]) == alle for K in (2, 3) for art in ("entartet", "modenaufgelöst"))
      and quad3_rest == [(6, 12)] and all(v["Kegel_min_Eigenwert"] > 0 and v["J_Rest"] < 1e-9 for v in b2.values()))
check("B2 [E] Wechselwirkung: Die Mitglieder mit phi = +-pi/2 bei theta = pi/2 sind nach dem Kriterium aus z144 nicht "
      "einkörperig (Linearitätsrest > 0) und sind im Zustand freier Quarter Symmetrien der Standardform des Paars. Das "
      "Paar trägt also die Wechselwirkung, die eine Naht mit Partner nicht trägt (z174 A7). Es wählt aber nichts aus: Winkel "
      "und Phase bleiben frei. Die Fermionparität je Naht erhält das Gatter nur für theta in {0, pi/2, pi} (z125: P_univ "
      "erzwingt die Vierteldrehung); die Standardform des Paars verlangt das nicht, P_univ ist in der Paarlesart nicht "
      "nativ",
      all(v["Linearitätsrest (phi = k pi/4)"][2] > 0.5 and v["Linearitätsrest (phi = k pi/4)"][6] > 0.5
          and v["Linearitätsrest (phi = k pi/4)"][0] == 0 and v["Linearitätsrest (phi = k pi/4)"][4] == 0
          and v["theta-Index mit Paritätserhaltung je Naht (phi = pi/2)"] == [0, 6, 12] for v in b2.values()))
befund("Teil B: Die zwei Lesarten leisten Verschiedenes. Partnerlesart (z174): wählt die Phase, trägt keine "
       "Wechselwirkung. Paarlesart (hier): trägt die Wechselwirkung, wählt weder Winkel noch Phase. Keine der beiden "
       "schließt die Zusammensetzung allein; P_std und P_ord bleiben", True)
say(f"   {zeit()}")

# =================================================================== Teil C
say("\nTeil C: arXiv:2611.34924 gegen den Stand nach z174 (gelesen: Abschnitte 1 bis 6, Anhänge D und E)")
c1 = {}
c1["(5.29)"] = float(max(abs(mp.e**eta_*(2**2)/(1 + mp.e**eta_)**2 - (1 - mp.tanh(eta_/2)**2))
                         for eta_ in (mp.mpf('0.3'), mp.mpf(1), mp.mpf('2.7'), mp.mpf(-1))))
c1["(5.76)"] = float(max(abs(mp.quad(lambda v: v*mp.cos(v*mp.log(uL/uR))/mp.sinh(PI*v), [0, 1, 4, 12, 40])
                             - uL*uR/(uL + uR)**2)
                         for uL, uR in ((mp.mpf(1), mp.mpf(1)), (mp.mpf(1), mp.mpf(3)), (mp.mpf('0.2'), mp.mpf(5)))))
c1["(5.75)"] = float(max(abs((1 - mp.exp(-2*bw))*mp.nsum(lambda n: n*mp.exp(-(2*n - 1)*bw), [1, mp.inf]) - 1/(2*mp.sinh(bw)))
                         for bw in (mp.mpf('0.4'), mp.mpf(1), mp.mpf('2.5'))))
cs = []
for _ in range(200):                                        # (D.49): Z_2,zyk/Z_1^2 >= 1/sum_s e^s d(s)
    s_ = rng.uniform(0, 6, size=5)
    k_ = [rng.uniform(0, 1, size=rng.integers(1, 6)) for _ in s_]
    z1 = sum(k.sum() for k in k_)
    z2 = sum(np.exp(-s)*(k**2).sum() for s, k in zip(s_, k_))
    cs.append(z2/z1**2 - 1/sum(np.exp(s)*len(k) for s, k in zip(s_, k_)))
c1["(D.49) kleinster Überschuss"] = float(min(cs))
ZAHLEN["C1_Formeln (Abweichung)"] = c1
check("C1 [E/N] Nachgerechnete Formeln der Arbeit: (5.29) e^eta G(-d, e^eta d)/G(-d, d) = 4 e^eta/(1 + e^eta)^2 = "
      "1 - tanh^2(eta/2); (5.76) int_0^inf nu cos(nu log(u_L/u_R))/sinh(pi nu) d nu = u_L u_R/(u_L + u_R)^2; (5.75) "
      "<a_L a_R> = 1/(2 sinh(b omega)) im Zustand sum_n e^(-n b omega) |n>|n*>; (D.49) Z_2,zyk/Z_1^2 >= 1/sum_s e^s d(s) "
      "(200 Zufallsbelegungen)",
      c1["(5.29)"] < 1e-40 and c1["(5.76)"] < 1e-12 and c1["(5.75)"] < 1e-30 and c1["(D.49) kleinster Überschuss"] >= 0)

c2 = {}
for K in (2, 3):
    m = modell(K, lambda q: np.pi*q)                       # jetzt wieder: Naht (A) und Tomita-Partner (B)
    d = m["d"]
    rho = np.diag(m["p"])
    x = rng.normal(size=(d, d)) + 1j*rng.normal(size=(d, d))
    erw = abs(m["Om"].conj() @ np.kron(x, np.eye(d)) @ m["Om"] - np.trace(rho @ x))
    s_op = -np.log(m["p"]) + np.log(m["p"].max())           # s = -log rho bis auf eine Konstante: pi Q
    qa = np.array([bin(i).count("1") for i in range(d)])
    boost = lambda G, s=s_op: np.exp(1j*G*s/(2*np.pi))      # die Arbeit: Boost der Rapidität Gamma wirkt als e^(i Gamma s/(2 pi))
    c2[K] = {"<Omega|x|Omega> - Tr(rho x)": float(erw), "s - pi Q": float(np.abs(s_op - np.pi*qa).max()),
             "W1: Boost pi = R": float(np.abs(boost(np.pi) - 1j**qa).max()),
             "W1: Boost 2 pi = R^2": float(np.abs(boost(2*np.pi) - (-1.0)**qa).max()),
             "W1: Boost 4 pi = 1": float(np.abs(boost(4*np.pi) - 1).max()),
             "W2: Boost 2 pi = R": float(np.abs(boost(2*np.pi, s_op/2) - 1j**qa).max()),
             "Gewicht von C: e^(-s/2) = Omega": float(np.abs(np.exp(-s_op/2)/np.linalg.norm(np.exp(-s_op/2))
                                                                - np.sqrt(m["p"])).max())}
ZAHLEN["C2_Wörterbuch"] = {str(K): v for K, v in c2.items()}
check("C2 [E] Wörterbuch im endlichen Modell (K = 2, 3). Die Arbeit [Lit]: Präparation C = U e^(-A/(8 G hbar)) P_D mit "
      "C^+ C = e^(-s), s = A/(4 G hbar); zwei Ränder tragen H (x) H*; ein Boost der Rapidität Gamma wirkt als "
      "e^(i Gamma s/(2 pi)). Die Naht [E]: Omega = vec(rho^(1/2)) hat die Gewichte e^(-s/2) von C mit s = pi Q, "
      "<Omega|x (x) 1|Omega> = Tr(rho x), und H (x) H* mit J ist die Standardform aus z130/z174. Wörterbuch W1 (rho der "
      "Arbeit = rho der Naht): Die Uhr R = i^Q ist der Boost um pi, R^2 = (-1)^F der um 2 pi, die Periode ist 4 pi. "
      "Wörterbuch W2 (rho der Naht = rho^2 der Arbeit, die Naht als zwei zyklisch verklebte Kopien): Dann wäre R der "
      "Boost um 2 pi",
      all(max(v.values()) < 1e-12 for v in c2.values()))

zeilen = Z174["zahlen"]["B6_spintabelle (2h -> Uhrexponent T, R = i^T)"]["X' = -1"]
tab = {int(k): v[0] for k, v in zeilen.items()}             # 2h -> T
w1 = {zh: Fr(T, 2) for zh, T in tab.items()}                # R = e^(i pi K): K = T/2 mod 2
w2 = {zh: Fr(T, 4) for zh, T in tab.items()}                # R = e^(2 pi i K): K = T/4 mod 1
w1_ok = all((2*K).denominator == 1 for K in w1.values())
w2_ok = all((2*K).denominator == 1 for K in w2.values())
g_ok = all((T - zh) % 2 == 0 for zh, T in tab.items())      # e^(2 pi i K) = (-1)^T gleich U = (-1)^(2h)
ZAHLEN["C3"] = {"W1: 2h -> Boostgewicht K mod 2": {str(k): str(v) for k, v in w1.items()},
                "W2: 2h -> K mod 1": {str(k): str(v) for k, v in w2.items()}}
check("C3 [E für die Arithmetik, L für die Lesart] Abgleich 2 pi gegen 4 pi (in z166 offen): Mit der Spintabelle der "
      "Verklebung (z174 B6, R = i^(2-2h)) gibt W1 die Boostgewichte K = 1 - h mod 2, also nur halbzahlige Werte (die 16: "
      "K = 1/2, ein Spinor), und G = U R^2 = +1 heißt dort: Die modulare Drehung um 2 pi ist die geometrische. W2 gibt für "
      "die 16 das Gewicht 1/4, das keine Darstellung der Lorentz-Gruppe hat. Zusammen mit A3 (die zyklische Verklebung ist "
      "anomal) scheidet W2 aus. Der KMS-Wert 4 pi des Viertelprozesses ist dann der Wert 2 pi der Arbeit in der "
      "Normierung, in der ein Quarter die Energie 1/4 und das Boostgewicht 1/2 hat (K = 2 H): der Parameter der Naht ist "
      "der doppelte Winkel, 4 pi ist die Periode der Spinoren",
      w1_ok and not w2_ok and g_ok and w1[1] == Fr(1, 2) and w2[1] == Fr(1, 4))

erw_tab = {zh: (2 - zh) % 4 for zh in range(-4, 5)}         # X' = i^T i^(2h) = -1
neutral = sorted(zh for zh, T in erw_tab.items() if T == 0)
np_daten = {"d_u phi": (Fr(1), 0), "psi": (Fr(1, 2), 1), "A_A": (Fr(0), 2), "sigma_AB": (Fr(1), 4)}      # (K, 2h)
np_gegen = {"phi": (Fr(0), 0), "d_u A_A": (Fr(1), 2), "h_AB": (Fr(0), 4)}
np_ok = all((2*K) % 4 == erw_tab[zh] for K, zh in np_daten.values())
np_nein = all((2*K) % 4 != erw_tab[zh] for K, zh in np_gegen.values())
produkt = {(a, b): ((erw_tab[a] + a) + (erw_tab[b] + b)) % 4 for a in (2, -2) for b in (2, -2)}          # Exponent von X'_1 X'_2
ZAHLEN["C4"] = {"2h -> T": {str(k): v for k, v in erw_tab.items()}, "neutral (T = 0) bei 2h": neutral,
                "Produkt zweier erlaubter Objekte (Exponent von X')": sorted(set(produkt.values()))}
check("C4 [E für die Arithmetik, L für die Lesart] Erweiterte Spintabelle (|h| <= 2): R = 1 (uhrneutral) gibt es im "
      f"verklebten Einzelobjekt nur bei ungeradem h (2h in {neutral}); h = +-2 verlangt wie h = 0 den Uhrwert R = -1. "
      "Ein uhrneutraler Spin-2-Zustand liegt nicht im Einzelobjekt. Die Regel X' = -1 ist nicht multiplikativ: Zwei "
      "erlaubte Vektoren (h = 1, R = 1) zusammen haben X' = +1, also h = 2 und R = 1 nur als Zustand zweier Objekte. Das "
      "passt zum Weg 'Doppelkopie' (P_diag, z112) und nicht zu einem Graviton als Anregung einer Naht. In W1 erfüllen "
      "genau die Nullkegeldaten d_u phi, psi, A_A und die Scherung sigma_AB (Boostgewicht 1, in der Arbeit eta-Faktor von "
      "(5.29)) die Regel e^(i pi (K + h)) = -1; phi, d_u A_A und h_AB erfüllen sie nicht. Die Scherung hat dann R = -1 "
      "als Boostgewicht, ist aber eichneutral: W1 und P_lift (Uhr = Zentrum von Spin(10)) zusammen lassen eine "
      "fundamentale Scherung nicht zu",
      neutral == [-2, 2] and erw_tab[4] == 2 and erw_tab[-4] == 2 and erw_tab[0] == 2
      and set(produkt.values()) == {0} and np_ok and np_nein)
befund("Teil C [Lit]: Die Arbeit hält G und Lambda als Eingaben, zeigt die Diskretheit der Fläche nur als Existenzaussage "
       "(endliche Dimension, keine Eigenwerte) und enthält keine Materie, keine Spinoren und keine Z4-Graduierung. Sie "
       "schließt P_G, P_curv und die Graviton-Klausel nicht", False)
befund("Teil C: Die Arbeit hilft an drei Stellen. (1) Sie ist eine unabhängige Konstruktion, in der die Standardform mit "
       "Partner (H (x) H*, Präparation mit halber Drehung) der Hilbertraum zweier Ränder ist. (2) Der Abgleich 2 pi gegen "
       "4 pi aus z166 ist entschieden (W1). (3) Ihr Unterschied zwischen getrennt und zyklisch geschlossenen Kopien ist "
       "genau der Test, der das Nahtpaar als getrennt summiert festlegt (A3)", True)
say(f"   {zeit()}")

# =================================================================== Teil D
say("\nTeil D: Kandidatentabelle, Suchprogramm und Abstand zur Komplettlösung")
(HERE / "z177_ergebnis.json").write_text(json.dumps(                            # die Belege, auf die das Suchprogramm zugreift
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": "vorläufig (Teil D läuft)"},
    indent=1, ensure_ascii=False, default=str))
import reconstruction_search as rs                                              # noqa: E402

uni = rs.unification_candidates()
dist = rs.complete_solution_distance()
objs, meta = rs.scenario_candidates()
for r in uni["candidates"]:
    say(f"      {r['verdict']:9s} {r['candidate']}  (würde schließen: {', '.join(r['would_close'])}; {r['source'][:4]})")
ZAHLEN["D_kandidaten"] = [{k: r[k] for k in ("candidate", "would_close", "source", "verdict")} for r in uni["candidates"]]
ZAHLEN["D_abstand"] = {"distance": dist["distance"], "open_total": dist["open_total"], "by_question": dist["by_question"],
                       "open_items": sorted(dist["open_items"]), "open_under_readings": uni["open_under_readings"]}
check("D1 [E] Kandidatentabelle im Suchprogramm (OriginLab, unification_candidates): elf Kandidaten für ein Objekt hinter "
      f"mehreren offenen Punkten, jeder mit einer gerechneten Prüfung aus z175, z176 oder z177. Urteile: {uni['counts']}. "
      "Kein Kandidat ist eine Herleitung; keiner ändert den Abstand",
      uni["counts"] == {"excluded": 8, "holds": 2, "reading": 1} and uni["closed_by_derivation"] == [])
impl = {(i["replaces"], tuple(i["by"])): i["holds"] for i in meta["implications"]}
check("D2 [E] Zuordnung im Suchprogramm: Die Implikation P_lift -> P_zentral (z175: Ordnung 4 und Vertauschen mit den "
      "Basiswechseln beider Deckhälften gibt die Uhr) ist eingetragen und gilt; P_zentral ist als Annahme benannt",
      impl.get(("P_lift", ("P_zentral",))) and "P_zentral" in rs.PREMISE_MEANING and "P_zentral" in dist["open_items"]
      and "P_lift" not in dist["open_items"])
check("D3 [E] Abstand zur Komplettlösung nach den drei Rechnungen: unverändert 1 fehlende Klausel (graviton), 1 freie "
      f"Zahl (P_G), 7 Annahmen ({sorted(x for x in dist['open_items'] if x.startswith(('P_', 'I_')))}), dazu K19 und K20: "
      f"zehn offene Punkte. Unter der Lesart 'Superauswahl der Paaruhr = Eichsumme' (z176) blieben "
      f"{uni['open_total_under_readings']}; die Lesart ist nicht hergeleitet und zählt nicht",
      dist["distance"] == {"missing_clauses": 1, "free_numbers": 1, "premises_beyond_axiom": 7}
      and dist["open_total"] == 10 and not dist["complete"] and uni["open_total_under_readings"] == 8
      and sorted(dist["open_items"]) == ["I_step", "P_G", "P_curv", "P_ord", "P_paar", "P_std", "P_zentral", "graviton"])
befund("Teil D: Keine Komplettlösung. Die drei Rechnungen schließen keinen Punkt; sie schließen sieben Wege mit "
       "Geltungsbereich aus und fassen drei Punkte enger (P_paar, P_lift -> P_zentral, K19/K20 unter einer Lesart)", True)

urteil = (
    "Drei Rechnungen parallel (z175 Zuordnung, z176 Ursprung, z177 Kandidaten und die Arbeit arXiv:2611.34924). Kein offener "
    "Punkt ist geschlossen; der Abstand bleibt bei zehn. Was sich ändert: (1) Das Nahtpaar. Normale Statistik besteht den "
    "Modultest genau für zwei Paare (N = 1..12), also ist die Zahl zwei keine Wahl mehr. Das Paar ist weder die "
    "Blattverdopplung von P1 (z176) noch eine zyklische Verklebung zweier Kopien (anomal); in P_paar bleibt, dass die "
    "Verklebung einseitig ist und über getrennt geschlossene Nähte summiert wird. (2) Zusammensetzung. Die Standardform des "
    "Paars trägt das wechselwirkende Gatter als Symmetrie, wählt aber weder Winkel noch Phase; die Partnerlesart wählt die "
    "Phase und trägt keine Wechselwirkung. Keine Lesart schließt die Zusammensetzung allein. (3) Die Arbeit. Sie liefert "
    "unabhängig die Standardform mit Partner als Hilbertraum zweier Ränder und entscheidet den Abgleich 2 pi gegen 4 pi: "
    "Die Uhr ist der Boost um pi, 4 pi die Periode der Spinoren. G und Lambda bleiben dort Eingaben, die Fläche ist nur als "
    "diskret erwiesen; P_G, P_curv und die Graviton-Klausel schließt sie nicht. (4) Graviton. Die Spintabelle der Verklebung "
    "hat uhrneutrale Zustände nur bei ungeradem h; ein neutraler Spin-2-Zustand existiert nur als Zustand zweier Objekte "
    "(Doppelkopie). (5) Zuordnung (z175): P_lift wird zu P_zentral, Sp R = 0 folgt nicht aus <T> = 0, I_step bleibt. "
    "(6) Ursprung (z176): K19 und K20 folgen aus P_std und A* nur unter der Lesart 'Superauswahl = Eichsumme'. Keine "
    "Komplettlösung.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z177_ergebnis.json").write_text(json.dumps({"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil},
                                                    indent=1, ensure_ascii=False, default=str))
(HERE / "ausgabe_z177.txt").write_text("\n".join(x for x in OUT if "Laufzeit bisher" not in x) + "\n")
say(f"{zeit()}")
