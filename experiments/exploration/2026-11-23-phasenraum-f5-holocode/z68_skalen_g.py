"""Z68: G und Skalen aus den Antworten des Originals (Lücke F der Vogelperspektive).

Frage: G ist bisher nur eine Einheit (8 pi G = 1/(c3 v_geo^2), c3 = 1/(8 pi), also G = v_geo^-2). Dimensionslose
Verhältnisse (a/l_P, M_Pl/m_e) brauchen eine skalentragende Antwort; Sakharovs induzierte Gravitation hängt am Regler
(Z27). Unter der Leitlinie "Rekonstruktion statt Auswahl" liefert der Originalzustand sein eigenes Spektrum, also seinen
eigenen Regler. Legt das eine dimensionslose Zahl fest?

Normierungen, vorab und vor jedem Vergleich festgelegt:
- Original: c(x) = 1/(1 + e^{4 pi x}), x = r + alpha, alpha = 1/4, 3/4 (z46, z58). Einteilchengenerator h = 4 pi x;
  voller Modulgenerator K = 4 pi L_0, L_0 = sum J_-m J_m + Q^2/2 + (alpha - 1/2) Q, Spektrum in pi N (z58).
- Physikalische Einheit: KMS bei beta_phys = L (z58), also H = K/L und Energien E = (K-Eigenwert)/L.
- Sakharov je skalarer Einheit (minimal gekoppelt, a_1 = R/6): 1/(16 pi G) = (1/(192 pi^2)) Int dt t^-2 rho(t)
  Tr e^{-t M^2}; Dirac = 2 Einheiten, Weyl = 1, konform gekoppelter Skalar = 0; l_P = sqrt(G) (hbar = c = 1).
- Spektralmoment: F_2 = Int_0^inf F(u) du, u = D^2/Lambda^2, Lambda = 1/L, |1/(16 pi G)| = F_2 Lambda^2/(96 pi^2) je Einheit.
- 2D: Gamma = (c/96 pi) Int sqrt(g) R Box^-1 R (Weyl-Anomalie -c R/(24 pi)); Spurkanal kappa_2 = beta^2 d^2_beta log Z.

Teil 1  Sakharov-Summen über die rekonstruierten Niveaus (eigener Regler statt äußerer Abschneidung).
Teil 2  Spurkanal kappa_2 und Tensorkanal c als Koeffizienten der induzierten 2D/4D-Wirkung.
Teil 3  beta = 4 pi, q = e^{-2 pi} und die Lesart Nahtskala ~ Planck (GRAV.GAPCORR.01).
Teil 4  Welche dimensionslose Zahl wird vorhergesagt? Bilanz.
"""
import json
import math
from itertools import accumulate
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.integrate import quad

HERE = Path(__file__).parent
mp.mp.dps = 40
OUT = []
CHECKS, BEFUNDE, ZAHLEN, STATUS = {}, {}, {}, {}


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val):
    BEFUNDE[name] = bool(val)
    say(f"   [Befund: {'ja' if val else 'nein'}] {name}")


PI = mp.pi
ALPHAS = (mp.mpf(1)/4, mp.mpf(3)/4)
Z58 = json.loads((HERE / "z58_ergebnis.json").read_text())
Z59 = json.loads((HERE / "z59_ergebnis.json").read_text())
C_REC = mp.mpf(Z59["zahlen"]["alpha=0.25"]["c"])
H_REC = mp.mpf(Z59["zahlen"]["alpha=0.25"]["h"])
KAPPA2_Z58 = mp.mpf(Z58["zahlen"]["energiekumulanten_K"]["kappa2"])
ZETA3 = mp.zeta(3)


def logZq(beta, a, QR=60, MR=None):
    """log Tr e^{-beta L_0} je Quarter (Thetasumme der Ladung mal Oszillatoren), L_0-Einheiten."""
    if MR is None:
        MR = int(120/float(beta)) + 40
    th = mp.fsum(mp.exp(-beta*(Qn**2/2 + (a - mp.mpf(1)/2)*Qn)) for Qn in range(-QR, QR + 1))
    osc = -mp.fsum(mp.log(1 - mp.exp(-beta*m)) for m in range(1, MR))
    return mp.log(th) + osc


# ================================================================== Teil 1
say("Teil 1: Sakharov-Summen über die rekonstruierten Niveaus (eigener Regler statt äußerer Abschneidung)")
say("  1a Einteilchenniveaus h = 4 pi x, x in Z + alpha (z58)")


def Zq(s, a):
    return (4*PI)**(-s)*(mp.zeta(s, a) + mp.zeta(s, 1 - a))


d_ = mp.mpf(10)**-15
res = [d_*Zq(1 + d_, a) for a in ALPHAS]
check("Zetafunktion der Einteilchenniveaus: einziger Pol bei s = 1 mit Residuum 1/(2 pi) je Quarter (Weyl-Gesetz, "
      "Spektraldimension 1 wie z60); bei s = 0, -1, -2 endlich",
      all(abs(r - 1/(2*PI)) < 1e-12 for r in res) and all(mp.isfinite(Zq(s, a)) for s in (0, -1, -2) for a in ALPHAS))
z0 = [Zq(0, a) for a in ALPHAS]
check("Z_alpha(0) = 0 je Quarter: die zetaregularisierte Determinante der Einteilchenniveaus hängt nicht von der Einheit L "
      "ab (keine Skalenanomalie, log det -> log det - Z(0) log L); der Stromturm 4 pi m hat dagegen zeta(0) = -1/2",
      all(abs(z) < mp.mpf(10)**-35 for z in z0) and abs(mp.zeta(0) + mp.mpf(1)/2) < mp.mpf(10)**-35)
evac = [-mp.zeta(-1, 1 - a) for a in ALPHAS]
e8vac = 8*mp.zeta(-1)/2
check("Casimirenergie (zetaregularisierte Dirac-See-Summe) -zeta(-1, 1 - alpha) = -1/96 je Quarter = h - c/24 mit c, h "
      "aus z59; E8-Naht 8 zeta(-1)/2 = -1/3 = -c/24: reglerfrei und gleich den rekonstruierten Virasoro-Daten",
      all(abs(e - (H_REC - C_REC/24)) < 1e-15 for e in evac) and abs(evac[0] + mp.mpf(1)/96) < 1e-30
      and abs(e8vac + mp.mpf(1)/3) < 1e-30)
ZAHLEN["einteilchen"] = {"residuum_s1": mp.nstr(res[0], 15), "Z(0)": [mp.nstr(z, 5) for z in z0],
                         "casimir_L0_je_quarter": mp.nstr(evac[0], 15), "casimir_phys": "-pi/(24 L) je Quarter"}
befund("Die reglerfreien Summen über die Einteilchenniveaus sind endlich, aber nur Zahlen mal Potenzen der einen Länge L "
       "(Casimirenergie -pi/(24 L) je Quarter); eine zweite Skale entsteht nicht", True)

say("  1b Voller Modulgenerator K = 4 pi L_0 (Niveaus pi n, Entartungen aus Theta x Oszillatoren)")
NMAX = 12000


def partitions(n):
    p = [1] + [0]*n
    for k in range(1, n + 1):
        s, j = 0, 1
        while True:
            g1 = j*(3*j - 1)//2
            if g1 > k:
                break
            sg = 1 if j % 2 else -1
            s += sg*p[k - g1]
            g2 = j*(3*j + 1)//2
            if g2 <= k:
                s += sg*p[k - g2]
            j += 1
        p[k] = s
    return p


PN = partitions(NMAX//4 + 1)


def degeneracies(a4):            # a4 = 4 alpha - 2: Niveau n = 2 Q^2 + a4 Q + 4 N (Einheit pi)
    d = [0]*(NMAX + 1)
    for Q in range(-100, 101):
        e = 2*Q*Q + a4*Q
        if 0 <= e <= NMAX:
            for N in range((NMAX - e)//4 + 1):
                d[e + 4*N] += PN[N]
    return d


DEG = degeneracies(-1)                                      # alpha = 1/4 (alpha = 3/4 ist das Spiegelbild Q -> -Q)
CUM = list(accumulate(DEG))
n1, n2 = 3000, 12000
slope = (math.log(CUM[n2]) - math.log(CUM[n1]))/(math.sqrt(n2) - math.sqrt(n1))


def loc_exp(n):
    return (math.log(CUM[int(1.1*n)]) - math.log(CUM[n]))/math.log(1.1)


ex1, ex2 = loc_exp(1000), loc_exp(10000)
ZAHLEN["entartung"] = {"steigung_ln_N_gegen_sqrt_n": slope, "cardy_pi/sqrt6": float(PI/mp.sqrt(6)),
                       "lokaler_exponent_n1000": ex1, "lokaler_exponent_n10000": ex2}
say(f"   ln N(n) gegen sqrt(n): Steigung {slope:.4f} (Cardy c = 1: pi/sqrt 6 = {float(PI/mp.sqrt(6)):.4f}); lokaler "
    f"Exponent d ln N/d ln n = {ex1:.1f} (n = 1000), {ex2:.1f} (n = 10000)")
check("Entartungen des vollen Modulspektrums wachsen wie exp(pi sqrt(n/6)) (Cardy, c = 1; Steigung auf 3 %), der lokale "
      "Exponent wächst wie sqrt(n) (Verhältnis sqrt 10 auf 10 %): Sum d_n (pi n)^-s divergiert für jedes s, die "
      "Zetafunktion des vollen Spektrums existiert nicht",
      abs(slope/float(PI/mp.sqrt(6)) - 1) < 0.03 and abs((ex2/ex1)/math.sqrt(10) - 1) < 0.10)
Zprod = mp.exp(logZq(4*PI, ALPHAS[0]))
Zlev = mp.fsum(DEG[n]*mp.exp(-PI*n) for n in range(0, 401))


def W(J):
    return mp.fsum(logZq(4*PI*(1 - J), a) for a in ALPHAS)


k2 = mp.diff(W, 0, 2)
Kmean = [mp.diff(lambda J: logZq(4*PI*(1 - J), a), 0, 1) for a in ALPHAS]
Sq = [logZq(4*PI, a) + km for a, km in zip(ALPHAS, Kmean)]
S_seam = mp.fsum(Sq)
check("am KMS-Punkt endlich: Tr e^{-K} aus der Niveausumme Sum d_n e^{-pi n} = Produktformel (relativ 1e-30); "
      "kappa_2 beider Quarter = z58 (0,79821)",
      abs(Zlev/Zprod - 1) < mp.mpf(10)**-30 and abs(k2 - KAPPA2_Z58) < 1e-9)
ZAHLEN["kms_punkt"] = {"Z_je_quarter": mp.nstr(Zprod, 15), "kappa2": mp.nstr(k2, 12),
                       "S_je_quarter": mp.nstr(Sq[0], 12), "S_naht": mp.nstr(S_seam, 12)}
say(f"   KMS-Punkt: Z = {mp.nstr(Zprod, 10)} je Quarter, kappa_2 = {mp.nstr(k2, 10)}, Entropie S = {mp.nstr(Sq[0], 8)} je "
    f"Quarter, {mp.nstr(S_seam, 8)} je Naht")
befund("Eine Sakharov-Summe über alle rekonstruierten Niveaus ist nur am KMS-Punkt endlich (dort ist e^{-K} spurklasse); "
       "dort liefert sie dimensionslose Zahlen (log Z, S, kappa_n), keinen Koeffizienten der Dimension Masse^2", True)

say("  1c Lesart KK: Einteilchenniveaus als 4D-Massen m = |h|/L (L = 1); Sakharov je skalarer Einheit")
T0 = 1e-3                     # unterhalb sind die Windungsterme < e^{-60}: Poisson-Volumenterm exakt


def S_f(t):                   # ein Quarter: Sum_{x in Z+1/4} e^{-t (4 pi x)^2} = Sum_{j>=0} e^{-t pi^2 (2j+1)^2}
    if t < T0:
        return math.sqrt(math.pi/t)/(4*math.pi)
    j = np.arange(int(math.sqrt(900/(t*math.pi**2))/2) + 3)
    return float(np.exp(-t*math.pi**2*(2*j + 1)**2).sum())


def S_b(t):                   # chiraler Stromturm 4 pi m, m >= 1
    if t < T0:
        return math.sqrt(math.pi/t)/(8*math.pi) - 0.5
    m = np.arange(1, int(math.sqrt(900/t)/(4*math.pi)) + 3)
    return float(np.exp(-t*(4*math.pi*m)**2).sum())


def bulk_f(t):
    return math.sqrt(math.pi/t)/(4*math.pi)


def bulk_b(t):
    return math.sqrt(math.pi/t)/(8*math.pi) - 0.5


def uint(f, lo, hi):          # Int dt t^-2 f(t) über u = ln t
    v, _ = quad(lambda u: math.exp(-u)*f(math.exp(u)), math.log(lo), math.log(hi), limit=800, epsabs=0, epsrel=1e-13)
    return v


j0 = np.arange(400)
direct_T0 = float(np.exp(-T0*math.pi**2*(2*j0 + 1)**2).sum())
TH = 60.0                     # oberhalb ist S = 0 bis auf e^{-600}: Rest -Int_TH^inf t^-2 Volumenterm analytisch
fin_f = uint(lambda t: S_f(t) - bulk_f(t), T0, TH) - (math.sqrt(math.pi)/(4*math.pi))*(2/3)*TH**-1.5
fin_b = uint(lambda t: S_b(t) - bulk_b(t), T0, TH) - (math.sqrt(math.pi)/(8*math.pi))*(2/3)*TH**-1.5 + 0.5/TH
exp_f, exp_b = -3*ZETA3/2, 8*ZETA3
say(f"   Endlichgrößenteil (Windungsterme, vor dem Faktor 1/(192 pi^2)): Fermionturm {fin_f:.10f} (erwartet -3 zeta(3)/2 = "
    f"{float(exp_f):.10f}), Stromturm {fin_b:.10f} (erwartet 8 zeta(3) = {float(exp_b):.10f})")
check("Endlichgrößenteil aus den direkten Niveausummen = Poisson-Windungssumme: Fermionturm (alpha = 1/4, 3/4) "
      "16 C(alpha) mit C = Sum cos(2 pi k alpha)/k^3 = -3 zeta(3)/32, chiraler Stromturm 8 zeta(3); bei t = 1e-3 stimmen "
      "direkte Summe und Volumenterm auf 1e-12 überein",
      abs(fin_f/float(exp_f) - 1) < 1e-7 and abs(fin_b/float(exp_b) - 1) < 1e-7
      and abs(direct_T0/bulk_f(T0) - 1) < 1e-12)
check("Vorzeichen: der reglerfreie Teil ist für die originalen Fermionniveaus negativ, für den Stromturm positiv -- die "
      "2D-Bosonisierung überlebt die KK-Lesart nicht; das Vorzeichen hängt an der Lesart, nicht am Original",
      fin_f < 0 < fin_b)
EPS = 1e-4
rho1 = lambda t, e: 1.0 if t >= e else 0.0
rho2 = lambda t, e: -math.expm1(-(t/e)**3)
b1 = uint(lambda t: rho1(t, EPS)*math.sqrt(1/t), EPS, 1e8) + 2/math.sqrt(1e8)
b2 = uint(lambda t: rho2(t, EPS)*math.sqrt(1/t), EPS*1e-12, 1e8) + 2/math.sqrt(1e8)
cl1, cl2 = (2/3)*EPS**-1.5, (2*math.sqrt(math.pi)/3)*EPS**-1.5
check("Volumenterm (Weyl-Gesetz der Niveaus) Int t^-5/2 rho: scharfer Regler (2/3) eps^-3/2, glatter Regler "
      "1 - e^{-(t/eps)^3} (2 sqrt(pi)/3) eps^-3/2 (numerisch auf 1e-8); Verhältnis sqrt(pi): reglerabhängig und nicht "
      "durch eine Superspur getilgt (er zählt dieselben Freiheitsgrade positiv)",
      abs(b1/cl1 - 1) < 1e-8 and abs(b2/cl2 - 1) < 1e-8 and abs(b2/b1 - math.sqrt(math.pi)) < 1e-7)
# zwei Welten: Regler auf der eigenen niedrigsten Stufe m_1 = pi/L, eps = 1/pi^2; gleiche Antworten, verschiedenes G
EPS0 = 1/math.pi**2
sig1 = uint(lambda t: rho1(t, EPS0)*S_f(t), EPS0, 60.0)
sig2 = uint(lambda t: rho2(t, EPS0)*S_f(t), EPS0*1e-6, 60.0)
lp = [math.sqrt(16*math.pi*2*s/(192*math.pi**2)) for s in (sig1, sig2)]
ZAHLEN["kk_lesart"] = {"endlich_fermion_vor_1/(192pi^2)": fin_f, "endlich_strom": fin_b,
                       "L^2/(16 pi G)_endlich_fermion_je_einheit": fin_f/(192*math.pi**2),
                       "L^2/(16 pi G)_endlich_strom_je_einheit": fin_b/(192*math.pi**2),
                       "L/l_P_zeta_strom_2_einheiten": math.sqrt(16*math.pi*2*fin_b/(192*math.pi**2)),
                       "eigene_stufe_eps=1/pi^2_L/l_P_scharf_glatt": lp}
say(f"   Regler auf der eigenen Stufe (eps = 1/pi^2, zwei Einheiten): L/l_P = {lp[0]:.4f} (scharf), {lp[1]:.4f} (glatt)")
check("zwei Welten: gleiche Antworten des Originals (c, h, kappa_n hängen nicht vom Regler ab), Regler auf der eigenen "
      "niedrigsten Stufe pi/L, nur die Form verschieden: L/l_P weicht um mehr als 10 % ab -- das eigene Spektrum legt die "
      "Reglerform nicht fest", abs(lp[1]/lp[0] - 1) > 0.10)
befund("Antwort auf (1): endlich ja; reglerfrei nur der Endlichgrößenteil, und der ist für die Fermionniveaus negativ "
       "(-zeta(3)/(128 pi^2 L^2) je Einheit und Quarter, kein positives G); der Sakharov-Hauptteil ist der Volumenterm des "
       "Weyl-Gesetzes und bleibt reglerabhängig, auch wenn die Reglerskala die eigene Stufe pi/L ist", True)

say("  1d Eigenes Gewicht als Abschneidefunktion (Spektralmoment F_2, vier vorab festgelegte Zuordnungen)")
fermi = lambda v: 1/(1 + mp.exp(v))
F2 = {"Fermi c(x) auf |D|": 2*mp.quad(lambda v: v*fermi(v), [0, mp.inf]),
      "Boltzmann e^{-h} auf |D|": 2*mp.quad(lambda v: v*mp.exp(-v), [0, mp.inf]),
      "Fermi auf D^2": mp.quad(fermi, [0, mp.inf]),
      "Boltzmann auf D^2": mp.quad(lambda v: mp.exp(-v), [0, mp.inf])}
F2_exact = [PI**2/6, mp.mpf(2), mp.log(2), mp.mpf(1)]
spread = max(F2.values())/min(F2.values())
ZAHLEN["F2_momente"] = {k: mp.nstr(v, 12) for k, v in F2.items()}
ZAHLEN["F2_L/l_P_zwei_einheiten"] = {k: mp.nstr(mp.sqrt(v/(3*PI)), 6) for k, v in F2.items()}
check("F_2 = pi^2/6, 2, ln 2, 1 (numerisch auf 1e-30); Spreizung 2/ln 2 = 2,885 > 2: Das Original legt die Gewichtsfunktion "
      "der modularen Energie fest, nicht ob sie auf |D| oder D^2 wirkt (Spektraldimension 1 statt 4, z60)",
      all(abs(v - e) < mp.mpf(10)**-30 for v, e in zip(F2.values(), F2_exact)) and spread > 2)

# ================================================================== Teil 2
say("\nTeil 2: Spur- und Tensorkanal als Koeffizienten der induzierten Wirkung (Normierungen vorab)")
say("  2a Tensorkanal in 2D: Polyakov-Wirkung mit c; der Einstein-Hilbert-Term hat keine Antwort")
NG = 64
kk = np.fft.fftfreq(NG, d=1/NG)
KX, KY = np.meshgrid(kk, kk, indexing="ij")
K2 = KX**2 + KY**2
rng = np.random.default_rng(68)
sig = np.real(np.fft.ifft2((rng.normal(size=(NG, NG)) + 1j*rng.normal(size=(NG, NG)))*np.exp(-K2/18)))
sig *= 0.4/np.abs(sig).max()
FT = np.fft.fft2


def lap(f):
    return np.real(np.fft.ifft2(-K2*FT(f)))


def invlap(f):
    fk, out = FT(f), np.zeros((NG, NG), dtype=complex)
    msk = K2 > 0
    out[msk] = -fk[msk]/K2[msk]
    return np.real(np.fft.ifft2(out))


dA = (2*np.pi/NG)**2
sgR = -2*lap(sig)                                   # sqrt(g) R für g = e^{2 sigma} delta
intR = sgR.sum()*dA
polI = (sgR*invlap(sgR)).sum()*dA                   # Int sqrt(g) R Box_g^-1 R (Box_g f = R <=> Delta f = sqrt(g) R)
gx = np.real(np.fft.ifft2(1j*KX*FT(sig)))
gy = np.real(np.fft.ifft2(1j*KY*FT(sig)))
grad2 = (gx**2 + gy**2).sum()*dA
check("Torus mit zufälligem Konformfaktor: Int sqrt(g) R = 0 (Einstein-Hilbert topologisch, keine Antwort, die eine "
      "Kopplung festlegen könnte) und Int sqrt(g) R Box^-1 R = -4 Int (d sigma)^2 (auf 1e-10): die induzierte Wirkung ist "
      "(c/96 pi) Int R Box^-1 R = -(c/24 pi) Int (d sigma)^2, ihr einziger lokaler Koeffizient ist das dimensionslose c",
      abs(intR) < 1e-10*np.abs(sgR).sum()*dA and abs(polI/(-4*grad2) - 1) < 1e-10 and grad2 > 0)
ZAHLEN["polyakov"] = {"c/(96 pi) Quarterquelle c=2": 2/(96*math.pi), "E8-Naht c=8": 8/(96*math.pi),
                      "Weyl-Anomalie": "-c R/(24 pi)"}

say("  2b Spurkanal: nullte, erste und zweite Antwort am KMS-Punkt beta = 4 pi (tau = i)")
odd = range(1, 401, 2)
prod_odd = mp.fprod(1 + mp.exp(-k*PI) for k in odd)
sum_odd = mp.fsum(k/(mp.exp(k*PI) + 1) for k in odd)
check("Ramanujan bei tau = i: Prod_{k ungerade}(1 + e^{-k pi}) = 2^{1/4} e^{-pi/24} (Klasseninvariante G_1 = 1) und "
      "Sum_{k ungerade} k/(e^{k pi} + 1) = 1/24 (auf 1e-35); die Einteilchenniveaus 4 pi |x| sind genau pi k, k ungerade",
      abs(prod_odd - mp.mpf(2)**(mp.mpf(1)/4)*mp.exp(-PI/24)) < mp.mpf(10)**-35
      and abs(sum_odd - mp.mpf(1)/24) < mp.mpf(10)**-35)
E_tot = [km + 4*PI*ev for km, ev in zip(Kmean, evac)]
S_off = [logZq(b, ALPHAS[0]) - b*mp.diff(lambda bb: logZq(bb, ALPHAS[0]), b) for b in (2*PI, 8*PI)]
check("erste Antwort: thermische Energie <K> = pi/24 je Quarter hebt die Casimirenergie 4 pi (-1/96) exakt auf "
      "(<T> = 0 am KMS-Punkt); nullte Antwort: Entropie S = ln 2/4 je Quarter (auf 1e-35), also ln 2/2 je Naht "
      "(halbes Bit); nicht für beliebiges beta (beta = 2 pi, 8 pi weichen ab)",
      all(abs(e) < mp.mpf(10)**-35 for e in E_tot) and all(abs(s - mp.log(2)/4) < mp.mpf(10)**-35 for s in Sq)
      and all(abs(s - mp.log(2)/4) > 0.1 for s in S_off))


def logZ_E8(beta):
    qq = mp.exp(-beta)
    E4 = 1 + 240*mp.fsum(mp.fsum(dd**3 for dd in range(1, nn + 1) if nn % dd == 0)*qq**nn for nn in range(1, 30))
    return mp.log(E4) + beta/3 - 8*mp.fsum(mp.log(1 - qq**nn) for nn in range(1, 60))


S_E8 = logZ_E8(2*PI) - 2*PI*mp.diff(logZ_E8, 2*PI)
check("E8-Naht bei tau = i: D log Z = 0, also S = log(E4/eta^8)(i) = log j(i)^{1/3} = ln 12 (auf 1e-25)",
      abs(S_E8 - mp.log(12)) < mp.mpf(10)**-25 and abs(mp.diff(logZ_E8, 2*PI)) < mp.mpf(10)**-25)
ZAHLEN["entropie_kms"] = {"je_quarter": "ln 2/4", "quarterquelle_je_naht": "ln 2/2 = " + mp.nstr(S_seam, 15),
                          "E8_naht": "ln 12 = " + mp.nstr(S_E8, 15), "Z_voll_je_quarter": "2^{1/4} = Weber f(i)"}


def kappa2_q(beta, a):
    return beta**2*mp.diff(lambda b: logZq(b, a), beta, 2)


cardy = [kappa2_q(mp.mpf(b), ALPHAS[0])*3*mp.mpf(b)/PI**2 for b in ("0.1", "0.2")]
k2_4pi = kappa2_q(4*PI, ALPHAS[0])
ratio4pi = k2_4pi/(PI/12)
k2_tw = {a_: kappa2_q(4*PI, mp.mpf(a_)) for a_ in ("0", "0.25", "0.5")}
ZAHLEN["kappa2_je_quarter"] = {"beta=4pi": mp.nstr(k2_4pi, 10), "cardy_pi_c/12": mp.nstr(PI/12, 10),
                               "verhaeltnis": mp.nstr(ratio4pi, 6), "cardy_test_beta=0.1,0.2": [mp.nstr(x, 15) for x in cardy],
                               "twist_alpha=0,1/4,1/2": [mp.nstr(v, 6) for v in k2_tw.values()]}
check("zweite Antwort: kappa_2(beta) 3 beta/pi^2 = c = 1 je Quarter im extensiven Grenzfall (beta = 0,1 und 0,2 auf "
      "1e-12, denn log Z = pi^2/(6 beta) + O(beta) bis auf e^{-4 pi^2/beta})", all(abs(x - 1) < 1e-12 for x in cardy))
check("bei beta = 4 pi ist kappa_2/(pi c/12) = " + mp.nstr(ratio4pi, 4) + ", und bei festem c = 1 hängt kappa_2 vom Twist "
      "ab (alpha = 0, 1/4, 1/2: Spreizung > 10): kappa_2 ist am Original eine globale Antwort der Ladungsnullmoden, keine "
      "lokale Kopplung", abs(ratio4pi - 1) > 0.3 and max(k2_tw.values())/min(k2_tw.values()) > 10)
befund("Antwort auf (2): Mit festen Normierungen ist c der einzige lokale Koeffizient der induzierten 2D-Wirkung "
       "(Polyakov, Weyl-Anomalie -c R/(24 pi)); kappa_2 ist die Hesse-Form der Modulfunktion log Z(tau) und fällt nur im "
       "extensiven Grenzfall mit pi^2 c/(3 beta) zusammen. Beide sind dimensionslos; der 2D-Einstein-Hilbert-Term ist "
       "topologisch. In 4D fehlt der Tensorkanal (z59); unter T sind die holomorphen Felder konform, und in jedem konform "
       "invarianten Schema ist der induzierte Einstein-Hilbert-Koeffizient null", True)

# ================================================================== Teil 3
say("\nTeil 3: beta = 4 pi, q = e^{-2 pi} und die Lesart Nahtskala ~ Planck (GRAV.GAPCORR.01)")
import sympy as sp                                   # noqa: E402
levels, ok_int = set(), True
for al in (sp.Rational(1, 4), sp.Rational(3, 4)):
    for Q in range(-6, 7):
        for N in range(0, 12):
            v = 4*(sp.Rational(Q*Q, 2) + (al - sp.Rational(1, 2))*Q) + 4*N
            ok_int &= bool(v.is_integer)
            levels.add(int(v))
nz = sorted(x for x in levels if x > 0)
window = [x for x in nz if math.pi*x <= 100]
check("alle Niveaus des Modulgenerators sind pi n mit n ganz; kleinste Stufe pi (Twist), Stromquant 4 pi, Paar 8 pi; mit "
      "beta_phys = L heißt das E = pi n/L, Temperatur 1/L",
      ok_int and nz[0] == 1 and 4 in nz and 8 in nz)
check("keine Hierarchie aus den Niveaus: unter allen Zuständen mit Gewicht e^{-K} >= e^{-100} ist das größte "
      "Energieverhältnis " + str(window[-1]) + " (< 100); M_Pl/m_e ~ 10^22 als Niveauverhältnis hätte das Gewicht "
      "e^{-pi 10^22}", window[-1]/nz[0] < 100 and window[-1] == int(100/math.pi))
rng2 = np.random.default_rng(6800)
lt = np.log(10**rng2.uniform(-30, -3, 4000))
dist = np.abs(lt - np.round(lt/(-2*np.pi))*(-2*np.pi))
ZAHLEN["q_turm_abdeckung"] = {"max_abstand_ln": float(dist.max()), "halbe_stufe_pi": math.pi}
befund("q = e^{-2 pi} tritt nur als Gewicht auf (Boltzmannfaktor je Paar, Toeplitz-Restmischung q^Tiefe), nicht als "
       "Niveauverhältnis; ein Turm q^n mit freier Tiefe n trifft jedes Verhältnis bis auf den Faktor e^pi = 23 (größter "
       f"Abstand im Logarithmus {dist.max():.3f} <= pi): ein solcher Treffer wäre kein Beleg", True)
befund("Antwort auf (3): beta = 4 pi und q legen nur dimensionslose Verhältnisse innerhalb der Nahtskala fest (Form "
       "tau = i, Stufen pi : 4 pi : 8 pi, T L = 1, Gewichte e^{-pi}, q, q^2). Die Lesart Nahtskala ~ Planck setzt "
       "L = lambda l_P mit offenem lambda; E_Twist/M_Pl = pi/lambda hängt genau an der offenen Zahl a/l_P", True)

# ================================================================== Teil 4
say("\nTeil 4: Welche dimensionslose Zahl wird vorhergesagt?")
cand = {"Spektralmoment F_2 (vier Zuordnungen)": [float(mp.sqrt(v/(3*PI))) for v in F2.values()],
        "Regler auf eigener Stufe pi/L (scharf, glatt)": lp,
        "Zeta-KK Stromturm (2 Einheiten)": [ZAHLEN["kk_lesart"]["L/l_P_zeta_strom_2_einheiten"]],
        "Jacobson-Zelle a^2 = 4 S l_P^2, S = ln 2/2": [math.sqrt(4*float(S_seam))]}
allv = [x for v in cand.values() for x in v]
ZAHLEN["kandidaten_L/l_P"] = cand
ZAHLEN["zeta_KK_fermion"] = "L^2/(16 pi G) = -zeta(3)/(128 pi^2) je Einheit und Quarter < 0: ausgeschlossen"
for k_, v_ in cand.items():
    say(f"   {k_}: L/l_P = {', '.join(f'{x:.3f}' for x in v_)}")
check("die vorab festgelegten Lesarten geben für L/l_P Werte über mehr als den Faktor 2 (" +
      f"{min(allv):.3f} bis {max(allv):.3f}" + "), die reglerfreie Fermionlesart sogar G < 0: keine eindeutige Zahl, "
      "also kein Datenvergleich", max(allv)/min(allv) > 2 and fin_f < 0)
befund("Skalensymmetrie: Kovarianz und alle Antworten hängen nur von modularen Größen ab (x, modulare Zeit, tau); die "
       "Dilatation L -> lambda L bildet eine Realisierung auf eine mit identischen Antworten ab, G/L^2 bleibt dabei eine "
       "freie Zahl. Sie wäre der Koeffizient einer Antwort der Dimension Masse^2; in 2D ist der einzige Kandidat "
       "topologisch (Teil 2a), in 4D fehlt der Tensorkanal (z59, z60)", True)
befund("Unter Prämisse D (Jacobson) ist der dimensionslose Inhalt von G jetzt exakt eine Antwort: S = ln 2/2 je Naht "
       "(E8-Naht ln 12); offen bleibt allein die Fläche je Naht, A = 4 S l_P^2 -- eine geometrische, keine spektrale Zahl. "
       "Unter Sakharov bleibt das Reglermoment offen", True)
befund("Antwort auf (4): vorhergesagt wird keine dimensionslose Zahl mit Messwert; M_Pl/m_e und a/l_P bleiben offen", True)

STATUS.update({
    "E": ["Zeta der Einteilchenniveaus: Pol nur bei s = 1, Z(0) = 0, Casimir -1/96 = h - c/24",
          "Cardy-Wachstum der Entartungen: keine Zetafunktion des vollen Spektrums",
          "KMS-Punkt: <K> + Casimir = 0, S = ln 2/4 je Quarter (Ramanujan G_1 = 1), E8-Naht S = ln 12",
          "Polyakov-Normierung, 2D-Einstein-Hilbert ohne Antwort", "kappa_2 -> pi^2 c/(3 beta) nur extensiv",
          "Niveaus pi n/L, keine Hierarchie"],
    "C": ["Endlichgrößenteil reglerfrei, Volumenterm reglerabhängig (Wärmekernrahmen)"],
    "L": ["KK-Lesart (Niveaus als 4D-Massen)", "Spektralmoment-Zuordnungen", "Jacobson-Zelle A = 4 S l_P^2",
          "Nahtskala ~ Planck"],
    "O": ["a/l_P", "M_Pl/m_e", "Fläche je Naht", "4D-Tensorkanal und Einstein-Hilbert-Koeffizient"]})
urteil = ("Das eigene Spektrum des Originals legt keine dimensionslose Skalenzahl fest. Sakharov-Summen über die "
          "rekonstruierten Niveaus sind endlich, aber nur am KMS-Punkt (volles Spektrum, Cardy-Wachstum) oder zeta-"
          "regularisiert (Einteilchen: Z(0) = 0, Casimir -1/96 = h - c/24); reglerfrei ist dort nur der Endlichgrößenteil, "
          "und der ist für die originalen Fermionniveaus negativ (-zeta(3)/(128 pi^2 L^2)), für den Stromturm positiv. Der "
          "Sakharov-Hauptteil ist der Volumenterm des Weyl-Gesetzes und bleibt reglerabhängig, auch mit der eigenen Stufe "
          "pi/L als Skala. c und kappa_2 sind dimensionslos: c ist der einzige lokale Koeffizient der 2D-Wirkung, kappa_2 die "
          "Hesse-Form der Modulfunktion; der 2D-Einstein-Hilbert-Term ist topologisch. Neu und exakt: Am KMS-Punkt hebt die "
          "thermische Energie die Casimirenergie auf, und die Entropie ist ln 2/4 je Quarter (ln 2/2 je Naht, E8-Naht ln 12). "
          "Unter Prämisse D ist damit der dimensionslose Inhalt von G eine Antwort; offen bleibt die Fläche je Naht. "
          "Vorhergesagt wird keine Zahl mit Messwert; a/l_P und M_Pl/m_e bleiben offen.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z68_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "status": STATUS, "urteil": urteil}, indent=1,
    ensure_ascii=False, default=str) + "\n")
