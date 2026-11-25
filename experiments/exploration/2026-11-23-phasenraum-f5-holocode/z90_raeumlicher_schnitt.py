"""Z90: Gravitation über einen räumlichen Schnitt -- Flächengesetz und nativer Regulator (Folgeauftrag F2, Stand H6h44).

Leitlinie 'Rekonstruktion statt Auswahl'. Firewalls: GRAV.NONCIRCULAR.01 (kein Schluss 'modular => Bisognano-Wichmann'; die
Kappenentropie wird direkt aus dem reduzierten Gaußzustand berechnet; Hislop-Longo wird nur in der erlaubten Reihenfolge benutzt:
zuerst das freie konforme Netz des Leiterfelds aus z83, dann der Satz), DIMENSION.UPLIFT.FIREWALL.01 (alles ist Kinematik eines
freien 3+1-Felds; keine wechselwirkende Schließung, keine Aussage über die Naht; die Fockquantisierung des Leiterfelds ist die
Standardannahme, ihre Herkunft aus vielen Nähten ist im Repo offen). Bestand: z83 (S2: nu = 1-Leiter = masseloser konformer
Skalar auf R x S^3, E = m, Vielfachheit m^2, Original x Partner = CFT bei beta' = 8 pi), z86 (S5: ln 2/2 je Naht, Typ I, keine
Flächennormierung; Empfehlung: Kappe von S^3 mit q^m als Regulator), z68/z78 (q = e^{-2 pi}, tau = i, KMS bei beta = 4 pi).
Ledger GRAV.NONLINEAR.01, GRAV.ENTROPY.EQUILIBRIUM.01, SEAM.AREACOEFF.01/02 (nur gelesen). Zehnte externe Analyse:
S = A/(48 pi eps^2) + ..., der Koeffizient hängt am Regulator, zwei Profile gleicher Nennlänge unterscheiden sich um 2.

Teil 1  Regulatorabhängigkeit und universelle Teile (konformer Skalar, 4D): Kegelwärmeleitungskern (1/12)(1/n - n) aus der
        Bildsumme, S = A I[f]/(48 pi) mit dem Profilfunktional I[f] = int f(s) s^-2 ds (symbolisch für sieben Profile, Faktor 2
        bei gleicher Nennlänge); universeller Logarithmus -1/90 = -4a dreifach (zeta(0) des konformen Laplace auf S^4 über
        Hurwitz und über a_2, CHM-Entropiedichte auf H^3 mal Logterm des Volumens); CHM-Flächenkoeffizient 1/(360 pi).
Teil 2  Kappe von S^3 im Leiterfeld mit nativen Regulatoren: Kanäle l mit Pöschl-Teller-Spektrum m >= l + 1, xi = 1/6 aus
        E = m; native Gewichte als Länge (G_eps = 1/(8 pi^2 (cosh eps - cos gamma))); scharfe native Abschneidung m <= M mit
        zwei Lokalisierungen (Gauß-Gegenbauer-Knoten = exakte Isometrie, Slepian-Teilung), Flächengesetz in M und Koeffizient;
        KMS-Gewichte (beta = 2 pi = tau = i, 4 pi, 8 pi); Schwellenabschneidung M(delta); Modulschnitt (Hislop-Longo nach dem
        Netz): kontinuierliches Modulspektrum (D_1), UV-endliche Dichte 1/(180 pi), Divergenz im H^3-Volumen, S/A_mod = a/pi.
Teil 3  Vergleich mit 1/(4G) in TFPT-Normierung (8 pi G = 1/c3 = 8 pi): nötige Längen je Regulator, induzierte Kopplung
        (Sakharov) gegen Verschränkung bei xi = 0 und beim nativen xi = 1/6, pi-Kommensurabilität 1/(48 pi) = c3/6,
        Dimensionstypisierung des Ledger-Satzes 'Flächenkoeffizient = Zentralladung'; was für G fehlt.
Teil 4  Urteil für F2 (E/C/L/O).
"""
import json
import time
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sp
from numpy.polynomial.legendre import leggauss
from scipy.linalg import eigh_tridiagonal
from scipy.special import gammaln, roots_gegenbauer

HERE = Path(__file__).parent
mp.mp.dps = 50
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
T0 = time.time()
PI = mp.pi
TOL = mp.mpf(10)**-40


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val):
    BEFUNDE[name] = bool(val)
    say(f"   [Befund: {'ja' if val else 'nein'}] {name}")


def ns(x, n=12):
    return mp.nstr(x, n)


# =================================================================== Teil 1
say("Teil 1: Regulatorabhängigkeit des Flächenkoeffizienten und universelle Teile (freier konformer Skalar, 4D)")
n_, N_, A_, I_, u_, R_, e_ = sp.symbols("n N A I u R epsilon", positive=True)
xi = sp.symbols("xi")

# --- 1.1 Kegel: Bildsumme der Orbifolds C/Z_N und Replikaentropie
bild = {N: mp.fsum(1/(4*mp.sin(PI*k/N)**2) for k in range(1, N)) for N in range(2, 41)}
ok_bild = all(abs(bild[N] - mp.mpf(N*N - 1)/12) < TOL for N in bild)
ok_bild &= sp.simplify((N_**2 - 1)/(12*N_) - (1/n_ - n_).subs(n_, 1/N_)/12) == 0
dTrK = lambda xi_: (1/n_ - n_)/12*A_/(4*sp.pi) - xi_*4*sp.pi*(1 - n_)*A_/(16*sp.pi**2)   # mal 1/s; int_Kegel R = 4 pi (1-n) A
dlnZ = lambda xi_: sp.Rational(1, 2)*I_*dTrK(xi_)                                      # (1/2) int ds/s f(s) (...)/s
S_kegel = lambda xi_: sp.simplify(sp.limit(dlnZ(xi_)/(1 - n_), n_, 1))
check("Kegel: Tr K(C_n) - n Tr K(R^2) = (1/12)(1/n - n) aus der Bildsumme der Orbifolds C/Z_N (Summe_k 1/(4 sin^2(pi k/N)) = "
      "(N^2 - 1)/12, N = 2..40 auf 1e-40, n = 1/N); mit dem Querfaktor A/(4 pi s) und ln Z_n = (1/2) int ds/s f(s) Tr K_n folgt "
      "S = A I[f]/(48 pi), I[f] = int_0^inf f(s) s^-2 ds (minimal gekoppelt)",
      ok_bild and sp.simplify(S_kegel(0) - A_*I_/(48*sp.pi)) == 0)

# --- 1.2 Profile: I[f] in Einheiten 1/eps^2 (s = eps^2 u)
PROFILE = {
    "scharf s >= eps^2": (None, 1),
    "scharf, Brownsche Konvention s >= eps^2/2 (sigma^2 = 2s)": (None, sp.Rational(1, 2)),
    "Q1 = e^{-eps^2/s}": (sp.exp(-1/u_), None),
    "Q2 = e^{-eps^2/s}(1 + eps^2/s)": (sp.exp(-1/u_)*(1 + 1/u_), None),
    "P2 = 1 - (1 + s/eps^2) e^{-s/eps^2}": (1 - (1 + u_)*sp.exp(-u_), None),
    "P3 = 1 - (1 + s/eps^2 + s^2/(2 eps^4)) e^{-s/eps^2}": (1 - (1 + u_ + u_**2/2)*sp.exp(-u_), None),
    "Pauli-Villars (1 - e^{-s/eps^2})^2 (Massen 1/eps, sqrt2/eps)": ((1 - sp.exp(-u_))**2, None),
}
I_wert, zul = {}, True
for k_, (f, unten) in PROFILE.items():
    if f is None:
        I_wert[k_] = sp.integrate(u_**-2, (u_, unten, sp.oo))
    else:
        I_wert[k_] = sp.simplify(sp.integrate(f/u_**2, (u_, 0, sp.oo)))
        zul &= sp.limit(f, u_, sp.oo) == 1 and sp.limit(f/u_, u_, 0) == 0                # Profil -> 1 im IR, f/s -> 0 im UV
soll = [1, 2, 1, 2, 1, sp.Rational(1, 2), 2*sp.log(2)]
ok_prof = zul and all(sp.simplify(sp.expand_log(v - w, force=True)) == 0 for v, w in zip(I_wert.values(), soll))
ZAHLEN["profile_I_mal_eps2"] = {k: str(v) for k, v in I_wert.items()}
ZAHLEN["profile_S_mal_eps2_durch_A"] = {k: ns(mp.mpf(str(sp.N(v, 60)))/(48*PI), 8) for k, v in I_wert.items()}
say("   I[f] eps^2: " + "; ".join(f"{k}: {v}" for k, v in I_wert.items()))
check("Profile (symbolisch): I[f] eps^2 = 1 (scharf, Q1, P2), 2 (Q2, Brownsche Konvention), 1/2 (P3), 2 ln 2 (Pauli-Villars); "
      "alle zulässig (f -> 1 für s -> inf, f = o(s) für s -> 0). Die externe Aussage stimmt: Q1/Q2 und P3/P2 haben dieselbe "
      "Nennlänge eps und unterscheiden sich um den Faktor 2, ebenso dasselbe scharfe Profil unter den Konventionen s = eps^2 "
      "und sigma^2 = 2s = eps^2. Der Flächenkoeffizient ist ein lineares Funktional des Profils, kein Datum des Felds",
      ok_prof and I_wert["Q2 = e^{-eps^2/s}(1 + eps^2/s)"]/I_wert["Q1 = e^{-eps^2/s}"] == 2
      and I_wert["P2 = 1 - (1 + s/eps^2) e^{-s/eps^2}"]/I_wert["P3 = 1 - (1 + s/eps^2 + s^2/(2 eps^4)) e^{-s/eps^2}"] == 2)
befund("zeta- und dimensionale Regularisierung setzen skalenfreie Potenzdivergenzen null: dort ist der Flächenkoeffizient 0. "
       "Zwischen 0 und jedem positiven Wert ist alles durch Wahl des Profils erreichbar", True)

# --- 1.3 universeller Logarithmus: dreifach
B = sp.bernoulli
zH = lambda m, a: -B(m + 1, a)/(m + 1)                           # Hurwitz zeta(-m, a)
a32 = sp.Rational(3, 2)
zeta0_hurwitz = sp.Rational(1, 3)*zH(3, a32) - sp.Rational(1, 12)*zH(1, a32) - sp.Rational(1, 96) + sp.Rational(1, 192)
# Eigenwerte (l+1)(l+2) = nu^2 - 1/4, Entartung (nu^3 - nu/4)/3, nu = l + 3/2; Polterme aus (1 - 1/(4 nu^2))^{-z}
ok_dl = all(sp.Rational((l + 1)*(l + 2)*(2*l + 3), 6) == ((l + a32)**3 - (l + a32)/4)/3 for l in range(20))
ok_ew = all((l*(l + 3) + 2) == (l + 1)*(l + 2) for l in range(20))           # -Laplace + R/6 = -Laplace + 2 auf S^4
zeta0_a2 = (sp.Rational(1, 180)*(24 - 36) + sp.Rational(1, 2)*(sp.Rational(1, 6) - sp.Rational(1, 6))**2*144)*(
    sp.Rational(8, 3)*sp.pi**2)/(16*sp.pi**2)                                  # Riem^2 = 24, Ric^2 = 36, R = 12, Vol = 8 pi^2/3
a_skalar = sp.Rational(1, 360)
x_ = sp.symbols("x", positive=True)
bose = sp.gamma(4)*sp.zeta(4)                                                  # int_0^inf x^3/(e^x - 1) dx
ok_bose = abs(mp.quad(lambda t: t**3/(mp.e**t - 1), [0, mp.inf]) - mp.mpf(str(sp.N(bose, 60)))) < mp.mpf(10)**-30
T_ = sp.symbols("T", positive=True)
f_frei = -T_**4/(2*sp.pi**2)*bose/3                                            # T int k^2/(2 pi^2) ln(1 - e^{-k/T}) dk
s_frei = -sp.diff(f_frei, T_)
s_mod = sp.simplify(s_frei.subs(T_, 1/(2*sp.pi)))
t_, r_ = sp.symbols("t rho", positive=True)
K_H3 = (4*sp.pi*t_)**sp.Rational(-3, 2)*(r_/sp.sinh(r_))*sp.exp(-t_ - r_**2/(4*t_))
ok_H3 = sp.simplify(sp.diff(K_H3, t_) - (sp.diff(K_H3, r_, 2) + 2*sp.coth(r_)*sp.diff(K_H3, r_))) == 0
ok_H3 &= sp.simplify(sp.limit(K_H3*sp.exp(t_), r_, 0) - (4*sp.pi*t_)**sp.Rational(-3, 2)) == 0
V_H3 = sp.simplify(4*sp.pi*sp.integrate(sp.sinh(u_)**2, (u_, 0, u_)))
S_CHM = sp.expand(sp.expand_log(sp.simplify((s_mod*V_H3).subs(u_, sp.log(2*R_/e_)).rewrite(sp.exp)), force=True))
flaeche_CHM = sp.simplify(S_CHM.coeff(e_, -2)/(4*sp.pi*R_**2))
log_CHM = sp.simplify(S_CHM.coeff(sp.log(R_)))                                     # Koeffizient von ln R (= von ln(R/eps))
ZAHLEN["universell"] = {"zeta(0) S^4 (Hurwitz)": str(zeta0_hurwitz), "zeta(0) S^4 (a_2)": str(sp.simplify(zeta0_a2)),
                        "-4a (a = 1/360)": str(-4*a_skalar), "s(T=1/2pi) auf H^3": str(s_mod),
                        "S_CHM": str(S_CHM), "CHM-Flächenkoeffizient je eps^-2": str(flaeche_CHM),
                        "CHM/Wärmeleitungskern (scharf)": str(sp.simplify(flaeche_CHM*48*sp.pi))}
say(f"   zeta(0)_S^4 = {zeta0_hurwitz} (Hurwitz) = {sp.simplify(zeta0_a2)} (a_2); s(1/2 pi) = {s_mod}; S_CHM = {S_CHM}")
check("universeller Logarithmus, dreifach und regulatorfrei: zeta(0) des konformen Laplace auf S^4 ist -1/90 über die "
      "Hurwitz-Fortsetzung (Eigenwerte (l+1)(l+2), Entartung (l+1)(l+2)(2l+3)/6) und über a_2 = (1/180)(Riem^2 - Ric^2) "
      "(xi = 1/6 tilgt R^2); CHM: der Wärmeleitungskern von -Laplace - 1 auf H^3 ist am Punkt der flache (4 pi t)^{-3/2} "
      "(Lösung der Wärmeleitungsgleichung symbolisch), also s(T = 1/2 pi) = (2 pi^2/45) T^3 = 1/(180 pi), und das Volumen "
      "pi (sinh 2u - 2u) mit u = ln(2R/eps) gibt den Logkoeffizienten -1/90 = -4a, a = 1/360",
      zeta0_hurwitz == sp.Rational(-1, 90) and sp.simplify(zeta0_a2 + sp.Rational(1, 90)) == 0 and ok_dl and ok_ew
      and ok_bose and ok_H3 and sp.simplify(s_mod - 1/(180*sp.pi)) == 0 and log_CHM == sp.Rational(-1, 90)
      and -4*a_skalar == sp.Rational(-1, 90))
check("CHM-Schnitt als weiterer Regulator: S = R^2/(90 eps^2) - (1/90) ln(2R/eps) - eps^2/(1440 R^2) exakt, also Flächen"
      "koeffizient 1/(360 pi) je eps^-2 gegen 1/(48 pi) beim scharfen Eigenzeitprofil (Verhältnis 2/15); unter eps -> lambda eps "
      "skaliert der Flächenterm mit lambda^-2, der Logterm verschiebt nur die Konstante",
      sp.simplify(flaeche_CHM - 1/(360*sp.pi)) == 0 and sp.simplify(flaeche_CHM*48*sp.pi - sp.Rational(2, 15)) == 0
      and sp.simplify(sp.expand_log(S_CHM - (R_**2/(90*e_**2) - sp.log(2*R_/e_)/90 - e_**2/(1440*R_**2)), force=True)) == 0)
befund("Teil 1: Die externe Analyse ist bestätigt. Fest ist nur der Logarithmus -1/90 = -4a; der Flächenkoeffizient ist "
       "I[f]/(48 pi), ein Funktional des Regulatorprofils (0 bis beliebig; Faktor 2 bei gleicher Nennlänge; CHM 2/15 des "
       "scharfen Werts). Ein Flächengesetz mit fester Zahl braucht also einen ausgezeichneten Regulator", True)
say(f"   (Laufzeit bisher {time.time() - T0:.1f} s)")

# =================================================================== Teil 2
say("\nTeil 2: Kappe von S^3 im Leiterfeld (z83) mit nativen Regulatoren")


# --- 2.1 das Leiterfeld je Kanal
def gegen_on(K, lam, x):
    """Orthonormale Gegenbauer-Polynome (Gewicht (1 - x^2)^(lam - 1/2)) vom Grad 0..K-1 an den Stellen x."""
    def logh(k):
        return (np.log(np.pi) + (1 - 2*lam)*np.log(2) + gammaln(k + 2*lam) - gammaln(k + 1) - np.log(k + lam)
                - 2*gammaln(lam))
    P = np.zeros((K, len(x)))
    P[0] = np.exp(-0.5*logh(0))
    pm1, p0 = np.zeros_like(x), P[0]
    for k in range(K - 1):
        a = 2*(k + lam)*np.exp(0.5*(logh(k) - logh(k + 1)))/(k + 1)
        b = (k + 2*lam - 1)*np.exp(0.5*(logh(k - 1) - logh(k + 1)))/(k + 1) if k > 0 else 0.0
        p1 = a*x*p0 - b*pm1
        P[k + 1] = p1
        pm1, p0 = p0, p1
    return P


def s_nu(nu):
    nu = nu[nu > 0.5 + 1e-14]
    return float(np.sum((nu + 0.5)*np.log(nu + 0.5) - (nu - 0.5)*np.log(nu - 0.5)))


def nu_sympl(X, P):
    L = np.linalg.cholesky(X)
    return np.sqrt(np.clip(np.linalg.eigvalsh(L.T @ P @ L), 0.25, None))


def faktoren(om, beta):
    if beta is None:
        return 1/(2*om), om/2
    c = 1/np.tanh(beta*om/2)
    return c/(2*om), om*c/2


def kanal_gauss(l, M):
    x, w = roots_gegenbauer(M - l, l + 1)
    return np.arccos(x), (gegen_on(M - l, l + 1, x)*np.sqrt(w)).T                 # O[i, k] orthogonal


def S_gauss(thetas, M, beta=None, komplement=False):
    S = np.zeros(len(thetas))
    for l in range(M):
        chi, O = kanal_gauss(l, M)
        fx, fp = faktoren(np.arange(M - l) + l + 1.0, beta)
        for j, th in enumerate(thetas):
            A = np.where((chi >= th) if komplement else (chi < th))[0]
            if len(A):
                S[j] += (2*l + 1)*s_nu(nu_sympl((O[A]*fx) @ O[A].T, (O[A]*fp) @ O[A].T))
    return S


def gram_kappe(l, M, th):
    t, wt = leggauss(2*(M - l) + l + 60)
    chi, wchi = 0.5*th*(t + 1), 0.5*th*wt
    U = gegen_on(M - l, l + 1, np.cos(chi))*np.sin(chi)**(l + 1)                 # u_k(chi), orthonormal auf (0, pi)
    return (U*wchi) @ U.T


def S_slepian(thetas, M, beta=None, schwelle=0.5):
    S = np.zeros(len(thetas))
    for l in range(M):
        fx, fp = faktoren(np.arange(M - l) + l + 1.0, beta)
        for j, th in enumerate(thetas):
            lam, V = np.linalg.eigh(gram_kappe(l, M, th))
            W = V[:, lam > schwelle]
            if W.shape[1]:
                S[j] += (2*l + 1)*s_nu(nu_sympl((W.T*fx) @ W, (W.T*fp) @ W))
    return S


nn_ = sp.symbols("n", integer=True, nonnegative=True)
xi_nativ = sp.solve(sp.Poly(sp.expand(nn_*(nn_ + 2) + 6*xi - (nn_ + 1)**2), nn_).coeffs(), xi)
iso = max(float(np.abs(kanal_gauss(l, 128)[1].T @ kanal_gauss(l, 128)[1] - np.eye(128 - l)).max()) for l in range(0, 128, 7))
voll = max(float(np.abs(gram_kappe(l, 40, np.pi) - np.eye(40 - l)).max()) for l in range(0, 40, 3))
m_, g_, eps_ = sp.symbols("m gamma eps", positive=True)
kern_ok = sp.limit(m_*sp.sin(m_*g_)/(2*sp.pi**2*sp.sin(g_)), g_, 0) == m_**2/(2*sp.pi**2)
G_eps = lambda e, g: 1/(8*PI**2*(mp.cosh(e) - mp.cos(g)))
kern_ok &= all(abs(mp.nsum(lambda m: mp.exp(-e*m)*mp.sin(m*g), [1, mp.inf])/(4*PI**2*mp.sin(g)) - G_eps(e, g)) < TOL
               for e, g in ((mp.mpf("0.3"), 1), (2*PI, mp.mpf("0.5")), (4*PI, 2)))
kern_ok &= sp.simplify(1/(8*sp.pi**2*(1 - sp.cos(g_))) - 1/(16*sp.pi**2*sp.sin(g_/2)**2)) == 0
check("[E] Leiterfeld je S^2-Kanal l: -u'' + l(l+1) u/sin^2 chi mit Frequenzen m = l+1, l+2, ... (Summe_l (2l+1) = m^2); "
      f"die ganzzahlige Leiter E = m erzwingt xi = {xi_nativ[xi]} (n(n+2) + 6 xi = (n+1)^2 für alle n). Native Abschneidung "
      f"m <= M exakt als Gauß-Gegenbauer-Knotengitter (Isometrie, Fehler {iso:.1e} bei M = 128) und als Kappen-Gram (volle "
      f"Sphäre = Einheit, {voll:.1e}); Projektorkern m sin(m gamma)/(2 pi^2 sin gamma), Vakuumkorrelator mit Dämpfung e^{{-eps m}} "
      "= 1/(8 pi^2 (cosh eps - cos gamma)) (1e-40), für eps -> 0 der flache 1/(4 pi^2 (2 sin(gamma/2))^2)",
      xi_nativ == {xi: sp.Rational(1, 6)} and iso < 1e-10 and voll < 1e-10 and kern_ok
      and all(sum(2*l + 1 for l in range(m)) == m*m for m in range(1, 30)))

# --- 2.2 native Gewichte als Länge
flach = {k: 2/(mp.cosh(v) - 1) for k, v in (("2 pi (q^m, tau = i)", 2*PI), ("4 pi (KMS Original)", 4*PI),
                                            ("8 pi (beta' aus z83)", 8*PI))}
ZAHLEN["native_daempfung_relative_variation_G_ueber_S3"] = {k: ns(v, 4) for k, v in flach.items()}
check("[E] native Gewichte als Regulatorlänge: e^{-eps m} ist die euklidische Zeitspaltung eps; nativ eps = 2 pi, 4 pi, 8 pi > pi "
      "(größter geodätischer Abstand auf S^3). Der so regulierte Korrelator ist über die ganze S^3 flach: relative Variation "
      f"2/(cosh eps - 1) = {', '.join(ns(v, 3) for v in flach.values())}. Ein Flächengesetz braucht eps << min(theta, pi - "
      "theta) <= pi/2; die nativen Gewichte liegen auf der falschen Seite (Krümmungsskala, nicht UV)",
      all(v < 0.01 for v in flach.values()) and 2*PI > PI)

# --- 2.3 scharfe native Abschneidung m <= M: Flächengesetz und Koeffizient je Lokalisierung
THETAS = [np.pi/6, np.pi/4, np.pi/3, np.pi/2]
MS = np.arange(32, 129, 12)
s2 = np.sin(np.array(THETAS))**2
SG = np.array([S_gauss(THETAS, M) for M in MS])
SS = np.array([S_slepian(THETAS, M) for M in MS])
say(f"   (Laufzeit bisher {time.time() - T0:.1f} s)")
fit = {}
for name, S in (("Gauss", SG), ("Slepian", SS)):
    for j, th in enumerate(THETAS):
        A3 = np.vstack([s2[j]*MS**2, MS, np.ones(len(MS))]).T
        c3_, *_ = np.linalg.lstsq(A3, S[:, j], rcond=None)
        A4 = np.vstack([s2[j]*MS**2, MS, np.log(MS), np.ones(len(MS))]).T
        c4_, *_ = np.linalg.lstsq(A4, S[:, j], rcond=None)
        fit[(name, j)] = (c3_, float(np.abs(A3 @ c3_ - S[:, j]).max()), c4_)
c2 = {name: np.array([fit[(name, j)][0][0] for j in range(4)]) for name in ("Gauss", "Slepian")}
streu = {name: float(np.ptp(v)/np.mean(v)) for name, v in c2.items()}
c2m = {name: float(np.mean(v)) for name, v in c2.items()}
delta_eff = {name: float(fit[(name, 3)][0][1]/(2*fit[(name, 3)][0][0])) for name in c2}
clog = {name: [float(fit[(name, j)][2][2]) for j in range(4)] for name in c2}
S_komp = S_gauss([np.pi/3], 64, komplement=True)[0]
S_kap = S_gauss([np.pi/3], 64)[0]
S_kompS = S_gauss([np.pi - np.pi/3], 64)[0]
ZAHLEN["kappe_scharfe_abschneidung"] = {
    "M": [int(m) for m in MS],
    "S_Gauss(theta=pi/6,pi/4,pi/3,pi/2)": {int(m): [round(float(v), 4) for v in SG[i]] for i, m in enumerate(MS)},
    "S_Slepian(theta=pi/6,pi/4,pi/3,pi/2)": {int(m): [round(float(v), 4) for v in SS[i]] for i, m in enumerate(MS)},
    "c2 je theta (S = c2 sin^2 theta M^2 + c1 M + c0)": {k: [round(float(x), 6) for x in v] for k, v in c2.items()},
    "c2 Mittel": {k: round(v, 6) for k, v in c2m.items()}, "theta-Streuung c2": {k: round(v, 4) for k, v in streu.items()},
    "S/A je eps^-2 mit eps = 1/M (c2/4 pi)": {k: round(v/(4*np.pi), 6) for k, v in c2m.items()},
    "Verhältnis zu 1/(48 pi) (12 c2)": {k: round(12*v, 4) for k, v in c2m.items()},
    "c1 (theta = pi/2)": {k: round(float(fit[(k, 3)][0][1]), 4) for k in c2},
    "effektive Verschiebung der Abschneidung delta = c1/(2 c2)": {k: round(v, 3) for k, v in delta_eff.items()},
    "Logkoeffizient aus 4-Parameter-Fit (theta = pi/6..pi/2)": {k: [round(x, 3) for x in v] for k, v in clog.items()},
    "max Fitrest (3 Parameter)": {k: round(max(fit[(k, j)][1] for j in range(4)), 3) for k in c2}}
say("   S/(sin^2 theta M^2) bei M = 128: Gauß " + ", ".join(f"{v:.5f}" for v in SG[-1]/(s2*MS[-1]**2)) + "; Slepian "
    + ", ".join(f"{v:.5f}" for v in SS[-1]/(s2*MS[-1]**2)))
say(f"   c2 (Fit über M = 32..128): Gauß {c2m['Gauss']:.5f} (Streuung {streu['Gauss']:.2%}), Slepian {c2m['Slepian']:.5f} "
    f"(Streuung {streu['Slepian']:.2%}); Logkoeffizienten frei gefittet: Gauß {[round(x, 2) for x in clog['Gauss']]}, "
    f"Slepian {[round(x, 2) for x in clog['Slepian']]}")
check("[E] Reinheit: Kappe und Komplement (Gauß-Knoten, theta = pi/3, M = 64) haben dieselbe Entropie, und die Kappe theta "
      "gleicht der Gegenkappe pi - theta (Knotensymmetrie x -> -x)",
      abs(S_komp - S_kap) < 1e-7*S_kap and abs(S_kompS - S_kap) < 1e-7*S_kap)
check("[E/C] native Abschneidung m <= M gibt ein Flächengesetz: S = c2 sin^2 theta M^2 + c1 M + c0 mit theta-unabhängigem c2 "
      f"(Streuung über theta = pi/6..pi/2: Gauß {streu['Gauss']:.1%}, Slepian {streu['Slepian']:.1%}), A = 4 pi sin^2 theta; "
      "bedingt nur durch die S2-Identifikation von z83 und die Fockquantisierung",
      streu["Gauss"] < 0.03 and streu["Slepian"] < 0.03 and all(fit[(k, j)][1] < 1.0 for k in c2 for j in range(4)))
check(f"[E] der Koeffizient hängt an der Lokalisierung derselben nativen Abschneidung: c2 = {c2m['Gauss']:.4f} (Gauß-Knoten) "
      f"gegen {c2m['Slepian']:.4f} (Slepian-Teilung), Verhältnis {c2m['Slepian']/c2m['Gauss']:.3f}; in Einheiten des "
      f"Wärmeleitungskerns bei eps = 1/M: {12*c2m['Gauss']:.3f} bzw. {12*c2m['Slepian']:.3f} x 1/(48 pi). Auch die Nennlänge "
      f"selbst ist unscharf: der lineare Term entspricht M -> M + delta mit delta = {delta_eff['Gauss']:+.2f} (Gauß) und "
      f"{delta_eff['Slepian']:+.2f} (Slepian)",
      c2m["Slepian"]/c2m["Gauss"] - 1 > 0.1 and delta_eff["Gauss"] > 0.5 and delta_eff["Slepian"] < 0)
befund("Der universelle Logarithmus -1/90 ist in diesen Daten nicht auflösbar: frei gefittete Logkoeffizienten streuen von "
       f"{min(clog['Gauss'] + clog['Slepian']):.1f} bis {max(clog['Gauss'] + clog['Slepian']):+.1f} (Knotensprünge, "
       "Entartung mit dem linearen Term); er bleibt die analytische Aussage aus Teil 1", True)

# --- 2.4 KMS-Gewichte (tau = i, Original, beta')
BETAS = {"2 pi (tau = i, q^m)": 2*np.pi, "4 pi (KMS Original)": 4*np.pi, "8 pi (beta' z83)": 8*np.pi}
kms = {}
for k_, b in BETAS.items():
    kms[k_] = {M: (S_gauss([np.pi/3], M, beta=b)[0] - S_gauss([np.pi/3], M)[0],
                   S_slepian([np.pi/3], M, beta=b)[0] - S_slepian([np.pi/3], M)[0]) for M in (32, 64, 96)}


def sB(x):
    if x <= 0 or x >= 1:
        return mp.mpf(0)
    nb = x/(1 - x)
    return (nb + 1)*mp.log(1 + nb) - nb*mp.log(nb)


S_th = {k: mp.fsum(m*m*sB(mp.exp(-mp.mpf(b)*m)) for m in range(1, 80)) for k, b in BETAS.items()}
S0_96 = S_gauss([np.pi/3], 96)[0]
ZAHLEN["kms_dS_kappe_pi/3"] = {k: {M: [f"{v[0]:.4e}", f"{v[1]:.4e}"] for M, v in d.items()} for k, d in kms.items()}
ZAHLEN["kms_S_thermisch_ganze_S3"] = {k: ns(v, 6) for k, v in S_th.items()}
konv = all(abs(d[96][i] - d[64][i]) < 0.05*abs(d[96][i]) for d in list(kms.values())[:2] for i in (0, 1))
lokal = all(abs(d[96][0]/d[96][1] - 1) < 0.03 for d in kms.values())
check("[E] KMS-Gewichte regulieren nicht: Delta S = S_beta - S_vak der Kappe theta = pi/3 ist UV-endlich (M = 64 -> 96 "
      "ändert < 5 %) und von der Lokalisierung unabhängig (Gauß/Slepian < 3 % bei M = 96): "
      + "; ".join(f"beta = {k}: {d[96][0]:.3e}" for k, d in kms.items())
      + f", gegen S_vak = {S0_96:.1f} bei M = 96. Der Flächenterm bleibt der der Abschneidung; die Gewichte e^{{-beta m}} "
      "geben nur einen endlichen thermischen Rest (ganze S^3: " + ", ".join(ns(v, 3) for v in S_th.values()) + ")",
      konv and lokal and all(abs(d[96][0]) < 1e-5*S0_96 for d in kms.values()))

# --- 2.5 native Schwellenabschneidung M(delta): q^M >= delta
schw = {}
for d in (1e-10, 1e-16, 1e-20, 1e-40):
    M = int(np.floor(np.log(1/d)/(2*np.pi)))
    sg = S_gauss(THETAS, M)
    schw[d] = (M, sg/(s2*M**2), float(sg[3]))
ZAHLEN["schwellenabschneidung"] = {f"{d:.0e}": {"M": v[0], "S/(sin^2 theta M^2), theta = pi/6, pi/4, pi/3, pi/2":
                                               [round(float(x), 4) for x in v[1]], "S(pi/2)": round(v[2], 4)}
                                   for d, v in schw.items()}
say("   Schwellenabschneidung q^M >= delta: " + "; ".join(f"delta = {d:.0e}: M = {v[0]}, S(pi/2) = {v[2]:.3f}, "
                                                        f"S/(s^2 M^2) = {[round(float(x), 3) for x in v[1]]}"
                                                        for d, v in schw.items()))
check("[E] native Schwelle q^M >= delta gibt M = ln(1/delta)/(2 pi) <= 14 auch bei delta = 1e-40; S(pi/2) wächst mit M(delta)^2 "
      f"({schw[1e-10][2]:.2f} bei M = 3 bis {schw[1e-40][2]:.1f} bei M = 14), und für M <= 7 streut S/(sin^2 theta M^2) über "
      "theta = pi/6..pi/2 um mehr als 10 % (kleine Kappen unter der Auflösung haben S = 0): kein festes S/A, die Zahl hängt an "
      "der willkürlichen Schwelle",
      all(v[0] <= 14 for v in schw.values()) and schw[1e-40][2]/schw[1e-10][2] > 10
      and all(np.ptp(v[1])/np.mean(v[1]) > 0.1 for v in schw.values() if v[0] <= 7))

# --- 2.6 Modulschnitt (Hislop-Longo nach dem Netz)
spek = {}
for N in (100, 400, 1600, 6400):
    k = np.arange(N - 1)
    ev = eigh_tridiagonal(np.zeros(N), np.sqrt((k + 1)*(k + 2))/2, eigvals_only=True)
    pos = np.sort(ev[ev > 0])
    spek[N] = (float(pos[1] - pos[0]), float((pos[1] - pos[0])*np.log(N)), float(ev.max()/N))
ZAHLEN["modulspektrum_D1"] = {N: {"Abstand bei 0": round(v[0], 4), "Abstand x ln N": round(v[1], 3), "max/N": round(v[2], 3)}
                             for N, v in spek.items()}
Ns = sorted(spek)
check("[E] Modulgenerator der Kappe (Hislop-Longo für das freie konforme Netz; hyperbolisches Element der radialen sl(2,R), "
      "die mit den Drehungen um den Kappenmittelpunkt vertauscht; s-Welle = D_1): die Abschneidungen von (L_1 + L_-1)/2 haben "
      "Eigenwertabstände bei 0 von " + ", ".join(f"{spek[N][0]:.3f}" for N in Ns) + " für N = 100..6400 (Abstand x ln N "
      "fast konstant, max/N -> 1): Spektrum R, kontinuierlich. Die diskreten q^m = e^{-2 pi m} sind Boltzmanngewichte von "
      "2 pi E (elliptisch, global), nicht des Modulhamiltonoperators der Kappe (Typ III_1, S unendlich)",
      all(spek[Ns[i + 1]][0] < spek[Ns[i]][0] for i in range(3)) and all(2.4 < spek[N][1] < 3.3 for N in Ns)
      and spek[6400][2] > 0.99)
dichte = mp.quad(lambda lam: lam**2/(2*PI**2)*sB(mp.exp(-2*PI*lam)), [0, 1, mp.inf])
dichte_3 = mp.quad(lambda lam: lam**2/(2*PI**2)*sB(mp.exp(-2*PI*lam)), [0, 1, 3])
s_mod_num = 1/(180*PI)
lim_SA = sp.limit(s_mod*V_H3/(4*sp.pi*sp.sinh(u_)**2), u_, sp.oo)
ZAHLEN["modulschnitt"] = {"s(1/2pi) aus der Plancherel-Dichte": ns(dichte, 20), "1/(180 pi)": ns(s_mod_num, 20),
                          "relativer Anteil der Modulenergien lambda > 3": ns(1 - dichte_3/dichte, 3),
                          "S/A_mod (u -> inf)": str(lim_SA), "a/pi": str(a_skalar/sp.pi)}
check("[E] die Divergenz sitzt nicht im Modulspektrum: die Entropiedichte bei T = 1/(2 pi) ist UV-endlich, int dlambda "
      "lambda^2/(2 pi^2) s_B(e^{-2 pi lambda}) = 1/(180 pi) (1e-30), Modulenergien lambda > 3 tragen nur "
      f"{ns(1 - dichte_3/dichte, 2)} davon (< 1e-4); divergent ist das H^3-Volumen am Schnitt. Ein Schnitt in Modulenergie "
      "(native Gewichte) reguliert daher nichts. Im Modulrahmen ist S/A_mod -> s/2 = 1/(360 pi) = a/pi fest (V/A -> 1/2 auf H^3)",
      abs(dichte - s_mod_num) < mp.mpf(10)**-30 and 0 < 1 - dichte_3/dichte < 1e-4
      and sp.simplify(lim_SA - a_skalar/sp.pi) == 0)
befund("Modulschnitt: Die einzige feste Zahl ist S/A_mod = a/pi, die Anomalie a im Modulrahmen. A_mod = 4 pi sinh^2 u_max ist "
       "die H^3-Fläche des Schnitts in Einheiten des Kappenradius, keine physikalische Fläche: bei nativ festem, "
       "dimensionslosem u_max ist S = (sinh 2u - 2u)/180 unabhängig von der Kappengröße (kein Flächengesetz); ein Flächengesetz "
       "entsteht erst mit u_max = ln(2R/eps), also wieder mit einer Länge eps (und der Konvention 2R/eps gegen R/eps: Faktor 4)",
       True)
def dez(x, fmt):
    return format(x, fmt).replace(".", ",")


C2G, C2S = dez(c2m["Gauss"], ".4f"), dez(c2m["Slepian"], ".4f")
PROZ = f"{100*(c2m['Slepian']/c2m['Gauss'] - 1):.0f} %"
KMS_T = ", ".join(dez(d[96][0], ".1e") for d in kms.values())
befund(f"Teil 2: Im rekonstruierten Leiterfeld gibt jede native Abschneidung m <= M ein Flächengesetz in M, aber mit einem "
       f"Koeffizienten, der an der Lokalisierung hängt (Gauß {C2G}, Slepian {C2S} je sin^2 theta M^2, {PROZ}), und mit "
       "unscharfer Nennlänge (M + delta). Die nativen Gewichte q^m, e^{-4 pi m}, e^{-8 pi m} sind keine UV-Regulatoren: als "
       "Länge liegen sie über dem Durchmesser der S^3, als KMS-Zustand ändern sie nur einen endlichen Rest, als Schwelle "
       "geben sie M <= 14 ohne Flächengesetz; der Modulschnitt fixiert nur a/pi in Modulfläche. S/A ist nicht fest", True)
say(f"   (Laufzeit bisher {time.time() - T0:.1f} s)")

# =================================================================== Teil 3
say("\nTeil 3: Vergleich mit 1/(4G) in TFPT-Normierung (GRAV.NONLINEAR.01: 8 pi G = 1/c3 = 8 pi)")
c3 = 1/(8*sp.pi)
Gs = sp.symbols("G", positive=True)
G_tfpt = sp.solve(sp.Eq(8*sp.pi*Gs, 1/c3), Gs)[0]
eta_tfpt = 1/(4*G_tfpt)
koeff = {"scharf/Q1/P2": sp.Integer(1)/(48*sp.pi), "Q2/Brown": sp.Integer(2)/(48*sp.pi),
         "P3": sp.Rational(1, 2)/(48*sp.pi), "Pauli-Villars": 2*sp.log(2)/(48*sp.pi), "CHM": 1/(360*sp.pi)}
eps2_c3 = {k: sp.simplify(v/eta_tfpt/c3) for k, v in koeff.items()}
ZAHLEN["noetiges_eps2_durch_c3_fuer_S/A=1/4"] = {k: str(v) for k, v in eps2_c3.items()}
werte = [float(v) for v in eps2_c3.values()]
check(f"Normierung: 8 pi G = 1/c3 gibt G = {G_tfpt}, 1/(4G) = {eta_tfpt} je v_geo-Fläche. S/A = 1/4 verlangt eps^2 = 4 G "
      "I[f] eps^2/(48 pi): eps^2/c3 = " + ", ".join(f"{v} ({k})" for k, v in eps2_c3.items()) + "; Streuung Faktor "
      f"{max(werte)/min(werte):.0f}, Pauli-Villars irrational (ln 2). Die nativen Abschneidungen geben eps in Einheiten des "
      "S^3-Radius, dessen Verhältnis zu v_geo das Repo nicht kennt: es gibt keinen regulatorfreien Vergleich",
      G_tfpt == 1 and eta_tfpt == sp.Rational(1, 4) and max(werte)/min(werte) > 10
      and sp.simplify(eps2_c3["Pauli-Villars"] - sp.Rational(4, 3)*sp.log(2)) == 0
      and eps2_c3["scharf/Q1/P2"] == sp.Rational(2, 3) and eps2_c3["CHM"] == sp.Rational(4, 45))
inv16piG = sp.Rational(1, 2)*(sp.Rational(1, 6) - xi)*I_/(16*sp.pi**2)          # Gamma_E = -(1/2) int ds/s f Tr K, a_1 = (1/6 - xi) R
inv4G = sp.simplify(4*sp.pi*inv16piG)
Sk = S_kegel(xi)
S_ent = S_kegel(0)                                                           # Hilbertraum: xi R phi^2 verschwindet am flachen Schnitt
check("[E] induzierte Kopplung gegen Entropie (Sakharov, Susskind-Uglum): 1/(4 G_ind) = (1 - 6 xi) I[f]/(48 pi) aus a_1 = "
      "(1/6 - xi) R, und die Kegelentropie mit Kontaktterm (int_Kegel R = 4 pi (1 - n) A) ist (1 - 6 xi) A I[f]/(48 pi): "
      "S_Kegel = A/(4 G_ind) für jedes Profil und jedes xi, der Regulator kürzt sich. Die Hilbertraum-Verschränkung ist am flachen "
      "Schnitt xi-unabhängig, A I[f]/(48 pi); sie gleicht A/(4 G_ind) nur bei xi = 0",
      sp.simplify(Sk - A_*inv4G) == 0 and sp.simplify(inv4G - (1 - 6*xi)*I_/(48*sp.pi)) == 0
      and sp.simplify(S_ent - A_*inv4G.subs(xi, 0)) == 0)
check("[E] beim nativen xi = 1/6 (erzwungen durch E = m) ist 1/G_ind = 0 und S_Kegel = 0 für jedes Profil, während die "
      f"Verschränkung der Kappe (Teil 2) einen Flächenterm c2 = {c2m['Gauss']:.4f} bzw. {c2m['Slepian']:.4f} trägt: "
      "S_ent/(A/4 G_ind) ist unendlich. Das Leiterfeld induziert keine Newtonkopplung; die regulatorfreie Gleichung "
      "S = A/(4G) gibt es für dieses Feld nicht",
      inv4G.subs(xi, xi_nativ[xi]) == 0 and Sk.subs(xi, xi_nativ[xi]) == 0 and min(c2m.values()) > 0.01)
check("pi-Kommensurabilität: 1/(48 pi) = c3/6 exakt, aber Q2 gibt c3/3, P3 c3/12, Pauli-Villars (ln 2/3) c3 und CHM "
      "(1/45) c3: die 'Relation' zu c3 folgt dem Profil, nicht dem Feld, und zählt nicht",
      sp.simplify(koeff["scharf/Q1/P2"] - c3/6) == 0 and sp.simplify(koeff["Q2/Brown"] - c3/3) == 0
      and sp.simplify(koeff["P3"] - c3/12) == 0 and sp.simplify(koeff["Pauli-Villars"] - sp.log(2)*c3/3) == 0
      and sp.simplify(koeff["CHM"] - c3/45) == 0)
befund("Dimensionstypisierung (DIMENSION.UPLIFT.FIREWALL.01, F2): Der Satz in GRAV.ENTROPY.EQUILIBRIUM.01 'der universelle "
       "Flächen-/Verschränkungskoeffizient einer CFT ist ihre Zentralladung (c = 8)' gilt nur in 1+1D (S = (c/3) ln L/eps, "
       "keine Fläche). In 3+1D ist der Flächenkoeffizient nicht universell (Teil 1: Profil, Teil 2: Lokalisierung); universell "
       "ist der Logkoeffizient -4a (a = 1/360 je reellem Skalar), und c = 8 ist die chirale Zentralladung der Naht. Der Satz ist "
       "für die 4D-Entropiedichte gescopt negativ (Ledger unverändert)", True)
befund("a/pi = 1/4 ist für freie Felder unmöglich (a ist rational, pi transzendent); der Vergleich wäre ohnehin sinnlos, weil "
       "A_mod keine physikalische Fläche ist", True)
befund("Was für G fehlt, genau: (i) eine native UV-Länge eps in Einheiten der Kappe und von l_P (die nativen Gewichte geben "
       "beta >= 2 pi, die Krümmungsskala); (ii) ein ausgezeichnetes Regulatorprofil bzw. eine ausgezeichnete Lokalisierung "
       f"(Faktor 2 bzw. {PROZ}); (iii) eine Bedeutung für 'S = A/4G': induziert geht das nur mit Materie mit a_1 != 0 im selben "
       "Regulator, das Leiterfeld hat a_1 = 0; (iv) alternativ eine native Zellzählung, etwa Nähte je Fläche mit ln 2/2 je Naht "
       f"(z86): dafür wäre die Dichte 1/(2 ln 2) = {ns(1/(2*mp.log(2)), 6)} je G-Fläche nötig, die nichts im Repo liefert", True)

# =================================================================== Teil 4
say("\nTeil 4: Urteil für F2")
STATUS = {
    "Regulatorabhängigkeit des Flächenkoeffizienten": ("E", "S = A I[f]/(48 pi); Faktor 2 bei gleicher Nennlänge (Q1/Q2, "
                                                           "P2/P3, Brownsche Konvention), Pauli-Villars 2 ln 2, CHM 2/15, "
                                                           "zeta 0: externe Analyse bestätigt"),
    "universeller Teil": ("E", "Logkoeffizient -1/90 = -4a, a = 1/360 (zeta(0) auf S^4 doppelt, CHM); in der nativen "
                               "Abschneidung numerisch nicht aufgelöst"),
    "Flächengesetz in nativer Abschneidung m <= M": ("E/C", f"c2 = {C2G} (Gauß) bzw. {C2S} (Slepian) je sin^2 theta M^2, "
                                                            "theta-unabhängig < 2 %, lokalisierungsabhängig "
                                                            f"{PROZ}; bedingt durch S2 (z83, C) und Fockquantisierung"),
    "native Gewichte q^m, KMS 4 pi, beta' 8 pi als Regulator": ("E (negativ)", "Länge > Durchmesser der S^3; KMS-Rest "
                                                                               f"UV-endlich ({KMS_T}); "
                                                                               "Schwelle M <= 14 ohne Flächengesetz"),
    "Modulschnitt": ("E", "Modulspektrum kontinuierlich (Typ III_1), Dichte 1/(180 pi) endlich, Divergenz im H^3-Volumen; "
                          "S/A_mod = a/pi fest, A_mod keine physikalische Fläche"),
    "induzierte Newtonkopplung des Leiterfelds": ("E (negativ)", "xi = 1/6 aus E = m, a_1 = 0, 1/G_ind = 0 und S_Kegel = 0 "
                                                                 "für jedes Profil; S_ent != 0"),
    "Vergleich mit 1/(4G) = 1/4": ("O", "nötiges eps^2/c3 in {4/45, 1/3, 2/3, (4/3) ln 2, 4/3}: ohne native Länge keine "
                                        "Bedeutung; 1/(48 pi) = c3/6 ist pi-Kommensurabilität"),
    "Ledger-Satz 'Flächenkoeffizient = Zentralladung'": ("E (negativ, 4D)", "nur 1+1D; in 4D universell ist a, nicht die "
                                                                             "Fläche"),
    "G über einen räumlichen Schnitt": ("O", "braucht native UV-Länge, Profil und eine Bedeutung (a_1 != 0 im selben "
                                             "Regulator oder native Zellzählung)"),
}
for k_, (s_, t_) in STATUS.items():
    say(f"   [{s_}] {k_}: {t_}")
ZAHLEN["status"] = STATUS
befund("F2 gegen seine Frage: Ist S/A mit den nativen Gewichten als Regulator fest und unabhängig von den Details? Nein, "
       "zweifach: die nativen Gewichte regulieren das UV nicht (Krümmungsskala), und jede native Abschneidung, die es tut, "
       "hat einen lokalisierungsabhängigen Koeffizienten. Der Weg zu G über die Kappe bleibt offen; die Bedingung dafür ist "
       "jetzt scharf: eine native Länge plus ein ausgezeichnetes Profil plus eine Bedeutung", True)
urteil = ("F2 ist in der gestellten Form gescopt negativ. Die zehnte externe Analyse stimmt: Für den freien konformen Skalar ist "
          "der Flächenkoeffizient I[f]/(48 pi), ein lineares Funktional des Regulatorprofils (Faktor 2 bei gleicher Nennlänge, "
          "Pauli-Villars 2 ln 2, CHM-Schnitt 2/15 des scharfen Werts, zeta-Regularisierung 0); regulatorfrei ist nur der "
          "Logarithmus -1/90 = -4a, dreifach bestätigt (zeta(0) des konformen Laplace auf S^4 über Hurwitz und a_2, "
          "CHM-Entropiedichte 1/(180 pi) auf H^3 mal Logterm des Volumens). Im Leiterfeld von z83 (E = m erzwingt xi = 1/6) gibt "
          "jede native Abschneidung m <= M ein echtes Flächengesetz S = c2 sin^2 theta M^2 + ..., theta-unabhängig auf 2 %, "
          f"aber c2 hängt an der Lokalisierung derselben Abschneidung: {C2G} mit den Gauß-Gegenbauer-Knoten (exakte Isometrie), "
          f"{C2S} mit der Slepian-Teilung ({PROZ}), und die Nennlänge ist selbst unscharf (M + delta, delta = "
          f"{dez(delta_eff['Gauss'], '+.1f')} bzw. {dez(delta_eff['Slepian'], '+.1f')}). "
          "Die nativen Gewichte sind keine UV-Regulatoren: als Zeitspaltung eps = 2 pi, 4 pi, 8 pi liegen sie über dem "
          "Durchmesser pi der S^3 (der regulierte Korrelator ist auf 0,75 % flach), als KMS-Zustände bei tau = i, beta = 4 pi und "
          f"beta' = 8 pi ändern sie die Kappe nur um einen UV-endlichen, lokalisierungsfreien Rest ({KMS_T}), als "
          "Schwelle q^M >= delta geben sie M <= 14 ohne Flächengesetz. Der Modulschnitt (Hislop-Longo nach dem Netz, "
          "GRAV.NONCIRCULAR.01 eingehalten) zeigt den Grund: das Modulspektrum der Kappe ist kontinuierlich (Typ III_1), die "
          "Entropiedichte bei T = 1/(2 pi) ist endlich, die Divergenz sitzt im H^3-Volumen am Schnitt, also im Ort, nicht in der "
          "Energie; fest ist dort nur S/A_mod = a/pi in Modulfläche, die keine physikalische Fläche ist. Für G heißt das: Ohne "
          "native UV-Länge und ausgezeichnetes Profil gibt es keinen Vergleich mit 1/(4G) = 1/4 (8 pi G = 1/c3); die nötigen "
          "eps^2/c3 streuen um den Faktor 15, und 1/(48 pi) = c3/6 ist bloße pi-Kommensurabilität. Der Sakharov-Weg, auf dem "
          "sich der Regulator kürzt (S_Kegel = A/(4 G_ind) für jedes Profil), ist für das Leiterfeld selbst verschlossen: bei "
          "xi = 1/6 ist a_1 = 0, also 1/G_ind = 0, während die Kappe verschränkt ist. Der Ledger-Satz 'universeller "
          "Flächenkoeffizient = Zentralladung c = 8' ist 1+1D-typisiert und für die 4D-Entropiedichte gescopt negativ. Offen und "
          "jetzt scharf gestellt: G über einen räumlichen Schnitt braucht eine native UV-Länge, ein ausgezeichnetes Profil und "
          "eine Bedeutung (Materie mit a_1 != 0 im selben Regulator oder eine native Zellzählung, z. B. Nähte je Fläche mit "
          "Dichte 1/(2 ln 2)).")
say("\nUrteil\n   " + urteil)
say(f"\nLaufzeit {time.time() - T0:.1f} s")
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z90_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
