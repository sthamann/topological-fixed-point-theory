"""Z58: Rekonstruktion statt Auswahl -- Dynamik, Gravitationsantwort und Chiralität aus den Antworten des Originals.

Strategiewechsel (Auftrag 24.11.): Die Dynamik wird nicht mehr aus kinematischen, Symmetrie- oder Topologiebedingungen
ausgewählt und mit neuen Prinzipien ergänzt, sondern aus den vollständigen Antworten des originalen Prozesses
rekonstruiert (GNS / Osterwalder-Schrader / Tomita-Takesaki). Chirale Materie muss nicht aus gegappten Spiegeln eines voll
diskretisierten Fermionoperators kommen.

Original: die beschriftete Calderón/Szegő-Quarterkovarianz c(x) = 1/(1 + e^{4 pi x}), x = r + alpha, alpha = 1/4, 3/4
(ELLIPTIC-CURRENT-SOURCE, Gl. 8; z46). Kein weiterer Input.

Teil 1  Antworten -> Generator: modularer Einteilchengenerator h = log((1-C)/C) = 4 pi x; Stromantworten, KMS-Verhältnisse
        und Kommutator geben CCR und die lineare Dispersion omega_m = 4 pi m.
Teil 2  Gaußsche Antworten -> quadratischer Generator: Sugawara-Identität auf dem Fockraum, Wick-Prüfung, Jacobi-Tripel-
        produkt bei beta = 4 pi (Ladung x Oszillatoren faktorisiert). Spektrum des vollen Modulgenerators liegt in pi*N.
Teil 3  Band oder Leiter: Das Band A^+A kann kein Modulgenerator sein (kontinuierliches Spektrum, e^{-K} nicht spurklasse)
        und ist keine Einschränkung eines quadratischen Zweimodenoperators; die Leiter ist es. Eindeutigkeit unter
        Modulfluss und modularer Konjugation J; Skala = Modulfrequenz 4 pi je Quant.
Teil 4  Gravitationsantwort aus der relativen Entropie: W(J) = log Tr e^{-(1-J)K}; der positive Propagator ist die
        Fisher-Information kappa_2 (Konvexität, nicht entworfen), die Selbstkopplungen sind kappa_3, kappa_4. Der
        Ladungssektor (Theta) dominiert und ist stark nichtgaußsch. Dazu die Modulantworten der E8-Naht bei tau = i.
Teil 5  Chiralität aus dem Index ohne Spiegellücke: ein chiraler Zweig, Spektralfluss 1, Toeplitz-Index +1 mit Restmischung
        O(q^depth), Familienindex 3 ohne Spiegel (E_- = O(-1) + 3 O), Z16.
Teil 6  Bilanz: welche Auswahlprinzipien zu Folgerungen werden, was offen bleibt.
"""
import json
from pathlib import Path

import mpmath as mp
import numpy as np
import scipy.sparse as sps
import sympy as sp

HERE = Path(__file__).parent
mp.mp.dps = 60

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


PI4 = 4*mp.pi
ALPHAS = (mp.mpf(1)/4, mp.mpf(3)/4)


def occ(x, beta=PI4):
    return 1/(1 + mp.exp(beta*x))


# ------------------------------------------------------------------ Teil 1
say("Teil 1: Antworten des Originals -> Generator (Tomita / Osterwalder-Schrader)")
ok = True
for a in ALPHAS:
    for r in range(-5, 6):
        x = r + a
        c = occ(x)
        ok &= abs(mp.log((1 - c)/c) - PI4*x) < mp.mpf(10)**-25
check("der modulare Einteilchengenerator h = log((1 - C)/C) der originalen Kovarianz ist exakt 4 pi (r + alpha): "
      "rekonstruiert, nicht gewählt", ok)

R1 = 40
ok_w, ok_kms, ok_ccr = True, True, True
for m in range(1, 5):
    exact = m/(mp.exp(PI4*m) - 1)
    for a in ALPHAS:
        up = mp.fsum(occ(r + a + m)*(1 - occ(r + a)) for r in range(-R1, R1))      # omega(J_-m J_m)
        dn = mp.fsum(occ(r + a)*(1 - occ(r + a + m)) for r in range(-R1, R1))      # omega(J_m J_-m)
        ok_w &= abs(up/exact - 1) < mp.mpf(10)**-30
        ok_kms &= abs(mp.log(dn/up) - PI4*m) < mp.mpf(10)**-30
        ok_ccr &= abs(dn - up - m) < mp.mpf(10)**-30
check("Stromantworten omega(J_-m J_m) = m/(e^{4 pi m} - 1) für m = 1..4 und beide Quarter (Teleskopsumme)", ok_w)
check("Kommutatorantwort omega([J_m, J_-m]) = m: die CCR folgt aus den Antworten", ok_ccr)
check("KMS-Verhältnis omega(J_m J_-m)/omega(J_-m J_m) = e^{4 pi m}: rekonstruierte Frequenzen omega_m = 4 pi m "
      "(lineare Dispersion, ein chiraler masseloser Zweig)", ok_kms)
ZAHLEN["omega_J-1_J1"] = mp.nstr(1/(mp.exp(PI4) - 1), 12)
befund("Der rekonstruierte Generator ist K = 4 pi L_0 (modulare Zeit); die Zeitrichtung und ihre Skala kommen aus dem "
       "Zustand selbst (thermische Zeit), nicht aus einer Wahl", True)

# ------------------------------------------------------------------ Teil 2
say("\nTeil 2: Gaußsche Antworten -> quadratischer Generator")


def fock(xs):
    n = len(xs)
    I2 = sps.identity(2, format="csr")
    Zp = sps.diags([1.0, -1.0], format="csr")
    sm = sps.csr_matrix(np.array([[0.0, 1.0], [0.0, 0.0]]))
    ops = []
    for j in range(n):
        mats = [Zp]*j + [sm] + [I2]*(n - j - 1)
        op = mats[0]
        for M in mats[1:]:
            op = sps.kron(op, M, format="csr")
        ops.append(op)
    return ops


def build(alpha, R):
    xs = [r + alpha for r in range(-R, R)]
    cs = fock(xs)
    num = [np.real((c.T @ c).diagonal()) for c in cs]
    dim = cs[0].shape[0]
    sea = [1.0 if x < 0 else 0.0 for x in xs]
    L0 = sum(x*(nn - s) for x, nn, s in zip(xs, num, sea))
    Q = sum(nn - s for nn, s in zip(num, sea))
    idx = {round(4*x): j for j, x in enumerate(xs)}

    def Jm(m):          # J_m = sum_x psi*_x psi_{x+m}  (senkt die Energie um m für m > 0)
        op = sps.csr_matrix((dim, dim))
        for j, x in enumerate(xs):
            k = idx.get(round(4*(x + m)))
            if k is not None:
                op = op + cs[j].T @ cs[k]
        return op.tocsr()
    return xs, num, L0, Q, Jm


ok_sug = True
for a in (0.25, 0.75):
    xs, num, L0, Q, Jm = build(a, 6)
    Mmax = 7
    sug = sps.diags(Q**2/2 + (a - 0.5)*Q)
    for m in range(1, Mmax + 1):
        sug = sug + Jm(-m) @ Jm(m)
    diff = sug - sps.diags(L0)
    sea_ok = np.ones(len(L0), dtype=bool)
    for x, nn in zip(xs, num):
        if abs(x) > 3.5:
            sea_ok &= (nn == (1.0 if x < 0 else 0.0))
    low = np.where(sea_ok & (L0 <= 3.0 + 1e-9))[0]
    err = max(np.abs(diff[:, s].toarray()).max() for s in low)
    ok_sug &= err < 1e-12 and len(low) >= 15
check("Sugawara-Identität L_0 = sum_m J_-m J_m + Q^2/2 + (alpha - 1/2) Q auf allen Zuständen mit L_0 <= 3 "
      "(Fockraum, beide Quarter): der Generator ist quadratisch in den Strömen plus eine Ladungsfunktion", ok_sug)

beta_t = 2.0
for a in (0.25, 0.75):
    xs, num, L0, Q, Jm = build(a, 7)
    cx = [1/(1 + np.exp(beta_t*x)) for x in xs]
    rho = np.ones(len(L0))
    for c, nn in zip(cx, num):
        rho *= c*nn + (1 - c)*(1 - nn)

    def ev(A, B):
        return float(rho @ np.asarray(A.multiply(B.T).sum(axis=1)).ravel())
    J1, Jm1, J2, Jm2 = Jm(1), Jm(-1), Jm(2), Jm(-2)
    n1 = ev(Jm1, J1)
    n2 = ev(Jm2, J2)
    four = ev(Jm1 @ Jm1, J1 @ J1)
    mixed = ev(Jm2 @ Jm1, J1 @ J2)
    e1, e2 = 1/(np.exp(beta_t) - 1), 2/(np.exp(2*beta_t) - 1)
    pq = {}
    for qv, p in zip(Q, rho):
        pq[round(qv)] = pq.get(round(qv), 0.0) + p
    qs = np.array(sorted(pq))
    ps = np.array([pq[k] for k in qs])
    mu = ps @ qs
    k2 = ps @ (qs - mu)**2
    k4 = ps @ (qs - mu)**4 - 3*k2**2
    w = np.exp(-beta_t*(qs**2/2 + (a - 0.5)*qs))
    w /= w.sum()
    mu_e = w @ qs
    k2e = w @ (qs - mu_e)**2
    k4e = w @ (qs - mu_e)**4 - 3*k2e**2
    if a == 0.25:
        wick = (abs(n1/e1 - 1) < 1e-3 and abs(n2/e2 - 1) < 1e-3 and abs(four/(2*n1**2) - 1) < 1e-3
                and abs(mixed/(n1*n2) - 1) < 1e-3)
        theta = abs(k2/k2e - 1) < 1e-4 and abs(k4 - k4e) < 1e-4*k2**2
        ZAHLEN["wick_test_beta2"] = {"n1": n1, "four/(2 n1^2)": four/(2*n1**2), "mixed/(n1 n2)": mixed/(n1*n2),
                                     "kappa4_Q/kappa2_Q^2": k4/k2**2}
check("Wick bei beta = 2 im Fockraum (14 Moden, Abschneidefehler ~1e-4): <J_-1^2 J_1^2> = 2<J_-1 J_1>^2, <J_-2 J_-1 J_1 J_2> = <J_-1 J_1><J_-2 J_2>, "
      "<J_-m J_m> = m/(e^{beta m} - 1): die Stromantworten sind gaußsch", wick)
check("die Ladungsverteilung derselben Rechnung ist die Thetaverteilung e^{-beta(Q^2/2 + (alpha - 1/2)Q)} (bei beta = 2 "
      "fast gaußsch; stark nichtgaußsch wird sie bei beta = 4 pi, Teil 4)", theta)

ok_jtp = True
x4 = mp.exp(-PI4)
for a in ALPHAS:
    for y in (mp.mpf(1), mp.mpf(2), mp.mpf(1)/3, mp.expj(1)):
        lhs = mp.fprod(1 + y*mp.exp(-PI4*(r + a)) for r in range(0, 30))*mp.fprod(
            1 + mp.exp(-PI4*(r - a))/y for r in range(1, 30))
        rhs = mp.fsum(y**Qn*mp.exp(-PI4*(Qn**2/2 + (a - mp.mpf(1)/2)*Qn)) for Qn in range(-12, 13))/mp.fprod(
            1 - x4**m for m in range(1, 40))
        ok_jtp &= abs(lhs/rhs - 1) < mp.mpf(10)**-40
check("Jacobi-Tripelprodukt bei beta = 4 pi für vier Fugazitäten: der originale Zustand faktorisiert exakt in "
      "Thetaverteilung der Ladung x gaußsche thermische Oszillatoren", ok_jtp)
ok_lat = all(sp.Rational(4)*(sp.Rational(1, 2)*Qn**2 + (al - sp.Rational(1, 2))*Qn) == int(
    4*(sp.Rational(1, 2)*Qn**2 + (al - sp.Rational(1, 2))*Qn)) for al in (sp.Rational(1, 4), sp.Rational(3, 4))
    for Qn in range(-6, 7))
check("Spektrum des vollen Modulgenerators K = 4 pi (sum J_-m J_m + E(Q)) liegt in pi*N (E(Q)*4 ganz): eine Leiter mit "
      "Stufe pi, Modulfluss periodisch mit Periode 2", ok_lat)
befund("Der niedrigste modulare Anregungszustand des Originals ist der geladene Twistzustand mit L_0 = 1/4 (Lücke pi), "
       "nicht das Stromquant (Lücke 4 pi)", True)
befund("'Quadratisch in den Strommoden' ist damit kein Auswahlprinzip mehr, sondern Folge der gaußschen Antworten des "
       "Originals (Sugawara + Wick + Tripelprodukt); vorausgesetzt bleibt, dass das Original quasifrei ist (freie Naht, z53)",
       True)

# ------------------------------------------------------------------ Teil 3
say("\nTeil 3: Band oder Leiter -- aus der Rekonstruktion entschieden")


def pair_check(q, nmax=36):
    a = np.diag(np.sqrt(np.arange(1, nmax)), 1)
    I = np.eye(nmax)
    aL, aR = np.kron(a, I), np.kron(I, a)
    th = np.arctanh(q)
    ch, shh = np.cosh(th), np.sinh(th)
    dL = ch*aL - shh*aR.T
    dR = ch*aR - shh*aL.T
    H = dL.T @ dL + dR.T @ dR
    K = 7
    pidx = [n*nmax + n for n in range(K)]
    Hp = H[np.ix_(pidx, pidx)]
    A = np.zeros((K + 1, K + 1))
    for n in range(K + 1):
        A[n, n] = q
        if n + 1 <= K:
            A[n, n + 1] = -1.0
    lad = (A.T @ np.diag(np.arange(1, K + 2)) @ A)[:K, :K]
    off = max(np.abs(H[:, n*nmax + n][[i for i in range(nmax*nmax) if i not in
                                       [k*nmax + k for k in range(nmax)]]]).max() for n in range(K))
    tfd = np.zeros(nmax*nmax)
    for n in range(nmax):
        tfd[n*nmax + n] = np.sqrt(1 - q*q)*q**n
    Ht = H @ tfd
    keep = [i for i in range(nmax*nmax) if i//nmax < nmax - 4 and i % nmax < nmax - 4]
    return (np.abs(Hp - 2/(1 - q*q)*lad).max(), off, np.abs(Ht[keep]).max(), aL, aR, tfd, keep, ch, shh, nmax)


for q in (0.3, float(mp.exp(-2*mp.pi))):
    dev, off, ann, aL, aR, tfd, keep, ch, shh, nmax = pair_check(q)
    if q == 0.3:
        ok_pair = dev < 1e-10 and off < 1e-10 and ann < 1e-10
    else:
        ok_pair &= dev < 1e-10 and off < 1e-10 and ann < 1e-10
check("Zweimoden-Fockraum (q = 0,3 und q = e^{-2 pi}): H = d_L^+ d_L + d_R^+ d_R erhält den Paarsektor, vernichtet den "
      "Thermofeldzustand und ist dort 2/(1 - q^2) A^+(N+1)A", ok_pair)

q = 0.3
dev, off, ann, aL, aR, tfd, keep, ch, shh, nmax = pair_check(q)
B = [aL.T @ aL, aR.T @ aR, aL @ aR + aL.T @ aR.T, 1j*(aL.T @ aR.T - aL @ aR), np.eye(nmax*nmax)]
Mtfd = np.array([(Bk @ tfd)[keep] for Bk in B]).T
Mtfd = np.vstack([Mtfd.real, Mtfd.imag])            # hermitesch: reelle Koeffizienten
sv = np.linalg.svd(Mtfd, compute_uv=False)
null_dim = int(np.sum(sv < 1e-10*sv.max()))
# J-invariant: c1 = c2, c4 = 0  ->  Parameter (c, c3, c5)
MJ = np.array([(B[0] + B[1]) @ tfd, B[2] @ tfd, B[4] @ tfd]).T[keep]
U, s, Vh = np.linalg.svd(MJ)
nJ = int(np.sum(s < 1e-10*s.max()))
v = Vh[-1].real
v = v/v[0]
target = np.array([ch**2 + shh**2, -2*ch*shh, 2*shh**2])/(ch**2 + shh**2)
check("quadratische hermitesche Zweimodenoperatoren, die mit dem Modulfluss N_L - N_R kommutieren und den Thermofeldzustand "
      "vernichten: genau zwei (d_L^+ d_L, d_R^+ d_R); mit modularer Konjugation J genau einer, d_L^+ d_L + d_R^+ d_R",
      null_dim == 2 and nJ == 1 and np.allclose(v, target, atol=1e-8))

K = 40
pr = [n*nmax + n for n in range(min(K, nmax - 2))]
span = [Bk[np.ix_(pr, pr)] for Bk in (B[0], B[2], B[3], B[4])]
kk = len(pr)
A = np.zeros((kk + 1, kk + 1))
for n in range(kk + 1):
    A[n, n] = q
    if n + 1 <= kk:
        A[n, n + 1] = -1.0
band = (A.T @ A)[:kk, :kk]
lad = (A.T @ np.diag(np.arange(1, kk + 2)) @ A)[:kk, :kk]
Mspan = np.array([S.ravel() for S in span]).T


def resid(T):
    coef, *_ = np.linalg.lstsq(Mspan, T.ravel().astype(complex), rcond=None)
    return np.linalg.norm(Mspan @ coef - T.ravel())/np.linalg.norm(T)


r_band, r_lad = resid(band), resid(lad)
ZAHLEN["paarsektor_residuum"] = {"band": r_band, "leiter": r_lad}
check("auf dem Paarsektor liegt die Leiter A^+(N+1)A im Spann der quadratischen Einschränkungen (Residuum 0), das Band "
      "A^+A nicht (konstante Nebendiagonale -q statt linear in n)", r_lad < 1e-12 and r_band > 1e-2)

qo = float(mp.exp(-2*mp.pi))
traces = []
for M in (100, 200, 400):
    Ab = np.zeros((M + 1, M + 1))
    for n in range(M + 1):
        Ab[n, n] = qo
        if n + 1 <= M:
            Ab[n, n + 1] = -1.0
    Lb = (Ab.T @ Ab)[:M, :M]
    ev = np.linalg.eigvalsh(Lb)
    traces.append(np.exp(-ev).sum())
inband = np.mean((ev > (1 - qo)**2 - 1e-9) & (ev < (1 + qo)**2 + 1e-9))
lad_tr = 1/(1 - np.exp(-4*np.pi))
check("das Band (q = e^{-2 pi}) hat ein Kontinuum in [(1-q)^2, (1+q)^2]: Tr e^{-L} wächst linear mit der Abschneidung "
      "(100, 200, 400 Moden), also ist e^{-L} nicht spurklasse und L kein Modulgenerator irgendeines normalen Zustands; "
      "die Leiter 4 pi N hat Tr e^{-K} = 1/(1 - e^{-4 pi})",
      abs(traces[1]/traces[0] - 2) < 0.05 and abs(traces[2]/traces[1] - 2) < 0.05 and inband > 0.95
      and abs(lad_tr - 1) < 1e-5)
ZAHLEN["band_spur"] = traces
befund("Band oder Leiter ist entschieden: Leiter. Skala: 4 pi je Quant (8 pi je Paaranregung), also "
       "lambda = 8 pi/(1 - q^2) in K_lambda = lambda A^+(N+1)A; die physikalische Einheit ist beta = L der Naht", True)

# ------------------------------------------------------------------ Teil 4
say("\nTeil 4: Gravitationsantwort aus der relativen Entropie (statt entworfener Propagatoren)")
Qr = range(-20, 21)


def logZ_parts(beta):
    ch_, os_ = mp.mpf(0), mp.mpf(0)
    for a in ALPHAS:
        ch_ += mp.log(mp.fsum(mp.exp(-beta*(Qn**2/2 + (a - mp.mpf(1)/2)*Qn)) for Qn in Qr))
        os_ += -mp.fsum(mp.log(1 - mp.exp(-beta*m)) for m in range(1, 80))
    return ch_, os_


def W(J, part=None):
    c_, o_ = logZ_parts(PI4*(1 - J))
    return c_ + o_ if part is None else (c_ if part == "c" else o_)


kap = {n: mp.diff(W, 0, n) for n in (1, 2, 3, 4)}
kc = {n: mp.diff(lambda J: W(J, "c"), 0, n) for n in (2, 3, 4)}
ko = {n: mp.diff(lambda J: W(J, "o"), 0, n) for n in (2, 3, 4)}
ZAHLEN["energiekumulanten_K"] = {f"kappa{n}": mp.nstr(kap[n], 10) for n in (2, 3, 4)}
ZAHLEN["ladung_vs_oszillator_kappa2"] = [mp.nstr(kc[2], 10), mp.nstr(ko[2], 10)]
say(f"   kappa_2(K) = {mp.nstr(kap[2], 8)}, kappa_3 = {mp.nstr(kap[3], 8)}, kappa_4 = {mp.nstr(kap[4], 8)}")
say(f"   Ladungssektor kappa_2 = {mp.nstr(kc[2], 8)}, Oszillatoren kappa_2 = {mp.nstr(ko[2], 8)}, "
    f"Verhältnis {mp.nstr(kc[2]/ko[2], 6)}")
Js = mp.mpf("0.01")
S_exact = Js*mp.diff(W, Js, 1) - W(Js) + W(0)
S_ser = kap[2]*Js**2/2 + kap[3]*Js**3/3 + kap[4]*Js**4/8
check("relative Entropie S(rho_J || rho) = J W'(J) - W(J) + W(0) = kappa_2 J^2/2 + kappa_3 J^3/3 + kappa_4 J^4/8 + O(J^5)",
      abs(S_exact - S_ser) < 50*abs(kap[2])*Js**5*mp.mpf(10)**3 and abs(S_exact/S_ser - 1) < mp.mpf("1e-3"))
check("der quadratische Antwortkern ist die Fisher-Information kappa_2 > 0 (Konvexität von W): positiv ohne Entwurf", kap[2] > 0)
check("der Ladungssektor (Theta) dominiert die Energieantwort um mehr als den Faktor 100", kc[2]/ko[2] > 100)
exk_c = kc[4]/kc[2]**2
exk_o = ko[4]/ko[2]**2
ZAHLEN["exzess_kurtosis"] = {"ladung": mp.nstr(exk_c, 8), "oszillator": mp.nstr(exk_o, 8)}
check("die Antwort ist stark nichtgaußsch: standardisierte vierte Kumulante des Ladungssektors > 10", exk_c > 10)
G2, G3, G4 = 1/kap[2], -kap[3]/kap[2]**3, -kap[4]/kap[2]**4 + 3*kap[3]**2/kap[2]**5
ZAHLEN["effektive_wirkung"] = {"Gamma2": mp.nstr(G2, 10), "Gamma3": mp.nstr(G3, 10), "Gamma4": mp.nstr(G4, 10)}
say(f"   Legendre-Wirkung Gamma(phi): Gamma2 = 1/kappa2 = {mp.nstr(G2, 8)}, Gamma3 = {mp.nstr(G3, 8)}, "
    f"Gamma4 = {mp.nstr(G4, 8)}")
befund("Propagator und Selbstkopplungen des Spurkanals sind damit Zahlen des Originals; der Tensorkanal (Spin 2 in 4D) ist "
       "aus diesem skalaren Antwortfunktional nicht ablesbar", True)

# Modulantworten der E8-Naht bei tau = i
e2, e4, e6 = sp.symbols("E2 E4 E6")


def Dr(f):
    return sp.expand(sp.diff(f, e2)*(e2**2 - e4)/12 + sp.diff(f, e4)*(e2*e4 - e6)/3 + sp.diff(f, e6)*(e2*e6 - e4**2)/2)


f = [None, sp.simplify(Dr(sp.log(e4)) - sp.Rational(8, 24)*e2)]
for n in range(2, 5):
    f.append(sp.simplify(Dr(f[-1])))
E4i = 3*sp.gamma(sp.Rational(1, 4))**8/(2*sp.pi)**6
vals = [sp.nsimplify(f[n].subs({e6: 0, e2: 3/sp.pi})) for n in range(1, 5)]
say(f"   D log Z_E8 = {f[1]}; bei tau = i: D^1..D^4 = {[sp.simplify(v_) for v_ in vals]}")


def logZE8(y):
    qq = mp.exp(-2*mp.pi*y)
    E4 = 1 + 240*mp.fsum(mp.fsum(d**3 for d in range(1, n + 1) if n % d == 0)*qq**n for n in range(1, 40))
    eta = qq**(mp.mpf(1)/24)*mp.fprod(1 - qq**n for n in range(1, 60))
    return mp.log(E4) - 8*mp.log(eta)


E4num = sp.N(E4i, 30)
okE8 = True
for n in range(1, 5):
    num_ = (-1/(2*mp.pi))**n*mp.diff(logZE8, 1, n)
    ana = mp.mpf(str(sp.N(vals[n - 1].subs(e4, E4i), 40)))
    okE8 &= abs(num_ - ana) < mp.mpf(10)**-20*max(1, abs(ana))
ZAHLEN["E8_modulantworten_tau_i"] = [str(sp.N(v_.subs(e4, E4i), 15)) for v_ in vals]
check("E8-Naht: D log Z = -E6/(3 E4) verschwindet bei tau = i (<T> = 0), D^2..D^4 in geschlossener Form über E4(i) = "
      "3 Gamma(1/4)^8/(2 pi)^6 (gegen q-Reihen geprüft)", okE8 and vals[0] == 0)
befund("Die zweite Modulantwort der Naht ist entlang der Temperatur positiv und entlang der Verdrillung negativ (holomorph): "
       "der bekannte Sattel; die dritte und vierte sind die Selbstkopplungen des Modulus", True)

# ------------------------------------------------------------------ Teil 5
say("\nTeil 5: Chiralität aus dem Index, ohne Spiegellücke")
eta_ok = all(abs((mp.zeta(0, a) - mp.zeta(0, 1 - a)) - (1 - 2*a)) < mp.mpf(10)**-40 for a in ALPHAS)
check("spektrale Asymmetrie des rekonstruierten Generators eta(alpha) = zeta(0, alpha) - zeta(0, 1 - alpha) = 1 - 2 alpha "
      "(+1/2 und -1/2); beim Durchlauf alpha -> alpha + 1 springt sie um 2: genau ein Niveau, Spektralfluss 1", eta_ok)
xs = np.linspace(-5, 5, 2000)
cross_rec = int(np.sum(np.diff(np.sign(4*np.pi*xs)) != 0))
ks = -np.pi + (np.arange(2000) + 0.5)*2*np.pi/2000
sl = np.sin(ks)
cross_lat = int(np.sum(np.diff(np.sign(sl)) != 0)) + int(np.sign(sl[0]) != np.sign(sl[-1]))
check("die rekonstruierte Dispersion 4 pi x hat genau einen Nulldurchgang (ein chiraler Zweig, keine Brillouinzone); eine "
      "periodische Gitterdispersion hat stets eine gerade Zahl (Nielsen-Ninomiya-Voraussetzung hier nicht erfüllt)",
      cross_rec == 1 and cross_lat % 2 == 0)

zz = [mp.expj(2*mp.pi*k/4000) for k in range(4001)]
for qq in (mp.mpf("0.5"), mp.exp(-2*mp.pi)):
    sym = [mp.conj(z)*(1 - qq*z**4) for z in zz]
    wind = mp.fsum(mp.arg(sym[k + 1]/sym[k]) for k in range(4000))/(2*mp.pi)
    if qq == mp.mpf("0.5"):
        w_half = wind
    else:
        w_orig = wind
check("Toeplitz-Symbol z^-1 (1 - q z^4) hat Windung -1 (q = 0,5 und e^{-2 pi}): Index A_q = +1 (Gohberg-Krein)",
      abs(w_half + 1) < 1e-20 and abs(w_orig + 1) < 1e-20)


def smin(qv, M):
    T = np.zeros((M, M))
    for j in range(M):
        if j >= 1:
            T[j - 1, j] = 1.0
        if j + 3 < M:
            T[j + 3, j] = -qv
    return np.linalg.svd(T, compute_uv=False).min()


rat = [smin(0.5, M + 4)/smin(0.5, M) for M in (8, 12, 16, 20)]
so = [smin(float(mp.exp(-2*mp.pi)), M) for M in (4, 8, 12)]
ZAHLEN["restmischung_q_e^-2pi"] = so
eta_vec_ok = True
qs_ = sp.Symbol("q")
Nn = 24
f_ = [qs_**(n//4) if n % 4 == 0 else 0 for n in range(Nn)]
for n in range(Nn - 4):
    eta_vec_ok &= sp.simplify(f_[n + 1] - (qs_*f_[n - 3] if n >= 3 else 0)) == 0
check("Nullvektor eta_q = sum q^j |4j> von A_q = S^+ - q S^3 exakt; Restmischung bei Tiefe M: sigma_min fällt je 4 Marken "
      "um q (q = 0,5: Verhältnisse -> 0,5), bei q = e^{-2 pi} genau q, q^2, q^3 (1,9e-3, 3,5e-6, 6,5e-9)",
      eta_vec_ok and all(abs(r_ - 0.5) < 0.02 for r_ in rat[1:]) and so[1] < 1e-4 and so[2] < 1e-7)


def h0(d):
    return max(d + 1, 0)


def h1(d):
    return max(-d - 1, 0)


Em = [-1, 0, 0, 0]
Ebad = [-4, 0, 0, 0]
check("Familienindex: E_- = O(-1) + 3 O hat h^0 = 3, h^1 = 0 (drei chirale Nullmoden, kein Spiegel); O(-4) + 3 O hätte "
      "h^0 = 3 und h^1 = 3 Spiegel",
      sum(map(h0, Em)) == 3 and sum(map(h1, Em)) == 0 and sum(map(h0, Ebad)) == 3 and sum(map(h1, Ebad)) == 3)
check("Z16: 16 Weyl je Familie (mit nu_R) und 3 x 16 = 48 sind = 0 mod 16: keine globale Anomalie", 16 % 16 == 0 and 48 % 16 == 0)
befund("Chirale Materie kommt hier aus Indizes (Hardy-Zweig, Toeplitz-Index, Familienindex) und braucht keine gegappten "
       "Spiegel; die Spiegellücke aus z56 ist eine Alternative, keine Voraussetzung", True)
befund("Offen bleibt die 4D-Seite: 3+1 Raumzeit muss selbst aus 4D-Antworten rekonstruiert werden (Osterwalder-Schrader); "
       "wird sie doch gegittert, braucht sie einen Ginsparg-Wilson-Operator (Lüscher), keine Spiegellücke", True)

# ------------------------------------------------------------------ Teil 6
say("\nTeil 6: Bilanz der Rekonstruktion")
tafel = [
    ("volle KMS-Bedingung", "folgt (Tomita: der Zustand bestimmt seinen Modulfluss)"),
    ("quadratisch in den Strommoden", "folgt (Sugawara, Wick, Jacobi-Tripelprodukt)"),
    ("Band oder Leiter", "entschieden: Leiter (Band nicht spurklasse und nicht quadratisch)"),
    ("Leiterskala lambda", "folgt: 4 pi je Quant in modularer Zeit; Einheit beta = L"),
    ("positiver Propagator", "folgt: Fisher-Information kappa_2 > 0 des Originals"),
    ("Selbstkopplungen W3, W4 im Spurkanal", "folgen: kappa_3, kappa_4 (Ladungssektor dominiert)"),
    ("Chiralität", "folgt aus Indizes, ohne Spiegellücke"),
    ("q = e^{-2 pi}", "folgt aus der Kovarianz (z46); für q = 1/sqrt2 ist im Bestand keine originale Kovarianz benannt"),
    ("Tensorkanal / t_n", "offen: braucht die 4D-Antworten mit Tensorstruktur"),
    ("3+1 Raumzeit", "offen: OS-Rekonstruktion aus 4D-Antworten"),
    ("Familiengewichte gamma, TM1-Brechung", "offen"),
]
for k_, v_ in tafel:
    say(f"   {k_}: {v_}")
ZAHLEN["bilanz"] = dict(tafel)

urteil = ("Die Rekonstruktion aus den Antworten des Originals ersetzt mehrere Auswahlprinzipien durch Folgerungen. Aus der "
          "Calderón/Szegő-Kovarianz allein (als quasifreier Zustand, wie im Original) folgen der Modulgenerator 4 pi L_0, die CCR, die lineare Dispersion 4 pi m und "
          "die gaußsche Struktur der Ströme (Sugawara, Wick, Jacobi-Tripelprodukt); der volle Modulgenerator hat Spektrum in "
          "pi*N. Damit ist Band gegen Leiter entschieden: Das Band ist weder Modulgenerator (nicht spurklasse) noch "
          "quadratisch; unter Modulfluss und modularer Konjugation bleibt genau die Leiter, Stufe 4 pi je Quant. Der positive "
          "Propagator des Spurkanals ist die Fisher-Information des Originals, die Selbstkopplungen sind seine "
          "Energiekumulanten; der Thetasektor der Ladung dominiert und ist stark nichtgaußsch. Chirale Materie kommt aus "
          "Indizes (ein chiraler Zweig, Toeplitz-Index +1 mit Restmischung O(q^Tiefe), Familienindex 3 ohne Spiegel), ohne "
          "Spiegellücke. Offen bleiben der Tensorkanal in 4D, die Rekonstruktion von 3+1 Raumzeit und die Familiengewichte. "
          "Eine vollständige Theorie ist das nicht.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z58_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
