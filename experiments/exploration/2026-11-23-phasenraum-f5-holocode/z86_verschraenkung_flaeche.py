"""Z86: Verschränkung der Purifikation, gleiche Uhr und Flächennormierung (Strategie S5, Stand H6h41).

Leitlinie 'Rekonstruktion statt Auswahl'. Firewalls: GRAV.NONCIRCULAR.01 (kein Schluss 'modular => Bisognano-Wichmann'; der
Schnitt Original|Partner wird nicht als Rindlerhorizont gelesen), DIMENSION.UPLIFT.FIREWALL.01 (Aussagen über die Naht bleiben
1D, die Orte aus z78 bleiben Kinematik). Bestand: z58 (K = 4 pi L_0, Ladung x Oszillatoren, kappa_2 = 0,798), z68 (S = ln 2/4
je Quarter, <T> = 0 am KMS-Punkt, E8-Naht ln 12, keine Skala), z78 (Purifikation Omega, q = e^{-2 pi}, kappa, Orte als
Relativzustände), z79; Ledger GRAV.NONLINEAR.01 und GRAV.ENTROPY.EQUILIBRIUM.01 (8 pi = 1/c3 = 2 pi/eta, eta = 1/4, unter
Jacobson- und BW-Prämissen, nur gelesen).

Teil 1  Verschränkung von Omega zwischen Original und Tomita-Partner = thermische Entropie des Originals bei beta = 4 pi:
        je Fermionmode, je Stromoszillator m (q_m = e^{-2 pi m}) und im Ladungssektor (Gewichte e^{-4 pi E_alpha(Q)}), beide
        Quarter; geschlossene Formen über theta_3(i), eta(i), eta(2i), eta(4i); Rényi-Entropien und S-Dualität n <-> 1/n;
        <K>, Kumulanten kappa_n und Kapazität.
Teil 2  Vergleich mit TFPT-Daten ohne Fit (1/c3 = 8 pi, c3, 48 c3^4, phi0, eta_B = 1/4, t = 1, ln 2, E8-Naht ln 12):
        abgeleitete Relationen, PSLQ, arithmetische Klassen, Kontrolle der 1/4, Zahlenspielbasislinie.
Teil 3  Erster Hauptsatz je Sektor mit derselben beta (Ladung, Oszillatoren, Fermionmoden, native Translationen e^{iP(X)}),
        Eötvös-Parameter der Sektoren, Kontrolle mit getrennten Temperaturen, zweite Ordnung und Linearitätsfenster.
Teil 4  Flächennormierung: Kandidaten (Moden mit Abschneidung, Nahtlänge, Zentralladung, Orte, Kopien, c3) und die Typfrage
        (Thermofeldschnitt einer Typ-I-Algebra gegen Horizont).
Teil 5  Urteil für S5 gegen Erfolgs- und Abbruchkriterium.
"""
import json
from collections import Counter
from pathlib import Path

import mpmath as mp
import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spla
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


def ns(x, n=12):
    return mp.nstr(x, n)


PI = mp.pi
BETA = 4*PI                                                     # KMS-Punkt des Originals (Quarter-Einheiten von L_0)
ALPHAS = (mp.mpf(1)/4, mp.mpf(3)/4)
QR, MR, RR = 40, 80, 60                                         # Ladungen |Q| <= 40, Oszillatoren m < 80, Fermionmoden |r| <= 60
LN2 = mp.log(2)
TOL = mp.mpf(10)**-45
TH3 = PI**(mp.mpf(1)/4)/mp.gamma(mp.mpf(3)/4)                   # theta_3(i)
A4 = TH3**4                                                     # theta_3(i)^4 = pi/Gamma(3/4)^4
ETA_I = mp.gamma(mp.mpf(1)/4)/(2*PI**(mp.mpf(3)/4))
ETA_2I = mp.gamma(mp.mpf(1)/4)/(2**(mp.mpf(11)/8)*PI**(mp.mpf(3)/4))
ETA_4I = mp.gamma(mp.mpf(1)/4)*(mp.sqrt(2) - 1)**(mp.mpf(1)/4)/(2**(mp.mpf(29)/16)*PI**(mp.mpf(3)/4))
Z58 = json.loads((HERE / "z58_ergebnis.json").read_text())


def E_lad(a, Q):
    return mp.mpf(Q)**2/2 + (a - mp.mpf(1)/2)*Q


def ladung(a, beta=BETA, mu=0):
    """Ladungssektor eines Quarters mit Gewichten e^{-beta E(Q) + mu Q}: (ln Z, <E>, <Q>, S, p)."""
    ws = {Q: mp.exp(-beta*E_lad(a, Q) + mu*Q) for Q in range(-QR, QR + 1)}
    Z = mp.fsum(ws.values())
    p = {Q: w/Z for Q, w in ws.items()}
    Em = mp.fsum(p[Q]*E_lad(a, Q) for Q in p)
    Qm = mp.fsum(p[Q]*Q for Q in p)
    S = -mp.fsum(v*mp.log(v) for v in p.values())
    return mp.log(Z), Em, Qm, S, p


def lnZ_lad(beta, a):
    return mp.log(mp.fsum(mp.exp(-beta*E_lad(a, Q)) for Q in range(-QR, QR + 1)))


def lnZ_osz(beta):
    return -mp.fsum(mp.log(1 - mp.exp(-beta*m)) for m in range(1, MR))


def s_bose(nu):
    return (nu + 1)*mp.log(nu + 1) - (nu*mp.log(nu) if nu > 0 else 0)


def s_mode(m, beta=BETA):
    """Verschränkung einer Oszillatormode m mit ihrem Partner: Schmidtgewichte (1 - q_m^2) q_m^{2n}, q_m = e^{-beta m/2}."""
    x = mp.exp(-beta*m)
    return -mp.log(1 - x) + beta*m*x/(1 - x)


def h_bin(c):
    return -c*mp.log(c) - (1 - c)*mp.log(1 - c)


def c_fermi(x, beta=BETA):
    return 1/(1 + mp.exp(beta*abs(x)))                          # h(c(x)) = h(c(|x|)): Teilchen und Löcher gleich


def eta_y(y):
    return mp.exp(-PI*y/12)*mp.qp(mp.exp(-2*PI*y))


def th(k, y):
    return mp.jtheta(k, 0, mp.exp(-PI*y))


def lnZv_dirac(y):
    """ln(theta_3/eta)(iy): voller (Casimir-inklusiver) Zustandssummenlogarithmus des komplexen NS-Fermions."""
    return mp.log(th(3, y)) - mp.log(eta_y(y))


def lnZv_e8(y):
    return mp.log((th(2, y)**8 + th(3, y)**8 + th(4, y)**8)/2) - 8*mp.log(eta_y(y))


def entropie_voll(f, y=1):
    return f(y) - y*mp.diff(f, y)


# =================================================================== Teil 1
say("Teil 1: Verschränkung von Omega zwischen Original und Tomita-Partner (beta = 4 pi)")

# --- 1.1 Spektrum: zwei Quarter = ein komplexes NS-Fermion bei tau = i
lev = {a: Counter(int(4*abs(r + a)) for r in range(-20, 20)) for a in (sp.Rational(1, 4), sp.Rational(3, 4))}
ung = Counter(range(1, 80, 2))
ok11 = all(lev[a] == ung for a in lev) and lev[sp.Rational(1, 4)] + lev[sp.Rational(3, 4)] == Counter(
    {k: 2 for k in range(1, 80, 2)})
beta_ns = [mp.log(1/c_fermi(mp.mpf(k)/4) - 1)/(mp.mpf(k)/2) for k in (1, 3, 5, 7)]    # NS-Energie 2|x| = k/2
ok11 &= all(abs(b - 2*PI) < TOL for b in beta_ns)
check("Spektrum: die modularen Einteilchenenergien 4 pi |x| eines Quarters sind pi l mit l ungerade, je einmal (Spektrum eines "
      "Majorana-NS-Fermions bei tau = i); beide Quarter zusammen geben jedes pi l zweimal, also Teilchen und Löcher eines "
      "komplexen NS-Fermions mit 2 pi k, k in N - 1/2: das Original ist spektral ein komplexes NS-Fermion bei tau = i "
      "(beta = L, KMS-Normierung 2 pi bezüglich L_0^NS = 2 L_0)", ok11)

# --- 1.2 Schmidt-Zerlegung der Purifikation auf dem Fockraum (Konstruktion wie z78)
qo = float(mp.exp(-2*PI))
kap2 = (1 - qo)/(1 + qo)
nbar = qo**2/(1 - qo**2)


def fock4(n):
    a = sps.diags(np.sqrt(np.arange(1, n)), 1, format="csr")
    Id = sps.identity(n, format="csr")
    ops = []
    for k in range(4):
        mats = [Id]*4
        mats[k] = a
        op = mats[0]
        for Mx in mats[1:]:
            op = sps.kron(op, Mx, format="csr")
        ops.append(op)
    return ops


NF = 10
a1, a2, b1, b2 = fock4(NF)
om = np.zeros([NF]*4)
for i in range(NF):
    for j in range(NF):
        om[i, j, i, j] = (1 - qo**2)*qo**(i + j)
om = om.reshape(-1).astype(complex)
dvec = (a1 - qo*b1.T) @ om                                       # Araki-Woods-Vernichter (z78)


def s_reduziert(st):
    sv = np.linalg.svd(st.reshape(NF*NF, NF*NF), compute_uv=False)
    p = sv**2
    p = p[p > 1e-300]/np.sum(p)
    return float(-np.sum(p*np.log(p)))


def sB_float(nu):
    return (nu + 1)*np.log1p(nu) - (nu*np.log(nu) if nu > 0 else 0.0)


S_fock0 = s_reduziert(om)
keep = [i for i in range(NF**4) if all((i // NF**k) % NF < NF - 2 for k in range(4))]
check("Fockraum (vier Moden a_1, a_2, b_1, b_2, je 10 Stufen, q = e^{-2 pi}): Omega wird von d_1 = (a_1 - q b_1^+)/sqrt(1 - q^2) "
      f"vernichtet und hat die Schmidtentropie {S_fock0:.6e} = 2 s(nbar) der Stromoszillatoren m = 1 beider Quarter "
      "(Verschränkung Original|Partner = thermische Entropie des Originals)",
      np.max(np.abs(dvec[keep])) < 1e-12 and abs(S_fock0/(2*sB_float(nbar)) - 1) < 1e-9)

# --- 1.3 Fermionmoden gegen Bosonisierung
SF_q = {a: mp.fsum(h_bin(c_fermi(r + a)) for r in range(-RR, RR)) for a in ALPHAS}
S_lad = {a: ladung(a)[3] for a in ALPHAS}
S_osz_q = mp.fsum(s_mode(m) for m in range(1, MR))
S_fermi = mp.fsum(SF_q.values())
S_bos = mp.fsum(S_lad.values()) + 2*S_osz_q
S_TOT = S_bos
check("Verschränkung von Omega: Summe der Fermionmoden-Entropien h(c(x)) beider Quarter = Ladungssektoren + Stromoszillatoren "
      "(Bosonisierung auf Entropieebene) = ln 2/2 je Naht und ln 2/4 je Quarter (auf 1e-45; z68 hatte die Summe)",
      abs(S_fermi - S_bos) < TOL and abs(S_bos - LN2/2) < TOL and all(abs(SF_q[a] - LN2/4) < TOL for a in ALPHAS)
      and all(abs(S_lad[a] + S_osz_q - LN2/4) < TOL for a in ALPHAS))

# --- 1.4 geschlossene Formen
E2 = lambda y: 1 - 24*mp.nsum(lambda n: n*mp.exp(-2*PI*y*n)/(1 - mp.exp(-2*PI*y*n)), [1, mp.inf])
E4 = lambda y: 1 + 240*mp.nsum(lambda n: n**3*mp.exp(-2*PI*y*n)/(1 - mp.exp(-2*PI*y*n)), [1, mp.inf])
ok_gamma = (abs(th(3, 1) - TH3) < TOL and abs(eta_y(1) - ETA_I) < TOL and abs(eta_y(2) - ETA_2I) < TOL
            and abs(eta_y(4) - ETA_4I) < TOL)
ok_eis = (abs(E2(1) - 3/PI) < TOL and abs(E2(2) - (3/(2*PI) + 3*A4/8)) < TOL and abs(E4(1) - 3*A4**2/4) < TOL
          and abs(E4(2) - 11*E4(1)/16) < TOL)
S_lad_cf = mp.log(TH3) - 5*LN2/8 + mp.mpf(1)/4 + PI*A4/16
S_osz_cf = -mp.log(ETA_2I) - mp.mpf(1)/4 - PI*A4/16
K_lad_cf = mp.mpf(1)/4 - PI/8 + PI*A4/16
K_osz_cf = PI/6 - mp.mpf(1)/4 - PI*A4/16
K_lad = {a: BETA*ladung(a)[1] for a in ALPHAS}
K_osz_q = mp.fsum(BETA*m*mp.exp(-BETA*m)/(1 - mp.exp(-BETA*m)) for m in range(1, MR))
check("theta_3(i) = pi^{1/4}/Gamma(3/4), eta(i) = Gamma(1/4)/(2 pi^{3/4}), eta(2i) = 2^{-3/8} eta(i), eta(4i) = Gamma(1/4) "
      "(sqrt2 - 1)^{1/4}/(2^{29/16} pi^{3/4}); E_2(i) = 3/pi, E_2(2i) = 3/(2 pi) + (3/8) theta_3(i)^4, E_4(2i) = (11/16) E_4(i), "
      "E_4(i) = (3/4) theta_3(i)^8 (gegen q-Reihen, 1e-45)", ok_gamma and ok_eis)
check("geschlossene Formen je Quarter: S_Ladung = ln theta_3(i) - (5/8) ln 2 + 1/4 + (pi/16) theta_3(i)^4, S_Osz = -ln eta(2i) "
      "- 1/4 - (pi/16) theta_3(i)^4 (Summe ln 2/4); <K_Ladung> = 1/4 - pi/8 + (pi/16) theta_3(i)^4, <K_Osz> = pi/6 - 1/4 - "
      "(pi/16) theta_3(i)^4 (Summe pi/24); für beide Quarter gleich",
      all(abs(S_lad[a] - S_lad_cf) < TOL and abs(K_lad[a] - K_lad_cf) < TOL for a in ALPHAS)
      and abs(S_osz_q - S_osz_cf) < TOL and abs(K_osz_q - K_osz_cf) < TOL and abs(S_lad_cf + S_osz_cf - LN2/4) < TOL)
anteil_lad = 2*S_lad_cf/S_TOT
S_osz_beide = 2*S_osz_cf
ZAHLEN["entropie"] = {"Naht (beide Quarter)": "ln 2/2 = " + ns(S_TOT, 20), "je Quarter": "ln 2/4",
                      "Ladung je Quarter": ns(S_lad_cf, 15), "Oszillatoren je Quarter": ns(S_osz_cf, 15),
                      "Oszillatoren beide Quarter": ns(S_osz_beide, 15), "Anteil Ladung": ns(anteil_lad, 12),
                      "Anteil Oszillatoren": ns(1 - anteil_lad, 6)}
moden = {f"m={m}": {"q_m": ns(mp.exp(-2*PI*m), 6), "s_m je Quarter": ns(s_mode(m), 8)} for m in range(1, 5)}
fmoden = {f"k={k}/2": ns(2*h_bin(c_fermi(mp.mpf(k)/4)), 8) for k in (1, 3, 5, 7)}
ZAHLEN["je_mode"] = {"Oszillatoren": moden, "NS-Fermionmoden 2 pi k (Teilchen + Loch)": fmoden}
say(f"   S = ln 2/2 = {ns(S_TOT, 15)}: Ladung {ns(2*S_lad_cf, 12)}, Oszillatoren {ns(S_osz_beide, 6)} (Anteil "
    f"{ns(1 - anteil_lad, 4)}); Mode m = 1: s_1 = {ns(s_mode(1), 6)} je Quarter, m = 2: {ns(s_mode(2), 4)}")
check("die Verschränkung sitzt im Ladungssektor: Oszillatoranteil < 3e-4 (beide Quarter 9,46e-5), getragen fast ganz von der "
      "untersten Fermionmode (2 pi k = pi, Twistzustand Q = +-1); die Stromoszillatoren, die in z78/z79 die Raumzeitkinematik "
      "tragen, sind mit ihren Partnern fast unverschränkt", 1 - anteil_lad < 3e-4 and S_osz_beide < 1e-4
      and 2*h_bin(c_fermi(mp.mpf(1)/4))/S_TOT > 0.98)

# --- 1.5 Rényi-Entropien
def lnspur_q(a, n, beta=BETA):
    """ln Tr rho^n eines Quarters (Ladung x Oszillatoren)."""
    osz = mp.fsum(n*mp.log(1 - mp.exp(-beta*m)) - mp.log(1 - mp.exp(-n*beta*m)) for m in range(1, MR))
    return lnZ_lad(n*beta, a) - n*lnZ_lad(beta, a) + osz


def S_n_spek(n):
    if n == 1:
        return S_TOT
    return mp.fsum(lnspur_q(a, n) for a in ALPHAS)/(1 - n)


def S_n_q(n, a=ALPHAS[0]):
    return LN2/4 if n == 1 else lnspur_q(a, n)/(1 - n)


def S_n_lad(n, a=ALPHAS[0]):
    if n == 1:
        return S_lad[a]
    return (lnZ_lad(n*BETA, a) - n*lnZ_lad(BETA, a))/(1 - n)


def S_n_cf(n):
    return (lnZv_dirac(n) - n*LN2/2)/(1 - n)


S2_cf = 7*LN2/8 - mp.log(1 + mp.sqrt(2))/2
Sh_cf = mp.log(1 + mp.sqrt(2)) - LN2/4
Sinf_cf = LN2/2 - PI/12
S_inf = mp.fsum(lnZ_lad(BETA, a) for a in ALPHAS) + 2*lnZ_osz(BETA)    # -ln lambda_max, Grundzustand K = 0 eindeutig
ok_r = (abs(S_n_spek(2) - S2_cf) < TOL and abs(S_n_spek(mp.mpf(1)/2) - Sh_cf) < TOL and abs(S_inf - Sinf_cf) < TOL
        and all(abs(S_n_spek(n) - S_n_cf(n)) < TOL for n in (3, mp.mpf(1)/3, 5)))
check("Rényi aus dem Spektrum = (ln(theta_3/eta)(in) - n ln sqrt2)/(1 - n): S_2 = (7/8) ln 2 - (1/2) ln(1 + sqrt2), "
      "S_{1/2} = ln(1 + sqrt2) - (1/4) ln 2, S_inf = ln 2/2 - pi/12 (n = 3, 1/3, 5 numerisch), auf 1e-45", ok_r)
dual = {}
for n in (2, 3, mp.mpf(3)/2, 5):
    dual[str(n)] = (S_n_spek(1/mp.mpf(n)) + n*S_n_spek(n) - (n + 1)*LN2/2,
                    S_n_q(1/mp.mpf(n)) + n*S_n_q(n) - (n + 1)*LN2/4,
                    S_n_lad(1/mp.mpf(n)) + n*S_n_lad(n) - (n + 1)*S_lad[ALPHAS[0]])
ZAHLEN["renyi"] = {"S_2": "(7/8) ln 2 - (1/2) ln(1+sqrt2) = " + ns(S2_cf), "S_1/2": "ln(1+sqrt2) - (1/4) ln 2 = " + ns(Sh_cf),
                   "S_inf": "ln 2/2 - pi/12 = " + ns(Sinf_cf), "S_3": ns(S_n_spek(3)),
                   "Dualitätsrest Ladung allein": {k: ns(v[2], 6) for k, v in dual.items()}}
check("S-Dualität der Rényi-Entropien: S_{1/n} + n S_n = (n + 1) S_1 für das Original (S_1 = ln 2/2) und für jedes Quarter "
      "(ln 2/4), n = 2, 3, 3/2, 5 (auf 1e-45); für den Ladungssektor allein gilt sie nicht (Rest 6e-4 bis 0,11): sie ist "
      "eine "
      "Eigenschaft des ganzen Quarters (Weber f(i) = 2^{1/4}, S-invariant), nicht der Bosonisierungsfaktoren",
      all(abs(v[0]) < TOL and abs(v[1]) < TOL and abs(v[2]) > 1e-4 for v in dual.values()))
S2_lad = -(lnZ_lad(2*BETA, ALPHAS[0]) - 2*lnZ_lad(BETA, ALPHAS[0]))
S2_osz = -(lnZ_osz(2*BETA) - 2*lnZ_osz(BETA))
check("Rényi-2 der Sektoren je Quarter: S_2(Ladung) = ln theta_3(i) = (1/4) ln pi - ln Gamma(3/4) (Reinheit 1/theta_3(i) = "
      "0,920), S_2(Osz) = ln eta(4i) - 2 ln eta(2i)", abs(S2_lad - mp.log(TH3)) < TOL
      and abs(S2_osz - (mp.log(ETA_4I) - 2*mp.log(ETA_2I))) < TOL)

# --- 1.6 Modulhamiltonoperator: Erwartung, Kumulanten, Kapazität
e_, x_, y_, z_ = sp.symbols("e x y z")
E4s = (x_**2 + y_**2 + z_**2)/2


def Dram(f):
    """q d/dq im Ramanujan-System: x, y, z = theta_2^4, theta_3^4, theta_4^4, e = E_2."""
    return sp.expand(sp.diff(f, e_)*(e_**2 - E4s)/12 + sp.diff(f, x_)*x_*(e_ + y_ + z_)/6
                     + sp.diff(f, y_)*y_*(e_ + x_ - z_)/6 + sp.diff(f, z_)*z_*(e_ - x_ - y_)/6)


aa = sp.symbols("a", positive=True)
ders = [(x_ - z_)/24]                                           # D ln(theta_3/eta) = (theta_2^4 - theta_4^4)/24
for _ in range(4):
    ders.append(Dram(ders[-1]))
sub_i = {e_: 3/sp.pi, x_: aa/2, y_: aa, z_: aa/2}
kap_sym = {n: sp.factor(sp.simplify((2*sp.pi)**n*ders[n - 1].subs(sub_i))) for n in range(1, 6)}
r_n = {n: sp.nsimplify(sp.simplify(kap_sym[n]/(sp.pi**2*aa**2))) for n in range(2, 6)}


def W(J):
    return mp.fsum(lnZ_lad(BETA*(1 - J), a) for a in ALPHAS) + 2*lnZ_osz(BETA*(1 - J))


kap = {n: mp.diff(W, 0, n) for n in (1, 2, 3, 4)}
casimir = 2*BETA*(-mp.mpf(1)/96)                               # 4 pi (h - c/24) je Quarter, z68
ok_k = abs(kap[1] - PI/12) < TOL and abs(kap[1] + casimir) < TOL and kap_sym[1] == 0
ok_k &= all(abs(kap[n] - mp.mpf(str(sp.N(r_n[n], 70)))*PI**2*A4**2) < mp.mpf(10)**-35 for n in (2, 3, 4))
ok_k &= all(abs(kap[n] - mp.mpf(Z58["zahlen"]["energiekumulanten_K"][f"kappa{n}"])) < 1e-9 for n in (2, 3, 4))
ZAHLEN["modulhamilton"] = {"<K>": "pi/12 = " + ns(kap[1]), "<K> + Casimir": ns(kap[1] + casimir, 3),
                           "<K_Ladung> beide": ns(2*K_lad_cf), "<K_Osz> beide": ns(2*K_osz_q, 8),
                           "kappa_n / (pi^2 theta_3(i)^8)": {n: str(r_n[n]) for n in r_n},
                           "kappa_2": "Gamma(1/4)^8/(384 pi^4) = " + ns(kap[2])}
say(f"   kappa_n = r_n pi^2 theta_3(i)^8 mit r_2..r_5 = {[str(r_n[n]) for n in range(2, 6)]} (Ramanujan-System bei tau = i)")
check("Modulhamiltonoperator: <K> = pi/12 hebt die Casimirenergie 2 x 4 pi (-1/96) auf (<K'> = 0, z68); Kumulanten in "
      "geschlossener Form kappa_2 = pi^2 theta_3(i)^8/24 = Gamma(1/4)^8/(384 pi^4), kappa_3 = 3 kappa_2, kappa_4 = 9 kappa_2 "
      "(E_2-Terme heben sich auf; gegen Ableitungen und die z58-Zahlen 0,7982/2,3946/7,1839 geprüft)",
      ok_k and abs(kap[2] - mp.gamma(mp.mpf(1)/4)**8/(384*PI**4)) < TOL and r_n[3] == 3*r_n[2] and r_n[4] == 9*r_n[2])
WE8 = lambda J: lnZv_e8(1 - J)
kE8 = {n: mp.diff(WE8, 0, n) for n in (1, 2, 3, 4)}
check("kappa_3 = 3 kappa_2 ist eine Folge der S-Dualität (ln Z_voll(y) = ln Z_voll(1/y) ist gerade in ln y, also verschwinden "
      "erste und dritte Ableitung in ln y): ebenso für die E8-Naht (kappa_2 = pi^2 theta_3(i)^8/2, kappa_3/kappa_2 = 3), "
      f"kappa_4/kappa_2 dagegen theorieabhängig (Original 9, E8 {ns(kE8[4]/kE8[2], 6)})",
      abs(kE8[1]) < mp.mpf(10)**-40 and abs(kE8[3]/kE8[2] - 3) < mp.mpf(10)**-35
      and abs(kE8[2] - PI**2*A4**2/2) < mp.mpf(10)**-35 and abs(kE8[4]/kE8[2] - 9) > 1)
C_E = kap[2]
ZAHLEN["kapazitaet"] = {"C_E = Var K": ns(C_E), "C_E/S": "pi^2 theta_3(i)^8/(12 ln 2) = " + ns(C_E/S_TOT)}

kl = mp.diff(lambda J: mp.fsum(lnZ_lad(BETA*(1 - J), a) for a in ALPHAS), 0, 2)
ko = mp.diff(lambda J: 2*lnZ_osz(BETA*(1 - J)), 0, 2)
kl_cf = 2*(mp.mpf(1)/8 + PI*A4/16)
ko_cf = 2*(PI**2*A4**2/48 - PI*A4/16 - mp.mpf(1)/8)
check("Fisher-Information je Sektor in geschlossener Form: Ladung 2(1/8 + (pi/16) theta_3(i)^4) = 0,79711, Oszillatoren "
      "2((pi^2/48) theta_3(i)^8 - (pi/16) theta_3(i)^4 - 1/8) = 1,1014e-3 (z58: 0,7971099 und 0,0011014; Verhältnis 724)",
      abs(kl - kl_cf) < mp.mpf(10)**-35 and abs(ko - ko_cf) < mp.mpf(10)**-35
      and abs(kl - mp.mpf(Z58["zahlen"]["ladung_vs_oszillator_kappa2"][0])) < 1e-9
      and abs(ko - mp.mpf(Z58["zahlen"]["ladung_vs_oszillator_kappa2"][1])) < 1e-12)
ZAHLEN["fisher_sektoren"] = {"Ladung": ns(kl_cf), "Oszillatoren": ns(ko_cf), "Verhältnis": ns(kl_cf/ko_cf, 6)}
befund("Am nativen Punkt ist die Verschränkung von Omega klein und innen: ein halbes Bit je Naht, zu 99,97 % im Ladungssektor "
       "(Uhr, Marken), 9,5e-5 in den Stromoszillatoren (Raumzeitkinematik). Die Kapazität C_E = Var K ist 2,30 S; im "
       "extensiven Grenzfall wäre sie S (Teil 4)", True)

# =================================================================== Teil 2
say("\nTeil 2: Vergleich mit TFPT-Daten ohne Fit")

# --- 2.1 abgeleitete Relationen
prod_q = mp.fprod(1 + mp.exp(-PI*l) for l in range(1, 400, 2))
Zv_q = prod_q*mp.exp(PI/24)
t2, t3, t4, et = th(2, 1), th(3, 1), th(4, 1), eta_y(1)
S_q_v = entropie_voll(lambda y: lnZv_dirac(y)/2)
S_d_v = entropie_voll(lnZv_dirac)
S_e8_v = entropie_voll(lnZv_e8)
ok21 = (abs(Zv_q - 2**(mp.mpf(1)/4)) < TOL and abs(t3/et - mp.sqrt(2)) < TOL and abs(t2*t3*t4 - 2*et**3) < TOL
        and abs(t2 - t4) < TOL and abs(mp.exp(lnZv_e8(1)) - 12) < mp.mpf(10)**-40)
ok21 &= all(abs(mp.diff(f, 1)) < mp.mpf(10)**-40 for f in (lnZv_dirac, lnZv_e8))
ok21 &= abs(S_q_v - LN2/4) < TOL and abs(S_d_v - LN2/2) < TOL and abs(S_e8_v - mp.log(12)) < mp.mpf(10)**-40
check("R1: Für jede S-invariante volle Zustandssumme ist d/dy ln Z_voll = 0 bei tau = i, also S = ln Z_voll(i): Quarter "
      "(Weber f(i) = 2^{1/4}) ln 2/4, Original (theta_3/eta = sqrt2) ln 2/2, E8-Naht (E_4/eta^8 = 12) ln 12", ok21)
check("R2: Herkunft der ln 2: theta_2 theta_3 theta_4 = 2 eta^3 (die 2 ist die Entartung des Ramond-Grundzustands, eine "
      "komplexe Nullmode) und theta_2(i) = theta_4(i) = 2^{-1/4} theta_3(i) (lemniskatischer Modul k^2 = 1/2) geben "
      "(theta_3/eta)(i)^3 = 2^{3/2}", abs((t3/et)**3 - 2**(mp.mpf(3)/2)) < TOL and abs(t2/t3 - 2**(-mp.mpf(1)/4)) < TOL)
gso = (1 + (t2/t3)**8 + (t4/t3)**8)/2
check("R3: als Zustandssumme ist die E8-Naht 16 Quarter (je ein Majorana-NS-Spektrum, Z_voll = 2^{1/4}) mit GSO-Projektion: "
      "Z_E8(i) = 16 x 3/4 mit dem GSO-Gewicht (1 + theta_2^8/theta_3^8 + theta_4^8/theta_3^8)/2 = 3/4 bei tau = i; also "
      "ln 12 = 16 (ln 2/4) + ln(3/4) (die Gleichsetzung der 16 Majoranas des Trägers mit 16 Quartern ist Lesart)",
      abs(gso - mp.mpf(3)/4) < TOL and abs(16*LN2/4 + mp.log(gso) - mp.log(12)) < TOL)
S2_e8 = -(lnZv_e8(2) - 2*lnZv_e8(1))
check("E8-Naht: Z_voll(2i) = 66 (E_4(2i) = (11/16) E_4(i), eta(2i)^8 = eta(i)^8/8), also S_2 = ln(24/11); S_inf = ln 12 - "
      "2 pi/3; Rényi-Dualität S_{1/2} + 2 S_2 = 3 ln 12", abs(mp.exp(lnZv_e8(2)) - 66) < mp.mpf(10)**-40
      and abs(S2_e8 - mp.log(mp.mpf(24)/11)) < TOL
      and abs(2*(lnZv_e8(mp.mpf(1)/2) - lnZv_e8(1)/2) + 2*S2_e8 - 3*mp.log(12)) < mp.mpf(10)**-40)
ZAHLEN["e8_naht"] = {"S": "ln 12", "S_2": "ln(24/11) = " + ns(S2_e8), "S_inf": "ln 12 - 2 pi/3 = " + ns(mp.log(12) - 2*PI/3),
                     "kappa_2": "pi^2 theta_3(i)^8/2 = " + ns(kE8[2]), "S/c": ns(mp.log(12)/8)}

# --- 2.2 PSLQ gegen die TFPT-Zahlen
c3 = 1/(8*PI)
ZIELE = {"1/c3 = 8 pi": 1/c3, "c3": c3, "48 c3^4": 48*c3**4, "phi0": 1/(6*PI) + 48*c3**4, "ln 2": LN2,
         "ln 12": mp.log(12)}
RATIONAL = {"eta_B = 1/4": mp.mpf(1)/4, "t = 1": mp.mpf(1)}
NATIV = {"S": S_TOT, "S je Quarter": LN2/4, "S_2": S2_cf, "S_1/2": Sh_cf, "S_inf": Sinf_cf, "<K>": kap[1], "kappa_2": kap[2],
         "C_E/S": C_E/S_TOT, "S_Ladung je Quarter": S_lad_cf, "S_Osz beide": S_osz_beide, "S_2 Ladung": S2_lad,
         "S_E8": S_e8_v, "S_2 E8": S2_e8, "S_inf E8": mp.log(12) - 2*PI/3}


def relation(vec):
    r = mp.pslq(vec, tol=mp.mpf(10)**-40, maxcoeff=10**4, maxsteps=10**5)
    return None if r is None or r[0] == 0 else r


gefunden = {}
for kn, xv in NATIV.items():
    for kt, tv in ZIELE.items():
        r = relation([xv, tv, 1])
        if r is not None:
            gefunden[(kn, kt)] = r
    r = relation([xv, 1])
    if r is not None:
        gefunden[(kn, "rational")] = r
erwartet = {("S", "ln 2"), ("S je Quarter", "ln 2"), ("<K>", "1/c3 = 8 pi"), ("S_E8", "ln 12")}
ZAHLEN["pslq"] = {f"{k[0]} ~ {k[1]}": str(v) for k, v in gefunden.items()}
say("   PSLQ (Koeffizienten <= 10^4, 60 Stellen): " + "; ".join(f"{k[0]} mit {k[1]}: {v}" for k, v in gefunden.items()))
check("PSLQ findet unter 14 nativen Zahlen x 8 TFPT-Zahlen (eta_B = 1/4 und t = 1 als Rationalitätstest) genau vier "
      "ganzzahlige Relationen: S = ln 2/2, S_Quarter = ln 2/4, "
      "S_E8 = ln 12 (abgeleitet, R1-R3) und <K> = pi/12 = (1/96)(1/c3); keine Relation mit c3, 48 c3^4, phi0, eta_B = 1/4, "
      "t = 1", set(gefunden) == erwartet)
befund("Die einzige Relation mit 1/c3 ist <K> = 1/(96 c3): pi-Kommensurabilität (Casimirenergie 2 pi c/24 gegen 8 pi) ohne "
       "Rechnung, die die 1/96 erklärt; sie zählt nicht. Von 1/c3 = 2 pi/eta (Ledger) ist der Faktor 2 pi nativ als "
       "KMS-Normierung bezüglich L_0^NS (tau = i, 1D-typisiert, keine BW-Aussage); eta = 1/4 hat in der Verschränkung von "
       "Omega kein Gegenstück", True)
befund("Arithmetische Klassen: S, S_2, S_{1/2}, S_E8, S_2(E8) sind Logarithmen algebraischer Zahlen; für algebraisches "
       "alpha > 0, alpha != 1 ist ln(alpha)/pi irrational (Gelfond-Schneider, e^pi transzendent), also ist keine davon ein "
       "rationales Vielfaches von 1/c3 oder c3. kappa_2 = Gamma(1/4)^8/(384 pi^4) ist über Q(pi) transzendent (Nesterenko "
       "1996: pi und Gamma(1/4) algebraisch unabhängig), steht also in keiner algebraischen Relation zu c3, 48 c3^4 oder phi0", True)

# --- 2.3 Kontrolle der 1/4 in S_Quarter = ln 2/4
S_twist = {}
for al in ("1/8", "1/4", "3/8", "1/2"):
    a_ = mp.mpf(sp.Rational(al))
    S_twist[al] = mp.fsum(h_bin(c_fermi(r + a_)) for r in range(-RR, RR))
rel_twist = {al: relation([S_twist[al], LN2, 1]) for al in S_twist}
ZAHLEN["twistkontrolle_S/ln2"] = {al: ns(v/LN2, 10) for al, v in S_twist.items()}
check("Kontrolle der 1/4: bei beta = 4 pi hat nur der native Twist alpha = 1/4 eine rationale Entropie in Einheiten ln 2 "
      "(S/ln 2 = 0,6628; 1/4; 0,0784; 0,0392 für alpha = 1/8, 1/4, 3/8, 1/2; PSLQ nur bei 1/4); kein Gesetz S = alpha ln 2: die "
      "1/4 kommt vom selbstdualen Punkt (Weber f(i) = 2^{1/4}), nicht aus mu_4. 'eta_B = 1/4 je Bit' wäre Zahlenspiel",
      rel_twist["1/4"] is not None and all(rel_twist[al] is None for al in ("1/8", "3/8", "1/2"))
      and all(abs(S_twist[al] - mp.mpf(sp.Rational(al))*LN2) > 0.01 for al in ("1/8", "3/8", "1/2")))
befund("t = 1 (z57) ist ein Formkoeffizient der 4D-Plebański-Kantenwirkung; die Verschränkung von Omega trägt keine "
       "Flächenvariable, also gibt es keinen Vergleich. Die native '1' am Original ist der selbstduale Punkt beta/L = 1, "
       "ein anderes Objekt", True)

# --- 2.4 Zahlenspielbasislinie
basis = [float(v) for v in ZIELE.values()] + [1.0]
rats = sorted({p/q for p in range(1, 13) for q in range(1, 13)})
fam = np.array([r*b**s for b in basis for s in (1, -1) for r in rats])


def bester(xv):
    d = np.abs(xv/fam - 1)
    return float(d.min()), int(d.argmin())


rng = np.random.default_rng(86)
zufall = 10**rng.uniform(-5, 2, 20000)
d_zuf = np.array([bester(v)[0] for v in zufall])
nicht_abgeleitet = ["S_2", "S_1/2", "S_inf", "kappa_2", "C_E/S", "S_Ladung je Quarter", "S_Osz beide", "S_2 Ladung", "S_2 E8",
                    "S_inf E8"]
pwerte = {}
for kn in nicht_abgeleitet:
    d, _ = bester(float(NATIV[kn]))
    pwerte[kn] = (d, float(np.mean(d_zuf <= d)))
ZAHLEN["basislinie"] = {k: {"bester relativer Abstand": f"{v[0]:.2e}", "Anteil Zufallszahlen ebenso gut": f"{v[1]:.3f}"}
                        for k, v in pwerte.items()}
say("   bester Fasttreffer p/q x Ziel^{+-1} (p, q <= 12) je native Zahl, Anteil der Zufallszahlen mit ebenso gutem Treffer: "
    + "; ".join(f"{k} {v[0]:.1e} ({v[1]:.2f})" for k, v in pwerte.items()))
say("   (S_Ladung liegt 2,7e-4 unter ln 2/4, weil S_Osz klein ist: strukturell erklärt, kein Zahlenspiel)")
check(f"Zahlenspielbasislinie: Die Familie p/q x Ziel^{{+-1}} ({len(fam)} Werte) trifft Zufallszahlen typisch auf "
      f"{np.median(d_zuf):.1e}; kein Fasttreffer der zehn nicht abgeleiteten nativen Zahlen ist besser als 1 % der "
      "Zufallszahlen: Zahlennähe ist hier kein Befund", min(v[1] for v in pwerte.values()) > 0.01)

# =================================================================== Teil 3
say("\nTeil 3: Erster Hauptsatz je Sektor und dieselbe Uhr")
# --- 3.1 beta je Sektor aus den KMS-Verhältnissen
betas = {}
for a in ALPHAS:
    p = ladung(a)[4]
    betas[f"Ladung alpha={ns(a, 2)}"] = [mp.log(p[Q]/p[Q + 1])/(E_lad(a, Q + 1) - E_lad(a, Q)) for Q in range(-3, 4)]
    betas[f"Fermionmoden alpha={ns(a, 2)}"] = [mp.log((1 - 1/(1 + mp.exp(BETA*(r + a))))*(1 + mp.exp(BETA*(r + a))))/(r + a)
                                               for r in range(-3, 4)]
for m in range(1, 7):
    nm = 1/(mp.exp(BETA*m) - 1)
    betas[f"Oszillator m={m}"] = [mp.log((1 + nm)/nm)/m]
eotvos = max(abs(b/BETA - 1) for v in betas.values() for b in v)
beta_kontroll = {"Ladung": BETA, "Oszillatoren": 2*PI}
eotvos_k = abs(beta_kontroll["Ladung"] - beta_kontroll["Oszillatoren"])/((beta_kontroll["Ladung"]
                                                                          + beta_kontroll["Oszillatoren"])/2)
check("dieselbe Uhr: beta = 4 pi aus den KMS-Verhältnissen aller Sektoren (Ladung beider Quarter, Q = -3..3; Fermionmoden; "
      f"Oszillatoren m = 1..6); Eötvös-Parameter max|beta_s/4 pi - 1| = {ns(eotvos, 3)} (Kontrolle mit Oszillatoren bei "
      "2 pi: 2/3)", eotvos < TOL and abs(eotvos_k - mp.mpf(2)/3) < TOL)

# --- 3.2 erster Hauptsatz delta S = 4 pi delta <L_0>
def S_E_lad(mu):
    _, Em, _, S, _ = ladung(ALPHAS[0], mu=mu)
    return S, Em


def S_E_osz(lam, m=1, beta=BETA):
    x = mp.exp(-(beta*m - lam))
    n = x/(1 - x)
    return s_bose(n), m*n


def S_E_fermi(eps, xm=mp.mpf(1)/4):
    c = c_fermi(xm) + eps
    return h_bin(c), xm*c


erster = {}
for name, fn in (("Ladung (Fugazität e^{mu Q})", S_E_lad), ("Oszillator m=1", S_E_osz),
                 ("Oszillator m=2", lambda l: S_E_osz(l, 2)), ("Fermionmode x=1/4", S_E_fermi)):
    dS = mp.diff(lambda t: fn(t)[0], 0)
    dE = mp.diff(lambda t: fn(t)[1], 0)
    eps = mp.mpf("1e-3")
    rest = (fn(eps)[0] - fn(0)[0]) - BETA*(fn(eps)[1] - fn(0)[1])
    erster[name] = (dS/dE, rest)
ZAHLEN["erster_hauptsatz"] = {k: {"dS/d<L_0>": ns(v[0], 20), "Rest bei 1e-3 (= -S_rel)": ns(v[1], 4)}
                              for k, v in erster.items()}
dSk = mp.diff(lambda t: S_E_osz(t, 1, 2*PI)[0], 0)/mp.diff(lambda t: S_E_osz(t, 1, 2*PI)[1], 0)
check("erster Hauptsatz delta S = 4 pi delta<L_0> je Sektor (Ladung über die Fugazität, Oszillatoren m = 1, 2, einzelne "
      "Fermionmode): dS/d<L_0> = 4 pi auf 1e-40 mit derselben beta; der Rest zweiter Ordnung ist -S_rel < 0; Kontrolle: mit "
      "Oszillatoren bei 2 pi gäbe derselbe Test 2 pi", all(abs(v[0] - BETA) < mp.mpf(10)**-40 and v[1] < 0
                                                          for v in erster.values()) and abs(dSk - 2*PI) < mp.mpf(10)**-40)

# --- 3.3 native Translationen e^{iP(X)} (Oszillatorsektor, z78)
PAULI = [np.eye(2), np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.diag([1., -1.])]


def S_X_gauss(lams, prec=False):
    """Verschränkung von |X> = e^{iP(X)} Omega: <a^+ a>_X = nbar + (kappa^2/4) X^2 (keine anomalen Terme), Eigenwerte von X."""
    if prec:
        q = mp.exp(-2*PI)
        k2, nb = (1 - q)/(1 + q), q**2/(1 - q**2)
        return mp.fsum(s_bose(nb + k2/4*l**2) for l in lams) - 2*s_bose(nb), k2/4*mp.fsum(l**2 for l in lams)
    return sum(sB_float(nbar + kap2/4*l**2) for l in lams)


kapf = np.sqrt(kap2)
ct = [kapf*(a1 + b1.T), kapf*(a2 + b2.T)]


def zustand_X(X):
    P = 0.5*sum(ct[i].conj().T @ ct[j]*X[i, j] for i in range(2) for j in range(2))
    return spla.expm_multiply(1j*P, om)


ok_g = True
for _ in range(3):
    X = sum(c*s for c, s in zip(rng.normal(size=4)*0.3, PAULI))
    ok_g &= abs(s_reduziert(zustand_X(X))/S_X_gauss(np.linalg.eigvalsh(X)) - 1) < 1e-7
lam_hat = [mp.mpf(v) for v in np.linalg.eigvalsh(0.7*PAULI[0] + 0.4*PAULI[1] - 0.5*PAULI[3])]
rest_X = {}
for eps in ("1e-4", "5e-5", "1e-3", "0.025"):
    dS, dN = S_X_gauss([mp.mpf(eps)*l for l in lam_hat], prec=True)
    rest_X[eps] = (dS - BETA*dN)/(BETA*dN)
ZAHLEN["translation_erster_hauptsatz_rest"] = {k: ns(v, 4) for k, v in rest_X.items()}
check("native Translationen |X> = e^{iP(X)} Omega (Fockraum = Gaußformel <a^+ a>_X = nbar + (kappa^2/4) X^2, drei zufällige X, "
      "1e-7): delta S = 4 pi delta<N> mit derselben beta, Rest O(X^2) relativ (Verhältnis 4 bei Halbierung; -5e-5 bei |X| = "
      "1e-4); das Linearitätsfenster ist |X| << q = e^{-2 pi} (Rest -0,27 bei |X| = 0,025)",
      ok_g and abs(rest_X["1e-4"]/rest_X["5e-5"] - 4) < 0.01 and abs(rest_X["1e-4"]) < 1e-4
      and abs(rest_X["0.025"]) > 0.2)

# --- 3.4 zweite Ordnung und Fenster
p_lad1 = ladung(ALPHAS[0])[4][1]
n_osz1 = 1/(mp.exp(BETA) - 1)
check("zweite Ordnung ist nicht universell: Fisher-Information Ladung/Oszillatoren = 724, Besetzung des ersten angeregten "
      f"Niveaus Ladung {ns(p_lad1, 3)} gegen Oszillator m = 1 {ns(n_osz1, 3)} (Faktor > 1e4): gleiche lineare Kopplung, "
      "verschiedene Steifigkeit und verschiedene Linearitätsfenster", kl_cf/ko_cf > 700 and p_lad1/n_osz1 > 1e4)
befund("Gleiche Uhr heißt: dS/d<L_0^s> = beta = 4 pi für jeden Sektor s exakt (Eötvös-Parameter 0). In jeder Gleichgewichts"
       "bedingung der Form delta S_Original + delta S_X = 0 tritt die Energie aller Sektoren mit demselben Gewicht auf: "
       "universelle Kopplung erster Ordnung an die modulare Uhr (nullter Hauptsatz). Nicht daraus folgt: (a) ein Flächenterm "
       "S_X = eta A, (b) BW (K = 4 pi L_0 ist die KMS-Drehung der Naht, kein Boost), (c) der Wert 1/c3 oder G, (d) eine "
       "Tensorstruktur, (e) Universalität zweiter Ordnung", True)

# =================================================================== Teil 4
say("\nTeil 4: Flächennormierung -- gibt es eine native Fläche mit festem S/Fläche?")
# --- 4.1 Moden mit Abschneidung
S_M = {M: 2*S_lad_cf + 2*mp.fsum(s_mode(m) for m in range(1, M + 1)) for M in (1, 2, 3)}
N_eps = {}
for e in (10, 20, 40):
    N_eps[e] = sum(1 for k in range(1, 400, 2) if 2*h_bin(c_fermi(mp.mpf(k)/4)) > mp.mpf(10)**-e)
ZAHLEN["moden_abschneidung"] = {"S - S_M (M=1,2,3)": [ns(S_TOT - v, 3) for v in S_M.values()],
                                "NS-Moden mit Entropie > 1e-10, 1e-20, 1e-40": list(N_eps.values())}
check("Moden: die Verschränkung ist UV-endlich (S - S_M = 6e-10, 2e-15, ... für M Oszillatormoden), S/M -> 0; die Zahl der "
      f"Moden mit Entropie > eps wächst wie ln(1/eps)/(2 pi) ({list(N_eps.values())} für eps = 1e-10, 1e-20, 1e-40): "
      "S/Modenzahl hängt an der Abschneidung, nicht fest", S_TOT - S_M[1] < 1e-9 and S_TOT - S_M[2] < 1e-14
      and N_eps[10] < N_eps[20] < N_eps[40])

# --- 4.2 Nahtlänge: Familie tau = i y, y = beta/L
fam_y = {}
for yv in ("0.05", "0.2", "1", "2"):
    y = mp.mpf(yv)
    Sy = entropie_voll(lnZv_dirac, y)
    Cy = y**2*mp.diff(lnZv_dirac, y, 2)
    fam_y[yv] = (y*Sy, Cy/Sy)
ZAHLEN["nahtlaenge_y=beta/L"] = {k: {"S/(L/beta)": ns(v[0], 10), "C_E/S": ns(v[1], 8)} for k, v in fam_y.items()}
check("Nahtlänge: S/(L/beta) = pi/6 (Cardy, c = 1) und C_E/S = 1 nur im extensiven Grenzfall (y = beta/L = 0,05 auf 1e-15); "
      "am nativen Punkt y = 1 ist S/(L/beta) = ln 2/2 = 0,662 x Cardy und C_E/S = 2,30, bei y = 2 0,054: kein festes "
      "Verhältnis, der native Punkt ist nicht extensiv; nativ ist nur y = 1 (tau = i), also ein Punkt und keine Familie von "
      "Flächen, über die eine Dichte normiert werden könnte", abs(fam_y["0.05"][0] - PI/6) < mp.mpf(10)**-15
      and abs(fam_y["0.05"][1] - 1) < mp.mpf(10)**-15 and abs(fam_y["1"][0] - LN2/2) < TOL
      and abs(fam_y["1"][1] - 2.3031510400738) < 1e-12 and fam_y["0.05"][0]/fam_y["2"][0] > 9)

# --- 4.3 Zentralladung / Majoranazahl bei tau = i
theorien = {"Majorana-NS (ein Quarter), c = 1/2": (lambda y: lnZv_dirac(y)/2, mp.mpf(1)/2),
            "Dirac-NS (Original), c = 1": (lnZv_dirac, mp.mpf(1)),
            "8 Dirac-NS (D8, 16 Quarter), c = 8": (lambda y: 8*lnZv_dirac(y), mp.mpf(8)),
            "E8-Naht, c = 8": (lnZv_e8, mp.mpf(8))}
proc = {}
for k_, (f, cc) in theorien.items():
    proc[k_] = (entropie_voll(f)/cc, mp.diff(f, 1, 2)/cc)
ZAHLEN["je_zentralladung_tau_i"] = {k: {"S/c": ns(v[0], 10), "kappa_2/c": ns(v[1], 10)} for k, v in proc.items()}
vals = list(proc.values())
check("Zentralladung als Fläche: S/c = ln 2/2 für alle freien NS-Fermionen (Additivität der Kopien), aber E8 hat bei gleichem "
      "c = 8 ln 12/8 = 0,311 statt ln 16/8 = 0,347 (GSO nimmt ln(4/3)); kappa_2/c ebenso (0,798 gegen 1,197): bei tau = i "
      "legt c weder S noch kappa_2 fest (das gilt nur extensiv)",
      all(abs(v[0] - LN2/2) < mp.mpf(10)**-40 for v in vals[:3]) and abs(vals[3][0] - mp.log(12)/8) < mp.mpf(10)**-40
      and abs(vals[3][1]/vals[2][1] - mp.mpf(3)/2) < mp.mpf(10)**-35)

# --- 4.4 Orte: Verschränkung der Relativzustände
const_gross = 2 + 4*np.log(kapf/2)
gross = {}
for name, Xg in (("zeitartig r I", PAULI[0]), ("raumartig r sigma_3", PAULI[3]), ("geboostet 0,6 I + 0,3 sigma_1",
                                                                                 0.6*PAULI[0] + 0.3*PAULI[1]),
                 ("lichtartig r (I + sigma_3)", PAULI[0] + PAULI[3])):
    lam_g = np.linalg.eigvalsh(1000*Xg)
    gross[name] = (S_X_gauss(lam_g), abs(np.linalg.det(1000*Xg)))
lam_b = np.linalg.eigvalsh(0.6*PAULI[0] + 0.3*PAULI[1])                    # det X = 0,27 (geboostet)
lam_r = np.linalg.eigvalsh(np.sqrt(0.27)*PAULI[0])                           # gleiches det X im Ruhesystem
S_b1, S_r1 = S_X_gauss(lam_b), S_X_gauss(lam_r)
S_null1, S_null0 = S_X_gauss(np.linalg.eigvalsh(PAULI[0] + PAULI[3])), S_X_gauss(np.array([0.0, 0.0]))
ZAHLEN["orte_r=1"] = {"geboostet det 0,27": f"{S_b1:.6f}", "Ruhesystem det 0,27": f"{S_r1:.6f}",
                      "lichtartig I + sigma_3 (det 0)": f"{S_null1:.6f}", "X = 0": f"{S_null0:.3e}"}
ZAHLEN["orte_r=1000"] = {k: {"S": f"{v[0]:.6f}", "S - 2 ln|det X|": (f"{v[0] - 2*np.log(v[1]):.6f}" if v[1] > 1e-9 else "-")}
                         for k, v in gross.items()}
check("Orte (z78): S(|X>) wächst für große X wie 2 ln|det X| + 2 + 4 ln(kappa/2) (zeit-, raumartig und geboostet bei r = "
      "1000 auf 1e-4), ist bei endlichem X keine Funktion von det X (bei |X| ~ 1 geboostet 0,65 gegen Ruhesystem 0,50 bei "
      "gleichem det X = 0,27; lichtartig S = 1,38 bei det X = 0) und für kleine X 2 pi kappa^2 (x_0^2 + |x|^2), die Euklidnorm "
      "im KMS-Ruhesystem: logarithmisch, kein Flächengesetz in s^2 = det X",
      all(abs(v[0] - 2*np.log(v[1]) - const_gross) < 1e-4 for k, v in gross.items() if v[1] > 1e-9)
      and abs(gross["zeitartig r I"][0] - gross["raumartig r sigma_3"][0]) < 1e-12 and gross["lichtartig r (I + sigma_3)"][1] < 1e-9
      and gross["zeitartig r I"][0]/1e6 < 1e-4 and abs(S_b1 - S_r1) > 0.1 and S_null1 - S_null0 > 1
      and abs(rest_X["1e-4"]) < 1e-4)

# --- 4.5 Typfrage und Dimension
lev_int = all((4*(sp.Rational(Q*Q, 2) + (al - sp.Rational(1, 2))*Q + N)).is_integer
              for al in (sp.Rational(1, 4), sp.Rational(3, 4)) for Q in range(-6, 7) for N in range(10))
spur = prod_q**2
check("Typfrage: spec K in pi N und Tr e^{-K} = sqrt2 e^{-pi/12} endlich (Produkt über ungerade l): der Originalfaktor ist Typ I, "
      "der Modulfluss periodisch (Periode 2), S endlich. Eine Horizontalgebra (Rindlerkeil) ist Typ III_1 mit kontinuierlichem "
      "Modulspektrum und flächendivergenter Entropie: der Schnitt Original|Partner ist ein Thermofeldschnitt, kein Horizont",
      lev_int and abs(spur - mp.sqrt(2)*mp.exp(-PI/12)) < TOL)


def S_4d(beta, N):
    n = np.arange(1, N + 1, dtype=float)
    x = np.exp(-beta*n)
    return float(np.sum(n**2*(-np.log1p(-x) + beta*n*x/(1 - x))))


sb = {b: S_4d(b, 6000)*b**3/(4*np.pi**4/45) for b in (0.2, 0.1, 0.05)}
S4_nativ = S_4d(4*np.pi, 40)
ZAHLEN["S2_lesart_4D"] = {"beta^3 S/(4 pi^4/45)": {str(b): f"{v:.8f}" for b, v in sb.items()}, "S bei beta = 4 pi": f"{S4_nativ:.6e}"}
check("S2-Lesart (Stufe m mit m^2 Zuständen = konformer Skalar auf R x S^3, z73): die Thermofeldentropie ist ein Volumengesetz, "
      "beta^3 S -> 4 pi^4/45 (Stefan-Boltzmann, Volumen 2 pi^2 R^3; Abweichung 4e-7 bei beta = 0,05), am nativen Punkt "
      f"{S4_nativ:.2e}: auch 4D gelesen entsteht keine Fläche", abs(sb[0.05] - 1) < 1e-6 and sb[0.2] < sb[0.1] < sb[0.05] < 1
      and S4_nativ < 1e-4)
befund("Kopien: N unabhängige Nähte tragen N ln 2/2 (Additivität). Das ist das einzige feste Verhältnis, eine Zählung ohne "
       "Geometrie; die Fläche je Naht (z68: A = 4 S l_P^2) liefert das Original nicht, und die E8-Naht bricht schon die "
       "Zählung (ln 12 statt 16 ln 2/4)", True)
befund("c3 als Flächeneinheit: S ist pi-frei, S c3^{+-1} ist nach Gelfond-Schneider irrational, und PSLQ findet keine "
       "Relation von S, S_n oder kappa_2 mit c3, 48 c3^4 oder phi0 (Teil 2); eine Fläche in Einheiten von c3 legt S/A nur mit "
       "einer neuen, nicht nativen Zahl fest", True)
befund("Eine Fläche entsteht aus einem Thermofeldschnitt nur, wenn er als Horizont gelesen wird (optisches Volumen ~ A/eps^2). "
       "Das verbietet GRAV.NONCIRCULAR.01, und es passt nicht zum Typ I mit diskretem Modulspektrum. Abbruchkriterium 2 von S5 "
       "greift damit für die Flächennormierung", True)

# =================================================================== Teil 5
say("\nTeil 5: Urteil für S5")
STATUS = {
    "Verschränkung von Omega": ("E", "ln 2/2 je Naht (z68), neu: Ladung/Oszillator-Aufteilung in geschlossener Form, "
                                     "Oszillatoranteil 2,7e-4, Rényi mit S-Dualität, kappa_n = r_n pi^2 theta_3(i)^8"),
    "Relationen zu TFPT-Daten": ("E", "ln 2 (Ramond-Entartung bei k^2 = 1/2), ln 12 = 16 ln 2/4 + ln(3/4) (E8 = 16 Quarter mit "
                                      "GSO); keine Relation zu 1/c3, c3, 48 c3^4, phi0, eta_B, t (Gelfond-Schneider, "
                                      "Nesterenko, PSLQ, Basislinie)"),
    "gleiche Uhr / universelle Kopplung erster Ordnung": ("E", "beta = 4 pi in allen Sektoren, Eötvös-Parameter 0, erster "
                                                               "Hauptsatz je Sektor, Translationen eingeschlossen"),
    "Universalität zweiter Ordnung": ("E (negativ)", "Fisher 724:1, Linearitätsfenster sektorabhängig"),
    "native Fläche mit festem S/Fläche": ("E (negativ)", "Moden, Nahtlänge, c, Orte, c3 nicht fest; Kopien nur als "
                                                       "Zählung"),
    "Schnitt Original|Partner als Horizont": ("E (negativ)", "Typ I, diskretes Modulspektrum; zudem GRAV.NONCIRCULAR.01"),
    "S als Majorana-Bindung (ln sqrt2)": ("L", "Zählung von Bindungen, keine Geometrie"),
    "Flächennormierung, G, eta = 1/4": ("O", "braucht einen räumlichen Schnitt in einem rekonstruierten 4D-Netz"),
}
for k_, (s_, t_) in STATUS.items():
    say(f"   [{s_}] {k_}: {t_}")
ZAHLEN["status"] = STATUS
befund("S5 gegen seine Kriterien: Erfolg (native Entropiedichte trifft 1/c3 und die Flächennormierung ohne BW) nicht erreicht; "
       "Abbruchkriterium 1 (Verhältnis nicht fest) für jede geometrische Flächenwahl erfüllt; Abbruchkriterium 2 (Horizontlesart "
       "nötig) für die Flächennormierung erfüllt. Erhalten bleibt die exakte gleiche Uhr als Randbedingung für jede "
       "Gravitationsquelle (S4)", True)
urteil = ("S5 ist in der gestellten Form gescopt verworfen. Die Verschränkung der Purifikation Omega zwischen Original und "
          "Tomita-Partner ist exakt und BW-frei berechenbar: ln 2/2 je Naht (z68), davon 99,97 % im Ladungssektor und nur "
          "9,5e-5 in den Stromoszillatoren, die die Raumzeitkinematik tragen; die Teilsummen haben geschlossene Formen in "
          "theta_3(i), eta(2i), die Rényi-Entropien erfüllen S_{1/n} + n S_n = (n+1) ln 2/2, und die Kumulanten des "
          "Modulhamiltonoperators sind kappa_n = r_n pi^2 theta_3(i)^8 mit kappa_3 = 3 kappa_2 aus der S-Dualität. Exakte "
          "Relationen zu TFPT-Daten gibt es nur mit ln 2 (Ramond-Entartung am lemniskatischen Punkt) und ln 12 = 16 (ln 2/4) "
          "+ ln(3/4) (E8-Naht = 16 Quarter mit GSO); zu 1/c3, c3, 48 c3^4, phi0, eta_B = 1/4 und t = 1 gibt es keine "
          "(Gelfond-Schneider, Nesterenko, PSLQ, Basislinie), <K> = 1/(96 c3) ist bloße pi-Kommensurabilität. Positiv und "
          "exakt ist die gleiche Uhr: beta = 4 pi in allen Sektoren (Eötvös-Parameter 0), der erste Hauptsatz gilt sektorweise "
          "mit derselben Kopplung, auch für die nativen Translationen e^{iP(X)} (Fenster |X| << e^{-2 pi}); das ist universelle "
          "Kopplung erster Ordnung an die modulare Uhr, aber ohne Flächenterm, ohne BW und ohne G, und zweiter Ordnung "
          "sektorabhängig (Fisher 724:1). Eine native Fläche mit festem S/Fläche gibt es nicht: Moden (abschneideabhängig), "
          "Nahtlänge (nur extensiv fest, der native Punkt hat C_E/S = 2,30), Zentralladung (E8 gegen D8 bei c = 8), Orte "
          "(2 ln|det X|) scheitern; fest ist nur die Zählung je Naht. Strukturgrund: Original|Partner ist ein Thermofeldschnitt "
          "einer Typ-I-Algebra mit diskretem Modulspektrum; seine Entropie ist thermisch (1D extensiv, in der S2-Lesart ein "
          "Volumengesetz), eine Fläche entstünde nur mit der verbotenen und hier unpassenden Horizontlesart. Die "
          "Flächennormierung braucht einen räumlichen Schnitt in einem rekonstruierten 4D-Netz.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z86_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
