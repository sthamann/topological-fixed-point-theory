"""Z137: Induzierte Gravitation -- ist die Gravitationskopplung (Prämisse P_G) im H6-Kandidaten eine Zählung statt eines
Parameters?

Leitlinie 'Rekonstruktion statt Auswahl'. Firewalls: DIMENSION.UPLIFT.FIREWALL.01 (die 4D-Lesart der Naht ist die radiale
Leiter von z83; ein lokales 3+1-Netz ist Prämisse, keine Herleitung), GRAV.NONCIRCULAR.01 (kein G wird vorausgesetzt),
Antinumerologie (eine passende Zahl ohne Mechanismus ist kein Ergebnis). Gelöste Blöcke (Flavour, Massen, alpha, 16, drei
Familien, E8) und H6 bleiben fest; der Ledger (COSMO.LAM.01, LAMBDA.BRANCH.01, Status 'Numerical') wird nur als
Vergleichsziel gelesen; nichts wird befördert, keine vollständige Theorie wird behauptet.
Bestand: z83 (Schalenleiter = masseloser konformer Skalar, Kern det(I - Z^+ W)^-1), z75/z99 (Flato-Fronsdal: zwei Leitern
geben die erhaltenen Ströme Delta = s + 2, Helizität +-2 nur als Komposita), z124 (Hypothese T: native Doppelspur bei n = 4
= konforme Supergravitation bei k = 1), z118 (Uhrspur = MacDowell-Mansouri, 64 je Index, Lambda = -12/R_Vak^2, G nur mit
dem gesetzten Gewicht 2 c_3), z98 (l = R_Vak/2, b = l^2/(64 pi G)), z27 (Sakharov gegen Jacobson: nur der log-Teil
universell), z121 (Kette = Destri-de-Vega-Backsteingitter), z68/z86 (Entropie ln 2/2 je Naht), z104 (positive Zeit =>
Lambda <= 0). Nicht wiederholt: Amplituden (z134) und Verschränkungsgleichgewicht (z136).
Statuskürzel: E exakt, C bedingt (Prämisse genannt), L Lesart bzw. Literaturformel, O offen.

Konventionen (durchgehend): euklidisch, d = 4, S_4 = 2 pi^2. <T_mn(x) T_rs(0)> = C_T I_mn,rs(x)/x^8 (Osborn-Petkou),
<T^m_m> = (c W^2 - a E_4)/(16 pi^2), c = pi^4 C_T/40, W^2 = Weyl-Quadrat (volle Kontraktion), E_4 = Riem^2 - 4 Ric^2 + R^2,
Int E_4 = 32 pi^2 chi. S_EH = -(1/(16 pi G)) Int (R - 2 Lambda), Lambda = -+3/l^2. L = ln(Abschneideskala/IR-Skala).

Teil 1  Der native Spannungstensor: erhaltene Delta = 4-Spin-2-Bilineare zweier Leitern aus dem nativen Kern, C_T aus der
        Ward-Identität; Gegenprobe im Ortsraum (Wick); Vergleich Skalar, Weyl, Maxwell; Negativkontrollen.
Teil 2  Induzierte Wirkung: a und c je Lesart (i) Einteilcheninhalt einer Naht, (ii) 3 x 16 Weyl, (iii) dazu Eichbosonen;
        universeller Weyl^2-Term (log), nicht universeller Einstein-Hilbert-Term (z27).
Teil 3  Vom Weyl-Koeffizienten zu G Lambda: Maldacena-Relation über Euler/renormiertes Volumen, geprüft an AdS4, S^4 und
        Hawking-Page; reine Zahlen l^2/G und G Lambda je Lesart; Vergleich mit z118 (64 je Index, Gewicht 2 c_3).
Teil 4  Ist die Kopplung eine Nahtzahl? Skalierung von Uhrspurgewicht, Entropie und C_T mit N (Tensorprodukt und Kette
        von z121); Vergleichsziel COSMO.LAM.01 / LAMBDA.BRANCH.01 (l^2/G = 32 pi^3 e^(2/alpha)).
Teil 5  Bilanz und Urteil zu P_G.
"""
import itertools
import json
from fractions import Fraction as Fr
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sp
from numpy.polynomial.legendre import leggauss
from scipy.linalg import block_diag

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
rng = np.random.default_rng(137)
mp.mp.dps = 40
PI = sp.pi
ALPHA_INV = mp.mpf("137.035999177")                                  # CODATA 2022, nur Vergleichsgröße (wie z98)


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


def lies(name):
    return json.loads((HERE / name).read_text())


def dl(a, b):
    return 1 if a == b else 0


Z27, Z98, Z118, Z121 = lies("z27_ergebnis.json"), lies("z98_ergebnis.json"), lies("z118_ergebnis.json"), lies("z121_ergebnis.json")

# =================================================================== Teil 1
say("Teil 1: Der native Spannungstensor und sein Zweipunktkoeffizient")
D4 = 4
S4 = 2*PI**2                                                         # Fläche der S^3

# --- 1.1 Der native Kern und die Leiter
xs, ys = sp.symbols("x0:4", real=True), sp.symbols("y0:4", real=True)
SIG = [sp.eye(2), sp.I*sp.Matrix([[0, 1], [1, 0]]), sp.I*sp.Matrix([[0, -sp.I], [sp.I, 0]]), sp.I*sp.Matrix([[1, 0], [0, -1]])]
Zm = sum((xs[i]*SIG[i] for i in range(4)), sp.zeros(2, 2))
Wm = sum((ys[i]*SIG[i] for i in range(4)), sp.zeros(2, 2))
xy, xx, yy = sum(a*b for a, b in zip(xs, ys)), sum(a*a for a in xs), sum(b*b for b in ys)
kern_rest = sp.expand((sp.eye(2) - Zm.H*Wm).det() - (1 - 2*xy + xx*yy))


def monome(deg):
    return [m for m in itertools.product(range(deg + 1), repeat=D4) if sum(m) == deg]


def einzelgram(Delta, maxdeg):
    """<P^b phi|P^a phi> = d_u^a d_v^b (1 - 2 u.v + u^2 v^2)^-Delta bei 0 (tiefster Zustand auf 1 normiert)."""
    u, v, t = sp.symbols("u0:4"), sp.symbols("v0:4"), sp.Symbol("t")
    uv, uu, vv = sum(a*b for a, b in zip(u, v)), sum(a*a for a in u), sum(b*b for b in v)
    ser = sp.series((1 - 2*t*uv + t**2*uu*vv)**(-Delta), t, 0, maxdeg + 1).removeO()
    G = {}
    for k in range(maxdeg + 1):
        for mon, c in sp.Poly(sp.expand(ser.coeff(t, k)), *u, *v).terms():
            a, b = mon[:4], mon[4:]
            G[(a, b)] = sp.factor(c)*sp.prod([sp.factorial(i) for i in a + b])
    return G


def e_(i):
    return tuple(1 if j == i else 0 for j in range(4))


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


NULL4 = (0, 0, 0, 0)


def basis2(k):
    """Zwei Leitern: Monome x1^a x2^b vom Gesamtgrad k (Zustände P_1^a P_2^b |phi>|phi>)."""
    return [(a, b) for da in range(k + 1) for a in monome(da) for b in monome(k - da)]


B1, B2, B3 = basis2(1), basis2(2), basis2(3)
I2x, I3x = {b: i for i, b in enumerate(B2)}, {b: i for i, b in enumerate(B3)}


def gram2(B, G1):
    return sp.Matrix(len(B), len(B), lambda i, j: G1.get((B[j][0], B[i][0]), 0)*G1.get((B[j][1], B[i][1]), 0))


def Pges(Bk, Bk1, mu):
    """Gesamtimpuls P_mu = P^(1)_mu + P^(2)_mu als Multiplikation mit x1_mu + x2_mu."""
    idx = {b: i for i, b in enumerate(Bk1)}
    M = sp.zeros(len(Bk1), len(Bk))
    for j, (a, b) in enumerate(Bk):
        M[idx[(add(a, e_(mu)), b)], j] += 1
        M[idx[(a, add(b, e_(mu)))], j] += 1
    return M


P12 = [Pges(B1, B2, m) for m in range(4)]
P23 = [Pges(B2, B3, m) for m in range(4)]


def Tvek(m, n, r):
    """x1_m x1_n + x2_m x2_n - r (x1_m x2_n + x1_n x2_m) als Vektor der Stufe 2."""
    v = sp.zeros(len(B2), 1)
    v[I2x[(add(e_(m), e_(n)), NULL4)]] += 1
    v[I2x[(NULL4, add(e_(m), e_(n)))]] += 1
    v[I2x[(e_(m), e_(n))]] -= r
    v[I2x[(e_(n), e_(m))]] -= r
    return v


def leiterdaten(Delta):
    G1 = einzelgram(Delta, 3)
    raenge = []
    for k in range(4):
        Mk = monome(k)
        raenge.append(sp.Matrix(len(Mk), len(Mk), lambda i, j: G1.get((Mk[j], Mk[i]), 0)).rank())
    G2, G3 = gram2(B2, G1), gram2(B3, G1)
    zeilen = sp.Matrix.vstack(*[P12[mu].T*G2 for mu in range(4)])   # <P_mu s|t> = 0 für alle s der Stufe 1: K_mu t = 0
    kern = zeilen.nullspace()
    Nmat = sp.Matrix.hstack(*kern)
    prim_dim = (Nmat.T*G2*Nmat).rank()                               # Primäre modulo Nullvektoren
    rs = sp.Symbol("r")
    r_wert = sp.solve(list(zeilen*Tvek(0, 1, rs)), rs)[rs]
    Tl = {}
    spur = sum((Tvek(a, a, r_wert) for a in range(4)), sp.zeros(len(B2), 1))/4
    for m in range(4):
        for n in range(4):
            Tl[(m, n)] = Tvek(m, n, r_wert) - (spur if m == n else sp.zeros(len(B2), 1))
    alle_prim = all(zeilen*Tl[k] == sp.zeros(zeilen.shape[0], 1) for k in sorted(Tl))
    rang_T = (sp.Matrix.hstack(*[Tl[k] for k in sorted(Tl)]).T*G2*sp.Matrix.hstack(*[Tl[k] for k in sorted(Tl)])).rank()
    n12 = (Tl[(0, 1)].T*G2*Tl[(0, 1)])[0]
    n00 = (Tl[(0, 0)].T*G2*Tl[(0, 0)])[0]
    div = [sum((P23[m]*Tl[(m, n)] for m in range(4)), sp.zeros(len(B3), 1)) for n in range(4)]
    ndiv = [(d_.T*G3*d_)[0] for d_ in div]
    return {"G1": G1, "G2": G2, "G3": G3, "raenge": raenge, "rang2": G2.rank(), "kern": len(kern), "prim": prim_dim,
            "r": r_wert, "T": Tl, "alle_prim": alle_prim, "rang_T": rang_T, "n12": n12, "n00": n00, "ndiv": ndiv,
            "zeilen": zeilen}


LD = {Dl: leiterdaten(Dl) for Dl in (sp.Integer(1), sp.Rational(5, 4), sp.Rational(3, 2), sp.Integer(2), sp.Integer(3))}
L1 = LD[sp.Integer(1)]
Ds = sp.Symbol("Delta", positive=True)
Gsym = einzelgram(Ds, 2)
norm_box = sp.factor(sum(Gsym[(add(e_(a), e_(a)), add(e_(b), e_(b)))] for a in range(4) for b in range(4)))   # <P^2 phi|P^2 phi>
ZAHLEN["leiter"] = {"det(I - Z^+ W) - (1 - 2 x.y + x^2 y^2)": str(kern_rest),
                    "Ränge der Stufen 0..3 (Delta = 1)": L1["raenge"], "Ränge (Delta = 3/2)": LD[sp.Rational(3, 2)]["raenge"],
                    "<P^2 phi|P^2 phi>(Delta)": str(norm_box)}
check(f"[E] Native Leiter (z83): der Kern det(I - Z^+ W)^-1 mit Z = x_0 + i x.sigma ist (1 - 2 x.y + x^2 y^2)^-1 (Rest "
      f"{kern_rest}), also der Zweipunktkern einer skalaren Primären mit Delta = 1. Gram der Zustände P^a|phi> aus seinen "
      f"Taylorkoeffizienten: Ränge der Stufen 0..3 = {L1['raenge']} = (k + 1)^2 (z83: Summe m^2 q^m), P^2|phi> ist "
      f"Nullvektor. Die Norm <P^2 phi|P^2 phi> = {norm_box} verschwindet für Delta > 0 nur bei Delta = 1: die nativen "
      f"Vielfachheiten legen Delta = 1 fest. Negativkontrolle Delta = 3/2: Ränge {LD[sp.Rational(3, 2)]['raenge']} (keine "
      f"Nullvektoren, lange Darstellung)",
      kern_rest == 0 and L1["raenge"] == [1, 4, 9, 16] and LD[sp.Rational(3, 2)]["raenge"] == [1, 4, 10, 20]
      and sp.solve(norm_box, Ds) == [1])

# --- 1.2 Die erhaltene Spin-2-Bilineare
sym_ok = all(L1["T"][(m, n)] == L1["T"][(n, m)] for m in range(4) for n in range(4))
tausch = sp.zeros(len(B2), len(B2))
for j, (a, b) in enumerate(B2):
    tausch[I2x[(b, a)], j] = 1
tausch_ok = all(tausch*L1["T"][k] == L1["T"][k] for k in sorted(L1["T"]))
ZAHLEN["spin2"] = {"Rang Stufe 2 (von 36)": L1["rang2"], "dim Kern von K": L1["kern"], "Primäre modulo Null": L1["prim"],
                   "r in x1x1 + x2x2 - r (x1x2 + x2x1)": str(L1["r"]), "Rang der T_mn": L1["rang_T"],
                   "n_00/n_01": str(L1["n00"]/L1["n12"]), "Norm der Divergenz": [str(x) for x in L1["ndiv"]]}
check(f"[E] Zwei native Leitern, Stufe 2 (E = 4): Gram-Rang {L1['rang2']} von 36 (zwei Nullvektoren P_i^2), Kern von "
      f"K_ges (orthogonal zu P_ges mal Stufe 1) Dimension {L1['kern']}, modulo Null {L1['prim']} = eine symmetrisch-"
      f"spurfreie Spin-2-Primäre (Rang der neun T_mn: {L1['rang_T']}). Die Primärbedingung legt die Bilineare fest: "
      f"T_mn ~ x1_m x1_n + x2_m x2_n - {L1['r']} (x1_m x2_n + x1_n x2_m) - Spur, also P_1 P_1 - 4 P_1 P_2 + P_2 P_2 "
      f"(Gegenbauer s = 2), d. h. T ~ (2/3) d phi d phi - (1/3) phi dd phi - Spur: der verbesserte Spannungstensor. "
      f"Symmetrisch im Leitertausch ({tausch_ok}), Verhältnis der Normen n_00/n_01 = {L1['n00']/L1['n12']} = "
      f"Pi_00,00/Pi_01,01 (Spin-2-Projektor)",
      L1["rang2"] == 34 and L1["kern"] == 11 and L1["prim"] == 9 and L1["r"] == 2 and L1["alle_prim"] and L1["rang_T"] == 9
      and sym_ok and tausch_ok and L1["n00"]/L1["n12"] == sp.Rational(3, 2))

# --- 1.3 Erhaltung und Negativkontrollen
verh = {str(Dl): sp.nsimplify(LD[Dl]["ndiv"][0]/LD[Dl]["n12"]) for Dl in LD}
gesetz = all(sp.simplify(LD[Dl]["ndiv"][0]/LD[Dl]["n12"] - 9*((2*Dl + 2) - 4)) == 0 for Dl in LD)
iso = all(len(set(LD[Dl]["ndiv"])) == 1 for Dl in LD)
t_unv = sp.zeros(len(B2), 1)                                         # unverbesserte Bilineare d_0 phi d_1 phi
t_unv[I2x[(e_(0), e_(1))]] += 1
t_unv[I2x[(e_(1), e_(0))]] += 1
t_nach = (P12[0]*sp.Matrix([1 if b in ((e_(1), NULL4), (NULL4, e_(1))) else 0 for b in B1]))   # (P_ges)_0 (P_ges)_1 |phi phi>
k_unv = (L1["zeilen"]*t_unv).norm() != 0
k_nach = (L1["zeilen"]*t_nach).norm() != 0
ZAHLEN["erhaltung"] = {"|div T|^2/|T_01|^2 je Delta": {k: str(v) for k, v in verh.items()},
                       "E der Spin-2-Primären je Delta": {str(Dl): str(2*Dl + 2) for Dl in LD},
                       "r je Delta": {str(Dl): str(LD[Dl]["r"]) for Dl in LD}}
check(f"[E] Erhaltung: die vier Zustände P_ges^m T_mn (Stufe 3) haben Norm {[str(x) for x in L1['ndiv']]}: die native "
      f"Bilineare ist erhalten und hat E = Delta = 4 = s + 2 (Flato-Fronsdal, z75/z99). Negativkontrolle über die "
      f"Leiterdimension: für Delta_phi = 5/4, 3/2, 2, 3 existiert die Spin-2-Primäre ebenfalls (r = "
      f"{[str(LD[Dl]['r']) for Dl in LD]}), liegt aber bei E = 2 Delta_phi + 2 = "
      f"{[str(2*Dl + 2) for Dl in LD]} und ist nicht erhalten: |div T|^2/|T_01|^2 = {verh} = 9 (Delta_T - 4) "
      f"({gesetz}, alle vier Komponenten gleich: {iso}). Eine nicht erhaltene Spin-2-Bilineare hat nie Delta = 4. "
      f"Zweite Kontrolle bei Delta_phi = 1: die unverbesserte Bilineare d phi d phi und der Nachkomme dd(phi^2) liegen "
      f"auch bei E = 4, werden von K aber nicht annulliert ({k_unv}, {k_nach})",
      all(x == 0 for x in L1["ndiv"]) and gesetz and iso and k_unv and k_nach
      and all(LD[Dl]["ndiv"][0] != 0 for Dl in LD if Dl != 1) and all(LD[Dl]["alle_prim"] for Dl in LD))

# --- 1.4 C_T aus der Ward-Identität (Zustandsbild)
t01 = L1["T"][(0, 1)]
yv1, yv2 = sp.symbols("p0:4"), sp.symbols("q0:4")
ueber = 0                                                            # <e^(P1.y1) e^(P2.y2) phi phi | t_01>
for j, (a, b) in enumerate(B2):
    if t01[j] == 0:
        continue
    for (a2, b2) in B2:
        g_ = L1["G1"].get((a, a2), 0)*L1["G1"].get((b, b2), 0)
        if g_ != 0:
            ueber += t01[j]*g_*sp.prod([yv1[i]**a2[i]/sp.factorial(a2[i]) for i in range(4)])*sp.prod(
                [yv2[i]**b2[i]/sp.factorial(b2[i]) for i in range(4)])
o0 = sp.simplify(ueber/((yv1[0] - yv2[0])*(yv1[1] - yv2[1])))
n0 = L1["n12"]
winkel = sp.Rational(3, 4)*S4                                        # Int dOmega n^m n^n (n_m n_n - delta_mn/4)
aW = 1/winkel                                                        # Ward: <T| phi phi> = -(d Delta/((d-1) S_d)) N (y12 y12 - Spur)
CT_zustand = sp.simplify(2*aW**2*(2*n0)/(2*o0)**2)                   # identische Felder: Norm und Überlapp je doppelt
CT_skalar = sp.Rational(D4, D4 - 1)/S4**2
c_aus_CT = sp.simplify(PI**4*CT_zustand/40)
ZAHLEN["C_T_zustand"] = {"Überlapp o_0 (Koeffizient von y12_0 y12_1)": str(o0), "Norm n_0": str(n0),
                         "Ward-Koeffizient d Delta/((d-1) S_d)": str(aW), "C_T": str(CT_zustand), "c = pi^4 C_T/40": str(c_aus_CT)}
check(f"[E] C_T aus der nativen Normierung: der Überlapp von T_01 mit dem Zweiteilchen-Kohärenzzustand ist "
      f"{o0} (y1 - y2)_0 (y1 - y2)_1 (hängt nur von y_12 ab, wie die Dreipunktstruktur verlangt), die Norm {n0} (N = 1). "
      f"Die Ward-Identität (T erzeugt die Dilatation mit Delta = 1; Winkelintegral {winkel}) normiert T unabhängig von "
      f"der Feldnormierung: C_T = 2 a_W^2 n/o^2 = {CT_zustand} = d/((d - 1) S_d^2) = {CT_skalar}, c = pi^4 C_T/40 = "
      f"{c_aus_CT}. Je reeller nativer Leiter also genau der Wert des freien reellen Skalars",
      o0.is_number and sp.simplify(CT_zustand - CT_skalar) == 0 and c_aus_CT == sp.Rational(1, 120)
      and sp.simplify(CT_zustand - 1/(3*PI**4)) == 0)

# --- 1.5 Gegenprobe im Ortsraum (Wick)
X = sp.symbols("X0:4", real=True)
Nn = sp.Symbol("N", positive=True)
Gp = Nn/sum(x*x for x in X)
_gc = {}


def gd(idx):
    k = tuple(sorted(idx))
    if k not in _gc:
        f = Gp
        for i in k:
            f = sp.diff(f, X[i])
        _gc[k] = f
    return _gc[k]


def Istr(m, n, r, s, pt):
    x2 = sum(p*p for p in pt)
    I_ = lambda a, b: dl(a, b) - sp.Rational(2*pt[a]*pt[b], x2)
    return sp.Rational(1, 2)*(I_(m, r)*I_(n, s) + I_(m, s)*I_(n, r)) - sp.Rational(1, 4)*dl(m, n)*dl(r, s)


def T_skalar(m, n, xi):
    """T_mn = d_m phi d_n phi - (1/2) delta (d phi)^2 - xi (d_m d_n - delta Box) phi^2 auf der Schale: (Koeff., A, B)."""
    terme = [(1 - 2*xi, (m,), (n,)), (-2*xi, (), (m, n))]
    terme += [((-sp.Rational(1, 2) + 2*xi)*dl(m, n), (a,), (a,)) for a in range(4)]
    return [t for t in terme if t[0] != 0]


def TT_skalar(m, n, r, s, pt, xi):
    tot = 0
    for c1, A, B in T_skalar(m, n, xi):
        for c2, C, Dd in T_skalar(r, s, xi):
            tot += c1*c2*(-1)**(len(C) + len(Dd))*(gd(A + C)*gd(B + Dd) + gd(A + Dd)*gd(B + C))
    return sp.nsimplify(sp.simplify(tot.subs(dict(zip(X, pt)))))


PUNKTE = [(1, 2, 3, 4), (2, -1, 1, 3)]
KOMP = [(0, 1, 0, 1), (0, 0, 0, 0), (0, 1, 2, 3), (0, 0, 1, 1)]


def ct_tabelle(fn, punkte=PUNKTE, komp=KOMP):
    out = []
    for pt in punkte:
        x2 = sum(p*p for p in pt)
        for cmp_ in komp:
            out.append(sp.simplify(fn(*cmp_, pt)/(Istr(*cmp_, pt)/x2**4)))
    return out


xi_k = sp.Rational(1, 6)
ct_s = ct_tabelle(lambda m, n, r, s, pt: TT_skalar(m, n, r, s, pt, xi_k))
ct_s0 = ct_tabelle(lambda m, n, r, s, pt: TT_skalar(m, n, r, s, pt, 0), punkte=PUNKTE[:1])
N_kan = 1/(4*PI**2)
CT_ort = sp.simplify(ct_s[0].subs(Nn, N_kan))

# Ward-Ladung im Ortsraum: Q = -Int_{S^3_r} dS n^m x^n <T_mn(x) phi(0) phi(y)> / <phi(0) phi(y)>
Yv = (3.0, 0.0, 0.0, 0.0)
gnum = {k: sp.lambdify(X, gd(k).subs(Nn, 1), "numpy") for k in [()] + [(a,) for a in range(4)]
        + [tuple(sorted((a, b))) for a in range(4) for b in range(a, 4)]}


def ladung(r_, xi):
    th, w = leggauss(120)
    th, w = (th + 1)*np.pi/2, w*np.pi/2
    tot = 0.0
    for t_, w_ in zip(th, w):
        nvec = (np.cos(t_), np.sin(t_), 0.0, 0.0)
        x_ = tuple(r_*c for c in nvec)
        xm = tuple(a - b for a, b in zip(x_, Yv))
        f = 0.0
        for m in (0, 1):
            for n in (0, 1):
                for c1, A, B in T_skalar(m, n, xi):
                    kA, kB = tuple(sorted(A)), tuple(sorted(B))
                    f += nvec[m]*nvec[n]*float(c1)*(gnum[kA](*x_)*gnum[kB](*xm) + gnum[kA](*xm)*gnum[kB](*x_))
        tot += w_*4*np.pi*np.sin(t_)**2*r_**4*f
    return -tot/gnum[()](*Yv)


Q_k = [ladung(r_, xi_k) for r_ in (0.5, 1.0, 2.0)]
Q_0 = [ladung(r_, 0) for r_ in (0.5, 1.0, 2.0)]
N_ward = 1/(4*np.pi**2)
ZAHLEN["C_T_ortsraum"] = {"<TT> x^8/I, xi = 1/6 (2 Punkte x 4 Komponenten)": sorted({str(v) for v in ct_s}),
                          "xi = 0": [str(v) for v in ct_s0], "C_T bei N = 1/(4 pi^2)": str(CT_ort),
                          "Ladung Q/(4 pi^2 N), xi = 1/6, r = 0,5/1/2": [float(round(q/(4*np.pi**2), 10)) for q in Q_k],
                          "|Ladung|, xi = 0": [float(round(abs(q)/(4*np.pi**2), 8)) for q in Q_0]}
check(f"[E] Gegenprobe im Ortsraum (Wick, sympy exakt an zwei Punkten, vier Komponenten): mit dem verbesserten "
      f"Spannungstensor (xi = 1/6) ist <T_mn(x) T_rs(0)> x^8/I_mn,rs = {sorted({str(v) for v in ct_s})} für alle acht "
      f"Auswertungen, also eine einzige Konstante C_T = (16/3) N^2; die Dilatationsladung Q = -Int_(S^3) n x T auf phi(0) "
      f"ist unabhängig vom Radius und gleich 4 pi^2 N Delta ({[float(round(q/(4*np.pi**2), 8)) for q in Q_k]} in Einheiten "
      f"4 pi^2 N bei r = 0,5, 1, 2): Delta = 1 verlangt N = 1/(4 pi^2), damit C_T = {CT_ort} wie im Zustandsbild. "
      f"Negativkontrolle xi = 0: die Verhältnisse {[str(sp.nsimplify(v/Nn**2)) for v in ct_s0]} N^2 sind keine gemeinsame "
      f"Konstante, und die Ladung ist {[float(round(abs(q)/(4*np.pi**2), 6)) for q in Q_0]} statt 1 (der unverbesserte "
      f"Tensor erzeugt nur den Bahnanteil x.d, nicht Delta)",
      len({str(v) for v in ct_s}) == 1 and sp.simplify(ct_s[0] - sp.Rational(16, 3)*Nn**2) == 0
      and max(abs(q*N_ward - 1) for q in Q_k) < 1e-9 and sp.simplify(CT_ort - CT_zustand) == 0
      and len({str(v) for v in ct_s0}) > 1 and max(abs(q*N_ward) for q in Q_0) < 1e-9)


# --- 1.6 Vergleich: Weyl-Fermion und Maxwell-Feld in derselben Konvention
def FF(m, n, r, s):
    """<F_mn(x) F_rs(0)> aus <A_m A_n> = delta_mn G (Feynman-Eichung; eichinvariant)."""
    return -(dl(n, s)*gd((m, r)) - dl(n, r)*gd((m, s)) - dl(m, s)*gd((n, r)) + dl(m, r)*gd((n, s)))


def T_maxwell(m, n):
    t = [(1, (m, a), (n, a)) for a in range(4)]
    if m == n:
        t += [(-sp.Rational(1, 4), (a, b), (a, b)) for a in range(4) for b in range(4)]
    return t


def TT_maxwell(m, n, r, s, pt):
    tot = 0
    for c1, A, B in T_maxwell(m, n):
        for c2, C, Dd in T_maxwell(r, s):
            tot += c1*c2*(FF(*A, *C)*FF(*B, *Dd) + FF(*A, *Dd)*FF(*B, *C))
    return sp.nsimplify(sp.simplify(tot.subs(dict(zip(X, pt)))))


s0_, sx_, sy_, sz_ = sp.eye(2), sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -sp.I], [sp.I, 0]]), sp.Matrix([[1, 0], [0, -1]])
kp = sp.kronecker_product
GAM = [kp(sx_, sx_), kp(sx_, sy_), kp(sx_, sz_), kp(sy_, s0_)]
G5 = GAM[0]*GAM[1]*GAM[2]*GAM[3]
cliff_ok = all(GAM[a]*GAM[b] + GAM[b]*GAM[a] == 2*dl(a, b)*sp.eye(4) for a in range(4) for b in range(4)) and G5*G5 == sp.eye(4)
Sprop = sp.zeros(4, 4)
for a in range(4):
    Sprop += -GAM[a]*sp.diff(Gp, X[a])                                # S = -gamma.d G = gamma.x/(2 pi^2 x^4) bei N = 1/(4 pi^2)
Sminus = Sprop.subs({x: -x for x in X}, simultaneous=True)


def dM(M, i):
    return M.applyfunc(lambda e: sp.diff(e, X[i]))


def TT_dirac(m, n, r, s, pt, proj=None):
    """T_mn = (1/4)[psi-quer gamma_m d_n psi - d_n psi-quer gamma_m psi + (m <-> n)]; eine Fermionschleife."""
    Pj = sp.eye(4) if proj is None else proj
    terme = lambda m_, n_: [(sg*sp.Rational(1, 4), a, b, w) for (a, b) in ((m_, n_), (n_, m_)) for sg, w in ((1, "psi"), (-1, "quer"))]
    tot = 0
    for c1, a1, b1, w1 in terme(m, n):
        for c2, a2, b2, w2 in terme(r, s):
            A, B = Sprop, Sminus                                     # A = <psi(x) psi-quer(y)>, B = <psi(y) psi-quer(x)>, z = x - y
            if w1 == "psi":
                A = dM(A, b1)
            else:
                B = dM(B, b1)
            if w2 == "psi":
                B = dM(B, b2)
            else:
                A = dM(A, b2)
            tot += c1*c2*(gam_spur(GAM[a1]*A*GAM[a2]*B*Pj))           # -1 (Schleife) mal -1 (Ableitung nach y)
    return sp.nsimplify(sp.simplify(tot.subs(dict(zip(X, pt)))))


def gam_spur(M):
    return M.trace()


ct_v = ct_tabelle(TT_maxwell, punkte=PUNKTE[:1], komp=[(0, 1, 0, 1), (0, 0, 0, 0), (0, 0, 1, 1)])
ct_d = ct_tabelle(TT_dirac, punkte=PUNKTE[:1], komp=[(0, 1, 0, 1), (0, 0, 1, 1)])
PL = (sp.eye(4) + G5)/2
ct_w = ct_tabelle(lambda m, n, r, s, pt: TT_dirac(m, n, r, s, pt, PL), punkte=PUNKTE[:1], komp=[(0, 1, 0, 1), (0, 0, 1, 1)])
einheit = sp.Rational(16, 3)*Nn**2
verh_frei = {"reeller Skalar": sp.simplify(ct_s[0]/einheit), "Weyl-Fermion": sp.simplify(ct_w[0]/einheit),
             "Dirac-Fermion": sp.simplify(ct_d[0]/einheit), "Maxwell-Feld": sp.simplify(ct_v[0]/einheit)}
C_FREI = {"skalar": Fr(1, 120), "weyl": Fr(1, 40), "vektor": Fr(1, 10)}
ZAHLEN["C_T_frei"] = {"C_T in Einheiten 1/(3 pi^4)": {k: str(v) for k, v in verh_frei.items()},
                      "c = pi^4 C_T/40": {k: str(sp.Rational(1, 120)*v) for k, v in verh_frei.items()},
                      "native Leiter": "1 (reeller Skalar)"}
check(f"[E] Freie Felder in derselben Konvention (Wick im Ortsraum, sympy exakt; Clifford-Relationen {cliff_ok}): "
      f"C_T in Einheiten des reellen Skalars 1/(3 pi^4): {{Skalar: {verh_frei['reeller Skalar']}, Weyl: "
      f"{verh_frei['Weyl-Fermion']}, Dirac: {verh_frei['Dirac-Fermion']}, Maxwell: {verh_frei['Maxwell-Feld']}}}, jeweils "
      f"eine Konstante über die geprüften Komponenten; also c = 1/120, 1/40, 1/10. Die native Leiter trägt 1 Einheit "
      f"(reeller Skalar), ein Weyl-Fermion 3, ein Eichboson 12",
      cliff_ok and len({str(v) for v in ct_v}) == 1 and len({str(v) for v in ct_d}) == 1 and len({str(v) for v in ct_w}) == 1
      and verh_frei == {"reeller Skalar": 1, "Weyl-Fermion": 3, "Dirac-Fermion": 6, "Maxwell-Feld": 12})
befund("E/C: Der native Spannungstensor existiert als erhaltene Delta = 4-Bilineare zweier Leitern mit C_T = 1/(3 pi^4) "
       "je reeller Leiter (c = 1/120), fest durch die Ward-Identität und die nativen Vielfachheiten (Delta = 1). Er ist "
       "ein zusammengesetzter Strom (z75: kein Graviton, Helizität +-2 nur als Komposita); die 4D-Lesart der Leiter als "
       "lokales Feld ist die Prämisse von z83/z90 (P_netz)", True)

# =================================================================== Teil 2
say("\nTeil 2: Die induzierte Wirkung für den nativen Inhalt")


# --- 2.1 a aus den Spektren auf S^4 (Zeta-Funktion), Einstein-Hilbert-Residuum daneben
def hurwitz_neg(n, a_):
    """zeta_H(-n, a) = -B_(n+1)(a)/(n + 1)."""
    return -sp.bernoulli(n + 1, a_)/(n + 1)


def zeta_s4(d3, d1, al2, nu0):
    """zeta(s) = Summe_(nu = nu0, nu0 + 1, ...) (d3 nu^3 + d1 nu) (nu^2 - al2)^-s: (zeta(0), Res_(s=1))."""
    z0 = d3*hurwitz_neg(3, nu0) + d1*hurwitz_neg(1, nu0) + al2*d1/2 + al2**2*d3/4
    return sp.nsimplify(z0), sp.nsimplify(d1/2 + al2*d3/2)


R_ = sp.Rational
SPEK = {"konformer Skalar": (R_(1, 3), -R_(1, 12), R_(1, 4), R_(3, 2)),      # lambda = (n+1)(n+2), Entartung (n+1)(n+2)(2n+3)/6
        "minimaler Skalar": (R_(1, 3), -R_(1, 12), R_(9, 4), R_(3, 2)),      # lambda = n(n+3)
        "Dirac D^2": (R_(4, 3), -R_(4, 3), R_(0), R_(2)),                    # lambda = (n+2)^2, Entartung 8 binom(n+3, 3)
        "transversaler Vektor": (R_(1), -R_(9, 4), R_(1, 4), R_(5, 2))}      # lambda = (n+1)(n+2), n >= 1, Entartung n(n+3)(2n+3)/2
ZS = {k: zeta_s4(*v) for k, v in SPEK.items()}
z_sk, z_min, z_dir, z_tv = ZS["konformer Skalar"], ZS["minimaler Skalar"], ZS["Dirac D^2"], ZS["transversaler Vektor"]
z_max0 = z_tv[0] + (z_min[0] - 1) - 2*z_min[0]                       # transversal + longitudinal (ohne Nullmode) - 2 Geister (lokal)
z_max1 = z_tv[1] + z_min[1] - 2*z_min[1]
A_FREI = {"skalar": Fr(-z_sk[0].p, 4*z_sk[0].q), "weyl": Fr(z_dir[0].p, 8*z_dir[0].q), "vektor": Fr(-z_max0.p, 4*z_max0.q)}
EH_EINH = {"konformer Skalar": z_sk[1]/z_min[1], "minimaler Skalar": z_min[1]/z_min[1], "Weyl": -z_dir[1]/(2*z_min[1]),
           "Vektor": z_max1/z_min[1]}


def waermespur(d3, d1, al2, nu0, t):
    tot, nu = mp.mpf(0), mp.mpf(nu0.p)/nu0.q
    al, a3, a1 = mp.mpf(al2.p)/al2.q, mp.mpf(d3.p)/d3.q, mp.mpf(d1.p)/d1.q
    while True:
        term = (a3*nu**3 + a1*nu)*mp.exp(-t*(nu**2 - al))
        tot += term
        if nu > 10 and abs(term) < mp.mpf(10)**-45:
            return tot
        nu += 1


def c0_numerisch(key):
    d3, d1, al2, nu0 = SPEK[key]
    f = lambda t: waermespur(d3, d1, al2, nu0, t) - (mp.mpf(d3.p)/d3.q/2)/t**2 - (mp.mpf(ZS[key][1].p)/ZS[key][1].q)/t
    t1 = mp.mpf("0.001")
    return 2*f(t1/2) - f(t1)                                         # Richardson gegen den O(t)-Term


c0_num = {k: c0_numerisch(k) for k in SPEK}
c0_abw = max(abs(c0_num[k] - mp.mpf(ZS[k][0].p)/ZS[k][0].q) for k in SPEK)
# lokale Wärmekernformel auf der Einheits-S^4: b_4 = (1/180)(Riem^2 - Ric^2), Riem^2 = 24, Ric^2 = 36, R = 12, Vol = 8 pi^2/3
aa, cc = sp.symbols("a c")
Riem2, Ric2, Rsk2 = sp.symbols("Riem2 Ric2 Rsk2")
W2_ = Riem2 - 2*Ric2 + Rsk2/3
E4_ = Riem2 - 4*Ric2 + Rsk2
lsg = sp.solve(sp.Poly(sp.expand(cc*W2_ - aa*E4_ - R_(1, 180)*(Riem2 - Ric2)), Riem2, Ric2).coeffs(), [aa, cc], dict=True)
lok_rest = sp.expand((cc*W2_ - aa*E4_ - R_(1, 180)*(Riem2 - Ric2)).subs(lsg[0])) if lsg else None
lok_s4 = R_(1, 180)*(24 - 36)*R_(8, 3)/16                            # Int b_4/(16 pi^2) auf S^4
ZAHLEN["anomalie_a"] = {"zeta(0) [Wärmekernkoeffizient mit Nullmoden]": {k: str(v[0]) for k, v in ZS.items()},
                        "Res_(s=1)": {k: str(v[1]) for k, v in ZS.items()}, "Maxwell lokal": str(z_max0),
                        "a": {k: str(v) for k, v in A_FREI.items()},
                        "numerisch (Wärmespur, Richardson) max. Abweichung": mp.nstr(c0_abw, 2),
                        "b_4 = (1/180)(Riem^2 - Ric^2) = c W^2 - a E_4": {"a": str(lsg[0][aa]), "c": str(lsg[0][cc])} if lsg else None}
check(f"[E] a aus den Spektren auf S^4 (Hurwitz-Zeta, exakt): zeta(0) = {z_sk[0]} (konformer Skalar), {z_dir[0]} (D^2 des "
      f"Dirac-Fermions), {z_max0} (Maxwell: transversal {z_tv[0]}, longitudinal ohne Nullmode, zwei Geister); mit "
      f"Int <T> = -4 a: a = {A_FREI['skalar']} (reeller Skalar), {A_FREI['weyl']} (Weyl), {A_FREI['vektor']} (Vektor). "
      f"Numerische Gegenprobe über die Wärmespur (Konstante der Entwicklung bei kleinem t, alle vier Spektren): max. "
      f"Abweichung {mp.nstr(c0_abw, 2)}. Lokale Gegenprobe für den Skalar: b_4 = (1/180)(Riem^2 - Ric^2) = c W^2 - a E_4 "
      f"gibt a = {lsg[0][aa] if lsg else None}, c = {lsg[0][cc] if lsg else None} (dasselbe c wie aus C_T in Teil 1, dasselbe a wie aus zeta(0); auf S^4 "
      f"{lok_s4}). Negativkontrolle: der minimal gekoppelte Skalar hat {z_min[0]}, kein konformes a",
      A_FREI == {"skalar": Fr(1, 360), "weyl": Fr(11, 720), "vektor": Fr(31, 180)} and c0_abw < mp.mpf(10)**-6
      and bool(lsg) and lok_rest == 0 and lsg[0][aa] == R_(1, 360) and lsg[0][cc] == R_(1, 120) and lok_s4 == z_sk[0]
      and z_min[0] != z_sk[0])


# --- 2.2 Die Lesarten des nativen Inhalts
def ac(ns=0, nh=0, nw=0, nv=0):
    """ns konforme Leiterskalare, nh minimal gekoppelte Skalare (Higgs), nw Weyl-Fermionen, nv Eichbosonen."""
    return ((ns + nh)*A_FREI["skalar"] + nw*A_FREI["weyl"] + nv*A_FREI["vektor"],
            (ns + nh)*C_FREI["skalar"] + nw*C_FREI["weyl"] + nv*C_FREI["vektor"])


LESARTEN = [("(i) eine native Leiter (reeller Skalar, z83)", dict(ns=1)),
            ("(i') Quarterpaar als C-Paar (komplexer Skalar)", dict(ns=2)),
            ("(ii) 3 x 16 Weyl-Fermionen", dict(nw=48)),
            ("(iii) (ii) + 12 Eichbosonen (SU(3) x SU(2) x U(1), wie z27)", dict(nw=48, nv=12)),
            ("(iii-H) (iii) + Higgs-Dublett (4 reelle Skalare)", dict(nh=4, nw=48, nv=12)),
            ("(iii-45) (ii) + 45 Eichbosonen (Spin(10))", dict(nw=48, nv=45)),
            ("(iii-248) (ii) + 248 Eichbosonen (E8)", dict(nw=48, nv=248))]
AC = {name: ac(**kw) for name, kw in LESARTEN}
HAUPT = [LESARTEN[0][0], LESARTEN[2][0], LESARTEN[3][0]]
EH_ZAHL = {name: int(kw.get("ns", 0)*EH_EINH["konformer Skalar"] + kw.get("nh", 0)*EH_EINH["minimaler Skalar"]
                     + kw.get("nw", 0)*EH_EINH["Weyl"] + kw.get("nv", 0)*EH_EINH["Vektor"])
           for name, kw in LESARTEN}                                 # Wärmekern: N_0 (1 - 6 xi) + N_Weyl - 4 N_1; Leiter xi = 1/6, Higgs minimal
ZAHLEN["lesarten"] = {name: {"a": str(AC[name][0]), "c": str(AC[name][1]), "c - a": str(AC[name][1] - AC[name][0]),
                             "Einstein-Hilbert-Zählung (Wärmekern)": EH_ZAHL[name]} for name, _ in LESARTEN}
say("   Lesart: (a, c), c - a, Einstein-Hilbert-Zählung im Wärmekernregler")
for name, _ in LESARTEN:
    say(f"      {name}: a = {AC[name][0]}, c = {AC[name][1]}, c - a = {AC[name][1] - AC[name][0]}, EH-Zählung {EH_ZAHL[name]}")
check(f"[E] Anomaliekoeffizienten je Lesart (Summen der freien Werte): (i) (a, c) = ({AC[HAUPT[0]][0]}, {AC[HAUPT[0]][1]}), (ii) "
      f"({AC[HAUPT[1]][0]}, {AC[HAUPT[1]][1]}), (iii) ({AC[HAUPT[2]][0]}, {AC[HAUPT[2]][1]}); c - a wechselt zwischen (ii) und "
      f"(iii) das Vorzeichen ({AC[HAUPT[1]][1] - AC[HAUPT[1]][0]}, {AC[HAUPT[2]][1] - AC[HAUPT[2]][0]}). Welche Lesart nativ ist, entscheidet die Rechnung nicht: (i) ist der eigene "
      f"Einteilcheninhalt einer Naht in der Lesart von z83 [C], (ii) und (iii) sind die gelösten Blöcke als 4D-Felder "
      f"gelesen [L]",
      AC[HAUPT[0]] == (Fr(1, 360), Fr(1, 120)) and AC[HAUPT[1]] == (Fr(11, 15), Fr(6, 5)) and AC[HAUPT[2]] == (Fr(14, 5), Fr(12, 5)))

# --- 2.3 Universeller Term und Einstein-Hilbert-Term
sig27 = Z27["zahlen"]["sigma_je_feld"]["skalar"]
sig_werte = [sig27[k] for k in sorted(sig27)]
r_alle = Z27["zahlen"]["alle_regeln"]["R_gesamt"]
alp_alle = Z27["zahlen"]["alle_regeln"]["a_lP_bei_R1"]
ZAHLEN["einstein_hilbert"] = {"Wärmekern-Einheiten (minimaler Skalar = 1)": {k: str(v) for k, v in EH_EINH.items()},
                              "z27 sigma = a^2/(4 G_ind), Skalar je Kopplungsregel": {k: round(sig27[k], 6) for k in sorted(sig27)},
                              "z27 R_gesamt über alle Regeln": [round(x, 4) for x in r_alle],
                              "z27 a/l_P bei R = 1": [round(x, 3) for x in alp_alle]}
check(f"[E/L] Universell ist allein der Logarithmus: Gamma_ind = -(L/(16 pi^2)) Int (c W^2 - a E_4), L = ln(Abschneide-/"
      f"IR-Skala), also d alpha_W/d ln mu = c/(16 pi^2) für S = alpha_W Int W^2. Der Einstein-Hilbert-Term hängt am "
      f"Regler: im Wärmekern (Res_(s=1) derselben S^4-Spektren, Einheiten des minimalen Skalars) zählen konformer Skalar, "
      f"Weyl, Vektor {EH_EINH['konformer Skalar']}, {EH_EINH['Weyl']}, {EH_EINH['Vektor']}: Lesart (i) gibt 0, (ii) "
      f"{EH_ZAHL[HAUPT[1]]}, (iii) {EH_ZAHL[HAUPT[2]]} (48 Weyl gegen 12 Vektoren heben sich in diesem Regler auf; nur eine "
      f"Regleraussage). Auf dem H6-Gitter (z27, nur gelesen) läuft der Skalar-Koeffizient sigma je Kopplungsregel von "
      f"{min(sig_werte):.3f} bis {max(sig_werte):.3f} (Vorzeichenwechsel), R_gesamt von {r_alle[0]:.2f} bis {r_alle[1]:.2f}, "
      f"a/l_P von {alp_alle[0]:.1f} bis {alp_alle[1]:.1f}",
      EH_EINH["konformer Skalar"] == 0 and EH_EINH["Weyl"] == 1 and EH_EINH["Vektor"] == -4 and EH_ZAHL[HAUPT[0]] == 0
      and EH_ZAHL[HAUPT[1]] == 48 and EH_ZAHL[HAUPT[2]] == 0 and min(sig_werte) < 0 < max(sig_werte)
      and max(sig_werte)/sig27["K"] > 10 and r_alle[0] < 0 < r_alle[1])
befund("E/L: Der induzierte Einstein-Hilbert-Koeffizient ist nicht festgelegt (Regler und Metrikkopplung, z27; im "
       "Wärmekern für (i) und (iii) sogar null). Fest ist nur c vor dem Weyl^2-Logarithmus", True)

# =================================================================== Teil 3
say("\nTeil 3: Vom Weyl-Koeffizienten zu G Lambda")

# --- 3.1 Punktweise Identität auf Einstein-Räumen (Schwarzschild-(A)dS, euklidisch)
tau_, r_s, th_s, ph_s = sp.symbols("tau r theta phi", real=True)
m_s, Lam = sp.symbols("m Lambda", real=True)
f_s = 1 - 2*m_s/r_s - Lam*r_s**2/3
koord = [tau_, r_s, th_s, ph_s]
gdn = sp.diag(f_s, 1/f_s, r_s**2, r_s**2*sp.sin(th_s)**2)
gup = gdn.inv()
Gam_ = [[[sum(gup[i, l]*(sp.diff(gdn[l, j], koord[k]) + sp.diff(gdn[l, k], koord[j]) - sp.diff(gdn[j, k], koord[l]))
              for l in range(4))/2 for k in range(4)] for j in range(4)] for i in range(4)]


def riem(i, j, k, l):
    """R^i_jkl."""
    out = sp.diff(Gam_[i][j][l], koord[k]) - sp.diff(Gam_[i][j][k], koord[l])
    out += sum(Gam_[i][k][p]*Gam_[p][j][l] - Gam_[i][l][p]*Gam_[p][j][k] for p in range(4))
    return out


Rudd = [[[[sp.simplify(riem(i, j, k, l)) for l in range(4)] for k in range(4)] for j in range(4)] for i in range(4)]
Ricd = sp.Matrix(4, 4, lambda j, l: sp.simplify(sum(Rudd[i][j][i][l] for i in range(4))))
einstein_rest = sp.simplify(Ricd - Lam*gdn)
Rdddd = [[[[sum(gdn[i, p]*Rudd[p][j][k][l] for p in range(4)) for l in range(4)] for k in range(4)] for j in range(4)] for i in range(4)]
Ruuuu_d = lambda i, j, k, l: Rudd[i][j][k][l]*gup[j, j]*gup[k, k]*gup[l, l]      # diagonale Metrik
Riem2_w = sp.simplify(sum(Rdddd[i][j][k][l]*Ruuuu_d(i, j, k, l) for i in range(4) for j in range(4) for k in range(4) for l in range(4)))
Ric2_w, Rsk_w = sp.simplify(16*Lam**2/4), 4*Lam
W2_w = sp.simplify(Riem2_w - 2*Ric2_w + Rsk_w**2/3)
E4_w = sp.simplify(Riem2_w - 4*Ric2_w + Rsk_w**2)
ident = sp.simplify(E4_w - W2_w - sp.Rational(8, 3)*Lam**2)
ZAHLEN["einstein_identitaet"] = {"Ric - Lambda g": str(einstein_rest.norm()), "W^2": str(W2_w), "E_4 - W^2": str(sp.simplify(E4_w - W2_w))}
check(f"[E] Punktweise auf Einstein-Räumen (sympy, Schwarzschild-(A)dS mit f = 1 - 2m/r - Lambda r^2/3: Ric = Lambda g, "
      f"Rest {einstein_rest.norm()}): W^2 = {W2_w} und E_4 - W^2 = {sp.simplify(E4_w - W2_w)} = 8 Lambda^2/3. Also "
      f"Int W^2 = 32 pi^2 chi - (8 Lambda^2/3) Vol und S_EH = -(Lambda/(8 pi G)) Vol = (3/(64 pi G Lambda)) (Int W^2 - "
      f"32 pi^2 chi): auf Einstein-Lösungen ist die Einstein-Wirkung die Weyl-Wirkung mit alpha_W = 3/(64 pi G Lambda) "
      f"plus Euler-Term (Maldacena; Anderson). Für AdS (Lambda = -3/l^2) alpha_W = -l^2/(64 pi G), für dS +l^2/(64 pi G)",
      einstein_rest == sp.zeros(4, 4) and sp.simplify(W2_w - 48*m_s**2/r_s**6) == 0 and ident == 0)

# --- 3.2 Bekannte Auswertungen: AdS4, S^4, Hawking-Page
rho, Ee, l_s, G_s, rp, Rc = sp.symbols("rho E l G r_p R_c", positive=True)
vol_h4 = sp.expand((2*PI**2*sp.integrate(sp.sinh(rho)**3, (rho, 0, rho))).rewrite(sp.exp).subs(sp.exp(rho), Ee))
vol_poly = sp.Poly(sp.expand(vol_h4*Ee**3), Ee)
vol_ren = vol_poly.coeff_monomial(Ee**3)                             # Konstante in e^rho; divergent sind nur e^(3 rho), e^rho
log_frei = (not vol_h4.has(rho)) and sorted(m[0] - 3 for m in vol_poly.monoms()) == [-3, -1, 0, 1, 3]
W2_ads = 32*PI**2*1 - 24*vol_ren                                     # l = 1, chi(B^4) = 1
S_ads = sp.simplify(3/(8*PI*G_s*l_s**2)*l_s**4*vol_ren)
vol_s4 = sp.Rational(8, 3)*PI**2
W2_s4 = 32*PI**2*2 - 24*vol_s4
S_s4 = sp.simplify(-3/(8*PI*G_s*l_s**2)*l_s**4*vol_s4)
# Hawking-Page: direkte regularisierte Wirkung (Hintergrundabzug), dann gegen top - k (l^2/G) Int W^2
f_bh = lambda r__, m__: 1 + r__**2/l_s**2 - 2*m__/r__
m_bh = sp.solve(f_bh(rp, m_s), m_s)[0]
beta_bh = sp.simplify(4*PI/sp.diff(f_bh(r_s, m_s), r_s).subs({r_s: rp, m_s: m_bh}))
V_bh = beta_bh*4*PI*(Rc**3 - rp**3)/3
V_ads = beta_bh*sp.sqrt(f_bh(Rc, m_bh)/f_bh(Rc, 0))*4*PI*Rc**3/3
I_hp = sp.simplify(sp.limit(3/(8*PI*G_s*l_s**2)*(V_bh - V_ads), Rc, sp.oo))
intW2_bh = sp.simplify(beta_bh*4*PI*sp.integrate(r_s**2*48*m_bh**2/r_s**6, (r_s, rp, sp.oo)))
k_s = sp.Symbol("k")
k_lsg = sp.solve(sp.simplify(PI*l_s**2/(2*G_s)*2 - k_s*l_s**2/G_s*intW2_bh - I_hp), k_s)
hp_lit = PI*rp**2*(l_s**2 - rp**2)/(G_s*(l_s**2 + 3*rp**2))
ZAHLEN["maldacena"] = {"Vol_ren(H^4), l = 1": str(vol_ren), "Int W^2 (AdS4)": str(sp.simplify(W2_ads)), "S_ren(AdS4)": str(S_ads),
                       "Int W^2 (S^4)": str(sp.simplify(W2_s4)), "S(S^4)": str(S_s4), "I_Hawking-Page": str(sp.factor(I_hp)),
                       "Int W^2 (AdS-Schwarzschild)": str(sp.factor(intW2_bh)), "k aus Hawking-Page": [str(x) for x in k_lsg],
                       "alpha_W": "3/(64 pi G Lambda)", "l^2/G": "64 pi |alpha_W|"}
check(f"[E] Auswertungen: renormiertes Volumen von H^4 (Konstante der Entwicklung in e^rho; nur die Potenzen e^(+-3 rho), e^(+-rho), 1, "
      f"kein Logarithmus: {log_frei}) = {vol_ren}, damit Int W^2 = 32 pi^2 - 24 Vol_ren = {sp.simplify(W2_ads)} (AdS4 ist konform flach) und "
      f"S_ren = {S_ads}; S^4: Int W^2 = {sp.simplify(W2_s4)}, S = {S_s4} = -S_dS. Mit W != 0: die direkt regularisierte "
      f"Wirkung des AdS-Schwarzschild-Lochs (Hintergrundabzug, Limes) ist {sp.factor(I_hp)} (Hawking-Page), und "
      f"I = (pi l^2/(2G)) chi - k (l^2/G) Int W^2 mit chi = 2 und Int W^2 = {sp.factor(intW2_bh)} verlangt genau k = "
      f"{k_lsg} = 1/(64 pi), für jedes r_+. Die reine Zahl ist also l^2/G = 64 pi |alpha_W|",
      vol_ren == sp.Rational(4, 3)*PI**2 and log_frei and sp.simplify(W2_ads) == 0 and sp.simplify(S_ads - PI*l_s**2/(2*G_s)) == 0
      and sp.simplify(W2_s4) == 0 and sp.simplify(S_s4 + PI*l_s**2/G_s) == 0 and sp.simplify(I_hp - hp_lit) == 0
      and len(k_lsg) == 1 and sp.simplify(k_lsg[0] - 1/(64*PI)) == 0)

# --- 3.3 Die reinen Zahlen je Lesart
Ls, cs = sp.symbols("L c", positive=True)
alpha_ind = -cs*Ls/(16*PI**2)                                        # rein induziert: nackter Weyl-Koeffizient null (P_ind)
GLam = sp.simplify(3/(64*PI*alpha_ind))                              # aus alpha_W = 3/(64 pi G Lambda)
l2G = sp.simplify(-3/GLam)                                           # Lambda = -3/l^2
RV = sp.Symbol("R_Vak", positive=True)
G_nativ = sp.simplify((RV/2)**2/l2G)                                 # l = R_Vak/2 (z98), Lambda = -12/R_Vak^2 (z118)
lam_konsistent = sp.simplify(-3/(RV/2)**2 + 12/RV**2) == 0
TAB = {}
for name, _ in LESARTEN:
    c_ = AC[name][1]
    TAB[name] = {"l^2/G je e-Faltung (L = 1)": float(4*c_/np.pi), "l^2/G bei L = 2 pi": float(8*c_),
                 "G Lambda mal L": float(-3*np.pi/(4*c_)), "G/R_Vak^2 mal L": float(np.pi/(16*c_)),
                 "l^2/G exakt": f"{4*c_} L/pi", "G Lambda exakt": f"-{Fr(3, 4)/c_} pi/L"}
ZAHLEN["reine_zahlen"] = {"alpha_W induziert": str(alpha_ind), "G Lambda": str(GLam), "l^2/G": str(l2G),
                          "G nativ (l = R_Vak/2)": str(G_nativ), "je Lesart": {k: {kk: (round(vv, 6) if isinstance(vv, float) else vv)
                                                                                 for kk, vv in v.items()} for k, v in TAB.items()}}
say("   Lesart: l^2/G = (4 c/pi) L, G Lambda = -3 pi/(4 c L)")
for name, _ in LESARTEN:
    say(f"      {name}: l^2/G = {TAB[name]['l^2/G exakt']} = {TAB[name]['l^2/G je e-Faltung (L = 1)']:.5f} L, "
        f"G Lambda = {TAB[name]['G Lambda exakt']} = {TAB[name]['G Lambda mal L']:.4f}/L")
check(f"[C] Reine Zahlen unter P_ind (nackter Weyl-Koeffizient null an der Nahtskala) und P_Neumann (Einschränkung auf "
      f"Einstein-Lösungen): alpha_W = {alpha_ind}, also G Lambda = {GLam} und l^2/G = {l2G}; mit l = R_Vak/2 und "
      f"Lambda = -12/R_Vak^2 (konsistent: {lam_konsistent}) G = {G_nativ}. Je e-Faltung: (i) l^2/G = "
      f"{TAB[HAUPT[0]]['l^2/G exakt']} = {TAB[HAUPT[0]]['l^2/G je e-Faltung (L = 1)']:.5f} L, (ii) "
      f"{TAB[HAUPT[1]]['l^2/G exakt']} = {TAB[HAUPT[1]]['l^2/G je e-Faltung (L = 1)']:.4f} L, (iii) "
      f"{TAB[HAUPT[2]]['l^2/G exakt']} = {TAB[HAUPT[2]]['l^2/G je e-Faltung (L = 1)']:.4f} L; G Lambda = "
      f"{TAB[HAUPT[0]]['G Lambda mal L']:.2f}/L, {TAB[HAUPT[1]]['G Lambda mal L']:.4f}/L, "
      f"{TAB[HAUPT[2]]['G Lambda mal L']:.4f}/L. L ist nicht nativ: eine Zahl entsteht erst mit einer Abschneideskala",
      sp.simplify(GLam + 3*PI/(4*cs*Ls)) == 0 and sp.simplify(l2G - 4*cs*Ls/PI) == 0 and lam_konsistent
      and abs(TAB[HAUPT[1]]["l^2/G je e-Faltung (L = 1)"] - 24/(5*np.pi)) < 1e-12
      and abs(TAB[HAUPT[2]]["G Lambda mal L"] + 5*np.pi/16) < 1e-12)
befund("E/C: Das Vorzeichen ist nicht frei. Mit c > 0 (Unitarität) und L > 0 ist der induzierte Weyl-Koeffizient negativ "
       "und sign(G Lambda) = sign(alpha_W) = -1: induzierte Gravitation aus unitärer Materie gibt bei G > 0 genau AdS. Das "
       "passt zum nativen Vorzeichen (z98, z104, z118) und widerspricht der Beobachtung (Lambda > 0); dS bräuchte einen "
       "positiven nackten Weyl-Koeffizienten (nicht induziert) oder c < 0", True)

# --- 3.4 Vergleich mit z118 (64 je Index, Gewicht 2 c_3)
c1_118 = sp.nsimplify(complex(Z118["zahlen"]["mm"]["c1 = L1/(R det e)"][0]).real)
c0_118 = sp.nsimplify(complex(Z118["zahlen"]["mm"]["c0/c1 (Euler)"]).real)*c1_118
c2_118 = sp.Rational(Z118["zahlen"]["mm"]["c2/c1"])*c1_118
lam_s = sp.Symbol("lambda", positive=True)
E4s = sp.Symbol("E4")
mm_dichte = c0_118*E4s + c1_118*lam_s*(4*(-3*lam_s)) + c2_118*lam_s**2      # Einstein-Raum: R = 4 Lambda, Lambda = -3 lambda
weyl_dichte = c0_118*(E4s - sp.Rational(8, 3)*(3*lam_s)**2)
mm_rest = sp.simplify(mm_dichte - weyl_dichte)
kap = sp.Symbol("kappa", positive=True)
l2G_118 = 64*PI*c0_118*kap                                          # b = kappa c_0 = l^2/(64 pi G) je Index
c3 = 1/(8*PI)
vl = sp.Symbol("vl", positive=True)                                  # v_geo l
kap_hand = sp.solve(sp.Eq(kap*c1_118/1, 2*c3/4*vl**2), kap)[0]       # kappa c_1 lambda = (2 c_3/4) v^2, lambda = 1/l^2
l2G_hand = sp.simplify(l2G_118.subs(kap, kap_hand))
e2a = mp.exp(2*ALPHA_INV)
ziel = 32*mp.pi**3*e2a                                               # l^2/G aus COSMO.LAM.01 (z96/z98), nur Vergleichsziel
kap_ziel = ziel/(64*mp.pi*int(c0_118))
kap_ind = {name: float(AC[name][1])/(16*np.pi**2*int(c0_118)) for name in HAUPT}
ZAHLEN["z118"] = {"c_0, c_1, c_2": [int(c0_118), int(c1_118), int(c2_118)], "MM-Dichte - c_0 W^2 auf Einstein-Räumen": str(mm_rest),
                  "l^2/G je Index": "1024 pi kappa", "kappa mit Gewicht 2 c_3": str(kap_hand), "l^2/G mit 2 c_3": str(l2G_hand),
                  "log10 kappa für das Ledgerziel": float(mp.log10(kap_ziel)), "kappa induziert je e-Faltung": {k: round(v, 8) for k, v in kap_ind.items()}}
check(f"[E] Brücke zu z118: mit dessen Zahlen (c_1 = {c1_118}, c_0 = {c0_118}, c_2 = {c2_118}) ist die "
      f"MacDowell-Mansouri-Dichte auf Einstein-Räumen genau c_0 W^2 (Rest {mm_rest}): die Uhrspur ist dort die "
      f"Weyl-Wirkung mit alpha_W = 16 kappa je Index, also l^2/G = 1024 pi kappa je Index (z98: b = l^2/(64 pi G)). Das "
      f"gesetzte Gewicht 2 c_3 heißt kappa = {kap_hand} und l^2/G = {l2G_hand} = (v_geo l)^2: keine reine Zahl, solange "
      f"v_geo l nicht bestimmt ist; das Ledgerziel 32 pi^3 e^(2/alpha) wäre kappa = 10^{float(mp.log10(kap_ziel)):.2f}. "
      f"Induziert wäre kappa = c L/(256 pi^2): {kap_ind[HAUPT[0]]:.2e} L, {kap_ind[HAUPT[1]]:.2e} L, "
      f"{kap_ind[HAUPT[2]]:.2e} L; ganzzahliges kappa = 1 verlangt c L = 256 pi^2 = {256*np.pi**2:.1f}",
      mm_rest == 0 and c0_118 == 16 and sp.simplify(l2G_hand - vl**2) == 0 and float(mp.log10(kap_ziel)) > 100
      and abs(kap_ind[HAUPT[1]] - 1.2/(256*np.pi**2)) < 1e-15)

# =================================================================== Teil 4
say("\nTeil 4: Ist die Kopplung eine Nahtzahl?")
QQ = float(np.exp(-2*np.pi))


# --- 4.1 Uhrspurgewicht auf N Nähten
def A_rechteck(n, q):
    """A = S^+ (1 - q S^4) auf C^n -> C^(n-1) (wie z118)."""
    S = np.zeros((n + 9, n + 9))
    for k in range(n + 8):
        S[k + 1, k] = 1
    return (S.T @ (np.eye(n + 9) - q*np.linalg.matrix_power(S, 4)))[:n - 1, :n]


def kern_daten(A, n_naht, n):
    """(dim Kern, dim Kokern, größtes Gewicht eines Kernvektors am Abschneiderand der Marken)."""
    u, sv, vh = np.linalg.svd(A)
    rang = int(np.sum(sv > 1e-10))
    kv = vh[rang:]
    rand = 0.0
    if len(kv):
        rand = max(float(np.sum(np.abs(v.reshape(n_naht, n)[:, -4:])**2)) for v in kv)
    return A.shape[1] - rang, A.shape[0] - rang, rand


n_m = 25
Aq = A_rechteck(n_m, QQ)
ind_tab, rand_tab = {}, {}
for Nn_ in range(1, 7):
    AN = block_diag(*[Aq]*Nn_)
    Us = np.linalg.qr(rng.normal(size=(Nn_, Nn_)) + 1j*rng.normal(size=(Nn_, Nn_)))[0]
    Ut = np.linalg.qr(rng.normal(size=(Nn_, Nn_)) + 1j*rng.normal(size=(Nn_, Nn_)))[0]
    K_ = np.zeros(AN.shape, complex)                                 # kompakte Kopplung der Nähte nahe der Marke k = 0
    for i in range(Nn_):
        for j in range(Nn_):
            K_[i*(n_m - 1):i*(n_m - 1) + 5, j*n_m:j*n_m + 5] = 0.3*rng.normal(size=(5, 5))
    AN2 = np.kron(Ut, np.eye(n_m - 1)) @ AN @ np.kron(Us, np.eye(n_m)).conj().T + K_
    k1, c1_, r1 = kern_daten(AN, Nn_, n_m)
    k2, c2_, r2 = kern_daten(AN2, Nn_, n_m)
    ind_tab[Nn_] = (k1 - c1_, k2 - c2_)
    rand_tab[Nn_] = max(r1, r2)
k_neg, c_neg, r_neg = kern_daten(A_rechteck(n_m, 2.0), 1, n_m)       # |q| > 1: Symbolwindung +3, keine normierbare Nullmode
n_p = 12
Ap = A_rechteck(n_p, QQ)
Dsharp = np.block([[np.kron(Ap, np.eye(n_p)), -np.kron(np.eye(n_p - 1), Ap.T)],
                   [np.kron(np.eye(n_p), Ap), np.kron(Ap.T, np.eye(n_p - 1))]])
u_, sv_, vh_ = np.linalg.svd(Dsharp)
rang_p = int(np.sum(sv_ > 1e-10))
zm1 = np.linalg.svd(Ap)[2][-1]
ueberl_p = float(abs(np.vdot(np.kron(zm1, zm1), vh_[-1][:n_p**2]))**2)
ind_prod = (Dsharp.shape[1] - rang_p) - (Dsharp.shape[0] - rang_p)
ZAHLEN["uhrspur_N"] = {"Index direkte Summe (frei, gekoppelt)": {str(k): list(v) for k, v in ind_tab.items()},
                       "Gewicht der Nullmoden am Abschneiderand < 1e-8": bool(max(rand_tab.values()) < 1e-8),
                       "Kontrolle q = 2: (Kern, Kokern, Randgewicht)": [k_neg, c_neg, round(r_neg, 6)],
                       "Index des graduierten Produkts zweier Nähte": ind_prod, "Überlapp seiner Nullmode mit zm x zm": round(ueberl_p, 10),
                       "Uhrspur": "64 N (direkte Summe) gegen 64 (graduiertes Produkt)"}
check(f"[E] Uhrspurgewicht auf N Nähten: im Tensorprodukt der Fockräume (gemeinsame CAR, z106) ist der Einteilchenoperator "
      f"die direkte Summe der A_q; ihr Index ist N = 1..6: {[v[0] for v in ind_tab.values()]}, unverändert unter "
      f"nahtmischenden Unitären und kompakter Kopplung an den Marken: {[v[1] for v in ind_tab.values()]} (echte Nullmoden: "
      f"Gewicht am Abschneiderand < 1e-8: {max(rand_tab.values()) < 1e-8}). Mit dem inneren Faktor von z118 ist die Uhrspur 64 N: "
      f"additiv. Kontrollen: q = 2 gibt nur eine Scheinnullmode am Abschneiderand (Gewicht {r_neg:.3f}); das graduierte "
      f"Produkt zweier Nähte (Dirac-Operator des Produktraums) hat Index {ind_prod} = 1 x 1 (Nullmode = zm x zm, Überlapp "
      f"{ueberl_p:.6f}): multiplikativ, zählt nicht. Additivität ist die Lesart 'direkte Summe unterscheidbarer Nähte' (P_ord)",
      all(v == (k, k) for k, v in ind_tab.items()) and max(rand_tab.values()) < 1e-8 and r_neg > 0.5 and ind_prod == 1
      and abs(ueberl_p - 1) < 1e-8)

# --- 4.2 Entropie des N-Naht-Zustands
S_quarter = mp.fsum(mp.log(1 + mp.exp(-k*mp.pi)) + k*mp.pi/(mp.exp(k*mp.pi) + 1) for k in range(1, 201, 2))
S_naht = 2*S_quarter


def swap_idx(L, x, y):
    idx = np.arange(2**L)
    bx, by = (idx >> x) & 1, (idx >> y) & 1
    return np.where(bx != by, idx ^ ((1 << x) | (1 << y)), idx)


def F_anw(V, L, u=1.0):
    """Backsteinschritt von z121: Gatter R(u) = (1 + i u P)/(1 + i u) auf (0,1), (2,3), ..., dann (1,2), ..., (L-1,0)."""
    a, b = 1/(1 + 1j*u), 1j*u/(1 + 1j*u)
    for x in list(range(0, L, 2)) + list(range(1, L, 2)):
        V = a*V + b*V[swap_idx(L, x, (x + 1) % L)]
    return V


def entropie(rho_):
    w = np.linalg.eigvalsh((rho_ + rho_.conj().T)/2)
    w = w[w > 1e-14]
    return float(-(w*np.log(w)).sum())


def teilspur_links(rho_, L, nl):
    """Reduktion auf die ersten nl Plätze (niedrige Bits)."""
    t = rho_.reshape(2**(L - nl), 2**nl, 2**(L - nl), 2**nl)
    return np.einsum("aibj->ij", t)


p_occ = 1/(1 + np.exp(-np.pi))                                      # Besetzung der tiefsten Quartermode bei beta = 4 pi
p_kon = 1/(1 + np.exp(-3*np.pi))                                    # Kontrolle: jede zweite Naht anders besetzt
h2 = lambda p: float(-(p*np.log(p) + (1 - p)*np.log(1 - p)))
ent_tab = {}
for L in (2, 4, 6, 8, 10):
    Fm = F_anw(np.eye(2**L, dtype=complex), L)
    unit = float(np.abs(Fm.conj().T @ Fm - np.eye(2**L)).max())
    zeile = {"unitär": unit}
    for lab, ps in (("gleich", [p_occ]*L), ("gestaffelt", [p_occ, p_kon]*(L//2))):
        rhoN = np.diag([ps[0], 1 - ps[0]]).astype(complex)
        for p in ps[1:]:
            rhoN = np.kron(np.diag([p, 1 - p]).astype(complex), rhoN)   # Platz j = Bit j
        rr = rhoN
        for _ in range(L if lab == "gleich" else 1):
            rr = Fm @ rr @ Fm.conj().T
        S_prod = sum(h2(p) for p in ps)
        zeile[lab] = {"|F rho F^+ - rho|": float(np.abs(rr - rhoN).max()), "S_ges/Summe": entropie(rr)/S_prod,
                      "S(Naht 0)/s": entropie(teilspur_links(rr, L, 1))/h2(ps[0]),
                      "S(N/2)/Summe": entropie(teilspur_links(rr, L, L//2))/sum(h2(p) for p in ps[:L//2])}
    ent_tab[L] = zeile
ZAHLEN["entropie_N"] = {"je Quarter": mp.nstr(S_quarter, 20), "je Naht": mp.nstr(S_naht, 20), "ln 2/2": mp.nstr(mp.log(2)/2, 20),
                        "Kette (N Schritte), gleiche Nähte: S_ges/Summe, S(N/2)/Summe": {
                            str(k): [round(v["gleich"]["S_ges/Summe"], 8), round(v["gleich"]["S(N/2)/Summe"], 8)] for k, v in ent_tab.items()},
                        "Kontrolle gestaffelt: S_ges/Summe, S(Naht 0)/s, S(N/2)/Summe": {
                            str(k): [round(v["gestaffelt"]["S_ges/Summe"], 8), round(v["gestaffelt"]["S(Naht 0)/s"], 5),
                                     round(v["gestaffelt"]["S(N/2)/Summe"], 5)] for k, v in ent_tab.items()}}
gleich_inv = max(v["gleich"]["|F rho F^+ - rho|"] for v in ent_tab.values())
check(f"[E] Entropie: je Quarter Summe über die Niveaus pi k (k ungerade, beta = 4 pi) = ln 2/4, je Naht ln 2/2 "
      f"(Abweichung < 1e-30; z68), im Produktzustand von N Nähten N ln 2/2. Die Kette von z121 (N Backsteinschritte, "
      f"N = 2..10, unitär auf den Originalen) lässt den Produktzustand gleicher Nähte exakt invariant (|F rho F^+ - rho| < "
      f"1e-12: {gleich_inv < 1e-12}; die Gatter sind a + b P, der Zustand ist tauschsymmetrisch): Gesamt- und alle "
      f"Blockentropien bleiben additiv, S(N/2)/((N/2) s_1) = "
      f"{[round(v['gleich']['S(N/2)/Summe'], 8) for v in ent_tab.values()]}. Kontrolle mit gestaffelten Nähten "
      f"(p = 1/(1 + e^-pi), 1/(1 + e^-3pi), ein Backsteinschritt): die Gesamtentropie bleibt die Summe "
      f"({[round(v['gestaffelt']['S_ges/Summe'], 8) for v in ent_tab.values()]}), die einer einzelnen Naht nicht "
      f"({[round(v['gestaffelt']['S(Naht 0)/s'], 3) for v in ent_tab.values()]} s): additiv ist dann nur die Summe",
      abs(S_quarter - mp.log(2)/4) < mp.mpf(10)**-30 and gleich_inv < 1e-12
      and all(v["unitär"] < 1e-12 and abs(v["gleich"]["S_ges/Summe"] - 1) < 1e-9 and abs(v["gleich"]["S(N/2)/Summe"] - 1) < 1e-9
              and abs(v["gestaffelt"]["S_ges/Summe"] - 1) < 1e-9 and abs(v["gestaffelt"]["S(Naht 0)/s"] - 1) > 1e-3
              for v in ent_tab.values()))

# --- 4.3 C_T und c auf N Nähten
t01n = np.array([float(x) for x in t01])
G2n = np.array(L1["G2"].tolist(), dtype=float)
tsw = np.array(tausch.tolist(), dtype=float)


def ct_arten(N):
    """Zweiteilchenraum von N unterscheidbaren Leitern (Artenpaare i <= j): Gram der T_ii, T_ij und von T_ges."""
    paare = [(i, j) for i in range(N) for j in range(i, N)]
    gram_ = np.zeros((len(paare), len(paare)))
    for a_, (i, j) in enumerate(paare):
        # bosonischer Fockraum: <(i,j,p)|(i,j,q)> = <p|q> + delta_ij <p|q vertauscht>
        gram_[a_, a_] = t01n @ G2n @ t01n + (t01n @ G2n @ tsw @ t01n if i == j else 0.0)
    diag = np.array([1.0 if i == j else 0.0 for (i, j) in paare])
    return len(paare), float(diag @ gram_ @ diag), paare, gram_


ct_N = {N: ct_arten(N) for N in (1, 2, 3, 4, 6)}
ct_verh = {N: ct_N[N][1]/ct_N[1][1] for N in ct_N}
anz_strom = {N: ct_N[N][0] for N in ct_N}
# Artentausch (das Gatter von z121 ist a + b P): T_ges invariant, T_11 - T_22 nicht
paare3, gram3 = ct_N[3][2], ct_N[3][3]
perm = {0: 1, 1: 0, 2: 2}
Pm = np.zeros((len(paare3), len(paare3)))
for a_, (i, j) in enumerate(paare3):
    Pm[paare3.index(tuple(sorted((perm[i], perm[j])))), a_] = 1
v_ges = np.array([1.0 if i == j else 0.0 for (i, j) in paare3])
v_diff = np.array([1.0 if (i, j) == (0, 0) else (-1.0 if (i, j) == (1, 1) else 0.0) for (i, j) in paare3])
inv_ges, inv_diff = float(np.abs(Pm @ v_ges - v_ges).max()), float(np.abs(Pm @ v_diff - v_diff).max())
# Uhr und Kette: einheitliche Uhr vertauscht mit F, eine Uhr auf einer einzelnen Naht nicht
Lc = 6
Fc = F_anw(np.eye(2**Lc, dtype=complex), Lc)
r1 = np.diag([1.0, 1j])
Runi, Rein = r1, r1
for _ in range(Lc - 1):
    Runi, Rein = np.kron(Runi, r1), np.kron(Rein, np.eye(2))
komm_uni, komm_ein = float(np.abs(Fc @ Runi - Runi @ Fc).max()), float(np.abs(Fc @ Rein - Rein @ Fc).max())


# Verschränkung zwischen den Nähten im Vakuum der Kette (stetig aus dem Heisenberg-Grundzustand verfolgt)
def sektor(L, nup):
    return np.array([s_ for s_ in range(2**L) if bin(s_).count("1") == nup])


def verschr(psi, L):
    s_ = np.linalg.svd(psi.reshape(2**(L - L//2), 2**(L//2)), compute_uv=False)**2
    s_ = s_[s_ > 1e-15]
    return float(-(s_*np.log(s_)).sum())


kette = {}
for L in (6, 8, 10, 12):
    Dd_, b_ = 2**L, sektor(L, L//2)
    Eb = np.zeros((Dd_, len(b_)))
    Eb[b_, np.arange(len(b_))] = 1
    Hh = sum(0.5*Eb[swap_idx(L, x, (x + 1) % L)] - 0.25*Eb for x in range(L))[b_]
    g0 = np.linalg.eigh(Hh)[1][:, 0]
    psi_h = np.zeros(Dd_)
    psi_h[b_] = g0
    prev, verfolgt = g0.astype(complex), True
    for u in np.linspace(0.05, 1.0, 20 if L < 12 else 10):
        Fm = F_anw(Eb.astype(complex), L, u)[b_]
        ev, Wv = np.linalg.eig(Fm)
        k_ = int(np.argmax(np.abs(np.linalg.solve(Wv, prev))))
        Qr = np.linalg.qr(Wv[:, np.abs(ev - ev[k_]) < 1e-8])[0]
        neu = Qr @ (Qr.conj().T @ prev)
        verfolgt &= bool(np.linalg.norm(neu) > 0.9)
        prev = neu/np.linalg.norm(neu)
    psi = np.zeros(Dd_, complex)
    psi[b_] = prev
    Fpsi = F_anw(psi[:, None], L)[:, 0]
    lam_ = np.vdot(psi, Fpsi)
    kette[L] = {"S_Heisenberg": verschr(psi_h, L), "S_Vakuum": verschr(psi, L), "Residuum": float(np.linalg.norm(Fpsi - lam_*psi)),
                "verfolgt": verfolgt}
Lk = np.array(sorted(kette))
c_vak = float(3*np.polyfit(np.log(Lk), [kette[L]["S_Vakuum"] for L in Lk], 1)[0])
c_hei = float(3*np.polyfit(np.log(Lk), [kette[L]["S_Heisenberg"] for L in Lk], 1)[0])
z121_S = Z121["zahlen"]["kontinuum"]["S(L/2) Vakuum"]
z121_abw = max(abs(kette[L]["S_Vakuum"] - z121_S[str(L)]) for L in Lk)
ZAHLEN["C_T_N"] = {"C_T(N)/C_T(1), unterscheidbare Leitern": {str(k): round(v, 10) for k, v in ct_verh.items()},
                   "Zahl der erhaltenen Spin-2-Ströme im Zweiteilchenraum": {str(k): v for k, v in anz_strom.items()},
                   "Artentausch: Änderung von T_ges, T_11 - T_22": [round(inv_ges, 10), round(inv_diff, 10)],
                   "[F, Uhr^(xN)] < 1e-12, [F, Uhr auf einer Naht]": [bool(komm_uni < 1e-12), round(komm_ein, 6)],
                   "Kette S(N/2): Vakuum": {str(L): round(kette[L]["S_Vakuum"], 5) for L in Lk},
                   "Kette S(N/2): Heisenberg": {str(L): round(kette[L]["S_Heisenberg"], 5) for L in Lk},
                   "c_2D aus (c/3) ln N: Vakuum, Heisenberg": [round(c_vak, 4), round(c_hei, 4)],
                   "Abweichung von z121 (S Vakuum) < 1e-4": bool(z121_abw < 1e-4)}
check(f"[E] C_T auf N Nähten: für N unterscheidbare Leitern (P_ord) enthält der Zweiteilchenraum N (N + 1)/2 = "
      f"{[anz_strom[N] for N in sorted(anz_strom)]} erhaltene Spin-2-Ströme T_ij (N = 1, 2, 3, 4, 6); der Spannungstensor "
      f"ist das Singulett T_ges = Summe T_ii mit C_T(N)/C_T(1) = {[round(ct_verh[N], 6) for N in sorted(ct_verh)]}: "
      f"additiv, c = N c_1 (für entkoppelte Arten nach Konstruktion). T_ges ist unter dem Artentausch invariant "
      f"({inv_ges < 1e-12}), T_11 - T_22 nicht (Änderung {inv_diff:g}); das Gatter a + b P von z121 ändert C_T also nicht. Unter P_ident "
      f"(N Quanten derselben Leiter) bleibt C_T = C_T(1): c gehört zum Netz, nicht zur Besetzung. Die einheitliche Uhr "
      f"vertauscht mit der Kette ({komm_uni < 1e-12}), eine Uhr auf einer einzelnen Naht nicht (Kommutator {komm_ein:.3f})",
      all(anz_strom[N] == N*(N + 1)//2 for N in anz_strom) and all(abs(ct_verh[N] - N) < 1e-12 for N in ct_verh)
      and inv_ges < 1e-12 and inv_diff > 0.5 and komm_uni < 1e-12 and komm_ein > 0.1)
check(f"[E/N] Die Kette selbst zählt nicht mit: ihr Vakuum (Floquet-Eigenzustand bei u = 1, stetig aus dem "
      f"Heisenberg-Grundzustand verfolgt: {all(kette[L]['verfolgt'] for L in Lk)}; Werte wie z121 auf 1e-4: "
      f"{z121_abw < 1e-4}) hat zwischen den Nähten die Verschränkung S(N/2) = "
      f"{[round(kette[L]['S_Vakuum'], 4) for L in Lk]} (N = 6, 8, 10, 12), also (c/3) ln N mit c_2D = {c_vak:.3f} "
      f"(Heisenberg {c_hei:.3f}): logarithmisch, nicht proportional zu N. Die 2D-Zentralladung der Kette ist eine Zahl "
      f"der Ordnung 1 für die ganze Kette und nicht der 4D-Koeffizient c; additiv in N ist nur der Leiterinhalt",
      all(kette[L]["verfolgt"] and kette[L]["Residuum"] < 1e-9 for L in Lk) and z121_abw < 1e-4 and 0.75 < c_vak < 1.4
      and 0.75 < c_hei < 1.2 and kette[12]["S_Vakuum"] < 0.5*12*np.log(2)/2)
befund("E/C: Unter P_ord (unterscheidbare, geordnete Nähte als direkte Summe) sind alle drei nativen Größen additiv: "
       "Uhrspur 64 N, Entropie N ln 2/2, C_T = N C_T(1); die Kette von z121 ändert keine davon. Unter P_ident oder im "
       "graduierten Produkt zählt nichts (C_T = C_T(1), Index 1). Dann ist l^2/G = (4 c_1/pi) L N: proportional zur "
       "Nahtzahl", True)

# --- 4.4 Was die Zählung verlangt: Vergleichsziel COSMO.LAM.01 / LAMBDA.BRANCH.01 (Status 'Numerical', nur gelesen)
log10 = lambda x: float(mp.log10(x))
ziel_SdS = mp.pi*ziel                                                # S_dS = pi l^2/G = 32 pi^4 e^(2/alpha)
ziel_z98 = Z98["zahlen"]["skalen"]["log10 l^2/G"]
pref = sp.simplify(48*c3**2)
rho_s, Mb = sp.symbols("rho Mbar", positive=True)
l2G_aus_rho = sp.simplify(24*PI/(pref*sp.exp(-2*sp.Symbol("ainv"))))  # rho/Mbar^4 = 24 pi G/l^2 (rho = 3/(8 pi G l^2), Mbar^-2 = 8 pi G)
EINH = [("native Leiter", Fr(1, 120)), ("ein Weyl-Fermion", Fr(1, 40)), ("eine 16", Fr(2, 5)), ("3 x 16", Fr(6, 5)),
        ("3 x 16 + 12 Eichbosonen", Fr(12, 5))]
NL = {name: 8*mp.pi**4*e2a/(mp.mpf(c_.numerator)/c_.denominator) for name, c_ in EINH}      # N L = 8 pi^4 e^(2/alpha)/c_1
s_je_naht = {name: 4*c_ for name, c_ in EINH}                         # S_dS/N = pi (l^2/G)/N = 4 c_1 L
L_fuer_ln2 = {name: float(mp.log(2)/2/(4*mp.mpf(c_.numerator)/c_.denominator)) for name, c_ in EINH}
N_jacobson = ziel_SdS/(mp.log(2)/2)
zwei_a = 2*ALPHA_INV
bits = {"2/(alpha ln 2)": zwei_a/mp.log(2), "(2/alpha)/(ln 2/2)": zwei_a/(mp.log(2)/2), "(2/alpha)/ln 12": zwei_a/mp.log(12)}
abst = {k: float(abs(v - mp.nint(v))) for k, v in bits.items()}
# Coleman-Weinberg-Faktor: d^3/ds^3 von (1/2) Int d^4k/(2 pi)^4 ln(k^2 + s) = 1/(32 pi^2 s) <-> s^2 ln s/(64 pi^2)
kk, ss = sp.symbols("k s", positive=True)
cw3 = sp.simplify(sp.integrate(2*PI**2*kk**3/(2*PI)**4/(kk**2 + ss)**3, (kk, 0, sp.oo)))
cw_koeff = sp.simplify(cw3/sp.diff(ss**2*sp.log(ss), ss, 3))
Nw = sp.Symbol("N_W", positive=True)
exp_ledger = sp.simplify(sp.diff(sp.log(24*PI/(Nw*c3**2)), Nw)*Nw)    # Ledgerform mit 48 als Zahl N_W gelesen
exp_ind = sp.simplify(sp.diff(sp.log(4*(Nw/40)*Ls/PI), Nw)*Nw)        # induziert: c = N_W/40
ZAHLEN["vergleichsziel"] = {"log10 e^(2/alpha)": log10(e2a), "log10 l^2/G (Ziel)": log10(ziel), "z98": ziel_z98,
                            "log10 S_dS (Ziel)": log10(ziel_SdS), "48 c_3^2": str(pref), "l^2/G aus rho_Lambda": str(l2G_aus_rho),
                            "log10 N L je Einheit": {k: round(log10(v), 4) for k, v in NL.items()},
                            "nötige Entropie je Naht S_dS/N": {k: f"{v} L" for k, v in s_je_naht.items()},
                            "L, damit S_dS/N = ln 2/2": {k: round(v, 6) for k, v in L_fuer_ln2.items()},
                            "log10 N bei ln 2/2 je Naht (Flächenzählung)": log10(N_jacobson),
                            "e^(2/alpha) als Dimension": {k: mp.nstr(v, 9) for k, v in bits.items()},
                            "Coleman-Weinberg-Koeffizient je reellem Freiheitsgrad": str(cw_koeff), "c_3^2": str(sp.simplify(c3**2)),
                            "Exponent d ln(l^2/G)/d ln N_W: Ledgerform, induziert": [str(exp_ledger), str(exp_ind)]}
check(f"[E/C] Vergleichsziel (Ledger 'Numerical', kein hergeleiteter Eingang): rho_Lambda/Mbar^4 = 48 c_3^2 e^(-2/alpha), "
      f"48 c_3^2 = {pref}, gleichwertig l^2/G = {l2G_aus_rho} = 32 pi^3 e^(2/alpha) = 10^{log10(ziel):.3f} (z98: "
      f"{ziel_z98:.3f}), S_dS = 32 pi^4 e^(2/alpha) = 10^{log10(ziel_SdS):.3f}. Mit l^2/G = (4 c_1/pi) L N verlangt es "
      f"N L = 8 pi^4 e^(2/alpha)/c_1: 10^{log10(NL['native Leiter']):.2f} (Leiter), 10^{log10(NL['ein Weyl-Fermion']):.2f} "
      f"(Weyl), 10^{log10(NL['3 x 16']):.2f} (3 x 16), 10^{log10(NL['3 x 16 + 12 Eichbosonen']):.2f} (mit Eichbosonen). "
      f"Eine Gleichung für das Produkt N L: jede Nahtzahl passt bei geeignetem L, eine Vorhersage ist das nicht. Die "
      f"nötige Entropie je Naht ist S_dS/N = 4 c_1 L (zielunabhängig); ln 2/2 je Naht (z68) hieße L = "
      f"{L_fuer_ln2['native Leiter']:.4f} = 15 ln 2 (Leiter) bzw. {L_fuer_ln2['3 x 16']:.4f} (3 x 16), ohne native "
      f"Herkunft. Reine Flächenzählung mit ln 2/2 je Naht: N = 10^{log10(N_jacobson):.2f}",
      pref == 3/(4*PI**2) and sp.simplify(l2G_aus_rho - 32*PI**3*sp.exp(2*sp.Symbol('ainv'))) == 0
      and abs(log10(ziel) - ziel_z98) < 1e-6 and abs(L_fuer_ln2["native Leiter"] - 15*np.log(2)) < 1e-12
      and 123 < log10(NL["native Leiter"]) < 125 and 121 < log10(NL["3 x 16"]) < 122.5)
check(f"[E] Kann e^(2/alpha) oder 48 c_3^2 aus einer nativen Zählung kommen? (1) e^(2/alpha) = 10^{log10(e2a):.3f} ist als "
      f"additive Zahl N keine Zählung mit Mechanismus; als Dimension wäre es 2^M mit M = {mp.nstr(bits['2/(alpha ln 2)'], 8)}, "
      f"bzw. {mp.nstr(bits['(2/alpha)/(ln 2/2)'], 8)} Nähte zu ln 2/2 oder {mp.nstr(bits['(2/alpha)/ln 12'], 8)} "
      f"E8-Nähte zu ln 12: keine ganze Zahl (Abstände {[round(v, 3) for v in abst.values()]}), und eine Dimension ist "
      f"multiplikativ, die drei nativen Größen sind additiv. (2) c_3^2 = 1/(64 pi^2) ist exakt der Einschleifenfaktor "
      f"des quartischen Vakuumterms je reellem Freiheitsgrad ({cw_koeff}); 48 c_3^2 wäre der Betrag für 48 reelle "
      f"bosonische Freiheitsgrade, 48 Weyl-Fermionen haben 96 mit negativem Vorzeichen: Faktor 2 (dieselbe "
      f"Dirac-/Majorana-Weyl-Zweideutigkeit wie in z118) und Vorzeichen fehlen, und der quartische Term ist der am "
      f"wenigsten universelle. (3) Liest man 48 als Zahl N_W der Weyl-Fermionen, so fällt l^2/G in der Ledgerform wie "
      f"N_W^{exp_ledger}, induziert wächst es wie N_W^{exp_ind}: entgegengesetzte Abhängigkeit vom Inhalt",
      all(v > 0.1 for v in abst.values()) and sp.simplify(cw_koeff - c3**2) == 0 and exp_ledger == -1 and exp_ind == 1)
GLam_obs = 3/ziel
befund(f"E/O: Unter der Zählung bleibt als Prämisse genau das Produkt N L (N nicht nativ festgelegt, L nicht nativ: es "
       f"gibt keine native Abschneide-zu-Vakuum-Skala). Der beobachtete Betrag G Lambda = 10^{log10(GLam_obs):.2f} verlangt "
       f"c_ges L = 8 pi^4 e^(2/alpha) = 10^{log10(8*mp.pi**4*e2a):.2f}; mit nativem Inhalt der Ordnung 1 und L der Ordnung "
       f"1 bis 300 ist l^2/G von der Ordnung 0,01 bis 10^3, also G von der Ordnung R_Vak^2: um rund 120 Größenordnungen "
       f"verfehlt. Das Vorzeichen ist nativ AdS (Lambda < 0), beobachtet ist Lambda > 0: Konflikt, den die Zählung nicht "
       f"behebt. Weder e^(2/alpha) noch 48 c_3^2 kommen aus einer nativen Zählung", True)

# =================================================================== Teil 5
say("\nTeil 5: Bilanz und Urteil")
bilanz = [
    "[E] nativer Spannungstensor: erhaltene Delta = 4-Bilineare P_1P_1 - 4 P_1P_2 + P_2P_2 zweier Leitern, C_T = 1/(3 pi^4), "
    "c = 1/120 je reeller Leiter (Ward-Identität; Delta = 1 aus den nativen Vielfachheiten); Komposit, kein Graviton",
    "[E] freie Vergleichswerte in derselben Konvention: C_T Skalar : Weyl : Maxwell = 1 : 3 : 12; a = 1/360, 11/720, 31/180",
    "[E/L] Lesarten: (i) (1/360, 1/120), (ii) (11/15, 6/5), (iii) (14/5, 12/5); welche nativ ist, bleibt Lesart",
    "[E/L] universell nur der Weyl^2-Logarithmus mit c; Einstein-Hilbert nicht fest (z27; Wärmekern: (i) 0, (ii) 48, (iii) 0)",
    "[E] Maldacena: alpha_W = 3/(64 pi G Lambda), l^2/G = 64 pi |alpha_W| (AdS4, S^4, Hawking-Page); z118s Uhrspur = 16 W^2 "
    "je Index auf Einstein-Räumen",
    "[C] induziert: l^2/G = (4 c/pi) L, G Lambda = -3 pi/(4 c L); Vorzeichen zwingend AdS; L nicht nativ",
    "[E/C] N Nähte: Uhrspur 64 N, Entropie N ln 2/2, C_T = N C_T(1) unter P_ord; Kette ändert nichts; unter P_ident oder im "
    "graduierten Produkt keine Zählung",
    "[E/O] Vergleichsziel 32 pi^3 e^(2/alpha): verlangt N L c_1 = 8 pi^4 e^(2/alpha); e^(2/alpha) und 48 c_3^2 nicht aus "
    "einer Zählung; entgegengesetzte Abhängigkeit vom Inhalt; Vorzeichen falsch",
]
for z in bilanz:
    say("   " + z)
urteil = ("Induzierte Gravitation ersetzt P_G bedingt durch eine Zählung, beseitigt die Prämisse aber nicht; nichts ist "
          "befördert, eine vollständige Theorie wird nicht behauptet. (1) Der native Spannungstensor existiert: Aus dem "
          "nativen Kern det(I - Z^+ W)^-1 folgt Delta = 1, die Spin-2-Primäre zweier Leitern ist P_1P_1 - 4 P_1P_2 + "
          "P_2P_2, erhalten und bei Delta = 4; die Ward-Identität gibt C_T = 1/(3 pi^4), c = 1/120 je reeller Leiter, im "
          "Ortsraum gegengerechnet (Skalar : Weyl : Maxwell = 1 : 3 : 12). Nicht erhaltene Bilineare liegen nie bei "
          "Delta = 4 (|div T|^2 = 9 (Delta_T - 4) |T|^2). (2) Induzierte Wirkung: (a, c) = (1/360, 1/120) für eine Leiter, "
          "(11/15, 6/5) für 3 x 16 Weyl, (14/5, 12/5) mit zwölf Eichbosonen. Universell ist nur Gamma = -(L/(16 pi^2)) Int "
          "(c W^2 - a E_4); der Einstein-Hilbert-Koeffizient ist reglerabhängig (z27: Vorzeichenwechsel je Metrikkopplung; "
          "Wärmekern: null für die Leiter und für 48 Weyl + 12 Vektoren) und wird nicht als fest ausgegeben. (3) Auf "
          "Einstein-Lösungen ist die Einstein-Wirkung die Weyl-Wirkung mit alpha_W = 3/(64 pi G Lambda) plus Euler "
          "(geprüft an AdS4, S^4 und Hawking-Page), dieselbe Zahl wie z98s b und z118s Uhrspur (16 W^2 je Index). "
          "Rein induziert folgt l^2/G = (4 c/pi) L und G Lambda = -3 pi/(4 c L): je e-Faltung 0,0106 (Leiter), 1,528 "
          "(3 x 16), 3,056 (mit Eichbosonen); G Lambda = -282,7/L, -1,963/L, -0,982/L. Das Vorzeichen ist zwingend AdS, "
          "in Übereinstimmung mit dem nativen Vorzeichen (z104, z118) und im Widerspruch zur Beobachtung. (4) Unter P_ord "
          "sind Uhrspur (64 N), Entropie (N ln 2/2) und C_T (N C_T(1)) additiv, und die Kette von z121 ändert keine der "
          "drei; ihre eigene 2D-Zentralladung (c nahe 1) zählt nicht mit. Dann ist l^2/G = (4 c_1/pi) L N: die Kopplung "
          "in Einheiten des Vakuumradius zählt Nähte. Unter P_ident oder im graduierten Produkt zählt nichts. Das "
          "Vergleichsziel l^2/G = 32 pi^3 e^(2/alpha) (Ledger 'Numerical') verlangt N L c_1 = 8 pi^4 e^(2/alpha), rund "
          "10^122 bis 10^124 für N L; die nötige Entropie je Naht wäre 4 c_1 L, und ln 2/2 passt nur mit einem nicht "
          "nativen L = 15 ln 2. Weder e^(2/alpha) noch 48 c_3^2 kommen aus einer Zählung: e^(2/alpha) ist keine ganze "
          "Dimension und wäre multiplikativ, 48 c_3^2 trifft den quartischen Einschleifenfaktor nur bis auf Faktor 2 und "
          "Vorzeichen, und l^2/G hängt in der Ledgerform umgekehrt vom Inhalt ab wie in der induzierten Form. "
          "Einordnung für P_G: (a) bedingt, mit Rest (c), nicht (b). Der universelle Teil macht aus der Kopplung das "
          "Produkt Inhalt mal Logarithmus, c_ges L mit c_ges = N c_1 unter P_ord; der eigentliche "
          "Einstein-Hilbert-Koeffizient bleibt frei. Aus dem einen freien Parameter werden zwei nicht native Größen (N "
          "und L), fest wird keine Zahl. Prämissen der Route: P_netz (die Leitern bzw. die gelösten Blöcke als lokales "
          "3+1-Netz mit T_mn), P_ord (unterscheidbare Nähte als direkte Summe), P_ind (kein nackter Weyl- und "
          "Einstein-Term an der Nahtskala), P_cut (Abschneideskala, L), P_dyn (die Metrik ist dynamisch und nicht nur "
          "Quelle) und P_Neumann (Maldacenas Randbedingung auf einem (A)dS-Hintergrund); Hypothese T (z124) ersetzt P_ind "
          "durch eine feste konforme Kopplung, deren Zahl hier nicht berechnet ist. Graviton: Der universelle Teil liefert "
          "nur die metrische Antwort <T T> ~ c p^4 ln p^2 ohne masselosen Pol; zwei Helizitäten entstehen erst unter P_dyn "
          "und P_Neumann um (A)dS, nicht nativ. Mit nativem Inhalt der Ordnung 1 ist G von der Ordnung R_Vak^2, rund 120 "
          "Größenordnungen neben dem Vergleichsziel und mit dem falschen Vorzeichen von Lambda.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z137_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil, "bilanz": bilanz}, indent=1,
    ensure_ascii=False, default=str) + "\n")
