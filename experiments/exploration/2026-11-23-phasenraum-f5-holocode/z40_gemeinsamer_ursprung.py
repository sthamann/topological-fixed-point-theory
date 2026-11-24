"""Z40: Gemeinsamer Ursprung -- Spin(10)-Niveaus, Eulerbündel, Eichspektrum, Polarisation, Matrixkern, Raumzeit, Familien, Spin 2.

Prüft die neunte externe Analyse ("H6h22: gemeinsame positive Matrixgeometrie und relative Eichstruktur") und den
Polarisationsteil der achten:
Teil 1  Stromniveaus komplexer chiraler Spin(10)-Multipletts: 2T = (2, 16, 56, 4, 4) für 10, 45, 120, 16, 16b; alle gerade,
        Produktregel l(R x S) = dim S l(R) + dim R l(S); Niveau 1 nur mit reeller (Majorana-)Struktur.
Teil 2  Eulerbündel: ein Operator auf W = H^0(O(1)), der alle Kerne der Auswertung erhält, ist skalar; H^0(O^2) = H^0(O(1)).
Teil 3  Adjungierter Laplace zu nabla_t = d + (1-t)[P, dP] auf S^2: Quadratur bestätigt den Block
        [[L + 2t^2, -2t sqrt L], [-2t sqrt L, L - 1 + t^2]] und den Eigenwert L - 1 + t^2; das J = 1-Triplett ist nur bei t = 1
        masselos, das J = 0-Singulett nur bei t = 0.
Teil 4  Polarisation: Schleifenintegral 1/(4 pi), Koeffizient 1/(8 pi), nabla_t P = t dP, Int tr(dP)^2 = 4 pi, Gamma_pol = 2 t^2;
        V(t) = a(1-t^2)^2 - 2 log t + 3t^2 hat für jedes a >= 0 genau ein Minimum in (0,1), bei a = 1 t = 1/sqrt 2;
        dort ist das Triplett massiv -- Konflikt, wenn Randverbindung und schwache Verbindung dieselbe sind.
Teil 5  Klebeklassen |s_5|^2 = 5/4, |s_3|^2 = 3/4; 5 k_5 + 3 k_3 = 8 hat genau die positive Lösung (1,1).
Teil 6  Matrixkern det(I - Z W^+)^-nu: Rang-1-Schnitt = Hardy-Kern; zweite Ordnung nu(nu+1)/2 h_2 + nu(nu-1)/2 e_2; Cayley
        det(d) det^-nu = nu(nu-1) det^(-nu-1); Wallach-Test mit der binären Ikosaedergruppe (negativ bei nu = 1/2); Gaußidentität;
        det(lambda lambda^+) = 0; det(lambda lambda^+ + mu mu^+) = |lambda_1 mu_2 - lambda_2 mu_1|^2; Spin-2-Feldgleichung;
        K_2 = 24/det A Sym^4(A^-1) (exakt diagonal, Monte Carlo allgemein).
Teil 7  Tetraedrische Nullbasis: m^2 = [(Sum p)^2 - 3 Sum p^2]/4; positive Einzelenergien können raumartig sein.
Teil 8  Lichtkegelkompression: T(a) = e^{-i a m^2/p+}/(1 + i a sigma^2/p+), Defekt |T(2) - T(1)^2| = 1/sqrt 20.
Teil 9  Familien: h^0(E_-) = 3, h^1 = 0 (E_- = O(-1) + 3 O), h^0(E_+) = 5; eine Generation Lambda^gerade C^5 mit
        V4-Hyperladung anomaliefrei, 4 Dubletts; drei Generationen = 48 Weyl-Felder.
Teil 10 Spin 2: Källén-Lehmann verlangt Q^2 C(Q^2) >= Z; e^{-Q^2/M^2}/Q^2 verletzt das (Querbezug GRAV.AMPLITUDE.01 u. a.);
        C = Z/Q^2 + r/(Q^2 + M^2) gibt C^-1 = Q^2/Z - r Q^4/(M^2 Z^2) + ...
"""
import itertools
import json
from fractions import Fraction as Fr
from pathlib import Path
import csv
from math import factorial

import numpy as np
import sympy as sp
from numpy.polynomial.legendre import leggauss
from scipy.special import eval_legendre

OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
HERE = Path(__file__).parent
REPO = HERE.parents[2]


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val):
    BEFUNDE[name] = bool(val)
    say(f"   [Befund: {'ja' if val else 'nein'}] {name}")


# ------------------------------------------------------------------ Teil 1
say("Teil 1: Stromniveaus komplexer chiraler Spin(10)-Multipletts")
E = [tuple(Fr(1) if k == i else Fr(0) for k in range(5)) for i in range(5)]
vec = [v for i in range(5) for v in (E[i], tuple(-c for c in E[i]))]
add = lambda *vs: tuple(sum(c) for c in zip(*vs))
lam = lambda k: [add(*c) for c in itertools.combinations(vec, k)]
half = Fr(1, 2)
spin = [tuple(half * e for e in eps) for eps in itertools.product((1, -1), repeat=5)]
S16 = [w for w in spin if sum(1 for e in w if e < 0) % 2 == 0]
S16b = [w for w in spin if sum(1 for e in w if e < 0) % 2 == 1]
alpha = (Fr(1), Fr(-1), Fr(0), Fr(0), Fr(0))
dot = lambda u, v: sum(a * b for a, b in zip(u, v))
lvl = lambda W: sum(dot(w, alpha) ** 2 for w in W) / 2
reps = {"10": vec, "45": lam(2), "120": lam(3), "210": lam(4), "16": S16, "16b": S16b}
lv = {k: lvl(W) for k, W in reps.items()}
say(f"   Dimensionen {[(k, len(W)) for k, W in reps.items()]}; Niveau (komplex) = Sum (w.alpha)^2/2 = {dict((k, str(v)) for k, v in lv.items())}")
check("komplexe Niveaus 10, 45, 120, 16, 16b = 2, 16, 56, 4, 4 (= 2T mit T(10) = 1); auch 210 gerade",
      [lv[k] for k in ("10", "45", "120", "16", "16b")] == [2, 16, 56, 4, 4] and lv["210"] % 2 == 0)
prod = lambda A, B: [add(a, b) for a in A for b in B]
pr = {"10x16": prod(vec, S16), "16x16": prod(S16, S16), "10x10": prod(vec, vec), "16x16b": prod(S16, S16b)}
pl = {k: lvl(W) for k, W in pr.items()}
check("Produktregel l(R x S) = dim S l(R) + dim R l(S), z. B. 10 x 16 = 16b + 144 mit l(144) = 72 - 4 = 68",
      pl["10x16"] == 16 * lv["10"] + 10 * lv["16"] and pl["16x16"] == 2 * 16 * lv["16"] and pl["10x10"] == 2 * 10 * lv["10"]
      and pl["10x16"] - lv["16b"] == 68)
check("alle Fundamentaldarstellungen gerade => jedes komplexe Spin(10)-Multiplett hat gerades Niveau (Darstellungsring)",
      all(v % 2 == 0 for v in lv.values()) and all(v % 2 == 0 for v in pl.values()))
check("Niveau 1 braucht eine reelle Struktur: zehn Majorana-Fermionen im reellen Vektor geben lvl(10)/2 = 1", lv["10"] / 2 == 1)
befund("Gewöhnliche komplexe Flussreduktion kann Spin(10) auf Niveau 1 nicht erzeugen, wenn Spin(10) auf dem Elternoperator wirkt", True)
ZAHLEN["spin10_niveaus"] = {k: int(v) for k, v in lv.items()}

# ------------------------------------------------------------------ Teil 2
say("\nTeil 2: Schwaches Dublett und Eulerbündel")
rng = np.random.default_rng(39)
pts = rng.normal(size=(3, 2)) + 1j * rng.normal(size=(3, 2))
# s(z) = s0 z0 + s1 z1; Kern von ev_x: (x1, -x0)
rows = []
for x0, x1 in pts:
    v = np.array([x1, -x0])
    # T v parallel zu v  <=>  (Tv)_0 v_1 - (Tv)_1 v_0 = 0, linear in T = (t00, t01, t10, t11)
    rows.append([v[0] * v[1], v[1] * v[1], -v[0] * v[0], -v[1] * v[0]])
Amat = np.array(rows)
null = 4 - np.linalg.matrix_rank(Amat, tol=1e-10)
sv = np.linalg.svd(Amat)[2][-1]
Tn = sv.reshape(2, 2) / sv[0]
check("ein Operator auf W = H^0(O(1)), der die Kerne der Auswertung an drei Punkten erhält, ist skalar",
      null == 1 and np.allclose(Tn, np.eye(2), atol=1e-8))
hd = lambda d: (max(d + 1, 0), max(-d - 1, 0))
check("Euler: 0 -> O(-1) -> O^2 -> O(1) -> 0 mit h^0, h^1(O(-1)) = 0 gibt H^0(O^2) = H^0(O(1)) = C^2",
      hd(-1) == (0, 0) and 2 * hd(0)[0] == hd(1)[0] == 2)
befund("Zwei geometrische Nullmoden tragen keine lokale SU(2); im vollen Rang-2-Bündel O^2 wirkt SU(2) fasertreu", True)

# ------------------------------------------------------------------ Teil 3
say("\nTeil 3: Adjungierter Laplace zu nabla_t = d + (1-t)[P, dP] auf der Einheitssphäre")
th, ph, tt = sp.symbols("theta phi t", real=True)
nv = sp.Matrix([sp.sin(th) * sp.cos(ph), sp.sin(th) * sp.sin(ph), sp.cos(th)])
sig = [sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -sp.I], [sp.I, 0]]), sp.Matrix([[1, 0], [0, -1]])]
P = (sp.eye(2) + sum((nv[i] * sig[i] for i in range(3)), sp.zeros(2))) / 2
dP = [sp.diff(P, th), sp.diff(P, ph)]
nab = [dP[m] + (1 - tt) * ((P * dP[m] - dP[m] * P) * P - P * (P * dP[m] - dP[m] * P)) for m in range(2)]
check("nabla_t P = t dP (wegen P dP P = 0)", all(sp.simplify(nab[m] - tt * dP[m]) == sp.zeros(2) for m in range(2)))
trdP = sp.simplify((dP[0] * dP[0]).trace() + (dP[1] * dP[1]).trace() / sp.sin(th) ** 2)
intdP = sp.integrate(sp.integrate(trdP * sp.sin(th), (th, 0, sp.pi)), (ph, 0, 2 * sp.pi))
check("Int_{S^2} tr(dP)^2 = 4 pi (konform invariant)", sp.simplify(intdP - 4 * sp.pi) == 0)

xg, wg = leggauss(60)
thq = np.arccos(xg)
phq = np.linspace(0, 2 * np.pi, 48, endpoint=False)
TH, PH = np.meshgrid(thq, phq, indexing="ij")
W = np.outer(wg, np.full(len(phq), 2 * np.pi / len(phq)))
nq = np.stack([np.sin(TH) * np.cos(PH), np.sin(TH) * np.sin(PH), np.cos(TH)])
eth = np.stack([np.cos(TH) * np.cos(PH), np.cos(TH) * np.sin(PH), -np.sin(TH)])
eph = np.stack([-np.sin(PH), np.cos(PH), 0 * PH])
dn = [eth, np.sin(TH) * eph]


def felder(J):
    """Vektorharmonische mit m = 0: n Y, grad Y, n x grad Y, samt Ableitungen nach theta und phi (analytisch)."""
    c = np.cos(TH)
    Y = eval_legendre(J, c)
    dY = -np.sin(TH) * (J * (c * eval_legendre(J, c) - eval_legendre(J - 1, c)) / (c**2 - 1)) if J > 0 else 0 * c
    # zweite Ableitung von Y nach theta: Legendre-Gleichung Y'' = -cot Y' - J(J+1) Y
    d2Y = -np.cos(TH) / np.sin(TH) * dY - J * (J + 1) * Y
    # d e_theta/d theta = -n, d e_theta/d phi = cos th e_phi, d e_phi/d phi = -(sin th n + cos th e_theta), d e_phi/d theta = 0
    b = [nq * Y, eth * dY, eph * dY]
    db_th = [eth * Y + nq * dY, -nq * dY + eth * d2Y, eph * d2Y]
    db_ph = [np.sin(TH) * eph * Y, np.cos(TH) * eph * dY, -(np.sin(TH) * nq + np.cos(TH) * eth) * dY]
    return b, db_th, db_ph


def block_quadratur(J, t):
    b, dth, dph = felder(J)
    ks = [0] if J == 0 else [0, 1, 2]
    A_th = np.cross(nq, dn[0], axis=0)
    A_ph = np.cross(nq, dn[1], axis=0)
    D = lambda k: (dth[k] - (1 - t) * np.cross(A_th, b[k], axis=0), dph[k] - (1 - t) * np.cross(A_ph, b[k], axis=0))
    Dk = [D(k) for k in ks]
    G = np.array([[np.sum(W * np.sum(b[k] * b[l], axis=0)) for l in ks] for k in ks])
    Q = np.array([[np.sum(W * (np.sum(Dk[i][0] * Dk[j][0], axis=0) + np.sum(Dk[i][1] * Dk[j][1], axis=0) / np.sin(TH) ** 2))
                   for j in range(len(ks))] for i in range(len(ks))])
    return np.sort(np.linalg.eigvals(np.linalg.solve(G, Q)).real)


def block_formel(J, t):
    L = J * (J + 1)
    if J == 0:
        return np.array([2 * t * t])
    Bm = np.array([[L + 2 * t * t, -2 * t * np.sqrt(L)], [-2 * t * np.sqrt(L), L - 1 + t * t]])
    return np.sort(np.concatenate([np.linalg.eigvalsh(Bm), [L - 1 + t * t]]))


abw = 0.0
for J in range(0, 5):
    for t in (0.0, 0.3, 1 / np.sqrt(2), 1.0):
        abw = max(abw, np.max(np.abs(block_quadratur(J, t) - block_formel(J, t))))
say(f"   größte Abweichung Quadratur - Blockformel (J = 0..4, t = 0, 0.3, 1/sqrt2, 1): {abw:.2e}")
check("Quadratur bestätigt den 2x2-Block [[L+2t^2, -2t sqrt L],[-2t sqrt L, L-1+t^2]] und den Eigenwert L-1+t^2", abw < 1e-8)
m2 = lambda t: (3 * (1 + t * t) - np.sqrt(1 + 34 * t * t + t**4)) / 2
check("J = 1: unterster Eigenwert [3(1+t^2) - sqrt(1+34t^2+t^4)]/2, null nur bei t = 1 (Triplett masselos; m^2 >= det/Spur ~ (1-t^2)^2)",
      abs(m2(1.0)) < 1e-12 and all(m2(t) > 0.9 * (1 - t * t) ** 2 / (1 + t * t) / 2 for t in np.linspace(0, 0.99, 100)) and abs(block_formel(1, 0.4)[0] - m2(0.4)) < 1e-12)
check("J = 0: Eigenwert 2t^2, null nur bei t = 0 (eine U(1)-Richtung)", abs(block_formel(0, 0.0)[0]) < 1e-15 and block_formel(0, 0.5)[0] > 0)
ZAHLEN["triplett_masse2"] = {f"{t:.4f}": float(m2(t)) for t in (0.0, 0.5, 1 / np.sqrt(2), 0.9, 1.0)}

# ------------------------------------------------------------------ Teil 4
say("\nTeil 4: Fermionisch induzierte Polarisation")
p_, m_ = sp.symbols("p m", positive=True)
loop = sp.integrate(2 * sp.pi * p_ * m_**2 / (p_**2 + m_**2) ** 2, (p_, 0, sp.oo)) / (2 * sp.pi) ** 2
check("Int d^2p/(2 pi)^2 m^2/(p^2+m^2)^2 = 1/(4 pi)", sp.simplify(loop - 1 / (4 * sp.pi)) == 0)
# Gamma = -1/2 log det(D^+ D), D^+ D = -d^2 + m^2 - m gamma.dU; zweite Ordnung -1/2 Tr(GVGV): tr(gamma^mu gamma^nu) = 2 delta
koeff = sp.Rational(1, 2) * sp.Rational(1, 2) * 2 * loop
check("Ableitungsentwicklung: Gamma_pol = (1/(8 pi)) Int tr(nabla(2P - I))^2 (Koeffizient 1/2 x 1/2 x 2 x 1/(4 pi))", sp.simplify(koeff - 1 / (8 * sp.pi)) == 0)
Gpol = sp.simplify(koeff * 4 * tt**2 * intdP)
check("auf der runden Projektorsphäre Gamma_pol = 2 t^2 pro massivem komplexem Dublett", sp.simplify(Gpol - 2 * tt**2) == 0)
a_ = sp.symbols("a", nonnegative=True)
V = a_ * (1 - tt**2) ** 2 - 2 * sp.log(tt) + 3 * tt**2
s_ = sp.symbols("s", positive=True)
quad = sp.expand(sp.diff(V, tt) * tt).subs(tt**4, s_**2).subs(tt**2, s_)
tstar = {}
for av in (0, sp.Rational(1, 2), 1, 2, 10):
    roots = [r for r in sp.solve(quad.subs(a_, av), s_) if r.is_positive]
    tstar[str(av)] = [float(sp.sqrt(r)) for r in roots]
say(f"   stationäre Punkte t* (positiv): {tstar}")
check("V(t) hat für a = 0, 1/2, 1, 2, 10 genau einen stationären Punkt, und der liegt in (0,1) (Minimum: V -> oo bei t -> 0, V'(1) = 4)",
      all(len(v) == 1 and 0 < v[0] < 1 for v in tstar.values()) and sp.diff(V, tt).subs(tt, 1) == 4)
check("bei a = 1: t* = 1/sqrt 2, t*^2 = 1/2", abs(tstar["1"][0] - 1 / np.sqrt(2)) < 1e-12)
mt = m2(1 / np.sqrt(2))
say(f"   Triplett-Massenquadrat bei t* = 1/sqrt 2: {mt:.6f} / R^2")
check("bei t* = 1/sqrt 2 ist das schwache Triplett massiv (m^2 R^2 = (4.5 - sqrt 18.25)/2 = 0.114)", abs(mt - (4.5 - np.sqrt(18.25)) / 2) < 1e-12 and mt > 0.1)
befund("Stabilisierter Kanal und masselose SU(2) schließen sich aus, wenn Randverbindung = schwache Verbindung (Y1/Z28: eigene Eichverbindung)", True)

# ------------------------------------------------------------------ Teil 5
say("\nTeil 5: Zentraler Kleber D5 + A3 und Stromniveaus")
s5 = [half] * 5
s3 = [half] * 3
check("Klebeklassen |s_5|^2 = 5/4, |s_3|^2 = 3/4, zusammen 2 (E8-Wurzel)", sum(c * c for c in s5) == Fr(5, 4) and sum(c * c for c in s3) == Fr(3, 4))
sol = [(k5, k3) for k5 in range(1, 20) for k3 in range(1, 20) if 5 * k5 + 3 * k3 == 8]
mod8 = [(k5, k3) for k5 in range(1, 9) for k3 in range(1, 9) if (5 * k5 + 3 * k3) % 8 == 0]
say(f"   5 k5 + 3 k3 = 8: {sol}; nur ganzzahlige Antwort (= 0 mod 8): {mod8}")
check("5 k_5 + 3 k_3 = 8 hat genau die positive Lösung (1,1), c = 5 + 3 = 8; Ganzzahligkeit allein lässt mehrere zu",
      sol == [(1, 1)] and len(mod8) > 1)

# ------------------------------------------------------------------ Teil 6
say("\nTeil 6: Matrixkern det(I - Z W^+)^-nu")
nu = sp.symbols("nu", positive=True)
z, w = sp.symbols("z w")
Zs = sp.Matrix([[z, 0], [0, 0]])
Ws = sp.Matrix([[w, 0], [0, 0]])
check("Rang-1-Schnitt: det(I - Z W^+) = 1 - z wbar, bei nu = 1 der Hardy-Kern", sp.simplify((sp.eye(2) - Zs * Ws.H).det() - (1 - z * sp.conjugate(w))) == 0)
eps = sp.symbols("epsilon")
Xs = sp.Matrix(2, 2, sp.symbols("x11 x12 x21 x22"))
ser = sp.series(sp.exp(-nu * sp.log((sp.eye(2) - eps * Xs).det())), eps, 0, 3).removeO()
c2 = sp.expand(ser.coeff(eps, 2))
p1s, p2s = Xs.trace(), (Xs * Xs).trace()
h2, e2 = (p1s**2 + p2s) / 2, (p1s**2 - p2s) / 2
check("zweite Ordnung = nu(nu+1)/2 h_2 + nu(nu-1)/2 e_2 (Sym^2- und det-Kanal, verallgemeinerte Pochhammer-Symbole)",
      sp.simplify(c2 - (nu * (nu + 1) / 2 * h2 + nu * (nu - 1) / 2 * e2)) == 0 and sp.simplify(e2 - Xs.det()) == 0)
As = sp.Matrix(2, 2, sp.symbols("A11 A12 A21 A22"))
dA = As.det()
cay = sp.diff(dA**(-nu), As[0, 0], As[1, 1]) - sp.diff(dA**(-nu), As[0, 1], As[1, 0])
check("Cayley: (d11 d22 - d12 d21) det(A)^-nu = nu(nu-1) det(A)^(-nu-1); masselose Wellengleichung genau bei nu = 1 (und 0)",
      sp.simplify(cay - nu * (nu - 1) * dA ** (-nu - 1)) == 0)


def binaer_ikosaeder():
    phi = (1 + 5**0.5) / 2
    q2m = lambda q: np.array([[q[0] + 1j * q[1], q[2] + 1j * q[3]], [-q[2] + 1j * q[3], q[0] - 1j * q[1]]])
    gens = [q2m(np.array([1, 1, 1, 1]) / 2), q2m(np.array([phi, 1 / phi, 1, 0]) / 2)]
    G = [np.eye(2, dtype=complex)]
    frontier = list(G)
    while frontier:
        new = []
        for g in frontier:
            for s in gens:
                h = s @ g
                if not any(np.allclose(h, k, atol=1e-9) for k in G):
                    G.append(h)
                    new.append(h)
        frontier = new
    return G


G2I = binaer_ikosaeder()
chi = lambda jj, g: np.sum(np.exp(1j * np.arange(-2 * jj, 2 * jj + 1, 2) * np.angle(np.linalg.eigvals(g)[0]))).real
mitt = [np.mean([chi(jj, g) for g in G2I]) for jj in range(1, 6)]
check("binäre Ikosaedergruppe: 120 Elemente, mittelt alle SU(2)-Charaktere j = 1..5 exakt weg", len(G2I) == 120 and max(abs(m) for m in mitt) < 1e-9)


def wallach_form(nuv, r=0.6, K=8):
    pts = [r * np.exp(2j * np.pi * k / K) * g for k in range(K) for g in G2I]
    c = np.array([np.exp(-2j * 2 * np.pi * k / K) for k in range(K) for _ in G2I]) / (K * len(G2I))
    Kmat = np.empty((len(pts), len(pts)), complex)
    for i, Zi in enumerate(pts):
        for jx, Wj in enumerate(pts):
            ev_ = np.linalg.eigvals(np.eye(2) - Zi @ Wj.conj().T)
            Kmat[i, jx] = np.prod(ev_ ** (-nuv))
    return float((c @ Kmat @ c.conj()).real), float(np.linalg.eigvalsh((Kmat + Kmat.conj().T) / 2).min())


wf = {nv_: wallach_form(nv_) for nv_ in (0.5, 1.0, 2.0)}
say(f"   Testform im (U(1)-Ladung 2, SU(2)-trivial)-Kanal und kleinster Gram-Eigenwert: {wf}")
check("Wallach: bei nu = 1/2 negative Testform und negativer Gram-Eigenwert, bei nu = 1 null (bis Aliasing), bei nu = 2 positiv; Testform = r^4 nu(nu-1)/2 bis Aliasing r^20",
      wf[0.5][0] < -1e-4 and wf[0.5][1] < 0 and abs(wf[1.0][0]) < 1e-6 and wf[2.0][0] > 1e-4
      and abs(wf[0.5][0] - 0.6**4 * (0.5 * -0.5 / 2)) < 1e-5 and abs(wf[2.0][0] - 0.6**4 * (2 * 1 / 2)) < 1e-4)
lam1, lam2, mu1, mu2 = sp.symbols("l1 l2 m1 m2")
lv_ = sp.Matrix([lam1, lam2])
mv_ = sp.Matrix([mu1, mu2])
Pm = lv_ * lv_.H
check("P = lambda lambda^+: det P = 0 (masselos), tr P = |lambda|^2 > 0 (Vorwärtskegel)", sp.simplify(Pm.det()) == 0)
check("det(lambda lambda^+ + mu mu^+) = |lambda_1 mu_2 - lambda_2 mu_1|^2: zwei freie Quanten nicht masselos, außer kollinear",
      sp.simplify((Pm + mv_ * mv_.H).det() - (lam1 * mu2 - lam2 * mu1) * sp.conjugate(lam1 * mu2 - lam2 * mu1)) == 0)
epsM = sp.Matrix([[0, 1], [-1, 0]])
check("Spin-2-Feldgleichung für phi_ABCD = lambda_A lambda_B lambda_C lambda_D e^{ixP}: Kontraktion eps^AB lambda_A lambda_B = 0",
      sp.simplify((lv_.T * epsM * lv_)[0]) == 0)
Ah = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
Ah = Ah @ Ah.conj().T + 0.5 * np.eye(2)
Acm = Ah + 0.6j * np.array([[0.4, 0.1], [0.1, -0.2]])
nS = 400000
zz = (rng.normal(size=(nS, 2)) + 1j * rng.normal(size=(nS, 2))) / np.sqrt(2)
# Importance Sampling mit dem hermiteschen Teil H: Int e^{-z^+ A z} = (1/det H) E_H[e^{-i z^+ K z}], |Integrand| <= 1
Hh = (Acm + Acm.conj().T) / 2
Kh = (Acm - Acm.conj().T) / 2j
zH = zz @ np.linalg.cholesky(np.linalg.inv(Hh)).T
est = np.mean(np.exp(-1j * np.einsum("ni,ij,nj->n", zH.conj(), Kh, zH))) / np.linalg.det(Hh).real
say(f"   Gaußidentität: MC {est:.6f}, 1/det A = {1 / np.linalg.det(Acm):.6f}")
check("Gaußidentität 1/det A = Int e^{-lambda^+ A lambda} d^4 lambda / pi^2 auch für komplexes A mit positivem hermiteschem Teil (MC)",
      abs(est - 1 / np.linalg.det(Acm)) / abs(1 / np.linalg.det(Acm)) < 1e-2 and abs(np.imag(1 / np.linalg.det(Acm))) > 1e-2 * abs(1 / np.linalg.det(Acm)))
Psym = np.zeros((16, 16))
for perm in itertools.permutations(range(4)):
    for idx in itertools.product(range(2), repeat=4):
        i = int("".join(map(str, idx)), 2)
        j = int("".join(str(idx[p]) for p in perm), 2)
        Psym[i, j] += 1 / 24
k4 = lambda M: np.kron(np.kron(M, M), np.kron(M, M))


def K2_formel(A):
    C = np.linalg.inv(A)
    return 24 / np.linalg.det(A) * k4(C) @ Psym


Ad = np.diag([1.3, 0.7])
Kd = K2_formel(Ad)
# Diagonalelement zum Tensorindex (1^k 2^(4-k)): <|lambda_1|^2k |lambda_2|^2(4-k)> = k!(4-k)!/(a_1^(k+1) a_2^(5-k)) (mit 1/det A)
diag_formel = [Kd[int("0" * k + "1" * (4 - k), 2), int("0" * k + "1" * (4 - k), 2)] for k in range(5)]
diag_exakt = [factorial(k) * factorial(4 - k) / (1.3 ** (k + 1) * 0.7 ** (5 - k)) for k in range(5)]
check("K_2 diagonal: Momente <|lambda_1|^2k |lambda_2|^2(4-k)> = k!(4-k)!/(a_1^(k+1) a_2^(5-k)) aus 24/det A Sym^4(A^-1)",
      np.allclose(diag_formel, diag_exakt, rtol=1e-12))
Cq = np.linalg.inv(Ah)
Lc = np.linalg.cholesky(Cq)
lamS = zz @ Lc.T
v4 = np.einsum("ni,nj,nk,nl->nijkl", lamS, lamS, lamS, lamS).reshape(nS, 16)
mc = (v4.T @ v4.conj()) / nS / np.linalg.det(Ah)
rel = np.linalg.norm(mc - K2_formel(Ah)) / np.linalg.norm(K2_formel(Ah))
say(f"   K_2 Monte Carlo ({nS} Stichproben) gegen 24/det A Sym^4(A^-1): relative Abweichung {rel:.4f}")
check("K_2 = 24/det A Sym^4(A^-1) für zufälliges positives A (Monte Carlo, < 3 %)", rel < 0.03)
befund("nu = 1 ist der kleinste positive Wallach-Punkt und genau der Punkt der masselosen Wellengleichung (Leiterdarstellung von U(2,2))", True)

# ------------------------------------------------------------------ Teil 7
say("\nTeil 7: Tetraedrische Nullbasis und gemeinsame Spektralbedingung")
Ee, px, py, pz = sp.symbols("E px py pz", real=True)
vt = [sp.Matrix(v) / sp.sqrt(3) for v in ((1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1))]
pv = sp.Matrix([px, py, pz])
pis = [Ee - (v.T * pv)[0] for v in vt]
check("Nullrichtungen (1, v_i) mit regulärem Tetraeder: m^2 = E^2 - |p|^2 = [(Sum p_i)^2 - 3 Sum p_i^2]/4",
      sp.simplify(sum(pis) ** 2 - 3 * sum(p**2 for p in pis) - 4 * (Ee**2 - px**2 - py**2 - pz**2)) == 0)
tab = {str(t): Fr(sum(t) ** 2 - 3 * sum(c * c for c in t), 4) for t in ((1, 1, 1, 1), (0, 1, 1, 1), (1, 0, 0, 0), (0, 0, 1, 1))}
say(f"   m^2 für positive Translationsenergien: {dict((k, str(v)) for k, v in tab.items())}")
check("Tabelle: (1,1,1,1) -> 1, (0,1,1,1) -> 0, (1,0,0,0) -> -1/2, (0,0,1,1) -> -1/2 (raumartig trotz positiver p_i)",
      list(tab.values()) == [1, 0, Fr(-1, 2), Fr(-1, 2)])
befund("Vier unabhängige positive Hardy-Translationen reichen nicht; nötig ist ein gemeinsamer Impuls im Kegel Herm(2)_+", True)

# ------------------------------------------------------------------ Teil 8
say("\nTeil 8: Lichtkegelkompression und Gedächtnis")
aa, pp, sg, mm, kk = sp.symbols("a p_plus sigma m k", positive=True)
Tint = sp.integrate(2 * sp.pi * kk / (sp.pi * sg**2) * sp.exp(-kk**2 / sg**2) * sp.exp(-sp.I * aa * kk**2 / pp), (kk, 0, sp.oo), conds="none")
Ta = sp.exp(-sp.I * aa * mm**2 / pp) * Tint
check("T(a) = <psi| e^{-i a P_-} |psi> = e^{-i a m^2/p+}/(1 + i a sigma^2/p+) für ein gaußsches Querprofil",
      sp.simplify(Ta - sp.exp(-sp.I * aa * mm**2 / pp) / (1 + sp.I * aa * sg**2 / pp)) == 0)
Tf = lambda av, mv=0.7: np.exp(-1j * av * mv**2) / (1 + 1j * av)
defekt = abs(Tf(2) - Tf(1) ** 2)
say(f"   |T(2) - T(1)^2| bei a = b = p+ = sigma = 1: {defekt:.11f} (1/sqrt 20 = {1 / np.sqrt(20):.11f})")
check("Kompositionsdefekt |T(2) - T(1)^2| = 1/sqrt 20, unabhängig von m (die volle Entwicklung ist eine Gruppe)",
      abs(defekt - 1 / np.sqrt(20)) < 1e-12 and abs(abs(Tf(2, 0.1) - Tf(1, 0.1) ** 2) - defekt) < 1e-12)

# ------------------------------------------------------------------ Teil 9
say("\nTeil 9: Familienoperator und eine Generation")
hP = lambda ds: tuple(sum(hd(d)[i] for d in ds) for i in range(2))
check("E_- = O(-1) + 3 O: h^0 = 3, h^1 = 0 (Index 3, keine Spiegel); E_+ = O(1) + 3 O: h^0 = 5",
      hP([-1, 0, 0, 0]) == (3, 0) and hP([1, 0, 0, 0]) == (5, 0))
y1 = [Fr(-1, 3)] * 3 + [Fr(1, 2)] * 2
t3c = [Fr(1, 2), Fr(-1, 2), Fr(0), Fr(0), Fr(0)]
hc = [Fr(1), Fr(1), Fr(-2), Fr(0), Fr(0)]
t3w = [Fr(0)] * 3 + [Fr(1, 2), Fr(-1, 2)]
zust = [S for k in (0, 2, 4) for S in itertools.combinations(range(5), k)]
Ysum = lambda S, v: sum(v[i] for i in S)
an = {
    "U1^3": sum(Ysum(S, y1) ** 3 for S in zust),
    "U1-grav": sum(Ysum(S, y1) for S in zust),
    "SU3^2 U1": sum(Ysum(S, y1) * Ysum(S, t3c) ** 2 for S in zust),
    "SU2^2 U1": sum(Ysum(S, y1) * Ysum(S, t3w) ** 2 for S in zust),
    "SU3^3": sum(Ysum(S, hc) ** 3 for S in zust),
}
dubletts = sum(1 for S in zust if len(set(S) & {3, 4}) == 1) // 2
say(f"   Anomalien pro Generation: {dict((k, str(v)) for k, v in an.items())}; SU(2)-Dubletts: {dubletts}")
check("eine Generation Lambda^gerade C^5 mit V4-Hyperladung (-1/3 x 3, 1/2 x 2): 16 Zustände, alle Anomalien 0, 4 Dubletts (gerade)",
      len(zust) == 16 and all(v == 0 for v in an.values()) and dubletts == 4)
check("drei Generationen: 3 x 16 = 48 linkschirale Weyl-Felder", 3 * len(zust) == 48)

# ------------------------------------------------------------------ Teil 10
say("\nTeil 10: Positivität des Spin-2-Propagators")
Qg = np.geomspace(1e-3, 1e3, 400)
ok_pos = True
for trial in range(20):
    s_k = rng.uniform(0.1, 50, size=5)
    w_k = rng.uniform(0, 3, size=5)
    Cq_ = 1.0 / Qg + np.sum(w_k[:, None] / (Qg[None, :] + s_k[:, None]), axis=0)
    ok_pos &= bool(np.all(Qg * Cq_ >= 1.0 - 1e-15))
check("positive Spektraldarstellung C = Z/Q^2 + Int drho/(Q^2+s) mit rho >= 0 erfüllt Q^2 C >= Z (20 Zufallsmaße)", ok_pos)
Cexp = np.exp(-Qg / 1.0) / Qg
check("exponentiell gedämpfter Propagator Z e^{-Q^2/M^2}/Q^2 verletzt Q^2 C >= Z für jedes Q^2 > 0", np.all(Qg * Cexp < 1.0))
Q2, Zs_, rs, Ms = sp.symbols("Q2 Z r M", positive=True)
Cinv = sp.series(1 / (Zs_ / Q2 + rs / (Q2 + Ms**2)), Q2, 0, 3).removeO()
check("C = Z/Q^2 + r/(Q^2+M^2): C^-1 = Q^2/Z - r Q^4/(M^2 Z^2) + O(Q^6)", sp.simplify(Cinv - (Q2 / Zs_ - rs * Q2**2 / (Ms**2 * Zs_**2))) == 0)
ledger = {r["claim_id"]: r for r in csv.DictReader(open(REPO / "verification" / "status_ledger.csv", newline="", encoding="utf-8"))}
quer = [k for k in ("GRAV.SPIN2.UNITARITY.01", "GRAV.KMS.HESSIAN.01", "QGAMB.IDG.01", "GRAV.AMPLITUDE.01", "GRAV.LOOP.FINITE.01",
                    "GRAV.NONLOCAL.01", "GRAV.SPIN2.EMERGENCE.01") if k in ledger]
amp = ledger.get("GRAV.AMPLITUDE.01", {}).get("claim", "")
check("Querbezug (nur lesend): sieben Ledgerzeilen nutzen den ganzen Formfaktor; GRAV.AMPLITUDE.01 nennt P(p^2) = e^{-p^2/M^2}/p^2",
      len(quer) == 7 and "e^{-p^2/M^2}/p^2" in amp)
befund("Die Ledger-Spin-2-Zeilen sichern Polfreiheit, nicht Spektralpositivität; als physischer positiver Propagator gelesen ist e^{-p^2/M^2}/p^2 ausgeschlossen", True)
ZAHLEN["ledger_querbezug"] = quer

urteil = ("Die neunte Analyse hält in allen geprüften Rechnungen. Komplexe chirale Spin(10)-Multipletts haben gerade Niveaus "
          "(2, 16, 56, 4, 4 für 10, 45, 120, 16, 16b; Produktregel), Niveau 1 braucht eine reelle Majorana-Struktur; zusammen "
          "mit Z39 (q-neutral sind nur 20 der 40 D5-Wurzeln) heißt das: Die Flussnaht trägt c = 8 und das E8-Gitter, aber "
          "Spin(10) auf Niveau 1 nicht als Symmetrie des Elternoperators. Zwei Schnitte von O(1) tragen keine lokale SU(2) "
          "(nur Skalare erhalten alle Auswertungskerne); im Eulerbündel O^2 schon. Die Quadratur bestätigt das Spektrum des "
          "adjungierten Laplace zu nabla_t: Das J = 1-Triplett ist nur bei t = 1 masselos. Die fermionisch induzierte "
          "Polarisation (Koeffizient 1/(8 pi), Gamma_pol = 2t^2) gibt V(t) = a(1-t^2)^2 - 2 log t + 3t^2 mit genau einem "
          "Minimum in (0,1), bei a = 1 t = 1/sqrt 2 -- dort ist das Triplett massiv (m^2 R^2 = 0.114). Randverbindung und "
          "schwache Verbindung müssen also verschieden sein. 5 k_5 + 3 k_3 = 8 hat genau (1,1). Der Matrixkern "
          "det(I - Z W^+)^-nu enthält den Hardy-Kern als Rang-1-Schnitt; der det-Kanal hat nu(nu-1)/2, die Cayley-Identität "
          "gibt nu(nu-1) det^(-nu-1), der Wallach-Test mit der binären Ikosaedergruppe ist bei nu = 1/2 negativ: nu = 1 ist der "
          "kleinste positive Punkt und genau der masselose. Tetraedrische Nullbasis, Lichtkegelkompression (Defekt 1/sqrt 20), "
          "Familienindex 3 mit anomaliefreier Generation und K_2 = 24/det A Sym^4(A^-1) rechnen sich nach. Spin 2: Eine positive "
          "Spektraldarstellung verlangt Q^2 C >= Z; der Formfaktor e^{-p^2/M^2}/p^2 aus GRAV.AMPLITUDE.01 und sechs weiteren "
          "Ledgerzeilen verletzt das. Diese Zeilen sichern Polfreiheit, nicht Positivität (Querbezug, Ledger unverändert). "
          "Offen bleiben die wechselwirkende Realisierung, G und die eindeutige Prozessregel Z[J,A,g,Phi].")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z40_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False, default=str) + "\n")
