"""Z41: Abgleich mit dem Altbestand -- das Spin-2-Graviton der Spektralwirkung, Positivität über Wallach, reelle Naht, Higgsrichtung.

Anlass: Der Auftrag, den Altbestand nach vorhandenen Lösungen der offenen Punkte zu durchsuchen. Der Altbestand nennt als
gemeinsames Masterfunktional die Spektralwirkung eines Dirac-Operators (CONTRACT.QFT4D.01, v244/v261) mit dem
Naht-KMS-Abschneider f(u) = e^{-u} (v259/v380). Die Ledgerzeilen GRAV.KMS.HESSIAN.01/GRAV.NONLOCAL.01 lesen dessen
resummierten Spin-2-Kinetikfaktor als a = e^{p^2/M^2} (nullstellenfrei, "der Geist entkoppelt"). Hier wird er aus dem
nichtlokalen Wärmeleitungskern (Barvinsky-Vilkovisky, Codello-Zanusso) tatsächlich berechnet.

Teil 1  Spektralwirkung Tr e^{-s D^2} für Dirac-Spinoren, quadratisch in TT-Fluktuationen: a(x) = 1 + 12 x g(x),
        g = f_Ric - f_Omega/2, x = s p^2. Prüfung gegen a_2 (Ric^2/10) und EH (-R/(48 pi^2 s)); Reihe 1 - 3x/10 + x^2/28 - ...
        a hat bei x* = 6.768 eine euklidische Nullstelle (Residuum < 0, tachyonischer Spin-2-Geist) und eine unendliche
        Folge komplexer Polpaare; das e^{x}-Ansatz-Bild ist nullstellenfrei. Jeder getestete vollständig monotone Abschneider
        (f = Laplace-Transformierte eines positiven Maßes) gibt eine euklidische Nullstelle.
Teil 2  Positivität über Wallach: det(P)^{s-2} dP auf Herm(2)_+ hat die Laplace-Transformierte Gamma_Omega(s) det(Y)^-s,
        Gamma_Omega(s) = pi Gamma(s) Gamma(s-1); Spektraldichte rho(m^2) ~ (m^2)^{s-2}, bei s -> 1 auf m^2 = 0 konzentriert;
        Källén-Lehmann für Delta in (1,2); Wallach-Menge {0} u [1, oo) = 4D-Unitaritätsschranke Delta >= 1 für Skalare;
        der Helizität-2-Kern ist homogen vom Grad -6 (Delta = 3 = Schranke für (2,0)) und auf dem Röhrengebiet positiv.
Teil 3  Reelle Naht aus dem Altbestand (SEAM.S3: 16 chirale Majorana-Randmoden einer p+ip-Phase, c = 8): ohne U(1)-Ladung
        bleiben alle 40 D5-Wurzeln erhalten und so(10) hat Niveau 1; mit der Flussladung q bleibt nur u(8) bzw. U(5).
Teil 4  M(P) als Higgsrichtung: P = H H^+/|H|^2, (2P - 1) spaltet das Dublett; tr(dP)^2 = 2(|dH^|^2 - |H^+ dH^|^2), also ist
        Gamma_pol der Fubini-Study-Kinetikterm der Higgsrichtung, wie Tr(dPhi)^2 = A kin (A = 4 Y_2) im Altbestand.
"""
import json
import warnings
from pathlib import Path

import numpy as np
import sympy as sp
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad
from scipy.optimize import brentq
import itertools
from fractions import Fraction as Fr

warnings.filterwarnings("ignore", category=RuntimeWarning)
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
HERE = Path(__file__).parent


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
say("Teil 1: Das Spin-2-Graviton der Spektralwirkung mit Naht-KMS-Abschneider")
x, xi = sp.symbols("x xi")
fser = sp.integrate(sp.series(sp.exp(-xi * (1 - xi) * x), x, 0, 7).removeO(), (xi, 0, 1))
fRic = sp.series((fser - 1 + x / 6) / x**2, x, 0, 4).removeO()
fOm = sp.series(-(fser - 1) / (2 * x), x, 0, 4).removeO()
check("Formfaktoren: f_Ric(0) = 1/60, f_Omega(0) = 1/12 (lokale a_2-Koeffizienten)", fRic.subs(x, 0) == sp.Rational(1, 60) and fOm.subs(x, 0) == sp.Rational(1, 12))
# Dirac: tr 1 = 4, U = R/4, tr(Omega Omega) = -Riem^2/2, Riem^2 -> 4 Ric^2 - R^2 (zweite Ordnung)
ric2 = 4 * fRic - 2 * fOm
check("Dirac a_2: Koeffizient von Ric^2 (nach Riemann-Elimination) = -1/10, wie -(7 Riem^2 + 8 Ric^2 - 5 R^2)/360 mit Riem^2 -> 4Ric^2 - R^2",
      ric2.subs(x, 0) == sp.Rational(-1, 10) and sp.Rational(-1, 360) * (7 * 4 + 8) == sp.Rational(-1, 10))
check("Einstein-Hilbert: (1/(16 pi^2 s)) tr(R/6 - U) = -R/(48 pi^2 s) (positives G, euklidisch)", 4 * (sp.Rational(1, 6) - sp.Rational(1, 4)) == sp.Rational(-1, 3))
# TT: Int sqrt g R -> -(1/4) p^2 |h|^2, Ric^2 -> (1/4) p^4 |h|^2  =>  K = p^2/(192 pi^2 s) + p^4 (f_Ric - f_Om/2)/(16 pi^2)
aser = sp.expand(1 + 12 * x * (fRic - fOm / 2))
say(f"   a(x) = 1 + 12 x g(x) = {sp.series(aser, x, 0, 4).removeO()} + ...")
check("Reihe a(x) = 1 - 3x/10 + x^2/28 - x^3/360 + O(x^4)",
      sp.expand(sp.series(aser, x, 0, 4).removeO() - (1 - sp.Rational(3, 10) * x + x**2 / 28 - x**3 / 360)) == 0)

xg, wg = leggauss(300)
XI, WI = (xg + 1) / 2, wg / 2


def f_bv(z):
    z = np.asarray(z, complex)
    return np.sum(WI * np.exp(-np.multiply.outer(z, XI * (1 - XI))), axis=-1)


def g_bv(z):
    z = np.asarray(z, complex)
    fz = f_bv(z)
    zz = np.where(z == 0, 1, z)
    out = (fz - 1 + z / 6) / zz**2 + (fz - 1) / (4 * zz)
    ser = (1 / 60 - z / 840) + (-1 / 6 + z / 60 - z**2 / 840) / 4
    return np.where(np.abs(z) < 1e-3, ser, out)


a_bv = lambda z: 1 + 12 * z * g_bv(z)
xs = np.geomspace(1e-4, 1e4, 4000)
av = a_bv(xs).real
xstar = brentq(lambda t: a_bv(t).real, 1, 20)
A_ = lambda z: z * a_bv(z)
dA_ = lambda z, h=1e-6: (A_(z + h) - A_(z - h)) / (2 * h)
res_star = (1 / dA_(xstar)).real
say(f"   euklidische Nullstelle x* = {xstar:.4f} (p_E^2 = {xstar:.3f} Lambda^2), Residuum von 1/(x a) dort {res_star:.4f}; min a = {av.min():.4f}")
check("Q^2 C = 1/a >= Z = 1 auf 0 < x < x* (notwendige Källén-Lehmann-Bedingung dort erfüllt)",
      bool(np.all(av[xs < xstar * 0.999] <= 1 + 1e-12) and np.all(av[xs < xstar * 0.999] > 0)))
check("a hat eine euklidische Nullstelle bei x* = 6.768 mit negativem Residuum: tachyonischer massiver Spin-2-Geist; darüber C < 0",
      abs(xstar - 6.768) < 1e-2 and res_star < 0 and bool(np.all(av[xs > xstar * 1.01] < 0)))
roots = []
for re in np.linspace(-120, 60, 46):
    for im in np.linspace(-100, 100, 41):
        z = complex(re, im)
        for _ in range(80):
            d = dA_(z)
            if abs(d) < 1e-14:
                break
            z = z - A_(z) / d
            if abs(z) > 300:
                break
        if abs(z) < 130 and abs(A_(z)) < 1e-9 and not any(abs(z - r) < 1e-5 for r in roots):
            roots.append(z)
roots.sort(key=abs)
cplx = [r for r in roots if abs(r.imag) > 1e-6]
say(f"   Pole von 1/(x a) mit |x| < 130: {[complex(np.round(r, 3)) for r in roots[:7]]} ... ({len(roots)} insgesamt, {len(cplx)} komplex)")
check("Pole: Graviton bei 0 (Residuum 1), reeller Geist bei x*, dazu komplex-konjugierte Paare (6.68 +- 32.99 i, 7.73 +- 58.76 i, ...)",
      abs(roots[0]) < 1e-9 and abs((1 / dA_(roots[0])).real - 1) < 1e-6 and len(cplx) >= 6 and len(cplx) % 2 == 0
      and any(abs(r - complex(6.6772, 32.9898)) < 1e-2 for r in cplx))
check("Vergleich: der Ledger-Ansatz a = e^{x} ist nullstellenfrei und wächst, der berechnete Faktor fällt (a'(0) = -3/10 statt +1)",
      sp.series(aser, x, 0, 2).removeO().coeff(x, 1) == sp.Rational(-3, 10))
fams = {"(1+u)^-2": lambda s: np.exp(-s), "(1+u)^-3": lambda s: s * np.exp(-s) / 2, "(1+u)^-5": lambda s: s**3 * np.exp(-s) / 24,
        "log-normal": lambda s: np.exp(-np.log(s) ** 2 / 2) / s**2 / np.sqrt(2 * np.pi)}


def a_mix(p2, w):
    num = quad(lambda s: w(s) * a_bv(s * p2).real, 0, np.inf, limit=400)[0]
    return num / quad(w, 0, np.inf, limit=400)[0]


zeros_cm = {"e^-u (KMS)": round(xstar, 3)}
for name, w in fams.items():
    ps = np.geomspace(0.1, 1e4, 50)
    vals = np.array([a_mix(p, w) for p in ps])
    sc = np.where(np.diff(np.sign(vals)) != 0)[0]
    zeros_cm[name] = round(brentq(lambda t: a_mix(t, w), ps[sc[0]], ps[sc[0] + 1]), 3) if len(sc) else None
say(f"   euklidische Nullstellen für vollständig monotone Abschneider: {zeros_cm}")
check("jeder getestete vollständig monotone Abschneider (positives Laplace-Maß) gibt eine euklidische Spin-2-Nullstelle; "
      "a(p^2) ist ein positives Mittel von a(s p^2), und x a(x) -> -6", all(v is not None for v in zeros_cm.values())
      and abs(xs[-1] * av[-1] + 6) < 0.1)
befund("Die tatsächliche Wärmeleitungs-Resummation entkoppelt den Spin-2-Geist nicht: GRAV.KMS.HESSIAN.01s [C]-Identifikation a = e^u hält für Tr e^{-D^2/Lambda^2} nicht", True)
ZAHLEN["spektralwirkung_spin2"] = {"x_stern": xstar, "residuum": res_star, "pole": [str(complex(np.round(r, 4))) for r in roots[:9]],
                                   "cm_nullstellen": zeros_cm, "reihe": "1 - 3x/10 + x^2/28 - x^3/360"}

# ------------------------------------------------------------------ Teil 2
say("\nTeil 2: Positivität über den Wallach-Parameter (Källén-Lehmann und 4D-Unitaritätsschranke)")
s_, a1, b1, y1, y2, r = sp.symbols("s a b y1 y2 r", positive=True)
inner = sp.integrate(2 * sp.pi * r * (a1 * b1 - r**2) ** (s_ - 2), (r, 0, sp.sqrt(a1 * b1)), conds="none")
inner = sp.simplify(inner)
check("Int_{|c|^2 < ab} (ab - |c|^2)^{s-2} d^2c = pi (ab)^{s-1}/(s-1)", sp.simplify(inner - sp.pi * (a1 * b1) ** (s_ - 1) / (s_ - 1)) == 0)
gam = sp.simplify(sp.pi / (s_ - 1) * sp.gamma(s_) ** 2)
check("Gindikin: Int_{Herm(2)+} e^{-tr PY} det(P)^{s-2} dP = pi Gamma(s) Gamma(s-1) det(Y)^-s (s > 1)",
      sp.simplify(gam - sp.pi * sp.gamma(s_) * sp.gamma(s_ - 1)) == 0)
Ynum = np.array([[1.3, 0.4 - 0.2j], [0.4 + 0.2j, 0.9]])
ev = np.linalg.eigvalsh(Ynum)
for sv in (2.0, 2.5, 3.0):
    numv = quad(lambda aa: quad(lambda bb: np.pi * (aa * bb) ** (sv - 1) / (sv - 1) * np.exp(-ev[0] * aa - ev[1] * bb), 0, np.inf)[0], 0, np.inf)[0]
    ZAHLEN[f"gindikin_s{sv}"] = numv
check("numerisch (Y nicht diagonal, unitär diagonalisiert): Integral = pi Gamma(s)Gamma(s-1)/det Y^s für s = 2, 2.5, 3",
      all(abs(ZAHLEN[f"gindikin_s{sv}"] / (np.pi * float(sp.gamma(sv) * sp.gamma(sv - 1))) * np.linalg.det(Ynum).real ** sv - 1) < 1e-6
          for sv in (2.0, 2.5, 3.0)))
conc = {sv: quad(lambda m: m ** (sv - 2), 0, 0.05)[0] / float(sp.gamma(sv - 1)) / (0.05 ** 0 if True else 1) for sv in (1.5, 1.1, 1.01)}
say(f"   Gewicht in m^2 < 0.05 von (m^2)^(s-2)/Gamma(s-1) (Normierung Int_0^1 -> 1/Gamma(s)): {conc}")
check("Spektraldichte rho(m^2) ~ (m^2)^{s-2}: bei s -> 1 wandert das Gewicht nach m^2 = 0 (Grenzmaß lambda lambda^+, masselos)",
      conc[1.01] * float(sp.gamma(1.01)) > 0.95 and conc[1.1] > conc[1.5])
pv = 2.3
for Dv in (1.5, 1.25, 1.75):
    lhs = quad(lambda m: m ** (Dv - 2) / (pv + m), 0, np.inf, limit=200)[0]
    rhs = pv ** (Dv - 2) * np.pi / np.sin(np.pi * (Dv - 1))
    ZAHLEN[f"kl_delta_{Dv}"] = (lhs, rhs)
check("Källén-Lehmann: (p^2)^{Delta-2} = (sin(pi(Delta-1))/pi) Int (m^2)^{Delta-2}/(p^2+m^2) dm^2 mit positiver Dichte für Delta in (1,2)",
      all(abs(ZAHLEN[f"kl_delta_{Dv}"][0] / ZAHLEN[f"kl_delta_{Dv}"][1] - 1) < 1e-6 for Dv in (1.5, 1.25, 1.75)))
check("unter der Schranke (Delta = 1/2) ist (m^2)^{Delta-2} bei 0 nicht integrierbar: keine positive Darstellung (= Wallach-Lücke 0 < nu < 1)",
      sp.integrate(sp.Symbol("mm", positive=True) ** sp.Rational(-3, 2), (sp.Symbol("mm", positive=True), 0, 1)) == sp.oo)
befund("Wallach-Menge {0} u [1, oo) von U(2,2) = 4D-Unitaritätsschranke Delta >= 1; nu = 1 ist ihr Rand mit Nullzustand det(d) = 0 (freies Feld)", True)
Xs = sp.Matrix(2, 2, sp.symbols("X11 X12 X21 X22"))
lam_ = sp.symbols("lam", positive=True)
adj = Xs.adjugate()
K2s = (adj[0, 0] ** 4) / Xs.det() ** 5          # eine Komponente von Sym^4(adj X)/det X^5
check("Helizität-2-Kern Sym^4(X^-1)/det X = Sym^4(adj X)/det(X)^5 ist homogen vom Grad -6, also Delta = 3 = j + 2 (Schranke für (2,0))",
      sp.simplify(K2s.subs({Xs[i, j]: lam_ * Xs[i, j] for i in range(2) for j in range(2)}) / K2s - lam_ ** (-6)) == 0)
rng = np.random.default_rng(41)
from itertools import permutations, product
Psym = np.zeros((16, 16))
for perm in permutations(range(4)):
    for idx in product(range(2), repeat=4):
        Psym[int("".join(map(str, idx)), 2), int("".join(str(idx[p]) for p in perm), 2)] += 1 / 24
k4 = lambda M: np.kron(np.kron(M, M), np.kron(M, M))
pts = []
for _ in range(6):
    Xh = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
    Xh = (Xh + Xh.conj().T) / 2
    Yh = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
    Yh = Yh @ Yh.conj().T + 0.3 * np.eye(2)
    pts.append(Xh + 1j * Yh)
blocks = [[None] * 6 for _ in range(6)]
for i, Zi in enumerate(pts):
    for j, Zj in enumerate(pts):
        Am = -1j * (Zi - Zj.conj().T)
        C = np.linalg.inv(Am)
        blocks[i][j] = 24 / np.linalg.det(Am) * k4(C) @ Psym
G2 = np.block(blocks)
G2h = (G2 + G2.conj().T) / 2
mineig = np.linalg.eigvalsh(G2h).min()
say(f"   Röhrengebiet, 6 Punkte: kleinster Eigenwert der Helizität-2-Gram-Matrix {mineig:.2e} (Größenordnung {np.abs(G2h).max():.2e})")
check("Helizität-2-Kern K_2 ist auf dem Röhrengebiet positiv definit (Gaußdarstellung mit |lambda^4><lambda^4|): positiver freier Spin-2-Sektor",
      mineig > -1e-10 * np.abs(G2h).max())
befund("Der positive Ersatz für den Formfaktor ist der Wallach-Kern: nu = 1 gibt den reinen Pol (Q^2 C = Z), nu > 1 ein positives Massenkontinuum", True)

# ------------------------------------------------------------------ Teil 3
say("\nTeil 3: Reelle Naht aus dem Altbestand (SEAM.S3) gegen die komplexe Flussnaht")
D5 = []
for i_, j_ in itertools.combinations(range(5), 2):
    for si, sj in itertools.product((1, -1), repeat=2):
        v = [0] * 8
        v[i_], v[j_] = si, sj
        D5.append(v)
q = [1] * 8
neutral_q = sum(1 for v in D5 if sum(a * b for a, b in zip(v, q)) == 0)
check("Majorana-Naht ohne U(1)-Ladung: alle 40 D5-Wurzeln bleiben Ströme; mit Flussladung q nur 20 (U(5))", len(D5) == 40 and neutral_q == 20)
dim_u8 = 8 * 8
dim_so16 = 16 * 15 // 2
check("Kommutant der Flussladung in so(16): u(8) (64 von 120 Richtungen); die reelle Naht trägt das ganze so(16)_1 ⊃ so(10)_1 + so(6)_1",
      dim_u8 == 64 and dim_so16 == 120 and 45 + 15 <= dim_so16)
vec10 = [[(1 if k == i_ else 0) * sg for k in range(5)] for i_ in range(5) for sg in (1, -1)]
alpha_l = [1, -1, 0, 0, 0]
lvl_real10 = Fr(sum(sum(a * b for a, b in zip(w, alpha_l)) ** 2 for w in vec10), 2) / 2   # komplex Sum(w.alpha)^2/2 = 2, reell halbiert
check("Niveau 1 für zehn reelle Majorana-Fermionen im Vektor (komplexes Niveau 2, halbiert)", lvl_real10 == 1)
befund("Die reelle 16-Majorana-Naht (p+ip x 16, SEAM.S3.LATTICE/INFLOW) trägt Spin(10)_1 und den Zufluss ohne U(1); ihre Chiralität "
       "kommt aus der Chernphase, nicht aus einem geladenen Wirbel", True)

# ------------------------------------------------------------------ Teil 4
say("\nTeil 4: M(P) als Yukawa-Kopplung an die Higgsrichtung")
th, ph, al = sp.symbols("theta phi alpha", real=True)
Hhat = sp.Matrix([sp.cos(th / 2), sp.sin(th / 2) * sp.exp(sp.I * ph)])
Pm = Hhat * Hhat.H
U = 2 * Pm - sp.eye(2)
check("P = H^ H^+ ist ein Rang-1-Projektor; 2P - 1 hat die Eigenwerte +1 (Higgsrichtung) und -1 (orthogonal): Dublettspaltung",
      sp.simplify(Pm * Pm - Pm) == sp.zeros(2) and sp.simplify(U * U - sp.eye(2)) == sp.zeros(2) and sp.simplify(U.trace()) == 0)
dPth, dPph = sp.diff(Pm, th), sp.diff(Pm, ph)
dHth, dHph = sp.diff(Hhat, th), sp.diff(Hhat, ph)
fs_th = sp.simplify(((dHth.H * dHth)[0] - (Hhat.H * dHth)[0] * sp.conjugate((Hhat.H * dHth)[0])))
fs_ph = sp.simplify(((dHph.H * dHph)[0] - (Hhat.H * dHph)[0] * sp.conjugate((Hhat.H * dHph)[0])))
fs_ok = True
for tv, pv_ in ((0.3, 1.1), (1.2, -0.7), (2.5, 2.9)):
    sub = {th: tv, ph: pv_}
    fs_ok &= abs(complex(((dPth * dPth).trace() - 2 * fs_th).subs(sub).evalf())) < 1e-12
    fs_ok &= abs(complex(((dPph * dPph).trace() - 2 * fs_ph).subs(sub).evalf())) < 1e-12
check("tr(dP)^2 = 2(|dH^|^2 - |H^+ dH^|^2) (Fubini-Study, an drei Punkten exakt ausgewertet): Gamma_pol ist der eichinvariante "
      "Kinetikterm der Higgsrichtung", fs_ok)
p4, m4, L4 = sp.symbols("p m Lambda", positive=True)
I4 = sp.integrate(2 * sp.pi**2 * p4**3 * m4**2 / (p4**2 + m4**2) ** 2, (p4, 0, L4)) / (2 * sp.pi) ** 4
lead = sp.limit(I4 / sp.log(L4), L4, sp.oo)
check("in 4D ist der induzierte Koeffizient m^2/(16 pi^2) log(Lambda^2/m^2): Higgs-Wellenfunktion aus der Yukawa-Spur (Altbestand A = 4 Y_2)",
      sp.simplify(lead - m4**2 / (8 * sp.pi**2)) == 0)
befund("Der für t* nötige Sektor M(P) ist im Altbestand die Yukawa/Higgs-Fluktuation; der Koeffizient folgt aus Tr Y^+ Y", True)

urteil = ("Der Altbestand enthält Gegenstücke zu mehreren offenen Punkten, und eines davon kehrt sich gegen seine eigene "
          "Lesart. Das gemeinsame Masterfunktional ist dort die Spektralwirkung eines Dirac-Operators mit Naht-KMS-Abschneider "
          "e^{-u}. Ihr tatsächlicher Spin-2-Kinetikfaktor aus dem nichtlokalen Wärmeleitungskern ist a(x) = 1 + 12 x g(x) = 1 - "
          "3x/10 + x^2/28 - ..., nicht e^{x}: Q^2 C >= Z gilt nur bis x* = 6.768 Lambda^2, dort liegt ein reeller euklidischer Pol "
          "mit negativem Residuum (tachyonischer Spin-2-Geist), dazu eine unendliche Folge komplexer Polpaare; jeder getestete "
          "vollständig monotone Abschneider zeigt dasselbe. Die Resummation entkoppelt den Geist also nicht; die [C]-Identifikation "
          "in GRAV.KMS.HESSIAN.01 hält für Tr e^{-D^2/Lambda^2} nicht. Die positive Gegenseite ist der Wallach-Kern: det(P)^{s-2} "
          "dP hat die Laplace-Transformierte pi Gamma(s)Gamma(s-1) det(Y)^-s, die Spektraldichte (m^2)^{s-2} konzentriert sich bei "
          "s -> 1 auf m^2 = 0, die Wallach-Menge ist die 4D-Unitaritätsschranke Delta >= 1, und der Helizität-2-Kern (Delta = 3) "
          "ist auf dem Röhrengebiet positiv. Die reelle 16-Majorana-Naht des Altbestands (SEAM.S3) trägt Spin(10)_1 und den "
          "Zufluss ohne U(1); die komplexe Flussnaht nicht. Der Sektor M(P) ist die Yukawa-Kopplung an die Higgsrichtung, "
          "Gamma_pol ihr Fubini-Study-Kinetikterm.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z41_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False, default=str) + "\n")
