"""Z72: Zeittranslation aus der periodischen Uhr? Parabolischer Generator, halbseitige modulare Inklusionen, Originalzustand.

Frage (VERDICT H6h34, offen (3)): Die globale Uhr K = 4 pi L_0 ist periodisch (z70), die lokale Modulzeit auf Nahtbögen
kontinuierlich (z71). Es fehlen lokale Translationsbeziehungen. Liefert sie (a) der parabolische Generator P_0 im
Vakuumsektor und/oder (b) eine halbseitige modulare Inklusion (hsm, Wiesbrock/Borchers) im Originalzustand (beta = 4 pi)?

Konventionen (explizit)
  Naht S^1, theta in [0, 2 pi); Moden e_x = e^{i x theta}/sqrt(2 pi), x in Z + alpha (NS-Vakuum alpha = 1/2, Original 1/4, 3/4).
  Einteilchengenerator eines Vektorfelds v d_theta auf Halbdichten: G_v = -i (v d_theta + v'/2); e^{-isG_v} schiebt
  Wellenpakete längs v vorwärts. L_n := G_{e^{-in theta}}: L_0 = -i d_theta, L_n e_x = (x - n/2) e_{x-n},
  [L_m, L_n] = (m - n) L_{m+n}, L_n^+ = L_{-n}.
  Diskrete Serie D_h: |n> = |h+n>, L_0|n> = (h+n)|n>, L_{-1}|n> = sqrt((n+1)(n+2h)) |n+1>, L_1 = L_{-1}^+.
  U(1)-Strom: h = 1 (J_{-m} Omega); Fermion (NS-Teilchen x = 1/2, 3/2, ...): h = 1/2.
  Reelle Erzeuger: elliptisch L_0; hyperbolisch H_1 = (L_1 + L_{-1})/2 = G_{cos}, H_2 = i(L_1 - L_{-1})/2 = G_{sin};
  parabolisch P_0 = L_0 + H_1 = G_{1 + cos} (Fixpunkt theta = pi), K_0 = L_0 - H_1 = G_{1 - cos} (Fixpunkt 0).
  Cayley y = tan(theta/2): P_0 = -i d_y (Translation), H_2 = Dilatation y d_y, K_0 = y^2 d_y, L_0 = (P_0 + K_0)/2.
  Allgemein P^(t*) = G_{1 - cos(theta - t*)} = L_0 - (e^{i t*} L_1 + e^{-i t*} L_{-1})/2, D^(t*) = G_{sin(theta - t*)}.
  Modular: Delta^{it} = e^{-itK}, K = -log Delta; quasifrei K = dGamma(k). Peschel: h_A = log((1 - rho_A)/rho_A) mit dem
  Einteilchen-Dichteoperator rho (<c^+(v) c(w)> = <w, rho v>); reiner Zustand: k_A = h_A (+) (-h_A'); gemischter Zustand:
  Purifizierung der gemischten Moden. dGamma(X) mit [X, Pi] = 0 ist (normalgeordnet) positiv genau dann, wenn X(1 - 2 Pi) >= 0.

Teil 1  sl(2,R) auf D_h: L_0 diskret {h, h+1, ...}; P_0 exakt Laguerre (Jacobi-Matrix), also absolut stetig [0, inf) mit
        Gamma-Spektralmass; H_2 absolut stetig auf R; [H_2, P_0] = i P_0 (Dilatation), keine eigene Skala.
Teil 2  Gitter-Chiralfermion (M Plätze, M Moden, SLAC-artig, ohne Verdopplung), Vakuum: Hislop-Longo, Wiesbrock-Differenz
        = parabolischer Generator, Borchers-Relation, affine Translationsfamilie, Positivität, Gamma-Gesetz; Kontrolle:
        streng geschachtelte Bögen geben einen elliptischen (getakteten) Generator.
Teil 3  Originalzustand beta = 4 pi: global periodisch (Takesaki: keine hsm), Bögen: Lie-Defekt, Translationsfamilie,
        Gamma-Gesetz; hochgenau (python-flint, 512 Bit) die hsm-Definition direkt (Leckage unter dem Modulfluss).
Teil 4  Urteil mit Dimensionstypisierung (GRAV.NONCIRCULAR.01, DIMENSION.UPLIFT.FIREWALL.01).
"""
import json
import time
from pathlib import Path

import numpy as np
import sympy as sp
from flint import acb, acb_mat, arb, ctx
from scipy.special import eval_genlaguerre, gammaln, jn_zeros, roots_genlaguerre

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
T0 = time.time()


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val):
    BEFUNDE[name] = bool(val)
    say(f"   [Befund: {'ja' if val else 'nein'}] {name}")


def r(v, n=4):
    return float(f"{float(v):.{n}g}")


# ================================================================== Teil 1
say("Teil 1: sl(2,R) auf der diskreten Serie D_h -- elliptisch, parabolisch, hyperbolisch")


def ds_sym(h, N):
    L0 = sp.diag(*[h + n for n in range(N)])
    Lm1 = sp.zeros(N, N)
    for n in range(N - 1):
        Lm1[n + 1, n] = sp.sqrt((n + 1)*(n + 2*h))
    return L0, Lm1.T, Lm1


ok_alg, ok_cas = True, True
for h in (sp.Rational(1, 2), sp.Integer(1)):
    N = 10
    L0, L1, Lm1 = ds_sym(h, N)
    inner = list(range(N - 1))
    c1 = (L1*Lm1 - Lm1*L1 - 2*L0)[inner, inner]
    c2 = (L0*L1 - L1*L0 + L1)
    c3 = (L0*Lm1 - Lm1*L0 - Lm1)
    ok_alg &= c1.applyfunc(sp.simplify).is_zero_matrix and c2.is_zero_matrix and c3.is_zero_matrix
    cas = (L0*L0 - (L1*Lm1 + Lm1*L1)/2)[inner, inner].applyfunc(sp.simplify)
    ok_cas &= cas == h*(h - 1)*sp.eye(N - 1)
check("Konventionen: auf D_h (h = 1/2, 1; exakt, sympy) gelten [L_1, L_-1] = 2 L_0, [L_0, L_+-1] = -+L_+-1 und die Casimir-"
      "Identität L_0^2 - (L_1 L_-1 + L_-1 L_1)/2 = h(h - 1); L_0 ist diagonal mit Spektrum {h, h+1, ...}, für den Strom "
      "(h = 1) also {1, 2, 3, ...}", ok_alg and ok_cas)

# Strukturkonstanten aus den Matrizen und Killing-Form der reellen Form (E_k = i X_k)
h = sp.Integer(1)
N = 8
L0, L1, Lm1 = ds_sym(h, N)
H1 = (L1 + Lm1)/2
H2 = sp.I*(L1 - Lm1)/2
X = [L0, H1, H2]
inner = list(range(N - 2))
coef = np.zeros((3, 3, 3), dtype=complex)
for a in range(3):
    for b in range(3):
        C = (X[a]*X[b] - X[b]*X[a])[inner, inner]
        A = sp.Matrix([[X[k][inner, inner][i, j] for k in range(3)] for i in range(len(inner)) for j in range(len(inner))])
        sol = A.solve_least_squares(sp.Matrix([C[i, j] for i in range(len(inner)) for j in range(len(inner))]))
        coef[a, b] = [complex(sp.N(s_)) for s_ in sol]
# reelle Basis E_k = i X_k: [E_a, E_b] = -[X_a, X_b] = sum_c (-coef_abc) X_c = sum_c (i coef_abc) E_c
f_ = 1j*coef
ad = [np.array([[f_[a, b, c] for b in range(3)] for c in range(3)]) for a in range(3)]


def killing(v):
    A = sum(v[k]*ad[k] for k in range(3))
    return float(np.real(np.trace(A @ A)))


kap = {"L0": killing([1, 0, 0]), "P0 = L0 + H1": killing([1, 1, 0]), "H2": killing([0, 0, 1])}
ZAHLEN["killing_form"] = {k_: r(v) for k_, v in kap.items()}
check(f"Klassifikation aus den gerechneten Strukturkonstanten ([L_0, H_1] = i H_2, [L_0, H_2] = -i H_1, [H_1, H_2] = -i L_0): "
      f"Killing-Form B(X, X) = {ZAHLEN['killing_form']}, also L_0 elliptisch (< 0), P_0 parabolisch (= 0), H_2 hyperbolisch "
      f"(> 0); allgemein B = 2(-a^2 + b^2 + c^2) für a L_0 + b H_1 + c H_2",
      np.allclose(np.imag(f_), 0, atol=1e-12) and kap["L0"] < -1 and abs(kap["P0 = L0 + H1"]) < 1e-10 and kap["H2"] > 1)

# exakte Spektraltheorie von P_0: 2 P_0 ist die Jacobi-Matrix der Laguerre-Gewichtsfunktion x^(2h-1) e^-x
n_, h_ = sp.symbols("n h", nonnegative=True)
a_ = 2*h_ - 1
diag_ok = sp.simplify(2*(h_ + n_) - (2*n_ + a_ + 1)) == 0
off_ok = sp.simplify((n_ + 1)*(n_ + 2*h_) - (n_ + 1)*(n_ + a_ + 1)) == 0
rec_err = 0.0
for hh in (0.5, 1.0, 1.5):
    a = 2*hh - 1
    for x in (0.3, 1.7, 4.2, 9.0):
        ell = [np.exp(0.5*(gammaln(k + 1) - gammaln(k + a + 1)))*eval_genlaguerre(k, a, x) for k in range(18)]
        psi = [(-1)**k*ell[k] for k in range(18)]
        for k in range(1, 16):
            lhs = x*psi[k]
            rhs = (2*k + a + 1)*psi[k] + np.sqrt((k + 1)*(k + a + 1))*psi[k + 1] + np.sqrt(k*(k + a))*psi[k - 1]
            rec_err = max(rec_err, abs(lhs - rhs)/(1 + abs(lhs)))
mom_ok = True
for hh in (sp.Rational(1, 2), sp.Integer(1), sp.Rational(3, 2)):
    N = 7
    L0s, L1s, Lm1s = ds_sym(hh, N)
    P0s = L0s + (L1s + Lm1s)/2
    Pk = sp.eye(N)
    for k in range(1, 9):
        Pk = Pk*P0s
        mom_ok &= sp.nsimplify(sp.simplify(Pk[0, 0])) == sp.rf(2*hh, k)/2**k
check("exakt: 2 P_0 auf D_h hat Diagonale 2n + a + 1 und Nebendiagonale sqrt((n+1)(n+a+1)) mit a = 2h - 1 (symbolisch), "
      "psi_n(s) = (-1)^n l_n^(a)(2s) (orthonormierte Laguerre-Funktionen) erfüllt P_0 psi = s psi (Rekursion, Fehler "
      f"{rec_err:.1e}), und <h|P_0^k|h> = (2h)_k / 2^k für k = 1..8, h = 1/2, 1, 3/2 (sympy): P_0 ist unitär äquivalent zur "
      "Multiplikation mit s auf L^2((0, inf), (2s)^(2h-1) e^(-2s) ds) -- Spektrum [0, inf), absolut stetig, einfach, "
      "Spektralmass des Grundzustands = Gamma(2h, 1/2)",
      diag_ok and off_ok and rec_err < 1e-10 and mom_ok)


def ds_num(hh, N):
    n = np.arange(N)
    L0 = np.diag(hh + n).astype(float)
    Lm1 = np.diag(np.sqrt((n[:-1] + 1)*(n[:-1] + 2*hh)), -1)
    return L0, Lm1.T.copy(), Lm1


trunc = {}
for hh in (1.0, 0.5):
    a = 2*hh - 1
    for N in (25, 100, 400):
        L0n, L1n, Lm1n = ds_num(hh, N)
        P0n = L0n + (L1n + Lm1n)/2
        ev = np.linalg.eigvalsh(P0n)
        zerr = 0.0
        if N <= 100:                                   # scipy liefert für N = 400 keine endlichen größten Nullstellen
            zl = np.sort(roots_genlaguerre(N, a)[0])/2
            zerr = float(np.max(np.abs(ev - zl)/zl))
        trunc[(hh, N)] = (zerr, float(ev[0]), int(np.sum(ev < 1)), float(N*ev[0]),
                          float(np.linalg.eigvalsh(L0n)[0]), int(np.sum(np.diag(L0n) < hh + 0.5)))
jz = jn_zeros(1, 1)[0]
ZAHLEN["trunkierung_h1"] = {N: {"min_eig_P0": r(trunc[(1.0, N)][1], 5), "anzahl_eig_P0_unter_1": trunc[(1.0, N)][2],
                                "N*min": r(trunc[(1.0, N)][3], 5), "min_eig_L0": trunc[(1.0, N)][4]} for N in (25, 100, 400)}
ZAHLEN["j_11^2/8"] = r(jz**2/8, 5)
say(f"   Trunkierung h = 1: {ZAHLEN['trunkierung_h1']}")
t1 = trunc[(1.0, 25)], trunc[(1.0, 100)], trunc[(1.0, 400)]
check("Trunkierung (N = 25, 100, 400): die Eigenwerte der N x N-Trunkierung von P_0 sind genau die Laguerre-Nullstellen/2 "
      "(Gauss-Quadratur, rel. Fehler < 1e-9 bei N = 25, 100); der kleinste fällt wie 1/N gegen 0 (N s_min -> j_{1,1}^2/8 = 1,8353 für h = 1), "
      "die Zahl der Eigenwerte unter 1 wächst wie sqrt(N) (Verdichtung); L_0 behält das Minimum h = 1 und genau einen "
      "Eigenwert unter 1,5",
      all(trunc[k][0] < 1e-9 for k in trunc) and t1[0][1] > t1[1][1] > t1[2][1] and abs(t1[2][3]/(jz**2/8) - 1) < 0.01
      and t1[0][2] < t1[1][2] < t1[2][2] and 1.6 < t1[2][2]/t1[1][2] < 2.4
      and all(trunc[(1.0, N)][4] == 1.0 and trunc[(1.0, N)][5] == 1 for N in (25, 100, 400)))


def char(Xm, t, v0=None):
    w, V = np.linalg.eigh(Xm)
    v = np.zeros(Xm.shape[0]) if v0 is None else v0
    if v0 is None:
        v[0] = 1.0
    amp = V.conj().T @ v
    return complex(np.sum(np.abs(amp)**2*np.exp(1j*t*w)))


N = 600
errs = []
for hh in (1.0, 0.5):
    L0n, L1n, Lm1n = ds_num(hh, N)
    P0n = L0n + (L1n + Lm1n)/2
    H2n = 1j*(L1n - Lm1n)/2
    for t in (0.5, 1.0, 2.0, 4.0):
        errs.append(abs(char(P0n, t) - (1 - 1j*t/2)**(-2*hh)))
        errs.append(abs(char(H2n, t) - np.cosh(t/2)**(-2*hh)))
        errs.append(abs(char(L0n, t) - np.exp(1j*hh*t)))
L0n, L1n, Lm1n = ds_num(1.0, N)
P0n = L0n + (L1n + Lm1n)/2
decayP = abs(char(P0n, 12.0))
returnL = abs(char(L0n, 2*np.pi) - 1)
ZAHLEN["charakteristische_funktionen"] = {"max_fehler": r(max(errs), 3), "|<e^{12 i P0}>| (h=1)": r(decayP, 4),
                                          "|<e^{2 pi i L0}> - 1|": r(returnL, 3)}
check("Spektralmasse des Grundzustands (N = 600, t = 0,5..4): <h|e^{itP_0}|h> = (1 - it/2)^{-2h} (Gamma, absolut stetig "
      "auf (0, inf)), <h|e^{itH_2}|h> = cosh(t/2)^{-2h} (absolut stetig auf R, Meixner-Pollaczek), <h|e^{itL_0}|h> = e^{iht} "
      f"(Punktmass); max. Fehler {max(errs):.1e}. Riemann-Lebesgue: |<e^{{12iP_0}}>| = {decayP:.4f} = (1 + 36)^-1 zerfällt, "
      "die Uhr kehrt bei t = 2 pi exakt zurück",
      max(errs) < 1e-8 and abs(decayP - 1/37) < 1e-6 and returnL < 1e-12)

# Dilatation: [H_2, P_0] = i P_0 exakt, e^{isH_2} P_0 e^{-isH_2} = e^{-s} P_0
L0s, L1s, Lm1s = ds_sym(sp.Integer(1), 9)
P0s = L0s + (L1s + Lm1s)/2
K0s = L0s - (L1s + Lm1s)/2
H2s = sp.I*(L1s - Lm1s)/2
inn = list(range(7))
dil_ok = ((H2s*P0s - P0s*H2s - sp.I*P0s)[inn, inn]).applyfunc(sp.simplify).is_zero_matrix
conf_ok = (L0s - (P0s + K0s)/2).is_zero_matrix
H2n = 1j*(L1n - Lm1n)/2
wH, VH = np.linalg.eigh(H2n)
derr = 0.0
for s_, t in ((0.7, 1.3), (-0.5, 2.0), (1.2, 0.8)):
    e0 = np.zeros(N)
    e0[0] = 1.0
    w0 = VH @ (np.exp(-1j*s_*wH)*(VH.conj().T @ e0))
    derr = max(derr, abs(char(P0n, t, w0) - (1 - 1j*t*np.exp(-s_)/2)**(-2)))
check("Dilatation: [H_2, P_0] = i P_0 exakt (sympy), also e^{isH_2} P_0 e^{-isH_2} = e^{-s} P_0; numerisch "
      f"<h|e^{{isH_2}} e^{{itP_0}} e^{{-isH_2}}|h> = (1 - i t e^-s/2)^-2 (Fehler {derr:.1e}); dazu L_0 = (P_0 + K_0)/2 exakt",
      dil_ok and conf_ok and derr < 1e-8)
befund("Teil 1 (Typ: eindimensional, eine Lichtstrahlkoordinate): Die Uhr L_0 ist elliptisch (Spektrum {1, 2, ...} auf dem "
       "Strom), P_0 parabolisch (absolut stetig [0, inf), positiv), H_2 hyperbolisch (absolut stetig R, erzeugt "
       "Dilatationen, die P_0 skalieren). Weil lambda P_0 zu P_0 unitär äquivalent ist, trägt die Translation keine eigene "
       "Skala; die Einheit kommt erst mit der Uhr L_0 = (P_0 + K_0)/2, also mit dem Kreis (beta = 4 pi). Der Bestand v438 "
       "(SEAM.EQUIV.BW.HSMI.01) prüfte nur die Positivität einer Trunkierung; die ist für Kompressionen positiver Operatoren "
       "automatisch und zeigt kein Kontinuum bei 0. Hier steht die exakte Spektraltheorie", True)
ZAHLEN["laufzeit_teil1_s"] = r(time.time() - T0, 3)

# ================================================================== Teil 2
say("\nTeil 2: Gitter-Chiralfermion im Vakuum -- Peschel-Modulgeneratoren, Wiesbrock, Borchers")


def lattice_state(M, alpha, beta, eps=1e-14):
    xs = np.arange(M) - M//2 + alpha
    th = 2*np.pi*np.arange(M)/M
    U = np.exp(1j*np.outer(th, xs))/np.sqrt(M)
    if beta is None:
        c = (xs < 0).astype(float)
        mixed = np.array([], dtype=int)
    else:
        c = 1/(1 + np.exp(np.clip(beta*xs, -700, 700)))
        mixed = np.where((c > eps) & (c < 1 - eps))[0]
        c = np.where(c <= eps, 0.0, np.where(c >= 1 - eps, 1.0, c))
    m = len(mixed)
    n = M + m
    G = np.zeros((n, n), complex)
    G[:M, :M] = (U*c) @ U.conj().T
    if m:
        cm = c[mixed]
        G[:M, M:] = U[:, mixed]*np.sqrt(cm*(1 - cm))
        G[M:, :M] = G[:M, M:].conj().T
        G[M:, M:] = np.diag(1 - cm)
    W = np.zeros((n, n), complex)
    W[:M, :M] = U
    W[M:, M:] = np.eye(m)
    return dict(M=M, xs=xs, th=th, U=U, c=c, G=G, W=W, labels=np.concatenate([xs, xs[mixed]]), m=m, alpha=alpha)


def hfun(nu, hmax):
    nu = np.clip(nu, 0, 1)
    with np.errstate(divide="ignore"):
        hv = np.log(1 - nu) - np.log(nu)
    return np.clip(hv, -hmax, hmax)


def K_mode(st, sites, hmax=34.0):
    """voller Einteilchen-Modulgenerator k_A = h_A (+) (-h_{A^c}) (A^c inkl. Purifizierer), in der Modenbasis"""
    G = st["G"]
    n = G.shape[0]
    A = np.unique(np.asarray(sites) % st["M"])
    Ac = np.setdiff1d(np.arange(n), A)
    k = np.zeros((n, n), complex)
    for R_, sg in ((A, 1.0), (Ac, -1.0)):
        w, V = np.linalg.eigh(G[np.ix_(R_, R_)])
        k[np.ix_(R_, R_)] = sg*(V*hfun(w, hmax)) @ V.conj().T
    return st["W"].conj().T @ k @ st["W"]


def mob(xs):
    M = len(xs)
    L0 = np.diag(xs).astype(complex)
    L1 = np.zeros((M, M), complex)
    Lm1 = np.zeros((M, M), complex)
    for k in range(M):
        if k >= 1:
            L1[k - 1, k] = xs[k] - 0.5
        if k + 1 < M:
            Lm1[k + 1, k] = xs[k] + 0.5
    return L0, L1, Lm1


def pad(Aop, n):
    out = np.zeros((n, n), complex)
    out[:Aop.shape[0], :Aop.shape[1]] = Aop
    return out


def arcs_sharing(M, js, fracs):
    """Bögen [js - f M, js) mit gemeinsamem rechtem Endpunkt t* = 2 pi (js - 1/2)/M; Cayley-Koordinate a = cot(delta/2)"""
    tstar = 2*np.pi*(js - 0.5)/M
    out = []
    for f in fracs:
        ja = js - int(round(f*M))
        delta = tstar - 2*np.pi*(ja - 0.5)/M
        out.append((np.arange(ja, js), 1/np.tan(delta/2)))
    return tstar, out


def analyse(st, js, fracs, nlo=2.0, hmax=34.0):
    M = st["M"]
    n = len(st["labels"])
    tstar, arcs = arcs_sharing(M, js, fracs)
    Ks = [K_mode(st, s_, hmax) for s_, _ in arcs]
    a = np.array([a_ for _, a_ in arcs])
    L0, L1, Lm1 = mob(st["xs"] - (st["alpha"] - 0.5))     # Referenz in NS-Etiketten (Eichung W = e^{-i(1/2-alpha) theta})
    Pst = pad(L0 - (np.exp(1j*tstar)*L1 + np.exp(-1j*tstar)*Lm1)/2, n)
    Dst = pad(1j*(np.exp(1j*tstar)*L1 - np.exp(-1j*tstar)*Lm1)/2, n)
    lo = np.where(np.abs(st["labels"]) < nlo)[0]
    lo = lo[lo < M]                                     # Systemblock (untere Schranke für Verletzungen)
    b = lambda Xm: Xm[np.ix_(lo, lo)]
    res = {}
    E = [2*np.pi*(-a_*Pst - Dst) for a_ in a]          # Hislop-Longo: K_A = 2 pi G_{v_A}
    res["hl"] = [float(np.abs(b(K - Ee)).max()/np.abs(b(Ee)).max()) for K, Ee in zip(Ks, E)]
    KJ, KI = Ks[0], Ks[1]
    fac = a[1] - a[0]
    PW = (KJ - KI)/(2*np.pi)
    res["fac"] = float(fac)
    res["wiesbrock"] = float(np.abs(b(PW - fac*Pst)).max()/(fac*np.abs(b(Pst)).max()))
    comm = KJ @ KI - KI @ KJ
    scale = np.abs(b(2*np.pi*(KJ - KI))).max()
    res["lie_pred"] = float(np.abs(b(comm + 2j*np.pi*(KJ - KI))).max()/scale)     # vorhergesagt: [K_J,K_I] = -2 pi i (K_J - K_I)
    res["lie_other"] = float(np.abs(b(comm - 2j*np.pi*(KJ - KI))).max()/scale)
    # affine Translationsfamilie K_A = X - 2 pi a_A Y (Kleinste Quadrate über alle Bögen)
    Bm = np.array([b(K).ravel() for K in Ks])
    Dm = np.vstack([np.ones(len(a)), -2*np.pi*a]).T
    coef, *_ = np.linalg.lstsq(Dm, Bm, rcond=None)
    resid = Bm - Dm @ coef
    res["familie"] = float(np.abs(resid).max()/np.abs(Bm - Bm.mean(0)).max())
    Y = coef[1].reshape(len(lo), len(lo))
    res["steigung_vs_P"] = float(np.abs(Y - b(Pst)).max()/np.abs(b(Pst)).max())
    # Fock-Positivität und Gamma-Gesetz
    Gm = st["W"].conj().T @ st["G"] @ st["W"]
    Xt = PW @ (np.eye(n) - 2*Gm)
    Xt = (Xt + Xt.conj().T)/2
    lc = np.where(np.abs(st["labels"]) < 3)[0]
    res["poscomp"] = float(np.linalg.eigvalsh(Xt[np.ix_(lc, lc)])[0])
    ev, V = np.linalg.eigh(Xt)
    xs = st["xs"]
    ip = int(np.argmin(np.where(xs > 0, xs, np.inf)))
    im = int(np.argmax(np.where(xs < 0, xs, -np.inf)))
    lap = {}
    for t in (1.0, 4.0, 16.0):
        Lp = float(np.sum(np.abs(V[ip])**2*np.exp(-t*ev)))
        Lm = float(np.sum(np.abs(V[im])**2*np.exp(-t*ev)))
        lap[t] = (Lp*(1 + t*fac/2), Lm*(1 + t*fac/2), Lp*Lm*(1 + t*fac/2)**2)
    res["gamma"] = lap
    return res


FR = (0.5, 0.4375, 0.375, 0.3125)
vac = {}
for M in (512, 1024):
    vac[M] = analyse(lattice_state(M, 0.5, None), M//2, FR)
v5 = vac[512]
ZAHLEN["vakuum_M512"] = {"hislop_longo_rel": [r(x, 3) for x in v5["hl"]], "wiesbrock_rel": r(v5["wiesbrock"], 3),
                         "fac=cot(dI/2)-cot(dJ/2)": r(v5["fac"], 5), "lie_vorhergesagt": r(v5["lie_pred"], 3),
                         "lie_andere_Orientierung": r(v5["lie_other"], 3), "familie_resid": r(v5["familie"], 3),
                         "steigung_vs_P": r(v5["steigung_vs_P"], 3), "fock_min_kompression": r(v5["poscomp"], 3),
                         "gamma_ratio": {t: [r(x, 6) for x in v] for t, v in v5["gamma"].items()}}
ZAHLEN["vakuum_M1024"] = {"hislop_longo_rel": [r(x, 3) for x in vac[1024]["hl"]], "lie_vorhergesagt": r(vac[1024]["lie_pred"], 3),
                          "familie_resid": r(vac[1024]["familie"], 3)}
say(f"   Vakuum M = 512: {ZAHLEN['vakuum_M512']}")
check("Hislop-Longo als Resultat aus dem Zustand: die Peschel-Modulgeneratoren von vier Bögen mit gemeinsamem Endpunkt t* "
      "(Längen pi, 7pi/8, 3pi/4, 5pi/8) sind auf den niedrigen Moden 2 pi G_{v_A}, v_A = Möbiusfeld mit Nullstellen in den "
      f"Endpunkten (rel. Abweichung <= {max(v5['hl']):.1e} bei M = 512)", max(v5["hl"]) < 5e-3)
check("Wiesbrock: (K_J - K_I)/2 pi = (cot(delta_I/2) - cot(delta_J/2)) P^(t*) mit dem parabolischen Generator "
      f"P^(t*) = G_{{1 - cos(theta - t*)}} (rel. Abweichung {v5['wiesbrock']:.1e}; Faktor {v5['fac']:.4f})",
      v5["wiesbrock"] < 5e-3)
check("Borchers/Wiesbrock-Algebra: [K_J, K_I] = -2 pi i (K_J - K_I) (vorhergesagte Orientierung für den gemeinsamen rechten "
      f"Endpunkt) mit rel. Defekt {v5['lie_pred']:.1e}; die andere Orientierung verfehlt ({v5['lie_other']:.2f}): die zwei "
      "Modulgeneratoren spannen die ax+b-Algebra", v5["lie_pred"] < 5e-3 and v5["lie_other"] > 0.3)
check("Translationsverhalten: die vier Modulgeneratoren liegen auf einer affinen Geraden K_A = X - 2 pi a_A Y in der Cayley-"
      f"Koordinate a_A = cot(delta_A/2) des freien Endpunkts (Residuum {v5['familie']:.1e}), und die Steigung ist der "
      f"parabolische Generator, Y = P^(t*) (rel. {v5['steigung_vs_P']:.1e}): e^{{-isP}} verschiebt den freien Endpunkt um s",
      v5["familie"] < 1e-2 and v5["steigung_vs_P"] < 1e-2)
g = v5["gamma"]
gam_dev = max(abs(x - 1) for v in g.values() for x in v)
check("Positivität und Gamma-Gesetz: die Fock-Anregung P_W(1 - 2 Pi) ist auf der Kompression (|x| < 3) positiv "
      f"(min {v5['poscomp']:.3f}); die Laplace-Transformierten der niedrigsten Teilchen- und Lochmode sind "
      "(1 + t fac/2)^-1 (h = 1/2, Teil 1) und ihr Produkt, der Teilchen-Loch-Zustand J_-1 Omega des Stroms, ist "
      f"(1 + t fac/2)^-2 (h = 1, Teil 1); max. rel. Abweichung {gam_dev:.1e} bei t = 1, 4, 16",
      v5["poscomp"] > -1e-3 and gam_dev < 1e-3)
check("Stabilität in M: bei M = 1024 bleiben Hislop-Longo-, Lie- und Familiendefekt auf dem Niveau von M = 512 (Grenze: "
      "Doppelgenauigkeit der Eigenwerte nahe 0 und 1, nicht die Gitterweite)",
      max(vac[1024]["hl"]) < 5e-3 and vac[1024]["lie_pred"] < 5e-3 and vac[1024]["familie"] < 1e-2)

# verdrehtes Vakuum alpha = 1/4 ist eichäquivalent (lokale Multiplikation e^{-i theta/4})
st14 = lattice_state(512, 0.25, None)
stNS = lattice_state(512, 0.5, None)
Wg = np.diag(np.exp(-1j*stNS["th"]/4))
gauge = float(np.abs(st14["G"] - Wg @ stNS["G"] @ Wg.conj().T).max())
check(f"verdrehtes Vakuum (alpha = 1/4): rho_alpha = W rho_NS W^+ mit der lokalen Eichung W = e^{{-i theta/4}} (Fehler "
      f"{gauge:.1e}); also sind alle Modulgeneratoren konjugiert und die hsm-Struktur gilt mit P_alpha = W P^(t*) W^+",
      gauge < 1e-12)

# Kontrolle: streng geschachtelte Bögen (kein gemeinsamer Endpunkt) geben einen elliptischen Generator
M = 512
stv = lattice_state(M, 0.5, None)
J_c, I_c = np.arange(M//4, 3*M//4), np.arange(3*M//8, 5*M//8)
KJc, KIc = K_mode(stv, J_c), K_mode(stv, I_c)
PWc = (KJc - KIc)/(2*np.pi)
th_ = np.linspace(0, 2*np.pi, 4096, endpoint=False)


def vfield(sites):
    ta = 2*np.pi*(sites[0] - 0.5)/M
    tb = 2*np.pi*(sites[-1] + 0.5)/M
    return 2*np.sin((th_ - ta)/2)*np.sin((tb - th_)/2)/np.sin((tb - ta)/2)


dv = vfield(J_c) - vfield(I_c)
a0, b0, c0 = np.mean(dv), 2*np.mean(dv*np.cos(th_)), 2*np.mean(dv*np.sin(th_))
gap_pred = np.sqrt(a0**2 - b0**2 - c0**2)/2
Gm = stv["W"].conj().T @ stv["G"] @ stv["W"]
Xc = PWc @ (np.eye(M) - 2*Gm)
Xc = (Xc + Xc.conj().T)/2
lc = np.where(np.abs(stv["xs"]) < 3)[0]
gap_num = float(np.linalg.eigvalsh(Xc[np.ix_(lc, lc)])[0])
lo = np.where(np.abs(stv["xs"]) < 2)[0]
comm_c = KJc @ KIc - KIc @ KJc
lie_c = min(np.abs((comm_c - e*2j*np.pi*(KJc - KIc))[np.ix_(lo, lo)]).max() for e in (1, -1)) \
    / np.abs((2*np.pi*(KJc - KIc))[np.ix_(lo, lo)]).max()
evc, Vc = np.linalg.eigh(Xc)
ip = int(np.argmin(np.where(stv["xs"] > 0, stv["xs"], np.inf)))
lap16 = float(np.sum(np.abs(Vc[ip])**2*np.exp(-16*evc)))
ZAHLEN["elliptische_kontrolle"] = {"v_J - v_I": [r(a0, 4), r(b0, 4), r(c0, 4)], "luecke_vorhergesagt": r(gap_pred, 4),
                                   "luecke_kompression": r(gap_num, 4), "lie_defekt_min_ueber_vorzeichen": r(lie_c, 3),
                                   "laplace_t16": r(lap16, 3)}
check(f"Kontrolle ohne gemeinsamen Endpunkt (J = (pi/2, 3pi/2), I = (3pi/4, 5pi/4)): v_J - v_I = {a0:.3f} + {b0:.3f} cos "
      f"theta ist elliptisch (a^2 > b^2 + c^2), K_J - K_I ist positiv mit Lücke (vorhergesagt {gap_pred:.3f}, Kompression "
      f"{gap_num:.3f}), die Laplace-Transformierte zerfällt exponentiell ({lap16:.1e} bei t = 16) und die Lie-Relation "
      f"verfehlt für beide Orientierungen ({lie_c:.2f}): Positivität allein macht keine Translation, sondern wieder eine Uhr",
      a0**2 > b0**2 + c0**2 and gap_num > 0.9*gap_pred and gap_num > 0.3 and lap16 < 1e-2 and lie_c > 0.3)
befund("Positivität von K_M - K_N ist allgemein (S_N in S_M gibt Delta_M <= Delta_N in Formen, log ist operatormonoton); "
       "erst der gemeinsame Endpunkt macht die Differenz parabolisch (Kontinuum ab 0, Translation) statt elliptisch "
       "(Lücke, periodisch)", True)
ZAHLEN["laufzeit_bis_teil2_s"] = r(time.time() - T0, 3)

# ================================================================== Teil 3
say("\nTeil 3: Originalzustand (KMS bei beta = 4 pi) -- global und auf Bögen")

# 3a global: Modulfluss des ganzen Kreises periodisch -> keine echte hsm (Takesaki)
per = {}
for al in (0.25, 0.75):
    xs = np.arange(-40, 40) + al
    per[al] = float(np.max(np.abs(np.exp(-2j*4*np.pi*xs) - 1)))
st = lattice_state(512, 0.25, 4*np.pi)
Mg = st["M"]
I_g = np.arange(Mg//2 - 3*Mg//8, Mg//2)
Xpos = st["U"] @ np.diag(st["xs"]) @ st["U"].conj().T
wX, VX = np.linalg.eigh(Xpos)
thg = st["th"]
fpk = np.exp(-(thg - (np.pi/4 + 1.0))**2/(2*0.2**2)).astype(complex)
fpk /= np.linalg.norm(fpk)
outI = ~np.isin(np.arange(Mg), I_g)
leak_rot = {}
for sgn in (1, -1):
    vals = []
    for t in (0.05, 0.1, 0.15):
        gvec = VX @ (np.exp(-1j*sgn*t*4*np.pi*wX)*(VX.conj().T @ fpk))
        vals.append(float(np.linalg.norm(gvec[outI])**2))
    leak_rot[sgn] = max(vals)
ZAHLEN["global"] = {"max|e^{-2iK}-1|": {a_: r(v, 3) for a_, v in per.items()}, "leckage_rotation_+-": {k_: r(v, 3) for k_, v in leak_rot.items()}}
check("global: spec(K) = 4 pi (Z + alpha) liegt in pi Z (alpha = 1/4: pi(4Z+1), 3/4: pi(4Z+3)), also e^{-2iK} = 1 "
      f"(Fehler {max(per.values()):.1e}): der Modulfluss der ganzen Naht ist periodisch und verschiebt ein Paket im Bogen I "
      f"in beiden Zeitrichtungen heraus (Leckage + {leak_rot[1]:.2f}, - {leak_rot[-1]:.2f})",
      max(per.values()) < 1e-11 and min(leak_rot.values()) > 0.1)
befund("Satz (Takesaki): Ist sigma^M periodisch und sigma_t(N) in N für t >= 0, dann auch für alle t, also gibt es einen "
       "zustandserhaltenden bedingten Erwartungswert M -> N; ist Omega zyklisch für N, folgt N = M. Die globale Algebra des "
       "Originals hat daher keine echte halbseitige modulare Teilalgebra: die periodische Uhr schließt globale "
       "Translationen strukturell aus (passt zu z70)", True)

# 3b Bögen in Doppelgenauigkeit (Systemblock = untere Schranke der Verletzung)
orig = {}
for key, (M, al, beta) in {"a1/4_M512": (512, 0.25, 4*np.pi), "a3/4_M512": (512, 0.75, 4*np.pi),
                           "a1/4_M1024": (1024, 0.25, 4*np.pi), "NS_b4pi_M512": (512, 0.5, 4*np.pi),
                           "a1/4_b8pi": (512, 0.25, 8*np.pi), "a1/4_b16pi": (512, 0.25, 16*np.pi),
                           "a1/4_b32pi": (512, 0.25, 32*np.pi)}.items():
    orig[key] = analyse(lattice_state(M, al, beta), M//2, FR)
o = orig["a1/4_M512"]
ZAHLEN["original_doppelt"] = {k_: {"lie": r(v["lie_pred"], 3), "familie": r(v["familie"], 3), "wiesbrock": r(v["wiesbrock"], 3),
                                   "fock_min": r(v["poscomp"], 2),
                                   "gamma_t16": [r(x, 4) for x in v["gamma"][16.0]]} for k_, v in orig.items()}
say(f"   Original (Systemblock): {ZAHLEN['original_doppelt']}")
floor_lie, floor_fam = v5["lie_pred"], v5["familie"]
gdev = max(abs(x - 1) for x in o["gamma"][16.0])
check(f"Original beta = 4 pi, alpha = 1/4 und 3/4: Lie-Defekt {o['lie_pred']:.1e} (> 10x Vakuum {floor_lie:.1e}), "
      f"Translationsfamilie {o['familie']:.1e} (> 4x Vakuum {floor_fam:.1e}), Gamma-Gesetz verletzt (bis {gdev:.2f}), "
      "unverändert bei M = 1024: keine Wiesbrock-Algebra, keine affine Translationsfamilie",
      o["lie_pred"] > 10*floor_lie and orig["a3/4_M512"]["lie_pred"] > 10*floor_lie and o["familie"] > 4*floor_fam
      and gdev > 0.05 and abs(orig["a1/4_M1024"]["lie_pred"]/o["lie_pred"] - 1) < 0.2)
check("Positivität bleibt (Kompressionsminimum >= -1e-3 in allen Originalfällen), wie es der allgemeine Satz verlangt: "
      "Positivität unterscheidet nicht, die Algebra schon",
      all(v["poscomp"] > -1e-3 for v in orig.values()))
check("beta-Abtastung (Doppelgenauigkeit): Lie-Defekt fällt monoton mit beta (4 pi, 8 pi, 16 pi, 32 pi) und erreicht bei "
      "32 pi den Vakuumboden; der Defekt ist thermisch, nicht von der Verdrehung alpha (NS bei 4 pi ebenso verletzt)",
      orig["a1/4_M512"]["lie_pred"] > orig["a1/4_b8pi"]["lie_pred"] > orig["a1/4_b16pi"]["lie_pred"]
      > orig["a1/4_b32pi"]["lie_pred"] and orig["a1/4_b32pi"]["lie_pred"] < 2*floor_lie
      and orig["NS_b4pi_M512"]["lie_pred"] > 10*floor_lie)

# 3c hochgenau: die hsm-Definition direkt -- Delta_J^{it} A(I) Delta_J^{-it} in A(I) für t >= 0?
say("   hochgenau (python-flint, approx. Eigenzerlegung): Leckage eines glatten Pakets aus I unter exp(-i t h_J)")


def hp_rows(M, an, ad, cfun, prec):
    ctx.prec = prec
    xs = [arb(k - M//2) + arb(an)/ad for k in range(M)]
    cs = [cfun(x) for x in xs]
    rows = {}
    for d in range(-(M - 1), M):
        s_ = acb(0)
        for k in range(M):
            s_ += cs[k]*acb(2*xs[k]*d/M).exp_pi_i()
        rows[d] = s_/M
    return rows


def c_vac(x):
    return arb(1) if float(x.mid()) < 0 else arb(0)


def c_therm(beta_over_pi, xmin=0.0, xmax=1e9):
    def cf(x):
        xf = abs(float(x.mid()))
        if xmin <= xf < xmax:
            return 1/(1 + (arb.pi()*beta_over_pi*x).exp())
        return c_vac(x)
    return cf


def hp_h(rows, sites, prec, H=200.0):
    ctx.prec = prec
    n = len(sites)
    A = acb_mat([[rows[j - k] for k in sites] for j in sites])
    E, R = A.eig(right=True, algorithm="approx")
    lo_ = arb(2)**(-int(H/0.6931) - 5)
    hs = []
    for e in E:
        nu = e.real
        if nu < lo_:
            nu = lo_
        if 1 - nu < lo_:
            nu = 1 - lo_
        hs.append(max(-H, min(H, float(((1 - nu)/nu).log().mid()))))
    Rn = np.array([[complex(float(R[i, j].real.mid()), float(R[i, j].imag.mid())) for j in range(n)] for i in range(n)])
    Q, _ = np.linalg.qr(Rn)
    return np.array(hs), Q


def hp_leak(M, an, ad, cfun, prec=512, wid=0.2, d0=1.0, ts=(0.05, 0.1, 0.15, 0.2, 0.25, 0.3), tback=(0.3, 0.45, 0.6)):
    rows = hp_rows(M, an, ad, cfun, prec)
    js = M//2
    J = list(range(js - M//2, js))
    jI0 = js - 3*M//8
    hJ, V = hp_h(rows, J, prec)
    th = 2*np.pi*np.array(J)/M
    inI = np.array([j >= jI0 for j in J])
    f = np.exp(-(th - (2*np.pi*(jI0 - 0.5)/M + d0))**2/(2*wid**2)).astype(complex)
    f /= np.linalg.norm(f)
    fh = V.conj().T @ f
    fwd = [float(np.linalg.norm((V @ (np.exp(-1j*t*hJ)*fh))[~inI])**2) for t in ts]
    bwd = [float(np.linalg.norm((V @ (np.exp(1j*t*hJ)*fh))[~inI])**2) for t in tback]
    return max(fwd), max(bwd), float(np.linalg.norm(f[~inI])**2)


LK = {}
for key, args in {"vak_NS_M256": (256, 1, 2, c_vac), "orig_a1/4_M256": (256, 1, 4, c_therm(4)),
                  "orig_a3/4_M256": (256, 3, 4, c_therm(4)), "a1/4_b8pi_M256": (256, 1, 4, c_therm(8)),
                  "a1/4_b16pi_M256": (256, 1, 4, c_therm(16)), "vak_NS_M128": (128, 1, 2, c_vac),
                  "vak_a1/4_M128": (128, 1, 4, c_vac), "orig_a1/4_M128": (128, 1, 4, c_therm(4)),
                  "NS_b4pi_M128": (128, 1, 2, c_therm(4)), "nur_|x|<1_M128": (128, 1, 4, c_therm(4, 0, 1)),
                  "nur_|x|>=1_M128": (128, 1, 4, c_therm(4, 1))}.items():
    LK[key] = hp_leak(*args)
for pr in (384, 640):
    LK[f"orig_a1/4_M128_{pr}bit"] = hp_leak(128, 1, 4, c_therm(4), prec=pr)
ZAHLEN["leckage_hochgenau"] = {k_: {"vorwaerts_max": r(v[0], 3), "rueckwaerts_max": r(v[1], 3), "start": r(v[2], 2)}
                               for k_, v in LK.items()}
say(f"   {ZAHLEN['leckage_hochgenau']}")
check("Vakuum (M = 256, 512 Bit): unter dem Modulfluss von J bleibt ein Gausspaket aus I für t > 0 in I (Leckage "
      f"{LK['vak_NS_M256'][0]:.1e}), für t < 0 verlässt es I ({LK['vak_NS_M256'][1]:.2f}): Delta_J^(it) A(I) Delta_J^(-it) in A(I) "
      "genau einseitig, also halbseitige modulare Inklusion; verdrehtes Vakuum identisch (Eichäquivalenz)",
      LK["vak_NS_M256"][0] < 1e-10 and LK["vak_NS_M256"][1] > 0.3 and LK["vak_NS_M128"][0] < 1e-10
      and abs(LK["vak_a1/4_M128"][0] - LK["vak_NS_M128"][0]) < 1e-12 and LK["vak_NS_M256"][2] < 1e-11)
lo4, lo4b = LK["orig_a1/4_M256"][0], LK["orig_a1/4_M128"][0]
check(f"Original beta = 4 pi: Vorwärtsleckage {lo4:.1e} (alpha = 1/4) und {LK['orig_a3/4_M256'][0]:.1e} (alpha = 3/4), "
      f"konvergiert in M (M = 128: {lo4b:.1e}) und unabhängig von der Rechengenauigkeit (384/512/640 Bit: "
      f"{LK['orig_a1/4_M128_384bit'][0]:.4e}/{lo4b:.4e}/{LK['orig_a1/4_M128_640bit'][0]:.4e}); das sind > 10^6 x Vakuum: "
      "die Inklusion A(I) in A(J) ist im Originalzustand keine halbseitige modulare Inklusion",
      lo4 > 1e6*max(LK["vak_NS_M256"][0], 1e-16) and lo4 > 1e-4 and abs(lo4/lo4b - 1) < 0.1
      and abs(LK["orig_a3/4_M256"][0]/lo4 - 1) < 0.05
      and all(abs(LK[f"orig_a1/4_M128_{pr}bit"][0]/lo4b - 1) < 1e-6 for pr in (384, 640)))
check("die Verletzung fällt mit beta, aber nur langsam (4 pi > 8 pi > 16 pi, bei 16 pi noch > 1e-5); sie ist nicht an die "
      f"Verdrehung gebunden (NS bei 4 pi: {LK['NS_b4pi_M128'][0]:.1e}) und stammt zur Hälfte aus Moden mit Besetzung "
      f"< e^-5pi ({LK['nur_|x|>=1_M128'][0]:.1e} gegen {LK['nur_|x|<1_M128'][0]:.1e} vom niedrigsten Paar): "
      "logarithmische Empfindlichkeit des Modulflusses",
      lo4 > LK["a1/4_b8pi_M256"][0] > LK["a1/4_b16pi_M256"][0] > 1e-5 and LK["NS_b4pi_M128"][0] > 0.5*lo4b
      and LK["nur_|x|>=1_M128"][0] > 0.3*lo4b)
befund("Wo es geht und wo nicht: Im Grundzustand der rekonstruierten Uhr (beta = inf, jedes alpha) bilden Bögen mit "
       "gemeinsamem Endpunkt halbseitige modulare Inklusionen; die Wiesbrock-Differenz ist der parabolische Generator. Im "
       "Originalzustand (beta = 4 pi) scheitert es global exakt (Periodizität, Takesaki) und auf Bögen numerisch "
       f"konvergent (Leckage {lo4:.1e}, Lie-Defekt {o['lie_pred']:.1e}, Familie {o['familie']:.1e}); Positivität bleibt, "
       "die ax+b-Algebra nicht", True)
ZAHLEN["laufzeit_gesamt_s"] = r(time.time() - T0, 3)

# ================================================================== Teil 4
say("\nTeil 4: Urteil und Dimensionstypisierung")
befund("(a) Vakuumsektor: ja. P^(t*) ist positiv, absolut stetig auf [0, inf), nicht periodisch und entsteht ohne neues "
       "Prinzip aus Modulgeneratoren des Grundzustands der Uhr (Wiesbrock: (K_J - K_I)/2 pi). Er wirkt auf demselben "
       "Fockraum, in dem das Original die Dichtematrix e^{-4 pi L_0}/Z ist. Wahl bleibt der Fixpunkt t* (Punkt im "
       "Unendlichen); kanonisch ist nur die Familie {P^(t*)}; die Uhr ist das Mittel zweier antipodaler Translationen, "
       "L_0 = (P^(pi) + P^(0))/2", True)
befund("(b) Originalzustand: nein. Er liefert keine halbseitige modulare Inklusion -- global strukturell nicht "
       "(periodisch), auf Bögen nicht (thermische Nichtlokalität des Modulflusses). Die Verbindung zum Vakuum ist ein "
       "Connes-Kozyklus (Dw_orig : Dw_vak)_t in A(J), nicht eine eigene Translation des Originals", True)
befund("Dimension (DIMENSION.UPLIFT.FIREWALL.01): alle Aussagen von z72 sind eindimensional (Kreis bzw. eine "
       "Lichtstrahlkoordinate, chirale Hälfte von 1+1); P^(t*) ist eine lichtartige Translation eines chiralen Netzes. "
       "Keine 3+1-Aussage. Unter der Himmelssphären-Lesart (T') wäre die Möbiusgruppe der Naht die Lorentzgruppe: P^(t*) "
       "wäre dort eine Nullrotation, der Bogen-Modulfluss ein Boost, L_0 eine Drehung -- keine Raumzeit-Zeittranslation", True)
befund("Nichtzirkularität (GRAV.NONCIRCULAR.01): Richtung ist Zustand -> Modulgeneratoren (Peschel) -> hsm (gerechnet) -> "
       "positive Translation (Satz von Wiesbrock/Borchers). Hislop-Longo dient nur als Prüfung des gerechneten Generators "
       "(RESULT, eindimensional), nicht als Annahme; es wird weder 'modular => Lorentz' noch BW in 3+1 gefolgert", True)

STATUS = {
    "Teil 1: L_0 elliptisch (diskret), P_0 parabolisch (absolut stetig [0, inf), Laguerre/Gamma), H_2 hyperbolisch, "
    "[H_2, P_0] = i P_0": "E",
    "Vakuum (beta = inf, jedes alpha): Bögen mit gemeinsamem Endpunkt sind hsm, (K_J - K_I)/2 pi = "
    "(cot(dI/2) - cot(dJ/2)) P^(t*)": "E (Satz Hislop-Longo/Wiesbrock; am Gitter aus dem Zustand nachgerechnet)",
    "streng geschachtelte Bögen: Differenz elliptisch, keine Translation": "E",
    "Positivität von K_J - K_I gilt allgemein und unterscheidet nicht": "E",
    "Original global: keine echte hsm (periodischer Modulfluss + Takesaki)": "E",
    "Original auf Bögen: keine hsm (Leckage, Lie-Defekt; konvergent in M und Rechengenauigkeit)": "E (numerisch)",
    "P^(t*) des Vakuumsektors als lichtartige Zeit der Naht": "L",
    "Himmelssphären-Lesart: P^(t*) ist Nullrotation, keine Raumzeit-Translation": "L (unter T')",
    "3+1-Zeittranslation, Skala von P^(t*), Wahl des Fixpunkts t*": "O",
}
ZAHLEN["status"] = STATUS
say("\nStatus (E/C/L/O)")
for k_, v_ in STATUS.items():
    say(f"   {v_}: {k_}")

urteil = ("Die fehlende lokale Translationsbeziehung existiert im Vakuumsektor der rekonstruierten Uhr, nicht im "
          "Originalzustand. Exakt: Auf dem Strom (h = 1) und dem Fermion (h = 1/2) ist L_0 elliptisch (diskret), "
          "P_0 = L_0 + (L_1 + L_-1)/2 parabolisch (unitär äquivalent zur Multiplikation mit s auf L^2((0, inf), Gamma-Mass), "
          "also absolut stetig [0, inf)) und H_2 hyperbolisch (Dilatation, e^{isH_2} P_0 e^{-isH_2} = e^{-s} P_0). Auf dem "
          "Gitter-Chiralfermion geben die Peschel-Modulgeneratoren zweier Bögen mit gemeinsamem Endpunkt im Grundzustand "
          "die Wiesbrock-Differenz (cot(delta_I/2) - cot(delta_J/2)) P^(t*), die ax+b-Relation und eine affine "
          f"Translationsfamilie (Defekte <= {max(max(v5['hl']), v5['lie_pred'], v5['familie']):.0e}, hochgenaue Leckage "
          f"{LK['vak_NS_M256'][0]:.0e}); ohne gemeinsamen Endpunkt wird die "
          "Differenz elliptisch (Lücke, wieder eine Uhr). Im Originalzustand beta = 4 pi scheitert die hsm global exakt "
          f"(periodischer Modulfluss, Takesaki) und auf Bögen konvergent (Leckage {lo4:.0e}, Lie-Defekt {o['lie_pred']:.0e}), "
          "obwohl die "
          "Positivität von K_J - K_I bleibt; der Defekt ist thermisch und fällt mit beta nur langsam. Alle Aussagen sind "
          "eindimensional (chirale Naht); P^(t*) ist eine lichtartige Translation ohne eigene Skala und unter der "
          "Himmelssphären-Lesart eine Nullrotation, keine 3+1-Zeittranslation.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)}); Laufzeit "
    f"{time.time() - T0:.0f} s")
(HERE / "z72_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
