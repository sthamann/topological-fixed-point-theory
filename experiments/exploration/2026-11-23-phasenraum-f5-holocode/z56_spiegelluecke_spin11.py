"""Z56: Prüfung des lokalen Spiegel-Hamiltonoperators H = -(J/4) Sum_A (B_A + B_A^+)^2 mit verborgener Spin(2k+1)-Algebra.

Teil 1  Spin(10): zehn symmetrische 16 x 16-Matrizen M_A mit M_A^+ M_B + M_B^+ M_A = 2 delta_AB (Blöcke von C Gamma_A).
Teil 2  Nach der Teilchen-Loch-Abbildung ist H = -J Sum_A X_A^2 mit X_A = c^+ L_(A,N) c, L_(A,N) = (i/2) Gamma_A Gamma_N, also
        H = J (C_2(Spin(2k)) - C_2(Spin(2k+1))). Exakt diagonalisiert für Spin(6) -> Spin(7) (256 Zustände) und Spin(8) ->
        Spin(9) (Halbfüllung 12 870 Zustände): Grundzustand eindeutig, E_0 = -k'(k' + N - 2) J mit k' = 2^(k-2), Lücke
        (N - 2) J, Fermionlücke aus Casimir-Differenzen.
Teil 3  Spin(11) ohne den 2^32-Fockraum: Polynomrealisierung der Darstellung (8,0,0,0,0) auf der Sphäre S^10 (Singulett =
        zonale Harmonische, X_A = -i(x_A d_N - x_N d_A)). Gegen die Fockzahlen für N = 7, 9 geeicht; für N = 11 folgen
        E_0 = -136 J, Lücke 9 J, Fermionlücke 29/2 J, <X_A X_B> = 68/5, die Vierpunktfunktion und E(h) bis h^4.
Teil 4  Starker Kopplungsgrenzfall: Die naive Gitterdispersion sum sin^2 k hat acht Knoten mit Chiralitäten 1+, 3-, 3+, 1-.
Teil 5  E8-Paarkanäle: der gekreuzte Operator auf End(16) = 1 + 45 + 210 hat die Eigenwerte 5/8, -3/8, 1/8.
"""
import json
from fractions import Fraction as Fr
from functools import reduce
from itertools import combinations, product
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigsh

HERE = Path(__file__).parent
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


s1 = np.array([[0, 1], [1, 0]], complex)
s2 = np.array([[0, -1j], [1j, 0]])
s3 = np.diag([1.0, -1.0]).astype(complex)
I2 = np.eye(2, dtype=complex)
kron = lambda *ms: reduce(np.kron, ms)


def gammas(k):
    out = []
    for j in range(k):
        for s in (s1, s2):
            out.append(kron(*([s3]*j + [s] + [I2]*(k - j - 1))))
    chi = kron(*([s3]*k))
    return out, chi


# ------------------------------------------------------------------ Teil 1
say("Teil 1: Spin(10)-Paarmatrizen")
G, chi = gammas(5)
plus = np.where(np.isclose(np.diag(chi).real, 1))[0]
found = None
for name, C in (("Gamma_odd", reduce(np.matmul, G[0::2])), ("Gamma_even", reduce(np.matmul, G[1::2]))):
    blocks = [(C @ g)[np.ix_(plus, plus)] for g in G]
    if all(np.linalg.norm(b) > 1e-9 for b in blocks):
        found = (name, blocks)
        break
name, M = found
sym = max(np.max(np.abs(m - m.T)) for m in M)
cliff = max(np.max(np.abs(M[a].conj().T @ M[b] + M[b].conj().T @ M[a] - 2*(a == b)*np.eye(16))) for a in range(10)
            for b in range(10))
say(f"   Ladungskonjugation aus dem Produkt {name}: zehn 16 x 16-Blöcke, Asymmetrie {sym:.1e}, Clifford-Fehler {cliff:.1e}")
check("zehn symmetrische 16 x 16-Matrizen M_A mit M_A^+ M_B + M_B^+ M_A = 2 delta_AB I existieren (Spin(10)-Vektorkanal "
      "in Sym^2(16))", sym < 1e-12 and cliff < 1e-12)


# ------------------------------------------------------------------ Teil 2
say("\nTeil 2: Exakte Diagonalisierung Spin(2k) -> Spin(2k+1)")


def generators(k):
    G, chi = gammas(k)
    LN = [0.5j*g @ chi for g in G]                                   # L_(A,N)
    return G, chi, LN


def sector(nmodes, nfill):
    states = [sum(1 << i for i in c) for c in combinations(range(nmodes), nfill)]
    return states, {s: i for i, s in enumerate(states)}


def bilinear(mat, nmodes, states, index):
    rows, cols, vals = [], [], []
    for a, b in zip(*np.nonzero(np.abs(mat) > 1e-12)):
        coef = mat[a, b]
        for j, s in enumerate(states):
            if not (s >> b) & 1:
                continue
            s2_ = s ^ (1 << b)
            sign = (-1)**bin(s & ((1 << b) - 1)).count("1")
            if (s2_ >> a) & 1:
                continue
            sign *= (-1)**bin(s2_ & ((1 << a) - 1)).count("1")
            s3_ = s2_ | (1 << a)
            rows.append(index[s3_])
            cols.append(j)
            vals.append(sign*coef)
    return csr_matrix((vals, (rows, cols)), shape=(len(states), len(states)))


def rep_casimir(lam, kind):
    """<lam, lam + 2 rho> for B_r (kind 'B') or D_r (kind 'D') with orthonormal weights."""
    r = len(lam)
    rho = [r - i - 0.5 for i in range(r)] if kind == "B" else [r - i - 1 for i in range(r)]
    return sum(l*(l + 2*p) for l, p in zip(lam, rho))


results = {}
for k in (3, 4):
    N = 2*k + 1
    G_, chi_, LN = generators(k)
    nm = 2**k
    kp = 2**(k - 2)
    levels = {}
    gs_vec = None
    for nf in range(nm + 1) if k == 3 else range(nm//2 - 2, nm//2 + 3):
        st, ix = sector(nm, nf)
        X = [bilinear(L, nm, st, ix) for L in LN]
        H = -sum(x @ x for x in X)
        if H.shape[0] <= 300:
            ev, vec = np.linalg.eigh(H.toarray())
        else:
            ev, vec = eigsh(H, k=6, which="SA")
            o = np.argsort(ev)
            ev, vec = ev[o], vec[:, o]
        levels[nf] = ev
        if nf == nm//2:
            gs_vec, Xh = vec[:, 0], X
    E0 = min(levels[nf][0] for nf in levels)
    nf0 = [nf for nf in levels if abs(levels[nf][0] - E0) < 1e-8]
    deg = sum(int(np.sum(np.abs(levels[nf] - E0) < 1e-8)) for nf in levels)
    second = min(e for nf in levels for e in levels[nf] if e > E0 + 1e-8)
    fermion = min(levels[nf][0] for nf in levels if abs(nf - nf0[0]) == 1) - E0
    mom = np.array([[np.vdot(gs_vec, Xh[a] @ (Xh[b] @ gs_vec)).real for b in range(2*k)] for a in range(2*k)])
    Cn = kp*(kp + N - 2)
    lam_f = [kp - 0.5] + [0.5]*(k - 1)
    pred_f = rep_casimir([0.5]*k, "D") + Cn - rep_casimir(lam_f, "B")
    results[k] = dict(E0=E0, deg=deg, gap=second - E0, fermion=fermion, x2=mom[0, 0], offdiag=np.max(np.abs(mom - np.diag(np.diag(mom)))))
    say(f"   Spin({2*k}) -> Spin({N}): E_0 = {E0:.6f} (Vorhersage {-Cn}), Füllung {nf0}, Entartung {deg}, Lücke {second - E0:.6f} "
        f"(Vorhersage {N - 2}), Fermionlücke {fermion:.6f} (Casimir {pred_f}), <X_A^2> = {mom[0, 0]:.6f} ({Cn/(N - 1):.6f})")
    check(f"Spin({2*k}) -> Spin({N}): eindeutiger Grundzustand bei Halbfüllung, E_0 = -k'(k'+N-2) J = {-Cn} J, Lücke "
          f"{N - 2} J, Fermionlücke {pred_f} J aus Casimir-Differenzen, <X_A X_B> = (C/(N-1)) delta_AB",
          abs(E0 + Cn) < 1e-8 and deg == 1 and abs(second - E0 - (N - 2)) < 1e-8 and abs(fermion - pred_f) < 1e-8
          and abs(mom[0, 0] - Cn/(N - 1)) < 1e-8 and results[k]["offdiag"] < 1e-8)
ZAHLEN["fock"] = {k: {kk: round(float(v), 9) for kk, v in d.items()} for k, d in results.items()}   # ARPACK-Rauschen gerundet

# ------------------------------------------------------------------ Teil 3
say("\nTeil 3: Spin(11) über die Polynomrealisierung")
x1, x2, z, s = sp.symbols("x1 x2 z s")


def zonal(N, kk):
    t = sp.symbols("t")
    r2 = x1**2 + x2**2 + z**2 + s
    poly = sp.Poly(sp.gegenbauer(kk, sp.Rational(N - 2, 2), t), t)
    out = 0
    for (deg,), c in poly.terms():
        out += c*z**deg*r2**sp.Rational(kk - deg, 2)
    return sp.expand(out)


def L(f, var):                                                    # L_(A,N) = x_A d_z - z d_x_A
    return sp.expand(var*sp.diff(f, z) - z*sp.diff(f, var))


def Lrest_sq(f, m):                                               # Sum over the m variables hidden in s of L_A^2
    g = sp.diff(f, z) - 2*z*sp.diff(f, s)
    return sp.expand(s*(sp.diff(g, z) - 2*z*sp.diff(g, s)) - m*z*g)


def inner(f, g, N):
    m = N - 3
    tot = 0
    for (a, b, c, d), coef in sp.Poly(sp.expand(f*g), x1, x2, z, s).terms():
        if a % 2 or b % 2 or c % 2:
            continue
        tot += coef*sp.gamma(sp.Rational(a + 1, 2))*sp.gamma(sp.Rational(b + 1, 2))*sp.gamma(sp.Rational(c + 1, 2)) \
            * sp.gamma(d + sp.Rational(m, 2))/sp.gamma(sp.Rational(m, 2))
    return sp.nsimplify(sp.simplify(tot))


def moments(N, kk):
    Y = zonal(N, kk)
    nrm = inner(Y, Y, N)
    X1Y, X2Y = L(Y, x1), L(Y, x2)
    x2m = inner(X1Y, X1Y, N)/nrm                                     # <X_1^2> = <L1 Y, L1 Y>
    x4 = inner(L(X1Y, x1), L(X1Y, x1), N)/nrm                        # <X_1^4>
    a_ = inner(L(L(Y, x1), x1), L(L(Y, x2), x2), N)/nrm              # <X_1 X_1 X_2 X_2> = alpha
    b_ = inner(L(X2Y, x1), L(X1Y, x2), N)/nrm                        # <X_2 X_1 X_2 X_1> = beta
    lap = -(L(L(Y, x1), x1) + L(L(Y, x2), x2) + Lrest_sq(Y, N - 3))
    cas = sp.simplify(sp.cancel(sp.expand(lap)/Y)) if Y != 0 else None
    return dict(x2=sp.nsimplify(x2m), x4=sp.nsimplify(x4), x12x21=sp.nsimplify(a_), x21x12=sp.nsimplify(b_), casimir=cas, Y=Y)


for N, kk in ((7, 2), (9, 4)):
    mo = moments(N, kk)
    say(f"   N = {N}, k' = {kk}: Sum_A X_A^2 Y = {mo['casimir']} Y, <X_1^2> = {mo['x2']}")
    check(f"N = {N}: Polynomrealisierung trifft die Fockzahlen (Sum X_A^2 = {kk*(kk + N - 2)}, <X_1^2> = {results[(N - 1)//2]['x2']:.4f})",
          mo["casimir"] == kk*(kk + N - 2) and abs(float(mo["x2"]) - results[(N - 1)//2]["x2"]) < 1e-9)
mo = moments(11, 8)
alpha, beta = mo["x12x21"], mo["x21x12"]
x4 = mo["x4"]
conn = sp.nsimplify(x4 - 3*mo["x2"]**2)
say(f"   N = 11, k' = 8: Sum X_A^2 = {mo['casimir']}, <X_1^2> = {mo['x2']}, <X_1^4> = {x4}, <X_1 X_1 X_2 X_2> = {alpha}, "
    f"<X_2 X_1 X_2 X_1> = {beta}, verbunden <X_1^4>_c = {conn}")
check("Spin(11): Casimir 136 (E_0 = -136 J), <X_A X_B> = 68/5 delta_AB", mo["casimir"] == 136 and mo["x2"] == sp.Rational(68, 5))
check("Vierpunktfunktion alpha delta_AB delta_CD + beta (delta_AC delta_BD + delta_AD delta_BC): alpha = 782/5, beta = 714/5, "
      "<X_1^4> = alpha + 2 beta = 442", alpha == sp.Rational(782, 5) and beta == sp.Rational(714, 5) and x4 == alpha + 2*beta == 442)
check("das verbundene vierte Moment ist -2822/25: die Antwort ist nicht gaußsch", conn == sp.Rational(-2822, 25))
# Spin(9)-Singulettkette exakt: Krylovraum von X_1, Phasen (-i)^n machen alle Matrizen reell
import mpmath as mpm
mpm.mp.dps = 60
Y = mo["Y"]
basis = [Y]
for _ in range(8):
    basis.append(L(basis[-1], x1))
g00 = inner(Y, Y, 11)
ratio = lambda e: sp.nsimplify(sp.simplify(e/g00))
Gm = sp.Matrix(9, 9, lambda i, j: ratio(inner(basis[i], basis[j], 11)))
Lm = sp.Matrix(9, 9, lambda i, j: ratio(inner(basis[i], L(basis[j], x1), 11)))
Sm = sp.Matrix(9, 9, lambda i, j: ratio(inner(basis[i], -(L(L(basis[j], x1), x1) + L(L(basis[j], x2), x2)
                                                         + Lrest_sq(basis[j], 8)), 11)))
ph = [(-sp.I)**n for n in range(9)]
Gp = sp.Matrix(9, 9, lambda i, j: sp.simplify(sp.conjugate(ph[i])*ph[j]*Gm[i, j]))
Sp_ = sp.Matrix(9, 9, lambda i, j: sp.simplify(sp.conjugate(ph[i])*ph[j]*Sm[i, j]))
Xp = sp.Matrix(9, 9, lambda i, j: sp.simplify(sp.conjugate(ph[i])*ph[j]*(-sp.I)*Lm[i, j]))
real = all(sp.im(v) == 0 for M_ in (Gp, Sp_, Xp) for v in M_)
C10 = sp.Matrix(136*sp.eye(9) - Gp.inv()*Sp_)
ev10 = sorted(sp.Matrix(C10).eigenvals().keys(), key=lambda v: float(v))
say(f"   Spin(9)-Singulettkette (exakt): C_2(Spin(10))-Eigenwerte {ev10}; Matrizen reell: {real}")
check("die Kette trägt j = 0..8 mit C_2(Spin(10)) = j(j+8) exakt: Lücke 9 J zum ersten Vektor",
      real and ev10 == [j*(j + 8) for j in range(9)])
toM = lambda A: mpm.matrix([[mpm.mpf(sp.Rational(A[i, j]).p)/mpm.mpf(sp.Rational(A[i, j]).q) for j in range(9)] for i in range(9)])
Gq, Sq, Xq = toM(Gp), toM(Sp_), toM(Xp)
Lc = mpm.cholesky(Gq)
Li = Lc**-1


def E_of(h):
    Hm = Li*(-Sq + h*Xq)*Li.T
    Hm = (Hm + Hm.T)/2
    return min(mpm.eigsy(Hm)[0]) + 136


hs = [mpm.mpf("0.002"), mpm.mpf("0.004"), mpm.mpf("0.006"), mpm.mpf("0.008")]
A_ = mpm.matrix([[h**2, h**4, h**6, h**8] for h in hs])
cof = mpm.lu_solve(A_, mpm.matrix([E_of(h) for h in hs]))
c2, c4 = float(cof[0]), float(cof[1])
say(f"   E(h) - E(0) = {c2:.12f} h^2 + {c4:.12f} h^4 (Angabe: {-68/45:.12f}, {17323/182250:.12f})")
check("E(h) - E(0) = -(68/45J) h^2 + (17323/(182250 J^3)) h^4: quadratische und quartische Antwort aus demselben Operator",
      abs(c2 + 68/45) < 1e-10 and abs(c4 - 17323/182250) < 1e-8)
lam_f = [7.5, 0.5, 0.5, 0.5, 0.5]
df = rep_casimir([0.5]*5, "D") + 136 - rep_casimir(lam_f, "B")
check("Fermionlücke aus den Casimir-Differenzen: C(16) + 136 - C(15/2, 1/2, 1/2, 1/2, 1/2) = 29/2 J", abs(df - 14.5) < 1e-12)
befund("Die Spin(10)-Zahlen E_0 = -136 J, Lücke 9 J, Fermionlücke 29/2 J und G(z) = z/(z^2 - Delta_f^2) folgen aus derselben "
       "Darstellungstheorie, die für Spin(7) und Spin(9) durch volle Diagonalisierung bestätigt ist; der 2^32-Fockraum selbst "
       "ist nicht diagonalisiert", True)

# ------------------------------------------------------------------ Teil 4
say("\nTeil 4: Starker Kopplungsgrenzfall -- naive Knoten")
corners = list(product((0, np.pi), repeat=3))
chis = [int(np.sign(np.prod(np.cos(c)))) for c in corners]
counts = {}
for c, ch in zip(corners, chis):
    npi = sum(1 for v in c if v > 0)
    counts[(npi, ch)] = counts.get((npi, ch), 0) + 1
say(f"   Knoten der Dispersion sqrt(sum sin^2 k): {len(corners)}; (Zahl der pi-Komponenten, Chiralität): {counts}")
check("acht Knoten mit Chiralitäten 1+ (0 pi), 3- (1 pi), 3+ (2 pi), 1- (3 pi), Summe null (Nielsen-Ninomiya): die starke "
      "Projektion lässt sieben Zusatzknoten stehen", counts == {(0, 1): 1, (1, -1): 3, (2, 1): 3, (3, -1): 1} and sum(chis) == 0)
befund("Die starke Projektion auf den Spiegel-Singulett ist kein Weg zur einzelnen Weylmode; die gesuchte Lösung liegt bei "
       "endlicher Kopplung oder in einem anderen lokalen Elternoperator (Rechnung des Elternraums aus der Analyse, hier nur "
       "die Knotenbilanz nachgeprüft)", True)

# ------------------------------------------------------------------ Teil 5
say("\nTeil 5: E8-Paarkanäle auf End(16)")
basisE = [np.zeros((16, 16), complex) for _ in range(256)]
for i in range(256):
    basisE[i][i//16, i % 16] = 1
Phi = np.array([(sum(M[a] @ E.T @ M[a].conj().T for a in range(10))/16).reshape(-1) for E in basisE]).T
ev = np.linalg.eigvals(Phi)
vals = sorted({round(v.real, 9) for v in ev})
mult = {v: int(np.sum(np.abs(ev.real - v) < 1e-7)) for v in vals}
say(f"   X -> (1/16) Sum_A M_A X^T M_A^+: Eigenwerte mit Vielfachheit {mult}")
check("der gekreuzte Vektorkanal hat genau drei Eigenwerte 5/8 (1), -3/8 (45), 1/8 (210): End(16) = 1 + 45 + 210",
      mult == {-0.375: 45, 0.125: 210, 0.625: 1})
befund("Die Kanalzahlen folgen aus Clifford-Kontraktionen ((10 - 2p)/16 mit Vorzeichen der Transposition); ihre Deutung als "
       "Kraftstärken braucht Hamiltonnormierung und Operatorordnung", True)

urteil = ("Die externe Spiegelkonstruktion hält. Zehn symmetrische Spin(10)-Paarmatrizen existieren; nach der "
          "Teilchen-Loch-Abbildung ist H = -J Sum X_A^2 = J (C_2(Spin(10)) - C_2(Spin(11))). Für Spin(7) und Spin(9) bestätigt die "
          "volle Diagonalisierung die Darstellungsformeln: eindeutiger Grundzustand bei Halbfüllung, E_0 = -k'(k'+N-2) J, Lücke "
          "(N-2) J, Fermionlücke aus Casimir-Differenzen. Mit derselben Theorie und einer exakten Polynomrealisierung folgen für "
          "Spin(11) E_0 = -136 J, Lücke 9 J, Fermionlücke 29/2 J, <X_A X_B> = 68/5, die Vierpunktfunktion mit 782/5 und 714/5, "
          "das nichtgaußsche verbundene Moment -2822/25 und E(h) = -(68/45) h^2 + (17323/182250) h^4. Der starke Grenzfall lässt "
          "acht naive Knoten (Chiralitäten 1+, 3-, 3+, 1-) stehen; der gekreuzte Kanal hat 5/8, -3/8, 1/8 auf 1 + 45 + 210. Die "
          "Wechselwirkung ist eine zusätzliche Regel, nicht aus der Naht abgeleitet, und die 4D-Einbettung bleibt offen.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z56_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False, default=str) + "\n")
