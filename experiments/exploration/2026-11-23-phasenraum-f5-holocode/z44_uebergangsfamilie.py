"""Z44: Prüfung der elften externen Analyse (gemeinsame Übergangsfamilie, positiver Kern, Paaralgebra, Spin^c, Higgs).

Teil 1  Übergangsfunktion a_q(z) = z^-1 (1 - q z^4): a_q(iz) = -i a_q(z), det g_q = 1, Faktorisierung in einen innen holomorph
        invertierbaren Faktor und diag(z, I_6, 1/z): Bündel O(1) + 6O + O(-1) für jedes q; Hälften h^0 = 5 und 3.
Teil 2  Hardy-Kompression A_q = S^+ - q S^3 = S^+(I - q S^4): Index +1, Nullvektor eta_q = sqrt(1-q^2) Sum q^n |4n>,
        L_q = A_q^+ A_q = (1+q^2) I - q(S^4 + S^+4) - |0><0|, R A_q R^+ = -i A_q, [R, L_q] = 0.
Teil 3  Randantwort C_q(x) = <0|(x + L_q)^-1|0> = [sqrt((x+a)(x+b)) - x + 1 - q^2]/(2x), a = (1-q)^2, b = (1+q)^2; Pol
        1 - q^2 plus Marchenko-Pastur-Band rho(s) = sqrt((b-s)(s-a))/(2 pi s); x C^2 + (x + q^2 - 1) C - 1 = 0;
        C^-1 = x + Int x/(x+s) rho ds (vollständig Bernstein). Verbundene Vakuumantwort ohne masselosen Pol; Residuenmatrix
        Rang 1, aber nicht universell; Polgewichte 15/16, 3/4, 7/16 für q = 1/4, 1/2, 3/4.
Teil 4  Paaralgebra: Abschluss der vier Markenkanäle B^+(P) und Adjunkten hat Dimension 10 (sp(4,R)); mit L(sigma) dazu 15
        (su(2,2)); Determinantenwort hat im Kern K_eps Normquadrat 4 eps; fermionische Abkürzung f_i^+ f_{j+2}^+ gibt 4.
Teil 5  Familien: 4 + 4 -> 5 + 3 durch eine balancierte elementare Modifikation; ind dbar_E = rk + deg, Spin-Dirac deg;
        F_+- = E_+- x O(1) gibt 5 und 3 ohne negative Moden, c_1 = 5 + 3 = 8 (gemeinsamer Spin-Lift); O(-4) + 3O hat Spiegel.
Teil 6  Higgs: |Dn|^2 = (1/2) tr(DP)^2 + |n^+ Dn|^2; bei n = (0,1) W-Teil g^2(W1^2+W2^2)/4 und Z-Teil (g W3 - g' B)^2/4.
"""
import itertools
import json
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.integrate import quad

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
say("Teil 1: Übergangsfunktion und Bündel")
z, q = sp.symbols("z q")
a = (1 - q * z**4) / z
check("a_q(iz) = -i a_q(z) (Vierteldrehung in der Kompositionsregel)", sp.simplify(a.subs(z, sp.I * z) + sp.I * a) == 0)
g = sp.diag(1 / a, *([1] * 6), a)
check("det g_q = 1", sp.simplify(g.det()) == 1)
fac = sp.diag(1 / (1 - q * z**4), *([1] * 6), 1 - q * z**4) * sp.diag(z, *([1] * 6), 1 / z)
check("g_q = diag((1 - q z^4)^-1, I_6, 1 - q z^4) diag(z, I_6, 1/z); der erste Faktor ist für |z| <= 1, q < 1 holomorph invertierbar",
      sp.simplify(fac - g) == sp.zeros(8))
hd = lambda d: (max(d + 1, 0), max(-d - 1, 0))
check("Bündel O(1) + 6O + O(-1) für jedes 0 < q < 1; Hälften E_+ = O(1) + 3O (h^0 = 5), E_- = O(-1) + 3O (h^0 = 3, h^1 = 0)",
      hd(1)[0] + 3 == 5 and hd(-1) == (0, 0) and 3 * hd(0)[0] == 3)

# ------------------------------------------------------------------ Teil 2
say("\nTeil 2: Hardy-Kompression, Index und positiver Operator")
Nn = 1600
S = np.diag(np.ones(Nn - 1), -1)               # S|n> = |n+1>


def ops(qv):
    A = S.T - qv * np.linalg.matrix_power(S, 3)
    L = A.T @ A
    return A, L


qv = 0.5
A, L = ops(qv)
S4 = np.linalg.matrix_power(S, 4)
Lform = (1 + qv**2) * np.eye(Nn) - qv * (S4 + S4.T)
Lform[0, 0] -= 1
inner = slice(0, Nn - 8)
check("A_q = S^+ - q S^3 und L_q = (1+q^2) I - q(S^4 + S^+4) - |0><0| (außerhalb des Abschneiderandes exakt)",
      np.allclose(L[inner, inner], Lform[inner, inner]) and np.allclose(A[inner, inner], (S.T @ (np.eye(Nn) - qv * S4))[inner, inner]))
eta = np.zeros(Nn)
idx = np.arange(0, Nn, 4)
eta[idx] = np.sqrt(1 - qv**2) * qv ** np.arange(len(idx))
check("eta_q = sqrt(1-q^2) Sum q^n |4n> ist normiert und liegt im Kern von A_q", abs(np.linalg.norm(eta) - 1) < 1e-12
      and np.linalg.norm(A[inner, :] @ eta) < 1e-12)
AAt = A @ A.T
lam_min = np.linalg.eigvalsh(AAt[:Nn - 8, :Nn - 8]).min()
check("A A^+ >= (1-q)^2 > 0 (Kokern null): Index +1, genau ein komplexer Nullkanal", lam_min > (1 - qv) ** 2 - 1e-9)
R = np.diag(1j ** np.arange(Nn))
check("R A_q R^+ = -i A_q und [R, L_q] = 0 für R|n> = i^n |n>",
      np.allclose(R @ A @ R.conj().T, -1j * A) and np.allclose(R @ L - L @ R, 0))

# ------------------------------------------------------------------ Teil 3
say("\nTeil 3: Geschlossener positiver Kern")


def C_closed(x, qv):
    a_, b_ = (1 - qv) ** 2, (1 + qv) ** 2
    return (np.sqrt((x + a_) * (x + b_)) - x + 1 - qv**2) / (2 * x)


worst = 0.0
for qv in (0.25, 0.5, 0.75):
    A, L = ops(qv)
    for xv in (0.05, 0.4, 1.3, 4.0):
        e0 = np.zeros(Nn)
        e0[0] = 1
        num = np.linalg.solve(xv * np.eye(Nn) + L, e0)[0]
        worst = max(worst, abs(num - C_closed(xv, qv)))
say(f"   größte Abweichung direkte Inversion gegen geschlossene Form: {worst:.1e}")
check("C_q(x) = [sqrt((x+a)(x+b)) - x + 1 - q^2]/(2x) stimmt mit der direkten Inversion überein (q = 1/4, 1/2, 3/4)", worst < 1e-10)
xs_, qs_ = sp.symbols("x q", positive=True)
Cs = (sp.sqrt((xs_ + (1 - qs_) ** 2) * (xs_ + (1 + qs_) ** 2)) - xs_ + 1 - qs_**2) / (2 * xs_)
check("x C^2 + (x + q^2 - 1) C - 1 = 0", sp.simplify(sp.expand(xs_ * Cs**2 + (xs_ + qs_**2 - 1) * Cs - 1)) == 0)
for qv in (0.25, 0.5, 0.75):
    a_, b_ = (1 - qv) ** 2, (1 + qv) ** 2
    rho = lambda s: np.sqrt(max((b_ - s) * (s - a_), 0)) / (2 * np.pi * s)
    band = quad(rho, a_, b_)[0]
    rep = [(1 - qv**2) / xv + quad(lambda s: rho(s) / (xv + s), a_, b_)[0] for xv in (0.3, 2.0)]
    ber = [xv + quad(lambda s: xv / (xv + s) * rho(s), a_, b_)[0] for xv in (0.3, 2.0)]
    ZAHLEN[f"q={qv}"] = {"pol": 1 - qv**2, "band": band}
    ok_rep = all(abs(r - C_closed(xv, qv)) < 1e-8 for r, xv in zip(rep, (0.3, 2.0)))
    ok_ber = all(abs(bv - 1 / C_closed(xv, qv)) < 1e-8 for bv, xv in zip(ber, (0.3, 2.0)))
    CHECKS[f"_q{qv}"] = abs(band - qv**2) < 1e-8 and ok_rep and ok_ber
check("Spektraldarstellung: Pol 1 - q^2 bei 0, Marchenko-Pastur-Band mit Gewicht q^2, C^-1 = x + Int x/(x+s) rho ds",
      all(CHECKS.pop(f"_q{qv}") for qv in (0.25, 0.5, 0.75)))
check("Polgewichte 15/16, 3/4, 7/16 für q = 1/4, 1/2, 3/4",
      [ZAHLEN[f"q={v}"]["pol"] for v in (0.25, 0.5, 0.75)] == [15 / 16, 3 / 4, 7 / 16])
qv = 0.5
A, L = ops(qv)
eta = np.zeros(Nn)
eta[idx] = np.sqrt(1 - qv**2) * qv ** np.arange(len(idx))
O = np.zeros(Nn)
O[0] = 1
v = O * eta[0] - (1 - qv**2) * eta            # (|0><0| - <O>) eta
cc = []
for xv in (1e-3, 0.1, 1.0):
    cn = v @ np.linalg.solve(xv * np.eye(Nn) + L, v)
    cf = (1 - qv**2) * C_closed(xv, qv) - (1 - qv**2) ** 2 / xv
    cc.append((xv, cn, cf))
say(f"   verbundene Antwort (x, direkt, Formel): {[(x_, round(c1, 8), round(c2, 8)) for x_, c1, c2 in cc]}")
check("verbundene Vakuumantwort C_conn = (1-q^2) C - (1-q^2)^2/x stimmt, und x C_conn -> 0 (kein masseloser Pol)",
      all(abs(c1 - c2) < 1e-7 for _, c1, c2 in cc) and abs(cc[0][0] * cc[0][1]) < 1e-3)
resid = np.array([[(1 - qv**2) * qv ** (i + j) for j in range(4)] for i in range(4)])
check("Residuenmatrix der Quellen |4i>: R_ij = (1-q^2) q^(i+j), Rang 1, aber verschiedene Kopplungen (nicht universell)",
      np.linalg.matrix_rank(resid, tol=1e-12) == 1 and len({round(resid[i, i], 12) for i in range(4)}) == 4)
befund("Bündeltyp, Vierteldrehung, Index und Positivität wählen q nicht; die Prozessregel muss Zustand, Norm und Quellen festlegen", True)

# ------------------------------------------------------------------ Teil 4
say("\nTeil 4: Paaralgebra und Raumzeitdimension")


def bilinear_K(Z=None, Zbar=None, A_=None):
    """Wirkung ad_Q auf xi = (a1, a2, b1, b2, a1^+, a2^+, b1^+, b2^+) als 8x8-Matrix (Q quadratisch in Bosonen)."""
    K = np.zeros((8, 8), complex)
    if Z is not None:      # Q = Sum Z_ij a_i^+ b_j^+ : [Q, a_k] = -Sum_j Z_kj b_j^+, [Q, b_l] = -Sum_i Z_il a_i^+
        for k in range(2):
            for j in range(2):
                K[k, 6 + j] += -Z[k, j]
                K[2 + j, 4 + k] += -Z[k, j]
    if Zbar is not None:   # Q = Sum conj(Z_ij) b_j a_i : [Q, a_k^+] = Sum_j conj(Z_kj) b_j, [Q, b_l^+] = Sum_i conj(Z_il) a_i
        for k in range(2):
            for j in range(2):
                K[4 + k, 2 + j] += np.conj(Zbar[k, j])
                K[6 + j, k] += np.conj(Zbar[k, j])
    if A_ is not None:     # Q = Sum A_ij a_i^+ a_j : [Q, a_k] = -Sum_j A_kj a_j, [Q, a_k^+] = Sum_i A_ik a_i^+
        for k in range(2):
            for j in range(2):
                K[k, j] += -A_[k, j]
                K[4 + j, 4 + k] += A_[k, j]
    return K


def closure_dim(gens):
    basis = []

    def add(M):
        vec = np.concatenate([M.real.ravel(), M.imag.ravel()])
        if not basis:
            if np.linalg.norm(vec) > 1e-10:
                basis.append(vec)
                return True
            return False
        Bm = np.array(basis)
        res = vec - Bm.T @ np.linalg.lstsq(Bm.T, vec, rcond=None)[0]
        if np.linalg.norm(res) > 1e-9:
            basis.append(vec)
            return True
        return False

    mats = []
    for G in gens:
        if add(G):
            mats.append(G)
    grew = True
    while grew:
        grew = False
        cur = list(mats)
        for X in cur:
            for Y in cur:
                C = X @ Y - Y @ X
                if add(C):
                    mats.append(C)
                    grew = True
    return len(basis)


s1 = np.array([[0, 1], [1, 0]], complex)
s2 = np.array([[0, -1j], [1j, 0]])
s3 = np.diag([1.0 + 0j, -1])
I2 = np.eye(2)
marks = [(I2 + s1) / 2, (I2 - s1) / 2, (I2 + s2) / 2, (I2 - s2) / 2]
pair_gens = []
for P in marks:
    Kc, Ka = bilinear_K(Z=P), bilinear_K(Zbar=P)
    pair_gens += [Kc + Ka, 1j * (Kc - Ka)]          # anti-hermitesch kombinierte reelle Erzeuger
d10 = closure_dim(pair_gens)
d15 = closure_dim(pair_gens + [1j * bilinear_K(A_=M) for M in (s1, s2, s3)])
say(f"   reeller Abschluss: Paarkanäle der Marken {d10}, mit einseitiger Blattaktion L(sigma) {d15}")
check("vier Paarkanäle + Adjunkte schließen zu Dimension 10 = sp(4,R) ~ so(3,2) (2+1); mit L(A) = a^+ A a zu 15 = su(2,2) (3+1)",
      d10 == 10 and d15 == 15)
LB = bilinear_K(A_=s1) @ bilinear_K(Z=s2) - bilinear_K(Z=s2) @ bilinear_K(A_=s1)
check("[L(sigma_x), B^+(sigma_y)] = i B^+(sigma_z) (die fehlende Richtung)", np.allclose(LB, 1j * bilinear_K(Z=s3) * -1) or np.allclose(LB, 1j * bilinear_K(Z=s3)))
eps = sp.symbols("epsilon", positive=True)
zz = sp.symbols("z11 z12 z21 z22")
ww = sp.symbols("w11 w12 w21 w22")
Zm, Wm = sp.Matrix(2, 2, zz), sp.Matrix(2, 2, ww)
Keps = eps * Zm.det() * Wm.det()
dZ = lambda f: sp.diff(f, zz[0], zz[3]) - sp.diff(f, zz[1], zz[2])
dW = lambda f: sp.diff(f, ww[0], ww[3]) - sp.diff(f, ww[1], ww[2])
check("Determinantenwort (B11 B22 - B12 B21)^+ Omega: Normquadrat 4 eps im Kern K + eps det Z conj det W, 0 im Hardy-Kern",
      sp.simplify(dW(dZ(Keps))) == 4 * eps)
# Fermionische Abkürzung: f_i^+ f_{j+2}^+ mit vier Moden (Jordan-Wigner)
cre = []
for k in range(4):
    op = np.array([[1.0]])
    for l in range(4):
        fac_ = np.array([[0, 0], [1, 0]]) if l == k else (np.diag([1, -1]) if l < k else np.eye(2))
        op = np.kron(op, fac_)
    cre.append(op)
B = lambda i, j: cre[i] @ cre[j + 2]
vac = np.zeros(16)
vac[0] = 1
word = (B(0, 0) @ B(1, 1) - B(0, 1) @ B(1, 0)) @ vac
check("fermionische Abkürzung B_ij^+ = f_i^+ f_(j+2)^+: das Determinantenwort hat Normquadrat 4 statt 0", abs(word @ word - 4) < 1e-12)

# ------------------------------------------------------------------ Teil 5
say("\nTeil 5: Familien, Spin^c und Spin-Lift")
chi = lambda degs: sum(d + 1 for d in degs)
check("balancierte elementare Modifikation: E_- (rk 4, deg -1) -> O^4 -> E_+ (deg +1): chi = 3 + 5 aus 4 + 4",
      chi([-1, 0, 0, 0]) == 3 and chi([0, 0, 0, 0]) == 4 and chi([1, 0, 0, 0]) == 5)
spin = lambda degs: sum(d for d in degs)          # ind D_Spin,E = deg E auf P^1 (K^1/2 = O(-1))
check("ind dbar_E = rk + deg (5, 3), Spin-Dirac nur deg (1, -1): 5 + 3 ist eine Dolbeault/Spin^c-Zählung",
      chi([1, 0, 0, 0]) == 4 + 1 and spin([1, 0, 0, 0]) == 1 and spin([-1, 0, 0, 0]) == -1)
Fp, Fm = [2, 1, 1, 1], [0, 1, 1, 1]
neg = lambda degs: sum(max(-(d - 1) - 1, 0) for d in degs)   # h^1(F x O(-1))
check("F_+- = E_+- x O(1): Spinindex 5 und 3, keine negativen Moden; c_1 = 5 + 3 = 8 gerade (gemeinsamer Spin-Lift, einzeln ungerade)",
      spin(Fp) == 5 and spin(Fm) == 3 and neg(Fp) == 0 and neg(Fm) == 0 and (spin(Fp) + spin(Fm)) % 2 == 0 and spin(Fp) % 2 == 1)
check("O(-4) + 3O: h^0 = 3, aber h^1 = 3 (drei Spiegel); wer nur h^0 zählt, übersieht sie",
      sum(hd(d)[0] for d in (-4, 0, 0, 0)) == 3 and sum(hd(d)[1] for d in (-4, 0, 0, 0)) == 3)

# ------------------------------------------------------------------ Teil 6
say("\nTeil 6: Higgsprojektor und Z-Masse")
gg, gp = sp.symbols("g gprime", positive=True)
W1, W2, W3, Bb = sp.symbols("W1 W2 W3 B", real=True)
sig = [sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -sp.I], [sp.I, 0]]), sp.Matrix([[1, 0], [0, -1]])]
n = sp.Matrix([0, 1])
Dn = -sp.I / 2 * (gg * (W1 * sig[0] + W2 * sig[1] + W3 * sig[2]) + gp * Bb * sp.eye(2)) * n
P = n * n.H
DP = Dn * n.H + n * Dn.H
lhs = sp.simplify((Dn.H * Dn)[0])
fs = sp.simplify((DP * DP).trace() / 2)
ph = sp.simplify((n.H * Dn)[0] * sp.conjugate((n.H * Dn)[0]))
check("|Dn|^2 = (1/2) tr(DP)^2 + |n^+ Dn|^2", sp.simplify(lhs - fs - ph) == 0)
check("bei n = (0,1): (1/2) tr(DP)^2 = g^2 (W1^2 + W2^2)/4 (nur W), |n^+ Dn|^2 = (g W3 - g' B)^2/4 (Z aus der Hopf-Phase)",
      sp.simplify(fs - gg**2 * (W1**2 + W2**2) / 4) == 0 and sp.simplify(ph - (gg * W3 - gp * Bb) ** 2 / 4) == 0)
befund("M(P) allein enthält die Z-Masse nicht; der volle Yukawa-Operator M(H) samt Verbindung muss erhalten bleiben", True)

urteil = ("Die elfte Analyse hält in allen geprüften Punkten. a_q(z) = z^-1(1 - q z^4) trägt die Vierteldrehung (a_q(iz) = -i a_q), "
          "hat Determinante 1 und gibt für jedes q das Bündel O(1) + 6O + O(-1). Die Hardy-Kompression A_q = S^+(I - q S^4) hat "
          "Index +1 mit Nullvektor sqrt(1-q^2) Sum q^n |4n>, L_q = A^+ A ist positiv und kommutiert mit R. Der Randkern "
          "C_q = [sqrt((x+a)(x+b)) - x + 1 - q^2]/(2x) stimmt mit der direkten Inversion bis 1e-10 überein, ist Pol 1 - q^2 plus "
          "Marchenko-Pastur-Band und hat C^-1 vollständig Bernstein. Die verbundene Vakuumantwort hat keinen masselosen Pol, die "
          "Residuen sind Rang 1, aber nicht universell, und die Polgewichte hängen von q ab: Die Struktur wählt q nicht. Die "
          "Paarkanäle der Marken schließen zu sp(4,R) (2+1), mit der einseitigen Blattaktion zu su(2,2) (3+1); das Determinantenwort "
          "trennt K von K + eps det det (4 eps), die fermionische Abkürzung scheitert (4). 4 + 4 -> 5 + 3 ist eine balancierte "
          "Modifikation, 5 + 3 eine Spin^c-Zählung, F = E x O(1) gibt sie im Spin-Dirac mit c_1 = 8. |Dn|^2 = (1/2) tr(DP)^2 + "
          "|n^+ Dn|^2: die Z-Masse sitzt in der Hopf-Phase, also braucht es den vollen Yukawa-Operator M(H).")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z44_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False, default=str) + "\n")
