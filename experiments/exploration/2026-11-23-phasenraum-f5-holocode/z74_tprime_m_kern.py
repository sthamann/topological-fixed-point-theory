"""Z74: T' aus Tomita, die Windungszahl m und die Eindeutigkeit des Kerns (Agentenauftrag nach H6h34).

Leitlinie 'Rekonstruktion statt Auswahl'. Firewalls: GRAV.NONCIRCULAR.01 (die Lorentzgruppe kommt hier aus dem natürlichen
Kegel, nicht aus dem Modulfluss; kein Schluss modular => Lorentz => BW) und DIMENSION.UPLIFT.FIREWALL.01 (alles ist Kinematik;
keine 3+1-Schließung). Bestand: z62 (Standardform), z63 (rho_m), z64/z70 (V-bar^*, eta), z40 (det-Kern, Wallach), z71.

Teil 1  T' aus der Standardform von (End V, rho_m): J aus der Polarzerlegung von S; die drei Realformen von (M_2(C), det)
        und welche M mit M' vertauscht; Lorentzalgebra des Kegels; Cliffordmodul und Spinhebung (Partner = V-bar^*);
        Tomita-Partner = J-Bild; Twistorform Sigma (Phase aus J, Vorzeichen aus P); Geraden L_x, x = adj(a); was bleibt.
Teil 2  Windungszahl m in {1, 2, 3, 4} gegen (a) Spinbrücke, (b) Neutralität, (c) Realstruktur/J, (d) Zeitorientierung.
        Satz: R = Delta^(-i/2) ist auf jedem Quartersektor eine Funktion der Ladung allein.
Teil 3  Kern: U(2) x U(2)-invariante positive holomorphe Kerne auf der 2x2-Matrixkugel; zwei native Auslesungen (Hardy-Kern
        auf dem Rang-1-Schnitt, Stromkern auf dem Skalarschnitt) erzwingen det(I - Z W^+)^-1; Cayley nach Minkowski.
Teil 4  Bilanz.
"""
import itertools
import json
from math import factorial
from pathlib import Path

import mpmath as mp
import numpy as np
import scipy.linalg as sla
import sympy as sp
from scipy.optimize import linprog

HERE = Path(__file__).parent
mp.mp.dps = 40
OUT = []
CHECKS, BEFUNDE, ZAHLEN, STATUS = {}, {}, {}, {}
rng = np.random.default_rng(74)


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val):
    BEFUNDE[name] = bool(val)
    say(f"   [Befund: {'ja' if val else 'nein'}] {name}")


def status(key, lab, text):
    STATUS[key] = {"status": lab, "aussage": text}
    say(f"   [{lab}] {key}: {text}")


def crand(*s):
    return rng.normal(size=s) + 1j * rng.normal(size=s)


def rsl2():
    g = crand(2, 2)
    return g / np.sqrt(np.linalg.det(g))


def nullraum(M, tol=1e-9):
    _, s, vh = np.linalg.svd(M)
    r = int(np.sum(s > tol * max(1.0, s.max())))
    return vh[r:].conj()


I2 = np.eye(2, dtype=complex)
SX = np.array([[0, 1], [1, 0]], dtype=complex)
SY = np.array([[0, -1j], [1j, 0]])
SZ = np.diag([1.0, -1.0]).astype(complex)
PAULI = [I2, SX, SY, SZ]
EPS = np.array([[0, 1], [-1, 0]], dtype=complex)
EINH = [np.outer(I2[i], I2[j]) for i in range(2) for j in range(2)]


def adj(x):
    return np.array([[x[1, 1], -x[0, 1]], [-x[1, 0], x[0, 0]]])


def rho_m(m):
    w = np.array([1.0, np.exp(-np.pi * m)])
    return np.diag(w / w.sum()).astype(complex)


def msqrt(h):
    w, v = np.linalg.eigh(h)
    return (v * np.sqrt(np.clip(w, 0, None))) @ v.conj().T


def reell(f):
    """reelle 8x8-Matrix einer reell-linearen (auch antilinearen) Abbildung von M_2(C); Koordinaten (Re a_ij, Im a_ij)."""
    cols = []
    for E in EINH + [1j * E for E in EINH]:
        F = f(E)
        cols.append(np.concatenate([F.real.ravel(), F.imag.ravel()]))
    return np.array(cols).T


def aus_reell(v):
    return (v[:4] + 1j * v[4:]).reshape(2, 2)


def det_bilinear(a, b):
    return (np.linalg.det(a + b) - np.linalg.det(a) - np.linalg.det(b)) / 2


# =================================================================== Teil 1
say("Teil 1: T' aus der Standardform von (End V, rho_m)")

# --- 1.1 J aus der Polarzerlegung von S, für alle m
A_ = crand(2, 2)
rho_zuf = A_ @ A_.conj().T
rho_zuf /= np.trace(rho_zuf).real
zustaende = {f"rho_{m}": rho_m(m) for m in (1, 2, 3, 4)}
zustaende["zufällig"] = rho_zuf
J_dag = reell(lambda a: a.conj().T)
ok_J, ok_D, ok_comm = True, True, True
for nm, rho in zustaende.items():
    xi = msqrt(rho)
    xi_i = np.linalg.inv(xi)
    S_hat = reell(lambda a: xi_i @ a.conj().T @ xi)                  # S(x xi) = x^+ xi
    D_hat = S_hat.T @ S_hat                                             # Delta = S^* S
    w, v = np.linalg.eigh((D_hat + D_hat.T) / 2)
    J_hat = S_hat @ (v / np.sqrt(w)) @ v.T                              # J = S Delta^(-1/2)
    ok_J &= np.allclose(J_hat, J_dag, atol=1e-8) and np.allclose(J_hat @ J_hat, np.eye(8), atol=1e-8)
    ok_D &= np.allclose(D_hat, reell(lambda a: rho @ a @ np.linalg.inv(rho)), atol=1e-6 * np.abs(D_hat).max())
    for _ in range(5):
        x = crand(2, 2)
        JLJ = J_hat @ reell(lambda a: x @ a) @ J_hat
        ok_comm &= np.allclose(JLJ, reell(lambda a: a @ x.conj().T), atol=1e-8)
check("für rho_1 .. rho_4 und einen zufälligen treuen Zustand gibt die Polarzerlegung von S (S x xi = x^+ xi) J a = a^+ und "
      "Delta a = rho a rho^-1; J ist zustandsunabhängig, J L_x J = R_(x^+) (J M J = M')", ok_J and ok_D and ok_comm)

# --- 1.2 Die drei Realformen von (M_2(C), det)
realformen = {
    "J_T: a -> a^+ (Tomita)": lambda a: a.conj().T,
    "J_K: a -> s_x conj(a) s_x (Kreis |z| = 1 durch die Marken)": lambda a: SX @ a.conj() @ SX,
    "J_A: a -> eps conj(a) eps^-1 (Antipoden)": lambda a: EPS @ a.conj() @ np.linalg.inv(EPS),
}
tab_rf = {}
for nm, f in realformen.items():
    Jh = reell(f)
    fix = nullraum(Jh - np.eye(8)).real
    basis = [aus_reell(v) for v in fix]
    G = np.array([[det_bilinear(a, b) for b in basis] for a in basis])
    ew = np.linalg.eigvalsh((G.real + G.real.T) / 2)
    x = crand(2, 2)
    JLJ = Jh @ reell(lambda a: x @ a) @ Jh
    E_links = all(np.allclose(JLJ @ reell(lambda a: a @ y), reell(lambda a: a @ y) @ JLJ) for y in EINH)
    E_rechts = all(np.allclose(JLJ @ reell(lambda a: y @ a), reell(lambda a: y @ a) @ JLJ) for y in EINH)
    tab_rf[nm] = {"reell_dim": len(basis), "det_reell": bool(np.abs(G.imag).max() < 1e-10),
                  "signatur": (int(np.sum(ew > 1e-9)), int(np.sum(ew < -1e-9))),
                  "J M J": "M'" if (E_rechts and not E_links) else ("M" if (E_links and not E_rechts) else "?")}
    say(f"   {nm}: {tab_rf[nm]}")
ZAHLEN["realformen"] = tab_rf
sig = [tab_rf[k]["signatur"] for k in realformen]
art = [tab_rf[k]["J M J"] for k in realformen]
check("die drei antilinearen Involutionen geben reelle Formen der Dimension 4 mit reellem det und den Signaturen (1,3), (2,2), "
      "(4,0); nur die Tomita-Konjugation vertauscht M und M', die beiden eigenen Realstrukturen der Naht erhalten M",
      sig == [(1, 3), (2, 2), (4, 0)] and art == ["M'", "M", "M"]
      and all(v["reell_dim"] == 4 and v["det_reell"] for v in tab_rf.values()))
befund("Lorentzsignatur heißt hier genau: Die Realstruktur vertauscht die Algebra mit ihrer Kommutante. Das ist die definierende "
       "Eigenschaft der Tomita-Konjugation; Realstrukturen, die M erhalten (Spiegelung am Markenkreis, Antipoden), geben (2,2) "
       "bzw. (4,0) wie in z64. 'Konjugiert' in V-bar^* ist also die Antilinearität von J", True)

# --- 1.3 Lorentzalgebra des natürlichen Kegels
G4 = np.array([[det_bilinear(a, b).real for b in PAULI] for a in PAULI])
rows = []
for k in range(16):
    A = np.zeros(16)
    A[k] = 1
    A = A.reshape(4, 4)
    rows.append((A.T @ G4 + G4 @ A).ravel())
Mco = np.column_stack(rows + [-G4.ravel()])
co_basis = nullraum(Mco).real
dim_co = co_basis.shape[0]


def pauli_koord(a):
    return np.array([np.trace(s.conj().T @ a).real / 2 for s in PAULI])


def A_von_X(X):
    return np.column_stack([pauli_koord(X @ s + s @ X.conj().T) for s in PAULI])


gl_basis = EINH + [1j * E for E in EINH]
AX = [A_von_X(X) for X in gl_basis]
rank_AX = np.linalg.matrix_rank(np.array([a.ravel() for a in AX]), tol=1e-9)
co_A = [v[:16] for v in co_basis]
rank_both = np.linalg.matrix_rank(np.array([a.ravel() for a in AX] + co_A), tol=1e-9)
ok_mod = True
for m in (1, 2, 3, 4):
    rho = rho_m(m)
    Xm = 1j * np.diag(np.log(np.diag(rho).real))
    Am = A_von_X(Xm)
    ok_mod &= np.allclose(Am.T @ G4 + G4 @ Am, 0) and np.allclose(Am @ pauli_koord(msqrt(rho)), 0)
    ok_mod &= np.all(np.linalg.eigvals(Am).real < 1e-9) and np.abs(np.linalg.eigvals(Am).imag).max() > 1
check(f"Kegelautomorphismen: co(det) auf Herm(2) hat Dimension {dim_co} (Lorentz 6 + Dilatation 1) und ist genau "
      f"{{a -> X a + a X^+}} (Rang {rank_AX}, gemeinsamer Rang {rank_both}); der Modulgenerator von rho_m (m = 1..4) liegt in "
      "so(1,3), hält sqrt(rho) fest und ist eine Drehung (rein imaginäre Eigenwerte): Lorentzgruppe aus dem Kegel, nicht aus dem "
      "Modulfluss", dim_co == 7 and rank_AX == 7 and rank_both == 7 and ok_mod)


# --- 1.4 Cliffordmodul und Spinhebung
def gam(a):
    Z0 = np.zeros((2, 2), dtype=complex)
    return np.block([[Z0, a], [adj(a), Z0]])


ok_cl = True
for _ in range(20):
    a, b = crand(2, 2), crand(2, 2)
    ok_cl &= np.allclose(gam(a) @ gam(b) + gam(b) @ gam(a), 2 * det_bilinear(a, b) * np.eye(4))
prods = [np.eye(4, dtype=complex)]
for L in range(1, 5):                                                   # bis Länge 4: das Volumenelement gehört dazu
    for w in itertools.product(EINH, repeat=L):
        P = np.eye(4, dtype=complex)
        for e in w:
            P = P @ gam(e.astype(complex))
        prods.append(P)
dim_cl = np.linalg.matrix_rank(np.array([P.ravel() for P in prods]), tol=1e-9)
check(f"gamma(a) = [[0, a], [adj a, 0]] auf V + V erfüllt gamma(a) gamma(b) + gamma(b) gamma(a) = 2 B(a, b) mit B = polarisierte "
      f"Determinante (ohne Wahl: adj ist kanonisch); die erzeugte Algebra ist M_4(C) (Dimension {dim_cl}), der Modul irreduzibel",
      ok_cl and dim_cl == 16)

ok_lift, ok_char, ok_ctrl = True, True, True
for _ in range(6):
    g = rsl2()
    rows = []
    for a in EINH:
        a = a.astype(complex)
        Ga, Gb = gam(a), gam(g @ a @ g.conj().T)
        rows.append(np.kron(Ga.T, np.eye(4)) - np.kron(np.eye(4), Gb))
    nul = nullraum(np.vstack(rows))
    ok_lift &= nul.shape[0] == 1
    S = nul[0].reshape(4, 4, order="F")
    A, B = S[:2, :2], S[2:, 2:]
    k = A.ravel()[np.argmax(np.abs(g.ravel()))] / g.ravel()[np.argmax(np.abs(g.ravel()))]
    ok_lift &= np.allclose(A, k * g) and np.allclose(B, k * np.linalg.inv(g.conj().T))
    ok_lift &= np.abs(S[:2, 2:]).max() < 1e-9 and np.abs(S[2:, :2]).max() < 1e-9
    trB = np.trace(B / k)
    ok_char &= abs(trB - np.conj(np.trace(g))) < 1e-9 and min(abs(trB - np.trace(g)), abs(trB + np.trace(g))) > 1e-3
    St = np.block([[g, np.zeros((2, 2))], [np.zeros((2, 2)), g]])
    ok_ctrl &= max(np.abs(St @ gam(a) - gam(g @ a @ g.conj().T) @ St).max() for a in EINH) > 1e-3
check("Spinhebung der Lorentzwirkung a -> g a g^+ (6 zufällige g in SL(2,C)): der Verflechter ist bis auf Skala eindeutig, "
      "S = diag(g, (g^+)^-1); der Nahtplatz trägt g (Linkswirkung von M auf V), der Partnerplatz (g^+)^-1 mit Charakter "
      "conj(tr g) != tr g, also V-bar^* und nicht V; Kontrolle: diag(g, g) verflicht nicht", ok_lift and ok_char and ok_ctrl)

# --- 1.5 Tomita-Partner = J-Bild
ok_vec, ok_eps = True, True
gs = [rsl2() for _ in range(4)]
for g in gs:
    a = crand(2, 2)
    ok_vec &= np.allclose((g @ a @ g.conj().T).ravel(order="F"), np.kron(g.conj(), g) @ a.ravel(order="F"))
    ok_eps &= np.allclose(EPS @ g.conj() @ np.linalg.inv(EPS), np.linalg.inv(g.conj().T))


def verflechter_dim(rep1, rep2, gs_):
    rows_ = [np.kron(rep1(g).T, np.eye(rep2(g).shape[0])) - np.kron(np.eye(rep1(g).shape[0]), rep2(g)) for g in gs_]
    return nullraum(np.vstack(rows_)).shape[0]


d_bar_bstern = verflechter_dim(lambda g: g.conj(), lambda g: np.linalg.inv(g.conj().T), gs)
d_bar_v = verflechter_dim(lambda g: g.conj(), lambda g: g, gs)
ok_ideal = np.allclose(np.array([[1, 0], [2, 0]]).conj().T, np.array([[1, 2], [0, 0]]))
check(f"H = V (x) V-bar: vec(g a g^+) = (conj g (x) g) vec a, der J-Partnerfaktor trägt conj g; J bildet das Linksideal M e "
      f"(Spalten, = V) auf das Rechtsideal e M ab; Hom(V-bar, V-bar^*) hat Dimension {d_bar_bstern} (eps-bar), "
      f"Hom(V-bar, V) Dimension {d_bar_v}: Das konjugierte Partnerdublett V-bar^* ist das J-Bild, dualisiert über eps-bar",
      ok_vec and ok_eps and d_bar_bstern == 1 and d_bar_v == 0 and ok_ideal)


# --- 1.6 Twistorform Sigma
def herm_basis(n):
    B = []
    for i in range(n):
        M = np.zeros((n, n), complex)
        M[i, i] = 1
        B.append(M)
    for i in range(n):
        for j in range(i + 1, n):
            M = np.zeros((n, n), complex)
            M[i, j] = M[j, i] = 1
            B.append(M)
            M = np.zeros((n, n), complex)
            M[i, j], M[j, i] = -1j, 1j
            B.append(M)
    return B


HB4 = herm_basis(4)


def S_von(g):
    return np.block([[g, np.zeros((2, 2))], [np.zeros((2, 2)), np.linalg.inv(g.conj().T)]])


rows = []
for g in gs:
    S = S_von(g)
    rows.append(np.array([np.concatenate([(S.conj().T @ B @ S - B).real.ravel(), (S.conj().T @ B @ S - B).imag.ravel()])
                          for B in HB4]).T)
M_inv = np.vstack(rows)
dim_inv = nullraum(M_inv).shape[0]
rows_n = []
for _ in range(12):
    X = crand(2, 2)
    x = X + X.conj().T
    u = crand(2)
    Zt = np.concatenate([u, x @ u])
    rows_n.append([np.real(Zt.conj() @ B @ Zt) for B in HB4])
nul_n = nullraum(np.vstack([M_inv, np.array(rows_n)])).real
dim_n = nul_n.shape[0]
Gs = sum(c * B for c, B in zip(nul_n[0], HB4))
Gs = Gs / np.abs(Gs).max()


def sigma_auf_gerade(G, x, n=30):
    vals = []
    for _ in range(n):
        u = crand(2)
        Zt = np.concatenate([u, x @ u])
        vals.append(np.real(Zt.conj() @ G @ Zt))
    return np.array(vals)


Y = crand(2, 2)
ypos = Y @ Y.conj().T + 0.1 * np.eye(2)
xr = (lambda X: (X + X.conj().T) / 2)(crand(2, 2))
if sigma_auf_gerade(Gs, xr + 1j * ypos).min() < 0:
    Gs = -Gs                                                             # Vorzeichen aus dem Kegel P
v_vor, v_rueck = sigma_auf_gerade(Gs, xr + 1j * ypos), sigma_auf_gerade(Gs, xr - 1j * ypos)
v_reell = sigma_auf_gerade(Gs, xr)
c_ist = Gs[2, 0]
ew_G = np.linalg.eigvalsh(Gs)
G_falsch = np.block([[np.zeros((2, 2)), I2], [I2, np.zeros((2, 2))]])    # Phase c reell statt imaginär
v_falsch = sigma_auf_gerade(G_falsch, xr)
ZAHLEN["twistorform"] = {"dim_invariant": dim_inv, "dim_mit_J": dim_n, "c": [float(c_ist.real), float(c_ist.imag)],
                         "eigenwerte": [float(e) for e in ew_G]}
check(f"Sigma: SL(2,C)-invariante Hermitesche Formen auf V + V-bar^* bilden einen {dim_inv}-dimensionalen Raum "
      f"[[0, conj(c) I], [c I, 0]]; dass die J-reellen Geraden null sind, lässt {dim_n} Dimension (c rein imaginär); das "
      f"Vorzeichen fixiert der Kegel (Vorwärtsröhre Im x in P: Sigma > 0, Rückwärtsröhre < 0); Signatur (2,2); Kontrolle: "
      f"reelles c macht die J-reellen Geraden nicht null (max |Sigma| = {np.abs(v_falsch).max():.2f})",
      dim_inv == 2 and dim_n == 1 and abs(c_ist.real) < 1e-9 and v_vor.min() > 0 and v_rueck.max() < 0
      and np.abs(v_reell).max() < 1e-9 and np.sum(ew_G > 0) == 2 and np.sum(ew_G < 0) == 2 and np.abs(v_falsch).max() > 1e-2)
ok_tr, ok_tc = True, True
for _ in range(10):
    b = (lambda X: X + X.conj().T)(crand(2, 2))
    T = np.block([[I2, np.zeros((2, 2))], [b, I2]])
    ok_tr &= np.allclose(T.conj().T @ Gs @ T, Gs)
    Ti = np.block([[I2, np.zeros((2, 2))], [1j * b, I2]])
    ok_tc &= not np.allclose(Ti.conj().T @ Gs @ Ti, Gs)
check("die J-reellen Translationen (b in Herm(2)) erhalten Sigma, die imaginären nicht: Die reellen Orte sind genau die "
      "Sigma-erhaltenden Verschiebungen der Gerade", ok_tr and ok_tc)

# --- 1.7 Geraden L_x, x = adj(a), und die Naht L_0 = P(V)
ok_adj, ok_meet = True, True
for _ in range(20):
    g, a = rsl2(), crand(2, 2)
    ok_adj &= np.allclose(adj(g @ a @ g.conj().T), np.linalg.inv(g.conj().T) @ adj(a) @ np.linalg.inv(g))
    ok_adj &= np.allclose(adj(a.conj().T), adj(a).conj().T) and abs(np.linalg.det(adj(a)) - np.linalg.det(a)) < 1e-10
    p = (lambda X: X @ X.conj().T)(crand(2, 2))
    ok_adj &= np.linalg.eigvalsh(adj(p)).min() > -1e-12
    u = crand(2)
    x = adj(a)
    up = g @ u
    ok_adj &= np.allclose(S_von(g) @ np.concatenate([u, x @ u]),
                          np.concatenate([up, (np.linalg.inv(g.conj().T) @ x @ np.linalg.inv(g)) @ up]))
    y1 = x + np.outer(crand(2), crand(2))                               # det(x - y1) = 0
    y2 = x + crand(2, 2)
    r1 = np.linalg.matrix_rank(np.block([[I2, I2], [x, y1]]), tol=1e-9)
    r2 = np.linalg.matrix_rank(np.block([[I2, I2], [x, y2]]), tol=1e-9)
    ok_meet &= r1 == 3 and r2 == 4
L0_null = np.abs(sigma_auf_gerade(Gs, np.zeros((2, 2)))).max() < 1e-12
check("x = adj(a) verflicht die GNS-Wirkung a -> g a g^+ mit x -> (g^+)^-1 x g^-1 auf Hom(V, V-bar^*) = H0(O(1) (x) V-bar^*), "
      "erhält J-Realität, Kegel und det; S(g) bildet L_x auf L_x' ab; L_x und L_y schneiden sich genau für det(x - y) = 0; "
      "die Naht L_0 = P(V) ist Sigma-null, ihre Lagen sind Hom(V, V-bar^*): Normalenbündel O(1) (x) V-bar^*",
      ok_adj and ok_meet and L0_null)


# --- 1.8 Was bleibt: die Identifikation I_pos
def rep_H(g):
    return np.kron(g.conj(), g)


def rep_Hom(rW):
    return lambda g: np.kron(np.linalg.inv(g).T, rW(g))


kand_W = {"V-bar^*": lambda g: np.linalg.inv(g.conj().T), "V-bar": lambda g: g.conj(), "V": lambda g: g,
          "C^2 (trivial)": lambda g: np.eye(2, dtype=complex)}
gsH = [rsl2() for _ in range(4)]
hom_dims = {k: verflechter_dim(rep_H, rep_Hom(rW), gsH) for k, rW in kand_W.items()}
hom_dims_dual = {k: verflechter_dim(lambda g: np.linalg.inv(rep_H(g)).T, rep_Hom(rW), gsH) for k, rW in kand_W.items()}
ZAHLEN["Hom_SL2(H, Hom(V,W))"] = hom_dims
ZAHLEN["Hom_SL2(H^*, Hom(V,W))"] = hom_dims_dual
check(f"äquivariante Abbildungen GNS-Raum -> Nahtlagen Hom(V, W): {hom_dims} (vom Dualraum H^*: {hom_dims_dual}); nur das "
      "konjugierte Partnerdublett (V-bar^* = V-bar über eps-bar) lässt die Nahtlagen als GNS-Vektoren lesen: T' <=> I_pos",
      hom_dims == {"V-bar^*": 1, "V-bar": 1, "V": 0, "C^2 (trivial)": 0} and hom_dims_dual == hom_dims)
befund("Die Standardform (H, J, Delta, P) ist für jedes W dieselbe: Sie kennt keine Einbettung der Naht in einen umgebenden Raum. "
       "Gegenmodelle mit denselben Tomita-Daten und N = O(1) (x) V, O(1) (x) C^2 oder O(1) + O(-1) (z61) existieren. T' folgt also "
       "nicht aus der Standardform allein; es ist äquivalent zur einen Identifikation I_pos", True)
say("")
status("Tomita-Inhalt von T'", "E", "Für jeden treuen Zustand liefert die Standardform kanonisch (bis auf positive Skala): "
       "(Herm(2), det) = R^(1,3) mit Zeitorientierung P; den Spinormodul V + V-bar^* (Partner = J-Bild, dual über eps-bar); "
       "Sigma der Signatur (2,2) mit Phase aus J und Vorzeichen aus P; die Geraden L_x, x = adj(a); L_0 = P(V) mit N = O(1) (x) "
       "V-bar^*; J-reelle Punkte = Herm(2)")
status("Warum konjugiert", "E", "Lorentzsignatur verlangt eine Realstruktur, die M mit M' vertauscht; das ist J. Realstrukturen, "
       "die M erhalten (die eigenen der Naht), geben (2,2) oder (4,0)")
status("T' selbst", "C", "T' <=> I_pos: Die Lagen der Naht sind äquivariant die (dualen) GNS-Vektoren ihres Dublettzustands "
       "(Naht = L_0 im projektiven Spinormodul ihrer eigenen Standardform). Nicht frei wählbar mehr, aber nicht abgeleitet")
status("I_pos aus Antworten", "O", "dass sich die Naht des Originals in Richtungen des GNS-Raums bewegt, rekonstruiert keine "
       "Antwort (vgl. z64 'Nahtlagen als Antworten')")

# =================================================================== Teil 2
say("\nTeil 2: Windungszahl m gegen Spinbrücke, Neutralität, Realstruktur und Zeitorientierung")
MARKEN = np.array([1, 1j, -1, -1j])


def R_auf_V(m):
    rho = rho_m(m)
    return np.diag(np.diag(rho) ** (-0.5j))                             # Delta^(-i/2) = Ad(rho^(-i/2))


def bahn(m):
    R = R_auf_V(m)
    pts, v = [], np.array([1.0, 1.0], complex)
    for _ in range(8):
        z = v[1] / v[0]
        if not any(abs(z - p) < 1e-9 for p in pts):
            pts.append(z)
        v = R @ v
    return pts


# (a) Spinbrücke
paarungen = [((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))]
tab_a = {}
for m in (1, 2, 3, 4):
    pts = bahn(m)
    n = len(pts)
    eta = None
    if n == 4:
        perm = [int(np.argmin(np.abs(np.array(pts) - (R_auf_V(m)[1, 1] / R_auf_V(m)[0, 0]) * p))) for p in pts]
        inv = [p for p in paarungen if sorted(tuple(sorted((perm[a], perm[b]))) for a, b in p) == sorted(
            tuple(sorted(x)) for x in p)]
        eta = [[[complex(np.round(pts[i], 9)) for i in x] for x in p] for p in inv]
    tab_a[m] = {"marken": n, "geschlecht": (n // 2 - 1) if n % 2 == 0 else None, "eta": eta,
                "a_erfuellt": bool(n == 4 and eta is not None and len(eta) == 1)}
ZAHLEN["(a)_spinbruecke"] = {str(k): {kk: (str(vv) if kk == "eta" else vv) for kk, vv in v.items()} for k, v in tab_a.items()}
set_a = {m for m, v in tab_a.items() if v["a_erfuellt"]}
ok_anti = all(any(abs(pr[0][0] + pr[0][1]) < 1e-9 for pr in [tab_a[m]["eta"][0]]) for m in set_a)
check(f"(a) R = Delta^(-i/2) auf P(V_m): Markenbahn der Länge {[tab_a[m]['marken'] for m in (1, 2, 3, 4)]} für m = 1..4; vier "
      f"getrennte Verzweigungspunkte (elliptische Doppelüberlagerung, Geschlecht 1) und genau eine R-invariante Paarung, die "
      f"antipodale (eta), nur für m in {sorted(set_a)}", set_a == {1, 3} and ok_anti)


# (b) Neutralität und der Satz 'R ist Ladungsdrehung'
def fock(alpha, rmin=-5, rmax=4):
    mo = [r + alpha for r in range(rmin, rmax + 1)]
    out = []
    for occ in itertools.product((0, 1), repeat=len(mo)):
        dl0 = sum(x for x, n in zip(mo, occ) if x > 0 and n) + sum(-x for x, n in zip(mo, occ) if x < 0 and not n)
        Q = sum(n for x, n in zip(mo, occ) if x > 0) - sum(1 - n for x, n in zip(mo, occ) if x < 0)
        par = (sum(n for x, n in zip(mo, occ) if x > 0) + sum(1 - n for x, n in zip(mo, occ) if x < 0)) % 2
        out.append((int(round(4 * dl0)), Q, par))
    return [s for s in out if s[0] <= 16]


sek = {al: fock(al) for al in (0.25, 0.75)}
ok_RQ = all((M - (2 * Q * Q + (4 * al - 2) * Q)) % 4 == 0 for al, st in sek.items() for M, Q, p in st)
ok_par = all(p == Q % 2 for st in sek.values() for M, Q, p in st)
einsektor = {m: {al: sorted(set((Q, p) for M, Q, p in st if M == m)) for al, st in sek.items()} for m in (1, 2, 3, 4)}
zweisektor = {}
for m in (1, 2, 3, 4):
    lst = set()
    for M1 in range(m + 1):
        for (Q1, p1) in set((Q, p) for M, Q, p in sek[0.25] if M == M1):
            for (Q2, p2) in set((Q, p) for M, Q, p in sek[0.75] if M == m - M1):
                lst.add((Q1, Q2, (p1 + p2) % 2))
    zweisektor[m] = sorted(lst)
set_b = {m for m in (1, 2, 3, 4) if (0, 0, 0) in zweisektor[m]}
ZAHLEN["(b)_ladungen_einsektor"] = {str(m): {str(a): v for a, v in d.items()} for m, d in einsektor.items()}
ZAHLEN["(b)_ladungen_zweisektor_(Q_1/4, Q_3/4, Parität)"] = {str(m): v for m, v in zweisektor.items()}
neutrale_paare = [(s, t) for st in sek.values() for s, t in itertools.combinations(st, 2) if s[1] == t[1]]
R_skalar = all((1j ** s[0]) / (1j ** t[0]) == 1 for s, t in neutrale_paare)
check(f"(b) Fockzustände beider Quartersektoren bis Marke 16: Marke M und Ladung Q erfüllen M = 2Q^2 + (4 alpha - 2) Q mod 4, "
      f"also R = i^M = Funktion von Q allein, und die Parität ist Q mod 2; Marke m = 1, 2, 3, 4 hat die Ladungen "
      f"(Q_1/4, Q_3/4, Parität) {zweisektor}: neutral nur m in {sorted(set_b)}; auf allen {len(neutrale_paare)} Paaren gleicher "
      "Ladung ist R skalar", ok_RQ and ok_par and set_b == {4} and R_skalar)
befund("Satz: Die Vierteldrehung R = Delta^(-i/2) ist auf jedem Quartersektor eine Ladungsdrehung (mal Fermionparität). Auf einem "
       "Dublett mit neutralen Matrixeinheiten wirkt sie skalar, die Marken fallen zusammen. (a) und (b) schließen sich für jedes "
       "m aus, nicht nur für m <= 4: Das erklärt den OriginLab-Befund 'kein neutrales Dublett trennt die Marken'", True)

# (c) Realstruktur und J
ok_cL = ok_J                                                            # J aus S: für alle rho_m dasselbe a -> a^+
vor_m = {}
for m in (1, 2, 3, 4):
    xi = msqrt(rho_m(m))
    vor_m[m] = bool(np.all(np.linalg.eigvalsh(adj(xi)) > 0)) and sigma_auf_gerade(Gs, 1j * adj(xi)).min() > 0
jq = lambda v: np.array([-np.conj(v[1]), np.conj(v[0])])                 # J mit Hardy-Form und eps: quaternionisch
sigma = lambda z: -1 / np.conj(z)
tab_c = {}
for m in (1, 2, 3, 4):
    pts = bahn(m)
    inv_set = all(any(abs(sigma(z) - p) < 1e-9 for p in pts) for z in pts)
    paar_eta = None
    if len(pts) == 4:
        paar_eta = all(abs(sigma(z) + z) < 1e-9 for z in pts)
    xi = msqrt(rho_m(m))
    a_, b_ = xi[0, 0].real, xi[1, 1].real
    n1 = np.outer([1, 1], [1, 1]).astype(complex)
    lam = np.trace(adj(xi) @ n1).real / np.linalg.det(xi).real
    n2 = lam * xi - n1
    w_, v_ = np.linalg.eigh(n2)
    vv = v_[:, -1]
    z_rest = vv[1] / vv[0]
    tab_c[m] = {"marken_J_invariant": bool(inv_set), "J_paarung_gleich_eta": paar_eta,
                "ruhesystem_antipode_von_1": [float(z_rest.real), float(z_rest.imag)],
                "formel_-(b/a)": float(-(b_ / a_)), "rapiditaet": float(np.arctanh((a_ - b_) / (a_ + b_)))}
ZAHLEN["(c)_realstruktur"] = {str(k): v for k, v in tab_c.items()}
set_cE = {m for m, v in tab_c.items() if v["marken_J_invariant"]}
set_eta = {m for m, v in tab_c.items() if v["J_paarung_gleich_eta"]}
su2 = [(lambda u: u / np.sqrt(np.linalg.det(u)))(sla.expm(1j * (X + X.conj().T))) for X in crand(3, 2, 2)]
ok_q = all(np.allclose(jq(jq(v)), -v) for v in crand(5, 2)) and all(
    np.allclose(jq(u @ v), u @ jq(v)) for u in su2 for v in crand(2, 2))
ok_q &= not np.allclose(jq(np.diag([2.0, 0.5]) @ np.array([1.0, 1j])), np.diag([2.0, 0.5]) @ jq(np.array([1.0, 1j])))
ok_rest = all(abs(tab_c[m]["ruhesystem_antipode_von_1"][0] - tab_c[m]["formel_-(b/a)"]) < 1e-9
              and abs(tab_c[m]["ruhesystem_antipode_von_1"][0] + 1) > 0.1 for m in (1, 2, 3, 4))
check(f"(c) Lorentzisch: J, Kegel und Sigma hängen nicht von m ab (Teil 1), jede Naht L_0 ist Sigma-null, keine Bedingung. "
      f"Euklidisch: J mit Hardy-Form und eps ist quaternionisch (j^2 = -1, SU(2)-äquivariant), auf P(V_m) die Antipodenabbildung; "
      f"sie erhält die Markenmenge nur für m in {sorted(set_cE)} (ein einzelner Punkt ist nie antipodal invariant) und induziert "
      f"eta nur für m in {sorted(set_eta)}; im Ruhesystem von sqrt(rho_m) geht die Marke 1 auf -e^(-pi m/2), für kein m >= 1 "
      "eine Marke", ok_cL and ok_q and set_cE == {1, 2, 3} and set_eta == {1, 3} and ok_rest)

# (d) Zeitorientierung
set_d = {m for m, v in vor_m.items() if v}
check(f"(d) sqrt(rho_m) liegt für m = 1..4 im Innern des Zukunftskegels, KMS bei beta = +1, und die Gerade zu i adj(sqrt rho_m) "
      f"liegt in PT^+ (Sigma > 0): erfüllt für m in {sorted(set_d)}, keine Bedingung", set_d == {1, 2, 3, 4})

# Uhrenquotient
uhr = {}
for m in (1, 2, 3, 4):
    deck = np.exp(2j * np.pi / m)
    uhr[m] = bool(all(any(abs(deck * z - p) < 1e-9 for p in MARKEN) for z in MARKEN))
set_uhr = {m for m, v in uhr.items() if v}
check(f"Wörterbuch z = S^m: die Deckgruppe mu_m erhält die Marken {{S^4 = 1}} genau für m | 4, also m in {sorted(set_uhr)}; für "
      "m = 3 identifiziert der Quotient Marken mit Nichtmarken", set_uhr == {1, 2, 4})
alle = {1, 2, 3, 4}
schnitt_abcd = set_a & set_b & alle & set_d
schnitt_ohne_b = set_a & set_d & set_cE
schnitt_uhr = schnitt_ohne_b & set_uhr
ZAHLEN["m_mengen"] = {"(a)": sorted(set_a), "(b)": sorted(set_b), "(c) lorentzisch": sorted(alle), "(c) euklidisch": sorted(set_cE),
                      "(d)": sorted(set_d), "Uhrenquotient m|4": sorted(set_uhr), "(a)(b)(c)(d)": sorted(schnitt_abcd),
                      "(a)(c)(d) ohne (b)": sorted(schnitt_ohne_b), "dazu m|4": sorted(schnitt_uhr)}
check(f"Bilanz m: (a) {sorted(set_a)}, (b) {sorted(set_b)}, (c) lorentzisch alle / euklidisch {sorted(set_cE)}, (d) "
      f"{sorted(set_d)}: alle vier zusammen {sorted(schnitt_abcd)} (unverträglich); ohne (b) {sorted(schnitt_ohne_b)} "
      f"(Gleichstand); mit dem Uhrenquotienten m | 4 {sorted(schnitt_uhr)}",
      schnitt_abcd == set() and schnitt_ohne_b == {1, 3} and schnitt_uhr == {1})
S4_rho1 = None
r1 = rho_m(1)
big = r1
for _ in range(3):
    big = np.kron(big, r1)
Bm = np.zeros((16, 5))
for idx in range(16):
    Bm[idx, bin(idx).count("1")] = 1
Bm = Bm / np.linalg.norm(Bm, axis=0)
S4 = Bm.T @ big.real @ Bm
e04 = np.diag([S4[0, 0], S4[4, 4]]) / (S4[0, 0] + S4[4, 4])
check("Verträglichkeit der Ebenen: Sym^4 rho_1 auf dem R-invarianten Dublett {1, zeta^4} des Trägers ist rho_4 (z63), und der "
      "Markenzustand 1 - zeta^4 liegt in diesem neutralen Dublett: m = 1 (Naht) und m = 4 (neutrales Trägerdublett) sind "
      "zwei Stufen, keine Konkurrenten", np.allclose(e04, rho_m(4).real))
say("")
status("(a) gegen (b)", "E", "unverträglich für jedes m: R = Delta^(-i/2) ist auf jedem Quartersektor eine Ladungsdrehung mal "
       "Parität, also auf neutralen Dubletts skalar")
status("(c), (d)", "E", "lorentzisch keine Bedingung (Tomita-Daten zustandsunabhängig, jedes treue rho_m zeitartig in der "
       "Zukunft); euklidisch (Antipoden) schließt m = 4 aus")
status("m = 1", "C", "aus (a) + (c) + (d) und dem Wörterbuch z = S^m als Quotient nach einer Uhrenuntergruppe (m | 4); "
       "Neutralität gilt dann für das Trägerdublett {1, zeta^4} (rho_4 = Sym^4 rho_1), nicht für den Nahtspinor")
status("Neutralität des Nahtspinors", "L", "als Forderung an den Spinor V selbst unverträglich mit der Spinbrücke; als "
       "Forderung an Observable (gerade Bilineare, Trägerdublett) erfüllt")

# =================================================================== Teil 3
say("\nTeil 3: Eindeutigkeit des Kerns auf der Matrixkugel und nu")


# --- 3.1 Native Auslesungen des Originals
def c_mp(x, beta):
    return 1 / (1 + mp.exp(beta * x))


def cb_mp(x, beta):
    return 1 / (1 + mp.exp(-beta * x))


def strom(m, beta, al=mp.mpf(1) / 4, R1=40):
    a = mp.fsum(c_mp(r + al + m, beta) * cb_mp(r + al, beta) for r in range(-R1, R1))       # omega(J_-m J_m)
    b = mp.fsum(c_mp(r + al, beta) * cb_mp(r + al + m, beta) for r in range(-R1, R1))       # omega(J_m J_-m)
    return a, b


B0 = 4 * mp.pi
ok_schw, ok_kms, ok_bose, ok_vak = True, True, True, True
for m in range(1, 6):
    a, b = strom(m, B0)
    ok_schw &= abs(b - a - m) < mp.mpf(10) ** -30
    ok_kms &= abs(a / b * mp.exp(4 * mp.pi * m) - 1) < mp.mpf(10) ** -25
    ok_bose &= abs(a - m / (mp.exp(4 * mp.pi * m) - 1)) < mp.mpf(10) ** -30 * m
    av, bv = strom(m, 100 * B0)
    ok_vak &= abs(bv - m) < mp.mpf(10) ** -30 and abs(av) < mp.mpf(10) ** -30
check("Stromauslesung aus den Fermisummen des Originals (m = 1..5): Schwinger-Term omega([J_m, J_-m]) = m (Niveau 1), KMS "
      "omega(J_-m J_m) = e^(-4 pi m) omega(J_m J_-m), also omega(J_-m J_m) = m/(e^(4 pi m) - 1); im Vakuumlimes "
      "<J_m J_-m> = m: Stromkern sum m t^(m-1) = (1 - t)^-2; die Fermiauslesung ist der Hardy-Kern sum t^k = (1 - t)^-1",
      ok_schw and ok_kms and ok_bose and ok_vak)

# --- 3.2 Faraut-Koranyi und Fock-Normierung
nu, t1, t2, tt = sp.symbols("nu t1 t2 t")


def schur2(m1, m2, x, y):
    return sp.expand(sp.cancel((x * y) ** m2 * (x ** (m1 - m2 + 1) - y ** (m1 - m2 + 1)) / (x - y)))


def hook(m1, m2):
    H = 1
    for j in range(1, m1 + 1):
        H *= m1 - j + 1 + (1 if j <= m2 else 0)
    for j in range(1, m2 + 1):
        H *= m2 - j + 1
    return H


def poch(x, n):
    r = sp.Integer(1)
    for i in range(n):
        r *= x + i
    return r


def bis_grad(expr, D):
    P = sp.Poly(sp.expand(expr), t1, t2)
    return sum(c * t1 ** i * t2 ** j for (i, j), c in P.terms() if i + j <= D)


DG = 6
lhs = bis_grad(sp.series(sp.series((1 - t1) ** (-nu) * (1 - t2) ** (-nu), t1, 0, DG + 1).removeO(), t2, 0, DG + 1).removeO(), DG)
rhs = sum(poch(nu, k - m2) * poch(nu - 1, m2) * schur2(k - m2, m2, t1, t2) / hook(k - m2, m2)
          for k in range(DG + 1) for m2 in range(k // 2 + 1))
fock_l = bis_grad(sp.series(sp.series(sp.exp(t1 + t2), t1, 0, DG + 1).removeO(), t2, 0, DG + 1).removeO(), DG)
fock_r = sum(schur2(k - m2, m2, t1, t2) / hook(k - m2, m2) for k in range(DG + 1) for m2 in range(k // 2 + 1))
check("Faraut-Koranyi bis Grad 6 symbolisch in nu: det(I - X)^-nu = sum_m (nu)_m1 (nu - 1)_m2 s_m(X)/H(m), und e^(tr X) = "
      "sum_m s_m(X)/H(m): K_m = s_m(Z W^+)/H(m) sind die Fischer-Reproduktionskerne der Schmid-Typen",
      sp.expand(lhs - rhs) == 0 and sp.expand(fock_l - fock_r) == 0)

# --- 3.3 U(2) x U(2) allein: ein Kegel, kein Kern
a1, a2, b1, b2 = sp.symbols("a1 a2 b1 b2")
v4 = [a1 * b1, a1 * b2, a2 * b1, a2 * b2]
ok_cauchy, zahl_typen = True, {}
for k in range(1, 7):
    hk = sp.expand(sum(sp.prod(c) for c in itertools.combinations_with_replacement(v4, k)))
    rhs_k = sp.expand(sum(schur2(k - m2, m2, a1, a2) * schur2(k - m2, m2, b1, b2) for m2 in range(k // 2 + 1)))
    ok_cauchy &= sp.expand(hk - rhs_k) == 0
    zahl_typen[k] = k // 2 + 1
ZAHLEN["schmid_typen_je_grad"] = zahl_typen
check(f"Cauchy: der Charakter von P^k(M_2) unter U(2) x U(2) ist sum s_m(a) s_m(b) (k = 1..6), also multiplizitätsfrei mit "
      f"{list(zahl_typen.values())} Typen: invariante positive Kerne bilden den Kegel {{sum c_m K_m, c_m >= 0}}, unendlich-"
      "dimensional; die det^-nu-Familie ist darin nur eine Kurve", ok_cauchy)

# --- 3.4 Die det^-nu-Familie: zwei Auslesungen, ein nu
r1k = [sp.simplify(poch(nu, k) / factorial(k)) for k in range(1, 9)]
sk = [sp.simplify(poch(2 * nu, k) / factorial(k)) for k in range(1, 9)]
sol_h = set(sp.solve(sp.Eq(r1k[0], 1), nu)) & set(sp.solve(sp.Eq(r1k[1], 1), nu))
sol_s = set(sp.solve(sp.Eq(sk[0], 2), nu)) & set(sp.solve(sp.Eq(sk[1], 3), nu))
ok_alle = all(sp.simplify(r1k[k].subs(nu, 1) - 1) == 0 and sp.simplify(sk[k].subs(nu, 1) - (k + 2)) == 0 for k in range(8))
wal = {}
for nv in (sp.Rational(1, 4), sp.Rational(1, 2), sp.Rational(3, 4), 1, sp.Rational(3, 2), 2, 3):
    koeff = [poch(sp.Integer(1) * nv, k - m2) * poch(nv - 1, m2) for k in range(9) for m2 in range(k // 2 + 1)]
    wal[str(nv)] = bool(min(koeff) >= 0)
ZAHLEN["wallach_positiv"] = wal
check(f"det(I - Z W^+)^-nu: Rang-1-Schnitt (1 - t)^-nu, Skalarschnitt Z = z I (1 - t)^-2nu; Hardy-Auslesung gibt nu in "
      f"{sorted(sol_h)}, Stromauslesung nu in {sorted(sol_s)}, beide bis Grad 8 erfüllt; Positivität (alle (nu)_m >= 0) gilt "
      f"nur auf der Wallach-Menge: {wal}",
      sol_h == {1} and sol_s == {1} and ok_alle
      and wal == {"1/4": False, "1/2": False, "3/4": False, "1": True, "3/2": True, "2": True, "3": True})

# --- 3.5 Stärker: Eindeutigkeit im ganzen invarianten Kegel (LP)
DM = 10
teile = [(k - m2, m2) for k in range(DM + 1) for m2 in range(k // 2 + 1)]
idx = {p: i for i, p in enumerate(teile)}
Aeq, beq = [], []
for k in range(DM + 1):
    row = np.zeros(len(teile))
    row[idx[(k, 0)]] = 1 / factorial(k)
    Aeq.append(row)
    beq.append(1.0)                                                     # Hardy: Koeffizient von t^k auf dem Rang-1-Schnitt
    row = np.zeros(len(teile))
    for m2 in range(k // 2 + 1):
        row[idx[(k - m2, m2)]] = (k - 2 * m2 + 1) / hook(k - m2, m2)
    Aeq.append(row)
    beq.append(k + 1.0)                                                 # Strom: Koeffizient von t^k auf dem Skalarschnitt
Aeq, beq = np.array(Aeq), np.array(beq)
bounds = [(0, None)] * len(teile)
ok_lp = True
for i, p in enumerate(teile):
    c = np.zeros(len(teile))
    c[i] = 1
    lo = linprog(c, A_eq=Aeq, b_eq=beq, bounds=bounds)
    hi = linprog(-c, A_eq=Aeq, b_eq=beq, bounds=bounds)
    ziel = float(factorial(p[0])) if p[1] == 0 else 0.0                 # (1)_m = m1! delta_(m2,0)
    ok_lp &= lo.status == 0 and hi.status == 0 and abs(lo.fun - ziel) < 1e-6 and abs(-hi.fun - ziel) < 1e-6
c = np.zeros(len(teile))
c[idx[(1, 1)]] = -1
ohne_strom = linprog(c, A_eq=Aeq[0::2], b_eq=beq[0::2], bounds=bounds)
c = np.zeros(len(teile))
c[idx[(2, 0)]] = 1
ohne_hardy = (linprog(c, A_eq=Aeq[1::2], b_eq=beq[1::2], bounds=bounds).fun,
              -linprog(-c, A_eq=Aeq[1::2], b_eq=beq[1::2], bounds=bounds).fun)
ohne_pos = nullraum(Aeq).shape[0]
ZAHLEN["kern_lp"] = {"teile": len(teile), "ohne_strom_status": int(ohne_strom.status), "ohne_hardy_c20": list(ohne_hardy),
                     "ohne_positivitaet_freie_richtungen": int(ohne_pos)}
check(f"im ganzen U(2) x U(2)-invarianten Kegel (alle {len(teile)} Typen bis Grad {DM}) legen Hardy- und Stromauslesung mit "
      "c_m >= 0 jeden Koeffizienten eindeutig fest: c_m = (1)_m, also K = det(I - Z W^+)^-1 (Sättigung: der Typ (k,0) allein "
      "gibt schon k + 1); Kontrollen: ohne Strom unbeschränkt (Status "
      f"{ohne_strom.status}), ohne Hardy c_(2,0) in [{ohne_hardy[0]:.1f}, {ohne_hardy[1]:.1f}], ohne Positivität {ohne_pos} "
      "freie Richtungen",
      ok_lp and ohne_strom.status == 3 and ohne_hardy[1] - ohne_hardy[0] > 1 and ohne_pos > 0)


# --- 3.6 SU(2,2)-Kovarianz folgt, sie wird nicht gebraucht
def su22():
    A = crand(2, 2)
    D = crand(2, 2)
    B = crand(2, 2)
    X = np.block([[(A - A.conj().T) / 2, B], [B.conj().T, (D - D.conj().T) / 2]]) * 0.4
    return sla.expm(X)


def kugel():
    Z = crand(2, 2)
    return Z / np.linalg.norm(Z, 2) * rng.uniform(0.2, 0.9)


def moeb(g, Z):
    a, b, c_, d = g[:2, :2], g[:2, 2:], g[2:, :2], g[2:, 2:]
    return (a @ Z + b) @ np.linalg.inv(c_ @ Z + d), np.linalg.det(c_ @ Z + d)


def K1(Z, W):
    return 1 / np.linalg.det(np.eye(2) - Z @ W.conj().T)


def Kst(Z, W, e=0.3):
    return K1(Z, W) + e * np.linalg.det(Z) * np.conj(np.linalg.det(W))


ok_cov, ok_kug, def_K1, def_Kst = True, True, [], []
eta_ = np.diag([1, 1, -1, -1])
for _ in range(10):
    g = su22()
    ok_cov &= np.allclose(g.conj().T @ eta_ @ g, eta_)
    Zs_ = [kugel() for _ in range(4)]
    im = [moeb(g, Z) for Z in Zs_]
    ok_kug &= all(np.linalg.norm(gz, 2) < 1 for gz, _ in im)
    for (Z, (gZ, jZ)), (W, (gW, jW)) in itertools.product(zip(Zs_, im), repeat=2):
        ok_cov &= np.isclose(K1(gZ, gW), jZ * K1(Z, W) * np.conj(jW))
    Rk = lambda K, i, j: K(im[i][0], im[j][0]) / K(Zs_[i], Zs_[j])
    for K, store in ((K1, def_K1), (Kst, def_Kst)):
        store.append(abs(Rk(K, 0, 1) * Rk(K, 2, 3) - Rk(K, 0, 3) * Rk(K, 2, 1)) / abs(Rk(K, 0, 1) * Rk(K, 2, 3)))
ZAHLEN["produktform_defekt"] = {"det^-1_max": float(max(def_K1)), "det^-1+0.3det_min": float(min(def_Kst))}
check("det(I - Z W^+)^-1 ist SU(2,2)-kovariant mit Multiplikator det(CZ + D) (10 zufällige g, Kugel invariant, Produktform-"
      f"defekt max {max(def_K1):.1e}); der ebenfalls positive, U(2) x U(2)-invariante Kern det^-1 + 0,3 det Z conj(det W) hat "
      f"keinen Multiplikator (Defekt min {min(def_Kst):.1e}): Kovarianz unter der konformen Gruppe ist Folge der Auslesungen, "
      "keine Prämisse", ok_cov and ok_kug and max(def_K1) < 1e-10 and min(def_Kst) > 1e-6)

# --- 3.7 Masselos: det(d_Z) und die Leiterdarstellung
zs = sp.symbols("z11 z12 z21 z22")
wb = sp.symbols("w11 w12 w21 w22")
Zm = sp.Matrix([[zs[0], zs[1]], [zs[2], zs[3]]])
Wd = sp.Matrix([[wb[0], wb[2]], [wb[1], wb[3]]])                        # W^+ (Einträge conj(w_ji))
Adet = sp.expand((sp.eye(2) - Zm * Wd).det())
Q1 = sp.expand(sp.diff(Adet, zs[0]) * sp.diff(Adet, zs[3]) - sp.diff(Adet, zs[1]) * sp.diff(Adet, zs[2]))
Q2 = sp.expand(sp.diff(Adet, zs[0], zs[3]) - sp.diff(Adet, zs[1], zs[2]))
dW = sp.expand(Wd.det())
ok_cay = sp.expand(Q1 - dW * Adet) == 0 and sp.expand(Q2 - 2 * dW) == 0
# det(d) A^-nu = nu(nu+1) A^(-nu-2) Q1 - nu A^(-nu-1) Q2 = nu(nu-1) det(W^+) A^(-nu-1)


def detd(f):
    return sp.expand(sp.diff(f, zs[0], zs[3]) - sp.diff(f, zs[1], zs[2]))


Wn = sp.Matrix([[sp.Rational(1, 3), sp.Rational(-1, 5)], [sp.Rational(2, 7), sp.Rational(1, 2)]])
X_ = Zm * Wn
hks = [sp.Integer(1), sp.expand(X_.trace())]
for k in range(2, 6):
    hks.append(sp.expand(X_.trace() * hks[-1] - X_.det() * hks[-2]))
ok_hk = all(detd(h) == 0 for h in hks) and detd(sp.expand(X_.det())) != 0
mono = {k: [sp.prod(c) for c in itertools.combinations_with_replacement(zs, k)] for k in range(6)}
kern_dim = {}
for k in range(2, 6):
    bild = [sp.Poly(detd(mn), *zs) for mn in mono[k]]
    basis_k2 = {mn: i for i, mn in enumerate(mono[k - 2])}
    M = sp.zeros(len(mono[k - 2]), len(mono[k]))
    for j, P in enumerate(bild):
        for mon, cf in P.terms():
            if cf == 0:
                continue
            term = sp.prod([z_ ** e for z_, e in zip(zs, mon)])
            M[basis_k2[term], j] = cf
    kern_dim[k] = len(mono[k]) - M.rank()
ZAHLEN["dim_ker_det_d"] = kern_dim
check(f"det(d_Z) det(I - Z W^+)^-nu = nu(nu - 1) det(W^+) det(I - Z W^+)^(-nu-1) (Polynomidentitäten Q1 = det(W^+) A, Q2 = "
      f"2 det(W^+)); ker det(d) auf P^k hat Dimension {kern_dim} = (k + 1)^2 = dim P_(k,0), die h_k(Z W^+) liegen darin, "
      "det(Z W^+) nicht: nu = 1 ist genau die Leiterdarstellung (nur Typen (k, 0))",
      ok_cay and ok_hk and all(kern_dim[k] == (k + 1) ** 2 for k in kern_dim))


# --- 3.8 Cayley nach Minkowski
def cayley(Z):
    return 1j * (np.eye(2) + Z) @ np.linalg.inv(np.eye(2) - Z)


ok_c1 = True
for _ in range(40):
    Z, W = kugel(), kugel()
    X, Yc = cayley(Z), cayley(W)
    ok_c1 &= np.linalg.eigvalsh((X - X.conj().T) / 2j).min() > 0
    lhs_ = np.eye(2) - Z @ W.conj().T
    rhs_ = -2j * np.linalg.inv(X + 1j * np.eye(2)) @ (X - Yc.conj().T) @ np.linalg.inv(Yc.conj().T - 1j * np.eye(2))
    ok_c1 &= np.allclose(lhs_, rhs_)
    ok_c1 &= np.isclose(K1(Z, W), -0.25 * np.linalg.det(X + 1j * np.eye(2)) * np.conj(np.linalg.det(Yc + 1j * np.eye(2)))
                        / np.linalg.det(X - Yc.conj().T))
    U = sla.expm(1j * (lambda A: (A + A.conj().T) / 2)(crand(2, 2)))
    Xu = cayley(U)
    ok_c1 &= np.allclose(Xu, Xu.conj().T)
xs, ys = sp.symbols("x0:4"), sp.symbols("y0:4")
quad = (xs[0] - ys[0]) ** 2 - sum((xs[i] - ys[i]) ** 2 for i in (1, 2, 3))
Xs = sp.Matrix([[xs[0] + xs[3], xs[1] - sp.I * xs[2]], [xs[1] + sp.I * xs[2], xs[0] - xs[3]]])
Ys = sp.Matrix([[ys[0] + ys[3], ys[1] - sp.I * ys[2]], [ys[1] + sp.I * ys[2], ys[0] - ys[3]]])
ok_q = sp.expand((Xs - Ys).det() - quad) == 0
fW = quad ** (-nu)
box = sp.diff(fW, xs[0], 2) - sum(sp.diff(fW, xs[i], 2) for i in (1, 2, 3))
box_fak = sp.simplify(box / quad ** (-nu - 1))
lam_, mu_ = sp.symbols("lambda mu")
Pn = sp.Matrix([[1, 0], [0, 0]])
nulllinie = sp.expand(((2 * sp.I) * sp.eye(2) + (lam_ - mu_) * Pn).det())
zeitachse = sp.expand(((lam_ - mu_) * sp.eye(2)).det())
check(f"Cayley X = i(I + Z)(I - Z)^-1: Kugel -> Vorwärtsröhre (Im X > 0), Shilov-Rand U(2) -> Herm(2) = R^(1,3), "
      f"I - Z W^+ = -2i (X + i)^-1 (X - Y^+) (Y^+ - i)^-1, also det(I - Z W^+)^-1 = -(1/4) det(X + i) conj det(Y + i) / "
      f"det(X - Y^+) mit det(X - Y^+) = (x - y)^2 (Minkowski); Box (x - y)^(-2 nu) = {box_fak} (x - y)^(-2 nu - 2), null "
      f"genau für nu = 0, 1; Nulllinie durch iI: det = {nulllinie} (Grad 1, Szegő), Zeitachse: det = {zeitachse} (Grad 2, Strom)",
      ok_c1 and ok_q and sp.simplify(box_fak - 4 * nu * (nu - 1)) == 0
      and sp.degree(nulllinie, lam_) == 1 and sp.degree(zeitachse, lam_) == 2)
say("")
status("Kern", "E", "Unter (K1) holomorph/antiholomorph auf der 2x2-Matrixkugel, (K2) positiv, (K3) U(2) x U(2)-invariant, "
       "(K4) Rang-1-Schnitt = Hardy-Kern (1 - t)^-1 und (K5) Skalarschnitt = Stromkern sum m t^(m-1) (Schwinger-Term m) gilt "
       "K = det(I - Z W^+)^-1, also nu = 1; SU(2,2)-Kovarianz und det(d) K = 0 folgen")
status("Notwendigkeit", "E", "(K4) allein, (K5) allein oder ohne (K2) ist der Kern nicht festgelegt (Kontrollen); U(2) x U(2) "
       "allein lässt einen unendlichdimensionalen Kegel")
status("Lesart", "C", "die Matrixkugel als Raum der Nahtlagen (Cayley-Bild der Vorwärtsröhre H_J + i P der Standardform) und die "
       "Zuordnung Rang-1-Schnitt = Lichtstrahl (Szegő), Skalarschnitt = Zeitachse (Strom) setzen T' bzw. I_pos voraus; die "
       "Dimension 4 ist Eingabe (K3), nicht Ausgabe (z71, DIMENSION.UPLIFT.FIREWALL.01)")

# =================================================================== Teil 4
say("\nTeil 4: Bilanz")
urteil = ("T' aus Tomita: teilweise. Die Standardform von (End V, rho) liefert für jeden treuen Zustand kanonisch (bis auf positive "
          "Skala) den ganzen Inhalt von T': Lorentzsignatur entsteht genau aus einer Realstruktur, die M mit M' vertauscht, also "
          "aus J (die eigenen Realstrukturen der Naht geben (2,2), (4,0)); der Spinormodul des J-reellen Raums (Herm(2), det) ist "
          "V + V-bar^*, der Partnerplatz trägt (g^+)^-1 und ist das J-Bild V-bar, dualisiert über eps-bar; die Twistorform hat "
          "die Phase aus J und das Vorzeichen aus dem Kegel; L_0 = P(V) mit N = O(1) (x) V-bar^*, J-reelle Punkte Herm(2). Nicht "
          "abgeleitet ist die eine Identifikation I_pos: Die Lagen der Naht sind äquivariant die (dualen) GNS-Vektoren ihres "
          "Dublettzustands. T' <=> I_pos; die Standardform allein entscheidet das Normalenbündel nicht. Windungszahl: (a) "
          "Spinbrücke verlangt m ungerade, (b) Neutralität m = 4, und beides schließt sich für jedes m aus, weil R = "
          "Delta^(-i/2) auf jedem Quartersektor eine Ladungsdrehung ist. (c) lorentzisch und (d) sind keine Bedingung, (c) "
          "euklidisch schließt m = 4 aus. Ohne (b) bleibt der Gleichstand {1, 3}; das Wörterbuch z = S^m als Quotient nach einer "
          "Uhrenuntergruppe (m | 4) gibt m = 1, und Neutralität gilt dann für das Trägerdublett {1, zeta^4} mit rho_4 = Sym^4 "
          "rho_1. Kern: Im ganzen U(2) x U(2)-invarianten positiven Kegel erzwingen Hardy-Kern (Rang-1-Schnitt) und Stromkern "
          "(Skalarschnitt, Schwinger-Term m) zusammen K = det(I - Z W^+)^-1, also nu = 1; jede Auslesung allein reicht nicht. "
          "Daraus folgen SU(2,2)-Kovarianz, die Leiterdarstellung und über Cayley die masselose Wellengleichung auf R^(1,3) "
          "(Kinematik, unter T').")
say("   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z74_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "status": STATUS, "zahlen": ZAHLEN, "urteil": urteil}, indent=1,
    ensure_ascii=False, default=str) + "\n")
