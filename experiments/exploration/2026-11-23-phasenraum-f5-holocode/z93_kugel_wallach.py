"""Z93: Die Kugel D(V, eta) -- Wallach-Parameter aus dem Ladungssektor (Folgeauftrag F5 nach H6h44; Prämisse P1 aus z84).

Leitlinie 'Rekonstruktion statt Auswahl'. Firewalls: DIMENSION.UPLIFT.FIREWALL.01 (jede Aussage ist typisiert: Ladungsgitter der
beiden Quarter, Spinbrücke H0(N) = V^* (x) eta, 2x2-Matrixkugel; keine 3+1-Aussage), GRAV.NONCIRCULAR.01 (die Ladungsuhr ist der
Modulfluss und wird nicht als Raumzeitzeit gelesen). Bestand: z58 (Ladung x Oszillatoren, Sugawara, Jacobi-Tripelprodukt),
z63/z79 (Nahtdublett V_1 = {Omega, psi*_1/4 Omega}, Quarterpaar = eta-Dublett, R = e^{2 pi i L_0}), z74/z40 (det-Kern,
Wallach-Menge {0, 1} u (1, inf), Hardy- und Stromauslesung), z76 (K/pi = Q_a - Q_b mod 4), z84 (eindeutige Markenform auf H0(N),
D(V, eta), T' punktweise, Prämissen P1/P2), z86 (Ladungssektor bei beta = 4 pi, theta_3(i)). Ledger nur gelesen.

Teil 1  Ladungsdaten auf der Kugel: R = i^(Q_a - Q_b) auf dem Ladungsgitter, Dreieckszahlen, Theta bei tau = i. Der R-Verflechter
        Ladungsbasis -> H0(N) (OriginLab _transport) ist bis auf Diagonalphasen und Skala eindeutig und unitär für die
        Markenform; Deck und iota in Ladungssprache. Der Ladungstorus T^2 enthält R = T(1/4, -1/4), wirkt in K, die Ladungsuhr ist
        nichtzentral und keine Kontraktion, die Dilatation liegt nicht im Torus.
Teil 2  Der Kern des Ladungssektors: kohärente Zustände Phi_u (u in C^2), K_ch = Theta_1/4 Theta_3/4 (ganz, positiv, R-invariant,
        eine Nullstelle je Zelle). Satz: auf keiner ladungsäquivarianten Einbettung in D(V, eta) ist K_ch ein Kugelkern
        det(I - Z W^+)^-nu (Liouville, Fourier-Träger Z^2 gegen Gewichtskegel, Nullstellen), Ausgleichsrechnung; Paarzustände
        als Kontrollen (Bosonen nu = 1, Fermionen kompakt, Gitter flach); radiale Lesarten; Gram-Tests nu = 1/4, 1/2, 1, 3/2, 2.
Teil 3  Woher nu = 1 kommt: Bosonisierung bei beta = 4 pi (Fermi-Propagator = Twist x K_ch x Oszillatorfaktor), Strom =
        Oszillatoren + Var Q; Entscheidung über P1, Folgen für T' und die Innen/Außen-Architektur.
Teil 4  Bilanz und Urteil zu F5.
"""
import itertools
import json
import math
import sys
import time
from pathlib import Path

import mpmath as mp
import numpy as np
import scipy.linalg as sla
import sympy as sp

HERE = Path(__file__).parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "research" / "originlab"))
from reconstruction_search import Z as ZS, _hcoords, _hmul, _transport  # noqa: E402

mp.mp.dps = 40
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
T0 = time.time()
rng = np.random.default_rng(93)
PI = mp.pi
BETA = 4*PI
QNOM = mp.exp(-BETA)                                              # q = e^{-4 pi}: Euklidischer Torus tau_E = 2i (Naht 2 pi, beta)
A14, A34 = mp.mpf(1)/4, mp.mpf(3)/4
TH3I = PI**(mp.mpf(1)/4)/mp.gamma(mp.mpf(3)/4)                    # theta_3(i)
I2 = np.eye(2, dtype=complex)


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


def npm(M):
    return np.array(sp.Matrix(M).evalf(30).tolist(), dtype=complex)


def mpm(A):
    return mp.matrix([[mp.mpc(complex(A[i, j])) for j in range(A.shape[1])] for i in range(A.shape[0])])


# =================================================================== Teil 1
say("Teil 1: Ladungsdaten auf der Kugel D(V, eta)")

# --- 1.1 Ladungsgitter der beiden Quarter
E_s = lambda a, Q: sp.Rational(Q*Q, 2) + (a - sp.Rational(1, 2))*Q
QA_S, QB_S = sp.Rational(1, 4), sp.Rational(3, 4)
ok_R = True
for qa, qb in itertools.product(range(-10, 11), repeat=2):
    E4 = 4*(E_s(QA_S, qa) + E_s(QB_S, qb))
    ok_R &= E4.is_integer and (E4 - (qa - qb)) % 4 == 0             # e^{2 pi i E} = i^{4E} = i^{Q_a - Q_b}
dreieck = sorted(int(4*E_s(QA_S, Q)) for Q in range(-30, 31))
ok_tri = dreieck == [k*(k + 1)//2 for k in range(61)] and all(E_s(QB_S, Q) == E_s(QA_S, -Q) for Q in range(-30, 31))


def E_m(a, Q):
    return mp.mpf(Q)**2/2 + (a - mp.mpf(1)/2)*Q


def Theta(a, w, QR=40):
    """Ladungssektor eines Quarters mit Fugazität: sum_Q e^{-4 pi E_a(Q)} e^{2 pi i Q w} (Gitter |Q| <= QR)."""
    return mp.fsum(mp.exp(-BETA*E_m(a, Q) + 2j*PI*Q*w) for Q in range(-QR, QR + 1))


def Theta_j(a, w):
    """dieselbe Funktion geschlossen: jtheta_3(pi w -+ i pi/2 | Nome e^{-2 pi}) (tau = 2i, Charakteristik -+1/4)."""
    return mp.jtheta(3, PI*w + (-1 if a == A14 else 1)*0.5j*PI, mp.exp(-2*PI))


TH0 = Theta(A14, 0)
geschl = abs(TH0 - mp.exp(PI/8)*2**(-mp.mpf(5)/8)*TH3I)
ws_t = [mp.mpc(0.3, 0.2), mp.mpc(-0.7, 0.9), mp.mpc(0.1, -0.4), mp.mpc(0.45, 1.7)]
jdef = max(abs(Theta(a, w) - Theta_j(a, w))/abs(Theta_j(a, w)) for a in (A14, A34) for w in ws_t)
ZAHLEN["ladungsgitter"] = {"4E_1/4(Q), Q = -3..3": [int(4*E_s(QA_S, Q)) for Q in range(-3, 4)],
                           "Theta_1/4(0)": mp.nstr(TH0, 20), "Defekt e^(pi/8) 2^(-5/8) theta_3(i)": r(geschl, 3),
                           "Defekt jtheta": r(jdef, 3)}
check(f"Ladungsgitter: R = e^(2 pi i L_0) = i^(Q_a - Q_b) für |Q_a|, |Q_b| <= 10 exakt (4E ganz und = Q_a - Q_b mod 4, z76); "
      f"4 E_1/4(Q) = Q(2Q - 1) durchläuft die Dreieckszahlen 0, 1, 3, 6, 10, ... je genau einmal (|Q| <= 30, Gauß), E_3/4(Q) = "
      f"E_1/4(-Q): die Ladungsniveaus wachsen quadratisch, eines je Stufe und Quarter; Theta_1/4(0) = Theta_3/4(0) = "
      f"e^(pi/8) 2^(-5/8) theta_3(i) (Defekt {float(geschl):.0e}), Fugazitätsform = jtheta_3 bei tau = 2i mit Charakteristik "
      f"-+1/4 (Defekt {float(jdef):.0e})", ok_R and ok_tri and geschl < mp.mpf(10)**-38 and jdef < mp.mpf(10)**-36)

# --- 1.2 der native Verflechter Ladungsbasis -> H0(N)
BASIS_N = [(ZS**k, sp.Integer(0)) for k in range(2)] + [(sp.Integer(0), ZS**k/(ZS**2 - 1)) for k in range(2)]
G_KAND = {"1": (sp.Integer(1), sp.Integer(0)), "z": (ZS, sp.Integer(0)), "1/z": (1/ZS, sp.Integer(0)),
          "F": (sp.Integer(0), 1/(ZS**2 - 1)), "zF": (sp.Integer(0), ZS/(ZS**2 - 1))}


def holo_transport(phi, psi, g):
    """Linearer Transport s -> g (s o sigma^-1) mit sigma^-1(z, w) = (phi(z), psi(z) w) (wie z84)."""
    cols = []
    for a, b in BASIS_N:
        img = (sp.cancel(a.subs(ZS, phi)), sp.cancel(b.subs(ZS, phi)*psi))
        c = _hcoords(_hmul(g, img), 1)
        if c is None:
            return None
        cols.append(c)
    return sp.Matrix(cols).T


R_S, DECK_S = _transport(1, -1), _transport(1, 0)
IOTA_S = next(M for M in (holo_transport(1/ZS, -sp.I/ZS**2, g) for g in G_KAND.values()) if M is not None)


def herm_basis(n):
    B = []
    for i in range(n):
        M = np.zeros((n, n), complex)
        M[i, i] = 1
        B.append(M)
    for i, j in itertools.combinations(range(n), 2):
        M = np.zeros((n, n), complex)
        M[i, j] = M[j, i] = 1
        B.append(M)
        M = np.zeros((n, n), complex)
        M[i, j], M[j, i] = -1j, 1j
        B.append(M)
    return B


HB = herm_basis(4)
zeilen = [np.array([np.concatenate([(A.conj().T @ B @ A - B).real.ravel(), (A.conj().T @ B @ A - B).imag.ravel()])
                    for B in HB]).T for A in (npm(R_S), npm(DECK_S), npm(IOTA_S))]
_, s_, vh_ = np.linalg.svd(np.vstack(zeilen))
inv_formen = vh_[np.sum(s_ > 1e-9):]
h_inv = sum(c*B for c, B in zip(inv_formen[0], HB))
h_inv = h_inv/h_inv[0, 0]

# Ladungsseite: V = {|Q_a = 0>, |Q_a = 1>} (Nahtdublett V_1, Fock-Ecke), eta = {|x = 1/4>, |x = -1/4>} (Quarterpaar, Einteilchen)
L_V, L_E = [sp.Integer(0), sp.Rational(1, 4)], [sp.Rational(1, 4), -sp.Rational(1, 4)]
Q_V, Q_E = [(0, 0), (1, 0)], [(1, 0), (0, 1)]                      # Ladungen (Q_a, Q_b) der vier Objekte
IDX = [(i, j) for j in range(2) for i in range(2)]                  # Index k = i + 2j (i: V^*, j: eta), wie z84
R_CH = sp.diag(*[sp.nsimplify(sp.exp(2*sp.pi*sp.I*(L_E[j] - L_V[i]))) for i, j in IDX])
P_ETA = sp.Matrix([[1, 1], [-sp.I, sp.I]])/sp.sqrt(2)               # x = 1/4 -> (1 - iF)/sqrt2, x = -1/4 -> (1 + iF)/sqrt2
U_S = sp.kronecker_product(P_ETA, sp.eye(2))                        # 1 <-> <Q_a = 0|, z <-> <Q_a = 1|
ok_U = ((U_S.H*U_S - sp.eye(4)).applyfunc(sp.simplify).is_zero_matrix
        and (U_S.inv()*R_S*U_S - R_CH).applyfunc(sp.simplify).is_zero_matrix)
U_np = npm(U_S)
h_ch = U_np.conj().T @ h_inv @ U_np
kommutant = 16 - np.linalg.matrix_rank(np.kron(np.eye(4), npm(R_CH)) - np.kron(npm(R_CH).T, np.eye(4)), tol=1e-9)
D_CH, I_CH = npm(U_S.inv()*DECK_S*U_S), npm(U_S.inv()*IOTA_S*U_S)
monomial = lambda M: all(np.sum(np.abs(M[k]) > 1e-12) == 1 for k in range(4))
ZAHLEN["verflechter"] = {"R_Ladung": str(R_CH.tolist()), "Deck_Ladung": np.round(D_CH, 6).tolist(),
                         "iota_Ladung": np.round(I_CH, 6).tolist(), "Kommutantendimension": int(kommutant),
                         "invariante Formen": len(inv_formen)}
check(f"nativer Verflechter: die Vierteldrehung R von H0(N) (OriginLab _transport, Basis 1, z, F, zF) ist unter U = "
      f"(1 -> <Q_a = 0|, z -> <Q_a = 1|) (x) (x = 1/4 -> (1 - iF)/sqrt2, x = -1/4 -> (1 + iF)/sqrt2) exakt die Ladungs-R "
      f"e^(2 pi i (L_eta - L_V)) auf Hom(V, eta) mit L_V = (0, 1/4), L_eta = (1/4, -1/4), Eigenwerte {{i, 1, -i, -1}}; die "
      f"Markengruppe lässt {len(inv_formen)} hermitesche Form zu, U^+ h U = I (Ladungs-Orthonormalbasen = Markenform); der "
      f"Kommutant von R hat Dimension {kommutant}: U ist bis auf Diagonalphasen und die Skala (P2) eindeutig",
      ok_U and len(inv_formen) == 1 and np.allclose(h_ch, np.eye(4), atol=1e-12) and kommutant == 4)

# --- 1.3 Ladungstorus und Ladungsuhr auf D(V, eta)
GEW = np.array([[Q_E[j][0] - Q_V[i][0], Q_E[j][1] - Q_V[i][1]] for i, j in IDX])     # Gewichte der Einträge |e_j><v_i|


def T_ch(u):
    """Ladungstorus e^{2 pi i (u_a Q_a + u_b Q_b)} auf Hom(V, eta) (Vektorform, Index k = i + 2j)."""
    return np.diag(np.exp(2j*np.pi*(GEW @ np.asarray(u))))


def T_mat(u, Zm):
    """dieselbe Wirkung in Matrixform Z (Zeilen eta, Spalten V): Z -> A(u) Z B(u)^-1."""
    A = np.diag(np.exp(2j*np.pi*(np.array(Q_E) @ np.asarray(u))))
    B = np.diag(np.exp(2j*np.pi*(np.array(Q_V) @ np.asarray(u))))
    return A @ Zm @ np.linalg.inv(B)


ok_RT = np.allclose(T_ch([0.25, -0.25]), npm(R_CH), atol=1e-14)
uhr_lin = all(E_s(QA_S, qa) + E_s(QB_S, qb) - sp.Rational(qa - qb, 4) == sp.Rational(qa*(qa - 1) + qb*(qb + 1), 2)
              for qa, qb in itertools.product(range(-10, 11), repeat=2))
uhr_ecken = [L_V[0], L_V[1], L_E[0], L_E[1]]                         # V: Fock-Energien, eta: Einteilchenenergien x
ok_ecken = all(L == sp.Rational(q[0] - q[1], 4) for L, q in zip(uhr_ecken, Q_V + Q_E))
loch_R = [sp.I**(qa - qb) for qa, qb in ((1, 0), (0, -1))]          # Lochlesart der Mode x = -1/4: Ladung (0, -1)
ok_see = loch_R == [sp.I, sp.I] and [sp.I**(qa - qb) for qa, qb in Q_E] == [sp.I, -sp.I]
uhr_ecken = [float(x) for x in uhr_ecken]
gew_uhr = GEW @ np.array([0.25, -0.25])
rang_zentrum = np.linalg.matrix_rank(np.column_stack([GEW, np.ones(4)]))
faktor = lambda y: np.exp(-2*np.pi*(GEW @ np.asarray(y)))           # |Faktor| der Einträge unter T(i y)
uhr_exp = {f"tau = {'+' if s > 0 else '-'}beta/2": np.round(faktor(np.array([0.25, -0.25])*s), 4).tolist()
           for s in (1, -1)}                                        # imaginäre Uhr e^{-tau L}: T(i y) mit y = tau/(2 pi) (1/4, -1/4)
streng = 0
for _ in range(400):
    y = rng.normal(size=2)
    if np.all(faktor(y) < 1 - 1e-12):
        streng += 1
Zt = 0.99*np.array([[0, 0], [0, 1]], complex)                       # Eintrag |x = -1/4><Q_a = 1| (Gewicht (-1, 1))
raus = np.linalg.norm(T_mat(1j*np.array([0.25, -0.25])/(2*np.pi), Zt), 2)
ZAHLEN["ladungstorus"] = {"Gewichte (Q_a, Q_b) der Einträge": GEW.tolist(), "Uhrgewichte L_eta - L_V": gew_uhr.tolist(),
                          "L_0 auf V (+) eta": uhr_ecken, "Rang [Torusalgebra | Zentrum]": int(rang_zentrum),
                          "Betragsfaktoren der imaginären Uhr": uhr_exp, "strikte Kontraktionen in 400 Richtungen": streng,
                          "||T(i tau) Z||, ||Z|| = 0.99, tau = 1": r(raus, 5)}
check(f"Ladungstorus auf D(V, eta): die Einträge |e_j><v_i| tragen die Ladungen {GEW.tolist()}, R = T(1/4, -1/4) exakt; auf "
      f"der Fock-Ecke V und den Einteilchenmoden eta ist die Uhr linear, L_0 = (Q_a - Q_b)/4 = {uhr_ecken}; auf dem Fock-Gitter "
      f"gilt L_0 = (Q_a - Q_b)/4 + [Q_a(Q_a - 1) + Q_b(Q_b + 1)]/2 (exakt, |Q| <= 10; der quadratische Teil ist ganz, für R "
      f"unsichtbar und verschwindet genau auf den Fock-Ecken Q_a in {{0, 1}}, Q_b in {{-1, 0}}). Das eta-Dublett mit R = (i, -i) "
      f"verlangt die Einteilchenlesart der Mode x = -1/4 (L_0 = -1/4, Dirac-See); die Lochlesart (Ladung (0, -1), L_0 = +1/4) "
      f"gäbe R = (i, i). Die Uhr wirkt in K mit Gewichten {gew_uhr.tolist()} (nichtzentral); die "
      f"Torusalgebra enthält das Zentrum von k (Dilatation) nicht (Rang {rang_zentrum}); die imaginäre Uhr ist für beide "
      f"Vorzeichen keine Kontraktion (Faktoren {uhr_exp}; ||T Z|| = {raus:.3f} > 1 bei ||Z|| = 0,99), und keine von 400 "
      f"Richtungen des komplexen Torus kontrahiert strikt (Gewicht 0 des Eintrags |x = 1/4><Q_a = 1|)",
      ok_RT and uhr_lin and ok_ecken and ok_see and rang_zentrum == 3 and streng == 0 and raus > 1
      and all(any(f > 1 for f in v) and any(f < 1 for f in v) for v in uhr_exp.values()))
iota_torus = []
for _ in range(4):
    u = rng.uniform(-0.5, 0.5, size=2)
    X = I_CH @ T_ch(u) @ np.linalg.inv(I_CH)
    iota_torus.append(bool(np.allclose(X, np.diag(np.diag(X)), atol=1e-10)))
iRi = (IOTA_S*R_S*IOTA_S.inv() + R_S**3*DECK_S).applyfunc(sp.simplify).is_zero_matrix     # iota R iota^-1 = -R^3 Deck
ii_deck = (IOTA_S*IOTA_S - DECK_S).applyfunc(sp.simplify).is_zero_matrix                  # iota^2 = Deck
ZAHLEN["markengruppe_ladung"] = {"iota R iota^-1 = -R^3 Deck": iRi, "iota^2 = Deck": ii_deck,
                                 "iota T(u) iota^-1 diagonal (vier zufällige u)": iota_torus}
befund(f"Die Markengruppe in Ladungssprache: R ist ein Element des Ladungstorus, Deck vertauscht die beiden Quartermoden x = +-1/4 "
       f"auf eta (monomial: {monomial(D_CH)}), iota ist Teilchen-Loch auf V (Omega <-> psi*_1/4 Omega) mal einer Hadamard-"
       f"Mischung der Quartermoden (monomial: {monomial(I_CH)}); exakt iota R iota^-1 = -R^3 Deck ({iRi}) und iota^2 = Deck "
       f"({ii_deck}). Die Markengruppe normalisiert also <R, Deck>, aber weder <R> allein noch den Ladungstorus (iota T(u) "
       f"iota^-1 diagonal für zufällige u: {iota_torus}): Der Torus ist ladungsnativ, nicht markennativ. Der Ladungssektor wirkt "
       "auf D(V, eta) nur durch einen kompakten Torus in K; die radiale Richtung (Dilatation, deren Spektrum nu + N den "
       "Wallach-Parameter trägt) erreicht er nicht",
       monomial(D_CH) and not monomial(I_CH) and not any(iota_torus) and iRi and ii_deck)
ZAHLEN["laufzeit_teil1_s"] = r(time.time() - T0, 3)

# =================================================================== Teil 2
say("\nTeil 2: Der Kern des Ladungssektors und sein Wallach-Parameter")


def herm_info(G, rel=mp.mpf(10)**-25):
    n = G.rows
    skala = max(abs(G[i, j]) for i in range(n) for j in range(n))
    defekt = max(abs(G[i, j] - mp.conj(G[j, i])) for i in range(n) for j in range(n))/skala
    E = mp.eighe((G + G.H)/2, eigvals_only=True)
    E = sorted(mp.re(E[k]) for k in range(n))
    emax = max(abs(e) for e in E)
    pos, neg = sum(1 for e in E if e > rel*emax), sum(1 for e in E if e < -rel*emax)
    return {"herm_defekt": defekt, "min": E[0], "max": E[-1], "inertia": (pos, neg, n - pos - neg)}


def kurz(info):
    return {"herm_defekt": r(info["herm_defekt"], 3), "min_eig": r(info["min"], 4), "max_eig": r(info["max"], 4),
            "inertia": info["inertia"]}


def K_ch(u, v):
    """<Phi_u|Phi_v> mit Phi_u = sum_Q e^{-2 pi E(Q)} e^{2 pi i Q.u} |Q_a, Q_b> (Amplitude e^{-beta E/2}, beta = 4 pi)."""
    return Theta_j(A14, v[0] - mp.conj(u[0]))*Theta_j(A34, v[1] - mp.conj(u[1]))


# --- 2.1 kohärente Zustände des Ladungsgitters
QRV = 12
GITTER = list(itertools.product(range(-QRV, QRV + 1), repeat=2))


def phi_vec(u):
    return [mp.exp(-2*PI*(E_m(A14, qa) + E_m(A34, qb)) + 2j*PI*(qa*u[0] + qb*u[1])) for qa, qb in GITTER]


US = [(mp.mpc(*rng.uniform(-0.5, 0.5, 2)*[1, 0.8]), mp.mpc(*rng.uniform(-0.5, 0.5, 2)*[1, 0.8])) for _ in range(24)]
ov_def = mp.mpf(0)
for u, v in [(US[0], US[1]), (US[2], US[3]), (US[4], US[4])]:
    pu, pv = phi_vec(u), phi_vec(v)
    ov = mp.fsum(mp.conj(a)*b for a, b in zip(pu, pv))
    ov_def = max(ov_def, abs(ov - K_ch(u, v))/abs(ov))
G_ch = mp.matrix(24, 24)
for i in range(24):
    for j in range(24):
        G_ch[i, j] = K_ch(US[i], US[j])
gi_ch = herm_info(G_ch)
RSH = (A14, -A14)
r_inv = max(abs(K_ch((u[0] + RSH[0], u[1] + RSH[1]), (v[0] + RSH[0], v[1] + RSH[1])) - K_ch(u, v))/abs(K_ch(u, v))
            for u, v in [(US[5], US[6]), (US[7], US[8]), (US[9], US[9])])
wq = mp.mpc(0.3, 0.2)
qp_def = max(abs(Theta_j(A14, wq + 1) - Theta_j(A14, wq)), abs(Theta_j(A14, wq + 2j) - mp.exp(PI - 2j*PI*wq)*Theta_j(A14, wq)))
mp.mp.dps = 25
dTh = lambda w: mp.fsum(2j*PI*Q*mp.exp(-BETA*E_m(A14, Q) + 2j*PI*Q*w) for Q in range(-30, 31))
ecken = [mp.mpc(0, 0), mp.mpc(1, 0), mp.mpc(1, 2), mp.mpc(0, 2), mp.mpc(0, 0)]
windung = mp.fsum(mp.quad(lambda t, a=a, b=b: dTh(a + t*(b - a))/Theta_j(A14, a + t*(b - a))*(b - a), [0, 0.5, 1])
                  for a, b in zip(ecken[:-1], ecken[1:]))/(2j*PI)
mp.mp.dps = 40
W0 = {A14: mp.mpc(0.5, 1.5), A34: mp.mpc(0.5, 0.5)}
null_def = max(abs(Theta_j(a, W0[a] + m + 2j*k))/abs(Theta_j(a, W0[a] + m + 2j*k + 0.25))
               for a in (A14, A34) for m in (0, 1) for k in (-1, 0, 1, 2))
ZAHLEN["koh_ladung"] = {"Überlapp gegen geschlossene Form": r(ov_def, 3), "Gram 24": kurz(gi_ch), "R-Invarianz": r(r_inv, 3),
                        "Quasiperiodizität": r(qp_def, 3), "Nullstellenzahl je Zelle [0,1] x [0,2i]": r(mp.re(windung), 12),
                        "Nullstellen": {"1/4": "1/2 + 3i/2 + Z + 2iZ", "3/4": "1/2 + i/2 + Z + 2iZ"},
                        "|Theta(Null)|/|Theta(Null + 1/4)| max": r(null_def, 3)}
check(f"kohärente Zustände Phi_u = sum_Q e^(-2 pi E(Q)) e^(2 pi i Q.u)|Q_a, Q_b> (u in C^2, normierbar für jedes u): K_ch(u, v) = "
      f"Theta_1/4(v_a - conj u_a) Theta_3/4(v_b - conj u_b) (Überlapp der Gittervektoren gegen geschlossene Form "
      f"{float(ov_def):.0e}); Gram von 24 Zuständen positiv definit {gi_ch['inertia']}; R-invariant, K_ch(u + r, v + r) = "
      f"K_ch(u, v) mit r = (1/4, -1/4) (Defekt {float(r_inv):.0e}); Theta(w + 1) = Theta(w), Theta(w + 2i) = e^(pi - 2 pi i w) "
      f"Theta(w) (Defekt {float(qp_def):.0e}); Argumentprinzip auf der Zelle [0, 1] x [0, 2i]: {mp.nstr(mp.re(windung), 12)} "
      f"Nullstelle, Nullgitter 1/2 + 3i/2 + Z + 2iZ (1/4) und 1/2 + i/2 + Z + 2iZ (3/4) (|Theta| relativ {float(null_def):.0e}): "
      "der Kern ist ganz, ohne Rand, mit einer Nullstelle je Periodenzelle",
      ov_def < mp.mpf(10)**-30 and gi_ch["inertia"] == (24, 0, 0) and r_inv < mp.mpf(10)**-35 and qp_def < mp.mpf(10)**-35
      and abs(windung - 1) < mp.mpf(10)**-15 and null_def < mp.mpf(10)**-35)


# --- 2.2 Satz: kein Kugelkern auf ladungsäquivarianten Einbettungen
def zuf_kugel(smin=0.15, smax=0.85):
    U = np.linalg.qr(rng.normal(size=(2, 2)) + 1j*rng.normal(size=(2, 2)))[0]
    V = np.linalg.qr(rng.normal(size=(2, 2)) + 1j*rng.normal(size=(2, 2)))[0]
    return U @ np.diag(rng.uniform(smin, smax, size=2)) @ V


# (a) Liouville: Normierbarkeitsbereiche
Z0s = []
while len(Z0s) < 6:
    Zc = zuf_kugel(0.2, 0.8)
    if np.abs(Zc).min() > 0.05:
        Z0s.append(Zc)
aus_kugel = min(np.linalg.norm(T_mat(np.array([0, -0.6j]), Zc), 2) for Zc in Z0s)
norm_ch = [mp.re(K_ch((mp.mpc(0, y), mp.mpc(0, -y)), (mp.mpc(0, y), mp.mpc(0, -y)))) for y in (1, 3, 5)]
norm_bos = [float(1/np.linalg.det(np.eye(2) - (s*I2).conj().T @ (s*I2)).real) for s in (0.9, 0.99, 0.999)]
check(f"Normierbarkeitsbereiche (Liouville): Phi_u ist für jedes u in C^2 normierbar (||Phi_(iy, -iy)||^2 = {[mp.nstr(x, 4) for x in norm_ch]} "
      f"für y = 1, 3, 5), der Parameterraum der Ladungszustände ist also ganz C^2; jede holomorphe Abbildung C^2 -> D(V, eta) "
      f"(beschränkt) ist konstant, und die äquivariante Bahnabbildung u -> T(u) Z_0 verlässt die Kugel (min ||T(0, -0,6i) Z_0|| = "
      f"{aus_kugel:.2f} > 1 für sechs zufällige Z_0). Bosonische Paarzustände e^(a^+ Z b^+)|0> sind dagegen genau auf der Kugel "
      f"normierbar (det(I - Z^+ Z)^-1 = {[r(x, 4) for x in norm_bos]} für Z = 0,9 / 0,99 / 0,999 I): D(V, eta) ist ihr "
      "natürlicher Parameterraum, nicht der des Ladungsgitters",
      aus_kugel > 1 and norm_bos[-1] > 1e5 and norm_bos[0] < norm_bos[1] < norm_bos[2]
      and all(mp.isfinite(x) and x > 0 for x in norm_ch))

# (b) Fourier-Träger
NF = 64
gw = np.array([(j/NF, k/NF) for j in range(NF) for k in range(NF)])


def fourier_traeger(werte):
    C = np.fft.fft2(werte.reshape(NF, NF))/NF**2                     # Koeffizient von e^{2 pi i (n_a w_a + n_b w_b)}
    aussen = inner = 0.0
    for na in range(-NF//2, NF//2):
        for nb in range(-NF//2, NF//2):
            c = abs(C[na % NF, nb % NF])
            if nb >= 0 and na + nb >= 0:
                inner = max(inner, c)
            else:
                aussen = max(aussen, c)
    return aussen, inner


kegel = {}
for nu in (0.5, 1, 1.5, 2):
    Zc = zuf_kugel(0.2, 0.55)
    werte = np.array([np.linalg.det(I2 - T_mat(w, Zc) @ Zc.conj().T)**(-nu) for w in gw])
    kegel[str(nu)] = fourier_traeger(werte)
ch_aussen = {"(0,-1)": mp.exp(-BETA*(E_m(A14, 0) + E_m(A34, -1))), "(-1,0)": mp.exp(-BETA*(E_m(A14, -1) + E_m(A34, 0))),
             "(0,-2)": mp.exp(-BETA*E_m(A34, -2)), "(0,-3)": mp.exp(-BETA*E_m(A34, -3))}
werte_ch = np.array([complex(Theta_j(A14, w[0])*Theta_j(A34, w[1])) for w in gw])
ch_fft = fourier_traeger(werte_ch)
ZAHLEN["fourier_traeger"] = {"Kugelkerne (max außerhalb, max innerhalb des Kegels)": {k: [r(v[0], 3), r(v[1], 3)]
                                                                                     for k, v in kegel.items()},
                             "Ladungskern außerhalb (exakt)": {k: mp.nstr(v, 6) for k, v in ch_aussen.items()},
                             "Ladungskern FFT (außen, innen)": [r(ch_fft[0], 4), r(ch_fft[1], 4)]}
check(f"Fourier-Träger auf der Ladungsbahn: jeder Kugelkern det(I - Z(w) Z_0^+)^-nu entlang u -> T(u) Z_0 hat nur Frequenzen im "
      f"Gewichtskegel S = {{n_b >= 0, n_a + n_b >= 0}}, erzeugt von den Eintragsladungen (1,0), (0,0), (0,1), (-1,1) (FFT 64 x 64: "
      f"außerhalb max {max(v[0] for v in kegel.values()):.0e}, innerhalb bis {min(v[1] for v in kegel.values()):.2f}, nu = 1/2, "
      f"1, 3/2, 2): ein Niedrigstgewichtsmodul. Der Ladungskern hat alle Ladungen (Q_a, Q_b) in Z^2 mit Gaußgewicht e^(-4 pi E), "
      f"z. B. (0,-1): {mp.nstr(ch_aussen['(0,-1)'], 4)} = e^-pi, (-1,0): {mp.nstr(ch_aussen['(-1,0)'], 4)}, (0,-2), (0,-3) != 0 "
      f"(FFT außen {ch_fft[0]:.4f}): sein Träger liegt in keiner Verschiebung k_0 + S, also ist K_ch für kein nu und keinen "
      "Multiplikator ein zurückgeholter Kugelkern (beidseitige Ladung gegen positive Gewichte)",
      max(v[0] for v in kegel.values()) < 1e-12 and min(v[1] for v in kegel.values()) > 1e-3
      and abs(ch_fft[0] - float(ch_aussen["(0,-1)"])) < 1e-12 and all(v > 0 for v in ch_aussen.values()))

KAND = {"1": [(0, 0), (1, 1), (-1, -1)], "i": [(1, 0), (0, -1)], "-i": [(-1, 0), (0, 1)]}   # R-Charakter i^(Q_a - Q_b), |Q| <= 1
RICHT = np.array([(np.cos(t), np.sin(t)) for t in np.linspace(0, 2*np.pi, 1441)])
zensus = {"eigentlicher Kegel": 0, "ganze Ebene": 0}
for v1, v2, e1, e2 in itertools.product(KAND["1"], KAND["i"], KAND["i"], KAND["-i"]):
    Wg = np.array([[e[0] - v[0], e[1] - v[1]] for v in (v1, v2) for e in (e1, e2)])
    eigentlich = bool(np.any(np.all(Wg @ RICHT.T >= -1e-12, axis=0)))
    zensus["eigentlicher Kegel" if eigentlich else "ganze Ebene"] += 1
    if (v1, v2, e1, e2) == ((0, 0), (1, 0), (1, 0), (0, 1)):
        zensus["nativ (Fock-Ecke, Quarterpaar)"] = "eigentlicher Kegel" if eigentlich else "ganze Ebene"
ZAHLEN["fourier_traeger"]["Zensus R-verträglicher Ladungszuordnungen (|Q| <= 1)"] = zensus

# (c) Ausgleichsrechnung log K_ch = c - nu log D und Kontrolle
gfit = np.array([(j/16, k/16) for j in range(16) for k in range(16)])
LF = np.array([complex(mp.log(Theta_j(A14, w[0])*Theta_j(A34, w[1]))) for w in gfit])
LF = LF - LF[0]
fits, kontrollen = [], []
for Zc in Z0s:
    LD = np.array([np.sum(np.log(np.linalg.eigvals(I2 - T_mat(w, Zc) @ Zc.conj().T))) for w in gfit])
    LD = LD - LD[0]
    nu_f = -np.vdot(LD, LF)/np.vdot(LD, LD)
    fits.append((nu_f, np.linalg.norm(LF + nu_f*LD)/np.linalg.norm(LF)))
    for nu0 in (1.0, 1.5):
        LS = -nu0*LD
        nu_c = -np.vdot(LD, LS)/np.vdot(LD, LD)
        kontrollen.append((abs(nu_c - nu0), np.linalg.norm(LS + nu_c*LD)/np.linalg.norm(LS)))
ZAHLEN["ausgleich"] = {"nu_fit": [r(f[0].real, 3) for f in fits], "Imaginärteil nu_fit": [r(f[0].imag, 2) for f in fits],
                       "relatives Residuum": [r(f[1], 3) for f in fits],
                       "Kontrolle max |nu_fit - nu| / Residuum": [r(max(k[0] for k in kontrollen), 2),
                                                                  r(max(k[1] for k in kontrollen), 2)]}
check(f"Ausgleichsrechnung auf der reellen Torusbahn (16 x 16 Punkte, sechs zufällige Z_0): log K_ch = c - nu log D lässt "
      f"relative Residuen {[r(f[1], 2) for f in fits]} mit nu_fit = {[r(f[0].real, 2) for f in fits]} (je Basispunkt anders); "
      f"Kontrolle mit synthetischem D^-nu (nu = 1, 3/2): nu exakt zurück (Abweichung {max(k[0] for k in kontrollen):.0e}, "
      f"Residuum {max(k[1] for k in kontrollen):.0e})",
      min(f[1] for f in fits) > 0.2 and max(k[0] for k in kontrollen) < 1e-10 and max(k[1] for k in kontrollen) < 1e-10)
befund("Satz (Nullstellen): Gilt K_ch(u, v) = conj f(u) f(v) det(I - Z(v) Z(u)^+)^-nu auf der reellen Bahn einer "
       "ladungsäquivarianten Einbettung Z(u) = T(u) Z_0, dann im Komplexen; der Quotient hängt nur von w = v - conj u ab und "
       "faktorisiert, also ist er e^(lambda.w) (Cauchy). K_ch verschwindet auf den Hyperebenen w_a in 1/2 + 3i/2 + Z + 2iZ; "
       "e^(lambda.w) D^-nu ist für nu >= 0 dort nullstellenfrei, für nu < 0 wäre das Laurentpolynom D in y_a = e^(2 pi i w_a) "
       "an unendlich vielen Stellen y_0 e^(-4 pi k) null, also D = 0, Widerspruch zu D(0) = det(I - Z_0 Z_0^+) > 0. Das gilt für "
       "jeden Basispunkt Z_0, jedes reelle nu und beliebige Gewichte der Einbettung; R-Fixpunkte Z_0 = kappa |x = 1/4><Q_a = 1| "
       f"geben eine einpunktige Bahn und einen konstanten Kern. Das Fourier-Argument hängt dagegen an der Zuordnung: von 24 "
       f"R-verträglichen Ladungszuordnungen mit |Q| <= 1 haben {zensus['eigentlicher Kegel']} einen eigentlichen Gewichtskegel "
       f"(darunter die native: {zensus['nativ (Fock-Ecke, Quarterpaar)']}), {zensus['ganze Ebene']} spannen die Ebene auf; "
       "für diese trägt nur der Nullstellensatz", True)

# (d) Paarzustände als Kontrollen: Bosonen, Fermionen, Gitter
def S_koef(W, r1, r2, c1, c2):
    s = 0j
    for k in range(min(r1, c1) + 1):
        n = (k, r1 - k, c1 - k, r2 - c1 + k)
        if min(n) < 0:
            continue
        s += W[0, 0]**n[0]*W[0, 1]**n[1]*W[1, 0]**n[2]*W[1, 1]**n[3]/math.prod(math.factorial(x) for x in n)
    return s


def boson_ueberlapp(Zm, Wm, Nmax=34):
    tot = 0j
    for N in range(Nmax + 1):
        for r1, c1 in itertools.product(range(N + 1), repeat=2):
            f = math.factorial(r1)*math.factorial(N - r1)*math.factorial(c1)*math.factorial(N - c1)
            tot += f*np.conj(S_koef(Zm, r1, N - r1, c1, N - c1))*S_koef(Wm, r1, N - r1, c1, N - c1)
    return tot


def fock_ops(n):
    a = np.array([[0, 1], [0, 0]], complex)
    zt = np.diag([1., -1.]).astype(complex)
    ops = []
    for k in range(n):
        mats = [zt]*k + [a] + [I2]*(n - k - 1)
        op = mats[0]
        for Mx in mats[1:]:
            op = np.kron(op, Mx)
        ops.append(op)
    return ops


FO = fock_ops(4)
VAK = np.zeros(16, complex)
VAK[0] = 1
fzust = lambda W: sla.expm(sum(W[i, j]*FO[i].conj().T @ FO[2 + j].conj().T for i in range(2) for j in range(2))) @ VAK
bos_def = fer_def = 0.0
for _ in range(3):
    Zm, Wm = zuf_kugel(0.05, 0.4), zuf_kugel(0.05, 0.4)
    bos_def = max(bos_def, abs(boson_ueberlapp(Zm, Wm)*np.linalg.det(I2 - Zm.conj().T @ Wm) - 1))
    fer_def = max(fer_def, abs(np.vdot(fzust(Zm), fzust(Wm)) - np.linalg.det(I2 + Zm.conj().T @ Wm)))
pois = max(abs(Theta_j(A14, w) - mp.exp(PI/8)/mp.sqrt(2)*mp.fsum(mp.exp(0.5j*PI*(w - m) - PI*(w - m)**2/2)
                                                                    for m in range(-60, 61))) for w in ws_t)
Zf, Wf = zuf_kugel(), zuf_kugel()
flach = [abs(np.linalg.det(I2 - Zf @ Wf.conj().T/n)**(-n) - np.exp(np.trace(Zf @ Wf.conj().T))) for n in (10, 100, 1000)]
ZAHLEN["paarzustaende"] = {"Bosonen |<Z|W> det(I - Z^+W) - 1|": r(bos_def, 3), "Fermionen |<Z|W> - det(I + Z^+W)|": r(fer_def, 3),
                           "Poisson Theta = periodisierte Gaußglocke": r(pois, 3),
                           "flacher Grenzfall |det(I - ZW^+/n)^-n - e^tr(ZW^+)|, n = 10, 100, 1000": [r(x, 3) for x in flach]}
check(f"Paarzustände über Hom(V, eta) als Kontrollen: bosonisch e^(sum Z_ij a_i^+ b_j^+)|0> (Fockraum bis Grad 34) gibt exakt "
      f"det(I - Z^+ W)^-1, nu = 1 (Defekt {bos_def:.0e}); fermionisch (16-dim. Fockraum) det(I + Z^+ W), der kompakte Partner "
      f"nu = -1 (Defekt {fer_def:.0e}); das Ladungsgitter gibt Theta = periodisierte Gaußglocke, Theta_1/4(w) = e^(pi/8)/sqrt2 "
      f"sum_m e^(i pi (w - m)/2 - pi (w - m)^2/2) (Poisson, Defekt {float(pois):.0e}), also einen flachen (Heisenberg-, Fock-) "
      f"Kern; in der Kugelfamilie ist das der Kontraktionslimes det(I - Z W^+/n)^-n -> e^(tr Z W^+) (Abweichung "
      f"{[r(x, 2) for x in flach]} für n = 10, 100, 1000), kein Punkt der Wallach-Menge",
      bos_def < 1e-12 and fer_def < 1e-12 and pois < mp.mpf(10)**-35 and flach[0] > flach[1] > flach[2])

# (e) radiale Lesarten (Ladungsuhr als radiale Zeit, Einheit 1/4)
NMAX = 60000                                                        # t^NMAX < e^-60 für 1 - t = 1e-3
TRI = [k*(k + 1)//2 for k in range(400) if k*(k + 1)//2 <= NMAX]
p1 = np.zeros(NMAX + 1)
p1[TRI] = 1
p2 = np.zeros(NMAX + 1)
for a_ in TRI:
    for b_ in TRI:
        if a_ + b_ <= NMAX:
            p2[a_ + b_] += 1


def c_eff(p, t):
    n = np.arange(len(p))
    w_ = p*np.exp(n*np.log(t))
    return (1 - t)*np.sum(n*w_)/t/np.sum(w_)


ceff = {f"1 - t = 1e-{k}": (r(c_eff(p1, 1 - 10.0**-k), 4), r(c_eff(p2, 1 - 10.0**-k), 4)) for k in (1, 2, 3)}
koef_dev = {"ein Quarter gegen (1/2)_n/n!": r(max(abs(p1[n] - math.gamma(n + 0.5)/(math.gamma(0.5)*math.factorial(n)))
                                                   for n in range(41)), 3),
            "beide Quarter gegen 1": int(max(abs(int(p2[n]) - 1) for n in range(41)))}
zaehl = {L: (int(np.sum(p2[:L + 1])), r(np.sum(p2[:L + 1])/L, 4), r(sum((n + 1)**2 for n in range(L + 1))/L**3, 4))
         for L in (100, 1000, 3000)}
ZAHLEN["radiale_lesarten"] = {"c_eff (ein Quarter, beide)": ceff, "Koeffizienten": koef_dev,
                              "r(n), n = 0..12": [int(x) for x in p2[:13]],
                              "Niveauzahl N(L), N(L)/L, Kugel nu = 1 sum (n+1)^2 / L^3": zaehl,
                              "nu_eff": {"ein Quarter": {"Rang-1-Schnitt": 0.5, "Skalarschnitt": 0.25},
                                         "beide Quarter": {"Rang-1-Schnitt": 1.0, "Skalarschnitt": 0.5}}}
check(f"radiale Lesart (Ladungsuhr als Dilatation, Stufen n = 4 L_0): ein Quarter sum_k t^(T_k) hat Koeffizienten 0/1 statt "
      f"(c)_n/n! und den Randexponenten c_eff -> 1/2 (Poisson), beide Quarter haben r(n) = #{{T_j + T_k = n}} = "
      f"{[int(x) for x in p2[:13]]} (schwankend, auch 0) und c_eff -> 1: {ceff}. Der Mittelwert passt zum Hardy-Exponenten "
      f"einer komplexen Dimension, die Koeffizienten zu keinem det^-nu-Schnitt (Rang-1: alle 1; Abweichung {koef_dev['beide Quarter gegen 1']}). "
      f"Niveauzahl N(L)/L -> pi/2 ({[zaehl[L][1] for L in zaehl]}) gegen die nu = 1-Kugel sum (n+1)^2 ~ L^3/3 "
      f"({[zaehl[L][2] for L in zaehl]}): linear gegen kubisch, keine uhrtreue Isometrie auf einen nu-Raum der Kugel; "
      "nu_eff hängt am gewählten Schnitt (1/4, 1/2 oder 1)",
      abs(ceff["1 - t = 1e-3"][0] - 0.5) < 0.02 and abs(ceff["1 - t = 1e-3"][1] - 1) < 0.03
      and koef_dev["beide Quarter gegen 1"] >= 1 and koef_dev["ein Quarter gegen (1/2)_n/n!"] > 0.3
      and abs(zaehl[3000][1] - np.pi/2) < 0.05 and abs(zaehl[3000][2] - 1/3) < 0.01)

# (f) Gram-Tests der Kugelkerne auf R-Bahnen in D(V, eta)
A_R, B_R = np.diag([1j, -1j]), np.diag([1, 1j])
R_mat = lambda Zm: A_R @ Zm @ np.linalg.inv(B_R)
PTS = []
for _ in range(6):
    Zc = zuf_kugel(0.15, 0.9)
    for k in range(4):
        PTS.append(Zc)
        Zc = R_mat(Zc)
PTS_mp = [mpm(P) for P in PTS]


def det2(A):
    return A[0, 0]*A[1, 1] - A[0, 1]*A[1, 0]


def kball(Zm, Wm, nu):
    """det(I - Z W^+)^-nu über die Eigenwerte von Z W^+ (|lambda| < 1, Hauptzweig stetig; wie z84)."""
    A = Zm*Wm.H
    tr, de = A[0, 0] + A[1, 1], det2(A)
    d = mp.sqrt(tr**2 - 4*de)
    return mp.exp(-nu*(mp.log(1 - (tr + d)/2) + mp.log(1 - (tr - d)/2)))


mp.mp.dps = 30
I2_MP = mp.eye(2)
grams = {}
for nu in ("0.25", "0.5", "1", "1.5", "2"):
    G = mp.matrix(24, 24)
    for i in range(24):
        for j in range(24):
            G[i, j] = kball(PTS_mp[i], PTS_mp[j], mp.mpf(nu))
    grams[nu] = herm_info(G, rel=mp.mpf(10)**-20)
rinv_k = max(abs(kball(mpm(R_mat(PTS[i])), mpm(R_mat(PTS[j])), mp.mpf("0.5")) - kball(PTS_mp[i], PTS_mp[j], mp.mpf("0.5")))
             for i, j in [(0, 5), (3, 17), (9, 22)])
mp.mp.dps = 40
ZAHLEN["gram_kugel"] = {k: kurz(v) for k, v in grams.items()} | {"R-Invarianz det-Kern": r(rinv_k, 3)}
say(f"   {ZAHLEN['gram_kugel']}")
check(f"Gram-Tests auf D(V, eta) (sechs zufällige Punkte mit ihren R-Bahnen, 24 Punkte, 30 Stellen): det(I - Z W^+)^-nu ist "
      f"für nu = 1, 3/2, 2 positiv definit ({grams['1']['inertia']}, {grams['1.5']['inertia']}, {grams['2']['inertia']}), für "
      f"nu = 1/4 und 1/2 indefinit ({grams['0.25']['inertia']}, {grams['0.5']['inertia']}); jeder nu-Kern ist exakt R-invariant "
      f"(Defekt {float(rinv_k):.0e}), R wählt also kein nu. Die radialen Ladungslesarten nu_eff = 1/4 und 1/2 liegen im "
      "indefiniten Bereich; nur die Lesart 'beide Quarter, Rang-1-Schnitt' trifft nu = 1, und dort stimmen die Koeffizienten nicht",
      all(grams[k]["inertia"] == (24, 0, 0) for k in ("1", "1.5", "2"))
      and all(grams[k]["inertia"][1] > 0 for k in ("0.25", "0.5")) and rinv_k < mp.mpf(10)**-25)
ZAHLEN["laufzeit_bis_teil2_s"] = r(time.time() - T0, 3)

# =================================================================== Teil 3
say("\nTeil 3: Woher nu = 1 kommt, Entscheidung über P1")


def G_fermi(a, u, R=200):
    """thermischer Fermi-Propagator des Quarters a: sum_{x in Z + a} e^{i x u}/(1 + e^{-4 pi x}), u = theta + i tau, 0 < tau < beta."""
    return mp.fsum(mp.exp(1j*(n + a)*u)/(1 + mp.exp(-BETA*(n + a))) for n in range(-R, R + 1))


def osz_faktor(zeta):
    """(q;q)^2/theta(zeta; q), theta(x; q) = (x;q)(q/x;q): Stromoszillatoren, einfacher Pol bei zeta = 1."""
    return mp.qp(QNOM, QNOM)**2/(mp.qp(zeta, QNOM)*mp.qp(QNOM/zeta, QNOM))


bos_split = mp.mpf(0)
for a in (A14, A34):
    for u in (mp.mpc(0.7, 1.3), mp.mpc(2.1, 5.0), mp.mpc(-1.0, 9.0), mp.mpc(3.0, 12.0)):
        rhs = mp.exp(1j*a*u)*(K_ch((0, 0), (u/(2*PI), 0)) if a == A14 else K_ch((0, 0), (0, u/(2*PI))))/K_ch((0, 0), (0, 0))
        rhs *= osz_faktor(mp.exp(1j*u))
        bos_split = max(bos_split, abs(G_fermi(a, u) - rhs)/abs(rhs))
eps = mp.mpf(10)**-15
residuum = abs((1 - mp.exp(1j*eps))*osz_faktor(mp.exp(1j*eps)) - 1)
varQ = {}
for a in (A14, A34):
    ps = {Q: mp.exp(-BETA*E_m(a, Q)) for Q in range(-40, 41)}
    Zs = mp.fsum(ps.values())
    mq = mp.fsum(Q*p for Q, p in ps.items())/Zs
    var_lad = mp.fsum(Q*Q*p for Q, p in ps.items())/Zs - mq**2
    var_fermi = mp.fsum(1/((1 + mp.exp(BETA*(n + a)))*(1 + mp.exp(-BETA*(n + a)))) for n in range(-80, 81))
    varQ[str(a)] = (var_lad, abs(var_lad - var_fermi))
ZAHLEN["bosonisierung"] = {"Defekt G = Twist x K_ch x Oszillator": r(bos_split, 3), "Residuum Oszillatorpol - 1": r(residuum, 3),
                           "Var Q (Ladung) = sum c(1 - c) (Fermi)": {k: [mp.nstr(v[0], 12), r(v[1], 3)] for k, v in varQ.items()}}
check(f"Bosonisierung bei beta = 4 pi exakt (Kronecker, 40 Stellen): der Fermi-Propagator jedes Quarters ist G_a(u) = e^(i a u) x "
      f"K_ch(0, u/2 pi)/K_ch(0, 0) x (q;q)^2/theta(e^(iu); q) (Defekt {float(bos_split):.0e}, a = 1/4, 3/4, vier Punkte): Twist x "
      f"kohärenter Ladungskern aus Teil 2 x Oszillatorfaktor. Der Pol sitzt allein im Oszillatorfaktor, (1 - zeta) x Faktor -> 1 "
      f"(Defekt {float(residuum):.0e}), das ist z74s Hardy-Auslesung (1 - t)^-1, also nu = 1 auf dem Rang-1-Schnitt; der "
      f"Ladungsfaktor ist ganz (Teil 2) und trägt nur die Charakteristik (den R-Twist). Strom: die Nullmode gibt die Konstante "
      f"Var Q = {mp.nstr(varQ[str(A14)][0], 8)} (Ladungstheta = Fermi-Varianz sum c(1 - c), Defekt "
      f"{float(max(v[1] for v in varQ.values())):.0e}), der Doppelpol kommt aus den Oszillatoren",
      bos_split < mp.mpf(10)**-35 and residuum < mp.mpf(10)**-12 and max(v[1] for v in varQ.values()) < mp.mpf(10)**-35)
befund("P1 entschieden: Der Ladungssektor ist der einzige Sektor, auf dem R nichttrivial wirkt, und er liefert auf D(V, eta) die "
       "hermitesche Form (Teil 1) und R (ein Torus-Element), aber keinen Kugelkern: Sein natürlicher Kern ist flach (Theta, "
       "Heisenberg), ganz und beidseitig geladen, auf keiner ladungsäquivarianten Einbettung ein det^-nu (Teil 2). Beide "
       "Auslesungen, die nu = 1 erzwingen (Hardy-Pol und Strom-Doppelpol, z74), sind Oszillatorfaktoren der Bosonisierung, und auf "
       "den Oszillatoren ist R trivial (z79). Der nu = 1-Kern auf D(V, eta) ist zulässig (R-invariant, positiv), aber nicht "
       "rekonstruiert: P1 ist eine Wahl, keine Folgerung. Der Ladungssektor wählt auch kein anderes nu", True)
befund("Folge für T': Die radiale Inversion rho(Z) = (Z^+)^-1 auf D(V, eta) erzeugt in jedem Randpunkt die Lorentzform (z84), "
       "ist aber nur dann OS-Spiegelung eines Zustands, wenn ein positiver Kern auf D(V, eta) physikalisch realisiert ist. Ohne "
       "Ladungskern bleibt T' punktweise mit freiem R kinematisch (Geometrie der Markenform); physikalisch ist T' nur auf der "
       "Oszillatorkugel, wo R trivial ist. Architektur: Innen = Ladungsgitter mit Marken, R und Spin-Z4, flacher Theta-Typ; "
       "Außen = Oszillatoren mit nu = 1-Leiter, SU(2,2) und Lorentzform. D(V, eta) gehört nach innen als Symmetriedatum "
       "(Form + endliche Markengruppe), nicht als Raumzeitkugel; das bestätigt S1 und die Lesart 'Marken innen' (z85)", True)
ZAHLEN["laufzeit_bis_teil3_s"] = r(time.time() - T0, 3)

# =================================================================== Teil 4
say("\nTeil 4: Bilanz und Urteil zu F5")
STATUS = {
    "Ladungsdaten auf H0(N): R = T(1/4, -1/4) im Ladungstorus; R-Verflechter bis auf Diagonalphasen und Skala eindeutig, "
    "Ladungs-Orthonormalbasen = Markenform": "E",
    "Markengruppe in Ladungssprache: Deck = Quartertausch auf eta, iota = Teilchen-Loch auf V (x) Hadamard auf eta, iota R "
    "iota^-1 = -R^3 Deck; normalisiert wird <R, Deck>, nicht der Ladungstorus": "E",
    "Ladungsuhr auf D(V, eta): linear (Q_a - Q_b)/4 auf Fock-Ecke und Einteilchenmoden, nichtzentral in K (Gewicht -1/2 aus der "
    "Dirac-See-Mode x = -1/4), keine Kontraktion; Dilatation nicht im Torus": "E",
    "kohärenter Ladungskern K_ch = Theta_1/4 Theta_3/4 (tau = 2i, Charakteristik -+1/4, = theta_3(i)-Klasse): ganz, positiv, "
    "R-invariant, eine Nullstelle je Zelle": "E",
    "K_ch ist auf keiner ladungsäquivarianten Einbettung in D(V, eta) ein Kugelkern det^-nu, für kein reelles nu "
    f"(Liouville, Fourier-Träger, Nullstellen; Ausgleich Residuum >= {min(f[1] for f in fits):.2f})": "E (negativ)",
    "Kerntyp nach Statistik: Bosonen-Paare nu = 1 (normierbar genau auf der Kugel), Fermionen-Paare kompakt (nu = -1), "
    "Ladungsgitter flach (Heisenberg, Kontraktionslimes)": "E",
    "radiale Lesarten des Ladungssektors: nu_eff = 1/4, 1/2 oder 1 je nach Schnitt, Koeffizienten arithmetisch, Niveauzahl "
    "linear statt kubisch": "E (keine native Auslesung)",
    "det^-nu auf D(V, eta): R-invariant für jedes nu, positiv für nu = 1, 3/2, 2, indefinit für 1/4, 1/2": "E",
    "Bosonisierung bei beta = 4 pi: Fermi-Propagator = Twist x K_ch x Oszillatorfaktor; Hardy-Pol und Strom-Doppelpol in den "
    "Oszillatoren, Ladung gibt Twist und Var Q": "E",
    "P1 (nu = 1-Kern gehört zu D(V, eta)) aus dem Ladungssektor": "E (negativ: nicht rekonstruiert, kein anderes nu gewählt)",
    "P1 als Prämisse": "O (zulässige Wahl; nach 'Rekonstruktion statt Auswahl' nicht nativ)",
    "T' punktweise auf D(V, eta) mit freiem R": "L (kinematisch; ohne realisierten positiven Kern keine OS-Spiegelung eines Zustands)",
    "Innen/Außen: Marken innen (flacher Ladungstyp), Raumzeit außen (nu = 1-Oszillatoren)": "C (gestützt; Gravitation, 3+1 unberührt)",
}
ZAHLEN["status"] = STATUS
for k_, v_ in STATUS.items():
    say(f"   [{v_}] {k_}")
z = ZAHLEN
urteil = (
    "F5 endet negativ für P1, mit positivem Strukturbefund. (1) Die Ladungsdaten landen nativ auf H0(N) = Hom(V, eta): V ist "
    "die Fock-Ecke {Q_a = 0, 1} des Quarters 1/4, eta das Quarterpaar x = +-1/4, der Verflechter zur OriginLab-Basis ist bis auf "
    "Diagonalphasen und die Skala P2 eindeutig, und die Ladungs-Orthonormalbasen geben genau die eindeutige Markenform von z84. "
    "R ist das Element T(1/4, -1/4) des Ladungstorus; die Ladungsuhr ist auf der Fock-Ecke und den Einteilchenmoden die lineare "
    "Ladungsdrehung (Q_a - Q_b)/4 (der quadratische Sugawara-Teil verschwindet dort), wirkt in K nichtzentral mit Gewichten "
    "(1/4, 0, -1/4, -1/2) (das negative Gewicht stammt von der Dirac-See-Mode x = -1/4, die das eta-Dublett erzwingt) und ist "
    "keine Kontraktion; die Dilatation liegt nicht im Torus, und die Markengruppe normalisiert <R, Deck>, nicht den Torus. "
    "(2) Der natürliche Kern des Ladungssektors, der Überlapp der "
    "kohärenten Gitterzustände, ist K_ch = Theta_1/4 Theta_3/4 (Theta bei tau = 2i mit Charakteristik -+1/4, Theta(0) = "
    "e^(pi/8) 2^(-5/8) theta_3(i)): positiv, R-invariant, ganz, mit einer Nullstelle je Zelle und Ladungsträger Z^2. Er ist auf "
    "keiner ladungsäquivarianten Einbettung in D(V, eta) ein Kugelkern det(I - Z W^+)^-nu, für kein reelles nu: Liouville (die "
    "Zustände sind auf ganz C^2 normierbar, D ist beschränkt), Fourier-Träger (Kugelkerne bleiben im Gewichtskegel, Theta ist "
    f"beidseitig geladen, (0,-1) mit e^-pi) und Nullstellen; die Ausgleichsrechnung lässt Residuen {z['ausgleich']['relatives Residuum']}. "
    "Kontrollen: bosonische Paarzustände geben exakt nu = 1 und sind genau auf der Kugel normierbar, fermionische den kompakten "
    "Partner, das Gitter den flachen Heisenberg-Kern (Kontraktionslimes nu -> unendlich). Radiale Lesarten geben nu_eff = 1/4, "
    "1/2 oder 1 je nach Schnitt, mit arithmetischen Koeffizienten und linearer statt kubischer Niveauzahl. Gram-Tests: det^-nu "
    "ist auf R-Bahnen positiv für nu = 1, 3/2, 2, indefinit für 1/4 und 1/2, und für jedes nu R-invariant. (3) Woher nu = 1 "
    "kommt, zeigt die exakte Bosonisierung bei beta = 4 pi: G_a = Twist x K_ch x Oszillatorfaktor, der Hardy-Pol (und der "
    "Strom-Doppelpol) sitzt in den Oszillatoren, wo R trivial ist; die Ladung trägt nur den Twist und die Konstante Var Q. P1 ist "
    "damit aus dem Ladungssektor nicht rekonstruierbar und auch kein anderes nu: der nu = 1-Kern auf D(V, eta) bleibt eine "
    "zulässige Wahl. T' punktweise mit freiem R ist kinematisch; physikalisch trägt nur die Oszillatorkugel (R trivial) die "
    "Lorentzform. Die Architektur wird schärfer: Marken innen im flachen Ladungstyp, Raumzeit außen in der nu = 1-Leiter, "
    "D(V, eta) als inneres Symmetriedatum. Nichts befördert; Ledger unverändert.")
say("\nUrteil\n   " + urteil)
ZAHLEN["laufzeit_gesamt_s"] = r(time.time() - T0, 3)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)}); Laufzeit "
    f"{time.time() - T0:.0f} s")
(HERE / "z93_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
