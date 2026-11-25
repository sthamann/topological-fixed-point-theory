"""Z84: Euklidisch zuerst -- native Reflexion und Osterwalder-Schrader (Strategie S3 aus H6h41; Fortsetzung von
WOIT.OS.TWISTOR.01, Ledger nur gelesen).

Leitlinie 'Rekonstruktion statt Auswahl'. Firewalls: GRAV.NONCIRCULAR.01 (die Zeit wird aus der Reflexionspositivität
rekonstruiert und nicht am Modulfluss abgelesen; wo die Rekonstruktion den Modulfluss zurückgibt, zählt das als Befund),
DIMENSION.UPLIFT.FIREWALL.01 (jede Aussage ist typisiert: 1D-Naht, Spinbrücke H0(N), 2x2-Matrixkugel, euklidischer
Schnitt). Bestand: z61 (Antipoden, quaternionische Form), z70 (eta, periodischer Modulfluss), z72 (parabolische Zeit im
Vakuumsektor), z74 (T' <=> I_pos, Kern nu = 1, Cayley), z79 (R auf H0(N) = V^* (x) eta euklidisch, tau), z83 (radiale
su(2,2)); Ledger v519/v522/v524/v565 (gelesen, nicht verändert).

Teil 1  Kandidaten theta. Exakter Transport auf H0(N) = L(R_1 + R_-1 + oo_+ + oo_-) (OriginLab _transport, hier um
        antiholomorphe Abbildungen erweitert) für alle zwölf Hebungen der markentreuen antiholomorphen Involutionen der
        Naht (Geraden z -> mu conj z, Schwarz-Spiegelung am Markenkreis = Teilchen-Loch-Spiegel, Antipoden).
        Satz A: jede geometrische Abbildung wirkt als Produkt auf V^* (x) eta (sie vertauscht nie die Faktoren).
        Satz B: keine faktortauschende (anti)lineare Abbildung normalisiert R (projektive Ordnungen 4 und 2); auf
        V (x) V-quer (Prämisse T') vertauscht J die Faktoren und zentralisiert R. Tomita-J, thermische Spiegelung und
        Teilchen-Loch-Spiegel bilden eine Klasse (KMS-Identität, Standardform).
Teil 2  Reflexionspositivität mit Gram-Matrizen (40 Stellen): thermische Schwingerfunktionen des Originals (Einteilchen,
        Strom, Zweiteilchen), nu-Kern unter radialer Inversion, Halbraumspiegelung, Weyl-Spinor mit Theta_n (externe
        Analyse), thermischer 4D-Kern bei beta = 4 pi. Kontrollen: nu = 1/2, Antipoden, Drehung statt Spiegelung,
        falsche Spinstruktur.
Teil 3  OS-Rekonstruktion: thermisch -> Liouville-Operator (zweiseitig, periodisch, also der Modulfluss); radial ->
        Gradoperator (diskret, Vielfachheiten m^2); Halbraum -> P_0 (stetig auf [0, inf), Gamma-Gesetz von z72).
Teil 4  T': Lorentz = Euklid o Inversion; die bi-unitäre Kugel der Spinbrücke; R wirkt frei auf dem Shilov-Rand und auf
        S^3; Theta_n ist die U(2) x U(2)-Bahn von z79s tau. T' gilt punktweise, nie markeninvariant.
Teil 5  Bilanz und Urteil zu S3, Bezug zu WOIT.OS.TWISTOR.01.
"""
import itertools
import json
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
rng = np.random.default_rng(84)
BETA = 4*mp.pi
I2 = np.eye(2, dtype=complex)
SIG = [np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]]), np.diag([1., -1.]).astype(complex)]


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


def realign_rang(M, tausch=False, tol=1e-9):
    """M auf V^* (x) eta, Index k = i + 2 j (i: V^*, j: eta). Produkt a (x) b <=> Rang 1 von R[(i,k),(j,l)];
    Faktortausch (M[i + 2j, k + 2l] = a[i,l] b[j,k], V^* <-> eta) <=> Rang 1 von R[(i,l),(j,k)]."""
    Rm = np.zeros((4, 4), complex)
    for i, j, k, l in itertools.product(range(2), repeat=4):
        if tausch:
            Rm[2*i + l, 2*j + k] = M[i + 2*j, k + 2*l]
        else:
            Rm[2*i + k, 2*j + l] = M[i + 2*j, k + 2*l]
    s = np.linalg.svd(Rm, compute_uv=False)
    return int(np.sum(s > tol*s[0]))


def faktoren(M):
    Rm = np.zeros((4, 4), complex)
    for i, j, k, l in itertools.product(range(2), repeat=4):
        Rm[2*i + k, 2*j + l] = M[i + 2*j, k + 2*l]
    u, s, vh = np.linalg.svd(Rm)
    return (u[:, 0]*np.sqrt(s[0])).reshape(2, 2), (vh[0]*np.sqrt(s[0])).reshape(2, 2)


def proj_ordnung2(g):
    """projektive Ordnung eines 2x2-Elements endlicher Ordnung aus dem Eigenwertverhältnis."""
    e = np.linalg.eigvals(g)
    q = e[0]/e[1]
    for k in range(1, 33):
        if abs(q**k - 1) < 1e-9:
            return k
    return None


def proportional(X, Y, tol=1e-9):
    k = np.unravel_index(np.argmax(np.abs(Y)), Y.shape)
    c = X[k]/Y[k]
    return bool(np.abs(X - c*Y).max() < tol*max(1.0, np.abs(X).max())), c


# =================================================================== Teil 1
say("Teil 1: Kandidaten theta und ihre Wirkung auf H0(N) = V^* (x) eta")


def conj_k(f):
    """f -> conj(f(conj z)): Koeffizienten konjugieren (Koeffizienten in Q(i))."""
    return sp.sympify(f).subs(sp.I, -sp.I)


BASIS_N = [(ZS**k, sp.Integer(0)) for k in range(2)] + [(sp.Integer(0), ZS**k/(ZS**2 - 1)) for k in range(2)]
G_KAND = {"1": (sp.Integer(1), sp.Integer(0)), "z": (ZS, sp.Integer(0)), "1/z": (1/ZS, sp.Integer(0)),
          "F": (sp.Integer(0), 1/(ZS**2 - 1)), "zF": (sp.Integer(0), ZS/(ZS**2 - 1))}


def anti_transport(phi, psi, g):
    """Antilinearer Transport s -> g conj(s o sigma) für die antiholomorphe Involution sigma(z, w) = (phi(conj z),
    psi(conj z) conj w) von E: w^2 = 1 - z^4. Rückgabe M mit A(v) = M conj(v), oder None, falls L(D) verlassen wird."""
    phib, psib = conj_k(phi), conj_k(psi)
    cols = []
    for a, b in BASIS_N:
        img = (sp.cancel(conj_k(a).subs(ZS, phib)), sp.cancel(conj_k(b).subs(ZS, phib)*psib))
        c = _hcoords(_hmul(g, img), 1)
        if c is None:
            return None
        cols.append(c)
    return sp.Matrix(cols).T


def holo_transport(phi, psi, g):
    """Linearer Transport s -> g (s o sigma^-1) mit sigma^-1(z, w) = (phi(z), psi(z) w)."""
    cols = []
    for a, b in BASIS_N:
        img = (sp.cancel(a.subs(ZS, phi)), sp.cancel(b.subs(ZS, phi)*psi))
        c = _hcoords(_hmul(g, img), 1)
        if c is None:
            return None
        cols.append(c)
    return sp.Matrix(cols).T


# --- 1.1 R und Deck (OriginLab), Kontrolle der eigenen Transportroutine
R_S = _transport(1, -1)
DECK_S = _transport(1, 0)
R_eig = holo_transport(-sp.I*ZS, sp.Integer(-1), G_KAND["F"])
R_np, DECK_np = npm(R_S), npm(DECK_S)
gV, gEta = faktoren(R_np)
ordV, ordEta = proj_ordnung2(gV), proj_ordnung2(gEta)
ZAHLEN["R_auf_H0N"] = {"matrix": str(R_S.tolist()), "projektive Ordnung V^*": ordV, "projektive Ordnung eta": ordEta}
check(f"R = OriginLab _transport(1, -1) auf H0(N) (Basis 1, z, F, zF) stimmt mit der eigenen Transportroutine überein, ist "
      f"ein Produkt g_V (x) g_eta (Realignment-Rang {realign_rang(R_np)}) mit projektiven Ordnungen {ordV} auf V^* = "
      f"span{{1, z}} und {ordEta} auf eta = span{{1, F}} (z79: Viertel auf V^*, Halb auf eta)",
      (R_eig - R_S).applyfunc(sp.simplify).is_zero_matrix and realign_rang(R_np) == 1 and (ordV, ordEta) == (4, 2))

# --- 1.2 Zensus der markentreuen antiholomorphen Involutionen
SPIEGEL = {}
for mu, nm in ((sp.Integer(1), "Gerade durch +-1"), (sp.I, "Diagonale arg = pi/4"), (sp.Integer(-1), "Gerade durch +-i"),
               (-sp.I, "Diagonale arg = -pi/4")):
    for eps in (1, -1):
        SPIEGEL[f"D: z -> {mu} conj z ({nm}), Hebung {eps:+d}"] = (mu*ZS, sp.Integer(eps), "D")
for s_, nm in ((1, "Schwarz am Markenkreis z -> 1/conj z (Teilchen-Loch, thermischer Schatten)"),
               (-1, "Antipoden z -> -1/conj z")):
    for eps in (1, -1):
        SPIEGEL[f"A: {nm}, Hebung {eps:+d}"] = (sp.Integer(s_)/ZS, eps*sp.I/ZS**2, "A")


def reelle_signatur(M):
    """Signatur von det(c_ij) auf dem reellen Schnitt {v = M conj v} (falls M conj(M) = I)."""
    Mr, Mi = M.real, M.imag
    A8 = np.block([[Mr - np.eye(4), Mi], [Mi, -Mr - np.eye(4)]])
    _, s, vh = np.linalg.svd(A8)
    nb = vh[np.sum(s > 1e-9):]
    if nb.shape[0] == 0:
        return "kein reeller Schnitt", 0
    basis = [v[:4] + 1j*v[4:] for v in nb]
    q = lambda u, v: (u[0]*v[3] + v[0]*u[3] - u[2]*v[1] - v[2]*u[1])/2
    B = np.array([[q(u, v) for v in basis] for u in basis])
    k = np.unravel_index(np.argmax(np.abs(B)), B.shape)
    ph = B[k]/abs(B[k])
    Bp = B/ph
    if np.abs(Bp.imag).max() > 1e-8*np.abs(Bp).max():
        return "komplex", len(basis)
    ew = np.linalg.eigvalsh((Bp.real + Bp.real.T)/2)
    p, n_ = int(np.sum(ew > 1e-9)), int(np.sum(ew < -1e-9))
    return f"({max(p, n_)},{min(p, n_)})", len(basis)


zensus, ok_inv, ok_g, ok_prod, ok_deck = {}, True, True, True, True
R3 = np.linalg.matrix_power(R_np, 3)
for nm, (phi, psi, fam) in SPIEGEL.items():
    ok_inv &= sp.cancel(psi**2*(1 - ZS**4) - (1 - phi**4)) == 0
    ok_inv &= sp.cancel(phi.subs(ZS, conj_k(phi)) - ZS) == 0 and sp.cancel(psi.subs(ZS, conj_k(phi))*conj_k(psi) - 1) == 0
    gs = {k: anti_transport(phi, psi, g) for k, g in G_KAND.items()}
    passend = [k for k, M in gs.items() if M is not None]
    ok_g &= len(passend) == 1
    M = npm(gs[passend[0]])
    A2 = M @ M.conj()
    sq = A2[0, 0]
    ist_skalar = np.allclose(A2, sq*np.eye(4), atol=1e-12) and abs(sq.imag) < 1e-12
    Mn = M/np.sqrt(abs(sq))
    sign = int(np.sign(sq.real))
    rel = "?"
    X = M @ R_np.conj() @ np.linalg.inv(M)
    for lab, Y in (("R", R_np), ("R Deck", R_np @ DECK_np), ("R^3", R3), ("R^3 Deck", R3 @ DECK_np)):
        if proportional(X, Y)[0]:
            rel = lab
    ok_deck &= proportional(M @ DECK_np.conj() @ np.linalg.inv(M), DECK_np)[0]
    pr, tr = realign_rang(M), realign_rang(M, tausch=True)
    ok_prod &= pr == 1
    sigt = reelle_signatur(Mn)[0] if sign > 0 else "kein reeller Schnitt"
    zensus[nm] = {"g": passend[0], "A^2": sign if ist_skalar else "nicht skalar", "Produktrang": pr, "Tauschrang": tr,
                  "R": rel, "reeller Schnitt": sigt, "Familie": fam}
    say(f"   {nm}: g = {passend[0]}, A^2 = {sign:+d}, Produktrang {pr}, {rel}, Schnitt {sigt}")
ZAHLEN["zensus_H0N"] = zensus
check("alle zwölf Hebungen sind antiholomorphe Involutionen von E (w^2 = 1 - z^4 exakt erhalten, sigma^2 = id), und für "
      "jede ist der Multiplikator g mit div g = sigma(D) - D erzwungen: genau einer der Kandidaten {1, z, 1/z, F, zF} hält "
      "H0(N) fest (Geraden durch Marken: 1; Diagonalen: F; Schwarz und Antipoden: z)", ok_inv and ok_g)
fam_rel = {fam: sorted({v["R"] for v in zensus.values() if v["Familie"] == fam}) for fam in ("D", "A")}
check(f"Satz A (Zensus): jede Hebung wirkt als Produkt a (x) b auf V^* (x) eta (Realignment-Rang 1), vertauscht die Faktoren "
      f"also nie, und normalisiert das Deck; A R A^-1 = {fam_rel}: Familie D invertiert die Vierteldrehung, Familie A "
      "zentralisiert sie auf der Sphäre (auf E tauscht sie die Spinorhebung w -> -w gegen die Hebung w -> w), wie v519 "
      "auf dem Gitter (Familie D invertiert die Uhr, Familie A zentralisiert sie)",
      ok_prod and ok_deck and set(fam_rel["D"]) <= {"R^3", "R^3 Deck"} and set(fam_rel["A"]) <= {"R", "R Deck"})
schnitte = sorted({v["reeller Schnitt"] for v in zensus.values()})
quadrate = {nm: v["A^2"] for nm, v in zensus.items()}
check(f"reelle Formen auf H0(N): die zwölf Hebungen geben nur die Signatur (2,2) (reell (x) reell) oder keinen reellen "
      f"Schnitt (A^2 = -1, quaternionisch (x) reell); gefunden {schnitte}. Weder die euklidische (4,0) noch die lorentzsche "
      "(1,3) Form entsteht; die Antipoden sind auf V^* quaternionisch, auf eta reell (A^2 = -1 für beide Hebungen), die "
      "Schwarz-Spiegelung reell (x) reell (2,2)",
      set(schnitte) <= {"(2,2)", "kein reeller Schnitt"} and all(quadrate[k] == -1 for k in quadrate if "Antipoden" in k)
      and all(zensus[k]["reeller Schnitt"] == "(2,2)" for k in zensus if "Schwarz" in k))

# --- 1.3 Satz B: faktortauschende Abbildungen und R
ord_k = {k: (proj_ordnung2(np.linalg.matrix_power(gV, k)), proj_ordnung2(np.linalg.matrix_power(gEta, k))) for k in range(8)}
tausch_ord = (ordEta, ordV)
q_V, q_eta = [np.linalg.eigvals(g) for g in (gV, gEta)]
verh = lambda e: {complex(np.round(e[0]/e[1], 9)), complex(np.round(e[1]/e[0], 9))}
R_V1 = np.diag([1, 1j])                                          # Nahtdublett V_1 (z79 2a)
AdR = np.kron(np.linalg.inv(R_V1).T, R_V1)                       # X -> R X R^-1 auf End V = V (x) V-quer^*
M_J = np.zeros((4, 4), complex)
for i, j in itertools.product(range(2), repeat=2):
    M_J[i + 2*j, j + 2*i] = 1                                    # J X = X^+ (antilinear: v -> M_J conj v)
J_zentral = proportional(M_J @ AdR.conj() @ np.linalg.inv(M_J), AdR)[0]
ZAHLEN["satz_B"] = {"Ordnungen (V^*, eta) von R^k": {k: v for k, v in ord_k.items()}, "getauscht": tausch_ord,
                    "Eigenwertverhältnisse V^*": sorted(map(str, verh(q_V))), "eta": sorted(map(str, verh(q_eta)))}
check(f"Satz B: eine faktortauschende (anti)lineare Abbildung theta konjugiert R = g_V (x) g_eta in ein Produkt mit den "
      f"projektiven Ordnungen {tausch_ord}; die Potenzen R^k haben {sorted(set(ord_k.values()))}, also normalisiert keine "
      "solche theta die Markengruppe <R>; ebenso gibt es keinen projektiven R-Verflechter eta -> V^* (Eigenwertverhältnisse "
      f"{sorted(map(str, verh(q_eta)))} gegen {sorted(map(str, verh(q_V)))}). Kontrolle unter der Prämisse T': auf End V = "
      f"V (x) V-quer^* hat Ad R die Ordnungen ({proj_ordnung2(R_V1)}, {proj_ordnung2(np.linalg.inv(R_V1).T)}), J: X -> X^+ "
      f"vertauscht die Faktoren (Tauschrang {realign_rang(M_J, tausch=True)}) und zentralisiert Ad R",
      tausch_ord not in set(ord_k.values()) and not (verh(q_V) & verh(q_eta)) and J_zentral
      and realign_rang(M_J, tausch=True) == 1 and realign_rang(AdR) == 1)
befund("Folge aus Satz A und B: Eine OS-Reflexion, die V^* und eta vertauscht (wie es die euklidische Orientierungsumkehr für "
       "die beiden Chiralitäten verlangt) und dabei die Marken respektiert, gibt es als konstante Abbildung auf H0(N) nicht; "
       "jede geometrische Spiegelung der Naht erhält die Faktoren. Auf V (x) V-quer^* (T' gesetzt) existiert sie: J", True)

# --- 1.4 Tomita-J, thermische Spiegelung und Teilchen-Loch-Spiegel: eine Klasse
say("   thermische Klasse (KMS bei beta = 4 pi)")


def c_f(x, beta=BETA):
    return 1/(1 + mp.exp(beta*x))


def cb_f(x, beta=BETA):
    return 1/(1 + mp.exp(-beta*x))                                 # 1 - c(x) ohne Auslöschung


def g_prop(x, tau, beta_c=BETA):
    return mp.exp(-tau*x)*cb_f(x, beta_c)


xs_t = [mp.mpf(k)/4 for k in (1, -1, 3, -3, 5, -5, 9, -9)]
taus_t = [mp.mpf(t) for t in ("0.3", "1.7", "3.1", "5.9", "9.0", "12.0")]
kms_def = max(abs(g_prop(x, BETA - t) - g_prop(-x, t))/abs(g_prop(-x, t)) for x in xs_t for t in taus_t)
kms_ctrl = min(max(abs(g_prop(x, BETA - t, 3*mp.pi) - g_prop(-x, t, 3*mp.pi))/abs(g_prop(-x, t, 3*mp.pi))
                   for t in taus_t) for x in xs_t)
ZAHLEN["kms_identitaet"] = {"rel_defekt_beta_4pi": r(kms_def, 3), "kontrolle_beta_3pi_min": r(kms_ctrl, 3)}
check(f"Teilchen-Loch = thermische Spiegelung: g_x(beta - tau) = g_-x(tau) mit g_x(tau) = e^(-tau x)(1 - c(x)) für "
      f"x in +-{{1/4, 3/4, 5/4, 9/4}} und sechs tau (rel. Defekt {float(kms_def):.0e}, 40 Stellen), äquivalent zu KMS und "
      f"c(-x) = 1 - c(x); Kontrolle: mit Besetzungen bei beta' = 3 pi und Spiegelung bei beta/2 = 2 pi verfehlt (min. "
      f"rel. Abweichung {float(kms_ctrl):.2f})", kms_def < mp.mpf(10)**-32 and kms_ctrl > 1e-2)


def fock_ops(n):
    a = np.array([[0, 1], [0, 0]], complex)
    zt = np.diag([1., -1.]).astype(complex)
    ops = []
    for k in range(n):
        mats = [zt]*k + [a] + [np.eye(2, dtype=complex)]*(n - k - 1)
        op = mats[0]
        for Mx in mats[1:]:
            op = np.kron(op, Mx)
        ops.append(op)
    return ops


MODEN = [0.25, -0.75, 1.25, 0.75, -0.25]                          # Quartersektoren 1/4 (drei Moden) und 3/4 (zwei)
psi = fock_ops(len(MODEN))
hdiag = np.real(np.diag(sum(x*p.conj().T @ p for x, p in zip(MODEN, psi))))
bf = 4*np.pi
eH = lambda u: np.diag(np.exp(-u*hdiag))                          # e^{-u H}
rho_f = eH(bf)/np.trace(eH(bf))
mono = psi + [p.conj().T for p in psi] + [psi[a_] @ psi[b_].conj().T for a_, b_ in ((0, 1), (2, 3), (1, 4))]
zuf_op = lambda: sum(complex(*rng.normal(size=2))*m for m in mono)
Phi = lambda A, tau: eH(bf/2 - tau) @ A @ eH(tau)
heis = lambda A, u: np.diag(np.exp(u*hdiag)) @ A @ eH(u)          # A(u) = e^{uH} A e^{-uH}
Rf = np.diag(np.exp(2j*np.pi*hdiag))
ok_sf = True
for _ in range(4):
    A, B = zuf_op(), zuf_op()
    tau, sig = rng.uniform(0.2, bf/2 - 0.2, size=2)
    s_ = rng.uniform(0.05, bf/2 - 0.2 - max(tau, sig))
    os_form = np.trace(Phi(A, tau).conj().T @ Phi(B, sig))
    eukl = np.trace(eH(bf) @ heis(A.conj().T, bf - tau) @ heis(B, sig))
    ok_sf &= abs(os_form - eukl) < 1e-10*abs(eukl)
    D_ = lambda X, u: np.diag(np.exp(-u*bf*hdiag)) @ X @ np.diag(np.exp(u*bf*hdiag))   # Delta^u X = rho^u X rho^-u
    ok_sf &= np.allclose(Phi(A, tau + s_), D_(Phi(A, tau), -s_/bf), rtol=1e-10, atol=1e-12*np.abs(Phi(A, tau)).max())
    ok_sf &= np.allclose(Phi(A, tau).conj().T, Phi(A.conj().T, bf/2 - tau), rtol=1e-10, atol=1e-14)
    X = zuf_op()
    AdRf = lambda Y: Rf @ Y @ Rf.conj().T
    ok_sf &= np.allclose(AdRf(X).conj().T, AdRf(X.conj().T)) and np.allclose(AdRf(X).conj(), Rf.conj().T @ X.conj() @ Rf)
    ok_sf &= np.allclose(AdRf(X), D_(X, -0.5j))                     # R = Delta^(-i/2)
check("Standardform (fünf Moden beider Quartersektoren, Fockraum 32, Hilbert-Schmidt 1024): die OS-Form der Spiegelung "
      "tau -> beta - tau ist das Standard-Skalarprodukt der Vektoren Phi(A, tau) = e^-(beta/2 - tau)H A e^-tau H; die "
      "Verschiebung tau -> tau + s ist Delta^(-s/beta); J (X -> X^+) bildet Phi(A, tau) auf Phi(A^+, beta/2 - tau) ab, ist "
      "also die Spiegelung an beta/4; J zentralisiert R = Delta^(-i/2) = Ad e^(2 pi i L_0), die modenerhaltende "
      "Konjugation K (Familie D) invertiert R", ok_sf)
befund("Eine Klasse: theta_cut (tau -> beta - tau, Spiegel tau = 0 und beta/2) definiert das Skalarprodukt (die S- bzw. "
       "+-Struktur), J ist die Spiegelung an beta/4 (bildet die Naht tau = 0 auf die Partnernaht tau = beta/2 ab), "
       "theta_cut o J = Drehung um beta/2 = Delta^(1/2) (wie v524/v565: theta_cut ~ Metrik, theta_perp ~ Theta_phys). Auf "
       "Einteilchenpropagatoren ist theta_cut der Teilchen-Loch-Spiegel x -> -x (KMS-Identität); sein Schatten auf der "
       "Naht (u = e^(i theta - tau) -> 1/conj u) ist die Schwarz-Spiegelung am Markenkreis: Produkt auf H0(N), Signatur "
       "(2,2), zentralisiert R. Kein Mitglied der Klasse vertauscht V^* und eta", True)
ZAHLEN["laufzeit_teil1_s"] = r(time.time() - T0, 3)

# =================================================================== Teil 2
say("\nTeil 2: Reflexionspositivität mit Gram-Matrizen (40 Stellen)")


def herm_info(G, rel=mp.mpf(10)**-28):
    n = G.rows
    skala = max(abs(G[i, j]) for i in range(n) for j in range(n))
    defekt = max(abs(G[i, j] - mp.conj(G[j, i])) for i in range(n) for j in range(n))/skala
    H = (G + G.H)/2
    E = mp.eighe(H, eigvals_only=True)
    E = sorted(mp.re(E[k]) for k in range(n))
    emax = max(abs(e) for e in E)
    pos, neg = sum(1 for e in E if e > rel*emax), sum(1 for e in E if e < -rel*emax)
    return {"herm_defekt": defekt, "min": E[0], "max": E[-1], "inertia": (pos, neg, n - pos - neg)}


def kurz(info):
    return {"herm_defekt": r(info["herm_defekt"], 3), "min_eig": r(info["min"], 4), "max_eig": r(info["max"], 4),
            "inertia": info["inertia"]}


# --- 2.1 thermische Schwingerfunktionen des Originals (Kreis tau in [0, beta), Naht theta in [0, 2 pi))
QA = (mp.mpf(1)/4, mp.mpf(3)/4)
MODE = {a: [(mp.mpf(k) + a, c_f(mp.mpf(k) + a), cb_f(mp.mpf(k) + a)) for k in range(-170, 170)] for a in QA}


def S_f(dth, dta, a, spin="A"):
    """<T psi(theta, tau) psi^*(0, 0)> des Quartersektors a; antiperiodisch (A, KMS) oder periodisch (P) fortgesetzt."""
    u = dth + 1j*dta
    if dta > 0:
        return mp.fsum(mp.exp(1j*x*u)*cb for x, c, cb in MODE[a])
    v = -mp.fsum(mp.exp(1j*x*u)*c for x, c, cb in MODE[a])
    return v if spin == "A" else -v


def gram_ferm(pts, a, spiegel, spin="A"):
    n = len(pts)
    G = mp.matrix(n, n)
    for i in range(n):
        thr, tar = spiegel(pts[i])
        for j in range(n):
            G[i, j] = -S_f(pts[j][0] - thr, pts[j][1] - tar, a, spin)     # <psi^*(theta p_i) psi(p_j)>, Fermi-Vorzeichen
    return G


refl = lambda p: (p[0], BETA - p[1])
drehung = lambda p: (p[0], p[1] + BETA/2)
P1 = [(mp.mpf(rng.uniform(0, 2*np.pi)), mp.mpf(rng.uniform(0.4, 2*np.pi - 0.4))) for _ in range(14)]
P_str = P1[:7] + [(p[0], BETA - p[1]) for p in P1[7:]]
th = {}
QL = {QA[0]: "1/4", QA[1]: "3/4"}
for a in QA:
    th[f"A_{QL[a]}"] = herm_info(gram_ferm(P1, a, refl))
th["P_1/4"] = herm_info(gram_ferm(P1, QA[0], refl, "P"))
th["Drehung_1/4"] = herm_info(gram_ferm(P1, QA[0], drehung))
th["gespreizt_1/4"] = herm_info(gram_ferm(P_str, QA[0], refl))
ZAHLEN["rp_thermisch_fermion"] = {k: kurz(v) for k, v in th.items()}
say(f"   {ZAHLEN['rp_thermisch_fermion']}")
n1 = len(P1)
check(f"thermische Spiegelung, Fermion: die OS-Gram-Matrix <psi^*(theta p_i) psi(p_j)> aus der euklidischen Zweipunktfunktion "
      f"(14 Punkte auf Naht x halbem Kreis, Quarter 1/4 und 3/4) ist hermitesch und positiv definit (Trägheit "
      f"{th['A_1/4']['inertia']}, {th['A_3/4']['inertia']}, min. Eigenwert {float(th['A_1/4']['min']):.1e}); Kontrollen: "
      f"periodische Fortsetzung (falsche Spinstruktur) negativ definit {th['P_1/4']['inertia']}, Drehung um beta/2 statt "
      f"Spiegelung nicht hermitesch (Defekt {float(th['Drehung_1/4']['herm_defekt']):.2f}), über den Spiegel gespreizte "
      f"Punkte indefinit {th['gespreizt_1/4']['inertia']}: RP wählt die antiperiodische Struktur (KMS)",
      all(th[f"A_{QL[a]}"]["inertia"] == (n1, 0, 0) and th[f"A_{QL[a]}"]["herm_defekt"] < 1e-30 for a in QA)
      and th["P_1/4"]["inertia"] == (0, n1, 0) and th["Drehung_1/4"]["herm_defekt"] > 1e-3
      and th["gespreizt_1/4"]["inertia"][0] > 0 and th["gespreizt_1/4"]["inertia"][1] > 0)


def S_J(dth, dta, spin="P"):
    """Stromzweipunktfunktion sum_m e^(i m theta - m tau) m/(1 - e^(-beta m)), 0 < tau < beta, periodisch (Boson)."""
    if dta < 0:
        dta += BETA
        sg = 1 if spin == "P" else -1
    else:
        sg = 1
    u = dth + 1j*dta
    return sg*mp.fsum(mp.exp(1j*m*u)*m/(1 - mp.exp(-BETA*m)) for m in range(-170, 171) if m != 0)


def gram_strom(pts, spin="P"):
    n = len(pts)
    G = mp.matrix(n, n)
    for i in range(n):
        thr, tar = refl(pts[i])
        for j in range(n):
            G[i, j] = S_J(pts[j][0] - thr, pts[j][1] - tar, spin)
    return G


G1 = gram_ferm(P1, QA[0], refl)
paare = [(0, 1), (2, 3), (4, 5), (6, 7), (8, 9), (10, 11), (12, 13), (0, 7), (3, 10), (5, 12)]
G2 = mp.matrix(len(paare), len(paare))
for k, (p, q) in enumerate(paare):
    for l, (p2, q2) in enumerate(paare):
        G2[k, l] = G1[p, p2]*G1[q, q2] - G1[p, q2]*G1[q, p2]
sj = {"P": herm_info(gram_strom(P1[:12])), "A": herm_info(gram_strom(P1[:12], "A"))}
s2 = herm_info(G2)
ZAHLEN["rp_thermisch_strom_zweiteilchen"] = {"Strom": kurz(sj["P"]), "Strom antiperiodisch": kurz(sj["A"]),
                                              "psi psi (Wick)": kurz(s2)}
check(f"Strom J (h = 1) mit periodischer Fortsetzung: Gram positiv definit {sj['P']['inertia']}, antiperiodisch negativ "
      f"definit {sj['A']['inertia']}; Zweiteilchen-Operatoren psi(p) psi(q) über Wick (Slater-Determinanten der "
      f"Einteilchenform, 10 Paare): positiv definit {s2['inertia']}. Die thermische Spiegelung ist auf allen geprüften "
      "Sektoren des Originals reflexionspositiv, und RP wählt Bose-periodisch / Fermi-antiperiodisch (Spin-Statistik)",
      sj["P"]["inertia"] == (12, 0, 0) and sj["A"]["inertia"] == (0, 12, 0) and s2["inertia"] == (len(paare), 0, 0))

# --- 2.2 nu-Kern auf der 2x2-Matrixkugel unter der radialen Inversion Z -> (Z^+)^-1
def mpm(A):
    return mp.matrix([[mp.mpc(complex(A[i, j])) for j in range(A.shape[1])] for i in range(A.shape[0])])


def zuf_kugel(smin=0.15, smax=0.85):
    U = np.linalg.qr(rng.normal(size=(2, 2)) + 1j*rng.normal(size=(2, 2)))[0]
    V = np.linalg.qr(rng.normal(size=(2, 2)) + 1j*rng.normal(size=(2, 2)))[0]
    return U @ np.diag(rng.uniform(smin, smax, size=2)) @ V


def det2(A):
    return A[0, 0]*A[1, 1] - A[0, 1]*A[1, 0]


def kball(Z, W, nu, vz=1):
    A = vz*(Z*W.H)
    tr, de = A[0, 0] + A[1, 1], det2(A)
    d = mp.sqrt(tr**2 - 4*de)
    return mp.exp(-nu*(mp.log(1 - (tr + d)/2) + mp.log(1 - (tr - d)/2)))


ZK = [mpm(zuf_kugel()) for _ in range(24)]
I2m = mp.eye(2)
os_def = mp.mpf(0)
for Zi, Zj in itertools.product(ZK[:8], repeat=2):
    thZ = mp.inverse(Zj.H)
    os_def = max(os_def, abs(kball(Zi, Zj, 1) - 1/(det2(thZ - Zi)*det2(Zj.H)))/abs(kball(Zi, Zj, 1)))


def gram_kugel(nu, vz=1, pts=ZK):
    n = len(pts)
    G = mp.matrix(n, n)
    for i in range(n):
        for j in range(n):
            G[i, j] = kball(pts[i], pts[j], nu, vz)
    return G


kg = {str(nu): herm_info(gram_kugel(mp.mpf(nu))) for nu in ("0.5", "1", "1.5")}
kg["Antipoden"] = herm_info(gram_kugel(mp.mpf(1), vz=-1))
ZAHLEN["rp_kugel"] = {k: kurz(v) for k, v in kg.items()} | {"os_identitaet_rel": r(os_def, 3)}
say(f"   {ZAHLEN['rp_kugel']}")
check(f"radiale Inversion theta(Z) = (Z^+)^-1 (bildet die Kugel auf ihr Äußeres ab, Fixmenge U(2)): det(I - Z W^+)^-1 = "
      f"det(theta W - Z)^-1 det(W^+)^-1, die Kugel-Gram ist also die OS-Form (rel. Defekt {float(os_def):.0e}); 24 Punkte: "
      f"nu = 1 und 3/2 positiv definit ({kg['1']['inertia']}, {kg['1.5']['inertia']}), nu = 1/2 (außerhalb der Wallach-Menge) "
      f"indefinit {kg['0.5']['inertia']} (min. {float(kg['0.5']['min']):.1e}), Antipoden Z -> -(Z^+)^-1 (fixpunktfrei) "
      f"indefinit {kg['Antipoden']['inertia']}",
      os_def < 1e-30 and kg["1"]["inertia"] == (24, 0, 0) and kg["1.5"]["inertia"] == (24, 0, 0)
      and kg["0.5"]["inertia"][1] > 0 and kg["Antipoden"]["inertia"][1] > 0 and kg["Antipoden"]["inertia"][0] > 0)

# --- 2.3 euklidischer Schnitt: Halbraumspiegelung, Weyl-Spinor mit Theta_n, thermischer 4D-Kern
PAULI_MP = [mpm(s) for s in SIG]


def Qm(v):
    """Q(v) = v_0 I + i v.sigma (Quaternion als 2x2-Matrix), mpmath."""
    return v[0]*I2m + 1j*(v[1]*PAULI_MP[0] + v[2]*PAULI_MP[1] + v[3]*PAULI_MP[2])


XH = [[mp.mpf(rng.uniform(0.2, 1.5))] + [mp.mpf(t) for t in rng.uniform(-1, 1, size=3)] for _ in range(20)]


def gram_halbraum(nu, pts=XH, t=0):
    n = len(pts)
    G = mp.matrix(n, n)
    for i in range(n):
        for j in range(n):
            G[i, j] = ((pts[i][0] + pts[j][0] + t)**2 + sum((pts[i][k] - pts[j][k])**2 for k in (1, 2, 3)))**(-nu)
    return G


hr = {str(nu): herm_info(gram_halbraum(mp.mpf(nu))) for nu in ("0.5", "1")}
GH1 = gram_halbraum(1)
cay_def = mp.mpf(0)
for i, j in itertools.product(range(8), repeat=2):
    Xi = 1j*XH[i][0]*I2m + sum(XH[i][k + 1]*PAULI_MP[k] for k in range(3))
    Xj = 1j*XH[j][0]*I2m + sum(XH[j][k + 1]*PAULI_MP[k] for k in range(3))
    Zi, Zj = mp.inverse(Xi + 1j*I2m)*(Xi - 1j*I2m), mp.inverse(Xj + 1j*I2m)*(Xj - 1j*I2m)
    soll = det2(Xi + 1j*I2m)*mp.conj(det2(Xj + 1j*I2m))*GH1[i, j]/4
    cay_def = max(cay_def, abs(kball(Zi, Zj, 1) - soll)/abs(soll))
ZAHLEN["rp_halbraum"] = {k: kurz(v) for k, v in hr.items()} | {"cayley_kongruenz_rel": r(cay_def, 3)}
check(f"Halbraumspiegelung x_0 -> -x_0 auf dem euklidischen Schnitt (20 Punkte, x_0 > 0): |x - theta y|^(-2 nu) für nu = 1 "
      f"positiv definit {hr['1']['inertia']}, nu = 1/2 (unter der Unitaritätsschranke) indefinit {hr['0.5']['inertia']}; "
      f"über Cayley ist sie dieselbe Spiegelung wie die radiale Inversion: det(I - Z_i Z_j^+)^-1 = det(X_i + i) conj "
      f"det(X_j + i)/4 |x_i - theta x_j|^-2 (rel. Defekt {float(cay_def):.0e}), also Kongruenz und gleiche Trägheit",
      hr["1"]["inertia"] == (20, 0, 0) and hr["0.5"]["inertia"][1] > 0 and cay_def < 1e-30)

# Theta_n der externen Analyse
G22 = np.diag([1., 1., -1., -1.])
TAU79 = np.block([[np.zeros((2, 2)), I2], [I2, np.zeros((2, 2))]])
Qn = lambda v: v[0]*I2 + 1j*sum(v[k + 1]*SIG[k] for k in range(3))
Theta = lambda v: np.block([[np.zeros((2, 2)), Qn(v)], [Qn(v).conj().T, np.zeros((2, 2))]])


def su2():
    u = np.linalg.qr(rng.normal(size=(2, 2)) + 1j*rng.normal(size=(2, 2)))[0]
    return u/np.sqrt(np.linalg.det(u))


ok_th = np.allclose(Theta([1, 0, 0, 0]), TAU79)
for _ in range(8):
    nv = rng.normal(size=4)
    nv /= np.linalg.norm(nv)
    T_ = Theta(nv)
    ok_th &= np.allclose(T_ @ T_, np.eye(4)) and np.allclose(T_.conj().T @ G22 @ T_, -G22)
    gL, gR = su2(), su2()
    Dg = sla.block_diag(gL, gR)
    Qs = gL @ Qn(nv) @ gR.conj().T
    n2 = np.array([Qs[0, 0].real, Qs[0, 1].imag, Qs[0, 1].real, Qs[0, 0].imag])      # Q = [[n0 + i n3, n2 + i n1], ...]
    ok_th &= abs(np.linalg.norm(n2) - 1) < 1e-12 and np.allclose(Dg @ T_ @ np.linalg.inv(Dg), Theta(n2))
check("Theta_n = [[0, Q(n)], [Q(n)^+, 0]] mit Q(n) = n_0 I + i n.sigma, n in S^3 (acht zufällige n): Theta_n^2 = I, "
      "antiisometrisch für die CCR-Form diag(1, 1, -1, -1) von z79, Spin(4)-kovariant (diag(g_L, g_R) Theta_n "
      "diag(g_L, g_R)^-1 = Theta_n' mit Q(n') = g_L Q(n) g_R^+), und Theta_(1,0,0,0) ist genau z79s tau (a_i <-> b_i^+): "
      "die Theta_n sind die U(2) x U(2)-Bahn von tau, die Clifford-Multiplikation mit der Normalen n", ok_th)

XW = []
while len(XW) < 12:
    v = rng.normal(size=4)
    v *= rng.uniform(0.1, 0.8)/np.linalg.norm(v)
    XW.append([mp.mpf(t) for t in v])
nrm = lambda v: mp.sqrt(sum(t**2 for t in v))
S_W = lambda x, y: Qm([x[k] - y[k] for k in range(4)])/nrm([x[k] - y[k] for k in range(4)])**4


def gram_weyl(vz=1):
    n = len(XW)
    G = mp.matrix(2*n, 2*n)
    for i in range(n):
        xi = XW[i]
        rr = nrm(xi)
        xs = [vz*t/rr**2 for t in xi]
        L = Qm([t/rr for t in xi]).H/rr**3
        for j in range(n):
            B = L*S_W(xs, XW[j])
            for a_, b_ in itertools.product(range(2), repeat=2):
                G[2*i + a_, 2*j + b_] = B[a_, b_]
    return G


GW = gram_weyl()
geschl = mp.mpf(0)
for i, j in itertools.product(range(5), repeat=2):
    A = I2m - Qm(XW[i]).H*Qm(XW[j])
    soll = A/det2(A)**2
    geschl = max(geschl, max(abs(GW[2*i + a_, 2*j + b_] - soll[a_, b_]) for a_ in range(2) for b_ in range(2)))
wy = {"radial": herm_info(GW), "Antipoden": herm_info(gram_weyl(-1))}
sk = {}
for lab, vz in (("radial", 1), ("Antipoden", -1)):
    n = len(XW)
    G = mp.matrix(n, n)
    for i in range(n):
        for j in range(n):
            dot = sum(XW[i][k]*XW[j][k] for k in range(4))
            G[i, j] = 1/(1 - vz*2*dot + nrm(XW[i])**2*nrm(XW[j])**2)
    sk[lab] = herm_info(G)
ZAHLEN["rp_weyl_skalar_radial"] = {"Weyl radial": kurz(wy["radial"]), "Weyl Antipoden": kurz(wy["Antipoden"]),
                                   "Skalar radial": kurz(sk["radial"]), "Skalar Antipoden": kurz(sk["Antipoden"]),
                                   "geschlossene_Form_Defekt": r(geschl, 3)}
say(f"   {ZAHLEN['rp_weyl_skalar_radial']}")
check(f"radial gespiegelter Weyl-Spinor (Theta psi(x) = |x|^-3 Q(x/|x|)^+ psi(x/|x|^2), Propagator Q(x - y)/|x - y|^4, "
      f"12 Punkte, 24 x 24): hermitesch, positiv semidefinit {wy['radial']['inertia']} (min. {float(wy['radial']['min']):.1e}),"
      f" Blöcke = (I - X_i^+ X_j)/det(I - X_i^+ X_j)^2 (Defekt {float(geschl):.0e}); skalarer Vakuumkern radial "
      f"{sk['radial']['inertia']}; Antipoden x -> -x/|x|^2 indefinit (Weyl {wy['Antipoden']['inertia']}, Skalar "
      f"{sk['Antipoden']['inertia']}): die Behauptungen der externen Analyse halten",
      wy["radial"]["inertia"][1] == 0 and wy["radial"]["herm_defekt"] < 1e-30 and geschl < 1e-30
      and sk["radial"]["inertia"] == (12, 0, 0) and wy["Antipoden"]["inertia"][0] > 0 and wy["Antipoden"]["inertia"][1] > 0
      and sk["Antipoden"]["inertia"][0] > 0 and sk["Antipoden"]["inertia"][1] > 0)

PT = []
for _ in range(16):
    v = rng.normal(size=4)
    PT.append((mp.mpf(rng.uniform(0.3, 2*np.pi - 0.3)), [mp.mpf(t) for t in v/np.linalg.norm(v)]))


def gram_zyl(rot=False):
    n = len(PT)
    G = mp.matrix(n, n)
    for i in range(n):
        for j in range(n):
            c = sum(PT[i][1][k]*PT[j][1][k] for k in range(4))
            t = (PT[j][0] - PT[i][0] + BETA/2) if rot else (PT[i][0] + PT[j][0])
            G[i, j] = mp.fsum(1/(2*(mp.cosh(t + k*BETA) - c)) for k in range(-8, 9))
    return G


zy = {"Spiegelung": herm_info(gram_zyl()), "Drehung": herm_info(gram_zyl(True))}
ZAHLEN["rp_thermisch_4d"] = {k: kurz(v) for k, v in zy.items()}
check(f"thermischer 4D-Kern bei beta = 4 pi auf S^1_beta x S^3 (konformer Skalar, Bildsumme von 1/(2(cosh tau - cos "
      f"gamma))), Spiegelung tau -> -tau = radiale Inversion, 16 Punkte: positiv definit {zy['Spiegelung']['inertia']}; "
      f"Kontrolle Drehung um beta/2: indefinit {zy['Drehung']['inertia']}",
      zy["Spiegelung"]["inertia"] == (16, 0, 0) and zy["Drehung"]["inertia"][1] > 0)
ZAHLEN["laufzeit_bis_teil2_s"] = r(time.time() - T0, 3)

# =================================================================== Teil 3
say("\nTeil 3: OS-Rekonstruktion für die positiven Kandidaten")


def gen_eig(G, T, rel=mp.mpf(10)**-30):
    """Verallgemeinerte Eigenwerte von (T, G) auf dem Bild von G (Nullvektoren = OS-Kern werden herausgeteilt)."""
    n = G.rows
    E, Q = mp.eighe((G + G.H)/2)
    emax = max(abs(E[k]) for k in range(n))
    idx = [k for k in range(n) if mp.re(E[k]) > rel*emax]
    m = len(idx)
    W = mp.matrix(n, m)
    for c_, k in enumerate(idx):
        for i in range(n):
            W[i, c_] = Q[i, k]/mp.sqrt(mp.re(E[k]))
    M = W.H*T*W
    lam = mp.eighe((M + M.H)/2, eigvals_only=True)
    return sorted(mp.re(lam[k]) for k in range(m)), m


# --- 3.1 thermischer Kreis des Originals: der Transfer ist Delta^(-s/beta)
MOD6 = [mp.mpf(k) + QA[0] for k in (-3, -2, -1, 0, 1, 2)]                 # Trunkierung |x| < 3 des Quarters 1/4
tau0, s_tr = BETA/8, BETA/8


def os_ein(gewicht, tau_i, tau_j, dth):
    return mp.fsum(mp.exp(1j*x*dth)*mp.exp(-x*(tau_i + tau_j))*gewicht(x) for x in MOD6)


ths = [2*mp.pi*k/6 for k in range(6)]
rec = {}
for lab, gw in (("thermisch", lambda x: cb_f(x)), ("Vakuumlimes", lambda x: mp.mpf(1) if x > 0 else mp.mpf(0))):
    G = mp.matrix(6, 6)
    T = mp.matrix(6, 6)
    for i in range(6):
        for j in range(6):
            G[i, j] = os_ein(gw, tau0, tau0, ths[j] - ths[i])
            T[i, j] = os_ein(gw, tau0, tau0 + s_tr, ths[j] - ths[i])
    lam, m = gen_eig(G, T)
    rec[lab] = (sorted(-mp.log(l_)/s_tr for l_ in lam), lam, m)
xs_rec, lam_th, _ = rec["thermisch"]
def_th = max(abs(a_ - b_) for a_, b_ in zip(xs_rec, sorted(MOD6)))
per = max(abs(mp.exp(-2j*(4*mp.pi*x)) - 1) for x in xs_rec)
ZAHLEN["os_thermisch"] = {"generator_thermisch": [r(x, 6) for x in xs_rec], "max_lambda": r(max(lam_th), 5),
                          "generator_vakuumlimes": [r(x, 6) for x in rec["Vakuumlimes"][0]],
                          "rang_vakuumlimes": rec["Vakuumlimes"][2], "max|e^(-2iK) - 1|": r(per, 3)}
say(f"   {ZAHLEN['os_thermisch']}")
check(f"thermische OS-Rekonstruktion (Quarter 1/4, Trunkierung auf sechs Moden, OS-Vektoren auf der Naht bei tau = beta/8, "
      f"Verschiebung s = beta/8): der rekonstruierte Transfer hat die Eigenwerte e^(-s x), der Generator also das Spektrum "
      f"{[r(x, 4) for x in xs_rec]} = Z + 1/4 beidseitig (Abweichung {float(def_th):.0e}); max. Eigenwert "
      f"{float(max(lam_th)):.1f} > 1, keine Kontraktion; e^(-2iK) = 1 mit K = 4 pi x (periodisch, z72): die Rekonstruktion "
      f"gibt genau den Modulfluss zurück. Vakuumlimes (1 - c -> theta(x)): Rang {rec['Vakuumlimes'][2]}, Generator "
      f"{[r(x, 4) for x in rec['Vakuumlimes'][0]]} positiv und diskret (L_0 der Naht, die Uhr)",
      def_th < mp.mpf(10)**-25 and max(lam_th) > 1 and min(xs_rec) < 0 < max(xs_rec) and per < mp.mpf(10)**-25
      and rec["Vakuumlimes"][2] == 3 and all(x > 0 for x in rec["Vakuumlimes"][0]))
befund("Abbruchkriterium von S3 für das Original: Die einzige native reflexionspositive Spiegelung des KMS-Zustands (die "
       "thermische Klasse aus Teil 1) rekonstruiert den Liouville-Operator, also den periodischen Modulfluss (Klein-Landau, "
       "Tomita); positive Zeit entsteht erst im Vakuumlimes, und dort als diskrete Uhr L_0", True)

# --- 3.2 radiale Inversion auf der Kugel (nu = 1): Gradoperator, diskret, Vielfachheiten m^2
nv0 = [mp.mpf(t) for t in rng.normal(size=4)]
U0 = mp.expj(mp.mpf("0.7"))*Qm([t/nrm(nv0) for t in nv0])         # exakt unitär (Einheitsquaternion mal Phase)
rr = mp.mpf("0.6")
Zr = rr*U0
spek_def = mp.mpf(0)
for t in (mp.mpf("0.5"), mp.mpf(1), mp.mpf(2)):
    lhs = kball(Zr, mp.exp(-t)*Zr, 1)/kball(Zr, Zr, 1)
    rhs = mp.fsum((n + 1)*(1 - rr**2)**2*rr**(2*n)*mp.exp(-n*t) for n in range(400))
    spek_def = max(spek_def, abs(lhs - rhs))
ZS30 = []
for _ in range(30):
    A = rng.normal(size=(2, 2)) + 1j*rng.normal(size=(2, 2))
    ZS30.append(mpm(0.01*A/np.linalg.norm(A, 2)))
t_r = mp.mpf(1)
with mp.workdps(70):
    G = gram_kugel(1, pts=ZS30)
    T = mp.matrix(30, 30)
    for i in range(30):
        for j in range(30):
            T[i, j] = kball(ZS30[i], mp.exp(-t_r)*ZS30[j], 1)
    lam_r, _ = gen_eig(G, T, rel=mp.mpf(10)**-60)
    D_rec = sorted(-mp.log(l_)/t_r for l_ in lam_r)
cl = [D_rec[:1], D_rec[1:5], D_rec[5:14]]
ok_cl = all(abs(v - n) < 0.01 for n, grp in enumerate(cl) for v in grp) and all(abs(v - 3) < 0.2 for v in D_rec[14:30])
ZAHLEN["os_radial"] = {"spektralgesetz_defekt": r(spek_def, 3), "D_rekonstruiert": [r(x, 4) for x in D_rec]}
say(f"   {ZAHLEN['os_radial']}")
check(f"radiale OS-Rekonstruktion des nu = 1-Kerns: Transfer = Dilatation Z -> e^-t Z; im OS-Vektor von Z = 0,6 U gilt "
      f"exakt <Phi, T(t) Phi>/<Phi, Phi> = sum_n (n+1)(1 - r^2)^2 r^(2n) e^(-nt) (Defekt {float(spek_def):.0e}), ein "
      f"reines Punktmaß; auf 30 kleinen Punkten (|Z| = 0,01) liegen die rekonstruierten Generatorwerte in den Clustern "
      f"0 (1x), 1 (4x), 2 (9x), 3 (16x, Trunkierungsrand): L_0 = D + 1 mit Spektrum {{1, 2, 3, ...}} und Vielfachheiten m^2 (z73/z83), e^(2 pi i D) = 1: "
      "positiv, aber elliptisch (eine Uhr)", spek_def < mp.mpf(10)**-30 and ok_cl)

# --- 3.3 dieselbe Spiegelung mit dem Halbraum-Transfer: P_0, stetig, Gamma-Gesetz
gam_def = mp.mpf(0)
for x0 in (mp.mpf(1), mp.mpf("0.4")):
    for t in (mp.mpf("0.5"), mp.mpf(1), mp.mpf(4)):
        lhs = gram_halbraum(1, pts=[[x0, 0, 0, 0]], t=t)[0, 0]/gram_halbraum(1, pts=[[x0, 0, 0, 0]])[0, 0]
        rhs = mp.quad(lambda E: (2*x0)**2*E*mp.exp(-2*x0*E)*mp.exp(-t*E), [0, mp.inf])
        gam_def = max(gam_def, abs(lhs - rhs))
z72_def = max(abs(gram_halbraum(1, pts=[[1, 0, 0, 0]], t=t)[0, 0]/gram_halbraum(1, pts=[[1, 0, 0, 0]])[0, 0]
                  - (1 + t/2)**-2) for t in (mp.mpf("0.5"), mp.mpf(3)))
par = {}
with mp.workdps(90):
    for N in (6, 10, 14):
        a_s = [mp.mpf("0.3")*mp.mpf(10)**(mp.mpf(k)/(N - 1)) for k in range(N)]
        pts = [[a, 0, 0, 0] for a in a_s]
        lam_p, m_ = gen_eig(gram_halbraum(1, pts=pts), gram_halbraum(1, pts=pts, t=mp.mpf(1)), rel=mp.mpf(10)**-80)
        par[N] = (min(lam_p), max(lam_p), m_)
ZAHLEN["os_halbraum"] = {"gamma_gesetz_defekt": r(gam_def, 3), "z72_h1_defekt": r(z72_def, 3),
                         "kompressionen": {N: {"min": r(v[0], 4), "1 - max": r(1 - v[1], 4), "rang": v[2]}
                                           for N, v in par.items()}}
say(f"   {ZAHLEN['os_halbraum']}")
check(f"Halbraum-Transfer x_0 -> x_0 + t derselben Spiegelung: <Phi_x, T(t) Phi_x>/<Phi_x, Phi_x> = (2 x_0)^2/(2 x_0 + t)^2 "
      f"ist die Laplace-Transformierte der stetigen Dichte (2 x_0)^2 E e^(-2 x_0 E) auf (0, inf) (Gamma(2), Defekt "
      f"{float(gam_def):.0e}); bei x_0 = 1, dem Cayley-Bild des Kugelzentrums, ist es (1 + t/2)^-2 = z72s Gesetz "
      f"(1 - it/2)^(-2h) für h = 1 (Defekt {float(z72_def):.0e}); die Kompressionen auf N = 6, 10, 14 Punkten der Zeitachse "
      f"sind Kontraktionen (Eigenwerte in (0, 1)), und 1 - max fällt ({[r(1 - par[N][1], 3) for N in (6, 10, 14)]}): "
      "P_0 >= 0 ohne Lücke, der parabolische Generator von z72",
      gam_def < mp.mpf(10)**-30 and z72_def < mp.mpf(10)**-30
      and all(0 < par[N][0] and par[N][1] < 1 for N in par) and 1 - par[6][1] > 1 - par[10][1] > 1 - par[14][1])
befund("Auf dem Vakuumkern nu = 1 ist die radiale Inversion dieselbe Spiegelung wie die Halbraumspiegelung (Cayley); der "
       "rekonstruierte Hilbertraum ist derselbe. Verschieden ist nur der Transfer: Dilatation um das Kugelzentrum (das "
       "KMS-Ruhesystem kappa^2 I von z79) gibt die Uhr L_0, Translation zu einem Randpunkt gibt P_0. Die Wahl dieses "
       "Randpunkts ist z72s t*; die Spiegelung allein legt ihn nicht fest", True)
ZAHLEN["laufzeit_bis_teil3_s"] = r(time.time() - T0, 3)

# =================================================================== Teil 4
say("\nTeil 4: Folgt T'? Lorentzform, Spinbrücken-Kugel und Theta_n")
adj = lambda X: np.array([[X[1, 1], -X[0, 1]], [-X[1, 0], X[0, 0]]])
rho_t = lambda X: X.conj().T                                       # Lorentz-Realstruktur in Röhrenkoordinaten
JA_t = lambda X: -adj(X).conj().T                                  # euklidische (antipodale) Realstruktur
Inv_t = lambda X: -adj(X)                                          # komplexe Zeitspiegelung (holomorphe Inversion)


def nullebenen_typ(F, X0):
    """alpha-Ebene {X0 + u v^T}: bleibt der Spaltenraum (Erhalt) oder der Zeilenraum (Tausch) des Bildes fest?"""
    u = rng.normal(size=2) + 1j*rng.normal(size=2)
    bilder = []
    for _ in range(3):
        v = rng.normal(size=2) + 1j*rng.normal(size=2)
        D = F(X0 + np.outer(u, v)) - F(X0)
        uu, ss, vh = np.linalg.svd(D)
        bilder.append((uu[:, 0], vh[0]))
    gleich = lambda a, b: abs(abs(np.vdot(a, b)) - 1) < 1e-9
    if all(gleich(bilder[0][0], b[0]) for b in bilder):
        return "erhält"
    if all(gleich(bilder[0][1], b[1]) for b in bilder):
        return "tauscht"
    return "?"


ok_wick = True
for _ in range(6):
    X = rng.normal(size=(2, 2)) + 1j*rng.normal(size=(2, 2))
    ok_wick &= np.allclose(JA_t(Inv_t(X)), rho_t(X)) and np.allclose(Inv_t(Inv_t(X)), X) and np.allclose(JA_t(JA_t(X)), X)
typen = {nm: nullebenen_typ(F, rng.normal(size=(2, 2)) + 1j*rng.normal(size=(2, 2)))
         for nm, F in (("rho (Lorentz)", rho_t), ("J_A (Euklid)", JA_t), ("I (Inversion)", Inv_t))}


def det_signatur(basis):
    q = lambda A, B: (np.linalg.det(A + B) - np.linalg.det(A) - np.linalg.det(B))/2
    Bm = np.array([[q(A, B) for B in basis] for A in basis])
    ew = np.linalg.eigvalsh(Bm.real)
    return (int(np.sum(ew > 1e-9)), int(np.sum(ew < -1e-9))), float(np.abs(Bm.imag).max())


sig_rho = det_signatur([I2] + SIG)
sig_JA = det_signatur([1j*I2] + SIG)
ok_fix = all(np.allclose(rho_t(B), B) for B in [I2] + SIG) and all(np.allclose(JA_t(B), B) for B in [1j*I2] + SIG)
cay = lambda Z: 1j*(np.eye(2) + Z) @ np.linalg.inv(np.eye(2) - Z)
ok_cay = all(np.allclose(cay(np.linalg.inv(Z.conj().T)), cay(Z).conj().T) for Z in [zuf_kugel() for _ in range(6)])
ZAHLEN["realformen_roehre"] = {"Nullebenen": typen, "Signatur Fix(rho)": sig_rho[0], "Signatur Fix(J_A)": sig_JA[0]}
check(f"Röhrenkoordinaten: rho(X) = X^+ (Lorentz, Fixmenge Herm(2), det-Signatur {sig_rho[0]}), J_A(X) = -(adj X)^+ "
      f"(Euklid, Fixmenge i x_0 + x.sigma, Signatur {sig_JA[0]}), I(X) = -adj X (holomorphe Zeitspiegelung x_0 -> -x_0); "
      f"exakt rho = J_A o I, alle drei involutiv; Nullebenen: {typen}; die radiale Inversion (Z^+)^-1 der Kugel ist unter "
      "Cayley genau rho. Lorentz = Euklid o OS-Spiegelung: wo die Euklidstruktur und die Spiegelung nativ sind, folgt die "
      "faktortauschende Lorentzform", ok_wick and ok_fix and ok_cay and sig_rho[0] in ((1, 3), (3, 1))
      and sig_JA[0] in ((0, 4), (4, 0)) and typen == {"rho (Lorentz)": "tauscht", "J_A (Euklid)": "erhält",
                                                     "I (Inversion)": "tauscht"})

# --- 4.2 bi-unitäre Struktur der Spinbrücke aus der Markengruppe D_4 = <R, iota>, iota: z -> 1/z
IOTA, g_iota = None, None
for k, g in G_KAND.items():
    M = holo_transport(1/ZS, -sp.I/ZS**2, g)
    if M is not None:
        IOTA, g_iota = M, k
        break
IOTA_np = npm(IOTA)


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
zeilen = []
for A in (R_np, DECK_np, IOTA_np):
    zeilen.append(np.array([np.concatenate([(A.conj().T @ B @ A - B).real.ravel(), (A.conj().T @ B @ A - B).imag.ravel()])
                            for B in HB]).T)
_, s_, vh_ = np.linalg.svd(np.vstack(zeilen))
inv_formen = vh_[np.sum(s_ > 1e-9):]
h_inv = sum(c*B for c, B in zip(inv_formen[0], HB)) if len(inv_formen) else None
ew_h = np.linalg.eigvalsh(h_inv) if h_inv is not None else np.array([0])
h_inv = h_inv*np.sign(ew_h[-1])
hV, hEta = faktoren(h_inv)
hV, hEta = hV/hV[0, 0], hEta/hEta[0, 0]                         # definite Faktoren, positiv normiert
check(f"die Markengruppe <R, iota, Deck> (iota: z -> 1/z, Hebung erzwungen mit g = {g_iota}) lässt auf H0(N) genau "
      f"{len(inv_formen)} hermitesche Form zu; sie ist definit (Eigenwerte {np.round(np.sort(np.abs(ew_h)), 6).tolist()}) "
      f"und ein Produkt h_V (x) h_eta (Realignment-Rang {realign_rang(h_inv)}): V^* und eta tragen native unitäre Strukturen, "
      "die Kugel D(V, eta) = {Z in Hom(V, eta): Z^+ Z < I} existiert nativ bis auf eine Dilatationsskala",
      IOTA is not None and len(inv_formen) == 1 and (np.all(ew_h > 1e-9) or np.all(ew_h < -1e-9))
      and realign_rang(h_inv) == 1)

# --- 4.3 R auf der Spinbrücken-Kugel: kommutiert mit rho, wirkt frei auf dem Shilov-Rand U(2)
wV, wE = sla.sqrtm(hV), sla.sqrtm(hEta)
gVo = wV @ gV @ np.linalg.inv(wV)                                # orthonormale Basen
gEo = wE @ gEta @ np.linalg.inv(wE)
dVo = wE @ faktoren(DECK_np)[1] @ np.linalg.inv(wE)
gVo, gEo, dVo = [g/np.sqrt(abs(np.linalg.det(g))) for g in (gVo, gEo, dVo)]     # Faktoren mit |det| = 1
R_kug = lambda Z, gv=gVo, ge=gEo: ge @ Z @ gv.T                   # v -> (g_V (x) g_eta) v als Z -> g_eta Z g_V^T
rho_k = lambda Z: np.linalg.inv(Z.conj().T)
ok_rk = all(np.allclose(g.conj().T @ g, np.eye(2)) for g in (gVo, gEo))
for _ in range(6):
    Z = zuf_kugel()
    ok_rk &= np.allclose(rho_k(R_kug(Z)), R_kug(rho_k(Z))) and abs(np.linalg.norm(R_kug(Z), 2) - np.linalg.norm(Z, 2)) < 1e-9


def verhaeltnisse(g):
    e = np.linalg.eigvals(g)
    return {complex(np.round(e[0]/e[1], 9)), complex(np.round(e[1]/e[0], 9))}


frei = {}
for k in (1, 2, 3):
    for mit_deck in (False, True):
        ge = np.linalg.matrix_power(gEo, k) @ (dVo if mit_deck else np.eye(2))
        gv = np.linalg.matrix_power(gVo, k)
        frei[f"R^{k}{' Deck' if mit_deck else ''}"] = not (verhaeltnisse(ge) & verhaeltnisse(np.linalg.inv(gv).T))
ok_dr = True
for _ in range(4):
    U = np.linalg.qr(rng.normal(size=(2, 2)) + 1j*rng.normal(size=(2, 2)))[0]
    dr = lambda dZ, U=U: -U @ dZ.conj().T @ U
    ok_dr &= nullebenen_typ(dr, np.zeros((2, 2))) == "tauscht"
    tang = [U @ (1j*B) for B in [I2] + SIG]
    ok_dr &= all(np.allclose(dr(T_), T_) for T_ in tang)
    sg, _ = det_signatur([T_/np.sqrt(np.linalg.det(U)) for T_ in tang])
    ok_dr &= sg in ((1, 3), (3, 1))
ZAHLEN["spinbruecke_kugel"] = {"R frei auf U(2)": frei}
check(f"auf D(V, eta) ist R = (g_V, g_eta) bi-unitär und kommutiert mit der radialen Inversion rho(Z) = (Z^+)^-1; das "
      f"Differential von rho in jedem Randpunkt U ist antilinear, faktortauschend und hat auf T_U U(2) die Lorentzsignatur "
      f"(1,3) (T' in U, mit eta = V über U); fixpunktfrei auf U(2) für jede Phasennormierung (Eigenwertverhältnisse von "
      f"g_eta und g_V verschieden): {frei}. R selbst hat keinen Fixpunkt; nur die Halbdrehung mit Deck R^2 Deck hat welche. "
      "Kein Randpunkt ist unter <R> fest, also gibt es keine markeninvariante Minkowski-Karte und keine konstante, "
      "markeninvariante Identifikation T'",
      ok_rk and ok_dr and all(v for k, v in frei.items() if k != "R^2 Deck") and not frei["R^2 Deck"])

# --- 4.4 Theta_n: nativ ist auf der Oszillatorseite genau tau = Theta_1
stab = []
for A in [np.diag([1j, 0]), np.diag([0, 1j]), np.array([[0, 1], [-1, 0]], complex), np.array([[0, 1j], [1j, 0]])]:
    for which in (0, 1):
        X = sla.block_diag(A if which == 0 else np.zeros((2, 2)), A if which == 1 else np.zeros((2, 2)))
        stab.append(X @ TAU79 - TAU79 @ X)
Mst = np.array([np.concatenate([S.real.ravel(), S.imag.ravel()]) for S in stab]).T
dim_stab = 8 - np.linalg.matrix_rank(Mst, tol=1e-9)
ok_frei_S3 = True
for _ in range(6):
    nv = rng.normal(size=4)
    nv /= np.linalg.norm(nv)
    Qs = gVo @ Qn(nv) @ np.linalg.inv(gEo)
    ok_frei_S3 &= not proportional(Qs, Qn(nv))[0]
check(f"Theta_n ist nativ nur als Bahn: der Stabilisator von tau = Theta_1 in u(2) + u(2) hat Dimension {dim_stab} (die "
      "Diagonale, dieselbe u(2), die z79 als Symmetrie des Oszillatorzustands findet), der Zustand wählt also n = 1; auf der "
      "markierten Spinbrücke bildet R jedes Theta_n auf Theta_(R n) mit R n != n ab (keine Fixpunkte auf S^3), die Familie "
      "ist R-kovariant, kein Mitglied R-invariant", dim_stab == 4 and ok_frei_S3)
befund("T' als Folge der radialen Spiegelung: ja punktweise, nein markeninvariant. Die radiale Inversion vertauscht in jedem "
       "Randpunkt U (Minkowskipunkt) die beiden Spinorfaktoren über U und erzeugt dort die Lorentzform; sie ist mit R "
       "verträglich, aber R verschiebt jeden Punkt. z79s Negativbefund betrifft die konstante Identifikation eta = V-quer^*, "
       "die nach Satz B unmöglich bleibt; die konforme (kompakte) Lorentzstruktur auf U(2) ist dagegen markenverträglich. "
       "Prämissen: (P1) der nu = 1-Kern gehört zur Spinbrücken-Kugel D(V, eta) (Identifikation mit den Auslesungen der "
       "Oszillatorkugel, z74), (P2) eine Dilatationsskala", True)
ZAHLEN["laufzeit_bis_teil4_s"] = r(time.time() - T0, 3)

# =================================================================== Teil 5
say("\nTeil 5: Bilanz, Urteil zu S3, Bezug zu WOIT.OS.TWISTOR.01")
befund("Bezug WOIT.OS.TWISTOR.01 (v519/v522/v524/v565, nur gelesen): (1) Die Familien stimmen überein: Spiegelungen an "
       "Geraden (Familie D) invertieren R, Schwarz-Spiegelung und Antipoden (Familie A) zentralisieren R. (2) Für das "
       "Original ist die native RP-Spiegelung aber die thermische Klasse, deren Nahtschatten Familie A ist (Schwarz am "
       "Markenkreis): Die euklidische Zeit des KMS-Zustands ist der radiale/thermische Kreis, R = Delta^(-i/2) ist reelle "
       "Modulzeit und wird zentralisiert. v519 liest dagegen den Nahtwinkel als euklidische Zeit (Präzision (iii)); beide "
       "sind reflexionspositiv und beantworten verschiedene Fragen. (3) v524s Befund 'kompakte euklidische Zeit "
       "rekonstruiert eine KMS-Darstellung' gilt hier exakt: OS auf dem thermischen Kreis = Tomita, Transfer = "
       "Delta^(-s/beta). (4) v565s Rollentrennung (Metrik vs. absteigender Spiegel) entspricht theta_cut vs. J. (5) Neu für "
       "Kill-Test (7) auf freier Ebene: Auf H0(N) ist keine markenverträgliche Abbildung faktortauschend (Satz B), und R "
       "hat auf dem kompaktifizierten Minkowskiraum U(2) keinen Fixpunkt. Keine Ledgeränderung", True)
STATUS = {
    "Satz A: geometrische Spiegelungen der Naht wirken auf H0(N) = V^* (x) eta als Produkt (nie Faktortausch)": "E",
    "reelle Formen auf H0(N) aus markentreuen Spiegelungen: nur (2,2) oder keine; weder (4,0) noch (1,3)": "E",
    "Satz B: keine faktortauschende (anti)lineare Abbildung normalisiert <R> (Ordnungen 4 und 2); unter T' tut es J": "E",
    "Tomita-J, thermische Spiegelung, Teilchen-Loch-Spiegel: eine Klasse; RP auf Fermion, Strom, Zweiteilchen": "E",
    "OS-Rekonstruktion des Originals (KMS, beta = 4 pi): Liouville-Operator, zweiseitig, periodisch = Modulfluss":
        "E (Abbruchkriterium von S3 für das Original)",
    "nu = 1-Kern: radiale Inversion = Halbraumspiegelung (Cayley), RP genau auf der Wallach-Menge (nu = 1/2 indefinit)": "E",
    "radialer Transfer: L_0 = D + 1, Spektrum {1, 2, 3, ...}, Vielfachheiten m^2 (diskret, Uhr)": "E",
    "Halbraum-Transfer: P_0 >= 0 stetig ohne Lücke, Gamma(2)-Gesetz = z72 (h = 1)": "E",
    "Theta_n (externe Analyse): Theta^2 = I, Spin(4)-kovariant, RP für Skalar, Weyl und thermischen 4D-Kern; Antipoden "
    "indefinit; nativ als U(2) x U(2)-Bahn von z79s tau": "E",
    "Lorentz = Euklid o Inversion (rho = J_A o I) mit Faktortausch": "E",
    "Spinbrücken-Kugel D(V, eta) mit nativer bi-unitärer Struktur; rho kommutiert mit R; T' punktweise in jedem Randpunkt":
        "C (Prämissen P1: nu = 1-Kern auf D(V, eta); P2: Dilatationsskala)",
    "markeninvariantes T' (konstante Identifikation eta = V-quer^*)": "E (negativ)",
    "markeninvariante Minkowski-Karte bzw. markeninvariantes P_0 (Randpunkt t*)": "E (negativ: R ohne Fixpunkt auf U(2))",
    "Wahl des Randpunkts t* für die stetige Zeit, Dynamik, Wechselwirkung, 3+1-Schließung": "O",
}
ZAHLEN["status"] = STATUS
for k_, v_ in STATUS.items():
    say(f"   [{v_}] {k_}")
z = ZAHLEN
urteil = (
    "S3 hält geteilt. (1) Native reflexionspositive Spiegelungen gibt es: die thermische Klasse des Originals (Tomita-J = "
    "Spiegelung an beta/4, theta_cut = Spiegelung an 0 und beta/2, auf Einteilchenpropagatoren der Teilchen-Loch-Spiegel "
    "x -> -x, KMS-Identität 3e-40) ist auf Fermion, Strom und Zweiteilchen positiv definit (40 Stellen), und RP wählt "
    "Fermi-antiperiodisch und Bose-periodisch; auf dem Vakuumkern nu = 1 ist die radiale Inversion (Z^+)^-1 die OS-Form "
    "(positiv für nu = 1, 3/2, indefinit für nu = 1/2 und für die Antipoden). (2) Die Rekonstruktion des Originals gibt "
    f"genau den Modulfluss zurück (Generator {z['os_thermisch']['generator_thermisch']}, zweiseitig, e^(-2iK) = 1): das "
    "Abbruchkriterium von S3 greift für den KMS-Zustand. Positive Zeit entsteht nur im Vakuumsektor: dieselbe Spiegelung "
    "gibt mit dem Dilatationstransfer die Uhr L_0 (Spektrum 1, 2, 3, ... mit Vielfachheiten 1, 4, 9, 16) und mit dem "
    "Halbraumtransfer P_0 (stetig auf [0, inf), Gamma(2)-Gesetz = z72 mit h = 1). Der Randpunkt t* bleibt Wahl. (3) T' als "
    "konstante, markeninvariante Aussage folgt nicht: Alle zwölf markentreuen antiholomorphen Involutionen der Naht wirken "
    "auf H0(N) = V^* (x) eta als Produkt und geben nur die Signatur (2,2) oder keinen reellen Schnitt; keine "
    "faktortauschende Abbildung normalisiert R (projektive Ordnungen 4 und 2), während J unter der Prämisse T' auf "
    "V (x) V-quer^* tauscht und R zentralisiert. (4) Die externe Theta_n ist richtig gerechnet (Theta^2 = I, Spin(4)-"
    "kovariant, RP für Skalar, Weyl-Spinor und den thermischen 4D-Kern, Antipoden indefinit) und nativ als U(2) x U(2)-"
    "Bahn von z79s tau; der Oszillatorzustand wählt tau = Theta_1. Auf der Spinbrücke trägt die Markengruppe eine "
    "eindeutige bi-unitäre Form, die Kugel D(V, eta) existiert, die radiale Inversion kommutiert mit R und erzeugt in "
    "jedem Randpunkt die Lorentzform mit Faktortausch (T' punktweise); R hat aber auf U(2) keinen Fixpunkt, also gibt es "
    "keine markeninvariante Minkowski-Karte und kein markeninvariantes T'. Die Marken wirken als fixpunktfreie konforme "
    "Abbildung, nicht als Lorentztransformation. Prämissen: nu = 1-Kern auf D(V, eta) (P1), eine Dilatationsskala (P2). "
    "Nichts befördert; Ledger unverändert.")
say("\nUrteil\n   " + urteil)
ZAHLEN["laufzeit_gesamt_s"] = r(time.time() - T0, 3)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)}); Laufzeit "
    f"{time.time() - T0:.0f} s")
(HERE / "z84_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
