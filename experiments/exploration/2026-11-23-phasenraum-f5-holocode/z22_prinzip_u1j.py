"""Z22: Prinzip U — die Deck-Z4 als Untergruppe einer inneren Drehung U(1)_J jedes Paars (psi, chi).

Hypothese U. J = exp(pi T/2) liegt in SO(2)_J = U(1)_J, die jedes Paar (psi, chi) dreht. xi_pm = (psi -+ i chi)/sqrt 2
tragen die Ladung +-1. Das Paarungsfeld Phi hat die feste Ladung -2 und koppelt an xi+ xi+ (plus h.c.). Das
Wilson-Glied W xi+ xi- hat Ladung 0 und ist erlaubt. In den eta-Variablen von Z7: eta+ windet wie bisher, eta-
hat die Majoranamasse m_minus = 0, das Wilson-Glied eta+ eta- bleibt.

Geprüft wird, ob U die vier offenen Entscheidungen von H6 zugleich löst.
 0. Ladungen aller Bilinearen; welche Yukawas U(1)_J erlaubt; ein Symmetriesatz für Phi* xi- xi-.
 3. Paarung: psi chi gegen xi+ xi+; Familienmassen im Z18-Kanal (q = f2^2 - f0^2).
 4. Gitter: Phasenpotential (wie Z9 Teil B, euklidisch wie Z19), [H, Q_J], Reflexionspositivität (wie Z19).
 I. Index der Naht mit lückenlosem xi-: Zählung am Kern (Z7-Methode), Kontrollmassen, Hybridisierung,
    Orientierung der Halbwindung und ihre Fermionenergie.
 1. Eine Naht: Holonomie als Viertelfluss, Minimum über Z4-wertige Schläuche, Z11-Parität.
 2. Ein Block: U(1)_J-Anomalien, Spin-Z4 (Z16), gerade unimodulare Klebung.
 D. Darstellungen: invariante Bilinearformen für die 16 von Spin(10), den Vektor von SO(16), SU(5).
"""
import itertools
import json
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sps
import sympy as sp
from scipy.linalg import eigh, orth
from scipy.sparse.linalg import eigsh

HERE = Path(__file__).parent
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
T0 = time.time()


def say(s=""):
    print(s, flush=True)


def check(name, cond):
    CHECKS[name] = bool(cond)
    return "ja" if cond else "NEIN"


def befund(name, val):
    BEFUNDE[name] = bool(val)
    return "ja" if val else "nein"


say("Z22: Prinzip U — Deck-Z4 als Untergruppe der inneren Drehung U(1)_J")
say("Masse m = 0,6, Wilson r = 1, Querschnitt mit Kernradius 3,5 wie Z7.\n")

# =====================================================================================================
say("0. Ladungen unter U(1)_J\n")
al = sp.symbols("alpha", real=True)
Tg = sp.Matrix([[0, 1], [-1, 0]])                     # exp(alpha T) dreht (psi, chi); J = exp(pi T/2)
Jm = sp.simplify((al * Tg).exp().subs(al, sp.pi / 2))
say(f"   J = exp(pi T/2) = {Jm.tolist()}: psi -> chi -> -psi wie Z7.")
vp = sp.Matrix([1, -sp.I]) / sp.sqrt(2)                # xi+ = vp^T (psi, chi)
vm = sp.Matrix([1, sp.I]) / sp.sqrt(2)
q_xi = [sp.simplify((Tg.T * v)[0] / v[0]) for v in (vp, vm)]
say(f"   T^T v = i q v: xi+ hat q = {sp.simplify(q_xi[0] / sp.I)}, xi- hat q = {sp.simplify(q_xi[1] / sp.I)}.")
# Bilineare psi_i B_ij psi_j, B symmetrisch; Wirkung delta B = T^T B + B T
b11, b22, b12 = sp.symbols("b11 b22 b12")
Bgen = sp.Matrix([[b11, b12], [b12, b22]])
dB = Tg.T * Bgen + Bgen * Tg
A = sp.Matrix([[sp.diff(e, v) for v in (b11, b22, b12)] for e in (dB[0, 0], dB[1, 1], dB[0, 1])])
ev = A.eigenvects()
charges = sorted([sp.simplify(e / sp.I) for e, mult, vecs in ev for _ in range(mult)], key=lambda x: float(x))
say(f"   Ladungen der drei Bilinearen (Eigenwerte von ad_T / i): {charges}")
Bpp, Bmm, Bpm = vp * vp.T, vm * vm.T, (vp * vm.T + vm * vp.T) / 2
ad = lambda B: sp.simplify(Tg.T * B + B * Tg)
ok_q = (sp.simplify(ad(Bpp) - 2 * sp.I * Bpp) == sp.zeros(2) and sp.simplify(ad(Bmm) + 2 * sp.I * Bmm) == sp.zeros(2)
        and ad(Bpm) == sp.zeros(2))
say(f"   xi+ xi+ hat +2, xi- xi- hat -2, xi+ xi- = (psi psi + chi chi)/2 = {sp.simplify(Bpm).tolist()} hat 0.")
check("U(1)_J-Ladungen der Bilinearen: xi+xi+ +2, xi-xi- -2, xi+xi- (Wilson) 0", ok_q and charges == [-2, 0, 2])
Bpc = sp.Matrix([[0, 1], [1, 0]]) / 2                   # psi chi
coef = sp.solve(sp.Matrix(Bpc - sp.Symbol("a") * Bpp - sp.Symbol("b") * Bmm - sp.Symbol("c") * Bpm),
                [sp.Symbol("a"), sp.Symbol("b"), sp.Symbol("c")], dict=True)[0]
a_, b_, c_ = [coef[sp.Symbol(s)] for s in "abc"]
say(f"   psi chi = ({a_}) xi+xi+ + ({b_}) xi-xi- + ({c_}) xi+xi-: Ladung +2 und -2 je zur Hälfte.")
not_cov = sp.Matrix([list(ad(Bpc)), list(Bpc)]).rank() == 2       # ad(B) nicht proportional zu B
check("psi chi (G) ist nicht U(1)_J-kovariant: Anteile +2 und -2 gleich groß",
      not_cov and sp.simplify(abs(a_) - abs(b_)) == 0 and c_ == 0)
# erlaubte Yukawas für ein komplexes Phi der Ladung c_Phi = -2
cPhi = -2
allowed_Phi = [nm for nm, qb in (("xi+xi+", 2), ("xi-xi-", -2), ("xi+xi-", 0)) if cPhi + qb == 0]
allowed_Phic = [nm for nm, qb in (("xi+xi+", 2), ("xi-xi-", -2), ("xi+xi-", 0)) if -cPhi + qb == 0]
say(f"   Phi (Ladung -2) darf an {allowed_Phi} koppeln, Phi* (Ladung +2) an {allowed_Phic}.")
say("   Die allgemeinste U(1)_J-invariante Kopplung ist y Phi xi+xi+ + y' Phi* xi-xi- + W xi+xi- + h.c.")
check("Phi koppelt nur an xi+xi+; U(1)_J lässt genau zwei Yukawas zu (Phi xi+xi+, Phi* xi-xi-)",
      allowed_Phi == ["xi+xi+"] and allowed_Phic == ["xi-xi-"])
befund("U(1)_J verbietet Phi* xi- xi- nicht (y' bleibt frei)", allowed_Phic == ["xi-xi-"])
# Satz für jede abelsche Symmetrie
qp, qm, qf = sp.symbols("q_p q_m q_Phi")
sol = sp.solve([qf + 2 * qp, qp + qm], [qf, qm], dict=True)[0]
thm = sp.simplify((-qf + 2 * qm).subs(sol))
say(f"   Satz: Erhält eine Ladung W xi+xi- (q_p + q_m = 0) und Phi xi+xi+ (q_Phi + 2 q_p = 0), so hat Phi* xi-xi-")
say(f"   die Ladung -q_Phi + 2 q_m = {thm}. Nichtabelsch gilt dasselbe: W verlangt xi- ~ konj(xi+), also")
say("   Phi* (x) Sym^2(xi-) = konj(Phi (x) Sym^2(xi+)). Keine Symmetrie, die Wilson zulässt, trennt y' von y.")
check("Satz: jede Ladung, die W xi+xi- und Phi xi+xi+ erhält, erhält Phi* xi-xi-", thm == 0)
befund("keine Symmetrie trennt y' von y, solange das Wilson-Glied erlaubt ist", thm == 0)

# vom Wilson-Glied induzierte xi- Masse (Schur-Komplement) im freien Gitter
m0 = 0.6
sx = np.array([[0, 1], [1, 0]], complex)
sy = np.array([[0, -1j], [1j, 0]])
sz = np.diag([1.0, -1.0]).astype(complex)
isy = 1j * sy


def bulk_bdg(k, Phi, mm, rw=1.0):
    hs = np.sin(k[0]) * sx + np.sin(k[1]) * sy + np.sin(k[2]) * sz
    W = rw * np.sum(1 - np.cos(k))
    M = np.array([[Phi, W], [W, mm]])
    D = np.kron(M, isy)
    hb = np.kron(np.eye(2), hs)
    return np.block([[hb, D], [D.conj().T, hb.T]]), W


phi0 = 0.7
dev, ratio = [], []
for kzz in (0.05, 0.1, 0.3, 0.8):
    Hb, W = bulk_bdg(np.array([0.0, 0.0, kzz]), m0 * np.exp(1j * phi0), 0.0)
    E = np.linalg.eigvalsh(Hb)
    smin = (np.sqrt(m0 ** 2 + 4 * W ** 2) - m0) / 2
    dev.append(abs(E[E > 0].min() - np.sqrt(np.sin(kzz) ** 2 + smin ** 2)))
    ratio.append(smin * m0 / W ** 2)
Meff = -(1.0 ** 2) / (m0 * np.exp(1j * phi0))
say(f"   Wilson induziert für xi- die Masse M_eff = -W^2/Phi (Schur): Phase {np.angle(Meff):+.3f} = pi - phi0 "
    f"({np.pi - phi0:+.3f}), also Phi*-artig.")
say(f"   Im Volumen gilt E_min(k_z) = sqrt(sin^2 k_z + s^2) mit s = leichter Singulärwert (Abweichung {max(dev):.1e});")
say(f"   s m / W^2 bei k_z = 0,05/0,1/0,3/0,8: {np.round(ratio, 4)}. xi- bleibt am Knoten masselos, s ~ k^4/m.")
check("freies Gitter: xi- am Knoten masselos, induzierte Masse W^2/m ~ k^4 mit Phase von Phi*",
      max(dev) < 1e-12 and abs(ratio[0] - 1) < 1e-3 and abs(np.angle(Meff) - (np.pi - phi0)) < 1e-12)
ZAHLEN["induzierte_xi_minus_masse_s_m_durch_W2"] = [float(x) for x in ratio]

# eta-Rahmen: Windungen bei Halbwindung nu = +-1/2 von Phi um den J-Schlauch
say("\n   Windungen im eta-Rahmen (eta_pm = e^{-+i phi/4} xi_pm, Holonomie i auf xi+):")
say("   Delta+ = Phi e^{i phi/2}: w+ = nu + 1/2. G: Delta- ~ Phi e^{-i phi/2}, w- = nu - 1/2. U mit y': Delta- ~ Phi* e^{-i phi/2},")
say("   w- = -nu - 1/2. Vorhergesagter Index = Summe der Windungen der massiven Arten (det M am Knoten).")
pred = {}
for nu in (sp.Rational(1, 2), -sp.Rational(1, 2)):
    wp = nu + sp.Rational(1, 2)
    tag = "ausgerichtet" if nu > 0 else "gegen"
    pred[f"G {tag}"] = (int(wp), int(nu - sp.Rational(1, 2)), int(wp + nu - sp.Rational(1, 2)))
    pred[f"U0 {tag}"] = (int(wp), None, int(wp))
    pred[f"Uy {tag}"] = (int(wp), int(-nu - sp.Rational(1, 2)), int(wp - nu - sp.Rational(1, 2)))
for k_, v_ in pred.items():
    say(f"      {k_:16s} w+ = {v_[0]:+d}, w- = {'--' if v_[1] is None else f'{v_[1]:+d}'}, Index {v_[2]:+d}")
say("   (ausgerichtet: Halbwindung nu = +1/2 im Sinn der Holonomie; gegen: nu = -1/2, gleiche Gradientenenergie nu^2)")

# =====================================================================================================
say("\n3. Paarung und Familienmassen im Z18-Kanal\n")
z = sp.symbols("z")
ident = sp.simplify(z ** 4 / (z ** 4 - 1) ** 2 - 1 / (z ** 4 - 1) ** 2 - 1 / (z ** 4 - 1))
say(f"   Z18-Identität q = f2^2 - f0^2: Rest {ident}. Kanäle aus Z18: <q,q> = Q, <q,f2^2> = Q/2, <q,f0^2> = -Q/2,")
say("   alle anderen Familienkanäle 0; Naht-Familie-Kanäle verschwinden (sqrt q deckungerade, q deckgerade).")
say("   Eichbedingung (Z18, eps = +1): Naht in xi+, Familien f0, f1, f2 in xi-.")
Qn = 1.0
Cch = np.zeros((4, 4), complex)                      # Basis (f0, f1, f2, Naht)
Cch[3, 3], Cch[2, 2], Cch[0, 0] = Qn, Qn / 2, -Qn / 2
ep, em = np.array([1, 0]), np.array([0, 1])          # Arten (xi+, xi-)


def phys_masses(Bsp):
    full = np.kron(Bsp, Cch)                          # (Art, Mode)
    seam = 0 * 4 + 3
    fam = [1 * 4 + j for j in range(3)]
    ms = abs(full[seam, seam])
    sv = np.sort(np.linalg.svd(full[np.ix_(fam, fam)], compute_uv=False))
    cross = np.abs(full[np.ix_([seam], fam)]).max()
    return ms, sv / ms, cross


B_G = 0.5j * (np.outer(ep, ep) - np.outer(em, em))
B_U = np.outer(ep, ep)
yp = 0.3
B_Uy = np.outer(ep, ep) + yp * np.outer(em, em)
mG, fG, xG = phys_masses(B_G)
mU, fU, xU = phys_masses(B_U)
mUy, fUy, xUy = phys_masses(B_Uy)
z18 = json.loads((HERE / "z18_ergebnis.json").read_text())["data"]["paarung"]["familienmassen_psichi"]
say(f"   G (psi chi):      Familienmassen / Nahtmasse = {np.round(fG, 6)} (Z18-Datei: {z18})")
say(f"   U, y' = 0:        {np.round(fU, 6)}")
say(f"   U, y' = {yp}:      {np.round(fUy, 6)} = |y'| x {{0, 1/2, 1/2}}")
check("Z18 reproduziert: psi chi gibt Familienmassen {0, 1/2, 1/2} der Nahtmasse",
      ident == 0 and np.allclose(fG, z18, atol=1e-12))
check("U mit y' = 0: alle drei Familien exakt masselos, Naht massiv", np.allclose(fU, 0) and mU > 0)
befund("mit y' != 0 erhalten die Familien |y'| x {0, 1/2, 1/2} der Nahtmasse", np.allclose(fUy, yp * np.array([0, .5, .5])))
say("   Das Profil von Phi* ist ebenfalls ~ q (reelle Charaktere von q), daher derselbe Kanal mit Faktor y'.")
ZAHLEN["familienmassen"] = {"G": fG.tolist(), "U_y0": fU.tolist(), f"U_y{yp}": fUy.tolist()}

# =====================================================================================================
say("\n4. Gitter: Phasenpotential, Erhaltung von Q_J, Reflexionspositivität\n")


def e_density_ham(phi, m=0.6, mm=0.6, conj_minus=False, n=12, rw=1.0):
    k = 2 * np.pi * (np.arange(n) + 0.5) / n
    Kx, Ky, Kz = [a.ravel() for a in np.meshgrid(k, k, k, indexing="ij")]
    hs = np.sin(Kx)[:, None, None] * sx + np.sin(Ky)[:, None, None] * sy + np.sin(Kz)[:, None, None] * sz
    W = rw * (3 - np.cos(Kx) - np.cos(Ky) - np.cos(Kz))
    Mm = np.zeros((len(W), 2, 2), complex)
    Mm[:, 0, 0] = m * np.exp(1j * phi)
    Mm[:, 1, 1] = mm * (np.exp(-1j * phi) if conj_minus else 1.0)
    Mm[:, 0, 1] = Mm[:, 1, 0] = W
    Dm = np.einsum("nij,ab->niajb", Mm, isy).reshape(-1, 4, 4)
    hb = np.einsum("ij,nab->niajb", np.eye(2), hs).reshape(-1, 4, 4)
    H = np.block([[hb, Dm], [np.conj(np.swapaxes(Dm, 1, 2)), np.swapaxes(hb, 1, 2)]])
    E = np.linalg.eigvalsh(H)
    return float(-0.5 * np.where(E > 0, E, 0).sum(axis=1).mean())


phis = np.linspace(0, 2 * np.pi, 9)
pot = {}
for lab, kw in (("G (Z9, m_minus = m)", dict(mm=0.6)), ("U (m_minus = 0)", dict(mm=0.0)),
                ("U mit y' = 0,6 (Phi*)", dict(mm=0.36, conj_minus=True))):
    v = np.array([e_density_ham(p, **kw) for p in phis])
    pot[lab] = v - v.min()
    say(f"   Hamiltonsch {lab:24s} V - V_min = " + " ".join(f"{x:.4f}" for x in pot[lab]) + f"  (Spannweite {np.ptp(v):.1e})")
spG, spU, spUy = [float(np.ptp(v)) for v in pot.values()]
check("Phasenpotential: Z9-Barriere 0,0882 reproduziert; U und U mit y' exakt flach (< 1e-12)",
      abs(spG - 0.0882) < 5e-4 and spU < 1e-12 and spUy < 1e-12)
ZAHLEN["phasenpotential_spannweite"] = {"G": spG, "U": spU, "U_y": spUy}

# euklidisch wie Z19 (W-Form)
gam = [np.array(x, complex) for x in (
    [[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]],
    [[0, 0, 0, -1j], [0, 0, -1j, 0], [0, 1j, 0, 0], [1j, 0, 0, 0]],
    [[0, 0, 0, -1], [0, 0, 1, 0], [0, 1, 0, 0], [-1, 0, 0, 0]],
    [[0, 0, -1j, 0], [0, 0, 0, 1j], [1j, 0, 0, 0], [0, -1j, 0, 0]])]
G4 = gam[3]
g5 = gam[0] @ gam[1] @ gam[2] @ gam[3]
I4 = np.eye(4)
P_p, P_m = (I4 + g5) / 2, (I4 - g5) / 2
k4 = 2 * np.pi * (np.arange(8) + 0.5) / 8
K8 = np.array(np.meshgrid(k4, k4, k4, k4, indexing="ij")).reshape(4, -1).T


def V_eukl(phi, mm, conj_minus=False, m=0.6):
    W = np.sum(1 - np.cos(K8), axis=1)
    Mk = np.zeros((len(K8), 2, 2), complex)
    Mk[:, 0, 0], Mk[:, 0, 1], Mk[:, 1, 0] = m * np.exp(1j * phi), W, W
    Mk[:, 1, 1] = mm * (np.exp(-1j * phi) if conj_minus else 1.0)
    Dn = 1j * np.einsum("nm,mab->nab", np.sin(K8), np.array(gam))
    Dt = np.einsum("ij,nab->niajb", np.eye(2), Dn).reshape(-1, 8, 8) + (
        np.einsum("nij,ab->niajb", Mk, P_p) + np.einsum("nij,ab->niajb", Mk.conj(), P_m)).reshape(-1, 8, 8)
    return -0.5 * np.linalg.slogdet(Dt)[1].mean()


vE_G = np.array([V_eukl(p, 0.6) for p in phis])
vE_U = np.array([V_eukl(p, 0.0) for p in phis])
say(f"   Euklidisch (Z19, W-Form, 8^4): G Spannweite {np.ptp(vE_G):.4f} (Z19: 0,083), U Spannweite {np.ptp(vE_U):.1e}")
check("Euklidisch: U flach, G mit Z19-Barriere", np.ptp(vE_U) < 1e-12 and abs(np.ptp(vE_G) - 0.083) < 2e-3)


# ---- Gitter mit Defekt: Builder aus Z7 mit eigener eta- Masse (vektorisiert, dünnbesetzt)
def geom(L):
    c = (L - 1) / 2.0
    xs, ys = np.meshgrid(np.arange(L), np.arange(L), indexing="ij")
    return np.arctan2(ys - c, xs - c), np.hypot(xs - c, ys - c)


def build(L, kz, n_plus, n_minus, flux, m=0.6, m_minus=0.6, conj_minus=False, phi0=0.0, rw=1.0):
    """Z7-build (eine Kopie). eta- Masse m_minus e^{i(n_minus phi + s phi0)}, s = -1 für die Phi*-Kopplung."""
    phi, rad = geom(L)
    prof = np.tanh(rad / 1.5)
    ns, nh = 2, L * L * 4
    I = lambda si, a, s: (si * ns + a) * 2 + s
    R, C, V = {"h": [], "z": [], "D": []}, {"h": [], "z": [], "D": []}, {"h": [], "z": [], "D": []}

    def add(key, r, c, v):
        r, c, v = np.broadcast_arrays(np.asarray(r).ravel(), np.asarray(c).ravel(), np.asarray(v).ravel())
        R[key].append(r)
        C[key].append(c)
        V[key].append(v.astype(complex))

    S = np.arange(L * L)
    ph, pr = phi.ravel(), prof.ravel()
    sg = -1.0 if conj_minus else 1.0
    for a in range(ns):
        for s in range(2):
            add("z", I(S, a, s), I(S, a, s), np.full(S.shape, np.sin(kz) * sz[s, s]))
        Md = m * pr * np.exp(1j * (n_plus * ph + phi0)) if a == 0 else m_minus * pr * np.exp(1j * (n_minus * ph + sg * phi0))
        add("D", I(S, a, 0), I(S, a, 1), Md * isy[0, 1])
        add("D", I(S, a, 1), I(S, a, 0), Md * isy[1, 0])
    Wos = rw * (3 - np.cos(kz))
    for (s, t) in ((0, 1), (1, 0)):
        add("D", I(S, 0, s), I(S, 1, t), np.full(S.shape, Wos * isy[s, t]))
        add("D", I(S, 1, s), I(S, 0, t), np.full(S.shape, Wos * isy[s, t]))
    for (dx, dy, sig) in ((1, 0, sx), (0, 1, sy)):
        x1, y1 = [u.ravel() for u in np.meshgrid(np.arange(L - dx), np.arange(L - dy), indexing="ij")]
        x2, y2 = x1 + dx, y1 + dy
        s1, s2 = x1 * L + y1, x2 * L + y2
        dphi = np.angle(np.exp(1j * (phi[x2, y2] - phi[x1, y1])))
        for a in range(ns):
            pph = np.exp(-1j * (flux if a == 0 else -flux) * dphi)
            for s in range(2):
                for t in range(2):
                    if sig[s, t] != 0:
                        val = sig[s, t] / (2j) * pph
                        add("h", I(s1, a, s), I(s2, a, t), val)
                        add("h", I(s2, a, t), I(s1, a, s), np.conj(val))
        for (aa, bb, q) in ((0, 1, +1), (1, 0, -1)):
            pw = -rw / 2 * np.exp(-1j * q * flux * dphi)
            for (s, t) in ((0, 1), (1, 0)):
                add("D", I(s1, aa, s), I(s2, bb, t), pw * isy[s, t])
                add("D", I(s2, bb, t), I(s1, aa, s), pw * isy[t, s])
    mk = lambda key: sps.coo_matrix((np.concatenate(V[key]), (np.concatenate(R[key]), np.concatenate(C[key]))),
                                    shape=(nh, nh)).tocsr()
    h, hz, D = mk("h"), mk("z"), mk("D")
    return sps.bmat([[h + hz, D], [D.conj().T, -(h - hz).T]]).tocsr(), rad.ravel()


def masks(L, rad, R=3.5):
    core = np.tile(np.repeat(rad <= R, 4), 2).astype(float)
    spec = np.tile(np.tile(np.repeat(np.arange(2), 2), L * L), 2)
    qJ = np.where(spec == 0, 1.0, -1.0)
    qJ[len(qJ) // 2:] *= -1                                  # Lochblock: -q
    return core, (spec == 1).astype(float), qJ


# Vergleich mit Z7s build
z7p = HERE / "z7_ursprung_naht.py"
ns7 = {"__file__": str(z7p)}
src7 = z7p.read_text()
exec(src7[:src7.index("L = 14")], ns7)
bdiff = max(np.abs(ns7["build"](*args)[0] - build(*args)[0].toarray()).max()
            for args in ((8, 0.37, 1, 0, 0.25), (7, -1.1, 1, 1, 0.0), (8, 0.2, 0, 0, 0.25)))
say(f"\n   Eigener dünnbesetzter Builder gegen Z7-build (m_minus = m): größte Differenz {bdiff:.1e}")
check("Builder identisch mit Z7-build", bdiff == 0.0)

Lc = 12
H0, radc = build(Lc, 0.4, 1, 0, 0.25, m=0.0, m_minus=0.0)
_, _, qJ = masks(Lc, radc)
comm0 = np.abs(H0.multiply(qJ[:, None]) - H0.multiply(qJ[None, :])).max()
alpha = 0.37
res = {}
for lab, kw in (("U", dict(n_minus=0, m_minus=0.0)), ("Uy", dict(n_minus=-1, m_minus=0.36, conj_minus=True)),
                ("G", dict(n_minus=0, m_minus=0.6))):
    Ha, _ = build(Lc, 0.4, n_plus=1, flux=0.25, phi0=0.3, **kw)
    Hb, _ = build(Lc, 0.4, n_plus=1, flux=0.25, phi0=0.3 + 2 * alpha, **kw)
    ph_ = np.exp(1j * alpha * qJ)
    Hrot = Ha.multiply(ph_[:, None]).multiply(np.conj(ph_)[None, :])
    res[lab] = float(np.abs(Hrot - Hb).max())
say(f"   Defektgitter {Lc}x{Lc}, Viertelfluss: max |[H, Q_J]| ohne Phi = {comm0:.1e}")
say(f"   e^(i a Q_J) H(phi0) e^(-i a Q_J) - H(phi0 + 2a): U {res['U']:.1e}, U mit y' {res['Uy']:.1e}, G {res['G']:.2f}")
check("Q_J erhalten ohne Phi; mit Phi dreht Q_J die Phase (U, U mit y'), bei G nicht",
      comm0 < 1e-14 and res["U"] < 1e-14 and res["Uy"] < 1e-14 and res["G"] > 0.1)
say("   Also ist das Potential der gleichförmigen Phase exakt flach: Drehungen sind unitär äquivalent.")
say("   Die Z9-Barriere stammt aus Phi xi-xi- (Gesamtladung -4); unter U fehlt dieser Term.")

# ---- Reflexionspositivität (OS-Form wie Z19, Paar aus zwei Weyl mit Majoranamassenmatrix)
say("\n   Reflexionspositivität (Link, OS-Form wie Z19, 4 Zeitschichten, 6 Raumimpulse):")
PS = [np.array(p) for p in ((0.3, 0.0, 0.0), (1.0, 0.5, 0.2), (2.5, 2.0, 0.3), (np.pi, 0.4, 0.0),
                            (2.8, 2.9, 3.0), (0.1, 0.05, 0.02))]


def D_pair(z, ps, Phi, mm, r=1.0):
    z = np.atleast_1d(np.asarray(z, complex))
    s4, c4 = (z - 1 / z) / (2j), (z + 1 / z) / 2
    ss, Wsp = np.sin(ps), np.sum(1 - np.cos(ps))
    kin = 1j * s4[:, None, None] * G4 + sum(1j * ss[i] * gam[i] for i in range(3))[None]
    n = len(z)
    if Phi is None:                                        # einzelnes Wilson-Dirac (Z19)
        return kin + (r * (Wsp + 1 - c4))[:, None, None] * I4[None]
    W = r * (Wsp + 1 - c4)
    M = np.zeros((n, 2, 2), complex)
    M[:, 0, 0], M[:, 1, 1], M[:, 0, 1], M[:, 1, 0] = Phi, mm, W, W
    Mb = M.copy()
    Mb[:, 0, 0], Mb[:, 1, 1] = np.conj(Phi), np.conj(mm)
    return (np.einsum("ij,nab->niajb", np.eye(2), kin) + np.einsum("nij,ab->niajb", M, P_p)
            + np.einsum("nij,ab->niajb", Mb, P_m)).reshape(n, 8, 8)


def os_test(Dfun, bs, nb=4, sigma=-1, N=4096):
    taus = np.arange(-(2 * nb + 1), 2 * nb + 2)
    p4 = 2 * np.pi * (np.arange(N) + 0.5) / N
    Si = np.linalg.inv(Dfun(np.exp(1j * p4)))
    S = dict(zip(taus, np.einsum("tn,nab->tab", np.exp(1j * np.outer(taus, p4)), Si) / N))
    G = np.kron(np.eye(bs // 4), G4)
    idx = range(1, 1 + nb)
    out = []
    for Am in (np.block([[sigma * G @ S[1 - s - t] for t in idx] for s in idx]),
               np.block([[-sigma * (S[s + t - 1] @ G).T for t in idx] for s in idx])):
        aher = np.max(np.abs(Am - Am.conj().T)) / np.max(np.abs(Am))
        Hm = (Am + Am.conj().T) / 2
        nrm = np.array([np.linalg.norm(Hm[bs * i:bs * i + bs, bs * i:bs * i + bs]) for i in range(nb)])
        d = np.repeat(np.sqrt(np.maximum(nrm, 1e-14 * nrm.max())), bs)
        e_ = np.linalg.eigvalsh(Hm / np.outer(d, d))
        out.append((aher, float(e_.min() / e_.max())))
    return max(a for a, _ in out), min(b for _, b in out)


rp = {}
for r in (1.0, 0.5, 2.0):
    rr_ = [os_test(lambda zz, p=p: D_pair(zz, p, None, None, r=r), 4) for p in PS]
    rp[f"Wilson einzeln r={r}"] = [min(b for _, b in rr_)]
TH = (0.0, np.pi / 2, np.pi, 2.0)
for lab, mmf in (("U", lambda t: 0.0), ("U y'=0,6", lambda t: 0.36 * np.exp(-1j * t)), ("G", lambda t: 0.6)):
    for r in ((1.0, 0.5, 2.0) if lab != "G" else (1.0, 2.0)):
        vals, ah = [], 0.0
        for t in TH:
            rr_ = [os_test(lambda zz, p=p: D_pair(zz, p, 0.6 * np.exp(1j * t), mmf(t), r=r), 8) for p in PS]
            vals.append(min(b for _, b in rr_))
            ah = max(ah, max(a for a, _ in rr_))
        rp[f"{lab} r={r}"] = vals
        rp[f"{lab} r={r} antiherm"] = [ah]
for k_, v_ in rp.items():
    if "antiherm" in k_:
        continue
    say(f"      {k_:22s} OS-Minimum über Phasen 0, pi/2, pi, 2: " + " ".join(f"{x:+.2e}" for x in v_))
TOL = -1e-9
ok_single = rp["Wilson einzeln r=1.0"][0] > TOL and rp["Wilson einzeln r=0.5"][0] > TOL and rp["Wilson einzeln r=2.0"][0] < -0.05
ok_U = all(x > TOL for k_ in ("U r=1.0", "U r=0.5", "U y'=0,6 r=1.0", "U y'=0,6 r=0.5") for x in rp[k_])
phase_indep = np.ptp(rp["U r=2.0"]) < 1e-9 and np.ptp(rp["U y'=0,6 r=2.0"]) < 1e-9 and np.ptp(rp["G r=2.0"]) > 0.1
herm = all(rp[k_][0] < 1e-12 for k_ in rp if "antiherm" in k_)
check("Z19 reproduziert: Wilson r <= 1 link-positiv, r = 2 nicht", ok_single)
check("U-Paarung mit Wilson r <= 1: link-positiv für alle Phasen (OS-Form hermitesch)", ok_U and herm)
check("OS-Spektrum unter U phasenunabhängig (auch r = 2), unter G nicht", phase_indep)
befund("G mit Wilson r <= 1 ist ebenfalls positiv: Positivität trennt U nicht von G, U entfernt die Barriere",
       all(x > TOL for x in rp["G r=1.0"]))
say("   Argument: Das Wilson-Glied ist ein Dirac-Wilson-Term für Psi = (xi+, xi-^c) mit Fermionzahl U(1)_J;")
say("   Phi steht nur auf einem Platz und ist unter Theta symmetrisch (Osterwalder-Seiler). U(1)_J ist")
say("   vektorartig und vertauscht mit Theta, deshalb hängt das OS-Spektrum nicht von der Phase ab.")
ZAHLEN["os_minimum"] = {k_: [float(x) for x in v_] for k_, v_ in rp.items()}

# =====================================================================================================
say("\nI. Index der Naht mit lückenlosem xi-\n")


def low_states(L, kz, kw, window, k=40, sigma=1e-4):
    H, rad = build(L, kz, **kw)
    if L <= 16:
        E, S = eigh(H.toarray(), subset_by_value=(-window, window), driver="evr")
    else:
        E, S = eigsh(H, k=k, sigma=sigma, which="LM")
        sel = np.abs(E) < window
        E, S = E[sel], S[:, sel]
    return E, S, rad, H


def core_branches(L, kw, R=3.5, window=0.15, dk=1e-3):
    E, S, rad, H = low_states(L, 0.0, kw, window)
    S = orth(S)
    core, minus, _ = masks(L, rad, R)
    cval, cvec = np.linalg.eigh(S.conj().T @ (core[:, None] * S))
    Hp, _ = build(L, dk, **kw)
    Hm, _ = build(L, -dk, **kw)
    dH = (Hp - Hm) / (2 * dk)
    cs, es = S @ cvec[:, cval > 0.5], S @ cvec[:, cval < 0.3]
    f = lambda X: np.real(np.einsum("ij,ij->j", X.conj(), dH @ X)) if X.shape[1] else np.array([])
    vc, ve = f(cs), f(es)
    return dict(net=int(np.sum(vc > 0.2) - np.sum(vc < -0.2)), net_rand=int(np.sum(ve > 0.2) - np.sum(ve < -0.2)),
                kern=[round(float(x), 3) for x in cval[cval > 0.5]], v=[round(float(x), 3) for x in vc],
                spalt=float(np.min(np.abs(E))) if len(E) else None)


cfg = {
    "G ausgerichtet": dict(n_plus=1, n_minus=0, flux=0.25, m_minus=0.6),
    "G gegen": dict(n_plus=0, n_minus=-1, flux=0.25, m_minus=0.6),
    "U0 ausgerichtet": dict(n_plus=1, n_minus=0, flux=0.25, m_minus=0.0),
    "U0 gegen": dict(n_plus=0, n_minus=0, flux=0.25, m_minus=0.0),
    "Uy ausgerichtet": dict(n_plus=1, n_minus=-1, flux=0.25, m_minus=0.6, conj_minus=True),
    "Uy gegen": dict(n_plus=0, n_minus=0, flux=0.25, m_minus=0.6, conj_minus=True),
}
counts = {}
say("   Netto chirale Kernzweige (Z7-Methode, Kernprojektor > 0,5, |v| > 0,2), L = 14 und 22:")
for k_, kw in cfg.items():
    r14, r22 = core_branches(14, kw), core_branches(22, kw)
    counts[k_] = (r14["net"], r22["net"])
    say(f"      {k_:16s} netto {r14['net']:+d}/{r22['net']:+d} (vorhergesagt {pred[k_][2]:+d}); Kern v = {r14['v']}; "
        f"Rand {r14['net_rand']:+d}")
check("Zählung wie vorhergesagt: G +1/-1, U0 +1/0, U mit y' 0/0",
      all(counts[k_] == (pred[k_][2], pred[k_][2]) for k_ in cfg))
befund("Gegenorientierung (nu = -1/2) unter U ohne Nahtmode, bei gleicher Gradientenenergie", counts["U0 gegen"] == (0, 0))
befund("mit y' != 0 ist die Naht nicht chiral (netto 0)", counts["Uy ausgerichtet"] == (0, 0))
ZAHLEN["zaehlung_L14_L22"] = {k_: list(v_) for k_, v_ in counts.items()}

say("\n   U0 ausgerichtet: Größe und Masse")
rob, spl = [], {}
for L in (14, 22, 30, 38):
    r_ = core_branches(L, cfg["U0 ausgerichtet"])
    rG = core_branches(L, cfg["G ausgerichtet"]) if L <= 22 else None
    rob.append(r_["net"])
    spl[L] = (r_["spalt"], rG["spalt"] if rG else None)
    say(f"      L = {L}: netto {r_['net']:+d}, Kernprojektor {r_['kern']}, v = {r_['v']}; Kern-Rand-Aufspaltung bei k_z = 0: "
        f"U0 {r_['spalt']:.4f}" + (f", G {rG['spalt']:.1e}" if rG else ""))
for m_ in (0.4, 0.9):
    kw = dict(cfg["U0 ausgerichtet"], m=m_)
    r_ = core_branches(22, kw)
    rob.append(r_["net"])
    say(f"      L = 22, m = {m_}: netto {r_['net']:+d}, Kernprojektor {r_['kern']}")
check("U0 ausgerichtet: netto eine chirale Mode für L = 14 ... 38 und m = 0,4/0,6/0,9", rob == [1] * 6)
pw_exp = -np.log(spl[38][0] / spl[14][0]) / np.log(38 / 14)
say(f"   Die Aufspaltung fällt unter U0 wie L^-{pw_exp:.2f} (Potenz), unter G exponentiell ({spl[14][1]:.1e} -> {spl[22][1]:.1e}).")
befund("Kern-Rand-Kopplung über das lückenlose xi- fällt nur als Potenz", pw_exp < 1.5 and spl[22][1] < 0.1 * spl[14][1])
ZAHLEN["aufspaltung_kz0"] = {str(k_): v_ for k_, v_ in spl.items()}

say("\n   Kontrollmassen für eta- (Grenzwert m_minus -> 0):")
ctrl = {}
for lab, conj in (("G-artig eps (bricht U(1)_J)", False), ("U-artig y' Phi* (U(1)_J-erlaubt)", True)):
    for e_ in (0.6, 0.3, 0.1, 0.03):
        kw = dict(n_plus=1, n_minus=(-1 if conj else 0), flux=0.25, m_minus=e_, conj_minus=conj)
        nets = [core_branches(L, kw)["net"] for L in (14, 30)]
        ctrl[(lab, e_)] = nets
    say(f"      {lab:34s} m_minus = 0,6/0,3/0,1/0,03: netto (L = 14, 30) "
        + "  ".join(f"{ctrl[(lab, e_)]}" for e_ in (0.6, 0.3, 0.1, 0.03)))
gl = ctrl
check("G-artige Kontrollmasse: netto 1 bis eps = 0,03; U-artiges y' = 0,6/0,3: netto 0",
      all(gl[("G-artig eps (bricht U(1)_J)", e_)] == [1, 1] for e_ in (0.6, 0.3, 0.1, 0.03))
      and all(gl[("U-artig y' Phi* (U(1)_J-erlaubt)", e_)] == [0, 0] for e_ in (0.6, 0.3)))
say("   Bei kleinem y' liegt die Gegenmode (xi-, v = -1) im Abstand ~ 1/(y' m) außerhalb des Kerns. Die Naht ist")
say("   nur bis zu diesem Abstand chiral. Der Grenzwert m_minus -> 0 hängt von der Richtung ab: G-artig 1, U-artig 0.")
say("   Am physischen Knoten ist det M = Delta+ Delta- - W(0)^2 = 0 für m_minus = 0: Der Index ist dort nicht")
say("   topologisch definiert; er gilt nur in der Klasse y' = 0.")
befund("Index unter U nicht topologisch geschützt (det M am Knoten = 0; U(1)_J-erlaubtes y' -> 0+ gibt 0)",
       gl[("U-artig y' Phi* (U(1)_J-erlaubt)", 0.3)] == [0, 0])

say("\n   Hybridisierung mit dem xi--Kontinuum (k_z = 0,3; Schwelle E_thr = sqrt(sin^2 k_z + s^2)):")
kz_ = 0.3
Wk = 1 - np.cos(kz_)
smin = (np.sqrt(m0 ** 2 + 4 * Wk ** 2) - m0) / 2
Ethr = np.sqrt(np.sin(kz_) ** 2 + smin ** 2)
say(f"      sin k_z = {np.sin(kz_):.5f}, E_thr = {Ethr:.5f} (Abstand {Ethr - np.sin(kz_):.1e})")
hyb = {}
for lab in ("U0 ausgerichtet", "G ausgerichtet"):
    for L in (14, 22, 30, 38):
        H, rad = build(L, kz_, **cfg[lab])
        E, Vv = eigsh(H, k=40, sigma=np.sin(kz_), which="LM")
        core, minus, _ = masks(L, rad)
        wc = np.real(np.einsum("ij,ij->j", Vv.conj(), core[:, None] * Vv))
        i_ = int(np.argmax(wc))
        rr = np.tile(np.repeat(rad, 4), 2)
        far = float(np.real(np.sum((rr > 6) * minus * np.abs(Vv[:, i_]) ** 2)))
        hyb[(lab, L)] = (float(E[i_]), float(wc[i_]), float(np.real(np.sum(minus * np.abs(Vv[:, i_]) ** 2))), far)
    say(f"      {lab}: " + "; ".join(f"L={L}: E={hyb[(lab, L)][0]:.5f}, Kern {hyb[(lab, L)][1]:.3f}, "
                                   f"xi- fern {hyb[(lab, L)][3]:.3f}" for L in (14, 22, 30, 38)))
dU = hyb[("U0 ausgerichtet", 14)][1] - hyb[("U0 ausgerichtet", 38)][1]
dG = abs(hyb[("G ausgerichtet", 14)][1] - hyb[("G ausgerichtet", 38)][1])
say(f"   Unter U0 sitzt die Nahtmode an der Kontinuumsschwelle; ihr Kerngewicht fällt um {dU:.3f} von L = 14 bis 38,")
say(f"   ihr fernes xi--Gewicht wächst. Unter G bleibt es konstant (Änderung {dG:.1e}).")
befund("U0: Nahtmode hybridisiert mit dem lückenlosen xi- (Kerngewicht fällt mit L, xi--Schwanz wächst)",
       dU > 0.08 and dG < 2e-3 and hyb[("U0 ausgerichtet", 38)][3] > 5 * hyb[("U0 ausgerichtet", 14)][3])
ZAHLEN["hybridisierung_kz0.3"] = {f"{a}|L={b}": v_ for (a, b), v_ in hyb.items()}

say("\n   Orientierung: Fermionische Vakuumenergie je Längeneinheit (k_z-Mittel über 8 Punkte)")


def egs(L, kw, nk=8, R=3.5):
    ks = 2 * np.pi * (np.arange(nk) + 0.5) / nk - np.pi
    tot = core_e = 0.0
    for kz in ks:
        H, rad = build(L, kz, **kw)
        E, Vv = np.linalg.eigh(H.toarray())
        pos = E > 0
        cm, _, _ = masks(L, rad, R)
        tot += -0.5 * E[pos].sum()
        core_e += -0.5 * (E[pos] * (np.abs(Vv[:, pos]) ** 2 * cm[:, None]).sum(0)).sum()
    return tot / nk, core_e / nk


dE = {}
for (L, m_) in ((12, 0.6), (14, 0.6), (12, 0.9)):
    a_ = egs(L, dict(cfg["U0 ausgerichtet"], m=m_))
    b_ = egs(L, dict(cfg["U0 gegen"], m=m_))
    dE[(L, m_)] = (a_[0] - b_[0], a_[1] - b_[1])
    say(f"      L = {L}, m = {m_}: E(ausgerichtet) - E(gegen) = {dE[(L, m_)][0]:+.5f} gesamt, {dE[(L, m_)][1]:+.5f} im Kern")
pref = all(v_[0] > 0 and v_[1] > 0 for v_ in dE.values())
say("   Beide Texturen haben dieselbe Gradientenenergie von Phi (nu^2 = 1/4). Die Fermionen heben die Entartung:")
say(f"   {'die Orientierung ohne Nahtmode ist billiger' if pref else 'die Orientierung mit Nahtmode ist billiger'} "
    f"(je Paar und Gitterlänge {dE[(14, 0.6)][0]:.4f} bei m = 0,6 und {dE[(12, 0.9)][0]:.4f} bei m = 0,9; 16 Paare: x 16).")
say("   Der Unterschied sitzt im Kern und hängt nicht von L ab.")
befund("Fermionenergie bevorzugt die Orientierung ohne Nahtmode", pref)
ZAHLEN["orientierung_dE"] = {f"L={a},m={b}": v_ for (a, b), v_ in dE.items()}

# =====================================================================================================
say("\n1. Eine Naht\n")
hol = {nm: np.exp(1j * np.pi / 2 * qq) for nm, qq in (("xi+", 1), ("xi-", -1), ("Phi", -2))}
lab4 = lambda x: {1: "1", 1j: "i", -1: "-1", -1j: "-i"}[complex(np.round(x.real, 9) + 0.0, np.round(x.imag, 9) + 0.0)]
say("   Deck-Holonomie als U(1)_J-Element exp(i pi Q_J / 2) (Viertelfluss): " +
    ", ".join(f"{nm} -> {lab4(v_)}" for nm, v_ in hol.items()))
Jn = np.array(Jm.tolist(), dtype=complex)
Bpp_n = np.array(Bpp.tolist(), dtype=complex)
check("Holonomie: xi+ -> i, xi- -> -i, Phi -> -1 (erzwungene Nullstelle); stimmt mit J^T B J = -B (Z20)",
      np.isclose(hol["xi+"], 1j) and np.isclose(hol["xi-"], -1j) and np.isclose(hol["Phi"], -1)
      and np.allclose(Jn.T @ Bpp_n @ Jn, -Bpp_n))
# Minimum über Z4-wertige Schläuche: Holonomien i^k, Produkt i; Phi-Kosten nu(k)^2 mit nu in Z + k/2
nu_min = {1: 0.5, 2: 0.0, 3: 0.5}
best = None
for n_t in range(1, 5):
    for ks_ in itertools.combinations_with_replacement((1, 2, 3), n_t):
        if sum(ks_) % 4 != 1:
            continue
        cost = sum(nu_min[k__] ** 2 for k__ in ks_) + 1e-3 * n_t     # Kernkosten je Schlauch
        if best is None or cost < best[0]:
            best = (cost, ks_)
odd_always = all(sum(1 for k__ in ks_ if k__ % 2) % 2 == 1 for n_t in range(1, 5)
                 for ks_ in itertools.combinations_with_replacement((1, 2, 3), n_t) if sum(ks_) % 4 == 1)
say(f"   Z4-wertige Schläuche mit Gesamtholonomie i (bis 4 Schläuche): Minimum {best[1]} mit Kosten {best[0] - 1e-3:.2f};")
say("   Schläuche mit Holonomie -1 erzwingen keine Nullstelle (nu = 0). Die Zahl der Schläuche mit Nullstelle ist")
say("   immer ungerade.")
check("Z4-Holonomien: minimal genau ein Schlauch (Holonomie i, nu = 1/2); Zahl der Nullstellen-Schläuche ungerade",
      best[1] == (1,) and odd_always)
split = [1 / (4 * n) for n in (1, 2, 4, 8)]
say(f"   Wäre die Holonomie kontinuierlich in U(1)_J: n Bruchschläuche mit je pi/(2n) kosten sum nu^2 = {split}.")
say("   Die Kosten fallen; ohne Z4-Quantisierung zerfließt der Schlauch. U(1)_J selbst quantisiert nicht; das tut")
say("   nur die Deck-Struktur (Strukturgruppe Z4, Eingabe aus Z20).")
befund("kontinuierliche U(1)_J quantisiert die Holonomie nicht (Teilung senkt die Kosten)", split[-1] < split[0])
say("   Geeichte U(1)_J: Mit einem Kondensat der Ladung 2 ist das Flussquant pi (Fermionholonomie -1). Ein")
say("   Viertelfluss wäre nicht stabil (Meißner). U muss als globale Symmetrie gelesen werden (Argument).")
befund("geeichte U(1)_J trüge keinen stabilen Viertelfluss (Argument, nicht gerechnet)", True)
# Z11-Parität: geschlossene Flächen in S^3 schneiden einen geschlossenen Schlauch gerade oft
rng = np.random.default_rng(22)
t = np.linspace(0, 2 * np.pi, 20001)
fib = np.stack([np.cos(t), np.sin(t), 0 * t, 0 * t], 1)            # Polfaser z1 = 0
cross = []
for _ in range(300):
    Qf = rng.normal(size=(4, 4))
    lin = rng.normal(size=4)
    c0 = rng.normal() * 0.5
    f = np.einsum("ni,ij,nj->n", fib, Qf, fib) + fib @ lin + c0          # Nullmenge = geschlossene Fläche in S^3
    cross.append(int(np.sum(np.sign(f[1:]) != np.sign(f[:-1]))))
say(f"   Z11-Parität: Die Polfaser schneidet 300 zufällige geschlossene Flächen in S^3 {sorted(set(cross))} mal: stets")
say("   gerade, weil H_2(S^3) = 0 (jede Fläche berandet). Die Paarregel von Z11 ist für einen geschlossenen Schlauch")
say("   in S^3 also von selbst erfüllt. Ungerade wäre nur der Durchstoß durch die Basissphäre; dort lebt Phi im")
say("   internen Zweig aber nicht.")
ZAHLEN["z11_schnittzahlen"] = sorted(set(cross))

# =====================================================================================================
say("\n2. Ein Block: legt U die Zahl N der 16er-Sätze fest?\n")
Nn = sp.symbols("N", positive=True, integer=True)
weyl = [(1, 16 * Nn), (-1, 16 * Nn)]                     # (U(1)_J-Ladung, Anzahl Weyl)
A3 = sp.simplify(sum(q ** 3 * n for q, n in weyl))
Ag = sp.simplify(sum(q * n for q, n in weyl))
TR = sp.symbols("T_R")
Amix = sp.simplify(sum(q * n / (16 * Nn) * TR * Nn for q, n in weyl))
n1 = 16 * Nn
n3 = 16 * Nn
nu_uv = sp.simplify(n1 - n3)
nu_ir = sp.Mod(-16 * Nn, 16)
say(f"   U(1)_J^3 = {A3}, U(1)_J-Gravitation = {Ag}, U(1)_J-G-G = {Amix} (Paare vektorartig: Psi = (xi+, xi-^c)).")
say(f"   Spin-Z4: xi+ Ladung 1, xi- Ladung 3: nu = n1 - n3 = {nu_uv}. Leichter Sektor allein (xi-): nu = -16N mod 16 = "
    f"{sp.simplify(nu_ir)}. Phi (Ladung 2) bricht Z4 -> Z2, also ist im Infraroten nichts anzupassen.")
glue = {n: (n / 4) % 2 == 0 for n in (8, 16, 24, 32)}
say(f"   Nahtgitter D_8N^+: Klebvektor (1/2)^(8N) hat Norm 2N, gerade für N = 1..4: {glue}. Jedes N gibt ein gerades")
say("   unimodulares Gitter (N = 1: E8; N = 2: D16^+, N = 3: D24^+ wie Z17). c_- = 8N ist stets in 8Z.")
befund("U(1)_J-Anomalien, Spin-Z4 und Klebung legen N nicht fest", A3 == 0 and Ag == 0 and Amix == 0 and nu_uv == 0
       and all(glue.values()))
befund("E8 auf der Naht bleibt eine eigene Forderung (N = 1)", True)
say("   Unter U koppeln alle Paare an Phi; die Zahl der windenden 16er-Sätze ist die Zahl der Sätze im Multiplett.")

# =====================================================================================================
say("\nD. Darstellungen: Phi xi+_a xi+_b braucht eine symmetrische invariante Form auf dem Index a\n")
px, py, pz = sx, sy, sz
I2 = np.eye(2)


def kronl(ms):
    out = ms[0]
    for mm_ in ms[1:]:
        out = np.kron(out, mm_)
    return out


Gm = []
for k_ in range(5):
    Gm.append(kronl([pz] * k_ + [px] + [I2] * (4 - k_)))
    Gm.append(kronl([pz] * k_ + [py] + [I2] * (4 - k_)))
G11 = kronl([pz] * 5)
cl = max(np.abs(Gm[i] @ Gm[j] + Gm[j] @ Gm[i] - 2 * (i == j) * np.eye(32)).max() for i in range(10) for j in range(10))
ip_ = np.where(np.diag(G11).real > 0)[0]
T16 = {(i, j): (0.5 * Gm[i] @ Gm[j])[np.ix_(ip_, ip_)] for i in range(10) for j in range(i + 1, 10)}
span45 = np.linalg.matrix_rank(np.array([v_.ravel() for v_ in T16.values()]))


def inv_space(gl_, gr_, sym=None, tol=1e-9):
    n_, m_ = gl_[0].shape[0], gr_[0].shape[0]
    rows = [np.kron(Aa.T, np.eye(m_)) + np.kron(np.eye(n_), Bb.T) for Aa, Bb in zip(gl_, gr_)]
    if sym is not None:
        P = np.zeros((n_ * m_, n_ * m_))
        for i in range(n_):
            for j in range(m_):
                P[i * m_ + j, j * n_ + i] = 1
        rows.append(np.eye(n_ * m_) - sym * P)
    _, s_, vh = np.linalg.svd(np.vstack(rows))
    return int(np.sum(s_ < tol) + (vh.shape[0] - len(s_))), vh[np.sum(s_ > tol):].conj()


adj = [T16[(i, i + 1)] for i in range(9)]
d_sym, _ = inv_space(adj, adj, sym=+1)
d_asym, _ = inv_space(adj, adj, sym=-1)
d_bar, _ = inv_space(adj, [a.conj() for a in adj])
hom = np.linalg.svd(np.vstack([np.kron(a.conj(), np.eye(16)) - np.kron(np.eye(16), a.T) for a in adj]), compute_uv=False)
d_hom = int(np.sum(hom < 1e-9))
so16 = []
for i in range(15):
    Aa = np.zeros((16, 16), complex)
    Aa[i, i + 1], Aa[i + 1, i] = 1, -1
    so16.append(Aa)
v_sym, _ = inv_space(so16, so16, sym=+1)
v_asym, _ = inv_space(so16, so16, sym=-1)
Xg = sum(T16[(2 * k_, 2 * k_ + 1)] for k_ in range(5))
keys = list(T16)
_, s_, vh = np.linalg.svd(np.array([(T16[k_] @ Xg - Xg @ T16[k_]).ravel() for k_ in keys]).T)
nul = vh[np.sum(s_ > 1e-9):].conj()
Kc = [sum(c * T16[k_] for c, k_ in zip(vec, keys)) for vec in nul]
comms = np.array([(Kc[a] @ Kc[b] - Kc[b] @ Kc[a]).ravel() for a in range(len(Kc)) for b in range(a + 1, len(Kc))])
_, s2, vh2 = np.linalg.svd(comms, full_matrices=False)
r5 = int(np.sum(s2 > 1e-9))
su5 = [v_.reshape(16, 16) for v_ in vh2[:r5]]
d5, n5 = inv_space(su5, su5, sym=+1)
rank5 = np.linalg.matrix_rank(n5[0].reshape(16, 16), tol=1e-9) if d5 else 0
du5, _ = inv_space(su5 + [Xg], su5 + [Xg], sym=+1)
say(f"   Spin(10): Clifford-Rest {cl:.0e}, 45 Erzeuger auf der 16 unabhängig: {span45 == 45}.")
say(f"   16 (x) 16: symmetrische Invarianten {d_sym}, antisymmetrische {d_asym}. 16 (x) 16-quer: {d_bar}. Hom(16, 16-quer): {d_hom}.")
say(f"   SO(16)-Vektor: symmetrische Invarianten {v_sym} (delta_ab), antisymmetrische {v_asym}.")
say(f"   u(5) = Kommutant von X in so(10): Dimension {len(Kc)}; su(5) = [u(5), u(5)]: {r5}. SU(5)-invariante")
say(f"   symmetrische Formen auf 16: {d5}, Rang {rank5}; mit U(1)_X: {du5}.")
check("Spin(10)-16: keine Invariante in 16 x 16 (weder Phi xi+xi+ als Singulett noch W xi+xi-); 16 x 16-quer: genau eine",
      cl == 0 and span45 == 45 and d_sym == 0 and d_asym == 0 and d_bar == 1 and d_hom == 0)
check("SO(16)-Vektor: genau eine symmetrische Invariante", v_sym == 1 and v_asym == 0)
check("SU(5)-invariante symmetrische Form auf 16: Dimension 1, Rang 1 (126-quer-Richtung)", len(Kc) == 25 and r5 == 24
      and d5 == 1 and rank5 == 1 and du5 == 0)
say("   Folgen:")
say("   - U braucht psi und chi in derselben Darstellung R (sonst gibt es keine G-äquivariante Drehung; Hom(16, 16-quer) = 0).")
say("   - Ein Singulett-Phi und das Wilson-Glied brauchen eine Invariante in R x R: R reell (SO(16)-Vektor).")
say("   - Dann trägt xi- dieselbe reelle R: Die Familien sind vektorartig, ihre Majoranamasse ist eichinvariant und")
say("     nur durch U(1)_J (also durch y' = 0) verboten. xi- kann nicht die chirale SM-Generation (16 von Spin(10)) sein.")
say("   - Mit R = 16 von Spin(10) müsste Phi in 126-quer liegen; seine SU(5)-Richtung paart 1 von 16 Komponenten:")
say("     Je 16er-Satz windet eine Mode statt 16, also keine c = 8 je Satz. Außerdem fehlte das Wilson-Glied.")
befund("Singulett-Phi verlangt reelle Darstellung; xi- dann vektorartig, nicht die chirale SM-Generation", d_sym == 0 and v_sym == 1)
befund("mit Spin(10)-16 und Phi in 126-quer windet nur 1 von 16 Komponenten", rank5 == 1)

# =====================================================================================================
entscheidungen = {
    "1_eine_naht": ("Lesart (wie Z20, verschärft): Die Holonomie ist ein Viertelfluss von U(1)_J; Phi erhält -1 und muss "
                    "verschwinden; mit Z4-wertigen Holonomien ist ein Schlauch minimal, verträglich mit Z11 (H_2(S^3) = 0). "
                    "Neu: U(1)_J quantisiert die Holonomie nicht (Z4-Deckstruktur bleibt Eingabe, U(1)_J global), und von den "
                    "zwei gleich teuren Halbwindungen trägt nur die ausgerichtete eine Nahtmode; die Fermionenergie bevorzugt "
                    "die andere."),
    "2_ein_block": ("offen: U(1)_J ist vektorartig, alle Anomalien und Spin-Z4 verschwinden für jedes N; D_8N^+ ist für "
                    "jedes N gerade unimodular. 'E8 auf der Naht' bleibt eine eigene Forderung."),
    "3_paarung": ("Lesart, auf eine Zahl reduziert: U schließt psi chi aus (nicht kovariant) und lässt für Phi nur xi+xi+ "
                  "zu; U(1)_J erlaubt aber Phi* xi-xi- (y'), und kein Symmetrieargument trennt y' von y, solange das "
                  "Wilson-Glied erlaubt ist. Bei y' = 0 sind die Familien exakt masselos, bei y' != 0 erhalten sie "
                  "|y'| x {0, 1/2, 1/2} und die Naht wird nicht chiral. Darstellungen: Singulett-Phi verlangt reelles R, "
                  "xi- ist dann vektorartig und nicht die SM-Generation."),
    "4_gitter": ("unter U hergeleitet: Wilson-Glied ist U(1)_J-neutral, das Phasenpotential exakt flach (Hamilton und "
                 "euklidisch), [H, Q_J] = 0, und die Paarung ist für Wilson r <= 1 link-reflexionspositiv für alle "
                 "Phasen. Keine Wand, keine Positivitätsverletzung."),
    "index": ("Lesart: Mit y' = 0 und ausgerichteter Halbwindung zählt der Kern netto eine chirale Majorana-Weyl-Mode je "
              "Paar (L = 14 ... 38, m = 0,4 ... 0,9). Sie ist nicht topologisch geschützt (det M am Knoten = 0): "
              "U(1)_J-erlaubtes y' > 0 gibt netto 0 jenseits ~1/(y' m). Sie hybridisiert mit dem lückenlosen xi- "
              "(Kerngewicht fällt mit L, Kern-Rand-Kopplung als Potenz)."),
}
eingaben = [
    "y' = 0: keine Kopplung Phi* xi- xi- (nicht durch U(1)_J geschützt; Holomorphie oder eine eigene U(1)_- nötig)",
    "Orientierung: Halbwindung von Phi mit der Deck-Holonomie ausgerichtet (Fermionenergie bevorzugt die Gegenrichtung)",
    "Deck-Holonomie Z4-wertig und physisch, U(1)_J global (aus Z20; U(1)_J quantisiert nicht)",
    "keine erzwungene Reziprozität der Konfiguration (Z20)",
    "reelle Darstellung R der Paare (SO(16)-Vektor); SM-Generation nicht aus xi-",
    "N = 1 ('E8 auf der Naht')",
]
urteil = ("Prinzip U löst Entscheidung 4 (flaches Potential, Wilson positiv, keine Wand) und schließt psi chi aus. Es löst "
          "Entscheidung 3 nicht ganz: U(1)_J erlaubt Phi* xi-xi-, und erst y' = 0 macht die Familien masselos und die Naht "
          "chiral; dann ist der Index eine Lesart ohne topologischen Schutz, mit Hybridisierung ins lückenlose xi-. "
          "Entscheidung 1 bleibt Lesart mit einer neuen Orientierungswahl, die die Fermionenergie gegen die Naht entscheidet. "
          "Entscheidung 2 bleibt offen. Die Darstellungen verlangen ein reelles R: xi- ist dann nicht die SM-Generation.")
say("\nUrteil je Entscheidung")
for k_, v_ in entscheidungen.items():
    say(f"   {k_}: {v_}")
say("\nVerbleibende Eingaben")
for e_ in eingaben:
    say(f"   - {e_}")
say(f"\nUrteil: {urteil}")
say(f"\n(Rechenzeit {time.time() - T0:.0f} s)")
say("\nZusammenfassung")
for n_, c_ in CHECKS.items():
    say(f"   {'bestanden' if c_ else 'NICHT bestanden'}: {n_}")
say("Befunde")
for n_, c_ in BEFUNDE.items():
    say(f"   {'ja' if c_ else 'nein'}: {n_}")


def _clean(o):
    if isinstance(o, dict):
        return {str(k): _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    return o


(HERE / "z22_ergebnis.json").write_text(json.dumps(_clean({
    "checks": CHECKS, "befunde": BEFUNDE, "urteil": urteil, "entscheidungen": entscheidungen,
    "eingaben": eingaben, "zahlen": ZAHLEN}), ensure_ascii=False, indent=1))
assert all(CHECKS.values())
