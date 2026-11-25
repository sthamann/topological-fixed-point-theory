"""Z136: Ist P_G ein freier Parameter, oder ist G die Einheit der einen nativen Länge und G/R_Vak^2 durch
Verschränkungsgleichgewicht nativ festgelegt?

Leitlinie 'Rekonstruktion statt Auswahl'. Firewalls: GRAV.NONCIRCULAR.01 (Bisognano-Wichmann / Hislop-Longo /
Casini-Huerta-Myers nur als RESULT, sonst als ASSUMPTION ausgewiesen; nie 'KMS => Lorentz => BW'),
DIMENSION.UPLIFT.FIREWALL.01 (Aussagen über die Naht bleiben 1D-typisiert; su(2,2) ist Kinematik einer Darstellung).
Gelöste Blöcke (Flavour, Massen, alpha, 16, E8) und H6 bleiben unberührt; die Ledgerwerte GRAV.ENTROPY.EQUILIBRIUM.01
(v358), GRAV.NONLINEAR.01 (v359), COSMO.LAM.01 / LAMBDA.BRANCH.01 (v7, v60, v274; Status 'Numerical', hier nur
Vergleichsziel) und der No-Unit-Satz (v153, v274, v364) werden nur gelesen; nichts wird befördert, keine vollständige
Theorie wird behauptet. Bestand: z72 (Wiesbrock, Peschel am Gitter), z83 (radiale su(2,2), E = L_0, beta = 4 pi, beta' =
8 pi), z86 (S = ln 2/2 je Naht, erster Hauptsatz je Sektor, keine native Fläche), z90 (Flächenkoeffizient reglerabhängig),
z98 (MM-Algebra nativ, sigma = +1, l = R_Vak/2), z104 (Kappe O_1, K_cap = 2 pi Q^-, Wick-Abbildung), z112 (L = kappa/(2 g^2)),
z118 (64 eps je Index, Gewicht 2 c_3), kernel-followups-20261120/check_6_gravity.py (Flächenzustand, G = l^2/(18 zeta(3) c_3)).
Statuskürzel: E exakt, C bedingt, L Lesart, O offen.

Teil 1  Nativer geometrischer Modulfluss (die Firewall): radiale sl(2,R) {E, Q^-, D} in su(2,2) auf dem Twistorraum,
        Konjugationsklassen, Doppelkegel (Zentrum, Radius, Bifurkationssphäre, Oberflächengravitation), Konjugation des
        Nahtgenerators K_Uhr = 4 pi E auf den Hislop-Longo-Generator K_cap = 2 pi Q^- (reell: unmöglich; komplex: Wick-
        Element), der Faktor 2 (Perioden auf Vektor, Twistor, Quarter; KMS-Wörterbuch), Peschel-Probe auf der halben Naht.
Teil 2  Erster Hauptsatz delta S = delta<K> im Fockraum der Naht (Kovarianz-, Quartik- und Möbiusstörungen),
        Negativkontrollen (beta = 2 pi, E allein, hyperbolischer Generator), zweite Ordnung (relative Entropie =
        Fisher-Information > 0), bosonische Leiter.
Teil 3  Native Entropie s_0 (Quarter, Naht, E8-Naht, Vielfachheiten, Leiter), Fläche der Bifurkationssphäre je Rahmen,
        Jacobson-Wert G/R_Vak^2, Überbestimmungstabelle aller Normierungen von G.
Teil 4  Einheit oder Parameter: No-Unit-Zählung, physikalische Bedeutung des nativen Werts.
Teil 5  Prämissen der Jacobson-Herleitung: nativ / bedingt / fehlend, mit Belegdatei.
"""
import itertools
import json
from pathlib import Path

import mpmath as mp
import numpy as np
import scipy.linalg as sla
import sympy as sp

HERE = Path(__file__).parent
REPO = HERE.parents[2]
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
rng = np.random.default_rng(136)
mp.mp.dps = 40


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


def prop(A, B):
    """A = lam B? (lam, relatives Residuum)"""
    a, b = np.ravel(A), np.ravel(B)
    lam = np.vdot(b, a)/np.vdot(b, b)
    return complex(lam), float(np.abs(a - lam*b).max()/np.abs(b).max())


def komm(A, B):
    return A @ B - B @ A


def g6(x, n=6):
    return f"{float(x):.{n}g}"


# ------------------------------------------------------------------ Twistordaten (wie z98, z104)
I2 = np.eye(2)
SIGMA = [np.eye(2, dtype=complex), np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]]),
         np.diag([1., -1.]).astype(complex)]
ETA = np.diag([1., -1., -1., -1.])
G4 = np.diag([1., 1., -1., -1.]).astype(complex)                        # CCR-Form auf T
uvec = [np.array([1, 0, -1, 0.]), np.array([0, 1, 0, -1.])]
vvec = [np.array([1, 0, 1, 0.]), np.array([0, 1, 0, 1.])]
NP = [-0.5*sum(s[i, j]*np.outer(uvec[i], vvec[j]) for i in range(2) for j in range(2)) for s in SIGMA]
NK = [-0.5*sum(s[i, j]*np.outer(vvec[i], uvec[j]) for i in range(2) for j in range(2)) for s in SIGMA]
NE = np.diag([-0.5, -0.5, 0.5, 0.5]).astype(complex)
MP = [1j*N for N in NP]
MKL = [ETA[m, m]*1j*NK[m] for m in range(4)]
ME = 1j*NE                                                              # Uhr E = L_0 (elliptisch)
MQ = (MP[0] - MKL[0])/2                                                  # Q^- = (P_0 - K_0)/2 (Kappe, hyperbolisch)
MD = komm(ME, MQ)                                                        # Dilatation
TMS_np = np.array([[1, 0, -1, 0], [0, 1, 0, -1], [1, 0, 1, 0], [0, 1, 0, 1]], dtype=complex)


def punkt(B):
    wc = TMS_np @ B
    if abs(np.linalg.det(wc[2:])) < 1e-9:
        return None
    xx = -1j*wc[:2] @ np.linalg.inv(wc[2:])
    return [complex(np.round(np.trace(xx @ s)/2, 9)) for s in SIGMA]


def eigenebenen(M):
    """Liste (Eigenwert, Gram-Eigenwerte der (2,2)-Form, Punkt x^mu oder None), nach Eigenwert sortiert."""
    w, V = np.linalg.eig(M)
    lams = []
    for lam in w:
        if not any(abs(lam - l2) < 1e-8 for l2 in lams):
            lams.append(lam)
    lams.sort(key=lambda z: (round(z.real, 8), round(z.imag, 8)))
    aus = []
    for lam in lams:
        B = np.linalg.qr(V[:, np.abs(w - lam) < 1e-8])[0]
        gram = B.conj().T @ G4 @ B
        aus.append((complex(np.round(lam, 9)), np.round(np.linalg.eigvalsh((gram + gram.conj().T)/2), 9).tolist(), punkt(B)))
    return aus


def vektor6(M):
    """so(4,2)-Vektordarstellung auf Bivektoren: I -> M I + I M^T (z98)."""
    paare = list(itertools.combinations(range(4), 2))
    Bs = []
    for i, j in paare:
        B = np.zeros((4, 4), complex)
        B[i, j], B[j, i] = 1, -1
        Bs.append(B)
    V = np.zeros((6, 6), complex)
    for k, B in enumerate(Bs):
        Y = M @ B + B @ M.T
        V[:, k] = [Y[i, j] for i, j in paare]
    return V


def verflechter(A, B):
    """Dimension von {g : B g = g A} und kleinster Singulärwert (vec spaltenweise)."""
    L = np.kron(np.eye(4), B) - np.kron(A.T, np.eye(4))
    s = np.linalg.svd(L, compute_uv=False)
    return int(np.sum(s < 1e-9)), float(s.min())


XS = sp.symbols("x0:4", real=True)
SIGS = [sp.eye(2), sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -sp.I], [sp.I, 0]]), sp.Matrix([[1, 0], [0, -1]])]
XM = sum((XS[k]*SIGS[k] for k in range(4)), sp.zeros(2))
TMS = sp.Matrix([[1, 0, -1, 0], [0, 1, 0, -1], [1, 0, 1, 0], [0, 1, 0, 1]])
TMSi = TMS.inv()


def zu_sympy(M):
    return sp.Matrix(4, 4, [sp.nsimplify(complex(x).real) + sp.I*sp.nsimplify(complex(x).imag) for x in np.ravel(M)])


def vektorfeld(Mz):
    """Physikalische Orientierung (P_0 = +d_0), wie z98/z104."""
    Mw = TMS*zu_sympy(Mz)*TMSi
    A, B, C, D = Mw[:2, :2], Mw[:2, 2:], Mw[2:, :2], Mw[2:, 2:]
    dx = A*XM - XM*D - sp.I*B - sp.I*XM*C*XM
    return [sp.expand(-(dx*SIGS[k]).trace()/2) for k in range(4)]


# =================================================================== Teil 1
say("Teil 1: Nativer geometrischer Modulfluss -- die Firewall")
tr2 = {n: float(np.trace(M @ M).real) for n, M in (("E", ME), ("Q", MQ), ("D", MD), ("P_0", MP[0]))}
in_su22 = all(np.abs(M.conj().T @ G4 + G4 @ M).max() < 1e-14 and abs(np.trace(M)) < 1e-14 for M in (ME, MQ, MD, MP[0]))
sl2 = (np.abs(komm(ME, MQ) - MD).max() < 1e-14 and np.abs(komm(MD, ME) - MQ).max() < 1e-14
       and np.abs(komm(MQ, MD) + ME).max() < 1e-14 and np.abs(ME - (MP[0] + MKL[0])/2).max() < 1e-14)
ebE, ebQ, ebD = eigenebenen(ME), eigenebenen(MQ), eigenebenen(MD)
spitzen = sorted(round(p[2][0].real, 9) for p in ebQ)
xiQ, xiE = vektorfeld(MQ), vektorfeld(ME)
t_s, r_s = sp.symbols("t r", real=True)
rr2 = XS[1]**2 + XS[2]**2 + XS[3]**2
xiQ_soll = [(1 - XS[0]**2 - rr2)/2] + [-XS[0]*XS[k] for k in (1, 2, 3)]
feld_ok = all(sp.expand(a - b) == 0 for a, b in zip(xiQ, xiQ_soll))
norm = lambda xi: sp.expand(xi[0]**2 - xi[1]**2 - xi[2]**2 - xi[3]**2)
nQ, nE = norm(xiQ), norm(xiE)
achse = {XS[2]: 0, XS[3]: 0}
nQ_r = sp.factor(nQ.subs(achse).subs({XS[0]: t_s, XS[1]: r_s}))
# Oberflächengravitation auf dem Zukunftshorizont t = 1 - r (Schnitt y = z = 0): d_mu (xi^2) = -2 kappa xi_mu
hor = {XS[0]: 1 - XS[1]}
xi_unten = [xiQ[0], -xiQ[1], -xiQ[2], -xiQ[3]]
kappa_werte = []
for m in (0, 1):
    grad = sp.diff(nQ, XS[m]).subs(achse).subs(hor)
    xi_m = xi_unten[m].subs(achse).subs(hor)
    kappa_werte.append(sp.simplify(-grad/(2*xi_m)))
bif = [sp.simplify(c.subs({XS[0]: 0, XS[1]: 1, XS[2]: 0, XS[3]: 0})) for c in xiQ]
w_chm = sp.simplify(xiQ[0].subs({XS[0]: 0, XS[1]: r_s, XS[2]: 0, XS[3]: 0}))
chm_int = sp.integrate(w_chm*4*sp.pi*r_s**2, (r_s, 0, 1))
nE_min = sp.factor(nE.subs(achse).subs({XS[0]: t_s, XS[1]: r_s}))
VE, VQ = vektor6(ME), vektor6(MQ)
evVE = sorted(np.round(np.linalg.eigvals(VE).imag, 8).tolist())
evVQ = sorted(np.round(np.linalg.eigvals(VQ).real, 8).tolist())
ZAHLEN["radial"] = {"tr M^2": tr2, "Eigenebenen E (Eigenwert, Gram, Punkt)": str(ebE), "Eigenebenen Q^-": str(ebQ),
                    "Spitzen": spitzen, "xi_Q": [str(c) for c in xiQ], "|xi_Q|^2 (t, r)": str(nQ_r),
                    "kappa (t-, r-Komponente)": [str(k) for k in kappa_werte], "CHM-Gewicht w(r)": str(w_chm),
                    "int_B w d^3x": str(chm_int), "Vektordarstellung E/i": evVE, "Vektordarstellung Q^-": evVQ}
check(f"[E] Radiale sl(2,R) in su(2,2) (Twistorraum, explizite 4x4-Matrizen): E, Q^-, D, P_0 sind G-antihermitesch und "
      f"spurfrei ({in_su22}), [E, Q^-] = D, [D, E] = Q^-, [Q^-, D] = -E, E = (P_0 + K_0)/2 ({sl2}); tr M^2 = {tr2}: E "
      f"elliptisch (+-i/2), Q^- und D hyperbolisch (+-1/2), P_0 parabolisch. Der Doppelkegel von Q^-: Null-Eigenebenen "
      f"(Gram {ebQ[0][1]}) mit den reellen Punkten x = ({spitzen[0]:+.0f}, 0, 0, 0), ({spitzen[1]:+.0f}, 0, 0, 0), also "
      f"Zentrum x = 0 (tau-Fixpunkt), Radius 1 = Im x_Vak = R_Vak (Vakuumpunkte x^0 = +-i aus den definiten Eigenebenen "
      f"von E: {[str(p[2][0]) for p in ebE]}). Feld xi_Q = ((1 - t^2 - r^2)/2, -t x) ({feld_ok}), |xi_Q|^2 = {nQ_r} "
      f"verschwindet genau auf dem Rand von O_1 = {{|t| + r < 1}}, xi_Q = 0 auf der Bifurkationssphäre t = 0, r = 1 "
      f"({bif}), Oberflächengravitation kappa = {kappa_werte[0]} (beide Komponenten gleich: {kappa_werte[0] == kappa_werte[1]}): "
      f"K = 2 pi Q^- ist genau die Hislop-Longo-/CHM-Normierung, mit dem Gewicht w(r) = {w_chm} und int_B w d^3x = "
      f"{chm_int} (v358: 4 pi R^4/15). Kontrolle: |xi_E|^2 = {nE_min} > 0, kein Horizont",
      in_su22 and sl2 and abs(tr2["E"] + 1) < 1e-12 and abs(tr2["Q"] - 1) < 1e-12 and abs(tr2["P_0"]) < 1e-12
      and spitzen == [-1.0, 1.0] and max(abs(x) for x in ebQ[0][1] + ebQ[1][1]) < 1e-9 and feld_ok
      and all(b == 0 for b in bif) and kappa_werte[0] == kappa_werte[1] and abs(kappa_werte[0]) == 1
      and w_chm == (1 - r_s**2)/2 and chm_int == 4*sp.pi/15 and evVE == [-1.0, 0.0, 0.0, 0.0, 0.0, 1.0]
      and evVQ == [-1.0, 0.0, 0.0, 0.0, 0.0, 1.0])

# --- 1.2 reelle Konjugation: unmöglich; Kontrollen
K_UHR, K_CAP = 4*np.pi*ME, 2*np.pi*MQ
vf = {str(c): verflechter(ME, c*MQ)[0] for c in (1, -1, 2, -2, 0.5, 4, 1j, -1j, 2j)}
spek_uhr = sorted(np.round(np.linalg.eigvals(K_UHR).imag/np.pi, 8).tolist())
spek_cap = sorted(np.round(np.linalg.eigvals(K_CAP).real/np.pi, 8).tolist())
gED = sla.expm(np.pi/2*ME)                                              # reelles Gruppenelement: Q^- -> D
qd_res = float(np.abs(gED @ MQ @ np.linalg.inv(gED) - MD).max())
gED_unit = float(np.abs(gED.conj().T @ G4 @ gED - G4).max())
a_ver, R_dil = np.array([0.3, 0.2, -0.1, 0.4]), 1.7
g_kegel = sla.expm(-sum(a_ver[m]*MP[m] for m in range(4))) @ sla.expm(np.log(R_dil)*MD)
Q_kegel = g_kegel @ MQ @ np.linalg.inv(g_kegel)
eb_kegel = eigenebenen(Q_kegel)
tips = sorted([[round(c.real, 8) for c in p[2]] for p in eb_kegel])
tips_soll = sorted([[round(a_ver[0] - R_dil, 8)] + [round(x, 8) for x in a_ver[1:]],
                    [round(a_ver[0] + R_dil, 8)] + [round(x, 8) for x in a_ver[1:]]])
g_kegel_unit = float(np.abs(g_kegel.conj().T @ G4 @ g_kegel - G4).max())
ZAHLEN["konjugation_reell"] = {"dim {g: g E = c Q^- g} je c": vf, "spec K_Uhr/(i pi)": spek_uhr, "spec K_cap/pi": spek_cap,
                               "Ad(e^(pi E/2)) Q^- = D, Residuum": qd_res, "Spitzen des konjugierten Kegels": tips,
                               "tr K_Uhr^2 / tr K_cap^2": float(np.trace(K_UHR @ K_UHR).real/np.trace(K_CAP @ K_CAP).real)}
check(f"[E, negativ] Reelle Konjugation: Der einseitige Modulgenerator des Nahtzustands ist K_Uhr = 4 pi E (z72, z83, "
      f"z104; auf der vollen Standardform ist der Tomita-Generator beta (L_0 - L-quer_0) = 2 beta h zentral, z83). K_Uhr "
      f"hat auf T das Spektrum i pi x {spek_uhr}, K_cap = 2 pi Q^- das Spektrum pi x {spek_cap}: Für kein reelles c gibt es "
      f"ein g in GL(4,C) mit g E g^-1 = c Q^- (Dimension des Verflechterraums {vf}; nur c = +-i erlaubt). In SU(2,2) und "
      f"selbst in SL(4,C) ist K_Uhr also zu keinem reellen Vielfachen eines Doppelkegelgenerators konjugiert "
      f"(tr K_Uhr^2/tr K_cap^2 = {ZAHLEN['konjugation_reell']['tr K_Uhr^2 / tr K_cap^2']:.1f}, Ad-invariant). "
      f"Kontrollen, dass der Test Konjugation findet: Ad(e^(pi E/2)) Q^- = D mit e^(pi E/2) in SU(2,2) (Residuen "
      f"{qd_res:.0e}, {gED_unit:.0e}); und jeder Doppelkegel ist nativ erreichbar: g = e^(-a.P) e^(ln R D) in SU(2,2) "
      f"({g_kegel_unit:.0e}) bringt Q^- auf den Generator des Kegels mit Zentrum a = {a_ver.tolist()} und Radius {R_dil} "
      f"(Spitzen {tips})",
      all(vf[str(c)] == 0 for c in (1, -1, 2, -2, 0.5, 4, 2j)) and vf["1j"] == 8 and vf["(-0-1j)"] == 8
      and spek_uhr == [-2.0, -2.0, 2.0, 2.0] and spek_cap == [-1.0, -1.0, 1.0, 1.0] and qd_res < 1e-12
      and gED_unit < 1e-12 and tips == tips_soll and g_kegel_unit < 1e-12)

# --- 1.3 komplexe Konjugation: das Wick-Element
W = sla.expm(1j*np.pi/2*MD)
Wi = np.linalg.inv(W)
wick_res = float(np.abs(W @ ME @ Wi - 1j*MQ).max())
wick_unit = float(np.abs(W.conj().T @ G4 @ W - G4).max())
wick_det = complex(np.linalg.det(W))
bild = []
for lam, gram, pkt in ebE:                                              # Eigenebenen von E -> Eigenebenen von Q^-
    wE, VEe = np.linalg.eig(ME)
    B = W @ np.linalg.qr(VEe[:, np.abs(wE - lam) < 1e-8])[0]
    bild.append((str(pkt[0]), str(punkt(np.linalg.qr(B)[0])[0])))
uhr_bild = float(np.abs(W @ K_UHR @ Wi - 2j*K_CAP).max())
ZAHLEN["wick"] = {"|W E W^-1 - i Q^-|": wick_res, "|W^+ G W - G|": wick_unit, "det W": str(np.round(wick_det, 12)),
                  "Vakuumpunkt x^0 -> Spitze x^0": bild, "|W K_Uhr W^-1 - 2 i K_cap|": uhr_bild}
check(f"[E] Komplexe Konjugation: Das einzige verbindende Element ist die Wick-Abbildung W = e^(i pi D/2) in SL(4,C) "
      f"(det W = {np.round(wick_det, 12).real:.0f}), W E W^-1 = i Q^- (Residuum {wick_res:.0e}); W liegt nicht in U(2,2) "
      f"(|W^+ G W - G| = {wick_unit:.2f}) und bildet die definiten Eigenebenen von E (Vakuumpunkte) auf die Nullebenen "
      f"von Q^- ab ({bild}: x^0 = +-i -> Spitzen t = -+1). Damit W K_Uhr W^-1 = 2 i K_cap (Residuum {uhr_bild:.0e}): "
      f"Zwischen dem nativen Modulgenerator und dem Hislop-Longo-Generator stehen ein Faktor i (Klassenwechsel "
      f"elliptisch -> hyperbolisch, keine Symmetrie) und ein Faktor 2",
      wick_res < 1e-12 and wick_unit > 0.5 and abs(wick_det - 1) < 1e-12 and uhr_bild < 1e-11
      and sorted(b[1] for b in bild) == sorted(str(complex(x, 0)) for x in (-1.0, 1.0)))

# --- 1.4 der Faktor 2
Id4, Id6 = np.eye(4), np.eye(6)
per = {"e^(2 pi E) auf T = -1": float(np.abs(sla.expm(2*np.pi*ME) + Id4).max()),
       "e^(4 pi E) auf T = +1": float(np.abs(sla.expm(4*np.pi*ME) - Id4).max()),
       "e^(2 pi E) auf R^(4,2) = +1": float(np.abs(sla.expm(2*np.pi*VE) - Id6).max()),
       "e^(2 pi i Q) auf T = -1": float(np.abs(sla.expm(2j*np.pi*MQ) + Id4).max()),
       "e^(4 pi i Q) auf T = +1": float(np.abs(sla.expm(4j*np.pi*MQ) - Id4).max()),
       "e^(2 pi i Q) auf R^(4,2) = +1": float(np.abs(sla.expm(2j*np.pi*VQ) - Id6).max()),
       "e^(pi i Q) auf R^(4,2) != 1": float(np.abs(sla.expm(1j*np.pi*VQ) - Id6).max()),
       "e^(4 pi i E) auf T nicht periodisch (Norm)": float(np.abs(sla.expm(4j*np.pi*ME)).max())}
woerterbuch = {"W e^(0,7 E) W^-1 = e^(0,7 i Q)": float(np.abs(W @ sla.expm(0.7*ME) @ Wi - sla.expm(0.7j*MQ)).max()),
               "W [Delta_Naht = e^(-4 pi E)] W^-1 = Kappenfluss bei t = 2":
                   float(np.abs(W @ sla.expm(4j*np.pi*ME) @ Wi - sla.expm(-4*np.pi*MQ)).max()/np.abs(sla.expm(-4*np.pi*MQ)).max()),
               "W [Uhr e^(2 pi i E) = Naht-Modulfluss bei t = 1/2] W^-1 = Delta_cap = e^(-2 pi Q)":
                   float(np.abs(W @ sla.expm(2*np.pi*ME) @ Wi - sla.expm(2j*np.pi*MQ)).max())}


def dh_matrix(h, N):
    """H_1 = (L_1 + L_-1)/2 auf D_h, <n+1|L_-1|n> = sqrt((n+1)(n+2h)) (z72, z104: Q^- auf D_h)."""
    n = np.arange(N - 1)
    off = np.sqrt((n + 1)*(n + 2*h))/2
    return np.diag(off, 1) + np.diag(off, -1)


sech = {}
for h in (0.5, 1.0):
    wq, Vq = np.linalg.eigh(dh_matrix(h, 600))
    fehler = max(abs(np.sum(np.abs(Vq[0])**2*np.exp(-1j*s*wq)) - np.cosh(s/2)**(-2*h)) for s in (0.5, 1.0, 2.0, 3.0))
    f_h = lambda s, h=h: mp.cosh(s/2)**(-2*h)
    forts2 = [complex(f_h(s - 2j*mp.pi)/f_h(s)) for s in (mp.mpf("0.3"), mp.mpf("1.1"))]
    forts4 = [complex(f_h(s - 4j*mp.pi)/f_h(s)) for s in (mp.mpf("0.3"), mp.mpf("1.1"))]
    sech[str(h)] = {"Fehler sech-Gesetz (N = 600)": float(fehler), "f(s - 2 pi i)/f(s)": [np.round(z, 12) for z in forts2],
                    "f(s - 4 pi i)/f(s)": [np.round(z, 12) for z in forts4]}
ALPHAS = (mp.mpf(1)/4, mp.mpf(3)/4)
ns_spek = all(mp.frac(2*(k + a)) == mp.mpf(1)/2 for a in ALPHAS for k in range(-5, 5))
quarter_per = [float(8*mp.pi*(k + a)/(2*mp.pi)) for a in ALPHAS for k in (0, 1)]       # 8 pi x in 2 pi Z?
n_4pi, n_2pi = 1/(1 + mp.exp(mp.pi)), 1/(1 + mp.exp(mp.pi/2))
ZAHLEN["faktor2"] = {"Perioden": per, "Wörterbuch": woerterbuch, "D_h": {k: {a: str(b) for a, b in v.items()} for k, v in sech.items()},
                     "2 x in Z + 1/2 (L_0^NS = 2 L_0)": bool(ns_spek), "8 pi x/(2 pi) für die Quartermoden": quarter_per,
                     "Besetzung der tiefsten Mode bei beta = 4 pi, 2 pi": [mp.nstr(n_4pi, 10), mp.nstr(n_2pi, 10)]}
check(f"[E] Der Faktor 2 (explizite Gruppenelemente): Die Uhr e^(theta E) hat auf Vektoren (R^(4,2), Raumzeit) die "
      f"reelle Periode 2 pi, auf Twistoren 4 pi (e^(2 pi E) = -1 auf T, zentral), auf den Quartermoden x in Z +- 1/4 der "
      f"Naht 8 pi (8 pi x/(2 pi) = {quarter_per}); die Kappe hat die imaginäre Periode 2 pi auf Vektoren und 4 pi auf "
      f"Twistoren (e^(2 pi i Q) = -1 auf T), Residuen <= {max(v for k, v in per.items() if 'nicht' not in k and '!=' not in k):.0e}. "
      f"In der unitären Darstellung D_h gilt <h|e^(-i s Q)|h> = cosh(s/2)^(-2h) (Fehler "
      f"{max(float(v['Fehler sech-Gesetz (N = 600)']) for v in sech.values()):.0e}), und die Fortsetzung s -> s - 2 pi i "
      f"gibt (-1)^(2h): +1 für h = 1 (Strom, Boson), -1 für h = 1/2 (Fermion), +1 erst bei 4 pi. Das ist die Spinor-"
      f"Halbwinkel-Aussage als Gruppenaussage: 4 pi = Spinorperiode, 2 pi = Vektorperiode. Als Erzeugerwahl: 4 pi L_0 = "
      f"2 pi L_0^NS mit L_0^NS = 2 L_0 (Spektrum Z + 1/2: {ns_spek}; tau = i). Wörterbuch unter W: Delta_cap = e^(-2 pi Q) "
      f"ist das Wick-Bild der vollen Uhrdrehung e^(2 pi i E) (Naht-Modulfluss bei t = 1/2), Delta_Naht = e^(-4 pi E) das "
      f"Wick-Bild des reellen Kappenflusses bei t = 2 (Residuen {max(woerterbuch.values()):.0e}). Als Zustandsaussage "
      f"bleibt es ein echter Unterschied: beta = 4 pi ist ein KMS-Datum des elliptischen E ohne imaginäre Periode "
      f"(|e^(4 pi i E)| = {per['e^(4 pi i E) auf T nicht periodisch (Norm)']:.0f} auf T), die fermionische KMS-Bedingung bei "
      f"4 pi verlangt Antiperiodizität bei 4 pi (h = 1/2 ist dort periodisch), und die Besetzung der tiefsten Mode wäre "
      f"{mp.nstr(n_4pi, 5)} (4 pi) gegen {mp.nstr(n_2pi, 5)} (2 pi)",
      max(v for k, v in per.items() if "nicht" not in k and "!=" not in k) < 1e-11 and per["e^(pi i Q) auf R^(4,2) != 1"] > 0.5
      and per["e^(4 pi i E) auf T nicht periodisch (Norm)"] > 100 and max(woerterbuch.values()) < 1e-9
      and max(float(v["Fehler sech-Gesetz (N = 600)"]) for v in sech.values()) < 1e-9
      and all(abs(z + 1) < 1e-12 for z in sech["0.5"]["f(s - 2 pi i)/f(s)"])
      and all(abs(z - 1) < 1e-12 for z in sech["0.5"]["f(s - 4 pi i)/f(s)"])
      and all(abs(z - 1) < 1e-12 for z in sech["1.0"]["f(s - 2 pi i)/f(s)"]) and ns_spek
      and all(abs(x - round(x)) < 1e-12 for x in quarter_per) and abs(n_4pi - n_2pi) > 0.1)


# --- 1.5 Peschel-Probe auf der halben Naht (1D-typisiert, wie z72)
def gitterzustand(M, alpha, beta, eps=1e-14):
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
    Wm = np.zeros((n, n), complex)
    Wm[:M, :M] = U
    Wm[M:, M:] = np.eye(m)
    return dict(M=M, xs=xs, G=G, W=Wm, labels=np.concatenate([xs, xs[mixed]]), alpha=alpha)


def hfun(nu, hmax=34.0):
    nu = np.clip(nu, 0, 1)
    with np.errstate(divide="ignore"):
        hv = np.log(1 - nu) - np.log(nu)
    return np.clip(hv, -hmax, hmax)


def halbe_naht(M, alpha, beta, nlo=3.0):
    """Peschel-Modulgenerator der halben Naht (-pi/2, pi/2) gegen 2 pi G_v (v Möbiusfeld, Nullstellen in den Enden)
    und gegen 4 pi L_0; dazu die Entropie des Schnitts."""
    st = gitterzustand(M, alpha, beta)
    G = st["G"]
    n = G.shape[0]
    A = np.unique(np.arange(-M//4, M//4) % M)
    Ac = np.setdiff1d(np.arange(n), A)
    k = np.zeros((n, n), complex)
    S_A = 0.0
    for R_, sg in ((A, 1.0), (Ac, -1.0)):
        w, V = np.linalg.eigh(G[np.ix_(R_, R_)])
        k[np.ix_(R_, R_)] = sg*(V*hfun(w)) @ V.conj().T
        if sg > 0:
            wc = np.clip(w, 1e-300, 1 - 1e-16)
            S_A = float(-(wc*np.log(wc) + (1 - wc)*np.log1p(-wc)).sum())
    K = st["W"].conj().T @ k @ st["W"]
    xs = st["xs"] - (alpha - 0.5)                                         # NS-Etiketten (Eichung, z72)
    L0 = np.diag(xs).astype(complex)
    L1 = np.zeros((M, M), complex)
    Lm1 = np.zeros((M, M), complex)
    for j in range(M):
        if j >= 1:
            L1[j - 1, j] = xs[j] - 0.5
        if j + 1 < M:
            Lm1[j + 1, j] = xs[j] + 0.5
    thm = -np.pi/M                                                       # Bogenmitte (Enden auf halben Plätzen)
    Gv = (np.exp(1j*thm)*L1 + np.exp(-1j*thm)*Lm1)/2                     # v = cos(theta - theta_m), Bogenlänge pi
    lo = np.where(np.abs(st["labels"][:M]) < nlo)[0]
    Kb = K[np.ix_(lo, lo)]
    lam_v, res_v = prop(Kb, 2*np.pi*Gv[np.ix_(lo, lo)])
    lam_E, res_E = prop(Kb, 4*np.pi*L0[np.ix_(lo, lo)])
    return {"lam (2 pi G_v)": lam_v.real, "Residuum (2 pi G_v)": res_v, "Residuum (4 pi L_0)": res_E, "S_A": S_A}


pes = {"Vakuum NS, M = 128": halbe_naht(128, 0.5, None), "Vakuum NS, M = 256": halbe_naht(256, 0.5, None),
       "Vakuum alpha = 1/4, M = 256": halbe_naht(256, 0.25, None),
       "Naht alpha = 1/4, beta = 4 pi, M = 128": halbe_naht(128, 0.25, 4*np.pi),
       "Naht alpha = 1/4, beta = 4 pi, M = 256": halbe_naht(256, 0.25, 4*np.pi),
       "Naht alpha = 3/4, beta = 4 pi, M = 256": halbe_naht(256, 0.75, 4*np.pi),
       "alpha = 1/4, beta = 32 pi, M = 256": halbe_naht(256, 0.25, 32*np.pi)}
ZAHLEN["peschel_halbe_naht"] = {k: {a: float(np.round(b, 10)) for a, b in v.items()} for k, v in pes.items()}
v256, n256 = pes["Vakuum NS, M = 256"], pes["Naht alpha = 1/4, beta = 4 pi, M = 256"]
check(f"[E numerisch, 1D-typisiert] Hislop-Longo als Resultat aus dem Zustand (Peschel, Gitter-Chiralfermion, wie z72): "
      f"Im Grundzustand der Uhr ist der Modulgenerator der halben Naht (-pi/2, pi/2) auf den niedrigen Moden lam x 2 pi "
      f"G_v mit dem Möbiusfeld v = cos theta (Q^- auf der Zentrumsweltlinie, z104): lam = {pes['Vakuum NS, M = 128']['lam (2 pi G_v)']:.5f} "
      f"(M = 128), {v256['lam (2 pi G_v)']:.5f} (M = 256), Residuum {pes['Vakuum NS, M = 128']['Residuum (2 pi G_v)']:.0e} -> "
      f"{v256['Residuum (2 pi G_v)']:.0e} (fällt mit M); für alpha = 1/4 "
      f"gleich ({pes['Vakuum alpha = 1/4, M = 256']['Residuum (2 pi G_v)']:.0e}). Im Nahtzustand (beta = 4 pi) ist er "
      f"nicht geometrisch: Residuum gegen 2 pi G_v {n256['Residuum (2 pi G_v)']:.2f} (alpha = 3/4: "
      f"{pes['Naht alpha = 3/4, beta = 4 pi, M = 256']['Residuum (2 pi G_v)']:.2f}), stabil in M "
      f"({pes['Naht alpha = 1/4, beta = 4 pi, M = 128']['Residuum (2 pi G_v)']:.2f} bei M = 128), und auch kein "
      f"Vielfaches der Uhr (Residuum gegen 4 pi L_0 {n256['Residuum (4 pi L_0)']:.2f}); erst für beta -> unendlich kehrt "
      f"die Geometrie zurück (beta = 32 pi: {pes['alpha = 1/4, beta = 32 pi, M = 256']['Residuum (2 pi G_v)']:.0e})",
      abs(v256["lam (2 pi G_v)"] - 1) < 1e-3 and abs(v256["lam (2 pi G_v)"] - 1) < abs(pes["Vakuum NS, M = 128"]["lam (2 pi G_v)"] - 1)
      and v256["Residuum (2 pi G_v)"] < 5e-4 and v256["Residuum (2 pi G_v)"] < pes["Vakuum NS, M = 128"]["Residuum (2 pi G_v)"]/2
      and pes["Vakuum alpha = 1/4, M = 256"]["Residuum (2 pi G_v)"] < 5e-4
      and n256["Residuum (2 pi G_v)"] > 0.2 and pes["Naht alpha = 3/4, beta = 4 pi, M = 256"]["Residuum (2 pi G_v)"] > 0.2
      and abs(n256["Residuum (2 pi G_v)"] - pes["Naht alpha = 1/4, beta = 4 pi, M = 128"]["Residuum (2 pi G_v)"]) < 0.01
      and n256["Residuum (4 pi L_0)"] > 0.1 and pes["alpha = 1/4, beta = 32 pi, M = 256"]["Residuum (2 pi G_v)"] < 5e-4)
befund("E/C: Nativ (RESULT) sind: die radiale sl(2,R) in su(2,2); der hyperbolische Erzeuger Q^- mit seinem Doppelkegel "
       "(Zentrum = tau-Fixpunkt, Radius R_Vak, kappa = 1, also 2 pi Q^- in Hislop-Longo-Normierung); die Konjugation aller "
       "Doppelkegel untereinander in SU(2,2); im Grundzustand der Uhr der geometrische Modulfluss der halben Naht (Peschel, "
       "1D). Nicht nativ: (a) der Modulgenerator des Nahtzustands ist 4 pi E (elliptisch) bzw. zentral und zu keinem "
       "Doppelkegelgenerator in SU(2,2) konjugiert; nur die Wick-Abbildung W = e^(i pi D/2) (kein Gruppenelement) verbindet "
       "sie, mit W K_Uhr W^-1 = 2 i K_cap; (b) der Faktor 2 ist als Gruppenaussage der Spinor-Halbwinkel (4 pi auf "
       "Twistoren, 2 pi auf Vektoren) und als Erzeugerwahl L_0^NS = 2 L_0, als Zustandsaussage aber ein echter "
       "Unterschied (kein gemeinsamer KMS-Zustand); (c) IMPORT nach GRAV.NONCIRCULAR.01 bleibt: die Lokalität eines 4D-"
       "Netzes (Kappenalgebra A(O_1), Reeh-Schlieder, Haag-Dualität), ohne die 'Modulfluss der Kappe = 2 pi Q^-' in 3+1 "
       "nur der Satz von Hislop-Longo für das freie Leiterfeld ist (z104: C nach dem Netz), und er gilt für das Vakuum, "
       "nicht für den Nahtzustand. QGEO.MODHAM.01 ist damit nicht geschlossen", True)

# =================================================================== Teil 2
say("\nTeil 2: Erster Hauptsatz delta S = delta<K> im Fockraum der Naht")
NM = 6
XQ = np.arange(-3, 3) + 0.25                                             # Quarter alpha = 1/4: x = k + 1/4
sz_, sm_, s0_ = np.diag([1., -1.]), np.array([[0, 1], [0, 0.]]), np.eye(2)


def kr(*m):
    out = m[0]
    for x in m[1:]:
        out = np.kron(out, x)
    return out


CA = [kr(*([sz_]*j + [sm_] + [s0_]*(NM - j - 1))) for j in range(NM)]    # Jordan-Wigner
CD = [a.T for a in CA]
DIMF = 2**NM
car = max(float(np.abs(CA[i] @ CD[j] + CD[j] @ CA[i] - (i == j)*np.eye(DIMF)).max()) for i in range(NM) for j in range(NM))


def dGamma(h):
    return sum(h[i, j]*(CD[i] @ CA[j]) for i in range(NM) for j in range(NM) if abs(h[i, j]) > 0)


def zustand(H):
    w, V = np.linalg.eigh(H)
    p = np.exp(-(w - w.min()))
    p /= p.sum()
    return (V*p) @ V.conj().T


def entropie(rho):
    p = np.linalg.eigvalsh(rho)
    p = p[p > 1e-300]
    return float(-(p*np.log(p)).sum())


def hb(n):
    return -(n*mp.log(n) + (1 - n)*mp.log(1 - n))


def s_mode(y):
    """Entropie einer Fermimode mit modularer Energie y (symmetrisch, stabil): ln(1 + e^-|y|) + |y|/(e^|y| + 1)."""
    y = abs(y)
    return mp.log(1 + mp.exp(-y)) + y/(mp.exp(y) + 1)


K_F = dGamma(np.diag(4*np.pi*XQ).astype(complex))
RHO = zustand(K_F)
S_F, K_F0 = entropie(RHO), float(np.trace(RHO @ K_F).real)
S_peschel = float(mp.fsum(s_mode(4*mp.pi*mp.mpf(x)) for x in XQ.tolist()))
P_F = np.diag(RHO).real
diag_ok = float(np.abs(RHO - np.diag(np.diag(RHO))).max())
check(f"[E] Fockraum eines Quarters (Moden x = k + 1/4, k = -3..2, Jordan-Wigner, dim {DIMF}; CAR-Residuum {car:.0e}): "
      f"rho = e^(-K)/Z mit K = 4 pi Summe x c_x^+ c_x ist diagonal ({diag_ok:.0e}), S(rho) = {S_F:.12f} = Summe der "
      f"Modenentropien (Peschel, mpmath) {S_peschel:.12f}; das andere Quarter (alpha = 3/4) ist das Teilchen-Loch-Bild "
      f"x -> -x mit gleichem Spektrum",
      car < 1e-14 and diag_ok < 1e-14 and abs(S_F - S_peschel) < 1e-12)

dh_list = []
for _ in range(4):
    Z_ = rng.normal(size=(NM, NM)) + 1j*rng.normal(size=(NM, NM))
    dh_list.append((Z_ + Z_.conj().T)/2)
X_GAUSS = [dGamma(dh) for dh in dh_list]                                 # Einteilchen-Kovarianzstörungen (delta h)
X_QUART = []
for _ in range(3):
    Xq = sum(rng.normal()*(CD[i] @ CA[i] @ CD[j] @ CA[j]) for i in range(NM) for j in range(i + 1, NM))
    for (i, j, k, l) in ((0, 1, 2, 3), (2, 3, 4, 5), (1, 2, 3, 4)):
        T_ = rng.normal()*(CD[i] @ CD[j] @ CA[k] @ CA[l])
        Xq = Xq + T_ + T_.conj().T
    X_QUART.append(Xq)
L1F = np.zeros((NM, NM))
for k in range(1, NM):
    L1F[k - 1, k] = XQ[k] - 0.5                                          # L_1 e_x = (x - 1/2) e_(x-1) (z72)
H1F = dGamma(((L1F + L1F.T)/2).astype(complex))                          # Möbius H_1 = (L_1 + L_-1)/2, hyperbolisch
EPS = (4e-3, 2e-3, 1e-3)


def variation(X, eps, Kop):
    r = zustand(K_F + eps*X)
    return entropie(r) - S_F, float(np.trace(r @ Kop).real - np.trace(RHO @ Kop).real)


def bkm_varianz(X):
    Pi, Pj = np.meshgrid(P_F, P_F, indexing="ij")
    with np.errstate(divide="ignore", invalid="ignore"):
        ker = np.where(np.abs(Pi - Pj) > 1e-13*np.maximum(Pi, Pj), (Pi - Pj)/(np.log(Pi) - np.log(Pj)), Pi)
    return float((np.abs(X)**2*ker).sum() - np.trace(RHO @ X).real**2)


erst = {}
for lab, Xs in (("Kovarianz (Gauß)", X_GAUSS), ("Quartik", X_QUART)):
    for nr, X in enumerate(Xs):
        zeile = []
        for eps in EPS:
            dS, dK = variation(X, eps, K_F)
            zeile.append(dS/dK - 1)
        erst[f"{lab} {nr}"] = zeile
max_abw = max(abs(z[-1]) for z in erst.values())
linear = all(1.9 < z[0]/z[1] < 2.1 and 1.9 < z[1]/z[2] < 2.1 for z in erst.values())
ZAHLEN["erster_hauptsatz"] = {k: [float(np.round(x, 10)) for x in v] for k, v in erst.items()}
check(f"[E] Erster Hauptsatz: Für rho_eps = e^(-(K + eps X))/Z mit vier Einteilchen-Kovarianzstörungen X = dGamma(delta h) "
      f"(hermitesch, nicht diagonal) und drei quartischen (nicht quasifreien) Störungen ist delta S/delta<K> - 1 = O(eps): "
      f"bei eps = 1e-3 höchstens {max_abw:.1e}, und die Abweichung halbiert sich mit eps (Verhältnisse in (1,9; 2,1): "
      f"{linear}); K = 4 pi L_0 ist der Modulhamiltonoperator, die Kopplung ist für alle Störungsarten dieselbe",
      max_abw < 0.02 and linear)

K_2PI, K_E = K_F/2, K_F/(4*np.pi)
K_HYP = 2*np.pi*H1F                                                      # 'Kappengenerator' 2 pi H_1 auf denselben Moden
neg = {"beta = 2 pi": [], "E allein": [], "delta<2 pi H_1>/delta S (Gauß)": [], "delta<2 pi H_1>/delta S (Quartik)": []}
for nr, X in enumerate(X_GAUSS + X_QUART):
    dS, dK2 = variation(X, 1e-3, K_2PI)
    neg["beta = 2 pi"].append(dS/dK2)
    neg["E allein"].append(dS/variation(X, 1e-3, K_E)[1])
    neg["delta<2 pi H_1>/delta S (Gauß)" if nr < len(X_GAUSS) else "delta<2 pi H_1>/delta S (Quartik)"].append(
        variation(X, 1e-3, K_HYP)[1]/dS)
hyp_g, hyp_q = neg["delta<2 pi H_1>/delta S (Gauß)"], neg["delta<2 pi H_1>/delta S (Quartik)"]
inv_hyp = float(np.abs(komm(RHO, K_HYP)).max())
ZAHLEN["erster_hauptsatz_negativ"] = {k: [float(np.round(x, 6)) for x in v] for k, v in neg.items()}
check(f"[E, Negativkontrolle] Falsches K scheitert: mit beta = 2 pi (Bisognano-Wichmann-Wert auf dem Erzeuger L_0) ist "
      f"delta S/delta<K> = {min(neg['beta = 2 pi']):.3f}..{max(neg['beta = 2 pi']):.3f} (2 statt 1), mit E allein "
      f"{min(neg['E allein']):.2f}..{max(neg['E allein']):.2f} (4 pi), und mit dem hyperbolischen Generator 2 pi H_1 "
      f"(Kappenform auf denselben Moden) ist delta<K>/delta S keine Konstante: {[float(np.round(x, 3)) for x in hyp_g]} für "
      f"die vier Kovarianzstörungen (beide Vorzeichen) und 0 für die quartischen (|.| <= {max(abs(x) for x in hyp_q):.0e}); "
      f"der Zustand ist unter ihm nicht invariant (|[rho, K]| = {inv_hyp:.2f})",
      all(abs(x - 2) < 0.06 for x in neg["beta = 2 pi"]) and all(abs(x - 4*np.pi) < 0.4 for x in neg["E allein"])
      and all(abs(x - 1) > 0.5 for x in neg["beta = 2 pi"] + neg["E allein"]) and all(abs(x - 1) > 0.3 for x in hyp_g)
      and min(hyp_g) < 0 < max(hyp_g) and max(abs(x) for x in hyp_q) < 1e-9 and inv_hyp > 1e-3)

zweit = {}
srel_min = np.inf
for lab, Xs in (("Kovarianz (Gauß)", X_GAUSS), ("Quartik", X_QUART)):
    for nr, X in enumerate(Xs):
        F = bkm_varianz(X)
        werte = []
        for eps in (2e-3, 1e-3, -1e-3):
            dS, dK = variation(X, eps, K_F)
            werte.append((dK - dS)/eps**2)
            srel_min = min(srel_min, dK - dS)
        zweit[f"{lab} {nr}"] = {"S_rel/eps^2": [float(np.round(x, 8)) for x in werte], "F/2 (BKM)": float(np.round(F/2, 8))}
zweit_ok = all(abs(v["S_rel/eps^2"][1]/v["F/2 (BKM)"] - 1) < 5e-3 and abs(v["S_rel/eps^2"][2]/v["F/2 (BKM)"] - 1) < 5e-3
               and v["F/2 (BKM)"] > 0 for v in zweit.values())
Pi_, Pj_ = np.meshgrid(P_F, P_F, indexing="ij")
with np.errstate(divide="ignore", invalid="ignore"):
    F_unit = float((np.abs(H1F)**2*(Pi_ - Pj_)*(np.log(Pi_) - np.log(Pj_))).sum())
unit = []
for eps in (4e-3, 2e-3, 1e-3):
    U_ = sla.expm(1j*eps*H1F)
    r_ = U_ @ RHO @ U_.conj().T
    unit.append((entropie(r_) - S_F, float(np.trace(r_ @ K_F).real - K_F0)/eps**2))
nx = [1/(1 + mp.exp(4*mp.pi*mp.mpf(x))) for x in XQ.tolist()]
hx = [4*mp.pi*mp.mpf(x) for x in XQ.tolist()]
kern_min = min([nx[i]*(1 - nx[i]) for i in range(NM)]
               + [(nx[j] - nx[i])/(hx[i] - hx[j]) for i in range(NM) for j in range(NM) if i != j])
ZAHLEN["zweite_ordnung"] = {"Störungen": zweit, "Möbius: (delta S, delta<K>/eps^2)": [(float(a), float(np.round(b, 10))) for a, b in unit],
                            "Möbius F/2": F_unit/2, "min S_rel": float(srel_min), "kleinster Fisher-Kern (Einteilchen)": mp.nstr(kern_min, 8)}
check(f"[E] Zweite Ordnung: delta<K> - delta S = S(rho_eps || rho) >= 0 für alle Störungen und beide Vorzeichen von eps "
      f"(Minimum {srel_min:.2e} > 0), und S_rel/eps^2 -> F/2 mit der Bogoliubov-Kubo-Mori-Fisher-Information F > 0 "
      f"(unabhängig aus dem Spektrum gerechnet; Abweichung < 0,5 %: {zweit_ok}). Native Möbiusstörung rho -> e^(i eps "
      f"dGamma(H_1)) rho e^(-i eps dGamma(H_1)): delta S = 0 ({max(abs(a) for a, _ in unit):.0e}), delta<K> = eps^2 x "
      f"{unit[-1][1]:.8f} = eps^2 F/2 ({F_unit/2:.8f}) > 0. Der Fisher-Kern der Einteilchenstörungen ist positiv definit "
      f"(kleinster Wert {mp.nstr(kern_min, 3)} > 0): Der Nahtzustand maximiert S bei festem <K>, die Bedingung, die "
      f"Jacobsons 'maximale Vakuumverschränkung' benutzt",
      srel_min > 0 and zweit_ok and max(abs(a) for a, _ in unit) < 1e-12 and abs(unit[-1][1]/(F_unit/2) - 1) < 1e-5
      and F_unit > 0 and kern_min > 0)

phi_falsch = []
for X in X_GAUSS + X_QUART:
    for eps in (1e-3, -1e-3):
        dS, dK2 = variation(X, eps, K_2PI)
        phi_falsch.append(dK2 - dS)
check(f"[E, Negativkontrolle] Mit beta = 2 pi ist delta<K'> - delta S nicht positiv: über dieselben Störungen und beide "
      f"Vorzeichen liegt es zwischen {min(phi_falsch):.2e} und {max(phi_falsch):.2e} ({sum(x < 0 for x in phi_falsch)} von "
      f"{len(phi_falsch)} negativ, linear in eps): rho ist kein stationärer Punkt von S - <K'>, also kein Gleichgewicht "
      f"zur Temperatur 1/(2 pi)",
      min(phi_falsch) < -1e-6 and max(phi_falsch) > 1e-6 and sum(x < 0 for x in phi_falsch) == len(phi_falsch)//2)


def gb(n):
    return mp.log1p(n) + n*mp.log1p(n) - n*mp.log(n)                    # (1 + n) ln(1 + n) - n ln n, stabil für kleine n


def leiter_S(beta, dn=None):
    """Thermische Entropie des Leiterfelds (z83: Stufe m mit m^2 Zuständen), optional mit Besetzungsstörung."""
    tot = mp.mpf(0)
    for m in range(1, 40):
        n = 1/(mp.exp(beta*m) - 1)
        if dn and m in dn:
            n += dn[m]
        tot += m**2*gb(n)
    return tot


bos = {}
for bl, beta in (("4 pi", 4*mp.pi), ("8 pi", 8*mp.pi)):
    abw = []
    for m in (1, 2, 3):
        n0 = 1/(mp.exp(beta*m) - 1)
        d = n0*mp.mpf("1e-8")
        dS = m**2*(gb(n0 + d) - gb(n0 - d))                                # nur die gestörte Stufe (die übrigen heben sich exakt auf)
        abw.append(abs(dS/(2*d*m**2*m) / beta - 1))                       # delta S = beta delta<E>, delta<E> = m^2 m delta n
    bos[bl] = float(max(abw))
n8 = 1/(mp.exp(8*mp.pi) - 1)
bos_kreuz = float(abs((gb(n8*(1 + mp.mpf("1e-8"))) - gb(n8*(1 - mp.mpf("1e-8"))))/(2*n8*mp.mpf("1e-8")*4*mp.pi) - 1))
ZAHLEN["leiter"] = {"|dS/(beta d<E>) - 1| bei beta": bos, "Kreuzkontrolle (Zustand 8 pi, K mit 4 pi)": bos_kreuz}
check(f"[E] Gleiche Uhr im bosonischen Leiterfeld (z83: Stufe E = m mit m^2 Zuständen): dS = beta d<E> je Stufe m = 1, 2, 3 "
      f"bei beta = 4 pi (einseitig) und beta' = 8 pi (Produktlesart), Abweichung <= {max(bos.values()):.0e}; Kreuzkontrolle: "
      f"der 8-pi-Zustand mit dem 4-pi-Generator verfehlt um {bos_kreuz:.2f} (Faktor 2)",
      max(bos.values()) < 1e-10 and abs(bos_kreuz - 1) < 0.05)
befund("E: Der erste Hauptsatz delta S = delta<K> gilt nativ, exakt in erster Ordnung, für Kovarianz-, Quartik- und "
       "Möbiusstörungen, mit K = 4 pi L_0 und positiver zweiter Ordnung (relative Entropie = Fisher-Information). Er ist "
       "aber der Hauptsatz des globalen elliptischen Generators (Thermofeldschnitt Original|Partner, z86), nicht der eines "
       "Doppelkegels: Es gibt nativ keine Zerlegung delta<K> = 2 pi int xi.T mit dem Kappenfeld xi_Q und keinen Flächenterm. "
       "Der Bisognano-Wichmann-Wert 2 pi auf demselben Erzeuger scheitert um den Faktor 2", True)

# =================================================================== Teil 3
say("\nTeil 3: Native Entropie s_0 und der Wert von G")


def s_quarter(alpha, beta=4*mp.pi, kmax=40):
    return mp.fsum(s_mode(beta*(k + alpha)) for k in range(-kmax, kmax))


sQ = [s_quarter(a) for a in ALPHAS]


def eta_fn(y):
    q = mp.exp(-2*mp.pi*y)
    return q**(mp.mpf(1)/24)*mp.nprod(lambda n: 1 - q**n, [1, mp.inf])


def E4_fn(y):
    q = mp.exp(-2*mp.pi*y)
    return 1 + 240*mp.nsum(lambda n: n**3*q**n/(1 - q**n), [1, mp.inf])


lnZ_E8 = lambda y: mp.log(E4_fn(y)/eta_fn(y)**8)
dlnZ = mp.diff(lnZ_E8, 1)
S_E8 = lnZ_E8(1) - dlnZ
S_leiter = {"4 pi": leiter_S(4*mp.pi), "8 pi": leiter_S(8*mp.pi)}
LN2 = mp.log(2)
S0 = [("Quarter (ein Majorana-NS-Spektrum, alpha = 1/4)", sQ[0], "E"),
      ("Naht (beide Quarter)", sQ[0] + sQ[1], "E"),
      ("Träger 5 (g_car = 5 Nähte)", 5*(sQ[0] + sQ[1]), "L"),
      ("16 Quarter ohne GSO (= c = 8 Nähte)", 16*sQ[0], "L"),
      ("E8-Naht (16 Quarter mit GSO)", S_E8, "E/L"),
      ("drei Familien: 3 E8-Nähte", 3*S_E8, "L"),
      ("drei Familien: 48 Quarter ohne GSO", 48*sQ[0], "L"),
      ("Leiterfeld (S2-Lesart), beta = 4 pi", S_leiter["4 pi"], "C"),
      ("Leiterfeld (S2-Lesart), beta' = 8 pi", S_leiter["8 pi"], "C")]
ZAHLEN["s0"] = {k: mp.nstr(v, 20) for k, v, _ in S0}
check(f"[E] Native Entropie (KMS beta = 4 pi über das Einteilchenspektrum 4 pi |x|, x in Z + alpha; Original|Partner): je "
      f"Quarter s = {mp.nstr(sQ[0], 15)} = ln 2/4 für alpha = 1/4 und 3/4 (Abweichung {mp.nstr(max(abs(s - LN2/4) for s in sQ), 3)}), "
      f"je Naht ln 2/2 = {mp.nstr(LN2/2, 12)}; E8-Naht: Z = E_4/eta^8 bei tau = i gleich {mp.nstr(mp.exp(lnZ_E8(1)), 12)}, "
      f"d ln Z/dy = {mp.nstr(abs(dlnZ), 2)} (S-Invarianz), also S = ln 12 = {mp.nstr(S_E8, 12)}; Leiterfeld (m^2 Zustände je "
      f"Stufe): {mp.nstr(S_leiter['4 pi'], 6)} (4 pi), {mp.nstr(S_leiter['8 pi'], 6)} (8 pi). Die Vielfachheiten (5, 16, 48, "
      f"drei E8-Nähte) sind Zählungen unabhängiger Kopien (Lesart), keine weitere Rechnung",
      max(abs(s - LN2/4) for s in sQ) < mp.mpf(10)**-35 and abs(S_E8 - mp.log(12)) < mp.mpf(10)**-30
      and abs(dlnZ) < mp.mpf(10)**-30 and abs(S_leiter["4 pi"] - mp.exp(-4*mp.pi)*(1 + 4*mp.pi)) < mp.mpf(10)**-8)

# --- 3.2 Fläche der Bifurkationssphäre je konformem Rahmen
rho_s = sp.symbols("rho", positive=True)
OM = {"flach (P_0-Zeit, Lambda = 0)": sp.Integer(1), "Einstein-Universum R x S^3 (Radius R_Vak)": (1 + rho_s**2)/2,
      "dS (Omega = 1 - det x, l = R_Vak/2; braucht P_cap)": 1 + rho_s**2,
      "AdS (Omega = 1 + det x, l = R_Vak/2; nativ, sigma = +1)": 1 - rho_s**2}
FLAECHE = {k: sp.simplify(4*sp.pi*rho_s**2/Om**2) for k, Om in OM.items()}
A_bif = {k: sp.limit(A, rho_s, 1, "-") for k, A in FLAECHE.items()}
dS_max = sp.solve(sp.diff(FLAECHE["dS (Omega = 1 - det x, l = R_Vak/2; braucht P_cap)"], rho_s), rho_s)
killing_dS = sp.simplify(sum(xiQ[m]*sp.diff(1 - XM.det(), XS[m]) for m in range(4))
                         - sum(sp.diff(xiQ[m], XS[m]) for m in range(4))*(1 - XM.det())/4)
killing_AdS = sp.simplify(sum(xiQ[m]*sp.diff(1 + XM.det(), XS[m]) for m in range(4))
                          - sum(sp.diff(xiQ[m], XS[m]) for m in range(4))*(1 + XM.det())/4)
ZAHLEN["flaeche"] = {k: str(v) for k, v in A_bif.items()}
ZAHLEN["flaeche"]["xi_Q Killing im dS-Rahmen (Rest)"] = str(killing_dS)
ZAHLEN["flaeche"]["xi_Q Killing im AdS-Rahmen (Rest)"] = str(killing_AdS)
check(f"[E] Die Fläche der Bifurkationssphäre (t = 0, r = R_Vak) ist nicht konform invariant; in Einheiten R_Vak^2 "
      f"(g = eta/Omega^2, sympy): flach {A_bif['flach (P_0-Zeit, Lambda = 0)']}, Einstein-Universum "
      f"{A_bif['Einstein-Universum R x S^3 (Radius R_Vak)']} (Äquator der S^3), dS "
      f"{A_bif['dS (Omega = 1 - det x, l = R_Vak/2; braucht P_cap)']} = 4 pi l^2 (kosmologischer Horizont, maximale Sphäre bei "
      f"rho = {dS_max}; xi_Q ist dort Killingfeld, Rest {killing_dS}), im nativen AdS-Rahmen "
      f"{A_bif['AdS (Omega = 1 + det x, l = R_Vak/2; nativ, sigma = +1)']}: Die Bifurkationssphäre liegt auf dem konformen Rand, "
      f"und xi_Q ist dort kein Killingfeld (Rest {sp.factor(killing_AdS)})",
      A_bif["flach (P_0-Zeit, Lambda = 0)"] == 4*sp.pi and A_bif["Einstein-Universum R x S^3 (Radius R_Vak)"] == 4*sp.pi
      and A_bif["dS (Omega = 1 - det x, l = R_Vak/2; braucht P_cap)"] == sp.pi and dS_max == [1]
      and A_bif["AdS (Omega = 1 + det x, l = R_Vak/2; nativ, sigma = +1)"] == sp.oo and killing_dS == 0 and killing_AdS != 0)

# --- 3.3 räumlicher Schnitt: keine endliche native Dichte
S_schnitt = {"Vakuum": [halbe_naht(M, 0.5, None)["S_A"] for M in (64, 128, 256)],
             "Naht (alpha = 1/4, beta = 4 pi)": [halbe_naht(M, 0.25, 4*np.pi)["S_A"] for M in (64, 128, 256)]}
zuw = {k: [v[1] - v[0], v[2] - v[1]] for k, v in S_schnitt.items()}
ZAHLEN["raeumlicher_schnitt"] = {"S_A (M = 64, 128, 256)": {k: [float(np.round(x, 8)) for x in v] for k, v in S_schnitt.items()},
                                 "Zuwachs je Verdopplung": {k: [float(np.round(x, 8)) for x in v] for k, v in zuw.items()},
                                 "(1/3) ln 2": float(mp.log(2)/3)}
check(f"[E numerisch, 1D-typisiert] Die Verschränkungsentropie eines räumlichen Schnitts (halbe Naht) ist nicht endlich: "
      f"S_A wächst je Verdopplung der Auflösung um {zuw['Vakuum'][1]:.4f} (Vakuum) bzw. "
      f"{zuw['Naht (alpha = 1/4, beta = 4 pi)'][1]:.4f} (Nahtzustand), also um (1/3) ln 2 = {float(mp.log(2)/3):.4f} "
      f"(logarithmische Divergenz des Gitterschnitts mit zwei Enden). Endlich ist nur die Thermofeldentropie Original|Partner (ln 2/2); eine Entropie 'der "
      f"Bifurkationssphäre' gibt es nativ nicht ohne UV-Länge (z86, z90)",
      all(abs(x - float(mp.log(2)/3)) < 5e-3 for v in zuw.values() for x in v))

# --- 3.4 Überbestimmungstabelle
C3 = 1/(8*mp.pi)
ALPHA_INV = mp.mpf("137.035999177")                                     # CODATA 2022, nur Größenordnung (wie z98, z104)
z118 = json.loads((HERE / "z118_ergebnis.json").read_text())
c1_z118 = complex(z118["zahlen"]["mm"]["c1 = L1/(R det e)"][0]).real      # 64 eps je Index
lam_s, Rv, Gs, ws, ls, bs = sp.symbols("lambda R_Vak G w l b", positive=True)
G_mm = sp.solve(sp.Eq(ws*sp.nsimplify(c1_z118)*4/Rv**2, 1/(16*sp.pi*Gs)), Gs)[0]         # w c1 lambda = 1/(16 pi G)
b_mm = sp.simplify((ls**2/(64*sp.pi*Gs)).subs({Gs: G_mm, ls: Rv/2}))                       # b = l^2/(64 pi G)
v2_z118 = sp.solve(sp.Eq(2*sp.Rational(1, 8)/sp.pi/4*lam_s, 1/(16*sp.pi*Gs)), lam_s)[0]   # v_geo^2 aus (2 c_3/4) v^2 = 1/(16 pi G)
vR = sp.sqrt(sp.simplify((v2_z118*Rv**2).subs(Gs, G_mm.subs(ws, sp.Rational(1, 4)/sp.pi))))  # v_geo R_Vak unter (a2)
TAB = []


def zeile(quelle, lesart, wert, status, bem):
    TAB.append({"Quelle": quelle, "Lesart": lesart, "G/R_Vak^2": wert, "Status": status, "Bemerkung": bem})


zeile("(a) z118", "(a1) Gewicht (2 c_3/4) v_geo^2 = 1/(16 pi G)", None, "E/C",
      "G = v_geo^-2; keine reine Zahl ohne v_geo R_Vak (zweite Länge)")
zeile("(a) z118", "(a2) reines Gewicht w = 1 (Index) auf 64 eps, lambda = 4/R_Vak^2", mp.mpf(1)/(4096*mp.pi), "C",
      "Spurnormierung von z118 (b = 16 w); Erzeugernormierung ändert um 16")
zeile("(a) z118", "(a2) reines Gewicht w = 2 c_3 auf 64 eps, lambda = 4/R_Vak^2", mp.mpf(1)/1024, "C",
      "entspricht v_geo R_Vak = 32; S_dS wäre 256 pi")
zeile("(a) z118", "(a2) reines Gewicht w = c_3 (Majorana-Weyl)", mp.mpf(1)/512, "C", "ein J-Summand")
zeile("(b) Ledger v358/v359", "8 pi G = 1/c_3 in Einheiten v_geo = 1", None, "E (Definition der Einheit)",
      "G v_geo^2 = 1; enthält R_Vak nicht")
TARGET = mp.exp(-2*ALPHA_INV)/(128*mp.pi**3)
zeile("(b') Ledger COSMO.LAM.01 / LAMBDA.BRANCH.01", "rho_Lambda/Mbar^4 = (3/4 pi^2) e^(-2/alpha), Lambda = +12/R_Vak^2",
      TARGET, "Numerical (Vergleichsziel)", "dS-Relation; nativ ist AdS")
ZETA3 = mp.zeta(3)


def s_flaeche(x):
    return mp.log(1 + mp.exp(-x)) + x/(mp.exp(x) + 1)


face_exakt = lambda ell: mp.pi/mp.nsum(lambda l: 4*(l + 1)*s_flaeche(2*ell*(l + 1)), [0, mp.inf])
zeile("(c) Flächenzustand", "G = l^2/(18 zeta(3) c_3), l frei", None, "E/O", "G/l^2 = " + mp.nstr(1/(18*ZETA3*C3), 6))
zeile("(c) Flächenzustand", "l = R_Vak/2 (MM-Länge), asymptotisch", (mp.mpf(1)/4)/(18*ZETA3*C3), "L",
      "exakte Sphärensumme: " + mp.nstr(face_exakt(mp.mpf(1)/2), 6))
zeile("(c) Flächenzustand", "2 l = beta R_Vak = 4 pi R_Vak (Kragenweite = KMS-Länge), exakte Summe", face_exakt(2*mp.pi), "L",
      "außerhalb des Flächengesetzes (l >> R); asymptotische Formel gäbe " + mp.nstr(4*mp.pi**2/(18*ZETA3*C3), 5))
for g2l, g2 in (("1", mp.mpf(1)), ("2", mp.mpf(2)), ("2 pi", 2*mp.pi), ("4 pi", 4*mp.pi)):
    for Ll, Lv in (("R_Vak/2", mp.mpf(1)/2), ("R_Vak", mp.mpf(1)), ("R_Vak/tanh(pi)", 1/mp.tanh(mp.pi))):
        zeile("(d) z112", f"L = kappa/(2 g^2) = {Ll}, g^2 = {g2l}", g2**2*Lv**2/(8*mp.pi), "C (P_L)", "G = g^4 L^2/(8 pi)")
for name, s0, st in S0:
    zeile("(e) Jacobson S = A/(4G)", f"s_0 = {name}; A = 4 pi R_Vak^2 (flach/ESU)", mp.pi/s0, "C (P_S)", "Lambda = 0 bzw. Rahmenwahl")
    zeile("(e) Jacobson S = A/(4G)", f"s_0 = {name}; A = pi R_Vak^2 (dS)", mp.pi/(4*s0), "C (P_S, P_cap)", "Lambda = +12/R_Vak^2")
zeile("(e) Jacobson S = A/(4G)", "nativer AdS-Rahmen: A = unendlich", None, "E negativ", "kein endliches G für endliches s_0")
say("   Überbestimmungstabelle (reine Zahl G/R_Vak^2; '--' = Umrechnung unmöglich):")
for z in TAB:
    wert = "--" if z["G/R_Vak^2"] is None else f"{mp.nstr(z['G/R_Vak^2'], 6)} (10^{float(mp.log10(z['G/R_Vak^2'])):+.2f})"
    say(f"     {z['Quelle']:<44} | {z['Lesart']:<86} | {wert:<28} | {z['Status']:<26} | {z['Bemerkung']}")
ZAHLEN["tabelle"] = [{k: (mp.nstr(v, 12) if isinstance(v, mp.mpf) else v) for k, v in z.items()} for z in TAB]
rho_L = 3/(4*mp.pi**2)*mp.exp(-2*ALPHA_INV)
target_wege = [rho_L/(96*mp.pi), 1/(4*32*mp.pi**3*mp.exp(2*ALPHA_INV))]
d_spec = [v for q in ("(d) z112",) for v in [z["G/R_Vak^2"] for z in TAB if z["Quelle"] == q]]
check(f"[E] Umrechnungen: (a) Mit w c_1 lambda = 1/(16 pi G), c_1 = {c1_z118:.0f} (z118_ergebnis.json) und lambda = "
      f"4/R_Vak^2 ist G = {G_mm} und b = l^2/(64 pi G) = {b_mm}; z118s eigene Fassung (2 c_3/4) v^2 entspricht dann "
      f"v_geo R_Vak = {vR}. (b') rho_Lambda/Mbar^4 = (3/4 pi^2) e^(-2/alpha) = {mp.nstr(rho_L, 5)} gibt G Lambda = "
      f"rho_Lambda/(8 pi Mbar^4) = {mp.nstr(rho_L/(8*mp.pi), 5)} und mit Lambda = 12/R_Vak^2 G/R_Vak^2 = e^(-2/alpha)/(128 "
      f"pi^3) = {mp.nstr(TARGET, 6)} = 10^{float(mp.log10(TARGET)):.2f} (zwei Wege, z98: l^2/G = 32 pi^3 e^(2/alpha)). (d) "
      f"z112: 10^{float(mp.log10(min(d_spec))):.2f} bis 10^{float(mp.log10(max(d_spec))):.2f} (dort: -2,00 bis 0,80)",
      G_mm == Rv**2/(4096*sp.pi*ws) and b_mm == 16*ws and vR == 32 and all(abs(t/TARGET - 1) < mp.mpf(10)**-30 for t in target_wege)
      and abs(float(mp.log10(min(d_spec))) + 2.0023) < 1e-3 and abs(float(mp.log10(max(d_spec))) - 0.80) < 5e-3
      and abs(1/(18*ZETA3*C3) - mp.mpf("1.1615618")) < 1e-6)

werte = [(z["Quelle"][:3], z["Lesart"], z["G/R_Vak^2"]) for z in TAB if z["G/R_Vak^2"] is not None]
paare = [(abs(float(mp.log(a[2]/b[2]))), a, b) for a, b in itertools.combinations(werte, 2) if a[0] != b[0]]
paare.sort(key=lambda p: (p[0], p[1][1], p[2][1]))
naechst, zweitn = paare[0], paare[1]
boltz = 4*mp.exp(-4*mp.pi)*(1 + 4*mp.pi)                                 # gemeinsamer Boltzmann-Leitterm der beiden Zeilen
S_face_kms = mp.pi/face_exakt(2*mp.pi)
boltz_abw = [float(abs(S_face_kms/boltz - 1)), float(abs(4*S_leiter["4 pi"]/boltz - 1))]
n_5proz = sum(p_[0] < np.log(1.05) for p_ in paare)
endl = [float(mp.log10(z[2])) for z in werte if z[0] != "(b'"]
jac_naht_dS = mp.pi/(4*(sQ[0] + sQ[1]))
jac_naht_fl = mp.pi/(sQ[0] + sQ[1])
fak_target = {name: mp.pi/(4*s0)/TARGET for name, s0, _ in S0}
S_dS_ziel = mp.pi/(4*TARGET)
N_naht = 256*mp.pi/(sQ[0] + sQ[1])
N_E8 = 256*mp.pi/S_E8
ZAHLEN["vergleich"] = {
    "nächstes Paar verschiedener Quellen": [naechst[1][1], naechst[2][1], float(np.round(np.exp(naechst[0]), 8))],
    "zweitnächstes Paar": [zweitn[1][1], zweitn[2][1], float(np.round(np.exp(zweitn[0]), 6))],
    "Paare verschiedener Quellen, davon innerhalb 5 %": [len(paare), int(n_5proz)],
    "nächstes Paar: Abweichung beider Entropien vom Boltzmann-Leitterm 4 e^(-4 pi)(1 + 4 pi)": boltz_abw,
    "Spanne der nativen Zeilen (log10)": [min(endl), max(endl)],
    "Jacobson (Naht, dS) / (a2, w = 2 c_3)": mp.nstr(jac_naht_dS*1024, 12),
    "= 256 pi/s_0 (Naht), nächste ganze Zahl, Abstand": [mp.nstr(N_naht, 12), int(mp.nint(N_naht)), mp.nstr(abs(N_naht - mp.nint(N_naht)), 6)],
    "256 pi/ln 12": mp.nstr(N_E8, 12),
    "Jacobson (dS) / Ledgerziel (b') je s_0 (log10)": {k: float(mp.log10(v)) for k, v in fak_target.items()},
    "S_dS des Ledgerziels = 32 pi^4 e^(2/alpha)": mp.nstr(S_dS_ziel, 8),
    "Nähte je Horizont für das Ledgerziel (s_0 = ln 2/2)": mp.nstr(S_dS_ziel/(sQ[0] + sQ[1]), 8),
    "Jacobson (flach) / Ledgerziel (Naht, log10)": float(mp.log10(jac_naht_fl/TARGET))}
check(f"[E] Vergleich: Keine zwei Normierungen verschiedener Herkunft stimmen exakt überein ({len(paare)} Paare). Das "
      f"nächste Paar ('{naechst[1][1]}' gegen '{naechst[2][1]}') unterscheidet sich um den Faktor "
      f"{np.exp(naechst[0]):.7f}; das ist keine Bestätigung, sondern derselbe Boltzmann-Leitterm: Beide Entropien sind bis "
      f"auf {max(boltz_abw):.0e} gleich 4 e^(-4 pi)(1 + 4 pi) (tiefste Stufe bei beta = 4 pi; die 4 ist einmal die "
      f"Vielfachheit der tiefsten Dirac-Stufe auf S^2, einmal das Flächenverhältnis flach/dS), beides Lesarten mit "
      f"G/R_Vak^2 ~ 1,7e4. Das zweitnächste Paar ('{zweitn[1][1]}' gegen '{zweitn[2][1]}') liegt bei "
      f"{np.exp(zweitn[0]):.4f}; {n_5proz} von {len(paare)} Paaren liegen innerhalb 5 %, wie bei so vielen O(1)-Zahlen zu "
      f"erwarten (Zahlennähe ist kein Befund, z86). Die nativen Zeilen streuen über 10^{min(endl):.2f} bis "
      f"10^{max(endl):.2f}. Jacobson (Naht, dS) gegen z118 mit w = 2 c_3: Faktor "
      f"{mp.nstr(jac_naht_dS*1024, 8)} = 256 pi/s_0 (E8-Naht: {mp.nstr(N_E8, 7)}); das ist keine ganze Zahl (Abstand "
      f"{mp.nstr(abs(N_naht - mp.nint(N_naht)), 3)} zur nächsten; pi/ln 2 ist irrational, weil e^pi transzendent ist), also "
      f"macht keine ganze Zahl von Nähten aus 2 c_3 die native Entropiedichte. Gegen das Ledgerziel (b'): Jacobson (Naht, "
      f"dS) liegt um 10^{float(mp.log10(fak_target['Naht (beide Quarter)'])):.2f} darüber (flach: "
      f"10^{float(mp.log10(jac_naht_fl/TARGET)):.2f}; je nach s_0 10^{min(float(mp.log10(v)) for k, v in fak_target.items() if 'Leiter' not in k):.2f} "
      f"bis 10^{max(float(mp.log10(v)) for k, v in fak_target.items() if 'Leiter' not in k):.2f}); der Faktor ist genau "
      f"S_dS(Ziel)/s_0 = 32 pi^4 e^(2/alpha)/s_0 = {mp.nstr(S_dS_ziel/(sQ[0] + sQ[1]), 5)}, die Zahl der Nähte je Horizont, "
      f"die das Ziel verlangte",
      naechst[0] > 1e-9 and max(boltz_abw) < 1e-4 and zweitn[0] > 0.01
      and abs(N_naht - mp.nint(N_naht)) > mp.mpf("0.01") and max(endl) - min(endl) > 5
      and abs(fak_target["Naht (beide Quarter)"]/(32*mp.pi**4*mp.exp(2*ALPHA_INV)/(sQ[0] + sQ[1])) - 1) < mp.mpf(10)**-30
      and float(mp.log10(fak_target["Naht (beide Quarter)"])) > 120)
befund("C/E negativ: Unter Jacobsons S = A/(4G) mit s_0 = Entropie des Nahtzustands (Prämisse P_S) ist G/R_Vak^2 = "
       "pi/s_0 (flach/ESU) bzw. pi/(4 s_0) (dS): 9,06 bzw. 2,27 für die Naht, 1,26 bzw. 0,316 für die E8-Naht, 0,378 "
       "bzw. 0,094 für 48 Quarter; im nativen AdS-Rahmen ist die Fläche unendlich und die Lesart leer. Keine zwei der "
       "fünf Normierungen stimmen exakt überein (die einzige enge Nähe, Flächenzustand mit Kragenweite 2 pi R_Vak gegen "
       "Leiterfeld im dS-Rahmen, ist derselbe Boltzmann-Leitterm zweier Lesarten). Das von Hand gesetzte 2 c_3 aus z118 ist in v_geo-Einheiten nur die "
       "Bekenstein-Dichte eta = 1/4 (Tautologie von G = v_geo^-2); als reines Gewicht auf der MM-Spur verlangte es "
       "S_dS = 256 pi, um den Faktor 256 pi/s_0 (irrational, mit konventionsabhängigem 256) über jeder nativen Entropie. "
       "Das Ledgerziel (b') ist eine dS-Relation und liegt 122 bis 123 Größenordnungen unter dem Jacobson-Wert; der "
       "Faktor 32 pi^4 e^(2/alpha)/s_0 hat nur die Bedeutung einer Zählforderung (Nähte je Horizont), keinen nativen "
       "Mechanismus (z112)", True)

# =================================================================== Teil 4
say("\nTeil 4: Einheit oder Parameter")
dimvek = {"H6-Gravitationssektor {G, R_Vak}": sp.Matrix([[2, 1]]),
          "mit Lambda {G, R_Vak, Lambda}, Lambda R_Vak^2 = -12 nativ": sp.Matrix([[2, 1, -2]]),
          "Ledger {G, v_geo}": sp.Matrix([[2, -1]])}
invarianten = {k: M.shape[1] - M.rank() for k, M in dimvek.items()}
kern = dimvek["H6-Gravitationssektor {G, R_Vak}"].nullspace()[0]
theta = np.log(1.7)
gD_ = sla.expm(theta*MD)
E_dil = gD_ @ ME @ np.linalg.inv(gD_)
vak_dil = sorted(round(abs(p[2][0].imag), 8) for p in eigenebenen(E_dil))
unit_dil = float(np.abs(gD_.conj().T @ G4 @ gD_ - G4).max())
ZAHLEN["no_unit"] = {"Invarianten (Größen minus Rang)": invarianten, "Kern {G, R_Vak}": [str(x) for x in kern],
                     "Im x_Vak nach Ad(e^(ln 1,7 D))": vak_dil}
check(f"[E] No-Unit-Zählung (v153: ein dimensionsloser Compiler liefert keine absolute Skala; genau eine Einheit ist "
      f"irreduzibel): Zahl der skaleninvarianten Kombinationen = Zahl der dimensionsbehafteten Größen minus Rang der "
      f"Dimensionsmatrix: {invarianten}. Im H6-Gravitationssektor bleibt genau eine reine Zahl, G^{kern[0]} R_Vak^{kern[1]} "
      f"(G/R_Vak^2); Lambda R_Vak^2 = -12 ist nativ fest (z98, z118) und fügt keine zweite hinzu (dann G Lambda = -12 "
      f"G/R_Vak^2). R_Vak ist selbst die Einheit der Kinematik: die native Dilatation Ad(e^(theta D)) in SU(2,2) "
      f"({unit_dil:.0e}) verschiebt die Vakuumpunkte nach x^0 = +-i x {vak_dil[0]} (theta = ln 1,7) und lässt alle "
      f"Relationen stehen. Der No-Unit-Satz verbraucht also R_Vak (oder G), nicht beide",
      invarianten["H6-Gravitationssektor {G, R_Vak}"] == 1 and invarianten["Ledger {G, v_geo}"] == 1
      and invarianten["mit Lambda {G, R_Vak, Lambda}, Lambda R_Vak^2 = -12 nativ"] == 2
      and [sp.nsimplify(x) for x in kern] == [sp.Rational(-1, 2), 1] and vak_dil == [1.7, 1.7] and unit_dil < 1e-12)

LP = {name: mp.sqrt(s0/mp.pi) for name, s0, _ in S0}                      # R_Vak/l_P = sqrt(R^2/G), flach/ESU
LP_dS = {name: mp.sqrt(4*s0/mp.pi) for name, s0, _ in S0}
GL_dS = {name: 3*mp.pi/s0 for name, s0, _ in S0}                          # G Lambda = 12 G/R^2, dS
GL_beob = 12*TARGET
naht_namen = [name for name, _, _ in S0 if "Leiter" not in name]
GL_spanne = [min(float(GL_dS[k]) for k in naht_namen), max(float(GL_dS[k]) for k in naht_namen)]     # dS-Lesart: 3 pi/s_0
fenster = max(float(v) for v in list(LP.values()) + list(LP_dS.values()))
ZAHLEN["bedeutung"] = {"R_Vak/l_P (flach/ESU)": {k: mp.nstr(v, 6) for k, v in LP.items()},
                       "R_Vak/l_P (dS)": {k: mp.nstr(v, 6) for k, v in LP_dS.items()},
                       "G Lambda (dS-Lesart)": {k: mp.nstr(v, 6) for k, v in GL_dS.items()},
                       "G Lambda des Ledgerziels (dS)": mp.nstr(GL_beob, 6), "R_Vak/l_P des Ledgerziels": mp.nstr(1/mp.sqrt(TARGET), 6)}
check(f"[E bedingt auf P_S] Physikalische Bedeutung des Jacobson-Werts: R_Vak/l_P = sqrt(R_Vak^2/G) liegt für alle "
      f"nativen s_0 zwischen {min(float(LP[k]) for k in naht_namen):.3g} und {fenster:.3g} (Leiterfeld-Lesart: bis hinab zu "
      f"{min(float(v) for v in LP.values()):.1e}; Naht: "
      f"{mp.nstr(LP['Naht (beide Quarter)'], 4)} flach, {mp.nstr(LP_dS['Naht (beide Quarter)'], 4)} dS; 48 Quarter: "
      f"{mp.nstr(LP['drei Familien: 48 Quarter ohne GSO'], 4)}, {mp.nstr(LP_dS['drei Familien: 48 Quarter ohne GSO'], 4)}): "
      f"ein planckgroßes Vakuum, |Lambda| G = 12 G/R_Vak^2 = 3 pi/s_0 zwischen {GL_spanne[0]:.3g} und {GL_spanne[1]:.3g} (dS-Lesart; Naht: "
      f"{mp.nstr(GL_dS['Naht (beide Quarter)'], 4)}) gegen {mp.nstr(GL_beob, 3)} im Ledgerziel (R_Vak/l_P = "
      f"{mp.nstr(1/mp.sqrt(TARGET), 3)}). Zugleich ist Jacobsons Bereich l_P << Ballradius << Krümmungsradius dann leer "
      f"(R_Vak/l_P <= {fenster:.2f}): Der bedingte Wert untergräbt die Herleitung, aus der er stammt",
      fenster < 4 and min(float(v) for v in LP.values()) > 1e-6 and float(GL_dS["Naht (beide Quarter)"]) > 1
      and float(mp.log10(GL_dS["Naht (beide Quarter)"]/GL_beob)) > 120)
befund("E/C: Entscheidung (ii). P_G ist keine bloße Einheit: Der No-Unit-Satz verbraucht genau eine Länge, H6 hat aber "
       "zwei dimensionsbehaftete Gravitationsgrößen (R_Vak aus dem Zustand, G aus der Wirkung), also bleibt genau eine "
       "reine Zahl N_G = G/R_Vak^2 (gleichwertig G Lambda = -12 N_G) nativ frei. (i) gälte nur, wenn R_Vak keine "
       "physikalische Länge wäre (dann ist aber Lambda = -12/R_Vak^2 aus z98/z118 keine Aussage über die kosmologische "
       "Konstante); (iii) gälte nur unter der Prämisse P_S, mit N_G = pi/s_0 bzw. pi/(4 s_0) von der Ordnung 0,1 bis 20. "
       "Minimale verbleibende Prämisse: P_eta, eine native, endliche Entropie je Fläche der Bifurkationssphäre in einem "
       "Rahmen endlicher Fläche; sie ist eine Zahl (eta R_Vak^2) und gleichwertig zu N_G = 1/(4 eta R_Vak^2)", True)
befund("L: Bedeutung. Unter P_S wäre R_Vak von der Ordnung der Planck-Länge und Lambda planckgroß und nativ negativ (AdS). "
       "Das widerspricht der Beobachtung (G Lambda ~ +2,8e-122) nur, wenn R_Vak als kosmologischer Krümmungsradius gelesen "
       "wird; als mikroskopische Zellgröße gelesen gibt es keinen direkten Widerspruch, aber dann fehlt die zweite Länge "
       "(Faktor e^(-1/alpha) in der Länge, z112) samt Vorzeichen, und Lambda = -12/R_Vak^2 ist nicht die beobachtete "
       "kosmologische Konstante", True)

# =================================================================== Teil 5
say("\nTeil 5: Prämissen der Jacobson-Herleitung")
PRAEMISSEN = [
    ("J1 geometrischer Modulfluss K = 2 pi x Boost des Doppelkegels",
     "bedingt (P_Netz) im Vakuumsektor; fehlend für den Nahtzustand",
     "Q^- und sein Kegel nativ [E]; Hislop-Longo nur nach dem freien Netz (z104), 1D als Resultat (z72, hier 1.5); "
     "Nahtgenerator 4 pi E elliptisch bzw. zentral, nur Wick-konjugiert (hier 1.2, 1.3)",
     ["experiments/exploration/2026-11-23-phasenraum-f5-holocode/z104_ds_modular.py",
      "experiments/exploration/2026-11-23-phasenraum-f5-holocode/z72_zeit_translation.py"]),
    ("J2 Unruh-/KMS-Normierung 2 pi",
     "nativ als Gruppenaussage; bedingt als Zustandsaussage",
     "2 pi = imaginäre Vektorperiode von Q^- = Wick-Bild der Uhrdrehung; beta = 4 pi der Naht ist ein anderes Datum "
     "(hier 1.4; z104)",
     ["experiments/exploration/2026-11-23-phasenraum-f5-holocode/z104_ds_modular.py"]),
    ("J3 endliche Flächendichte eta (S_UV = eta A, UV-endlich, universell)",
     "fehlend (gescopt negativ)",
     "Original|Partner ist Typ I, thermisch, ohne Fläche (z86); räumlicher Schnitt divergiert bzw. ist reglerabhängig "
     "(z90, hier 3.3); der Flächenzustand hat endliches eta, aber freie Weite l (check_6_gravity)",
     ["experiments/exploration/2026-11-23-phasenraum-f5-holocode/z86_verschraenkung_flaeche.py",
      "experiments/exploration/2026-11-23-phasenraum-f5-holocode/z90_raeumlicher_schnitt.py",
      "experiments/theory-contracts/kernel-followups-20261120/check_6_gravity.py"]),
    ("J4 erster Hauptsatz delta S = delta<K> und maximale Entropie (zweite Ordnung)",
     "nativ für K = 4 pi L_0 (global); nicht für die Kappe",
     "hier Teil 2; z86 (gleiche Uhr je Sektor)",
     ["experiments/exploration/2026-11-23-phasenraum-f5-holocode/z86_verschraenkung_flaeche.py"]),
    ("J5 lokaler Spannungstensor, delta<K> = 2 pi int xi.T",
     "bedingt (P_Netz)",
     "xi_Q und das CHM-Gewicht (1 - r^2)/2 sind nativ (hier 1.1); ein lokales T_ab gibt es nur im freien Leiterfeld "
     "nach dem Netz (z83, z90); v358 führt CHM/BW als ASSUMPTION",
     ["experiments/exploration/2026-11-23-phasenraum-f5-holocode/z83_radiale_quantisierung.py",
      "verification/v358_grav_entropy_equilibrium.py"]),
    ("J6 konforme Materie (sonst Variante bei festem Volumen)",
     "nativ in der Kinematik, bedingt als Feld",
     "Leiterfeld = konform gekoppelter Skalar xi = 1/6, Delta = 1 (z83, z104); die Variante bei festem Volumen braucht "
     "eine dynamische Metrik (v359: C)",
     ["experiments/exploration/2026-11-23-phasenraum-f5-holocode/z83_radiale_quantisierung.py",
      "verification/v359_grav_nonlinear_einstein.py"]),
    ("J7 kleine Bälle an jedem Punkt, alle Zeitrichtungen (l_P << Radius << Krümmungsradius)",
     "Kovarianz nativ; Skalentrennung fehlend",
     "alle Doppelkegel sind in SU(2,2) konjugiert (hier 1.2), aber es gibt nur eine Länge; unter P_S ist R_Vak/l_P <= 3,3 "
     "(hier Teil 4)",
     ["experiments/exploration/2026-11-23-phasenraum-f5-holocode/z112_hodges_gravitation.py"]),
    ("J8 dynamische Metrik (Flächenvariation delta A aus der Krümmung)",
     "fehlend",
     "nativ nur das MM-Vakuum F = 0, flaches Vierbein (z98); Wirkung und Gewicht sind Prämissen (z118: P_R, 2 c_3)",
     ["experiments/exploration/2026-11-23-phasenraum-f5-holocode/z98_kopplung_mm.py",
      "experiments/exploration/2026-11-23-phasenraum-f5-holocode/z118_indexspur_wirkung.py"]),
    ("J9 lokales 3+1-Netz (Lokalität, Reeh-Schlieder, Haag-Dualität)",
     "fehlend (Import; DIMENSION.UPLIFT.FIREWALL.01, GRAV.NONCIRCULAR.01)",
     "su(2,2) ist Kinematik einer Darstellung; ein wechselwirkendes 4D-Netz ist nicht rekonstruiert (z83, z86)",
     ["experiments/exploration/2026-11-23-phasenraum-f5-holocode/z86_verschraenkung_flaeche.py"]),
    ("J10 Rahmen endlicher Horizontfläche",
     "fehlend im nativen AdS-Rahmen; bedingt (P_cap) für dS",
     "Fläche 4 pi (flach/ESU), pi (dS), unendlich (AdS) in R_Vak^2 (hier 3.2); positive Zeit <=> Lambda <= 0 (z104)",
     ["experiments/exploration/2026-11-23-phasenraum-f5-holocode/z104_ds_modular.py"]),
]
for name, status, beleg, _ in PRAEMISSEN:
    say(f"   {name}\n      Status: {status}\n      Beleg: {beleg}")
dateien = sorted({d for p in PRAEMISSEN for d in p[3]})
fehlend_dat = [d for d in dateien if not (REPO / d).is_file()]
n_nativ = sum(p[1].startswith("nativ") for p in PRAEMISSEN)
n_fehl = sum(p[1].startswith("fehlend") for p in PRAEMISSEN)
ZAHLEN["praemissen"] = [{"Prämisse": p[0], "Status": p[1], "Beleg": p[2], "Dateien": p[3]} for p in PRAEMISSEN]
check(f"[E] Prämissenliste: {len(PRAEMISSEN)} Prämissen, davon {n_nativ} (teilweise) nativ, {n_fehl} fehlend, "
      f"{len(PRAEMISSEN) - n_nativ - n_fehl} bedingt oder gemischt; alle {len(dateien)} Belegdateien liegen im Repository "
      f"(fehlend: {fehlend_dat})", not fehlend_dat and len(PRAEMISSEN) == 10 and n_fehl >= 3)
befund("C/O: Von Jacobsons Prämissen sind nativ nur die kinematischen (Kegel, Kovarianz, 2 pi als Gruppenaussage, "
       "konformes Gewicht, erster Hauptsatz des globalen Generators). Es fehlen die drei, an denen der Wert von G hängt: "
       "eine endliche Flächendichte (J3), eine dynamische Metrik (J8) und das lokale 3+1-Netz (J9); dazu der Rahmen "
       "endlicher Fläche (J10) und die Skalentrennung (J7). Die Ledgerkette v358/v359 bleibt ASSUMPTION im Sinn von "
       "GRAV.NONCIRCULAR.01; H6 liefert das fehlende Stück QGEO.MODHAM.01 nicht", True)

# =================================================================== Bilanz und Urteil
say("\nBilanz")
bilanz = [
    "[E] radiale sl(2,R) in su(2,2); Kegel von Q^-: Zentrum tau-Fixpunkt, Radius R_Vak, kappa = 1, K_cap = 2 pi Q^- "
    "(Hislop-Longo-Normierung), CHM-Gewicht (1 - r^2)/2",
    "[E, negativ] K_Uhr = 4 pi E ist in SU(2,2) (und SL(4,C)) zu keinem reellen Vielfachen von Q^- konjugiert; nur Wick "
    "W = e^(i pi D/2): W K_Uhr W^-1 = 2 i K_cap",
    "[E] Faktor 2: Gruppenaussage ja (4 pi Twistor, 2 pi Vektor, 8 pi Quarter; (-1)^(2h)), Erzeugerwahl ja (L_0^NS = 2 "
    "L_0), Zustandsaussage nein",
    f"[E numerisch, 1D] Peschel: Vakuum geometrisch (lam -> 1), Nahtzustand nicht (Residuum {n256['Residuum (2 pi G_v)']:.2f})",
    "[E] erster Hauptsatz und Fisher-Positivität nativ für K = 4 pi L_0; beta = 2 pi, E allein, 2 pi H_1 scheitern",
    "[E] s_0: ln 2/4 je Quarter, ln 2/2 je Naht, ln 12 E8-Naht; räumlicher Schnitt divergiert",
    "[E] Fläche der Bifurkationssphäre: 4 pi (flach/ESU), pi (dS), unendlich (AdS, nativ) in R_Vak^2",
    "[C] Jacobson unter P_S: G/R_Vak^2 = pi/s_0 bzw. pi/(4 s_0); keine zwei Normierungen stimmen überein",
    "[E/C] P_G: Entscheidung (ii), eine reine Zahl N_G = G/R_Vak^2 bleibt frei; Minimalprämisse P_eta",
    "[C/O] Jacobson-Prämissen: J3, J8, J9 fehlen, J10 und J7 im nativen Rahmen nicht erfüllt",
]
for z in bilanz:
    say("   " + z)
urteil = ("P_G ist nativ ein freier Parameter, genauer genau eine reine Zahl N_G = G/R_Vak^2 (gleichwertig G Lambda = -12 "
          "N_G); das Verschränkungsgleichgewicht legt sie nicht nativ fest. Nichts ist befördert. (1) Firewall: In der "
          "nativen su(2,2) ist der Doppelkegelgenerator vorhanden und richtig normiert (Q^- mit Zentrum tau-Fixpunkt, "
          "Radius R_Vak, Oberflächengravitation 1, K_cap = 2 pi Q^-, CHM-Gewicht (1 - r^2)/2), und alle Doppelkegel sind "
          "in SU(2,2) konjugiert. Der Modulgenerator des Nahtzustands ist aber K_Uhr = 4 pi E (einseitig, elliptisch) "
          "bzw. zentral (Tomita auf der Standardform) und zu keinem reellen Vielfachen eines Doppelkegelgenerators "
          "konjugiert, weder in SU(2,2) noch in SL(4,C); allein die Wick-Abbildung W = e^(i pi D/2), kein Gruppenelement, "
          "gibt W K_Uhr W^-1 = 2 i K_cap. Der Faktor 2 ist als Gruppenaussage der Spinor-Halbwinkel (Periode 4 pi auf "
          "Twistoren, 2 pi auf Vektoren, 8 pi auf Quartermoden) und als Erzeugerwahl L_0^NS = 2 L_0 (tau = i), als "
          "Zustandsaussage ein echter Unterschied: beta = 4 pi ist ein KMS-Datum des elliptischen Erzeugers ohne "
          "imaginäre Periode. Als Resultat aus dem Zustand ist der geometrische Modulfluss nur im Grundzustand der Uhr "
          "belegt (Peschel auf der halben Naht, 1D); im Nahtzustand ist er nicht geometrisch. Import bleibt die "
          "Lokalität eines 4D-Netzes; QGEO.MODHAM.01 ist nicht geschlossen. (2) Erster Hauptsatz: delta S = delta<K> gilt "
          "nativ und exakt in erster Ordnung für Kovarianz-, Quartik- und Möbiusstörungen, mit positiver zweiter Ordnung "
          "(relative Entropie = Fisher-Information), aber für den globalen Generator 4 pi L_0; beta = 2 pi scheitert um "
          "den Faktor 2, ein hyperbolischer Generator ganz. (3) Wert: s_0 = ln 2/4 je Quarter, ln 2/2 je Naht, ln 12 für "
          "die E8-Naht. Die Fläche der Bifurkationssphäre hängt am konformen Rahmen: 4 pi R_Vak^2 (flach, Einstein-"
          "Universum), pi R_Vak^2 (dS, braucht P_cap), unendlich im nativen AdS-Rahmen. Unter der Prämisse P_S (s_0 ist "
          "die Entropie der Bifurkationssphäre) wäre G/R_Vak^2 = pi/s_0 = 9,06 bzw. pi/(4 s_0) = 2,27 für die Naht (E8-"
          "Naht 1,26 bzw. 0,316). Keine zwei der vorhandenen Normierungen stimmen überein: z118 als reines Gewicht gibt "
          "1/(4096 pi w) (w = 2 c_3: 1/1024), z112 10^-2,0 bis 10^0,8, der Flächenzustand hängt an der freien Weite, der "
          "Ledgerwert 8 pi G = 1/c_3 definiert nur die Einheit, und das Ledgerziel e^(-2/alpha)/(128 pi^3) = 10^-122,63 "
          "(Status Numerical, dS-Relation) liegt 122 bis 123 Größenordnungen darunter; der Faktor 32 pi^4 e^(2/alpha)/s_0 "
          "ist eine Zählforderung (Nähte je Horizont) ohne nativen Mechanismus. Das 2 c_3 aus z118 ist nicht die native "
          "Entropiedichte: Es verlangte S_dS = 256 pi, um den irrationalen Faktor 256 pi/s_0 über jeder nativen Entropie. "
          "(4) Einheit oder Parameter: Der No-Unit-Satz verbraucht eine Länge, H6 hat zwei (R_Vak, sqrt G); es bleibt die "
          "Zahl N_G. Entscheidung (ii); (iii) nur unter P_S, und dann mit planckgroßem Vakuum (R_Vak/l_P <= 3,3), "
          "negativem Lambda im nativen Rahmen und leerem Gültigkeitsbereich der Jacobson-Herleitung. Minimalprämisse: "
          "P_eta, eine native endliche Entropie je Fläche der Bifurkationssphäre in einem Rahmen endlicher Fläche. (5) Von "
          "Jacobsons Prämissen fehlen die endliche Flächendichte, die dynamische Metrik und das lokale 3+1-Netz; die "
          "Ledgerkette v358/v359 bleibt ASSUMPTION nach GRAV.NONCIRCULAR.01.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z136_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil, "bilanz": bilanz}, indent=1,
    ensure_ascii=False, default=str) + "\n")
