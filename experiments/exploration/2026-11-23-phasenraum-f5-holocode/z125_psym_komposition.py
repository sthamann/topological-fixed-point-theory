"""Z125: P_sym aus kohärenter Komposition? Identifikation, Leckfreiheit allgemein, Nativität, Energieroute zu S(n).

Folgeprüfung zu z120 Teil 3b (Leckfreiheit erzwingt die Vierteldrehung). Leitlinie 'Rekonstruktion statt Auswahl'; H6 und die
gelösten Bausteine bleiben unberührt, nichts wird befördert. Statuskürzel: E exakt, C bedingt, L Lesart, O offen.
Bestand: z120 (Transfergatter U = B^+ e^(-i phi H_Q) B, leckfrei nur bei theta in {0, pi/2, pi}, Energie-Eindeutigkeit über
Transzendenz), z94 (P_sym <=> Vierteldrehung A -> B -> -A innerlich <=> theta = pi/2 zwischen zwei Binäroperationen in der
irreduziblen 2D-Darstellung; Tausch immer innerlich; P_clock <=> P_tn), z80/z87 (Saat rho_1 + Spin-Z4 + CAR, Cl_2, tau = i),
z107 (P_tn aus A*, A* impliziert P_sym nicht: theta = pi/3 spannungsfrei mit j = 3905.8; S(n) ab n = 3 arithmetisch, Kern
2n - 3), z93/z110 (Ladungssektor ist Teil des Zustands), z105 (native Zweinahtterme je Naht gerade), z106 (gradierte und
ungradierte Komposition haben dieselbe gerade Algebra), z108 (erhaltene Einzeldynamik heißt Produktzustand).

Teil 1  Identifikation: Auf dem Einanregungssektor zweier Nähte ist A = -iR (Nahtparität, Uhr bis auf i) und B = B_theta^+ A
        B_theta; {A, B} = 2 cos theta. Leckformel der Zweianregungskomposition |1/2 {A,B}| |1/2 [A,B]| |v_P - v_D|; gleichwertige
        Formen ([U, (-1)^F_A] = 0, [R, B^+ R B] = 0 auf N = 2).
Teil 2  Verallgemeinerung: Paarsatz für ein beliebiges diagonales Phasengesetz und verschiedene Transferwinkel (symbolisch),
        Energieleiter mit Amplituden, Vorspannung, beliebiges phi, zwei Arten mit allgemeinen Ladungen (Austauschphase
        phi q_k q_l), drei Nähte (paarweise; zyklisch; F_3; zufällig; Kriterium im adjungierten Raum).
Teil 3  Nativität: 'eine Anregung je Naht' als Paritäts-, Uhr- und Energiesektor; Klassifikation aller Z4-Gradierungen
        i^(a Q_A + b Q_B) als Superauswahl; Artenzahl drei; drei Nähte unter Nahtparität; Tomita-Paaruhr Delta^(-i/2);
        Einzel- gegen Doppeltransfer; Verträglichkeit mit z107 und z94.
Teil 4  Energieroute zu S(n): Inhalt des Nahttorus der Kandidaten (Modul- gegen geometrische Energie), S-Defekt (A*(a) auf
        den Kandidaten), Restklassenkern mit Energieidentität für die Fenster n = 2..5 (z107-Tupel), Ganzzahligkeit (reelle
        Relaxation mit einer Energiegleichung).
Teil 5  Bilanz, Prämissen, Urteil.
"""
import itertools
import json
import time
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path

import flint
import mpmath as mp
import numpy as np
import sympy as sp
from scipy.linalg import expm, schur
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import csr_matrix, lil_matrix, vstack
from scipy.stats import unitary_group

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
rng = np.random.default_rng(125)
mp.mp.dps = 50
T_START = time.time()
YV = mp.exp(-mp.pi)
yS = sp.symbols("y")


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val):
    BEFUNDE[name] = bool(val)
    say(f"   [Befund: {'ja' if val else 'nein'}] {name}")


def jw(n):
    I_, Zp, sm = np.eye(2), np.diag([1.0, -1.0]), np.array([[0.0, 1.0], [0.0, 0.0]])
    out = []
    for j in range(n):
        m_ = [Zp]*j + [sm] + [I_]*(n - j - 1)
        op = m_[0]
        for M in m_[1:]:
            op = np.kron(op, M)
        out.append(op)
    return out


def nahtmodell(nahte, arten, ladungen=None):
    """Moden (Naht i, Art a) in Jordan-Wigner-Reihenfolge naht-major. Rückgabe: Vernichter, Nahtzahlen n_i und Nahtladungen
    Q_i (Ladung je Art) als Diagonalen."""
    c = jw(nahte*arten)
    num = [np.diag(x.T @ x).real.round().astype(int) for x in c]
    lad = ladungen or [1]*arten
    n_naht = [sum(num[i*arten + a] for a in range(arten)) for i in range(nahte)]
    q_naht = [sum(lad[a]*num[i*arten + a] for a in range(arten)) for i in range(nahte)]
    return c, n_naht, q_naht


def transfer(c, nahte, arten, Wn, t=None):
    """Zweitquantisierter Transfer sum_(i,j,a) t_a Wn[i,j] c_(i a)^+ c_(j a) (artenerhaltend, Wn hermitesch auf den Nähten)."""
    H = np.zeros(c[0].shape, dtype=complex)
    for i in range(nahte):
        for j in range(nahte):
            if abs(Wn[i, j]) > 0:
                for a in range(arten):
                    H += (1 if t is None else t[a])*Wn[i, j]*(c[i*arten + a].T @ c[j*arten + a])
    return H


def log_unitaer(u):
    """h mit u = exp(-i h) (Schur-Zerlegung, normal)."""
    T, Z = schur(u, output="complex")
    return Z @ np.diag(-np.angle(np.diag(T))) @ Z.conj().T


def gatter_sektor(H, Vdiag, idx, th):
    """U = B^+ V B auf dem Teilchenzahlsektor idx, B = exp(-i th H)."""
    Hs = H[np.ix_(idx, idx)]
    B = expm(-1j*th*Hs)
    return B.conj().T @ np.diag(Vdiag[idx]) @ B


def leck_austausch(U, phys_lok):
    """Leckage ||(1 - P) U P|| und kohärenter Austausch (nichtdiagonaler Teil auf dem physikalischen Raum, Besetzungsbasis)."""
    rest = [k for k in range(U.shape[0]) if k not in set(phys_lok)]
    leck = np.linalg.norm(U[np.ix_(rest, phys_lok)], 2) if rest else 0.0
    Up = U[np.ix_(phys_lok, phys_lok)]
    return leck, np.linalg.norm(Up - np.diag(np.diag(Up)), 2), Up


def scan(nahte, arten, gesetz, ladungen=None, t=None, Wn=None, thetas=None, tol=1e-9):
    """Leckfreie Winkel und leckfreie Winkel mit kohärentem Austausch (eine Anregung je Naht)."""
    c, nN, qN = nahtmodell(nahte, arten, ladungen)
    tot = sum(nN)
    idx = np.where(tot == nahte)[0]
    phys_lok = [k for k, s_ in enumerate(idx) if all(n_[s_] == 1 for n_ in nN)]
    if Wn is None:
        Wn = np.zeros((nahte, nahte))
        Wn[0, 1] = Wn[1, 0] = 0.5
    H = transfer(c, nahte, arten, Wn, t)
    Vd = np.array([gesetz(np.array([q[s_] for q in qN]), np.array([n_[s_] for n_ in nN]), s_) for s_ in range(2**(nahte*arten))])
    lf, lfa = [], []
    for th in thetas:
        l_, a_, _ = leck_austausch(gatter_sektor(H, Vd, idx, th), phys_lok)
        if l_ < tol:
            lf.append(round(float(th/np.pi), 4))
            if a_ > 1e-6:
                lfa.append(round(float(th/np.pi), 4))
    return lf, lfa


def implementierbar(bild_A, bild_B, A, B):
    """wie z94: gibt es U != 0 mit U A = bild_A U, U B = bild_B U (Intertwinergleichungen)?"""
    zeilen = []
    for i in range(2):
        for j in range(2):
            U = np.zeros((2, 2), dtype=complex)
            U[i, j] = 1
            zeilen.append(np.concatenate([(U @ A - bild_A @ U).ravel(), (U @ B - bild_B @ U).ravel()]))
    sv = np.linalg.svd(np.array(zeilen).T, compute_uv=False)
    return bool(np.sum(sv < 1e-10) > 0)


def j_marken(th):
    """wie z94: j-Invariante der vier Marken {1, -1, e^(i th), -e^(i th)}."""
    z = [1, -1, np.exp(1j*th), -np.exp(1j*th)]
    lam = (z[0] - z[2])*(z[1] - z[3])/((z[0] - z[3])*(z[1] - z[2]))
    return float((256*(1 - lam + lam**2)**3/(lam**2*(1 - lam)**2)).real)


TH181 = np.linspace(0, np.pi, 181)

# =================================================================== Teil 1
say("Teil 1: Identifikation -- die Transferdrehung trägt die Binäroperationen von z94")
c4, nN4, qN4 = nahtmodell(2, 2)
QA, QB = qN4
NT = QA + QB
Jx = transfer(c4, 2, 2, np.array([[0, 0.5], [0.5, 0]]))
Rd = 1j**((QA - QB) % 4)
R = np.diag(Rd)
HQd = (QA**2 + QB**2)/2
PAd = (-1.0)**QA
i1, i2 = np.where(NT == 1)[0], np.where(NT == 2)[0]
I4 = np.eye(len(i1))

# 1.1 Einanregungssektor
A1 = (-1j*R)[np.ix_(i1, i1)]
binaer, abw_anti = True, 0.0
for th in np.linspace(0, np.pi, 13):
    Bt = expm(-1j*th*Jx[np.ix_(i1, i1)])
    B1 = Bt.conj().T @ A1 @ Bt
    binaer &= np.allclose(A1 @ A1, I4) and np.allclose(B1 @ B1, I4)
    abw_anti = max(abw_anti, np.linalg.norm(A1 @ B1 + B1 @ A1 - 2*np.cos(th)*I4))
ein = [k for k, s_ in enumerate(i1) if (c4[1].T @ c4[1])[s_, s_] == 0 and (c4[3].T @ c4[3])[s_, s_] == 0]  # nur Art 0
A2 = A1[np.ix_(ein, ein)]
J2 = Jx[np.ix_(i1, i1)][np.ix_(ein, ein)]
winkel = np.unique(np.round(np.append(np.linspace(0.05, np.pi - 0.05, 61), np.pi/2), 12))
impl = {"Tausch A<->B": [], "Vorzeichen A->-A": [], "Viertel A->B->-A": []}
for th in winkel:
    Bt = expm(-1j*th*J2)
    B2 = Bt.conj().T @ A2 @ Bt
    for name, (bA, bB) in (("Tausch A<->B", (B2, A2)), ("Vorzeichen A->-A", (-A2, B2)), ("Viertel A->B->-A", (B2, -A2))):
        if implementierbar(bA, bB, A2, B2):
            impl[name].append(round(float(th/np.pi), 4))
selbst_viertel = []
for th in TH181:
    Bt = expm(-1j*th*Jx[np.ix_(i1, i1)])
    B1 = Bt.conj().T @ A1 @ Bt
    if np.allclose(Bt.conj().T @ B1 @ Bt, -A1, atol=1e-10):
        selbst_viertel.append(round(float(th/np.pi), 4))
Bt = expm(-1j*np.pi/2*Jx[np.ix_(i1, i1)])
B1 = Bt.conj().T @ A1 @ Bt
Gam = -1j*A1 @ B1
gam_abw = min(np.linalg.norm(Gam - s_*2*Jx[np.ix_(i1, i1)]) for s_ in (1, -1))
jm = {"pi/2": j_marken(np.pi/2), "pi/3": j_marken(np.pi/3)}
ZAHLEN["identifikation"] = {"max |{A,B} - 2 cos theta|": float(abw_anti), "z94-Intertwiner (theta/pi)":
                            {k: (v if len(v) < 5 else f"{len(v)} von {len(winkel)}") for k, v in impl.items()},
                            "B_theta implementiert A -> B -> -A bei theta/pi": selbst_viertel,
                            "|Gamma -+ 2 J_x| bei pi/2": float(gam_abw), "j(Marken)": jm}
check(f"[E] Einanregungssektor zweier Nähte (z120-Modell, zwei Arten als Vielfachheit): A = -iR = Nahtparität (die Uhr R = "
      f"i^(Q_A - Q_B) bis auf i) und B = B_theta^+ A B_theta sind Binäroperationen ({binaer}) mit {{A, B}} = 2 cos theta "
      f"(13 Winkel, Abweichung {abw_anti:.0e}): der Transferwinkel ist der Winkel von z94. z94-Intertwinertest auf der "
      f"irreduziblen 2D-Darstellung (eine Art): Tausch {len(impl['Tausch A<->B'])} von {len(winkel)}, Vorzeichen "
      f"{impl['Vorzeichen A->-A']}, Viertel {impl['Viertel A->B->-A']}; B_theta selbst implementiert A -> B -> -A genau bei "
      f"theta/pi = {selbst_viertel}; Gamma = -iAB ist der Transfergenerator 2 J_x (Abweichung {gam_abw:.0e}): die Energieachse "
      f"von z94 ist das Hüpfen, der Modulfluss der Saat die Transferdrehung; Marken j = {jm['pi/2']:.1f} (pi/2) bzw. "
      f"{jm['pi/3']:.1f} (pi/3)",
      binaer and abw_anti < 1e-12 and len(impl["Tausch A<->B"]) == len(winkel) and impl["Vorzeichen A->-A"] == [0.5]
      and impl["Viertel A->B->-A"] == [0.5] and selbst_viertel == [0.5] and gam_abw < 1e-12
      and abs(jm["pi/2"] - 1728) < 1e-6 and abs(jm["pi/3"] - 3905.8) < 0.1)

# 1.2 Leckformel auf N = 2
phys2 = [k for k, s_ in enumerate(i2) if QA[s_] == 1 and QB[s_] == 1]
phis = [np.pi/5, np.pi/3, np.pi/2, 2*np.pi/3, 1.0, 2.0, np.pi]
abw_formel, abw_2th, abw_halb = 0.0, 0.0, 0.0
for phi in phis:
    Vd = np.exp(-1j*phi*HQd)
    for th in np.linspace(0, np.pi, 37):
        l_, _, _ = leck_austausch(gatter_sektor(Jx, Vd, i2, th), phys2)
        fak = abs(1 - np.exp(-1j*phi))
        abw_formel = max(abw_formel, abs(l_ - abs(np.cos(th))*abs(np.sin(th))*fak))
        abw_2th = max(abw_2th, abs(l_ - abs(np.cos(2*th))*abs(np.sin(2*th))*fak))
        abw_halb = max(abw_halb, abs(l_ - abs(np.cos(th/2))*abs(np.sin(th/2))*fak))
Vd_h = np.exp(-1j*np.pi/2*HQd)
lf_12 = [round(float(th/np.pi), 4) for th in TH181 if leck_austausch(gatter_sektor(Jx, Vd_h, i2, th), phys2)[0] < 1e-10]
leck_pi3 = leck_austausch(gatter_sektor(Jx, Vd_h, i2, np.pi/3), phys2)[0]
ZAHLEN["leckformel"] = {"max Abweichung |cos||sin||1 - e^-i phi|": abw_formel, "Kontrolle 2 theta": abw_2th,
                        "Kontrolle theta/2": abw_halb, "leckfrei (phi = pi/2) theta/pi": lf_12, "Leckage bei pi/3": leck_pi3}
check(f"[E] Leckformel der Zweianregungskomposition U = B_theta^+ e^(-i phi H_Q) B_theta (37 Winkel x 7 Phasen phi): "
      f"||(1 - P) U P|| = |1/2 {{A,B}}| |1/2 [A,B]| |v_P - v_D| = |cos theta| |sin theta| |1 - e^(-i phi)| (Abweichung "
      f"{abw_formel:.0e}); Kontrolle: mit 2 theta bzw. theta/2 statt theta weicht sie um {abw_2th:.3f} bzw. {abw_halb:.3f} ab "
      f"(der Winkel ist derselbe, kein Vielfaches). Leckfrei (phi = pi/2, 181 Winkel) genau bei theta/pi = {lf_12}: A und B "
      f"vertauschen (kein Transfer) oder antivertauschen (P_sym); bei pi/3 Leckage {leck_pi3:.6f}",
      abw_formel < 1e-12 and abw_2th > 0.1 and abw_halb > 0.1 and lf_12 == [0.0, 0.5, 1.0])

# 1.3 gleichwertige Formen
PA = np.diag(PAd)
form_PA, form_R, form_RR, anti_N1 = [], [], [], None
for th in TH181:
    B = expm(-1j*th*Jx)
    U = B.conj().T @ np.diag(Vd_h) @ B
    if np.linalg.norm(U @ PA - PA @ U) < 1e-9:
        form_PA.append(round(float(th/np.pi), 4))
    if np.linalg.norm(U @ R - R @ U) < 1e-9:
        form_R.append(round(float(th/np.pi), 4))
    RB = (B.conj().T @ R @ B)
    if np.linalg.norm((R @ RB - RB @ R)[np.ix_(i2, i2)]) < 1e-9:
        form_RR.append(round(float(th/np.pi), 4))
B = expm(-1j*np.pi/2*Jx)
RB = B.conj().T @ R @ B
anti_N1 = np.linalg.norm((R @ RB + RB @ R)[np.ix_(i1, i1)])
komm_N2 = np.linalg.norm((R @ RB - RB @ R)[np.ix_(i2, i2)])
check(f"[E] gleichwertige Formen auf dem ganzen Fockraum (16 Zustände, 181 Winkel): [U, (-1)^F_A] = 0 bei theta/pi = {form_PA}, "
      f"[U, R] = 0 bei {form_R}, [R, B^+ R B] = 0 auf N = 2 bei {form_RR}. Bei pi/2 antivertauschen R und B^+ R B auf N = 1 "
      f"(Rest {anti_N1:.0e}) und vertauschen auf N = 2 (Rest {komm_N2:.0e}): Spinor gegen Vektor der Transfer-SU(2), "
      "R^2 = (-1)^N; die Nahtparität ist die Halbdrehung um die Nahtachse",
      form_PA == form_R == form_RR == [0.0, 0.5, 1.0] and anti_N1 < 1e-12 and komm_N2 < 1e-12)
befund("E: Die Identifikation ist wörtlich. Der Einanregungssektor zweier Nähte ist die irreduzible 2D-Darstellung von z94; "
       "A ist die Nahtparität (-iR), B ihr Bild unter der kohärenten Übertragung, Gamma = -iAB der Transfergenerator, der "
       "Saat-Modulfluss die Transferdrehung. z120s Leckfreiheit ist die Aussage '|1/2 {A,B}| |1/2 [A,B]| = 0': A und B "
       "vertauschen oder antivertauschen; nichttrivial bleibt nur {A, B} = 0, also P_sym mit tau = i (j = 1728)", True)
say(f"   (Laufzeit Teil 1: {time.time() - T_START:.1f} s)")

# =================================================================== Teil 2
say("\nTeil 2: Verallgemeinerung der Leckfreiheit")
# 2.1 Paarsatz (zwei unterscheidbare Anregungen k, l; Basis |AA>, |AB>, |BA>, |BB>, erste Marke = Naht von k)
vs = sp.symbols("v_AA v_AB v_BA v_BB")


def b_sym(t):
    return sp.Matrix([[sp.cos(t/2), -sp.I*sp.sin(t/2)], [-sp.I*sp.sin(t/2), sp.cos(t/2)]])


def paar_system(tk, tl):
    W = sp.kronecker_product(b_sym(tk), b_sym(tl))
    U = W.H*sp.diag(*vs)*W
    leck = [U[i, j] for i in (0, 3) for j in (1, 2)]
    M = sp.Matrix([[sp.nsimplify(sp.simplify(sp.diff(L, v))) for v in vs] for L in leck])
    o = sp.Matrix([[sp.nsimplify(sp.simplify(sp.diff(U[2, 1], v))) for v in vs]])
    return M, o, U


faelle = [(sp.pi/2, sp.pi/2), (sp.pi/2, 3*sp.pi/2), (sp.pi/2, sp.pi/3), (sp.pi/3, sp.pi/3), (sp.pi/2, sp.pi),
          (sp.pi/3, 2*sp.pi/3), (sp.pi/4, sp.pi/4)]
paar_tab = {}
for tk, tl in faelle:
    M, o, U = paar_system(tk, tl)
    ns = M.nullspace()
    austausch = any(sp.simplify((o*n)[0]) != 0 for n in ns)
    paar_tab[f"({tk}, {tl})"] = {"Lösungsraum": len(ns), "kohärenter Austausch": austausch}
M, o, U = paar_system(sp.pi/2, sp.pi/2)
ns = M.nullspace()
symm = len(ns) == 2 and all(sp.simplify(n[0] - n[3]) == 0 and sp.simplify(n[1] - n[2]) == 0 for n in ns)
vP, vD = sp.symbols("v_P v_D")
Up = sp.simplify(U.subs({vs[0]: vD, vs[3]: vD, vs[1]: vP, vs[2]: vP})[1:3, 1:3])
Xs = sp.Matrix([[0, 1], [1, 0]])
gatter_ok = sp.simplify(Up - (vD*(sp.eye(2) + Xs)/2 + vP*(sp.eye(2) - Xs)/2)) == sp.zeros(2, 2)
# geschlossene Form numerisch an 20 Zufallspunkten
sx = np.array([[0, 1], [1, 0]], dtype=complex)
abw_geschl = 0.0
for _ in range(20):
    tk_, tl_ = rng.uniform(0, 2*np.pi, 2)
    v_ = np.exp(1j*rng.uniform(0, 2*np.pi, 4))
    W = np.kron(expm(-0.5j*tk_*sx), expm(-0.5j*tl_*sx))
    Un = W.conj().T @ np.diag(v_) @ W
    ck, sk, cl, sl = np.cos(tk_/2), np.sin(tk_/2), np.cos(tl_/2), np.sin(tl_/2)
    vAA, vAB, vBA, vBB = v_
    geschl = [-1j*(ck**2*(vAA - vAB) + sk**2*(vBA - vBB))*sl*cl, -1j*(cl**2*(vAA - vBA) + sl**2*(vAB - vBB))*sk*ck,
              1j*(sl**2*(vAA - vBA) + cl**2*(vAB - vBB))*sk*ck, 1j*(sk**2*(vAA - vAB) + ck**2*(vBA - vBB))*sl*cl]
    abw_geschl = max(abw_geschl, max(abs(Un[i, j] - g_) for (i, j), g_ in zip(((0, 1), (0, 2), (3, 1), (3, 2)), geschl)))
# Querprobe gegen den Fockraum (Arten k = 0, l = 1 als unterscheidbare Anregungen)
idx_q = [s_ for s_ in i2 if sum((c4[m].T @ c4[m])[s_, s_] for m in (0, 2)) == 1]          # je ein Teilchen jeder Art
pos = {s_: ("A" if (c4[0].T @ c4[0])[s_, s_] else "B") + ("A" if (c4[1].T @ c4[1])[s_, s_] else "B") for s_ in idx_q}
v_ = np.exp(1j*rng.uniform(0, 2*np.pi, 4))
vmap = dict(zip(("AA", "AB", "BA", "BB"), v_))
Vf = np.ones(16, dtype=complex)
for s_, p_ in pos.items():
    Vf[s_] = vmap[p_]
Bf = expm(-1j*1.1*Jx)
Uf = (Bf.conj().T @ np.diag(Vf) @ Bf)
W = np.kron(expm(-0.55j*sx), expm(-0.55j*sx))
Un = W.conj().T @ np.diag(v_) @ W
ordn = {"AA": 0, "AB": 1, "BA": 2, "BB": 3}
quer = max(abs(abs(Uf[s1, s2]) - abs(Un[ordn[pos[s1]], ordn[pos[s2]]])) for s1 in idx_q for s2 in idx_q)
ZAHLEN["paarsatz"] = {"Fälle (theta_k, theta_l)": paar_tab, "pi/2, pi/2: Lösungsraum nahtsymmetrisch": symm,
                      "Gatter v_D P_+ + v_P P_-": gatter_ok, "geschlossene Form Abweichung": abw_geschl,
                      "Querprobe Fockraum |U|": quer}
aus = [k for k, v in paar_tab.items() if v["kohärenter Austausch"]]
check(f"[E] Paarsatz (symbolisch, beliebiges diagonales Phasengesetz v_AA, v_AB, v_BA, v_BB, Transferwinkel theta_k, theta_l): "
      f"die vier Leckamplituden sind <AA|U|AB> = -i[c_k^2 (v_AA - v_AB) + s_k^2 (v_BA - v_BB)] s_l c_l usw. (c = cos theta/2, "
      f"s = sin theta/2; 20 Zufallspunkte, Abweichung {abw_geschl:.0e}; Fockraum-Querprobe {quer:.0e}). Kohärenter Austausch "
      f"bei Leckfreiheit nur in den Fällen {aus} von {len(faelle)}; bei (pi/2, pi/2) ist der Lösungsraum genau v_AA = v_BB, "
      f"v_AB = v_BA ({symm}) und das Gatter v_D P_+ + v_P P_- ({gatter_ok}). Satz (aus den geschlossenen Formen: für "
      "sin theta_k sin theta_l != 0 zwei 2x2-Systeme mit Determinanten cos theta_k und cos theta_l; für sin = 0 bleibt U "
      "diagonal): leckfreier nichttrivialer Austausch <=> beide Winkel = pi/2 mod pi und nahtsymmetrisches Gesetz mit "
      "v_D != v_P",
      abw_geschl < 1e-12 and quer < 1e-12 and symm and gatter_ok
      and aus == ["(pi/2, pi/2)", "(pi/2, 3*pi/2)"])

# 2.2 Leiter, Amplituden, Vorspannung, beliebiges phi
eps_l = [0.25, 1.25, 2.25]
leiter = {}
for phi in (np.pi/5, np.pi/2, 2.0, 2*np.pi):
    g_ = lambda q, n_, s_, phi=phi: np.exp(-1j*phi*np.sum(q**2)/2)
    leiter[f"phi = {phi:.4f}"] = scan(2, 3, g_, thetas=TH181)
lad_e = {}
c6, nN6, qN6 = nahtmodell(2, 3)
num6 = [np.diag(x.T @ x).real for x in c6]
for phi in (np.pi/5, np.pi/2):
    g_ = lambda q, n_, s_, phi=phi: np.exp(-1j*(phi*np.sum(q**2)/2 + 0.7*sum(eps_l[a]*(num6[0*3 + a][s_] + num6[3 + a][s_])
                                                                              for a in range(3))))
    lad_e[f"phi = {phi:.4f}"] = scan(2, 3, g_, thetas=TH181)
TH361 = np.linspace(0, 2*np.pi, 361)
amp = {str(t_): scan(2, 3, lambda q, n_, s_: np.exp(-1j*np.pi/2*np.sum(q**2)/2), t=t_, thetas=TH361)
       for t_ in ([1, 3, 5], [1, 2, 1], [1, 1, 2])}
c4b, nN4b, _ = nahtmodell(2, 2)
vorsp = {}
for bias in (0.3, 1.0, np.pi):
    g_ = lambda q, n_, s_, bias=bias: np.exp(-1j*(np.pi/2*np.sum(q**2)/2 + bias*n_[0]))
    vorsp[f"{bias:.4f}"] = scan(2, 2, g_, thetas=TH181)
ZAHLEN["leiter_amplituden_vorspannung"] = {"Leiter (3 Niveaus) leckfrei/Austausch": leiter, "mit Niveauenergien": lad_e,
                                           "Amplituden t (theta in [0, 2 pi])": amp, "Vorspannung": vorsp}
check(f"[E] Energieleiter (zwei Nähte, drei Niveaus x = 1/4, 5/4, 9/4, Fockraum 64): leckfrei mit kohärentem Austausch für "
      f"phi = pi/5, pi/2, 2 genau bei theta/pi = {leiter['phi = 1.5708'][1]} (auch mit nahtsymmetrischen Niveauenergien: "
      f"{lad_e['phi = 1.5708'][1]}), bei phi = 2 pi nie ({leiter['phi = 6.2832'][1]}, Kontrolle); Amplituden t_k: (1,3,5) "
      f"{amp['[1, 3, 5]'][1]}, (1,2,1) {amp['[1, 2, 1]'][1]}, (1,1,2) {amp['[1, 1, 2]'][1]} (theta in [0, 2 pi]): jedes "
      f"ausgetauschte Niveau braucht theta t_k = pi/2 mod pi. Vorspannung (Gesetz nicht nahtsymmetrisch) 0,3 und 1: "
      f"{vorsp['0.3000'][1]}, {vorsp['1.0000'][1]} (keine Komposition), bei pi wieder {vorsp['3.1416'][1]} (pi ist "
      "nahtsymmetrisch mod 2 pi)",
      all(v[1] == [0.5] for k, v in leiter.items() if k != "phi = 6.2832") and leiter["phi = 6.2832"][1] == []
      and all(v[1] == [0.5] for v in lad_e.values()) and amp["[1, 3, 5]"][1] == [0.5, 1.5] and amp["[1, 2, 1]"][1] == []
      and amp["[1, 1, 2]"][1] == [] and vorsp["0.3000"][1] == [] and vorsp["1.0000"][1] == [] and vorsp["3.1416"][1] == [0.5])

# 2.3 zwei Arten mit allgemeinen Ladungen
lad_tab, phase_ok = {}, True
for q in ([1, 1], [1, -1], [1, 3], [3, -1], [1, 2], [2, 2]):
    for kappa in (0, 1):
        for phi in (np.pi/2, 1.0):
            g_ = lambda qq, n_, s_, phi=phi, kappa=kappa: np.exp(-1j*phi*np.sum(qq**2)/2)*(1j**((kappa*(qq[0] - qq[1])) % 4))
            lad_tab[f"q = {q}, kappa = {kappa}, phi = {phi:.4f}"] = scan(2, 2, g_, ladungen=q, thetas=TH181)[1]
    # Austauschphase bei pi/2 ohne Uhr
    cq, nq, qq_ = nahtmodell(2, 2, q)
    Hq = transfer(cq, 2, 2, np.array([[0, 0.5], [0.5, 0]]))
    iq = np.where(nq[0] + nq[1] == 2)[0]
    paarzust = [k for k, s_ in enumerate(iq) if nq[0][s_] == 1 and nq[1][s_] == 1
                and (cq[0].T @ cq[0])[s_, s_] != (cq[2].T @ cq[2])[s_, s_]]           # verschiedene Arten in A und B
    for phi in (np.pi/2, 1.0):
        Vq = np.exp(-1j*phi*(qq_[0]**2 + qq_[1]**2)/2)
        Uq = gatter_sektor(Hq, Vq, iq, np.pi/2)[np.ix_(paarzust, paarzust)]
        ew = np.linalg.eigvals(Uq)
        vP_ = np.exp(-1j*phi*(q[0]**2 + q[1]**2)/2)
        quot = sorted(ew/vP_, key=lambda z: abs(z - 1))
        phase_ok &= abs(quot[0] - 1) < 1e-10 and abs(quot[1] - np.exp(-1j*phi*q[0]*q[1])) < 1e-10
erwart = {}
for k_ in lad_tab:
    if "q = [1, 2], kappa = 1" in k_ or ("q = [2, 2]" in k_ and "phi = 1.5708" in k_):
        erwart[k_] = []
    else:
        erwart[k_] = [0.5]
ZAHLEN["ladungen"] = {"leckfreier Austausch theta/pi": lad_tab}
check(f"[E] zwei Arten mit Ladungen q_k, q_l (Gesetz e^(-i phi (Q_A^2 + Q_B^2)/2) R^kappa, 24 Fälle, 181 Winkel): leckfreier "
      f"Austausch genau bei theta = pi/2 für (1,1), (1,-1), (1,3), (3,-1) mit und ohne Uhr; (1,2) ohne Uhr bei pi/2, mit Uhr "
      f"nie (gemischte Parität: i^(q_k - q_l) != i^(q_l - q_k), das Gesetz ist nicht nahtsymmetrisch); (2,2) bei phi = pi/2 "
      f"nie (Austauschphase 4 phi = 2 pi trivial), bei phi = 1 bei pi/2. Austauschphase v_D/v_P = e^(-i phi q_k q_l) aus den "
      f"Gatter-Eigenwerten ({phase_ok}): bilinear in den Ladungen, die Z4-Wechselphase i^(Q_A Q_B) von z120 allgemein",
      lad_tab == erwart and phase_ok)

# 2.4 drei Nähte
c3, nN3, qN3 = nahtmodell(3, 2)
tot3 = sum(nN3)
idx3 = np.where(tot3 == 3)[0]
phys3 = [k for k, s_ in enumerate(idx3) if all(n_[s_] == 1 for n_ in nN3)]
V3 = {phi: np.exp(-1j*phi*sum(n_**2 for n_ in nN3)/2) for phi in (np.pi/2, 1.0)}


def gatter3(u, phi):
    H = transfer(c3, 3, 2, log_unitaer(u))
    Hs = H[np.ix_(idx3, idx3)]
    B = expm(-1j*Hs)
    return B.conj().T @ np.diag(V3[phi][idx3]) @ B


Wab = np.zeros((3, 3))
Wab[0, 1] = Wab[1, 0] = 0.5
hz = np.roll(np.eye(3), 1, 0) + np.roll(np.eye(3), -1, 0)
paarw = [round(float(th/np.pi), 4) for th in TH181
         if (lambda r: r[0] < 1e-9 and r[1] > 1e-6)(leck_austausch(gatter3(expm(-1j*th*Wab), np.pi/2), phys3))]
zykl, zykl_had = [], True
for th in np.linspace(0, 2*np.pi/3, 241)[:-1]:
    u = expm(-1j*th*hz)
    l_, a_, _ = leck_austausch(gatter3(u, np.pi/2), phys3)
    if l_ < 1e-9 and a_ > 1e-6:
        zykl.append(round(float(th/np.pi), 4))
        zykl_had &= np.allclose(np.abs(u)**2, 1/3, atol=1e-10)
w3 = np.exp(2j*np.pi/3)
F3 = np.array([[w3**(j*k) for k in range(3)] for j in range(3)])/np.sqrt(3)
f3_tab = {}
for phi in (np.pi/2, 1.0):
    l_, a_, Up3 = leck_austausch(gatter3(F3, phi), phys3)
    ew = np.linalg.eigvals(Up3)
    n_P = int(np.sum(np.abs(ew - np.exp(-1.5j*phi)) < 1e-9))
    n_D = int(np.sum(np.abs(ew - np.exp(-2.5j*phi)) < 1e-9))
    f3_tab[f"{phi:.4f}"] = {"Leckage": float(l_), "Austausch": float(a_), "Eigenwerte v_P": n_P, "Eigenwerte v_D": n_D}
# S = 3/2: v_P-Eigenraum = artensymmetrische Zustände
l_, a_, Up3 = leck_austausch(gatter3(F3, np.pi/2), phys3)
arten3 = [tuple(int((c3[2*i + 1].T @ c3[2*i + 1])[idx3[k], idx3[k]]) for i in range(3)) for k in phys3]   # Art je Naht
Sz_tot = np.diag([sum(0.5 - a for a in ar) for ar in arten3])
ew, ev = np.linalg.eig(Up3)
vP_raum = ev[:, np.abs(ew - np.exp(-1.5j*np.pi/2)) < 1e-9]
# Projektor auf total symmetrische Artenzustände: Symmetrisator über Permutationen der Nähte (mit Fermionvorzeichen aus JW)
def perm_op(p):
    """Nahtpermutation als Fockoperator c_(i a) -> c_(p(i) a); auf 'eine je Naht' = sgn(p) x Artenumbenennung."""
    u = np.zeros((3, 3))
    for i in range(3):
        u[p[i], i] = 1
    H = transfer(c3, 3, 2, log_unitaer(u))
    return expm(-1j*H[np.ix_(idx3, idx3)])[np.ix_(phys3, phys3)]


def sgn(p):
    return (-1)**sum(1 for i in range(len(p)) for j in range(i + 1, len(p)) if p[i] > p[j])


Sym = sum(sgn(p)*perm_op(p) for p in itertools.permutations(range(3)))/6     # Projektor auf artensymmetrisch (S = 3/2)
rg_sym = int(np.linalg.matrix_rank(Sym, tol=1e-9))
proj = vP_raum @ np.linalg.pinv(vP_raum)
sym_rel = float(np.linalg.norm(proj - Sym)) if rg_sym == 4 else 1.0
zuf_leck = []
for _ in range(40):
    zuf_leck.append(leck_austausch(gatter3(unitary_group.rvs(3, random_state=rng), np.pi/2), phys3)[0])


def adj_krit(u):
    """[P_h, P_h'] im Raum der 3x3-Matrizen: h = spurfreie Diagonale, h' = u^+ h u."""
    D1, D2 = np.diag([1, -1, 0]).astype(complex), np.diag([1, 1, -2]).astype(complex)/np.sqrt(3)
    def proj_(mats):
        Mv = np.array([m.ravel() for m in mats]).T
        Q_, _ = np.linalg.qr(Mv)
        return Q_ @ Q_.conj().T
    Ph = proj_([D1/np.sqrt(2), D2/np.sqrt(2)])
    Phs = proj_([u.conj().T @ D1 @ u, u.conj().T @ D2 @ u])
    return np.linalg.norm(Ph @ Phs - Phs @ Ph)


krit = []
for u in [F3, expm(-1j*np.pi/2*Wab), expm(-1j*np.pi/3*Wab), expm(-1j*2*np.pi/9*hz)] + \
        [unitary_group.rvs(3, random_state=rng) for _ in range(10)]:
    krit.append((leck_austausch(gatter3(u, 1.0), phys3)[0] < 1e-9, adj_krit(u) < 1e-9))
ZAHLEN["drei_naehte"] = {"paarweise theta/pi": paarw, "zyklisch theta/pi": zykl, "zyklisch Hadamard": zykl_had,
                         "F_3": f3_tab, "zufällig min Leckage": float(min(zuf_leck)),
                         "Kriterium [P_h, P_h'] = 0 <=> leckfrei": all(a == b for a, b in krit)}
check(f"[E] drei Nähte (je zwei Arten, N = 3, eine Anregung je Naht, Gesetz e^(-i phi sum Q^2/2)): paarweiser Transfer A-B "
      f"mit Austausch genau bei theta/pi = {paarw}; zyklischer Dreinahttransfer e^(-i theta h_zykl) leckfrei mit Austausch bei "
      f"theta/pi = {zykl} (= 2/9, 4/9), dort |u_ij|^2 = 1/3 ({zykl_had}): komplexe Hadamard-Matrix (Klasse F_3); F_3 selbst "
      f"leckfrei (phi = pi/2: Leckage {f3_tab['1.5708']['Leckage']:.0e}), Gatter v_P auf {f3_tab['1.5708']['Eigenwerte v_P']} "
      f"und v_D auf {f3_tab['1.5708']['Eigenwerte v_D']} Zuständen (v_P-Raum = artensymmetrisch S = 3/2: Abstand "
      f"{sym_rel:.0e}); 40 zufällige U(3): Leckage >= {min(zuf_leck):.3f}; leckfrei <=> [P_h, P_(u^+ h u)] = 0 im Raum der "
      f"Matrizen (14 Fälle: {all(a == b for a, b in krit)}), also u monomial, Paar-Viertel mit Zuschauer oder |u_ij|^2 = 1/3",
      paarw == [0.5] and zykl == [0.2222, 0.4444] and zykl_had and f3_tab["1.5708"]["Leckage"] < 1e-12
      and f3_tab["1.5708"]["Eigenwerte v_P"] == 4 and f3_tab["1.5708"]["Eigenwerte v_D"] == 4 and sym_rel < 1e-9
      and min(zuf_leck) > 0.05 and all(a == b for a, b in krit) and sum(a for a, b in krit) == 3)
befund("E/C: Für zwei Nähte ist die Vierteldrehung der einzige nichttriviale leckfreie Transfer, für jede Leiter, Artenzahl, "
       "Ladung und jede Phase phi außerhalb 2 pi Z, sofern das Gesetz nahtsymmetrisch ist (sonst gibt es keinen); mit "
       "niveauabhängigen Amplituden braucht jedes ausgetauschte Niveau die Vierteldrehung. Drei Nähte: paarweise wieder pi/2, "
       "echte Dreinahttransfers lassen zusätzlich die Hadamard-Klasse F_3 zu (zyklisch theta = 2 pi/9): ein Gegenbeispiel zum "
       "wörtlichen Winkel, nicht zum Prinzip. Gemeinsam ist die gegenseitige Unverzerrtheit von Nahtbasis und übertragener "
       "Nahtbasis; für zwei Nähte ist das {A, B} = 0", True)
say(f"   (Laufzeit bis Teil 2: {time.time() - T_START:.1f} s)")

# =================================================================== Teil 3
say("\nTeil 3: Ist Leckfreiheit nativ?")
# 3.1 Sektoren
sek_par = sorted(k for k, s_ in enumerate(i2) if PAd[s_] == -1)
sek_uhr = sorted(k for k, s_ in enumerate(i2) if abs(Rd[s_] - 1) < 1e-12)
sek_E = sorted(k for k, s_ in enumerate(i2) if abs(HQd[s_] - 1) < 1e-12)
sek_ok = sek_par == sek_uhr == sek_E == sorted(phys2)
trennt = {}
for q in ([1, -1], [1, 3], [2, 2]):
    cq, nq, qq_ = nahtmodell(2, 2, q)
    iq = [s_ for s_ in np.where(nq[0] + nq[1] == 2)[0]                     # gemischter Artenblock (je eine Anregung
          if (cq[0].T @ cq[0])[s_, s_] + (cq[2].T @ cq[2])[s_, s_] == 1]     # jeder Art); Transfer erhält die Arten
    ph = [s_ for s_ in iq if nq[0][s_] == 1 and nq[1][s_] == 1]
    rest_ = [s_ for s_ in iq if s_ not in ph]
    par = lambda s_: (-1)**nq[0][s_]
    uhr = lambda s_: 1j**((qq_[0][s_] - qq_[1][s_]) % 4)
    trennt[str(q)] = {"Parität": len({par(s_) for s_ in ph}) == 1 and not ({par(s_) for s_ in ph} & {par(s_) for s_ in rest_}),
                      "Paaruhr": len({uhr(s_) for s_ in ph}) == 1 and not ({uhr(s_) for s_ in ph} & {uhr(s_) for s_ in rest_})}
gew_doppel = float(mp.exp(-4*mp.pi*(2 - 1)))
ZAHLEN["sektoren"] = {"N = 2: Parität = Uhr = H_Q-Niveau = eine je Naht": sek_ok, "Trennung je Ladung": trennt,
                      "KMS-Gewicht (2,0) relativ (1,1) bei beta = 4 pi, E = Q^2/2": gew_doppel}
check(f"[E] 'Eine Anregung je Naht' auf N = 2: dieselbe Zustandsmenge wie ungerade Parität jeder Naht, Paaruhr R = +1 und "
      f"Ladungsenergie H_Q = 1 ({sek_ok}). Mit Ladungen (1,-1) und (1,3) trennt im gemischten Artenblock (der Transfer "
      f"erhält die Arten; gleiche Arten sind Transfer-Singuletts) die Nahtparität weiter "
      f"({trennt['[1, -1]']['Parität']}, {trennt['[1, 3]']['Parität']}), die Paaruhr auch ({trennt['[1, -1]']['Paaruhr']}, "
      f"{trennt['[1, 3]']['Paaruhr']}); bei geraden Ladungen (2,2) trennt die Paaruhr nicht ({trennt['[2, 2]']['Paaruhr']}), "
      f"die Parität schon ({trennt['[2, 2]']['Parität']}). Der native Zustand besetzt den Doppelsektor (z93: der Ladungssektor "
      f"ist Teil des Zustands; Gewicht e^(-4 pi) = {gew_doppel:.3e} > 0): Leckfreiheit muss auf allen Sektoren gelten",
      sek_ok and trennt["[1, -1]"]["Parität"] and trennt["[1, 3]"]["Parität"] and trennt["[1, -1]"]["Paaruhr"]
      and trennt["[1, 3]"]["Paaruhr"] and not trennt["[2, 2]"]["Paaruhr"] and trennt["[2, 2]"]["Parität"] and gew_doppel > 0)

# 3.2 Klassifikation der Z4-Gradierungen als Superauswahl
U_th = {}
for th in TH181:
    B = expm(-1j*th*Jx)
    U_th[th] = B.conj().T @ np.diag(Vd_h) @ B
klassen = {}
for a in range(4):
    for b in range(4):
        G = np.diag(1j**((a*QA + b*QB) % 4))
        erlaubt = [round(float(th/np.pi), 4) for th in TH181 if np.linalg.norm(U_th[th] @ G - G @ U_th[th]) < 1e-9]
        klassen[(a, b)] = ("leer" if len(erlaubt) == len(TH181) else "erzwingt" if erlaubt == [0.0, 0.5, 1.0]
                           else "verbietet" if erlaubt == [0.0, 1.0] else str(erlaubt))
QAm = np.diag(QA.astype(float))
erlaubt_Q = [round(float(th/np.pi), 4) for th in TH181 if np.linalg.norm(U_th[th] @ QAm - QAm @ U_th[th]) < 1e-9]
s20 = int(np.where((QA == 2) & (QB == 0))[0][0])
s02 = int(np.where((QA == 0) & (QB == 2))[0][0])
mischung = abs(U_th[TH181[90]][s02, s20])
erzw = sorted(k for k, v in klassen.items() if v == "erzwingt")
leer = sorted(k for k, v in klassen.items() if v == "leer")
verb = sorted(k for k, v in klassen.items() if v == "verbietet")
ZAHLEN["superauswahl_z4"] = {"erzwingt pi/2": [str(k) for k in erzw], "leer": [str(k) for k in leer],
                             "verbietet jede Komposition": [str(k) for k in verb], "Nahtladung Q_A": erlaubt_Q,
                             "|<0,2|U|2,0>| bei pi/2": mischung}
check(f"[E] Superauswahl durch Z4-Gradierungen i^(a Q_A + b Q_B) auf dem ganzen Fockraum (16 x 181 Winkel, phi = pi/2): "
      f"die Vierteldrehung erzwingen genau (a, b) = {erzw} (b - a = 2 mod 4: Nahtparitäten (-1)^F_A, (-1)^F_B und Paaruhr "
      f"i^(+-(Q_A - Q_B)) = i^(-+N) (-1)^F_A); nichts fordern die {len(leer)} Gradierungen mit a = b (Gesamtladung, "
      f"Gesamtparität, Produktuhr i^N); jede nichttriviale Komposition verbieten die übrigen {len(verb)} (darunter die "
      f"Nahtuhr i^Q_A) und die Nahtladung Q_A selbst ({erlaubt_Q}): Bei pi/2 mischt U (2,0) und (0,2) mit |sin(phi/2)| = "
      f"{mischung:.4f}, verträglich mit z108 (erhaltene Einzeldynamik heißt Produkt)",
      erzw == [(0, 2), (1, 3), (2, 0), (3, 1)] and leer == [(0, 0), (1, 1), (2, 2), (3, 3)] and len(verb) == 8
      and erlaubt_Q == [0.0, 1.0] and abs(mischung - np.sin(np.pi/4)) < 1e-12)

# 3.3 drei Arten: Parität bleibt, Ladungsenergie nicht
c6, nN6, qN6 = nahtmodell(2, 3)
QA6, QB6 = qN6
N6 = QA6 + QB6
Jx6 = transfer(c6, 2, 3, np.array([[0, 0.5], [0.5, 0]]))
PA6 = np.diag((-1.0)**QA6)
par6 = {}
for trial in range(3):
    g_ = rng.normal(size=4)
    h_ = rng.normal()
    f_ = np.array([g_[qa] + g_[qb] + h_*qa*qb for qa, qb in zip(QA6, QB6)])
    V6 = np.diag(np.exp(-1j*f_))
    for th in (np.pi/2, np.pi/3, 0.4):
        B = expm(-1j*th*Jx6)
        U6 = B.conj().T @ V6 @ B
        par6[(trial, round(th, 4))] = float(np.linalg.norm(U6 @ PA6 - PA6 @ U6))
B = expm(-1j*np.pi/2*Jx6)
U6 = B.conj().T @ np.diag(np.exp(-1j*np.pi/2*(QA6**2 + QB6**2)/2)) @ B
HQ6 = np.diag((QA6**2 + QB6**2)/2.0)
E_verl = {N: float(np.linalg.norm((U6 @ HQ6 - HQ6 @ U6)[np.ix_(np.where(N6 == N)[0], np.where(N6 == N)[0])])) for N in range(7)}
ZAHLEN["drei_arten"] = {"||[U, (-1)^F_A]|| (Gesetz, theta)": {str(k): v for k, v in par6.items()},
                        "||[U, H_Q]|| je N bei pi/2": E_verl}
check(f"[E] Nahtparität allgemein (zwei Nähte, drei Arten, Fockraum 64, drei zufällige nahtsymmetrische Gesetze "
      f"g(Q_A) + g(Q_B) + h Q_A Q_B): bei pi/2 ist [U, (-1)^F_A] = 0 (max {max(v for k, v in par6.items() if k[1] == 1.5708):.0e}), "
      f"bei pi/3 und 0,4 nicht (min {min(v for k, v in par6.items() if k[1] != 1.5708):.3f}): (-1)^F_A ~ e^(i pi J_z) und U = "
      f"f(J_y^2) ist unter der Halbdrehung invariant. Die Ladungsenergie (z105s Kriterium 'vertauscht mit E') bleibt dagegen "
      f"nur auf N <= 2 erhalten, auf N = 3 (Transferspin 3/2) nicht ({E_verl[3]:.3f}): 'auf der Schale' ist nur auf dem "
      "Zweianregungssektor gleichwertig, die native Form ist die Parität",
      max(v for k, v in par6.items() if k[1] == 1.5708) < 1e-12 and min(v for k, v in par6.items() if k[1] != 1.5708) > 0.1
      and E_verl[3] > 0.1 and all(E_verl[N] < 1e-12 for N in (0, 1, 2, 4, 5, 6)))

# 3.4 drei Nähte unter Nahtparität
H3full = {}
QN3 = [np.diag(q.astype(float)) for q in qN3]
V3full = np.diag(np.exp(-1j*np.pi/2*sum(q**2 for q in qN3)/2))
par3 = {}
for name, u in (("Paar pi/2", expm(-1j*np.pi/2*Wab)), ("Paar pi/3", expm(-1j*np.pi/3*Wab)), ("F_3", F3),
                ("zyklisch 2 pi/9", expm(-1j*2*np.pi/9*hz))):
    B = expm(-1j*transfer(c3, 3, 2, log_unitaer(u)))
    U3 = B.conj().T @ V3full @ B
    par3[name] = [round(float(np.linalg.norm(U3 @ np.diag((-1.0)**q) - np.diag((-1.0)**q) @ U3)), 6) for q in qN3]
ZAHLEN["drei_naehte_paritaet"] = par3
check(f"[E] drei Nähte (Fockraum 64) unter Nahtparität: Paar-Vierteldrehung erhält alle drei Paritäten ({par3['Paar pi/2']}), "
      f"Paar pi/3 nicht ({par3['Paar pi/3']}), F_3 und zyklisches Hadamard verletzen alle drei ({par3['F_3']}, "
      f"{par3['zyklisch 2 pi/9']}; sie mischen (2,0,1), (1,2,0), (0,1,2)): die Nahtparität wählt die paarweise Komposition "
      "und damit tau = i statt der Z3-Struktur",
      max(par3["Paar pi/2"]) < 1e-12 and min(par3["Paar pi/3"][:2]) > 0.1 and min(par3["F_3"]) > 0.1
      and min(par3["zyklisch 2 pi/9"]) > 0.1)

# 3.5 Tomita-Paaruhr; Einzel- gegen Doppeltransfer
pA = np.array([np.exp(-np.pi*q) for q in (0, 1, 1, 2)])
pA /= pA.sum()
QA4 = np.array([0, 1, 1, 2])
Psi = np.zeros(16, dtype=complex)
for k in range(4):
    Psi[4*k + k] = np.sqrt(pA[k])
bx = []
for i in range(16):
    m_ = np.zeros((4, 4), dtype=complex)
    m_[i // 4, i % 4] = 1
    bx.append(m_)
T_ = np.array([np.kron(x, np.eye(4)) @ Psi for x in bx]).T
Td = np.array([np.kron(x.conj().T, np.eye(4)) @ Psi for x in bx]).T
L_ = Td @ np.linalg.inv(T_).conj()                                   # S v = L conj(v)
Delta = (L_.conj().T @ L_).conj()
tom_abw = np.linalg.norm(Delta - np.kron(np.diag(pA), np.diag(1/pA)))
uhr_t = np.exp(-0.5j*np.log(np.diag(Delta).real))
QAf, QBf = np.repeat(QA4, 4), np.tile(QA4, 4)
rel = float(np.ptp(np.angle(uhr_t*1j**(-(QAf - QBf) % 4))))
prod = float(np.ptp(np.angle(uhr_t*1j**(-(QAf + QBf) % 4))))
einzel = c4[0].T @ c4[2]                                             # c_(A,0)^+ c_(B,0)
doppel = c4[0].T @ c4[1].T @ c4[3] @ c4[2]                           # z105-Typ: zwei Teilchen B -> A
Pges = np.diag((-1.0)**NT)
e_ges = np.linalg.norm(einzel @ Pges - Pges @ einzel)
e_A = np.linalg.norm(einzel @ PA + PA @ einzel)
d_A = np.linalg.norm(doppel @ PA - PA @ doppel)
ZAHLEN["tomita_transfer"] = {"|Delta - rho_A (x) rho_B^-1|": float(tom_abw), "Phasenspanne Delta^(-i/2) i^-(Q_A-Q_B)": rel,
                             "Phasenspanne Delta^(-i/2) i^-(Q_A+Q_B)": prod, "Einzeltransfer gesamtgerade": float(e_ges),
                             "Einzeltransfer antivertauscht mit (-1)^F_A": float(e_A), "Doppeltransfer vertauscht": float(d_A)}
check(f"[E] Tomita: Naht A (zwei Moden, K_A = pi Q_A) mit Partner B in der Reinigung; der Tomita-Operator gibt Delta = "
      f"rho_A (x) rho_B^-1 (Abweichung {tom_abw:.0e}) und Delta^(-i/2) ~ i^(Q_A - Q_B) (Phasenspanne {rel:.0e}), nicht "
      f"i^(Q_A + Q_B) (Spanne {prod:.3f}): die erzwingende Paaruhr ist die Moduluhr des Original-Partner-Paars. Einzeltransfer "
      f"c_A^+ c_B: gesamtgerade (Rest {e_ges:.0e}), aber ungerade in der Nahtparität (Antikommutator-Rest {e_A:.0e}); "
      f"z105-Doppeltransfer: nahtgerade (Rest {d_A:.0e}). Die globale Univalenz (z106: dieselbe gerade Algebra) lässt "
      "Einzeltransfers zu; die Nahtparität ist eine eigene Prämisse, die z105s Kopplungen erfüllen",
      tom_abw < 1e-12 and rel < 1e-12 and prod > 0.1 and e_ges < 1e-12 and e_A < 1e-12 and d_A < 1e-12)

# 3.6 Verträglichkeit mit z107 und z94
phi_blind = {f"{phi:.4f}": [round(float(th/np.pi), 4) for th in TH181
                            if leck_austausch(gatter_sektor(Jx, np.exp(-1j*phi*HQd), i2, th), phys2)[0] < 1e-10]
             for phi in (np.pi/5, np.pi/3, np.pi/2, 2.0, np.pi)}
z107_pi3 = leck_austausch(gatter_sektor(Jx, Vd_h, i2, np.pi/3), phys2)[0]
ZAHLEN["vertraeglichkeit"] = {"leckfrei je phi": phi_blind, "Leckage theta = pi/3 (z107 spannungsfrei)": z107_pi3,
                              "j(pi/2)": jm["pi/2"], "j(pi/3)": jm["pi/3"]}
check(f"[E] Verträglichkeit: z107s spannungsfreier Fall theta = pi/3 (A* gilt, j = {jm['pi/3']:.1f}) ist als Komposition "
      f"undicht (Leckage {z107_pi3:.6f} = |cos||sin||1 - e^(-i pi/2)|): A* und Komposition sind unabhängig. Die leckfreie Menge "
      f"ist für phi = pi/5, pi/3, pi/2, 2, pi dieselbe ({set(map(tuple, phi_blind.values()))}) und hängt vom Zustand nicht ab: "
      "die Komposition legt theta fest, nicht beta oder g, gibt also P_sym, nicht P_tn/P_clock (z94: unabhängig)",
      abs(z107_pi3 - np.cos(np.pi/3)*np.sin(np.pi/3)*np.sqrt(2)) < 1e-12
      and all(v == [0.0, 0.5, 1.0] for v in phi_blind.values()))
befund("C: Unter den geprüften Superauswahloptionen (alle Z4-Gradierungen, Nahtladung; der native Zustand besetzt alle "
       "Ladungssektoren, z93) ist Leckfreiheit genau unter der Nahtparität P_univ nativ: physikalische Kompositionen erhalten "
       "die Fermionparität jeder Naht (für zwei Nähte mit ungeraden Ladungen gleichwertig: die Paaruhr i^(Q_A - Q_B) = "
       "Delta^(-i/2) des Original-Partner-Paars ist superselektiert; das sind z108s Bedingungen 'Uhr und Parität erhalten'). "
       "Dann gilt: nichttriviale kohärente Komposition + P_univ => nahtsymmetrisches Gesetz (Vorzeichensymmetrie A -> -A), "
       "paarweise Komposition, theta = pi/2 => {A, B} = 0 auf dem Einanregungssektor => P_sym (Cl_2, tau = i). P_univ liegt "
       "zwischen den Rahmenoptionen: feinere Superauswahl (Nahtladung, Nahtuhr) verbietet jede Komposition (z108, z93), "
       "gröbere (Gesamtladung, Gesamtparität, Produktuhr) fordert nichts; aus der globalen Univalenz folgt sie nicht (z106), "
       "z105s native Zweinahtterme erfüllen sie", True)
say(f"   (Laufzeit bis Teil 3: {time.time() - T_START:.1f} s)")

# =================================================================== Teil 4
say("\nTeil 4: Energieroute zu S(n)")
# 4.1 Gauß
spann = mp.nsum(lambda k: (2*k + 1)*YV**(2*k + 1)/(1 + YV**(2*k + 1)), [0, mp.inf])


def Zvoll(t):
    t = mp.mpf(t)
    return mp.exp(mp.pi*t/12)*mp.nprod(lambda k: (1 + mp.exp(-mp.pi*t*(2*k + 1)))**2, [0, mp.inf])


sdef_gauss = max(abs(Zvoll(t) - Zvoll(1/mp.mpf(t)))/Zvoll(t) for t in ("0.5", "0.8", "2", "3"))
ZAHLEN["gauss"] = {"sum m y^m/(1 + y^m) - 1/24": mp.nstr(spann - mp.mpf(1)/24, 5), "S-Defekt Inhalt": mp.nstr(sdef_gauss, 3)}
check(f"[E] Gaußzustand: mittlere Modulenergie sum_(m ungerade) m y^m/(1 + y^m) = 1/24 bei y = e^-pi (Abweichung "
      f"{mp.nstr(spann - mp.mpf(1)/24, 3)}), also <T> ~ E_0 + <e>/2 = 0 mit dem Casimir E_0 = -1/24; der Inhalt "
      f"q^(-1/24) prod (1 + q^(m/2))^2 ist S-invariant (Defekt {mp.nstr(sdef_gauss, 2)}): A*(a) und A*(b) gelten (z107)",
      abs(spann - mp.mpf(1)/24) < mp.mpf(10)**-45 and sdef_gauss < mp.mpf(10)**-40)

# 4.2 Inhalt der Kandidaten: Modul- gegen geometrische Energie
z94 = json.loads((HERE / "z94_ergebnis.json").read_text())
fz = z94["zahlen"]["aussenbehauptung_gleiche_saaten"]["f(n) für N = 6"]
C6 = list(itertools.product((0, 1), repeat=6))
H6 = {c_: sum(c_) + 4*fz.get("".join(map(str, c_)), 0) for c_ in C6}
Zk = mp.fsum(YV**H6[c_] for c_ in C6)
E_mod_k = mp.fsum(H6[c_]*YV**H6[c_] for c_ in C6)/Zk
E_geo_k = mp.fsum(sum(c_)*YV**H6[c_] for c_ in C6)/Zk
E_g = 6*YV/(1 + YV)
Zs_k = sp.expand(sum(yS**H6[c_] for c_ in C6))
Z0s_k = sp.expand((1 + yS)**6)
Dk = sp.expand(sp.diff(Zs_k, yS)*yS*Z0s_k - Zs_k*yS*sp.diff(Z0s_k, yS))
Dk_int = all(cc.is_integer for cc in sp.Poly(Dk, yS).coeffs())
casimir_noetig = -mp.mpf(1)/24 + (E_g - E_mod_k)/2
ZAHLEN["kandidat_z94"] = {"<e> modular": mp.nstr(E_mod_k, 12), "<e> geometrisch (Einteilchen)": mp.nstr(E_geo_k, 12),
                          "<e> Gauß": mp.nstr(E_g, 12), "D(y) ganzzahlig, Grad": [Dk_int, sp.degree(Dk, yS)],
                          "Casimir, der <T> = 0 erzwingen würde": mp.nstr(casimir_noetig, 10)}
check(f"[E] Kandidat aus z94 (sechs Saaten, y^(|n| + 4 f(n)), gleiche Einzelbesetzungen): geometrische (Einteilchen-) Energie "
      f"{mp.nstr(E_geo_k, 10)} = Gauß {mp.nstr(E_g, 10)}, also blind; Modulenergie des eigenen Gesetzestorus "
      f"{mp.nstr(E_mod_k, 10)}: mit dem zustandsunabhängigen Casimir -1/24 ist <T> ~ ({mp.nstr(E_mod_k, 6)} - "
      f"{mp.nstr(E_g, 6)})/2 != 0, A*(b) scheitert. Das Polynom D(y) = y Z' Z_0 - Z' y Z_0' ist ganzzahlig ({Dk_int}, Grad "
      f"{sp.degree(Dk, yS)}) und nicht null, also (Gelfond) D(e^-pi) != 0. Ein pro Zustand gewählter Casimir "
      f"({mp.nstr(casimir_noetig, 8)} statt -1/24) machte A*(b) leer",
      abs(E_geo_k - E_g) < mp.mpf(10)**-40 and abs(E_mod_k - E_g) > mp.mpf("0.1") and Dk_int and Dk != 0)


# 4.3 S-Defekt (A*(a) auf den Kandidaten)
def s_defekt(Zp, Z0, ts):
    """min_c max_t |g(t) - 2 pi c (t - 1/t)| mit g = log R(t) - log R(1/t), R = Z'/Z_0 bei y = e^(-pi t)."""
    def logR(t):
        yy = mp.exp(-mp.pi*t)
        return mp.log(mp.fsum(w*yy**e for e, w in Zp.items())/mp.fsum(w*yy**e for e, w in Z0.items()))
    g = [logR(t) - logR(1/t) for t in ts]
    a = [2*mp.pi*(t - 1/t) for t in ts]
    kand = [gi/ai for gi, ai in zip(g, a)]
    for i in range(len(ts)):
        for j in range(i + 1, len(ts)):
            for sg in (1, -1):
                if a[i] - sg*a[j] != 0:
                    kand.append((g[i] - sg*g[j])/(a[i] - sg*a[j]))
    best = min(kand, key=lambda cc: max(abs(gi - cc*ai) for gi, ai in zip(g, a)))
    return best, max(abs(gi - best*ai) for gi, ai in zip(g, a))


TS = [mp.mpf(x) for x in ("0.4", "0.55", "0.7", "0.85", "1.2", "1.5", "2", "2.5")]
Zk_d = {e: mp.mpf(v) for e, v in Counter(H6.values()).items()}
Z0_d = {k: mp.mpf(sp.binomial(6, k)) for k in range(7)}
sd_k = s_defekt(Zk_d, Z0_d, TS)
sd_shift = s_defekt({e + 4: w for e, w in Z0_d.items()}, Z0_d, TS)
sd_orig = s_defekt(Z0_d, Z0_d, TS)


def eps_von(ms):
    return [1 if m % 4 == 1 else -1 for m in ms]


def paarbloecke(ms):
    """wie z94/z107: Q_- -Randblöcke (Exponenten je Ladung) und mu = N/D gekürzt."""
    eps = eps_von(ms)
    rand = defaultdict(list)
    for st in itertools.product(((0, 0), (1, 1), (1, -1), (2, 0)), repeat=len(ms)):
        E = sum(m*nu for m, (nu, cc) in zip(ms, st))
        qm = sum(e*nu for e, (nu, cc) in zip(eps, st))
        rand[qm].append(E)
    mu = sum(e*2*yS**m/(1 + yS**m) for e, m in zip(eps, ms))
    N, D = sp.fraction(sp.cancel(sp.together(mu)))
    cf = lambda P: np.array([int(cc) for cc in reversed(sp.Poly(sp.expand(P), yS).all_coeffs())], dtype=np.int64)
    return {k: sorted(v) for k, v in rand.items()}, cf(N), cf(D), sp.expand(N), sp.expand(D)


def aufbau(ms, lo, hi):
    """wie z107: reelles Q_- -Randproblem (Blocksummen, Ladungsidentität D M' = N Z', Randblock fest)."""
    rand, Nc, Dc, _, _ = paarbloecke(ms)
    keys = sorted(rand)
    var = [(q, e) for q in keys for e in range(lo, hi + 1) if (e - q) % 4 == 0]
    idx = {v: i for i, v in enumerate(var)}
    G = max(len(Nc), len(Dc)) - 1
    A = lil_matrix((len(keys) + (hi + G - lo + 1), len(var)))
    lb, ub, r = [], [], 0
    for q in keys:
        for e in range(lo, hi + 1):
            if (q, e) in idx:
                A[r, idx[(q, e)]] = 1
        lb.append(len(rand[q]))
        ub.append(len(rand[q]))
        r += 1
    for e2 in range(lo, hi + G + 1):
        for (q, e), i in idx.items():
            d, val = e2 - e, 0
            if 0 <= d < len(Dc):
                val += q*Dc[d]
            if 0 <= d < len(Nc):
                val -= Nc[d]
            if val:
                A[r, i] = val
        lb.append(0)
        ub.append(0)
        r += 1
    rand_q = keys[0]
    fix = lil_matrix((1, len(var)))
    fix[0, idx[(rand_q, rand[rand_q][0])]] = 1
    return vstack([csr_matrix(A[:r]), csr_matrix(fix)]), lb + [1], ub + [1], var, rand


def energie_lp(ms, lo, hi, vorzeichen):
    """Extremum (vorzeichen = +1: max, -1: min) des Energiefunktionals f(w) = sum (e - E_0) y^e w über die reelle Relaxation
    (Ziel mit y^(e - lo) skaliert); Rückgabe (f in 50 Stellen, L1(Z' - Z), Spektrum Z')."""
    A, lb, ub, var, rand = aufbau(ms, lo, hi)
    Z0 = Counter(e for v in rand.values() for e in v)
    E0 = mp.fsum(e*cc*YV**e for e, cc in Z0.items())/mp.fsum(cc*YV**e for e, cc in Z0.items())
    yv = float(YV)
    cobj = np.array([-vorzeichen*(e - float(E0))*yv**(e - lo) for (q, e) in var])
    res = milp(cobj, constraints=LinearConstraint(A, lb, ub), integrality=np.zeros(len(var)),
               bounds=Bounds(np.zeros(len(var)), np.full(len(var), np.inf)))
    w = res.x
    f = mp.fsum((e - E0)*mp.mpf(float(x))*YV**e for (q, e), x in zip(var, w) if abs(x) > 1e-12)
    Zp = Counter()
    for (q, e), x in zip(var, w):
        Zp[e] += x
    L1 = sum(abs(Zp[e] - Z0.get(e, 0)) for e in set(Zp) | set(Z0))
    return f, L1, {e: mp.mpf(float(x)) for e, x in Zp.items() if abs(x) > 1e-12}, {e: mp.mpf(v) for e, v in Z0.items()}


f135_max, L135_max, Zp135, Z0135 = energie_lp((1, 3, 5), 0, 2*9 + 24, +1)
sd_reell = s_defekt(Zp135, Z0135, TS)
ZAHLEN["s_defekt"] = {"Kandidat z94": [mp.nstr(sd_k[0], 6), mp.nstr(sd_k[1], 6)],
                      "Verschiebung y^4 Z": [mp.nstr(sd_shift[0], 6), mp.nstr(sd_shift[1], 3)],
                      "Original": [mp.nstr(sd_orig[0], 3), mp.nstr(sd_orig[1], 3)],
                      "reelle Ecke (1,3,5)": [mp.nstr(sd_reell[0], 6), mp.nstr(sd_reell[1], 6), round(L135_max, 3)]}
check(f"[E] A*(a) auf den Kandidaten (Inhalt = Gaußschwanz x Fensterquotient R = Z'/Z, S-Invarianz heißt log R(t) - "
      f"log R(1/t) = 2 pi c (t - 1/t) für alle t; Minimax über c, acht Werte t): Kandidat z94 Defekt {mp.nstr(sd_k[1], 4)}, "
      f"reelle Nicht-Spektrenmischung (1,3,5) (L1 = {L135_max:.1f}) Defekt {mp.nstr(sd_reell[1], 4)}; Kontrollen: globale "
      f"Verschiebung y^4 Z Defekt {mp.nstr(sd_shift[1], 2)} bei c = {mp.nstr(sd_shift[0], 4)} (Casimir gleicht aus), Original "
      f"{mp.nstr(sd_orig[1], 2)}. Lemma: f(t) = e^(-2 pi c t) R(e^(-pi t)) ist meromorph in t und (bis auf den Faktor) "
      "2i-periodisch; f(t) = f(1/t) macht f bei unendlich meromorph, also rational, also konstant: R = y^(4k), "
      "Isospektralität, für jedes n und ohne Ganzzahligkeit",
      sd_k[1] > mp.mpf("1e-3") and sd_reell[1] > mp.mpf("1e-5") and L135_max > 1 and sd_shift[1] < mp.mpf("1e-30")
      and abs(sd_shift[0] + 2) < mp.mpf("1e-30") and sd_orig[1] < mp.mpf("1e-40"))


# 4.4 Restklassenkern mit Energieidentität (Fenster n = 2..5 wie z107)
def kern_mit_energie(ms, Yval):
    """Unbekannte P_q(Y) je Ladung q und c; Ladungsidentitäten (vier Restklassen, wie z107), Energieidentität Z'_r = c Z_r
    (aus (Z'/Z)' = 0) und verborgene Richtungen (Z'_r = 0, c = 0)."""
    rand, _, _, Ne, De = paarbloecke(ms)
    qs = sorted(rand)
    nq = len(qs)
    Yv = Fraction(Yval)
    lad = [[Fraction(0)]*(nq + 1) for _ in range(4)]
    for iq, q in enumerate(qs):
        for (e,), cc in sp.Poly(sp.expand((q*De - Ne)*yS**(q % 4)), yS).terms():
            lad[e % 4][iq] += Fraction(int(cc))*Yv**(e//4)
    ener = [[Fraction(0)]*(nq + 1) for _ in range(4)]
    verb_ = [[Fraction(0)]*(nq + 1) for _ in range(4)]
    for iq, q in enumerate(qs):
        r = q % 4
        ener[r][iq] += 1
        ener[r][nq] -= sum(Yv**((E - r)//4) for E in rand[q])
        verb_[r][iq] += 1
    cfix = [[Fraction(0)]*nq + [Fraction(1)]]

    def kern(rows):
        M = flint.fmpq_mat(len(rows), nq + 1, [flint.fmpq(x.numerator, x.denominator) for r_ in rows for x in r_])
        return nq + 1 - M.rank()
    return kern(lad) - 1, kern(lad + ener), kern(lad + verb_ + cfix)


tupel = {2: list(itertools.combinations(range(1, 14, 2), 2)), 3: list(itertools.combinations(range(1, 14, 2), 3)),
         4: list(itertools.combinations(range(1, 12, 2), 4)), 5: [(1, 3, 5, 7, 9), (1, 3, 5, 7, 11), (3, 5, 7, 9, 11)]}
kern_tab = {}
for n, liste in tupel.items():
    kern_tab[n] = Counter(kern_mit_energie(ms, Fraction(3, 7)) for ms in liste)
ZAHLEN["restklassenkern_energie"] = {str(n): {"Tupel": sum(v.values()), "(Ladung, mit Energie, verborgen)": [list(k) for k in v]}
                                     for n, v in kern_tab.items()}
kern_ok = all(len(v) == 1 and list(v)[0][0] == 2*n - 3 and list(v)[0][1] == 1 + list(v)[0][2]
              for n, v in kern_tab.items()) and [list(kern_tab[n])[0][2] for n in (2, 3, 4, 5)] == [0, 0, 1, 3]
check(f"[E] Fenster n = 2..5 (die {sum(sum(v.values()) for v in kern_tab.values())} Tupel aus z107, exakt bei Y = 3/7): "
      f"Restklassenkern der Ladungsidentitäten 2n - 3 = "
      f"{[list(kern_tab[n])[0][0] for n in (2, 3, 4, 5)]}; mit der Energieidentität (Z'/Z)' = 0 bleiben "
      f"{[list(kern_tab[n])[0][1] for n in (2, 3, 4, 5)]} = 1 + verborgen, verborgen = "
      f"{[list(kern_tab[n])[0][2] for n in (2, 3, 4, 5)]} (Umverteilungen innerhalb der Blöcke, die Z' nicht ändern). Die "
      "Energie tilgt also jede spektral sichtbare Richtung außer dem Original; die verborgenen nimmt der Max-Entropie-Schritt "
      "(gleiche Einteilchendaten + isospektral => rho' = rho_0)",
      kern_ok)

# 4.5 Ganzzahligkeit in der Variante B2
lp_tab = {}
for ms in ((1, 3, 5), (3, 5, 7), (1, 3, 7), (1, 5, 9)):
    fmax, Lmax, _, _ = energie_lp(ms, 0, 2*sum(ms) + 24, +1)
    fmin, Lmin, _, _ = energie_lp(ms, -24, 2*sum(ms) + 24, -1)
    lp_tab[str(ms)] = {"max f (Fenster e >= 0)": mp.nstr(fmax, 6), "L1 bei max": round(Lmax, 3),
                       "min f (Fenster -24)": mp.nstr(fmin, 6), "L1 bei min": round(Lmin, 3)}
    lp_tab[str(ms)]["_roh"] = (fmax, Lmax, fmin, Lmin)
nicht_starr = [ms for ms, v in lp_tab.items() if v["_roh"][0] > mp.mpf("1e-6") and v["_roh"][1] > 1
               and v["_roh"][2] < -1 and v["_roh"][3] > 1]
einseitig = [ms for ms, v in lp_tab.items() if abs(v["_roh"][0]) < mp.mpf("1e-12") and v["_roh"][1] < 1e-6
             and v["_roh"][2] < -1 and v["_roh"][3] > 1]
starr = [ms for ms, v in lp_tab.items() if v["_roh"][1] < 1e-6 and v["_roh"][3] < 1e-6]
for v in lp_tab.values():
    del v["_roh"]
ZAHLEN["reell_mit_energie"] = {"Tabelle": lp_tab, "zweiseitig (Gegenbeispiel zu B2 reell)": nicht_starr,
                               "einseitig (max f = 0 am Original)": einseitig, "reell starr": starr}
check(f"[E] Variante B2 braucht Ganzzahligkeit: in der reellen Relaxation von z107 nimmt das Energiefunktional f(w) = "
      f"sum (e - E_0) y^e w auf Nicht-Spektrenmischungen beide Vorzeichen an (max/min: "
      f"{', '.join(f'{k}: {v['max f (Fenster e >= 0)']} / {v['min f (Fenster -24)']}' for k, v in lp_tab.items() if k in nicht_starr)}); "
      f"die Lösungsmenge ist konvex, also gibt es reelle nicht isospektrale Lösungen mit exakt gleicher Energie für {nicht_starr}. "
      f"Einseitig ist {einseitig} (Nicht-Spektrenmischungen nur mit negativen Exponenten und f < 0, Maximum f = 0 am Original), "
      f"reell starr {starr} (Kontrolle). Mit ganzzahligen Vielfachheiten wird die eine Gleichung über Gelfond zur "
      "Polynomidentität (wie beim Kandidaten aus z94)",
      nicht_starr == ["(1, 3, 5)", "(3, 5, 7)"] and einseitig == ["(1, 3, 7)"] and starr == ["(1, 5, 9)"])
befund("E/C: S(n) folgt für jedes n aus A*, wenn A* auf den tatsächlichen Modulinhalt wirkt. Nötig ist: (T1) der Inhalt des "
       "Nahttorus ist die Modul-Zustandssumme Tr e^(-t K') des tatsächlichen Zustands (Gesetzestorus = Modultorus wie z94/z107, "
       "K' durch die Uhr normiert: ganze Exponenten = Q_- mod 4, also dasselbe tau_Gesetz), nicht die geometrische "
       "Einteilchenenergie (die ist blind); (T2) entweder B1: A*(a) gilt für diesen Inhalt (S-Invarianz mit irgendeinem "
       "Casimir) -- dann ist jeder Fensterkandidat mit Gaußschwanz isospektral, ohne Ganzzahligkeit; oder B2: der Casimir ist "
       "die zustandsunabhängige Vakuumgröße -1/24 und A*(b) gilt -- dann gleiche Modulenergie, mit ganzzahligen Vielfachheiten "
       "über Gelfond Z' = Z, für jedes n. Ein pro Zustand gewählter Casimir macht A*(b) leer; reelle Vielfachheiten genügen B2 "
       "nicht. Danach schließt Max-Entropie (z94)", True)
say(f"   (Laufzeit bis Teil 4: {time.time() - T_START:.1f} s)")

# =================================================================== Teil 5
say("\nTeil 5: Bilanz")
STATUS = {
    "Identifikation Transfer = z94-Binäroperationen": ("E", "A = -iR (Nahtparität), B = B_theta^+ A B_theta, {A, B} = 2 cos "
                                                           "theta, Gamma = 2 J_x; Leckage |1/2{A,B}||1/2[A,B]||v_P - v_D|"),
    "Paarsatz zwei Nähte": ("E", "leckfreier Austausch <=> theta_k = theta_l = pi/2 mod pi und nahtsymmetrisches Gesetz, "
                                 "v_D != v_P; Gatter v_D P_+ + v_P P_-; Leiter, Arten, Ladungen, phi"),
    "drei Nähte": ("E (Gegenbeispiel zum Winkel)", "paarweise pi/2; Dreinahttransfer: Hadamard-Klasse F_3 (zyklisch 2 pi/9)"),
    "Nativität der Leckfreiheit": ("C", "genau unter P_univ (Nahtparität / Paaruhr Delta^(-i/2)); feiner verbietet, gröber "
                                        "fordert nichts; nicht aus globaler Univalenz (z106)"),
    "P_sym aus Komposition": ("C unter P_univ", "nichttriviale kohärente Komposition + P_univ => paarweise, theta = pi/2, "
                                               "nahtsymmetrisches Gesetz => P_sym"),
    "Verträglichkeit z107, z94": ("E", "theta = pi/3 spannungsfrei, aber undicht; Komposition phi- und zustandsblind: kein "
                                       "P_tn/P_clock"),
    "S(n) über A*": ("E unter A* auf dem Modulinhalt", "B1 (S-Invarianz): jedes n, ohne Ganzzahligkeit; B2 (fester Casimir, "
                                                       "<T> = 0): jedes n mit Ganzzahligkeit (Gelfond)"),
    "S(n) ohne A*": ("O (unverändert, z107)", "arithmetisch ab n = 3"),
    "K19 Existenz der Naht": ("L", "nicht hergeleitet; die Komposition setzt zwei Nähte (Original und Partner) voraus"),
}
for k_, (s_, t_) in STATUS.items():
    say(f"   [{s_}] {k_}: {t_}")
ZAHLEN["status"] = STATUS
PRAEMISSEN = {
    "A* (Nahtprämisse)": "(a) S-Invarianz des Inhalts, (b) <T> = 0; gibt P_tn (z107) und, auf den Modulinhalt des "
                         "tatsächlichen Zustands angewandt, S(n) für jedes n",
    "Existenz einer nichttrivialen kohärenten Komposition": "zwei Nähte tauschen Anregungen kohärent aus (U nicht diagonal)",
    "P_univ (neu benannt)": "physikalische Kompositionen erhalten die Fermionparität jeder Naht; für zwei Nähte = Paaruhr "
                            "i^(Q_A - Q_B) = Delta^(-i/2) superselektiert",
    "Rahmen": "Hardy-Naht, Viertelzählung, CAR; die Ladungssuperauswahl je Naht gilt für die Einzelnaht, nicht für die "
              "Komposition (z108)",
}
ZAHLEN["praemissen"] = PRAEMISSEN
check("[C] Auftragskriterien: Identifikation wörtlich (Teil 1), Leckfreiheit für zwei Nähte allgemein und für drei Nähte mit "
      "benanntem Gegenbeispiel (Teil 2), Nativität auf eine benannte Superauswahl P_univ zurückgeführt und mit z107/z94 "
      "verträglich (Teil 3), Energieroute mit präzisen Annahmen über den Torusinhalt (Teil 4)",
      all(v for k, v in CHECKS.items()))

urteil = ("Der Ursprung reduziert sich in gescopter Form auf A* und die Existenz einer nichttrivialen kohärenten Komposition, "
          "mit einer benannten Superauswahlprämisse P_univ. (1) Die Identifikation ist wörtlich: Auf dem Einanregungssektor "
          "zweier Nähte (der irreduziblen 2D-Darstellung von z94) ist A = -iR die Nahtparität und B = B_theta^+ A B_theta ihr "
          "Bild unter der Übertragung; {A, B} = 2 cos theta, Gamma = -iAB ist der Transfergenerator 2 J_x, der Saat-Modulfluss "
          "die Transferdrehung, und B_theta implementiert A -> B -> -A genau bei pi/2. Die Leckage der Zweianregungskomposition "
          "ist |1/2 {A,B}| |1/2 [A,B]| |v_P - v_D| = |cos theta sin theta| |1 - e^(-i phi)|: z120s Leckfreiheit heißt "
          "wörtlich, dass A und B vertauschen (kein Transfer) oder antivertauschen (P_sym). (2) Allgemein für zwei Nähte "
          "(symbolischer Paarsatz, Fockraum-Querproben): leckfreier kohärenter Austausch genau dann, wenn jeder ausgetauschte "
          "Modus eine Vierteldrehung erfährt und das Phasengesetz nahtsymmetrisch ist (v_AA = v_BB, v_AB = v_BA, v_D != v_P); "
          "das Gatter ist v_D P_+ + v_P P_-, bei quadratischer Ladungsenergie mit Austauschphase phi q_k q_l. Leiter, Arten, "
          "allgemeine Ladungen und jedes phi außerhalb 2 pi Z ändern nichts; ohne Nahtsymmetrie (Vorspannung, gemischte "
          "Ladungsparität mit Uhr) gibt es keine Komposition, mit niveauabhängigen Amplituden braucht jedes Niveau die "
          "Vierteldrehung. Drei Nähte sind ein Gegenbeispiel zum wörtlichen Winkel: neben der paarweisen Vierteldrehung ist "
          "die Hadamard-Klasse F_3 leckfrei (zyklisch bei theta = 2 pi/9), mit Gatter v_P P_(S=3/2) + v_D P_(S=1/2); das "
          "gemeinsame Prinzip ist die gegenseitige Unverzerrtheit von Nahtbasis und übertragener Basis. (3) Unter den "
          "geprüften Superauswahloptionen ist die Leckfreiheit genau unter P_univ nativ: Physikalische Kompositionen erhalten "
          "die Fermionparität jeder Naht, für zwei Nähte mit ungeraden Ladungen gleichwertig die Paaruhr i^(Q_A - Q_B), die "
          "der Tomita-Operator als Delta^(-i/2) des Original-Partner-Paars liefert. Unter allen Z4-Gradierungen "
          "i^(a Q_A + b Q_B) erzwingen genau die mit b - a = 2 mod 4 die Vierteldrehung; "
          "feinere Superauswahl (Nahtladung, Nahtuhr) verbietet jede Komposition, weil der native Zustand den Doppelsektor "
          "besetzt (z93) und U dort (2,0) mit (0,2) mischt, verträglich mit z108; gröbere fordert nichts. Die Parität bleibt "
          "für jede Artenzahl erhalten, die Ladungsenergie (z105s 'auf der Schale') nur auf dem Zweianregungssektor. Unter "
          "P_univ verletzen F_3 und das zyklische Hadamard die Nahtparitäten, die Paar-Vierteldrehung nicht: P_univ wählt "
          "paarweise Komposition und tau = i. Also: nichttriviale kohärente Komposition + P_univ => nahtsymmetrisches Gesetz "
          "und {A, B} = 0 => P_sym. P_univ folgt nicht aus der globalen Univalenz (Einzeltransfers sind gesamtgerade, z106), "
          "wird aber von z105s nativen Zweinahttermen erfüllt. Verträglich mit z107 (theta = pi/3 ist spannungsfrei, aber "
          "undicht) und z94 (die Komposition ist phi- und zustandsblind, gibt also P_sym, nicht P_clock). (4) Energieroute: "
          "S(n) folgt für jedes n aus A*, wenn A* auf den Modulinhalt Tr e^(-t K') des tatsächlichen Zustands wirkt (nicht "
          "auf die geometrische Einteilchenenergie, die für Kandidaten mit gleichen Rändern blind ist). Variante B1 (S-Invarianz "
          "des Kandidateninhalts): Der Fensterquotient muss ein Monom sein, jeder Kandidat ist isospektral, ohne "
          "Ganzzahligkeit (Defekte 0,05 beim z94-Kandidaten, Kontrollen 0). Variante B2 (zustandsunabhängiger Casimir -1/24 "
          "und <T> = 0): gleiche Modulenergie, mit ganzzahligen Vielfachheiten über Gelfond Z' = Z; in den Fenstern n = 2..5 "
          "sinkt der Restklassenkern von 2n - 3 auf 1 + verborgen (0, 0, 1, 3), die verborgenen Richtungen nimmt Max-Entropie. "
          "Reell genügt B2 nicht ((1,3,5), (3,5,7)). Damit ist z107s Hindernis unter A* umgangen. K19 bleibt L: Die Komposition "
          "setzt zwei Nähte voraus. Nichts befördert, Ledger unverändert.")
ZAHLEN["laufzeit_s"] = round(time.time() - T_START, 1)
say("\nUrteil\n   " + urteil)
say(f"\nLaufzeit {time.time() - T_START:.0f} s")
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z125_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
(HERE / "ausgabe_z125.txt").write_text("\n".join(OUT) + "\n", encoding="utf-8")
