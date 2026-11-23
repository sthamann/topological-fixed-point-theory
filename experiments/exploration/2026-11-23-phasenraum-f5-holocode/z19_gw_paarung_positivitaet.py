"""Z19: Ginsparg-Wilson-Paarung, Reflexionspositivität und Gravitation bei endlichem Gitterabstand.

Offene Verbindung 3 aus H6h10. H6 paart in der eta-Basis zwei Weyl-Fermionen: eta_+ trägt eine windende
Majoranamasse m e^{i phi}, eta_- eine feste Masse m. In 4D ist ein Weyl-Fermion mit Majoranamasse ein
Majorana-Fermion lambda, lambda-quer = lambda^T C. Z7/Z9 entfernen die Verdoppler mit dem Wilson-Glied
W eta_+ eta_-, das die Paarungsphase bricht (Barriere 0,088 je Zelle). Hier:

 1. Overlap-Fassung der Paarung (euklidisch, 4D, rho = 1, Wilson-Kern r = 1).
    Symmetrische Lüscher-Drehung für Majoranafelder: delta lambda = i alpha X lambda, X = g5 (1 - D/2).
    Vier Massenformen für eta_+ (eta_- bleibt D + m (1 - D/2)):
      W   Wilson-Majorana wie Z7 (euklidisch): M(k) = [[m e^{i phi}, W], [W, m]], W = sum (1 - cos k).
      X   kovariant: m (1 - D/2) e^{i phi X} (exakt aus der Drehung, lokal, aber X nicht quantisiert).
      Xq  quantisiert: m (1 - D/2) [Pi_+ e^{i phi} Pi_+ + Pi_- e^{-i phi} Pi_-], Pi = (1 +- sign X)/2 (nichtlokal).
      GJ  Gerhold-Jansen/Z16-Form P_+ Phi P^_+ + P_- Phi* P^_-, auf den Majoranateil projiziert (lokal).
    (a) gleichförmiges Phasenpotential, (b) Index des Defekts im 2D-Querschnitt, (c) Wand: Energie einer
    Windung gegen gleichförmig, Spannung des Schnitts.
 2. Reflexionspositivität (Link und Platz) für naiv, Wilson r, Overlap; Spektraldarstellung des Overlap
    (Pole und Schnitt); gedrehte Massen (Phasen theta), wie sie der windende Defekt durchläuft.
 3. Flächendichte kappa der Verschränkung (Methode wie Z15) für chirale Regulatoren; a/l_P = 2 sqrt(kappa).

Zitiert, nicht geprüft: Lu und Shao (laut externem Bericht), Craig und Witten 2026, Lüscher 1977,
Osterwalder-Seiler 1978, Menotti-Pelissetto 1987, Gerhold-Jansen 2007.
"""
import json
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
LINES, CHECKS, BEFUNDE, ZAHLEN = [], {}, {}, {}
T0 = time.time()


def say(s=""):
    print(s, flush=True)
    LINES.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    return "ja" if cond else "NEIN"


def befund(name, val):
    BEFUNDE[name] = bool(val)
    return "ja" if val else "nein"


# ------------------------------------------------------------------------------------------- Grundlagen
gam = [np.array(x, complex) for x in (
    [[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]],
    [[0, 0, 0, -1j], [0, 0, -1j, 0], [0, 1j, 0, 0], [1j, 0, 0, 0]],
    [[0, 0, 0, -1], [0, 0, 1, 0], [0, 1, 0, 0], [-1, 0, 0, 0]],
    [[0, 0, -1j, 0], [0, 0, 0, 1j], [1j, 0, 0, 0], [0, -1j, 0, 0]])]
G4 = gam[3]
g5 = gam[0] @ gam[1] @ gam[2] @ gam[3]
I4 = np.eye(4)
C = gam[0] @ gam[2]                      # C g C^-1 = -g^T, C^T = -C (numerisch gesucht)
Ci = np.linalg.inv(C)
P_p, P_m = (I4 + g5) / 2, (I4 - g5) / 2
SIG = -1j * gam[2] @ gam[3]              # Chiralität auf der Weltfläche (x3, x4) des Defekts
M_PAAR, RHO = 0.6, 1.0


def herm_exp(H, t):
    w, U = np.linalg.eigh(H)
    return (U * np.exp(1j * t * w[..., None, :])) @ np.swapaxes(U, -1, -2).conj()


def kgrid4(n):
    k = 2 * np.pi * (np.arange(n) + 0.5) / n
    return np.array(np.meshgrid(k, k, k, k, indexing="ij")).reshape(4, -1).T


def D_naive(K):
    return 1j * np.einsum("nm,mab->nab", np.sin(K), np.array(gam))


def D_overlap(K, rho=RHO, r=1.0):
    w = r * np.sum(1 - np.cos(K), axis=1) - rho
    f = np.sum(np.sin(K) ** 2, axis=1) + w ** 2
    return rho * (I4 + (D_naive(K) + w[:, None, None] * I4) / np.sqrt(f)[:, None, None])


def antisym(Dk, Dmk):
    """Majoranateil: (D(k) + C^-1 D(-k)^T C)/2; nur er trägt zum Pfaffian bei."""
    return 0.5 * (Dk + Ci @ np.swapaxes(Dmk, -1, -2) @ C)


def logabsdet(M):
    return np.linalg.slogdet(M)[1]


say("Z19: Ginsparg-Wilson-Paarung, Reflexionspositivität und Gravitation bei endlichem Gitterabstand")
say(f"Paarungsmasse m = {M_PAAR}, Overlap mit rho = {RHO} und Wilson-Kern r = 1. Majoranaanteil über C = g1 g3.\n")

# =========================================================================================== Teil 1
say("Teil 1: Overlap-Fassung der H6-Paarung")
say("1.0 Grundlagen der Majoranadrehung")
K8 = kgrid4(8)
D8, D8m = D_overlap(K8), D_overlap(-K8)
gw = np.max(np.abs(g5 @ D8 + D8 @ g5 - D8 @ g5 @ D8))
cc = np.max(np.abs(C @ D8 @ Ci - np.swapaxes(D8m, 1, 2)))
X8 = g5 @ (I4 - D8 / 2)
xh = np.max(np.abs(X8 - np.conj(np.swapaxes(X8, 1, 2))))
x2 = np.max(np.abs(X8 @ X8 - (I4 - np.conj(np.swapaxes(D8, 1, 2)) @ D8 / 4)))
say(f"   Ginsparg-Wilson-Relation: Rest {gw:.1e}. C D(k) C^-1 = D(-k)^T: Rest {cc:.1e}.")
say(f"   X = g5 (1 - D/2) ist hermitesch (Rest {xh:.1e}) mit X^2 = 1 - D^+D/4 (Rest {x2:.1e}).")
check("Overlap erfüllt GW und Majoranabedingung C D C^-1 = D^T", gw < 1e-12 and cc < 1e-12)
check("X hermitesch, X^2 = 1 - D^+D/4", xh < 1e-12 and x2 < 1e-12)
say("   Die symmetrische Drehung lambda -> e^{i alpha X} lambda lässt die Majoranawirkung exakt invariant.")
say("   Für Dirac-Felder gibt es die asymmetrische Drehung mit g5^ = g5 (1 - D), g5^2 = 1 (Z16). Für Majorana")
say("   sind lambda und lambda-quer verbunden; übrig bleibt nur X, und X hat kein ganzzahliges Spektrum.")
absX = np.abs(np.linalg.eigvalsh(X8))
J4dev = np.abs(np.exp(2j * np.pi * absX) - 1).max(axis=1)
frac = float(np.mean(J4dev > 0.01))
say(f"   |X| auf dem Gitter 8^4 (Mittelpunkte): {absX.min():.3f} bis {absX.max():.3f}; |X| -> 1 an der physischen Ecke,")
say("   |X| -> 0 an den Verdopplerecken (dort D = 2).")
say(f"   Die Gitter-Vierteldrehung J = e^{{i pi X/2}} hat J^4 = e^{{2 pi i X}} != 1: max |J^4 - 1| = {J4dev.max():.3f},")
say(f"   auf {100 * frac:.0f} % der Gitterimpulse 8^4 über 0,01. J^4 = 1 gilt nur für k -> 0.")
befund("J der Majorana-Overlapdrehung hat nicht Ordnung 4 (X nicht quantisiert)", J4dev.max() > 0.5)
ZAHLEN["X_spektrum"] = [float(absX.min()), float(absX.max())]
ZAHLEN["J4_abweichung_max"] = float(J4dev.max())
ZAHLEN["J4_anteil_ueber_0.01"] = frac

# Lokalität von X gegen sign X: Richtungsgrenzwerte an einer Verdopplerecke
corner = np.array([np.pi, 0.0, 0.0, 0.0])


def at(k):
    X = g5 @ (I4 - D_overlap(np.atleast_2d(k)) / 2)
    w, U = np.linalg.eigh(X)
    return X[0], ((U * np.sign(w)[:, None, :]) @ np.conj(np.swapaxes(U, 1, 2)))[0]


jump = {}
for delta in (1e-2, 1e-4):
    vals = [at(corner + delta * e) for e in np.eye(4)]
    jump[delta] = (max(np.linalg.norm(a[0] - b[0]) for a in vals for b in vals),
                   max(np.linalg.norm(a[1] - b[1]) for a in vals for b in vals))
say("   Lokalität: Annäherung an die Verdopplerecke (pi,0,0,0) aus den vier Achsenrichtungen, Abstand 1e-2/1e-4:")
say(f"   größte Differenz von X {jump[1e-2][0]:.1e}/{jump[1e-4][0]:.1e} (stetig), von sign X {jump[1e-2][1]:.2f}/{jump[1e-4][1]:.2f}.")
say("   sign X springt an den Verdopplerecken: Sein Ortsraumkern fällt nur wie eine Potenz, es ist nichtlokal.")
say("   X ist analytisch (exponentiell lokal), weil D_W - rho keine Nullstelle hat.")
check("X stetig an den Verdopplerecken, sign X springt", jump[1e-4][0] < 1e-3 and jump[1e-4][1] > 0.5)
befund("Quantisierte Majoranadrehung (sign X) ist nichtlokal", jump[1e-4][1] > 0.5)
ZAHLEN["lokalitaet_sprung"] = {str(d): [float(v[0]), float(v[1])] for d, v in jump.items()}

# ------------------------------------------------------------------------------------------- 1a
say("\n1a Gleichförmiges Phasenpotential V(phi) = -1/2 <ln |det D_eff|> je Gitterzelle (8^4 Impulse)")


def V_uniform(form, phi, K=K8, m=M_PAAR):
    Km = -K
    if form == "W":
        W = np.sum(1 - np.cos(K), axis=1)
        Mk = np.zeros((len(K), 2, 2), complex)
        Mk[:, 0, 0], Mk[:, 1, 1], Mk[:, 0, 1], Mk[:, 1, 0] = m * np.exp(1j * phi), m, W, W
        Dt = np.einsum("ij,nab->niajb", np.eye(2), D_naive(K)).reshape(-1, 8, 8) + (
            np.einsum("nij,ab->niajb", Mk, P_p) + np.einsum("nij,ab->niajb", Mk.conj(), P_m)).reshape(-1, 8, 8)
        return -0.5 * logabsdet(Dt).mean(), 0.0
    Dk, Dmk = D_overlap(K), D_overlap(Km)

    def mass(D):
        X = g5 @ (I4 - D / 2)
        if form == "X":
            return m * (I4 - D / 2) @ herm_exp(X, phi)
        if form == "Xq":
            w, U = np.linalg.eigh(X)
            return m * (I4 - D / 2) @ ((U * np.exp(1j * phi * np.sign(w))[:, None, :]) @ np.conj(np.swapaxes(U, 1, 2)))
        if form in ("GJ", "GJ_Dirac"):
            return m * (np.cos(phi) * I4 + 1j * np.sin(phi) * g5) @ (I4 - D / 2)
    Mk, Mmk = mass(Dk), mass(Dmk)
    full = Dk + Mk
    eff = full if form == "GJ_Dirac" else antisym(full, Dmk + Mmk)
    base = logabsdet(Dk + m * (I4 - Dk / 2)).mean()          # eta_-, phasenunabhängig
    fac = 1.0 if form == "GJ_Dirac" else 0.5
    return -fac * (logabsdet(eff).mean() + (0 if form == "GJ_Dirac" else base)), float(np.max(np.abs(eff - full)))


phis = np.linspace(0, 2 * np.pi, 9)
pot, asy = {}, {}
for form in ("W", "X", "Xq", "GJ", "GJ_Dirac"):
    vals = [V_uniform(form, p) for p in phis]
    v = np.array([a for a, _ in vals])
    pot[form] = v - v.min()
    asy[form] = max(b for _, b in vals)
names = {"W": "Wilson-Majorana (Z7, euklidisch)", "X": "Overlap, kovariant e^{i phi X}",
         "Xq": "Overlap, quantisiert sign X", "GJ": "Overlap, GJ-Form auf Majorana projiziert",
         "GJ_Dirac": "Overlap, GJ-Form als Dirac (Z16)"}
for form in pot:
    say(f"   {names[form]:42s} V - V_min = " + " ".join(f"{x:.4f}" for x in pot[form])
        + f"  (Spannweite {np.ptp(pot[form]):.1e})")
say(f"   phi = 0, pi/4, ..., 2 pi. Majoranaprojektion ändert die X-Form um {asy['X']:.1e} (automatisch konsistent),")
say(f"   die GJ-Form um {asy['GJ']:.2f}: Die Dirac-Form ist keine Majoranamasse.")
check("X-Masse ist von selbst eine Majoranamasse (C-antisymmetrisch)", asy["X"] < 1e-12)
check("Wilson-Majorana: Barriere, Minimum bei pi, Periode 2 pi", np.ptp(pot["W"]) > 0.05 and abs(phis[pot["W"].argmin()] - np.pi) < 1e-9
      and abs(pot["W"][0] - pot["W"][-1]) < 1e-10)
check("X- und Xq-Form flach (Spannweite < 1e-12)", np.ptp(pot["X"]) < 1e-12 and np.ptp(pot["Xq"]) < 1e-12)
check("Z16 reproduziert: Dirac-GJ flach", np.ptp(pot["GJ_Dirac"]) < 1e-12)
befund("Drehidentität aus Z16 überträgt sich nicht auf die Majorana-GJ-Form", np.ptp(pot["GJ"]) > 1e-3)
say(f"   Die projizierte GJ-Form behält eine Barriere {np.ptp(pot['GJ']):.4f} mit Periode pi (Minima 0 und pi):")
say(f"   etwa {np.ptp(pot['W']) / np.ptp(pot['GJ']):.0f} mal kleiner als Wilson, aber nicht null.")

# Hamiltonsche Gegenprobe wie Z9
sx = np.array([[0, 1], [1, 0]], complex)
sy = np.array([[0, -1j], [1j, 0]])
sz = np.diag([1.0, -1.0]).astype(complex)
isy = 1j * sy


def e_density_ham(phi0, m=M_PAAR, n=12):
    k = 2 * np.pi * (np.arange(n) + 0.5) / n
    Kx, Ky, Kz = [a.ravel() for a in np.meshgrid(k, k, k, indexing="ij")]
    hs = np.sin(Kx)[:, None, None] * sx + np.sin(Ky)[:, None, None] * sy + np.sin(Kz)[:, None, None] * sz
    W = 3 - np.cos(Kx) - np.cos(Ky) - np.cos(Kz)
    Mm = np.zeros((len(W), 2, 2), complex)
    Mm[:, 0, 0], Mm[:, 1, 1], Mm[:, 0, 1], Mm[:, 1, 0] = m * np.exp(1j * phi0), m, W, W
    Dm = np.einsum("nij,ab->niajb", Mm, isy).reshape(-1, 4, 4)
    hb = np.einsum("ij,nab->niajb", np.eye(2), hs).reshape(-1, 4, 4)
    H = np.block([[hb, Dm], [np.conj(np.swapaxes(Dm, 1, 2)), np.swapaxes(hb, 1, 2)]])
    E = np.linalg.eigvalsh(H)
    return float(-0.5 * np.where(E > 0, E, 0).sum(axis=1).mean())


vh = np.array([e_density_ham(p) for p in (0.0, np.pi)])
say(f"   Gegenprobe Hamiltonsch (Z9): Barriere {vh[0] - vh[1]:.4f} (Z9: 0,0882); euklidisch {np.ptp(pot['W']):.4f}.")
check("Z9-Barriere reproduziert (0,0882)", abs(vh[0] - vh[1] - 0.0882) < 0.0005)
ZAHLEN["phasenpotential"] = {f: [float(x) for x in pot[f]] for f in pot}
ZAHLEN["barriere_hamilton_z9"] = float(vh[0] - vh[1])

# ------------------------------------------------------------------------------------------- 1b
say("\n1b Index des verdrillten Defekts: 2D-Querschnitt L x L (offen) bei k3 = k4 = 0")
L2 = 16
nsite = L2 * L2


def transverse(L, r=1.0):
    idx = lambda x, y: x * L + y
    n = L * L
    N = np.zeros((4 * n, 4 * n), complex)
    Wop = np.zeros((n, n))
    for x in range(L):
        for y in range(L):
            i = idx(x, y)
            Wop[i, i] += 2 * r
            for mu, (dx, dy) in enumerate(((1, 0), (0, 1))):
                if x + dx < L and y + dy < L:
                    j = idx(x + dx, y + dy)
                    N[4 * i:4 * i + 4, 4 * j:4 * j + 4] += gam[mu] / 2
                    N[4 * j:4 * j + 4, 4 * i:4 * i + 4] -= gam[mu] / 2
                    Wop[i, j] -= r / 2
                    Wop[j, i] -= r / 2
    return N, Wop


Nt, Wt = transverse(L2)
G5b = np.kron(np.eye(nsite), g5)
Hw = G5b @ (Nt + np.kron(Wt - RHO * np.eye(nsite), I4))
Hw = (Hw + Hw.conj().T) / 2
wH, UH = np.linalg.eigh(Hw)
Dq = RHO * (np.eye(4 * nsite) + G5b @ ((UH * np.sign(wH)) @ UH.conj().T))
Ib = np.eye(4 * nsite)
Cb, Cbi = np.kron(np.eye(nsite), C), np.kron(np.eye(nsite), Ci)
Tq = lambda A: Cbi @ A.T @ Cb
gwq = np.max(np.abs(G5b @ Dq + Dq @ G5b - Dq @ G5b @ Dq))
say(f"   Overlap auf dem offenen Querschnitt: GW-Rest {gwq:.1e}, kleinster |H_W| = {np.abs(wH).min():.3f}.")
cq = (L2 - 1) / 2
xs, ys = np.meshgrid(np.arange(L2), np.arange(L2), indexing="ij")
phi_q = np.arctan2(ys - cq, xs - cq).ravel()
rad = np.hypot(xs - cq, ys - cq).ravel()
prof = np.tanh(rad / 1.5)
Xq_ = G5b @ (Ib - Dq / 2)
Xq_ = (Xq_ + Xq_.conj().T) / 2
wX, UX = np.linalg.eigh(Xq_)
qX = np.where(np.abs(wX) > 1e-9, np.sign(wX), 0.0)
Pip, Pim = (UX * ((1 + qX) / 2)) @ UX.conj().T, (UX * ((1 - qX) / 2)) @ UX.conj().T
g5hat = G5b @ (Ib - Dq)
Phat_p, Phat_m = (Ib + g5hat) / 2, (Ib - g5hat) / 2
PPb, PMb = np.kron(np.eye(nsite), P_p), np.kron(np.eye(nsite), P_m)
Sq = np.kron(np.diag(np.sqrt(prof)), I4)


def E_ph(n):
    return np.kron(np.diag(prof * np.exp(1j * n * phi_q)), I4)


def op_defect(form, n):
    m = M_PAAR
    if form == "X":
        Phi = np.kron(np.diag(n * phi_q), I4)
        K = herm_exp((Phi @ Xq_ + Xq_ @ Phi) / 2, 0.5)
        return Dq + m * Tq(K) @ Sq @ (Ib - Dq / 2) @ Sq @ K
    if form == "Xq":
        E = E_ph(n)
        D = Dq + m * (Ib - Dq / 2) @ (Pip @ E @ Pip + Pim @ E.conj() @ Pim)
        return 0.5 * (D + Tq(D))
    if form == "GJ":
        E = E_ph(n)
        D = Dq + m * (PPb @ E @ Phat_p + PMb @ E.conj() @ Phat_m)
        return 0.5 * (D + Tq(D))


def op_wilson_pair(npl, nmi):
    m = M_PAAR
    Mp = m * (E_ph(npl) @ PPb + E_ph(npl).conj() @ PMb)
    Mm = m * (E_ph(nmi) @ PPb + E_ph(nmi).conj() @ PMb)
    Wb = np.kron(Wt, I4)
    D = np.block([[Nt + Mp, Wb], [Wb, Nt + Mm]])
    perm = np.array([[f * 4 * nsite + 4 * s + a for f in range(2) for a in range(4)] for s in range(nsite)]).ravel()
    return D[np.ix_(perm, perm)]


def chiral_count(Deff, ns, R=3.5, tol=0.25):
    """Kernnahe rechte Nullvektoren je Sektor von Sigma34; netto = n_+ - n_- (Majorana-Weyl-Moden)."""
    wS, US = np.linalg.eigh(SIG)
    res, lows = {}, {}
    core = np.repeat(rad <= R, 4 * ns)
    for s in (1, -1):
        B = np.kron(np.eye(nsite * ns), US[:, wS * s > 0])
        _, sv, vh = np.linalg.svd(B.conj().T @ Deff @ B)
        V = (B @ vh.conj().T)[:, sv < tol]
        res[s] = int(np.sum(np.sum(np.abs(V[core]) ** 2, axis=0) > 0.5))
        lows[s] = float(np.sort(sv)[0])
    return res[1] - res[-1], res, lows


index_tab = {}
for form in ("X", "Xq", "GJ"):
    row = {}
    for n in (1, 0, 2, -1):
        net, res, lows = chiral_count(op_defect(form, n), 1)
        row[n] = net
    index_tab[form] = row
    say(f"   {names[form]:42s} Windung 1/0/2/-1: netto {row[1]}/{row[0]}/{row[2]}/{row[-1]} chirale MW-Moden am Kern")
wp = {}
for pair in ((1, 0), (1, 1), (0, 0)):
    net, res, lows = chiral_count(op_wilson_pair(*pair), 2)
    wp[pair] = net
say(f"   Wilson-Paar wie Z7 (eta_+, eta_-) = (1,0)/(1,1)/(0,0): netto {wp[(1, 0)]}/{wp[(1, 1)]}/{wp[(0, 0)]}")
ok_idx = all(index_tab[f] == {1: 1, 0: 0, 2: 2, -1: -1} for f in index_tab)
check("Index = Windung für X, Xq, GJ (1, 0, 2, -1)", ok_idx)
check("Wilson-Paar wie Z7: 1, 2, 0", wp == {(1, 0): 1, (1, 1): 2, (0, 0): 0})
say("   Verdrillter Defekt (eta_+ windet, eta_- nicht): eine chirale Majorana-Weyl-Mode je Paar, für jede")
say("   chirale Form. Die Verdoppler tragen beim Overlap nichts bei (D = 2 dort, die Masse fällt heraus).")
# Schnittdiagnose der X-Form
Dx1 = op_defect("X", 1)
wS, US = np.linalg.eigh(SIG)
B = np.kron(np.eye(nsite), US[:, wS > 0])
u_, sv_, vh_ = np.linalg.svd(B.conj().T @ Dx1 @ B)
V_ = B @ vh_.conj().T
cutm = np.repeat((np.abs(np.abs(phi_q) - np.pi) < 0.35) & (rad > 2), 4)
order = np.argsort(sv_)
cut_w = [float(np.sum(np.abs(V_[cutm, j]) ** 2)) for j in order[:2]]
say(f"   X-Form: Die Kernmode trägt {cut_w[0]:.2f} ihres Gewichts auf dem Schnitt (phi = pi); der nächste Zustand")
say(f"   (Singulärwert {sv_[order[1]]:.2f} statt Lücke ~0,57) liegt zu {cut_w[1]:.2f} auf dem Schnitt: Der Schnitt bindet Zustände.")
befund("X-Form: Verzweigungsschnitt ist physikalisch sichtbar (gebundene Zustände)", cut_w[1] > 0.3)
ZAHLEN["index"] = {f: {str(k): v for k, v in index_tab[f].items()} for f in index_tab}
ZAHLEN["index_wilson_paar"] = {str(k): v for k, v in wp.items()}

# ------------------------------------------------------------------------------------------- 1c
say("\n1c Wand: freie Energie einer Phasenwindung entlang x1 (Ring der Länge L, Querimpulse k2,k3,k4)")


def mixed_ops(L, kp):
    k1 = 2 * np.pi * (np.arange(L) + 0.5) / L
    F = np.exp(1j * np.outer(np.arange(L), k1)) / np.sqrt(L)
    sp, Wp = np.sin(kp), np.sum(1 - np.cos(kp))
    Dk, Nk, Wk = [], [], []
    for q in k1:
        kin = 1j * np.sin(q) * gam[0] + sum(1j * sp[i] * gam[i + 1] for i in range(3))
        w = Wp + 1 - np.cos(q)
        f = np.sin(q) ** 2 + np.sum(sp ** 2) + (w - RHO) ** 2
        Dk.append(RHO * (I4 + (kin + (w - RHO) * I4) / np.sqrt(f)))
        Nk.append(kin)
        Wk.append(w)
    tox = lambda Bl: np.einsum("xk,kab,yk->xayb", F, np.array(Bl), F.conj()).reshape(4 * L, 4 * L)
    return tox(Dk), tox(Nk), (F * np.array(Wk)) @ F.conj().T


def F_config(form, L, phi, nperp):
    """-1/2 <ln|det D_eff>> über Querimpulse, summiert über den Ring (eta_+; bei W das ganze Paar)."""
    kk = 2 * np.pi * (np.arange(nperp) + 0.5) / nperp
    I = np.eye(4 * L)
    Cb_, Cbi_ = np.kron(np.eye(L), C), np.kron(np.eye(L), Ci)
    T_ = lambda A: Cbi_ @ A.T @ Cb_
    Phi = np.kron(np.diag(phi), I4)
    E = np.kron(np.diag(np.exp(1j * phi)), I4)
    G5_ = np.kron(np.eye(L), g5)
    PP_, PM_ = np.kron(np.eye(L), P_p), np.kron(np.eye(L), P_m)
    acc = 0.0
    for a in kk:
        for b in kk:
            for c in kk:
                kp = np.array([a, b, c])
                if form == "W":
                    _, N, Wx = mixed_ops(L, kp)
                    Wb = np.kron(Wx, I4)
                    Mp = M_PAAR * (E @ PP_ + E.conj() @ PM_)
                    acc += logabsdet(np.block([[N + Mp, Wb], [Wb, N + M_PAAR * I]]))
                    continue
                o = {}
                for sg in (1, -1):
                    D, _, _ = mixed_ops(L, sg * kp)
                    X = G5_ @ (I - D / 2)
                    o[sg] = (D, (X + X.conj().T) / 2)
                if form == "X":
                    Xp = lambda s: (Phi @ o[s][1] + o[s][1] @ Phi) / 2
                    acc += logabsdet(o[1][0] + M_PAAR * T_(herm_exp(Xp(-1), 0.5)) @ (I - o[1][0] / 2) @ herm_exp(Xp(1), 0.5))
                    continue

                def mass(s):
                    D, X = o[s]
                    if form == "Xq":
                        w, U = np.linalg.eigh(X)
                        q = np.where(np.abs(w) > 1e-9, np.sign(w), 0.0)
                        Pp, Pm = (U * ((1 + q) / 2)) @ U.conj().T, (U * ((1 - q) / 2)) @ U.conj().T
                        return M_PAAR * (I - D / 2) @ (Pp @ E @ Pp + Pm @ E.conj() @ Pm)
                    gh = G5_ @ (I - D)
                    return M_PAAR * (PP_ @ E @ (I + gh) / 2 + PM_ @ E.conj() @ (I - gh) / 2)
                acc += logabsdet(0.5 * ((o[1][0] + mass(1)) + T_(o[-1][0] + mass(-1))))
    return -0.5 * acc / nperp ** 3


ref_phase = {"W": np.pi, "X": 0.0, "Xq": 0.0, "GJ": 0.0}
NP = 6
spread = {f: {} for f in ref_phase}
t1 = time.time()
for f in ref_phase:
    for L in (8, 16, 32):
        lin = 2 * np.pi * (np.arange(L) + 0.5) / L
        base = F_config(f, L, np.full(L, ref_phase[f]), NP)
        spread[f][L] = F_config(f, L, (lin + ref_phase[f]) % (2 * np.pi) if f != "X" else lin, NP) - base
for f in spread:
    s = spread[f]
    say(f"   {names[f]:42s} gleichmäßig verteilte Windung: dF(L=8/16/32) = {s[8]:.4f}/{s[16]:.4f}/{s[32]:.4f};"
        f" L dF = {8 * s[8]:.3f}/{16 * s[16]:.3f}/{32 * s[32]:.3f}")
say(f"   (Querimpulse {NP}^3, je Einheit Wandquerschnitt; Rechenzeit {time.time() - t1:.0f} s)")
say("   Wilson und GJ: dF wächst mit L (Potentialkosten): die Windung zieht sich zu einer Wand zusammen.")
say("   Xq: L dF ist konstant, dF ~ 1/L, reine Gradientenenergie: keine Wand.")
say("   X: dF bleibt bei L -> unendlich endlich: Die Kosten sitzen am Sprung phi = 2 pi -> 0, dem Schnitt.")
say("   Ihre Höhe hängt davon ab, zwischen welchen Zweigwerten der Sprung liegt (e^{2 pi i X} != 1).")
check("Wilson: Windungskosten wachsen linear mit L", 1.9 < spread["W"][32] / spread["W"][16] < 2.1)
check("Xq: Windungskosten fallen wie 1/L (keine Wand)", abs(32 * spread["Xq"][32] / (16 * spread["Xq"][16]) - 1) < 0.1)
befund("GJ-Form: Wand bleibt (Kosten wachsen mit L)", spread["GJ"][32] / spread["GJ"][16] > 1.6)
# Schnittspannung der X-Form: reiner Schnitt (phi = 0 | 2 pi, physikalisch dieselbe Masse)
Lc = 12
cut_conf = np.where(np.arange(Lc) < Lc // 2, 0.0, 2 * np.pi)
sig_cut = {}
for npp in (6, 8, 10):
    sig_cut[npp] = (F_config("X", Lc, cut_conf, npp) - F_config("X", Lc, np.zeros(Lc), npp)) / 2
say(f"   X-Form, reiner Schnitt (phi = 0 links, 2 pi rechts; e^{{i phi}} überall 1): Spannung je Schnitt "
    + "/".join(f"{sig_cut[n]:.5f}" for n in sig_cut) + " (Querimpulse 6^3/8^3/10^3).")
befund("X-Form: Schnitt zwischen gleichen Vakua kostet Spannung", sig_cut[10] > 1e-3)
# Wandspannungen Wilson und GJ über Kinkbreiten
Lk = 32
x = np.arange(Lk) + 0.5


def step(xx, x0, w):
    s = 1 / (1 + np.exp(-(xx - x0) / w))
    s0, s1 = 1 / (1 + np.exp(x0 / w)), 1 / (1 + np.exp(-(Lk - x0) / w))
    return (s - s0) / (s1 - s0)


kink = {"W": {}, "GJ": {}}
bW = F_config("W", Lk, np.full(Lk, np.pi), NP)
bG = F_config("GJ", Lk, np.zeros(Lk), NP)
for w in (0.7, 1.5, 3.0, 5.0):
    kink["W"][w] = F_config("W", Lk, np.pi + 2 * np.pi * step(x, Lk / 2, w), NP) - bW
    kink["GJ"][w] = F_config("GJ", Lk, np.pi * step(x, Lk / 4, w) + np.pi * step(x, 3 * Lk / 4, w), NP) - bG
sW, sG = min(kink["W"].values()), min(kink["GJ"].values()) / 2
say("   Kinkbreiten 0,7/1,5/3/5 (L = 32; bei 0,7 springt die Phase höchstens 2,2 je Link, Windung noch definiert):")
say("   Wilson dF = " + "/".join(f"{v:.4f}" for v in kink["W"].values())
    + "; GJ (zwei pi-Kinks) = " + "/".join(f"{v:.4f}" for v in kink["GJ"].values()))
say(f"   Die Wände sind gitterscharf (dF fällt bis zur kleinsten Breite). Obere Schranken je Einheit Wandquerschnitt:")
say(f"   Wilson {sW:.4f} (eine Wand je Windung), GJ {sG:.4f} je pi-Wand (zwei je Windung), X-Schnitt {sig_cut[10]:.4f}.")
say("   Xq: keine Wand (nichtlokal).")
ZAHLEN["wand"] = {"verteilt": {f: {str(L): float(v) for L, v in spread[f].items()} for f in spread},
                  "sigma_W": float(sW), "sigma_GJ_je_pi_wand": float(sG),
                  "sigma_X_schnitt": {str(k): float(v) for k, v in sig_cut.items()},
                  "kink": {f: {str(w): float(v) for w, v in kink[f].items()} for f in kink}}
say("   Ergebnis 1c: Lokal und eindeutig (W, GJ) gibt es eine Wand; lokal und flach (X) wird der Schnitt")
say("   selbst zur Wand (kleinere Spannung); wandfrei ist nur die nichtlokale Form Xq.")

# =========================================================================================== Teil 2
say("\nTeil 2: Reflexionspositivität (OS) bei endlichem Gitterabstand, freies Feld, fester Raumimpuls p")
say("   OS-Form für lineare Observablen F = sum_t psi-quer(t) f_t (t >= 1, Link) bzw. t >= 0 (Platz):")
say("   M_st = -g4 S(1 - s - t) (Link), -g4 S(-s - t) (Platz); dazu die psi-Seite. RP <=> M >= 0.")
say("   Das Vorzeichen ist eine Konvention der Form; es ist einmal an Wilson r = 1 festgelegt.")


def D_time(kind, z, ps, m=0.0, r=1.0, rho=1.0, th=0.0):
    z = np.atleast_1d(np.asarray(z, complex))
    s4, c4 = (z - 1 / z) / (2j), (z + 1 / z) / 2
    ss, Wsp = np.sin(ps), np.sum(1 - np.cos(ps))
    kin = 1j * s4[:, None, None] * G4 + sum(1j * ss[i] * gam[i] for i in range(3))[None]
    if kind == "naiv":
        return kin + m * I4[None]
    if kind == "wilson":
        return kin + (r * (Wsp + 1 - c4))[:, None, None] * I4[None] + m * (np.cos(th) * I4 + 1j * np.sin(th) * g5)[None]
    w = r * (Wsp + 1 - c4) - rho
    f = s4 ** 2 + np.sum(ss ** 2) + w ** 2
    D = rho * (I4[None] + (kin + w[:, None, None] * I4[None]) / np.sqrt(f)[:, None, None])
    if kind == "overlap":
        return (1 - m / (2 * rho)) * D + m * I4[None]
    X = g5[None] @ (I4[None] - D / (2 * rho))
    if kind == "ov_X":
        return D + m * (I4[None] - D / (2 * rho)) @ herm_exp((X + np.conj(np.swapaxes(X, 1, 2))) / 2, th)
    if kind == "ov_Xq":
        w_, U_ = np.linalg.eigh(X)
        R = (U_ * np.exp(1j * th * np.sign(w_))[:, None, :]) @ np.conj(np.swapaxes(U_, 1, 2))
        return D + m * (I4[None] - D / (2 * rho)) @ R
    if kind == "ov_GJ":
        return D + m * (np.cos(th) * I4 + 1j * np.sin(th) * g5)[None] @ (I4[None] - D / (2 * rho))


def S_tau(kind, ps, taus, N=4096, **kw):
    p4 = 2 * np.pi * (np.arange(N) + 0.5) / N
    Si = np.linalg.inv(D_time(kind, np.exp(1j * p4), ps, **kw))
    return np.einsum("tn,nab->tab", np.exp(1j * np.outer(taus, p4)), Si) / N


def os_test(kind, ps, nb=4, link=True, sigma=-1, **kw):
    off = 1 if link else 0
    taus = np.arange(-(2 * nb + 1), 2 * nb + 2)
    S = dict(zip(taus, S_tau(kind, ps, taus, **kw)))
    idx = range(off, off + nb)
    out = []
    for A in (np.block([[sigma * G4 @ S[off - s - t] for t in idx] for s in idx]),
              np.block([[-sigma * (S[s + t - off] @ G4).T for t in idx] for s in idx])):
        aher = np.max(np.abs(A - A.conj().T)) / np.max(np.abs(A))
        H = (A + A.conj().T) / 2
        nrm = np.array([np.linalg.norm(H[4 * i:4 * i + 4, 4 * i:4 * i + 4]) for i in range(nb)])
        d = np.repeat(np.sqrt(np.maximum(nrm, 1e-14 * nrm.max())), 4)      # Skala je Zeitschicht
        ev = np.linalg.eigvalsh(H / np.outer(d, d))
        out.append((aher, float(ev.min() / ev.max())))
    return max(a for a, _ in out), min(b for _, b in out)


PS = [np.array(p) for p in ((0.3, 0.0, 0.0), (1.0, 0.5, 0.2), (2.5, 2.0, 0.3), (np.pi, 0.4, 0.0),
                            (2.8, 2.9, 3.0), (0.1, 0.05, 0.02))]
regs = [("naiv", dict(m=0.3)), ("wilson", dict(r=1.0)), ("wilson", dict(r=0.5)), ("wilson", dict(r=2.0)),
        ("overlap", dict(rho=1.0)), ("overlap", dict(rho=1.0, m=0.3)), ("overlap", dict(rho=0.5)),
        ("overlap", dict(rho=1.5))]
say("\n2a Normierte OS-Matrix (4 Zeitschichten, 6 Raumimpulse): kleinster Eigenwert, Link / Platz")
rp = {}
for kind, kw in regs:
    lk = [os_test(kind, p, link=True, **kw) for p in PS]
    st = [os_test(kind, p, link=False, **kw) for p in PS]
    lab = kind + " " + ",".join(f"{k}={v}" for k, v in kw.items())
    rp[lab] = {"link_min": min(b for _, b in lk), "platz_min": min(b for _, b in st),
               "antiherm": max(a for a, _ in lk + st)}
    say(f"   {lab:26s} Link {rp[lab]['link_min']:+.2e}   Platz {rp[lab]['platz_min']:+.2e}   (antihermitesch {rp[lab]['antiherm']:.0e})")
flip = os_test("wilson", PS[1], sigma=+1, r=1.0)[1]
TOL = -1e-9
check("Wilson r = 1: Link-RP, Platz nicht; Gegenvorzeichen verletzt", rp["wilson r=1.0"]["link_min"] > TOL
      and rp["wilson r=1.0"]["platz_min"] < -0.1 and flip < -0.1)
check("naiv und Wilson r = 0,5 Link-RP; Wilson r = 2 verletzt", rp["naiv m=0.3"]["link_min"] > TOL
      and rp["wilson r=0.5"]["link_min"] > TOL and rp["wilson r=2.0"]["link_min"] < -0.05)
ov_ok = all(rp[l]["link_min"] > TOL for l in rp if l.startswith("overlap"))
check("Overlap bei theta = 0 (rho 0,5/1/1,5, m 0/0,3): Link-OS-Matrix positiv", ov_ok)
say("   Wilson r <= 1 und naiv sind Link-reflexionspositiv, r = 2 nicht; Platzreflexion verletzt alle")
say("   (Literatur: Lüscher 1977, Osterwalder-Seiler 1978, Menotti-Pelissetto 1987; zitiert). Naiv: T hat")
say("   negative Eigenwerte, T^2 > 0 (wie Z6).")

say("\n2b Spektraldarstellung des Overlap-Propagators in w = e^{-E}: Pole und Schnitt (0, z*)")
say("   D_W^+ D_W = f(z) = A - B'(z + 1/z) ist skalar; sqrt(f) hat Verzweigungspunkte bei 0 und z*.")
say("   Zwischen ihnen liegt ein Kontinuum von Transfereigenwerten (Zustände auf Gitterskala).")


def spectral(ps, m=0.0, rho=1.0, ngl=240):
    ss, Wsp = np.sin(ps), np.sum(1 - np.cos(ps))
    Bp = Wsp + 1 - rho
    cst = (1 + np.sum(ss ** 2) + Bp ** 2) / (2 * Bp)
    zst = cst - np.sign(cst) * np.sqrt(cst ** 2 - 1)
    Sinv = lambda w: D_time("overlap", 1 / np.asarray(w, complex), ps, m=m, rho=rho)
    # Schnittgewicht rho(lambda) = [S(1/(l - i0)) - S(1/(l + i0))]/(2 pi i), Substitution l = z* sin^2(pi u/2)
    u, wu = np.polynomial.legendre.leggauss(ngl)
    u = (u + 1) / 2
    wu = wu / 2
    lam = zst * np.sin(np.pi * u / 2) ** 2
    jac = zst * np.pi * np.sin(np.pi * u / 2) * np.cos(np.pi * u / 2)
    eps = 1e-12 * max(abs(zst), 1e-3)
    rho_c = (np.linalg.inv(Sinv(lam - 1j * eps)) - np.linalg.inv(Sinv(lam + 1j * eps))) / (2j * np.pi)
    Acut = -G4[None] @ rho_c
    Hc = (Acut + np.conj(np.swapaxes(Acut, 1, 2))) / 2
    ev = np.linalg.eigvalsh(Hc)
    cut_min = float((ev.min(axis=1) / np.maximum(ev.max(axis=1), 1e-300)).min())
    # Pole auf der reellen Achse zwischen z* und 1 (und für z* < 0 auch in (0, 1))
    grid = np.concatenate([np.linspace(max(zst, 0.0) + 1e-6, 0.999, 4000), np.linspace(-0.999, min(zst, 0.0) - 1e-6, 4000)])
    smin = np.linalg.svd(Sinv(grid), compute_uv=False)[:, -1]
    poles = []
    for i in list(range(1, 3999)) + list(range(4001, 7999)):
        if smin[i] < smin[i - 1] and smin[i] < smin[i + 1] and smin[i] < 1e-2:
            a, b = grid[i - 1], grid[i + 1]
            for _ in range(60):
                c1, c2 = a + (b - a) / 3, b - (b - a) / 3
                f1 = np.linalg.svd(Sinv(c1), compute_uv=False)[0, -1]
                f2 = np.linalg.svd(Sinv(c2), compute_uv=False)[0, -1]
                a, b = (a, c2) if f1 < f2 else (c1, b)
            poles.append((a + b) / 2)
    res = []
    for w0 in poles:
        rr = 0.3 * min(abs(w0 - zst), 1 - abs(w0), abs(w0), 0.05)
        th = 2 * np.pi * (np.arange(64) + 0.5) / 64
        wc = w0 + rr * np.exp(1j * th)
        Sv = np.linalg.inv(Sinv(wc))
        res.append(np.einsum("n,nab->ab", 1j * rr * np.exp(1j * th), Sv) / 64 / 1j)
    pole_min = 1.0
    for R in res:
        A = -G4 @ R
        H = (A + A.conj().T) / 2
        e = np.linalg.eigvalsh(H)
        pole_min = min(pole_min, e.min() / e.max())
    # Rekonstruktion von G(-tau) gegen FFT
    taus = np.arange(1, 7)
    G_fft = S_tau("overlap", ps, -taus, m=m, rho=rho)
    rec = []
    for t in taus:
        g = sum(w0 ** (t - 1) * R for w0, R in zip(poles, res))
        g = g + np.sign(zst) * np.einsum("n,nab->ab", wu * jac * lam ** (t - 1), rho_c)
        rec.append(g)
    err = float(np.max(np.abs(np.array(rec) - G_fft)) / np.max(np.abs(G_fft)))
    return dict(zstar=float(zst), E_cut=float(-np.log(abs(zst))), poles=[float(p) for p in poles],
                pole_min=float(pole_min), cut_min=cut_min, rek_fehler=err)


spec = {}
ok_spec = True
for rho_, m_ in ((1.0, 0.0), (1.0, 0.3), (0.5, 0.0), (1.5, 0.0)):
    for p in PS[:5]:
        sp_ = spectral(p, m=m_, rho=rho_)
        spec[f"rho={rho_},m={m_},p=({p[0]:.2f},{p[1]:.2f},{p[2]:.2f})"] = sp_
        ok_spec &= sp_["rek_fehler"] < 1e-6 and sp_["cut_min"] > -1e-9 and sp_["pole_min"] > -1e-9
for key in list(spec)[:5] + list(spec)[10:11] + list(spec)[15:16]:
    s_ = spec[key]
    say(f"   {key:32s} z* = {s_['zstar']:+.3f} (E_Schnitt {s_['E_cut']:.2f}), Pol {np.round(s_['poles'], 3)},"
        f" Gewichte >= 0: {min(s_['pole_min'], s_['cut_min']) > -1e-9}, Rekonstruktion {s_['rek_fehler']:.0e}")
check("Overlap-Spektrum: Pole + Schnitt rekonstruieren G(tau), alle Gewichte positiv-semidefinit", ok_spec)
nopole = [k for k in spec if not spec[k]["poles"]]
say(f"   Bei {len(nopole)} von {len(spec)} Fällen (Impulse nahe den Verdopplern) fehlt der Teilchenpol ganz: Dort gibt")
say("   es nur das Kontinuum. Das Overlap hat keine Verdopplerteilchen, sondern Kontinuumszustände auf Gitterskala.")
befund("Overlap: Teilchenpol fehlt nahe den Verdopplerimpulsen (nur Kontinuum)", len(nopole) > 0)
ecut = [spec[k]["E_cut"] for k in spec if k.startswith("rho=1.0,m=0.0")]
say(f"   Ergebnis: Das freie Overlap (theta = 0) ist Link-reflexionspositiv: Pol mit positivem Residuum,")
say(f"   dazu ein Kontinuum mit positivem Gewicht ab E ~ {min(ecut):.2f} (je nach p, gegen unendlich für p -> 0).")
say("   Die Behauptung des Berichts (Lu-Shao, Overlap-Weylfermionen) ist für das vektorartige freie Overlap")
say("   und damit für die Majoranapaarung bei theta = 0 hier nicht reproduziert. Sie mag die chirale")
say("   Projektion mit P^ = (1 - g5^)/2 betreffen, die unter Theta nicht invariant ist; das ist nicht geprüft.")
befund("Lu-Shao-Behauptung für freies Overlap bei theta = 0 reproduziert", not (ov_ok and ok_spec))
befund("Overlap-Transfermatrix hat kontinuierliches Spektrum (unendlich dimensionaler lokaler Raum)",
       all(abs(spec[k]["zstar"]) > 1e-6 for k in spec))

say("\n2c Gedrehte Massen: der windende Defekt durchläuft alle Phasen theta (Dirac-Form; Majorana = psi-Seite)")
rot = {}
for kind, ths in (("wilson", (0.0, np.pi / 2, np.pi)), ("ov_X", (0.0, 0.3, np.pi / 2, np.pi, 2 * np.pi)),
                  ("ov_Xq", (np.pi / 4, np.pi / 2, np.pi)), ("ov_GJ", (np.pi / 2,))):
    for th in ths:
        res_ = [os_test(kind, p, nb=3, m=0.3, th=th) for p in PS[:4]]
        rot[f"{kind} theta={th:.2f}"] = {"antiherm": max(a for a, _ in res_), "link_min": min(b for _, b in res_)}
for k, v in rot.items():
    say(f"   {k:22s} Link-OS min {v['link_min']:+.2e}, antihermitescher Anteil {v['antiherm']:.0e}")
small_p = [os_test("ov_X", np.array([0.3, 0.1, 0.0]) * sc, nb=3, m=0.3, th=np.pi / 2)[1] for sc in (1.0, 0.25, 0.0625)]
say("   X-Form, theta = pi/2, p = (0,3; 0,1; 0) x 1 / 1/4 / 1/16: OS min " + "/".join(f"{v:+.3f}" for v in small_p)
    + ": Die Verletzung bleibt für p -> 0 und sitzt auf der ersten Zeitschicht (Gitterskala).")
wil_rot_ok = all(rot[k]["link_min"] > TOL for k in rot if k.startswith("wilson"))
check("Wilson mit gedrehter Masse e^{i theta g5}: RP für alle theta", wil_rot_ok)
flat_viol = (rot["ov_X theta=1.57"]["link_min"] < -1e-3 and rot["ov_X theta=3.14"]["link_min"] < -1e-3
             and rot["ov_Xq theta=1.57"]["link_min"] < -1e-4 and rot["ov_X theta=6.28"]["link_min"] < -1e-3)
befund("Flache Overlapformen (X, Xq) verletzen RP bei theta != 0 (Xq: != 0, pi)", flat_viol)
befund("GJ-Form bei theta = pi/2 nicht reflexionssymmetrisch (OS-Form nicht hermitesch)", rot["ov_GJ theta=1.57"]["antiherm"] > 1e-3)
befund("RP-Verletzung verschwindet für p -> 0", min(small_p) > -1e-3)
say("   Ergebnis 2: Die gedrehten Overlapwirkungen bleiben reflexionssymmetrisch (OS-Form hermitesch), sind")
say("   aber nicht positiv. Die exakte Drehung e^{i theta X} ist nicht ultralokal in der Zeit; sie mischt")
say("   Zeitschichten über den Spiegel hinweg. theta = 0 ist positiv, gedrehte Vakua nicht. Wilson ist für jede Phase")
say("   positiv, hat aber die Barriere. Für H6 heißt das: entweder Wand (Wilson) oder Positivität verletzt")
say("   auf Gitterskala (flache Overlapformen). Die euklidisch-lorentzsche Rekonstruktion (KMS-Zustand, Z4/Z6)")
say("   steht mit den flachen Formen bei endlichem a nicht zur Verfügung; ein Hamiltongitter (Z7) hat sie")
say("   automatisch, verliert aber die Gleichwertigkeit von Zeit- und Raumachsen aus Z6.")
ZAHLEN["rp"] = rp
ZAHLEN["rp_spektral"] = spec
ZAHLEN["rp_gedreht"] = rot
ZAHLEN["rp_X_kleines_p"] = small_p

# =========================================================================================== Teil 3
say("\nTeil 3: Gravitation unter D, kappa für chirale Regulatoren (Kettenmethode wie Z15)")


def ent_f(nu):
    nu = nu[(nu > 1e-14) & (nu < 1 - 1e-14)]
    return float(-np.sum(nu * np.log(nu) + (1 - nu) * np.log(1 - nu)))


def S_chain(dx_q, s, B, ring):
    q = 2 * np.pi * (np.arange(ring) + 0.5) / ring
    dx, dy, dz = dx_q(q), np.full(ring, s), B(q)
    E = np.sqrt(dx ** 2 + dy ** 2 + dz ** 2)
    blocks = {(0, 0): 0.5 * (1 - dz / E), (1, 1): 0.5 * (1 + dz / E),
              (0, 1): -0.5 * (dx - 1j * dy) / E, (1, 0): -0.5 * (dx + 1j * dy) / E}
    half = ring // 2
    ph = np.exp(1j * np.outer(np.arange(-half + 1, half), q)) / ring
    idx = np.arange(half)[:, None] - np.arange(half)[None, :] + half - 1
    Cm = np.zeros((half, 2, half, 2), complex)
    for (a, b), v in blocks.items():
        Cm[:, a, :, b] = (ph @ v)[idx]
    return ent_f(np.linalg.eigvalsh(Cm.reshape(2 * half, 2 * half))) / 2


saw = lambda q: np.where(q > np.pi, q - 2 * np.pi, q)


def kappa(kind, m, ring, nk, r=1.0):
    k = 2 * np.pi * (np.arange(nk // 2) + 0.5) / nk
    tot = wsum = 0.0
    for a in range(nk // 2):
        for b in range(a, nk // 2):
            w = 1.0 if a == b else 2.0
            if kind == "wilson":
                s = np.hypot(np.sin(k[a]), np.sin(k[b]))
                Bq = lambda q: m + r * (2 - np.cos(k[a]) - np.cos(k[b])) + r * (1 - np.cos(q))
                tot += w * S_chain(np.sin, s, Bq, ring)
            else:
                s = np.hypot(k[a], k[b])                       # SLAC: Impuls selbst, |k| < pi
                tot += w * S_chain(saw, s, lambda q: m + 0 * q, ring)
            wsum += w
    return 2 * tot / wsum


z16 = json.loads((HERE / "z16_ergebnis.json").read_text())
t3 = time.time()
kw1 = kappa("wilson", 0.0, 160, 48)
say(f"   Wilson r = 1, masselos (eigener Code, 48^2 Querimpulse, Ring 160): kappa = {kw1:.4f} je Dirac (Z15/Z16: 0,2209).")
check("Wilson kappa(r = 1) reproduziert", abs(kw1 - 0.2209) < 0.002)
kw_m = {ring: kappa("wilson", 0.5, ring, 16) for ring in (100, 200, 400)}
ks_m = {ring: kappa("slac", 0.5, ring, 16) for ring in (100, 200, 400)}
say("   Masse 0,5, Ring 100/200/400 (16^2): Wilson " + "/".join(f"{v:.4f}" for v in kw_m.values())
    + "; SLAC " + "/".join(f"{v:.4f}" for v in ks_m.values()))
dS1, dS2 = ks_m[200] - ks_m[100], ks_m[400] - ks_m[200]
say(f"   SLAC wächst je Verdopplung um {dS1:.4f} und {dS2:.4f}: logarithmisch, also ohne Flächengesetz. Grund: der")
say("   Sprung der Dispersion am Zonenrand (q = +-pi) wirkt wie eine Fermifläche. Wilson ist konvergiert.")
check("Wilson konvergiert, SLAC wächst logarithmisch", abs(kw_m[400] - kw_m[200]) < 1e-4 and dS1 > 0.01 and abs(dS2 / dS1 - 1) < 0.2)
befund("SLAC (chiral, ohne Verdoppler, lokal in k) hat kein endliches kappa", dS2 > 0.01)
say(f"   (Rechenzeit Teil 3: {time.time() - t3:.0f} s)")
kap_r = {float(r): v for r, v in z16["kappa_wilson"].items()}
say("   Wilson-Band aus Z16 (je Dirac): " + ", ".join(f"r = {r:g}: {v:.3f}" for r, v in kap_r.items())
    + "; nur 24 leichte Dirac: a/l_P = " + ", ".join(f"{2 * np.sqrt(24 * v):.2f}" for v in kap_r.values()))
say("   Overlap: Die OS-Rekonstruktion bei theta = 0 existiert (Teil 2), aber jede Impulsschicht trägt ein")
say("   Kontinuum von Zuständen (Schnitt). Der Gitterhilbertraum ist kein Produkt endlicher Platzräume;")
say("   kappa als Entropie je Gitterfläche ist ohne weitere Abschneidewahl nicht definiert. Bei gedrehter")
say("   Phase fehlt sogar die Positivität (Teil 2c).")
beta = G4
alph = [1j * G4 @ gam[i] for i in range(3)]
rng = np.random.default_rng(19)
qerr, qspec = 0.0, []
for _ in range(20):
    kv = rng.uniform(-np.pi, np.pi, 3)
    sv = np.sin(kv)
    Hk = sum(sv[i] * alph[i] for i in range(3)) + beta * np.sum(1 - np.cos(kv))
    E = np.sqrt(np.sum(sv ** 2) + np.sum(1 - np.cos(kv)) ** 2)
    Pp, Pm = (I4 + Hk / E) / 2, (I4 - Hk / E) / 2
    Q = Pp @ g5 @ Pp + Pm @ g5 @ Pm
    qerr = max(qerr, np.max(np.abs(Q @ Hk - Hk @ Q)), np.max(np.abs(np.sort(np.abs(np.linalg.eigvalsh(Q))) - np.linalg.norm(sv) / E)))
    qspec.append(np.linalg.norm(sv) / E)
say(f"   Hamiltonsch (Wilson-Dirac, r = 1): Q5 = P_+ g5 P_+ + P_- g5 P_- vertauscht mit H, Eigenwerte +-|sin k|/E")
say(f"   (Rest {qerr:.1e}), hier zwischen {min(qspec):.2f} und {max(qspec):.2f}: exakt, aber nicht quantisiert, wie X.")
check("Hamiltonsche Chiralladung Q5: [Q5, H] = 0, Spektrum |s|/E", qerr < 1e-12)
say("   Hamiltonsch mit exakter, aber nicht quantisierter Chiralladung (bandprojiziertes g5, wie X) bleibt der")
say("   Wilson-Hamiltonoperator selbst: kappa = kappa_Wilson(r), die Chiralität legt r nicht fest.")
befund("Chiraler Regulator legt a/l_P fest", False)
ZAHLEN["kappa"] = {"wilson_r1": kw1, "wilson_m0.5": {str(k): v for k, v in kw_m.items()},
                   "slac_m0.5": {str(k): v for k, v in ks_m.items()}, "wilson_band_z16": z16["kappa_wilson"],
                   "a_lP_24_dirac": {str(r): float(2 * np.sqrt(24 * v)) for r, v in kap_r.items()}}

# =========================================================================================== Urteil
URTEIL = {
    "1a_phasenpotential": f"hergeleitet: flach (Spannweite < 1e-12) nur mit der Majoranadrehung X (lokal, nicht "
                          f"quantisiert, J^4 != 1) oder sign X (quantisiert, nichtlokal); die lokale eindeutige GJ-Form "
                          f"behält {np.ptp(pot['GJ']):.4f} je Zelle mit Periode pi; Wilson {np.ptp(pot['W']):.3f} euklidisch, "
                          f"{vh[0] - vh[1]:.4f} hamiltonsch",
    "1b_index": "hergeleitet: eine chirale Majorana-Weyl-Mode je Paar (Index = Windung von eta_+) für X, Xq, GJ und "
                "Wilson; Verdoppler tragen beim Overlap nicht bei; AB-Viertelfluss nicht im Overlap enthalten",
    "1c_wand": f"gescheitert für lokale Formen, Lesart für Xq: wandfrei nur mit nichtlokaler Drehung (dF ~ 1/L); "
               f"lokal-flach (X) wird der Schnitt zur Wand (Spannung {sig_cut[10]:.4f} gegen Wilson <= {sW:.3f} je "
               f"Einheit Wandquerschnitt), lokal-eindeutig (GJ) bleiben zwei pi-Wände (je <= {sG:.4f})",
    "2_positivitaet": "hergeleitet (Zielkonflikt): Wilson r <= 1 Link-RP für jede Massenphase; freies Overlap bei "
                      "theta = 0 Link-RP (Pol + positives Kontinuum), Lu-Shao dafür nicht reproduziert; flache gedrehte "
                      "Formen (X, Xq) verletzen RP auf Gitterskala, auch für p -> 0; GJ bei pi/2 nicht reflexionssymmetrisch. "
                      "Für H6: Wand oder Positivität",
    "3_gravitation": f"gescheitert: kein chiraler Regulator macht kappa bestimmter; Overlap ohne endlichen lokalen "
                     f"Hilbertraum, SLAC ohne Flächengesetz (+{dS2:.3f} je Verdopplung), Hamilton-chiral = Wilson(r); "
                     f"a/l_P bleibt regulatorabhängig (Wilson r = 1: kappa = {kw1:.4f} je Dirac)",
}
say("\nUrteil")
for k, v in URTEIL.items():
    say(f"   {k}: {v}")
say(f"\nPrüfungen: {sum(CHECKS.values())}/{len(CHECKS)} bestanden. Befunde: " +
    ", ".join(f"{k} = {'ja' if v else 'nein'}" for k, v in BEFUNDE.items()))
say(f"Laufzeit {time.time() - T0:.0f} s.")
(HERE / "ausgabe_z19.txt").write_text("\n".join(LINES) + "\n")
(HERE / "z19_ergebnis.json").write_text(json.dumps({"checks": CHECKS, "befunde": BEFUNDE, "urteil": URTEIL,
                                                    "zahlen": ZAHLEN}, indent=1, ensure_ascii=False, default=float))
assert all(CHECKS.values())
