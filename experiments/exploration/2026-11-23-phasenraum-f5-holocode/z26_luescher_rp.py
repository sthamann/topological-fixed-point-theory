"""Z26: Reflexionspositivität für Lüschers projizierte Weylfelder (Overlap, Ginsparg-Wilson), frei, abelsch.

Lücke aus Z19/Z21 (S3, Querbefunde 6.b/6.c): Entscheidung D4 in H6 lautet "Wand mit Positivität (G mit zwei
abgestimmten Gegentermen) oder Positivitätsverletzung ohne Wand (flacher GW-Zweig)". Z19 hat offen gelassen, ob
Lüschers projizierte Weylfelder reflexionspositiv (RP) sind. Hier, bei endlichem Gitterabstand (a = 1):

    D = 1 + g5 sign(g5 (D_W - m0)),  g5 D + D g5 = D g5 D,  g5^ = g5 (1 - D),
    Weylfeld psi = P^_- psi, psi-quer = psi-quer P_+,  Propagator S_W = P^_- D^-1 P_+ = D^-1 P_+.

 0. Exakte Identitäten (GW): D^-1 = 1/2 + (chiral ungerade), also S_W = P_+/2 delta + P_- D^-1 P_+;
    Verflechter (1 - D/2) g5^ = g5 (1 - D/2); mit Lüscher-Paarung (GJ-Form, Phase theta)
    S_theta = (1 - (m/2) e^{i theta g5}) G_chi(theta) + 1/2,  G_chi = <chi psi-quer>,  chi = (1 - D/2) psi.
 1. Bildet die OS-Reflexion Theta den projizierten Raum auf sich ab? (Link- und Platzreflexion, R D R = D^+)
 2. RP-Matrizen M_ij = <Theta(F_i) F_j> für die Weylfelder: lineare Sektoren (psi-quer, psi) im Impulsraum
    (m0, L_t, L_s, Holonomie, Abstand zur Ebene), im Ortsraum (ganze positive Halbzeit, auch mit nicht
    konstantem, Theta-symmetrischem U(1)-Hintergrund), Bilineare samt Eins über das Wick-Theorem, dazu die
    Spektraldarstellung (Pole + Schnitt) für L_t = unendlich. Kalibriert an Wilson (r = 1 positiv, r = 2 nicht,
    Gegenvorzeichen nicht) und am vollen Overlap-Dirac.
 3. Gedrehte Phase: zwei Lüscher-Weylfelder (L, R) mit Paarung m e^{i theta}; theta ist genau die GW-chirale
    Drehung. Zwei Reflexionen: Theta_D auf (psi, psi-quer) wie Z19, Theta_W auf (chi, psi-quer).
 4. 4D: Impulsraum-Scans, Spektraldarstellung, kleine Ortsraumgitter.
 5. Einordnung ins Z24-Trilemma: Flachheit, Index der Naht (Z19-Querschnitt), Nichtlokalität der Variablen,
    in denen Theta_W ortslokal wirkt (D_c = 2D/(2 - D)), Majoranadrehung X aus Z19.
 6. H6s G selbst: psi, chi Weylfelder gleicher Chiralität, Paarung Phi psi chi, J: psi -> chi -> -psi (J^2 = -1).
    Diraclesart Psi = (psi, chi^c); Gitter-J: Psi -> g5 C Psi-quer^T, Psi-quer -> Psi^T C^-1 g5. Erhält J Lüschers
    Projektoren und die Overlap-Wirkung? Gibt es eine exakte (ggf. nichtlokale) J, und verträgt sie sich mit Theta_W?

Zitiert, nicht geprüft: Lüscher 1977/1998, Osterwalder-Seiler 1978, Neuberger 1998, Furman-Shamir 1995.
"""
import json
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
T0 = time.time()
TOL = -1e-9


def say(s=""):
    print(s, flush=True)


def check(name, ok):
    CHECKS[name] = bool(ok)
    say(f"   [{'ok' if ok else 'FEHLER'}] {name}")


def befund(name, val):
    BEFUNDE[name] = bool(val)
    say(f"   [Befund: {'ja' if val else 'nein'}] {name}")


# =========================================================================================== Grundlagen
s1 = np.array([[0, 1], [1, 0]], complex)
s2 = np.array([[0, -1j], [1j, 0]])
s3 = np.diag([1.0, -1.0]).astype(complex)
gam4 = [np.array(x, complex) for x in (                          # Z19-Konvention, Zeit = gam4[3]
    [[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]],
    [[0, 0, 0, -1j], [0, 0, -1j, 0], [0, 1j, 0, 0], [1j, 0, 0, 0]],
    [[0, 0, 0, -1], [0, 0, 1, 0], [0, 1, 0, 0], [-1, 0, 0, 0]],
    [[0, 0, -1j, 0], [0, 0, 0, 1j], [1j, 0, 0, 0], [0, -1j, 0, 0]])]
DIM = {2: dict(gam=[s2, s1], g5=s3),                             # 2D: g_x = s2, g_t = s1, g5 = s3 = -i g_t g_x
       4: dict(gam=gam4, g5=gam4[0] @ gam4[1] @ gam4[2] @ gam4[3])}
for _d in DIM:
    _g5 = DIM[_d]["g5"]
    _w, _U = np.linalg.eigh(_g5)
    DIM[_d]["I"] = np.eye(len(_g5), dtype=complex)
    DIM[_d]["gt"] = DIM[_d]["gam"][-1]
    DIM[_d]["E+"] = _U[:, _w > 0]
    DIM[_d]["E-"] = _U[:, _w < 0]
    DIM[_d]["Ef"] = np.eye(len(_g5), dtype=complex)


def rot(dim, th):
    g5 = DIM[dim]["g5"]
    return np.cos(th) * DIM[dim]["I"] + 1j * np.sin(th) * g5


# ------------------------------------------------------------------------------- Impulsraum (fester Raumimpuls)
def ops_k(dim, z, ps, m0=1.0, r=1.0):
    """D_W und Overlap D bei z = e^{i p_t} (auch komplex, analytisch fortgesetzt) und Raumimpuls ps."""
    gam, I = DIM[dim]["gam"], DIM[dim]["I"]
    z = np.atleast_1d(np.asarray(z, complex))
    ps = np.atleast_1d(np.asarray(ps, float))
    st, ct = (z - 1 / z) / 2j, (z + 1 / z) / 2
    ss, Wsp = np.sin(ps), np.sum(1 - np.cos(ps))
    kin = 1j * st[:, None, None] * gam[-1] + sum(1j * ss[i] * gam[i] for i in range(dim - 1))[None]
    W = r * (Wsp + 1 - ct)
    DW = kin + W[:, None, None] * I
    f = st ** 2 + np.sum(ss ** 2) + (W - m0) ** 2
    D = I + (kin + (W - m0)[:, None, None] * I) / np.sqrt(f)[:, None, None]
    return DW, D


def inv_op(kind, dim, z, ps, m0=1.0, m=0.0, th=0.0, r=1.0):
    I = DIM[dim]["I"]
    DW, D = ops_k(dim, z, ps, m0, r)
    if kind == "wilson":
        return DW + m * rot(dim, th)
    if kind == "weyl":
        return D
    return D + m * rot(dim, th) @ (I - D / 2)


def prop(kind, dim, z, ps, m0=1.0, m=0.0, th=0.0, r=1.0):
    """wilson: (D_W + m e^{i th g5})^-1; ov: Overlap mit Lüscher-Paarung (GJ); chi: (1 - D/2) S_ov; weyl: D^-1 P_+."""
    I, g5 = DIM[dim]["I"], DIM[dim]["g5"]
    DW, D = ops_k(dim, z, ps, m0, r)
    if kind == "wilson":
        return np.linalg.inv(DW + m * rot(dim, th))
    if kind == "weyl":
        return np.linalg.inv(D) @ ((I + g5) / 2)
    S = np.linalg.inv(D + m * rot(dim, th) @ (I - D / 2))
    return S if kind == "ov" else (I - D / 2) @ S


def G_tau(kind, dim, ps, taus, Lt, **kw):
    """Zeitpropagator G(tau) bei antiperiodischer Zeit L_t (Vorzeichen symmetrisch zur Ebene t = 1/2)."""
    pt = 2 * np.pi * (np.arange(Lt) + 0.5) / Lt
    S = prop(kind, dim, np.exp(1j * pt), ps, **kw)
    ph = np.exp(1j * np.outer(np.asarray(sorted(taus)), pt))
    vals = np.einsum("tn,nab->tab", ph, S) / Lt
    out = {int(t): v for t, v in zip(sorted(taus), vals)}
    out["rausch"] = 1e-15 * float(np.max(np.abs(S)))
    return out


def theta_t(link):
    return (lambda s: 1 - s) if link else (lambda s: -s)


def os_form(G, dim, T, link, side, sigma, E):
    """psi-quer-Seite: sigma g_t G(theta s - t); psi-Seite: -sigma (G(t - theta s) g_t)^T; komprimiert mit E."""
    gt, th = DIM[dim]["gt"], theta_t(link)
    k, n = E.shape[1], len(T)
    A = np.zeros((n * k, n * k), complex)
    for i, s in enumerate(T):
        for j, t in enumerate(T):
            if side == "bar":
                blk = E.conj().T @ (sigma * gt @ G[th(s) - t]) @ E
            else:
                blk = (E.conj().T @ (-sigma * G[t - th(s)] @ gt) @ E).T
            A[i * k:(i + 1) * k, j * k:(j + 1) * k] = blk
    return A


def assess(A, k):
    """(antihermitescher Anteil relativ, kleinster/größter Eigenwert der je Block normierten hermiteschen Form)."""
    aher = float(np.max(np.abs(A - A.conj().T)) / max(np.max(np.abs(A)), 1e-300))
    H = (A + A.conj().T) / 2
    n = H.shape[0] // k
    nrm = np.array([np.linalg.norm(H[k * i:k * i + k, k * i:k * i + k]) for i in range(n)])
    d = np.repeat(np.sqrt(np.maximum(nrm, 1e-14 * nrm.max())), k)
    ev = np.linalg.eigvalsh(H / np.outer(d, d))
    return aher, float(ev.min() / np.abs(ev).max()), float(nrm.min()), float(np.linalg.eigvalsh(H).min())


def taus_for(T, link):
    th = theta_t(link)
    return sorted({th(s) - t for s in T for t in T} | {t - th(s) for s in T for t in T})


RAUSCH = {"rauschfaelle": 0, "max_roh_durch_rausch": 0.0, "min_normiert": 0.0, "roh_dazu": 0.0}


def os_scan(kind, dim, ps, Lt, link, forms, Ts, sigma=-1, **kw):
    """forms: Liste (Seite, E); Ts: Liste von Zeitschichtmengen. Gibt Liste (aher, min) zurück. Blöcke unter
    1e3 x Rauschgrenze der Fouriersumme (1e-15 max|G(p_t)|) werden ausgelassen. Ein negativer normierter Eigenwert
    zählt nur, wenn auch der rohe Eigenwert unter -1e3 x Rauschgrenze liegt; sonst ist er Rundung (gezählt, auf 0)."""
    allt = sorted(set().union(*[taus_for(T, link) for T in Ts]))
    G = G_tau(kind, dim, ps, allt, Lt, **kw)
    out = []
    for T in Ts:
        for side, E in forms:
            a, b, nmin, raw = assess(os_form(G, dim, T, link, side, sigma, E), E.shape[1])
            if nmin > 1e3 * G["rausch"]:
                if b < 0 and raw > -1e3 * G["rausch"]:
                    RAUSCH["rauschfaelle"] += 1
                    RAUSCH["max_roh_durch_rausch"] = max(RAUSCH["max_roh_durch_rausch"], -raw / G["rausch"])
                    if b < RAUSCH["min_normiert"]:
                        RAUSCH["min_normiert"], RAUSCH["roh_dazu"] = b, raw
                    b = 0.0
                out.append((a, b))
    return out


def spectral(kind, dim, ps, side, E, sigma=-1, ngl=240, nrec=6, **kw):
    """Spektraldarstellung für L_t = unendlich: G(-n) (psi-quer-Seite) bzw. G(+n) (psi-Seite)
    = Summe_Pole w0^(n-1) Res + sign(z*) Int_0^z* l^(n-1) rho(l) dl. RP <=> alle Gewichte (komprimiert) >= 0."""
    gt = DIM[dim]["gt"]
    ps = np.atleast_1d(np.asarray(ps, float))
    m0, r = kw.get("m0", 1.0), kw.get("r", 1.0)
    zz = (lambda w: 1 / np.asarray(w, complex)) if side == "bar" else (lambda w: np.asarray(w, complex))
    F = lambda w: prop(kind, dim, zz(w), ps, **kw)
    Fi = lambda w: inv_op(kind, dim, zz(w), ps, **kw)
    ss, Wsp = np.sin(ps), np.sum(1 - np.cos(ps))
    Bp = r * (Wsp + 1) - m0
    if kind == "wilson" or abs(Bp) < 1e-12:
        zst = 0.0
    else:
        cst = (1 + np.sum(ss ** 2) + Bp ** 2) / (2 * Bp)
        zst = float(cst - np.sign(cst) * np.sqrt(cst ** 2 - 1))
    u, wu = np.polynomial.legendre.leggauss(ngl)
    u, wu = (u + 1) / 2, wu / 2
    lam = zst * np.sin(np.pi * u / 2) ** 2
    jac = zst * np.pi * np.sin(np.pi * u / 2) * np.cos(np.pi * u / 2)
    eps = 1e-12 * max(abs(zst), 1e-3)
    nsp = len(DIM[dim]["g5"])
    rho_c = (F(lam - 1j * eps) - F(lam + 1j * eps)) / (2j * np.pi) if zst != 0 else np.zeros((ngl, nsp, nsp))

    def wmat(Wt):
        if side == "bar":
            return E.conj().T @ (sigma * gt @ Wt) @ E
        return (E.conj().T @ (-sigma * Wt @ gt) @ E).T
    ngrid = 4000
    grid = np.concatenate([np.linspace(max(zst, 0) + 1e-6, 0.999, ngrid), np.linspace(-0.999, min(zst, 0) - 1e-6, ngrid)])
    smin = np.linalg.svd(Fi(grid), compute_uv=False)[:, -1]
    poles = []
    for i in list(range(1, ngrid - 1)) + list(range(ngrid + 1, 2 * ngrid - 1)):
        if smin[i] < smin[i - 1] and smin[i] < smin[i + 1] and smin[i] < 1e-2:
            a, b = grid[i - 1], grid[i + 1]
            for _ in range(80):
                c1, c2 = a + (b - a) / 3, b - (b - a) / 3
                f1 = np.linalg.svd(Fi(c1), compute_uv=False)[0, -1]
                f2 = np.linalg.svd(Fi(c2), compute_uv=False)[0, -1]
                a, b = (a, c2) if f1 < f2 else (c1, b)
            poles.append(float((a + b) / 2))
    res = []
    for w0 in poles:
        rr = 0.3 * min(abs(w0 - zst), 1 - abs(w0), abs(w0), 0.05)
        thc = 2 * np.pi * (np.arange(64) + 0.5) / 64
        res.append(np.einsum("n,nab->ab", rr * np.exp(1j * thc), F(w0 + rr * np.exp(1j * thc))) / 64)
    Wc = [wmat(x) for x in rho_c] if zst != 0 else []
    Wp = [wmat(R) for R in res]
    scale = max([np.abs(np.linalg.eigvalsh((A + A.conj().T) / 2)).max() for A in Wp] +
                [abs(zst) * np.abs(np.linalg.eigvalsh((A + A.conj().T) / 2)).max() for A in Wc] + [1e-300])
    wmin, aher = 0.0, 0.0
    for A, fac in [(A, 1.0) for A in Wp] + [(A, abs(zst)) for A in Wc]:
        wmin = min(wmin, fac * np.linalg.eigvalsh((A + A.conj().T) / 2).min() / scale)
        aher = max(aher, fac * np.max(np.abs(A - A.conj().T)) / scale)
    taus = np.arange(1, nrec + 1)
    Gf = G_tau(kind, dim, ps, list(-taus) + list(taus), 4096, **kw)
    sc = max(np.abs(Gf[int(t)]).max() for t in list(-taus) + list(taus))
    err = 0.0
    for t in taus:
        g = sum(w0 ** (t - 1) * R for w0, R in zip(poles, res))
        if zst != 0:
            g = g + np.sign(zst) * np.einsum("n,nab->ab", wu * jac * lam ** (t - 1), rho_c)
        err = max(err, float(np.abs(g - Gf[int(-t if side == "bar" else t)]).max() / sc))
    return dict(zst=zst, poles=poles, err=err, wmin=float(wmin), aher=float(aher))


# ------------------------------------------------------------------------------------- Ortsraum
def build_DW(dim, Lt, Ls, r=1.0, phi=0.0, bg=0.0):
    """Wilson-Dirac auf L_t x L_s^(d-1), antiperiodisch in der Zeit (Vorzeichen auf dem Link L_t/2 -> L_t/2+1,
    symmetrisch zu t -> 1 - t), Raumholonomie phi, in 2D optional ein glatter, flussfreier Hintergrund, der unter
    der Linkreflexion symmetrisch ist (Raumlinks gerade in t - 1/2, Zeitlinks ungerade)."""
    gam, I = DIM[dim]["gam"], DIM[dim]["I"]
    ns = len(I)
    shape = (Lt,) + (Ls,) * (dim - 1)
    V = int(np.prod(shape))
    co = np.array(list(np.ndindex(*shape)))
    D = np.zeros((ns * V, ns * V), complex)
    Dv = D.reshape(V, ns, V, ns)
    ii = np.arange(V)
    Dv[ii, :, ii, :] += r * dim * I
    t, x = co[:, 0], co[:, 1]
    for mu in range(dim):
        e = np.zeros(dim, int)
        e[mu] = 1
        c2 = (co + e) % np.array(shape)
        jj = np.ravel_multi_index(tuple(c2.T), shape)
        if mu == 0:
            g = gam[-1]
            a = 0.7 * bg * np.sin(2 * np.pi * x / Ls) * np.sin(2 * np.pi * t / Lt)
            U = np.exp(1j * a) * np.where(t == Lt // 2, -1.0, 1.0)
        else:
            g = gam[mu - 1]
            a = phi / Ls + (bg * np.cos(2 * np.pi * x / Ls) * np.cos(2 * np.pi * (t - 0.5) / Lt) if mu == 1 else 0.0)
            U = np.exp(1j * a) * np.ones(V)
        Dv[ii, :, jj, :] += -(r * I - g)[None] / 2 * U[:, None, None]
        Dv[jj, :, ii, :] += -(r * I + g)[None] / 2 * np.conj(U)[:, None, None]
    return D, co


def overlap(DW, dim, m0=1.0):
    n = DW.shape[0]
    G5 = np.kron(np.eye(n // len(DIM[dim]["g5"])), DIM[dim]["g5"])
    H = G5 @ (DW - m0 * np.eye(n))
    H = (H + H.conj().T) / 2
    w, U = np.linalg.eigh(H)
    return np.eye(n) + G5 @ ((U * np.sign(w)) @ U.conj().T), G5, float(np.abs(w).min())


def refl(dim, Lt, Ls, link=True):
    shape = (Lt,) + (Ls,) * (dim - 1)
    V = int(np.prod(shape))
    ns = len(DIM[dim]["g5"])
    co = np.array(list(np.ndindex(*shape)))
    c2 = co.copy()
    c2[:, 0] = ((1 - co[:, 0]) if link else (-co[:, 0])) % Lt
    jj = np.ravel_multi_index(tuple(c2.T), shape)
    R = np.zeros((ns * V, ns * V), complex)
    R.reshape(V, ns, V, ns)[np.arange(V), :, jj, :] = DIM[dim]["gt"]
    return R


def half_basis(dim, co, Lt, P, slices=None):
    """Isometrie auf {positive Halbzeit} x Spinorraum P (Spalten); slices: Auswahl von Zeitschichten."""
    sl = set(range(1, Lt // 2 + 1)) if slices is None else set(slices)
    sites = [i for i, c in enumerate(co) if c[0] in sl]
    ns = len(DIM[dim]["g5"])
    k = P.shape[1]
    E = np.zeros((ns * len(co), len(sites) * k), complex)
    for a, i in enumerate(sites):
        E[ns * i:ns * i + ns, k * a:k * a + k] = P
    return E, sites


def os_real(G, R, E, side, sigma=-1):
    if side == "bar":
        return sigma * E.conj().T @ (R @ G) @ E
    return (E.conj().T @ (-sigma * G @ R) @ E).T


say("Z26: Reflexionspositivität für Lüschers projizierte Weylfelder (Overlap, GW), frei, abelsch, a = 1")
say("Konventionen: D = 1 + g5 sign(g5 (D_W - m0)), r = 1; Weylfeld psi = P^_- psi, psi-quer = psi-quer P_+;")
say("OS-Reflexion Theta psi-quer(t) = -g_t psi(theta t), Theta psi(t) = -psi-quer(theta t) g_t (antilinear,")
say("Ordnung umkehrend), Link theta t = 1 - t, Platz theta t = -t. Das Vorzeichen ist an Wilson r = 1 geeicht.\n")

# =========================================================================================== Teil 0
say("Teil 0: exakte Identitäten (2D-Ortsraum 16 x 16, Holonomie 0,7, mit und ohne Theta-symmetrischen Hintergrund)")
id_rest = {}
REAL2 = {}
for bg in (0.0, 0.8):
    DW, co = build_DW(2, 16, 16, phi=0.7, bg=bg)
    D, G5, hmin = overlap(DW, 2)
    n = D.shape[0]
    Ib = np.eye(n)
    R = refl(2, 16, 16)
    PP, PM = (Ib + G5) / 2, (Ib - G5) / 2
    Di = np.linalg.inv(D)
    g5h = G5 @ (Ib - D)
    m_, th_ = 0.3, 1.1
    Eth = np.kron(np.eye(n // 2), rot(2, th_))
    S_th = np.linalg.inv(D + m_ * Eth @ (Ib - D / 2))
    Gchi = (Ib - D / 2) @ S_th
    id_rest[bg] = {
        "GW": float(np.abs(G5 @ D + D @ G5 - D @ G5 @ D).max()),
        "RDWR=DW+": float(np.abs(R @ DW @ R - DW.conj().T).max()),
        "RDR=D+": float(np.abs(R @ D @ R - D.conj().T).max()),
        "D^-1-1/2 ungerade": float(np.abs(G5 @ (Di - Ib / 2) + (Di - Ib / 2) @ G5).max()),
        "S_W=P+/2+P-D^-1P+": float(np.abs(Di @ PP - (PP / 2 + PM @ Di @ PP)).max()),
        "Verflechter": float(np.abs((Ib - D / 2) @ g5h - G5 @ (Ib - D / 2)).max()),
        "R g5^ R = -(1-D) g5": float(np.abs(R @ g5h @ R + (Ib - D) @ G5).max()),
        "S_theta-Zerlegung": float(np.abs(S_th - ((Ib - m_ / 2 * np.kron(np.eye(n // 2), rot(2, th_))) @ Gchi + Ib / 2)).max()),
    }
    REAL2[bg] = dict(DW=DW, D=D, G5=G5, R=R, co=co, hmin=hmin)
    say(f"   Hintergrund {bg}: kleinster |H_W| = {hmin:.3f}; Reste " + ", ".join(f"{k} {v:.1e}" for k, v in id_rest[bg].items()))
check("GW, R D_W R = D_W^+ und R D R = D^+ (auch mit Theta-symmetrischem Hintergrund)",
      all(v["GW"] < 1e-12 and v["RDWR=DW+"] < 1e-12 and v["RDR=D+"] < 1e-12 for v in id_rest.values()))
check("D^-1 = 1/2 + chiral ungerade, also S_W = P_+/2 delta + P_- D^-1 P_+ (Rest < 1e-12)",
      all(v["D^-1-1/2 ungerade"] < 1e-12 and v["S_W=P+/2+P-D^-1P+"] < 1e-12 for v in id_rest.values()))
check("Verflechter (1 - D/2) g5^ = g5 (1 - D/2) und S_theta = (1 - (m/2) e^{i theta g5}) G_chi + 1/2",
      all(v["Verflechter"] < 1e-12 and v["S_theta-Zerlegung"] < 1e-12 for v in id_rest.values()))
say("   Folge: Der Weylpropagator ist ein reiner Kontaktterm P_+/2 plus der chirale Block des Overlap-Propagators.")
say("   Die P_+-Komponente des Weylfelds, P_+ psi = (1/2) P_+ D psi, ist ein Bewegungsgleichungsfeld: Ihre")
say("   Korrelation mit psi-quer ist (1/2) delta. Die GW-Drehung e^{i a g5^} wirkt auf chi = (1 - D/2) psi ortslokal als e^{i a g5}.")
ZAHLEN["identitaeten_rest_2d"] = {str(k): v for k, v in id_rest.items()}

# =========================================================================================== Teil 1
say("\nTeil 1: Bildet Theta den projizierten Raum auf sich ab?")
r0 = REAL2[0.0]
D, G5, R = r0["D"], r0["G5"], r0["R"]
n = D.shape[0]
Ib = np.eye(n)
Phm = (Ib - G5 @ (Ib - D)) / 2
PP, PM = (Ib + G5) / 2, (Ib - G5) / 2
dev_bar = float(np.linalg.norm(R @ Phm @ R - PP, 2))
dev_psi = float(np.linalg.norm(R @ PP @ R - Phm, 2))
mirror = float(np.abs(R @ Phm @ R - (Ib + (Ib - D) @ G5) / 2).max())
dnorm = float(np.linalg.norm(D, 2)) / 2
say(f"   R P^_- R = (1 + (1 - D) g5)/2 (gespiegelter Hut, Rest {mirror:.1e}), nicht P_+:")
say(f"   ||R P^_- R - P_+|| = {dev_bar:.4f}, ||R P_+ R - P^_-|| = {dev_psi:.4f}, ||D||/2 = {dnorm:.4f} (16 x 16).")
pk = {}
for p in ((0.05, 0.0), (0.2, 0.1), (1.0, 0.5), (np.pi - 0.05, 0.0), (np.pi - 0.05, np.pi - 0.05)):
    _, Dk = ops_k(2, np.exp(1j * p[0]), np.array([p[1]]))
    pk[f"({p[0]:.2f},{p[1]:.2f})"] = float(np.linalg.norm(Dk[0], 2) / 2)
say("   Impulsaufgelöst ||D(p)||/2 bei (p_t, p_x) = " + ", ".join(f"{k}: {v:.3f}" for k, v in pk.items()))
say("   Die Abweichung ist O(a p) bei kleinen Impulsen und O(1) an den Verdopplerecken.")
check("Theta bildet (P^_-, P_+) auf (P_-, gespiegelter Hut) ab; Abweichung = ||D||/2 (Rest < 1e-12)",
      mirror < 1e-12 and abs(dev_bar - dnorm) < 1e-10 and abs(dev_psi - dnorm) < 1e-10)
befund("Theta bildet Lüschers projizierten Raum auf sich ab", dev_bar < 1e-8)
Di = np.linalg.inv(D)
SW = Di @ PP
cont = PP @ SW
diag_c = np.array([cont[2 * i, 2 * i] for i in range(n // 2)])
off = cont - np.diag(np.diag(cont))
say(f"   Kontaktfeld: <P_+psi(x) psi-quer(y)> hat Diagonale {diag_c.real.min():.6f} ... {diag_c.real.max():.6f},")
say(f"   außerhalb der Diagonale max {np.abs(off).max():.1e}. Theta(P_+ psi) = psi-quer P_- g_t ist im Weylsektor null,")
say("   P_+ psi selbst nicht (Selbstkontraktion 1/2): Theta ist auf der vollen Komponentenalgebra nicht involutiv.")
say("   Abgeschlossen ist die Ankeralgebra {psi-quer P_+, chi_L = P_- psi = P_-(1 - D/2) psi}: Theta(psi-quer P_+) = -g_t chi_L.")
check("P_+-Komponente des Weylfelds ist reiner Kontakt (Diagonale 1/2, sonst < 1e-12)",
      np.abs(diag_c - 0.5).max() < 1e-12 and np.abs(off).max() < 1e-12)
befund("Theta ist auf der vollen Komponentenalgebra der Weylfelder involutiv", False)
ZAHLEN["theta_abweichung"] = {"R_Phat_R_minus_P+": dev_bar, "D_norm_halbe": dnorm, "impulsaufgeloest": pk}

# =========================================================================================== Teil 2
say("\nTeil 2: RP-Matrizen der Weylfelder, masselos")
say("2a Kalibrierung (2D, Impulsraum, L_t = 4096, Raumimpulse 0,3/1,0/2,5/3,1, 4 Zeitschichten an der Ebene)")
PS2 = [np.array([p]) for p in (0.3, 1.0, 2.5, 3.1)]
Ef2, Ep2, Em2 = DIM[2]["Ef"], DIM[2]["E+"], DIM[2]["E-"]
cal = {}
for lab, kind, kw, fb, fp in (("Wilson r=1", "wilson", dict(r=1.0, m=0.2), Ef2, Ef2),
                              ("Wilson r=2", "wilson", dict(r=2.0, m=0.2), Ef2, Ef2),
                              ("Overlap m0=0,5", "ov", dict(m0=0.5), Ef2, Ef2),
                              ("Overlap m0=1", "ov", dict(m0=1.0), Ef2, Ef2),
                              ("Overlap m0=1,5", "ov", dict(m0=1.5), Ef2, Ef2),
                              ("Overlap m0=1, m=0,3", "ov", dict(m0=1.0, m=0.3), Ef2, Ef2),
                              ("Weyl m0=1", "weyl", dict(m0=1.0), Ep2, Em2)):
    row = {}
    for link in (True, False):
        for sg in (-1, 1):
            vals = []
            for p in PS2:
                T = [1, 2, 3, 4] if link else [0, 1, 2, 3]
                vals += os_scan(kind, 2, p, 4096, link, [("bar", fb), ("psi", fp)], [T], sigma=sg, **kw)
            row[("Link" if link else "Platz") + ("" if sg == -1 else " Gegenvorzeichen")] = (
                max(a for a, _ in vals), min(b for _, b in vals))
    cal[lab] = row
    say(f"   {lab:22s} Link {row['Link'][1]:+.2e}  Platz {row['Platz'][1]:+.2e}  Link mit Gegenvorzeichen "
        f"{row['Link Gegenvorzeichen'][1]:+.2e}  (antihermitesch max {max(v[0] for v in row.values()):.0e})")
check("Wilson r = 1 Link-RP; r = 2, Gegenvorzeichen und Platzreflexion verletzt",
      cal["Wilson r=1"]["Link"][1] > TOL and cal["Wilson r=2"]["Link"][1] < -0.05
      and cal["Wilson r=1"]["Link Gegenvorzeichen"][1] < -0.5 and cal["Wilson r=1"]["Platz"][1] < -0.05)
check("Voller Overlap-Dirac (m0 = 0,5/1/1,5; m = 0 und 0,3) Link-RP in 2D",
      all(cal[k]["Link"][1] > TOL for k in cal if k.startswith("Overlap")))
befund("Platzreflexion (naiv, Felder bei t = 0 fest) ist für Wilson, Overlap und Weyl verletzt",
       all(cal[k]["Platz"][1] < -0.05 for k in cal))
ZAHLEN["kalibrierung_2d"] = {k: {kk: list(vv) for kk, vv in v.items()} for k, v in cal.items()}

say("\n2b Lüscher-Weyl, lineare Sektoren (psi-quer P_+, Anker P_- psi), Link, 2D-Impulsraum")
say("   Scan: m0, L_t, L_s mit Holonomie phi (p = (2 pi n + phi)/L_s), Zeitschichten {d, .., d+nb-1} mit d = 1, 2, 3.")
weyl_tab = {}
t1 = time.time()
r_vor = dict(RAUSCH)
for m0 in (0.5, 1.0, 1.5, 1.9):
    for Lt in (8, 16, 32, 64, 4096):
        worst, nt = 1.0, 0
        Ts = [list(range(d, d + min(4, Lt // 2 - d + 1))) for d in (1, 2, 3) if d <= Lt // 2]
        for Ls in (8, 16, 32, 64):
            for phi in (0.0, 0.37 * np.pi, np.pi):
                for nn in range(Ls):
                    p = np.array([(2 * np.pi * nn + phi) / Ls])
                    res = os_scan("weyl", 2, p, Lt, True, [("bar", Ep2), ("psi", Em2)], Ts, m0=m0)
                    nt += len(res)
                    worst = min([worst] + [b for _, b in res])
        weyl_tab[(m0, Lt)] = (worst, nt)
for m0 in (0.5, 1.0, 1.5, 1.9):
    say(f"   m0 = {m0}: kleinster normierter Eigenwert je L_t = 8/16/32/64/4096: "
        + "/".join(f"{weyl_tab[(m0, Lt)][0]:+.1e}" for Lt in (8, 16, 32, 64, 4096))
        + f"  ({sum(weyl_tab[(m0, Lt)][1] for Lt in (8, 16, 32, 64, 4096))} Formen)")
n_rausch = RAUSCH["rauschfaelle"] - r_vor["rauschfaelle"]
n_formen = sum(v[1] for v in weyl_tab.values())
say(f"   {n_rausch} von {n_formen} Formen haben einen negativen normierten Eigenwert, der roh höchstens "
    f"{RAUSCH['max_roh_durch_rausch']:.2f} x Rauschgrenze (1e-15 max|G|) groß ist: Rundung bei Einträgen ~1e-11")
say(f"   (normiert wurden daraus bis {RAUSCH['min_normiert']:+.1e}, roh {RAUSCH['roh_dazu']:+.1e}). Signifikant negativ: "
    f"{'keine' if min(v[0] for v in weyl_tab.values()) >= 0 else 'ja'}.")
say(f"   (Rechenzeit {time.time() - t1:.0f} s)")
weyl_min = min(v[0] for v in weyl_tab.values())
check("Lüscher-Weyl Link-RP: alle linearen Formen >= 0 (über Rundung) für alle m0, L_t, L_s, Holonomien, Abstände",
      weyl_min > TOL and RAUSCH["max_roh_durch_rausch"] < 1e3)
ZAHLEN["weyl_2d_scan"] = {f"m0={k[0]},Lt={k[1]}": v[0] for k, v in weyl_tab.items()}
ZAHLEN["weyl_2d_scan_rundung"] = {"formen": n_formen, "nur_rundung_negativ": n_rausch,
                                  "max_roh_durch_rauschgrenze": RAUSCH["max_roh_durch_rausch"]}

say("\n2c Spektraldarstellung des Weylblocks (L_t = unendlich): Pole und Schnitt, Gewichte komprimiert")
spec_w = {}
for m0 in (0.5, 1.0, 1.5):
    for p in (-3.0, -2.0, -1.0, -0.3, -0.05, 0.05, 0.3, 1.0, 2.0, 2.5, 3.1):
        a = spectral("weyl", 2, p, "bar", Ep2, m0=m0)
        b = spectral("weyl", 2, p, "psi", Em2, m0=m0)
        spec_w[f"m0={m0},p={p}"] = dict(zst=a["zst"], pole_bar=a["poles"], err=max(a["err"], b["err"]),
                                         wmin=min(a["wmin"], b["wmin"]), aher=max(a["aher"], b["aher"]))
for key in ("m0=1.0,p=0.3", "m0=1.0,p=-0.3", "m0=1.0,p=2.5", "m0=0.5,p=1.0", "m0=1.5,p=-2.0"):
    s_ = spec_w[key]
    say(f"   {key:15s} z* = {s_['zst']:+.3f}, Pole {np.round(s_['pole_bar'], 3)}, Rekonstruktion {s_['err']:.0e},"
        f" kleinstes Gewicht {s_['wmin']:+.1e}, antihermitesch {s_['aher']:.0e}")
sw_ok = all(v["err"] < 1e-8 and v["wmin"] > TOL and v["aher"] < 1e-8 for v in spec_w.values())
say(f"   {len(spec_w)} Fälle: Rekonstruktion max {max(v['err'] for v in spec_w.values()):.0e}, kleinstes Gewicht "
    f"{min(v['wmin'] for v in spec_w.values()):+.1e}. Für p > 0 trägt die psi-quer-Seite nur den Schnitt (Gitterskala),")
say("   der Pol des Weylteilchens erscheint auf der psi-Seite (und umgekehrt): eine chirale Mode.")
check("Spektralgewichte des Weylblocks positiv-semidefinit, Rekonstruktion < 1e-8 (also RP für alle Abstände)", sw_ok)
ZAHLEN["weyl_2d_spektral"] = spec_w

say("\n2d Ortsraum 2D: ganze positive Halbzeit (t = 1 .. L_t/2, alle x), psi-quer P_+ und Anker P_- psi")
real_tab = {}
for L in (8, 12, 16, 24):
    for phi, bg in ((0.0, 0.0), (0.7, 0.0), (0.7, 0.8)):
        DW, co = build_DW(2, L, L, phi=phi, bg=bg)
        D, G5, hmin = overlap(DW, 2)
        R = refl(2, L, L)
        n = D.shape[0]
        SW = np.linalg.inv(D) @ ((np.eye(n) + G5) / 2)
        Eb, _ = half_basis(2, co, L, Ep2)
        Ea, _ = half_basis(2, co, L, Em2)
        Efu, _ = half_basis(2, co, L, Ef2)
        rb = assess(os_real(SW, R, Eb, "bar"), 1)
        ra = assess(os_real(SW, R, Ea, "psi"), 1)
        rf = assess(os_real(SW, R, Efu, "psi"), 2)
        real_tab[(L, phi, bg)] = dict(bar=rb[:2], anker=ra[:2], voll_psi=rf[:2], hmin=hmin)
for (L, phi, bg), v in real_tab.items():
    say(f"   L = {L:2d}, phi = {phi}, Hintergrund {bg}: psi-quer {v['bar'][1]:+.1e} (antiherm. {v['bar'][0]:.0e}), "
        f"Anker {v['anker'][1]:+.1e}, volle psi-Komponenten {v['voll_psi'][1]:+.1e} (antiherm. {v['voll_psi'][0]:.0e})")
check("Ortsraum 2D: lineare Weyl-OS-Matrizen auf der ganzen Halbzeit hermitesch und >= 0 (L = 8 .. 24, mit Hintergrund)",
      all(v["bar"][1] > TOL and v["anker"][1] > TOL and v["voll_psi"][1] > TOL and v["bar"][0] < 1e-10
          and v["anker"][0] < 1e-10 for v in real_tab.values()))
say("   Lineare Sektoren sehen den Kontaktteil nicht (Link: theta x != y), daher auch die vollen Komponenten >= 0.")
ZAHLEN["weyl_2d_ortsraum"] = {f"L={k[0]},phi={k[1]},bg={k[2]}": v for k, v in real_tab.items()}

say("\n2e Bilineare Monome psi-quer(x) e * u^+ psi(y) samt Eins (Ladung 0), Wick mit S_W, Zeitschichten t = 1, 2")


def bilinear(G, R, Eb, Ep_):
    """M_FG = <Theta F><G> - [u'^+ G(y', theta y) g_t u][e^+ g_t G(theta x, x') e'] mit Eins; sigma sigma' = +1."""
    a = -(Eb.conj().T @ (R @ G @ R) @ Ep_)                      # <Theta F>, F = (i, j)
    b = -(Ep_.conj().T @ G @ Eb).T                              # <G>
    Mb0 = Eb.conj().T @ (R @ G) @ Eb                            # (i, i')
    Mp0 = Ep_.conj().T @ (G @ R) @ Ep_                          # (j', j)
    conn = -np.einsum("ik,lj->ijkl", Mb0, Mp0).reshape(a.size, a.size)
    M = np.zeros((a.size + 1, a.size + 1), complex)
    M[0, 0] = 1.0
    M[0, 1:] = b.ravel()
    M[1:, 0] = a.ravel()
    M[1:, 1:] = np.outer(a.ravel(), b.ravel()) + conn
    kron_ok = float(np.abs(conn - np.kron(-Mb0, Mp0.T)).max())
    asym = float(np.abs(a - b.conj()).max())
    return M, asym, kron_ok


bil = {}
for L, bg in ((8, 0.0), (12, 0.0), (12, 0.8)):
    DW, co = build_DW(2, L, L, phi=0.7, bg=bg)
    D, G5, _ = overlap(DW, 2)
    R = refl(2, L, L)
    n = D.shape[0]
    SW = np.linalg.inv(D) @ ((np.eye(n) + G5) / 2)
    Eb, _ = half_basis(2, co, L, Ep2, slices=(1, 2))
    for lab, P in (("voll (mit Kontaktkomponente P_+ psi)", Ef2), ("Anker (P_- psi)", Em2)):
        Ep_, _ = half_basis(2, co, L, P, slices=(1, 2))
        M, asym, kr = bilinear(SW, R, Eb, Ep_)
        aher, mn = assess(M, 1)[:2]
        bil[(L, bg, lab)] = dict(dim=M.shape[0], asym=asym, aher=aher, min=mn, kron=kr)
        say(f"   L = {L:2d}, Hintergrund {bg}, {lab:36s}: {M.shape[0]:4d} Monome, |<Theta F> - <F>*| max {asym:.3f},"
            f" antiherm. {aher:.1e}, min {mn:+.1e}")
full_asym = [v["asym"] for k, v in bil.items() if k[2].startswith("voll")]
anc = [v for k, v in bil.items() if k[2].startswith("Anker")]
check("Wick-Zerlegung: verbundener Teil = Tensorprodukt der linearen Formen (Rest < 1e-12)",
      all(v["kron"] < 1e-12 for v in bil.values()))
check("Volle Komponenten: Kontaktasymmetrie |<Theta F> - <F>*| = 1/2 exakt (negatives Ergebnis festgestellt)",
      all(abs(x - 0.5) < 1e-10 for x in full_asym))
check("Ankeralgebra: Bilineare mit Eins hermitesch und >= 0", all(v["aher"] < 1e-10 and v["min"] > TOL for v in anc))
befund("Volle Komponentenalgebra der Weylfelder ist RP (Bilineare mit Eins)", False)
say("   Die Verletzung sitzt nur in Komposita psi-quer P_+ psi am selben Punkt (Bewegungsgleichungsfeld, <.> = -1/2,")
say("   gespiegelt 0); sie hängt nicht von L, m0 oder dem Hintergrund ab und erreicht keine getrennten Punkte.")
ZAHLEN["bilineare_2d"] = {f"L={k[0]},bg={k[1]},{k[2]}": v for k, v in bil.items()}

# =========================================================================================== Teil 3
say("\nTeil 3: Gedrehte Phase. Zwei Lüscher-Weylfelder psi_L = P^_- psi, psi_R = P^_+ psi mit Paarung")
say("   m (e^{i theta} psi-quer_L psi_R + e^{-i theta} psi-quer_R psi_L) = psi-quer m e^{i theta g5} (1 - D/2) psi (GJ-Form, Z16/Z19).")
say("   theta ist genau die GW-chirale Drehung. Theta_D wirkt auf (psi, psi-quer) wie Z19, Theta_W auf (chi, psi-quer),")
say("   chi = (1 - D/2) psi; Theta_W(psi-quer P_+) = -g_t chi_L, Theta_W(psi-quer P_-) = -g_t chi_R.")
rot_tab = {}
for m in (0.1, 0.3, 0.6):
    for th in (0.0, 0.3, np.pi / 2, np.pi, 4.0):
        for Lt in (16, 4096):
            for kind in ("ov", "chi"):
                vals = []
                for p in PS2 + [np.array([0.05]), np.array([-1.3])]:
                    vals += os_scan(kind, 2, p, Lt, True, [("bar", Ef2), ("psi", Ef2)], [[1, 2, 3], [2, 3, 4]], m0=1.0, m=m, th=th)
                rot_tab[(m, round(th, 3), Lt, kind)] = (max(a for a, _ in vals), min(b for _, b in vals))
for m in (0.1, 0.3, 0.6):
    say(f"   m = {m}: theta = 0/0,3/pi/2/pi/4,0 (L_t = 4096)")
    for kind, lab in (("ov", "Theta_D"), ("chi", "Theta_W")):
        say(f"      {lab}: min " + "/".join(f"{rot_tab[(m, round(th, 3), 4096, kind)][1]:+.1e}" for th in (0.0, 0.3, np.pi / 2, np.pi, 4.0))
            + "; antiherm. " + "/".join(f"{rot_tab[(m, round(th, 3), 4096, kind)][0]:.1e}" for th in (0.0, 0.3, np.pi / 2, np.pi, 4.0)))
stab = {Lt: rot_tab[(0.3, round(np.pi / 2, 3), Lt, "ov")] for Lt in (16, 4096)}
say(f"   Theta_D bei m = 0,3, theta = pi/2: L_t = 16 {stab[16][1]:+.2e} (antiherm. {stab[16][0]:.3f}), "
    f"L_t = 4096 {stab[4096][1]:+.2e} ({stab[4096][0]:.3f}): stabil in L, also kein Endlichkeitsartefakt.")
ah = np.array([rot_tab[(m, round(np.pi / 2, 3), 4096, "ov")][0] for m in (0.1, 0.3, 0.6)])
slope_m = float(np.polyfit(np.log([0.1, 0.3, 0.6]), np.log(ah), 1)[0])
say(f"   Antihermitescher Anteil bei pi/2 gegen m = 0,1/0,3/0,6: {ah[0]:.3f}/{ah[1]:.3f}/{ah[2]:.3f}, Exponent {slope_m:.2f}:")
say("   linear in a m. Grund (exakt): M_D = B M_W mit B = 1 - (m/2) e^{-i theta g5}, einseitig, ortslokal.")
check("Theta_W: positiv für alle theta, m, L_t (Formen hermitesch)",
      all(v[1] > TOL and v[0] < 1e-9 for k, v in rot_tab.items() if k[3] == "chi"))
check("Theta_D: positiv bei theta = 0 und pi, bei theta = 0,3/pi/2/4,0 nicht hermitesch (negatives Ergebnis)",
      all(rot_tab[(m, 0.0, Lt, 'ov')][1] > TOL and rot_tab[(m, round(np.pi, 3), Lt, 'ov')][1] > TOL for m in (0.1, 0.3, 0.6) for Lt in (16, 4096))
      and all(rot_tab[(m, round(th, 3), Lt, 'ov')][0] > 1e-2 for m in (0.1, 0.3, 0.6) for th in (0.3, np.pi / 2, 4.0) for Lt in (16, 4096)))
check("Theta_D-Verletzung linear in m (Exponent 1 auf 0,15) und stabil in L_t (16 gegen 4096 auf 5 %)",
      abs(slope_m - 1) < 0.15 and abs(stab[16][0] / stab[4096][0] - 1) < 0.05)
befund("Lüscher-gepaarte Weylfelder bei gedrehter Phase RP unter Theta_D (Reflexion auf psi)", False)
befund("Lüscher-gepaarte Weylfelder bei gedrehter Phase RP unter Theta_W (Reflexion auf chi = (1 - D/2) psi)", True)
# exakte Struktur M_D = B M_W und theta-Unabhängigkeit der Theta_W-Spektren
p = np.array([0.7])
T = [1, 2, 3]
Gd = G_tau("ov", 2, p, taus_for(T, True), 4096, m0=1.0, m=0.3, th=1.1)
Gw = G_tau("chi", 2, p, taus_for(T, True), 4096, m0=1.0, m=0.3, th=1.1)
Gw0 = G_tau("chi", 2, p, taus_for(T, True), 4096, m0=1.0, m=0.3, th=0.0)
Md = os_form(Gd, 2, T, True, "bar", -1, Ef2)
Mw = os_form(Gw, 2, T, True, "bar", -1, Ef2)
Mw0 = os_form(Gw0, 2, T, True, "bar", -1, Ef2)
Bm = np.kron(np.eye(3), np.eye(2) - 0.15 * rot(2, -1.1))
rest_B = float(np.abs(Md - Bm @ Mw).max())
rest_U = float(np.abs(np.linalg.eigvalsh((Mw + Mw.conj().T) / 2) - np.linalg.eigvalsh((Mw0 + Mw0.conj().T) / 2)).max())
say(f"   Gegenprobe p = 0,7, theta = 1,1: |M_D - B M_W| = {rest_B:.1e}; Spektrum von M_W(theta) gegen M_W(0): {rest_U:.1e}.")
check("Exakt: M_D = B M_W und Spektrum von M_W unabhängig von theta (Reste < 1e-12)", rest_B < 1e-12 and rest_U < 1e-12)
spec_r = {}
for kind in ("chi", "ov"):
    for th in (0.0, 1.2, np.pi / 2):
        for p in (0.3, -0.3, 1.0, 2.5, -2.0):
            a = spectral(kind, 2, p, "bar", Ef2, m0=1.0, m=0.3, th=th)
            b = spectral(kind, 2, p, "psi", Ef2, m0=1.0, m=0.3, th=th)
            spec_r[f"{kind},theta={th:.2f},p={p}"] = dict(err=max(a["err"], b["err"]), wmin=min(a["wmin"], b["wmin"]),
                                                          aher=max(a["aher"], b["aher"]), poles=a["poles"])
sr_chi = [v for k, v in spec_r.items() if k.startswith("chi")]
sr_ov = [v for k, v in spec_r.items() if k.startswith("ov") and "theta=0.00" not in k]
say(f"   Spektral (L_t = unendlich), m = 0,3: Theta_W kleinstes Gewicht {min(v['wmin'] for v in sr_chi):+.1e}, antiherm. "
    f"{max(v['aher'] for v in sr_chi):.0e}; Theta_D bei theta != 0 kleinstes Gewicht {min(v['wmin'] for v in sr_ov):+.1e},"
    f" antiherm. {max(v['aher'] for v in sr_ov):.2f}; Rekonstruktion max {max(v['err'] for v in spec_r.values()):.0e}.")
check("Spektral: Theta_W-Gewichte >= 0 bei jedem theta; Theta_D-Gewichte bei theta != 0 nicht hermitesch",
      all(v["wmin"] > TOL and v["aher"] < 1e-8 and v["err"] < 1e-8 for v in sr_chi)
      and all(v["aher"] > 1e-2 for v in sr_ov))
# Ortsraum mit Hintergrund
rr, symm = {}, {}
for L, bg in ((12, 0.0), (16, 0.8)):
    DW, co = build_DW(2, L, L, phi=0.7, bg=bg)
    D, G5, _ = overlap(DW, 2)
    R = refl(2, L, L)
    n = D.shape[0]
    Efu, _ = half_basis(2, co, L, Ef2)
    for th in (0.0, np.pi / 2):
        Eth = np.kron(np.eye(n // 2), rot(2, th))
        S = np.linalg.inv(D + 0.3 * Eth @ (np.eye(n) - D / 2))
        Gc = (np.eye(n) - D / 2) @ S
        for lab, G in (("Theta_D", S), ("Theta_W", Gc)):
            rb = assess(os_real(G, R, Efu, "bar"), 2)
            rp_ = assess(os_real(G, R, Efu, "psi"), 2)
            rr[(L, bg, round(th, 3), lab)] = (max(rb[0], rp_[0]), min(rb[1], rp_[1]))
            symm[(L, bg, round(th, 3), lab)] = float(np.abs(R @ G @ R - G.conj().T).max())
for k, v in rr.items():
    say(f"   Ortsraum L = {k[0]}, Hintergrund {k[1]}, theta = {k[2]}, {k[3]}: min {v[1]:+.1e}, antiherm. {v[0]:.1e}")
say("   Reflexionssymmetrie des Zustands auch am selben Punkt, max |R G R - G^+|: " + ", ".join(
    f"L={k[0]},bg={k[1]},theta={k[2]},{k[3]}: {v:.1e}" for k, v in symm.items()))
say("   Für einen quasifreien Zustand reichen diese Symmetrie und die linearen Formen für alle Monome (Wick, Teil 2e).")
check("Theta_W-Zustand reflexionssymmetrisch an allen Punkten (< 1e-12), Theta_D-Zustand bei pi/2 nicht (> 1e-2)",
      all(v < 1e-12 for k, v in symm.items() if k[3] == "Theta_W")
      and all(v > 1e-2 for k, v in symm.items() if k[3] == "Theta_D" and k[2] > 0))
check("Ortsraum 2D (ganze Halbzeit, mit Hintergrund): Theta_W >= 0 bei 0 und pi/2; Theta_D bei pi/2 nicht hermitesch",
      all(v[1] > TOL and v[0] < 1e-9 for k, v in rr.items() if k[3] == "Theta_W")
      and all(v[0] > 1e-2 for k, v in rr.items() if k[3] == "Theta_D" and k[2] > 0))
ZAHLEN["gedreht_2d"] = {f"m={k[0]},theta={k[1]},Lt={k[2]},{'Theta_D' if k[3] == 'ov' else 'Theta_W'}": list(v)
                        for k, v in rot_tab.items()}
ZAHLEN["gedreht_2d_antiherm_gegen_m"] = {"m": [0.1, 0.3, 0.6], "antiherm": ah.tolist(), "exponent": slope_m}
ZAHLEN["gedreht_2d_spektral"] = spec_r
ZAHLEN["gedreht_2d_ortsraum"] = {f"L={k[0]},bg={k[1]},theta={k[2]},{k[3]}": list(v) for k, v in rr.items()}
ZAHLEN["gedreht_2d_reflexionssymmetrie"] = {f"L={k[0]},bg={k[1]},theta={k[2]},{k[3]}": v for k, v in symm.items()}
# Majoranadrehung X aus Z19
K = 2 * np.pi * (np.arange(64) + 0.5) / 64
absX = []
for px in K:
    _, Dk = ops_k(2, np.exp(1j * K), np.array([px]))
    X = DIM[2]["g5"][None] @ (DIM[2]["I"][None] - Dk / 2)
    absX.append(np.abs(np.linalg.eigvalsh((X + np.conj(np.swapaxes(X, 1, 2))) / 2)))
absX = np.array(absX)
say(f"   Zum Vergleich die Majoranadrehung X = g5 (1 - D/2) aus Z19: |X| zwischen {absX.min():.3f} und {absX.max():.3f}")
say("   (X^2 = 1 - D^+D/4). Ein ortslokaler Verflechter K X = Gamma K mit Gamma^2 = 1 müsste alle Eigenwerte |X| != 1")
say("   vernichten, also fast die ganze Brillouinzone: Für die flachen Majoranaformen (X, Xq) gibt es kein Theta_W.")
check("Majoranadrehung X nicht quantisiert (|X| < 0,5 auf Teilen der Zone), also kein ortslokaler Verflechter",
      absX.min() < 0.5 and absX.max() > 0.99)
ZAHLEN["majorana_X_spektrum_2d"] = [float(absX.min()), float(absX.max())]

# =========================================================================================== Teil 4
say("\nTeil 4: 4D")
Ef4, Ep4, Em4 = DIM[4]["Ef"], DIM[4]["E+"], DIM[4]["E-"]
PS4 = [np.array(p) for p in ((0.3, 0.0, 0.0), (1.0, 0.5, 0.2), (2.5, 2.0, 0.3), (np.pi, 0.4, 0.0),
                             (2.8, 2.9, 3.0), (0.1, 0.05, 0.02))]
say("4a Kalibrierung und Z19-Gegenprobe (L_t = 4096, Z19-Impulse)")
c4 = {}
for lab, kind, kw in (("Wilson r=1", "wilson", dict(r=1.0, m=0.3)), ("Overlap m0=1", "ov", dict(m0=1.0)),
                      ("Overlap m0=1, m=0,3", "ov", dict(m0=1.0, m=0.3))):
    vals = []
    for p in PS4:
        vals += os_scan(kind, 4, p, 4096, True, [("bar", Ef4), ("psi", Ef4)], [[1, 2, 3, 4]], **kw)
    c4[lab] = (max(a for a, _ in vals), min(b for _, b in vals))
    say(f"   {lab:22s} Link min {c4[lab][1]:+.2e}, antiherm. {c4[lab][0]:.0e}")
z19v = []
for p in PS4[:4]:
    z19v += os_scan("ov", 4, p, 4096, True, [("bar", Ef4), ("psi", Ef4)], [[1, 2, 3]], m0=1.0, m=0.3, th=np.pi / 2)
z19_ah, z19_min = max(a for a, _ in z19v), min(b for _, b in z19v)
say(f"   Z19 ov_GJ theta = pi/2 nachgerechnet: min {z19_min:+.3e} (Z19: -2,86e-3), antiherm. {z19_ah:.2f} (Z19: 3e-1)")
check("4D: Wilson r = 1 und Overlap Link-RP; Z19-Wert ov_GJ (pi/2) reproduziert",
      c4["Wilson r=1"][1] > TOL and c4["Overlap m0=1"][1] > TOL and c4["Overlap m0=1, m=0,3"][1] > TOL
      and abs(z19_min + 2.86e-3) < 2e-5 and 0.2 < z19_ah < 0.4)
say("\n4b Lüscher-Weyl 4D, lineare Sektoren (P_+ bzw. Anker P_-, je zwei Komponenten), Link")
w4 = {}
t1 = time.time()
grid4 = [np.array([(2 * np.pi * a + ph) / 4, (2 * np.pi * b + ph) / 4, (2 * np.pi * c + ph) / 4])
         for ph in (0.0, 0.37 * np.pi) for a in range(4) for b in range(4) for c in range(4)]
for m0 in (0.5, 1.0, 1.5):
    for Lt in (8, 16, 4096):
        worst = 1.0
        Ts = [[1, 2, 3, 4], [2, 3, 4]] if Lt > 8 else [[1, 2, 3, 4]]
        for p in grid4 + PS4:
            res = os_scan("weyl", 4, p, Lt, True, [("bar", Ep4), ("psi", Em4)], Ts, m0=m0)
            worst = min([worst] + [b for _, b in res])
        w4[(m0, Lt)] = worst
    say(f"   m0 = {m0}: kleinster normierter Eigenwert L_t = 8/16/4096: " + "/".join(f"{w4[(m0, Lt)]:+.1e}" for Lt in (8, 16, 4096))
        + f" ({len(grid4) + len(PS4)} Raumimpulse, Holonomie 0 und 0,37 pi)")
say(f"   (Rechenzeit {time.time() - t1:.0f} s)")
check("4D Lüscher-Weyl Link-RP für m0 = 0,5/1/1,5, L_t = 8/16/4096", min(w4.values()) > TOL)
ZAHLEN["weyl_4d_scan"] = {f"m0={k[0]},Lt={k[1]}": v for k, v in w4.items()}
rot4 = {}
for th in (0.0, np.pi / 2, 4.0):
    for kind in ("ov", "chi"):
        vals = []
        for p in PS4:
            vals += os_scan(kind, 4, p, 4096, True, [("bar", Ef4), ("psi", Ef4)], [[1, 2, 3]], m0=1.0, m=0.3, th=th)
        rot4[(round(th, 3), kind)] = (max(a for a, _ in vals), min(b for _, b in vals))
say("   Gedrehte Phase 4D (m = 0,3), theta = 0/pi/2/4,0: Theta_D min " + "/".join(f"{rot4[(round(t, 3), 'ov')][1]:+.1e}" for t in (0.0, np.pi / 2, 4.0))
    + "; Theta_W min " + "/".join(f"{rot4[(round(t, 3), 'chi')][1]:+.1e}" for t in (0.0, np.pi / 2, 4.0)))
check("4D gedrehte Phase: Theta_W >= 0 und hermitesch für alle theta; Theta_D bei pi/2 und 4,0 nicht hermitesch",
      all(v[1] > TOL and v[0] < 1e-9 for k, v in rot4.items() if k[1] == "chi")
      and rot4[(round(np.pi / 2, 3), "ov")][0] > 1e-2 and rot4[(4.0, "ov")][0] > 1e-2)
ZAHLEN["gedreht_4d"] = {f"theta={k[0]},{'Theta_D' if k[1] == 'ov' else 'Theta_W'}": list(v) for k, v in rot4.items()}
sp4 = {}
for p in PS4[:5]:
    a = spectral("weyl", 4, p, "bar", Ep4, m0=1.0)
    b = spectral("weyl", 4, p, "psi", Em4, m0=1.0)
    c = spectral("chi", 4, p, "bar", Ef4, m0=1.0, m=0.3, th=1.2)
    sp4[f"p=({p[0]:.2f},{p[1]:.2f},{p[2]:.2f})"] = dict(err=max(a["err"], b["err"], c["err"]),
                                                        wmin=min(a["wmin"], b["wmin"], c["wmin"]),
                                                        aher=max(a["aher"], b["aher"], c["aher"]))
say(f"   Spektral 4D (Weylblock, Theta_W mit theta = 1,2): Rekonstruktion max {max(v['err'] for v in sp4.values()):.0e}, "
    f"kleinstes Gewicht {min(v['wmin'] for v in sp4.values()):+.1e}")
check("4D Spektralgewichte (Weyl, Theta_W gedreht) >= 0, Rekonstruktion < 1e-8",
      all(v["err"] < 1e-8 and v["wmin"] > TOL and v["aher"] < 1e-8 for v in sp4.values()))
ZAHLEN["spektral_4d"] = sp4

say("\n4c Kleine 4D-Ortsraumgitter L_t x 4^3 (Holonomie 0,7), ganze positive Halbzeit")
r4 = {}
t1 = time.time()
for Lt in (4, 8):
    DW, co = build_DW(4, Lt, 4, phi=0.7)
    D, G5, hmin = overlap(DW, 4)
    R = refl(4, Lt, 4)
    n = D.shape[0]
    Ib = np.eye(n)
    PP = (Ib + G5) / 2
    Di = np.linalg.inv(D)
    SW = Di @ PP
    gw = float(np.abs(G5 @ D + D @ G5 - D @ G5 @ D).max())
    rdr = float(np.abs(R @ D @ R - D.conj().T).max())
    cont = float(np.abs(SW - (PP / 2 + (Ib - PP) @ Di @ PP)).max())
    Eb, _ = half_basis(4, co, Lt, Ep4)
    Ea, _ = half_basis(4, co, Lt, Em4)
    Efu, _ = half_basis(4, co, Lt, Ef4)
    rb = assess(os_real(SW, R, Eb, "bar"), 2)
    ra = assess(os_real(SW, R, Ea, "psi"), 2)
    Eth = np.kron(np.eye(n // 4), rot(4, np.pi / 2))
    S = np.linalg.inv(D + 0.3 * Eth @ (Ib - D / 2))
    Gc = (Ib - D / 2) @ S
    rD = assess(os_real(S, R, Efu, "bar"), 4)
    rW = assess(os_real(Gc, R, Efu, "bar"), 4)
    r4[Lt] = dict(GW=gw, RDR=rdr, kontakt=cont, hmin=hmin, weyl_bar=rb[:2], weyl_anker=ra[:2], ThetaD_pi2=rD[:2], ThetaW_pi2=rW[:2])
    say(f"   L_t = {Lt}: GW {gw:.0e}, RDR-D^+ {rdr:.0e}, Kontaktzerlegung {cont:.0e}; Weyl psi-quer {rb[1]:+.1e}, Anker {ra[1]:+.1e};"
        f" gedreht pi/2: Theta_D antiherm. {rD[0]:.2f}, Theta_W min {rW[1]:+.1e} (antiherm. {rW[0]:.0e})")
say(f"   (Rechenzeit {time.time() - t1:.0f} s)")
check("4D Ortsraum: Identitäten exakt, Weyl-Formen >= 0, Theta_W gedreht >= 0, Theta_D gedreht nicht hermitesch",
      all(v["GW"] < 1e-11 and v["RDR"] < 1e-11 and v["kontakt"] < 1e-11 and v["weyl_bar"][1] > TOL
          and v["weyl_anker"][1] > TOL and v["ThetaW_pi2"][1] > TOL and v["ThetaW_pi2"][0] < 1e-9
          and v["ThetaD_pi2"][0] > 1e-2 for v in r4.values()))
ZAHLEN["ortsraum_4d"] = {str(k): v for k, v in r4.items()}

# =========================================================================================== Teil 5
say("\nTeil 5: Einordnung ins Z24-Trilemma (geschützte Naht / positiver Transfer / flache Phase ohne Abstimmung)")
K16 = 2 * np.pi * (np.arange(16) + 0.5) / 16
Vgj, Vwi = [], []
for th in np.linspace(0, 2 * np.pi, 9):
    a = b = 0.0
    for px in K16:
        DWk, Dk = ops_k(2, np.exp(1j * K16), np.array([px]))
        a += np.linalg.slogdet(Dk + 0.3 * rot(2, th)[None] @ (DIM[2]["I"][None] - Dk / 2))[1].sum()
        b += np.linalg.slogdet(DWk + 0.3 * rot(2, th)[None])[1].sum()
    Vgj.append(-a / 256)
    Vwi.append(-b / 256)
sp_gj, sp_wi = float(np.ptp(Vgj)), float(np.ptp(Vwi))
say(f"   5a Phasenpotential (2D, 16^2, m = 0,3): GJ-Paarung Spannweite {sp_gj:.1e} (flach, exakte GW-Drehung);")
say(f"      Wilson-Dirac mit e^{{i theta g5}}: Spannweite {sp_wi:.3f} (Barriere).")
check("GJ-gepaarte Lüscher-Weylfelder: Phasenpotential flach (< 1e-12), Wilson nicht", sp_gj < 1e-12 and sp_wi > 0.1)
# 5b Index (Z19-Querschnitt, 4D, k3 = k4 = 0), Dirac-GJ
I4, g54 = DIM[4]["I"], DIM[4]["g5"]
P4p, P4m = (I4 + g54) / 2, (I4 - g54) / 2
SIG = -1j * gam4[2] @ gam4[3]
L2 = 16
nsite = L2 * L2


def transverse(L, r=1.0):
    idx = lambda x, y: x * L + y
    N = np.zeros((4 * L * L, 4 * L * L), complex)
    Wop = np.zeros((L * L, L * L))
    for x in range(L):
        for y in range(L):
            i = idx(x, y)
            Wop[i, i] += 2 * r
            for mu, (dx, dy) in enumerate(((1, 0), (0, 1))):
                if x + dx < L and y + dy < L:
                    j = idx(x + dx, y + dy)
                    N[4 * i:4 * i + 4, 4 * j:4 * j + 4] += gam4[mu] / 2
                    N[4 * j:4 * j + 4, 4 * i:4 * i + 4] -= gam4[mu] / 2
                    Wop[i, j] -= r / 2
                    Wop[j, i] -= r / 2
    return N, Wop


Nt, Wt = transverse(L2)
G5b = np.kron(np.eye(nsite), g54)
Hw = G5b @ (Nt + np.kron(Wt - np.eye(nsite), I4))
Hw = (Hw + Hw.conj().T) / 2
wH, UH = np.linalg.eigh(Hw)
Dq = np.eye(4 * nsite) + G5b @ ((UH * np.sign(wH)) @ UH.conj().T)
Ibq = np.eye(4 * nsite)
cq = (L2 - 1) / 2
xs, ys = np.meshgrid(np.arange(L2), np.arange(L2), indexing="ij")
phi_q = np.arctan2(ys - cq, xs - cq).ravel()
rad = np.hypot(xs - cq, ys - cq).ravel()
prof = np.tanh(rad / 1.5)
g5hq = G5b @ (Ibq - Dq)
Phq_p, Phq_m = (Ibq + g5hq) / 2, (Ibq - g5hq) / 2
PPq, PMq = np.kron(np.eye(nsite), P4p), np.kron(np.eye(nsite), P4m)


def chiral_count(Deff, R=3.5, tol=0.25):
    wS, US = np.linalg.eigh(SIG)
    res = {}
    core = np.repeat(rad <= R, 4)
    for s in (1, -1):
        B = np.kron(np.eye(nsite), US[:, wS * s > 0])
        _, sv, vh = np.linalg.svd(B.conj().T @ Deff @ B)
        V = (B @ vh.conj().T)[:, sv < tol]
        res[s] = int(np.sum(np.sum(np.abs(V[core]) ** 2, axis=0) > 0.5))
    return res[1] - res[-1]


idx_gj = {}
for nw in (1, 0, 2, -1):
    Eph = np.kron(np.diag(prof * np.exp(1j * nw * phi_q)), I4)
    idx_gj[nw] = chiral_count(Dq + 0.6 * (PPq @ Eph @ Phq_p + PMq @ Eph.conj() @ Phq_m))
say(f"   5b Index des GJ-Wirbels im Z19-Querschnitt (16 x 16, k3 = k4 = 0), Windung 1/0/2/-1: netto "
    f"{idx_gj[1]}/{idx_gj[0]}/{idx_gj[2]}/{idx_gj[-1]} chirale Moden auf der Naht.")
check("GJ-Wirbel: Index = Windung (1, 0, 2, -1), also geschützte Naht", idx_gj == {1: 1, 0: 0, 2: 2, -1: -1})
# 5c Nichtlokalität der Variablen chi: D_c = D (1 - D/2)^-1 (Wirkung psi-quer (D_c + m e^{i theta g5}) chi)
Lk = 1024
pk_ = 2 * np.pi * (np.arange(Lk) + 0.5) / Lk
xs_ = np.arange(Lk)
tails = {}
for m0 in (0.5, 1.0, 1.5):
    Kc = np.zeros((Lk, 2, 2), complex)
    Kd = np.zeros((Lk, 2, 2), complex)
    for j, px in enumerate(pk_):
        _, Dk = ops_k(2, np.exp(1j * pk_), np.array([px]), m0=m0)
        Kc[j] = (Dk @ np.linalg.inv(DIM[2]["I"][None] - Dk / 2)).mean(axis=0)
        Kd[j] = Dk.mean(axis=0)
    F = np.exp(1j * np.outer(pk_, xs_)) / Lk
    nc = np.linalg.norm(np.einsum("jx,jab->xab", F, Kc), axis=(1, 2))
    nd = np.linalg.norm(np.einsum("jx,jab->xab", F, Kd), axis=(1, 2))
    sel = np.array([16, 32, 64, 128])
    tails[m0] = dict(Dc_exponent=float(np.polyfit(np.log(sel), np.log(nc[sel]), 1)[0]),
                     D_rate=float(np.polyfit(np.arange(4, 20), np.log(nd[4:20]), 1)[0]),
                     Dc_bei_64=float(nc[64]), D_bei_16=float(nd[16]))
say("   5c Theta_W wirkt ortslokal auf chi = (1 - D/2) psi; dort lautet die Wirkung psi-quer (D_c + m e^{i theta g5}) chi mit")
say("      D_c = D (1 - D/2)^-1, {g5, D_c} = 0 exakt: eine ortslokale U(1)_A. D_c hat Pole an den Verdopplerecken (D = 2):")
for m0, v in tails.items():
    say(f"      m0 = {m0}: |D_c(x)| ~ x^{v['Dc_exponent']:.3f} (x = 16 .. 128), |D(x)| ~ exp({v['D_rate']:.3f} x)")
dc_anti = []
for p in ((0.4, 1.1), (2.9, 0.2), (1.7, 3.0)):
    _, Dk = ops_k(2, np.exp(1j * p[0]), np.array([p[1]]))
    Dc = Dk[0] @ np.linalg.inv(DIM[2]["I"] - Dk[0] / 2)
    dc_anti.append(float(np.abs(DIM[2]["g5"] @ Dc + Dc @ DIM[2]["g5"]).max() + np.abs(Dc + Dc.conj().T).max()))
check("D_c chiral exakt ({g5, D_c} = 0, antihermitesch) und potenzartig nichtlokal (Exponent -1 auf 0,02 in 2D)",
      max(dc_anti) < 1e-12 and all(abs(v["Dc_exponent"] + 1) < 0.02 and v["D_rate"] < -0.3 for v in tails.values()))
say("      Damit umgeht Theta_W den Z24-Satz (exakte ortslokale U(1) => Index 0) nur über diese Nichtlokalität:")
say("      An den Verdopplerecken ist D_c unendlich, die Masse fällt dort heraus, nur die physische Ecke windet.")
befund("Positive Transfer mit flacher Phase und geschützter Naht ohne Abstimmung (Theta_W, nichtultralokal)", True)
befund("Theta_W ist in den Variablen, auf denen es ortslokal wirkt, mit lokaler Wirkung verträglich", False)
ZAHLEN["einordnung"] = {"V_spannweite_GJ": sp_gj, "V_spannweite_Wilson_twisted": sp_wi,
                        "index_GJ_wirbel": {str(k): v for k, v in idx_gj.items()}, "D_c_schwanz": {str(k): v for k, v in tails.items()}}

# =========================================================================================== Teil 6
say("\nTeil 6: Gilt der vierte Ausweg für H6s G-Paarung selbst?")
say("   G: psi, chi Weylfelder gleicher Chiralität, Paarung Phi psi chi, J: psi -> chi -> -psi (J^2 = -1), Phi -> -Phi.")
say("   Diraclesart Psi = (psi, chi^c): psi = P^_- Psi (Lüscher, mit Hut), chi = Psi-quer P_- (ohne Hut), Phi psi chi + h.c.")
say("   = Psi-quer m e^{i theta g5} (1 - D/2) Psi ist die GJ-Paarung aus Teil 3. Gitter-J: Psi -> g5 C Psi-quer^T,")
say("   Psi-quer -> Psi^T C^-1 g5 mit C = g1 g3 (Z19). Auf K in Psi-quer K Psi wirkt J als K -> -g5 (C K^T C^-1) g5.")
C4 = gam4[0] @ gam4[2]
C4i = np.linalg.inv(C4)
g54, gt4 = DIM[4]["g5"], DIM[4]["gt"]
Z4 = np.zeros((4, 4))
cprop = max([float(np.abs(C4 @ g @ C4i + g.T).max()) for g in gam4]
            + [float(np.abs(C4.T + C4).max()), float(np.abs(C4 @ g54.T @ C4i - g54).max())])
Jm = np.block([[Z4, g54 @ C4], [g54.T @ np.linalg.inv(C4.T), Z4]])
J2 = float(np.abs(Jm @ Jm + np.eye(8)).max())
Um = DIM[4]["E-"]
Vimg = C4.T @ g54.T @ Um.conj()                                  # J(u^+ psi) = psi-quer v, v = C^T g5^T u*
jmap = float(np.abs((DIM[4]["I"] - g54) / 2 @ Vimg - Vimg).max())
al = 0.7
Ra = rot(4, al)
Ua = np.block([[Ra, Z4], [Z4, Ra.T]])                             # U(1)_A: Psi -> e^{i a g5} Psi, Psi-quer -> Psi-quer e^{i a g5}
jcomm = float(np.abs(Jm @ Ua - Ua @ Jm).max())
A_JT = -gt4 @ C4.T @ g54.T                                       # erst J, dann Theta_W (auf u^+ chi)
A_TJ = -np.conj(C4i @ g54 @ gt4)                                 # erst Theta_W, dann J
tjcomm = float(np.abs(A_JT - A_TJ).max())
say(f"   6a Algebra: C-Eigenschaften Rest {cprop:.0e}; J^2 + 1 = {J2:.0e}; J bildet die psi-Komponenten (P_- Psi) auf")
say(f"      psi-quer P_- = chi ab (Rest {jmap:.0e}); [J, U(1)_A] = {jcomm:.0e}; Theta_W J = J Theta_W (Rest {tjcomm:.0e}).")
check("6a: J^2 = -1 = (-1)^F, J: psi -> chi, [J, U(1)_A] = 0 und [Theta_W, J] = 0 exakt",
      cprop < 1e-14 and J2 < 1e-14 and jmap < 1e-14 and jcomm < 1e-14 and tjcomm < 1e-14)
# 6b Lüschers Originalvariablen, 4D-Ortsraum 4 x 4^3 ohne Holonomie (J kehrt die Vektorladung von Psi um)
DW, co = build_DW(4, 4, 4, phi=0.0)
D, G5, _ = overlap(DW, 4)
n = D.shape[0]
Ib = np.eye(n)
Cb, Cbi = np.kron(np.eye(n // 4), C4), np.kron(np.eye(n // 4), C4i)
Phm = (Ib - G5 @ (Ib - D)) / 2
PMb = (Ib - G5) / 2
Qimg = G5 @ Cb @ Phm.T @ Cbi @ G5
q_hat = float(np.abs(Qimg - Phm).max())
q_dev = float(np.linalg.norm(Qimg - PMb, 2))
d_half = float(np.linalg.norm(D, 2) / 2)
KC = lambda K: Cb @ K.T @ Cbi
Eth = np.kron(np.eye(n // 4), rot(4, 1.1))
Kd = lambda mm: D + mm * Eth @ (Ib - D / 2)
Dc = D @ np.linalg.inv(Ib - D / 2)
Kw = lambda mm: Dc + mm * Eth
res_kin = float(np.abs(-G5 @ KC(D) @ G5 + D.conj().T).max())
res_DDd = float(np.abs(D + D.conj().T - D.conj().T @ D).max())
res_Kd = float(np.abs(Kd(0.3) + G5 @ KC(Kd(-0.3)) @ G5).max() / np.abs(Kd(0.3)).max())
res_Kw = float(np.abs(Kw(0.3) + G5 @ KC(Kw(-0.3)) @ G5).max() / np.abs(Kw(0.3)).max())
Sd = lambda mm: np.linalg.inv(Kd(mm))
Gw_ = lambda mm: (Ib - D / 2) @ Sd(mm)
cov = lambda G: float(np.abs(-G5 @ Cb @ G(0.3).T @ Cbi @ G5 - G(-0.3)).max() / np.abs(G(-0.3)).max())
cov_d, cov_w = cov(Sd), cov(Gw_)
PPb = (Ib + G5) / 2
q_anchor = float(np.abs(G5 @ Cb @ PMb.T @ Cbi @ G5 - PMb).max())
say("   6b Lüschers Originalvariablen (4D-Ortsraum 4 x 4^3, m = 0,3, theta = 1,1):")
say(f"      J bildet den Hut-Raum P^_- Psi auf psi-quer Q mit Q = P^_- (rechtsseitiger Hut, Rest {q_hat:.0e}) ab, nicht auf")
say(f"      psi-quer P_-: ||Q - P_-|| = {q_dev:.4f} = ||D||/2 = {d_half:.4f}. Wirkung: D -> -D^+ (Rest {res_kin:.0e}), der Bruch ist")
say(f"      D + D^+ = D^+D (Rest {res_DDd:.0e}), der hermitesche Wilson-artige Teil des Overlap ist J-ungerade. Relativer Rest")
say(f"      der Wirkung mit Paarung {res_Kd:.2f}, des Korrelators {cov_d:.2f} (Gitterskala, Verdoppler).")
pj = {}
for lab, (pt_, ps_) in (("(0,05;0,02,0,0)", (0.05, (0.02, 0, 0))), ("(0,2;0,1,0,0)", (0.2, (0.1, 0, 0))),
                        ("(1;0,5,0,2,0)", (1.0, (0.5, 0.2, 0))), ("(pi-0,05;0,03,0,0)", (np.pi - 0.05, (0.03, 0, 0)))):
    ps_ = np.array(ps_, float)
    _, Dp = ops_k(4, np.exp(1j * pt_), ps_)
    _, Dmp = ops_k(4, np.exp(-1j * pt_), -ps_)
    Dp, Dmp = Dp[0], Dmp[0]
    E1 = rot(4, 1.1)
    Kp = Dp + 0.3 * E1 @ (DIM[4]["I"] - Dp / 2)
    Kmc = C4 @ (Dmp - 0.3 * E1 @ (DIM[4]["I"] - Dmp / 2)).T @ C4i
    Dcp = Dp @ np.linalg.inv(DIM[4]["I"] - Dp / 2)
    Dcm = Dmp @ np.linalg.inv(DIM[4]["I"] - Dmp / 2)
    Kwp, Kwmc = Dcp + 0.3 * E1, C4 @ (Dcm - 0.3 * E1).T @ C4i
    absp = float(np.sqrt(pt_ ** 2 + np.sum(ps_ ** 2)))
    pj[lab] = dict(abs_p=absp, D_halbe=float(np.linalg.norm(Dp, 2) / 2),
                   J_D_rest=float(np.linalg.norm(Kp + g54 @ Kmc @ g54, 2) / np.linalg.norm(Kp, 2)),
                   J_W_rest=float(np.linalg.norm(Kwp + g54 @ Kwmc @ g54, 2) / np.linalg.norm(Kwp, 2)))
say("      Impulsaufgelöst (4D, m = 0,3): " + "; ".join(
    f"{k}: ||D||/2 = {v['D_halbe']:.3f}, J_D-Rest {v['J_D_rest']:.3f}, J_W-Rest {v['J_W_rest']:.0e}" for k, v in pj.items()))
small = [v for v in pj.values() if v["abs_p"] < 0.3]
check("6b: J in Lüschers Variablen um O(a) gebrochen: Q = rechtsseitiger Hut, ||Q - P_-|| = ||D||/2, Wirkung D -> -D^+,"
      " Korrelator nicht kovariant (negatives Ergebnis)",
      q_hat < 1e-12 and abs(q_dev - d_half) < 1e-10 and res_kin < 1e-12 and res_DDd < 1e-12 and res_Kd > 0.5 and cov_d > 0.5)
check("6b: J-Bruch bei kleinen Impulsen O(a p): ||D||/2 < |p| und J_D-Rest < 2 |p| + m/2 für |p| < 0,3",
      all(v["D_halbe"] < v["abs_p"] and v["J_D_rest"] < 2 * v["abs_p"] + 0.15 for v in small))
say("   6c Ankervariablen chi = (1 - D/2) Psi, Psi-quer: J_W: chi -> g5 C Psi-quer^T, Psi-quer -> chi^T C^-1 g5.")
say(f"      Wirkung psi-quer (D_c + m e^{{i theta g5}}) chi: J_W-Rest {res_Kw:.0e} (D_c ist antihermitesch und chiral, also J-gerade;")
say(f"      die Masse ist J-ungerade, Phi -> -Phi). Korrelator kovariant (Rest {cov_w:.0e}); P_- chi -> psi-quer P_- exakt (Rest {q_anchor:.0e}).")
say("      Die Verdoppler entfernt in diesen Variablen der Pol von D_c (J-gerade, U(1)_A-gerade), keine Wilson-Majoranamasse.")
# 6c Nichtlokalität von J_W in Psi-Variablen: Psi -> (1 - D/2)^-1 g5 C Psi-quer^T; (1 - D/2)^-1 - 1 = D_c/2
Lq = 48
pq = 2 * np.pi * (np.arange(Lq) + 0.5) / Lq
P1, P2, P3 = [a.ravel() for a in np.meshgrid(pq, pq, pq, indexing="ij")]
acc_c = np.zeros((Lq, 4, 4), complex)
acc_d = np.zeros((Lq, 4, 4), complex)
I4m = DIM[4]["I"]
for j, pt_ in enumerate(pq):
    st, ct = np.sin(pt_), np.cos(pt_)
    ssq = [np.sin(P1), np.sin(P2), np.sin(P3)]
    Wq = 4 - np.cos(P1) - np.cos(P2) - np.cos(P3) - ct
    kin = 1j * st * gam4[3][None] + sum(1j * ssq[i][:, None, None] * gam4[i][None] for i in range(3))
    fq = st ** 2 + ssq[0] ** 2 + ssq[1] ** 2 + ssq[2] ** 2 + (Wq - 1.0) ** 2
    Dq4 = I4m[None] + (kin + (Wq - 1.0)[:, None, None] * I4m[None]) / np.sqrt(fq)[:, None, None]
    acc_c[j] = np.linalg.inv(I4m[None] - Dq4 / 2).mean(axis=0) - I4m
    acc_d[j] = Dq4.mean(axis=0)
Fq = np.exp(1j * np.outer(pq, np.arange(Lq))) / Lq
ncq = np.linalg.norm(np.einsum("jx,jab->xab", Fq, acc_c), axis=(1, 2))
ndq = np.linalg.norm(np.einsum("jx,jab->xab", Fq, acc_d), axis=(1, 2))
odd = np.array([9, 11, 13, 15])
exp4 = float(np.polyfit(np.log(odd), np.log(ncq[odd]), 1)[0])
rate4 = float(np.polyfit(np.arange(6, 12), np.log(ndq[6:12]), 1)[0])
say(f"      In Psi-Variablen ist J_W nichtlokal: Kern von (1 - D/2)^-1 - 1 = D_c/2 entlang der Zeitachse (4D, 48^4, m0 = 1)")
say(f"      |K(x)| ~ x^{exp4:.2f} (x = 9 .. 15 ungerade; Erwartung -(d-1) = -3), Overlap |D(x)| ~ exp({rate4:.2f} x).")
check("6c: J_W in Ankervariablen exakt (Wirkung, Korrelator, Projektoren; Reste < 1e-12, impulsaufgelöst < 1e-12)",
      res_Kw < 1e-12 and cov_w < 1e-12 and q_anchor < 1e-14 and all(v["J_W_rest"] < 1e-12 for v in pj.values()))
check("6c: J_W in Psi-Variablen potenzartig nichtlokal (4D-Exponent zwischen -3,6 und -2,8), Overlap exponentiell lokal",
      -3.6 < exp4 < -2.8 and rate4 < -0.5)
say("   6d Majorana-Lesart: psi und chi je als Majoranafeld mit eigenem Overlap. Dann ist J eine SO(2)-Drehung zweier")
say("      gleicher Felder, lokal und exakt; U(1)_A gibt es dort nur als e^{i a X}, X = g5 (1 - D/2) nicht quantisiert")
say(f"      (|X| ab {absX.min():.3f}, Teil 3), also kein Theta_W, und Z19s Verletzung bei gedrehter Phase bleibt. Die Ankerform")
say("      zeigt aber: Eine exakte J erzwingt X nicht. In (chi, psi-quer) sind J (J^2 = -1), U(1)_A (e^{i a g5}, quantisiert)")
say("      und Theta_W zugleich exakt und ortslokal; Theta_W-Positivität (Teil 3, 4) gilt für denselben Zustand.")
befund("Gitter-J in Lüschers Originalvariablen exakt", False)
befund("Exakte Gitter-J mit J^2 = -1, verträglich mit Theta_W und U(1)_A (Ankerform)", True)
befund("Exakte J erzwingt die Majoranaform X", False)
befund("Vierter Ausweg gilt für H6s G selbst (J exakt, gleicher Preis: nichtlokale Wirkung in den Ankervariablen)", True)
ZAHLEN["teil6_J"] = {"J2_rest": J2, "Q_minus_P-": q_dev, "D_halbe": d_half, "J_D_rest_wirkung": res_Kd,
                     "J_D_rest_korrelator": cov_d, "J_W_rest_wirkung": res_Kw, "J_W_rest_korrelator": cov_w,
                     "impulsaufgeloest": pj, "J_W_kern_exponent_4d": exp4, "D_rate_4d": rate4,
                     "Theta_W_J_kommutator": tjcomm}

# =========================================================================================== Urteil
URTEIL = (
    "Fall (a), mit genau benanntem Preis. Für freie Lüscher-Weylfelder (Overlap, abelsch, a = 1, 2D und 4D) gilt Link-RP: "
    f"Alle linearen OS-Formen (psi-quer P_+, Anker P_- psi) sind >= 0 für m0 = 0,5 bis 1,9, L_t = 8 bis 4096, L_s = 8 bis 64, "
    f"Holonomien und Abstände 1 bis 3 zur Ebene (kein signifikant negativer Eigenwert; negative Werte nur in Rundung, roh "
    f"<= {RAUSCH['max_roh_durch_rausch']:.1f} x 1e-15 max|G|), auf der ganzen Halbzeit im "
    "Ortsraum bis 24 x 24 auch mit nicht konstantem Theta-symmetrischem Hintergrund, und die Spektralgewichte (Pol plus "
    "Schnitt, L_t = unendlich) sind positiv-semidefinit. Theta bildet Lüschers Raum (P^_-, P_+) aber nicht auf sich ab, "
    "sondern auf (P_-, gespiegelter Hut), Abweichung ||D||/2 (O(ap), an den Verdopplern 1). Exakt ist S_W = P_+/2 delta "
    "+ P_- D^-1 P_+: Die ganze Abweichung ist Kontakt. Mit der vollen Komponentenalgebra verletzt der Kontakt die "
    "Reflexionssymmetrie der Komposita psi-quer P_+ psi um genau 1/2. Die Ankeralgebra {psi-quer P_+, P_-(1 - D/2) psi} ist "
    "abgeschlossen und RP. Bei gedrehter Phase (Lüscher-Paarung L-R = GW-Drehung, die G-Paarung psi chi aus Z24) verletzt "
    f"die Reflexion auf psi (Theta_D, wie Z19) RP: nicht hermitesch, {ah[1]:.2f} bei m = 0,3, linear in a m, stabil in L. "
    "Die Ursache ist exakt lokal, M_D = (1 - (m/2) e^{-i theta g5}) M_W. Die Reflexion auf chi = (1 - D/2) psi (Theta_W) ist "
    "bei jedem theta positiv, weil die GW-Drehung auf chi ortslokal als e^{i theta g5} wirkt. Mit flacher Phase (Spannweite "
    f"{sp_gj:.0e}) und Index = Windung gibt das den vierten Ausweg aus dem Z24-Trilemma: geschützt, flach, positiv, ohne "
    "Abstimmung, aber nicht ultralokal. Die Variablen chi tragen eine potenzartig nichtlokale Wirkung (D_c ~ 1/x in 2D, Pole "
    "an den Verdopplern), und der Transfer hat ein Kontinuum je Impuls (Z19). Für die flachen Majoranaformen X und Xq aus Z19 "
    "gibt es kein solches Theta_W (X nicht quantisiert); dort bleibt die Verletzung. Für H6s G selbst (psi, chi gleicher "
    "Chiralität, J: psi -> chi -> -psi) gilt der vierte Ausweg ebenfalls: In Lüschers Originalvariablen ist J um O(a) gebrochen "
    f"(||Q - P_-|| = ||D||/2 = {d_half:.3f}, Wirkung D -> -D^+), in den Ankervariablen (chi, psi-quer) sind J mit J^2 = -1, "
    "U(1)_A und Theta_W zugleich exakt und ortslokal, zum selben Preis (dort nichtlokale Wirkung; J in psi-Variablen "
    f"~ x^{exp4:.1f} in 4D); die Majorana-Lesart mit lokaler J erzwingt dagegen X ohne Theta_W. Grenzen: frei, fester Hintergrund, keine "
    "Eichfeldintegration und keine komplexe Weyl-Determinante; 4D nur Impulsraum und Gitter bis 8 x 4^3."
)
say("\nUrteil")
say("   " + URTEIL)
say(f"\nPrüfungen: {sum(CHECKS.values())}/{len(CHECKS)} bestanden. Befunde: "
    + ", ".join(f"{k} = {'ja' if v else 'nein'}" for k, v in BEFUNDE.items()))
say(f"Laufzeit {time.time() - T0:.0f} s.")
(HERE / "z26_ergebnis.json").write_text(json.dumps({"checks": CHECKS, "befunde": BEFUNDE, "urteil": URTEIL,
                                                    "zahlen": ZAHLEN}, indent=1, ensure_ascii=False, default=float))
assert all(CHECKS.values()), [k for k, v in CHECKS.items() if not v]
