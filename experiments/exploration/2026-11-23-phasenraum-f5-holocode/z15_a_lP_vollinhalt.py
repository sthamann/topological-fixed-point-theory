"""Z15: a/l_P mit dem vollen Feldinhalt von H6 (Folgeauftrag 5 zu Prämisse D, Anschluss an Z2).

Unter D gilt 8 pi G = 2 pi/eta mit eta = kappa_gesamt/a^2, also a/l_P = 2 sqrt(kappa_gesamt). Z2 hat nur
Fermionen gerechnet. Hier kommen hinzu:
 1. freies reelles Skalarfeld auf Z^3 (masselos und massiv), Dimensionsreduktion auf harmonische Ketten
    je Querimpuls, Infrarot über kleine Massen und Ringgröße kontrolliert;
 2. das U(1)-Gitterphoton: je Querimpuls exakt eine TE-Kette (Plätze) und eine TM-Kette (x-Links,
    phi = E_x/M, pi = Plakette b). Zentrumswahl: ohne Zentrum, elektrisches Zentrum (destillierbar)
    und der Shannon-Anteil des ganzzahligen Randflusses (Kante);
 3. Wilson- und naive Fermionen, massiv und masselos, samt Reproduktion von Z2;
 4. Szenarien (SM mit nu_R, schweres reelles 16-Plett, 12 Eichbosonen, 4 Higgs) mit Band und Empfindlichkeit.
"""
import json
import time
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.linalg import toeplitz
from scipy.special import ellipk

HERE = Path(__file__).parent
LINES, CHECKS, BEFUNDE = [], {}, {}
T0 = time.time()


def say(s=""):
    print(s)
    LINES.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    return "ja" if cond else "NEIN"


def befund(name, val):
    BEFUNDE[name] = bool(val)
    return "ja" if val else "nein"


# ---------------------------------------------------------------------------------------------------
# Bosonische Ketten
# ---------------------------------------------------------------------------------------------------
def ent_bos(nu):
    nu = nu[nu > 0.5 + 1e-12]
    return float(np.sum((nu + 0.5) * np.log(nu + 0.5) - (nu - 0.5) * np.log(nu - 0.5)))


def chain_xp(M, ring):
    """Grundzustandskorrelatoren x_d = <phi_0 phi_d>, p_d = <pi_0 pi_d> der Kette mit Masse M (antiperiodischer Ring)."""
    q = 2 * np.pi * (np.arange(ring) + 0.5) / ring
    w = np.sqrt(M * M + 4 * np.sin(q / 2) ** 2)
    cosm = np.cos(np.outer(np.arange(ring // 2), q))
    return cosm @ (1 / (2 * w)) / ring, cosm @ (w / 2) / ring


def S_boson(M, ring, center=False):
    """Entropie je Schnitt einer halben Ringkette. center=True: die beiden Randplätze tragen nur phi
    (elektrisches Zentrum der TM-Kette); zurück kommt der destillierbare Teil S(X_bedingt, P)."""
    x, p = chain_xp(M, ring)
    L = ring // 2
    X, P = toeplitz(x), toeplitz(p)
    if center:
        c, inner = [0, L - 1], np.arange(1, L - 1)
        Xic = X[np.ix_(inner, c)]
        X = X[np.ix_(inner, inner)] - Xic @ np.linalg.solve(X[np.ix_(c, c)], Xic.T)
        P = P[np.ix_(inner, inner)]
    Lc = np.linalg.cholesky(X)
    ev = np.linalg.eigvalsh(Lc.T @ P @ Lc)
    return ent_bos(np.sqrt(np.clip(ev, 0.25, None))) / 2


M_LO, M_TOP = 0.01, 3.05


def boson_table(base=400, center=False):
    Ms = np.concatenate([np.geomspace(M_LO, 0.1, 30, endpoint=False), np.linspace(0.1, M_TOP, 150)])
    S = np.array([S_boson(M, max(base, 2 * int(base / 800 * 10 / M)), center) for M in Ms])
    return Ms, S


def interp(table, ir="steigung"):
    Ms, S = table
    cs = CubicSpline(np.log(Ms), S)
    slope = (S[1] - S[0]) / (np.log(Ms[1]) - np.log(Ms[0]))

    def f(M):
        M = np.asarray(M, float)
        out = np.empty_like(M)
        lo = M < Ms[0]
        out[~lo] = cs(np.log(M[~lo]))
        out[lo] = S[0] + (slope * (np.log(M[lo]) - np.log(Ms[0])) if ir == "steigung" else 0.0)
        return out
    return f, slope


def kgrid(nk):
    k = 2 * np.pi * (np.arange(nk) + 0.5) / nk
    K1, K2 = np.meshgrid(k, k)
    return 4 * np.sin(K1 / 2) ** 2 + 4 * np.sin(K2 / 2) ** 2          # M^2 ohne Masse


def kappa_boson(f, m=0.0, nk=480):
    return float(f(np.sqrt(kgrid(nk) + m * m)).mean())


# ---------------------------------------------------------------------------------------------------
# Fermionische Ketten: h(q) = sx sin q + sy s + sz (B + r (1 - cos q)), zwei Kopien je Dirac-Feld
# ---------------------------------------------------------------------------------------------------
def ent_f(nu):
    nu = nu[(nu > 1e-14) & (nu < 1 - 1e-14)]
    return float(-np.sum(nu * np.log(nu) + (1 - nu) * np.log(1 - nu)))


def S_fermion(s, B, r, ring=200):
    q = 2 * np.pi * (np.arange(ring) + 0.5) / ring
    dx, dy, dz = np.sin(q), np.full(ring, s), B + r * (1 - np.cos(q))
    E = np.sqrt(dx ** 2 + dy ** 2 + dz ** 2)
    blocks = {(0, 0): 0.5 * (1 - dz / E), (1, 1): 0.5 * (1 + dz / E),
              (0, 1): -0.5 * (dx - 1j * dy) / E, (1, 0): -0.5 * (dx + 1j * dy) / E}
    half = ring // 2
    ph = np.exp(1j * np.outer(np.arange(-half + 1, half), q)) / ring
    idx = np.arange(half)[:, None] - np.arange(half)[None, :] + half - 1
    C = np.zeros((half, 2, half, 2), complex)
    for (a, b), v in blocks.items():
        C[:, a, :, b] = (ph @ v)[idx]
    return ent_f(np.linalg.eigvalsh(C.reshape(2 * half, 2 * half))) / 2


def kappa_dirac(m, r, nk=128, ring=200):
    """Direkte 2D-Integration über Querimpulse (Symmetrie k -> -k und k1 <-> k2)."""
    k = 2 * np.pi * (np.arange(nk // 2) + 0.5) / nk
    tot = wsum = 0.0
    for a in range(nk // 2):
        for b in range(a, nk // 2):
            w = 1.0 if a == b else 2.0
            s = np.hypot(np.sin(k[a]), np.sin(k[b]))
            B = m + r * (2 - np.cos(k[a]) - np.cos(k[b]))
            tot += w * S_fermion(s, B, r, ring)
            wsum += w
    return 2 * tot / wsum


def kappa_naive_table(masses, ring=200, nk=480):
    """Z2-Weg: Tabelle in mu, mu^2 = sin^2 k1 + sin^2 k2 + m^2."""
    mus = np.unique(np.concatenate([np.geomspace(1e-6, 0.05, 300), np.linspace(0.05, np.sqrt(3.0), 500)]))
    tab = np.array([S_fermion(mu, 0.0, 0.0, ring) for mu in mus])
    k = 2 * np.pi * (np.arange(nk) + 0.5) / nk
    K1, K2 = np.meshgrid(k, k)
    out = {}
    for m in masses:
        mu = np.sqrt(np.sin(K1) ** 2 + np.sin(K2) ** 2 + m * m)
        out[m] = 2 * float(np.interp(np.log(np.maximum(mu, 1e-6)).ravel(), np.log(mus), tab).mean())
    return out


def H_dg(sig, d=0.0):
    """Shannon-Entropie der diskreten Gaußverteilung p_n ~ exp(-(n-d)^2/(2 sig^2)) auf Z."""
    K = int(np.ceil(12 * sig)) + 5
    n = np.arange(-K, K + 1)
    w = np.exp(-(n - d) ** 2 / (2 * sig * sig))
    p = w / w.sum()
    p = p[p > 0]
    return float(-np.sum(p * np.log(p)))


def H_dg_mittel(sig):
    return float(np.mean([H_dg(sig, d) for d in np.linspace(0, 1, 64, endpoint=False)]))


# ===================================================================================================
say("Z15: Gitterabstand gegen Plancklänge mit dem vollen Feldinhalt von H6")
say("Prämisse D: a/l_P = 2 sqrt(kappa_gesamt), kappa in nat je Gitterfläche, ebener Schnitt in Z^3.")
say("Methode: je Querimpuls k zerfällt jedes freie Feld in 1D-Ketten. kappa = Mittel über k der Kettenentropie.\n")

# --------------------------------------------------------------------------- Teil 1: Skalar
say("Teil 1: freies reelles Skalarfeld")
say("   Kette je k: H = 1/2 sum [pi^2 + (phi_{n+1} - phi_n)^2 + M^2 phi^2], M^2 = m^2 + 4 sin^2(k1/2) + 4 sin^2(k2/2).")
tabs = {base: boson_table(base) for base in (200, 400, 800)}
f400, slope = interp(tabs[400])
x_s, _ = chain_xp(1.0, 2000)
say(f"   Kette bei kleiner Masse: S = -c ln M + const mit c = {-slope:.4f} (freies Boson: 1/6 = {1 / 6:.4f}).")
say(f"   Steigung wie für c = 1: {check('Skalarkette: Steigung -1/6 bei kleiner Masse', abs(slope + 1 / 6) < 0.003)}")
say("   Konvergenz des masselosen kappa_s (Ringgröße, Querimpulse):")
conv_s = {}
for base in (200, 400, 800):
    f, _ = interp(tabs[base])
    for nk in (240, 480, 960):
        if base != 400 and nk != 480:
            continue
        conv_s[(base, nk)] = kappa_boson(f, 0.0, nk)
        say(f"   Ring {base:4d}, Querimpulse {nk:4d}^2: kappa_s = {conv_s[(base, nk)]:.7f}")
ks0 = conv_s[(400, 960)]
spread_ring = max(abs(conv_s[(b, 480)] - conv_s[(400, 480)]) for b in (200, 800)) / ks0
spread_nk = max(abs(conv_s[(400, n)] - ks0) for n in (240, 480)) / ks0
say(f"   Ring konvergiert (rel. < 1e-4): {check('Skalar: Konvergenz im Ring', spread_ring < 1e-4)} ({spread_ring:.1e});"
    f" Querimpulse: {check('Skalar: Konvergenz in Querimpulsen', spread_nk < 1e-4)} ({spread_nk:.1e})")
f_const, _ = interp(tabs[400], ir="konstant")
d_ir = abs(kappa_boson(f_const, 0.0, 960) - ks0)
say(f"   Infrarot: Fortsetzung unter M = {M_LO} mit Steigung gegen konstant ändert kappa_s um {d_ir:.1e}. "
    f"Vernachlässigbar: {check('Skalar: Infrarotfortsetzung unerheblich', d_ir < 1e-5)}")
ks_m = {m: kappa_boson(f400, m, 960) for m in (0.003, 0.01, 0.03, 0.1, 0.3, 0.5, 1.0)}
say("   Kleine Massen als Infrarotregler: " + ", ".join(f"m = {m}: {v:.6f}" for m, v in ks_m.items()))
mono = all(a > b for a, b in zip([ks0] + list(ks_m.values()), list(ks_m.values())))
say(f"   m -> 0 stetig (|kappa_s(0.003) - kappa_s(0)| < 1e-4) und monoton: "
    f"{check('Skalar: masseloser Limes stetig und monoton', abs(ks_m[0.003] - ks0) < 1e-4 and mono)}")
say(f"   Ergebnis: kappa_s = {ks0:.5f} nat je Gitterfläche je reellem Skalar (masselos).")
say(f"   Zum Vergleich Srednicki (radiales Gitter): 0,30/(4 pi) = {0.30 / (4 * np.pi):.4f}. Gleiche Größenordnung, anderes Gitter.\n")

# --------------------------------------------------------------------------- Teil 2: Photon
say("Teil 2: U(1)-Gitterphoton (Kogut-Susskind, schwache Kopplung, e = g E_int)")
say("   Je Querimpuls zerfällt die eichinvariante Algebra exakt in zwei Ketten mit derselben Masse M = |s_perp|:")
say("   TE: a = Querlinks senkrecht zu s_perp, auf Plätzen. Sie ist genau die Skalarkette.")
say("   TM: phi = E_x/M auf den x-Links, pi = Plakette b = grad_x A_par - M A_x. Gauß eliminiert E_par.")
say("   Auch TM ist eine harmonische Kette mit Masse M, nur auf den Linkpositionen n + 1/2.")


def maxwell_mode(k1, k2, N=40):
    """Volles Gittermaxwell für festen Querimpuls auf einem periodischen x-Ring: Links A_x(n), A_1(n), A_2(n),
    Plaketten B_12(n), B_x1(n+1/2), B_x2(n+1/2). Vergleich mit TE- und TM-Kette."""
    e1, e2 = np.exp(1j * k1) - 1, np.exp(1j * k2) - 1
    eye, zero = np.eye(N), np.zeros((N, N))
    D = np.roll(eye, -1, axis=0) - eye
    curl = np.block([[zero, -e2 * eye, e1 * eye], [-e1 * eye, D, zero], [-e2 * eye, zero, D]])
    grad = np.vstack([D, e1 * eye, e2 * eye])
    lam, V = np.linalg.eigh(curl.conj().T @ curl)
    phys = lam > 1e-9
    M2 = abs(e1) ** 2 + abs(e2) ** 2
    q = 2 * np.pi * np.arange(N) / N
    w2 = M2 + 4 * np.sin(q / 2) ** 2
    spec_err = float(np.max(np.abs(np.sort(lam[phys]) - np.sort(np.concatenate([w2, w2])))))
    sq = np.where(phys, np.sqrt(np.abs(lam)), 0.0)
    ee = [float(np.real(0.5 * np.sum(sq * V[0, :] * V[d, :].conj()))) for d in (0, 1, 2)]
    ee_chain = [float(M2 * np.mean(np.cos(q * d) / (2 * np.sqrt(w2)))) for d in (0, 1, 2)]
    return float(np.max(np.abs(curl @ grad))), int((~phys).sum()), spec_err, ee, ee_chain


mx_ok = True
for kk in ((0.7, 1.9), (2.5, 0.3)):
    cg, nz, se, ee, eec = maxwell_mode(*kk)
    rel = max(abs(a - b) / abs(b) for a, b in zip(ee, eec))
    mx_ok &= cg < 1e-12 and nz == 40 and se < 1e-10 and rel < 1e-10
    say(f"   Volles Gittermaxwell bei k = {kk}: rot grad = {cg:.0e}, {nz} Eichmoden, Spektrum = 2 Ketten bis {se:.0e},")
    say(f"   <E_x(0) E_x(d)> d = 0,1,2: " + ", ".join(f"{a:.8f}/{b:.8f}" for a, b in zip(ee, eec)) + " (Maxwell/TM-Kette)")
say(f"   TE + TM ist exakt das Gitterphoton, phi = E_x/M stimmt: "
    f"{check('Photon: volles Gittermaxwell = TE-Kette + TM-Kette (Spektrum und Flusskorrelator)', mx_ok)}")
x00_ring = x_s[0]
x00_ell = ellipk(4 / 5) / (np.pi * np.sqrt(5))
say(f"   X_00(M=1) Ring gegen elliptisches Integral: {x00_ring:.12f} / {x00_ell:.12f}: "
    f"{check('Photon: X_00 Ring gleich K(m)/(pi sqrt(M^2+4))', abs(x00_ring - x00_ell) < 1e-10)}")

tabs_tm = {base: boson_table(base, center=True) for base in (200, 400, 800)}
ftm, _ = interp(tabs_tm[400])
ktm = {base: kappa_boson(interp(tabs_tm[base])[0], 0.0, 480) for base in (200, 400, 800)}
ktm0 = kappa_boson(ftm, 0.0, 960)
spread_tm = max(abs(ktm[b] - ktm[400]) for b in (200, 800)) / ktm[400]
say(f"   Elektrisches Zentrum: Die x-Links über dem Schnitt tragen nur E_x. Bedingt auf diesen Fluss bleibt")
say(f"   der destillierbare TM-Anteil kappa_TM = {ktm0:.6f} (Ring 200/400/800: "
    + "/".join(f"{ktm[b]:.6f}" for b in (200, 400, 800)) + ").")
say(f"   Konvergiert: {check('Photon: destillierbarer TM-Anteil konvergiert im Ring', spread_tm < 1e-3)}")
say(f"   Bedingen senkt die Entropie bei jedem M: "
    f"{check('Photon: S_dest <= S_Kette für alle M', bool(np.all(tabs_tm[400][1] <= tabs[400][1] + 1e-12)))}")
# Randfluss (Kante)
M2f = kgrid(960)
C0 = M2f * ellipk(4 / (M2f + 4)) / (np.pi * np.sqrt(M2f + 4))
C0_marg, C0_geo = float(C0.mean()), float(np.exp(np.log(C0).mean()))
h_e = 0.5 * np.log(2 * np.pi * np.e) + 0.5 * np.log(C0_geo)
n3 = 64
p3 = 2 * np.pi * (np.arange(n3) + 0.5) / n3
w3 = np.sqrt(sum(4 * np.sin(P / 2) ** 2 for P in np.meshgrid(p3, p3, p3, indexing="ij")))
say(f"   Randfluss je Link in e-Einheiten: Varianz {C0_marg:.5f} (3D: <omega>/3 = {w3.mean() / 3:.5f}),"
    f" Innovationsvarianz {C0_geo:.5f} (Szegő).")
say(f"   Differentielle Entropie h_e = 1/2 ln(2 pi e) + 1/2 <ln C0(k)> = {h_e:.5f} nat je Link.")
say("   Ganzzahliger Fluss E_int = e/g: Shannon-Anteil s(g) ~ h_e - ln g (wächst ohne Schranke für g -> 0).")
gauss_ok = abs(H_dg_mittel(4.0) - 0.5 * np.log(2 * np.pi * np.e * 16)) < 1e-6
say(f"   Diskrete Entropie geht für breite Verteilung in die Gaußformel über: "
    f"{check('Kante: diskrete Entropie -> Gaußgrenze (sigma = 4)', gauss_ok)}")
G_LIST = (0.5, 0.55, 0.75, 1.0, 1.5)
edge = {}
for g in G_LIST:
    edge[g] = {"innovation": H_dg_mittel(np.sqrt(C0_geo) / g), "marginal": H_dg(np.sqrt(C0_marg) / g),
               "gauss": float(h_e - np.log(g))}
say("   Kante je Eichboson je Gitterfläche (Innovationsnäherung / Einzellink als obere Schranke):")
for g in G_LIST:
    say(f"   g = {g:4.2f}: s = {edge[g]['innovation']:.4f} / {edge[g]['marginal']:.4f}")
gam_nc = 2 * ks0
gam_ez = ks0 + ktm0
say(f"   Photon ohne Zentrum (Randlink ganz einer Seite oder keiner): 2 kappa_s = {gam_nc:.5f}.")
say(f"   Photon, elektrisches Zentrum, destillierbar: kappa_s + kappa_TM = {gam_ez:.5f}.")
say(f"   Photon, elektrisches Zentrum mit Kante bei g = 0,55: {gam_ez + edge[0.55]['innovation']:.4f}.")
say(f"   Die Kante ist {edge[0.55]['innovation'] / gam_ez:.0f} mal größer als der destillierbare Teil.")
say("   Mehrdeutigkeit: Die Wahl des Zentrums ist Teil der Definition. Im Wärmeleitungs- oder Kegelschema")
say("   (Kabat) erscheint derselbe Kantenbeitrag als negativer Kontaktterm (Vektor: 2 - 6 = -4 Skalareinheiten);")
say("   Donnelly-Wall zeigen, dass dieser Term die Randmoden sind. Auf dem kompakten Gitter ist er positiv.")
say("   Nichtabelsch: hier als Kopien von U(1) gezählt (führende Ordnung). In der Formel H(R) + <ln dim R>")
say("   hat ln(1/g) den Koeffizienten (dim G + rang G)/2 statt dim G: SU(3) 5 statt 8, SU(2) 2 statt 3.")
say("   Das gibt eine untere Variante mit 8 statt 12 effektiven Kopien. U(1)_Y mit Ladungsquant 1/6 hätte")
say(f"   zusätzlich ln 6 = {np.log(6):.3f} je Link.\n")

# --------------------------------------------------------------------------- Teil 3: Fermionen
say("Teil 3: Fermionen")
s1 = S_fermion(1.0, 0.0, 0.0)
s1r = S_fermion(0.6, 0.8, 0.0)
say(f"   Kette bei mu = 1: {s1:.10f}; gedreht (s, B) = (0,6; 0,8): {s1r:.10f}; ln2/2 = {np.log(2) / 2:.10f}")
say(f"   Kettencode stimmt: {check('Fermionkette: ln2/2 bei mu = 1 und Drehinvarianz', abs(s1 - np.log(2) / 2) < 1e-9 and abs(s1r - s1) < 1e-9)}")
z2 = json.loads((HERE / "z2_ergebnis.json").read_text())
kn_tab = kappa_naive_table((0.0, 0.3, 0.5, 1.0))
say(f"   Naiv über die mu-Tabelle (Weg von Z2): kappa = {kn_tab[0.0]:.6f}; Z2: {z2['kappa_naive_dirac']:.6f}")
say(f"   Z2 reproduziert: {check('Fermion: Z2-Wert reproduziert (Tabellenweg)', abs(kn_tab[0.0] - z2['kappa_naive_dirac']) < 1e-6)}")
kn_dir = kappa_dirac(0.0, 0.0, 128)
say(f"   Unabhängig über direkte 2D-Integration (128^2): {kn_dir:.6f}: "
    f"{check('Fermion: direkte 2D-Integration trifft Z2 (< 2e-4)', abs(kn_dir - z2['kappa_naive_dirac']) < 2e-4)}")
kn03_dir = kappa_dirac(0.3, 0.0, 64)
say(f"   Masse 0,3 direkt / Tabelle: {kn03_dir:.6f} / {kn_tab[0.3]:.6f}: "
    f"{check('Fermion: massiv naiv, beide Wege gleich (< 2e-4)', abs(kn03_dir - kn_tab[0.3]) < 2e-4)}")
say("   Gestaffelte Fermionen: Die Spindiagonalisierung ist platzlokal und unitär. Daher ist ein naives Feld")
say(f"   exakt vier gestaffelte Felder, je zwei Geschmäcker: kappa je Geschmack = kappa_naiv/8 = {kn_tab[0.0] / 8:.6f}.")
say("   Wilson: h = sum alpha_j sin p_j + beta (m + r sum (1 - cos p_j)); je k zwei Kopien einer Kette mit d = (sin q, s, B + r(1-cos q)).")
t1 = time.time()
kw = {0.0: kappa_dirac(0.0, 1.0, 128)}
kw64 = kappa_dirac(0.0, 1.0, 64)
kw64_r400 = kappa_dirac(0.0, 1.0, 64, ring=400)
for m in (0.3, 0.5, 1.0):
    kw[m] = kappa_dirac(m, 1.0, 64)
kw_r = {0.5: {0.0: kappa_dirac(0.0, 0.5, 64), 0.5: kappa_dirac(0.5, 0.5, 64)},
        2.0: {0.0: kappa_dirac(0.0, 2.0, 64), 0.5: kappa_dirac(0.5, 2.0, 64)}}
say(f"   Wilson r = 1, masselos: {kw[0.0]:.6f} (Querimpulse 64^2: {kw64:.6f}, Ring 400: {kw64_r400:.6f})")
conv_w = abs(kw[0.0] - kw64) / kw[0.0] < 1e-3 and abs(kw64_r400 - kw64) / kw64 < 1e-4
say(f"   Konvergiert: {check('Wilson: Konvergenz in Querimpulsen und Ring', conv_w)}")
say("   Massiv (je Dirac-Feld): " + ", ".join(f"m = {m}: naiv {kn_tab[m]:.4f}, Wilson {kw[m]:.4f}" for m in (0.3, 0.5, 1.0)))
mono_f = all(kn_tab[a] > kn_tab[b] and kw[a] > kw[b] for a, b in ((0.0, 0.3), (0.3, 0.5), (0.5, 1.0)))
say(f"   Monoton fallend mit der Masse: {check('Fermion: kappa fällt monoton mit der Masse', mono_f)}")
say(f"   Wilson-Parameter, masselos: r = 0,5: {kw_r[0.5][0.0]:.4f}; r = 1: {kw[0.0]:.4f}; r = 2: {kw_r[2.0][0.0]:.4f}")
say(f"   Ein Wilson-Geschmack trägt {kw[0.0] / (kn_tab[0.0] / 8):.2f} mal so viel wie ein gestaffelter Geschmack.")
say("   Die Verdoppler werden schwer, verschwinden aber nicht aus der Flächenentropie.")
say(f"   (Rechenzeit Wilson {time.time() - t1:.0f} s)\n")

# --------------------------------------------------------------------------- Teil 4: Szenarien
say("Teil 4: Feldinhalt und Szenarien")
say("   Leicht: 3 Generationen mit nu_R = 48 Weyl = 24 Dirac (masselos an der Gitterskala).")
say("   Schwer: reelles 16-Plett der Naht = 16 Majorana-Weyl = 8 Dirac-Äquivalente mit Masse m (0,3 bis 1),")
say("   dazu das komplexe Paarungsfeld Phi als 2 reelle Skalare mit derselben Masse (Prämisse G).")
say("   Eichfelder: 12 (8 + 3 + 1), je eine Photonkopie. Higgs: 4 reelle Skalare, masselos.")
say("   Die 16 MW-Moden der Naht leben auf einer Linie im Raum. Sie tragen je Durchstoßpunkt durch die Ebene")
say("   höchstens (c_L + c_R)/12 ln(xi) = (2/3) ln(xi), also O(1) je Punkt und nichts je Fläche. Nur ein Nahtgas")
say("   auf Gitterskala (Fluktuationsweg aus Z9) könnte je Fläche beitragen; das ist hier nicht gerechnet.")
say("   Gravitonen: keine eigenen Gitterfreiheitsgrade, Gravitation ist unter D Zustandsgleichung.")


def ferm_light(scheme):
    return {"gestaffelt": kn_tab[0.0] / 8, "wilson": kw[0.0], "wilson_r0.5": kw_r[0.5][0.0],
            "wilson_r2": kw_r[2.0][0.0], "naiv": kn_tab[0.0]}[scheme]


def ferm_heavy(scheme, m):
    if scheme == "gestaffelt":
        return kn_tab[m] / 8
    if scheme == "naiv":
        return kn_tab[m]
    if scheme == "wilson":
        return kw[m]
    r = 0.5 if scheme == "wilson_r0.5" else 2.0
    return kw_r[r][m] if m in kw_r[r] else float("nan")


def gauge(choice):
    kind = choice[0]
    if kind == "ohne_zentrum":
        return 12 * gam_nc
    if kind == "ez_dest":
        return 12 * gam_ez
    if kind == "kabat":
        return 12 * (-4 * ks0)
    _, g, neff, mode, u1 = choice                      # ("kante", g, n_eff, innovation|marginal, ln6 ja/nein)
    return 12 * gam_ez + neff * edge[g][mode] + (np.log(6) if u1 else 0.0)


def total(ferm, heavy, gch, scalars=True):
    k = 24 * ferm_light(ferm) + gauge(gch)
    if heavy is not None:
        k += 8 * ferm_heavy(ferm, heavy) + 2 * ks_m[heavy]
    if scalars:
        k += 4 * ks0
    return float(k)


def a_lP(k):
    return 2 * np.sqrt(k) if k > 0 else float("nan")


KANTE_REF = ("kante", 0.55, 12, "innovation", False)
CENTRAL = ("wilson", 0.5, KANTE_REF)
say("\n   Bausteine (nat je Gitterfläche):")
say(f"   24 leichte Dirac: gestaffelt {24 * ferm_light('gestaffelt'):.3f}, Wilson r=1 {24 * kw[0.0]:.3f}, naiv {24 * kn_tab[0.0]:.3f}")
say("   16-Plett + Phi bei m = 0,3/0,5/1: Wilson " + "/".join(f"{8 * kw[m] + 2 * ks_m[m]:.3f}" for m in (0.3, 0.5, 1.0))
    + ", naiv " + "/".join(f"{8 * kn_tab[m] + 2 * ks_m[m]:.3f}" for m in (0.3, 0.5, 1.0)))
say(f"   4 Higgs: {4 * ks0:.3f}. 12 Eichbosonen: ohne Zentrum {12 * gam_nc:.3f}, destillierbar {12 * gam_ez:.3f},"
    f" Kante bei g = 0,55 zusätzlich {12 * edge[0.55]['innovation']:.2f}")

say("\n   a/l_P (Zeilen: Fermionschema und schwerer Sektor; Spalten: Eichalgebra):")
cols = [("ohne Zentrum", ("ohne_zentrum",)), ("EZ dest.", ("ez_dest",)),
        ("Kante g=1", ("kante", 1.0, 12, "innovation", False)),
        ("Kante g=0,55", KANTE_REF), ("Kante g=0,5", ("kante", 0.5, 12, "innovation", False))]
say("   " + " " * 28 + "".join(f"{c:>14s}" for c, _ in cols))
table_rows = []
for ferm in ("gestaffelt", "wilson", "naiv"):
    for heavy in (None, 0.5):
        lab = f"{ferm}, schwer {'aus' if heavy is None else 'm=0,5'}"
        vals = [total(ferm, heavy, gc) for _, gc in cols]
        table_rows.append({"fermionen": ferm, "schwer": heavy, "kappa": dict(zip([c for c, _ in cols], vals)),
                           "a_lP": dict(zip([c for c, _ in cols], [a_lP(v) for v in vals]))})
        say(f"   {lab:28s}" + "".join(f"{a_lP(v):14.2f}" for v in vals))

# Band über alle Schemata
ferms = ("gestaffelt", "wilson", "wilson_r0.5", "wilson_r2", "naiv")
heavies = (None, 0.3, 0.5, 1.0)
gauges = [("ohne_zentrum",), ("ez_dest",)] + [("kante", g, n, "innovation", False) for g in (0.5, 0.55, 0.75, 1.0)
                                              for n in (8, 12)] + [("kante", 0.5, 12, "marginal", True)]
allk = []
for f_ in ferms:
    for h_ in heavies:
        if h_ is not None and np.isnan(ferm_heavy(f_, h_)):
            continue
        for g_ in gauges:
            allk.append((total(f_, h_, g_), f_, h_, g_))
kmin, kmax = min(allk), max(allk)
kc = total(*CENTRAL)
say(f"\n   H6-eigenes Schema (Wilson r = 1 wie in Z7-Z10, schwer bei m = 0,5, erweiterter Hilbertraum, g = 0,55):")
say(f"   kappa = {kc:.2f}, a/l_P = {a_lP(kc):.2f}.")
core = [total("wilson", h_, ("kante", g, n, "innovation", False)) for h_ in (0.3, 0.5, 1.0)
        for g in (0.5, 0.55, 0.75, 1.0) for n in (8, 12)]
say(f"   Kernband (Wilson r = 1, schwer ein, Kante ein mit g in [0,5; 1], 8 bis 12 Kopien): "
    f"a/l_P = {a_lP(min(core)):.2f} bis {a_lP(max(core)):.2f}")
say(f"   Gesamtband über alle Schemata: kappa = {kmin[0]:.2f} bis {kmax[0]:.2f}, a/l_P = {a_lP(kmin[0]):.2f} bis {a_lP(kmax[0]):.2f}")
def lab(t):
    return f"{t[1]}, schwer {'aus' if t[2] is None else t[2]}, Eich {' '.join(str(x) for x in t[3])}"


say(f"   Minimum: {lab(kmin)}. Maximum: {lab(kmax)}.")
k_kabat = total("wilson", 0.5, ("kabat",))
k_kabat_min = total("gestaffelt", None, ("kabat",))
say(f"   Außerhalb des Gitterschemas (Kabat-Vorzeichen, Vektor = -4 kappa_s): kappa = {k_kabat:.2f}, a/l_P = {a_lP(k_kabat):.2f};")
say(f"   mit gestaffelten Fermionen ohne schweren Sektor kappa = {k_kabat_min:.2f}, a/l_P = {a_lP(k_kabat_min):.2f}.")

# Empfindlichkeit
say("\n   Empfindlichkeit: je eine Wahl variiert, Rest wie H6-eigenes Schema (a/l_P von bis):")
axes = {
    "Kante (aus / obere Schranke g=0,5 + ln6)": [total("wilson", 0.5, ("ez_dest",)), total("wilson", 0.5, ("ohne_zentrum",)),
                                                total("wilson", 0.5, ("kante", 0.5, 12, "marginal", True))],
    "Verdoppler (gestaffelt / Wilson r / naiv)": [total(f_, 0.5, KANTE_REF) for f_ in ferms],
    "Kopplung g (0,5 bis 1,0)": [total("wilson", 0.5, ("kante", g, 12, "innovation", False)) for g in (0.5, 0.75, 1.0)],
    "nichtabelsche Zählung (8 / 12 Kopien)": [total("wilson", 0.5, ("kante", 0.55, n, "innovation", False)) for n in (8, 12)],
    "schwerer Sektor (aus / m = 0,3 bis 1)": [total("wilson", h_, KANTE_REF) for h_ in heavies],
    "Higgs und Phi (aus / ein)": [total("wilson", 0.5, KANTE_REF, scalars=False) - 2 * ks_m[0.5], kc],
}
sens = {}
for name, vals in axes.items():
    lo_, hi_ = min(vals + [kc]), max(vals + [kc])
    sens[name] = {"kappa": [lo_, hi_], "a_lP": [a_lP(lo_), a_lP(hi_)], "breite_a_lP": a_lP(hi_) - a_lP(lo_)}
for name, v in sorted(sens.items(), key=lambda t: -t[1]["breite_a_lP"]):
    say(f"   {name:45s} {v['a_lP'][0]:6.2f} bis {v['a_lP'][1]:6.2f}  (Breite {v['breite_a_lP']:.2f})")
order = sorted(sens, key=lambda n: -sens[n]["breite_a_lP"])

# SM-Kopplung an der Gitterskala (Referenz für g)
inv_lP = 1.2209e19
mu_c = inv_lP / a_lP(kc)
L = np.log(mu_c / 91.1876)
ainv = np.array([59.0, 29.6, 8.47]) - np.array([41 / 10, -19 / 6, -7.0]) / (2 * np.pi) * L
gs = np.sqrt(4 * np.pi / ainv)
say(f"\n   Gitterskala im H6-Schema: 1/a = {mu_c:.2e} GeV. SM-Einschleifenkopplungen dort (Eingabe MZ-Werte):")
say(f"   g1 (GUT-normiert) = {gs[0]:.3f}, g2 = {gs[1]:.3f}, g3 = {gs[2]:.3f}. Daher g = 0,55 als Referenz.")
say("   H6 selbst legt g an der Gitterskala nicht fest (induzierte Kopplung aus Fermionschleifen, Y1).\n")

# Befunde
say("Befunde (nicht erzwungen):")
say(f"   Gitter in Planckeinheiten möglich (kappa_min <= 1/4): {befund('Gitter in Planckeinheiten möglich', kmin[0] <= 0.25)}")
say(f"   a/l_P eindeutig (Bandbreite unter 20 %): {befund('a/l_P eindeutig festgelegt', a_lP(kmax[0]) / a_lP(kmin[0]) < 1.2)}")
dom = 12 * edge[0.55]["innovation"] > kc - 12 * edge[0.55]["innovation"]
say(f"   Kante dominiert im H6-Schema (mehr als alles andere zusammen): {befund('Kante dominiert kappa im H6-Schema', dom)}")
say(f"   Wilson-Geschmack gleich gestaffeltem Geschmack: {befund('Verdopplerentfernung schemaunabhängig', abs(kw[0.0] / (kn_tab[0.0] / 8) - 1) < 0.1)}")
say(f"   Destillierbarer TM-Anteil unter 5 % von kappa_s: {befund('TM destillierbar vernachlässigbar', ktm0 / ks0 < 0.05)}")
same_sign = np.sign(gauge(("kabat",))) == np.sign(gauge(KANTE_REF))
say(f"   Vorzeichen des Eichbeitrags schemaunabhängig: {befund('Vorzeichen Eichbeitrag schemaunabhängig', same_sign)}"
    " (Kabat-Schema gegen Gitter)")
say(f"   Größte Empfindlichkeit: {order[0]}; dann {order[1]}.\n")

say("Urteil:")
say(f"   a/l_P ist in H6 keine bestimmte Zahl. Robust ist nur die untere Schranke a/l_P >= {a_lP(kmin[0]):.1f} in jedem")
say(f"   Gitterschema. Mit Kabats Vorzeichen (außerhalb des Gitters) fiele sie auf {a_lP(k_kabat_min):.1f}; 1 wird nie erreicht.")
say(f"   Im H6-eigenen Schema liegt a/l_P bei etwa {a_lP(kc):.0f} (Kernband {a_lP(min(core)):.1f} bis {a_lP(max(core)):.1f}),")
say(f"   über alle Schemata bei {a_lP(kmin[0]):.1f} bis {a_lP(kmax[0]):.1f}.")
say("   Den Ausschlag gibt die Kante der Eichfelder: die Shannon-Entropie des ganzzahligen Randflusses,")
say("   rund 1,8 nat je Link und Eichboson bei g = 0,55, mit ln(1/g) wachsend. Danach kommen die Verdoppler.")
say("   Für eine Zahl braucht D drei Festlegungen: die Zentrumswahl als Teil von eta, g an der Gitterskala")
say("   aus H6 selbst und die chirale Gitterformulierung (Wilson r, gestaffelt, Domänenwand).")
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)}). Laufzeit {time.time() - T0:.0f} s.")

(HERE / "ausgabe_z15.txt").write_text("\n".join(LINES) + "\n")
(HERE / "z15_ergebnis.json").write_text(json.dumps({
    "kappa_skalar_masselos": ks0, "kappa_skalar_massiv": {str(m): v for m, v in ks_m.items()},
    "kappa_skalar_konvergenz": {f"ring{b}_nk{n}": v for (b, n), v in conv_s.items()},
    "kappa_TM_destillierbar": ktm0, "kappa_photon_ohne_zentrum": gam_nc, "kappa_photon_ez_destillierbar": gam_ez,
    "randfluss": {"varianz_e": C0_marg, "innovationsvarianz_e": C0_geo, "h_e": h_e},
    "kante_je_eichboson": {str(g): v for g, v in edge.items()},
    "kappa_dirac_naiv": {str(m): v for m, v in kn_tab.items()}, "kappa_dirac_naiv_direkt_128": kn_dir,
    "kappa_dirac_gestaffelt_je_geschmack": kn_tab[0.0] / 8,
    "kappa_dirac_wilson_r1": {str(m): v for m, v in kw.items()},
    "kappa_dirac_wilson_r": {str(r): {str(m): v for m, v in d.items()} for r, d in kw_r.items()},
    "szenario_tabelle": table_rows,
    "h6_schema": {"kappa": kc, "a_lP": a_lP(kc), "kernband_a_lP": [a_lP(min(core)), a_lP(max(core))]},
    "gesamtband": {"kappa": [kmin[0], kmax[0]], "a_lP": [a_lP(kmin[0]), a_lP(kmax[0])],
                   "min_schema": [kmin[1], kmin[2], list(kmin[3])], "max_schema": [kmax[1], kmax[2], list(kmax[3])]},
    "kabat_variante": {"kappa": k_kabat, "a_lP": a_lP(k_kabat), "kappa_min": k_kabat_min, "a_lP_min": a_lP(k_kabat_min)},
    "empfindlichkeit": sens, "sm_kopplungen_an_gitterskala": {"mu_GeV": mu_c, "g1_g2_g3": gs.tolist()},
    "checks": CHECKS, "befunde": BEFUNDE}, indent=1, ensure_ascii=False, default=float))
assert all(CHECKS.values())
