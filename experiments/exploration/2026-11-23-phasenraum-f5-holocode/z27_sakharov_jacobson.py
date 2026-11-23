"""Z27: Sakharov gegen Jacobson im selben Gitterregler (z21 S5, Prämisse D).

Unter D ist G durch die Verschränkungsdichte festgelegt: 1/(4G) = kappa/a^2 (Z2, Z15). Derselbe Feldinhalt
induziert über eine Schleife einen Einstein-Hilbert-Term 1/(16 pi G_ind). Gefragt ist das Verhältnis
R = (kappa/a^2) / (1/(4 G_ind)) je Feldtyp und für den vollen H6-Inhalt, im selben Regler wie Z2/Z15
(Hamiltonsches Gitter Z^3, kontinuierliche Zeit, Wilson r = 1).

Weg (begründet in Teil 1). Sakharov-Seite als statische Antwort der Grundzustandsenergie auf eine schwache,
langsam variierende räumliche TT-Metrik gamma_11 = -gamma_22 = eps cos(k z), eingebaut als Modulation der
Hüpf-, Plaketten- und Platzgewichte (sqrt(h) h^ij usw.). Für statische Metriken ist E_0 = -(1/16 pi G) Int sqrt(h) R^(3)
+ ..., und für diese TT-Mode gilt E/V = eps^2 k^2/(64 pi G) + O(k^0, k^4). Also sigma := a^2/(4 G_ind) = 16 pi e_2,
e_2 = k^2-Koeffizient von E^(2)/(V eps^2). Für freie Felder ist das eine Blase (zweite Ordnung) plus Kontaktterme,
im Impulsraum exakt. Die konische Methode wird nicht benutzt: Auf dem Gitter ist der Kegel die Replika-
Konstruktion selbst und gibt kappa zurück; der Vergleich wäre zirkulär.

Teil 0: Kontinuumskalibrierung (Wärmekern): Kegel über die Bildmethode, a_1 aus S^2-Wärmekernen, Kabats Kontaktterm.
Teil 1: Gitterformeln, geprüft gegen direkte Diagonalisierung und unabhängige Quadratur.
Teil 2: Universeller Teil auf dem Gitter: m^2 ln m beider Seiten (Skalar, naives Dirac).
Teil 3: Nicht-universeller Teil je Feld und Kopplungsregel; Kopplungsmehrdeutigkeit.
Teil 4: Eichfelder mit und ohne Randterm.
Teil 5: Voller H6-Inhalt (3 x 16 mit nu_R, 12 Eichbosonen, Higgs, schwere Eltern aus Z23 getrennt).
"""
import json
import textwrap
import time
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import CubicSpline
from scipy.linalg import toeplitz
from scipy.optimize import brentq

np.seterr(all="ignore")
HERE = Path(__file__).parent
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
T0 = time.time()
Z15 = json.loads((HERE / "z15_ergebnis.json").read_text())
Z2 = json.loads((HERE / "z2_ergebnis.json").read_text())


def say(s=""):
    print(s, flush=True)


def check(name, cond):
    CHECKS[name] = bool(cond)
    return "ja" if cond else "NEIN"


def befund(name, val):
    BEFUNDE[name] = bool(val)
    return "ja" if val else "nein"


def f4(x):
    return f"{x:.4f}"


# ===================================================================================================
say("Z27: Sakharov gegen Jacobson im selben Gitterregler (Prämisse D, z21 S5)")
say("Unter D: 1/(4G) = kappa/a^2, also a/l_P = 2 sqrt(kappa). Sakharov: derselbe Inhalt induziert 1/(16 pi G_ind).")
say("Gemessen wird sigma := a^2/(4 G_ind) in nat je Gitterfläche (gleiche Einheit wie kappa) und R = kappa/sigma.\n")

# ---------------------------------------------------------------------------------------------------
# Teil 0: Kontinuum
# ---------------------------------------------------------------------------------------------------
say("Teil 0: Kontinuumskalibrierung mit Wärmekernregler (Eigenzeit s >= eps)")
say("   Kegel mit Winkel 2 pi alpha, alpha = 1/n, über die Bildmethode (Z_n-Orbifold):")
say("   Tr K_Kegel - alpha Tr K_flach = (1/n) sum_j tr(Drehung_j) / (4 sin^2(pi j/n)) je Querfläche A/(4 pi s).")
ns = range(2, 10)
img_s = max(abs(sum(1 / np.sin(np.pi * j / n) ** 2 for j in range(1, n)) - (n * n - 1) / 3) for n in ns)
img_f = max(abs(sum((-1) ** j * np.cos(np.pi * j / n) / np.sin(np.pi * j / n) ** 2 for j in range(1, n)) + (n * n - 1) / 6)
            for n in (3, 5, 7, 9))
img_v = max(abs(sum(2 * np.cos(2 * np.pi * j / n) / (4 * np.sin(np.pi * j / n) ** 2) for j in range(1, n)) / n
                - ((n * n - 1) / (6 * n) - (n - 1) / n)) for n in ns)
say(f"   Skalar: sum 1/sin^2 = (n^2-1)/3, also (1/alpha - alpha)/12 (Abweichung {img_s:.0e}).")
say(f"   Dirac (Drehung um 2 pi = -1, tr S = 4 cos(theta/2), n ungerade): -(1/alpha - alpha)/6 (Abweichung {img_f:.0e});")
say("   mit dem Fermionvorzeichen in ln Z also 2 Skalareinheiten.")
say(f"   Vektor minus 2 Geister (tr = 2 cos theta): 2 (1/alpha - alpha)/12 - (1 - alpha) (Abweichung {img_v:.0e}).")
say("   Bei alpha -> 1 ist eine Skalareinheit (1 - alpha)/6, der Zusatz -(1 - alpha) also -6 Einheiten (Kabat).")
say(f"   Bildmethode exakt: {check('Kontinuum: Kegel-Bildmethode Skalar/Dirac/Vektor exakt (Kabat -6)', max(img_s, img_f, img_v) < 1e-10)}")
s_hk = 1e-3
ll = np.arange(0, 6000)
hk_s = float(np.sum((2 * ll + 1) * np.exp(-s_hk * ll * (ll + 1))) - 1 / s_hk)
hk_f = float(np.sum(4 * ll[1:] * np.exp(-s_hk * ll[1:] ** 2)) - 2 / s_hk)
hk_v = float(2 * np.sum((2 * ll[1:] + 1) * np.exp(-s_hk * ll[1:] * (ll[1:] + 1))) - 2 / s_hk)
say("   Glatter a_1-Koeffizient auf der Einheits-S^2 (R = 2, Fläche 4 pi), Tr e^{-s Delta} - tr(1)/s -> Int a_1/(4 pi):")
say(f"   Skalar {hk_s:.5f} (R/6: 1/3), Spinor {hk_f:.5f} (tr(1/6 - 1/4) R: -1/3), 1-Form {hk_v:.5f} (tr(R/6 - Ric): -4/3)")
check("Kontinuum: a_1 auf S^2 (Skalar 1/3, Spinor -1/3, 1-Form -4/3)",
      abs(hk_s - 1 / 3) < 2e-4 and abs(hk_f + 1 / 3) < 2e-4 and abs(hk_v + 4 / 3) < 5e-4)
say("   In 3+1 D: a_1(Vektor) - 2 a_1(Geist) = (4/6 - 1) R - (2/6) R = -(2/3) R = -4 Skalareinheiten.")
say("   Gleicher Regler: S/A = (Kegelkoeffizient)/(48 pi eps), 1/(4 G_ind) = (a_1 in R/6-Einheiten)/(48 pi eps). Damit")
say("   R_Skalar = 1, R_Dirac = 2/2 = 1, R_Vektor = 2/(-4) = -1/2 ohne Randmoden und (2 - 6)/(-4) = 1 mit Kabats")
say("   Kontaktterm (Donnelly-Wall: das sind die Randmoden, im Wärmekernschema mit negativem Vorzeichen).")
cont_R = {"skalar": 1.0, "dirac": 1.0, "vektor_ohne_rand": -0.5, "vektor_mit_kontakt": 1.0}
say("   Die Kontinuumsidentität gilt also für Spin 0 und 1/2 und für Eichfelder nur mit dem Kontaktterm.\n")

# ---------------------------------------------------------------------------------------------------
# Teil 1: Gitterformeln
# ---------------------------------------------------------------------------------------------------
say("Teil 1: Sakharov-Seite auf dem Gitter, statische TT-Antwort")
say("   Störung gamma_11 = -gamma_22 = eps cos(k z), k entlang z (transversal, spurfrei, diagonal: koppelt an die")
say("   x- und y-Gewichte). Kontinuumsgewichte: Skalar-Link i: sqrt(h) h^ii, Platz: pi^2/sqrt(h), m^2 sqrt(h);")
say("   Photon E_i und B_i: h_ii/sqrt(h); Dirac (kanonisches chi = h^(1/4) psi): naives Hüpfen e_i^i = h_ii^(-1/2),")
say("   Wilsonterm h_ii^(-c_W) (c_W = 1: kovarianter Laplace; 1/2: wie das naive Hüpfen; 0: metrikblind).")
say("   E^(2) = <V_2> - sum |<n|V_1|0>|^2/(E_n - E_0). Ultralokale <V_2> hängen nicht von k ab.")
say("   Blase: Bosonen -(1/16) <sum |Y/sqrt(w w') - X sqrt(w w')|^2/(w + w')>, Dirac -(1/2) <[g^2(1 + d.d') - 2(d.g)(d'.g)]/(E + E')>.")
say("   Ort der Gitterobjekte in z: ganzzahlig (x-,y-Links, Plätze, xy-Plaketten), halbzahlig (z-Links, xz-/yz-Plaketten).")
say("   Kopplungsregeln (wo die Metrik sitzt):")
say("     K  stetig an der Objektmitte abgetastet (kanonisch);")
say("     P  auf den Gitterplätzen, halbzahlige Objekte linear interpoliert (Faktor cos(k/2));")
say("     D  auf den Würfelmitten, ganzzahlige Objekte interpoliert (Faktor cos(k/2)).")
say("   Kontaktterme: kov (Gewichtsfunktion bis O(gamma^2)) oder no (Metrik koppelt an O - <O>, Vakuum lokal abgezogen).\n")

DSTEP = 2e-3


def nodes(lo=1e-4, n_gl=8, both=False):
    pts = [0.0] + list(np.geomspace(lo, 0.5, 14)) + list(np.linspace(0.5, np.pi / 2 if both else np.pi, 4)[1:])
    if both:
        pts = pts + [np.pi - x for x in pts]
    pts = np.array(sorted(set(np.round(pts, 14))))
    x, w = leggauss(n_gl)
    X = np.concatenate([(b - a) / 2 * x + (a + b) / 2 for a, b in zip(pts[:-1], pts[1:])])
    W = np.concatenate([(b - a) / 2 * w for a, b in zip(pts[:-1], pts[1:])])
    return X, W / np.pi


def bz(fun, n_gl=8, lo=1e-4, both=False):
    """Mittel über die Brillouinzone; der Integrand ist in jedem p_i gerade, also [0, pi]^3."""
    x, w = nodes(lo, n_gl, both)
    P2, P3 = np.meshgrid(x, x, indexing="ij")
    W23 = np.outer(w, w)
    out = None
    for i in range(len(x)):
        r = fun(x[i], P2, P3)
        r = np.array([np.sum(W23 * q) for q in r]) if isinstance(r, tuple) else np.sum(W23 * r)
        out = w[i] * r if out is None else out + w[i] * r
    return out


def k2(fk):
    """Punktweiser k^2-Koeffizient (halbe zweite Differenz; f(k) + f(-k) ist in p_3 gerade)."""
    return lambda a, b, c: (fk(a, b, c, DSTEP) + fk(a, b, c, -DSTEP) - 2 * fk(a, b, c, 0.0)) / (2 * DSTEP * DSTEP)


ONE = lambda k: 1.0 + 0.0 * k
COS = lambda k: np.cos(k / 2)
FACT = {"K": (ONE, ONE), "P": (ONE, COS), "D": (COS, ONE)}


# ---------------- Skalar
def scal_bub(p1, p2, p3, k, m, ai, ah):
    v = -4 * np.sin(p1 / 2) ** 2 + 4 * np.sin(p2 / 2) ** 2
    s12 = 4 * np.sin(p1 / 2) ** 2 + 4 * np.sin(p2 / 2) ** 2 + m * m
    w = np.sqrt(s12 + 4 * np.sin(p3 / 2) ** 2)
    wp = np.sqrt(s12 + 4 * np.sin((p3 + k) / 2) ** 2)
    return -(1 / 16) * (ai(k) * v) ** 2 / (w * wp * (w + wp))


def scal_dia(p1, p2, p3, m):
    s1, s2, s3 = 4 * np.sin(p1 / 2) ** 2, 4 * np.sin(p2 / 2) ** 2, 4 * np.sin(p3 / 2) ** 2
    w = np.sqrt(s1 + s2 + s3 + m * m)
    # <gamma^2> = 1/2; c2: pi^2/sqrt h +1/2, m^2 sqrt h -1/2, x/y-Link +1/2, z-Link -1/2
    return 0.5 * (w / 8 - m * m / (8 * w) + (s1 + s2) / (8 * w)), 0.5 * (-s3 / (8 * w)), w / 2


# ---------------- Photon (Kogut-Susskind, nichtkompakt, schwache Kopplung)
DD = np.diag([1.0, -1.0, 0.0])


def _cross(v):
    Z = np.zeros(v.shape[1:] + (3, 3))
    Z[..., 0, 1], Z[..., 0, 2], Z[..., 1, 0] = -v[2], v[1], v[2]
    Z[..., 1, 2], Z[..., 2, 0], Z[..., 2, 1] = -v[0], -v[1], v[0]
    return Z


def phot_bub(p1, p2, p3, k, m, ai, ah):
    s = np.stack(np.broadcast_arrays(2 * np.sin(p1 / 2), 2 * np.sin(p2 / 2), 2 * np.sin(p3 / 2)))
    sp = np.stack(np.broadcast_arrays(2 * np.sin(p1 / 2), 2 * np.sin(p2 / 2), 2 * np.sin((p3 + k) / 2)))
    w, wp = np.sqrt((s ** 2).sum(0)), np.sqrt((sp ** 2).sum(0))
    sq = np.sqrt(w * wp)[..., None, None]
    Z = ah(k) * np.swapaxes(_cross(s), -1, -2) @ DD @ _cross(sp) / sq - ai(k) * DD * sq
    sh, shp = np.moveaxis(s / w, 0, -1), np.moveaxis(sp / wp, 0, -1)
    PT = np.eye(3) - sh[..., :, None] * sh[..., None, :]
    PTp = np.eye(3) - shp[..., :, None] * shp[..., None, :]
    val = np.einsum("...ij,...jk,...lk,...il->...", PT, Z, PTp, Z)
    return np.nan_to_num(-(1 / 16) * val / (w + wp))


def phot_dia(p1, p2, p3, m=0.0):
    s = np.stack(np.broadcast_arrays(2 * np.sin(p1 / 2), 2 * np.sin(p2 / 2), 2 * np.sin(p3 / 2)))
    w = np.sqrt((s ** 2).sum(0))
    O = np.nan_to_num(sum((w / 4) * (1 - s[j] ** 2 / w ** 2) for j in range(3)))   # sum_j <E_j^2/2> = sum_j <B_j^2/2>
    return 0.25 * O, 0.25 * O, w                     # ganz: E1, E2, B3; halb: E3, B1, B2 (c2 = 1/2 alle)


# ---------------- Wilson-Dirac (4x4, alpha_1..3, beta)
def wil_bub(p1, p2, p3, k, m, ai, ah, r=1.0, cW=1.0):
    def dvec(q3):
        return np.stack(np.broadcast_arrays(np.sin(p1), np.sin(p2), np.sin(q3),
                                            m + r * (3 - np.cos(p1) - np.cos(p2) - np.cos(q3))))
    d, dp = dvec(p3), dvec(p3 + k)
    E, Ep = np.sqrt((d ** 2).sum(0)), np.sqrt((dp ** 2).sum(0))
    g = ai(k) * np.stack(np.broadcast_arrays(-0.5 * np.sin(p1), 0.5 * np.sin(p2), 0 * p3, cW * r * (np.cos(p1) - np.cos(p2))))
    g2, dd = (g ** 2).sum(0), (d * dp).sum(0) / (E * Ep)
    return np.nan_to_num(-0.5 * (g2 * (1 + dd) - 2 * ((d * g).sum(0) / E) * ((dp * g).sum(0) / Ep)) / (E + Ep))


def wil_dia(p1, p2, p3, m, r=1.0, cW=1.0):
    d = np.stack(np.broadcast_arrays(np.sin(p1), np.sin(p2), np.sin(p3), m + r * (3 - np.cos(p1) - np.cos(p2) - np.cos(p3))))
    E = np.sqrt((d ** 2).sum(0))
    dia = 0.5 * sum((3 / 8) * (-2 * np.sin(q) ** 2 / E) + cW * (cW + 1) / 2 * (-2 * d[3] * r * (1 - np.cos(q)) / E)
                    for q in (p1, p2))
    return np.nan_to_num(dia), 0 * E, -2 * E


# ---------------- direkte Diagonalisierung auf L x L x Nz (periodisch)
LD, NZ = 4, 8
SITES = [(x, y, z) for x in range(LD) for y in range(LD) for z in range(NZ)]
SIDX = {s: i for i, s in enumerate(SITES)}
NS = len(SITES)


def _sh(s, j):
    t = list(s)
    t[j] = (t[j] + 1) % (LD, LD, NZ)[j]
    return tuple(t)


def gam_at(zc, half, rule, eps, k):
    g = lambda z: eps * np.cos(k * z)
    if (rule == "P" and half) or (rule == "D" and not half):
        return 0.5 * (g(zc - 0.5) + g(zc + 0.5))
    return g(zc)


def direct_scalar(eps, k, m, rule):
    B, A = np.zeros((NS, NS)), np.zeros(NS)
    for s in SITES:
        a, z = SIDX[s], s[2]
        g = gam_at(z, False, rule, eps, k)
        A[a] = 1 / np.sqrt(1 - g * g)
        B[a, a] += np.sqrt(1 - g * g) * m * m
        for j in range(3):
            gl = gam_at(z + 0.5, True, rule, eps, k) if j == 2 else g
            h = (1 + gl, 1 - gl, 1.0)
            wl = np.sqrt(h[0] * h[1] * h[2]) / h[j]
            b = SIDX[_sh(s, j)]
            B[a, a] += wl
            B[b, b] += wl
            B[a, b] -= wl
            B[b, a] -= wl
    return 0.5 * np.sum(np.sqrt(np.linalg.eigvalsh(np.sqrt(A)[:, None] * B * np.sqrt(A)[None, :])))


_sx, _sy, _sz = np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.array([[1, 0], [0, -1]])
ALPHA = [np.kron(_sx, q) for q in (_sx, _sy, _sz)]
BETA = np.kron(_sz, np.eye(2))


def direct_wilson(eps, k, m, rule, r=1.0, cW=1.0):
    H = np.zeros((4 * NS, 4 * NS), complex)
    for s in SITES:
        a = SIDX[s]
        g = gam_at(s[2], False, rule, eps, k)
        H[4 * a:4 * a + 4, 4 * a:4 * a + 4] += m * BETA
        for j in range(3):
            h = (1 + g, 1 - g, 1.0)[j]
            wN, wW = h ** -0.5, h ** -cW
            b = SIDX[_sh(s, j)]
            T = wN * 0.5j * ALPHA[j] - wW * 0.5 * r * BETA
            H[4 * b:4 * b + 4, 4 * a:4 * a + 4] += T
            H[4 * a:4 * a + 4, 4 * b:4 * b + 4] += T.conj().T
            H[4 * a:4 * a + 4, 4 * a:4 * a + 4] += wW * r * BETA
    ev = np.linalg.eigvalsh(H)
    return float(np.sum(ev[ev < 0]))


PLAQ = [(s, j, a, b) for s in SITES for j, (a, b) in enumerate(((1, 2), (2, 0), (0, 1)))]
CURL = np.zeros((len(PLAQ), 3 * NS))
for _i, (_s, _j, _a, _b) in enumerate(PLAQ):
    CURL[_i, 3 * SIDX[_sh(_s, _a)] + _b] += 1
    CURL[_i, 3 * SIDX[_s] + _b] -= 1
    CURL[_i, 3 * SIDX[_sh(_s, _b)] + _a] -= 1
    CURL[_i, 3 * SIDX[_s] + _a] += 1


def direct_photon(eps, k, m, rule):
    wE, wB = np.zeros(3 * NS), np.zeros(len(PLAQ))
    for s in SITES:
        for j in range(3):
            g = gam_at(s[2] + (0.5 if j == 2 else 0.0), j == 2, rule, eps, k)
            h = (1 + g, 1 - g, 1.0)
            wE[3 * SIDX[s] + j] = h[j] / np.sqrt(h[0] * h[1] * h[2])
    for i, (s, j, a, b) in enumerate(PLAQ):
        half = 2 in (a, b)
        g = gam_at(s[2] + (0.5 if half else 0.0), half, rule, eps, k)
        h = (1 + g, 1 - g, 1.0)
        wB[i] = np.sqrt(h[0] * h[1] * h[2]) / (h[a] * h[b])
    M = (np.sqrt(wE)[:, None] * CURL.T) @ (wB[:, None] * CURL) * np.sqrt(wE)[None, :]
    ev = np.linalg.eigvalsh(M)
    nz = NS - 1 + 3                                  # Eichmoden + drei flache Holonomien
    return 0.5 * np.sum(np.sqrt(np.clip(ev[nz:], 0, None)))


def formula_finite(bub, dia, k, m, rule, kontakt="kov"):
    pa = 2 * np.pi * np.arange(LD) / LD
    pz = 2 * np.pi * np.arange(NZ) / NZ
    P1, P2, P3 = np.meshgrid(pa, pa, pz, indexing="ij")
    ai, ah = FACT[rule]
    dI, dH, _ = dia(P1, P2, P3, m)
    val = bub(P1, P2, P3, k, m, ai, ah).mean()
    if kontakt == "kov":
        val += ai(k) ** 2 * np.mean(dI) + ah(k) ** 2 * np.mean(dH)
    return float(val)


say("   Prüfung gegen direkte Diagonalisierung (Gitter 4 x 4 x 8, eps = 2e-3, k = 2 pi n/8), je Platz:")
dmax = 0.0
eps_d = 2e-3
for lab, dfun, bub, dia, m, rules in (("Skalar m=0,3", direct_scalar, scal_bub, scal_dia, 0.3, ("K", "P", "D")),
                                     ("Wilson m=0,2", direct_wilson, wil_bub, wil_dia, 0.2, ("K", "D")),
                                     ("Photon", direct_photon, phot_bub, phot_dia, 0.0, ("K", "P", "D"))):
    parts = []
    for rule in rules:
        for n in (1, 2):
            k = 2 * np.pi * n / NZ
            v = [dfun(e, k, m, rule) for e in (eps_d, -eps_d, 0.0)]
            direct = (v[0] + v[1] - 2 * v[2]) / (2 * eps_d * eps_d) / NS
            form = formula_finite(bub, dia, k, m, rule)
            dmax = max(dmax, abs(direct - form))
            parts.append(f"{rule}{n}: {direct:.6f}/{form:.6f}")
    say(f"   {lab}: " + ", ".join(parts))
say(f"   Größte Abweichung {dmax:.1e} (endliches eps). Formeln samt Kontakttermen und Regeln stimmen: "
    f"{check('Gitter: Blase + Kontaktterme = direkte Diagonalisierung (Skalar, Wilson, Photon; Regeln K, P, D)', dmax < 2e-5)}")


def uniform_e2(bub, m, N=96, ns_=(1, 2, 3, 4, 5)):
    p = 2 * np.pi * (np.arange(N) + 0.5) / N
    P1, P2, P3 = np.meshgrid(p, p, p, indexing="ij")
    ks = np.array([2 * np.pi * n / N for n in ns_])
    e = np.array([bub(P1, P2, P3, kk, m, ONE, ONE).mean() for kk in ks])
    A = np.stack([np.ones_like(ks), ks ** 2, ks ** 4, ks ** 6], 1)
    return float(np.linalg.lstsq(A, e, rcond=None)[0][1])


u_s = uniform_e2(scal_bub, 0.2)
g_s = bz(k2(lambda a, b, c, k: scal_bub(a, b, c, k, 0.2, ONE, ONE)), lo=1e-2)
u_w = uniform_e2(wil_bub, 0.2, N=64)
g_w = bz(k2(lambda a, b, c, k: wil_bub(a, b, c, k, 0.2, ONE, ONE)), lo=1e-2)
say(f"   Unabhängige Quadratur (gleichmäßiges Gitter, Anpassung in k) gegen Gauß-Legendre mit Punktdifferenz, m = 0,2:")
say(f"   Skalar e_2 = {u_s:.7f}/{g_s:.7f}, Wilson e_2 = {u_w:.7f}/{g_w:.7f}")
say(f"   Übereinstimmung: {check('Quadratur: gleichmäßiges Gitter = Gauß-Legendre für e_2 (rel. < 1e-3)', abs(u_s / g_s - 1) < 1e-3 and abs(u_w / g_w - 1) < 1e-3)}")
say("   Begründung des Weges: Die statische Antwort ist im selben Regler wie kappa definiert (Hamiltonsches Gitter,")
say("   gleiche Felder, gleiche Hüpfterme); sie misst den Koeffizienten von sqrt(h) R^(3). Zeitableitungsterme")
say("   (K_ij K^ij), die ein Gitter mit eigenem Koeffizienten tragen kann, sind nicht gemessen.\n")

# ---------------------------------------------------------------------------------------------------
# Teil 2: universeller Teil
# ---------------------------------------------------------------------------------------------------
say("Teil 2: Universeller Teil auf dem Gitter (m^2 ln m, reglerunabhängig)")
say("   Kontinuum (Wärmekern): 1/(4 G_ind) und S/A enthalten je Skalar (m^2/(24 pi)) ln m, je Dirac das Doppelte.")
say("   Die statische TT-Blase gibt im Kontinuum denselben Koeffizienten: e_2 enthält m^2 ln m/(384 pi^2) (analytisch:")
say("   Winkelmittel 4/105 und 4/15, ln-Koeffizienten 7/2 und 5/2, zusammen -1/(12 pi^2) mal -1/32).")


def ent_bos(nu):
    nu = nu[nu > 0.5 + 1e-12]
    return float(np.sum((nu + 0.5) * np.log(nu + 0.5) - (nu - 0.5) * np.log(nu - 0.5)))


def S_boson(M, ring):
    q = 2 * np.pi * (np.arange(ring) + 0.5) / ring
    w = np.sqrt(M * M + 4 * np.sin(q / 2) ** 2)
    cosm = np.cos(np.outer(np.arange(ring // 2), q))
    X, P = toeplitz(cosm @ (1 / (2 * w)) / ring), toeplitz(cosm @ (w / 2) / ring)
    Lc = np.linalg.cholesky(X)
    return ent_bos(np.sqrt(np.clip(np.linalg.eigvalsh(Lc.T @ P @ Lc), 0.25, None))) / 2


def ent_f(nu):
    nu = nu[(nu > 1e-14) & (nu < 1 - 1e-14)]
    return float(-np.sum(nu * np.log(nu) + (1 - nu) * np.log(1 - nu)))


def S_fermion(s, B, r, ring=200):
    """Kette d = (sin q, s, B + r(1 - cos q)) wie in Z15."""
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


MS = np.array([0.025, 0.035, 0.05, 0.07, 0.1, 0.14, 0.2])
t1 = time.time()
Mtab = np.concatenate([np.geomspace(0.012, 0.3, 40, endpoint=False), np.linspace(0.3, 3.05, 60)])
cs_b = CubicSpline(np.log(Mtab), [S_boson(M, max(300, 2 * int(20 / M))) for M in Mtab])
mutab = np.concatenate([np.geomspace(0.012, 0.3, 36, endpoint=False), np.linspace(0.3, np.sqrt(2.3), 40)])
cs_f = CubicSpline(np.log(mutab), [S_fermion(mu, 0.0, 0.0, max(200, 2 * int(12 / mu))) for mu in mutab])


def kappa_s_of(m):
    x, w = nodes(lo=min(1e-3, m / 10), n_gl=10)
    K1, K2 = np.meshgrid(x, x, indexing="ij")
    return float(np.sum(np.outer(w, w) * cs_b(np.log(np.sqrt(m * m + 4 * np.sin(K1 / 2) ** 2 + 4 * np.sin(K2 / 2) ** 2)))))


def kappa_naive_of(m):
    x, w = nodes(lo=min(1e-3, m / 10), n_gl=10, both=True)
    K1, K2 = np.meshgrid(x, x, indexing="ij")
    return 2 * float(np.sum(np.outer(w, w) * cs_f(np.log(np.sqrt(m * m + np.sin(K1) ** 2 + np.sin(K2) ** 2)))))


def fit_log(y, with_const, odd=False):
    cols = ([np.ones_like(MS)] if with_const else []) + [MS ** 2 * np.log(MS), MS ** 2, MS ** 4 * np.log(MS), MS ** 4]
    c = np.linalg.lstsq(np.stack(cols, 1), y, rcond=None)[0]
    return (c[0], c[1]) if with_const else (None, c[0])


ks_m = np.array([kappa_s_of(m) for m in MS])
kn_m = np.array([kappa_naive_of(m) for m in MS])
ss_m = np.array([16 * np.pi * bz(k2(lambda a, b, c, k: scal_bub(a, b, c, k, m, ONE, ONE)), lo=min(1e-4, m / 20)) for m in MS])
sn_m = np.array([16 * np.pi * bz(k2(lambda a, b, c, k: wil_bub(a, b, c, k, m, ONE, ONE, r=0.0, cW=0.0)),
                                 lo=min(1e-4, m / 20), both=True) for m in MS])
sig_s0 = 16 * np.pi * bz(k2(lambda a, b, c, k: scal_bub(a, b, c, k, 0.0, ONE, ONE)))
sig_n0 = 16 * np.pi * bz(k2(lambda a, b, c, k: wil_bub(a, b, c, k, 0.0, ONE, ONE, r=0.0, cW=0.0)), both=True)
k_s0_fit, A_ent_s = fit_log(ks_m, True)
k_n0_fit, A_ent_n = fit_log(kn_m, True)
_, A_sak_s = fit_log(ss_m - sig_s0, False)
_, A_sak_n = fit_log(sn_m - sig_n0, False)
tgt_s, tgt_n = 1 / (24 * np.pi), 8 / (12 * np.pi)
say(f"   Anpassung an m = {', '.join(str(m) for m in MS)} mit m^2 ln m, m^2, m^4 ln m, m^4 (Verschränkung mit Konstante):")
say(f"   Skalar:      Verschränkung {A_ent_s:.6f}, Sakharov {A_sak_s:.6f}, Kontinuum 1/(24 pi) = {tgt_s:.6f}")
say(f"   naives Dirac (8 Geschmäcker): Verschränkung {A_ent_n:.5f}, Sakharov {A_sak_n:.5f}, Kontinuum 8/(12 pi) = {tgt_n:.5f}")
say(f"   Konstanten: kappa_s(0) = {k_s0_fit:.6f} (Z15 {Z15['kappa_skalar_masselos']:.6f}), "
    f"kappa_naiv(0) = {k_n0_fit:.6f} (Z2 {Z2['kappa_naive_dirac']:.6f})")
check("Universell Skalar: m^2 ln m beider Seiten = 1/(24 pi) (je < 2 %)",
      abs(A_ent_s / tgt_s - 1) < 0.02 and abs(A_sak_s / tgt_s - 1) < 0.02)
check("Universell Dirac (naiv, 8 Geschmäcker): m^2 ln m beider Seiten = 8/(12 pi) (je < 3 %)",
      abs(A_ent_n / tgt_n - 1) < 0.03 and abs(A_sak_n / tgt_n - 1) < 0.03)
check("Reproduktion: kappa_s(0) wie Z15 und kappa_naiv(0) wie Z2 (< 2e-5 bzw. < 2e-4)",
      abs(k_s0_fit - Z15["kappa_skalar_masselos"]) < 2e-5 and abs(k_n0_fit - Z2["kappa_naive_dirac"]) < 2e-4)
say(f"   Identität im universellen Teil auf dem Gitter erfüllt: Skalar "
    f"{'ja' if CHECKS['Universell Skalar: m^2 ln m beider Seiten = 1/(24 pi) (je < 2 %)'] else 'NEIN'}, Dirac "
    f"{'ja' if CHECKS['Universell Dirac (naiv, 8 Geschmäcker): m^2 ln m beider Seiten = 8/(12 pi) (je < 3 %)'] else 'NEIN'}"
    f" (R_universell Skalar {A_ent_s / A_sak_s:.3f}, Dirac {A_ent_n / A_sak_n:.3f}).")
say(f"   Damit ist die Normierung sigma = 16 pi e_2 auf dem Gitter bestätigt ({time.time() - t1:.0f} s).")
say("   Der nicht-universelle Teil (Ordnung 1/a^2) ist damit nicht festgelegt; um ihn geht es unter D.\n")
ZAHLEN["kalibrierung_universell"] = {"skalar_verschraenkung": A_ent_s, "skalar_sakharov": A_sak_s, "ziel_skalar": tgt_s,
                                     "naivdirac_verschraenkung": A_ent_n, "naivdirac_sakharov": A_sak_n, "ziel_naivdirac": tgt_n,
                                     "kappa_s0_fit": k_s0_fit, "kappa_naiv0_fit": k_n0_fit}
ZAHLEN["kalibrierung_kontinuum"] = {"bildmethode_max_abw": max(img_s, img_f, img_v), "S2_a1": [hk_s, hk_f, hk_v], "R_kontinuum": cont_R}

# ---------------------------------------------------------------------------------------------------
# Teil 3: nicht-universell, je Feld und Regel
# ---------------------------------------------------------------------------------------------------
say("Teil 3: Nicht-universeller Teil: sigma = a^2/(4 G_ind) je Feld und Kopplungsregel")
t1 = time.time()


def sigma_table(bub, dia, m, both=False, n_gl=8, has_half_vertex=False):
    dI, dH, rho = bz(lambda a, b, c: dia(a, b, c, m), n_gl=n_gl, both=both)
    e0 = bz(lambda a, b, c: bub(a, b, c, 0.0, m, ONE, ONE), n_gl=n_gl, both=both)
    eb = {}
    for rule in ("K", "P", "D"):
        if rule == "P" and not has_half_vertex:
            eb["P"] = eb["K"]
            continue
        ai, ah = FACT[rule]
        eb[rule] = bz(k2(lambda a, b, c, k: bub(a, b, c, k, m, ai, ah)), n_gl=n_gl, both=both)
    sig = {"K": 16 * np.pi * eb["K"],
           "P_kov": 16 * np.pi * (eb["P"] - 0.25 * dH), "P_no": 16 * np.pi * eb["P"],
           "D_kov": 16 * np.pi * (eb["D"] - 0.25 * dI), "D_no": 16 * np.pi * eb["D"]}
    return {"sigma": {k_: float(v) for k_, v in sig.items()}, "e0_blase": float(e0), "kontakt_ganz": float(dI),
            "kontakt_halb": float(dH), "rho": float(rho), "e0_LV": float(e0 + dI + dH + rho / 4)}


RULES = ("K", "P_kov", "P_no", "D_kov", "D_no")
T = {}
T["skalar"] = sigma_table(scal_bub, scal_dia, 0.0)
T["skalar_m0.5"] = sigma_table(scal_bub, scal_dia, 0.5)
for cW in (0.0, 0.5, 1.0):
    T[f"wilson_cW{cW}"] = sigma_table(lambda a, b, c, k, m, ai, ah, cW=cW: wil_bub(a, b, c, k, m, ai, ah, 1.0, cW),
                                      lambda a, b, c, m, cW=cW: wil_dia(a, b, c, m, 1.0, cW), 0.0)
    T[f"wilson_cW{cW}_m0.5"] = sigma_table(lambda a, b, c, k, m, ai, ah, cW=cW: wil_bub(a, b, c, k, m, ai, ah, 1.0, cW),
                                           lambda a, b, c, m, cW=cW: wil_dia(a, b, c, m, 1.0, cW), 0.5)
T["photon"] = sigma_table(phot_bub, phot_dia, 0.0, n_gl=6, has_half_vertex=True)
ph8 = 16 * np.pi * bz(k2(lambda a, b, c, k: phot_bub(a, b, c, k, 0.0, ONE, ONE)), n_gl=8)
sw10 = 16 * np.pi * bz(k2(lambda a, b, c, k: wil_bub(a, b, c, k, 0.0, ONE, ONE)), n_gl=10, lo=1e-5)
conv_ok = abs(ph8 / T["photon"]["sigma"]["K"] - 1) < 1e-4 and abs(sw10 / T["wilson_cW1.0"]["sigma"]["K"] - 1) < 1e-5
say(f"   Konvergenz Gauß-Legendre: Photon K {T['photon']['sigma']['K']:.6f} (6) / {ph8:.6f} (8 Knoten je Abschnitt); "
    f"Wilson K {T['wilson_cW1.0']['sigma']['K']:.7f} / {sw10:.7f}: {check('Quadratur konvergiert (Photon 1e-4, Wilson 1e-5)', conv_ok)}")
sig_r = {r_: 16 * np.pi * bz(k2(lambda a, b, c, k, r_=r_: wil_bub(a, b, c, k, 0.0, ONE, ONE, r=r_, cW=1.0))) for r_ in (0.5, 2.0)}

kap = {"skalar": Z15["kappa_skalar_masselos"], "skalar_m0.5": Z15["kappa_skalar_massiv"]["0.5"],
       "wilson": Z15["kappa_dirac_wilson_r1"]["0.0"], "wilson_m0.5": Z15["kappa_dirac_wilson_r1"]["0.5"],
       "naiv": Z2["kappa_naive_dirac"], "photon_ohne_zentrum": Z15["kappa_photon_ohne_zentrum"],
       "photon_ez_dest": Z15["kappa_photon_ez_destillierbar"]}


# Wilson-kappa aus Z15 grob nachrechnen (48^2 Querimpulse, Z15-Kette)
def kappa_wilson(m, r, nk=48, ring=200):
    k = 2 * np.pi * (np.arange(nk // 2) + 0.5) / nk
    tot = wsum = 0.0
    for a in range(nk // 2):
        for b in range(a, nk // 2):
            w = 1.0 if a == b else 2.0
            tot += w * S_fermion(np.hypot(np.sin(k[a]), np.sin(k[b])), m + r * (2 - np.cos(k[a]) - np.cos(k[b])), r, ring)
            wsum += w
    return 2 * tot / wsum


kw48 = kappa_wilson(0.0, 1.0)
say(f"   kappa_Wilson(r = 1) aus Z15 {kap['wilson']:.5f}, grob nachgerechnet (48^2) {kw48:.5f}: "
    f"{check('kappa_Wilson aus Z15 nachgerechnet (< 1e-3)', abs(kw48 - kap['wilson']) < 1e-3)}")
s_naiv_K = sig_n0

say(f"\n   sigma = a^2/(4 G_ind) in nat je Gitterfläche (Spalten: Regel und Kontaktterm):")
say("   " + " " * 26 + "".join(f"{r_:>10s}" for r_ in RULES) + "     kappa")
rows = [("Skalar m=0", "skalar", kap["skalar"]), ("Skalar m=0,5", "skalar_m0.5", kap["skalar_m0.5"]),
        ("Wilson-Dirac c_W=1", "wilson_cW1.0", kap["wilson"]), ("Wilson-Dirac c_W=1/2", "wilson_cW0.5", kap["wilson"]),
        ("Wilson-Dirac c_W=0", "wilson_cW0.0", kap["wilson"]), ("Wilson c_W=1 m=0,5", "wilson_cW1.0_m0.5", kap["wilson_m0.5"]),
        ("Photon", "photon", kap["photon_ohne_zentrum"])]
for lab, key, kv in rows:
    say(f"   {lab:26s}" + "".join(f"{T[key]['sigma'][r_]:10.4f}" for r_ in RULES) + f"{kv:10.4f}")
say(f"   naives Dirac (8 Geschmäcker), K: sigma = {s_naiv_K:.4f}, kappa = {kap['naiv']:.4f}")
say(f"   Wilson r = 0,5 / 2 (c_W = 1, K): sigma = {sig_r[0.5]:.4f} / {sig_r[2.0]:.4f}; kappa (Z15) "
    f"{Z15['kappa_dirac_wilson_r']['0.5']['0.0']:.4f} / {Z15['kappa_dirac_wilson_r']['2.0']['0.0']:.4f}")
say(f"   (Rechenzeit Teil 3 {time.time() - t1:.0f} s)")

R_field = {}
for lab, key, kv in rows[:5]:
    R_field[lab] = {r_: kv / T[key]["sigma"][r_] for r_ in RULES}
R_field["naives Dirac"] = {"K": kap["naiv"] / s_naiv_K}
R_field["Wilson r=0,5"] = {"K": Z15["kappa_dirac_wilson_r"]["0.5"]["0.0"] / sig_r[0.5]}
R_field["Wilson r=2"] = {"K": Z15["kappa_dirac_wilson_r"]["2.0"]["0.0"] / sig_r[2.0]}
say("\n   R = kappa/sigma je Feld:")
for lab in ("Skalar m=0", "Wilson-Dirac c_W=1", "Wilson-Dirac c_W=1/2", "Wilson-Dirac c_W=0"):
    say(f"   {lab:26s}" + "".join(f"{R_field[lab][r_]:10.3f}" for r_ in RULES))
say(f"   naives Dirac K: {R_field['naives Dirac']['K']:.3f}; Wilson r = 0,5 / 2 (K): "
    f"{R_field['Wilson r=0,5']['K']:.3f} / {R_field['Wilson r=2']['K']:.3f}")

# Bruch der Identität in der kanonischen Kopplung, als xi ausgedrückt
phi2 = bz(lambda a, b, c: np.nan_to_num(1 / (2 * np.sqrt(4 * np.sin(a / 2) ** 2 + 4 * np.sin(b / 2) ** 2 + 4 * np.sin(c / 2) ** 2))))
xi_star = (T["skalar"]["sigma"]["K"] - kap["skalar"]) / (2 * np.pi * phi2)
say(f"\n   Gitterbruch der Identität in der kanonischen Kopplung K (Kontinuum: R = 1):")
say(f"   Skalar R = {R_field['Skalar m=0']['K']:.3f}, naives Dirac R = {R_field['naives Dirac']['K']:.3f}, "
    f"Wilson R = {R_field['Wilson-Dirac c_W=1']['K']:.3f} (c_W = 1) bis {R_field['Wilson-Dirac c_W=0']['K']:.3f} (c_W = 0).")
say(f"   Ein Term xi R phi^2 verschiebt sigma um -2 pi xi <phi^2>, <phi^2> = {phi2:.5f} (Kontinuum mit Eigenzeitregler:")
say(f"   dieselbe Formel gibt die bekannte (1/6 - xi)-Abhängigkeit). R_Skalar = 1 verlangt xi = {xi_star:.4f} statt 0.")
check("Gitter bricht die Identität in K für Skalar und naives Dirac (R zwischen 1,2 und 2): negativer Befund festgestellt",
      1.2 < R_field["Skalar m=0"]["K"] < 2 and 1.2 < R_field["naives Dirac"]["K"] < 2)

# Kopplungsmehrdeutigkeit
spread = {lab: (min(T[key]["sigma"].values()), max(T[key]["sigma"].values())) for lab, key, _ in rows}
say("\n   Kopplungsmehrdeutigkeit: Der k^0-Teil der Antwort (Vakuum-Schermodul) ist Ordnung 1/a^4, der gesuchte k^2-Teil")
say("   Ordnung 1/a^2. Jede Kopplungsänderung der Ordnung a^2 d^2 (Ort der Metrik, Interpolation, Kontaktterme)")
say("   verschiebt sigma daher um Größen derselben Ordnung wie sigma selbst:")
for lab, key, kv in rows:
    d = T[key]
    say(f"   {lab:22s} e0_Blase = {d['e0_blase']:8.4f}, Kontakt ganz/halb = {d['kontakt_ganz']:8.4f}/{d['kontakt_halb']:8.4f},"
        f" e_2(K) = {d['sigma']['K'] / (16 * np.pi):.5f}")
wil_all = [v for c_ in ("0.0", "0.5", "1.0") for v in T[f"wilson_cW{c_}"]["sigma"].values()]
k_edge055 = kap["photon_ez_dest"] + Z15["kante_je_eichboson"]["0.55"]["innovation"]
brackets = {"Skalar": min(T["skalar"]["sigma"].values()) < kap["skalar"] < max(T["skalar"]["sigma"].values()),
            "Wilson": min(wil_all) < kap["wilson"] < max(wil_all),
            "Photon ohne Rand": min(T["photon"]["sigma"].values()) < kap["photon_ohne_zentrum"] < max(T["photon"]["sigma"].values()),
            "Photon mit Rand g=0,55": min(T["photon"]["sigma"].values()) < k_edge055 < max(T["photon"]["sigma"].values())}
dev_vs_spread = {"Skalar": (abs(T["skalar"]["sigma"]["K"] - kap["skalar"]), spread["Skalar m=0"]),
                 "Wilson c_W=1": (abs(T["wilson_cW1.0"]["sigma"]["K"] - kap["wilson"]), spread["Wilson-Dirac c_W=1"]),
                 "Photon mit Rand": (abs(T["photon"]["sigma"]["K"] - k_edge055), spread["Photon"])}
say(f"   Spannweite sigma je Feld: " + "; ".join(f"{lab} {a_:.3f} bis {b_:.3f}" for lab, (a_, b_) in spread.items() if "m=0,5" not in lab))
say("   Umschließen die natürlichen Regeln kappa? " + ", ".join(f"{k_} {'ja' if v else 'nein'}" for k_, v in brackets.items()))
say(f"   R = 1 liegt je Feld innerhalb der Kopplungsmehrdeutigkeit: "
    f"{befund('R = 1 liegt je Feld innerhalb der Kopplungsmehrdeutigkeit', all(brackets.values()))}")
say("   Abweichung |sigma_K - kappa| gegen Spannweite der Regeln: "
    + ", ".join(f"{k_} {d_:.3f} gegen {b_[1] - b_[0]:.3f}" for k_, (d_, b_) in dev_vs_spread.items()))
say(f"   Mehrdeutigkeit größer als die Abweichung von R = 1 in K, für jedes Feld: "
    f"{befund('Kopplungsmehrdeutigkeit übersteigt die Abweichung von R = 1', all(b_[1] - b_[0] > d_ for d_, b_ in dev_vs_spread.values()))}")
lv = {key: T[key]["e0_LV"] for key in ("skalar", "wilson_cW1.0", "photon")}
say(f"   Gradientenfreier TT-Term nach Abzug der kosmologischen Konstanten (e0_LV, kovariantes Vakuum: 0): "
    + ", ".join(f"{k_} {v:.3f}" for k_, v in lv.items()) + ".")
say(f"   Das Hamiltonsche Gittervakuum ist nicht lorentzinvariant (Druck != -rho); als induzierte Gravitation gelesen")
say(f"   ist e0_LV ein TT-Massenterm an der Gitterskala: {befund('Gitter erzeugt gradientenfreien TT-Term (Gravitonmasse an der Gitterskala)', all(abs(v) > 0.1 for v in lv.values()))}\n")
ratio_pg = T["photon"]["sigma"]["K"] / T["skalar"]["sigma"]["K"]
# Kontinuum, scharfer räumlicher Impulsschnitt: Winkelmittel der Blase je |p| (k^0 und k^2-Dichte bei |p| = 1)
_x, _w = leggauss(40)
_ph = np.linspace(0, 2 * np.pi, 81)[:-1]
_CT, _PH = np.meshgrid(_x, _ph, indexing="ij")
_W = np.outer(_w, np.ones_like(_ph)) / (2 * len(_ph))
_n = np.stack([np.sqrt(1 - _CT ** 2) * np.cos(_PH), np.sqrt(1 - _CT ** 2) * np.sin(_PH), _CT])


def cont_bub(field, k):
    s_, sp_ = _n.copy(), _n.copy()
    sp_[2] = sp_[2] + k
    w_, wp_ = np.linalg.norm(s_, axis=0), np.linalg.norm(sp_, axis=0)
    if field == "skalar":
        return -(1 / 16) * (-s_[0] ** 2 + s_[1] ** 2) ** 2 / (w_ * wp_ * (w_ + wp_))
    sq = np.sqrt(w_ * wp_)[..., None, None]
    Z = np.swapaxes(_cross(s_), -1, -2) @ DD @ _cross(sp_) / sq - DD * sq
    sh, shp = np.moveaxis(s_ / w_, 0, -1), np.moveaxis(sp_ / wp_, 0, -1)
    PT, PTp = np.eye(3) - sh[..., :, None] * sh[..., None, :], np.eye(3) - shp[..., :, None] * shp[..., None, :]
    return -(1 / 16) * np.einsum("...ij,...jk,...lk,...il->...", PT, Z, PTp, Z) / (w_ + wp_)


cont = {f: (float(np.sum(_W * cont_bub(f, 0.0))),
            float(np.sum(_W * (cont_bub(f, 1e-3) + cont_bub(f, -1e-3) - 2 * cont_bub(f, 0.0)))) / 2e-6) for f in ("skalar", "photon")}
say(f"   Kontinuum mit scharfem räumlichem Impulsschnitt (Winkelmittel bei |p| = 1): Skalar k^0 {cont['skalar'][0]:.6f}"
    f" (-1/120), k^2 {cont['skalar'][1]:.7f} (11/3360 = {11 / 3360:.7f}); Photon k^0 {cont['photon'][0]:.4f}, k^2 {cont['photon'][1]:.7f} (3/280).")
check("Kontinuum-TT-Blase: Skalar -1/120 und 11/3360 analytisch, Photon k^2 = 3/280",
      abs(cont["skalar"][0] + 1 / 120) < 1e-9 and abs(cont["skalar"][1] - 11 / 3360) < 1e-6 and abs(cont["photon"][1] - 3 / 280) < 1e-6)
say(f"   Reglerabhängigkeit des Lambda^2-Koeffizienten, Photon/Skalar: Wärmekern -4, scharfer räumlicher Schnitt")
say(f"   {cont['photon'][1] / cont['skalar'][1]:.3f} (= 36/11), Gitter K {ratio_pg:.1f}.")

# ---------------------------------------------------------------------------------------------------
# Teil 4: Eichfelder mit und ohne Rand
# ---------------------------------------------------------------------------------------------------
say("\nTeil 4: Eichfelder mit und ohne Randterm")
C0_geo = Z15["randfluss"]["innovationsvarianz_e"]


def H_dg(sig, d=0.0):
    K = int(np.ceil(12 * sig)) + 5
    n = np.arange(-K, K + 1)
    w = np.exp(-(n - d) ** 2 / (2 * sig * sig))
    p = w / w.sum()
    p = p[p > 0]
    return float(-np.sum(p * np.log(p)))


def edge(g):
    return float(np.mean([H_dg(np.sqrt(C0_geo) / g, d) for d in np.linspace(0, 1, 64, endpoint=False)]))


edge_ok = abs(edge(0.55) - Z15["kante_je_eichboson"]["0.55"]["innovation"]) < 1e-9
say(f"   Randterm (Shannon des ganzzahligen Randflusses) aus Z15 nachgerechnet, g = 0,55: {edge(0.55):.4f}: "
    f"{check('Randterm s(g) wie Z15 (g = 0,55)', edge_ok)}")
k_dest = kap["photon_ez_dest"]
say("   Sakharov-Seite: Bei schwacher Kopplung ist das kompakte Photon Gaußsch plus O(g^2); g skaliert aus der")
say("   Einschleifenantwort heraus. Der Randterm wächst dagegen wie ln(1/g): R_Eich(g) kann nur bei einem g gleich 1 sein.")
say("   Photon je Eichboson, R = kappa/sigma:")
say("   " + " " * 30 + "".join(f"{r_:>10s}" for r_ in RULES))
cases = [("ohne Zentrum (2 kappa_s)", kap["photon_ohne_zentrum"]), ("EZ destillierbar", k_dest),
         ("EZ + Rand g=0,55", k_dest + edge(0.55)), ("EZ + Rand g=1", k_dest + edge(1.0))]
R_gauge = {}
for lab, kv in cases:
    R_gauge[lab] = {r_: kv / T["photon"]["sigma"][r_] for r_ in RULES}
    say(f"   {lab:30s}" + "".join(f"{R_gauge[lab][r_]:10.3f}" for r_ in RULES))
noedge_far = all(abs(R_gauge[lab][r_]) < 0.1 for lab in ("ohne Zentrum (2 kappa_s)", "EZ destillierbar") for r_ in RULES)
say(f"   Ohne Randterm verfehlen Eichfelder R = 1 in allen fünf Regeln um mehr als Faktor 10: "
    f"{check('Eichfeld ohne Randterm: |R| < 0,1 in allen Regeln (negativer Befund)', noedge_far)}")
say("   Unter D ist der Randterm (Zentrumswahl) daher notwendig, wie Donnelly-Wall im Kontinuum.")
gstar = {}
for r_ in RULES:
    sg = T["photon"]["sigma"][r_]
    need = sg - k_dest
    if need <= 0.05:
        gstar[r_] = None
        continue
    gstar[r_] = brentq(lambda g: edge(g) - need, 0.05, 20.0) if need < edge(0.05) else None
say("   g, bei dem R_Eich = 1 (Rand = sigma - kappa_dest): "
    + ", ".join(f"{r_}: {'keins (sigma < kappa_dest)' if v is None else f'{v:.3f}'}" for r_, v in gstar.items()))
say(f"   Vorzeichen: Der Gitterrand ist positiv (Shannon), Kabats Kontaktterm im Wärmekern negativ. Das Gitterphoton")
say(f"   induziert in K, P_no und D_no positives 1/G (anders als -4 im Wärmekern), dort kann der Rand ausgleichen;")
say(f"   in P_kov und D_kov ist sigma_Photon < 0 und kein g hilft.\n")

# ---------------------------------------------------------------------------------------------------
# Teil 5: voller H6-Inhalt
# ---------------------------------------------------------------------------------------------------
say("Teil 5: Voller H6-Inhalt (Schema wie Z15: Wilson r = 1, schwer bei m = 0,5, 12 Eichbosonen, 4 Higgs)")
say("   Leicht: 3 x 16 mit nu_R = 24 Dirac masselos. Schwer (Z23-Eltern, getrennt): 16 Majorana-Weyl = 8 Dirac")
say("   bei m = 0,5 und Phi als 2 reelle Skalare bei m = 0,5. Eich: 12 Kopien, Rand mit n_eff = 12 (oder 8, nichtabelsch).")


def kappa_parts(g=0.55, neff=12):
    return {"leicht": 24 * kap["wilson"], "schwer": 8 * kap["wilson_m0.5"] + 2 * kap["skalar_m0.5"],
            "higgs": 4 * kap["skalar"], "eich": 12 * k_dest + neff * edge(g)}


def sigma_parts(rule, cW):
    return {"leicht": 24 * T[f"wilson_cW{cW}"]["sigma"][rule],
            "schwer": 8 * T[f"wilson_cW{cW}_m0.5"]["sigma"][rule] + 2 * T["skalar_m0.5"]["sigma"][rule],
            "higgs": 4 * T["skalar"]["sigma"][rule], "eich": 12 * T["photon"]["sigma"][rule]}


kp = kappa_parts()
k_tot = sum(kp.values())
say(f"   kappa (H6-Schema, g = 0,55): " + ", ".join(f"{k_} {v:.3f}" for k_, v in kp.items()) + f"; gesamt {k_tot:.3f}"
    f" (Z15: {Z15['h6_schema']['kappa']:.3f})")
check("kappa_gesamt des H6-Schemas wie Z15 (< 1e-6)", abs(k_tot - Z15["h6_schema"]["kappa"]) < 1e-6)
say("\n   sigma_gesamt und R_gesamt (g = 0,55; n_eff = 12 / 8):")
say(f"   {'Regel':8s}{'c_W':>6s}{'leicht':>9s}{'schwer':>9s}{'Higgs':>8s}{'Eich':>9s}{'sigma':>9s}{'R(12)':>8s}{'R(8)':>8s}"
    f"{'g*(12)':>9s}{'g*(8)':>8s}{'a/l_P|R=1':>11s}")
TOT = []
k8 = sum(kappa_parts(neff=8).values())
EDGE_MAX = edge(0.02)
for rule in RULES:
    for cW in (1.0, 0.5, 0.0):
        sp = sigma_parts(rule, cW)
        st = sum(sp.values())
        gs = []
        for neff in (12, 8):
            base = sum(kappa_parts(neff=neff).values()) - neff * edge(0.55)
            need = (st - base) / neff
            gs.append(brentq(lambda g: edge(g) - need, 0.02, 50.0) if 0.05 < need < EDGE_MAX else (0.0 if need >= EDGE_MAX else None))
        row = {"regel": rule, "c_W": cW, "sigma_teile": sp, "sigma": st, "R12": k_tot / st, "R8": k8 / st,
               "g_stern_12": gs[0], "g_stern_8": gs[1], "a_lP_bei_R1": 2 * np.sqrt(st) if st > 0 else None}
        TOT.append(row)
        fmt = lambda v: "   -  " if v is None else ("<0.02" if v == 0.0 else f"{v:6.3f}")
        say(f"   {rule:8s}{cW:6.1f}" + "".join(f"{sp[k_]:9.3f}" if k_ != "higgs" else f"{sp[k_]:8.3f}" for k_ in sp)
            + f"{st:9.3f}{row['R12']:8.3f}{row['R8']:8.3f}{fmt(gs[0]):>9s}{fmt(gs[1]):>8s}"
            + (f"{row['a_lP_bei_R1']:11.2f}" if st > 0 else "         -"))
canon = [r for r in TOT if r["regel"] == "K"]
c1 = [r for r in canon if r["c_W"] == 1.0][0]
say(f"\n   Schwere Eltern getrennt (K, c_W = 1): kappa {kp['schwer']:.3f}, sigma {c1['sigma_teile']['schwer']:.3f}, "
    f"R_schwer = {kp['schwer'] / c1['sigma_teile']['schwer']:.3f}; Anteil an sigma {c1['sigma_teile']['schwer'] / c1['sigma']:.1%}.")
Rk = [r["R12"] for r in canon] + [r["R8"] for r in canon]
say(f"   Kanonische Kopplung K: R_gesamt = {min(Rk):.2f} bis {max(Rk):.2f} (c_W und n_eff variiert), bei g = 0,55.")
gk = [r["g_stern_12"] for r in canon] + [r["g_stern_8"] for r in canon]
gk = [g for g in gk if g]
ak = [r["a_lP_bei_R1"] for r in canon]
say(f"   R = 1 verlangt in K g(1/a) = {min(gk):.2f} bis {max(gk):.2f} (SM-Wert dort nach Z15 etwa 0,5 bis 0,6) und gibt dann")
say(f"   a/l_P = 2 sqrt(sigma) = {min(ak):.2f} bis {max(ak):.2f}, unabhängig von g und Zentrumswahl.")
allR = [r["R12"] for r in TOT] + [r["R8"] for r in TOT]
pos = [r["a_lP_bei_R1"] for r in TOT if r["a_lP_bei_R1"] is not None]
say(f"   Über alle Regeln: R_gesamt von {min(allR):.2f} bis {max(allR):.2f}; a/l_P bei R = 1 von {min(pos):.1f} bis {max(pos):.1f}"
    f" (wo sigma > 0); in {sum(r['sigma'] <= 0 for r in TOT)} von {len(TOT)} Kombinationen ist sigma <= 0.")
nar = [r for r in TOT if r["regel"] in ("K", "P_no", "D_no")]
nR = [r["R12"] for r in nar] + [r["R8"] for r in nar]
nA = [r["a_lP_bei_R1"] for r in nar]
say(f"   Engeres Band mit lokal abgezogenem Vakuum (K, P_no, D_no): R_gesamt {min(nR):.2f} bis {max(nR):.2f},"
    f" a/l_P bei R = 1 von {min(nA):.1f} bis {max(nA):.1f}.")
say("   Die kov-Varianten zeigen, was geschieht, wenn die Vakuumenergie der Gitterobjekte an Metrikgradienten koppelt.")
k_band = any(r["R12"] > 1 for r in TOT) and any(0 < r["R12"] < 1 for r in TOT)
check("Gesamtbilanz berechnet: R_gesamt über die natürlichen Regeln oberhalb und unterhalb von 1", k_band)
befund("R_gesamt = 1 in der kanonischen Kopplung K bei g = 0,55 (auf 20 %)", any(abs(r) < 1.2 and r > 0.8 for r in Rk))
befund("R_gesamt != 1 robust über alle natürlichen Kopplungen", not k_band)
befund("a/l_P aus Sakharov = Jacobson als Vorhersage (Spanne unter Faktor 1,5 über alle Regeln)", max(pos) / min(pos) < 1.5)
befund("a/l_P bei R = 1 in K enger als die Z15-Spanne 3,3 bis 14,3", min(ak) > 3.3 and max(ak) < 14.3)
R_ohne_schwer = (k_tot - kp["schwer"]) / (c1["sigma"] - c1["sigma_teile"]["schwer"])
say(f"   Ohne schwere Eltern (K, c_W = 1): R = {R_ohne_schwer:.3f} statt {c1['R12']:.3f}.")
befund("schwere Eltern verschieben R_gesamt in K um mehr als 10 %", abs(c1["R12"] / R_ohne_schwer - 1) > 0.1)

# ---------------------------------------------------------------------------------------------------
# Einordnung
# ---------------------------------------------------------------------------------------------------
say("\nEinordnung")
say("   Prämisse D: Weder widerlegt noch bestätigt. Der universelle Teil der Identität gilt auf dem Gitter (Teil 2).")
say("   Der Teil, der a/l_P festlegt, ist auf der Sakharov-Seite nicht durch den Feldinhalt bestimmt: Er hängt an der")
say("   Kopplung des Gitters an die Metrik in Ordnung a^2 d^2, weil das Gittervakuum einen Schermodul der Ordnung 1/a^4")
say(f"   trägt. Mit der kanonischen Kopplung K übersteigt die Verschränkung die induzierte Flächendichte (R_ges {min(Rk):.2f} bis {max(Rk):.2f});")
say("   die übrigen natürlichen Regeln überdecken R = 1 von beiden Seiten. Konsistenz von D mit H6-Inhalt verlangt eine")
say("   Zusatzannahme: entweder eine Kopplungsregel, die R = 1 ergibt, oder einen nackten Term (kappa - sigma)/a^2 aus")
say("   nicht-materiellen Freiheitsgraden, die H6 nicht hat.")
say("   a/l_P: bleibt unter D die Definition 2 sqrt(kappa). Fordert man zusätzlich R = 1, folgt a/l_P = 2 sqrt(sigma) ohne")
say(f"   Rand- und g-Abhängigkeit, aber mit der Kopplungsregel als neuer Eingabe: in K {min(ak):.1f} bis {max(ak):.1f},")
say(f"   mit lokal abgezogenem Vakuum {min(nA):.1f} bis {max(nA):.1f}, mit Vakuumkontakttermen bis {max(pos):.0f}.")
say("   T7: Die fehlende Festlegung ist eine Gitterfassung der Diffeomorphismus-Wardidentität (wie T_mu nu an h_mu nu")
say("   koppelt, bis O(a^2 d^2)), dieselbe Struktur, die ein masseloser Spin-2-Pol braucht. Dazu kommt der")
say("   gradientenfreie TT-Term des Gittervakuums: Ein über Sakharov dynamisches Metrikfeld bekäme eine Masse an der")
say("   Gitterskala, die wie Lambda abgestimmt werden müsste. Unter D (kein Graviton) stört das nicht, zeigt aber, dass")
say("   der Schritt von D zu einem Graviton (T7) ohne neues Prinzip nicht aus dem Feldinhalt folgt.")
say("   Grenzen: nur statischer Koeffizient von sqrt(h) R^(3); TT-Mode '+' mit k entlang einer Gitterachse (die kubische")
say("   Gruppe erlaubt einen zweiten Koeffizienten für 'x'); freie Felder (Einschleifen); nichtabelsche Eichfelder als")
say("   U(1)-Kopien; Rand in der Innovationsnäherung aus Z15; Kopplungsfamilie mit drei Ortsregeln, zwei")
say("   Kontakttermen und c_W in {0, 1/2, 1}, nicht erschöpfend.")

def de(x, nd=2):
    return f"{x:.{nd}f}".replace(".", ",")


urteil = (
    f"Sakharov gegen Jacobson ist auf dem H6-Gitter nicht entscheidbar, mit benennbarem Grund. Der universelle Teil der "
    f"Identität gilt im selben Regler (m^2 ln m: Skalar {de(A_ent_s, 4)}/{de(A_sak_s, 4)}, naives Dirac {de(A_ent_n, 3)}/"
    f"{de(A_sak_n, 3)} gegen 1/(24 pi) bzw. 8/(12 pi)). Der nicht-universelle Flächenkoeffizient a^2/(4 G_ind), der a/l_P "
    f"festlegt, hängt aber um mehr als eine Größenordnung und im Vorzeichen davon ab, wie das Hüpfen in Ordnung a^2 d^2 an "
    f"die Metrik koppelt, weil das Gittervakuum einen Schermodul der Ordnung 1/a^4 trägt. In der kanonischen Kopplung bricht "
    f"das Gitter die Identität: R = {de(R_field['Skalar m=0']['K'])} (Skalar), {de(R_field['naives Dirac']['K'])} (naives "
    f"Dirac), {de(R_field['Wilson-Dirac c_W=1']['K'])} bis {de(R_field['Wilson-Dirac c_W=0']['K'])} (Wilson, je nach c_W), für "
    f"den H6-Inhalt R_ges = {de(min(Rk))} bis {de(max(Rk))} bei g = 0,55. Eichfelder ohne Randterm verfehlen R = 1 in jeder "
    f"Regel um mehr als Faktor 10; unter D ist der Gitterrand also notwendig. Über die natürlichen Kopplungen liegt R = 1 "
    f"innerhalb der Mehrdeutigkeit. D ist daher nicht widerlegt, a/l_P wird aber keine Vorhersage: Fordert man R = 1, gilt "
    f"a/l_P = 2 sqrt(sigma) = {de(min(ak), 1)} bis {de(max(ak), 1)} in der kanonischen Kopplung (dort mit g(1/a) = "
    f"{de(min(gk))} bis {de(max(gk))} statt etwa 0,55), {de(min(nA), 1)} bis {de(max(nA), 1)} mit lokal abgezogenem Vakuum "
    f"und bis {max(pos):.0f}, wenn die Vakuumenergie an Metrikgradienten koppelt. Es fehlt eine Gitter-Wardidentität für die "
    f"Metrikkopplung; dieselbe Struktur braucht T7 für einen Spin-2-Pol.")
say("\nUrteil:")
for line in textwrap.wrap(urteil, 112):
    say("   " + line)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)}). Laufzeit {time.time() - T0:.0f} s.")

ZAHLEN.update({
    "sigma_je_feld": {k_: v["sigma"] for k_, v in T.items()},
    "vakuum_je_feld": {k_: {"e0_blase": v["e0_blase"], "kontakt_ganz": v["kontakt_ganz"], "kontakt_halb": v["kontakt_halb"],
                            "rho": v["rho"], "e0_LV": v["e0_LV"]} for k_, v in T.items()},
    "sigma_naiv_dirac_K": s_naiv_K, "sigma_wilson_r_K": {str(k_): v for k_, v in sig_r.items()},
    "kappa_je_feld": kap, "kappa_wilson_48_nachgerechnet": kw48,
    "R_je_feld": R_field, "R_eich": R_gauge, "g_stern_eich_je_regel": gstar,
    "xi_stern_skalar": xi_star, "phi2_skalar": phi2, "photon_zu_skalar_K": ratio_pg,
    "kappa_h6_schema": k_tot, "kappa_h6_neff8": k8, "gesamt": TOT,
    "kanonisch_K": {"R_gesamt": [min(Rk), max(Rk)], "g_stern": [min(gk), max(gk)], "a_lP_bei_R1": [min(ak), max(ak)]},
    "alle_regeln": {"R_gesamt": [min(allR), max(allR)], "a_lP_bei_R1": [min(pos), max(pos)]},
    "vakuum_lokal_abgezogen": {"R_gesamt": [min(nR), max(nR)], "a_lP_bei_R1": [min(nA), max(nA)]},
    "kontinuum_scharfer_schnitt": {"skalar_k0_k2": list(cont["skalar"]), "photon_k0_k2": list(cont["photon"])},
    "laufzeit_s": time.time() - T0,
})
(HERE / "z27_ergebnis.json").write_text(json.dumps({"checks": CHECKS, "befunde": BEFUNDE, "urteil": urteil, "zahlen": ZAHLEN},
                                                   indent=1, ensure_ascii=False, default=float))
assert all(CHECKS.values())
