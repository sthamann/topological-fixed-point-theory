"""Z13: Gravitative Anomalie der Naht und ihr Zufluss aus dem Volumen (Folgeaufgabe 3 zu H6h8).

Gittermodell aus Z7 (verdrillter Defekt eines Paarungsmassenfelds, Defektlinie entlang z, k_z erhalten).

A. Casimir-Impuls. Auf einem Ring der Länge L_z trägt eine chirale CFT den Grundzustandsimpuls
   P = (2 pi / L_z)(h - c/24) für Rechtsläufer (Linksläufer mit umgekehrtem Vorzeichen). Antiperiodisch
   (NS) ist h = 0, periodisch (R) ist h das Twistgewicht. Gemessen wird der Impuls der besetzten
   BdG-Zustände, projiziert auf Kern, Ring (lückenhaftes Volumen) und äußeren Rand:
     - Gitterimpulsdichte sum_z (1/2i)(c_z^+ c_{z+1} - h.c.) = sum_k sin k n_k,
     - Impulspolarisation <T_A> (Tu, Zhang, Qi 2013): Verschiebung nur der Region A um eine Zelle.
   Aus dem Koeffizienten von 1/L_z folgt die chirale Zentralladung c_A = c(+z) - c(-z) je Region.
B. Thermische Kontrolle. Der Energiestrom entlang z einer chiralen CFT ist J_E = (pi/12) c T^2.
C. Volumen. Eine Wand, über die die Massenphase einmal um 2 pi windet (ohne Naht, periodisch), ist ein
   2D System. Seine BdG-Chernzahl C ist die thermische Hallantwort, c_Wand = C/2. Das ist der
   Koeffizient des gravitativen theta-Terms im Volumen, also der Zufluss je Windung (Callan-Harvey).
D. Bezug zu I4 = (c/24) p1 mit c = 8 für 16 Paare.
"""
import json
import os
import time
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import get_context
from pathlib import Path

import numpy as np
from scipy.linalg import eigh, svdvals

HERE = Path(__file__).resolve().parent
_src = (HERE / "z7_ursprung_naht.py").read_text()
_ns = {"__file__": str(HERE / "z7_ursprung_naht.py")}
exec(_src[:_src.index("L = 14")], _ns)
build = _ns["build"]
MASS = _ns["MASS"]

L = 12                        # Querschnitt L x L wie in Z9
R_KERN, R_RING = 3.0, 5.0     # Kern r <= 3, Ring 3 < r <= 5, Rand r > 5
M_MAIN = 1.0                  # Masse; Kern-Rand-Hybridisierung bei k_z = 0 dann < 1e-4 (bei m = 0.6: 8e-3)
CONFIGS = {"verdrillt": (1, 0, 0.25), "gewöhnlich": (1, 1, 0.0), "ohne Windung": (0, 0, 0.25)}
LZ_NS = [8, 10, 12, 14, 16, 20, 24, 28, 32, 40, 48, 64, 80, 96]
LZ_R = [9, 11, 13, 15, 17, 21, 25, 29, 33, 41, 49, 65, 81, 97]
LZ_FIT_MIN = 16
TS = np.array([0.015, 0.02, 0.025, 0.03, 0.04, 0.05])
KTH = np.linspace(0.0, 0.64, 161)

_W = {}


# ---------------------------------------------------------------- Arbeiterprozesse (Defektgeometrie)
def _init(Lc, cfg, m):
    MASS.clear()
    MASS["m"] = m
    H0, (rad, ns) = build(Lc, 0.0, *cfg)
    Hh, _ = build(Lc, np.pi / 2, *cfg)
    Hp, _ = build(Lc, np.pi, *cfg)
    MASS.clear()
    A, C = (H0 + Hp) / 2, (H0 - Hp) / 2
    B = Hh - A
    site_rows = lambda mask: np.tile(np.repeat(mask, ns * 2), 2)
    regions = {"Kern": rad <= R_KERN, "Ring": (rad > R_KERN) & (rad <= R_RING), "Rand": rad > R_RING,
               "Innen": rad <= R_RING}
    _W.update(A=A, B=B, C=C, rad=rad, ns=ns, L=Lc, rows={k: site_rows(v) for k, v in regions.items()})


def _H(k):
    return _W["A"] + np.sin(k) * _W["B"] + np.cos(k) * _W["C"]


def _point(k):
    """Besetzte BdG-Zustände bei k > 0: Überschussbesetzung je Gitterplatz und Phasen für <T_A>."""
    E, V = eigh(_H(k), driver="evr")
    Vo = V[:, E < 0]
    Lc, ns = _W["L"], _W["ns"]
    dens = (np.abs(Vo) ** 2).sum(axis=1)
    f_site = dens.reshape(2, Lc * Lc, ns * 2).sum(axis=(0, 2)) - 2 * ns
    targ = {}
    for name, rows in _W["rows"].items():
        nu = svdvals(Vo[rows]) ** 2
        n_a = rows.sum() // 2
        targ[name] = float(np.sum(np.angle(1 - nu + nu * np.exp(1j * k))) - k * n_a)
    return f_site, targ, float(np.min(np.abs(E)))


def _thermal_point(k):
    """Energien und regionsprojizierte Geschwindigkeiten <n|Pi_A dH/dk|n> der nullnahen Zustände."""
    E, V = eigh(_H(k), driver="evr")
    sel = np.abs(E) < 0.75
    Vs, Es = V[:, sel], E[sel]
    dV = (np.cos(k) * _W["B"] - np.sin(k) * _W["C"]) @ Vs
    out = {}
    for name, rows in _W["rows"].items():
        out[name] = np.real(np.einsum("ij,ij->j", Vs[rows].conj(), dV[rows]))
    return Es, out


def run_defect(cfg, m, lz_ns, lz_r, thermal, workers):
    ks_ns = sorted({round(2 * np.pi * (n + 0.5) / Lz, 12) for Lz in lz_ns for n in range(Lz // 2)})
    ks_r = sorted({round(2 * np.pi * n / Lz, 12) for Lz in lz_r for n in range(1, (Lz - 1) // 2 + 1)})
    ks = sorted(set(ks_ns) | set(ks_r))
    with ProcessPoolExecutor(max_workers=workers, mp_context=get_context("spawn"), initializer=_init,
                             initargs=(L, cfg, m)) as ex:
        res = dict(zip(ks, ex.map(_point, ks, chunksize=4)))
        th = list(ex.map(_thermal_point, KTH, chunksize=4)) if thermal else None
    return res, th


# ---------------------------------------------------------------- Auswertung
def grid(Lz, sector):
    if sector == "NS":
        return [round(2 * np.pi * (n + 0.5) / Lz, 12) for n in range(Lz // 2)]
    return [round(2 * np.pi * n / Lz, 12) for n in range(1, (Lz - 1) // 2 + 1)]


def series(res, mask, lzs, sector, kind):
    out = []
    for Lz in lzs:
        tot = 0.0
        for k in grid(Lz, sector):
            f = res[k][0][mask].sum()
            tot += {"sin": np.sin(k) * f, "naiv": k * f}[kind]
        out.append(tot)
    return np.array(out)


def series_T(res, name, lzs, sector):
    return np.array([sum(res[k][1][name] for k in grid(Lz, sector)) for Lz in lzs])


def lsq(Ls, P, powers):
    X = np.array([[float(Lz) ** p for p in powers] for Lz in Ls])
    return np.linalg.lstsq(X, P, rcond=None)[0]


def kappa(Ls, P):
    """Koeffizient von 1/L_z mit Unsicherheit aus Variation von Fitmodell und Fitfenster."""
    Ls, P = np.asarray(Ls), np.asarray(P)
    sel = Ls >= LZ_FIT_MIN
    main = lsq(Ls[sel], P[sel], (1, -1, -3))[1]
    alts = [lsq(Ls[sel], P[sel], (1, -1))[1], lsq(Ls[sel], P[sel], (1, -1, -3, -5))[1],
            lsq(Ls[Ls >= 24], P[Ls >= 24], (1, -1, -3))[1], lsq(Ls, P, (1, -1, -3, -5))[1]]
    const = lsq(Ls[sel], P[sel], (1, 0, -1, -3))[1]
    return main, max(abs(a - main) for a in alts) + 1e-6, const


def c_from(res, mask_or_name, kind="sin"):
    if kind == "T":
        Pns, Pr = series_T(res, mask_or_name, LZ_NS, "NS"), series_T(res, mask_or_name, LZ_R, "R")
    else:
        Pns, Pr = series(res, mask_or_name, LZ_NS, "NS", kind), series(res, mask_or_name, LZ_R, "R", kind)
    kn, en, bn = kappa(LZ_NS, Pns)
    kr, er, br = kappa(LZ_R, Pr)
    c, dc = -24 * kn / (2 * np.pi), 24 * en / (2 * np.pi)
    hr, dh = kr / (2 * np.pi) + c / 24, er / (2 * np.pi) + dc / 24
    return {"c": float(c), "dc": float(dc), "h_R": float(hr), "dh_R": float(dh), "kappa_NS": float(kn),
            "kappa_R": float(kr), "R_zu_NS": float(kr / kn) if abs(kn) > 1e-9 else None,
            "konst_NS": float(bn), "konst_R": float(br)}


def jump_of(res, mask):
    """Sprung f_A(0-) - f_A(0+) = -2 f_A(0+) der projizierten Besetzung (Gewicht der chiralen Moden in A),
    unabhängig vom Casimir-Fit aus den zwei kleinsten k linear nach k -> 0+ extrapoliert."""
    k1, k2 = sorted(res)[:2]
    f1, f2 = -2 * res[k1][0][mask].sum(), -2 * res[k2][0][mask].sum()
    return float(f1 + (f1 - f2) * k1 / (k2 - k1))


def thermal_c(th, name):
    dk = KTH[1] - KTH[0]
    wts = np.full(len(KTH), dk)
    wts[0] = dk / 2
    wts[-1] = dk / 2
    J = []
    for T in TS:
        tot = 0.0
        for w, (E, v) in zip(wts, th):
            occ = 0.5 * (1 - np.tanh(E / (2 * T)))
            tot += w * np.sum(E * v[name] * (occ - (E < 0)))
        J.append(tot / (2 * np.pi))       # Faktor 2 (k < 0 gespiegelt) mal 1/2 (BdG)
    J = np.array(J)
    X = np.stack([np.ones_like(TS), TS ** 2, TS ** 4], axis=1)
    coef = np.linalg.lstsq(X, J, rcond=None)[0]
    X2 = np.stack([np.ones_like(TS[1:]), TS[1:] ** 2], axis=1)
    coef2 = np.linalg.lstsq(X2, J[1:], rcond=None)[0]
    c = coef[1] * 12 / np.pi
    return float(c), float(abs(coef2[1] * 12 / np.pi - c)) + 1e-6


# ---------------------------------------------------------------- Volumen: Wand ohne Naht
SX = np.array([[0, 1], [1, 0]], complex)
SY = np.array([[0, -1j], [1j, 0]])
SZ = np.array([[1, 0], [0, -1]], complex)
ISY = 1j * SY


def slab_H(ky, kz, Nx, n_plus, n_minus, flux, m, rw=1.0, wall=None):
    """x periodisch (Nx Plätze), k_y und k_z erhalten. Massenphase theta(x) windet einmal um 2 pi:
    gleichförmig oder als schmale Wand der Breite `wall` um theta = pi (Minimum aus Z9)."""
    x = np.arange(Nx)
    if wall is None:
        theta = 2 * np.pi * x / Nx
    else:
        theta = np.pi + np.pi * (1 + np.tanh((x - Nx / 2 + 0.5) / wall))
    ns = 2
    n = Nx * ns * 2
    idx = lambda xx, a, s: (xx * ns + a) * 2 + s
    hk = np.zeros((n, n), complex)
    hm = np.zeros((n, n), complex)
    D = np.zeros((n, n), complex)
    Aphase = 2 * np.pi * flux / Nx
    wind = (n_plus, n_minus)
    for xx in range(Nx):
        for a in range(ns):
            b0 = idx(xx, a, 0)
            hk[b0:b0 + 2, b0:b0 + 2] += np.sin(ky) * SY + np.sin(kz) * SZ
            hm[b0:b0 + 2, b0:b0 + 2] += -np.sin(ky) * SY - np.sin(kz) * SZ
            D[b0:b0 + 2, b0:b0 + 2] += m * np.exp(1j * wind[a] * theta[xx]) * ISY
        Wos = rw * (3 - np.cos(ky) - np.cos(kz))
        p0, q0 = idx(xx, 0, 0), idx(xx, 1, 0)
        D[p0:p0 + 2, q0:q0 + 2] += Wos * ISY
        D[q0:q0 + 2, p0:p0 + 2] += Wos * ISY
        x2 = (xx + 1) % Nx
        for a in range(ns):
            q = 1 if a == 0 else -1
            val = SX / (2j) * np.exp(-1j * q * Aphase)
            i0, j0 = idx(xx, a, 0), idx(x2, a, 0)
            for H_ in (hk, hm):
                H_[i0:i0 + 2, j0:j0 + 2] += val
                H_[j0:j0 + 2, i0:i0 + 2] += val.conj().T
        for (aa, bb, q) in ((0, 1, +1), (1, 0, -1)):
            pw = -rw / 2 * np.exp(-1j * q * Aphase)
            i0, j0 = idx(xx, aa, 0), idx(x2, bb, 0)
            D[i0:i0 + 2, j0:j0 + 2] += pw * ISY
            D[j0:j0 + 2, i0:i0 + 2] += pw * ISY          # D^T = -D, weil ISY^T = -ISY
    return np.block([[hk, D], [D.conj().T, -hm.T]])


def wall_chern(Nx, n_plus, n_minus, flux, m, Nk=24, wall=None):
    ks = 2 * np.pi * np.arange(Nk) / Nk
    U = {}
    gap = np.inf
    for i, ky in enumerate(ks):
        for j, kz in enumerate(ks):
            E, V = eigh(slab_H(ky, kz, Nx, n_plus, n_minus, flux, m, wall=wall))
            U[i, j] = V[:, E < 0]
            gap = min(gap, np.min(np.abs(E)))
    F = 0.0
    for i in range(Nk):
        for j in range(Nk):
            a, b = U[i, j], U[(i + 1) % Nk, j]
            c_, d = U[(i + 1) % Nk, (j + 1) % Nk], U[i, (j + 1) % Nk]
            loop = (np.linalg.det(a.conj().T @ b) * np.linalg.det(b.conj().T @ c_)
                    * np.linalg.det(c_.conj().T @ d) * np.linalg.det(d.conj().T @ a))
            F += np.angle(loop)
    return F / (2 * np.pi), float(gap)


# ---------------------------------------------------------------- Hauptprogramm
def main():
    t0 = time.time()
    lines = []

    def say(s=""):
        print(s, flush=True)
        lines.append(s)

    checks, befunde, data = {}, {}, {}
    workers = max(2, min(12, (os.cpu_count() or 4) - 2))

    say("Z13: Gravitative Anomalie der Naht und ihr Zufluss\n")
    say(f"Gitter aus Z7: Querschnitt {L} x {L}, m = {M_MAIN}, r_w = 1, k_z erhalten.")
    say(f"Regionen: Kern r <= {R_KERN}, Ring {R_KERN} < r <= {R_RING} (lückenhaftes Volumen), Rand r > {R_RING}.")
    say("Innen = Kern + Ring.\n")

    say("A. Casimir-Impuls auf einem Ring der Länge L_z\n")
    say("Für einen Rechtsläufer gilt P = (2 pi / L_z)(h - c/24). Antiperiodisch (NS): h = 0.")
    say("Periodisch (R): h = 1/16 je Majorana-Weyl Mode. Gemessen wird P_A(L_z) je Region A.")
    say("Die BdG-Verdopplung ist berücksichtigt: P = sum_{k>0} sin k (Psi_k^+ Psi_k - N).")
    say(f"NS: L_z = {LZ_NS}; R: L_z = {LZ_R}; Fit P = a L_z + kappa/L_z + d/L_z^3 ab L_z >= {LZ_FIT_MIN}.")
    say("Unsicherheit: größte Abweichung unter anderen Fitmodellen und Fitfenstern.\n")

    results = {}
    for name, cfg in CONFIGS.items():
        t1 = time.time()
        res, th = run_defect(cfg, M_MAIN, LZ_NS, LZ_R, True, workers)
        rad = np.hypot(*np.meshgrid(np.arange(L) - (L - 1) / 2, np.arange(L) - (L - 1) / 2, indexing="ij")).ravel()
        masks = {"Kern": rad <= R_KERN, "Ring": (rad > R_KERN) & (rad <= R_RING), "Rand": rad > R_RING,
                 "Innen": rad <= R_RING, "Alle": np.ones_like(rad, bool)}
        out = {"sin": {r: c_from(res, mk, "sin") for r, mk in masks.items()},
               "T": {r: c_from(res, r, "T") for r in ("Kern", "Ring", "Rand", "Innen")},
               "naiv": {r: c_from(res, mk, "naiv") for r, mk in masks.items()},
               "thermisch": {r: thermal_c(th, r) for r in ("Kern", "Ring", "Rand", "Innen")}}
        jump = {r: jump_of(res, mk) for r, mk in masks.items()}
        k1 = min(res)
        out["sprung"] = jump
        out["min_gap"] = float(min(v[2] for v in res.values()))
        # radiales Profil c(r) für Scheiben r <= R
        out["profil"] = {f"{R:.1f}": c_from(res, rad <= R + 1e-9, "sin")["c"] for R in (1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0)}
        results[name] = out
        say(f"{name} (n+, n-, Fluss) = {cfg}   [{len(res)} k-Punkte, {time.time() - t1:.0f} s]")
        say(f"   kleinste |E| auf den Gittern: {out['min_gap']:.4f}; k_min = {k1:.4f}")
        for r in ("Kern", "Ring", "Rand", "Innen"):
            s, T, th_c = out["sin"][r], out["T"][r], out["thermisch"][r]
            say(f"   {r:5s}: c = {s['c']:+.4f} ± {s['dc']:.4f} (Impulsdichte), {T['c']:+.4f} ± {T['dc']:.4f} "
                f"(<T_A>), {th_c[0]:+.4f} ± {th_c[1]:.4f} (J_E); Sprung/2 = {jump[r] / 2:+.4f}")
            say(f"          h_R = {s['h_R']:+.4f} ± {s['dh_R']:.4f}; kappa_R/kappa_NS = "
                f"{s['R_zu_NS'] if s['R_zu_NS'] is None else round(s['R_zu_NS'], 3)}; Konstante im NS-Fit {s['konst_NS']:+.1e}")
        say(f"   Profil c(r <= R): " + ", ".join(f"R = {k}: {v:+.3f}" for k, v in out["profil"].items()))
        say(f"   Summe aller Regionen (Impulsdichte): c = {out['sin']['Alle']['c']:+.1e} (folgt aus der"
            " Teilchen-Loch-Symmetrie)\n")

    tw, ov, nw = results["verdrillt"], results["gewöhnlich"], results["ohne Windung"]
    checks["verdrillt: c_Kern = 1/2 (Casimir NS)"] = abs(tw["sin"]["Kern"]["c"] - 0.5) < 0.02
    checks["verdrillt: c_Innen = 1/2 (Kern + Ring)"] = abs(tw["sin"]["Innen"]["c"] - 0.5) < 0.01
    checks["verdrillt: Ring trägt nur den Ausläufer (|c| < 0.02)"] = abs(tw["sin"]["Ring"]["c"]) < 0.02
    checks["verdrillt: Twistgewicht h_R = 1/16"] = abs(tw["sin"]["Innen"]["h_R"] - 1 / 16) < 0.003
    checks["verdrillt: R/NS = -2 (h_R = c/8)"] = abs(tw["sin"]["Innen"]["R_zu_NS"] + 2) < 0.05
    checks["gewöhnlicher Wirbel: c_Innen = 1, h_R = 1/8"] = (abs(ov["sin"]["Innen"]["c"] - 1) < 0.02
                                                         and abs(ov["sin"]["Innen"]["h_R"] - 1 / 8) < 0.006)
    checks["ohne Windung: c = 0 in Kern, Ring und Rand"] = all(abs(nw["sin"][r]["c"]) < 0.01
                                                             for r in ("Kern", "Ring", "Rand"))
    checks["Casimir-Koeffizient = halber Sprung der Kernbesetzung"] = all(
        abs(res_["sin"]["Innen"]["c"] - res_["sprung"]["Innen"] / 2) < 0.01 for res_ in results.values())
    checks["Impulspolarisation <T_A> gibt dasselbe c"] = all(
        abs(res_["T"][r]["c"] - res_["sin"][r]["c"]) < 0.01 for res_ in results.values() for r in ("Kern", "Innen"))
    checks["Energiestrom (pi/12) c T^2 gibt dasselbe c"] = all(
        abs(res_["thermisch"]["Innen"][0] - res_["sin"]["Innen"]["c"]) < 0.02 for res_ in results.values())
    checks["kein konstanter Term im Casimir-Fit"] = all(abs(res_["sin"]["Innen"]["konst_NS"]) < 1e-4
                                                     for res_ in results.values())

    say("Lesart von A.")
    say(f"   Der verdrillte Defekt trägt c = {tw['sin']['Innen']['c']:.4f} ± {tw['sin']['Innen']['dc']:.4f}, also eine")
    say("   Majorana-Weyl Mode (c = 1/2). Der Rand trägt das Negative. Der Ring trägt nur den Ausläufer der")
    say(f"   Kernmode ({tw['sin']['Ring']['c']:.4f}; beim gewöhnlichen Wirbel {ov['sin']['Ring']['c']:.4f}), keine eigene Chiralität.")
    say(f"   Im R-Sektor ist h_R = {tw['sin']['Innen']['h_R']:.4f} (erwartet 1/16 = 0.0625). Bei 16 Paaren ist das")
    say("   h_R = 1, das Twistgewicht aus Z10 (Punkt 5), jetzt aus dem Gitter.")
    say("   Der gewöhnliche Wirbel gibt c = 1 und h_R = 1/8. Ohne Windung ist die projizierte Besetzung für")
    say("   jedes k exakt k-unabhängig (Symmetrie des reellen Massenfelds); dort ist alles null.")
    say("   Der Koeffizient ist der halbe Sprung der Kernbesetzung bei k_z = 0. Für freie Fermionen misst die")
    say("   Casimir-Rechnung also die Netto-Chiralität; neu sind die Normierung 1/24 und das Twistgewicht.\n")

    say("Kontrolle: naive Gewichtung mit k statt sin k (k im Intervall (-pi, pi])")
    for name in CONFIGS:
        s = results[name]["naiv"]
        say(f"   {name:13s}: c_Kern = {s['Kern']['c']:+.4f}, c_Ring = {s['Ring']['c']:+.4f}, c_Rand = {s['Rand']['c']:+.4f}")
    naive_dev = max(abs(results[n]["naiv"][r]["c"] - results[n]["sin"][r]["c"]) for n in CONFIGS
                    for r in ("Kern", "Ring", "Rand"))
    say(f"   Größte Abweichung von der Impulsdichte: {naive_dev:.4f}. Der Sprung der Gewichtung bei k = pi gibt")
    say("   einen nicht universellen 1/L_z-Beitrag; sin k und <T_A> sind am Zonenrand glatt.\n")
    befunde["naive k-Gewichtung weicht ab (Zonenrand)"] = bool(naive_dev > 0.02)

    say("Robustheit (verdrillt): andere Masse")
    rob = {}
    for m in (0.8, 0.6):
        res, _ = run_defect(CONFIGS["verdrillt"], m, LZ_NS, LZ_R, False, workers)
        rad = np.hypot(*np.meshgrid(np.arange(L) - (L - 1) / 2, np.arange(L) - (L - 1) / 2, indexing="ij")).ravel()
        r_ = c_from(res, rad <= R_RING, "sin")
        r_["sprung_halb"] = jump_of(res, rad <= R_RING) / 2
        prof = {R: c_from(res, rad <= R + 1e-9, "sin")["c"] for R in (3.0, 3.5, 4.0, 4.5, 5.0)}
        r_["plateau"] = max(prof.values())
        rob[m] = r_
        say(f"   m = {m}: c_Innen = {r_['c']:+.4f} ± {r_['dc']:.4f}, Sprung/2 = {r_['sprung_halb']:+.4f}, "
            f"h_R = {r_['h_R']:+.4f} ± {r_['dh_R']:.4f}, R/NS = {r_['R_zu_NS']:.3f}, Konstante {r_['konst_NS']:+.1e}")
        say("            Profil c(r <= R): " + ", ".join(f"R = {R}: {v:+.3f}" for R, v in prof.items()))
    checks["robust bei m = 0.8 (c = Sprung/2, R/NS = -2, Plateau 1/2)"] = (
        abs(rob[0.8]["c"] - rob[0.8]["sprung_halb"]) < 0.01 and abs(rob[0.8]["R_zu_NS"] + 2) < 0.05
        and abs(rob[0.8]["plateau"] - 0.5) < 0.03)
    befunde["m = 0.6 (Z7-Standard): Casimir-Fit gestört (R/NS weicht um > 0,05 von -2 ab)"] = bool(
        abs(rob[0.6]["R_zu_NS"] + 2) > 0.05)
    befunde["m = 0.6: Plateau von c(r) verfehlt 1/2 um mehr als 0,02"] = bool(abs(rob[0.6]["plateau"] - 0.5) > 0.02)
    say("   Bei kleinerem m ist die Kernmode breiter (Länge ~ 1/m) und überlappt im Querschnitt 12 mit der")
    say("   Randmode; c(r) erreicht das Plateau 1/2 nicht mehr ganz. Bei m = 0.6 kommt die Hybridisierung von")
    say("   Kern- und Randmode bei k_z = 0 (8e-3) hinzu; sie verschmiert den Sprung, sobald 2 pi / L_z klein wird,")
    say("   und R/NS weicht von -2 ab. Deshalb m = 1 als Hauptrechnung; die Aussage ist c = Sprung/2.\n")

    say("B. Volumen: Wand ohne Naht (Zuflusskoeffizient)\n")
    say("Die Massenphase windet entlang einer periodischen Richtung x einmal um 2 pi, y und z sind periodisch.")
    say("Es gibt keine Naht und keinen Rand. Die BdG-Chernzahl C in der (k_y, k_z)-Ebene ist die thermische")
    say("Hallantwort des Volumens auf den Phasengradienten: kappa_xy = (pi/6) T c_Wand mit c_Wand = C/2.")
    say("Nach Callan-Harvey ist das der Zufluss je Windung auf die Naht, an der die Wand endet.")
    wall = {}
    wall_cases = [("verdrillt (1, 0), Viertelfluss", (1, 0, 0.25), None),
                  ("verdrillt (1, 0), ohne Fluss", (1, 0, 0.0), None),
                  ("verdrillt (1, 0), schmale Wand um theta = pi", (1, 0, 0.25), 1.2),
                  ("gewöhnlich (1, 1)", (1, 1, 0.0), None),
                  ("ohne Windung (0, 0)", (0, 0, 0.25), None),
                  ("doppelt (2, 0)", (2, 0, 0.25), None),
                  ("umgekehrt (-1, 0)", (-1, 0, -0.25), None)]
    Nx = 16
    for label, (a, b, fl), wd in wall_cases:
        C, gap = wall_chern(Nx, a, b, fl, M_MAIN, wall=wd)
        wall[label] = {"C": float(C), "gap": gap}
        say(f"   {label:45s}: C = {C:+.4f}, Lücke {gap:.3f}, c_Wand = {C / 2:+.3f}")
    Cs = {k: round(v["C"]) for k, v in wall.items()}
    checks["Wand: C ganzzahlig und Lücke offen"] = all(abs(v["C"] - round(v["C"])) < 1e-6 and v["gap"] > 0.05
                                                     for v in wall.values())
    checks["Wand: |C| = n+ + n- (1, 1, 1, 2, 0, 2, 1)"] = [abs(Cs[k]) for k, _, _ in wall_cases] == [1, 1, 1, 2, 0, 2, 1]
    checks["Wand: Vorzeichen folgt der Windungsrichtung"] = (
        np.sign(Cs["verdrillt (1, 0), Viertelfluss"]) == -np.sign(Cs["umgekehrt (-1, 0)"])
        and np.sign(Cs["gewöhnlich (1, 1)"]) == np.sign(Cs["verdrillt (1, 0), Viertelfluss"]))
    c_wall = abs(wall["verdrillt (1, 0), Viertelfluss"]["C"]) / 2
    checks["Zufluss je Windung = c der Naht (verdrillt und gewöhnlich)"] = (
        abs(c_wall - tw["sin"]["Innen"]["c"]) < 0.01
        and abs(abs(wall["gewöhnlich (1, 1)"]["C"]) / 2 - ov["sin"]["Innen"]["c"]) < 0.02)
    say(f"   Der Zufluss je Windung ist c_Wand = {c_wall:.3f}, die Naht trägt c = {tw['sin']['Innen']['c']:.4f}.")
    say("   Der Viertelfluss (J) ändert C nicht; die schmale Wand um das Minimum theta = pi aus Z9 auch nicht.")
    say("   Die Wand ist der Verzweigungsschnitt aus Z9. Sie trägt eine 2D chirale Phase mit c = 1/2 je Paar.")
    say("   Ihre beiden Enden sind die Naht und der äußere Rand, mit c = +1/2 und -1/2 wie in A.\n")

    say("C. Bezug zu I4 = (c/24) p1\n")
    c_pair = tw["sin"]["Innen"]["c"]
    say(f"   Je Paar: c = {c_pair:.4f} (Naht) und c_Wand = {c_wall:.3f} (Volumen). I4 je Paar = p1/48.")
    say(f"   16 entkoppelte Paare (Blöcke additiv): c = 16 x {c_pair:.4f} = {16 * c_pair:.3f}, Wand C = 16, c_Wand = 8.")
    say(f"   Twistgewicht 16 x {tw['sin']['Innen']['h_R']:.4f} = {16 * tw['sin']['Innen']['h_R']:.3f}.")
    say("   I4 = (8/24) p1 = p1/3 auf der Naht wird vom gravitativen theta-Term des Volumens ausgeglichen:")
    say("   S = (c_Wand/24) int theta p1. Weil d(d theta) = 2 pi delta_Naht, gibt die Abstiegsgleichung für eine")
    say("   Windung genau (c_Wand/24) p1 auf der Naht mit umgekehrtem Vorzeichen. Der Koeffizient c_Wand ist hier")
    say("   gerechnet, nicht gesetzt. theta -> theta + 2 pi ändert S um 2 pi c_Wand sigma/8; das ist die Phase aus")
    say("   Z5, trivial auf allen orientierten 4-Mannigfaltigkeiten genau für c_Wand in 8Z.")
    say("   Der Schnitt trägt damit c = 8. Lesart (nicht gerechnet): 16 Kopien der c = 1/2 Phase (nu = 16 im")
    say("   16-fach-Weg) sind die E8-Phase mit trivialem Fermion, also die Phase, die Z5 für F angenommen hat.\n")

    say("Was nicht gezeigt ist")
    say("   - Wechselwirkungen und ein dynamisches Massenfeld; das Feld ist vorgegeben.")
    say("   - Der 4D Koeffizient von theta p1 selbst. Gerechnet ist seine 2D Reduktion (Chernzahl der Wand),")
    say("     der Rest ist die Abstiegsgleichung.")
    say("   - Ein ortsaufgelöster Energie-Impuls-Fluss in gekrümmtem Hintergrund (Gitter mit z-abhängiger")
    say("     Geschwindigkeit). Das braucht ein Gitter ohne k_z-Erhaltung und ist hier zu groß.")
    say("   - Gemischte Anomalien (SO(16) mit Gravitation).")
    befunde["4D Koeffizient von theta p1 nicht direkt gerechnet (nur 2D Reduktion)"] = True
    befunde["ortsaufgelöster gravitativer Fluss in Krümmung nicht gerechnet"] = True

    urteil = ("Gravitative Anomalie der Naht: im freien Gittermodell gerechnet (c = 1/2 und h_R = 1/16 je Paar). "
              "Zufluss: Koeffizient aus dem Volumen gerechnet und gleich c der Naht; der Fluss selbst in "
              "gekrümmtem Hintergrund ist nicht gerechnet. Also teilweise gerechnet, kein Widerspruch.")
    say("\nUrteil")
    say("   Gravitative Anomalie der Naht: im freien Gittermodell gerechnet, c = 1/2 und h_R = 1/16 je Paar.")
    say("   Zufluss: Der Koeffizient kommt aus dem Volumen (Wand ohne Naht) und ist gleich c der Naht.")
    say("   Der ortsaufgelöste Fluss in gekrümmtem Hintergrund ist nicht gerechnet. Also teilweise gerechnet.")

    say("\nZusammenfassung")
    for n, c in checks.items():
        say(f"   {'bestanden' if c else 'NICHT bestanden'}: {n}")
    for n, c in befunde.items():
        say(f"   Befund: {n}: {'ja' if c else 'nein'}")
    say(f"\nLaufzeit {time.time() - t0:.0f} s")

    def slim(res_):
        return {"sin": res_["sin"], "T": res_["T"], "naiv": {r: res_["naiv"][r]["c"] for r in res_["naiv"]},
                "thermisch": {r: {"c": v[0], "dc": v[1]} for r, v in res_["thermisch"].items()},
                "sprung": res_["sprung"], "profil": res_["profil"], "min_gap": res_["min_gap"]}

    data = {"checks": {k: bool(v) for k, v in checks.items()}, "befunde": befunde, "urteil": urteil,
            "parameter": {"L": L, "m": M_MAIN, "R_Kern": R_KERN, "R_Ring": R_RING, "LZ_NS": LZ_NS, "LZ_R": LZ_R,
                          "LZ_fit_min": LZ_FIT_MIN, "T": TS.tolist(), "Nx_Wand": Nx},
            "defekt": {k: slim(v) for k, v in results.items()},
            "robust": {str(k): v for k, v in rob.items()},
            "wand": wall,
            "I4": {"c_je_Paar": c_pair, "c_Wand_je_Windung": c_wall, "c_16": 16 * c_pair,
                   "h_R_16": 16 * tw["sin"]["Innen"]["h_R"]}}
    (HERE / "z13_ergebnis.json").write_text(json.dumps(data, ensure_ascii=False, indent=1))
    (HERE / "ausgabe_z13.txt").write_text("\n".join(lines) + "\n")
    assert all(checks.values())


if __name__ == "__main__":
    main()
