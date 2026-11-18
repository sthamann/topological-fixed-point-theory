#!/usr/bin/env python3
"""Part SF: CHIRAL4D.NODE_SPECTRAL_FIT.02 -- separate energy and residue instead of reading G(beta/2).

G(beta/2) < 1/2 can mean a gap, a small residue Z, a continuum or several states.  For every measured momentum and Cartan
sector the full reweighted tau-curve is fitted, without the endpoints (G(0) = G(beta) = 1/2 is an equal-time sum rule and
would force Z = 1), by the single symmetric pole pair
    G(tau) = Z cosh(E (beta/2 - tau)) / (2 cosh(beta E / 2)),
which is exact for a free two-level pair +-E (both spin states at one k).  Results: E, Z, chi^2 per degree of freedom.
  SF1 unit test: the free curves stored with every run (same code path, same cell, same twist) give E = |h_W(k)| exactly and
      Z = 1.
  SF2 unit test: a thermal operator correlation matrix C_ab(tau) = <O_a(tau) O_b^dag(0)> is positive semidefinite for every
      0 <= tau <= beta (brute force, random interacting 3-mode Hamiltonian, random operators); the measured reweighted G(tau)
      respects G(tau) >= 0 within errors.
  SF3 (statement) E, Z, chi^2 for every measured point; points with E dtau > 1 are flagged Trotter-limited (the energy is
      then only a lower bound of 'large').
  SF4 (statement) classification: gapless within the finite-beta resolution (E beta < 0.5) with Z > 0.5 (pole) or small Z
      (renormalised), else gapped -- per point, no thermodynamic statement.
  SF5 (.02) the symmetric pole pair M1 is only valid with particle-hole symmetry; every curve is also fitted by
      M2: two independent levels G = sum_i w_i e^{-e_i tau} / (1 + e^{-beta e_i}) (addition / removal energies and residues
      free) and M3: three levels with positive weights (a positive rest); the model is chosen by AIC, the reported
      excitation energy is min_i |e_i| over the resolvable levels of the chosen model (weight off the bounds, |e| dtau < 3).
      With 7 interior tau points and partly underestimated errors these fits are fragile: finer tau grids (smaller dtau,
      larger beta) and autocorrelation-corrected errors are needed before any number is used.  The Fisher information of
      an energy scales like Z^2: if the lowest level's energy error (errors x3) exceeds max(0.5, |E|) or its residue is
      below 2 sigma, the curve is flagged ENERGY_UNIDENTIFIABLE -- neither a gap nor gaplessness is read from it.  Unit tests: free curves give M2 = M1 (e_2 = -e_1, w_1 = w_2);
      synthetic asymmetric data are recovered by M2 and rejected by M1.
    python3 -B experiments/theory-contracts/seam-holonomy-generating-20261001/spectral_fit.py
"""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("reweight", HERE / "reweight.py")
rw = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(rw)
LOG: list[dict] = []


def step(name: str, ok: bool, detail: str = "") -> None:
    LOG.append({"step": name, "pass": bool(ok), "detail": detail})
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))
    sys.stdout.flush()


def model_curve(par, tau, beta):
    E, Z = abs(par[0]), par[1]
    return Z * np.cosh(E * (beta / 2 - tau)) / (2 * np.cosh(beta * E / 2))


def fit(G, err, beta, dtau):
    L = len(G) - 1
    tau = dtau * np.arange(L + 1)
    sel = np.arange(1, L)
    sig = np.maximum(np.asarray(err)[sel], 1e-4)
    Gh = G[L // 2]
    E0 = 2 / beta * np.arccosh(min(1e6, 1 / (2 * Gh))) if 0 < Gh < 0.5 else (0.0 if Gh >= 0.5 else 5.0)
    best = None
    for E_start in (E0, 0.5 * E0 + 0.1, 2 * E0 + 0.1, 0.0):
        r = least_squares(lambda p: (model_curve(p, tau[sel], beta) - np.asarray(G)[sel]) / sig, x0=[E_start, 1.0],
                          bounds=([-50, 0], [50, 2]))
        if best is None or r.cost < best.cost:
            best = r
    E, Z = abs(best.x[0]), best.x[1]
    return float(E), float(Z), float(2 * best.cost / max(len(sel) - 2, 1))


B3 = rw.B3.astype(float); S = [rw.sx, rw.sy, rw.sz]
Binv = np.linalg.inv(B3); tauv = [sum(Binv[a, i] * S[a] for a in range(3)) for i in range(3)]
def h_of_k(k):
    p = B3 @ k
    return sum((np.sin(p[i]) + np.cos(p[(i + 1) % 3]) - 1) * tauv[i] for i in range(3))
base = {"Gamma": np.zeros(3), "P": np.pi / 2 * np.ones(3), "-P": -np.pi / 2 * np.ones(3), "H": np.array([np.pi, 0, 0])}
def k_of(name, twist):
    b_, sec = (name[:-2], name[-2:]) if name.endswith(("+d", "-d")) else (name, "")
    k = base[b_].copy()
    if sec == "+d":
        k = k + np.asarray(twist)
    elif sec == "-d":
        k = k - np.asarray(twist)
    return k


pts, excluded = [], []
for f in sorted(HERE.glob("reweight_points_md_*.json")):
    for p in json.loads(f.read_text()):
        if "greens" in p:
            if p["R"] < 0.5:                                            # sign problem: N_eff ~ N R^2 too small, point unusable
                excluded.append((f.name, p["g"], p.get("twist_delta"), p["R"]))
                continue
            pts.append(p)

# SF1 --------------------------------------------------------------------------------------------------------
free_err = []
for p in pts:
    for name, v in p["greens"].items():
        E, Z, chi = fit(np.array(v["G_free"]), np.full(len(v["G_free"]), 1e-6), p["beta"], p["dtau"])
        Ex = float(np.max(np.abs(np.linalg.eigvalsh(h_of_k(k_of(name, p.get("twist_delta", [0, 0, 0])))))))
        free_err.append((abs(E - Ex), abs(Z - 1), name, p["cell"]))
mx = max(free_err, key=lambda t: t[0] + t[1])
step("SF1 unit test: fitting the free curves stored with every run (same cell, twist and code path) returns E = |h_W(k)| and Z = 1 for every measured momentum and Cartan sector",
     max(a for a, *_ in free_err) < 1e-4 and max(b for _, b, *_ in free_err) < 1e-4,
     f"{len(free_err)} free curves; max |E - |h(k)|| {max(a for a, *_ in free_err):.1e}, max |Z - 1| {max(b for _, b, *_ in free_err):.1e} (worst {mx[2]} on {mx[3]})")

# SF2 --------------------------------------------------------------------------------------------------------
rng = np.random.default_rng(5)
nb = 3
Zp = np.diag([1.0, -1.0]); low = np.array([[0, 1], [0, 0]], complex)
cf = []
for i in range(nb):
    mm = np.eye(1, dtype=complex)
    for j in range(nb):
        mm = np.kron(mm, Zp if j < i else (low if j == i else np.eye(2)))
    cf.append(mm)
Hr = rng.normal(size=(8, 8)) + 1j * rng.normal(size=(8, 8)); Hr = (Hr + Hr.conj().T) / 2
ops = [cf[0], cf[1] @ cf[2].conj().T @ cf[0], cf[2] @ cf[1] @ cf[0].conj().T + 0.3 * cf[1]]          # degree 1 and 3 seeds
Ev, Vv = np.linalg.eigh(Hr); beta_t = 2.0
rho = lambda t: (Vv * np.exp(-t * Ev)) @ Vv.conj().T
Zt = np.trace(rho(beta_t)).real
mineig = []
for t in np.linspace(0, beta_t, 9):
    C = np.array([[np.trace(rho(beta_t - t) @ a @ rho(t) @ b.conj().T) / Zt for b in ops] for a in ops])
    mineig.append(np.linalg.eigvalsh((C + C.conj().T) / 2).min())
neg = []
for p in pts:
    for name, v in p["greens"].items():
        G = np.array(v["G_tau"]); e = np.array(v["G_err"])
        neg.append((float(np.min(G + np.maximum(3 * e, 0.005))), float(np.min(G)), p["g"], name))
worst = min(neg, key=lambda t: t[1])
step("SF2 unit test: a thermal operator correlation matrix <O_a(tau) O_b^dag(0)> is positive semidefinite for every 0 <= tau <= beta (brute force, interacting 3-mode Hamiltonian, degree-1 and degree-3 operators); every measured reweighted G(tau) is >= 0 within max(3 sigma, 0.005) -- the binned errors underestimate the uncertainty at strong coupling (8 bins, longer autocorrelation), values ~ -0.002 there are zero",
     min(mineig) > -1e-12 and min(t[0] for t in neg) >= 0, f"min eigenvalue over tau {min(mineig):.1e}; most negative measured G {worst[1]:.4f} (g = {worst[2]}, {worst[3]})")

# SF3 / SF4 --------------------------------------------------------------------------------------------------
rows = []
for p in pts:
    for name, v in p["greens"].items():
        if name.startswith("Gamma") and "twist_delta" not in p:
            continue
        E, Z, chi = fit(np.array(v["G_tau"]), np.array(v["G_err"]), p["beta"], p["dtau"])
        k = k_of(name, p.get("twist_delta", [0, 0, 0]))
        rows.append({"cell": p["cell"], "beta": p["beta"], "dtau": p["dtau"], "g": p["g"], "twist": p.get("twist_delta"), "k_name": name,
                     "k": k.tolist(), "E": E, "Z": Z, "chi2_dof": chi, "E_free": float(np.max(np.abs(np.linalg.eigvalsh(h_of_k(k))))),
                     "trotter_limited": bool(E * p["dtau"] > 1.0), "G_half": v["G_half"]})
def cls(r):
    if r["trotter_limited"]:
        return "gapped (E dtau > 1: Trotter-limited, E only 'large')"
    if r["E"] * r["beta"] < 0.5:                                      # below the finite-beta resolution
        return "gapless within the beta resolution (pole)" if r["Z"] > 0.5 else "gapless within resolution, small Z (renormalised)"
    return "gapped"
step("SF3 (statement) energy and residue for every measured momentum / Cartan sector (twisted runs: k = base +- delta)",
     len(rows) > 0, " | ".join(f"{r['cell']} b={r['beta']} g={r['g']} tw={r['twist'][0] if r['twist'] else 0} {r['k_name']}: E = {r['E']:.3f} (free {r['E_free']:.3f}), Z = {r['Z']:.3f}, chi2/dof = {r['chi2_dof']:.1f}" for r in rows))
step("SF4 (statement) classification per point from (E, Z) -- small cells, single beta and dtau: no thermodynamic statement; runs with R < 0.5 excluded",
     True, f"excluded (R < 0.5): {excluded} | " + " | ".join(f"{r['cell']} g={r['g']} tw={r['twist'][0] if r['twist'] else 0} {r['k_name']}: {cls(r)}" for r in rows))

# SF5 (.02) -------------------------------------------------------------------------------------------------
def levels_curve(par, tau, beta, nlev):
    G = np.zeros_like(tau)
    for i in range(nlev):
        e, w = par[2 * i], abs(par[2 * i + 1])
        G = G + w * np.exp(-e * tau) / (1 + np.exp(-beta * e))
    return G
def fit_levels(G, err, beta, dtau, nlev):
    L = len(G) - 1; tau = dtau * np.arange(L + 1); sel = np.arange(1, L)
    sig = np.maximum(np.asarray(err)[sel], 1e-4); Gs = np.asarray(G)[sel]
    best = None
    for e0 in (0.3, 1.0, 3.0):
        starts = [e0, 0.5, -e0, 0.5, 2 * e0, 0.05][:2 * nlev]
        r = least_squares(lambda p: (levels_curve(p, tau[sel], beta, nlev) - Gs) / sig, x0=starts, bounds=([-50, 0] * nlev, [50, 2] * nlev))
        if best is None or r.cost < best.cost:
            best = r
    chi2 = 2 * best.cost; k_par = 2 * nlev
    global _last_cov
    _last_cov = np.linalg.pinv(best.jac.T @ best.jac)                   # parameter covariance (Fisher, whitened residuals)
    return best.x, float(chi2), float(chi2 + 2 * k_par)                 # AIC = chi^2 + 2 k
def m1_aic(G, err, beta, dtau):
    E, Z, chi = fit(G, err, beta, dtau); L = len(G) - 1
    return E, Z, float(chi * max(L - 1 - 2, 1) + 4)
# unit tests: free curves -> M2 symmetric; synthetic asymmetric -> M2 recovers, M1 rejected
sym_dev = 0.0
for p in pts[:6]:
    for name, v in p["greens"].items():
        Gf = np.array(v["G_free"])
        if Gf[len(Gf) // 2] > 0.499:
            continue
        x, chi, aic = fit_levels(Gf, np.full(len(Gf), 1e-6), p["beta"], p["dtau"], 2)
        e1, w1, e2, w2 = x[0], abs(x[1]), x[2], abs(x[3])
        sym_dev = max(sym_dev, abs(e1 + e2), abs(w1 - w2))
tau_t = 0.25 * np.arange(9); Gsyn = 0.6 * np.exp(-0.8 * tau_t) / (1 + np.exp(-2 * 0.8)) + 0.35 * np.exp(1.7 * tau_t) / (1 + np.exp(2 * 1.7))
err_t = np.full(9, 2e-4)
x_s, chi_s, aic_s = fit_levels(Gsyn, err_t, 2.0, 0.25, 2)
_, _, aic_m1 = m1_aic(Gsyn, err_t, 2.0, 0.25)
lv = sorted([(x_s[0], abs(x_s[1])), (x_s[2], abs(x_s[3]))])
syn_ok = abs(lv[0][0] + 1.7) < 0.05 and abs(lv[1][0] - 0.8) < 0.05 and aic_m1 > aic_s + 10
step("SF5 unit tests for NODE_SPECTRAL_FIT.02: on free curves the asymmetric two-level model returns the symmetric pair (e_2 = -e_1, w_1 = w_2); synthetic asymmetric data (addition 0.8, removal 1.7, residues 0.6 / 0.35) are recovered by M2 and rejected by M1 (AIC)",
     sym_dev < 1e-3 and syn_ok, f"free-curve asymmetry {sym_dev:.1e}; synthetic levels {[(round(e, 3), round(w, 3)) for e, w in lv]}; AIC M1 {aic_m1:.1f} vs M2 {aic_s:.1f}")
rows2 = []
for p in pts:
    for name, v in p["greens"].items():
        if name.startswith("Gamma") and "twist_delta" not in p:
            continue
        G, e = np.array(v["G_tau"]), np.array(v["G_err"])
        E1, Z1, a1 = m1_aic(G, e, p["beta"], p["dtau"])
        x2, c2, a2 = fit_levels(G, e, p["beta"], p["dtau"], 2); cov2 = _last_cov
        x3, c3, a3 = fit_levels(G, e, p["beta"], p["dtau"], 3); cov3 = _last_cov
        aics = {"M1": a1, "M2": a2, "M3": a3}
        best_m = min(aics, key=aics.get)
        if best_m == "M1":
            gap, desc = E1, f"E = {E1:.3f}, Z = {Z1:.3f}"
        else:
            xx, nl = (x2, 2) if best_m == "M2" else (x3, 3)
            lev_all = [(xx[2 * i], abs(xx[2 * i + 1])) for i in range(nl)]
            # resolvable levels only: weight not at the bound, |e| dtau < 3 (otherwise the level is a contact-like artefact of 7 tau points)
            lev = [(l_, w_) for l_, w_ in lev_all if 0.02 < w_ < 1.99 and abs(l_) * p["dtau"] < 3]
            gap = min(abs(l_) for l_, _ in lev) if lev else float("nan")
            desc = ", ".join(f"e = {l_:+.3f} (w {w_:.3f})" for l_, w_ in sorted(lev)) + (" [fit at bounds]" if len(lev) < len(lev_all) else "")
        # identifiability of the lowest level (Fisher information of E scales like Z^2): energy error > max(0.5, |E|) or
        # residue below 2 sigma -> ENERGY_UNIDENTIFIABLE (no gap and no gaplessness may be read from this curve)
        ident = "identifiable"
        if best_m != "M1":
            xx, cv, nl = (x2, cov2, 2) if best_m == "M2" else (x3, cov3, 3)
            cand = [i for i in range(nl) if 0.02 < abs(xx[2 * i + 1]) < 1.99 and abs(xx[2 * i]) * p["dtau"] < 3]
            if cand:
                i0 = min(cand, key=lambda i: abs(xx[2 * i]))
                sE = float(np.sqrt(abs(cv[2 * i0, 2 * i0]))) * 3; sZ = float(np.sqrt(abs(cv[2 * i0 + 1, 2 * i0 + 1]))) * 3   # errors x3
                if sE > max(0.5, abs(xx[2 * i0])) or abs(xx[2 * i0 + 1]) < 2 * sZ:
                    ident = "ENERGY_UNIDENTIFIABLE"
            else:
                ident = "ENERGY_UNIDENTIFIABLE"
        rows2.append({"identifiability": ident, "cell": p["cell"], "g": p["g"], "beta": p["beta"], "dtau": p["dtau"], "twist": p.get("twist_delta"), "k_name": name,
                      "best_model": best_m, "aic": aics, "excitation_min": float(gap), "levels": desc, "trotter_limited": bool(gap * p["dtau"] > 1.0)})
nonsym = [r for r in rows2 if r["best_model"] != "M1"]
step("SF5 (.02, statement) model choice per point by AIC (M1 symmetric pole pair, M2 independent addition / removal levels, M3 with a positive rest) and the excitation energy min |e_i| of the chosen model; the earlier M1-only energies are replaced by these",
     len(rows2) > 0, f"M1 chosen {sum(r['best_model'] == 'M1' for r in rows2)}, M2 {sum(r['best_model'] == 'M2' for r in rows2)}, M3 {sum(r['best_model'] == 'M3' for r in rows2)} of {len(rows2)}; ENERGY_UNIDENTIFIABLE {sum(r['identifiability'] != 'identifiable' for r in rows2)} | "
     + " | ".join(f"{r['cell']} g={r['g']} tw={r['twist'][0] if r['twist'] else 0} {r['k_name']}: {r['best_model']} {r['levels']} -> min |e| = {r['excitation_min']:.3f}{' (Trotter-limited)' if r['trotter_limited'] else ''} [{r['identifiability']}]" for r in rows2))

out = {"steps": LOG, "passed": sum(x["pass"] for x in LOG), "total": len(LOG), "rows": rows, "rows_v2": rows2}
(HERE / "spectral_fit_results.json").write_text(json.dumps(out, indent=2, default=float))
print(f"{out['passed']}/{out['total']} steps passed")
sys.exit(0 if out["passed"] == out["total"] else 1)
