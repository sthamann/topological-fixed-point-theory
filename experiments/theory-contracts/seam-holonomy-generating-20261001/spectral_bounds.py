#!/usr/bin/env python3
"""Part SB: CHIRAL4D.SPECTRAL_WEIGHT_BOUNDS.01 -- model-free bounds on low-energy spectral weight, and honest errors.

(1) Errors.  Binned errors underestimate the uncertainty (two independent chains differed by 5 binned sigma).  The
    Gamma method (U. Wolff, Comput. Phys. Commun. 156 (2004) 143) with automatic windowing gives the integrated
    autocorrelation time tau_int and the error of a function of means (here the reweighted ratio <O cos theta>/<cos theta>).
(2) Bounds.  G(tau) = int A(w) e^{-w tau} / (1 + e^{-beta w}) dw with A >= 0 and int A = 1.  For a window Omega the
    linear programs  min / max  W_Omega = int_{|w| < Omega} A  subject to |G_model(tau_j) - G(tau_j)| <= k sigma_j
    give the range of low-energy weight compatible with the data -- no pole model.  W_min(Omega) > 0 forces
    low-energy weight; W_max(Omega) ~ 0 for a fixed Omega > 0 is evidence for a gap.  At small beta the kernel cannot
    resolve energies below ~ 1/beta, so the bounds are informative only for Omega >~ 1/beta.
  SB1 unit test (Gamma method): an AR(1) series with known tau_int = (1 + rho) / (2 (1 - rho)) is recovered.
  SB2 unit test (bounds): exact free curves -- a gapless node (E = 0) gives W_min(Omega) = 1 for Omega >= 0.25, a gapped
      level pair (E = 1.41) gives W_max(0.5) small; the bounds bracket the true weight.
  SB3 (statement) bounds for the measured central points (errors inflated x3 as long as no Gamma-method error exists).
  SB4 matrix bounds over the 2 x 2 spin block (scalar LPs over spin directions, conservative): W_max^{any spin} and
      W_min^{some spin}; unit test with a gapless level in only one spin direction.
  SB5 (statement) matrix bounds U_Omega = W_max^{any spin}, L_Omega = W_min^{some spin} on the confirmation and residue-series
      runs for the preregistered fixed windows 0.1 / 0.25 / 0.5 / 1 (0.1 registered 2026-10-03T08:41:49 before these data).
    python3 -B experiments/theory-contracts/seam-holonomy-generating-20261001/spectral_bounds.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

HERE = Path(__file__).resolve().parent
LOG: list[dict] = []


def step(name: str, ok: bool, detail: str = "") -> None:
    LOG.append({"step": name, "pass": bool(ok), "detail": detail})
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))
    sys.stdout.flush()


# ---------------------------------------------------------------------------------------------- Gamma method
def uwerr(series, f=lambda m: m[0], S=1.5):
    """Wolff's Gamma method for f(<a_1>, ..., <a_n>); series shape (n_obs, N).  Returns (value, error, tau_int)."""
    a = np.atleast_2d(np.asarray(series, float)); n, N = a.shape
    m = a.mean(1)
    eps = 1e-7 * np.maximum(np.abs(m), 1e-10)
    grad = np.array([(f(m + eps[i] * np.eye(n)[i]) - f(m - eps[i] * np.eye(n)[i])) / (2 * eps[i]) for i in range(n)])
    proj = grad @ (a - m[:, None])
    G0 = np.mean(proj * proj)
    if G0 <= 0:
        return float(f(m)), 0.0, 0.5
    gam = [G0]; W = 0
    for t in range(1, N // 2):
        gam.append(np.mean(proj[:-t] * proj[t:]))
        tau = 0.5 + np.sum(gam[1:]) / G0
        if tau <= 0:
            tau = 1e-6
        g_ = np.exp(-t / (S * tau)) - S * tau / np.sqrt(t * N) if tau > 0 else -1
        if g_ < 0:
            W = t; break
    W = W or N // 2 - 1
    tau = 0.5 + np.sum(gam[1:W + 1]) / G0
    var = 2 * tau * G0 / N
    return float(f(m)), float(np.sqrt(max(var, 0))), float(tau)


def ratio_error(O, theta):
    """reweighted <O cos th> / <cos th> (imaginary parts negligible for R ~ 1) with Gamma-method error."""
    O = np.asarray(O); c = np.cos(theta)
    return uwerr(np.vstack([O * c, c]), f=lambda m: m[0] / m[1])


# SB1 --------------------------------------------------------------------------------------------------------
rng = np.random.default_rng(13)
res1 = []
for rho in (0.5, 0.8):
    x = np.zeros(40000)
    for t in range(1, len(x)):
        x[t] = rho * x[t - 1] + rng.normal()
    _, err, tau = uwerr(x)
    exact = (1 + rho) / (2 * (1 - rho))
    res1.append((rho, tau, exact, err, np.sqrt(2 * exact / len(x)) * np.std(x)))
step("SB1 unit test: the Gamma method recovers tau_int = (1 + rho) / (2 (1 - rho)) of AR(1) series and the corresponding error",
     all(abs(t - e) / e < 0.2 and abs(er - ee) / ee < 0.2 for _, t, e, er, ee in res1),
     "; ".join(f"rho = {r}: tau_int {t:.2f} (exact {e:.2f}), error {er:.4f} (expected {ee:.4f})" for r, t, e, er, ee in res1))


# ---------------------------------------------------------------------------------------------- LP bounds
def bounds(G, err, beta, dtau, Om, k_sig=3.0, wmax=14.0, nw=561):
    L = len(G) - 1
    tau = dtau * np.arange(L + 1)
    w = np.linspace(-wmax, wmax, nw)
    Kmat = np.exp(-np.outer(tau, w)) / (1 + np.exp(-beta * w))[None, :]
    sel = np.arange(0, L + 1)
    sig = np.maximum(k_sig * np.asarray(err), 1e-4)
    sig[0] = sig[-1] = max(sig[0], 1e-3)
    A_ub = np.vstack([Kmat[sel], -Kmat[sel]]); b_ub = np.concatenate([np.asarray(G)[sel] + sig[sel], -(np.asarray(G)[sel] - sig[sel])])
    A_eq = np.ones((1, nw)); b_eq = [1.0]
    c = (np.abs(w) < Om).astype(float)
    lo = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq, bounds=(0, None), method="highs")
    hi = linprog(-c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq, bounds=(0, None), method="highs")
    if lo.status != 0 or hi.status != 0:
        return np.nan, np.nan
    return float(lo.fun), float(-hi.fun)


def free_curve(E, beta, dtau, L):
    tau = dtau * np.arange(L + 1)
    return 0.5 * (np.exp(-E * tau) / (1 + np.exp(-beta * E)) + np.exp(E * tau) / (1 + np.exp(beta * E)))


# SB2 --------------------------------------------------------------------------------------------------------
beta, dtau, L = 2.0, 0.25, 8
G0 = free_curve(0.0, beta, dtau, L); G1 = free_curve(1.41, beta, dtau, L)
e_small = np.full(L + 1, 2e-4)
b_gapless = {Om: bounds(G0, e_small, beta, dtau, Om) for Om in (0.25, 0.5, 1.0)}
b_gapped = {Om: bounds(G1, e_small, beta, dtau, Om) for Om in (0.5, 1.0, 2.0)}
ok2 = b_gapless[0.25][0] > 0.9 and b_gapless[0.5][0] > 0.98 and b_gapped[0.5][1] < 0.2 and b_gapped[2.0][0] > 0.8
step("SB2 unit test (beta = 2, 9 tau points): a free gapless node forces W_min(0.25) > 0.9 and W_min(0.5) > 0.98; a free level pair at E = 1.41 allows at most W_max(0.5) < 0.2 and forces W_min(2) > 0.8 -- the bounds discriminate gapless from gapped at this beta (finite-beta smearing ~ 1/beta)",
     ok2, f"gapless: { {k: (round(a, 3), round(b, 3)) for k, (a, b) in b_gapless.items()} }; gapped: { {k: (round(a, 3), round(b, 3)) for k, (a, b) in b_gapped.items()} }")

# SB3 --------------------------------------------------------------------------------------------------------
pts = []
for f in sorted(HERE.glob("reweight_points_md_*.json")):
    for p in json.loads(f.read_text()):
        if "greens" in p and p["R"] >= 0.5:
            pts.append((f.name, p))
OMEGAS = (0.1, 0.25, 0.5, 1.0)                                          # predeclared fixed windows (t = 1 units); 0.1 added 2026-10-03 BEFORE the beta = 4 untwisted / 32-site data were evaluated
rows = []
for fname, p in pts:
    if p["g"] not in (0.01, 0.015, 0.018, 0.021, 0.024, 0.03, 0.1) or "multitwist" in p and False:
        continue
    for name, v in p["greens"].items():
        if not name.startswith("P"):
            continue
        G, e = np.array(v["G_tau"]), np.array(v["G_err"])
        for k_sig in (3.0, 5.0, 8.0, 12.0, 20.0):                       # smallest error inflation compatible with A >= 0, int A = 1
            bb = {Om: bounds(G, e, p["beta"], p["dtau"], Om, k_sig=k_sig) for Om in OMEGAS}
            if not np.isnan(bb[0.5][0]):
                break
        rows.append({"file": fname, "cell": p["cell"], "g": p["g"], "beta": p["beta"], "twist": p.get("twist_delta"), "k": name, "G_half": v["G_half"],
                     "W_0.1": bb[0.1], "W_0.25": bb[0.25], "W_0.5": bb[0.5], "W_1.0": bb[1.0], "k_sigma": k_sig})
def fmt(r):
    tw = r["twist"][0] if r["twist"] else 0
    return f"{r['cell']} g={r['g']} tw={tw} {r['k']}: G(b/2) = {r['G_half']:.3f}, W(0.1) in [{r['W_0.1'][0]:.2f}, {r['W_0.1'][1]:.2f}], W(0.25) in [{r['W_0.25'][0]:.2f}, {r['W_0.25'][1]:.2f}], W(0.5) in [{r['W_0.5'][0]:.2f}, {r['W_0.5'][1]:.2f}], W(1.0) in [{r['W_1.0'][0]:.2f}, {r['W_1.0'][1]:.2f}] (errors x{r['k_sigma']:g})"
step("SB3 (statement) model-free low-energy weight intervals for the measured P-sector curves: errors inflated by the smallest factor (>= 3) that makes the data compatible with A >= 0 and int A = 1 (a factor > 3 is itself evidence of underestimated errors); W_min > 0 forces low-energy weight, W_max ~ 0 indicates a gap at that momentum",
     len(rows) > 0, " | ".join(fmt(r) for r in rows))

# SB4 matrix (2 x 2 spin block) bounds via scalar LPs over a grid of spin directions v (no SDP solver available):
#   W_max^{any spin} = max over v and feasible A of v^dag W_Omega v  <=  max_v LP_max(v)      (matrix positivity relaxed: valid upper bound)
#   W_min^{some spin} = min over feasible A of lambda_max(W_Omega) >= max_v LP_min(v)       (lambda_max >= v^dag W v: valid lower bound)
def spin_dirs(n_th=7, n_ph=12):
    out = []
    for th in np.linspace(0, np.pi, n_th):
        for ph in (np.linspace(0, 2 * np.pi, n_ph, endpoint=False) if 0 < th < np.pi else [0.0]):
            out.append(np.array([np.cos(th / 2), np.exp(1j * ph) * np.sin(th / 2)]))
    return out
def matrix_bounds(Cs, err, beta, dtau, Om, k_sig=3.0):
    """Cs: (L+1, 2, 2) spin-resolved correlator (per spin component normalised: trace/2 = scalar G); err: proxy errors of the scalar curve.
    The smallest error inflation (>= k_sig) that makes EVERY spin direction feasible is used; infeasible -> (nan, nan, k)."""
    dirs = spin_dirs()
    for k_ in (k_sig, 5.0, 8.0, 12.0, 20.0):
        vals = []
        for v in dirs:
            Gv = np.real(np.einsum("i,tij,j->t", v.conj(), Cs, v))
            vals.append(bounds(Gv, err, beta, dtau, Om, k_sig=k_))
        if not any(np.isnan(a_) for a_, _ in vals):
            return max(a_ for a_, _ in vals), max(b_ for _, b_ in vals), k_
    return np.nan, np.nan, np.nan
# unit test: only spin direction up has a level at 0, down is gapped at 1.41 -> spin-averaged bound blurs, matrix bound finds it
tau_ = dtau * np.arange(L + 1)
Cup = free_curve(0.0, beta, dtau, L); Cdn = free_curve(1.41, beta, dtau, L)
Cs_test = np.zeros((L + 1, 2, 2)); Cs_test[:, 0, 0] = Cup; Cs_test[:, 1, 1] = Cdn
e_small = np.full(L + 1, 2e-4)
avg = bounds((Cup + Cdn) / 2, e_small, beta, dtau, 0.25)
mb = matrix_bounds(Cs_test, e_small, beta, dtau, 0.25)[:2]
Cg = np.zeros((L + 1, 2, 2)); Cg[:, 0, 0] = Cdn; Cg[:, 1, 1] = Cdn
mg = matrix_bounds(Cg, e_small, beta, dtau, 0.5)[:2]
step("SB4 matrix bounds over the 2 x 2 spin block (scalar LPs over a grid of spin directions; conservative): a gapless level in only one spin direction is found -- W_min^{some spin}(0.25) ~ 1 although the spin-averaged bound only forces ~ 1/2 -- and a block gapped in every direction gives a small W_max^{any spin}(0.5)",
     mb[0] > 0.9 and avg[0] < 0.6 and mg[1] < 0.2, f"one-direction gapless: spin-averaged W(0.25) in [{avg[0]:.2f}, {avg[1]:.2f}], matrix W_min^some = {mb[0]:.2f}, W_max^any = {mb[1]:.2f}; fully gapped: W_max^any(0.5) = {mg[1]:.3f}")
mrows = []
for fname, p in pts:                                                     # pts already excludes R < 0.5
    if "cf" not in fname and "_sus" not in fname and "_gl" not in fname and "_g24" not in fname and "_g27" not in fname and fname not in ("reweight_points_md_cubic16_ds3.json", "reweight_points_md_cubic16_ds4.json"):   # ds3 / ds4: beta = 2 points of the residue series
        continue
    for name, v in p["greens"].items():
        if "C_spin_re" not in v or not name.startswith("P"):
            continue
        Cs = np.array(v["C_spin_re"]) + 1j * np.array(v["C_spin_im"])
        e = np.array(v["G_err"])
        res_ = {Om: matrix_bounds(Cs, e, p["beta"], p["dtau"], Om) for Om in OMEGAS}
        mrows.append({"file": fname, "g": p["g"], "beta": p["beta"], "R": p["R"], "twist": p.get("twist_delta"), "k": name, "W": res_})
step("SB5 (statement) matrix spin-block bounds L_Omega = W_min^{some spin} / U_Omega = W_max^{any spin} for the preregistered fixed windows 0.1 / 0.25 / 0.5 / 1 on the confirmation and residue-series runs (errors: scalar-curve errors x3 as proxy for the components)",
     True, " | ".join(f"{r['file']} g={r['g']} beta={r['beta']} R={r['R']:.3f} {r['k']}: " + ", ".join(f"W({Om}) in [{a:.2f}, {b:.2f}] (x{k:g})" for Om, (a, b, k) in r["W"].items()) for r in mrows) if mrows else "confirmation runs not finished yet")

out = {"steps": LOG, "passed": sum(x["pass"] for x in LOG), "total": len(LOG), "rows": rows, "matrix_rows": mrows}
(HERE / "spectral_bounds_results.json").write_text(json.dumps(out, indent=2, default=float))
print(f"{out['passed']}/{out['total']} steps passed")
sys.exit(0 if out["passed"] == out["total"] else 1)
