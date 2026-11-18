#!/usr/bin/env python3
"""Part ZB: NODE.ZERO.BRACKET.01 -- a sign change of a signed thermal spectral line.

If the same identifiable spectral line has e(k1) < 0 < e(k2) with residues Z(k1), Z(k2) > Z_min, it must cross omega = 0 in
between.  Status: THERMAL_SPECTRAL_BRACKET, NOT 'mirror gapless': a thermal excited -> excited transition can pass through
omega = 0 while the ground-sector gap stays finite.  A bracket becomes ground-connected only if its residue survives
cooling (CHIRAL4D.GROUNDLINE.RESIDUE.01: Z_beta -> Z_0 > 0) and is then tested by the active sector trace.  Per point the
lowest resolvable level of the two- / three-level fit (independent addition and removal energies, positive weights) is
taken with its sign; errors are the fit errors inflated x3 (binned errors are underestimated).  A bracket requires both
|e| > 2 sigma_e and both Z > Z_min = 0.3.  A bracket is then bisected by new twist points.
  ZB1 unit test: a synthetic signed branch e(q) = -0.2 + 0.6 q (residue 0.7, plus a high level) sampled at q = 0, 0.15, 0.3,
      0.45, 0.6 with realistic noise is bracketed on either side of its zero q = 1/3.
  ZB2 unit test: a gapped branch e(q) = -0.35 - 0.3 q^2 gives no bracket.
  ZB3 (statement) the measured g = 0.03 double star (16 sites, beta = 2): signed lowest levels, residues and brackets.
    python3 -B experiments/theory-contracts/seam-holonomy-generating-20261001/zero_bracket.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
LOG: list[dict] = []
Z_MIN, SIG_INFL = 0.3, 3.0


def step(name: str, ok: bool, detail: str = "") -> None:
    LOG.append({"step": name, "pass": bool(ok), "detail": detail})
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))
    sys.stdout.flush()


def levels_curve(par, tau, beta, nlev):
    G = np.zeros_like(tau)
    for i in range(nlev):
        e, w = par[2 * i], abs(par[2 * i + 1])
        G = G + w * np.exp(-e * tau) / (1 + np.exp(-beta * e))
    return G


def lowest_signed(G, err, beta, dtau):
    """lowest resolvable level (sign, weight, sigma_e) from the better (AIC) of the two- and three-level fits."""
    L = len(G) - 1; tau = dtau * np.arange(L + 1); sel = np.arange(1, L)
    sig = np.maximum(np.asarray(err)[sel], 1e-4); Gs = np.asarray(G)[sel]
    best = None
    for nlev in (2, 3):
        for e0 in (0.2, 0.8, 2.5):
            x0 = [e0, 0.5, -e0, 0.5, 2 * e0, 0.05][:2 * nlev]
            r = least_squares(lambda p: (levels_curve(p, tau[sel], beta, nlev) - Gs) / sig, x0=x0, bounds=([-50, 0] * nlev, [50, 2] * nlev))
            aic = 2 * r.cost + 4 * nlev
            if best is None or aic < best[0]:
                best = (aic, r, nlev)
    _, r, nlev = best
    cov = np.linalg.pinv(r.jac.T @ r.jac)
    cand = [i for i in range(nlev) if 0.02 < abs(r.x[2 * i + 1]) < 1.99 and abs(r.x[2 * i]) * dtau < 3]
    if not cand:
        return None
    i0 = min(cand, key=lambda i: abs(r.x[2 * i]))
    return float(r.x[2 * i0]), float(abs(r.x[2 * i0 + 1])), float(SIG_INFL * np.sqrt(abs(cov[2 * i0, 2 * i0]))), float(SIG_INFL * np.sqrt(abs(cov[2 * i0 + 1, 2 * i0 + 1])))


def brackets(points):
    """points: list of (label, q, e, w, sig_e, sig_w); returns pairs with a significant sign change and residues > Z_MIN."""
    good = [p for p in points if p[2] is not None and abs(p[2]) > 2 * p[4] and p[3] - 2 * p[5] > Z_MIN]
    out = []
    for i in range(len(good)):
        for j in range(i + 1, len(good)):
            if good[i][2] * good[j][2] < 0:
                out.append((good[i][0], good[j][0], good[i][2], good[j][2]))
    return out, len(good)


def synth_curve(e_low, w_low, beta=2.0, dtau=0.25, L=8, noise=0.002, rng=None):
    tau = dtau * np.arange(L + 1)
    G = w_low * np.exp(-e_low * tau) / (1 + np.exp(-beta * e_low)) + (1 - w_low) * np.exp(-2.5 * tau) / (1 + np.exp(-beta * 2.5))
    return G + noise * rng.normal(size=L + 1), np.full(L + 1, noise)


rng = np.random.default_rng(41)
qs = [0.0, 0.15, 0.3, 0.45, 0.6]
pts1 = []
for q in qs:
    G, e = synth_curve(-0.2 + 0.6 * q, 0.7, rng=rng)
    lo = lowest_signed(G, e, 2.0, 0.25)
    pts1.append((f"q={q}", q, *(lo if lo else (None, 0, 1, 1))))
br1, ng1 = brackets(pts1)
step("ZB1 unit test: a synthetic signed branch e(q) = -0.2 + 0.6 q (residue 0.7) is bracketed by two significant levels of opposite sign on either side of its zero q = 1/3",
     any(min(float(a[2:]), float(b[2:])) < 1 / 3 < max(float(a[2:]), float(b[2:])) for a, b, *_ in br1),
     f"fitted (q, e, Z): {[(p[1], round(p[2], 3) if p[2] is not None else None, round(p[3], 2)) for p in pts1]}; brackets {[(a, b) for a, b, *_ in br1]}")
pts2 = []
for q in qs:
    G, e = synth_curve(-0.35 - 0.3 * q ** 2, 0.7, rng=rng)
    lo = lowest_signed(G, e, 2.0, 0.25)
    pts2.append((f"q={q}", q, *(lo if lo else (None, 0, 1, 1))))
br2, ng2 = brackets(pts2)
step("ZB2 unit test: a gapped branch e(q) = -0.35 - 0.3 q^2 gives no bracket",
     len(br2) == 0, f"fitted e: {[round(p[2], 3) if p[2] is not None else None for p in pts2]}; brackets {br2}")

files = sorted(HERE.glob("reweight_points_md_cubic16_s13_*.json"))
if files:
    pts = []
    for f in files:
        p = json.loads(f.read_text())[-1]
        cands = [("centre:" + f.name.split("_s13_")[1][:-5], np.zeros(3), p["greens"]["P+d"])]
        for t in p.get("multitwist", []):
            cands.append((f"{f.name.split('_s13_')[1][:-5]}:{np.round(np.array(t['twist']) - 0.3, 3).tolist()}", np.array(t["twist"]) - 0.3, t["greens"]["P+d"]))
        for lab, q, v in cands:
            lo = lowest_signed(np.array(v["G_tau"]), np.array(v["G_err"]), p["beta"], p["dtau"])
            pts.append((lab, q, *(lo if lo else (None, 0, 1, 1))))
    br, ng = brackets(pts)
    step("ZB3 (statement) measured g = 0.03 double star (16 sites, beta = 2): signed lowest levels with residues; a bracket would be a THERMAL_SPECTRAL_BRACKET (not yet ground-connected), its absence is no gap proof (beta = 2 resolution)",
         True, f"{len(pts)} points, {ng} significant (|e| > 2 sigma, Z > {Z_MIN}); brackets: {br if br else 'none'} | "
               + "; ".join(f"{lab}: e = {e_:+.3f}+-{se:.3f}, Z = {w_:.2f}" if e_ is not None else f"{lab}: unresolved" for lab, q, e_, w_, se, sw in pts))
else:
    step("ZB3 (statement) measured double star", True, "runs not finished yet")

out = {"steps": LOG, "passed": sum(x["pass"] for x in LOG), "total": len(LOG)}
(HERE / "zero_bracket_results.json").write_text(json.dumps(out, indent=2, default=float))
print(f"{out['passed']}/{out['total']} steps passed")
sys.exit(0 if out["passed"] == out["total"] else 1)
