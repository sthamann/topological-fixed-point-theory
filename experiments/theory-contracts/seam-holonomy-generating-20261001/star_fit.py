#!/usr/bin/env python3
"""Part SR: NODE.DISPERSION.13PT.01 -- local dispersion fit around a candidate node with the JOINT covariance.

Data: the double star at g = 0.03 (16 sites, beta = 2), centre k0 = P + 0.3 (1,1,1), radii h = 0.15 and 0.075, directions
+-e1, +-e2, +-e3, +-(e1+e2), +-(e1+e3), +-(e2+e3) (13 points per radius, 25 distinct momenta); three independent anchor
chains per radius, each carrying the centre and four targets by multi-twist reweighting (reweight_points_md_cubic16_s13_*).
Observable: E^2(k) from G(beta/2) = 1/(2 cosh(beta E/2)) (single pole with Z = 1 -- the fit is CONDITIONAL on that; Z is
checked separately by spectral_fit.py / spectral_bounds.py).  Model: E^2 = Delta^2 + (k - k*)^T Gm (k - k*), 10 parameters.
Errors: delete-one-block jackknife within each chain (block length >= 2 tau_int from the Gamma method), all quantities of a
chain jointly (targets of one chain are correlated through the shared HS configurations); chains are independent.
Search and confirmation are separated: the fit locates k*; confirmation must use NEW independent chains at the fixed k*.
  SR1 synthetic, exactly gapless (Delta^2 = 0) with a strong common-mode noise shared by all points of a chain: with the
      joint covariance the Delta^2 interval contains 0 at the nominal rate (diagonal-error coverage reported).
  SR2 synthetic gapped (Delta = 0.6): recovered within errors.
  SR3 the measured double star (when present): k*, Delta^2 with joint-covariance error, velocity metric, chi^2/dof, the
      13-point design rank (10), and the per-chain centre values (multi-seed check).
    python3 -B experiments/theory-contracts/seam-holonomy-generating-20261001/star_fit.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
LOG: list[dict] = []
_sb = (HERE / "spectral_bounds.py").read_text()
exec(compile(_sb.split("# ---------------------------------------------------------------------------------------------- Gamma method")[1]
             .split("# SB1 ---")[0], "gamma", "exec"))


def step(name: str, ok: bool, detail: str = "") -> None:
    LOG.append({"step": name, "pass": bool(ok), "detail": detail})
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))
    sys.stdout.flush()


BETA = 2.0
e1 = np.array([1, -1, 0]) / np.sqrt(2); e2 = np.array([1, 1, -2]) / np.sqrt(6); e3 = np.ones(3) / np.sqrt(3)
DIRS = [e1, -e1, e2, -e2, e3, -e3, e1 + e2, -(e1 + e2), e1 + e3, -(e1 + e3), e2 + e3, -(e2 + e3)]


def E2_of_G(G):
    """E^2 from G(beta/2) = 1/(2 cosh(beta E/2)); continued to E^2 < 0 for G > 1/2 (statistical overshoot)."""
    x = 1 / (2 * np.asarray(G, float))
    return np.where(x >= 1, (2 / BETA * np.arccosh(np.maximum(x, 1))) ** 2, -(2 / BETA * np.arccos(np.minimum(x, 1))) ** 2)


def model(par, q):
    ks, d2 = par[:3], par[3]
    L = np.zeros((3, 3)); L[np.tril_indices(3)] = par[4:10]
    Gm = L @ L.T                                                        # positive semidefinite velocity metric
    dq = q - ks[None, :]
    return d2 + np.einsum("ni,ij,nj->n", dq, Gm, dq)


def design_rank(q):
    """rank of the linearised design of E^2 = c + b.q + q^T A q (10 coefficients)."""
    X = np.column_stack([np.ones(len(q)), q, q[:, 0] ** 2, q[:, 1] ** 2, q[:, 2] ** 2, q[:, 0] * q[:, 1], q[:, 0] * q[:, 2], q[:, 1] * q[:, 2]])
    return np.linalg.matrix_rank(X, tol=1e-10)


def gls_fit(q, y, C):
    Lc = np.linalg.cholesky(C + 1e-14 * np.eye(len(y)))
    best = None
    for d0 in (0.0, 0.3, 1.0):
        x0 = np.concatenate([np.zeros(3), [d0], [1, 0, 1, 0, 0, 1]])
        r = least_squares(lambda p: np.linalg.solve(Lc, model(p, q) - y), x0)
        if best is None or r.cost < best.cost:
            best = r
    J = best.jac                                                       # Jacobian of the WHITENED residuals already
    cov_p = np.linalg.pinv(J.T @ J)
    return best.x, cov_p, 2 * best.cost / max(len(y) - 10, 1)


# SR1 / SR2 synthetic ---------------------------------------------------------------------------------------------
rng = np.random.default_rng(31)
def synthetic(d2_true, common=0.15, indiv=0.02):
    q, y, chains = [], [], []
    kst = np.array([0.02, -0.03, 0.04]); Gt = np.diag([1.0, 1.3, 0.8])
    qs = [np.zeros(3)] + [h * v for h in (0.15, 0.075) for v in DIRS]
    for ci in range(6):                                                 # six chains, each centre + four targets
        idx = [0] + list(range(1 + 4 * (ci % 3) + 12 * (ci // 3), 1 + 4 * (ci % 3) + 12 * (ci // 3) + 4))
        cm = common * rng.normal()                                       # common-mode noise of this chain
        for i in idx:
            qq = qs[i]; dq = qq - kst
            q.append(qq); y.append(d2_true + dq @ Gt @ dq + cm + indiv * rng.normal()); chains.append(ci)
    q, y, chains = np.array(q), np.array(y), np.array(chains)
    C_full = np.zeros((len(y), len(y))); C_diag = np.zeros_like(C_full)
    for ci in range(6):
        sel = np.where(chains == ci)[0]
        C_full[np.ix_(sel, sel)] = common ** 2
    C_full += indiv ** 2 * np.eye(len(y))
    C_diag = (common ** 2 + indiv ** 2) * np.eye(len(y))                # correct individual error bars, correlations ignored
    return q, y, C_full, C_diag
cover_full, cover_diag, n_rep = 0, 0, 60
for rep in range(n_rep):
    q, y, Cf, Cd = synthetic(0.0)
    pf, cf, _ = gls_fit(q, y, Cf); pd, cd, _ = gls_fit(q, y, Cd)
    cover_full += abs(pf[3]) <= 2 * np.sqrt(cf[3, 3])
    cover_diag += abs(pd[3]) <= 2 * np.sqrt(cd[3, 3])
step("SR1 synthetic exactly gapless double star (Delta^2 = 0) with a common-mode noise shared by the points of each chain: the joint-covariance fit contains Delta^2 = 0 in its 2 sigma interval at the nominal rate (the size of the under-coverage of a diagonal-error fit depends on the correlation structure; here it is moderate -- a covariance-free fit is not accepted in any case)",
     cover_full / n_rep > 0.9, f"2-sigma coverage of Delta^2 = 0: joint covariance {cover_full}/{n_rep}, diagonal errors {cover_diag}/{n_rep}; 13-point design rank {design_rank(np.array([np.zeros(3)] + [0.15 * v for v in DIRS]))} (7-point star: {design_rank(np.array([np.zeros(3)] + [0.15 * v for v in DIRS[:6]]))})")
q, y, Cf, Cd = synthetic(0.36, common=0.02)
pg, cg, _ = gls_fit(q, y, Cf)
step("SR2 synthetic gapped double star (Delta = 0.6, Delta^2 = 0.36) is recovered by the joint-covariance fit",
     abs(pg[3] - 0.36) < 3 * np.sqrt(cg[3, 3]) + 0.02, f"Delta^2 = {pg[3]:.3f} +- {np.sqrt(cg[3, 3]):.3f}")

# SR3 measured ------------------------------------------------------------------------------------------------------
files = sorted(HERE.glob("reweight_points_md_cubic16_s13_*.json"))
if files:
    rows, center_vals = [], []
    qs_all, ys_all, blocks = [], [], []
    for f in files:
        p = json.loads(f.read_text())[-1]
        ser = p.get("series")
        if ser is None:
            continue
        th = np.array(ser["theta"]); Ga = np.array(ser["G_half"]["P+d"])
        _, _, tau = ratio_error(Ga, th)
        nmeas = len(th); bl = max(1, int(np.ceil(2 * tau))); nb = max(4, nmeas // bl); edges = np.array_split(np.arange(nmeas), nb)
        k_center = np.array(p["twist_delta"]) - 0.3 * np.ones(3)          # q relative to k0 = P + 0.3 (1,1,1)
        quantities = [("center", np.zeros(3), Ga, th, np.zeros(nmeas))]
        for t, ts in zip(p["multitwist"], ser["targets"]):
            quantities.append(("target", np.array(t["twist"]) - 0.3 * np.ones(3), np.array(ts["G_half"]["P+d"]), np.array(ts["theta"]), np.array(ts["log_ratio"])))
        def estimates(mask):
            out = []
            for kind, qv, Gs, ths, lr in quantities:
                w = np.exp(lr[mask] - lr[mask].max()) * np.cos(ths[mask])
                out.append(np.sum(w * Gs[mask]) / np.sum(w))
            return np.array(out)
        full = estimates(np.ones(nmeas, bool))
        jk = np.array([estimates(np.array([i not in set(e) for i in range(nmeas)])) for e in edges])
        Cj = (nb - 1) * np.cov(E2_of_G(jk).T, bias=True)                 # delete-one-block jackknife covariance of the chain
        qs_all.append(np.array([qv for _, qv, *_ in quantities])); ys_all.append(E2_of_G(full)); blocks.append(Cj)
        center_vals.append((f.name, full[0], float(np.sqrt(max((nb - 1) * np.var(jk[:, 0]), 0))), tau))
    q = np.vstack(qs_all); y = np.concatenate(ys_all)
    C = np.zeros((len(y), len(y))); o = 0
    for Cj in blocks:
        n_ = Cj.shape[0]; C[o:o + n_, o:o + n_] = Cj + 1e-8 * np.eye(n_); o += n_
    pf, cf, chi = gls_fit(q, y, C)
    pdg, cdg, _ = gls_fit(q, y, np.diag(np.diag(C)))
    rk = design_rank(np.unique(np.round(q, 6), axis=0))
    step("SR3 (statement) measured double star at g = 0.03 (16 sites, beta = 2; E^2 from G(beta/2), conditional on Z = 1): joint-covariance fit of E^2 = Delta^2 + (k - k*)^T G (k - k*); multi-seed centre values; search only -- confirmation needs new chains at the fixed k*",
         True, ("FIT NOT DETERMINED (flat landscape within errors and an ill-conditioned jackknife covariance with ~10 blocks per chain for 5 quantities): " if (not np.isfinite(chi) or chi > 10) else "") + f"{len(files)} chains, {len(y)} values, distinct momenta {len(np.unique(np.round(q, 6), axis=0))}, design rank {rk}; "
               f"k* - k0 = {np.round(pf[:3], 3).tolist()}, Delta^2 = {pf[3]:.3f} +- {np.sqrt(cf[3, 3]):.3f} (diagonal-error fit: {pdg[3]:.3f} +- {np.sqrt(cdg[3, 3]):.3f}), chi^2/dof = {chi:.2f}; "
               f"centre per chain: " + ", ".join(f"{n.split('_s13_')[1][:-5]}: {g:.3f}+-{e:.3f} (tau_int {t:.1f})" for n, g, e, t in center_vals))
else:
    step("SR3 (statement) measured double star", True, "runs not finished yet")

out = {"steps": LOG, "passed": sum(x["pass"] for x in LOG), "total": len(LOG)}
(HERE / "star_fit_results.json").write_text(json.dumps(out, indent=2, default=float))
print(f"{out['passed']}/{out['total']} steps passed")
sys.exit(0 if out["passed"] == out["total"] else 1)
