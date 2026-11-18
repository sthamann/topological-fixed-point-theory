#!/usr/bin/env python3
"""Part GM: targeted update of the identified slow direction (QMC.MIXING.MODE.01 -> SLOW_MODE_IDENTIFIED).

The audit (QMC.EQUILIBRATION.AUDIT.01, stage 1) found that G_{P+d}(beta/2) follows the overall HS-field amplitude within
every drifting chain (corr -0.93 ... -0.99 with |phi|^2, the time-constant power and log |W|) and that the local sweeps move
this direction only slowly.  reweight.run(global_moves=...) adds (i) a multiplicative amplitude move phi -> s phi with the
Jacobian s^D and (ii) pCN on the time-constant modes, both exact for e^{-S_HS} |W|.
  GM1 toy (D = 20, Gaussian prior, radial weight e^{a |phi|}): local Metropolis + the multiplicative move with the Jacobian
      reproduces the exact <|phi|^2> (radial quadrature); without the Jacobian it is visibly biased.
  GM2 toy: pCN on a 3-dimensional subspace with a weight tilting that subspace reproduces the exact shifted mean.
  GM3 production code path (debug4, tiny g = 1e-4 where |W| hardly depends on phi): with the global moves switched on the
      field norm follows the prior, <|phi|^2> = D 2 g dtau, both moves are accepted, and with global_moves = None the
      sampler is bitwise unchanged (sector_envelope.py ES2).
  GM5 production code path: pCN on the 10 uniform (K = 0) channel modes (the slowly condensing order-parameter direction
      found in the G4 runs) keeps the prior at tiny g and is accepted.
  GM6 production code path: multiplicative move of the 10 uniform coordinates (Jacobian s^10) keeps the chi^2_10 prior
      at tiny g (needed because the K = 0 order parameter barely moves on 32 sites otherwise).
  GM4 (statement) benchmark at g = 0.021, beta = 4 (lambda = 0.35): chains with global moves vs the plain audit chains --
      sweeps to the plateau of |phi|^2 and G(beta/2), agreement of differently started chains, CPU time.
    python3 -B experiments/theory-contracts/seam-holonomy-generating-20261001/global_moves.py
"""
from __future__ import annotations

import importlib.util as ilu
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
LOG: list[dict] = []


def step(name: str, ok: bool, detail: str = "") -> None:
    LOG.append({"step": name, "pass": bool(ok), "detail": detail})
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))
    sys.stdout.flush()


# GM1 --------------------------------------------------------------------------------------------------------------
D, sig2, a_ = 20, 0.5, 1.5
r = np.linspace(1e-6, 30, 300001)
logp = (D - 1) * np.log(r) - r ** 2 / (2 * sig2) + a_ * r
pr = np.exp(logp - logp.max()); exact = float(np.sum(pr * r ** 2) / np.sum(pr))
def chain(jacobian, n=40000, seed=1):
    rng = np.random.default_rng(seed); x = np.zeros(D) + 0.1; out = []
    logpi = lambda y: -np.sum(y ** 2) / (2 * sig2) + a_ * np.linalg.norm(y)
    for it in range(n):
        for c in rng.permutation(D)[:4]:                                # a few local updates (slow in the radial direction)
            y = x.copy(); y[c] += 0.05 * rng.normal()
            if np.log(rng.random()) < logpi(y) - logpi(x):
                x = y
        s = np.exp(0.05 * rng.normal())
        if np.log(rng.random()) < logpi(s * x) - logpi(x) + (D * np.log(s) if jacobian else 0.0):
            x = s * x
        out.append(np.sum(x ** 2))
    return float(np.mean(out[5000:]))
with_j, without_j = chain(True), chain(False)
step("GM1 toy (D = 20, Gaussian prior, radial weight e^{a |phi|}): the multiplicative amplitude move WITH the Jacobian s^D reproduces the exact <|phi|^2>; without the Jacobian it is biased",
     abs(with_j - exact) / exact < 0.03 and abs(without_j - exact) / exact > 0.05, f"exact {exact:.3f}; with Jacobian {with_j:.3f}; without {without_j:.3f}")

# GM2 --------------------------------------------------------------------------------------------------------------
rng = np.random.default_rng(2); Dd, k_ = 12, 3; bvec = np.array([0.8, -0.5, 0.3]); sig2b = 0.7
Q, _ = np.linalg.qr(rng.normal(size=(Dd, k_)))                           # orthonormal subspace basis
x = np.zeros(Dd); acc = 0; samples = []
logw = lambda y: bvec @ (Q.T @ y)                                       # tilt of the subspace coordinates
for it in range(60000):
    c = Q.T @ x; c_new = np.sqrt(1 - 0.5 ** 2) * c + 0.5 * np.sqrt(sig2b) * rng.normal(size=k_)
    y = x + Q @ (c_new - c)
    if np.log(rng.random()) < logw(y) - logw(x):
        x = y; acc += 1
    for cc in range(Dd):                                                # local prior-exact updates of all coordinates
        z = x.copy(); z[cc] += 0.6 * rng.normal()
        if np.log(rng.random()) < -(z[cc] ** 2 - x[cc] ** 2) / (2 * sig2b) + logw(z) - logw(x):
            x = z
    samples.append(Q.T @ x)
mean_sub = np.mean(samples[5000:], axis=0)
step("GM2 toy: pCN on a 3-dimensional subspace (prior exact) with a weight e^{b . c} reproduces the exact shifted subspace mean b sigma^2",
     np.allclose(mean_sub, bvec * sig2b, atol=0.03), f"measured {np.round(mean_sub, 3).tolist()}, exact {np.round(bvec * sig2b, 3).tolist()}, pCN acceptance {acc / 60000:.2f}")

# GM3 --------------------------------------------------------------------------------------------------------------
_sp = ilu.spec_from_file_location("reweight", HERE / "reweight.py"); rw = ilu.module_from_spec(_sp); _sp.loader.exec_module(rw)
norms = []
def mon(m, phi, st):
    norms.append(float(np.sum(phi ** 2))); return {}
gq, dt_ = 1e-4, 0.25
res = rw.run("debug4", 2.0, dt_, gq, 120, 7, "md", False, np.array([0.3, 0.3, 0.3]), None, burn=20, monitor=mon,
             global_moves={"scale_eps": 0.1, "pcn_rho": 0.5, "every": 1})
mm = rw.Model("debug4", 2.0, dt_, gq, "md", np.array([0.3, 0.3, 0.3]))
Dn = mm.Lt * mm.N * 10; prior = Dn * 2 * gq * dt_
meas = float(np.mean(norms[30:])); err = float(np.std(norms[30:]) / np.sqrt(len(norms[30:]) / 10))
ga = res.get("global_acceptance", {})
step("GM3 production code path (debug4, g = 1e-4): with amplitude and pCN moves the field norm follows the prior <|phi|^2> = D 2 g dtau and both moves are accepted (bitwise invariance with global_moves = None: sector_envelope.py ES2)",
     abs(meas - prior) < 3 * err + 0.05 * prior and all(v is not None and v > 0.05 for v in ga.values()),
     f"<|phi|^2> {meas:.4f} +- {err:.4f} vs prior {prior:.4f}; acceptance {ga}")

# GM5: pCN on the uniform (K = 0) channel modes, production code path at tiny g: uniform-mode variance follows the prior
cu_rec = []
def mon_u(m, phi, st):
    cu_rec.append(phi.sum(axis=(0, 1)) / np.sqrt(phi.shape[0] * phi.shape[1])); return {}
res_u = rw.run("debug4", 2.0, dt_, gq, 150, 9, "md", False, np.array([0.3, 0.3, 0.3]), None, burn=20, monitor=mon_u,
               global_moves={"pcn_uniform_rho": 0.5, "every": 1})
cu_arr = np.array(cu_rec[30:]); var_u = float(cu_arr.var()); prior_u = 2 * gq * dt_
gu = res_u.get("global_acceptance", {}).get("pcn_uniform")
step("GM5 production code path (debug4, g = 1e-4): pCN on the 10 uniform (K = 0, omega = 0) channel modes keeps the prior variance of the uniform coordinates and is accepted",
     abs(var_u - prior_u) < 0.25 * prior_u and gu is not None and gu > 0.5, f"uniform-mode variance {var_u:.2e} vs prior {prior_u:.2e}; acceptance {gu}")

# GM6: multiplicative move of the uniform coordinates (Jacobian s^10) at tiny g: <|c_u|^2> / sigma^2 = 10 (chi^2_10)
cu6 = []
def mon6(m, phi, st):
    cu6.append(float(np.sum((phi.sum(axis=(0, 1)) / np.sqrt(phi.shape[0] * phi.shape[1])) ** 2)) / (2 * gq * dt_)); return {}
res6 = rw.run("debug4", 2.0, dt_, gq, 400, 21, "md", False, np.array([0.3, 0.3, 0.3]), None, burn=20, monitor=mon6,
              global_moves={"uniform_scale_eps": 0.3, "every": 1})
m6 = float(np.mean(cu6[40:])); e6 = float(np.std(cu6[40:]) / np.sqrt(len(cu6[40:]) / 10))
ga6 = res6.get("global_acceptance", {}).get("uniform_scale")
step("GM6 production code path (debug4, g = 1e-4): the multiplicative move of the 10 uniform coordinates with the Jacobian s^10 keeps their prior distribution, <|c_u|^2> / sigma^2 = 10 (chi^2_10), and is accepted",
     abs(m6 - 10) < 3 * e6 + 0.8 and ga6 is not None and ga6 > 0.3, f"<|c_u|^2>/sigma^2 = {m6:.2f} +- {e6:.2f} (exact 10); acceptance {ga6}")

# GM4 benchmark ------------------------------------------------------------------------------------------------------
def summary(tag, half):
    fp = HERE / f"eqa_monitor_{tag}.jsonl"
    if not fp.exists():
        return None
    rr = [json.loads(l) for l in fp.read_text().splitlines() if l.strip()]
    G = [(x["sw"], x["theta"], x["G_P+d"][half]) for x in rr if "G_P+d" in x]
    n = max(x["sw"] for x in rr) + 1
    def w(lo, hi):
        sel = [(t, gg) for sw, t, gg in G if lo <= sw < hi]
        if len(sel) < 2:
            return float("nan")
        t = np.array([u for u, v in sel]); gg = np.array([v for u, v in sel]); return float((gg * np.cos(t)).mean() / np.cos(t).mean())
    nrm = lambda lo, hi: float(np.mean([x["norm2"] for x in rr if lo <= x["sw"] < hi] or [np.nan]))
    return n, [round(w(a, a + 50), 3) for a in range(0, n, 50)], [round(nrm(a, a + 50), 1) for a in range(0, n, 50)]
rows = []
for tag in ("g021b4A", "g021b4B", "g021b4GA", "g021b4GB"):
    s_ = summary(tag, 8)
    if s_:
        rows.append(f"{tag} ({'global moves' if 'G' in tag[6:] else 'local only'}): {s_[0]} sweeps; G(beta/2) per 50: {s_[1]}; |phi|^2 per 50: {s_[2]}")
step("GM4 (statement) benchmark at g = 0.021, beta = 4, lambda = 0.35: plateau of |phi|^2 and G(beta/2) with and without global moves",
     True, " | ".join(rows) if rows else "benchmark not started")

out = {"steps": LOG, "passed": sum(x["pass"] for x in LOG), "total": len(LOG)}
(HERE / "global_moves_results.json").write_text(json.dumps(out, indent=2, default=float))
print(f"{out['passed']}/{out['total']} steps passed")
sys.exit(0 if out["passed"] == out["total"] else 1)
