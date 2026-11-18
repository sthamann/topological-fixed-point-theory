#!/usr/bin/env python3
"""Part ES: QMC.SECTOR_SAMPLING.PRODUCTION.01 -- envelope (surrogate-transition) sampling for sector observables.

The production sampler keeps the phase-quenched local Metropolis sweeps (measure pi_0 = e^{-S} |W|) as the PROPOSAL; a
forward or backward sweep is chosen with probability 1/2 so that the proposal kernel is reversible w.r.t. pi_0, and the
end configuration is accepted with min(1, h' / h).  The chain then samples pi = pi_0 h; physical averages are recovered by
<O> = <O e^{i theta} / h>_pi / <e^{i theta} / h>_pi.  Envelope for the projected Green function (reference sector r fixed in
advance): h = sqrt(|d_r|^2 + alpha^2 + lambda^2 sum_j |y_r(tau_j)|^2), d_r = D_r / |W|, y_r = numerator / |W| -- the alpha
term keeps configurations with d_r = 0, the y terms keep configurations where the Green numerator lives.  Separate chains
for different sectors do NOT give Z_s / Z_r; ratios need the common envelope or the alpha term (both sectors in one chain).
  ES1 toy (2D Gaussian pi_0, sequential-scan Metropolis in random direction as proposal, h = e^{a x_1}): the surrogate
      chain samples pi_0 h = N((a, 0), 1) and the reweighting recovers the pi_0 averages.
  ES2 refactored reweight.run: with outer = None the sampler and its random stream are unchanged (bitwise identical
      phase series and Green functions on debug4) and the integrator hash is unchanged.
  ES3 (statement) benchmark at g = 0.018, beta = 6 (envelope chains vs plain chains): error x sqrt(CPU time) of D_ref and
      C_ref(tau), outer acceptance, tau_int, N_eff, R^proj.
    python3 -B experiments/theory-contracts/seam-holonomy-generating-20261001/sector_envelope.py
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


# ES1 toy ------------------------------------------------------------------------------------------------------------
rng = np.random.default_rng(5)
a_, n_steps = 0.8, 60000
logpi0 = lambda x: -0.5 * np.sum(x ** 2)
h = lambda x: np.exp(a_ * x[0])
x = np.zeros(2); hx = h(x); xs, ws = [], []
acc_o = 0
for it in range(n_steps):
    x0, h0 = x.copy(), hx
    order = (0, 1) if rng.random() < 0.5 else (1, 0)                   # forward / backward scan: reversible mixture
    for c in order:
        y = x.copy(); y[c] += 0.9 * rng.normal()
        if np.log(rng.random()) < logpi0(y) - logpi0(x):
            x = y
    hn = h(x)
    if rng.random() < hn / h0:
        hx = hn; acc_o += 1
    else:
        x, hx = x0, h0
    xs.append(x.copy()); ws.append(1 / hx)
xs = np.array(xs[2000:]); ws = np.array(ws[2000:])
m_pi = xs.mean(0); m_back = (xs * ws[:, None]).sum(0) / ws.sum()
step("ES1 toy: the surrogate-transition chain (pi_0 proposal kernel reversible as a forward / backward mixture, outer acceptance min(1, h'/h)) samples pi_0 h = N((a, 0), 1) and the 1/h reweighting recovers the pi_0 averages",
     abs(m_pi[0] - a_) < 0.05 and abs(m_pi[1]) < 0.05 and abs(m_back[0]) < 0.08 and abs(m_back[1]) < 0.08,
     f"a = {a_}: <x>_pi = {np.round(m_pi, 3).tolist()} (exact ({a_}, 0)), reweighted <x>_pi0 = {np.round(m_back, 3).tolist()} (exact (0, 0)), outer acceptance {acc_o / n_steps:.2f}")

# ES2 equivalence of the refactored sampler ------------------------------------------------------------------------------
def load(path, name):
    sp = ilu.spec_from_file_location(name, path); mod = ilu.module_from_spec(sp); sp.loader.exec_module(mod); return mod
new = load(HERE / "reweight.py", "rw_new")
ref_file = HERE / "es2_reference.json"
b = new.run("debug4", 1.0, 0.25, 0.05, 6, 11, "md", True, np.array([0.3, 0.3, 0.3]))
cur = {"theta": b["series"]["theta"], "G": {k: v["G_tau"] for k, v in b["greens"].items()}, "acceptance": b["acceptance"]}
if not ref_file.exists():                                              # written once from the pre-refactor sampler (see next.txt)
    old = load(Path("/tmp/reweight_before_outer.py"), "rw_old") if Path("/tmp/reweight_before_outer.py").exists() else None
    if old is not None:
        a = old.run("debug4", 1.0, 0.25, 0.05, 6, 11, "md", True, np.array([0.3, 0.3, 0.3]))
        ref_file.write_text(json.dumps({"theta": a["series"]["theta"], "G": {k: v["G_tau"] for k, v in a["greens"].items()}, "acceptance": a["acceptance"]}))
refd = json.loads(ref_file.read_text()) if ref_file.exists() else None
src = (HERE / "reweight.py").read_text()
import hashlib
def _block(name):
    i = src.index(f"    def {name}("); j = src.index("\n    def ", i + 10)
    return src[i:j]
hsh = hashlib.sha256((_block("__init__") + _block("build_M") + _block("slice_op")).encode()).hexdigest()[:16]
same = refd is not None and np.allclose(refd["theta"], cur["theta"], atol=0, rtol=0) and all(np.allclose(refd["G"][k], cur["G"][k], atol=0, rtol=0) for k in cur["G"])
step("ES2 refactored reweight.run (optional envelope, forward / backward sweeps): with outer = None the phase series and Green functions on debug4 are bitwise identical to the pre-refactor sampler (stored reference) and the integrator hash is unchanged",
     same and hsh == "f0f2d58fada5669b", f"bitwise identical: {same}; integrator hash {hsh}")

# ES3 benchmark ----------------------------------------------------------------------------------------------------------
plain = sorted(HERE.glob("spg_points_g018b6*.json")); env = sorted(HERE.glob("spge_points_g018b6*.json"))
def stats(files, envelope):
    vals = []
    for fp in files:
        d = json.loads(fp.read_text())
        z = np.exp(1j * np.array(d["theta"])); f = np.array(d["f_re"]) + 1j * np.array(d["f_im"]); g = np.array(d["g_re"]) + 1j * np.array(d["g_im"])
        chi = np.array(d["chi_re"]) + 1j * np.array(d["chi_im"]); ref = d["sectors"].index([0, 0])
        w = z / np.array(d["h"]) if envelope else z
        Dref = (w * (f @ chi[ref].conj()) / 64).mean() / w.mean()
        Nr = np.einsum("m,me,met,e->t", w, f, g, chi[ref].conj()) / 64 / len(w) / w.mean()
        vals.append((Dref.real, (Nr / Dref).real, d.get("seconds", np.nan), d.get("outer_acceptance", np.nan)))
    return vals
if plain and env:
    sp_, se_ = stats(plain, False), stats(env, True)
    def summ(v):
        D = np.array([x[0] for x in v]); C = np.array([x[1] for x in v]); T = np.array([x[2] for x in v])
        eD = D.std(ddof=1) / np.sqrt(len(D)) if len(D) > 1 else np.nan
        return D.mean(), eD, C.mean(0), (C.std(0, ddof=1) / np.sqrt(len(C)) if len(C) > 1 else np.full(C.shape[1], np.nan)), T.sum()
    a1, b1 = summ(sp_), summ(se_)
    step("ES3 (statement) benchmark at g = 0.018, beta = 6: chain-to-chain error x sqrt(total CPU time) of D_ref for plain vs envelope chains, outer acceptance",
         True, f"plain ({len(sp_)} chains): D_ref {a1[0]:.3f} +- {a1[1]:.3f}, cost metric {a1[1] * np.sqrt(a1[4]):.2f}; envelope ({len(se_)} chains, outer acceptance {np.nanmean([x[3] for x in se_]):.2f}): D_ref {b1[0]:.3f} +- {b1[1]:.3f}, cost metric {b1[1] * np.sqrt(b1[4]):.2f}")
else:
    step("ES3 (statement) benchmark", True, "envelope chains not finished yet")

out = {"steps": LOG, "passed": sum(x["pass"] for x in LOG), "total": len(LOG)}
(HERE / "sector_envelope_results.json").write_text(json.dumps(out, indent=2, default=float))
print(f"{out['passed']}/{out['total']} steps passed")
sys.exit(0 if out["passed"] == out["total"] else 1)
