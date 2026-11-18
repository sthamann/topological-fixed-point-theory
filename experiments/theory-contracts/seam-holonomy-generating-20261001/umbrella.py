#!/usr/bin/env python3
"""Part UM: QMC.UMBRELLA.AMPLITUDE.01 -- free-energy profile along the slow coordinate |phi|^2 at the two-branch point.

At g = 0.015, beta = 6, lambda = 0.32 two branches persist (|phi|^2 ~ 36-37 with G_{P+d}(beta/2) ~ 0.32, ~ 43-44 with ~ 0.20)
even with the global amplitude move; tempering in |W|^s would need ~ 20 replicas (the branches differ by Delta log |W| ~ 300
while log |W| fluctuates by 2-30 inside a branch).  The slow coordinate is known, so umbrella windows
U_i = kappa / 2 (|phi|^2 - c_i)^2 (reweight.run(bias=...), with local sweeps + global amplitude / pCN moves) map the
phase-quenched free energy F(|phi|^2) = -log P(|phi|^2); WHAM removes the bias.  Which minimum is lower decides the branch;
the barrier height tells whether the two-branch behaviour is a first-order-like coexistence or slow drift.
  UM1 WHAM unit test: synthetic double-well F(n) sampled in harmonic windows is recovered (minimum depths and positions).
  UM2 (statement) the measured profile F(|phi|^2) at g = 0.015, beta = 6 with the window statistics and G(beta/2) per bin.
  UM4 order-parameter windows (bias kind 'uniform', X = |M|^2 / N) through the production code path at tiny g (debug4,
      g = 1e-4): the umbrella mean force -kappa (<X> - c) reproduces the analytic prior force F'(X) = N/2 - 4/X of
      X = chi^2_10 / N (G4 order-parameter umbrella validated).
  UM5 (statement) G4 order-parameter profile per cell at g = 0.021, beta = 4: umbrella windows in X = |M|^2 / N (kappa 1),
      umbrella integration of the mean forces F'(<X>) = -kappa (<X> - c) (no window overlap needed), position X* of the
      minimum of F(X) per cell; rule (fixed 2026-10-06 before the window data): X* > 0 and X*_32 / X*_16 within [0.7, 1.4]
      -> FAIL_ORDERED (long-range uniform pairing order); X*_32 / X*_16 < 0.6 (minimum moves towards 0 ~ 1/N) -> SHORT_RANGE.
  UM3 (statement) the two minima and the barrier of the phase-corrected F_phys(A) = F_pq(A) - log <cos theta>_A (binwise),
      Delta F = F_strong - F_weak and the barrier height with block-bootstrap errors (blocks of 5 sweeps, 200 resamples),
      and the number of samples supporting each overlap between consecutive windows.
    python3 -B umbrella.py run <beta> <g> <lambda> <center> <kappa> <sweeps> <seed> <tag> [start.npy]   (env EQA_GLOBAL)
    python3 -B umbrella.py
"""
from __future__ import annotations

import importlib.util as ilu
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
LOG: list[dict] = []


def step(name: str, ok: bool, detail: str = "") -> None:
    LOG.append({"step": name, "pass": bool(ok), "detail": detail})
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))
    sys.stdout.flush()


if len(sys.argv) > 1 and sys.argv[1] == "run":
    _sp = ilu.spec_from_file_location("reweight", HERE / "reweight.py"); rw = ilu.module_from_spec(_sp); _sp.loader.exec_module(rw)
    beta, g_, lam, c, kap, sweeps, seed, tag = float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6]), int(sys.argv[7]), int(sys.argv[8]), sys.argv[9]
    kind = os.environ.get("UM_KIND", "norm2"); cell = os.environ.get("UM_CELL", "cubic16")   # 'uniform': window in X = |M|^2 / N
    start = np.load(sys.argv[10]) if len(sys.argv) > 10 else None
    gm = None
    if os.environ.get("EQA_GLOBAL"):
        vals = [float(x) for x in os.environ["EQA_GLOBAL"].split(",")]   # "eps,rho[,rho_uniform]"
        gm = {"scale_eps": vals[0], "pcn_rho": vals[1], "pcn_uniform_rho": vals[2] if len(vals) > 2 else 0.0,
              "uniform_scale_eps": vals[3] if len(vals) > 3 else 0.0, "every": 1}
    mon = HERE / f"um_monitor_{tag}.jsonl"; mon.write_text("")
    def monitor(m, phi, state):
        Mv = phi.sum(axis=(0, 1)); X = float(np.sum(Mv ** 2)) / (phi.shape[0] * phi.shape[1] * 2 * g_ * 0.25) / phi.shape[1]
        rec = {"sw": state["sw"], "theta": state["theta"], "logabsW": state["logabsW"], "norm2": float(np.sum(phi ** 2)), "X": X}
        if state["sw"] % 4 == 0:
            G = m.greens_k(phi, m.momenta()); rec["G_P+d"] = np.real(G["P+d"]).tolist()
        with mon.open("a") as fh:
            fh.write(json.dumps(rec) + "\n")
        return {}
    t0 = time.time()
    res = rw.run(cell, beta, 0.25, g_, sweeps, seed, "md", False, np.array([lam] * 3), None, burn=20, monitor=monitor,
                 save_phi=str(HERE / f"um_phi_{tag}.npy"), phi0=start, global_moves=gm, bias={"center": c, "kappa": kap, "kind": kind})
    meta = {"tag": tag, "cell": cell, "kind": kind, "beta": beta, "g": g_, "lambda": lam, "center": c, "kappa": kap, "sweeps": sweeps, "R": res["R"], "seconds": time.time() - t0,
            "global_acceptance": res.get("global_acceptance")}
    with (HERE / "um_runs.jsonl").open("a") as fh:
        fh.write(json.dumps(meta) + "\n")
    print(json.dumps(meta)); sys.exit(0)


def wham(samples, centers, kappa, bins):
    """samples: list of arrays of n per window; returns bin centres, F (min 0), window free energies."""
    H = np.array([np.histogram(x, bins=bins)[0] for x in samples], float); Nw = H.sum(1)
    kap_l = np.broadcast_to(np.asarray(kappa, float), (len(centers),))   # one kappa per window (scalar allowed)
    nb = 0.5 * (bins[1:] + bins[:-1]); U = np.array([0.5 * k_ * (nb - c) ** 2 for c, k_ in zip(centers, kap_l)])
    f = np.zeros(len(centers))
    for _ in range(5000):
        denom = np.sum(Nw[:, None] * np.exp(f[:, None] - U), axis=0)
        P = H.sum(0) / np.maximum(denom, 1e-300)
        f_new = -np.log(np.maximum(np.sum(P[None, :] * np.exp(-U), axis=1), 1e-300))
        f_new -= f_new[0]
        if np.max(np.abs(f_new - f)) < 1e-10:
            f = f_new; break
        f = f_new
    with np.errstate(divide="ignore"):
        F = -np.log(P / P.max())
    return nb, F, f


# UM1 --------------------------------------------------------------------------------------------------------------
rng = np.random.default_rng(4)
Ftrue = lambda n: 0.15 * (n - 36) ** 2 * (n - 44) ** 2 / 16 - 1.0 * (n > 40) * np.exp(-((n - 44) ** 2) / 4)   # two minima, the right one deeper
centers = np.arange(32, 47, 1.0); kap = 1.0; samp = []
for c in centers:
    x = c; xs = []
    for it in range(20000):
        y = x + 0.4 * rng.normal()
        if np.log(rng.random()) < -(Ftrue(y) + 0.5 * kap * (y - c) ** 2) + (Ftrue(x) + 0.5 * kap * (x - c) ** 2):
            x = y
        xs.append(x)
    samp.append(np.array(xs[1000:]))
bins = np.linspace(31, 48, 69)
nb, F, _ = wham(samp, centers, kap, bins)
ok = np.isfinite(F)
Fe = Ftrue(nb) - Ftrue(nb[ok]).min()
iL = np.argmin(np.where((nb < 40) & ok, F, np.inf)); iR = np.argmin(np.where((nb > 40) & ok, F, np.inf))
dF_meas = F[iR] - F[iL]; dF_true = Ftrue(nb[iR]) - Ftrue(nb[iL])
step("UM1 WHAM unit test: a synthetic double-well free energy sampled in harmonic umbrella windows is recovered (positions of both minima, their difference within 0.2)",
     abs(nb[iL] - 36) < 0.6 and abs(nb[iR] - 44) < 0.6 and abs(dF_meas - dF_true) < 0.2, f"minima at {nb[iL]:.2f} / {nb[iR]:.2f} (true 36 / ~44), F_right - F_left measured {dF_meas:+.2f}, true {dF_true:+.2f}")

# UM4 order-parameter umbrella at tiny g -------------------------------------------------------------------------------
_sp4 = ilu.spec_from_file_location("reweight", HERE / "reweight.py"); rw4 = ilu.module_from_spec(_sp4); _sp4.loader.exec_module(rw4)
m4 = rw4.Model("debug4", 2.0, 0.25, 1e-4, "md", np.array([0.3, 0.3, 0.3])); N4 = m4.N
cache4 = HERE / "umbrella_um4_cache.json"                               # the 18-minute check is rerun only with UM_RETEST=1
mf4 = json.loads(cache4.read_text()) if (cache4.exists() and not os.environ.get("UM_RETEST")) else []
for c4, sd in (() if mf4 else ((1.0, 11), (2.0, 12), (3.0, 13), (4.0, 14))):
    Xs = []
    def mon4(m, phi, st):
        Mv = phi.sum(axis=(0, 1)); Xs.append(float(np.sum(Mv ** 2)) / (phi.shape[0] * phi.shape[1] * 2 * 1e-4 * 0.25) / phi.shape[1]); return {}
    rw4.run("debug4", 2.0, 0.25, 1e-4, 400, sd, "md", False, np.array([0.3, 0.3, 0.3]), None, burn=20, monitor=mon4,
            global_moves={"scale_eps": 0.05, "pcn_rho": 0.0, "pcn_uniform_rho": 0.6, "every": 1}, bias={"center": c4, "kappa": 4.0, "kind": "uniform"})
    Xm = float(np.mean(Xs[40:])); mf4.append((Xm, -4.0 * (Xm - c4), N4 / 2 - 4 / Xm))
cache4.write_text(json.dumps(mf4))
err4 = max(abs(a - b) for _, a, b in mf4)
step("UM4 order-parameter umbrella (bias on X = |M|^2 / N) through the production code at tiny g (debug4): measured mean force -kappa (<X> - c) matches the analytic prior force N/2 - 4/X of X = chi^2_10 / N",
     err4 < 0.6, f"N = {N4}; (<X>, measured F', analytic F'): {[(round(a, 2), round(b, 2), round(c_, 2)) for a, b, c_ in mf4]}; max deviation {err4:.2f}")

# UM2 measured --------------------------------------------------------------------------------------------------------
files = sorted(HERE.glob("um_monitor_g015b6*.jsonl"))
runs = [json.loads(l) for l in (HERE / "um_runs.jsonl").read_text().splitlines()] if (HERE / "um_runs.jsonl").exists() else []
if files:
    meta = {}
    for r_ in runs:
        meta[r_["tag"]] = r_
    samp, cents, Gs, kaps, samp_recs = [], [], [], [], []
    for fp in files:
        tag = fp.stem.split("_")[-1]; rr = [json.loads(l) for l in fp.read_text().splitlines() if l.strip()]
        if tag not in meta or len(rr) < 40:
            continue
        rr = [x for x in rr if x["sw"] >= 20]
        samp.append(np.array([x["norm2"] for x in rr])); cents.append(meta[tag]["center"]); kaps.append(meta[tag]["kappa"]); samp_recs.append(rr)
        Gs += [(x["norm2"], x["theta"], x["G_P+d"][len(x["G_P+d"]) // 2]) for x in rr if "G_P+d" in x]
    if samp:
        lo, hi = min(x.min() for x in samp), max(x.max() for x in samp)
        bins = np.linspace(lo, hi, 41); nb, F, fw = wham(samp, cents, kaps, bins)
        # binwise phase: the conditional mean <e^{i theta}>_A does not depend on the bias (a function of A only), so all
        # windows are pooled per bin; F_phys(A) = F_pq(A) - log <cos theta>_A (the imaginary part averages to zero)
        allA = np.concatenate([np.array([x["norm2"] for x in rr_]) for rr_ in samp_recs]); allT = np.concatenate([np.array([x["theta"] for x in rr_]) for rr_ in samp_recs])
        Rb, Cb = [], []
        for i in range(len(nb)):
            sel = (allA >= bins[i]) & (allA < bins[i + 1])
            z_ = np.exp(1j * allT[sel]) if sel.sum() >= 3 else None
            Rb.append(float(abs(z_.mean())) if z_ is not None else np.nan); Cb.append(float(z_.mean().real) if z_ is not None else np.nan)
        Fphys = np.where(np.array(Cb) > 0, F - np.log(np.maximum(np.array(Cb), 1e-12)), np.nan); Fphys = Fphys - np.nanmin(Fphys)
        Garr = np.array(Gs); gb = []
        for i in range(len(nb)):
            sel = (Garr[:, 0] >= bins[i]) & (Garr[:, 0] < bins[i + 1])
            gb.append(float((Garr[sel, 2] * np.cos(Garr[sel, 1])).sum() / np.cos(Garr[sel, 1]).sum()) if sel.sum() >= 3 else np.nan)
        prof = "; ".join(f"{n:.1f}: F_pq {f_:.2f}, F_phys {fp_:.2f}, R {r_:.2f}, G {g:.3f}" for n, f_, fp_, r_, g in zip(nb[::3], F[::3], Fphys[::3], np.array(Rb)[::3], gb[::3]))
        # overlap chain: consecutive windows (sorted by mean) must share histogram support
        ws = sorted(zip([x.mean() for x in samp], [x.min() for x in samp], [x.max() for x in samp]))
        gaps = [(round(ws[i][2], 2), round(ws[i + 1][1], 2)) for i in range(len(ws) - 1) if ws[i][2] < ws[i + 1][1]]
        # umbrella integration (robust where WHAM overlap is missing): mean force F'(<n>_i) = -kappa_i (<n>_i - c_i)
        mf = sorted((float(x.mean()), -k_ * (float(x.mean()) - c_)) for x, c_, k_ in zip(samp, cents, kaps))
        ui = "; ".join(f"<n> {a:.2f}: F' {b:+.2f}" for a, b in mf)
        step("UM2 (statement) phase-quenched free-energy profile F(|phi|^2) at g = 0.015, beta = 6, lambda = 0.32 from umbrella windows (WHAM) with G(beta/2) per bin, plus the umbrella-integration mean forces",
             True, f"{len(samp)} windows (centres {sorted(cents)}, kappa {sorted(set(kaps))}); WHAM profile: {prof}; global minimum at |phi|^2 = {nb[np.nanargmin(F)]:.1f} (F_pq) / {nb[np.nanargmin(Fphys)]:.1f} (F_phys); overlap gaps between consecutive windows: {gaps if gaps else 'none'}; mean forces: {ui}")
    else:
        step("UM2 (statement) profile", True, "windows not finished")
else:
    step("UM2 (statement) profile", True, "umbrella runs not started")

def profile_phys(samp_list, recs_list, cents_, kaps_, bins):
    nb, F, _ = wham(samp_list, cents_, kaps_, bins)
    allA = np.concatenate([np.array([x["norm2"] for x in rr_]) for rr_ in recs_list]); allT = np.concatenate([np.array([x["theta"] for x in rr_]) for rr_ in recs_list])
    Cb = []
    for i in range(len(nb)):
        sel = (allA >= bins[i]) & (allA < bins[i + 1])
        Cb.append(float(np.cos(allT[sel]).mean()) if sel.sum() >= 3 else np.nan)
    Fp = np.where(np.array(Cb) > 0, F - np.log(np.maximum(np.array(Cb), 1e-12)), np.nan)
    return nb, Fp - np.nanmin(Fp)
def minima(nb, Fp):
    ok = np.isfinite(Fp)
    if ok.sum() < 5:
        return None
    x, y = nb[ok], Fp[ok]
    # barrier = highest point between the global minimum on the right half and the minimum on the left
    iS = int(np.argmin(np.where(x > 39.5, y, np.inf))); iW = int(np.argmin(np.where(x < 39.0, y, np.inf)))
    if not (iW < iS):
        return None
    iB = iW + int(np.argmax(y[iW:iS + 1]))
    return x[iW], y[iW], x[iS], y[iS], x[iB], y[iB]
if files and samp:
    lo, hi = min(x.min() for x in samp), max(x.max() for x in samp); bins = np.linspace(lo, hi, 41)
    nb, Fp = profile_phys(samp, samp_recs, cents, kaps, bins)
    mres = minima(nb, Fp)
    rngb = np.random.default_rng(17); boots = []
    for b_ in range(200):
        s2, r2 = [], []
        for x_, rr_ in zip(samp, samp_recs):
            nbk = max(1, len(x_) // 5); idx = np.concatenate([np.arange(k * 5, min(len(x_), k * 5 + 5)) for k in rngb.integers(0, nbk, nbk)])
            s2.append(x_[idx]); r2.append([rr_[i] for i in idx])
        nb2, Fp2 = profile_phys(s2, r2, cents, kaps, bins); m2 = minima(nb2, Fp2)
        if m2:
            boots.append((m2[3] - m2[1], m2[5] - m2[1]))
    ws = sorted(zip([x.mean() for x in samp], samp))
    supp = []
    for i in range(len(ws) - 1):
        a_, b_ = ws[i][1], ws[i + 1][1]; lo_, hi_ = b_.min(), a_.max()
        supp.append(int(((a_ >= lo_) & (a_ <= hi_)).sum() + ((b_ >= lo_) & (b_ <= hi_)).sum()) if hi_ >= lo_ else 0)
    if mres:
        bt = np.array(boots)
        dF, dFe = mres[3] - mres[1], (bt[:, 0].std() if len(bt) > 10 else np.nan)
        bar, bare = mres[5] - mres[1], (bt[:, 1].std() if len(bt) > 10 else np.nan)
        step("UM3 (statement) phase-corrected minima, barrier and Delta F = F_strong - F_weak at g = 0.015, beta = 6 (block bootstrap), overlap support",
             True, f"weak minimum at |phi|^2 = {mres[0]:.1f}, strong minimum at {mres[2]:.1f}, barrier at {mres[4]:.1f}; Delta F = {dF:+.2f} +- {dFe:.2f}, barrier above the weak minimum {bar:.2f} +- {bare:.2f} ({len(bt)} bootstrap samples); samples in the overlaps of consecutive windows: {supp}")
    else:
        step("UM3 (statement) minima / Delta F", True, "two minima not resolved yet")

prof5, Xstar = [], {}
runs5 = [json.loads(l) for l in (HERE / "um_runs.jsonl").read_text().splitlines()] if (HERE / "um_runs.jsonl").exists() else []
meta5 = {r_["tag"]: r_ for r_ in runs5 if r_.get("kind") == "uniform"}
for cellkey in ("g021X16", "g021X32"):
    pts = []
    for fp in sorted(HERE.glob(f"um_monitor_{cellkey}c*.jsonl")):
        tag = fp.stem.split("_")[-1]; rr = [json.loads(l) for l in fp.read_text().splitlines() if l.strip()]
        if tag not in meta5 or len(rr) < 30:
            continue
        late = [x for x in rr if x["sw"] >= max(20, len(rr) // 3)]
        Xm = float(np.mean([x["X"] for x in late])); Xs = float(np.std([x["X"] for x in late]))
        z = np.exp(1j * np.array([x["theta"] for x in late]))
        pts.append((Xm, -meta5[tag]["kappa"] * (Xm - meta5[tag]["center"]), meta5[tag]["center"], Xs, float(abs(z.mean()))))
    if len(pts) >= 3:
        pts.sort()
        Xv = np.array([p_[0] for p_ in pts]); Fp = np.array([p_[1] for p_ in pts])
        F = np.concatenate([[0.0], np.cumsum(0.5 * (Fp[1:] + Fp[:-1]) * np.diff(Xv))])
        i0 = int(np.argmin(F)); Xstar[cellkey] = float(Xv[i0])
        prof5.append(f"{cellkey}: " + "; ".join(f"c {c_:g}: <X> {x_:.2f} (sd {s_:.2f}), F' {f_:+.2f}, R {r_:.2f}" for x_, f_, c_, s_, r_ in pts)
                     + f" -> F(X) by integration: " + ", ".join(f"{x_:.1f}: {f_:+.1f}" for x_, f_ in zip(Xv, F)) + f"; minimum at X* = {Xv[i0]:.2f}")
dec5 = "UNDECIDED"
if "g021X16" in Xstar and "g021X32" in Xstar and Xstar["g021X16"] > 0:
    rr5 = Xstar["g021X32"] / Xstar["g021X16"]
    dec5 = "FAIL_ORDERED" if 0.7 <= rr5 <= 1.4 else ("SHORT_RANGE" if rr5 < 0.6 else "UNDECIDED")
    prof5.append(f"X*_32 / X*_16 = {rr5:.2f} -> {dec5}")
step("UM5 (statement) G4 order-parameter profile F(X), X = |M|^2 / N, per cell (umbrella integration), minimum X* and the 16 / 32 comparison",
     True, " | ".join(prof5) if prof5 else "order-parameter windows not finished")

out = {"steps": LOG, "passed": sum(x["pass"] for x in LOG), "total": len(LOG), "g4_umbrella_decision": dec5}
(HERE / "umbrella_results.json").write_text(json.dumps(out, indent=2, default=float))
print(f"{out['passed']}/{out['total']} steps passed")
sys.exit(0 if out["passed"] == out["total"] else 1)
