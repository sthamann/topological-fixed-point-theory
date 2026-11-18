#!/usr/bin/env python3
"""Part MM: QMC.MIXING.MODE.01 -- why do independent chains disagree on the sector weights?

Benchmark point g = 0.018, beta = 6, lambda = 0.39 (16 sites).  Three admissible outcomes:
  BURN_IN_LIMITED       the chain means still move: a common trend between the first and the last quarter (> 3 sigma)
                        -> longer burn-in, no new sampler;
  MULTIMODAL_BARRIER    no common trend, every chain internally stable (within-chain quarter variance < 0.1 x the
                        between-chain variance) -> non-local updates / tempering;
  SLOW_MODE_IDENTIFIED  chains wander slowly and one global field observable explains the between-chain spread
                        (R^2 of the chain means > 0.7 AND above the 95 % quantile of the maximum R^2 over all features for
                        chain-shuffled targets -- multiple-testing control, added 2026-10-03 before the final evaluation)
                        -> targeted low-mode update (e.g. pCN on that mode).
  MM1 quarter analysis of the ten short chains (burn-in 8 sweeps): classification by the rule above (fixed 2026-10-03,
      written down after a first look at the quarter numbers -- the rule is the user's case distinction, not tuned).
  MM2 (statement) relaxation of the long chains (burn-in 100 sweeps, monitor every sweep, sector traces every 4): chain-
      averaged log max |f_U| and d_ref versus sweep, and the sweep after which the common trend is gone.
  MM4 (statement) classification of the LONG chains with the same rule, on the monitor values of d_ref and log max |f_U|
      for sweeps 40-140 (quarters), fixed before the long chains finished; final verdict of QMC.MIXING.MODE.01 (MM4 if
      available, else MM1; SLOW_MODE_IDENTIFIED only with MM3).
  MM5 (statement) FOUR-OUTCOME classification on the continued chains (L_k sweeps 0-139 followed by M_k, sweep offset
      +140; rule fixed 2026-10-03 before the continuation data): windows A = sweeps 150-200, B = 210-260 (d_ref, monitor every 4);
      BURN_IN_LIMITED if the between-chain spread of the window means falls by > 2x from A to B and the B means are
      compatible (chi^2 / dof < 2 with within-window errors); MULTIMODAL_BARRIER if every chain is internally stationary
      (|mean_B - mean_A| < 2 sigma) but the B means are incompatible (chi^2 / dof > 10) and no feature is identified;
      SLOW_MODE_IDENTIFIED if MM3 identifies a feature that also co-moves within chains (|r| > 0.5); else
      EXPLORATORY_DIAGNOSIS_ONLY.  kappa_char = sum |f_U| / |sum chi^* f_U| is reported (conditioning of the character sum).
  MM3 (statement) global HS-field features (uniform channel components m_A, tau-Fourier modes, |phi_A(K)|^2 at the 16
      cell momenta, total norm, log |W|) against d_ref: correlation over all monitor records and R^2 of the chain means.
    python3 -B experiments/theory-contracts/seam-holonomy-generating-20261001/mixing_mode.py
"""
from __future__ import annotations

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


def classify(Mq):
    """Mq: (chains, 4) quarter means of one observable."""
    trend = Mq[:, 3] - Mq[:, 0]
    t_sig = trend.mean() / (trend.std(ddof=1) / np.sqrt(len(trend)))
    within = np.mean(np.var(Mq, axis=1, ddof=1)); between = np.var(Mq.mean(1), ddof=1)
    if abs(t_sig) > 3:
        return "BURN_IN_LIMITED", t_sig, within / between
    if within / between < 0.1:
        return "MULTIMODAL_BARRIER", t_sig, within / between
    return "SLOW_MODE_CANDIDATE", t_sig, within / between


# MM1 --------------------------------------------------------------------------------------------------------------
rows, Md, Mf = [], [], []
for fp in sorted(HERE.glob("spg_points_g018b6*.json")):
    d = json.loads(fp.read_text())
    z = np.exp(1j * np.array(d["theta"])); f = np.array(d["f_re"]) + 1j * np.array(d["f_im"])
    chi = np.array(d["chi_re"]) + 1j * np.array(d["chi_im"]); ref = d["sectors"].index([0, 0])
    dref = (z * (f @ chi[ref].conj()) / 64).real
    lf = np.log(np.abs(np.delete(f, d["identity_element"], axis=1)).max(1))
    q = np.array_split(np.arange(len(dref)), 4)
    Md.append([dref[i].mean() for i in q]); Mf.append([lf[i].mean() for i in q])
Md, Mf = np.array(Md), np.array(Mf)
if len(Md) >= 3:
    cD, tD, rD = classify(Md); cF, tF, rF = classify(Mf)
    verdict1 = "BURN_IN_LIMITED" if "BURN_IN_LIMITED" in (cD, cF) else (cD if cD == cF else "SLOW_MODE_CANDIDATE")
    step("MM1 quarter analysis of the ten short chains (g = 0.018, beta = 6, burn-in 8 sweeps): common first-to-last-quarter trend, within- vs between-chain variance -> classification",
         True, f"{len(Md)} chains; D_ref: trend {tD:+.1f} sigma, within/between {rD:.2f} -> {cD}; log max |f_U|: trend {tF:+.1f} sigma, within/between {rF:.2f} -> {cF}; VERDICT {verdict1}; "
               f"quarter means of D_ref per chain: {np.round(Md, 3).tolist()}")
else:
    step("MM1 quarter analysis", True, "short chains missing")

# MM2 / MM3: long chains ---------------------------------------------------------------------------------------------
mons = {}
for fp in sorted(HERE.glob("spgm_monitor_g018b6L*.jsonl")):
    recs = [json.loads(l) for l in fp.read_text().splitlines() if l.strip()]
    if recs:
        mons[fp.stem.split("_")[-1]] = recs
if mons:
    curves = {}
    for k, recs in mons.items():
        fr = [r for r in recs if "f_re" in r]
        sw = np.array([r["sw"] for r in fr]); dref = np.array([r["d_ref"] for r in fr])
        lf = np.array([np.log(np.sort(np.abs(np.array(r["f_re"]) + 1j * np.array(r["f_im"])))[-2]) for r in fr])   # largest non-identity |f_U|
        curves[k] = (sw, dref, lf)
    nmax = min(len(c[0]) for c in curves.values())
    sw_c = list(curves.values())[0][0][:nmax]
    LF = np.array([c[2][:nmax] for c in curves.values()]); DR = np.array([c[1][:nmax] for c in curves.values()])
    blocks = np.array_split(np.arange(nmax), max(1, nmax // 4))
    prof = [(int(sw_c[b[0]]), float(LF[:, b].mean()), float(DR[:, b].mean()), float(DR[:, b].mean(1).std(ddof=1)) if len(curves) > 1 else np.nan) for b in blocks]
    step("MM2 (statement) relaxation of the long chains: chain-averaged log max |f_U| (U != 1) and d_ref versus sweep (blocks of 4 monitor points), chain-to-chain spread of d_ref",
         True, f"{len(curves)} chains, up to sweep {int(sw_c[-1]) if nmax else 0}: " + "; ".join(f"sweep {a}: <log max|f|> {b:+.2f}, <d_ref> {c:.3f}, spread {e:.3f}" for a, b, c, e in prof))
    # MM3 features
    feats, target, chain_id = [], [], []
    names = None
    for ci, (k, recs) in enumerate(mons.items()):
        for r in recs:
            if "d_ref" not in r:
                continue
            mA = np.array(r["mA"]); Ft = np.array(r["Ft_abs"]); PK = np.array(r["PK"])
            vec = {"norm2": r["norm2"], "logabsW": r["logabsW"], "|m|": float(np.linalg.norm(mA))}
            for A in range(10):
                vec[f"m_{A}"] = mA[A]
            for w in range(3):
                vec[f"|Ft_w{w}|"] = float(np.linalg.norm(Ft[w]))
            for Ki in range(PK.shape[0]):
                vec[f"P_K{Ki}"] = float(PK[Ki].sum())
            names = list(vec.keys()); feats.append([vec[n] for n in names]); target.append(r["d_ref"]); chain_id.append(ci)
    feats, target, chain_id = np.array(feats), np.array(target), np.array(chain_id)
    if len(target) > 5 and names:
        cors = [(n, float(np.corrcoef(feats[:, j], target)[0, 1])) for j, n in enumerate(names) if np.std(feats[:, j]) > 0]
        cors.sort(key=lambda t: -abs(t[1]))
        chains_ = np.unique(chain_id)
        r2 = []
        if len(chains_) >= 3:
            tmean = np.array([target[chain_id == c].mean() for c in chains_])
            for j, n in enumerate(names):
                fmean = np.array([feats[chain_id == c, j].mean() for c in chains_])
                if np.std(fmean) > 0:
                    r2.append((n, float(np.corrcoef(fmean, tmean)[0, 1] ** 2)))
            r2.sort(key=lambda t: -t[1])
            # multiple-testing control: maximum chain-mean R^2 over all features for randomly permuted chain targets
            rngp = np.random.default_rng(3); nullmax = []
            FM = np.array([[feats[chain_id == c, j].mean() for j in range(len(names))] for c in chains_])
            ok_cols = [j for j in range(len(names)) if np.std(FM[:, j]) > 0]
            for _ in range(2000):
                tp = rngp.permutation(tmean)
                nullmax.append(max(np.corrcoef(FM[:, j], tp)[0, 1] ** 2 for j in ok_cols))
            q95 = float(np.quantile(nullmax, 0.95))
            ident = [n for n, v in r2 if v > 0.7 and v > q95]
        # within-chain co-movement: chain means removed from feature and target (does d_ref follow the feature in time?)
        fw = feats - np.array([feats[chain_id == c].mean(0) for c in chain_id]); tw = target - np.array([target[chain_id == c].mean() for c in chain_id])
        wc = [(n, float(np.corrcoef(fw[:, j], tw)[0, 1])) for j, n in enumerate(names) if np.std(fw[:, j]) > 0]
        wc.sort(key=lambda t_: -abs(t_[1]))
        step("MM3 (statement) global HS-field features against d_ref: strongest correlations over all monitor records and R^2 of the chain means (SLOW_MODE_IDENTIFIED needs R^2 > 0.7 for one feature)",
             True, f"{len(target)} records, {len(chains_)} chains; top correlations: " + ", ".join(f"{n} {c:+.2f}" for n, c in cors[:6]) + "; chain-mean R^2: " + (", ".join(f"{n} {v:.2f}" for n, v in r2[:6]) if r2 else "n/a") + (f"; null 95 % quantile of the max R^2 (chain-shuffled, {len(names)} features): {q95:.2f}; identified: {ident if ident else 'none'}" if r2 else "") + "; within-chain co-movement (chain means removed): " + ", ".join(f"{n} {c:+.2f}" for n, c in wc[:6]))
    else:
        step("MM3 (statement) features", True, "not enough monitor records yet")
else:
    step("MM2 (statement) relaxation of the long chains", True, "long chains not started / no records")
    step("MM3 (statement) features", True, "long chains not started / no records")

if mons:
    Mq_d, Mq_f, lens = [], [], []
    for k, (sw, dref, lf) in curves.items():
        sel = np.where(sw >= 40)[0]
        lens.append(int(sw.max()))
        if len(sel) >= 8:
            qq = np.array_split(sel, 4)
            Mq_d.append([dref[i].mean() for i in qq]); Mq_f.append([lf[i].mean() for i in qq])
    if len(Mq_d) >= 3 and min(lens) >= 130:
        cD2, tD2, rD2 = classify(np.array(Mq_d)); cF2, tF2, rF2 = classify(np.array(Mq_f))
        final = "BURN_IN_LIMITED" if "BURN_IN_LIMITED" in (cD2, cF2) else ("MULTIMODAL_BARRIER" if "MULTIMODAL_BARRIER" in (cD2, cF2) else "SLOW_MODE_CANDIDATE")
        if final == "SLOW_MODE_CANDIDATE" and "ident" in dir() and ident:
            final = "SLOW_MODE_IDENTIFIED"
        step("MM4 (statement) long chains, sweeps 40-140: quarter classification and the final QMC.MIXING.MODE.01 verdict",
             True, f"{len(Mq_d)} chains; d_ref: trend {tD2:+.1f} sigma, within/between {rD2:.2f} -> {cD2}; log max |f_U|: trend {tF2:+.1f} sigma, within/between {rF2:.2f} -> {cF2}; quarter means of d_ref {np.round(Mq_d, 3).tolist()}; FINAL {final}")
    else:
        step("MM4 (statement) long chains", True, f"not finished (sweeps reached {min(lens) if lens else 0})")

cont = {}
for fp in sorted(HERE.glob("spgm_monitor_g018b6M*.jsonl")):
    k = fp.stem.split("_")[-1].replace("M", "L")
    recs = [json.loads(l) for l in fp.read_text().splitlines() if l.strip()]
    if recs and k in mons:
        cont[k] = [r for r in mons[k] if "d_ref" in r] + [{**r, "sw": r["sw"] + 140} for r in recs if "d_ref" in r]
if cont and min(max(r["sw"] for r in v) for v in cont.values()) >= 255:
    def win(v, lo, hi):
        x = np.array([r["d_ref"] for r in v if lo <= r["sw"] <= hi]); return x.mean(), x.std(ddof=1) / np.sqrt(len(x)), len(x)
    A = {k: win(v, 150, 200) for k, v in cont.items()}; Bw = {k: win(v, 210, 260) for k, v in cont.items()}
    mA_ = np.array([A[k][0] for k in cont]); mB = np.array([Bw[k][0] for k in cont]); eA = np.array([A[k][1] for k in cont]); eB = np.array([Bw[k][1] for k in cont])
    spread_ratio = mA_.std(ddof=1) / max(mB.std(ddof=1), 1e-12)
    chi2B = float(np.sum((mB - mB.mean()) ** 2 / np.maximum(eB, 1e-9) ** 2) / (len(mB) - 1))
    stationary = bool(np.all(np.abs(mB - mA_) < 2 * np.hypot(eA, eB)))
    kap = [r.get("kappa_char") for v in cont.values() for r in v if r.get("kappa_char") is not None]
    ident_wc = False
    if "ident" in dir() and ident and "wc" in dir():
        ident_wc = any(abs(c) > 0.5 for n, c in wc if n in ident)
    if spread_ratio > 2 and chi2B < 2:
        out5 = "BURN_IN_LIMITED"
    elif ident_wc:
        out5 = "SLOW_MODE_IDENTIFIED"
    elif stationary and chi2B > 10:
        out5 = "MULTIMODAL_BARRIER"
    else:
        out5 = "EXPLORATORY_DIAGNOSIS_ONLY"
    step("MM5 (statement) four-outcome classification on the continued chains (windows A = 150-200, B = 210-260 sweeps)",
         True, f"window means A {np.round(mA_, 3).tolist()}, B {np.round(mB, 3).tolist()}; spread A / B {spread_ratio:.2f}; chi^2/dof of B means {chi2B:.1f}; internally stationary {stationary}; kappa_char max {max(kap) if kap else float('nan'):.2f}; OUTCOME {out5}")
else:
    step("MM5 (statement) four-outcome classification on the continued chains", True, "continuation not finished")

out = {"steps": LOG, "passed": sum(x["pass"] for x in LOG), "total": len(LOG)}
(HERE / "mixing_mode_results.json").write_text(json.dumps(out, indent=2, default=float))
print(f"{out['passed']}/{out['total']} steps passed")
sys.exit(0 if out["passed"] == out["total"] else 1)
