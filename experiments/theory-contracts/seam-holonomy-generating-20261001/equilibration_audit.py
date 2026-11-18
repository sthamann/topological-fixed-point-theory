#!/usr/bin/env python3
"""Part EA: QMC.EQUILIBRATION.AUDIT.01 -- which g >= 0.015 QMC points are actually in equilibrium?

Round 56 showed that 8-10 burn-in sweeps are far too few at g >= 0.018 (G_{P+d}(beta/2) at g = 0.018, beta = 6 drifted from
0.46 to 0.39 with the burn-in).  Every audit point (g, beta, lambda) gets >= 2 chains with DIFFERENT starts (A: fresh prior
draw; B: a configuration from an equilibrated nearby point), the unchanged sampler (step adapted during the first 20 sweeps,
then fixed), and per-sweep monitors: phase, log |W|, field low modes (uniform channel components, tau-Fourier modes, |phi(K)|^2,
alternating imaginary-time mode, norm), every 4 sweeps the full G(tau) at P +- delta, every 8 sweeps the 64 sector traces
(D_ref).  Chains are extended in blocks of 50 sweeps until the analysis releases the point.
  EQUILIBRATED (rule fixed 2026-10-04 before the audit data): with windows of 50 sweeps, (i) no late drift: for every chain
  |mean(W_last) - mean(W_prev)| < 2 sigma for G(beta/2) and log |W|; (ii) chains compatible: chi^2 / dof < 2 of the last-window
  means of G(beta/2) (Gamma-method errors); (iii) last-window phase R > 0.5 in every chain; (iv) no chain dominates: dropping
  any one chain moves the pooled G(beta/2) by < 2 sigma; (v) >= 150 sweeps per chain.  Otherwise NOT_EQUILIBRATED (extend).
  EA1 (statement) per audit point: window means per chain, the rule (i)-(v), the verdict, and for released points the
      equilibrium G(tau) (phase-reweighted, pooled over chains and the last window) next to the earlier short-burn-in value.
    python3 -B equilibration_audit.py run <beta> <g> <lambda> <sweeps> <seed> <tag> [start_phi.npy]
    python3 -B equilibration_audit.py                      (analysis)
"""
from __future__ import annotations

import importlib.util as ilu
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
LOG: list[dict] = []
_sp = ilu.spec_from_file_location("reweight", HERE / "reweight.py"); rw = ilu.module_from_spec(_sp); _sp.loader.exec_module(rw)


def step(name: str, ok: bool, detail: str = "") -> None:
    LOG.append({"step": name, "pass": bool(ok), "detail": detail})
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))
    sys.stdout.flush()


def load_spg_parts():
    src = (HERE / "spg_production.py").read_text()
    ns = {"__name__": "spg_parts", "__file__": str(HERE / "spg_production.py")}
    pre = src.split("# ------------------------------------------------------------------------------------------------ SPP1 brute force")[0]
    mid = src.split("# ------------------------------------------------------------------------------------------------ production operator")[1].split("# ------------------------------------------------------------------------------------------------ CLI: production run")[0]
    feat = "def field_features" + src.split("def field_features")[1].split('if len(sys.argv) > 1 and sys.argv[1] == "runlong":')[0]
    exec(compile(pre.replace('if len(sys.argv) > 1', 'if False'), "spg_pre", "exec"), ns)
    exec(compile(mid, "spg_mid", "exec"), ns); exec(compile(feat, "spg_feat", "exec"), ns)
    return ns


if len(sys.argv) > 1 and sys.argv[1] == "run":
    beta, g_, lam, sweeps, seed, tag = float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6]), sys.argv[7]
    start = np.load(sys.argv[8]) if len(sys.argv) > 8 else None
    import os
    ns = load_spg_parts()
    twist = np.array([lam] * 3)
    m0 = rw.Model("cubic16", beta, 0.25, g_, "md", twist)
    if start is not None and start.shape[0] != m0.Lt:                  # nearby point with another beta: take the first Lt slices
        start = start[:m0.Lt]
    holder = {}
    mon = HERE / f"eqa_monitor_{tag}.jsonl"
    if not (start is not None and len(sys.argv) > 9 and sys.argv[9] == "append"):
        mon.write_text("")
    offset = int(sys.argv[10]) if len(sys.argv) > 10 else 0
    def monitor(m, phi, state):
        if "sec" not in holder:
            holder["sec"] = ns["Sector"](m); holder["ref"] = holder["sec"].sectors.index((0, 0))
        sec = holder["sec"]
        rec = {"sw": state["sw"] + offset, "theta": state["theta"], "logabsW": state["logabsW"], **ns["field_features"](sec, phi)}
        if state["sw"] % 4 == 0:
            G = m.greens_k(phi, m.momenta())
            rec["G_P+d"] = np.real(G["P+d"]).tolist(); rec["G_P-d"] = np.real(G["P-d"]).tolist()
            rec["G_imag_max"] = float(max(np.abs(np.imag(G["P+d"])).max(), np.abs(np.imag(G["P-d"])).max()))
        if state["sw"] % 8 == 0:
            f = ns["f_only"](sec, phi)
            rec["d_ref"] = complex(np.exp(1j * state["theta"]) * (f @ sec.chi[holder["ref"]].conj()) / 64).real
            rec["kappa_char"] = float(np.sum(np.abs(f)) / max(abs(f @ sec.chi[holder["ref"]].conj()), 1e-300))
            if os.environ.get("EQA_FULL_F"):                               # all 64 sector traces (active-sector preparation)
                rec["f_re"] = f.real.tolist(); rec["f_im"] = f.imag.tolist()
        with mon.open("a") as fh:
            fh.write(json.dumps(rec) + "\n")
        return {}
    t0 = time.time()
    import os
    gm = None
    if os.environ.get("EQA_GLOBAL"):                                    # "eps,rho": targeted amplitude / low-mode update (global_moves.py)
        vals = [float(x) for x in os.environ["EQA_GLOBAL"].split(",")]   # "eps,rho[,rho_uniform]"
        gm = {"scale_eps": vals[0], "pcn_rho": vals[1], "pcn_uniform_rho": vals[2] if len(vals) > 2 else 0.0,
              "uniform_scale_eps": vals[3] if len(vals) > 3 else 0.0, "every": 1}
    res = rw.run("cubic16", beta, 0.25, g_, sweeps, seed, "md", False, twist, None, burn=20, monitor=monitor, monitor_every=1,
                 save_phi=str(HERE / f"eqa_phi_{tag}.npy"), phi0=start, global_moves=gm)
    meta = {k: res[k] for k in ("beta", "g", "sweeps", "seed", "R", "R_err", "acceptance")}
    meta.update({"lambda": lam, "tag": tag, "seconds": time.time() - t0, "start": sys.argv[8] if len(sys.argv) > 8 else "prior", "offset": offset,
                 "global_moves": gm, "global_acceptance": res.get("global_acceptance")})
    with (HERE / "eqa_runs.jsonl").open("a") as fh:
        fh.write(json.dumps(meta) + "\n")
    print(json.dumps(meta))
    sys.exit(0)


# ------------------------------------------------------------------------------------------------ analysis
_sb = (HERE / "spectral_bounds.py").read_text()
exec(compile(_sb.split("# ---------------------------------------------------------------------------------------------- Gamma method")[1].split("# SB1 ---")[0], "gamma", "exec"))
POINTS = {"g015b6": (0.015, 6.0, 0.32), "g018b4": (0.018, 4.0, 0.39), "g021b4": (0.021, 4.0, 0.35), "g021b6": (0.021, 6.0, 0.35),
          "g024b4": (0.024, 4.0, 0.37), "g018b2": (0.018, 2.0, 0.39), "g018b6": (0.018, 6.0, 0.39),
          "g015b6G": (0.015, 6.0, 0.32), "g021b4G": (0.021, 4.0, 0.35), "g018b4G": (0.018, 4.0, 0.39)}   # ...G: chains with the global amplitude / pCN moves, evaluated as their own ensemble
SHORT = {"g018b4": ("reweight_points_md_cubic16_gl4c.json", "P+d"), "g021b4": ("reweight_points_md_cubic16_gl4d.json", "P+d"),
         "g021b6": ("reweight_points_md_cubic16_gl6a.json", "P+d"), "g024b4": ("reweight_points_md_cubic16_g24c.json", "P+d"),
         "g018b2": ("reweight_points_md_cubic16_ds3.json", "P+d")}
rows = []
for key, (g_, beta, lam) in POINTS.items():
    files = [fp for fp in sorted(HERE.glob(f"eqa_monitor_{key}*.jsonl")) if fp.stem.split("_")[-1][len(key):] in ("A", "B", "C", "D")]   # plain and global-move chains kept apart
    chains = {}
    for fp in files:
        recs = [json.loads(l) for l in fp.read_text().splitlines() if l.strip()]
        if recs:
            chains[fp.stem.split("_")[-1]] = recs
    if len(chains) < 1:
        continue
    L = int(round(beta / 0.25)); half = L // 2
    short_key = key.rstrip("G")
    per = {}
    for cname, recs in chains.items():
        Gr = [(r["sw"], r["theta"], r["G_P+d"][half]) for r in recs if "G_P+d" in r]
        Wr = [(r["sw"], r["logabsW"], r["theta"]) for r in recs]
        nsw = max(r["sw"] for r in recs) + 1
        per[cname] = (np.array(Gr), np.array(Wr), nsw, recs)
    nmin = min(v[2] for v in per.values())
    def window(arr, lo, hi, col):
        sel = (arr[:, 0] >= lo) & (arr[:, 0] < hi)
        return arr[sel, col], arr[sel, 1] if col != 1 else None
    w_last = (nmin - 50, nmin); w_prev = (nmin - 100, nmin - 50)
    crit = {"i": True, "ii": None, "iii": True, "iv": None, "v": nmin >= 150}
    lastG, lastE, lines = [], [], []
    for cname, (Gr, Wr, nsw, recs) in per.items():
        def wstats(lo, hi):
            sG = (Gr[:, 0] >= lo) & (Gr[:, 0] < hi); sW = (Wr[:, 0] >= lo) & (Wr[:, 0] < hi)
            if sG.sum() < 3 or sW.sum() < 3:
                return None
            th = Gr[sG, 1]; Gv = Gr[sG, 2]
            gm, ge, _ = uwerr(np.vstack([Gv * np.cos(th), np.cos(th)]), f=lambda mm: mm[0] / mm[1])
            wm, we, _ = uwerr(Wr[sW, 1])
            R = float(abs(np.exp(1j * Wr[sW, 2]).mean()))
            return gm, ge, wm, we, R
        a, b = wstats(*w_prev), wstats(*w_last)
        if a is None or b is None:
            crit["i"] = False; continue
        drift_G = abs(b[0] - a[0]) / max(np.hypot(a[1], b[1]), 1e-9); drift_W = abs(b[2] - a[2]) / max(np.hypot(a[3], b[3]), 1e-9)
        if drift_G >= 2 or drift_W >= 2:
            crit["i"] = False
        if b[4] <= 0.5:
            crit["iii"] = False
        lastG.append(b[0]); lastE.append(b[1])
        dr = [r["d_ref"] for r in recs if "d_ref" in r and r["sw"] >= w_last[0]]
        lines.append(f"{cname}: G(beta/2) prev {a[0]:.3f}+-{a[1]:.3f} -> last {b[0]:.3f}+-{b[1]:.3f} (drift {drift_G:.1f} sigma), log|W| drift {drift_W:.1f} sigma, R_last {b[4]:.2f}, d_ref_last {np.mean(dr) if dr else float('nan'):.3f}")
    if len(lastG) >= 2:
        lastG, lastE = np.array(lastG), np.array(lastE)
        chi2 = float(np.sum((lastG - lastG.mean()) ** 2 / np.maximum(lastE, 1e-9) ** 2) / (len(lastG) - 1)); crit["ii"] = chi2 < 2
        pooled = lastG.mean(); perr = max(lastG.std(ddof=1) / np.sqrt(len(lastG)), np.sqrt(np.sum(lastE ** 2)) / len(lastG))
        crit["iv"] = all(abs(np.delete(lastG, i).mean() - pooled) < 2 * perr for i in range(len(lastG)))
    else:
        chi2, pooled, perr = np.nan, (lastG[0] if lastG else np.nan), (lastE[0] if lastE else np.nan)
    verdict = "EQUILIBRATED" if all(v is True for v in crit.values()) and len(per) >= 2 else "NOT_EQUILIBRATED"
    short = ""
    if short_key in SHORT and (HERE / SHORT[short_key][0]).exists():
        sp = json.loads((HERE / SHORT[short_key][0]).read_text())[-1]["greens"][SHORT[short_key][1]]
        short = f"; earlier short-burn-in G(beta/2) {sp['G_half']:.3f}"
    rows.append(f"{key} (g = {g_}, beta = {beta:g}, lambda = {lam}; {len(per)} chains, {nmin} sweeps): {verdict} -- criteria {crit}, chi^2/dof {chi2:.2f}, pooled last-window G(beta/2) {pooled:.3f} +- {perr:.3f}{short} || " + " | ".join(lines))
step("EA1 (statement) QMC.EQUILIBRATION.AUDIT.01 per audit point: window rule (i)-(v), verdict, pooled last-window G(beta/2) and the earlier short-burn-in value",
     True, " #### ".join(rows) if rows else "audit runs not started")

out = {"steps": LOG, "passed": sum(x["pass"] for x in LOG), "total": len(LOG)}
(HERE / "equilibration_audit_results.json").write_text(json.dumps(out, indent=2, default=float))
print(f"{out['passed']}/{out['total']} steps passed")
sys.exit(0 if out["passed"] == out["total"] else 1)
