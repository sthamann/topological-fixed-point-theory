#!/usr/bin/env python3
"""Part VT: CHIRAL4D.PAIRING.ORDER.VOLUME.01 -- is the high-amplitude state of Theta_1 ordered or symmetric (Gate B, G4)?

In equilibrium (QMC.EQUILIBRATION.AUDIT.01 with the global amplitude moves) the HS field at g = 0.021, beta = 4 sits at
|phi|^2 ~ 2.3 x prior with the time-constant pairing power enhanced 3-6 x over many momenta and channels, and the mirror
line is strongly suppressed (G_{P+d}(beta/2) = 0.039).  Ordered or symmetric?  The time-constant mode power per momentum
    P_A(K) = |sum_{x, l} e^{-i K x} phi_{l, A}(x)|^2 / (N L)
is intensive (N independent) for short-range pairing correlations and grows ~ N for long-range order at some K* (in units of
the prior 2 g dtau).  Same g, beta, lambda and sampler (local sweeps + global amplitude / pCN moves) on 8, 16 and 32 sites.
  VT3 (statement) G4 block decision (rule fixed 2026-10-05 before the 32-site extension data): per cell and chain, on
      the last two blocks of 25 sweeps: (A) stationarity -- |phi|^2, S_0/N, S_max(K != 0)/N and the channel anisotropy agree
      within 15 %; (B) seed consistency -- the two chains of a cell agree within 30 % on S_0/N and S_max/N; (C) R_eq =
      |<e^{i theta}>| on the stationary part only, with a block-jackknife lower bound R_eq - 2 sigma > 0.5.  All three on 8,
      16, 32 and S/N falling by > 1.6 per doubling -> G4_PASS_CANDIDATE; R_eq clearly < 0.5 on 32 -> SOLVER_LIMIT (no further
      sweeps; SOLVER.FEASIBILITY.01); otherwise CONTINUE.
  VT2 (statement) channel-resolved G4 per cell (same definitions on 8 / 16 / 32 sites): S_0/N at K = 0 (uniform pairing:
      U(1) / Spin(10) order parameter), S_max/N over K != 0 (finite-momentum pairing: translation breaking), channel
      anisotropy max_A / mean_A of the K-summed power (a preferred Spin(10) direction); the density channel is not
      monitored here (open).
  VT1 (statement) per cell: late-window |phi|^2 / prior, max_{K, A} P_A(K) / prior, the volume-normalised order parameter
      S_max / N = max_{K, A} P_A(K) / (N prior) (constant for long-range order, ~ 1/N for short range), the median, the max over K of the
      channel sum, G_{P+d}(beta/2) and R; ORDER_CANDIDATE if the maximal mode power grows ~ N between the cells (ratio of
      max power per site-doubling > 1.6 both times), SHORT_RANGE if it stays flat (< 1.3), else UNDECIDED.  Rule fixed
      2026-10-05 before the 8- and 32-site data; GUARD added after the first evaluation (which wrongly used 4-sweep 32-site
      chains): a cell enters only if each of its chains has >= 100 sweeps and its max-mode power is stationary (last two
      windows of 25 sweeps within 25 %); otherwise the verdict is INSUFFICIENT.
    python3 -B volume_test.py run <cell> <beta> <g> <lambda> <sweeps> <seed> <tag> [start.npy]   (env EQA_GLOBAL="eps,rho")
    python3 -B volume_test.py                                                                     (analysis)
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
_sp = ilu.spec_from_file_location("reweight", HERE / "reweight.py"); rw = ilu.module_from_spec(_sp); _sp.loader.exec_module(rw)


def step(name: str, ok: bool, detail: str = "") -> None:
    LOG.append({"step": name, "pass": bool(ok), "detail": detail})
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))
    sys.stdout.flush()


def cell_momenta(m):
    """all N momenta of the cell (distinct characters e^{-i K x} on the site vectors)."""
    X = np.array(m.sites, float); Rv = rw.CELLS[m.cell].astype(float)
    chars, Ks = [], []
    rng = range(0, 16)
    for nv in np.ndindex(16, 16, 16):
        K = 2 * np.pi * np.linalg.solve(Rv, np.array(nv, float))
        ch = np.round(np.exp(-1j * X @ K), 8)
        key = tuple(np.round(ch.real, 6)) + tuple(np.round(ch.imag, 6))
        if key not in chars:
            chars.append(key); Ks.append(K)
            if len(Ks) == m.N:
                break
    return np.array(Ks)


if len(sys.argv) > 1 and sys.argv[1] == "run":
    cell, beta, g_, lam, sweeps, seed, tag = sys.argv[2], float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]), int(sys.argv[6]), int(sys.argv[7]), sys.argv[8]
    start = np.load(sys.argv[9]) if len(sys.argv) > 9 else None
    twist = np.array([lam] * 3)
    gm = None
    if os.environ.get("EQA_GLOBAL"):
        vals = [float(x) for x in os.environ["EQA_GLOBAL"].split(",")]   # "eps,rho[,rho_uniform]"
        gm = {"scale_eps": vals[0], "pcn_rho": vals[1], "pcn_uniform_rho": vals[2] if len(vals) > 2 else 0.0,
              "uniform_scale_eps": vals[3] if len(vals) > 3 else 0.0, "every": 1}
    g_every = int(os.environ.get("VT_G_EVERY", "4"))
    mon = HERE / f"vt_monitor_{tag}.jsonl"
    offset = int(os.environ.get("VT_OFFSET", "0"))
    if offset == 0:
        mon.write_text("")
    holder = {}
    def monitor(m, phi, state):
        if "K" not in holder:
            holder["K"] = cell_momenta(m); holder["ph"] = np.exp(-1j * holder["K"] @ np.array(m.sites, float).T)
        Phi = phi.sum(axis=0)                                            # (N, 10)
        PK = np.abs(holder["ph"] @ Phi) ** 2 / (m.N * m.Lt)               # (N_K, 10)
        rec = {"sw": state["sw"] + offset, "theta": state["theta"], "logabsW": state["logabsW"], "norm2": float(np.sum(phi ** 2)), "PK": PK.tolist()}
        if state["sw"] % g_every == 0:
            G = m.greens_k(phi, m.momenta())
            rec["G_P+d"] = np.real(G["P+d"]).tolist(); rec["G_P-d"] = np.real(G["P-d"]).tolist()
        with mon.open("a") as fh:
            fh.write(json.dumps(rec) + "\n")
        return {}
    t0 = time.time()
    res = rw.run(cell, beta, 0.25, g_, sweeps, seed, "md", False, twist, None, burn=20, monitor=monitor, monitor_every=1,
                 save_phi=str(HERE / f"vt_phi_{tag}.npy"), phi0=start, global_moves=gm)
    meta = {"cell": cell, "beta": beta, "g": g_, "lambda": lam, "sweeps": sweeps, "seed": seed, "tag": tag, "R": res["R"], "acceptance": res["acceptance"],
            "global_acceptance": res.get("global_acceptance"), "seconds": time.time() - t0, "N": res["N"]}
    with (HERE / "vt_runs.jsonl").open("a") as fh:
        fh.write(json.dumps(meta) + "\n")
    print(json.dumps(meta))
    sys.exit(0)


# ------------------------------------------------------------------------------------------------ analysis
def late_stats(fp, g, nlast=40):
    rr = [json.loads(l) for l in fp.read_text().splitlines() if l.strip()]
    if not rr:
        return None
    n = max(x["sw"] for x in rr) + 1
    late = [x for x in rr if x["sw"] >= n - nlast]
    prior = 2 * g * 0.25
    PK = np.array([x["PK"] for x in late]).mean(0) / prior
    Nsites = PK.shape[0]; L = None
    G = [(x["theta"], x["G_P+d"]) for x in late if "G_P+d" in x]
    if G:
        th = np.array([t for t, _ in G]); L = len(G[0][1]) - 1
        Gh = float((np.array([v[L // 2] for _, v in G]) * np.cos(th)).mean() / np.cos(th).mean())
    else:
        Gh = np.nan
    Dn = (L or 16) * Nsites * 10
    Kmom = np.array([x["PK"] for x in late])
    k0 = 0
    S0 = float(PK[k0].max()); Sfin = float(np.delete(PK, k0, axis=0).max()); chan = PK.sum(0); aniso = float(chan.max() / chan.mean())
    mx_w = [np.mean([np.max(x["PK"]) for x in rr if a <= x["sw"] < a + 25] or [np.nan]) for a in (n - 50, n - 25)]
    stationary = bool(n >= 100 and np.isfinite(mx_w).all() and abs(mx_w[1] - mx_w[0]) <= 0.25 * mx_w[0])
    return {"S0": S0, "Sfin": Sfin, "aniso": aniso, "stationary": stationary, "n": n, "N": Nsites, "norm_ratio": float(np.mean([x["norm2"] for x in late]) / (Dn * prior)), "max_mode": float(PK.max()), "median_mode": float(np.median(PK)),
            "max_K_sum": float(PK.sum(1).max() / 10), "G_half": Gh, "R": float(abs(np.exp(1j * np.array([x["theta"] for x in late])).mean()))}
rows, per_cell = [], {}
for fp in sorted(HERE.glob("vt_monitor_*.jsonl")) + sorted(HERE.glob("eqa_monitor_g021b4G*.jsonl")):
    st = late_stats(fp, 0.021)
    if st is None:
        continue
    tag = fp.stem.split("_")[-1]
    if st["stationary"]:
        per_cell.setdefault(st["N"], []).append(st)
    rows.append(f"{tag} (N = {st['N']}, {st['n']} sweeps, max-mode stationary {st['stationary']}): |phi|^2/prior {st['norm_ratio']:.2f}, max mode {st['max_mode']:.1f}, S_max/N {st['max_mode'] / st['N']:.2f}, median {st['median_mode']:.2f}, max over K of the channel mean {st['max_K_sum']:.2f}, G(beta/2) {st['G_half']:.3f}, R {st['R']:.2f}")
verdict = "INSUFFICIENT"
Ns = sorted(per_cell)
if len(Ns) >= 3:
    mx = [np.mean([s["max_mode"] for s in per_cell[N]]) for N in Ns]
    ratios = [mx[i + 1] / mx[i] for i in range(len(mx) - 1)]
    verdict = "ORDER_CANDIDATE" if all(r_ > 1.6 for r_ in ratios) else ("SHORT_RANGE" if all(r_ < 1.3 for r_ in ratios) else "UNDECIDED")
    rows.append(f"max-mode power per cell {dict(zip(Ns, np.round(mx, 1)))}, S_max/N per cell {dict(zip(Ns, np.round(np.array(mx) / np.array(Ns), 2)))}, ratios per doubling {np.round(ratios, 2).tolist()} -> {verdict}")
step("VT1 (statement) pairing-order volume test at g = 0.021, beta = 4, lambda = 0.35 (8 / 16 / 32 sites, global moves): intensive vs extensive maximal time-constant mode power",
     True, " | ".join(rows) if rows else "volume runs not started")

rows2 = []
for N_ in sorted({st["N"] for st in [late_stats(fp, 0.021) for fp in sorted(HERE.glob("vt_monitor_*.jsonl")) + sorted(HERE.glob("eqa_monitor_g021b4G*.jsonl"))] if st}):
    sts = [st for st in [late_stats(fp, 0.021) for fp in sorted(HERE.glob("vt_monitor_*.jsonl")) + sorted(HERE.glob("eqa_monitor_g021b4G*.jsonl"))] if st and st["N"] == N_]
    rows2.append(f"N = {N_}: S_0/N " + "/".join(f"{st['S0'] / N_:.2f}" for st in sts) + ", S_max(K != 0)/N " + "/".join(f"{st['Sfin'] / N_:.2f}" for st in sts)
                 + ", channel anisotropy " + "/".join(f"{st['aniso']:.2f}" for st in sts) + " (stationary " + "/".join(str(st["stationary"]) for st in sts) + ")")
step("VT2 (statement) channel-resolved G4: uniform pairing S_0/N, finite-momentum pairing S_max/N, Spin(10) channel anisotropy per cell (density channel not monitored)",
     True, " | ".join(rows2) if rows2 else "no data")

def block_stats(fp, g):
    rr = [json.loads(l) for l in fp.read_text().splitlines() if l.strip()]
    if not rr:
        return None
    n = max(x["sw"] for x in rr) + 1; prior = 2 * g * 0.25
    def blk(lo, hi):
        b = [x for x in rr if lo <= x["sw"] < hi]
        PK = np.array([x["PK"] for x in b]).mean(0) / prior; Nn = PK.shape[0]
        return {"norm": float(np.mean([x["norm2"] for x in b])), "S0": float(PK[0].max() / Nn), "Sf": float(np.delete(PK, 0, axis=0).max() / Nn),
                "an": float(PK.sum(0).max() / PK.sum(0).mean()), "N": Nn}
    b1, b2 = blk(n - 50, n - 25), blk(n - 25, n)
    th = np.array([x["theta"] for x in rr if x["sw"] >= n - 50]); z = np.exp(1j * th)
    R = float(abs(z.mean())); parts = np.array_split(np.arange(len(z)), 5)
    Rj = [abs(np.delete(z, p_).mean()) for p_ in parts]; Re = float(np.sqrt(4 / 5 * np.sum((np.array(Rj) - np.mean(Rj)) ** 2)))
    stat = all(abs(b2[k] - b1[k]) <= 0.15 * abs(b1[k]) for k in ("norm", "S0", "Sf", "an"))
    return {"n": n, "N": b2["N"], "b": b2, "stationary": stat, "R": R, "R_err": Re}
cells = {}
for fp in sorted(HERE.glob("vt_monitor_*.jsonl")) + sorted(HERE.glob("eqa_monitor_g021b4G*.jsonl")):
    st = block_stats(fp, 0.021)
    if st and st["n"] >= 100:
        cells.setdefault(st["N"], []).append((fp.stem.split("_")[-1], st))
rows3, okA, okB, okC, Rlow = [], True, True, True, False
for N_, lst in sorted(cells.items()):
    sA = all(st["stationary"] for _, st in lst)
    sB = len(lst) >= 2 and all(abs(lst[0][1]["b"][k] - lst[1][1]["b"][k]) <= 0.30 * max(lst[0][1]["b"][k], lst[1][1]["b"][k]) for k in ("S0", "Sf"))
    sC = all(st["R"] - 2 * st["R_err"] > 0.5 for _, st in lst)
    okA &= sA; okB &= sB; okC &= sC
    if N_ == 32 and all(st["R"] + 2 * st["R_err"] < 0.5 for _, st in lst):
        Rlow = True
    rows3.append(f"N = {N_}: " + "; ".join(f"{t}: S0/N {st['b']['S0']:.2f}, Sf/N {st['b']['Sf']:.2f}, aniso {st['b']['an']:.2f}, R_eq {st['R']:.2f} +- {st['R_err']:.2f}, stationary {st['stationary']}" for t, st in lst)
                 + f" -> A {sA}, B {sB}, C {sC}")
Nsorted = sorted(cells)
falling = len(Nsorted) >= 3 and all(np.mean([st["b"]["Sf"] for _, st in cells[Nsorted[i]]]) / np.mean([st["b"]["Sf"] for _, st in cells[Nsorted[i + 1]]]) > 1.6 for i in range(len(Nsorted) - 1))
dec = "SOLVER_LIMIT" if Rlow else ("G4_PASS_CANDIDATE" if (okA and okB and okC and falling and len(Nsorted) >= 3) else "CONTINUE")
step("VT3 (statement) G4 block decision: (A) stationarity, (B) seed consistency, (C) R_eq on the stationary part, S/N falling per doubling",
     True, " | ".join(rows3) + f" || S_max(K != 0)/N falling > 1.6 per doubling: {falling} -> DECISION {dec}")

out = {"steps": LOG, "passed": sum(x["pass"] for x in LOG), "total": len(LOG), "verdict": verdict, "g4_decision": dec}
(HERE / "volume_test_results.json").write_text(json.dumps(out, indent=2, default=float))
print(f"{out['passed']}/{out['total']} steps passed")
sys.exit(0 if out["passed"] == out["total"] else 1)
