#!/usr/bin/env python3
"""Part PT: QMC.TEMPERING.01 -- parallel tempering in the fermion-weight power for the two-branch point g = 0.015, beta = 6.

QMC.EQUILIBRATION.AUDIT.01 found two branches at g = 0.015, beta = 6 that do not meet in 500 local sweeps nor with the global
amplitude move (G_{P+d}(beta/2) 0.32 at |phi|^2 ~ 36-37 vs 0.20 at ~ 43-44).  Ladder p_s(phi) ~ e^{-S_HS} |W(phi)|^s, s = s_0 < ...
< s_{K-1} = 1: s = 0 is the Gaussian prior (no barrier), s = 1 the physical phase-quenched measure.  Per round every replica
does one local sweep (acceptance -dS + s dlog|W|) and one global amplitude move (Jacobian s^D, weight |W|^s); then neighbour
swaps with min(1, exp((s_i - s_j)(l_j - l_i))), l = log |W|.  Only s = 1 is measured; its phase is continued along its own
moves and rebuilt from phi = 0 after it receives a configuration by a swap.
  PT1 toy (bimodal 1D target e^{-x^2/2} w(x)^s with two well separated wells): plain Metropolis at s = 1 stays in its
      start well, the tempered ladder reproduces the exact well weights at s = 1.
  PT2 (statement) g = 0.015, beta = 6, lambda = 0.32: swap rates, round trips of the replica labels, the s = 1 history of
      |phi|^2 and G(beta/2), and which branch carries the weight.
    python3 -B tempering.py run <beta> <g> <lambda> <rounds> <seed> <tag> <s0,s1,...> <start_0.npy|prior> ...
    python3 -B tempering.py                            (toy check + analysis)
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
wrap = lambda a: (a + np.pi) % (2 * np.pi) - np.pi


def step(name: str, ok: bool, detail: str = "") -> None:
    LOG.append({"step": name, "pass": bool(ok), "detail": detail})
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))
    sys.stdout.flush()


if len(sys.argv) > 1 and sys.argv[1] == "run":
    _sp = ilu.spec_from_file_location("reweight", HERE / "reweight.py"); rw = ilu.module_from_spec(_sp); _sp.loader.exec_module(rw)
    beta, g_, lam, rounds, seed, tag = float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6]), sys.argv[7]
    svals = [float(x) for x in sys.argv[8].split(",")]; starts = sys.argv[9:9 + len(svals)]
    dtau = 0.25; m = rw.Model("cubic16", beta, dtau, g_, "md", np.array([lam] * 3))
    Lt, N = m.Lt, m.N; nn = 2 * m.n; unit = (np.eye(nn, dtype=complex), np.ones(nn), np.eye(nn, dtype=complex))
    rng = np.random.default_rng(seed); sig = np.sqrt(2 * g_ * dtau)
    reps = []
    for i, s in enumerate(svals):
        phi = np.load(starts[i]).reshape(Lt, N, 10).copy() if starts[i] != "prior" else sig * rng.normal(size=(Lt, N, 10))
        S = [m.slice_op(phi[l]) for l in range(Lt)]
        la, ph, _ = m.logdet(S)
        reps.append({"phi": phi, "S": S, "la": la, "ph": ph, "label": i, "step": 0.5 * sig / np.sqrt(N)})
    def track_from_zero(phi):
        """phase of W(phi) continued along t phi from t = 0 (W > 0 there)."""
        t_, dt, th, prev = 0.0, 0.125, 0.0, m.logdet([m.slice_op(0 * phi[l]) for l in range(Lt)])[1]
        while t_ < 1 - 1e-12:
            dt = min(dt, 1 - t_); cur = m.logdet([m.slice_op((t_ + dt) * phi[l]) for l in range(Lt)])[1]
            inc = wrap(cur - prev)
            if abs(inc) > np.pi / 4 and dt > 1e-4:
                dt /= 2; continue
            th += inc; prev = cur; t_ += dt; dt *= 1.5
        return th, prev
    top = len(svals) - 1
    theta, ph_ref = track_from_zero(reps[top]["phi"])                    # tracked phase (2 theta = arg det) of the s = 1 replica
    mon = HERE / f"pt_monitor_{tag}.jsonl"; mon.write_text("")
    acc_sw = np.zeros(len(svals) - 1); try_sw = np.zeros(len(svals) - 1)
    t0 = time.time()
    for rd in range(rounds):
        for i, rep in enumerate(reps):
            s = svals[i]; phi, S = rep["phi"], rep["S"]
            Y = m.sweep_cache(S); X = unit
            la_s, ph_s = m.logdet_split(S[0], unit, Y[0])
            if i == top:
                theta += wrap(ph_s - rep["ph"])
            rep["la"], rep["ph"] = la_s, ph_s
            acc = 0
            for l in range(Lt):
                new = phi[l] + rep["step"] * rng.normal(size=(N, 10))
                dS = (np.sum(new ** 2) - np.sum(phi[l] ** 2)) / (4 * g_ * dtau)
                Sl = m.slice_op(new)
                la2, ph2 = m.logdet_split(Sl, X, Y[l])
                inc = wrap(ph2 - rep["ph"])
                if i == top and abs(inc) > np.pi / 2:                        # continuity refinement only where the phase is used
                    inc_tot, ph_cur = 0.0, rep["ph"]
                    for k in range(1, 9):
                        phm = m.logdet_split(m.slice_op(phi[l] + (new - phi[l]) * k / 8), X, Y[l])[1]
                        inc_tot += wrap(phm - ph_cur); ph_cur = phm
                    inc = inc_tot
                if np.log(rng.random()) < -dS + s * 0.5 * (la2 - rep["la"]):
                    phi[l] = new; S[l] = Sl; rep["la"], rep["ph"] = la2, ph2; acc += 1
                    if i == top:
                        theta += inc
                X = m._udv_left(S[l].T, X)
            if rd < 10:                                                     # adapt the local step during the first rounds only
                rep["step"] *= 1.3 if acc / Lt > 0.6 else (0.7 if acc / Lt < 0.35 else 1.0)
            # global amplitude move (Jacobian s^D, tempered weight)
            sf = np.exp(0.02 * rng.normal()); n2 = float(np.sum(phi ** 2))
            Snew = [m.slice_op(sf * phi[l]) for l in range(Lt)]
            la_n, ph_n, _ = m.logdet(Snew)
            lw = phi.size * np.log(sf) - (sf ** 2 - 1) * n2 / (4 * g_ * dtau) + s * 0.5 * (la_n - rep["la"])
            if np.log(rng.random()) < lw:
                if i == top:                                                # phase along the straight path
                    t_, dt, inc_tot, prev = 0.0, 0.5, 0.0, rep["ph"]
                    while t_ < 1 - 1e-12:
                        dt = min(dt, 1 - t_); cur = m.logdet([m.slice_op((1 + (t_ + dt) * (sf - 1)) * phi[l]) for l in range(Lt)])[1]
                        inc = wrap(cur - prev)
                        if abs(inc) > np.pi / 4 and dt > 1e-4:
                            dt /= 2; continue
                        inc_tot += inc; prev = cur; t_ += dt; dt *= 1.5
                    theta += inc_tot
                rep["phi"] = sf * phi; rep["S"] = Snew; rep["la"], rep["ph"] = la_n, ph_n
        # neighbour swaps (alternating even / odd pairs)
        swapped_top = False
        for i in range(rd % 2, len(svals) - 1, 2):
            li, lj = 0.5 * reps[i]["la"], 0.5 * reps[i + 1]["la"]
            try_sw[i] += 1
            if np.log(rng.random()) < (svals[i] - svals[i + 1]) * (lj - li):
                for key in ("phi", "S", "la", "ph", "label"):
                    reps[i][key], reps[i + 1][key] = reps[i + 1][key], reps[i][key]
                acc_sw[i] += 1
                if i + 1 == top:
                    swapped_top = True
        if swapped_top:
            theta, ph_ref = track_from_zero(reps[top]["phi"]); reps[top]["ph"] = ph_ref
        rec = {"round": rd, "labels": [r_["label"] for r_ in reps], "norm2": [float(np.sum(r_["phi"] ** 2)) for r_ in reps],
               "logabsW": [0.5 * r_["la"] for r_ in reps], "theta_top": theta / 2, "swap_acc": (acc_sw / np.maximum(try_sw, 1)).tolist()}
        if rd % 4 == 0:
            G = m.greens_k(reps[top]["phi"], m.momenta())
            rec["G_P+d"] = np.real(G["P+d"]).tolist()
        with mon.open("a") as fh:
            fh.write(json.dumps(rec) + "\n")
    for i, rep in enumerate(reps):
        np.save(HERE / f"pt_phi_{tag}_s{i}.npy", rep["phi"])
    print(json.dumps({"tag": tag, "rounds": rounds, "seconds": time.time() - t0, "swap_acc": (acc_sw / np.maximum(try_sw, 1)).tolist()}))
    sys.exit(0)


# PT1 toy ------------------------------------------------------------------------------------------------------------
rng = np.random.default_rng(9)
logw = lambda x: np.logaddexp(-((x - 4.0) ** 2) / 0.08 + np.log(0.3), -((x + 4.0) ** 2) / 0.08 + np.log(0.7))   # two narrow wells (weights 0.3 / 0.7)
logprior = lambda x: -x ** 2 / 50.0
xs = np.linspace(-12, 12, 200001); dens = np.exp(logprior(xs) + logw(xs)); exact_right = float(dens[xs > 0].sum() / dens.sum())
def pt_toy(n=40000, ladder=(0.0, 0.1, 0.25, 0.5, 1.0)):
    x = [4.0] * len(ladder); right = []
    for it in range(n):
        for i, s in enumerate(ladder):
            y = x[i] + 0.3 * rng.normal()
            if np.log(rng.random()) < logprior(y) - logprior(x[i]) + s * (logw(y) - logw(x[i])):
                x[i] = y
        for i in range(it % 2, len(ladder) - 1, 2):
            if np.log(rng.random()) < (ladder[i] - ladder[i + 1]) * (logw(x[i + 1]) - logw(x[i])):
                x[i], x[i + 1] = x[i + 1], x[i]
        right.append(x[-1] > 0)
    return float(np.mean(right[2000:]))
def plain_toy(n=40000):
    x = 4.0; right = []
    for it in range(n):
        y = x + 0.3 * rng.normal()
        if np.log(rng.random()) < logprior(y) - logprior(x) + logw(y) - logw(x):
            x = y
        right.append(x > 0)
    return float(np.mean(right[2000:]))
pt_r, pl_r = pt_toy(), plain_toy()
step("PT1 toy (two narrow wells of weights 0.3 / 0.7 under a broad prior): plain Metropolis at s = 1 stays in its start well, the tempered ladder p_s ~ prior x w^s with neighbour swaps reproduces the exact well weight at s = 1",
     abs(pt_r - exact_right) < 0.05 and abs(pl_r - exact_right) > 0.3, f"exact weight of the right well {exact_right:.3f}; tempering {pt_r:.3f}; plain {pl_r:.3f}")

# PT2 analysis -------------------------------------------------------------------------------------------------------
rows = []
for fp in sorted(HERE.glob("pt_monitor_*.jsonl")):
    rr = [json.loads(l) for l in fp.read_text().splitlines() if l.strip()]
    if not rr:
        continue
    labels = np.array([x["labels"] for x in rr]); top_lab = labels[:, -1]
    n2_top = np.array([x["norm2"][-1] for x in rr])
    G = [(x["theta_top"], x["G_P+d"][len(x["G_P+d"]) // 2]) for x in rr if "G_P+d" in x]
    Gh = float(np.sum([gg * np.cos(t) for t, gg in G]) / np.sum([np.cos(t) for t, _ in G])) if G else np.nan
    # round trips: a label travelling from the bottom (s_0) to the top and back
    trips = 0
    for lab in range(labels.shape[1]):
        pos = np.array([list(row).index(lab) for row in labels]); state, = [None]
        for p_ in pos:
            if p_ == 0 and state == "top":
                trips += 1; state = "bottom"
            elif p_ == 0:
                state = "bottom"
            elif p_ == labels.shape[1] - 1 and state == "bottom":
                state = "top"
    rows.append(f"{fp.stem[11:]}: {len(rr)} rounds, swap acceptance {np.round(rr[-1]['swap_acc'], 2).tolist()}, round trips {trips}, s = 1 labels visited {sorted(set(top_lab.tolist()))}, "
                f"|phi|^2 at s = 1: first 10 {n2_top[:10].mean():.1f}, last 20 {n2_top[-20:].mean():.1f} (range {n2_top.min():.1f}-{n2_top.max():.1f}), G(beta/2) at s = 1 (reweighted) {Gh:.3f}")
step("PT2 (statement) tempering at g = 0.015, beta = 6: swap rates, round trips, s = 1 history", True, " | ".join(rows) if rows else "tempering runs not started")

out = {"steps": LOG, "passed": sum(x["pass"] for x in LOG), "total": len(LOG)}
(HERE / "tempering_results.json").write_text(json.dumps(out, indent=2, default=float))
print(f"{out['passed']}/{out['total']} steps passed")
sys.exit(0 if out["passed"] == out["total"] else 1)
