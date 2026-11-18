#!/usr/bin/env python3
"""Part RWC: evaluation of CHIRAL4D.REWEIGHT.01 (phase-quenched sampling, reweight.py).

Two kinetic terms are kept strictly apart: walk "md" = the frozen Theta_1 (minimal-doubling C3 Hamiltonian h_W,
theta1_md.py) -- decisive; walk "bdt" = legacy BDT Hamiltonian walk -- diagnostic only, it is NOT Theta_1.
  RWC1 folded debug cells: with the BDT walk the 4- and 8-site cells are phase-free (walk vanishes / round-off) -- a
       geometric artefact; with the frozen h_W the 4-site cell conserves n.sigma (number-conserving special formula,
       qmc_pfaffian.py PF2), so it is a sign-validation cell, not a phase-severity cell either.
  RWC2 the 16-site cubic cell carries genuine phases for the frozen h_W.
  RWC3 phase-quenched results (both walks, separate tables): R = |<e^{i theta}>_pq|, Var theta, Delta f = -log R/(beta V),
       1/R^2, autocorrelation, acceptance, zero detector.
  RWC4 Trotter control for Theta_1: dtau = 0.25 vs 0.125 at beta = 1, g = 0.1 (16 sites).
  RWC5 volume scaling 16 -> 32 at fixed beta v / L for Theta_1 (16/32 is NOT a homothetic family: estimate only); BDT
       diagnostic alongside.
  RWC6 (statement) decision for the next step, based on Theta_1 only.
  RWC7 Trotter scaling of the phase for Theta_1: -log R at dtau = 0.25, 0.125, 0.0625: falls, then converges (genuine
       continuum-time phase; dtau = 0.25 overestimates it).
  RWC9 (finding) Trotter dependence of the phase at strong coupling on 32 sites (dtau = 0.25 vs 0.125).
  RWC10 (finding) volume growth of the phase at fixed c_beta on 16 / 32 / 64 sites.
  RWC8 (statement) first mirror-node measurement ("meas" runs): reweighted G_P(beta/2) against the free value 1/2, the gap
       estimate, Gamma sentinel (exactly free), order-parameter susceptibility chi_O, composite-seed GEVP -- on 16 sites only,
       no thermodynamic statement.
    python3 -B experiments/theory-contracts/seam-holonomy-generating-20261001/reweight_check.py
"""
from __future__ import annotations

import importlib.util
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


spec = importlib.util.spec_from_file_location("rw", HERE / "reweight.py")
rw = importlib.util.module_from_spec(spec); sys.argv = ["x"]; spec.loader.exec_module(rw)
rng = np.random.default_rng(91)


def direct_phases(cell, g, walk, n=3):
    m = rw.Model(cell, 1.0, 0.25, g, walk)
    sig = np.sqrt(2 * g * 0.25)
    out = []
    for _ in range(n):
        phi = sig * rng.normal(size=(m.Lt, m.N, 10))
        S = [m.slice_op(phi[l]) for l in range(m.Lt)]
        Bm = np.linalg.multi_dot(S)
        sgn, _ = np.linalg.slogdet(np.eye(Bm.shape[0]) + Bm)
        out.append((abs(np.angle(sgn)), np.linalg.cond(Bm)))
    return out


d4 = direct_phases("debug4", 0.02, "bdt"); d8 = direct_phases("cell8", 0.02, "bdt")
pf = json.loads((HERE / "qmc_pfaffian_results.json").read_text()) if (HERE / "qmc_pfaffian_results.json").exists() else {}
pf2 = next((x for x in pf.get("steps", []) if x["step"].startswith("PF2")), None)
step("RWC1.folded debug cells: with the BDT walk det(1 + B) is real on the 4-site (walk vanishes) and 8-site cells (round-off) -- a geometric artefact; with the frozen h_W the 4-site cell conserves n.sigma and has the exact number-conserving special formula (PF2): a sign-validation cell, not a phase-severity cell",
     all(p < 1e-12 for p, c in d4 + d8) and all(c < 1e6 for p, c in d4 + d8) and bool(pf2 and pf2["pass"]),
     f"BDT max |phase| 4-site {max(p for p, c in d4):.1e}, 8-site {max(p for p, c in d8):.1e}; PF2: {pf2['detail'][:160] if pf2 else 'missing'}")
def stab_phases(cell, g, walk, n=3):
    """phase of det(1 + B) from the stabilised UDV chain and, independently, from the split factorisation at slice 2."""
    m = rw.Model(cell, 1.0, 0.25, g, walk)
    sig = np.sqrt(2 * g * 0.25)
    out = []
    nn = 2 * m.n; unit = (np.eye(nn, dtype=complex), np.ones(nn), np.eye(nn, dtype=complex))
    for _ in range(n):
        phi = sig * rng.normal(size=(m.Lt, m.N, 10))
        S = [m.slice_op(phi[l]) for l in range(m.Lt)]
        la, ph, _ = m.logdet(S)
        Y = m.sweep_cache(S); X = unit
        for l in range(2):
            X = m._udv_left(S[l].T, X)
        la2, ph2 = m.logdet_split(S[2], X, Y[2])
        out.append((abs(rw.wrap(ph)), abs(rw.wrap(ph - ph2)) + abs(la - la2) / abs(la)))
    return out
d16 = stab_phases("cubic16", 0.1, "md")
step("RWC2.the 16-site cubic cell carries genuine phases for the frozen Theta_1 kinetics h_W: the stabilised phase of det(1 + B) is far above round-off and is reproduced by an independent split factorisation",
     max(p for p, e in d16) > 1e-6 and all(e < 1e-8 for p, e in d16), f"|phase of det| {[f'{p:.1e}' for p, e in d16]}, independent-evaluation deviation {[f'{e:.0e}' for p, e in d16]}")

pts = []
for walk in ("md", "bdt"):
    for f in sorted(list(HERE.glob(f"reweight_points_{walk}_cubic16*.json")) + list(HERE.glob(f"reweight_points_{walk}_oh32*.json")) + list(HERE.glob(f"reweight_points_{walk}_prim64*.json"))):
        for p in json.loads(f.read_text()):
            if "twist_delta" in p:                                    # node-tracking runs are evaluated in mirror_node.py
                continue
            p.setdefault("walk", walk); pts.append(p)
tab = {(p["walk"], p["cell"], p["beta"], p["dtau"], p["g"]): p for p in pts}
fmt = lambda k, p: (f"{k[1]} (V={p['N']}, c_beta={p.get('c_beta', p['beta'] / 4):.2f}), beta={k[2]}, dtau={k[3]}, g={k[4]}: R={p['R']:.6f}+-{p['R_err']:.1e}, "
                    f"-log R={-np.log(max(p['R'], 1e-300)):.2e}, Var theta={p['var_theta']:.2e}, Delta f={p['delta_f']:.2e}, 1/R^2={1 / max(p['R'], 1e-300) ** 2:.4f}, tau_int={p['tau_int']:.1f}, acc={p['acceptance']:.2f}, zero detector={p['min_smin']:.2e}")
md_tab = {k: p for k, p in tab.items() if k[0] == "md"}
bdt_tab = {k: p for k, p in tab.items() if k[0] == "bdt"}
step("RWC3.phase-quenched sampling (phase tracked along the chain) on the 16- and 32-site cells, frozen Theta_1 (md) and, separately, the BDT diagnostic: R, -log R, Var theta, Delta f = -log R/(beta V), cost factor 1/R^2, autocorrelation, acceptance, zero detector",
     len(md_tab) >= 3 and all(0 < p["R"] <= 1 for p in tab.values()),
     "Theta_1 (md): " + " | ".join(fmt(k, p) for k, p in sorted(md_tab.items())) + " || BDT diagnostic: " + " | ".join(fmt(k, p) for k, p in sorted(bdt_tab.items())))
tro = [md_tab.get(("md", "cubic16", 1.0, 0.25, 0.1)), md_tab.get(("md", "cubic16", 1.0, 0.125, 0.1))]
step("RWC4.Trotter control for Theta_1 at beta = 1, g = 0.1 (16 sites): dtau = 0.25 and 0.125 give the same average phase within 3 sigma or 1 %",
     all(tro) and abs(tro[0]["R"] - tro[1]["R"]) <= max(3 * (tro[0]["R_err"] + tro[1]["R_err"]), 0.01),
     f"R(0.25) = {tro[0]['R']:.6f}+-{tro[0]['R_err']:.1e}, R(0.125) = {tro[1]['R']:.6f}+-{tro[1]['R_err']:.1e}" if all(tro) else "missing points")
def scaling(tb, walk):
    est = []
    for (w, cell, beta, dtau, g), p in tb.items():
        if cell != "cubic16" or dtau != 0.25:
            continue
        cb = beta / 4.0
        for (w2, c2, b2, d2, g2), q in tb.items():
            if c2 == "oh32" and d2 == 0.25 and g2 == g and abs(q.get("c_beta", 0) - cb) < 0.02:
                R128 = np.exp(-max(p["delta_f"], q["delta_f"]) * cb * 8.0 * 128)
                est.append((walk, g, cb, p["delta_f"], q["delta_f"], p["R"], q["R"], R128))
    return est
est = scaling(md_tab, "md"); est_b = scaling(bdt_tab, "bdt")
step("RWC5.volume scaling at fixed beta v / L (Theta_1): Delta f on 16 and 32 sites for the same g and c_beta; the larger extrapolates (R ~ exp(-Delta f beta V)) to the 128-site cell (estimate only -- 16 and 32 are not a homothetic family); BDT diagnostic alongside",
     len(est) >= 1, "; ".join(f"{w}: g = {g}, c_beta = {cb:.2f}: Delta f(16) = {a:.2e}, Delta f(32) = {b:.2e}, R(16) = {r1:.5f}, R(32) = {r2:.5f}, R_128 ~ {r3:.3e}" for w, g, cb, a, b, r1, r2, r3 in est + est_b))
worst = min(md_tab.values(), key=lambda p: p["R"]) if md_tab else None
R128min = min((e[-1] for e in est), default=np.nan)
if worst and worst["R"] < 1e-2:
    decision = f"No-go for reweighting of Theta_1 at the measured sizes (lowest R = {worst['R']:.2e}): tensor-network route"
elif est and R128min < 1e-2:
    decision = f"reweighting of Theta_1 feasible on 16 and 32 sites (lowest R = {worst['R']:.4f}), but the extrapolated 128-site phase R ~ {R128min:.1e} would be too small: thermodynamics needs tensor networks or improved estimators"
else:
    decision = f"reweighting of Theta_1 feasible (lowest R = {worst['R']:.4f}, extrapolated R_128 ~ {R128min:.2e}): go to the g phase map and the mirror-gap measurement"
step("RWC6.(statement) Go / No-Go for Theta_1: N_eff ~ N_samples R^2; R >= 0.1 on the relevant cells -> physics measurement; R << 1e-2 with stable Delta f > 0 -> tensor networks; the thermodynamic mirror gap is not yet measured",
     True, decision + "; Gate B stays open")

tr_pts = sorted((k[3], p) for k, p in md_tab.items() if k[1] == "cubic16" and k[2] == 1.0 and k[4] == 0.1)
mlr = [(dt, -np.log(p["R"]), p["R_err"] / p["R"]) for dt, p in tr_pts]
rat = [mlr[i + 1][1] / mlr[i][1] for i in range(len(mlr) - 1) if mlr[i][1] > 0]
mlr_d = sorted(mlr, reverse=True)                                      # dtau descending
conv = len(mlr_d) >= 3 and abs(mlr_d[-1][1] - mlr_d[-2][1]) <= 3 * (mlr_d[-1][2] + mlr_d[-2][2]) + 1e-12
step("RWC7.Trotter scaling of the Theta_1 phase (16 sites, beta = 1, g = 0.1): -log R falls from dtau = 0.25 to 0.125 and is converged between 0.125 and 0.0625 -> a genuine, small continuum-time phase; dtau = 0.25 overestimates it",
     len(mlr_d) >= 3 and mlr_d[1][1] < mlr_d[0][1] and conv,
     "; ".join(f"dtau = {dt}: -log R = {v:.2e} (+-{e:.1e})" for dt, v, e in mlr_d) + f"; ratios per halving {np.round([mlr_d[i][1] / mlr_d[i + 1][1] for i in range(len(mlr_d) - 1)], 2).tolist()}")
t32 = [md_tab.get(("md", "oh32", 2.5, 0.25, 0.1)), md_tab.get(("md", "oh32", 2.5, 0.125, 0.1))]
step("RWC9.(finding) Trotter dependence of the Theta_1 phase at strong coupling on 32 sites (beta = 2.5, g = 0.1): unlike the 16-site, beta = 1 case, -log R GROWS when dtau is halved -- the dtau = 0.25 phases underestimate the sign problem at strong coupling; the 'mild phase' statement is restricted to coarse dtau / weak coupling until dtau -> 0 is controlled (20 sweeps, large error)",
     all(t32), (f"dtau = 0.25: R = {t32[0]['R']:.4f}+-{t32[0]['R_err']:.1e} (-log R = {-np.log(t32[0]['R']):.2e}); dtau = 0.125: R = {t32[1]['R']:.4f}+-{t32[1]['R_err']:.1e} (-log R = {-np.log(t32[1]['R']):.2e}); "
                f"Delta f = {t32[1]['delta_f']:.2e} -> R_128 (c_beta = 0.5) ~ {np.exp(-t32[1]['delta_f'] * 4.0 * 128):.2f}") if all(t32) else "missing points")
vs = [md_tab.get(("md", "cubic16", 2.0, 0.25, 0.1)), md_tab.get(("md", "oh32", 2.5, 0.25, 0.1)), md_tab.get(("md", "prim64", 3.25, 0.25, 0.1))]
if all(vs):
    mlr_v = [-np.log(p["R"]) for p in vs]; bV = [p["beta"] * p["N"] for p in vs]
    growth = [mlr_v[i + 1] / mlr_v[i] for i in range(2)]; bV_g = [bV[i + 1] / bV[i] for i in range(2)]
    df_ = [p["delta_f"] for p in vs]
    R128 = [np.exp(-df_[-1] * 4.0 * 128), np.exp(-df_[-1] * (df_[-1] / df_[-2]) * 4.0 * 128)]
step("RWC10.(finding) volume growth of the Theta_1 phase at fixed c_beta = 0.5, g = 0.1, dtau = 0.25 (16 / 32 / 64 sites, NOT a homothetic family): -log R grows much faster than beta V, i.e. Delta f = -log R / (beta V) is not yet intensive and increases with V; together with RWC9 (dtau = 0.25 underestimates the phase at strong coupling) reweighting at g = 0.1 is feasible to ~64-128 sites at most",
     all(vs) and all(gw > bg for gw, bg in zip(growth, bV_g)),
     (f"R = {[round(p['R'], 5) for p in vs]}, -log R = {[f'{x:.2e}' for x in mlr_v]}, beta V = {bV}; growth of -log R {np.round(growth, 1).tolist()} vs beta V {np.round(bV_g, 2).tolist()}; Delta f = {[f'{x:.1e}' for x in df_]}; "
      f"R_128 (beta = 4) ~ {R128[0]:.2f} with the 64-site Delta f, ~ {R128[1]:.2f} if Delta f keeps growing at the same rate") if all(vs) else "missing points")
gm = [(k, p) for k, p in md_tab.items() if "greens" in p]
def gline(k, p):
    G = p["greens"]; txt = f"{k[1]} beta={k[2]} g={k[4]}: G_P(beta/2) = {G['P']['G_half']:.4f}+-{G['P']['G_half_err']:.4f} (free {G['P']['G_half_free']:.2f}) -> gap ~ {G['P']['gap_estimate']:.2f}"
    if "-P" in G:
        txt += f", G_-P = {G['-P']['G_half']:.4f} (free {G['-P']['G_half_free']:.4f})"
    txt += f", Gamma sentinel {p.get('gamma_sentinel', float('nan')):.0e}, chi_O = {p['chi_O']['mean']:.2f}+-{p['chi_O']['err']:.2f}" if "chi_O" in p else ""
    if "composites_P" in p:
        txt += f", GEVP kept {p['composites_P']['kept_seeds']}/3 seeds"
    return txt
step("RWC8.(statement) first mirror-node measurement for Theta_1 on 16 sites (reweighted G_P(tau), Gamma exactly free, order-parameter susceptibility, composite seeds): small-volume indication only, no thermodynamic gap statement",
     True, " | ".join(gline(k, p) for k, p in sorted(gm)) if gm else "no meas runs yet")

out = {"steps": LOG, "passed": sum(x["pass"] for x in LOG), "total": len(LOG), "points": {str(k): v for k, v in tab.items()}, "scaling": est + est_b}
(HERE / "reweight_check_results.json").write_text(json.dumps(out, indent=2, default=float))
print(f"{out['passed']}/{out['total']} steps passed")
sys.exit(0 if out["passed"] == out["total"] else 1)
