#!/usr/bin/env python3
"""Part GL: CHIRAL4D.GROUNDLINE.RESIDUE.01 -- does a low spectral line belong to the ground sector?

A line at omega ~ 0 in the thermal Green function can be an excited -> excited transition: its residue then needs a
thermally populated initial state and vanishes on cooling, Z_beta ~ A e^{-beta delta}.  A line connected to the ground
sector keeps a finite residue, Z_beta -> Z_0 > 0.  For every identified line (same cell, g, twist and momentum sector)
the series (beta, e_beta, Z_beta, R) is fitted by
    M_G: Z = Z_0 + A e^{-beta delta}  (Z_0 >= 0)      against      M_T: Z = A e^{-beta delta};
status GROUND_CONNECTED_<n>BETA_CANDIDATE (n = number of temperatures) only if Z_0 - 2 sigma > 0 (M_G) and the residue at
the coldest point exceeds 0.1 by 2 sigma -- a candidate, not a proof: a finite set of temperatures cannot exclude a slow
decay (GL4); promotion to MIRROR_GROUND_SECTOR_LIGHT needs the candidate + a low start-sector projected Green gap
(SPG_LOW_GAP) + a low direct sector-trace gap (SECTOR_LOW_GAP);
THERMAL if the residue at the coldest point is compatible with 0 and the decay rate delta of M_T is significant; otherwise
UNDECIDED (a nearly flat M_T fit with tiny delta is not 'thermal').  Fewer than three temperatures -> INSUFFICIENT.  Runs with R < 0.5 are excluded.
  GL1 unit test: a thermal line Z = 0.4 e^{-1.5 beta} sampled at beta = 2, 4, 6 is classified THERMAL; a ground line
      Z = 0.55 + 0.3 e^{-beta} is classified GROUND_CONNECTED_3BETA_CANDIDATE.
  GL4 counter-case (why the status is only a candidate): Z = 0.6 e^{-beta} + 0.45 e^{-0.03 beta} has Z_inf = 0 but at
      beta = 2, 4, 6 an apparent plateau Z_0 ~ 0.37 and is classified as a 3-beta candidate.
  GL2 (statement) the measured lines; for UNDECIDED the diagnostic reason (residue resolution at the coldest point, phase
      R, beta lever) is reported -- it does not change the classification.
  GL3 (statement) one table per coupling for the predeclared branch (Q+ sector at P + lambda (1,1,1)): beta, e_beta, Z_beta, R
      and the matrix windows L_Omega / U_Omega for Omega = 0.1, 0.25, 0.5, 1 (from spectral_bounds_results.json, SB5).
    python3 -B experiments/theory-contracts/seam-holonomy-generating-20261001/groundline_residue.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
LOG: list[dict] = []
_zb = (HERE / "zero_bracket.py").read_text()
exec(compile(_zb.split("def levels_curve")[0].split("LOG: list[dict] = []")[1] + "def levels_curve" + _zb.split("def levels_curve")[1].split("def brackets")[0], "zb", "exec"))


def step(name: str, ok: bool, detail: str = "") -> None:
    LOG.append({"step": name, "pass": bool(ok), "detail": detail})
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))
    sys.stdout.flush()


def classify(betas, Z, sZ):
    betas, Z, sZ = map(np.asarray, (betas, Z, sZ))
    if len(np.unique(betas)) < 3:
        return "INSUFFICIENT", {}
    w = 1 / np.maximum(sZ, 1e-3)
    rG = least_squares(lambda p: (abs(p[0]) + abs(p[1]) * np.exp(-betas * abs(p[2])) - Z) * w, [Z.min(), Z.max() - Z.min() + 1e-3, 0.5])
    rT = least_squares(lambda p: (abs(p[0]) * np.exp(-betas * abs(p[1])) - Z) * w, [Z.max() * np.exp(betas.min() * 0.5), 0.5])
    aG, aT = 2 * rG.cost + 6, 2 * rT.cost + 4
    covG = np.linalg.pinv(rG.jac.T @ rG.jac); covT = np.linalg.pinv(rT.jac.T @ rT.jac)
    Z0, sZ0 = abs(rG.x[0]), float(np.sqrt(abs(covG[0, 0])))
    dT, sdT = abs(rT.x[1]), float(np.sqrt(abs(covT[1, 1])))
    i_max = int(np.argmax(betas))
    info = {"Z0": Z0, "sigma_Z0": sZ0, "delta_T": dT, "sigma_delta_T": sdT, "AIC_G": aG, "AIC_T": aT}
    # physical rule: a thermal line needs a significant decay rate and a vanishing residue at the coldest point; a
    # ground-connected line keeps a finite residue (a nearly flat M_T fit with tiny delta is equivalent to Z_0 > 0)
    if Z0 - 2 * sZ0 > 0 and Z[i_max] - 2 * sZ[i_max] > 0.1:
        return f"GROUND_CONNECTED_{len(np.unique(betas))}BETA_CANDIDATE", info
    if Z[i_max] - 2 * sZ[i_max] <= 0.05 and dT - 2 * sdT > 0:
        return "THERMAL", info
    return "UNDECIDED", info


b_ = np.array([2.0, 4.0, 6.0])
rng = np.random.default_rng(3)
Zt = 0.4 * np.exp(-1.5 * b_) + 0.002 * rng.normal(size=3); Zg = 0.55 + 0.3 * np.exp(-b_) + 0.01 * rng.normal(size=3)
ct, _ = classify(b_, Zt, np.full(3, 0.003)); cg, ig = classify(b_, Zg, np.full(3, 0.015))
step("GL1 unit test: a thermal line Z = 0.4 e^{-1.5 beta} (beta = 2, 4, 6) is classified THERMAL and a ground line Z = 0.55 + 0.3 e^{-beta} GROUND_CONNECTED_3BETA_CANDIDATE",
     ct == "THERMAL" and cg == "GROUND_CONNECTED_3BETA_CANDIDATE", f"thermal -> {ct}; ground -> {cg} (Z_0 = {ig.get('Z0', np.nan):.3f} +- {ig.get('sigma_Z0', np.nan):.3f})")

Zs = 0.6 * np.exp(-b_) + 0.45 * np.exp(-0.03 * b_)
cs, is_ = classify(b_, Zs, np.full(3, 0.01))
step("GL4 counter-case: Z = 0.6 e^{-beta} + 0.45 e^{-0.03 beta} (Z_inf = 0) looks like a plateau at beta = 2, 4, 6 and is classified as a 3-beta candidate -- the status names the evidence level and is promoted only by the projected Green function and the direct sector trace",
     cs == "GROUND_CONNECTED_3BETA_CANDIDATE", f"Z(2, 4, 6) = {np.round(Zs, 3).tolist()} -> {cs}, apparent Z_0 = {is_.get('Z0', np.nan):.3f} +- {is_.get('sigma_Z0', np.nan):.3f}, true Z_inf = 0")

lines: dict = {}
for f in sorted(HERE.glob("reweight_points_md_*.json")):
    for p in json.loads(f.read_text()):
        if "greens" not in p or p["R"] < 0.5:
            continue
        tw = tuple(np.round(p["twist_delta"], 3)) if "twist_delta" in p else None
        for name, v in p["greens"].items():
            if not name.startswith("P"):
                continue
            lo = lowest_signed(np.array(v["G_tau"]), np.array(v["G_err"]), p["beta"], p["dtau"])
            if lo is None:
                continue
            e, Z, se, sZ = lo
            lines.setdefault((p["cell"], p["g"], tw, name), []).append((p["beta"], e, Z, se, sZ, p["R"], f.name))
rows = []
for key, pts in sorted(lines.items(), key=lambda t: (t[0][0], t[0][1], str(t[0][2]), t[0][3])):
    if len({b for b, *_ in pts}) < 2:
        continue
    pts.sort()
    st, info = classify([q[0] for q in pts], [q[2] for q in pts], [q[4] for q in pts])
    if st == "UNDECIDED":
        cold = max(pts); why = []
        if cold[4] > 0.15:
            why.append(f"residue at beta = {cold[0]:g} unresolved (sigma_Z = {cold[4]:.2f})")
        if min(q[5] for q in pts) < 0.8:
            why.append(f"phase R_min = {min(q[5] for q in pts):.2f}")
        if max(q[0] for q in pts) / min(q[0] for q in pts) < 3:
            why.append("beta lever < 3")
        if not why:
            why.append("Z_0 interval reaches 0 while the coldest residue is not compatible with 0")
        info["undecided_reason"] = "; ".join(why)
    rows.append({"cell": key[0], "g": key[1], "twist": key[2], "k": key[3], "status": st, "info": info,
                 "series": [{"beta": b, "e": e, "Z": Z, "sigma_e": se, "sigma_Z": sZ, "R": R, "file": fn} for b, e, Z, se, sZ, R, fn in pts]})
step("GL2 (statement) residue series of the identified lines across beta (runs with R < 0.5 excluded; several chains at one beta are listed separately)",
     True, " | ".join(f"{r['cell']} g={r['g']} tw={r['twist'][0] if r['twist'] else 0} {r['k']}: {r['status']} " + ", ".join(f"beta {s['beta']:g}: e = {s['e']:+.3f}+-{s['sigma_e']:.3f}, Z = {s['Z']:.2f}+-{s['sigma_Z']:.2f} (R {s['R']:.3f})" for s in r["series"]) + (f" [{r['info']['undecided_reason']}]" if "undecided_reason" in r["info"] else "") for r in rows) if rows else "no line with two or more temperatures")

sbf = HERE / "spectral_bounds_results.json"
mrows = json.loads(sbf.read_text()).get("matrix_rows", []) if sbf.exists() else []
tables = []
for r in rows:
    if r["k"] != "P+d" or r["cell"] != "cubic16" or r["g"] not in (0.018, 0.021, 0.024):
        continue
    lines_ = []
    for srow in r["series"]:
        m = [x for x in mrows if x["file"] == srow["file"] and x["k"] == "P+d"]
        W = m[0]["W"] if m else {}
        lines_.append(f"beta {srow['beta']:g} ({srow['file'].split('cubic16_')[1][:-5]}): e = {srow['e']:+.3f}+-{srow['sigma_e']:.3f}, Z = {srow['Z']:.2f}+-{srow['sigma_Z']:.2f}, R = {srow['R']:.3f}; "
                      + (", ".join(f"[L, U]({Om}) = [{v[0]:.2f}, {v[1]:.2f}]" for Om, v in W.items()) if W else "matrix windows n/a"))
    tables.append({"g": r["g"], "twist": r["twist"], "status": r["status"], "lines": lines_})
step("GL3 (statement) per-coupling table of the predeclared branch (Q+ sector, cubic16): beta, e_beta, Z_beta, R and the matrix windows L_Omega = W_min^{some spin}, U_Omega = W_max^{any spin} for Omega = 0.1 / 0.25 / 0.5 / 1",
     True, " || ".join(f"g = {t['g']} (lambda = {t['twist'][0] if t['twist'] else 0}): {t['status']} -- " + " | ".join(t["lines"]) for t in tables) if tables else "no table yet")

out = {"steps": LOG, "passed": sum(x["pass"] for x in LOG), "total": len(LOG), "rows": rows, "tables": tables}
(HERE / "groundline_residue_results.json").write_text(json.dumps(out, indent=2, default=float))
print(f"{out['passed']}/{out['total']} steps passed")
sys.exit(0 if out["passed"] == out["total"] else 1)
