"""KI1 control required by the experiment package: the same relational net with an ORDINARY
exact consistency message (factor message of 'four symbols sum to zero'), given the same
knowledge as the TFPT variant.  gen_xor knows the xor rule, gen_z4 knows the Z4 rule.
Both are trained on both worlds, so each is once well informed and once misinformed."""
import json, collections, sys
import torch
import worlds as W, models as M, e1_three_way as E1, e2_local_vs_global as E2

def job_grid(a):
    return E1.job(a)

def job_lat(a):
    return E2.job(a)

if __name__ == "__main__":
    mods = ["gen_xor", "gen_z4"]
    jobs = [(m, g, s, sd) for g in E1.GROUPS for s in E1.SIZES for m in mods for sd in E1.SEEDS]
    res = M.run_parallel(job_grid, jobs, nproc=10)
    t1 = collections.defaultdict(dict)
    for g in E1.GROUPS:
        for s in E1.SIZES:
            for m in mods:
                rs = [r for r in res if (r["model"], r["group"], r["size"]) == (m, g, s)]
                t1[f"{g} N={s}"][m] = {k: M.agg([r[k] for r in rs]) for k in ("id", "ood_combo", "ood_hard")}
    res2 = M.run_parallel(job_lat, [(m, g, s) for g in E2.GROUPS for m in mods for s in E2.SEEDS], nproc=8)
    t2 = collections.defaultdict(dict)
    for g in E2.GROUPS:
        for m in mods:
            rs = [r for r in res2 if (r["model"], r["group"]) == (m, g)]
            t2[g][m] = {k: M.agg([r[k] for r in rs]) for k in rs[0] if k.startswith("n") and not k.endswith("depth")}
    json.dump(dict(grid=t1, lattice=t2), open("out/ki1_generic_control.json", "w"), indent=1)
    print(json.dumps(dict(grid=t1, lattice=t2))[:6000])
