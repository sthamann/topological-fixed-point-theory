"""E1 - the three-way comparison of the chat (section 5), plus two controls.

Task: complete a partially hidden 4x4 world so that every active four-cell rule holds.
Rules come from five families (rows, columns, blocks, two diagonal families).
Measured: data efficiency, unseen rule combinations, harder puzzles (more hidden cells).
Groups: 'xor' (native to the quartic code) and 'z4' (same shape, wrong group).
"""
import json, sys, collections
import numpy as np, torch
import worlds as W, models as M

SIZES, SEEDS, STEPS = [30, 300, 3000], [0, 1, 2], 3000
BIG = 30000          # extra run: does the plain transformer get there with 10x the data?
MODELS = ["transformer", "pair", "plain", "tfpt", "rand"]
GROUPS = ["xor", "z4"]
if "--quick" in sys.argv:
    SIZES, SEEDS, STEPS, BIG = [300], [0], 200, 600


def data(group):
    mk = lambda name, n, combos, ks, seed: W.cached(
        f"e1_{group}_{name}", lambda: W.make_grid_set(n, combos, ks, group, seed))
    return dict(train=mk(f"train{BIG}", BIG, W.TRAIN_COMBOS, [1, 2, 3, 4, 5], 1),
                id=mk("id", 1000, W.TRAIN_COMBOS, [1, 2, 3, 4, 5], 2),
                ood_combo=mk("oodc", 1000, W.OOD_COMBOS, [1, 2, 3, 4, 5], 3),
                ood_hard=mk("oodh", 1000, W.TRAIN_COMBOS + W.OOD_COMBOS, [6, 7], 4))


def job(a):
    model, group, size, seed = a
    torch.set_num_threads(2); torch.manual_seed(seed)
    d = data(group)
    tr = {k: v[:size] for k, v in d["train"].items()}
    quads = torch.as_tensor(W.QUADS16)
    if model == "transformer":
        net, kw = M.Transformer(), dict(use_rule=True)
        M.train(net, tr, quads, (5 if size == BIG else 2) * STEPS, lr=5e-4, seed=seed, **kw)
    else:
        net, kw = M.RelNet(model, seed=seed), dict(use_active=True)
        M.train(net, tr, quads, STEPS, seed=seed, **kw)
    out = dict(model=model, group=group, size=size, seed=seed, params=M.n_params(net))
    for split in ("id", "ood_combo", "ood_hard"):
        out[split] = M.evaluate(net, d[split], quads, **kw)
    if model != "transformer":   # does more thinking time help on the hard set?
        out["ood_hard_T16"] = M.evaluate(net, d["ood_hard"], quads, T_steps=16, **kw)
    return out


if __name__ == "__main__":
    for g in GROUPS:
        data(g)
    jobs = [(m, g, s, sd) for g in GROUPS for s in SIZES for m in MODELS for sd in SEEDS]
    jobs += [("transformer", g, BIG, 0) for g in GROUPS]
    res = M.run_parallel(job, jobs, nproc=12)
    big = {r["group"]: {k: round(r[k], 4) for k in ("id", "ood_combo", "ood_hard")} for r in res if r["size"] == BIG}
    res = [r for r in res if r["size"] != BIG]
    table = collections.defaultdict(dict)
    for g in GROUPS:
        for s in SIZES:
            for m in MODELS:
                rs = [r for r in res if (r["model"], r["group"], r["size"]) == (m, g, s)]
                table[f"{g} N={s}"][m] = {k: M.agg([r[k] for r in rs])
                                          for k in ("id", "ood_combo", "ood_hard", "ood_hard_T16") if k in rs[0]}
    params = {m: next(r["params"] for r in res if r["model"] == m) for m in MODELS}
    json.dump(dict(params=params, table=table, transformer_10x_data=big, raw=res), open("out/e1_three_way.json", "w"), indent=1)
    print("params", params, "\ntransformer with 10x data:", big)
    for k, row in table.items():
        print(k)
        for m, v in row.items():
            print(f"   {m:12s} " + "  ".join(f"{kk}={vv[0]:.3f}±{vv[1]:.3f}" for kk, vv in v.items()))
