"""E2 - "emergent intelligence without a global attention mechanism" (chat 3.1).

A lattice of four-mark cells: every 2x2 window of an n x n grid must sum to zero.
Local model = weight-shared cellular updates on these windows only (no attention, no
positions).  Global model = transformer.  Train on 6x6, test on 6x6, 12x12, 16x16,
with random givens and with 'chain' puzzles where inference must cross the whole grid.
"""
import json, sys, os, collections
import numpy as np, torch
import worlds as W, models as M

SEEDS, STEPS, NTRAIN, NTEST = [0, 1, 2], 2500, 4000, 400
MODELS, GROUPS = ["transformer", "plain", "tfpt", "rand"], ["xor", "z4"]
TESTS = [(6, "random"), (6, "chain"), (12, "random"), (12, "chain"), (16, "random"), (16, "chain")]
if "--quick" in sys.argv:
    SEEDS, STEPS, NTRAIN, NTEST = [0], 150, 500, 100


def lat(name, m, n, group, seed, kind):
    return W.cached(f"lat_{group}_{name}_{n}_{kind}_{m}", lambda: W.make_lattice_set(m, n, group, seed, kind))


def ckpt(model, group, seed):
    return f"out/cache/e2_{model}_{group}_{seed}.pt"


def build(model, seed):
    return M.Transformer(pos="sin2d", n_rules=0) if model == "transformer" else M.RelNet(model, seed=seed)


def job(a):
    model, group, seed = a
    torch.set_num_threads(2); torch.manual_seed(seed)
    tr = lat("train", NTRAIN, 6, group, 10, "mix")
    net = build(model, seed)
    M.train(net, tr, torch.as_tensor(W.lattice_quads(6)), STEPS, T_steps=12, seed=seed)
    torch.save(net.state_dict(), ckpt(model, group, seed))
    out = dict(model=model, group=group, seed=seed, params=M.n_params(net))
    for n, kind in TESTS:
        d = lat("test", NTEST, n, group, 20 + n, kind)
        q = torch.as_tensor(W.lattice_quads(n))
        out[f"n{n}_{kind}_T12"] = M.evaluate(net, d, q, T_steps=12, bs=100)
        if model != "transformer":
            out[f"n{n}_{kind}_T4n"] = M.evaluate(net, d, q, T_steps=4 * n, bs=100)
        out[f"n{n}_{kind}_depth"] = float(d["depth"].mean())
    return out


if __name__ == "__main__":
    for g in GROUPS:
        lat("train", NTRAIN, 6, g, 10, "mix")
        for n, kind in TESTS:
            lat("test", NTEST, n, g, 20 + n, kind)
    res = M.run_parallel(job, [(m, g, s) for g in GROUPS for m in MODELS for s in SEEDS], nproc=8)
    table = collections.defaultdict(dict)
    for g in GROUPS:
        for m in MODELS:
            rs = [r for r in res if (r["model"], r["group"]) == (m, g)]
            keys = [k for k in rs[0] if k.startswith("n") and not k.endswith("depth")]
            table[g][m] = {k: M.agg([r[k] for r in rs]) for k in keys}
    depth = {k: res[0][k] for k in res[0] if k.endswith("depth")}
    json.dump(dict(table=table, mean_depth=depth, raw=res), open("out/e2_local_vs_global.json", "w"), indent=1)
    print("mean propagation depth", depth)
    for g in GROUPS:
        print(g)
        for m, v in table[g].items():
            print(f"  {m:12s} " + " ".join(f"{k[1:]}={x[0]:.2f}" for k, x in v.items()))
