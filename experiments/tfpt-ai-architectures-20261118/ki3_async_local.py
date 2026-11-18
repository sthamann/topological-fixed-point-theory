"""KI3 - local self-organisation under delayed updates and loss of internal state.

Uses the trained lattice models of E2 (and the generic-control models if present), no retraining.
  sync            every cell updates in every round, T = 4n rounds
  delayed, same rounds     20 % of the cells skip each round (random), T rounds
  delayed, same updates    20 % skip, rounds extended to T / 0.8 so the number of cell updates is equal
  state loss      after T/2 rounds half of the cells lose their internal state (reset to the input)
Also: the arithmetic of communication (neighbour reads per puzzle) for local cells against
global attention, and whether distant information could have arrived (rounds against depth).
"""
import json, os, collections
import numpy as np, torch
import worlds as W, models as M, e2_local_vs_global as E2

torch.set_num_threads(4)
MODELS = ["plain", "tfpt", "rand", "gen_xor", "gen_z4"]
TESTS = [(6, "random"), (12, "random"), (12, "chain")]

@torch.no_grad()
def run(net, d, q, **kw):
    x, known = M.tens(d, ["x", "known"])
    ok = []
    for s in range(0, len(x), 100):
        torch.manual_seed(s)
        o = net(x[s:s + 100], known[s:s + 100], q, **kw)
        ok.append(M.solved(o[-1], x[s:s + 100], known[s:s + 100]))
    return torch.cat(ok).float().mean().item()

res = collections.defaultdict(lambda: collections.defaultdict(dict))
for g in E2.GROUPS:
    for m in MODELS:
        paths = [E2.ckpt(m, g, s) for s in E2.SEEDS if os.path.exists(E2.ckpt(m, g, s))]
        if not paths:
            continue
        for n, kind in TESTS:
            d = E2.lat("test", E2.NTEST, n, g, 20 + n, kind)
            q = torch.as_tensor(W.lattice_quads(n))
            Tn = 4 * n
            rows = collections.defaultdict(list)
            for p in paths:
                net = M.RelNet(m, seed=0); net.load_state_dict(torch.load(p)); net.eval()
                rows["sync"].append(run(net, d, q, steps=Tn))
                rows["20% delayed, same rounds"].append(run(net, d, q, steps=Tn, drop=0.2))
                rows["20% delayed, same number of updates"].append(run(net, d, q, steps=int(round(Tn / 0.8)), drop=0.2))
                rows["half of the state lost once"].append(run(net, d, q, steps=Tn, wipe_at=Tn // 2))
            res[g][m][f"n={n} {kind}"] = {k: [round(float(np.mean(v)), 3), round(float(np.std(v)), 3)] for k, v in rows.items()}
            res[g][m]["seeds"] = len(paths)
        print(g, m, json.dumps(res[g][m]), flush=True)
comm = {}
for n in (6, 12, 16):
    Q = (n - 1) ** 2
    comm[f"n={n}"] = dict(local_reads_per_round=Q * 4 * 3, local_reads_total_T_4n=Q * 4 * 3 * 4 * n,
                          global_attention_reads_4_layers=4 * (n * n) ** 2,
                          minimum_rounds_for_information_to_cross=n - 1)
json.dump(dict(results=res, communication=comm), open("out/ki3_async_local.json", "w"), indent=1)
