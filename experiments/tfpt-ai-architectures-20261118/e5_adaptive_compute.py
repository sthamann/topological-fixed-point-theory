"""E5 - "adaptive compute time instead of fixed depth" (chat 3.4).

The trained local networks of E2 are run for up to T_max steps and stopped by a rule:
  fixed       always T steps
  stable      stop when the decoded world did not change for 2 steps (heuristic)
  certificate stop when every four-mark cell of the decoded world has zero defect.
              For xor worlds this defect is exactly the TFPT stabiliser check 1 - 4<Q_P>.
A halt is 'false' when the rule stopped on a wrong world.  Second part: the parameter-free
P-contraction dynamics, to see how computation time scales with inference depth.
"""
import json, os, sys, collections
import numpy as np, torch, torch.nn.functional as F
import worlds as W, models as M, tfpt_ops as T
import e2_local_vs_global as E2

torch.set_num_threads(8)
res = {}


def defect(pred, quads, group):
    v = pred[:, quads]
    s = v[..., 0] ^ v[..., 1] ^ v[..., 2] ^ v[..., 3] if group == "xor" else v.sum(-1) % 4
    return (s != 0).sum(1)


@torch.no_grad()
def trace(net, d, quads, group, tmax):
    x, known = M.tens(d, ["x", "known"])
    preds = []
    for s in range(0, len(x), 100):
        outs = net(x[s:s + 100], known[s:s + 100], quads, steps=tmax)
        preds.append(torch.stack([torch.where(known[s:s + 100], x[s:s + 100], o.argmax(-1)) for o in outs]))
    P = torch.cat(preds, 1)                                  # T,B,N
    correct = (P == x[None]).all(-1)                         # T,B
    cert = torch.stack([defect(p, quads, group) == 0 for p in P])
    same = torch.zeros_like(correct)
    same[2:] = (P[2:] == P[1:-1]).all(-1) & (P[1:-1] == P[:-2]).all(-1)
    return correct, cert, same


def halt(correct, flag, tmax):
    """first step where flag is true (else tmax); returns accuracy, mean steps, false-halt rate."""
    Tn, B = flag.shape
    idx = torch.where(flag.any(0), flag.float().argmax(0), torch.full((B,), Tn - 1))
    ok = correct[idx, torch.arange(B)]
    halted = flag.any(0)
    return dict(acc=round(ok.float().mean().item(), 4), steps=round((idx + 1).float().mean().item(), 2),
                false_halt=round((halted & ~ok).float().mean().item(), 4))


for group in ("xor", "z4"):
    for model in ("plain", "tfpt"):
        for n, kind in ((6, "chain"), (12, "random"), (12, "chain"), (16, "chain")):
            rows = collections.defaultdict(list)
            for seed in E2.SEEDS:
                path = E2.ckpt(model, group, seed)
                if not os.path.exists(path):
                    continue
                net = E2.build(model, seed); net.load_state_dict(torch.load(path)); net.eval()
                d = E2.lat("test", E2.NTEST, n, group, 20 + n, kind)
                quads = torch.as_tensor(W.lattice_quads(n))
                tmax = 6 * n
                c, cert, same = trace(net, d, quads, group, tmax)
                rows["fixed T=12"].append(dict(acc=c[11].float().mean().item(), steps=12, false_halt=0))
                rows[f"fixed T={tmax}"].append(dict(acc=c[-1].float().mean().item(), steps=tmax, false_halt=0))
                rows["stable"].append(halt(c, same, tmax))
                rows["certificate"].append(halt(c, cert, tmax))
                dep = torch.as_tensor(d["depth"])
                idx = torch.where(cert.any(0), cert.float().argmax(0) + 1, torch.full_like(dep, tmax))
                solved = cert.any(0)
                if solved.sum() > 5:
                    rows["_steps_per_depth"].append(dict(v=(idx[solved].float() / dep[solved].float()).mean().item()))
            key = f"{group} {model} n={n} {kind}"
            res[key] = {r: {k: round(float(np.mean([a[k] for a in v])), 4) for k in v[0]} for r, v in rows.items()}
            print(key, res[key], flush=True)

# ---- parameter-free dynamics: each unknown cell takes the P-contraction of its windows
ops = T.torch_ops()
@torch.no_grad()
def free_dynamics(d, n, tmax):
    x, known = M.tens(d, ["x", "known"])
    quads = torch.as_tensor(W.lattice_quads(n)); idx = quads.reshape(-1)
    xo = F.one_hot(x, 4).float()
    p = torch.where(known[..., None], xo, torch.full_like(xo, 0.25))
    first = torch.full((len(x),), -1)
    for t in range(1, tmax + 1):
        mu, _ = T.code_messages(ops["V4"], torch.sqrt(p[:, quads]))
        agg = torch.zeros_like(p).index_add_(1, idx, mu.reshape(len(x), -1, 4))
        newp = agg / agg.sum(-1, keepdim=True).clamp(min=1e-12)
        newp = torch.where(agg.sum(-1, keepdim=True) > 1e-12, newp, p)
        p = torch.where(known[..., None], xo, newp)
        pred = p.argmax(-1)
        done = (defect(pred, quads, "xor") == 0) & (p.max(-1).values > 0.99).all(1)
        first[(first < 0) & done] = t
    ok = (p.argmax(-1) == x).all(1)
    return first, ok

for n, kind in ((6, "chain"), (12, "chain"), (16, "chain"), (24, "chain"), (12, "random")):
    d = W.cached(f"e5_free_{n}_{kind}", lambda: W.make_lattice_set(200, n, "xor", 90 + n, kind))
    first, ok = free_dynamics(d, n, 4 * n)
    dep = torch.as_tensor(d["depth"]).float()
    res[f"parameter-free n={n} {kind}"] = dict(
        acc=round(ok.float().mean().item(), 4), mean_depth=round(dep.mean().item(), 2),
        mean_steps=round(first.float().mean().item(), 2),
        steps_equal_depth=round((first == dep).float().mean().item(), 4),
        corr_steps_depth=round(float(np.corrcoef(first.numpy(), dep.numpy())[0, 1]), 4) if dep.std() > 0 else None)
    print(f"parameter-free n={n} {kind}", res[f"parameter-free n={n} {kind}"], flush=True)

json.dump(res, open("out/e5_adaptive_compute.json", "w"), indent=1)
