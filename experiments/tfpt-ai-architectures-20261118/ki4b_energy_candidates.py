"""KI4 (b) - eight candidates through the same iteration, four energies, fixed-point bookkeeping.

Energies on every four-mark cell (same optimiser, same step size, same number of steps;
Adam removes the overall scale):
  learned     a generic small network trained on examples of valid and invalid cells
              (it is told nothing about the rule; 2000 gradient steps)
  constraint  the ordinary explicit constraint energy 1 - Pr[xor = 0]
  I - P       the quartic code projector on amplitudes (same ground space as H_mom)
  H_mom       the protection energy I - 20 M4
For every puzzle eight random starts.  Recorded per energy:
  candidates correct, puzzles with at least one correct candidate,
  lowest-energy candidate correct (energy used as selector),
  invalid fixed points (iteration settled on a state that breaks a rule),
  cycles / not settled, and
  'low energy hides a wrong state': the selected state is wrong AND has lower energy than
  the true solution.
"""
import json
import numpy as np, torch, torch.nn as nn, torch.nn.functional as F
import worlds as W, tfpt_ops as T, e4_energy_reasoning as E4

torch.set_num_threads(4)
torch.manual_seed(0)
ops = T.torch_ops()

# ---- generic learned energy from examples of single cells
tuples = torch.tensor([(a, b, c, d) for a in range(4) for b in range(4) for c in range(4) for d in range(4)])
valid = ((tuples[:, 0] ^ tuples[:, 1] ^ tuples[:, 2] ^ tuples[:, 3]) == 0).float()
net = nn.Sequential(nn.Linear(16, 64), nn.ReLU(), nn.Linear(64, 64), nn.ReLU(), nn.Linear(64, 1))
opt = torch.optim.Adam(net.parameters(), lr=2e-3)
for _ in range(2000):
    i = torch.randint(0, 256, (256,))
    soft = F.softmax(6 * F.one_hot(tuples[i], 4).float() + torch.randn(256, 4, 4), -1)     # noisy versions of the examples
    loss = F.binary_cross_entropy_with_logits(net(soft.reshape(256, 16)).squeeze(-1), 1 - valid[i])
    opt.zero_grad(); loss.backward(); opt.step()
for q in net.parameters():
    q.requires_grad_(False)
acc = ((net(F.one_hot(tuples, 4).float().reshape(256, 16)).squeeze(-1) > 0).float() == 1 - valid).float().mean().item()

ENERGY = {
    "learned": lambda p: torch.sigmoid(net(p.reshape(*p.shape[:-2], 16)).squeeze(-1)),
    "constraint": E4.ENERGY["Qc"],
    "I - P": E4.ENERGY["Pamp"],
    "H_mom": E4.ENERGY["Hmom"],
}

def run(d, quads, name, active=None, steps=400, cands=8):
    x = torch.as_tensor(d["x"]).long(); known = torch.as_tensor(d["known"]).bool()
    quads = torch.as_tensor(quads); B = len(x)
    act = torch.ones(B, len(quads)) if active is None else torch.as_tensor(active).float()
    xo = F.one_hot(x, 4).float()
    en = ENERGY[name]
    e_true = (en(xo[:, quads]) * act).sum(1)
    out = []
    for c in range(cands):
        torch.manual_seed(100 + c)
        th = (0.5 * torch.randn(B, x.shape[1], 4)).requires_grad_(True)
        opt = torch.optim.Adam([th], lr=0.3)
        hist = []
        for t in range(steps):
            p = torch.where(known[..., None], xo, F.softmax(th, -1))
            E = (en(p[:, quads]) * act).sum(1)
            opt.zero_grad(); E.sum().backward(); opt.step()
            if t >= steps - 40 and t % 10 == 9:
                hist.append(torch.where(known, x, th.argmax(-1)).clone())
        pred = hist[-1]
        settled = torch.stack([(h == pred).all(1) for h in hist]).all(0)
        v = pred[:, quads]
        broken = (((v[..., 0] ^ v[..., 1] ^ v[..., 2] ^ v[..., 3]) != 0).float() * act).sum(1) > 0
        with torch.no_grad():
            e_fin = (en(F.one_hot(pred, 4).float()[:, quads]) * act).sum(1)
        out.append(dict(ok=(pred == x).all(1), settled=settled, broken=broken, e=E.detach(), e_hard=e_fin))
    ok = torch.stack([o["ok"] for o in out]); e = torch.stack([o["e"] for o in out])
    sel = e.argmin(0); ar = torch.arange(B)
    settled = torch.stack([o["settled"] for o in out]); broken = torch.stack([o["broken"] for o in out])
    e_hard = torch.stack([o["e_hard"] for o in out])
    r = lambda t: round(t.float().mean().item(), 4)
    return {"candidates correct": r(ok), "puzzles with a correct candidate": r(ok.any(0)),
            "lowest-energy candidate correct": r(ok[sel, ar]),
            "invalid fixed points": r(settled & broken), "not settled": r(~settled),
            "low energy hides a wrong state": r(~ok[sel, ar] & (e_hard[sel, ar] < e_true - 1e-6))}

res = {"learned_energy_accuracy_on_the_256_cell_states": round(acc, 4)}
sets = [("4x4 worlds, 1-5 hidden", W.cached("e4_grid_1", lambda: None), W.QUADS16, True),
        ("4x4 worlds, 6-7 hidden", W.cached("e4_grid_6", lambda: None), W.QUADS16, True),
        ("lattice 6x6, random givens", W.cached("e4_lat_6_random_300", lambda: None), W.lattice_quads(6), False),
        ("lattice 6x6, chain", W.cached("e4_lat_6_chain_300", lambda: None), W.lattice_quads(6), False)]
for name, d, quads, use_act in sets:
    d = {k: v[:150] for k, v in d.items()}
    res[name] = {e: run(d, quads, e, d["active"] if use_act else None) for e in ENERGY}
    print(name, json.dumps(res[name]), flush=True)
json.dump(res, open("out/ki4b_energy_candidates.json", "w"), indent=1)
