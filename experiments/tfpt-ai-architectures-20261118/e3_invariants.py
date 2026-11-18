"""E3 - "learning through invariants instead of surface patterns" (chat 3.2).

Input: four registers z_i in R^4 (unit vectors).  The five quartic code coordinates
f = V^T (z1 z2 z3 z4) are invariant under register permutations and the Pauli group.
Three hidden targets:  inv (a function of f only), perm (permutation invariant only),
raw (no symmetry).  Students differ only in the features they see.  Metric: test R^2.
The point is to measure both sides: the gain when the symmetry fits and the price when not.
"""
import json, sys, collections
import numpy as np, torch, torch.nn as nn
import tfpt_ops as T, models as M

SIZES, SEEDS, STEPS = [100, 1000, 10000], [0, 1, 2], 2000
if "--quick" in sys.argv:
    SIZES, SEEDS, STEPS = [1000], [0], 300
TARGETS = ["inv", "perm", "raw"]
STUDENTS = ["raw", "deepsets", "tfpt5", "raw+tfpt5", "rand5", "raw+rand5"]


def feats(z, kind):
    V4 = T.torch_ops(kind)["V4"]
    f = T.code_feat(V4, z)
    return f * (8.0 if kind == "tfpt" else 16.0)


def sample(n, g):
    z = torch.randn(n, 4, 4, generator=g)
    return z / z.norm(dim=-1, keepdim=True)


def teacher(target, tseed):
    g = torch.Generator().manual_seed(100 + tseed)
    def rnd(i, o):
        return torch.randn(i, o, generator=g) / np.sqrt(i)
    if target == "inv":
        A, B, c = rnd(5, 32) * 2, rnd(32, 32) * 2, rnd(32, 1)
        return lambda z: (torch.tanh(torch.tanh(feats(z, "tfpt") @ A) @ B) @ c).squeeze(-1)
    if target == "perm":
        A, B, c = rnd(4, 32) * 3, rnd(32, 32) * 2, rnd(32, 1)
        return lambda z: (torch.tanh(torch.tanh(z @ A).sum(1) @ B) @ c).squeeze(-1)
    A, B, c = rnd(16, 32) * 3, rnd(32, 32) * 2, rnd(32, 1)
    return lambda z: (torch.tanh(torch.tanh(z.reshape(len(z), 16) @ A) @ B) @ c).squeeze(-1)


class Student(nn.Module):
    def __init__(self, kind):
        super().__init__()
        self.kind = kind
        if kind == "deepsets":
            self.phi, self.rho = M.mlp(4, 96, 96), M.mlp(96, 96, 1)
        else:
            d = {"raw": 16, "tfpt5": 5, "rand5": 5, "raw+tfpt5": 21, "raw+rand5": 21}[kind]
            self.net = nn.Sequential(nn.Linear(d, 128), nn.GELU(), nn.Linear(128, 128), nn.GELU(), nn.Linear(128, 1))

    def forward(self, z):
        if self.kind == "deepsets":
            return self.rho(self.phi(z).sum(1)).squeeze(-1)
        parts = []
        if "raw" in self.kind: parts.append(z.reshape(len(z), 16))
        if "tfpt5" in self.kind: parts.append(feats(z, "tfpt"))
        if "rand5" in self.kind: parts.append(feats(z, "rand"))
        return self.net(torch.cat(parts, -1)).squeeze(-1)


def job(a):
    target, student, size, seed = a
    torch.set_num_threads(1); torch.manual_seed(seed)
    g = torch.Generator().manual_seed(seed)
    t = teacher(target, 0)
    ztr, zte = sample(size, g), sample(5000, g)
    big = t(sample(20000, g)); mu, sd = big.mean(), big.std()
    ytr, yte = (t(ztr) - mu) / sd, (t(zte) - mu) / sd
    net = Student(student)
    opt = torch.optim.AdamW(net.parameters(), lr=2e-3, weight_decay=1e-4)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, STEPS)
    for _ in range(STEPS):
        i = torch.randint(0, size, (min(128, size),), generator=g)
        loss = ((net(ztr[i]) - ytr[i]) ** 2).mean()
        opt.zero_grad(); loss.backward(); opt.step(); sched.step()
    with torch.no_grad():
        r2 = 1 - ((net(zte) - yte) ** 2).mean().item() / yte.var().item()
    return dict(target=target, student=student, size=size, seed=seed, r2=r2)


if __name__ == "__main__":
    jobs = [(t, s, n, sd) for t in TARGETS for s in STUDENTS for n in SIZES for sd in SEEDS]
    res = M.run_parallel(job, jobs, nproc=10)
    table = collections.defaultdict(dict)
    for t in TARGETS:
        for s in STUDENTS:
            table[t][s] = {n: M.agg([r["r2"] for r in res if (r["target"], r["student"], r["size"]) == (t, s, n)])
                           for n in SIZES}
    json.dump(dict(table=table, raw=res), open("out/e3_invariants.json", "w"), indent=1)
    for t in TARGETS:
        print("target:", t)
        for s in STUDENTS:
            print(f"   {s:10s} " + "  ".join(f"N={n}: R2={v[0]:+.3f}" for n, v in table[t][s].items()))
