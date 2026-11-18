"""Three model classes for the comparison described in the chat (section 5):

  Transformer        ordinary network, global attention, one forward pass
  RelNet(plain)      relational net: learnable *local* transitions on four-cell hyperedges,
                     weight-shared and iterated
  RelNet(pair)       same, but hyperedges replaced by their pairwise clique expansion
  RelNet(tfpt)       plain + the fixed TFPT operators per four-mark cell: rank-5 code
                     coordinates V^T(z1 z2 z3 z4) and the parameter-free transition d<P>/dz_i
  RelNet(rand)       identical to tfpt but with a random rank-5 isometry (specificity control)
"""
import math, time
import numpy as np, torch, torch.nn as nn, torch.nn.functional as F
import tfpt_ops as T


def mlp(i, h, o):
    return nn.Sequential(nn.Linear(i, h), nn.GELU(), nn.Linear(h, o))


class RelNet(nn.Module):
    def __init__(self, mode="plain", d=64, seed=0):
        super().__init__()
        self.mode, self.d = mode, d
        self.enc = mlp(5, d, d)
        self.phi = mlp(d, d, d)
        self.msg = mlp(2 * d, d, d)
        self.gru = nn.GRUCell(d, d)
        self.read = nn.Linear(d, 4)
        if mode in ("tfpt", "rand"):
            self.V4 = nn.Parameter(T.torch_ops(mode, seed)["V4"], requires_grad=False)
            self.tf = mlp(4 + 5, d, d)
        if mode in ("gen_xor", "gen_z4"):
            # ordinary, equally informed control: exact factor message of the constraint
            # "the four symbols sum to zero" in the stated group; no TFPT object involved.
            G = torch.zeros(4, 4, 4, 4)
            for a in range(4):
                for b in range(4):
                    for c in range(4):
                        G[a, b, c, (a ^ b ^ c) if mode == "gen_xor" else (-(a + b + c)) % 4] = 1.0
            self.G4 = nn.Parameter(G, requires_grad=False)
            self.tf = mlp(4 + 5, d, d)
        if mode == "pair":
            self.pmsg = mlp(2 * d, d, d)

    def forward(self, x, known, quads, active=None, steps=8, drop=0.0, wipe_at=None):
        B, N = x.shape
        xo = F.one_hot(x, 4).float() * known[..., None]
        kf = known.float()
        h = self.enc(torch.cat([xo, kf[..., None]], -1))
        Q = quads.shape[0]
        act = torch.ones(B, Q, device=x.device) if active is None else active.float()
        idx = quads.reshape(-1)
        deg = torch.zeros(B, N, device=x.device).index_add_(1, idx, act.repeat_interleave(4, 1)).clamp(min=1)
        if self.mode == "pair":
            A = torch.zeros(B, N, N, device=x.device)
            for a in range(4):
                for b in range(4):
                    if a != b:
                        A[:, quads[:, a], quads[:, b]] = torch.maximum(A[:, quads[:, a], quads[:, b]], act)
            pdeg = A.sum(-1, keepdim=True).clamp(min=1)
        outs = []
        for _ in range(steps):
            if self.mode == "pair":
                hi = h[:, :, None, :].expand(B, N, N, self.d)
                m = self.pmsg(torch.cat([hi, hi.transpose(1, 2)], -1))
                agg = (m * A[..., None]).sum(2) / pdeg
            else:
                hq = h[:, quads]                                   # B,Q,4,d
                ph = self.phi(hq)
                m = self.msg(torch.cat([hq, ph.sum(2, keepdim=True) - ph], -1))
                if self.mode in ("tfpt", "rand"):
                    p = torch.where(known[..., None], xo, F.softmax(self.read(h), -1))
                    z = torch.sqrt(p + 1e-9)[:, quads]             # amplitudes, B,Q,4,4
                    mu, f = T.code_messages(self.V4, z)
                    mu = mu / (mu.norm(dim=-1, keepdim=True) + 1e-6)
                    m = m + self.tf(torch.cat([mu, f[:, :, None, :].expand(-1, -1, 4, -1)], -1))
                if self.mode in ("gen_xor", "gen_z4"):
                    p = torch.where(known[..., None], xo, F.softmax(self.read(h), -1))[:, quads]
                    a, b, c, d4 = (p[..., i, :] for i in range(4))
                    G = self.G4
                    mu = torch.stack([torch.einsum("abcd,...b,...c,...d->...a", G, b, c, d4),
                                      torch.einsum("abcd,...a,...c,...d->...b", G, a, c, d4),
                                      torch.einsum("abcd,...a,...b,...d->...c", G, a, b, d4),
                                      torch.einsum("abcd,...a,...b,...c->...d", G, a, b, c)], -2)
                    sat = (mu[..., 0, :] * a).sum(-1, keepdim=True)
                    mu = mu / (mu.norm(dim=-1, keepdim=True) + 1e-6)
                    f = torch.cat([sat, torch.zeros_like(sat).expand(-1, -1, 4)], -1)
                    m = m + self.tf(torch.cat([mu, f[:, :, None, :].expand(-1, -1, 4, -1)], -1))
                m = m * act[:, :, None, None]
                agg = torch.zeros(B, N, self.d, device=x.device).index_add_(1, idx, m.reshape(B, Q * 4, self.d))
                agg = agg / deg[..., None]
            hn = self.gru(agg.reshape(B * N, self.d), h.reshape(B * N, self.d)).reshape(B, N, self.d)
            if drop > 0:      # KI3: a random share of cells misses this update (asynchronous / delayed)
                hn = torch.where(torch.rand(B, N, 1, device=x.device) < drop, h, hn)
            h = hn
            if wipe_at is not None and len(outs) == wipe_at:   # KI3: erase half of the internal state once
                h = torch.where(torch.rand(B, N, 1, device=x.device) < 0.5, self.enc(torch.cat([xo, kf[..., None]], -1)), h)
            outs.append(self.read(h))
        return outs


class Transformer(nn.Module):
    """Global attention baseline. pos='learned' (4x4 worlds, gets the rule bits) or 'sin2d'."""
    def __init__(self, d=96, layers=4, heads=4, n_rules=5, pos="learned", n_cells=16):
        super().__init__()
        self.d, self.pos = d, pos
        self.inp = nn.Linear(5, d)
        self.rule = nn.Linear(n_rules, d) if n_rules else None
        if pos == "learned":
            self.pe = nn.Parameter(torch.randn(n_cells, d) * 0.02)
        layer = nn.TransformerEncoderLayer(d, heads, 2 * d, dropout=0.0, batch_first=True, norm_first=True)
        self.tr = nn.TransformerEncoder(layer, layers, enable_nested_tensor=False)
        self.read = nn.Linear(d, 4)

    def _sin(self, N, dev):
        n = int(round(math.sqrt(N)))
        r, c = torch.meshgrid(torch.arange(n, device=dev), torch.arange(n, device=dev), indexing="ij")
        f = torch.exp(-math.log(100.0) * torch.arange(self.d // 4, device=dev) / (self.d // 4))
        pe = [fn(v.reshape(-1, 1) * f) for v in (r, c) for fn in (torch.sin, torch.cos)]
        return torch.cat(pe, -1)

    def forward(self, x, known, quads=None, active=None, steps=None, rule=None):
        xo = F.one_hot(x, 4).float() * known[..., None]
        h = self.inp(torch.cat([xo, known.float()[..., None]], -1))
        h = h + (self.pe if self.pos == "learned" else self._sin(x.shape[1], x.device))
        if rule is not None and self.rule is not None:
            h = h + self.rule(rule.float())[:, None, :]
        return [self.read(self.tr(h))]


def n_params(m):
    return sum(p.numel() for p in m.parameters() if p.requires_grad)


def tens(d, keys):
    out = []
    for k in keys:
        a = torch.as_tensor(d[k])
        out.append(a.long() if k == "x" else a.bool())
    return out


def solved(logits, x, known):
    pred = torch.where(known, x, logits.argmax(-1))
    return (pred == x).all(1)


def train(model, data, quads, steps, bs=64, lr=2e-3, T_steps=8, use_rule=False, use_active=False, seed=0):
    g = torch.Generator().manual_seed(seed)
    x, known = tens(data, ["x", "known"])
    rule = torch.as_tensor(data["rule"]) if use_rule else None
    active = torch.as_tensor(data["active"]) if use_active else None
    opt = torch.optim.Adam(model.parameters(), lr=lr)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, steps)
    n = len(x)
    model.train()
    for _ in range(steps):
        i = torch.randint(0, n, (min(bs, n),), generator=g)
        kw = dict(steps=T_steps)
        if use_rule: kw["rule"] = rule[i]
        if use_active: kw["active"] = active[i]
        outs = model(x[i], known[i], quads, **kw)
        unk = ~known[i]
        loss = sum(F.cross_entropy(o[unk], x[i][unk]) for o in outs) / len(outs)
        opt.zero_grad(); loss.backward()
        nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        opt.step(); sched.step()
    return loss.item()


@torch.no_grad()
def evaluate(model, data, quads, T_steps=8, use_rule=False, use_active=False, bs=250):
    model.eval()
    x, known = tens(data, ["x", "known"])
    ok = []
    for s in range(0, len(x), bs):
        kw = dict(steps=T_steps)
        if use_rule: kw["rule"] = torch.as_tensor(data["rule"][s:s + bs])
        if use_active: kw["active"] = torch.as_tensor(data["active"][s:s + bs])
        outs = model(x[s:s + bs], known[s:s + bs], quads, **kw)
        ok.append(solved(outs[-1], x[s:s + bs], known[s:s + bs]))
    return torch.cat(ok).float().mean().item()


def run_parallel(fn, jobs, nproc=10):
    import multiprocessing as mp
    ctx = mp.get_context("spawn")
    t0 = time.time()
    with ctx.Pool(nproc) as pool:
        res = []
        for k, r in enumerate(pool.imap_unordered(fn, jobs)):
            res.append(r)
            if (k + 1) % max(1, len(jobs) // 10) == 0:
                print(f"  {k + 1}/{len(jobs)} jobs  {time.time() - t0:.0f}s", flush=True)
    return res


def agg(vals):
    a = np.array(vals, dtype=float)
    return [round(float(a.mean()), 4), round(float(a.std()), 4), len(a)]
