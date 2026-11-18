"""KI4c and KI7 - the two tests the chat asked for after its own toy run.

KI4c  "Use the real TFPT matrices, against random and learnable subspaces, without making the
      data fit."  256-dimensional inputs, 40 training examples, ridge regression (the chat's
      protocol).  Feature maps: the real rank-5 code projector V, a random rank-5 projector,
      a rank-5 projector learned from the inputs (PCA), a rank-5 projector learned with labels
      (partial least squares), and all 256 dimensions.  Data sets:
        aligned     signal placed inside the TFPT five-space (function control, must win)
        generic     signal in a random five-space
        symmetric   signal in a random five-space of Sym^4 (shares the permutation symmetry only)
        digits      handwritten digits upsampled to 16x16 = 256 pixels, read as four registers
                    of 4 values; nothing is aligned to TFPT
KI7   "Discrete reflection dynamics for long computation chains."  A recurrent step made of
      reflections U = I - 2 d d^T and a shift, against an ordinary tanh recurrent step and a
      GRU: how well is an input from T steps ago still recoverable, and how do gradients scale.
"""
import json
import numpy as np, torch, torch.nn as nn
from sklearn.linear_model import Ridge
from sklearn.cross_decomposition import PLSRegression
from sklearn.datasets import load_digits
import tfpt_ops as T

rng = np.random.default_rng(3)
o = T.build()
V = o["V"]
w, Us = np.linalg.eigh(o["S"]); Bsym = Us[:, w > 0.5]
res = {}

def run(make, trials=40, n=40):
    out = {k: [] for k in ("tfpt_rank5", "random_rank5", "pca_rank5", "pls_rank5", "full_256")}
    for t in range(trials):
        X, y = make(np.random.default_rng(1000 + t), n + 1000)
        Xtr, ytr, Xte, yte = X[:n], y[:n], X[n:], y[n:]
        mu = Xtr.mean(0)
        Rr = np.linalg.qr(np.random.default_rng(5000 + t).standard_normal((256, 5)))[0]
        pca = np.linalg.svd(Xtr - mu, full_matrices=False)[2][:5].T
        for name, F in (("tfpt_rank5", V), ("random_rank5", Rr), ("pca_rank5", pca), ("full_256", np.eye(256))):
            m = Ridge(alpha=10).fit((Xtr - mu) @ F, ytr)
            out[name].append(float(((m.predict((Xte - mu) @ F) - yte) ** 2).mean() / yte.var()))
        pls = PLSRegression(n_components=5).fit(Xtr, ytr)
        out["pls_rank5"].append(float(((pls.predict(Xte).ravel() - yte) ** 2).mean() / yte.var()))
    return {k: [round(float(np.mean(v)), 3), round(float(np.std(v)), 3)] for k, v in out.items()}

def synth(basis_fn):
    def make(g, m):
        X = g.standard_normal((m, 256))
        Bm = basis_fn(g)
        wv = g.standard_normal(5); wv /= np.linalg.norm(wv)
        return X, X @ Bm @ wv + 0.3 * g.standard_normal(m)
    return make

dg = load_digits()
img = dg.images / 16.0
up = np.kron(img, np.ones((1, 2, 2))).reshape(len(img), 256)       # 16x16, index = (r1 r2 c1 c2) registers
def digits(g, m):
    idx = g.permutation(len(up))[:m]
    X = up[idx] - up.mean(0)
    target = g.integers(10)
    return X + 0.0, (dg.target[idx] == target).astype(float)

res["KI4c"] = {
    "aligned to TFPT five-space (function control)": run(synth(lambda g: V)),
    "signal in a random five-space": run(synth(lambda g: np.linalg.qr(g.standard_normal((256, 5)))[0])),
    "signal in a random five-space of Sym^4": run(synth(lambda g: Bsym @ np.linalg.qr(g.standard_normal((35, 5)))[0])),
    "handwritten digits, one-vs-rest": run(digits),
    "metric": "test MSE / variance of target (1.0 = no better than predicting the mean), mean and std over 40 trials",
}
print(json.dumps(res["KI4c"], indent=1), flush=True)

# ------------------------------------------------------------------ KI7
torch.manual_seed(0)
d = 64
class Reflect(nn.Module):
    """h <- shift( prod_k (I - 2 d_k d_k^T) h ) + B x ; norm-preserving on h."""
    def __init__(self, k=4):
        super().__init__()
        self.D = nn.Parameter(torch.randn(k, d)); self.B = nn.Linear(8, d)
    def forward(self, h, x):
        for v in self.D:
            v = v / v.norm()
            h = h - 2 * (h @ v)[:, None] * v
        return torch.roll(h, 1, dims=1) + self.B(x)
class Tanh(nn.Module):
    def __init__(self):
        super().__init__()
        self.Wm = nn.Linear(d, d); self.B = nn.Linear(8, d)
    def forward(self, h, x):
        return torch.tanh(self.Wm(h) + self.B(x))
class Gru(nn.Module):
    def __init__(self):
        super().__init__()
        self.c = nn.GRUCell(8, d)
    def forward(self, h, x):
        return self.c(x, h)

def recall_task(cell, Tlen, steps=600, seed=0):
    """see one of 8 symbols at step 0, then Tlen blank steps, then name it."""
    torch.manual_seed(seed)
    read = nn.Linear(d, 8)
    opt = torch.optim.Adam(list(cell.parameters()) + list(read.parameters()), lr=3e-3)
    for _ in range(steps):
        y = torch.randint(0, 8, (64,))
        h = torch.zeros(64, d)
        h = cell(h, nn.functional.one_hot(y, 8).float())
        for _ in range(Tlen):
            h = cell(h, torch.zeros(64, 8))
        loss = nn.functional.cross_entropy(read(h), y)
        opt.zero_grad(); loss.backward()
        nn.utils.clip_grad_norm_(cell.parameters(), 1.0); opt.step()
    with torch.no_grad():
        y = torch.randint(0, 8, (2000,))
        h = cell(torch.zeros(2000, d), nn.functional.one_hot(y, 8).float())
        for _ in range(Tlen):
            h = cell(h, torch.zeros(2000, 8))
        return float((read(h).argmax(1) == y).float().mean())

def grad_norm(cell, Tlen):
    """|d h_T / d h_0| for a random direction, untrained cell."""
    torch.manual_seed(1)
    h0 = torch.randn(1, d, requires_grad=True)
    h = h0
    for _ in range(Tlen):
        h = cell(h, torch.zeros(1, 8))
    g = torch.autograd.grad(h.sum(), h0)[0]
    return float(g.norm())

torch.set_num_threads(4)
ki7 = {}
for name, mk in (("reflection + shift", Reflect), ("tanh recurrent", Tanh), ("GRU", Gru)):
    ki7[name] = dict(
        recall_accuracy={f"T={Tl}": round(float(np.mean([recall_task(mk(), Tl, seed=s) for s in (0, 1)])), 3) for Tl in (20, 100, 300)},
        gradient_norm_untrained={f"T={Tl}": float(f"{grad_norm(mk(), Tl):.3e}") for Tl in (20, 100, 300, 1000)})
    print(name, ki7[name], flush=True)
ki7["chance_level"] = 0.125
ki7["note"] = ("Norm-preserving recurrences built from Householder reflections are a known architecture class "
               "(unitary / orthogonal RNNs, Mhammedi et al. 2017); the TFPT pulse operator is one instance of it.")
res["KI7"] = ki7
json.dump(res, open("out/ki4c_ki7_projector_reflection.json", "w"), indent=1)
