"""KI6 - when does E8 symmetry help a learner?

Inputs are tuples of E8 roots (D8+ construction, 240 roots).  Four ways to read them:
  raw        the coordinates as they are
  augment    raw, but every training batch is moved by a random Weyl group element
  gram       all inner products of the tuple (general geometry, invariant under ALL rotations)
  e8         the tuple moved into the E8 fundamental chamber (invariant under exactly W(E8))
  d8         the same with the subgroup W(D8) (permutations and even sign changes)
Tasks:
  A  three roots a, b, c: is a + b - c again a root?   (decided by inner products)
  B  four mutually orthogonal roots: is (a+b+c+d)/2 a lattice point?  E8 has two kinds of
     such quadruples with identical inner products, so general geometry cannot tell them apart.
  C  symmetry broken by a fixed background direction h: is h.(a + b - c) > 0 ?
     variants: h hidden from the invariant models, h given as an input that is moved along,
     and the subgroup that keeps h fixed.
Splits: fresh representatives (new tuples) and, for A, held-out orbit families (inner-product
classes never seen in training).  Invariance is checked on all eight simple reflections.
"""
import itertools, json
import numpy as np, torch, torch.nn as nn
import tfpt_ops as T

torch.set_num_threads(2)
rng = np.random.default_rng(0)
R = T.e8_roots()
res = {}

def simple_roots(roots):
    rho = np.random.default_rng(7).standard_normal(roots.shape[1])
    pos = roots[roots @ rho > 0]
    key = {tuple(np.round(2 * p).astype(int)) for p in pos}
    out = [p for p in pos if not any(tuple(np.round(2 * (p - q)).astype(int)) in key for q in pos)]
    return np.array(out)

S_E8 = simple_roots(R)
S_D8 = simple_roots(R[:112])
h = np.zeros(8); h[0] = 1.0
S_D7 = simple_roots(R[:112][np.abs(R[:112] @ h) < 1e-9])            # keeps h fixed
assert len(S_E8) == 8 and len(S_D8) == 8 and len(S_D7) == 7
W = np.array([1.0, 0.371, 0.1193, 0.0431, 0.0149])

def canon(X, simple):
    """X (n, k, 8): move every tuple by one group element into the closed fundamental chamber."""
    X = X.copy()
    w = W[:X.shape[1]]
    for _ in range(400):
        v = np.einsum("k,nkd->nd", w, X)
        ip = v @ simple.T
        j = ip.argmin(1)
        bad = ip[np.arange(len(X)), j] < -1e-9
        if not bad.any():
            break
        a = simple[j[bad]]
        X[bad] -= np.einsum("nkd,nd->nk", X[bad], a)[:, :, None] * a[:, None, :]
    return X

def weyl(X, g, n=40):
    X = X.copy()
    for _ in range(n):
        a = R[g.integers(240, size=len(X))]
        X -= np.einsum("nkd,nd->nk", X, a)[:, :, None] * a[:, None, :]
    return X

def gram(X):
    k = X.shape[1]
    iu = np.triu_indices(k, 1)
    return np.einsum("nid,njd->nij", X, X)[:, iu[0], iu[1]]

def in_lattice(v):
    f = np.abs(v - np.rint(v)).max(1) < 1e-9
    hf = np.abs(v - 0.5 - np.rint(v - 0.5)).max(1) < 1e-9
    return (f | hf) & (np.abs(v.sum(1) / 2 - np.rint(v.sum(1) / 2)) < 1e-9)

# ------------------------------------------------------------------ data
def task_A(n, g):
    X, y = [], []
    while len(X) < n:
        t = R[g.integers(240, size=(4 * n, 3))]
        lab = np.isclose(((t[:, 0] + t[:, 1] - t[:, 2]) ** 2).sum(1), 2)
        pos, neg = np.flatnonzero(lab), np.flatnonzero(~lab)
        m = min(len(pos), len(neg), n // 2 + 1)
        X += list(t[pos[:m]]) + list(t[neg[:m]]); y += [1] * m + [0] * m
    p = g.permutation(len(X))[:n]
    return np.array(X)[p], np.array(y)[p]

def task_B(n, g):
    X, y = [], []
    G = np.rint(R @ R.T).astype(int)
    while len(X) < n:
        a = g.integers(240)
        c1 = np.flatnonzero(G[a] == 0); b = g.choice(c1)
        c2 = c1[G[b][c1] == 0]; c = g.choice(c2)
        c3 = c2[G[c][c2] == 0]
        lab = in_lattice((R[a] + R[b] + R[c] + R[c3]) / 2)
        want = len(X) % 2
        cand = c3[lab == bool(want)]
        if len(cand):
            d = g.choice(cand)
            X.append(R[[a, b, c, d]][g.permutation(4)]); y.append(want)
    return np.array(X), np.array(y)

def task_C(n, g):
    X = R[g.integers(240, size=(3 * n, 3))]
    s = (X[:, 0] + X[:, 1] - X[:, 2]) @ h
    X, s = X[np.abs(s) > 1e-9][:n], s[np.abs(s) > 1e-9][:n]
    return X, (s > 0).astype(int)

def with_h(X):
    return np.concatenate([X, np.broadcast_to(h, (len(X), 1, 8))], 1)

FEATS = {
    "raw": lambda X: X.reshape(len(X), -1),
    "gram (all rotations)": gram,
    "e8 chamber": lambda X: canon(X, S_E8).reshape(len(X), -1),
    "d8 chamber (subgroup)": lambda X: canon(X, S_D8).reshape(len(X), -1),
}

def fit(Ftr, ytr, Fte_list, seed, aug=None, steps=1500):
    torch.manual_seed(seed)
    net = nn.Sequential(nn.Linear(Ftr.shape[1], 128), nn.ReLU(), nn.Linear(128, 128), nn.ReLU(), nn.Linear(128, 2))
    opt = torch.optim.Adam(net.parameters(), lr=2e-3, weight_decay=1e-4)
    Xt, yt = torch.tensor(Ftr, dtype=torch.float32), torch.tensor(ytr)
    g = np.random.default_rng(seed)
    for _ in range(steps):
        idx = torch.randint(0, len(Xt), (min(128, len(Xt)),))
        xb = Xt[idx]
        if aug is not None:
            xb = torch.tensor(aug(idx.numpy(), g), dtype=torch.float32)
        loss = nn.functional.cross_entropy(net(xb), yt[idx])
        opt.zero_grad(); loss.backward(); opt.step()
    with torch.no_grad():
        return net, [float((net(torch.tensor(F, dtype=torch.float32)).argmax(1) == torch.tensor(y)).float().mean()) for F, y in Fte_list]

def invariance_defect(net, feat, X):
    """mean change of the class-1 probability under the eight simple reflections of E8."""
    with torch.no_grad():
        p0 = torch.softmax(net(torch.tensor(feat(X), dtype=torch.float32)), 1)[:, 1]
        d = []
        for a in S_E8:
            Xr = X - np.einsum("nkd,d->nk", X, a)[:, :, None] * a
            d.append(float((torch.softmax(net(torch.tensor(feat(Xr), dtype=torch.float32)), 1)[:, 1] - p0).abs().mean()))
    return round(float(np.mean(d)), 4)

SIZES, SEEDS = [30, 100, 300, 1000, 3000], [0, 1, 2]
def run(task, feats, tests, extra_aug=True):
    out = {}
    for name, feat in list(feats.items()) + ([("raw + Weyl augmentation", None)] if extra_aug else []):
        row = {}
        for n in SIZES:
            accs, defect = [], []
            for s in SEEDS:
                Xtr, ytr = task(n, np.random.default_rng(100 + s))
                if feat is None:
                    f = FEATS["raw"]
                    net, a = fit(f(Xtr), ytr, [(f(X), y) for X, y in tests], s,
                                 aug=lambda idx, g: weyl(Xtr[idx], g).reshape(len(idx), -1))
                else:
                    f = feat
                    net, a = fit(f(Xtr), ytr, [(f(X), y) for X, y in tests], s)
                accs.append(a); defect.append(invariance_defect(net, f, tests[0][0][:300]))
            m = np.mean(accs, 0)
            row[f"N={n}"] = [round(float(x), 3) for x in m] + [round(float(np.mean(defect)), 4)]
        out[name] = row
        print(name, row, flush=True)
    return out

gt = np.random.default_rng(999)
# --- A: fresh tuples and held-out inner-product classes
XA, yA = task_A(4000, gt)
cls = [tuple(r) for r in np.rint(gram(XA)).astype(int)]
classes = sorted(set(cls))
held = {c for i, c in enumerate(classes) if i % 4 == 0}
mask = np.array([c in held for c in cls])
def task_A_train(n, g):
    X, y = task_A(3 * n + 50, g)
    keep = np.array([tuple(r) not in held for r in np.rint(gram(X)).astype(int)])
    return X[keep][:n], y[keep][:n]
print("A: inner-product classes", len(classes), "held out", len(held), flush=True)
res["A  is a+b-c a root (decided by inner products)"] = dict(
    columns="accuracy on new tuples of seen classes, accuracy on held-out classes, invariance defect",
    n_classes=len(classes), n_held_out=len(held),
    results=run(task_A_train, FEATS, [(XA[~mask], yA[~mask]), (XA[mask], yA[mask])]))
# --- B
XB, yB = task_B(3000, gt)
res["B  which kind of orthogonal quadruple (same inner products)"] = dict(
    columns="accuracy on new quadruples, invariance defect", results=run(task_B, FEATS, [(XB, yB)]))
# --- C
XC, yC = task_C(3000, gt)
featsC = {
    "raw": FEATS["raw"],
    "gram, h hidden": gram,
    "e8 chamber, h hidden": FEATS["e8 chamber"],
    "gram, h given": lambda X: gram(with_h(X)),
    "e8 chamber, h given and moved along": lambda X: canon(with_h(X), S_E8).reshape(len(X), -1),
    "chamber of the subgroup that fixes h (D7)": lambda X: canon(X, S_D7).reshape(len(X), -1),
}
res["C  background direction breaks the symmetry"] = dict(
    columns="accuracy on new triples, invariance defect (non-zero is CORRECT here: the task is not invariant)",
    results=run(task_C, featsC, [(XC, yC)], extra_aug=False))
json.dump(res, open("out/ki6_e8_equivariance.json", "w"), indent=1)
