"""E7 - "relational world models with consistency checking" (chat priority 1).

A world model imagines long rollouts.  World: periodic n x n lattice whose 2x2 windows are
four-mark cells (xor defect must vanish).  Actions shift a whole row or column, which keeps
the world valid.  The imagined transition is imperfect: after every step each cell is
corrupted with probability eps (a controlled stand-in for the error of a learned model;
no network is trained here).  Three ways to imagine H steps ahead:
  none     trust the model
  detect   compute the stabiliser defects; count how many wrong worlds are noticed
  repair   local stable-state search: a cell whose windows agree on the same nonzero defect
           value is corrected by it (the classical syndrome of 4 Q_P on overlapping cells)
Control: 'repair_z4' applies the same procedure with the wrong group.
"""
import json, sys
import numpy as np

n, B = 8, (400 if "--quick" in sys.argv else 4000)
rng = np.random.default_rng(0)
R, C = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")
# window q=(r,c) covers (r,c),(r,c+1),(r+1,c),(r+1,c+1) periodically
WIN = np.stack([R * n + C, R * n + (C + 1) % n, ((R + 1) % n) * n + C, ((R + 1) % n) * n + (C + 1) % n], -1).reshape(-1, 4)
# cell i=(r,c) lies in the windows (r,c),(r,c-1),(r-1,c),(r-1,c-1)
CELLWIN = np.stack([R * n + C, R * n + (C - 1) % n, ((R - 1) % n) * n + C, ((R - 1) % n) * n + (C - 1) % n], -1).reshape(-1, 4)


def syn(x, group):
    v = x[:, WIN]
    return np.bitwise_xor.reduce(v, -1) if group == "xor" else v.sum(-1) % 4


def repair(x, group, sweeps=8):
    x = x.copy()
    for _ in range(sweeps):
        s = syn(x, group)
        if not s.any():
            break
        cs = s[:, CELLWIN]                                   # B,N,4 defect values around each cell
        for val in (1, 2, 3):
            hit = (cs == val).sum(-1) >= 3
            if group == "xor":
                x = np.where(hit, x ^ val, x)
            else:
                x = np.where(hit, (x - val) % 4, x)
    return x


def rollout(eps, H, mode):
    a, b = rng.integers(0, 4, (B, n)), rng.integers(0, 4, (B, n))
    true = (a[:, :, None] ^ b[:, None, :]).reshape(B, -1)
    img = true.copy()
    noticed = wrong = 0
    out = {}
    for t in range(1, H + 1):
        k, v = rng.integers(0, 2 * n, B), rng.integers(1, 4, B)
        shift = np.zeros((B, n, n), dtype=np.int64)
        rows = np.arange(B)
        isrow = k < n
        shift[rows[isrow], k[isrow], :] = v[isrow, None]
        shift[rows[~isrow], :, k[~isrow] - n] = v[~isrow, None]
        shift = shift.reshape(B, -1)
        true, img = true ^ shift, img ^ shift
        err = rng.random(img.shape) < eps
        img = np.where(err, img ^ rng.integers(1, 4, img.shape), img)
        if mode == "detect":
            bad = (img != true).any(1)
            wrong += bad.sum(); noticed += (bad & syn(img, "xor").any(1)).sum()
        if mode == "repair":
            img = repair(img, "xor")
        if mode == "repair_z4":
            img = repair(img, "z4")
        if t in (10, 50, 200):
            out[f"exact_after_{t}"] = round(float((img == true).all(1).mean()), 4)
    if mode == "detect":
        out["wrong_worlds_noticed"] = round(float(noticed / max(wrong, 1)), 5)
    return out


res = {"setup": dict(n=n, cells=n * n, windows=len(WIN), rollouts=B)}
for eps in (0.001, 0.003, 0.01, 0.03):
    res[f"eps={eps}"] = {m: rollout(eps, 200, m) for m in ("none", "detect", "repair", "repair_z4")}
    print(f"eps={eps}", json.dumps(res[f"eps={eps}"]), flush=True)
json.dump(res, open("out/e7_world_model_consistency.json", "w"), indent=1)
