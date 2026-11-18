"""E6 - "E8-based neural representations" (chat priority 3).

What E8 concretely offers a network is its lattice: the densest packing and best known
quantiser in 8 dimensions.  Three measurements at equal point density (unit covolume):
  (a) quantisation error (normalised second moment) of Z^8, D8, E8
  (b) discrete memory in a noisy recurrent loop: a latent token is perturbed and
      re-snapped to the lattice every step; how much noise survives 1000 steps?
  (c) a small autoencoder on handwritten digits with a 16-d latent (two 8-d blocks)
      quantised by each lattice (straight-through), plus an unquantised reference.
"""
import json, sys
import numpy as np, torch, torch.nn as nn
import tfpt_ops as T

LATS = ["Z8", "D8", "E8"]
rng = np.random.default_rng(0)
res = {}
quick = "--quick" in sys.argv

# (a) normalised second moment G = E|x - Q(x)|^2 / 8 at unit covolume
x = rng.uniform(-6, 6, (400_000 if quick else 4_000_000, 8))
res["second_moment_G"] = {l: round(float(((x - T.nearest(x, l)) ** 2).sum(1).mean() / 8), 5) for l in LATS}
res["second_moment_reference"] = dict(Z8=round(1 / 12, 5), D8=0.07591, E8=0.07168)

# (b) per-step symbol error under Gaussian noise, and survival over a 1000-step loop
origin = {"Z8": np.zeros(8), "D8": np.zeros(8), "E8": np.zeros(8)}
sig = [0.10, 0.12, 0.14, 0.16, 0.18, 0.20]
n = 300_000 if quick else 3_000_000
perr = {l: {} for l in LATS}
for s in sig:
    noise = rng.standard_normal((n, 8)) * s
    for l in LATS:
        perr[l][s] = float((np.abs(T.nearest(noise, l)).sum(1) > 1e-9).mean())
res["per_step_error"] = {l: {str(s): round(p, 7) for s, p in v.items()} for l, v in perr.items()}
res["survival_1000_steps_with_resnap"] = {l: {str(s): round((1 - p) ** 1000, 4) for s, p in v.items()}
                                          for l, v in perr.items()}
# the same loop without any snapping: noise accumulates as sqrt(steps)
acc = rng.standard_normal((20000, 8)) * 0.14 * np.sqrt(1000)
res["survival_1000_steps_without_resnap_sigma0.14"] = {
    l: round(float((np.abs(T.nearest(acc, l)).sum(1) < 1e-9).mean()), 4) for l in LATS}
res["min_distance_unit_covolume"] = dict(Z8=1.0, D8=round(float(np.sqrt(2) * 2 ** (-1 / 8)), 4),
                                         E8=round(float(np.sqrt(2)), 4))

# (c) lattice-quantised autoencoder on digits
from sklearn.datasets import load_digits
X = load_digits().data.astype(np.float32) / 16.0
perm = np.random.default_rng(1).permutation(len(X))
Xtr, Xte = torch.tensor(X[perm[:1400]]), torch.tensor(X[perm[1400:]])

def quant(z, lat, s):
    if lat == "none":
        return z
    zz = z.detach().numpy().reshape(-1, 8) / s
    q = torch.tensor(T.nearest(zz, lat).reshape(z.shape) * s, dtype=torch.float32)
    return z + (q - z).detach()

def run_ae(lat, s, seed, steps):
    torch.manual_seed(seed)
    enc = nn.Sequential(nn.Linear(64, 128), nn.GELU(), nn.Linear(128, 16))
    dec = nn.Sequential(nn.Linear(16, 128), nn.GELU(), nn.Linear(128, 64))
    opt = torch.optim.Adam(list(enc.parameters()) + list(dec.parameters()), lr=2e-3)
    for _ in range(steps):
        i = torch.randint(0, len(Xtr), (128,))
        z = torch.tanh(enc(Xtr[i])) * 2
        zq = quant(z, lat, s)
        loss = ((dec(zq) - Xtr[i]) ** 2).mean() + 0.25 * ((z - zq.detach()) ** 2).mean()
        opt.zero_grad(); loss.backward(); opt.step()
    with torch.no_grad():
        zq = quant(torch.tanh(enc(Xte)) * 2, lat, s)
        mse = ((dec(zq) - Xte) ** 2).mean().item()
        codes = len({tuple(np.round(r, 4)) for r in zq.numpy()})
    return mse, codes

torch.set_num_threads(4)
steps, seeds = (300, [0]) if quick else (3000, [0, 1, 2])
ae = {}
for s in (2.0, 1.0, 0.5):
    for lat in LATS:
        r = [run_ae(lat, s, sd, steps) for sd in seeds]
        ae[f"{lat} scale={s}"] = dict(test_mse=round(float(np.mean([a for a, _ in r])), 5),
                                      std=round(float(np.std([a for a, _ in r])), 5),
                                      distinct_codes_on_397_test=int(np.mean([b for _, b in r])))
        print(f"{lat} scale={s}", ae[f"{lat} scale={s}"], flush=True)
r = [run_ae("none", 1.0, sd, steps) for sd in seeds]
ae["unquantised"] = dict(test_mse=round(float(np.mean([a for a, _ in r])), 5), std=round(float(np.std([a for a, _ in r])), 5))
res["digits_autoencoder"] = ae
json.dump(res, open("out/e6_e8_representations.json", "w"), indent=1)
print(json.dumps({k: v for k, v in res.items() if k != "digits_autoencoder"}, indent=1))
