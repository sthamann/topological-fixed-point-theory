"""E0 - smallest decisive structural test, before any training.

Question: what do the literal TFPT operators (H_mom, P, Q_P) do on the states a
classical network can actually hold, i.e. product states / definite 4-tuples?
"""
import itertools, json, collections
import numpy as np, torch
import tfpt_ops as T

o = T.build()
res = {"identities": "60 rays; ranks S,Q,P = 35,16,5; 40 M4 = S+P; spec(H_mom) = 0^5, (1/2)^30, 1^221"}

# 1. classical shadow on the 256 definite 4-tuples
def pattern(w):
    c = sorted(collections.Counter(w).values(), reverse=True)
    return {(4,): "xxxx", (3, 1): "xxxy", (2, 2): "xxyy", (2, 1, 1): "xxyz", (1, 1, 1, 1): "xyzw"}[tuple(c)]

rows = collections.defaultdict(list)
for idx, w in enumerate(itertools.product(range(4), repeat=4)):
    xor0 = (w[0] ^ w[1] ^ w[2] ^ w[3]) == 0
    rows[pattern(w)].append((xor0, o["S"][idx, idx], o["P"][idx, idx], o["Q"][idx, idx], o["Hmom"][idx, idx]))
tab = {}
for k, v in rows.items():
    a = np.array(v, dtype=float)
    assert np.allclose(a, a[0]), k
    tab[k] = dict(n=len(v), xor_zero=bool(a[0, 0]), S=round(a[0, 1], 5), P=round(a[0, 2], 5),
                  Q_P=round(a[0, 3], 5), H_mom=round(a[0, 4], 5), one_minus_P=round(1 - a[0, 2], 5))
res["definite_tuples"] = tab
supp = [(w[0] ^ w[1] ^ w[2] ^ w[3]) == 0 for w in itertools.product(range(4), repeat=4)]
res["Q_diag_equals_quarter_xor_indicator"] = bool(np.allclose(np.diag(o["Q"]), 0.25 * np.array(supp)))
res["P_support_equals_xor_zero"] = bool(np.array_equal(np.diag(o["P"]) > 1e-12, np.array(supp)))
res["B_Bt_diag_equals_xor_indicator"] = bool(np.allclose(np.diag(o["B"] @ o["B"].T), np.array(supp)))

# 2. parameter-free transition: contraction of P with three definite registers
ops = T.torch_ops()
ok = 0
for a, b, c in itertools.product(range(4), repeat=3):
    z = torch.zeros(4, 4); z[0, a] = z[1, b] = z[2, c] = 1; z[3] = 0.5
    mu, _ = T.code_messages(ops["V4"], z)
    ok += int(mu[3].argmax().item() == (a ^ b ^ c) and (mu[3] > 1e-9).sum().item() == 1)
res["P_contraction_completes_xor"] = f"{ok}/64"

# 3. H_mom on product states: minimum, and what the minimisers are
torch.manual_seed(0)
n = 2000
def prod_energy(zraw, tied):
    z = zraw / zraw.norm(dim=-1, keepdim=True)
    z4 = z.expand(-1, 4, -1) if tied else z
    return 1 - T.exp_20M4(ops["stab"], z4)

out = {}
for tied in (True, False):
    zr = torch.randn(n, 1 if tied else 4, 4, dtype=torch.complex64, requires_grad=True)
    e0 = prod_energy(zr, tied).detach()
    opt = torch.optim.Adam([zr], lr=0.05)
    for _ in range(1500):
        opt.zero_grad(); e = prod_energy(zr, tied); e.sum().backward(); opt.step()
    e = prod_energy(zr, tied).detach()
    z = (zr / zr.norm(dim=-1, keepdim=True)).detach()
    ov = (torch.einsum("la,nia->nil", ops["stab"].conj(), z).abs() ** 2).max(-1).values.min(-1).values
    out["tied z^x4" if tied else "independent z1..z4"] = dict(
        random_start_mean=round(e0.mean().item(), 4), min_found=round(e.min().item(), 5),
        frac_at_min=round((e < e.min() + 1e-3).float().mean().item(), 3),
        frac_all_registers_stabilizer=round((ov > 0.999).float().mean().item(), 3))
res["H_mom_on_product_states"] = out
res["H_mom_stabilizer_copy_value"] = round(float(1 - T.exp_20M4(ops["stab"], ops["stab"][0].expand(4, 4))), 5)

json.dump(res, open("out/e0_structure.json", "w"), indent=1)
print(json.dumps(res, indent=1))
