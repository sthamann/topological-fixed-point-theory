"""E4 - "reasoning as a search for stable states" (chat 3.3), with ZERO learned parameters.

Unknown cells carry a distribution; we descend a fixed energy summed over all four-mark
cells and read off the stable state.  Compared energies, all taken from the TFPT operators:
  Qc    1 - Pr[xor = 0]            classical shadow of the Pauli stabiliser projector 4 Q_P
  Qamp  1 - 4 <Q_P>                the same projector on amplitudes sqrt(p)
  Pamp  1 - <P>                    rank-5 quartic code projector on amplitudes
  Hmom  1 - 20 <M4> = <H_mom>      the literal protection energy
and for Sudoku-4 the fifth code vector c4 (all four symbols distinct):
  ADc   1 - Pr[all different]      classical shadow of 24 |c4><c4|
  c4amp 1 - 24 <c4|psi>^2
"""
import json, sys, itertools, collections
import numpy as np, torch, torch.nn.functional as F
import worlds as W, tfpt_ops as T

torch.set_num_threads(8)
ops = T.torch_ops()
STEPS = 400
NS = 60 if "--quick" in sys.argv else 300

ENERGY = {
    "Qc": lambda p: 1 - T.prob_xor0(ops["had"], p),
    "Qamp": lambda p: 1 - 4 * T.exp_Q(ops["pauli"], torch.sqrt(p + 1e-8)),
    "Pamp": lambda p: 1 - T.exp_P(ops["V4"], torch.sqrt(p + 1e-8)),
    "Hmom": lambda p: 1 - T.exp_20M4(ops["stab"], torch.sqrt(p + 1e-8)),
    "ADc": lambda p: 1 - T.prob_alldiff(p),
    "c4amp": lambda p: 1 - 24 * T.code_feat(ops["V4"], torch.sqrt(p + 1e-8))[..., 4] ** 2,
}


def descend(d, quads, energy, active=None, steps=STEPS, lr=0.3, seed=0):
    noise = 0.01 if energy.endswith("+noise") else 0.0   # exact uniform start vs. tiny symmetry breaking
    energy = energy.replace("+noise", "")
    """Returns (solved fraction, median steps to first solution among solved)."""
    torch.manual_seed(seed)
    x = torch.as_tensor(d["x"]).long(); known = torch.as_tensor(d["known"]).bool()
    quads = torch.as_tensor(quads)
    act = torch.ones(len(x), len(quads)) if active is None else torch.as_tensor(active).float()
    xo = F.one_hot(x, 4).float()
    th = (noise * torch.randn(*x.shape, 4)).requires_grad_(True)
    opt = torch.optim.Adam([th], lr=lr)
    first = torch.full((len(x),), -1)
    for t in range(1, steps + 1):
        p = torch.where(known[..., None], xo, F.softmax(th, -1))
        E = (ENERGY[energy](p[:, quads]) * act).sum()
        opt.zero_grad(); E.backward(); opt.step()
        if t % 5 == 0:
            ok = (torch.where(known, x, th.argmax(-1)) == x).all(1)
            first[(first < 0) & ok] = t
    ok = (torch.where(known, x, th.argmax(-1)) == x).all(1)
    st = first[ok & (first > 0)]
    return round(ok.float().mean().item(), 4), (int(st.median()) if len(st) else None)


# ---- Sudoku 4x4 -------------------------------------------------------------
SUD_QUADS = W.FAM[:3].reshape(12, 4)
def sudoku_grids():
    out = []
    def rec(g):
        if len(g) == 16:
            out.append(list(g)); return
        i = len(g)
        for v in range(4):
            if all(g[j] != v for q in SUD_QUADS if i in q for j in q if j < i):
                rec(g + [v])
    rec([])
    return np.array(out)

def count_sol(g, known, limit=2):
    cells = [i for i in range(16) if not known[i]]
    g = g.copy(); cnt = 0
    peers = {i: sorted({j for q in SUD_QUADS if i in q for j in q if j != i}) for i in range(16)}
    kn = known.copy()
    def rec(k):
        nonlocal cnt
        if cnt >= limit: return
        if k == len(cells):
            cnt += 1; return
        i = cells[k]
        for v in range(4):
            if all(not kn[j] or g[j] != v for j in peers[i]):
                g[i] = v; kn[i] = True; rec(k + 1); kn[i] = False
    rec(0)
    return cnt

def sudoku_set(m, ks, seed):
    rng = np.random.default_rng(seed); G = sudoku_grids(); assert len(G) == 288
    X, K = [], []
    while len(X) < m:
        g = G[rng.integers(288)]; k = int(rng.choice(ks))
        known = np.ones(16, bool); known[rng.choice(16, k, replace=False)] = False
        if count_sol(g, known) == 1:
            X.append(g); K.append(known)
    return dict(x=np.array(X), known=np.array(K))


if __name__ == "__main__":
    res = collections.OrderedDict()
    xor_en = ["Qc", "Qamp", "Pamp", "Hmom", "Qc+noise", "Qamp+noise"]
    # (a) 4x4 rule-combination worlds, any combination, zero training
    for name, ks in (("grid4 k=1-5", [1, 2, 3, 4, 5]), ("grid4 k=6-7", [6, 7])):
        d = W.cached(f"e4_grid_{ks[0]}", lambda: W.make_grid_set(NS, W.TRAIN_COMBOS + W.OOD_COMBOS, ks, "xor", 40 + ks[0]))
        res[name] = {e: descend(d, W.QUADS16, e, d["active"]) for e in xor_en}
        print(name, res[name], flush=True)
    # (b) plaquette lattices of growing size, incl. puzzles that local propagation cannot solve
    for n in (6, 12, 24):
        for kind in ("random", "chain", "global"):
            m = NS if n < 24 else NS // 3
            d = W.cached(f"e4_lat_{n}_{kind}_{m}", lambda: W.make_lattice_set(m, n, "xor", 50 + n, kind))
            name = f"lattice n={n} {kind}"
            res[name] = {e: descend(d, W.lattice_quads(n), e, steps=STEPS if n < 24 else 2 * STEPS) for e in xor_en}
            res[name]["unit_propagation"] = [0.0 if kind == "global" else 1.0, int(np.median(d["depth"]))]
            print(name, res[name], flush=True)
    # (c) wrong group: the same energies on Z4 worlds must fail
    d = W.cached("e4_lat_z4", lambda: W.make_lattice_set(NS, 6, "z4", 77, "random"))
    res["lattice n=6 random, Z4 worlds (mismatch)"] = {e: descend(d, W.lattice_quads(6), e) for e in ("Qc", "Pamp")}
    print(res["lattice n=6 random, Z4 worlds (mismatch)"], flush=True)
    # (d) Sudoku 4x4 through the fifth code vector
    for name, ks in (("sudoku4 k=4-7", [4, 5, 6, 7]), ("sudoku4 k=8-11", [8, 9, 10, 11])):
        d = W.cached(f"e4_sud_{ks[0]}", lambda: sudoku_set(NS, ks, 60 + ks[0]))
        res[name] = {e: descend(d, SUD_QUADS, e) for e in ("ADc", "c4amp", "ADc+noise", "Qc", "Hmom")}
        print(name, res[name], flush=True)
    json.dump(res, open("out/e4_energy_reasoning.json", "w"), indent=1)
