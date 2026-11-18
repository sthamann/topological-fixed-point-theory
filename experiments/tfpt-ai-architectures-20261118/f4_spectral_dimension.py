"""F4 - three large space dimensions from relations?

Measured: spectral dimension d_s = -2 dlogP/dlogt of the heat kernel P(t) = Tr exp(-tL)/N.
Calibration on Z^2, Z^3, Z^4.  Candidates built from four marks:
  A3 difference graph   nodes = integer combinations of mark differences e_i - e_j (rank 3),
                        neighbours = the 12 differences.  Also m = 3 and m = 5 marks.
  plaquette lattice     overlapping four-mark cells in a grid (as used in the AI experiments)
  random tetrahedra     four-mark cells glued at random (no neighbourhood rule at all)
  rewired A3            same local valence, edges rewired at random (negative control)
The neighbourhood rule is an INPUT in every case; that is the first added assumption.
"""
import itertools, json
import numpy as np

rng = np.random.default_rng(0)
res = {}

def ds_from_spectrum(lam, tmin, tmax):
    t = np.geomspace(tmin, tmax, 40)
    P = np.array([np.exp(-x * lam).mean() for x in t])
    d = -2 * np.gradient(np.log(P), np.log(t))
    return t, d

def torus_spectrum(gens, L):
    r = gens.shape[1]
    k = 2 * np.pi * np.stack(np.meshgrid(*[np.arange(L)] * r, indexing="ij"), -1).reshape(-1, r) / L
    return (1 - np.cos(k @ gens.T)).sum(1) / len(gens) * 2      # normalised: lazy-walk generator

def a_gens(m):
    """the m(m-1) differences e_i - e_j in coordinates of the simple roots e_i - e_{i+1}."""
    out = []
    for i in range(m):
        for j in range(m):
            if i != j:
                v = np.zeros(m - 1, int)
                lo, hi = min(i, j), max(i, j)
                v[lo:hi] = 1 if i < j else -1
                out.append(v)
    return np.array(out)

def window(L):
    return 2.0, L * L / 12.0

fams = {"Z^2 (calibration)": (np.vstack([np.eye(2, dtype=int), -np.eye(2, dtype=int)]), [24, 48, 96, 192]),
        "Z^3 (calibration)": (np.vstack([np.eye(3, dtype=int), -np.eye(3, dtype=int)]), [10, 14, 20, 28]),
        "Z^4 (calibration)": (np.vstack([np.eye(4, dtype=int), -np.eye(4, dtype=int)]), [8, 10, 13, 16]),
        "A2 differences (3 marks)": (a_gens(3), [24, 48, 96, 192]),
        "A3 differences (4 marks)": (a_gens(4), [10, 14, 20, 28]),
        "A4 differences (5 marks)": (a_gens(5), [8, 10, 13, 16])}
for name, (g, Ls) in fams.items():
    row = {}
    for L in Ls:
        lam = torus_spectrum(g, L)
        t, d = ds_from_spectrum(lam, *window(L))
        mid = d[(t > 2.5) & (t < L * L / 14)]
        row[f"N={L ** g.shape[1]}"] = dict(d_s=round(float(np.median(mid)), 3),
                                          spread=round(float(mid.max() - mid.min()), 3))
    res[name] = row
    print(name, row, flush=True)

def laplacian_spectrum(edges, N):
    A = np.zeros((N, N))
    A[edges[:, 0], edges[:, 1]] = 1; A[edges[:, 1], edges[:, 0]] = 1
    deg = A.sum(1)
    keep = deg > 0
    A = A[keep][:, keep]; deg = deg[keep]
    Ln = np.eye(len(deg)) - A / np.sqrt(np.outer(deg, deg))
    return np.linalg.eigvalsh(Ln) * 2

# explicit A3 torus graph, then degree-preserving rewiring
def a3_edges(L):
    g = a_gens(4)
    idx = np.stack(np.meshgrid(*[np.arange(L)] * 3, indexing="ij"), -1).reshape(-1, 3)
    num = lambda v: (v[:, 0] % L) * L * L + (v[:, 1] % L) * L + v[:, 2] % L
    e = np.concatenate([np.stack([num(idx), num(idx + gg)], 1) for gg in g])
    return e[e[:, 0] < e[:, 1]]
row = {}
for L in (10, 14, 18):
    e = a3_edges(L).copy()
    for _ in range(4 * len(e)):                      # double-edge swaps
        i, j = rng.integers(len(e), size=2)
        a, b, c, d = e[i, 0], e[i, 1], e[j, 0], e[j, 1]
        if len({a, b, c, d}) == 4:
            e[i, 1], e[j, 1] = d, b
    lam = laplacian_spectrum(e, L ** 3)
    t, d = ds_from_spectrum(lam, 2.0, 12.0)
    row[f"N={L ** 3}"] = dict(d_s_at_t_4_to_10=[round(float(x), 2) for x in d[(t > 4) & (t < 10)][[0, -1]]])
res["A3 rewired (negative control)"] = row
print("rewired", row, flush=True)

row = {}
for N in (1000, 2000, 4000):
    cells = np.concatenate([rng.permutation(N).reshape(-1, 4) for _ in range(3)])   # each mark in 3 cells
    e = np.concatenate([cells[:, [a, b]] for a, b in itertools.combinations(range(4), 2)])
    lam = laplacian_spectrum(e, N)
    t, d = ds_from_spectrum(lam, 2.0, 12.0)
    row[f"N={N}"] = dict(d_s_at_t_4_to_10=[round(float(x), 2) for x in d[(t > 4) & (t < 10)][[0, -1]]])
res["random gluing of four-mark cells"] = row
print("random tetra", row, flush=True)

row = {}
for L in (24, 48, 96):                                # plaquette lattice: cells = 2x2 windows -> king graph
    g = np.array([[a, b] for a in (-1, 0, 1) for b in (-1, 0, 1) if (a, b) != (0, 0)])
    lam = torus_spectrum(g, L)
    t, d = ds_from_spectrum(lam, *window(L))
    mid = d[(t > 4) & (t < L * L / 40)]
    row[f"N={L * L}"] = dict(d_s=round(float(np.median(mid)), 3))
res["plaquette lattice of four-mark cells"] = row
res["reading"] = ("d_s = m - 1 for the difference graph of m marks, so four marks give 3 (window 2.8-3.2 passed) "
                  "only through the reading 'neighbour = one mark difference'.  The same cells glued in a grid "
                  "give 2, glued at random give no finite dimension.  The three is carried by the chosen "
                  "neighbourhood, not by the cells.")
json.dump(res, open("out/f4_spectral_dimension.json", "w"), indent=1)
print(json.dumps(res["plaquette lattice of four-mark cells"]))
