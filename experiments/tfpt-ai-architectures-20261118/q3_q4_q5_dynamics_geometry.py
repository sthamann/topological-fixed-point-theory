"""Q3, Q4, Q5 - pulse dynamics, ten equiangular lines, distance from entanglement.

Q3  Two quartic cells with H = H_cell(L) + H_cell(R) + g D_A(L) D_B(R).  Source claims: after
    T0 = pi / sqrt(Delta^2 + 256 g^2) the population outside the code returns to zero; with
    g/Delta = sqrt(2m-1) / (16 (m-1)) and m cycles the active block gets an exact minus sign;
    the largest leak on the way is (2m-1)/m^2.  Also: sensitivity to a miscalibrated coupling,
    and the pulse reflection exp(-i pi H/2) J0 = W for H = 2 D D^+.
Q4  The ten equiangular lines of the pair measurement: C^2 = 9I, the six Petersen sign frames,
    optimality (Welch bound), comparison with freely optimised and random frames, and use as
    a superposition memory (10 items in 5 dimensions).
Q5  Distance from entanglement: logical cells on a ring or a line with the exchange -K on
    neighbours.  Mutual information between cells against graph distance.
"""
import itertools, json
import numpy as np
import scipy.linalg as sla
from scipy.sparse.linalg import LinearOperator, eigsh
import tfpt_ops as T

rng = np.random.default_rng(2)
o = T.build()
V, P, S, pauli = o["V"], o["P"], o["S"], o["pauli"]
res = {}

def on(A, r):
    ops = [np.eye(4)] * 4; ops[r] = A
    return T.kron4(*ops)
def D_of(A):
    return sum(on(A, r) for r in range(4)).real if np.abs(A.imag).max() < 1e-12 else sum(on(A, r) for r in range(4))

# ---------------------------------------------------------------- Q3
Delta = 1.0
Hcell = Delta * (2 * np.eye(256) - P - S)
DA, DB = D_of(pauli[1]), D_of(pauli[15])                 # IX on the left cell, ZZ on the right
def cell_basis(D):
    Bq = np.linalg.qr(np.hstack([V, D @ V]))[0]
    rank = np.linalg.matrix_rank(np.hstack([V, D @ V]), tol=1e-9)
    Bq = sla.orth(np.hstack([V, D @ V]))
    return Bq, rank
BL, rL = cell_basis(DA); BR, rR = cell_basis(DB)
inv_err = max(np.abs(M @ B - B @ (B.T @ M @ B)).max() for B, M in ((BL, Hcell), (BL, DA), (BR, Hcell), (BR, DB)))
hL, dL, hR, dR = BL.T @ Hcell @ BL, BL.T @ DA @ BL, BR.T @ Hcell @ BR, BR.T @ DB @ BR
cL, cR = BL.T @ V, BR.T @ V                               # code inside the invariant spaces
Pc = np.kron(cL @ cL.T, cR @ cR.T)
def evolve(g, t):
    H = np.kron(hL, np.eye(len(hR))) + np.kron(np.eye(len(hL)), hR) + g * np.kron(dL, dR)
    return sla.expm(-1j * t * H)
def max_leak(U):
    Vc = np.kron(cL, cR)
    L = Vc.T @ U.conj().T @ (np.eye(len(Pc)) - Pc) @ U @ Vc
    return float(np.linalg.eigvalsh((L + L.conj().T) / 2).max())
q3 = {"invariant_subspace_dim_per_cell": [int(rL), int(rR)], "invariance_error": float(inv_err)}
for m in (2, 3, 5):
    g = Delta * np.sqrt(2 * m - 1) / (16 * (m - 1))
    T0 = np.pi / np.sqrt(Delta ** 2 + 256 * g ** 2)
    Ttot = m * T0
    U = evolve(g, Ttot)
    Vc = np.kron(cL, cR)
    Ulog = Vc.T @ U @ Vc
    ph = np.linalg.eigvals(Ulog)
    ts = np.linspace(0, Ttot, 400)
    q3[f"m={m}"] = dict(
        g_over_Delta=round(float(g), 6), total_time_vs_formula=[round(float(Ttot), 8), round(float((m - 1) * np.pi / Delta), 8)],
        leak_at_end=float(f"{max_leak(U):.2e}"), leak_after_one_cycle=float(f"{max_leak(evolve(g, T0)):.2e}"),
        max_leak_on_the_way=round(max(max_leak(evolve(g, t)) for t in ts), 6), formula_2m_minus_1_over_m2=round((2 * m - 1) / m ** 2, 6),
        logical_eigenphases={"+1": int((np.abs(ph - 1) < 1e-6).sum()), "-1": int((np.abs(ph + 1) < 1e-6).sum())},
        leak_at_end_if_g_off_by_1_percent=float(f"{max_leak(evolve(1.01 * g, Ttot)):.2e}"),
        leak_at_end_if_g_off_by_5_percent=float(f"{max_leak(evolve(1.05 * g, Ttot)):.2e}"))
# pulse reflection with the documented overlap s = sqrt(2/5): J0^+ W = s I
s = np.sqrt(2 / 5)
Qm = np.linalg.qr(rng.standard_normal((75, 10)))[0]
J0, Jp = Qm[:, :5], Qm[:, 5:]
Wev = s * J0 + np.sqrt(1 - s * s) * Jp
Dm = (J0 - Wev) / np.sqrt(2 * (1 - s))
Hp = 2 * Dm @ Dm.T
ev = np.round(np.linalg.eigvalsh(Hp), 8)
q3["pulse_reflection"] = dict(spectrum={"0": int((ev == 0).sum()), "2": int((ev == 2).sum())},
                              error_of_exp_minus_i_pi_H_over_2_J0_equals_W=float(np.abs(sla.expm(-1j * np.pi * Hp / 2) @ J0 - Wev).max()),
                              scope="any pair of isometries with J0^+ W = s I; the actual encoder W_ev was not rebuilt")
res["Q3"] = q3

# ---------------------------------------------------------------- Q4
F_num = np.array([[4, 4, 4, 4, 0], [0, 8, 0, 0, 8], [4, -4, 4, -4, 0], [0, 0, 8, 0, 8], [0, 0, 0, 8, 8],
                  [0, 0, 8, 0, -8], [0, 0, 0, 8, -8], [4, 4, -4, -4, 0], [0, 8, 0, 0, -8], [4, -4, -4, 4, 0]], float)
G = np.array([4, 12, 12, 12, 24.0])
W10 = F_num / 4 / np.sqrt(G)
Gam = W10 @ W10.T
C = 6 * Gam - 3 * np.eye(10)
frames = [sg for sg in itertools.product((1, -1), repeat=9) if np.allclose(C @ np.array((1,) + sg), 3 * np.array((1,) + sg))]
pet = []
for sg in frames:
    s_ = np.array((1,) + sg); A = (np.ones((10, 10)) - np.eye(10) - np.diag(s_) @ C @ np.diag(s_)) / 2
    pet.append(bool(np.allclose(A.sum(1), 3) and np.allclose(A @ A, 2 * np.eye(10) + np.ones((10, 10)) - A)))
unit = W10 * np.sqrt(2)
def coherence(Fm):
    Gm = np.abs(Fm @ Fm.T); np.fill_diagonal(Gm, 0); return float(Gm.max())
def optimise():
    X = rng.standard_normal((10, 5))
    lr = 0.05
    for it in range(4000):
        X /= np.linalg.norm(X, axis=1, keepdims=True)
        Gm = X @ X.T; np.fill_diagonal(Gm, 0)
        grad = 4 * (Gm ** 7) @ X                       # minimise sum of 8th powers of overlaps
        X -= lr * grad
    X /= np.linalg.norm(X, axis=1, keepdims=True)
    return X
opt = [coherence(optimise()) for _ in range(20)]
rnd = [coherence((lambda X: X / np.linalg.norm(X, axis=1, keepdims=True))(rng.standard_normal((10, 5)))) for _ in range(200)]
def memory(Fm, sigma, trials=20000):
    i = rng.integers(10, size=trials)
    y = Fm[i] + sigma * rng.standard_normal((trials, 5))
    return float((np.abs(y @ Fm.T).argmax(1) == i).mean())
mem = {}
for sigma in (0.2, 0.35, 0.5):
    rf = [memory((lambda X: X / np.linalg.norm(X, axis=1, keepdims=True))(rng.standard_normal((10, 5))), sigma, 4000) for _ in range(50)]
    mem[f"noise={sigma}"] = dict(tfpt_lines=round(memory(unit, sigma), 4), optimised_frame=round(memory(optimise(), sigma), 4),
                                 random_frames_mean=round(float(np.mean(rf)), 4))
res["Q4"] = dict(W_is_isometry=bool(np.allclose(W10.T @ W10, np.eye(5))),
                 off_diagonal_abs_values=sorted({round(float(abs(x)), 6) for x in Gam[~np.eye(10, dtype=bool)]}),
                 C_squared_equals_9I=bool(np.allclose(C @ C, 9 * np.eye(10))), sign_frames_with_Cs_3s=len(frames),
                 all_are_Petersen_graphs=all(pet), coherence_tfpt=round(coherence(unit), 6),
                 welch_bound=round(float(np.sqrt((10 - 5) / (5 * 9))), 6),
                 coherence_freely_optimised_20_starts=[round(min(opt), 4), round(float(np.median(opt)), 4), round(max(opt), 4)],
                 coherence_random_frames_mean=round(float(np.mean(rnd)), 4), recall_accuracy_10_items_in_5_dims=mem)

# ---------------------------------------------------------------- Q5
R = []
for A in pauli[1:]:
    D = D_of(A)
    R.append((V.T @ D @ D @ V).real / 16)
K4 = sum(np.kron(r, r) for r in R).reshape(5, 5, 5, 5)
def ground(n, edges):
    shape = (5,) * n
    def mv(v):
        x = v.reshape(shape); out = np.zeros_like(x)
        for i, j in edges:
            y = np.tensordot(K4, x, axes=([2, 3], [i, j]))       # result axes: (i', j', rest)
            out -= np.moveaxis(y, [0, 1], [i, j])
        return out.reshape(-1)
    w, v = eigsh(LinearOperator((5 ** n, 5 ** n), matvec=mv, dtype=float), k=3, which="SA")
    order = np.argsort(w)
    return w[order], v[:, order[0]].reshape(shape)
def entropy(rho):
    w = np.linalg.eigvalsh(rho); w = w[w > 1e-12]
    return float(-(w * np.log2(w)).sum())
def rdm(psi, keep):
    n = psi.ndim
    rest = [a for a in range(n) if a not in keep]
    m = np.transpose(psi, list(keep) + rest).reshape(5 ** len(keep), -1)
    return m @ m.T
q5 = {}
for name, n, edges in (("ring of 6 cells", 6, [(i, (i + 1) % 6) for i in range(6)]),
                       ("line of 6 cells", 6, [(i, i + 1) for i in range(5)]),
                       ("ring of 7 cells", 7, [(i, (i + 1) % 7) for i in range(7)])):
    w, psi = ground(n, edges)
    mi = {}
    for i, j in itertools.combinations(range(n), 2):
        if "ring" in name:
            d = min(abs(i - j), n - abs(i - j))
        else:
            d = abs(i - j)
        mi.setdefault(d, []).append(entropy(rdm(psi, [i])) + entropy(rdm(psi, [j])) - entropy(rdm(psi, [i, j])))
    q5[name] = dict(lowest_energies=[round(float(x), 5) for x in w],
                    mutual_information_bits_by_graph_distance={str(d): round(float(np.mean(v)), 5) for d, v in sorted(mi.items())},
                    monotone_decreasing=all(np.mean(mi[d]) > np.mean(mi[d + 1]) for d in sorted(mi)[:-1]))
q5["reading"] = ("Mutual information is largest between neighbours and falls with graph distance on the rings; on the open "
                 "line the two ends are slightly more correlated than distance 4, and the odd ring has a degenerate ground "
                 "state.  What is read back is the distance of the coupling graph that was put in; these sizes say nothing "
                 "about dimension.")
res["Q5"] = q5
json.dump(res, open("out/q3_q4_q5_dynamics_geometry.json", "w"), indent=1)
print(json.dumps(res, indent=1))
