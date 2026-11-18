"""F5 - time, quantum state and entropy arrow: the parts that can be computed without R113.

(a) Restricted readers.  Code state P/5 versus symmetric state S/35 seen through k of the four
    registers: relative entropy and trace distance.  When does a local reader lose the information?
(b) Global unitary control.  A single physical pair pulse exp(-i t A1 A2) on a code state: the
    global entropy stays zero; the leak out of the code is compared with (8/9) sin^2 t.
(c) Reflection positivity.  Two logical cells with exchange K = sum_A R_A x R_A.  Gram form
    G_ij = <Omega| Theta(O_i) O_j |Omega> for H = -K (the documented negative coupling) and
    H = +K (control).  A negative Gram mode rules out an Osterwalder-Schrader reconstruction.
The R113 windows, the actual time map and the Lorentzian part are NOT computed here.
"""
import itertools, json
import numpy as np
import tfpt_ops as T

o = T.build()
V, P, S = o["V"], o["P"], o["S"]
paulis = o["pauli"]
res = {}

def ptrace_keep(rho, k):
    r = rho.reshape([4] * 8)
    for _ in range(4 - k):
        n = r.ndim // 2
        r = np.trace(r, axis1=n - 1, axis2=2 * n - 1)
    return r.reshape(4 ** k, 4 ** k)

def rel_ent(a, b):
    wa, va = np.linalg.eigh(a); wb, vb = np.linalg.eigh(b)
    la = np.where(wa > 1e-12, np.log2(np.clip(wa, 1e-300, None)), 0.0)
    if (wb < 1e-12).any():
        supp = va[:, wa > 1e-12]
        if np.abs(vb[:, wb < 1e-12].T @ supp).max() > 1e-8:
            return float("inf")
    lb = np.where(wb > 1e-12, np.log2(np.clip(wb, 1e-300, None)), 0.0)
    return float((wa * la).sum() - np.trace(a @ (vb * lb) @ vb.T).real)

readers = {}
for k in (1, 2, 3, 4):
    a, b = ptrace_keep(P / 5, k), ptrace_keep(S / 35, k)
    readers[f"{k} registers"] = dict(relative_entropy_bits=round(rel_ent(a, b), 10),
                                     trace_distance=round(float(np.abs(np.linalg.eigvalsh(a - b)).sum() / 2), 10))
res["a_readers_code_vs_symmetric"] = readers
res["a_log2_7"] = round(float(np.log2(7)), 10)

# (b) one physical pair pulse
def on(A, r):
    ops = [np.eye(4)] * 4; ops[r] = A
    return T.kron4(*ops)
A = paulis[1]                                         # IX
h = (on(A, 0) @ on(A, 1)).real
leak = {}
for t in (0.2, 0.5, np.pi / 2):
    U = np.cos(t) * np.eye(256) - 1j * np.sin(t) * h
    L = V.T @ U.conj().T @ (np.eye(256) - P) @ U @ V
    leak[f"t={t:.3f}"] = dict(max_leak=round(float(np.linalg.eigvalsh((L + L.conj().T) / 2).max()), 10),
                              formula_8_9_sin2=round(float(8 / 9 * np.sin(t) ** 2), 10))
psi = V[:, 4]
U = np.cos(0.5) * np.eye(256) - 1j * np.sin(0.5) * h
out = U @ psi
rho = np.outer(out, out.conj())
def vn(r):
    w = np.linalg.eigvalsh(r); w = w[w > 1e-12]
    return float(-(w * np.log2(w)).sum())
res["b_pair_pulse"] = dict(leak=leak, global_entropy_after_pulse=round(vn(rho), 10),
                           one_register_entropy_before=round(vn(ptrace_keep(np.outer(psi, psi), 1)), 6),
                           one_register_entropy_after=round(vn(ptrace_keep(rho, 1)), 6))

# (c) response projectors and reflection positivity of the two-cell exchange
R = []
for A in paulis[1:]:
    D = sum(on(A, r) for r in range(4))
    assert np.abs(V.T @ D @ V).max() < 1e-9
    R.append((V.T @ D @ D @ V).real / 16)
ranks = sorted({int(np.linalg.matrix_rank(r, tol=1e-9)) for r in R})
idem = max(float(np.abs(r @ r - r).max()) for r in R)
K = sum(np.kron(r, r) for r in R)
basis = [np.zeros((5, 5)) for _ in range(25)]
for i in range(25):
    basis[i][i // 5, i % 5] = 1
om5 = np.eye(5).reshape(25) / np.sqrt(5)
rp = {}
for sign, name in ((-1, "H = -K (documented sign)"), (+1, "H = +K (control)")):
    w, v = np.linalg.eigh(sign * K)
    deg = int((w < w[0] + 1e-9).sum())
    Om = v[:, 0]
    G = np.array([[Om @ np.kron(bj, bi.conj()) @ Om for bj in basis] for bi in basis])   # O_j left, Theta(O_i) right
    G = (G + G.T) / 2
    rp[name] = dict(ground_energy=round(float(w[0]), 6), ground_degeneracy=deg, gap=round(float(w[deg] - w[0]), 6),
                    overlap_with_Omega5=round(float(abs(Om @ om5)), 6),
                    min_Gram_eigenvalue=round(float(np.linalg.eigvalsh(G).min()), 8),
                    max_Gram_eigenvalue=round(float(np.linalg.eigvalsh(G).max()), 8))
res["c_two_cell_exchange"] = dict(response_projector_ranks=ranks, idempotence_error=idem,
                                  K_spectrum=sorted({round(float(x), 6) for x in np.linalg.eigvalsh(K)}), **rp)
res["not_computed"] = "R113 weighted windows, the actual time map T(s+t)=T(s)T(t), causal cone and Lorentzian signature"
json.dump(res, open("out/f5_readers_time_positivity.json", "w"), indent=1)
print(json.dumps(res, indent=1))
