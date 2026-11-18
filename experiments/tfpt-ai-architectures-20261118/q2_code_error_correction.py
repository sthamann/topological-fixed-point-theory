"""Q2 (and KI4a) - is the quartic five-space a quantum error-correcting code, and is that special?

Knill-Laflamme test: an error set is correctable iff V^+ E_a^+ E_b V is a multiple of the identity.
  one known lost register   <=> every operator on ONE register is trivial on the code
  two known lost registers  <=> every operator on TWO registers is trivial on the code
  one unknown Pauli error   detectable <=> weight-1 trivial, correctable <=> weight-2 trivial
Compared codes: the TFPT quartic code, random 5-dim subspaces of the same symmetric space
Sym^4 (rank 35), random 5-dim subspaces of the full space, and for scale the standard
[[4,2,2]] qubit code is quoted in the text only.
Then an explicit recovery after losing a register, with the decoder of the source document.
"""
import itertools, json
import numpy as np
import tfpt_ops as T

rng = np.random.default_rng(1)
o = T.build()
V, S, pauli = o["V"], o["S"], o["pauli"]
res = {}

def on(ops):
    full = [np.eye(4)] * 4
    for r, A in ops.items():
        full[r] = A
    return T.kron4(*full)

W1 = [on({r: A}) for r in range(4) for A in pauli[1:]]
W2 = [on({r: A, s: B}) for r, s in itertools.combinations(range(4), 2) for A in pauli[1:] for B in pauli[1:]]

def kl_violation(Vc, ops):
    worst = 0.0
    for E in ops:
        M = Vc.conj().T @ E @ Vc
        worst = max(worst, float(np.abs(M - np.trace(M) / 5 * np.eye(5)).max()))
    return worst

def readout_ranks(Vc):
    out = []
    for k in (1, 2, 3):
        rows = []
        for i in range(25):
            X = np.zeros((5, 5)); X[i // 5, i % 5] = 1
            rho = (Vc @ X @ Vc.conj().T).reshape([4] * 8)
            for _ in range(4 - k):
                n = rho.ndim // 2
                rho = np.trace(rho, axis1=n - 1, axis2=2 * n - 1)
            rows.append(rho.reshape(-1))
        out.append(int(np.linalg.matrix_rank(np.array(rows), tol=1e-9)))
    return out

def recovery_fidelity(Vc, trials=200):
    """lose register 0, apply R(sigma) = sum_a W_a^+ sigma W_a with W_a = 2 V_a; mean and worst fidelity."""
    Va = Vc.reshape(4, 64, 5)
    f = []
    for _ in range(trials):
        psi = rng.standard_normal(5) + 1j * rng.standard_normal(5); psi /= np.linalg.norm(psi)
        rho3 = sum(Va[a] @ np.outer(psi, psi.conj()) @ Va[a].conj().T for a in range(4))
        rec = sum((2 * Va[a]).conj().T @ rho3 @ (2 * Va[a]) for a in range(4))
        rec = rec / np.trace(rec)
        f.append(float((psi.conj() @ rec @ psi).real))
    return round(float(np.mean(f)), 6), round(float(np.min(f)), 6)

def random_code(space):
    if space == "sym":
        w, U = np.linalg.eigh(S)
        Bs = U[:, w > 0.5]                                    # 256 x 35
        Q = np.linalg.qr(rng.standard_normal((35, 5)))[0]
        return Bs @ Q
    return np.linalg.qr(rng.standard_normal((256, 5)))[0]

tf = dict(kl_violation_weight1=kl_violation(V, W1), kl_violation_weight2=round(kl_violation(V, W2), 6),
          readout_ranks_1_2_3_registers=readout_ranks(V), recovery_after_one_lost_register=recovery_fidelity(V))
tf["distance"] = 2 if tf["kl_violation_weight1"] < 1e-9 < tf["kl_violation_weight2"] else None
res["TFPT quartic code"] = tf
for space in ("sym", "full"):
    k1, k2, rec = [], [], []
    for _ in range(20):
        Vc = random_code(space)
        k1.append(kl_violation(Vc, W1)); k2.append(kl_violation(Vc, W2)); rec.append(recovery_fidelity(Vc, 50)[0])
    res[f"random rank-5 code in {'Sym^4' if space == 'sym' else 'full space'} (20 samples)"] = dict(
        kl_violation_weight1=[round(float(np.min(k1)), 4), round(float(np.mean(k1)), 4)],
        kl_violation_weight2=[round(float(np.min(k2)), 4), round(float(np.mean(k2)), 4)],
        recovery_after_one_lost_register_mean=round(float(np.mean(rec)), 4))

# detection in a noisy memory: independent depolarising noise p on every register, then project on the code
def detect(Vc, p, trials=3000):
    acc = kept = 0.0
    for _ in range(trials):
        psi = rng.standard_normal(5) + 1j * rng.standard_normal(5); psi /= np.linalg.norm(psi)
        st = Vc @ psi
        ops = {r: pauli[rng.integers(1, 16)] for r in range(4) if rng.random() < p}
        if ops:
            st = on(ops) @ st
        back = Vc.conj().T @ st
        w = float(np.linalg.norm(back) ** 2)
        kept += w
        if w > 1e-12:
            acc += w * abs(psi.conj() @ back) ** 2 / w
    return round(kept / trials, 4), round(acc / kept, 5)
res["noisy_memory_with_postselection"] = {
    f"p={p}": dict(zip(("accepted", "fidelity_when_accepted"), detect(V, p)),
                   unprotected_fidelity_4_registers=round((1 - p) ** 4, 5)) for p in (0.01, 0.03, 0.1)}
res["reading"] = ("The quartic code detects any single-register error and corrects one known lost register exactly "
                  "(distance 2), as the source document states; two lost registers or one error at an unknown "
                  "place are not correctable.  Random rank-5 subspaces, also inside Sym^4, do not have this property, "
                  "so within this comparison the code is special.  In absolute terms it is a small detection code: "
                  "8 physical qubits for log2(5) = 2.32 logical qubits at distance 2, the same distance as the "
                  "textbook [[4,2,2]] code, which stores 2 qubits in 4.")
json.dump(res, open("out/q2_code_error_correction.json", "w"), indent=1)
print(json.dumps(res, indent=1))
