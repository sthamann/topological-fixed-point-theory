"""Many-body level. The Fock state of the original polarization on the two normal halves, and its mirror.

d modes per normal half. Modes 1..d span the half G = +1, modes d+1..2d the half G = -1 (Jordan-Wigner).
A_+ is the algebra of the first half; it is the tensor factor M_(2^d) x 1.
The polarization is P_C = (1/2) [[1, U*], [U, 1]] with a unitary U, as forced by the Lagrangian property.
"""
import itertools
import numpy as np
from common import Checks, randc

c = Checks("fock")
rng = np.random.default_rng(20261120)
norm = lambda A: float(np.abs(A).max())
I2, Z2 = np.eye(2), np.diag([1.0, -1.0])
low = np.array([[0, 1], [0, 0]], complex)            # annihilates |1>

def kron(*ops):
    out = np.eye(1)
    for o in ops: out = np.kron(out, o)
    return out

def run(d, seed):
    rng = np.random.default_rng(seed)
    n = 2 * d
    a = [kron(*([Z2] * k + [low] + [I2] * (n - k - 1))) for k in range(n)]
    ad = [x.conj().T for x in a]
    N = sum(ad[k] @ a[k] for k in range(n))
    par = np.diag(np.exp(1j * np.pi * np.diag(N).real)).real
    U = np.linalg.qr(randc(rng, d, d))[0]
    vac = np.zeros(2 ** n, complex); vac[0] = 1
    omega = vac
    for k in range(d):
        bk = (ad[k] + sum(U[j, k] * ad[d + j] for j in range(d))) / np.sqrt(2)
        omega = bk @ omega
    omega /= np.linalg.norm(omega)
    ev = lambda A: np.vdot(omega, A @ omega)
    out = {}
    # two-point function equals the polarization
    PC = 0.5 * np.block([[np.eye(d), U.conj().T], [U, np.eye(d)]])
    two = np.array([[ev(ad[q] @ a[p]) for q in range(n)] for p in range(n)])
    out["two_point"] = max(norm(two - PC), norm(two - PC.T), norm(two - PC.conj()))
    out["two_point_exact"] = min(norm(two - PC), norm(two - PC.T), norm(two - PC.conj()))
    # reduced state of the first half
    psi = omega.reshape(2 ** d, 2 ** d)
    rho = psi @ psi.conj().T
    out["flat"] = norm(rho - np.eye(2 ** d) / 2 ** d)
    # the second half, tested on all monomials
    worst = 0.0
    for bits in itertools.product((0, 1, 2, 3), repeat=d):
        op = np.eye(2 ** n, dtype=complex)
        tr = 1.0
        for k, b in enumerate(bits):
            m = [np.eye(2 ** n), a[d + k], ad[d + k], ad[d + k] @ a[d + k]][b]
            op = op @ m
            tr *= [1.0, 0.0, 0.0, 0.5][b]
        worst = max(worst, abs(ev(op) - tr))
    out["flat_second"] = worst
    # modular conjugation of (A_+, omega): J (X omega) = X^dagger omega
    units = []
    for i in range(2 ** d):
        for j in range(2 ** d):
            E = np.zeros((2 ** d, 2 ** d), complex); E[i, j] = 1
            units.append(np.kron(E, np.eye(2 ** d)))
    V = np.array([X @ omega for X in units]).T
    W = np.array([X.conj().T @ omega for X in units]).T
    out["cyclic_rank"] = int(np.linalg.matrix_rank(V))
    Jm = W @ np.linalg.inv(V.conj())                  # J = Jm o conj
    J = lambda x: Jm @ x.conj()
    out["J_unitary"] = norm(Jm.conj().T @ Jm - np.eye(2 ** n))
    out["J_square"] = norm(Jm @ Jm.conj() - np.eye(2 ** n))
    out["J_fixes_omega"] = float(np.linalg.norm(J(omega) - omega))
    jop = lambda F: Jm @ F.conj() @ Jm.conj().T      # J F J as a matrix (J^-1 = J)
    # reflected form against the trace
    worst = 0.0
    sample = [units[i] for i in rng.choice(len(units), size=min(len(units), 40), replace=False)]
    for F in sample:
        for Gm in sample:
            worst = max(worst, abs(ev(jop(F) @ Gm) - np.trace(F[::2 ** d, ::2 ** d].conj().T @ Gm[::2 ** d, ::2 ** d]) / 2 ** d))
    out["form_is_trace"] = worst
    out["reference"] = abs(ev(jop(np.eye(2 ** n)) @ np.eye(2 ** n)) - 1)
    # image of the first half lies in its commutant
    out["commutant"] = max(norm(jop(a[k]) @ X - X @ jop(a[k])) for k in range(d) for X in units[:: max(1, len(units) // 16)])
    # one-particle content of the mirror: J a(e_k) J against Klein factor times fields of the second half
    cand = {}
    kleins = {"1": np.eye(2 ** n), "parity": par}
    for kn, K in kleins.items():
        for fn, fields in (("a", a), ("a*", ad)):
            basis = [(K @ fields[d + j]).reshape(-1) for j in range(d)]
            A = np.array(basis).T
            coef, res = [], 0.0
            for k in range(d):
                t = jop(a[k]).reshape(-1)
                x = np.linalg.lstsq(A, t, rcond=None)[0]
                res = max(res, float(np.linalg.norm(A @ x - t)))
                coef.append(x)
            cand[f"{kn}.{fn}"] = (res, np.array(coef).T)
    best = min(cand, key=lambda k: cand[k][0])
    out["mirror_form"] = best
    out["mirror_residual"] = cand[best][0]
    L = cand[best][1]                                   # J a(e_k) J = Klein * sum_j L[j,k] field(f_j)
    out["mirror_matrix_unitary"] = norm(L.conj().T @ L - np.eye(d))
    out["mirror_matrix_vs_U"] = {"U": norm(np.abs(L) - np.abs(U)), "phase_U": min(norm(L - z * U) for z in (1, -1, 1j, -1j)),
                                 "phase_conjU": min(norm(L - z * U.conj()) for z in (1, -1, 1j, -1j))}
    out["other_residuals"] = {k: round(v[0], 6) for k, v in cand.items()}
    # the alternative mirror without the adjoint: second quantization of the anti-unitary  T = i_C o conj
    T = np.block([[np.zeros((d, d)), U.conj().T], [U, np.zeros((d, d))]])
    herm = []
    for _ in range(60):
        X = randc(rng, 2 ** d, 2 ** d); herm.append(X + X.conj().T)
    # theta(a(e_k)) = a(T conj e_k) = sum_j U[j,k] a_{d+j}; extended anti-linearly and multiplicatively on matrix units
    def theta(Fsmall):
        # express F in the monomial basis of the first half, map generators, conjugate coefficients
        res = np.zeros((2 ** n, 2 ** n), complex)
        for bits, mono_first, mono_second in MONO:
            coeff = np.trace(DUAL[bits].conj().T @ Fsmall)
            res += np.conj(coeff) * mono_second
        return res
    img = [sum(U[j, k] * a[d + j] for j in range(d)) for k in range(d)]
    imgd = [x.conj().T for x in img]
    MONO, DUAL = [], {}
    smalla = [kron(*([Z2] * k + [low] + [I2] * (d - k - 1))) for k in range(d)]
    for bits in itertools.product((0, 1, 2, 3), repeat=d):
        f, s = np.eye(2 ** d, dtype=complex), np.eye(2 ** n, dtype=complex)
        for k, b in enumerate(bits):
            f = f @ [np.eye(2 ** d), smalla[k], smalla[k].conj().T, smalla[k].conj().T @ smalla[k]][b]
            s = s @ [np.eye(2 ** n), img[k], imgd[k], imgd[k] @ img[k]][b]
        MONO.append((bits, f, s))
    gram = np.array([[np.trace(x[1].conj().T @ y[1]) for y in MONO] for x in MONO])
    ginv = np.linalg.inv(gram)
    for i, (bits, f, s) in enumerate(MONO):
        DUAL[bits] = sum(ginv[j, i].conj() * MONO[j][1] for j in range(len(MONO)))
    vals = [ev(theta(X) @ np.kron(X, np.eye(2 ** d))).real for X in herm]
    out["alt_mirror_values"] = (min(vals), max(vals))
    return out

res = {d: run(d, 20261120 + d) for d in (1, 2, 3)}
for d, o in res.items():
    c.ok(f"d = {d}: the Fock state of P_C has two-point function P_C", o["two_point_exact"] < 1e-9, f"{o['two_point_exact']:.1e}")
    c.ok(f"d = {d}: its restriction to the first normal half is the flat density 1/2^d", o["flat"] < 1e-12)
    c.ok(f"d = {d}: its restriction to the second normal half is the trace on all {4**d} monomials", o["flat_second"] < 1e-12)
    c.ok(f"d = {d}: omega is cyclic for A_+ (rank {4**d}); J(X omega) = X* omega is anti-unitary, J^2 = 1, J omega = omega",
         o["cyclic_rank"] == 4 ** d and o["J_unitary"] < 1e-9 and o["J_square"] < 1e-9 and o["J_fixes_omega"] < 1e-9)
    c.ok(f"d = {d}: reflected form <omega, (J F J) G omega> = Tr(F* G)/2^d, reference 1", o["form_is_trace"] < 1e-9 and o["reference"] < 1e-12)
    c.ok(f"d = {d}: J A_+ J commutes with A_+", o["commutant"] < 1e-9)
    c.ok(f"d = {d}: J a(e_k) J = parity x a*(L e_k) on the second half, L unitary and equal to U up to one phase",
         o["mirror_form"] == "parity.a*" and o["mirror_residual"] < 1e-9 and o["mirror_matrix_unitary"] < 1e-9
         and min(o["mirror_matrix_vs_U"]["phase_U"], o["mirror_matrix_vs_U"]["phase_conjU"]) < 1e-9,
         f"{o['mirror_form']}, residuals {o['other_residuals']}, match {o['mirror_matrix_vs_U']}")
for d in (2, 3):
    lo, hi = res[d]["alt_mirror_values"]
    c.ok(f"d = {d}: the side exchange without the adjoint (a -> a) is not reflection positive: both signs occur", lo < -1e-6 and hi > 1e-6,
         f"values from {lo:.4f} to {hi:.4f}")
c.done({str(d): {k: (v if not isinstance(v, np.ndarray) else v.tolist()) for k, v in o.items()} for d, o in res.items()})
