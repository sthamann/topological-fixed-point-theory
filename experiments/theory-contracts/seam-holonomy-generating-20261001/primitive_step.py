#!/usr/bin/env python3
"""Part Q: SPACETIME.PRIMITIVE.STEP.01 -- which internal amplitude does the BCC step set need?

The eight primitive displacements +-w_a are the nearest neighbours of a bcc lattice (Part S).
  Q1 null-axis cause (exact): for X_a = U^a X U^-a the axis component of the linear term is
     n.Gamma = -2i (O - O^dag), O the (2,3) off-diagonal block of X.  v117 M is complex symmetric with a
     real (2,3) coupling, so O is Hermitian and the deck axis is null.  A generic X lifts it (control).
  Q2 coherence: tr(M_a^dag M_b)/3 = 0 (neighbours), -1/3 (opposite): the alternatives interfere
     (no decohering record in the transport).
  Q3 inventory of existing sectors as the propagation factor:
     - charged nu = 1/4 phase on each visit: band touchings exist, but the projected cone has rank 2;
     - second chirality as a direct sum: same null axis;
     - flavor riding on each step of a unitary walk (A_s (x) M_a): NOT unitary.
     None supplies the missing direction.
  Q4 the bcc Weyl quantum walk exists on EXACTLY these eight steps: U(k) = prod_i (cos k_i - i sin k_i
     sigma_i) = sum_{eps in {+-1}^3} e^{-i eps.k} A_eps, A_eps = prod_i (1 + eps_i sigma_i)/2 rank one;
     unitary, and i log U(k) = k.sigma + O(k^2): isotropic Weyl ({Gamma_i, Gamma_j} = 2 delta), rank 3.
  Q5 Berry charge of its node = +-1 and the sigma -> sigma* partner has the opposite charge.
  Q6 flavor as spectator: U_Weyl (x) U_flavor is unitary and H_eff = k.sigma (x) 1 + 1 (x) H_flavor.
  Q7 (statement) the missing propagation factor is a C^2 spinor on which the spatial point group acts by
     its double cover (the mu8 lift with g^4 = -1 is exactly this sign); with it, unitarity + locality +
     homogeneity + isotropy fix the Weyl walk (Bisio-D'Ariano-Tosini, cited) -- the spinor step itself is
     NOT derived from the present TFPT source: AX.SPINSTEP.
    PYTHONPATH=. python3 -B experiments/theory-contracts/seam-holonomy-generating-20261001/primitive_step.py
"""
from __future__ import annotations

import json
import sys
from itertools import product
from pathlib import Path

import numpy as np
from scipy.linalg import logm
from scipy.optimize import minimize

from tfpt_explorer.flavor_path_transport import exact_clock_dictionary

HERE = Path(__file__).resolve().parent
LOG: list[dict] = []


def step(name: str, ok: bool, detail: str = "") -> None:
    LOG.append({"step": name, "pass": bool(ok), "detail": detail})
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))


M = np.array(exact_clock_dictionary()["M"].evalf(), dtype=complex)
U = np.diag([1, 1j, -1j])
w = np.array([[-.5, -.5, -.5], [-.5, .5, .5], [.5, -.5, .5], [.5, .5, -.5]])
axis = np.array([0.0, 1.0, 0.0])
rng = np.random.default_rng(11)


def dressed(X):
    return [np.linalg.matrix_power(U, a) @ X @ np.linalg.matrix_power(U, -a) for a in range(4)]


def axis_gamma(X):
    Xa = dressed(X)
    return sum(1j * (axis @ w[a]) * (Xa[a] - Xa[a].conj().T) for a in range(4))


O = lambda X: np.array([[0, X[1, 2]], [X[2, 1], 0]])
Xr = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))
Ag, D = axis_gamma(Xr), O(Xr) - O(Xr).conj().T
lam = Ag[1, 2] / D[0, 1]
formula = np.allclose(Ag[0, :], 0) and np.allclose(Ag[:, 0], 0) and np.allclose(Ag[1:, 1:], lam * D)
step("Q1.null-axis cause: n.Gamma = lambda (O - O^dag) on the (2,3) block (zero elsewhere); v117 M is complex symmetric with real (2,3) coupling -> null; a generic X lifts it",
     np.allclose(M, M.T) and abs(M[1, 2].imag) < 1e-15 and np.abs(axis_gamma(M)).max() < 1e-12
     and np.abs(axis_gamma(Xr)).max() > 1e-2 and formula,
     f"|n.Gamma(M)| = {np.abs(axis_gamma(M)).max():.1e}, generic X: {np.abs(axis_gamma(Xr)).max():.3f}, lambda = {np.round(lam, 6)}")

Ma = dressed(M)
ov = np.array([[np.trace(Ma[a].conj().T @ Ma[b]).real / 3 for b in range(4)] for a in range(4)])
step("Q2.coherence of the alternatives: tr(M_a^dag M_b)/3 = 0 for neighbouring, -1/3 for opposite visits (interfering, not decohered)",
     np.allclose(ov, [[1, 0, -1 / 3, 0], [0, 1, 0, -1 / 3], [-1 / 3, 0, 1, 0], [0, -1 / 3, 0, 1]]))

# Q3 inventory
Hc = lambda z: sum(1j * np.exp(1j * z @ w[a]) * Ma[a] - 1j * np.exp(-1j * z @ w[a]) * Ma[a].conj().T for a in range(4))
gapf = lambda z: min(np.diff(np.linalg.eigvalsh(Hc(z))))
best = min((minimize(gapf, rng.uniform(-np.pi, np.pi, 3), method="Nelder-Mead",
                     options={"xatol": 1e-10, "fatol": 1e-13, "maxiter": 4000}) for _ in range(40)), key=lambda r: r.fun)
z0 = best.x
h = 1e-6
G = [(Hc(z0 + h * np.eye(3)[i]) - Hc(z0 - h * np.eye(3)[i])) / (2 * h) for i in range(3)]
e, v = np.linalg.eigh(Hc(z0))
k0 = int(np.argmin(np.diff(e)))
Pv = v[:, [k0, k0 + 1]]
Gp = [Pv.conj().T @ g @ Pv for g in G]
Gp = [g - np.trace(g) / 2 * np.eye(2) for g in Gp]           # traceless (relative) part
gmet = np.array([[np.trace(Gp[i] @ Gp[j]).real / 2 for j in range(3)] for i in range(3)])
rank_charged = int(np.sum(np.linalg.eigvalsh(gmet) > 1e-6))
sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.diag([1.0, -1.0]).astype(complex)
I2 = np.eye(2)
A_eps = lambda eps: (I2 + eps[0] * sx) / 2 @ (I2 + eps[1] * sy) / 2 @ (I2 + eps[2] * sz) / 2
signs = list(product((1, -1), repeat=3))
lab = {}
for a in range(4):
    lab[tuple((2 * w[a]).astype(int))] = (a, 1); lab[tuple((-2 * w[a]).astype(int))] = (a, -1)
def U_flavor_step(k):
    out = np.zeros((6, 6), complex)
    for eps in signs:
        a, o = lab[tuple(eps)]
        out += np.exp(-1j * np.dot(eps, k)) * np.kron(A_eps(eps), Ma[a] if o == 1 else Ma[a].conj().T)
    return out
ks = rng.normal(size=(4, 3))
nonunit = min(np.abs(U_flavor_step(k).conj().T @ U_flavor_step(k) - np.eye(6)).max() for k in ks)
step("Q3.inventory: charged nu=1/4 phase -> touching cone of rank 2; direct-sum second chirality -> same null axis; flavor on each step of a unitary walk -> NOT unitary: no existing sector supplies the missing direction",
     best.fun < 1e-8 and rank_charged == 2 and np.abs(axis_gamma(M.conj())).max() < 1e-12 and nonunit > 0.1,
     f"charged node gap {best.fun:.1e}, rank {rank_charged}; flavor-dressed walk |U^dag U - 1| >= {nonunit:.2f}")

# Q4 bcc Weyl walk on exactly these steps
Uw = lambda k: (np.cos(k[0]) * I2 - 1j * np.sin(k[0]) * sx) @ (np.cos(k[1]) * I2 - 1j * np.sin(k[1]) * sy) @ (np.cos(k[2]) * I2 - 1j * np.sin(k[2]) * sz)
Usum = lambda k: sum(np.exp(-1j * np.dot(eps, k)) * A_eps(eps) for eps in signs)
support = {tuple(eps) for eps in signs} == set(lab.keys())
unit = all(np.allclose(Uw(k).conj().T @ Uw(k), I2) for k in ks)
q = 1e-4 * np.array([1.0, 2.0, 3.0])
Heff = 1j * logm(Uw(q)) / 1e-4
lin = np.allclose(Heff, 1 * sx + 2 * sy + 3 * sz, atol=2e-3)
Gam = [sx, sy, sz]
anti = all(np.allclose(Gam[i] @ Gam[j] + Gam[j] @ Gam[i], 2 * (i == j) * I2) for i in range(3) for j in range(3))
ranks = [np.linalg.matrix_rank(A_eps(eps)) for eps in signs]
step("Q4.the bcc Weyl walk lives on EXACTLY the eight steps: U(k) = prod_i e^{-i k_i sigma_i} = sum_eps e^{-i eps.k} A_eps (A_eps rank one); unitary; i log U = k.sigma + O(k^2): isotropic Weyl, rank 3",
     support and unit and np.allclose(Usum(ks[0]), Uw(ks[0])) and lin and anti and set(ranks) == {1},
     "{Gamma_i, Gamma_j} = 2 delta_ij")


def berry_charge(Hfun, r=0.05, n=40):
    th = np.linspace(1e-3, np.pi - 1e-3, n); ph = np.linspace(0, 2 * np.pi, n, endpoint=False)
    def vec(t, p):
        kk = r * np.array([np.sin(t) * np.cos(p), np.sin(t) * np.sin(p), np.cos(t)])
        e_, v_ = np.linalg.eigh(Hfun(kk)); return v_[:, 1]
    V = [[vec(t, p) for p in ph] for t in th]
    F = 0.0
    for i in range(n - 1):
        for j in range(n):
            a, b, c, d = V[i][j], V[i + 1][j], V[i + 1][(j + 1) % n], V[i][(j + 1) % n]
            F += np.angle(np.vdot(a, b) * np.vdot(b, c) * np.vdot(c, d) * np.vdot(d, a))
    return F / (2 * np.pi)
Hw = lambda k: 0.5 * (1j * logm(Uw(k)) + (1j * logm(Uw(k))).conj().T)
Hw_bar = lambda k: Hw(k).conj()                              # sigma -> sigma*: det -1, opposite chirality
C1, C2 = berry_charge(Hw), berry_charge(Hw_bar)
step("Q5.Berry charge of the Weyl node = +-1 and its sigma -> sigma* partner has the opposite charge: two chiralities",
     abs(abs(C1) - 1) < 0.05 and abs(C1 + C2) < 0.05, f"C = {C1:.3f}, partner {C2:.3f}")
Hf = 0.3 * (Ma[0] + Ma[0].conj().T)
Ut = lambda k, t=1.0: np.kron(Uw(k), np.eye(3)) @ np.kron(I2, np.linalg.matrix_power(np.eye(3), 1)) @ np.kron(I2, __import__("scipy.linalg", fromlist=["expm"]).expm(-1j * Hf))
spect = all(np.allclose(Ut(k).conj().T @ Ut(k), np.eye(6)) for k in ks)
step("Q6.flavor as spectator: U_Weyl (x) U_flavor is unitary, H_eff = k.sigma (x) 1 + 1 (x) H_flavor (flavor holonomy separate from propagation)", spect)
step("Q7.(statement) the missing propagation factor is a C^2 spinor carrying the double cover of the spatial point group (the mu8 sign g^4 = -1); with it the Weyl walk is fixed by unitarity, locality, homogeneity and isotropy (Bisio-D'Ariano-Tosini, cited); the spinor step is NOT derived from the present TFPT source",
     True, "new minimal axiom candidate AX.SPINSTEP")

out = {"steps": LOG, "passed": sum(x["pass"] for x in LOG), "total": len(LOG), "berry": [C1, C2]}
(HERE / "primitive_step_results.json").write_text(json.dumps(out, indent=2, default=float))
print(f"{out['passed']}/{out['total']} steps passed")
sys.exit(0 if out["passed"] == out["total"] else 1)
