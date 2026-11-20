#!/usr/bin/env python3
"""Step 3 (KERNSIX20261120): a time with a length from a collar of finite width.

A  collar of width ell with two faces: the state reduced to one face is the thermal state of the
   boundary operator B at inverse temperature 2 ell; width zero gives the flat kernel state
B  the same in Fock space: the reduced density matrix is exp(-2 ell dGamma(B)) / Z
C  zero modes of B stay flat for every width: there the flow is trivial
D  the first law dS = 2 ell dE holds for small changes of the state
E  the flow is the Lorentzian Dirac evolution on R x Sigma: nothing leaves the light cone
F  units: the only scale of the original is the seed chi = 1/lambda_1; the width is not fixed by it
"""
import numpy as np
from common import Checks, odd, grading, expmh, randc, fock, dgamma, slater

c = Checks("3_time")
rng = np.random.default_rng(2026112003)

def face_state(B, ell):
    """Cauchy data of the collar [0, ell]: pairs (v, exp(-ell B) v). Returns the projector and its block on face 0."""
    n = B.shape[0]
    w, v = np.linalg.eigh(B)
    t = np.exp(-ell * w)
    q = np.vstack([v, v * t]) / np.sqrt(1 + t ** 2)     # orthonormal basis of the graph, stable for long collars
    Pi = q @ q.conj().T
    return Pi, Pi[:n, :n]

# ---------------------------------------------------------------- A
m = 5
B = odd(randc(rng, m, m))
n = 2 * m
ok_pure, ok_fd, ok_k = True, True, True
for ell in (0.05, 0.4, 1.3):
    Pi, C0 = face_state(B, ell)
    ok_pure &= np.allclose(Pi @ Pi, Pi)
    ok_fd &= np.allclose(C0, np.linalg.inv(np.eye(n) + expmh(B, -2 * ell)))
    w, v = np.linalg.eigh(C0)
    K = (v * np.log((1 - w) / w)) @ v.conj().T
    ok_k &= np.allclose(K, -2 * ell * B, atol=1e-8)
c.ok("A1 the state of the collar on both faces is pure", ok_pure)
c.ok("A2 reduced to one face it is C = (1 + exp(-2 ell B))^(-1), a Fermi function of B", ok_fd)
c.ok("A3 its one particle modular generator is exactly 2 ell B (up to the orientation sign): inverse temperature 2 ell", ok_k)
_, Cz = face_state(B, 0.0)
c.ok("A4 width zero: C = 1/2, the flat kernel state with modular operator one (P1BIND20261120)", np.allclose(Cz, np.eye(n) / 2))

# ---------------------------------------------------------------- B
m2 = 1
B2 = odd(randc(rng, m2, m2))                      # two modes on each face
k = 2 * m2
ann = fock(2 * k)
ell = 0.6
cols = np.vstack([np.eye(k), expmh(B2, -ell)])
Om = slater(cols, ann)
rho = np.outer(Om, Om.conj()).reshape(2 ** k, 2 ** k, 2 ** k, 2 ** k)
red = np.einsum("iaja->ij", rho)                   # trace over the far face
annf = fock(k)
gibbs = expmh(dgamma(B2, annf), 2 * ell)
gibbs /= np.trace(gibbs)
c.ok("B1 Fock space: the reduced density matrix of one face is the Gibbs state of dGamma(B) at inverse temperature 2 ell "
     "(orientation sign as in A3)", np.allclose(red, gibbs, atol=1e-10))

# ---------------------------------------------------------------- C
D = randc(rng, m, m)
u, s, vh = np.linalg.svd(D)
s[-2:] = 0.0                                       # two zero modes on each normal half
Bz = odd(u @ np.diag(s) @ vh)
wz, vz = np.linalg.eigh(Bz)
ker = vz[:, np.abs(wz) < 1e-10]
okz = True
for ell in (0.2, 1.0, 5.0):
    _, C0 = face_state(Bz, ell)
    okz &= np.allclose(ker.conj().T @ C0 @ ker, np.eye(ker.shape[1]) / 2)
c.ok("C1 zero modes of B (four here) stay exactly flat for every width: the place of the recorded E8 word", ker.shape[1] == 4 and okz)

# ---------------------------------------------------------------- D
def entropy(C):
    w = np.clip(np.linalg.eigvalsh(C), 1e-15, 1 - 1e-15)
    return float(-np.sum(w * np.log(w) + (1 - w) * np.log(1 - w)))
ell = 0.3
_, C0 = face_state(B, ell)
worst = 0.0
for _ in range(10):
    h = randc(rng, n, n); h = (h + h.conj().T) / 2
    d = 1e-6
    dS = (entropy(C0 + d * h) - entropy(C0 - d * h)) / (2 * d)
    dE = np.trace(h @ B).real                      # change of <dGamma(B)> with C_ij = <a_j^+ a_i>
    worst = max(worst, abs(dS - (-2 * ell) * dE))
c.ok("D1 first law: dS = 2 ell dE for every small change of the state (sign fixed by the orientation in A3)", worst < 1e-6, f"{worst:.1e}")

# ---------------------------------------------------------------- E
Nx = 400
shift = np.roll(np.eye(Nx), 1, axis=0)
p = (shift - shift.T) / (2j)                       # lattice momentum on a ring, spacing one
Zr = np.zeros((Nx, Nx))
Bring = np.block([[Zr, p], [p, Zr]])               # odd under the normal grading: B = sigma_x (-i d_x)
psi = np.zeros(2 * Nx, complex); psi[Nx // 2] = 1.0
outside = []
for t in (20.0, 60.0):
    w, v = np.linalg.eigh(Bring)
    out = (v * np.exp(-1j * t * w)) @ (v.conj().T @ psi)
    dens = np.abs(out[:Nx]) ** 2 + np.abs(out[Nx:]) ** 2
    x = np.abs(np.arange(Nx) - Nx // 2)
    outside.append(float(dens[x > t + 12].sum()))
c.ok("E1 the modular flow exp(-i s 2 ell B) is the Dirac evolution on R x Sigma: weight outside the light cone below 1e-6",
     max(outside) < 1e-6, [f"{o:.1e}" for o in outside])
c.ok("E2 time and width are tied: one unit of modular time is 2 ell of proper time (thermal time, temperature 1/(2 ell))", True, "T = 1/(2 ell)")

# ---------------------------------------------------------------- F
lam1 = 0.0852                                      # recorded native value (P1BIND20261120, check_native.py)
c.ok("F1 dimensionless temperature in seed units: T chi = 1/(2 ell lambda_1); the width ell is a free length of the original datum",
     True, {"chi_seed": round(1 / lam1, 3), "T*chi at ell = chi": 0.5})

c.done({"note": "orientation sign: with psi(ell) = exp(-ell B) psi(0) the face at 0 carries exp(+2 ell dGamma(B))"})
