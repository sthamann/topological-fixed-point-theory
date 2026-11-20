#!/usr/bin/env python3
"""Step 1 (KERNSIX20261120): does the primitive action of the original produce the Fock state of P_C?

A  collar transfer: the boundary values of solutions of D = gamma_n (d_n + B) form a Lagrangian that
   converges to the positive spectral space of B, whatever the interior
B  second quantisation: the Euclidean fermion integral of the same collar gives a pure state whose
   two-point function is that projector, also for an interior state that is not Gaussian
C  the primitive action S = tr f(D_rel/chi) - tr f(D_ref/chi) + (i pi/2) Delta eta is minus the
   logarithm of the relative fermion determinant exactly for the profile f = -log|x|
D  the determinant line of a polarised block is its Fock vacuum line
"""
import numpy as np
from common import Checks, odd, grading, spectral, expmh, projector, randc, fock, dgamma, slater, twopoint

c = Checks("1_state")
rng = np.random.default_rng(2026112001)

# ---------------------------------------------------------------- A
m = 6
B = odd(randc(rng, m, m))
G = grading(m)
lam1 = np.min(np.abs(np.linalg.eigvalsh(B)))
Ppos = spectral(B, +1)
c.ok("A1 B is self-adjoint and odd under the normal grading G", np.allclose(B, B.conj().T) and np.allclose(G @ B + B @ G, 0))
def transport(L, ell, step=0.25):
    """exp(ell B) L, orthonormalised after every short step (stable for long collars)."""
    T = expmh(B, step)
    q, _ = np.linalg.qr(L)
    for _ in range(int(round(ell / step))):
        q, _ = np.linalg.qr(T @ q)
    return q @ q.conj().T

ells = (8.0, 16.0, 24.0)
interiors = [randc(rng, 2 * m, m) for _ in range(5)]   # arbitrary interior conditions of half dimension
devs = []
for ell in ells:
    worst = 0.0
    for L0 in interiors:
        worst = max(worst, np.linalg.norm(transport(L0, ell) - Ppos, 2))
    devs.append(worst)
rate = [np.log(devs[i] / devs[i + 1]) / (ells[i + 1] - ells[i]) for i in range(2)]   # per unit length
c.ok("A2 boundary values of collar solutions converge to the positive spectral projector of B for every interior condition",
     devs[-1] < 1e-3 * devs[0] and devs[-1] < 1e-3, [f"{d:.2e}" for d in devs])
c.ok("A3 the convergence rate is twice the smallest eigenvalue of |B| (the seed scale of the original)",
     abs(rate[1] - 2 * lam1) < 0.05 * 2 * lam1, (round(rate[1], 4), round(2 * lam1, 4)))
c.ok("A4 the limit obeys the Lagrange relation G P G = 1 - P", np.allclose(G @ Ppos @ G, np.eye(2 * m) - Ppos))

# ---------------------------------------------------------------- B
m2 = 2
n = 2 * m2
B2 = odd(randc(rng, m2, m2))
ann = fock(n)
H2 = dgamma(B2, ann)
P2 = spectral(B2, +1)
L0 = randc(rng, n, m2)
ell = 0.7
v = expmh(H2, ell) @ slater(L0, ann)
v /= np.linalg.norm(v)
C = twopoint(v, ann)
c.ok("B1 Gaussian interior state, finite collar: the two-point function is the projector onto the transported space",
     np.allclose(C, projector(expmh(B2, ell) @ L0), atol=1e-9))
c.ok("B2 that boundary state is pure for every collar length (C^2 = C)", np.allclose(C @ C, C, atol=1e-9))
N = sum(a.conj().T @ a for a in ann)
sector = np.isclose(np.diag(N).real, m2)
worst = 0.0
for _ in range(5):
    w = np.zeros(2 ** n, complex)
    w[sector] = randc(rng, int(sector.sum()))            # generic interior state, not Gaussian
    w = expmh(H2, 40.0 / np.min(np.abs(np.linalg.eigvalsh(B2)))) @ w
    w /= np.linalg.norm(w)
    worst = max(worst, np.linalg.norm(twopoint(w, ann) - P2, 2))
c.ok("B3 long collar: also a non-Gaussian interior state ends in the Fock state of the positive spectral projector",
     worst < 1e-8, f"{worst:.1e}")

# ---------------------------------------------------------------- C
def herm(k):
    a = randc(rng, k, k)
    return (a + a.conj().T) / 2

def eta(D):
    w = np.linalg.eigvalsh(D)
    return int(np.sum(w > 0) - np.sum(w < 0))

def trf(D, chi, f):
    return np.sum(f(np.linalg.eigvalsh(D) / chi))

k = 9
ok_log, ok_chi, bad_other = True, True, 0
for _ in range(20):
    Dref, Drel = herm(k), herm(k)
    ratio = np.linalg.det(Drel) / np.linalg.det(Dref)
    flog = lambda x: -np.log(np.abs(x))
    S = [trf(Drel, chi, flog) - trf(Dref, chi, flog) + 0.5j * np.pi * (eta(Drel) - eta(Dref)) for chi in (0.3, 1.0, 7.0)]
    ok_log &= abs(np.exp(-S[1]) - ratio) < 1e-8 * abs(ratio)
    ok_chi &= abs(S[0] - S[1]) < 1e-9 and abs(S[2] - S[1]) < 1e-9
    fg = lambda x: np.exp(-x ** 2)
    Sg = trf(Drel, 1.0, fg) - trf(Dref, 1.0, fg) + 0.5j * np.pi * (eta(Drel) - eta(Dref))
    bad_other += abs(np.exp(-Sg) - ratio) > 1e-3 * abs(ratio)
c.ok("C1 with f = -log|x| the primitive action is exactly minus the logarithm of det(D_rel)/det(D_ref), phase included",
     ok_log)
c.ok("C2 for that profile the seed scale chi drops out of the action", ok_chi)
c.ok("C3 the phase term (i pi/2) Delta eta alone already equals the phase of the relative fermion determinant",
     all(abs(np.exp(-0.5j * np.pi * (eta(a) - eta(b))) - np.sign(np.linalg.det(a).real * np.linalg.det(b).real)) < 1e-9
         for a, b in [(herm(k), herm(k)) for _ in range(20)]))
c.ok("C4 a Gaussian profile does not give the determinant (the modulus part depends on the declared profile)", bad_other == 20, bad_other)

# ---------------------------------------------------------------- D
E = randc(rng, n, m2)
g = randc(rng, m2, m2)
q, _ = np.linalg.qr(E)
def wedge(cols):
    vv = np.zeros(2 ** n, complex); vv[0] = 1.0
    for kk in range(cols.shape[1]):
        vv = sum(cols[i, kk] * ann[i].conj().T for i in range(n)) @ vv
    return vv
c.ok("D1 the Fock vacuum of a polarised block transforms with the determinant: it spans the determinant line of the block",
     np.allclose(wedge(E @ g), np.linalg.det(g) * wedge(E)))

c.done({"lambda1": lam1, "deviations": devs})
