#!/usr/bin/env python3
"""KERNORIG20261120 B: the reconstruction of the original applied to its own collar operator.

The original defines the form (F, G) = <Theta F, G> on observables at positive Euclidean time,
divides by its null space and reads the time semigroup exp(-t H), H >= 0.  Here for the free field
of D = gamma_n (d_tau + B) on the line times the seam (one particle level and a small Fock space).
B1  Green function of d_tau + B
B2  reflected kernel on positive times: exp(-(t + t') |B|) on one spectral half, positive
B3  the null space of the form is the other half: the reconstructed space is the range of the
    spectral projector, i.e. the polarisation
B4  the time semigroup is exp(-t |B|): generator |B| >= lambda_1 > 0 (a mass gap in seed units)
B5  equal time limit of the two point function = the projector (the Fock state of step 1)
B6  Fock space: that state is the ground state of the normal ordered dGamma(B), gap lambda_1
B7  the face state of a collar of width l (KERNSIX step 3) tends to this vacuum for l -> infinity
B8  Planck coefficient of the original, M^2 = chi^2/(2 pi^2) = 32 c3^2 chi^2, against the area law
    of KERNSIX step 6: both agree only for l chi = (3/4) sqrt(zeta(3)); this is a condition, not a result
"""
import numpy as np
from scipy.special import zeta
from common import Checks, odd, grading, spectral, expmh, randc, fock, dgamma, slater, twopoint

c = Checks("b_os")
rng = np.random.default_rng(2026112012)
m = 4
B = odd(randc(rng, m, m))
n = 2 * m
G = grading(m)
w, v = np.linalg.eigh(B)
Pp, Pm = spectral(B, +1), spectral(B, -1)
absB = (v * np.abs(w)) @ v.conj().T
lam1 = np.min(np.abs(w))

def S(tau):
    """Green function of d_tau + B on the line, decaying in both directions."""
    return expmh(B, -tau) @ Pp if tau > 0 else -expmh(B, -tau) @ Pm

h = 1e-6
jump = S(h) - S(-h)
dS = (S(0.5 + h) - S(0.5 - h)) / (2 * h)
c.ok("B1 (d_tau + B) S = 0 away from zero and S jumps by one at zero", np.allclose(jump, np.eye(n), atol=1e-5) and np.allclose(dS + B @ S(0.5), 0, atol=1e-5))

times = [0.1, 0.35, 0.8, 1.7]
K = np.block([[-S(-(a + b)) for b in times] for a in times])          # reflected kernel on positive times
ev = np.linalg.eigvalsh((K + K.conj().T) / 2)
ok_form = np.allclose(K, K.conj().T) and ev.min() > -1e-10
ok_shape = all(np.allclose(-S(-(a + b)), expmh(absB, -(a + b)) @ Pm) for a in times for b in times)
c.ok("B2 the reflected kernel is exp(-(t + t') |B|) on one spectral half and positive semidefinite", ok_form and ok_shape, f"min eigenvalue {ev.min():.1e}")
rank = int(np.sum(ev > 1e-9))
c.ok("B3 its rank is half the boundary dimension: the null space is the other spectral half, the quotient is the polarisation",
     rank == m and np.allclose(K @ np.kron(np.eye(len(times)), Pp), 0, atol=1e-10), rank)
c.ok("B4 time translation acts as exp(-t |B|) with generator |B| >= lambda_1 > 0", np.allclose(-S(-(0.3 + 0.5)), expmh(absB, -0.5) @ (-S(-0.3))) and lam1 > 0, round(float(lam1), 4))
c.ok("B5 equal time limits of the two point function are the two spectral projectors; the normal grading exchanges them",
     np.allclose(S(1e-9), Pp, atol=1e-7) and np.allclose(-S(-1e-9), Pm, atol=1e-7) and np.allclose(G @ Pp @ G, Pm))

m2 = 2
B2 = odd(randc(rng, m2, m2))
ann = fock(2 * m2)
H = dgamma(B2, ann)
ew, evec = np.linalg.eigh(H)
gs = evec[:, 0]
C = twopoint(gs, ann)
c.ok("B6 Fock space: the ground state of dGamma(B) has the spectral projector as two point function; the gap is lambda_1",
     np.allclose(C, spectral(B2, -1), atol=1e-9) and abs((ew[1] - ew[0]) - np.min(np.abs(np.linalg.eigvalsh(B2)))) < 1e-9)

devs = []
for ell in (1.0, 4.0, 16.0):
    Cf = (v / (1 + np.exp(-2 * ell * w))) @ v.conj().T
    devs.append(float(np.linalg.norm(Cf - Pp, 2)))
c.ok("B7 the face state (1 + exp(-2 l B))^(-1) of a collar tends to this vacuum as l grows: the original sits at temperature zero",
     devs[2] < 1e-3 and devs[0] > devs[1] > devs[2], [f"{d:.1e}" for d in devs])

c3 = 1 / (8 * np.pi)
lchi = np.sqrt(18 * zeta(3) / 32)
c.ok("B8 1/(2 pi^2) = 32 c3^2; equal to the area law coefficient 18 zeta(3) c3^2 / l^2 only if l chi = (3/4) sqrt(zeta(3)) = 0.8223",
     abs(1 / (2 * np.pi ** 2) - 32 * c3 ** 2) < 1e-15 and abs(lchi - 0.75 * np.sqrt(zeta(3))) < 1e-12 and abs(lchi - 0.8223) < 1e-4, round(float(lchi), 5))
c.done({"lambda1": float(lam1)})
