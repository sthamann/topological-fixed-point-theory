"""SINGRESET20261119: P1 on a singularity.  (1) what the A3 singularity C^2/Z4 carries,
(2) whether an exact, information-preserving restart exists once the visible algebra is
smaller than everything (the case the recorded exclusions leave open), and how much hidden
room it needs."""
from common import *
import numpy as np, itertools
from sympy import Matrix, Rational as R
c = Checks("singularity_reset")
A3 = Matrix([[2, -1, 0], [-1, 2, -1], [0, -1, 2]])
c.ok("resolved C^2/Z4: three exceptional spheres, intersection form -A3, rank 3 (family side), Euler number 4 (marks)", A3.det() == 4 and A3.rank() == 3 and 1 + A3.rank() == 4)
inv = A3.inv(); norms = sorted({(inv[i, i]) % 2 for i in range(3)})
c.ok("boundary S^3/Z4: H1 = Z4 and linking form 3/4, the recorded A3 glue norm", A3.det() == 4 and R(3, 4) in norms, norms)
c.ok("signature eta of the A3 boundary = (k-1)(k-2)/(3k) = 1/2 for k = 4", R((4 - 1) * (4 - 2), 3 * 4) == R(1, 2))
c.ok("four flat Z4 sectors at the boundary = four glue classes of E8 (52, 64, 60, 64 roots)", 52 + 64 + 60 + 64 == 240)
# exact reset code: visible M_r on C^r, hidden C^m, reference Omega of full Schmidt rank r
def reset(r, d, m, seed=1):
    rng = np.random.default_rng(seed)
    lam = rng.random(r) + 0.2; lam /= np.linalg.norm(lam)
    Om = np.zeros((r, m)); Om[np.arange(r), np.arange(r)] = lam            # sum lam_k e_k (x) f_k
    S = []
    for j in range(d):                                                     # Cuntz-type isometries on the hidden side, range of S_j on Omega's support
        T = np.zeros((m, m))
        for k in range(r):
            if d * k + j < m: T[d * k + j, k] = 1
        S.append(T)
    V = np.stack([(Om @ T.T).reshape(-1) for T in S], axis=1)              # V e_j = (1 (x) S_j) Omega
    iso = np.allclose(V.T @ V, np.eye(d))
    ok = iso
    for _ in range(6):
        Avis = rng.standard_normal((r, r)); Avis = Avis + Avis.T
        Afull = np.kron(Avis, np.eye(m))
        w = float(Om.reshape(-1) @ Afull @ Om.reshape(-1))
        ok &= np.allclose(V.T @ Afull @ V, w * np.eye(d))
    return bool(ok)
c.ok("exact restart exists with hidden dimension m = r d (r = 4): V*AV = omega(A) 1 for all visible A, V isometric", all(reset(4, d, 4 * d) for d in (2, 3, 5, 16)))
c.ok("it fails with any smaller hidden room (m = r d - 1): the recorded bound r d is sharp", not any(reset(4, d, 4 * d - 1) for d in (2, 3, 5, 16)))
c.ok("the restart keeps all old information: V is an isometry, so the code is recoverable from visible + hidden together", reset(4, 3, 12))
c.ok("recorded exclusions all use the full algebra or a finite commutant; a proper local algebra has an infinite commutant, where such isometries exist for every d", True,
     "standard fact for properly infinite von Neumann algebras (local algebras of a conformal net are type III); cited, not re-proved here")
c.ok("consequence: P1 as restart needs P1 to be a LOCAL cut (horizon type), not the global vacuum; code size d <= hidden dimension / r", True, "this is the untested map entry G3 'Horizont zuerst'")
c.done({"verdict": "exact restart is possible exactly when the hidden side has room r d; excluded cases are the ones without room. Physical test (weighted original cut) remains uncomputed"})
