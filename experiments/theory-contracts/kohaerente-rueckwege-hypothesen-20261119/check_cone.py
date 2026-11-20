"""RETURNCONE20261119 (hypotheses 6 and 5): the comparison model of the proposal.
H = [[A, B^+], [B, D]],  A = sin kx sx + sin ky sy + g^2/Delta,  D = Delta,  B = g + v sin kz sz,
Delta = 2, g = v = 1. Claimed: doubly degenerate pole at k = 0, velocities 0.8 sx, 0.8 sy, -0.8 sz,
spatial form of rank three, rank two without the third coupling."""
from common import *
import sympy as sp
import numpy as np

c = Checks("cone")
kx, ky, kz, z = sp.symbols("k_x k_y k_z z", real=True)
s0 = sp.eye(2); sx = sp.Matrix([[0, 1], [1, 0]]); sy = sp.Matrix([[0, -sp.I], [sp.I, 0]]); sz = sp.diag(1, -1)
D_, g_, v_ = 2, 1, 1
def model(Bz):
    A = sp.sin(kx) * sx + sp.sin(ky) * sy + sp.Rational(g_**2, D_) * s0
    B = g_ * s0 + v_ * sp.sin(kz) * Bz
    return sp.Matrix(sp.BlockMatrix([[A, B.H], [B, D_ * s0]]))
H = model(sz)
H0 = H.subs({kx: 0, ky: 0, kz: 0})
ev = H0.eigenvals()
c.ok("k = 0: eigenvalues 0 (twice) and 5/2 (twice)", ev == {0: 2, sp.Rational(5, 2): 2}, ev)
basis = sp.Matrix([[2, 0], [0, 2], [-1, 0], [0, -1]]) / sp.sqrt(5)      # kernel of H0, orthonormal
c.ok("kernel: (2 chi, -chi)/sqrt 5, weight 4/5 on the retained block", H0 * basis == sp.zeros(4, 2) and basis.H * basis == sp.eye(2))
def velocities(Hm):
    return [sp.simplify(basis.H * sp.diff(Hm, q).subs({kx: 0, ky: 0, kz: 0}) * basis) for q in (kx, ky, kz)]
V = velocities(H)
c.ok("velocities 0.8 sx, 0.8 sy, -0.8 sz as claimed", V == [sp.Rational(4, 5) * sx, sp.Rational(4, 5) * sy, -sp.Rational(4, 5) * sz])
G = sp.Matrix(3, 3, lambda i, j: sp.trace(V[i] * V[j]) / 2)
c.ok("spatial form G = 0.64 delta, rank three", G == sp.Rational(16, 25) * sp.eye(3))
V2 = velocities(H.subs(kz, 0))
c.ok("without the third coupling the rank is two", sp.Matrix(3, 3, lambda i, j: sp.trace(V2[i] * V2[j]) / 2).rank() == 2)
anti = lambda Vs: all((Vs[i] * Vs[j] + Vs[j] * Vs[i]) == sp.zeros(2) for i in range(3) for j in range(i))
c.ok("the three velocities anticommute: an isotropic Weyl cone omega^2 = 0.64 |k|^2", anti(V))
handed = sp.simplify((V[0] * V[1] * V[2] / sp.I)[0, 0])
c.ok("its handedness is opposite to the cone sigma.k: V_x V_y V_z = -i 0.512", handed == -sp.Rational(64, 125), handed)
Hn = sp.lambdify((kx, ky, kz), H, "numpy")
kk = np.array([3e-4, -2e-4, 5e-4])
low = sorted(np.linalg.eigvalsh(np.array(Hn(*kk), dtype=complex)), key=abs)[:2]
c.ok("numerical eigenvalues of the full operator confirm the linear cone", abs(abs(low[0]) - 0.8 * np.linalg.norm(kk)) < 1e-6 and abs(low[0] + low[1]) < 1e-6, low)

# --- control 1: the rank diagnostic alone does not certify a cone -------------
Hd = model(s0)                                            # third coupling without a Pauli factor
Vd = velocities(Hd)
Gd = sp.Matrix(3, 3, lambda i, j: sp.trace(Vd[i] * Vd[j]) / 2)
c.ok("control: with B = g + v sin kz (no sigma_z) the form (1/2) tr V_i V_j still has rank three",
     Gd.rank() == 3 and Vd[2] == -sp.Rational(4, 5) * s0)
c.ok("but V_z commutes with V_x, V_y: dispersion omega = -0.8 kz +- 0.8 |k_perp|, a drift, no third cone direction", not anti(Vd))
lowd = sorted(np.linalg.eigvalsh(np.array(sp.lambdify((kx, ky, kz), Hd, "numpy")(*kk), dtype=complex)), key=abs)[:2]
pred = sorted([-0.8 * kk[2] + 0.8 * np.hypot(kk[0], kk[1]), -0.8 * kk[2] - 0.8 * np.hypot(kk[0], kk[1])], key=abs)
c.ok("numerical check of the drift dispersion", np.allclose(sorted(lowd), sorted(pred), atol=1e-6), lowd)
c.ok("consequence: the test must demand the Clifford relation {V_i, V_j} = 2 G_ij, and the sigma_z in B is the neutral two-block put in by hand", True)

# --- control 2: back-reaction reveals, it does not create ---------------------
c.ok("the full 4 x 4 operator depends on k_z before any reduction (dH/dk_z != 0)", sp.diff(H, kz).subs(kz, 0) != sp.zeros(4))
direct = [sp.diff(H[:2, :2], q).subs({kx: 0, ky: 0, kz: 0}) for q in (kx, ky, kz)]
c.ok("the direct block alone has no k_z dependence: rank two", direct[2] == sp.zeros(2))

# --- hypothesis 5: the dynamical embedding given by the return term ------------
# coarse (retained) amplitude vs fine (full) one: iota_z psi = (psi, (z - D)^-1 B psi)
B0 = g_ * s0
iota = lambda zz: sp.Matrix.vstack(s0, B0 / (sp.sympify(zz) - D_))
Zren = lambda zz: sp.simplify((iota(zz).H * iota(zz))[0, 0])
c.ok("embedding weight Z(z) = 1 + g^2/(z - Delta)^2; pole normalisation 1/Z(0) = 4/5 is the factor in the velocities",
     sp.simplify(Zren(z) - 1 - sp.Integer(1) / (z - 2)**2) == 0 and 1 / Zren(0) == sp.Rational(4, 5))
c.ok("the embedding depends on z: (d/dz) iota != 0, so one fixed embedding reproduces the full amplitude at one energy only",
     sp.diff(iota(z), z) != sp.zeros(4, 2))
Sigma = lambda zz: sp.simplify((B0.H * B0 / (sp.sympify(zz) - D_))[0, 0])
c.ok("at k = 0 the exact retained equation z - g^2/Delta - Sigma(z) = 0 has both roots 0 and 5/2: the reduction keeps the full spectrum only with the z-dependent return term",
     sp.solve(sp.Eq(z - sp.Rational(1, 2) - Sigma(z), 0), z) == [0, sp.Rational(5, 2)])
c.done({"V": [str(m.tolist()) for m in V], "G": str(G.tolist()), "drift_control_G_rank": int(Gd.rank())})
