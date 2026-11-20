"""One-particle level. What the original boundary datum fixes about state, marker and cut.

Source facts used (read, not computed), boundary kernel source section 001:
  * the datum contains a self-adjoint tangential operator B_Sigma and the collar form D_+ = gamma_n (d_n + B_Sigma);
  * the polarization is the Calderon projector P_C of the interior problem, i_C = 2 P_C - 1;
  * the seed response is 1 / lambda_1^+(|B_Sigma|).
Convention: gamma_n = i G with G = sigma_z x 1 the normal grading. Interior n <= 0, boundary at n = 0.
"""
import numpy as np
from common import Checks, TOL, odd, grading, spectral, expmh, projector, randc

c = Checks("polarization")
rng = np.random.default_rng(20261120)
m = 6
G = grading(m)
gam = 1j * G
norm = lambda A: float(np.abs(A).max())

# 1 symmetry of the collar operator needs an odd B
D0 = randc(rng, m, m)
B = odd(D0)
c.ok("B odd under the normal grading: gamma_n B + B gamma_n = 0", norm(gam @ B + B @ gam) < TOL)
Bbad = B + G @ np.diag(rng.normal(size=2 * m)).astype(complex)
u0, v0 = randc(rng, 2 * m), randc(rng, 2 * m)
flux = lambda Bx: np.vdot(gam @ (-Bx @ u0), v0) + np.vdot(gam @ u0, -Bx @ v0)
c.ok("the boundary form <gamma_n u, v> is conserved along the normal for odd B and not for B with an even part",
     abs(flux(B)) < TOL and abs(flux(Bbad)) > 1e-3, f"odd {abs(flux(B)):.1e}, with even part {abs(flux(Bbad)):.3f}")

# 2 layered, non-product interior: B(n) piecewise constant, all layers odd; decaying solutions as n -> -infinity
layers = [odd(randc(rng, m, m)) for _ in range(40)]
h = 0.05
start = spectral(layers[0], -1)               # d_n u = -B u decays towards n -> -infinity on the negative part
w, v = np.linalg.eigh(start)
cols = v[:, w > 0.5]
for Bk in layers:
    cols = expmh(Bk, -h) @ cols
    cols, _ = np.linalg.qr(cols)
PC = projector(cols)
iC = 2 * PC - np.eye(2 * m)
c.ok("Cauchy data of the layered interior: a projector of half rank", norm(PC @ PC - PC) < TOL and round(np.trace(PC).real) == m)
c.ok("it is not the spectral projector of the boundary operator (the interior is not a product)",
     norm(PC - spectral(layers[-1], -1)) > 1e-2, f"distance {norm(PC - spectral(layers[-1], -1)):.3f}")
c.ok("Lagrangian: gamma_n maps the Cauchy data space onto its orthogonal complement, G P_C G = 1 - P_C",
     norm(G @ PC @ G - (np.eye(2 * m) - PC)) < 1e-7, f"defect {norm(G @ PC @ G - (np.eye(2 * m) - PC)):.1e}")
c.ok("hence polarization and normal grading anticommute: i_C G + G i_C = 0", norm(iC @ G + G @ iC) < 1e-7)

# 3 the state on the two normal halves
a, d, b = PC[:m, :m], PC[m:, m:], PC[m:, :m]
U = 2 * b
c.ok("restricted to either normal half the two-point function is exactly 1/2: the flat (tracial) state",
     norm(a - np.eye(m) / 2) < 1e-7 and norm(d - np.eye(m) / 2) < 1e-7)
c.ok("the correlation between the halves is a unitary U / 2: the state is pure and maximally correlated",
     norm(U.conj().T @ U - np.eye(m)) < 1e-6 and norm(iC[m:, :m] - U) < 1e-7)
flat = lambda P: max(norm(P[:m, :m] - np.eye(m) / 2), norm(P[m:, m:] - np.eye(m) / 2))
others = [flat(projector(randc(rng, 2 * m, m))) for _ in range(200)]
c.ok("200 random polarizations of half rank: none is flat on the halves", min(others) > 1e-2, f"smallest defect {min(others):.3f}")
conv = []
for _ in range(50):
    Ur = np.linalg.qr(randc(rng, m, m))[0]
    P = 0.5 * np.block([[np.eye(m), Ur.conj().T], [Ur, np.eye(m)]])
    conv.append(max(norm(P @ P - P), norm(G @ P @ G - (np.eye(2 * m) - P))))
c.ok("converse: every projector that is flat on both halves is Lagrangian (50 random unitaries)", max(conv) < TOL)
Beven = B + 0.7 * G @ odd(randc(rng, m, m)) @ G @ G + 0.7 * np.diag(rng.normal(size=2 * m))
c.ok("control: the positive half of an operator with an even part is not flat", flat(spectral(Beven, 1)) > 1e-2,
     f"defect {flat(spectral(Beven, 1)):.3f}")

# 4 the two sides of the double
PCop = projector(gam @ cols)
c.ok("Cauchy data of the reflected side M_+^op is gamma_n (Ran P_C) = ker P_C", norm(PCop - (np.eye(2 * m) - PC)) < 1e-7)
c.ok("the pointwise lift of the reflection, G, exchanges the two Cauchy data spaces and is not i_C",
     norm(G @ PC @ G - PCop) < 1e-7 and norm(G - iC) > 0.5 and norm(G + iC) > 0.5)
c.ok("an eigenprojector of the pointwise lift G is fixed by G, so it is never a Cauchy data space",
     norm(G @ ((np.eye(2 * m) + G) / 2) @ G - (np.eye(2 * m) + G) / 2) < TOL)

# 5 marker and cut: canonical marks on every line of the polarization
Y0 = 1j * G @ iC
c.ok("Y0 = i G i_C is a self-adjoint involution anticommuting with G and with i_C",
     norm(Y0 - Y0.conj().T) < 1e-7 and norm(Y0 @ Y0 - np.eye(2 * m)) < 1e-6 and norm(Y0 @ G + G @ Y0) < 1e-7 and norm(Y0 @ iC + iC @ Y0) < 1e-6)
worst = 0.0
for eta, P in ((1, PC), (-1, np.eye(2 * m) - PC)):
    for _ in range(40):
        psi = P @ randc(rng, 2 * m); psi /= np.linalg.norm(psi)
        ev = lambda A: np.vdot(psi, A @ psi)
        th = rng.uniform(0, 2 * np.pi)
        Um, Ym = np.cos(th) * G + np.sin(th) * Y0, -np.sin(th) * G + np.cos(th) * Y0
        worst = max(worst, abs(ev(G)), abs(ev(Y0)), abs(ev(G @ Y0) - 1j * eta), abs(ev(Um)), abs(ev(Ym)),
                    abs(ev(Um @ Ym) - 1j * eta), norm(Um @ Ym @ Um @ Ym + np.eye(2 * m)))
c.ok("on every line of the side eta: Z_0 = 1, Z_U = Z_Y = 0, Z_UY = i eta, ABAB = -1; unchanged under rotating (U, Y) about i_C",
     worst < 1e-6, f"largest deviation {worst:.1e}")
c.ok("the mixed mark is the polarization itself: G Y0 = i i_C", norm(G @ Y0 - 1j * iC) < 1e-6)
ok = 0.0
for _ in range(30):
    cm = np.block([[randc(rng, m, m), np.zeros((m, m))], [np.zeros((m, m)), randc(rng, m, m)]])
    cm = cm + cm.conj().T
    for eta, P in ((1, PC), (-1, np.eye(2 * m) - PC)):
        ok = max(ok, abs(np.trace(P @ cm).real / m - np.trace(cm).real / (2 * m)),
                 abs(np.trace(P @ cm @ Y0) / m - 1j * eta * np.trace(G @ cm) / (2 * m)), abs(np.trace(P @ Y0)))
c.ok("flat state of one side, marker c commuting with G: Z_c = <c>, Z_Y = 0, Z_cY = i eta <G c> (chiral average of the marker)",
     ok < 1e-6, f"largest deviation {ok:.1e}")

# 6 product collar: polarization = spectral half, cut = phase of the tangential part, seed line
Bs = odd(D0)
Ps = spectral(Bs, -1)
tang = gam @ Bs
wt, vt = np.linalg.eigh(tang)
phase = (vt * np.sign(wt)) @ vt.conj().T
c.ok("product collar: the tangential part gamma_n B of D_+ is self-adjoint and Y0 is minus its phase",
     norm(tang - tang.conj().T) < TOL and norm(1j * G @ (2 * Ps - np.eye(2 * m)) + phase) < 1e-7)
w, v = np.linalg.eigh(Bs)
lam1 = min(abs(w))
seed = v[:, np.argmin(np.abs(w + lam1))]
c.ok("the seed line (eigenvalue of smallest modulus on the polarized side) lies in Ran P_C and carries the canonical marks",
     abs(np.vdot(seed, Ps @ seed) - 1) < 1e-8 and abs(np.vdot(seed, G @ (1j * G @ (2 * Ps - np.eye(2 * m))) @ seed) - 1j) < 1e-6)

# 7 round S^2: spectrum +-(k+1) with multiplicity 2(k+1)
ks = range(0, 6)
blocks = []
for k in ks:
    for _ in range(2 * (k + 1)):
        blocks.append(k + 1)
Dd = np.diag(blocks).astype(complex)
Bsph = odd(Dd @ np.linalg.qr(randc(rng, len(blocks), len(blocks)))[0])
ws = np.sort(np.linalg.eigvalsh(Bsph))
mult = {k + 1: int(np.sum(np.abs(ws - (k + 1)) < 1e-8)) for k in ks}
c.ok("S^2 model: eigenvalues +-(k+1), each sign with multiplicity 2(k+1); lambda_1^+ = 1",
     all(mult[k + 1] == 2 * (k + 1) for k in ks) and abs(min(abs(ws)) - 1) < 1e-9, mult)
M = len(blocks)
m = M
c.ok("S^2 model: positive half flat on both normal halves", flat(spectral(Bsph, 1)) < 1e-8)
c.done({"modes": 2 * 6, "layers": len(layers)})
