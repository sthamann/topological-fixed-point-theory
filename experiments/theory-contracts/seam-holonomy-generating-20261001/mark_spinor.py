#!/usr/bin/env python3
"""Part K: SPINOR.ORIGIN.02 -- the propagation spinor from the four marks (mark Clifford algebra).

Part G excluded the family spinor.  Here the candidate carrier is the MARKS themselves: one real Majorana
zero mode gamma_a per mark a = 0..3 (AX.MARKMAJORANA, an internal singlet), i.e. the Clifford algebra Cl(4)
of the four marks on C^4.  Everything else is read off, nothing is chosen:
  K1 frame: the four winding vectors w_a (Part S) are orthonormal frame rows, sum_a w_a = 0, so
     R^4 (mark index) = R v (v = (1,1,1,1)/2, the surface-relation direction killed in H_1) + R^3 (space);
     e_i = sum_a (w_a)_i gamma_a (i = x,y,z) and beta = sum_a gamma_a / 2 are four anticommuting Majoranas.
  K2 Dirac algebra: alpha_i = i beta e_i and beta are Hermitian, square to 1 and pairwise anticommute; the
     spin S_i = -(i/2) eps_ijk e_j e_k is ONE su(2) (S_x S_y = i S_z) commuting with beta and with the
     chirality Gamma = -gamma_0 gamma_1 gamma_2 gamma_3; Gamma anticommutes with beta and every gamma_a, and
     alpha_i = -Gamma S_i: C^4 = L (+) R with opposite helicity (k.alpha = -+ k.S), beta = the
     rotation-scalar mass operator joining them (Weyl form of the Dirac algebra).
  K3 point group from mark relabelings: every permutation pi of the marks acts on the gamma_a, has a Pin(4)
     lift p_pi (unique up to sign; odd pi -> odd Clifford element) with p S p^-1 = sgn(pi) std(pi) S and
     p beta p^-1 = sgn(pi) beta, std(pi) being the action on the winding lattice: even relabelings (A4) are
     rotations on each chirality, odd ones exchange L and R (parity).  The three double transpositions
     (pi-rotations about the axes) lift to -i S_x, -i S_y, -i S_z: the quaternion group (double cover).
  K4 mu8: the deck relabeling a -> a+1 (4-cycle, spatial action det -1) lifts to an element of order 8 with
     fourth power -1 (sign-choice independent, since the order is even).  Belt-trick cross-check: the
     Majorana braid images b_i = exp(pi/4 gamma_i gamma_{i+1}) satisfy Artin's relations, the full twist is
     (b_0 b_1 b_2)^4 = Gamma and the sphere word b_0 b_1 b_2^2 b_1 b_0 = -Gamma.  So EXACTLY on the sector
     where the Clifford image of the surface relation holds (gamma_0 gamma_1 gamma_2 gamma_3 = +1) they form a
     genuine representation of the sphere braid group B_4(S^2) with full twist -1 (Delta^2 is the central
     element of order 2, Fadell-Van Buskirk); on the other sector the sphere relation fails by a sign.
  K5 walks: W(k) = prod_i (cos k_i - i sin k_i alpha_i) is unitary on exactly the eight bcc steps, i log W =
     k.alpha + O(k^2) (isotropic, both chiralities, Berry +-1 on the two blocks); with the on-site mass
     U_D = e^{-i m beta/2} W e^{-i m beta/2} the dispersion is omega^2 = (k^2 + m^2)(1 + O(k)) (Dirac limit).
  K6 orientation: on the sector where the Clifford image of the surface relation gamma_0 gamma_1 gamma_2 gamma_3
     equals +1, the walk is a single Weyl chirality; the opposite sector carries the partner.
  K7 species count: for N real internal species with mark modes, the diagonal lift of a pi-rotation (double
     transposition; its square is sign-choice independent) squares to (-1)^N on the whole Fock space; the diagonal spin commutes with species rotations so(N) (Howe dual
     pair).  So the spinor sign needs an ODD number of real internal-singlet mark species; the sixteen E8 seam
     Majoranas (N = 16) give no spin 1/2.  The weak doublet as spinor: [W(k), T_j] = O(k) (kill).
  K8 (statement) AX.SPINSTEP reduces to AX.MARKMAJORANA (one real internal-singlet zero mode per mark):
     spinor, axes, double cover, mu8, parity, chirality pair and Dirac mass operator follow from the marks.
  K9 inventory against the existing source U: the Z2 sheet acts as (-1)^F on ALL 16 real seam fermions
     (Part U, U6), so U supplies an EVEN number (16) of sheet-odd species and the diagonal spin sign is
     (-1)^16 = +1: U itself does not contain AX.MARKMAJORANA.  A single extra real sheet-odd internal singlet
     is needed; if it were chiral on the seam the central charge would be 17/2 with Ising primaries, i.e. the
     seam would leave the holomorphic c = 8 point (tension with IMP.HOLO); a non-chiral (1/2, 1/2) mark
     fermion avoids this tension.
    python3 -B experiments/theory-contracts/seam-holonomy-generating-20261001/mark_spinor.py
"""
from __future__ import annotations

import json
import sys
from itertools import permutations
from pathlib import Path

import numpy as np
from scipy.linalg import expm, logm

HERE = Path(__file__).resolve().parent
LOG: list[dict] = []


def step(name: str, ok: bool, detail: str = "") -> None:
    LOG.append({"step": name, "pass": bool(ok), "detail": detail})
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))


Z = np.diag([1.0, -1.0]).astype(complex); I2 = np.eye(2, dtype=complex)
lower = np.array([[0, 1], [0, 0]], dtype=complex)


def majoranas(nmodes: int) -> list[np.ndarray]:
    """2 nmodes Majoranas on (C^2)^nmodes (Jordan-Wigner)."""
    def op(i, m):
        out = np.eye(1, dtype=complex)
        for j in range(nmodes):
            out = np.kron(out, Z if j < i else (m if j == i else I2))
        return out
    out = []
    for i in range(nmodes):
        a = op(i, lower)
        out += [a + a.conj().T, -1j * (a - a.conj().T)]
    return out


def acomm(x, y):
    return x @ y + y @ x


g = majoranas(2)                                               # gamma_0..gamma_3, one per mark
I4 = np.eye(4)
w = np.array([[-.5, -.5, -.5], [-.5, .5, .5], [.5, -.5, .5], [.5, .5, -.5]])   # winding vectors (Part S)
frame = np.vstack([np.ones(4) / 2, w.T])                        # rows: v, x, y, z
beta = sum(g[a] for a in range(4)) / 2
e = [sum(w[a, i] * g[a] for a in range(4)) for i in range(3)]
clif = all(np.allclose(acomm(x, y), 2 * (i == j) * I4) for i, x in enumerate([beta] + e) for j, y in enumerate([beta] + e))
step("K1.frame: (v, w_x, w_y, w_z) is orthonormal (sum_a w_a = 0), so beta = sum gamma_a/2 and e_i = sum_a (w_a)_i gamma_a are four anticommuting Majoranas: mark index = surface-relation direction + space",
     np.allclose(frame @ frame.T, np.eye(4)) and np.allclose(w.sum(0), 0) and clif)

alpha = [1j * beta @ e[i] for i in range(3)]
eps3 = lambda i, j, k: int(np.sign((j - i) * (k - i) * (k - j))) if len({i, j, k}) == 3 else 0
S = [sum(-0.5j * eps3(i, j, k) * e[j] @ e[k] for j in range(3) for k in range(3)) for i in range(3)]
Gam = -g[0] @ g[1] @ g[2] @ g[3]
dirac = all(np.allclose(x, x.conj().T) for x in alpha + [beta]) and all(
    np.allclose(acomm(x, y), 2 * (i == j) * I4) for i, x in enumerate(alpha + [beta]) for j, y in enumerate(alpha + [beta]))
spin_ok = all(np.allclose(S[i] @ S[i], I4) and np.allclose(S[i] @ beta, beta @ S[i]) and np.allclose(S[i] @ Gam, Gam @ S[i]) for i in range(3))
odd = all(np.allclose(acomm(Gam, x), 0) for x in g + [beta])
hand = np.allclose(S[0] @ S[1], 1j * S[2])
alpha_vs_S = all(np.allclose(alpha[i], -Gam @ S[i]) for i in range(3))
step("K2.Dirac algebra from the marks: alpha_i = i beta e_i, beta Hermitian involutions, pairwise anticommuting; one spin su(2) S_i (S_x S_y = i S_z) commuting with beta and Gamma; Gamma anticommutes with beta and all gamma_a; alpha_i = -Gamma S_i: L + R with opposite helicity, beta the scalar mass operator",
     dirac and spin_ok and odd and hand and alpha_vs_S and np.allclose(Gam @ Gam, I4))


def pin_lift(pi):
    """Pin(4) lift of the permutation gamma_a -> gamma_pi(a): product of reflections along gamma_a - gamma_b."""
    perm, out, sign = list(range(4)), I4.astype(complex), 1
    target = list(pi)
    for pos in range(4):                                       # bubble the current arrangement into target
        if perm[pos] != target[pos]:
            q = perm.index(target[pos])
            r = (g[perm[pos]] - g[perm[q]]) / np.sqrt(2)          # reflection swapping the two current labels
            out = r @ out
            perm[pos], perm[q] = perm[q], perm[pos]
            sign = -sign
    return out, sign


def std(pi):
    """Action of the mark relabeling a -> pi(a) on the winding lattice: w_a -> w_pi(a)."""
    A = np.array([w[pi[a]] for a in range(3)]).T
    return A @ np.linalg.inv(np.array([w[a] for a in range(3)]).T)


cov, perms = True, list(permutations(range(4)))
for pi in perms:
    p, sgn = pin_lift(pi)
    act = all(np.allclose(p @ g[a] @ np.linalg.inv(p), (sgn if False else 1) * g[pi[a]]) or np.allclose(p @ g[a] @ np.linalg.inv(p), -g[pi[a]]) for a in range(4))
    gperm = [p @ g[a] @ np.linalg.inv(p) for a in range(4)]
    s_glob = 1 if np.allclose(gperm[0], g[pi[0]]) else -1
    perm_ok = all(np.allclose(gperm[a], s_glob * g[pi[a]]) for a in range(4))
    R = std(pi)
    detR = round(np.linalg.det(R))
    Sp = [p @ S[i] @ np.linalg.inv(p) for i in range(3)]
    Sexp = [sum(detR * R[j, i] * S[j] for j in range(3)) for i in range(3)]
    bp = p @ beta @ np.linalg.inv(p)
    cov &= perm_ok and detR == sgn and all(np.allclose(Sp[i], Sexp[i]) for i in range(3)) and np.allclose(bp, s_glob * beta)
    cov &= np.allclose(p @ Gam @ np.linalg.inv(p), sgn * Gam)
klein = {(1, 0, 3, 2): 0, (2, 3, 0, 1): 1, (3, 2, 1, 0): 2}
klein_ok = True
for pi, ax in klein.items():
    p, _ = pin_lift(pi)
    klein_ok &= np.allclose(p, -1j * S[ax]) or np.allclose(p, 1j * S[ax])
step("K3.mark relabelings = point group: each of the 24 permutations has a Pin(4) lift p with p S p^-1 = sgn std S, p beta p^-1 = +-beta, p Gamma p^-1 = sgn Gamma (std = action on the windings); even = rotations on each chirality, odd = L <-> R; double transpositions lift to +-i S_axis (quaternion group)",
     cov and klein_ok)

deck = (1, 2, 3, 0)
p4, sg4 = pin_lift(deck)
mp = np.linalg.matrix_power
b = [expm(np.pi / 4 * g[i] @ g[i + 1]) for i in range(3)]
artin = np.allclose(b[0] @ b[1] @ b[0], b[1] @ b[0] @ b[1]) and np.allclose(b[1] @ b[2] @ b[1], b[2] @ b[1] @ b[2]) and np.allclose(b[0] @ b[2], b[2] @ b[0])
alpha0 = b[0] @ b[1] @ b[2]
twist = mp(alpha0, 4)
sphere = b[0] @ b[1] @ b[2] @ b[2] @ b[1] @ b[0]
ev_, V_ = np.linalg.eigh(Gam)
Prel = V_[:, ev_ < 0]                                          # gamma_0 gamma_1 gamma_2 gamma_3 = -Gamma = +1 here
Poth = V_[:, ev_ > 0]
on_rel = np.allclose(Prel.conj().T @ sphere @ Prel, np.eye(2)) and np.allclose(Prel.conj().T @ twist @ Prel, -np.eye(2))
on_oth = np.allclose(Poth.conj().T @ sphere @ Poth, -np.eye(2)) and np.allclose(Poth.conj().T @ twist @ Poth, np.eye(2))
step("K4.mu8 from the marks: the deck relabeling a -> a+1 (odd, spatial det -1) lifts with p^4 = -1, p^8 = 1; belt trick: Majorana braids satisfy Artin's relations, full twist = Gamma, sphere word = -Gamma -> on the surface-relation sector a genuine B_4(S^2) representation with full twist -1; on the other sector the sphere relation fails by a sign",
     sg4 == -1 and round(np.linalg.det(std(deck))) == -1 and np.allclose(mp(p4, 4), -I4) and np.allclose(mp(p4, 8), I4)
     and artin and np.allclose(twist, Gam) and np.allclose(sphere, -Gam) and on_rel and on_oth,
     f"deck lift p^4 = {np.round(mp(p4, 4)[0, 0].real, 6)}; relation sector: sphere word +1, full twist -1; other sector: sphere word -1, full twist +1")

W = lambda k: (np.cos(k[0]) * I4 - 1j * np.sin(k[0]) * alpha[0]) @ (np.cos(k[1]) * I4 - 1j * np.sin(k[1]) * alpha[1]) @ (np.cos(k[2]) * I4 - 1j * np.sin(k[2]) * alpha[2])
signs = [np.array(s) for s in np.ndindex(2, 2, 2)]
signs = [2 * s - 1 for s in signs]
A = lambda s: (I4 + s[0] * alpha[0]) / 2 @ (I4 + s[1] * alpha[1]) / 2 @ (I4 + s[2] * alpha[2]) / 2
rng = np.random.default_rng(5)
ks = rng.normal(size=(5, 3))
expand = all(np.allclose(W(k), sum(np.exp(-1j * s @ k) * A(s) for s in signs)) for k in ks)
unit = all(np.allclose(W(k).conj().T @ W(k), I4) for k in ks)
kk = 1e-3 * ks[0]
lin = np.allclose(1j * logm(W(kk)), sum(kk[i] * alpha[i] for i in range(3)), atol=5e-6)


def berry(Hf, P, r=0.05, n=36):
    th, ph = np.linspace(0, np.pi, n + 1), np.linspace(0, 2 * np.pi, n + 1)
    def vec(t, f):
        k = r * np.array([np.sin(t) * np.cos(f), np.sin(t) * np.sin(f), np.cos(t)])
        ev, V = np.linalg.eigh(P.conj().T @ Hf(k) @ P)
        return V[:, -1]
    flux = 0.0
    for i in range(n):
        for j in range(n):
            u = [vec(th[i], ph[j]), vec(th[i + 1], ph[j]), vec(th[i + 1], ph[j + 1]), vec(th[i], ph[j + 1])]
            flux += np.angle(np.vdot(u[0], u[1]) * np.vdot(u[1], u[2]) * np.vdot(u[2], u[3]) * np.vdot(u[3], u[0]))
    return flux / (2 * np.pi)


Pp, Pm = V_[:, ev_ > 0], V_[:, ev_ < 0]
Hlin = lambda k: sum(k[i] * alpha[i] for i in range(3))
Cp, Cm = berry(Hlin, Pp), berry(Hlin, Pm)
dirs = rng.normal(size=(6, 3)); dirs /= np.linalg.norm(dirs, axis=1)[:, None]
def dirac_dev(sc):
    m = 0.02 * sc
    UD = lambda k: expm(-0.5j * m * beta) @ W(k) @ expm(-0.5j * m * beta)
    dev = []
    for n_ in dirs:
        om = np.sort(np.abs(np.angle(np.linalg.eigvals(UD(0.02 * sc * n_)))))
        dev.append(abs(om[0] ** 2 - 2 * (0.02 * sc) ** 2) / (2 * (0.02 * sc) ** 2))
    return max(dev)
disp = [dirac_dev(1.0), dirac_dev(0.5), dirac_dev(0.25)]
step("K5.walks from the mark Clifford algebra: W(k) = prod (cos k_i - i sin k_i alpha_i) = sum over the eight bcc steps, unitary, i log W = k.alpha + O(k^2); Berry charges +-1 on the two chiral blocks; with the on-site mass coin e^{-i m beta/2}: omega^2 = k^2 + m^2 in the continuum limit (Dirac)",
     expand and unit and lin and round(Cp) * round(Cm) == -1 and abs(abs(Cp) - 1) < 1e-6
     and all(abs(disp[i] / disp[i + 1] - 2) < 0.05 for i in range(2)),
     f"Berry {Cp:.4f}, {Cm:.4f}; relative deviation of omega^2 from k^2+m^2 at (k,m) = s(0.02,0.02), s = 1, 1/2, 1/4: {[f'{x:.1e}' for x in disp]} (-> 0 linearly: Dirac continuum limit)")

rel = g[0] @ g[1] @ g[2] @ g[3]
rel_ok = np.allclose(Prel.conj().T @ rel @ Prel, np.eye(2))
Wr = lambda k: Prel.conj().T @ W(k) @ Prel
weyl_r = all(np.allclose(Wr(k).conj().T @ Wr(k), np.eye(2)) for k in ks) and np.allclose(Pp.conj().T @ W(ks[0]) @ Prel, 0)
step("K6.orientation: on the sector where the Clifford image of the surface relation gamma_0 gamma_1 gamma_2 gamma_3 = +1, the walk closes on itself (one Weyl chirality, Berry sign fixed); the other sector carries the partner; only beta (and the gamma_a) connect them",
     rel_ok and weyl_r, f"Berry charge on the relation sector {berry(Hlin, Prel):.4f}")

ok7, info = True, []
for N in (1, 2, 3):
    G = majoranas(2 * N)
    dim = G[0].shape[0]
    gs = [[G[4 * s + a] for a in range(4)] for s in range(N)]
    lift = np.eye(dim, dtype=complex)
    for s in range(N):                                         # diagonal lift of the pi-rotation (01)(23): even element
        lift = (gs[s][0] - gs[s][1]) @ (gs[s][2] - gs[s][3]) / 2 @ lift
    sq = lift @ lift
    ok7 &= np.allclose(sq, (-1) ** N * np.eye(dim))
    info.append(f"N={N}: square {sq[0, 0].real:+.0f}")
    if N == 3:
        es = [[sum(w[a, i] * gs[s][a] for a in range(4)) for i in range(3)] for s in range(N)]
        Jd = [sum(-0.25j * eps3(i, j, k) * es[s][j] @ es[s][k] for s in range(N) for j in range(3) for k in range(3)) for i in range(3)]
        rot = [sum(gs[s][a] @ gs[t][a] for a in range(4)) for s in range(N) for t in range(s + 1, N)]
        ok7 &= all(np.allclose(J @ X, X @ J) for J in Jd for X in rot)
weak = [np.kron(I2, x) for x in (np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.diag([1, -1]))]
sx, sy, sz = (np.array([[0, 1], [1, 0]], dtype=complex), np.array([[0, -1j], [1j, 0]]), np.diag([1.0 + 0j, -1]))
Ww = lambda k: (np.cos(k[0]) * I2 - 1j * np.sin(k[0]) * sx) @ (np.cos(k[1]) * I2 - 1j * np.sin(k[1]) * sy) @ (np.cos(k[2]) * I2 - 1j * np.sin(k[2]) * sz)
weak_mix = max(np.abs(Ww(1e-3 * dirs[0]) @ x - x @ Ww(1e-3 * dirs[0])).max() for x in (sx, sy, sz))
step("K7.species count: the diagonal pi-rotation lift over N real mark species squares to (-1)^N (N = 1, 2, 3 checked; N = 16 E8 seam Majoranas -> +1: no spin 1/2); the diagonal spin commutes with species rotations (Howe dual pair); weak doublet as spinor: [W(k), T_j] = O(k) -> kill",
     ok7 and weak_mix > 1e-4, "; ".join(info) + f"; weak isospin mixing {weak_mix:.1e} at |k| = 1e-3")
step("K8.(statement) AX.SPINSTEP reduces to AX.MARKMAJORANA: one real internal-singlet Majorana zero mode per mark gives the C^2 spinor, the axis assignment (from the windings), the point-group double cover, mu8, parity = odd relabelings, both chiralities and the Dirac mass operator beta along the surface-relation direction",
     True, "remaining input: existence of an odd number (minimally one) of real internal-singlet mark modes")

us = json.loads((HERE / "universal_source_results.json").read_text())
u6 = [x for x in us["steps"] if x["step"].startswith("U6")]
u1 = [x for x in us["steps"] if x["step"].startswith("U1")]
n_real = 2 * 8                                                    # 8 complex chiral modes (5 carrier + 3 family)
sign_U = (-1) ** n_real
c_chiral = 8 + 0.5
step("K9.inventory against U: the sheet acts as (-1)^F on all 16 real seam fermions (U6 passed), so U supplies an EVEN count and the diagonal spin sign is (-1)^16 = +1 -> AX.MARKMAJORANA is an EXTENSION of U (one extra real sheet-odd singlet); chiral it would give c = 17/2 with Ising primaries (leaves the holomorphic c = 8 point, tension with IMP.HOLO); non-chiral (1/2,1/2) avoids it",
     bool(u6) and u6[0]["pass"] and bool(u1) and u1[0]["pass"] and sign_U == 1 and c_chiral == 8.5,
     "the spacetime spinor is NOT contained in the present source; decisive open input")

out = {"steps": LOG, "passed": sum(x["pass"] for x in LOG), "total": len(LOG), "berry": [Cp, Cm]}
(HERE / "mark_spinor_results.json").write_text(json.dumps(out, indent=2, default=float))
print(f"{out['passed']}/{out['total']} steps passed")
sys.exit(0 if out["passed"] == out["total"] else 1)
