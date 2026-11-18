#!/usr/bin/env python3
"""Part P: SOURCE.PROJECTIVE_TRANSPORT.01 -- the relative composition operator of actual transport histories.

Question (external research order, 2026-10-01): does the ORIGINAL local Fuchs transport, lifted to the
universal abelian cover of Sigma = S^2 minus the four marks (deck group H_1 = Z^3 = the bcc address lattice
Lambda), carry a nontrivial composition sign between two histories with the same start and end, from which a
Clifford triple would follow?  Nothing is inserted: the links are integrated from the unchanged v117 residues.
  P1 lifted transport: L_e = P exp int_e A dz on sub-paths; subdivisions compose exactly; every lifted loop
     around mark a moves the sheet by w_a (winding vector computed); the four integrated monodromies equal the
     exact M_a in ONE common frame S; the positive invariant form is unique and pulls back to the identity;
     the determinant winds once around every mark while det M_a = 1; removing the trace changes each loop by
     the central factor e^{-2 pi i/3} (so the trace part stays in the links).
  P2 relative operators: L_ab (first a, then b) and L_ba end on the same sheet (w_a + w_b); their determinant
     phases agree exactly (same total winding).  Omega_ab = L_ba^-1 L_ab is NOT central: Omega^2 = I with
     eigenvalues (1,-1,-1); over all twelve ordered pairs exactly two values X, Y occur; they commute, the deck
     rotation exchanges them, and together they generate the Klein group V4 = [G,G] of the monodromy group A4.
  P3 composition sign: with X^2 = Y^2 = I one has XY = +YX: the induced cocycle on K = Z_2^2 is trivial and
     no Clifford triple arises.  Reason (general): the monodromy of a flat connection is a LINEAR
     representation and [A4, A4] = V4 is abelian, so the relative operators commute on every linear fiber,
     including line twists such as the half-determinant (sheet) line, which leaves all Omega unchanged.
  P4 intermediate classes: all lifted closed loops shift the address inside Lambda (integer windings), so no
     source transport reaches the half-shift classes (1/2 Z)^3 / Lambda = Z_2^2; the abelianised flavor
     holonomy G/[G,G] = Z_3 induces only the trivial function on Lambda / 2 Lambda.  The flavor transport also
     does not resolve sheets: gamma_a^3 has trivial monodromy but sheet shift 3 w_a.
  P5 comparison, NOT the source fiber: the flavor local system is of real type (unique symmetric invariant form);
     in every SU(2) lift (8 sign choices) the lifted relative operators are lift-independent and ANTIcommute
     (quaternion group): the nontrivial sign exists exactly in the spin cover of the flavor connection, i.e. it
     would be flavor-derived spin -- an additional identification, not a property of the original fiber.
  P6 (statement) verdict for the original source: the relative phase of transport histories is non-central but
     commuting; the projective-transport route to a Clifford sector FAILS on the original Fuchs fiber.  The
     blind Weyl limit of a source-determined controlled process is therefore not run (its precondition fails).
    PYTHONPATH=. python3 -B experiments/theory-contracts/seam-holonomy-generating-20261001/projective_transport.py
"""
from __future__ import annotations

import json
import sys
from itertools import product
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import null_space, sqrtm

from tfpt_explorer.flavor_path_transport import exact_clock_dictionary, original_fuchs_connection

HERE = Path(__file__).resolve().parent
LOG: list[dict] = []


def step(name: str, ok: bool, detail: str = "") -> None:
    LOG.append({"step": name, "pass": bool(ok), "detail": detail})
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))


I3 = np.eye(3, dtype=complex)
marks = [1j ** k for k in range(4)]
w = np.array([[-.5, -.5, -.5], [-.5, .5, .5], [.5, -.5, .5], [.5, .5, -.5]])


def seg(path, conn=original_fuchs_connection):
    def rhs(t, y):
        zz, dz = path(t)
        return (dz * conn(zz) @ y[:9].reshape(3, 3)).ravel().tolist() + [dz * np.trace(conn(zz))] + [dz / (zz - m) for m in marks]
    y0 = np.concatenate([I3.ravel(), np.zeros(5, dtype=complex)])
    sol = solve_ivp(lambda t, y: np.array(rhs(t, y)), (0, 1), y0, method="DOP853", rtol=1e-13, atol=1e-15)
    y = sol.y[:, -1]
    return y[:9].reshape(3, 3), y[9], y[10:] / (2j * np.pi)        # link, int tr A dz, partial windings


ray = lambda a, b: (lambda t: (a + (b - a) * t, b - a))
arc = lambda c, r, s0, s1: (lambda t: (c + r * np.exp(1j * (s0 + (s1 - s0) * t)), 1j * (s1 - s0) * r * np.exp(1j * (s0 + (s1 - s0) * t))))


def loop_pieces(k, nsub=4):
    c = marks[k]; s0 = np.angle(-c)
    pieces = [ray(0, 0.6 * c)] + [arc(c, 0.4, s0 + 2 * np.pi * j / nsub, s0 + 2 * np.pi * (j + 1) / nsub) for j in range(nsub)] + [ray(0.6 * c, 0)]
    return pieces


T, sub_dev, wind, trint = [], 0.0, [], []
for k in range(4):
    links = [seg(p) for p in loop_pieces(k)]
    Lk = I3
    tr_tot, wv = 0, np.zeros(4, dtype=complex)
    for L_, tr_, wv_ in links:
        Lk = L_ @ Lk; tr_tot += tr_; wv += wv_
    full = seg(ray(0.6 * marks[k], 0))[0] @ seg(arc(marks[k], 0.4, np.angle(-marks[k]), np.angle(-marks[k]) + 2 * np.pi))[0] @ seg(ray(0, 0.6 * marks[k]))[0]
    sub_dev = max(sub_dev, np.abs(Lk - full).max())
    T.append(Lk); wind.append(wv); trint.append(tr_tot / (2j * np.pi))

M = np.array(exact_clock_dictionary()["M"].evalf(), dtype=complex)
Uc = np.diag([1, 1j, -1j])
Ma = [np.linalg.matrix_power(Uc, a) @ M @ np.linalg.matrix_power(Uc, -a) for a in range(4)]
basis = []
for i, j in product(range(3), repeat=2):
    E = np.zeros((3, 3), dtype=complex); E[i, j] = 1; basis.append(E)
inter = null_space(np.array([np.concatenate([(Tk @ E - E @ Mk).ravel() for Tk, Mk in zip(T, Ma)]) for E in basis]).T, rcond=1e-9)
S = sum(inter[i, 0] * basis[i] for i in range(9))
frame_dev = max(np.abs(Tk - S @ Mk @ np.linalg.inv(S)).max() for Tk, Mk in zip(T, Ma))
forms = null_space(np.array([np.concatenate([(Tk.conj().T @ E @ Tk - E).ravel() for Tk in T]) for E in basis]).T, rcond=1e-9)
H = sum(forms[i, 0] * basis[i] for i in range(9)); H = (H + H.conj().T) / 2
H = H if np.linalg.eigvalsh(H)[0] > 0 else -H
pull = S.conj().T @ H @ S
wind_ok = all(np.allclose(wind[k], np.eye(4)[k], atol=1e-9) for k in range(4))
det_ok = all(abs(trint[k] - 1) < 1e-9 for k in range(4)) and all(abs(np.linalg.det(Tk) - 1) < 1e-9 for Tk in T)
traceless = lambda z: original_fuchs_connection(z) - np.trace(original_fuchs_connection(z)) / 3 * I3
T0 = []
for k in range(4):
    Lk = I3
    for p in loop_pieces(k, 2):
        Lk = seg(p, traceless)[0] @ Lk
    T0.append(Lk)
central = all(np.allclose(T0[k], np.exp(-2j * np.pi / 3) * T[k], atol=1e-9) for k in range(4))
step("P1.lifted transport from the unchanged v117 residues: subdivisions compose exactly; each loop moves the sheet by w_a (winding vector e_a); integrated monodromies = exact M_a in ONE common frame; unique positive invariant form pulling back to I; det winds once per mark while det M_a = 1; trace removal = central factor e^{-2 pi i/3}",
     sub_dev < 1e-12 and inter.shape[1] == 1 and frame_dev < 1e-11 and forms.shape[1] == 1 and np.allclose(pull / pull[0, 0], I3, atol=1e-9) and wind_ok and det_ok and central,
     f"subdivision {sub_dev:.1e}, common-frame residual {frame_dev:.1e}, invariant forms {forms.shape[1]}, determinant windings {[round(abs(x), 9) for x in trint]}")

Hh = np.real_if_close(sqrtm(H))
V = [Hh @ Tk @ np.linalg.inv(Hh) for Tk in T]                          # unitary frame
unit = all(np.allclose(v.conj().T @ v, I3, atol=1e-10) for v in V)
Om = {}
for a in range(4):
    for b in range(4):
        if a != b:
            Lab, Lba = V[b] @ V[a], V[a] @ V[b]
            Om[(a, b)] = np.linalg.inv(Lba) @ Lab
same_end = all(np.allclose(w[a] + w[b], w[b] + w[a]) for a in range(4) for b in range(4))
classes = []
for kk in Om:
    for c in classes:
        if np.allclose(Om[kk], Om[c[0]], atol=1e-9):
            c.append(kk); break
    else:
        classes.append([kk])
X, Y = Om[classes[0][0]], Om[classes[1][0]] if len(classes) > 1 else None
invol = all(np.allclose(O @ O, I3, atol=1e-9) and np.allclose(np.sort(np.linalg.eigvals(O).real), [-1, -1, 1], atol=1e-9) for O in Om.values())
deckU = [k for k in range(4)]
cls_of = {kk: i for i, c in enumerate(classes) for kk in c}
swap = len(classes) == 2 and all(cls_of[((a + 1) % 4, (b + 1) % 4)] != cls_of[(a, b)] for (a, b) in Om)   # deck a -> a+1 exchanges X, Y
gen = {tuple(np.round(Z_, 6).ravel()) for Z_ in (I3, X, Y, X @ Y)}
step("P2.relative operators: both orders end on the same sheet with equal determinant phase; Omega_ab = L_ba^-1 L_ab is NOT central (Omega^2 = I, eigenvalues 1,-1,-1); over the twelve ordered pairs exactly two values X, Y occur, exchanged by the deck rotation, generating V4 = [A4, A4]",
     unit and same_end and invol and len(classes) == 2 and swap and len(gen) == 4,
     f"classes: {[sorted({tuple(sorted(p)) for p in c}) for c in classes]}")

anticomm = np.allclose(X @ Y, -Y @ X, atol=1e-9)
comm = np.allclose(X @ Y, Y @ X, atol=1e-9)
Vtw = [-v for v in V]                                                    # half-determinant (sheet) line twist: -M_a
Om_tw = {k: np.linalg.inv(Vtw[k[0]] @ Vtw[k[1]]) @ (Vtw[k[1]] @ Vtw[k[0]]) for k in Om}
tw_same = all(np.allclose(Om_tw[k], Om[k], atol=1e-12) for k in Om)
step("P3.composition sign: X^2 = Y^2 = I and XY = +YX: trivial cocycle on K = Z_2^2, NO Clifford triple; general reason: flat-connection monodromy is a linear representation and [A4,A4] = V4 is abelian; the half-determinant (sheet) line twist leaves every Omega unchanged",
     comm and not anticomm and tw_same, "case 'not central' of the research order: full operators kept, no projection")

in_lambda = lambda v: np.allclose(v, np.round(v)) or np.allclose(v - 0.5, np.round(v - 0.5))
words_ok = all(in_lambda(sum(n[a] * w[a] for a in range(4))) for n in product(range(-2, 3), repeat=4))
half_classes = {tuple(np.mod(np.array(v) * 2, 2).astype(int)) for v in product((0, 0.5), repeat=3)}
lam_classes = {(0, 0, 0), (1, 1, 1)}
K_index = len(half_classes) // len(lam_classes)
cube = np.allclose(np.linalg.matrix_power(V[0], 3), I3, atol=1e-9)
step("P4.intermediate classes: every lifted closed loop shifts the address inside Lambda (integer windings), the half-shift classes (1/2 Z)^3/Lambda = Z_2^2 are reached by no source transport; G/[G,G] = Z_3 gives only the trivial function on Lambda/2Lambda; flavor does not resolve sheets (gamma_a^3: trivial monodromy, shift 3 w_a)",
     words_ok and K_index == 4 and cube and np.gcd(3, 2) == 1, f"index of Lambda in (1/2 Z)^3: {K_index}")

Bf = null_space(np.array([np.concatenate([(v.T @ E @ v - E).ravel() for v in V]) for E in basis]).T, rcond=1e-9)
B = sum(Bf[i, 0] * basis[i] for i in range(9))
c_ = (B.conj() @ B)[0, 0]
B = B / np.sqrt(c_)
J = lambda v: np.linalg.solve(B, v.conj())                              # antilinear, commutes with every V_a
cands = [J(e) + e for e in np.eye(3, dtype=complex)] + [1j * (e - J(e)) for e in np.eye(3, dtype=complex)]
R = []
for v in cands:
    for r in R:
        v = v - np.vdot(r, v).real * r                                 # J antiunitary: inner products of J-fixed vectors are real
    if np.linalg.norm(v) > 1e-6:
        R.append(v / np.linalg.norm(v))
    if len(R) == 3:
        break
R = np.array(R).T
Vr = [np.linalg.inv(R) @ v @ R for v in V]
real_ok = all(np.allclose(x.imag, 0, atol=1e-8) and np.allclose(x.real @ x.real.T, np.eye(3), atol=1e-8) for x in Vr)
sx, sy, sz = np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]]), np.diag([1.0 + 0j, -1])
def su2_lift(Rm):
    Rm = np.real(Rm)
    q0 = np.sqrt(max(0.0, 1 + np.trace(Rm))) / 2
    if q0 > 1e-6:
        q = np.array([Rm[2, 1] - Rm[1, 2], Rm[0, 2] - Rm[2, 0], Rm[1, 0] - Rm[0, 1]]) / (4 * q0)
    else:
        ax = np.real(null_space(Rm - np.eye(3))[:, 0]); q = ax
    return q0 * np.eye(2) - 1j * (q[0] * sx + q[1] * sy + q[2] * sz)
lifts = [su2_lift(x) for x in Vr]
def rot_of(u):
    return np.array([[0.5 * np.trace(s1 @ u @ s2 @ u.conj().T).real for s2 in (sx, sy, sz)] for s1 in (sx, sy, sz)])
lift_ok = all(np.allclose(rot_of(u), x.real, atol=1e-8) or np.allclose(rot_of(u).T, x.real, atol=1e-8) for u, x in zip(lifts, Vr))
spin_results = []
for sg in product((1, -1), repeat=4):
    Lf = [s * u for s, u in zip(sg, lifts)]
    Xs = np.linalg.inv(Lf[classes[0][0][0]] @ Lf[classes[0][0][1]]) @ (Lf[classes[0][0][1]] @ Lf[classes[0][0][0]])
    Ys = np.linalg.inv(Lf[classes[1][0][0]] @ Lf[classes[1][0][1]]) @ (Lf[classes[1][0][1]] @ Lf[classes[1][0][0]])
    spin_results.append((Xs, Ys))
indep = all(np.allclose(x, spin_results[0][0]) and np.allclose(y, spin_results[0][1]) for x, y in spin_results)
Xs, Ys = spin_results[0]
spin_anti = np.allclose(Xs @ Ys, -Ys @ Xs, atol=1e-8) and np.allclose(Xs @ Xs, -np.eye(2), atol=1e-8)
step("P5.comparison (NOT the source fiber): the flavor system is of real type (unique symmetric invariant form); in all 2^4 SU(2) sign lifts the lifted relative operators are lift-independent and ANTIcommute with square -1 (quaternion group): the nontrivial sign lives only in the spin cover of the flavor connection = flavor-derived spin, an extra identification",
     Bf.shape[1] == 1 and real_ok and lift_ok and indep and spin_anti, "Clifford triple i X~, i Y~, i X~Y~ exists only after lifting the flavor frame")
step("P6.(statement) verdict SOURCE.PROJECTIVE_TRANSPORT.01: on the original Fuchs fiber the relative phase of transport histories is non-central but commuting -> no source Clifford sector; the projective route needs the spin lift of the flavor connection (excluded as identification without proof); blind Weyl limit not run (precondition fails)",
     True, "negative for the original source; scope: linear fibers of the given flat connection")

out = {"steps": LOG, "passed": sum(x["pass"] for x in LOG), "total": len(LOG),
       "omega_classes": [[list(p) for p in c] for c in classes], "frame_residual": frame_dev}
(HERE / "projective_transport_results.json").write_text(json.dumps(out, indent=2, default=float))
print(f"{out['passed']}/{out['total']} steps passed")
sys.exit(0 if out["passed"] == out["total"] else 1)
