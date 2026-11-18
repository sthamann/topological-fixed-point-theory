#!/usr/bin/env python3
"""Part W: SEAM.NEUTRAL.WEYL.01 -- charged access refuted, neutral current channel tested.

All microscopic numbers use the ORIGINAL QWZ strip blocks (all 16 bands) of
microscopic-charged-car-limit, p_j = 2 pi (j - 1/4)/N, D_N = N h/(2 pi).
  W1 exact: charged positive access D_nu(-2,-1)_ii = -(1/2pi)(1/(1+nu) + 1/(2-nu)) = -24/(35 pi) at
     nu = 1/4, geometric M(-2)-M(-1) = -1/(2 pi); difference -13/(70 pi) != 0 (external claim).
  W2 microscopic charged values [K_N(-2) - K_N(-1)]_ii (full bands) reproduce the external table and
     tend to -24/(35 pi): the direct charged identification is REFUTED (reference-independent).
  W3 neutral channel: rho_{-n} = sum_x e^{2 pi i n x/N} Psi_x^dag Psi_x on the top edge creates
     particle-hole pairs with total weight -> |n| concentrated at energy |n| (holonomy cancels).
  W4 neutral Weyl function: with the access weights |a_n|^2 = 1/(2 pi |n|) (top edge for one sign of n,
     bottom edge = opposite chirality for the other), K(-2) - K(-1) tends to the GEOMETRIC values
     (-1/(2 pi) on the diagonal, the cosine kernel off it) -- full dynamical agreement.
  W5 the weight 1/sqrt|n| is the unique isometric identification of the one-current space with
     L^2_0(S^1) intertwining H and |D| (J_{-n} Omega / sqrt n orthonormal); no local linear neutral
     field d^k phi has it (their weights are n^(2k-1)): the geometric access is the H^{1/2} boundary
     trace (Newton-Wigner) localisation -- canonical, but not a local field insertion.
    PYTHONPATH=. python3 -B experiments/theory-contracts/seam-holonomy-generating-20261001/neutral_weyl.py
"""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import sympy as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
LOG: list[dict] = []


def step(name: str, ok: bool, detail: str = "") -> None:
    LOG.append({"step": name, "pass": bool(ok), "detail": detail})
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))


# W1 exact telescoping -----------------------------------------------------------
nu = sp.Rational(1, 4)
Dnu = -(1 / (2 * sp.pi)) * (1 / (1 + nu) + 1 / (2 - nu))
m = sp.symbols("m", integer=True, nonnegative=True)
Dgeo = sp.summation(1 / (m + 3) - 1 / (m + 2), (m, 0, sp.oo)) / sp.pi
step("W1.exact: charged D_{1/4}(-2,-1)_ii = -24/(35 pi), geometric = -1/(2 pi), difference -13/(70 pi)",
     sp.simplify(Dnu + sp.Rational(24, 35) / sp.pi) == 0 and sp.simplify(Dgeo + 1 / (2 * sp.pi)) == 0
     and sp.simplify(Dnu - Dgeo + sp.Rational(13, 70) / sp.pi) == 0)

# microscopic strip blocks ---------------------------------------------------------
spec = importlib.util.spec_from_file_location(
    "car_limit", ROOT / "experiments/theory-contracts/microscopic-charged-car-limit/checker.py")
car = importlib.util.module_from_spec(spec)
spec.loader.exec_module(car)
r_top, r_bot = car.row_spinor("top"), car.row_spinor("bottom")
_, source = car.inherited()                                   # validates the source pins


def blocks(N):
    labels = list(range(-N // 2 + 1, N // 2 + 1))
    ev, vec = [], []
    for j in labels:
        h = np.asarray(source.strip_at_momentum(2 * np.pi * (j - .25) / N, 8), dtype=complex)  # original block, all p
        e, v = np.linalg.eigh(h)
        ev.append(e); vec.append(v)
    return labels, np.array(ev), np.array(vec)


res = {}
for N in (128, 512, 2048):
    labels, ev, vec = blocks(N)
    A = N / (2 * np.pi) * np.abs(ev)                         # positive comparison |D_N|
    w_top = np.abs(np.einsum("i,jib->jb", r_top.conj(), vec)) ** 2
    w_bot = np.abs(np.einsum("i,jib->jb", r_bot.conj(), vec)) ** 2
    charged = ((w_top / (A + 2)) - (w_top / (A + 1))).sum() / (2 * np.pi)
    res[N] = dict(labels=labels, ev=ev, w_top=w_top, w_bot=w_bot, charged=charged)
ext = {128: -0.2220908092, 512: -0.2196902728, 2048: -0.2186719860}
step("W2.microscopic charged [K_N(-2)-K_N(-1)]_ii (all 16 bands) reproduces the external table and approaches -24/(35 pi), not -1/(2 pi)",
     all(abs(res[N]["charged"] - ext[N]) < 1e-8 for N in ext)
     and abs(res[2048]["charged"] - float(-sp.Rational(24, 35) / sp.pi)) < 1e-3
     and abs(res[2048]["charged"] + 1 / (2 * np.pi)) > 5e-2,
     ", ".join(f"N={N}: {res[N]['charged']:.10f}" for N in ext) + f"; limit {float(-sp.Rational(24, 35) / sp.pi):.10f}")


def neutral(N, edge, nmax=None, zs=(-2.0, -1.0)):
    """Particle-hole spectral data of rho_{-n} on one edge: weight S(n), weight near |n|, resolvents."""
    d = res[N]
    ev, w = d["ev"], d["w_top"] if edge == "top" else d["w_bot"]
    E = N / (2 * np.pi) * ev
    occ = ev < 0
    nmax = nmax or N // 2 - 1
    out = {}
    for n in range(-nmax, nmax + 1):
        if n == 0:
            continue
        j2 = np.roll(np.arange(N), -n)                      # label j -> j + n (periodic labels)
        wo = np.where(occ, w, 0.0)                          # hole in occupied band at j
        we = np.where(~occ, w, 0.0)[j2]                     # particle in empty band at j + n
        dE = E[j2][:, None, :] - E[:, :, None]               # (j, b_occ, b_emp)
        amp = wo[:, :, None] * we[:, None, :]
        low = dE < N / 16                                   # below the bulk gap scale N*gap/(2 pi)
        S = amp[low].sum()
        near = amp[np.abs(dE - abs(n)) < 0.5].sum()
        out[n] = dict(S=S, near=near, bulk=amp[~low].sum(), R={z: (amp / (dE - z)).sum() for z in zs})
    return out


Nn = 512
top, bot = neutral(Nn, "top"), neutral(Nn, "bottom")
edge_sign = -1 if sum(top[n]["near"] for n in range(-8, 0)) > sum(top[n]["near"] for n in range(1, 9)) else 1
chk = [(n, top[edge_sign * n]["S"], top[edge_sign * n]["near"]) for n in (1, 2, 4, 8)]
step("W3.neutral top-edge density channel: the low-energy weight S(n) -> |n|, concentrated at energy |n| (holonomy cancels); bulk interband weight sits at E ~ N gap",
     all(abs(S / n - 1) < 0.05 and near / S > 0.95 for n, S, near in chk),
     "; ".join(f"n={n}: S/n={S / n:.4f}, at E=n {near / S:.4f}" for n, S, near in chk))


def neutral_weyl(theta, N, top, bot, sgn):
    """(1/2pi) sum_n e^{i n theta} (1/|n|) R_n: top edge for the edge sign, bottom edge for the other."""
    tot = 0.0
    for n in range(1, N // 2 - 1):
        Rt = top[sgn * n]["R"]; Rb = bot[-sgn * n]["R"]
        tot += (np.exp(1j * n * theta) * (Rt[-2.0] - Rt[-1.0]) + np.exp(-1j * n * theta) * (Rb[-2.0] - Rb[-1.0])) / n
    return tot / (2 * np.pi)


geo = lambda th: sum(np.cos(n * th) * (1 / (n + 2) - 1 / (n + 1)) for n in range(1, 200000)) / np.pi
thetas = (0.0, np.pi / 2, np.pi)
gvals = {th: geo(th) for th in thetas}
conv = {}
for NN in (128, 512, 2048):
    tp, bt = neutral(NN, "top"), neutral(NN, "bottom")
    sg = -1 if sum(tp[n]["near"] for n in range(-4, 0)) > sum(tp[n]["near"] for n in range(1, 5)) else 1
    vv = {th: neutral_weyl(th, NN, tp, bt, sg) for th in thetas}
    conv[NN] = (max(abs(vv[t].real - gvals[t]) for t in thetas), max(abs(vv[t].imag) for t in thetas), vv)
vals = conv[2048][2]
devs = [conv[NN][0] for NN in (128, 512, 2048)]
step("W4.neutral Weyl function from the microscopic source: K(-2)-K(-1) -> geometric (-1/2pi diagonal, cosine kernel at pi/2, pi), real; deviation shrinks with N",
     devs[2] < devs[1] < devs[0] and devs[2] < 5e-4 and conv[2048][1] < 1e-3,
     "max dev " + ", ".join(f"N={NN}: {d:.1e}" for NN, d in zip((128, 512, 2048), devs))
     + "; N=2048 " + "; ".join(f"th={t:.3f}: {vals[t].real:.6f} vs {gvals[t]:.6f}" for t in thetas))

# W5 canonical weight vs local fields ---------------------------------------------------
K = 9
occ = list(range(1, K + 1))                                  # chiral Fock model: modes j, occupied j >= 1
def J_minus(n):                                             # J_{-n} Omega = sum_j a_j^dag a_{j+n} Omega
    return [(j, j + n) for j in range(1 - n, 1)]            # create at j <= 0, annihilate at j + n >= 1
norms = [len(J_minus(n)) for n in range(1, 7)]
energies = {((0.25 - j) - (0.25 - (j + n))) for n in range(1, 7) for j, _ in J_minus(n)} - set(range(1, 7))
local_weights = {k: [n ** (2 * k - 1) for n in range(1, 6)] for k in range(0, 4)}
step("W5.||J_{-n}Omega||^2 = n at energy exactly n; the isometry e_n -> J_{-n}Omega/sqrt n needs |a_n|^2 = 1/n; no local d^k phi (weights n^(2k-1)) has weight 1",
     norms == [1, 2, 3, 4, 5, 6] and not energies and all(w != [1] * 5 for w in local_weights.values()),
     "the geometric seam access = H^{1/2} boundary trace / Newton-Wigner localisation (canonical, non-local in the current algebra)")

out = {"steps": LOG, "passed": sum(x["pass"] for x in LOG), "total": len(LOG),
       "charged": {N: res[N]["charged"] for N in res},
       "neutral_weyl": {f"{t:.6f}": [vals[t].real, gvals[t]] for t in vals},
       "neutral_weyl_convergence": {NN: conv[NN][0] for NN in conv}}
(HERE / "neutral_weyl_results.json").write_text(json.dumps(out, indent=2, default=float))
print(f"{out['passed']}/{out['total']} steps passed")
sys.exit(0 if out["passed"] == out["total"] else 1)
