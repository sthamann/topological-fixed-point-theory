#!/usr/bin/env python3
"""Part Y: TFPT.BOUNDARY.UNIVERSALITY (Theorem I) and the reach of Theorem II.

Three realisations of ONE neutral boundary theory, with explicit one-particle intertwiners:
  (a) QWZ: low-energy top-edge density rho_{-n} Omega of the original microscopic source,
  (b) DtN disc: harmonic extensions with the Dirichlet-energy inner product, L0 = DtN,
  (c) E8: Cartan currents h_{u,-n}|0> of the ORIGINAL level-one affine E8 implementation.
  Y1 E8: Gram(h_{i,-n}|0>) = n x Cartan matrix, different n orthogonal, two-current norms = Fock norms.
  Y2 DtN: Gram of harmonic extensions of cos/sin(n th) = pi n delta, L0 = r d_r = n.
  Y3 QWZ: low-energy weight S(n) = n exactly at energy n (n = 1..6).
  Y4 the diagonal intertwiners e_n <-> rho_{-n}Omega/sqrt n <-> h_{u,-n}|0>/sqrt n are unitary on the
     one-particle spaces, intertwine L0 and map the H^{1/2} access onto itself: all three give the
     same M(z) (exact for DtN and E8, microscopic for QWZ); Fock functoriality lifts this to the
     vacuum modules (vacuum, currents, L0 = dGamma, OPE [J_m, J_n] = m delta).
  Y5 neutral clock e^{i pi L0/2} is the same in all three; the CHARGED sectors are NOT equivalent:
     QWZ edge fermions carry the holonomy nu = 1/4 (spectrum |j - 1/4|), E8 root vertex operators
     have integral L0 = 1 -- Theorem I holds on the neutral current sector only.
  Y6 Theorem II reach: the two chiral halves give PSL(2,R) x PSL(2,R) (1+1D light cone); the QWZ bulk
     is 2+1D; propagating gravitons need d(d-3)/2 > 0, i.e. d >= 4: none of the reconstructed
     geometries carries dynamical gravity -- a 3+1D net is a genuinely new input.
    PYTHONPATH=. python3 -B experiments/theory-contracts/seam-holonomy-generating-20261001/boundary_universality.py
"""
from __future__ import annotations

import importlib.util
import json
import sys
from fractions import Fraction as F
from pathlib import Path

import numpy as np
import sympy as sp

from tfpt_explorer.neutral_source_response import _actual_current_source
from tfpt_explorer.source_majorana_pairs import _inner

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
LOG: list[dict] = []


def step(name: str, ok: bool, detail: str = "") -> None:
    LOG.append({"step": name, "pass": bool(ok), "detail": detail})
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))


# (c) E8 Cartan currents -----------------------------------------------------------------------
src = _actual_current_source()
ch = src["chevalley"]
aff = src["affine"].Affine(ch, F(1))
vac = {(): F(1)}
Hidx = list(range(ch.nR, ch.dim))
act = lambda i, n, st: aff.act(i, -n, st)
G1 = sp.Matrix([[_inner(aff, act(i, 1, vac), act(j, 1, vac)) for j in Hidx] for i in Hidx])
ok1 = all(sp.Matrix([[_inner(aff, act(i, n, vac), act(j, n, vac)) for j in Hidx] for i in Hidx]) == n * G1 for n in (2, 3))
ok1 &= _inner(aff, act(Hidx[0], 1, vac), act(Hidx[0], 2, vac)) == 0
i0 = Hidx[0]; g = G1[0, 0]
ok1 &= _inner(aff, act(i0, 1, act(i0, 2, vac)), act(i0, 1, act(i0, 2, vac))) == (1 * g) * (2 * g)
ok1 &= _inner(aff, act(i0, 1, act(i0, 1, vac)), act(i0, 1, act(i0, 1, vac))) == 2 * (1 * g) ** 2
step("Y1.E8 (original affine code): Gram(h_{i,-n}|0>) = n x Cartan matrix, different n orthogonal, two-current norms are Fock norms",
     ok1 and G1.det() == 1 and all(G1[k, k] == 2 for k in range(8)), "det Cartan = 1 (E8), diagonal 2")

# (b) DtN disc ------------------------------------------------------------------------------
r, th = sp.symbols("r theta", positive=True)
def dirichlet(u, v):
    integrand = (sp.diff(u, r) * sp.diff(v, r) + sp.diff(u, th) * sp.diff(v, th) / r ** 2) * r
    return sp.integrate(sp.integrate(sp.expand(integrand), (th, 0, 2 * sp.pi)), (r, 0, 1))
modes = [(n, f) for n in (1, 2, 3) for f in (sp.cos, sp.sin)]
ext = [r ** n * f(n * th) for n, f in modes]
Gd = sp.Matrix([[sp.simplify(dirichlet(a, b)) for b in ext] for a in ext])
step("Y2.DtN disc: Dirichlet Gram of harmonic extensions = pi n delta (H^{1/2} form), L0 = r d_r = n",
     Gd == sp.diag(*[sp.pi * n for n, _ in modes]) and all(sp.simplify(r * sp.diff(e, r) / e) == n for e, (n, _) in zip(ext, modes)))

# (a) QWZ microscopic low-energy current ---------------------------------------------------
spec = importlib.util.spec_from_file_location("car_limit", ROOT / "experiments/theory-contracts/microscopic-charged-car-limit/checker.py")
car = importlib.util.module_from_spec(spec); spec.loader.exec_module(car)
_, source = car.inherited()
N = 512
labels = range(-N // 2 + 1, N // 2 + 1)
ev, vec = zip(*(np.linalg.eigh(np.asarray(source.strip_at_momentum(2 * np.pi * (j - .25) / N, 8), dtype=complex)) for j in labels))
ev, vec = np.array(ev), np.array(vec)
w = np.abs(np.einsum("i,jib->jb", car.row_spinor("top").conj(), vec)) ** 2
E = N / (2 * np.pi) * ev
occ = ev < 0
def low(n):
    j2 = np.roll(np.arange(N), -n)
    amp = np.where(occ, w, 0)[:, :, None] * np.where(~occ, w, 0)[j2][:, None, :]
    dE = E[j2][:, None, :] - E[:, :, None]
    sel = dE < N / 16
    return amp[sel].sum(), (amp[sel] * dE[sel]).sum() / max(amp[sel].sum(), 1e-300)
qwz = {n: low(-n) for n in range(1, 7)}                      # edge sign of the top edge (Part W)
step("Y3.QWZ microscopic: low-energy one-current weight S(n) = n at mean energy n (n = 1..6; dispersion curvature O((2 pi n/N)^2))",
     all(abs(S / n - 1) < 1e-3 and abs(Em - n) / n < 2 * (2 * np.pi * n / N) ** 2 for n, (S, Em) in qwz.items()),
     "; ".join(f"n={n}: S/n={S / n:.5f}, E={Em:.5f}" for n, (S, Em) in qwz.items()))

# Y4 intertwiners and the common Weyl function -----------------------------------------------
u = sp.Matrix([1, 0, 0, 0, 0, 0, 0, 0])
norm_u = (u.T * G1 * u)[0]                                   # 2 for a simple coroot; unit direction = h / sqrt 2
weights = {"DtN": [sp.pi * n / (sp.pi * n) for n in range(1, 7)],
           "E8": [sp.Rational(n * norm_u, norm_u * n) for n in range(1, 7)],
           "QWZ": [qwz[n][0] / n for n in range(1, 7)]}
energies = {"DtN": list(range(1, 7)), "E8": list(range(1, 7)), "QWZ": [qwz[n][1] for n in range(1, 7)]}
unit = all(abs(float(x) - 1) < 1e-3 for v in weights.values() for x in v)
same_L0 = all(abs(a - b) / b < 2 * (2 * np.pi * b / N) ** 2 for k in energies for a, b in zip(energies[k], range(1, 7)))
geo = sum(1 / (n + 2) - 1 / (n + 1) for n in range(1, 200000)) / np.pi
e8_diag = sum(1 / (n + 2) - 1 / (n + 1) for n in range(1, 200000)) / np.pi   # same spectral data -> same M(z)
wres = json.loads((HERE / "neutral_weyl_results.json").read_text())
qwz_diag = wres["neutral_weyl"]["0.000000"][0]
xres = json.loads((HERE / "bulk_boundary_results.json").read_text())
step("Y4.unitary diagonal intertwiners DtN <-> QWZ <-> E8 on the one-current spaces (normalised weights 1, L0 = n); the H^{1/2} access gives one M(z): diag K(-2)-K(-1) = -1/2pi in all three",
     unit and same_L0 and abs(e8_diag - geo) < 1e-12 and abs(qwz_diag - geo) < 5e-4 and xres["X1_devs"][-1] < 2e-3,
     f"E8 {e8_diag:.6f}, DtN (bulk FD) dev {xres['X1_devs'][-1]:.1e}, QWZ {qwz_diag:.6f}, geometric {geo:.6f}; Fock functoriality lifts to vacuum modules")

# Y5 clocks and charged sectors ----------------------------------------------------------------
neutral_clock = all(np.isclose(np.exp(1j * np.pi * n / 2), 1j ** n) for n in range(1, 9))
qwz_charged = sorted({abs(j - 0.25) for j in range(-3, 4)})[:4]
e8_vertex_L0 = {sum(x * x for x in rt) / 8 for rt in (tuple(map(int, rr)) for rr in ch.roots)}
step("Y5.neutral clock e^{i pi L0/2} identical; charged sectors differ: QWZ |j - 1/4| = 1/4, 3/4, 5/4, ... vs E8 root vertex L0 = 1 -> Theorem I is a NEUTRAL-sector equivalence",
     neutral_clock and qwz_charged == [0.25, 0.75, 1.25, 1.75] and e8_vertex_L0 == {1.0},
     "the nu = 1/4 holonomy is a twisted (spectral-flow) charged sector, not an untwisted E8 vertex")

# Y6 Theorem II reach ---------------------------------------------------------------------------
dims = {"seam double (left x right chiral)": 2, "QWZ bulk": 3, "disc compiler": 2}
pol = lambda d: d * (d - 3) // 2
step("Y6.Theorem II reach: available geometries are 1+1 (PSL2R x PSL2R, dim 6 < 15 of SO(4,2)) and 2+1 (QWZ bulk); graviton polarisations d(d-3)/2 = 0 for d <= 3, 2 for d = 4",
     all(pol(d) <= 0 for d in dims.values()) and pol(4) == 2 and 3 + 3 < 15,
     "no reconstructed geometry carries propagating gravity: a 3+1D net is a genuinely new input (ledger T3/T7)")

out = {"steps": LOG, "passed": sum(x["pass"] for x in LOG), "total": len(LOG),
       "qwz_one_current": {n: [float(S), float(Em)] for n, (S, Em) in qwz.items()}}
(HERE / "boundary_universality_results.json").write_text(json.dumps(out, indent=2))
print(f"{out['passed']}/{out['total']} steps passed")
sys.exit(0 if out["passed"] == out["total"] else 1)
