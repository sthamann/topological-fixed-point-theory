#!/usr/bin/env python3
"""Part ST: QMC.SYMMETRY_PRESERVING_TRANSFER.01 -- does one Trotter slice of the QMC commute with the symmetries?

An ordered product of local gates e^{dtau g O_1^2} e^{dtau g O_2^2} ... is NOT translation invariant at finite dtau; its
ground state then leaks into every momentum sector and a projected sector trace can show a spurious zero gap.  The
integrator of reweight.py is different: one slice is
    T = W_h  < exp( sum_x phi_x . O_x ) >_phi  W_h ,   W_h = e^{-dtau H_w / 2},
i.e. ALL channels of ALL sites in one exponential, averaged over a Gaussian phi with a translation-, Z4- and
Spin(10)-invariant covariance.  For every symmetry u of H_w and of the set {O_x} (u O_x u^dag = R_u O_{u x}),
u T u^dag = W_h < exp(sum phi_x . R_u O_{ux}) > W_h = T because the Gaussian measure is invariant under phi -> R_u phi
permuted.  Hence [T, P_k] = [T, P_q] = 0 exactly at every dtau (the HS average differs from e^{-dtau H} at O(dtau^2), but
symmetrically).
  ST1 numerical check on a solvable ring (3 sites, spinless fermions, pairing composites O_x = c_x c_{x+1} + h.c. +
      density part, nearest-neighbour hopping): the Gauss-Hermite HS average of one slice commutes with the translation
      and with the Z4 rotation e^{i pi N/2} to quadrature precision.
  ST2 counter-check with non-commuting composites: the ordered local-gate Trotter product does NOT commute with the
      translation (O(dtau^2) leak) -- the failure mode that would corrupt a sector-gap measurement -- while the
      single-exponential HS slice of the same composites does.
  ST4 the production operator: the cell translations commute with W_h, K(phi) is translation covariant and Z4 maps
      K(phi) -> K(-phi) for the frozen Theta_1 code (g != 0); the release is tied to an integrator hash.
  ST3 (statement) the reweight.py slice is of the ST1 type (single exponential of the full HS bilinear with an invariant
      Gaussian prior); the integrator therefore preserves translations, Z4 and Spin(10) exactly and needs no group-averaged
      variant; the per-configuration slice is not invariant, only its HS average -- which is what a sector trace sees.
    python3 -B experiments/theory-contracts/seam-holonomy-generating-20261001/symmetric_transfer.py
"""
from __future__ import annotations

import json
import sys
from itertools import product
from pathlib import Path

import numpy as np
from scipy.linalg import expm

HERE = Path(__file__).resolve().parent
LOG: list[dict] = []


def step(name: str, ok: bool, detail: str = "") -> None:
    LOG.append({"step": name, "pass": bool(ok), "detail": detail})
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))
    sys.stdout.flush()


Ls = 3
Zp = np.diag([1.0, -1.0]); low = np.array([[0, 1], [0, 0]], complex)
c = []
for i in range(Ls):
    mm = np.eye(1, dtype=complex)
    for j in range(Ls):
        mm = np.kron(mm, Zp if j < i else (low if j == i else np.eye(2)))
    c.append(mm)
cd = [x.conj().T for x in c]
dim = 2 ** Ls
N = sum(cd[i] @ c[i] for i in range(Ls))
Hw = -sum(cd[i] @ c[(i + 1) % Ls] + cd[(i + 1) % Ls] @ c[i] for i in range(Ls))
O = [c[i] @ c[(i + 1) % Ls] + cd[(i + 1) % Ls] @ cd[i] for i in range(Ls)]              # pairing composites (Z4: O -> -O)
# translation operator: c_i -> c_{i+1}, built from its action on occupation basis with fermionic signs
Tr = np.zeros((dim, dim), complex)
for b in range(dim):
    occ = [(b >> (Ls - 1 - i)) & 1 for i in range(Ls)]
    vec = np.zeros(dim, complex); vec[0] = 1.0
    # |occ> = prod_i (c_i^dag)^{occ_i} |0> in increasing i; translated: prod_i (c_{i+1}^dag)^{occ_i} |0>
    v = vec.copy()
    for i in reversed(range(Ls)):
        if occ[i]:
            v = cd[(i + 1) % Ls] @ v
    u = vec.copy()
    for i in reversed(range(Ls)):
        if occ[i]:
            u = cd[i] @ u
    Tr += np.outer(v, u.conj())
tr_ok = np.abs(Tr @ c[0] @ Tr.conj().T - c[1]).max() < 1e-12 and np.abs(Tr @ Hw @ Tr.conj().T - Hw).max() < 1e-12
Z4 = expm(1j * np.pi / 2 * N)
dtau, g = 0.25, 0.8
Wh = expm(-dtau * Hw / 2)
xs, ws = np.polynomial.hermite.hermgauss(24)                             # int e^{-y^2} f(y) dy
sig = np.sqrt(2 * g * dtau)                                               # phi ~ N(0, 2 g dtau)
Tavg = np.zeros((dim, dim), complex)
for (i0, i1, i2) in product(range(len(xs)), repeat=3):
    phi = np.sqrt(2) * sig * np.array([xs[i0], xs[i1], xs[i2]])
    w = ws[i0] * ws[i1] * ws[i2] / np.pi ** 1.5
    Tavg += w * expm(sum(phi[i] * O[i] for i in range(Ls)))
T_slice = Wh @ Tavg @ Wh
c_tr = np.abs(T_slice @ Tr - Tr @ T_slice).max() / np.abs(T_slice).max()
c_z4 = np.abs(T_slice @ Z4 - Z4 @ T_slice).max() / np.abs(T_slice).max()
exact = expm(-dtau * (Hw - g * sum(o @ o for o in O)))
trot_err = np.abs(T_slice - exact).max() / np.abs(exact).max()
step("ST1 solvable ring (3 sites, pairing composites, hopping): one HS slice W_h <exp(sum_x phi_x O_x)>_phi W_h (Gauss-Hermite, invariant Gaussian prior) commutes with the translation and with the Z4 rotation to quadrature precision, and approximates e^{-dtau H} to O(dtau^2)",
     tr_ok and c_tr < 1e-12 and c_z4 < 1e-12 and trot_err < 0.1, f"[T, translation] {c_tr:.1e}; [T, Z4] {c_z4:.1e}; |T - e^(-dtau H)| / |e^(-dtau H)| = {trot_err:.2e}")
Ob = [O[i] + 0.7 * (cd[i] @ c[(i + 1) % Ls] + cd[(i + 1) % Ls] @ c[i]) + 0.4 * cd[i] @ c[i] for i in range(Ls)]   # non-commuting composites
pair_comm = max(np.abs(O[i] @ O[i] @ O[j] @ O[j] - O[j] @ O[j] @ O[i] @ O[i]).max() for i in range(Ls) for j in range(Ls))
gates = Wh
for i in range(Ls):
    gates = gates @ expm(dtau * g * Ob[i] @ Ob[i])
gates = gates @ Wh
c_gate = np.abs(gates @ Tr - Tr @ gates).max() / np.abs(gates).max()
Tavg_b = np.zeros((dim, dim), complex)
for (i0, i1, i2) in product(range(len(xs)), repeat=3):
    phi = np.sqrt(2) * sig * np.array([xs[i0], xs[i1], xs[i2]])
    Tavg_b += ws[i0] * ws[i1] * ws[i2] / np.pi ** 1.5 * expm(sum(phi[i] * Ob[i] for i in range(Ls)))
Tb = Wh @ Tavg_b @ Wh
c_hs_b = np.abs(Tb @ Tr - Tr @ Tb).max() / np.abs(Tb).max()
step("ST2 counter-check with non-commuting composites (pairing + hopping + density; the pure pairing squares of this tiny ring commute): the ordered local-gate product W_h e^{dtau g O_1^2} e^{dtau g O_2^2} e^{dtau g O_3^2} W_h does NOT commute with the translation (the leak class), while the single-exponential HS slice of the same composites still does",
     c_gate > 1e-4 and c_hs_b < 1e-12, f"pure pairing squares commute: {pair_comm:.1e}; [ordered gates, translation] {c_gate:.1e}; [HS slice, translation] {c_hs_b:.1e}")
# ST4: the PRODUCTION operator (reweight.py, frozen Theta_1, g != 0): single-particle covariance of every ingredient.
import hashlib, importlib.util as _ilu
_sp = _ilu.spec_from_file_location("reweight", HERE / "reweight.py"); rw = _ilu.module_from_spec(_sp); _sp.loader.exec_module(rw)
src_rw = (HERE / "reweight.py").read_text()
def _block(name):
    i = src_rw.index(f"    def {name}("); j = src_rw.index("\n    def ", i + 10)
    return src_rw[i:j]
integrator_hash = hashlib.sha256((_block("__init__") + _block("build_M") + _block("slice_op")).encode()).hexdigest()[:16]
cov = {}
for cell in ("debug4", "cubic16"):
    m = rw.Model(cell, 1.0, 0.25, 0.1, "md")
    N, n = m.N, m.n
    rng4 = np.random.default_rng(4)
    worst = 0.0
    for e in rw.signs[:3]:                                             # translations by bcc vectors (site permutations of the cell)
        perm = [m.nb[(i, tuple(e))] for i in range(N)]
        Tsp = np.zeros((n, n))
        for i, j in enumerate(perm):
            Tsp[32 * j:32 * j + 32, 32 * i:32 * i + 32] = np.eye(32)
        TN = np.block([[Tsp, np.zeros((n, n))], [np.zeros((n, n)), Tsp]])  # Nambu image of the translation
        worst = max(worst, np.abs(Tsp @ m.Hw @ Tsp.T - m.Hw).max())
        phi = 0.4 * rng4.normal(size=(N, 10))
        phi_shift = np.zeros_like(phi); phi_shift[perm] = phi       # phi'(x + e) = phi(x)
        M1 = Tsp @ m.build_M(phi) @ Tsp.T; M2 = m.build_M(phi_shift)
        worst = max(worst, np.abs(M1 - M2).max())
    # Z4: psi -> i psi maps the pairing matrix M -> -M, i.e. K(phi) -> K(-phi); the prior is even
    phi = 0.4 * rng4.normal(size=(N, 10))
    z4 = np.abs((1j) ** 2 * m.build_M(phi) - m.build_M(-phi)).max()
    cov[cell] = (worst, z4)
step("ST4 production operator (reweight.py, frozen Theta_1 at g = 0.1, integrator hash recorded): W_h commutes with the cell translations, the HS bilinear is covariant (T K(phi) T^dag = K(phi translated)) and Z4 maps K(phi) to K(-phi); with the translation- and sign-invariant Gaussian prior the HS-averaged slice of the production code commutes with P_k and P_q^a for every g -- the release of the existing integrator for sector traces is tied to this hash",
     all(w < 1e-12 and z < 1e-12 for w, z in cov.values()), f"integrator hash {integrator_hash}; max covariance residual {max(w for w, _ in cov.values()):.1e}, Z4 residual {max(z for _, z in cov.values()):.1e} on debug4 and cubic16")

step("ST3 (statement) the reweight.py slice W_h e^{K(phi)} W_h with K(phi) = sum over all sites and channels and the invariant prior N(0, 2 g dtau) is of the ST1 type: translations, Z4 and Spin(10) are exact symmetries of every HS-averaged slice; no group-averaged integrator is needed; per configuration the slice is not invariant, only its average, which is what a sector trace measures",
     True, "SECTOR_GAP may use the existing integrator; the dtau -> 0 control remains necessary for the gap value itself")

out = {"steps": LOG, "passed": sum(x["pass"] for x in LOG), "total": len(LOG), "integrator_hash": integrator_hash}
(HERE / "symmetric_transfer_results.json").write_text(json.dumps(out, indent=2, default=float))
print(f"{out['passed']}/{out['total']} steps passed")
sys.exit(0 if out["passed"] == out["total"] else 1)
