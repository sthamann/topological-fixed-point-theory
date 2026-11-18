#!/usr/bin/env python3
"""Part TW: TWIST.COVARIANCE.01 -- the Spin(10)-Cartan flavour twist is a flat background connection of Theta_1.

Twisted cell: psi(x + R) = exp(2 i (R.delta) Q) psi(x), Q = Sigma_12 = 1/2 - n_0 on the 16 (reweight.py, twist = delta).
The Q = +1/2 flavours then carry momenta k + delta, the Q = -1/2 flavours k - delta.  Six explicit tests:
  TW1 delta -> 0: the twisted code path converges to the untwisted Theta_1 (Hamiltonian and pairing matrix, O(delta)).
  TW2 delta -> -delta: a Spin(10) element U (rotation by pi in the (Gamma_2, Gamma_3) plane, U Q U^dag = -Q) maps the
      twisted cell at delta onto the one at -delta: H_w(-delta) = W H_w(delta) W^dag and M(-delta, phi) = W^T M(delta, R phi) W
      with an orthogonal channel map R -> the two Cartan sectors are exchanged, G_+(delta) = G_-(-delta).
  TW3 sum_A O_A^2 is holonomy covariant: the one-loop mean field of the twisted cell is translation covariant and, for each
      Cartan sector, equals h_W(k +- delta) + f(k +- delta)^2 Sigma^rho_+- at every allowed momentum.
  TW4 no channel is renormalised: ||M_{x,A}||_F is site- and twist-independent for every channel A.
  TW5 a twist by a full superlattice period is trivial (identical H_w and M), i.e. the background is a flat connection.
  TW6 at delta != 0 the sign-exact Pfaffian weight agrees with the stabilised UDV determinant (modulus and phase) and the
      analytic HMC force agrees with central finite differences of the Pfaffian log W (twisted debug4 cell).
Twist runs are for k_star, Delta_min and topology only; symmetry-breaking diagnostics must use untwisted runs (the twist
breaks Spin(10) explicitly to the commutant of Q).
    python3 -B experiments/theory-contracts/seam-holonomy-generating-20261001/twist_covariance.py
"""
from __future__ import annotations

import importlib.util
import json
import sys
from itertools import product
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("reweight", HERE / "reweight.py")
rw = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(rw)
_src = (HERE / "qmc_pfaffian.py").read_text()
exec(compile(_src.split("# ---------------------------------------------------------------------------------------- Pfaffian Gaussian algebra")[1]
             .split("# ---------------------------------------------------------------------------------------- PF1 brute force")[0], "pfaffian_algebra", "exec"))
_mn = (HERE / "mirror_node.py").read_text()
exec(compile("def tau_x(m):" + _mn.split("def tau_x(m):")[1].split("# ---------------------------------------------------------------------------------------- MN1 brute force")[0], "hfb", "exec"))
LOG: list[dict] = []
wrap = lambda a: (a + np.pi) % (2 * np.pi) - np.pi


def step(name: str, ok: bool, detail: str = "") -> None:
    LOG.append({"step": name, "pass": bool(ok), "detail": detail})
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))
    sys.stdout.flush()


rng = np.random.default_rng(17)
S = [rw.sx, rw.sy, rw.sz]
d0 = 0.17 * np.array([1.0, 1.0, 1.0])

# TW1 -------------------------------------------------------------------------------------------------------
m_u = rw.Model("cubic16", 0.5, 0.25, 0.1, "md")
phr = 0.3 * rng.normal(size=(m_u.N, 10))
dev = []
for eps in (1e-3, 1e-5):
    m_e = rw.Model("cubic16", 0.5, 0.25, 0.1, "md", twist=eps * np.ones(3))
    dev.append(max(np.abs(m_e.Hw - m_u.Hw).max(), np.abs(m_e.build_M(phr) - m_u.build_M(phr)).max()))
step("TW1 delta -> 0: the twisted code path converges linearly to the untwisted Theta_1 (Hamiltonian and pairing matrix)",
     dev[1] < 1e-4 and 50 < dev[0] / dev[1] < 200, f"deviation at delta = 1e-3 / 1e-5: {dev[0]:.2e} / {dev[1]:.2e}")

# TW2 -------------------------------------------------------------------------------------------------------
Zq = np.diag([1.0, -1.0]).astype(complex); lowc = np.array([[0, 1], [0, 0]], dtype=complex)
c5 = []
for j in range(5):
    mats = [Zq] * j + [lowc] + [np.eye(2, dtype=complex)] * (4 - j)
    o = mats[0]
    for m_ in mats[1:]:
        o = np.kron(o, m_)
    c5.append(o)
Gm = []
for j in range(5):
    Gm += [c5[j] + c5[j].conj().T, -1j * (c5[j] - c5[j].conj().T)]
even = [b for b in range(32) if bin(b).count("1") % 2 == 0]
Pe = np.eye(32)[:, even]
U16 = Pe.T @ (Gm[1] @ Gm[2]) @ Pe                                       # pi rotation in the (Gamma_2, Gamma_3) plane
Qm = np.diag(rw.QF)
flipQ = np.abs(U16 @ Qm @ U16.conj().T + Qm).max()
unit16 = np.abs(U16 @ U16.conj().T - np.eye(16)).max()
BA = [MA[:16, 16:] for MA in rw.MA]
basis = np.array([B_.ravel() for B_ in BA]).T
Rch = np.array([np.linalg.lstsq(basis, (U16.T @ B_ @ U16).ravel(), rcond=None)[0] for B_ in BA]).real   # U^T B_A U = sum_B R_AB B_B
m_p = rw.Model("cubic16", 0.5, 0.25, 0.1, "md", twist=d0); m_m = rw.Model("cubic16", 0.5, 0.25, 0.1, "md", twist=-d0)
W = np.kron(np.eye(m_p.N), np.kron(np.eye(2), U16))
h_err = np.abs(m_m.Hw - W @ m_p.Hw @ W.conj().T).max()
# psi -> W psi: psi^T M(delta, phi) psi = psi'^T W^-T M W^-1 psi'; with U^T B_A U = R_AB B_B the field maps by R
cands = {"R phi": phr @ Rch, "R^T phi": phr @ Rch.T}
m_err = min(np.abs(m_m.build_M(phr) - W.conj() @ m_p.build_M(v) @ W.conj().T).max() for v in cands.values())
step("TW2 delta -> -delta: a Spin(10) element (pi rotation in the (Gamma_2, Gamma_3) plane, U Q U^dag = -Q) maps the twisted cell at delta onto the one at -delta, for the kinetic term and for the pairing matrix with an orthogonal channel map -> the two Cartan sectors are exchanged, G_+(delta) = G_-(-delta)",
     flipQ < 1e-12 and unit16 < 1e-12 and np.allclose(Rch @ Rch.T, np.eye(10), atol=1e-12) and h_err < 1e-12 and m_err < 1e-12,
     f"U Q U^dag + Q = {flipQ:.1e}; H_w(-delta) vs W H_w(delta) W^dag {h_err:.1e}; pairing matrix {m_err:.1e}; channel map orthogonal")

# TW3 -------------------------------------------------------------------------------------------------------
B3 = rw.B3.astype(float); Binv = np.linalg.inv(B3)
tauv = [sum(Binv[a, i] * S[a] for a in range(3)) for i in range(3)]
def h_of_k(k):
    p = B3 @ k
    return sum((np.sin(p[i]) + np.cos(p[(i + 1) % 3]) - 1) * tauv[i] for i in range(3))
f_of_k = lambda k: 1 - np.cos(k[0]) * np.cos(k[1]) * np.cos(k[2])
Kloc = [np.block([[np.zeros((32, 32)), 2 * MA.conj().T], [2 * MA, np.zeros((32, 32))]]) for MA in rw.MA]
def all_momenta(Rv):
    G_ = 2 * np.pi * np.linalg.inv(Rv.astype(float)).T
    fcc_ = lambda k: np.allclose(np.cos(np.array([k @ np.array(s_) for s_ in product((-1, 1), repeat=3)])), 1, atol=1e-9)
    Nn = int(round(abs(np.linalg.det(Rv.astype(float))) / 4)); ks_ = []
    for mm_ in product(range(-Nn, Nn + 1), repeat=3):
        k = np.array(mm_) @ G_
        if not any(fcc_(k - q) for q in ks_):
            ks_.append(k)
        if len(ks_) == Nn:
            break
    return ks_
beta = 2.0
E_, U_ = np.linalg.eigh(m_p.Hw); G0 = (U_ * (1 / (1 + np.exp(-beta * E_)))) @ U_.conj().T
Sig = []
for i in range(m_p.N):
    Gr = m_p.Ftw[i] @ G0 @ m_p.Ftw[i].conj().T; Gh = m_p.Ftw[i].conj() @ (np.eye(m_p.n) - G0).T @ m_p.Ftw[i].T
    Sig.append(hfb_sigma(Kloc, np.block([[Gr, np.zeros((32, 32))], [np.zeros((32, 32)), Gh]]), 1.0)[:32, :32])
site_dev = max(np.abs(Sg - Sig[0]).max() for Sg in Sig)
gq = 0.02
H1 = m_p.Hw + gq * sum(m_p.Ftw[i].conj().T @ Sig[0] @ m_p.Ftw[i] for i in range(m_p.N))
e_cell = np.sort(np.linalg.eigvalsh(H1))
pred = []
for q in (0.5, -0.5):
    sel = rw.QF == q
    Sq = Sig[0].reshape(2, 16, 2, 16)[:, sel][:, :, :, sel]             # (2, 8, 2, 8) Cartan-sector block
    for k in all_momenta(rw.CELLS["cubic16"]):
        kq = k + 2 * q * d0
        Hk = np.kron(h_of_k(kq), np.eye(8)) + gq * f_of_k(kq) ** 2 * Sq.reshape(16, 16)
        pred += list(np.linalg.eigvalsh(Hk))
spec3 = np.abs(e_cell - np.sort(pred)).max()
step("TW3 sum_A O_A^2 is holonomy covariant: on the twisted cell the one-loop mean field is site-independent and its spectrum equals, Cartan sector by sector, h_W(k +- delta) + f(k +- delta)^2 Sigma^rho_+- at every allowed momentum",
     site_dev < 1e-10 and spec3 < 1e-10, f"site dependence of Sigma^rho {site_dev:.1e}; twisted one-loop spectrum vs momentum-space prediction {spec3:.1e}")

# TW4 -------------------------------------------------------------------------------------------------------
norms = lambda m: np.array([[np.linalg.norm(m.Ftw[i].T @ rw.MA[A] @ m.Ftw[i]) if m.twist is not None else
                             np.linalg.norm(np.kron(m.FF[i], rw.MA[A])) for A in range(10)] for i in range(m.N)])
n_u, n_p = norms(m_u), norms(m_p)
step("TW4 no channel is renormalised by the twist: ||M_{x,A}||_F is the same for every site, every channel and both twists",
     np.abs(n_p - n_u).max() < 1e-12 and np.abs(n_u - n_u[0, 0]).max() < 1e-12, f"max deviation {np.abs(n_p - n_u).max():.1e}; common norm {n_u[0, 0]:.4f}")

# TW5 -------------------------------------------------------------------------------------------------------
dfull = 2 * np.pi * np.linalg.inv(rw.CELLS["cubic16"].astype(float)) @ np.array([1.0, 0.0, 0.0])
m_f = rw.Model("cubic16", 0.5, 0.25, 0.1, "md", twist=dfull + d0); m_d = m_p
full_err = max(np.abs(m_f.Hw - m_d.Hw).max(), np.abs(m_f.build_M(phr) - m_d.build_M(phr)).max())
m_t = rw.Model("cubic16", 0.5, 0.25, 0.1, "md", twist=dfull)
triv_err = max(np.abs(m_t.Hw - m_u.Hw).max(), np.abs(m_t.build_M(phr) - m_u.build_M(phr)).max())
step("TW5 flat connection: a twist by a full superlattice period is trivial (identical H_w and pairing matrix, alone and on top of a twist)",
     full_err < 1e-12 and triv_err < 1e-12, f"period twist alone {triv_err:.1e}; period + delta vs delta {full_err:.1e}")

# TW6 -------------------------------------------------------------------------------------------------------
MAd = np.array([x.conj().T for x in rw.MA])
def dlogW_tw(m, phi):
    """stabilised symmetric HMC force (qmc_pfaffian.py PF5) for a general (twisted) pairing matrix M = sum_i Ftw_i^T (phi_i.MA) Ftw_i."""
    Ss = [m.slice_op(phi[l]) for l in range(m.Lt)]; L = len(Ss); nn = 2 * m.n; I = np.eye(nn, dtype=complex)
    out = np.zeros((L, m.N, 10), complex)
    for l in range(L):
        M = m.build_M(phi[l]); K = np.block([[np.zeros_like(M), 2 * M.conj().T], [2 * M, np.zeros_like(M)]])
        lam, U = np.linalg.eigh(K)
        Eh = (U * np.exp(lam / 2)) @ U.conj().T
        udv = (I, np.ones(nn), I)
        for Mx in [Eh, m.Wh] + [Ss[k] for k in range(l - 1, -1, -1)] + [Ss[k] for k in range(L - 1, l, -1)] + [m.Wh, Eh]:
            udv = m._udv_left(Mx, udv)
        Uq, Dq, Vq = udv
        Db, Ds = np.maximum(Dq, 1.0), np.minimum(Dq, 1.0)
        Vi = np.linalg.inv(Vq)
        Gt = Vi @ np.linalg.solve((Uq.conj().T @ Vi) / Db[:, None] + np.diag(Ds), Uq.conj().T / Db[:, None])
        dl = (lam[:, None] - lam[None, :]) / 2
        with np.errstate(divide="ignore", invalid="ignore"):
            Sg = np.where(np.abs(dl) > 1e-10, np.sinh(dl) / np.where(np.abs(dl) > 1e-10, dl, 1), 1.0)
        Z = U @ ((U.conj().T @ (I - Gt) @ U) * Sg) @ U.conj().T
        n = m.n
        Z12, Z21 = Z[:n, n:], Z[n:, :n]
        for i in range(m.N):
            Y12 = m.Ftw[i] @ Z12 @ m.Ftw[i].T                           # Tr(Z_12 dM), dM = Ftw^T MA Ftw
            Y21 = m.Ftw[i].conj() @ Z21 @ m.Ftw[i].conj().T             # Tr(Z_21 dM^dag), dM^dag = Ftw^dag MA^dag Ftw^*
            out[l, i] = np.einsum("ab,Aba->A", Y12, rw.MA) + np.einsum("ab,Aba->A", Y21, MAd)
    return out
m4 = rw.Model("debug4", 1.0, 0.25, 0.5, "md", twist=0.13 * np.ones(3))
pf_dev = []
for _ in range(6):
    phi = np.sqrt(2 * 0.5 * 0.25) * rng.normal(size=(m4.Lt, m4.N, 10))
    Ss = [m4.slice_op(phi[l]) for l in range(m4.Lt)]
    lp = logW_pf(Ss); la, ph, _ = m4.logdet(Ss)
    pf_dev.append((abs(2 * lp.real - la) / abs(la), abs(wrap(2 * lp.imag - ph)), abs(wrap(lp.imag))))
m16t = rw.Model("cubic16", 1.0, 0.25, 0.1, "md", twist=d0)             # >= 3 slices: genuinely complex weights
phi16 = np.sqrt(2 * 0.1 * 0.25) * rng.normal(size=(m16t.Lt, m16t.N, 10))
S16 = [m16t.slice_op(phi16[l]) for l in range(m16t.Lt)]
lp16 = logW_pf(S16); la16, ph16, _ = m16t.logdet(S16)
pf16 = (abs(2 * lp16.real - la16) / abs(la16), abs(wrap(2 * lp16.imag - ph16)), abs(wrap(lp16.imag)))
phi = np.sqrt(2 * 0.5 * 0.25) * rng.normal(size=(m4.Lt, m4.N, 10))
an = dlogW_tw(m4, phi)
fd_err = []
for t_ in range(10):
    l, i, A = t_ % m4.Lt, (3 * t_) % m4.N, (7 * t_) % 10
    h = 1e-5
    pp = phi.copy(); pp[l, i, A] += h; pm = phi.copy(); pm[l, i, A] -= h
    lpp = logW_pf([m4.slice_op(pp[k]) for k in range(m4.Lt)]); lpm = logW_pf([m4.slice_op(pm[k]) for k in range(m4.Lt)])
    fd = ((lpp.real - lpm.real) + 1j * wrap(lpp.imag - lpm.imag)) / (2 * h)
    fd_err.append(abs(fd - an[l, i, A]) / max(abs(an[l, i, A]), 1e-3))
step("TW6 at delta != 0 (twisted debug4, g = 0.5, and a twisted 16-site configuration with complex weight): the sign-exact Pfaffian weight agrees with the stabilised UDV determinant in modulus and phase, and the analytic HMC force with the twisted pairing matrix agrees with central finite differences of the Pfaffian log W",
     max(a for a, _, _ in pf_dev) < 1e-9 and max(b for _, b, _ in pf_dev) < 1e-8 and max(fd_err) < 1e-5 and pf16[0] < 1e-9 and pf16[1] < 1e-8 and pf16[2] > 1e-6,
     f"debug4: rel |2 log|W| - UDV| {max(a for a, _, _ in pf_dev):.1e}, phase {max(b for _, b, _ in pf_dev):.1e} (weights real there, |arg W| {max(c for _, _, c in pf_dev):.1e}); twisted cubic16 with a genuinely complex weight (|arg W| = {pf16[2]:.2e}): rel {pf16[0]:.1e}, phase {pf16[1]:.1e}; force vs FD {max(fd_err):.1e}")

out = {"steps": LOG, "passed": sum(x["pass"] for x in LOG), "total": len(LOG)}
(HERE / "twist_covariance_results.json").write_text(json.dumps(out, indent=2, default=float))
print(f"{out['passed']}/{out['total']} steps passed")
sys.exit(0 if out["passed"] == out["total"] else 1)
