#!/usr/bin/env python3
"""Part GD: GAUGE.DYNAMIC.U1.01 -- dynamical U(1) links on the smallest bcc plaquette (separate from Gate B).

Rhombus plaquette v0 = 0, v1 = eps1, v2 = eps1 + eps2, v3 = eps2 (eps1 = (1,1,1), eps2 = (1,-1,-1)); a two-component Weyl
spinor per vertex with the walk hopping blocks h_eps (Part HS); link Hilbert spaces truncated to E = -E_max..E_max with
the raising operator U; H = sum_links (psi_y^dag h_eps U_l psi_x + h.c.) + (g^2/2) sum E_l^2 - (1/(2 g^2)) (U_box + h.c.).
Gauss operators G_x = (div E)_x + n_x - q_x.
  GD1 exact Gauss law: [H, G_x] = 0 and [G_x, G_y] = 0 for every vertex (truncation does not spoil it).
  GD2 process invariance: EVERY primitive gate (each link hopping, each electric term, the plaquette term) commutes with
      every G_x, so each single gate exp(-i tau h_gate) keeps a physical state physical (leakage 0), not only the total H.
  GD3 physical subspace = common kernel (diagonal in the occupation x electric basis); static charges +-1 at adjacent
      vs opposite corners: at strong coupling the energy difference approaches g^2/2 per extra flux link (string).
  GD4 truncation: E_max = 1 and 2 give the same strong-coupling static potential to O(1/g^4).
  GD5 (statement) Gauss law and per-gate invariance are exact for dynamical U(1) links; NOT shown: massless photons in
      the infinite-volume limit, SU(2)/SU(3), the joint test with the filtered interaction (Gate B + C).
    python3 -B experiments/theory-contracts/seam-holonomy-generating-20261001/gauge_dynamic.py
"""
from __future__ import annotations

import json
import sys
from itertools import product
from pathlib import Path

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh, expm_multiply

HERE = Path(__file__).resolve().parent
LOG: list[dict] = []


def step(name: str, ok: bool, detail: str = "") -> None:
    LOG.append({"step": name, "pass": bool(ok), "detail": detail})
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))
    sys.stdout.flush()


sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.diag([1.0 + 0j, -1]); I2 = np.eye(2)
def A_eps(e):
    return (I2 + e[0] * sx) / 2 @ (I2 + e[1] * sy) / 2 @ (I2 + e[2] * sz) / 2
def h_of(e):
    e = np.array(e)
    return 0.5j * (A_eps(e) - A_eps(-e).conj().T)
e1, e2 = (1, 1, 1), (1, -1, -1)
links = [(0, 1, e1), (1, 2, e2), (3, 2, e1), (0, 3, e2)]          # (from, to, direction)
nv, nmode = 4, 8                                                  # 4 vertices x 2 spinor components


def build(Emax, g2, t=1.0):
    d = 2 * Emax + 1
    Zp = sp.csr_matrix(np.diag([1.0, -1.0])); Ii = sp.identity(2, format="csr")
    low = sp.csr_matrix(np.array([[0, 1], [0, 0]], dtype=complex))
    def fop(i):
        out = sp.identity(1, format="csr")
        for j in range(nmode):
            out = sp.kron(out, Zp if j < i else (low if j == i else Ii), format="csr")
        return out
    c = [fop(i) for i in range(nmode)]
    nF = 2 ** nmode
    Ed = sp.diags(np.arange(-Emax, Emax + 1).astype(float)).tocsr()
    Ur = sp.diags(np.ones(d - 1), -1, shape=(d, d)).tocsr().astype(complex)       # |E> -> |E+1>, truncated
    Il = sp.identity(d, format="csr")
    def emb_link(op, l):
        mats = [Il] * 4
        mats[l] = op
        out = mats[0]
        for m_ in mats[1:]:
            out = sp.kron(out, m_, format="csr")
        return sp.kron(sp.identity(nF, format="csr"), out, format="csr")
    def emb_f(op):
        return sp.kron(op, sp.identity(d ** 4, format="csr"), format="csr")
    E = [emb_link(Ed, l) for l in range(4)]
    U = [emb_link(Ur, l) for l in range(4)]
    cf = [emb_f(x) for x in c]
    n_v = [cf[2 * v].conj().T @ cf[2 * v] + cf[2 * v + 1].conj().T @ cf[2 * v + 1] for v in range(nv)]
    hops = []
    for l, (x, y, e) in enumerate(links):
        he = h_of(e)
        T = sum(he[a, b] * cf[2 * y + a].conj().T @ U[l] @ cf[2 * x + b] for a in range(2) for b in range(2))
        hops.append(t * (T + T.conj().T))
    elec = [0.5 * g2 * E[l] @ E[l] for l in range(4)]
    Ubox = U[0] @ U[1] @ U[2].conj().T @ U[3].conj().T
    plaq = -(0.5 / g2) * (Ubox + Ubox.conj().T)
    H = sum(hops) + sum(elec) + plaq
    divE = [sum(E[l] for l, (x, y, e) in enumerate(links) if x == v) - sum(E[l] for l, (x, y, e) in enumerate(links) if y == v) for v in range(nv)]
    return dict(H=H, hops=hops, elec=elec, plaq=plaq, divE=divE, n_v=n_v, E=E, dim=nF * d ** 4)


def comm_norm(A, B):
    C = A @ B - B @ A
    return 0.0 if C.nnz == 0 else float(np.abs(C.data).max())


M = build(1, 2.0)
G = [M["divE"][v] + M["n_v"][v] for v in range(nv)]
c1 = max(comm_norm(M["H"], Gx) for Gx in G)
c2 = max(comm_norm(G[a], G[b]) for a in range(nv) for b in range(nv))
step("GD1.exact Gauss law with dynamical truncated U(1) links (E_max = 1) and two-component Weyl hoppings: [H, G_x] = 0 and [G_x, G_y] = 0 for every vertex",
     c1 < 1e-12 and c2 < 1e-12, f"max |[H, G_x]| = {c1:.1e}, max |[G_x, G_y]| = {c2:.1e}; Hilbert dimension {M['dim']}")

gates = M["hops"] + M["elec"] + [M["plaq"]]
cg = max(comm_norm(gt, Gx) for gt in gates for Gx in G)
qbg = np.array([1, 0, 0, 0])                                         # one static background charge: one dynamical fermion + flux strings
diagG = np.array([Gx.diagonal().real for Gx in G]) - qbg[:, None]
phys = np.where(np.all(np.abs(diagG) < 1e-12, axis=0))[0]
rng = np.random.default_rng(71)
psi = np.zeros(M["dim"], complex)
psi[phys] = rng.normal(size=len(phys)) + 1j * rng.normal(size=len(phys)); psi /= np.linalg.norm(psi)
leak = 0.0
mask = np.ones(M["dim"], bool); mask[phys] = False
for gt in gates:
    out = expm_multiply(-1j * 0.37 * gt, psi)
    leak = max(leak, np.linalg.norm(out[mask]))
step("GD2.process invariance: every primitive gate (four link hoppings, four electric terms, the plaquette term) commutes with every G_x, so each single gate keeps a physical state physical",
     cg < 1e-12 and leak < 1e-12 and len(phys) > 3, f"max |[gate, G_x]| = {cg:.1e}; sector with one static charge (one fermion + flux strings), physical dimension {len(phys)}; max leakage under single gates {leak:.1e}")


def static_energy(Emax, g2, qcharges):
    Mq = build(Emax, g2)
    Gq = np.array([(Mq["divE"][v] + Mq["n_v"][v]).diagonal().real - qcharges[v] for v in range(nv)])
    ph = np.where(np.all(np.abs(Gq) < 1e-12, axis=0))[0]
    Hp = Mq["H"][ph][:, ph]
    if len(ph) < 3:
        return np.min(np.linalg.eigvalsh(Hp.toarray()))
    return eigsh(Hp, k=1, which="SA")[0][0]
pot = {}
for Emax in (1, 2):
    for g2 in (4.0, 8.0):
        adj = static_energy(Emax, g2, [1, -1, 0, 0])
        opp = static_energy(Emax, g2, [1, 0, -1, 0])
        vac = static_energy(Emax, g2, [0, 0, 0, 0])
        pot[(Emax, g2)] = (opp - vac, adj - vac)
ratios = {k: (v[0] - v[1]) / (k[1] / 2) for k, v in pot.items()}
step("GD3.static charges +-1 at adjacent vs opposite corners (fermion-free physical sector): at strong coupling the extra flux link costs ~ g^2/2 (string), the ratio tends to 1 as g^2 grows",
     all(v[0] > v[1] for v in pot.values()) and abs(ratios[(2, 8.0)] - 1) < abs(ratios[(2, 4.0)] - 1) + 1e-9 and abs(ratios[(2, 8.0)] - 1) < 0.2,
     "; ".join(f"E_max={k[0]}, g^2={k[1]}: (opp, adj) = ({v[0]:.3f}, {v[1]:.3f}), ratio {ratios[k]:.3f}" for k, v in pot.items()))
step("GD4.truncation: E_max = 1 and 2 give the same strong-coupling static potential within O(1/g^4)",
     abs(pot[(1, 8.0)][0] - pot[(2, 8.0)][0]) < 0.05 * pot[(2, 8.0)][0] and abs(pot[(1, 8.0)][1] - pot[(2, 8.0)][1]) < 0.05 * pot[(2, 8.0)][1],
     f"g^2 = 8: opposite {pot[(1, 8.0)][0]:.4f} vs {pot[(2, 8.0)][0]:.4f}, adjacent {pot[(1, 8.0)][1]:.4f} vs {pot[(2, 8.0)][1]:.4f}")
step("GD5.(statement) Gauss law and per-gate invariance are exact for dynamical U(1) links on the bcc plaquette; NOT shown: massless photons in the infinite-volume limit, SU(2)/SU(3), the joint Gate B + C test",
     True, "Gate C: first exact checkpoint only")

out = {"steps": LOG, "passed": sum(x["pass"] for x in LOG), "total": len(LOG), "static_potential": {str(k): v for k, v in pot.items()}}
(HERE / "gauge_dynamic_results.json").write_text(json.dumps(out, indent=2, default=float))
print(f"{out['passed']}/{out['total']} steps passed")
sys.exit(0 if out["passed"] == out["total"] else 1)
