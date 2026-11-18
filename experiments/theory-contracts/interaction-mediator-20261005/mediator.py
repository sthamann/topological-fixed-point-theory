#!/usr/bin/env python3
"""SOURCE.INTERACTION.MEDIATOR.01 -- can the Theta_1 vertex -g sum_A O_A^2 come from eliminating a heavy source mode?

  MD1 elimination mechanism (exact check): H = H_0 + Delta sum_A b_A^dag b_A + lambda sum_A (b_A + b_A^dag) O_A gives, in the
      mediator vacuum, H_eff = H_0 - (lambda^2/Delta) sum_A O_A^2 + O(1/Delta): for a single channel with H_0 = 0 the
      completed square is exact; with noncommuting channels and hopping the error falls like 1/Delta at fixed
      g = lambda^2/Delta (toy fermion system, truncated oscillators, truncation converged).
  MD2 mediator dictionary on V = 16 (x) 2: the invariant symmetric bilinear form B(u, v) = <z, u ^ v> on Lambda^2 V (z the
      unique Delta N = 4 invariant) restricted to each channel decides whether a REAL mediator in that channel also produces the
      Delta N = 4 term: B is nondegenerate on (10,1) AND on (120,3), zero on (126,1).  A complex (charge-2) mediator always gives
      only the U(1)-invariant |P_R|^2.
  MD3 consequence: Theta_1's ray sum_A O_A^2 = sum_A (P_A + P_A^dag)^2 is exactly ONE real (10,1) mediator: a heavy real
      Spin(10) vector, SU(2)-spin singlet, Z4-odd (Phi -> -Phi under psi -> i psi), coupled on site to the filtered pair
      operator.  The QMC's Hubbard-Stratonovich field has the same quantum numbers but is a mathematical rewriting of the
      quartic, not evidence of a physical mode.
  MD4 source inventory: E8 roots split under D5 x D3 into (45,1) + (1,15) + (10,6) + (16,4) + (16b,4b) (40+5, 12+3, 60, 128); the
      only bosonic Spin(10)-vector currents are the 60 (10,6) roots, which carry a family SO(6) vector (no singlet under
      SO(6), SU(4) or the family SO(3)) and are weight-1 (massless) seam currents: no elementary heavy family-singlet real 10.
  MD5 (statement) SOURCE.INTERACTION.MEDIATOR.01 = MEDIATOR_SPECIFIED_NOT_IN_ELEMENTARY_SOURCE.
    python3 -B experiments/theory-contracts/interaction-mediator-20261005/mediator.py
"""
from __future__ import annotations

import sys
from itertools import combinations, product
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "minimal-runtime-reconstruction-20261005"))
from ur_lib import Log, jordan_wigner  # noqa: E402

log = Log()
rng = np.random.default_rng(17)

# ---------------------------------------------------------------- MD1 elimination mechanism
nf = 4
a = jordan_wigner(nf)
ad = [x.conj().T for x in a]


def pair_op(M):
    P = sum(M[i, j] * a[i] @ a[j] for i in range(nf) for j in range(nf))
    return P + P.conj().T


Ms = []
for _ in range(2):
    m = rng.normal(size=(nf, nf)); Ms.append(m - m.T)
O = [pair_op(m) for m in Ms]
O = [o / np.linalg.norm(o, 2) for o in O]                                    # unit operator norm
hop = rng.normal(size=(nf, nf)); hop = 0.1 * (hop + hop.T)
H0 = sum(hop[i, j] * ad[i] @ a[j] for i in range(nf) for j in range(nf))


def boson(nmax):
    b = np.diag(np.sqrt(np.arange(1, nmax)), 1).astype(complex)
    return b


def full_low(H0, Os, lam, Delta, nmax, k=6):
    dF = H0.shape[0]
    b = boson(nmax); Ib = np.eye(nmax)
    nb = len(Os)
    def emb_b(op, A):
        mats = [Ib] * nb; mats[A] = op
        out = np.array([[1.0 + 0j]])
        for m_ in mats:
            out = np.kron(out, m_)
        return out
    Db = nmax ** nb
    H = np.kron(H0, np.eye(Db))
    for A, OA in enumerate(Os):
        H = H + Delta * np.kron(np.eye(dF), emb_b(b.conj().T @ b, A)) + lam * np.kron(OA, emb_b(b + b.conj().T, A))
    return np.linalg.eigvalsh(H)[:k]


g = 0.2
errs, single_exact = [], None
for Delta in (10.0, 20.0, 40.0):
    lam = np.sqrt(g * Delta)
    Heff = H0 - g * sum(o @ o for o in O)
    ev_eff = np.linalg.eigvalsh(Heff)[:6]
    ev_full = full_low(H0, O, lam, Delta, 7)
    errs.append(np.abs(ev_full - ev_eff).max())
conv = abs(full_low(H0, O, np.sqrt(g * 20), 20.0, 9)[0] - full_low(H0, O, np.sqrt(g * 20), 20.0, 7)[0])
ev1 = full_low(np.zeros_like(H0), O[:1], np.sqrt(g * 10), 10.0, 14, k=4)
ev1_eff = np.linalg.eigvalsh(-g * O[0] @ O[0])[:4]
single_exact = np.abs(ev1 - ev1_eff).max()
log.step("MD1 elimination: a heavy oscillator per channel coupled lambda (b + b^dag) O_A gives -(lambda^2/Delta) sum O_A^2; single channel "
         "with H_0 = 0 exact (completed square); with two noncommuting channels plus hopping the error falls ~1/Delta at fixed g",
         single_exact < 1e-6 and errs[0] > errs[1] > errs[2] and 1.6 < errs[0] / errs[1] < 2.6 and 1.6 < errs[1] / errs[2] < 2.6
         and conv < 1e-8 and errs[2] < 1e-2, f"single-channel error {single_exact:.1e}; errors at Delta = 10, 20, 40: {np.round(errs, 4)}; truncation {conv:.1e}")

# ---------------------------------------------------------------- Spin(10) x SU(2) on V = 16 (x) 2 (same construction as quartic_vertex)
X = np.array([[0, 1], [1, 0]], dtype=complex); Y = np.array([[0, -1j], [1j, 0]]); Zm = np.diag([1, -1]).astype(complex)
I2 = np.eye(2, dtype=complex)


def kron(*ms):
    out = np.array([[1.0 + 0j]])
    for m_ in ms:
        out = np.kron(out, m_)
    return out


Gam = []
for j in range(5):
    Gam += [kron(*[Zm] * j, X, *[I2] * (4 - j)), kron(*[Zm] * j, Y, *[I2] * (4 - j))]
chir = kron(*[Zm] * 5)
plus = [i for i in range(32) if chir[i, i].real > 0]
gens10 = [(0.5j * Gam[p] @ Gam[q])[np.ix_(plus, plus)] for p, q in combinations(range(10), 2)]
gens = [np.kron(gg, I2) for gg in gens10] + [np.kron(np.eye(16), s / 2) for s in (X, Y, Zm)]
weights = np.array([[np.real(np.kron((0.5j * Gam[2 * j] @ Gam[2 * j + 1])[np.ix_(plus, plus)], I2)[i, i]) for j in range(5)]
                    + [np.real(np.kron(np.eye(16), Zm / 2)[i, i])] for i in range(32)])
cols_nz = [[np.nonzero(np.abs(gg[:, c]) > 1e-12)[0] for c in range(32)] for gg in gens]


def act(gi, vec):
    gg = gens[gi]; out = {}
    for t, cf in vec.items():
        for slot, c in enumerate(t):
            for m_ in cols_nz[gi][c]:
                if m_ in t and m_ != c:
                    continue
                arr = list(t); arr[slot] = int(m_)
                sgn = 1
                for i in range(len(arr)):
                    for j in range(len(arr) - 1 - i):
                        if arr[j] > arr[j + 1]:
                            arr[j], arr[j + 1] = arr[j + 1], arr[j]; sgn = -sgn
                key = tuple(arr)
                out[key] = out.get(key, 0) + cf * gg[m_, c] * sgn
    return {k: v for k, v in out.items() if abs(v) > 1e-13}


b2 = list(combinations(range(32), 2)); i2 = {t: i for i, t in enumerate(b2)}
G2 = []
for gi in range(len(gens)):
    M = np.zeros((len(b2), len(b2)), dtype=complex)
    for i, t in enumerate(b2):
        for kk, vv in act(gi, {t: 1.0}).items():
            M[i2[kk], i] += vv
    G2.append(M)
Ctot = sum(m_ @ m_ for m_ in G2[:45]) + 1000 * sum(m_ @ m_ for m_ in G2[45:])
w2, v2 = np.linalg.eigh(Ctot)
channels = {"(10,1)": 9.0, "(126,1)": 25.0, "(120,3)": 2021.0}
zero4 = [t for t in combinations(range(32), 4) if np.allclose(weights[list(t)].sum(0), 0)]
iz = {t: i for i, t in enumerate(zero4)}
C4 = np.zeros((len(zero4), len(zero4)), dtype=complex)
for i, t in enumerate(zero4):
    acc = {}
    for gi in range(len(gens)):
        for kk, vv in act(gi, act(gi, {t: 1.0})).items():
            acc[kk] = acc.get(kk, 0) + vv
    for kk, vv in acc.items():
        if kk in iz:
            C4[iz[kk], i] += vv
zvec = np.linalg.eigh((C4 + C4.conj().T) / 2)[1][:, 0]
Z4 = np.array(zero4); I_, J_, K_, L_ = Z4.T


def as_matrix(vec):
    M = np.zeros((32, 32), dtype=complex)
    for (i, j), c in zip(b2, vec):
        M[i, j], M[j, i] = c, -c
    return M


def B_matrix(vecs):
    Vs = np.array([as_matrix(x) for x in vecs])
    T1 = lambda p, q: Vs[:, p, q]                                       # (n_vec, 240)
    terms = [(I_, J_, K_, L_, 1), (I_, K_, J_, L_, -1), (I_, L_, J_, K_, 1), (J_, K_, I_, L_, 1), (J_, L_, I_, K_, -1), (K_, L_, I_, J_, 1)]
    zc = zvec.conj()
    Bm = np.zeros((len(vecs), len(vecs)), dtype=complex)
    for p, q, r, s, sg in terms:
        Bm += sg * (T1(p, q) * zc[None, :]) @ T1(r, s).T
    return Bm


ranks = {}
for name, ev in channels.items():
    vecs = v2[:, np.isclose(w2, ev, atol=1e-6)].T
    Bm = B_matrix(vecs)
    ranks[name] = (vecs.shape[0], int(np.linalg.matrix_rank(Bm, tol=1e-8)), float(np.abs(Bm - Bm.T).max()))
Ball = B_matrix(v2.T)
cross_ok = True
idx = {name: np.where(np.isclose(w2, ev, atol=1e-6))[0] for name, ev in channels.items()}
for n1, n2 in combinations(channels, 2):
    cross_ok &= np.abs(Ball[np.ix_(idx[n1], idx[n2])]).max() < 1e-9
log.step("MD2 mediator dictionary: the invariant bilinear B(u,v) = <z, u ^ v> on Lambda^2 V is symmetric, block-diagonal in the channels, "
         "nondegenerate on (10,1) and on (120,3), zero on (126,1) (complex rep): real 10 and real (120,3) mediators both feed the "
         "Delta N = 4 term, 126 mediators cannot",
         ranks["(10,1)"][1] == 10 and ranks["(120,3)"][1] == 360 and ranks["(126,1)"][1] == 0 and cross_ok and all(r[2] < 1e-9 for r in ranks.values()),
         "; ".join(f"{k}: dim {v[0]}, rank B {v[1]}" for k, v in ranks.items()))
pins = {"B_ranks": {k: v[1] for k, v in ranks.items()}}
log.step("MD3 Theta_1's ray sum_A O_A^2 = sum_A (P_A + P_A^dag)^2 is exactly one REAL (10,1) mediator (|P_10|^2 plus the Delta N = 4 "
         "term through B_(10,1)); a complex (charge-2) 10 gives only sum {P_A, P_A^dag}; a real (120,3) adds |P_(120,3)|^2, so it cannot "
         "reproduce the ray; (126,1) mediators give only |P_126|^2; required "
         "quantum numbers: heavy real Spin(10) vector, spin singlet, Z4-odd, on-site coupling to the filtered pair operator (the QMC "
         "Hubbard-Stratonovich field has these quantum numbers but is only a rewriting, not a physical mode)", ranks["(10,1)"][1] == 10 and ranks["(126,1)"][1] == 0,
         f"real (120,3) mediator feeds Delta N = 4: {ranks['(120,3)'][1] > 0}")

# ---------------------------------------------------------------- MD4 source inventory from E8
roots = []
for i, j in combinations(range(8), 2):
    for si, sj in product((1, -1), repeat=2):
        r = np.zeros(8); r[i], r[j] = si, sj; roots.append(r)
for signs in product((0.5, -0.5), repeat=8):
    if sum(1 for s in signs if s < 0) % 2 == 0:
        roots.append(np.array(signs))
roots = np.array(roots)
car, famc = roots[:, :5], roots[:, 5:]
is_int = np.all(np.abs(roots - np.round(roots)) < 1e-9, axis=1)
blocks = {
    "(45,1) roots": int(np.sum(is_int & (np.abs(famc).sum(1) == 0))),
    "(1,15) roots": int(np.sum(is_int & (np.abs(car).sum(1) == 0))),
    "(10,6)": int(np.sum(is_int & (np.abs(car).sum(1) == 1) & (np.abs(famc).sum(1) == 1))),
    "spinor": int(np.sum(~is_int)),
}
vec10_6 = roots[is_int & (np.abs(car).sum(1) == 1) & (np.abs(famc).sum(1) == 1)]
fam_dirs = sorted({tuple(r[5:]) for r in vec10_6})
w = np.array([[-.5, -.5, -.5], [-.5, .5, .5], [.5, -.5, .5], [.5, .5, -.5]])
Rs = []
for perm in [(1, 2, 0, 3), (1, 0, 3, 2), (2, 0, 1, 3)]:
    P = np.zeros((4, 4)); P[list(perm), range(4)] = 1; Rs.append(w.T @ P @ w)
six = [np.kron(np.eye(2), R) for R in Rs]                                   # family SO(3) on the 6 real family Majoranas
inv_dim = 6 - np.linalg.matrix_rank(np.vstack([S - np.eye(6) for S in six]))
log.step("MD4 source inventory: E8 roots split under D5 x D3 as 40 (+5 Cartan) = (45,1), 12 (+3) = (1,15), 60 = (10,6), 128 spinor = "
         "(16,4)+(16b,4b); the only bosonic Spin(10)-vector currents are (10,6): six family directions (an SO(6) vector) with no "
         "invariant vector under the family SO(3) -> no family-singlet elementary real 10; they are weight-1 (massless) seam currents, "
         "not heavy modes",
         len(roots) == 240 and blocks == {"(45,1) roots": 40, "(1,15) roots": 12, "(10,6)": 60, "spinor": 128} and len(fam_dirs) == 6
         and inv_dim == 0, f"{blocks}; family directions {len(fam_dirs)}; SO(3)-invariant family vectors {inv_dim}")
log.step("MD5 (statement) SOURCE.INTERACTION.MEDIATOR.01 = MEDIATOR_SPECIFIED_NOT_IN_ELEMENTARY_SOURCE: the elimination mechanism is "
         "exact (g = lambda^2/Delta), Theta_1's vertex needs one heavy real (10,1), spin-singlet, Z4-odd mode; the elementary source has "
         "none (only massless family-vector (10,6) currents); composites or an extended source remain", True)

log.dump(HERE / "mediator_results.json", {"research_id": "SOURCE.INTERACTION.MEDIATOR.01",
                                          "verdict": "MEDIATOR_SPECIFIED_NOT_IN_ELEMENTARY_SOURCE", "pins": pins})
