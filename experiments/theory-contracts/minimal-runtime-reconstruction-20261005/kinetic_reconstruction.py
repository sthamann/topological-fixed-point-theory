#!/usr/bin/env python3
"""SOURCE.KINETIC.RECONSTRUCTION.01 -- which kinetic structure does the source force, and where does Theta_1 sit?

Precondition: the residual quotient of the source carries no Weyl coin (SOURCE.MARK_CLIFFORD.01).  The kinetic
reconstruction therefore runs CONDITIONALLY on a coin C^2 with the source-derived H_1 Gram (AX.WEYLCOIN), with the
mark rotation group T = A4 acting through its spin lift.  Single-chirality Weyl Hamiltonians h(k) = eps(k) + d(k).sigma
are covariant iff d(R k) = R d(k).  Nothing is fitted to Theta_1; Theta_1 is compared only at the end.
  KR1 the points of the bcc Brillouin torus fixed by all of T are exactly Gamma, H, P, -P; T acts irreducibly on R^3,
      so EVERY T-covariant d vanishes at all four: at least four band-touching points -- minimal doubling is impossible
      with the full mark rotation symmetry (also for unitary walks: U(k*) commutes with 2T, hence is scalar).
  KR2 numerical witness: random T-covariant fields (range <= 2) vanish at the four points to machine precision; random
      C3-covariant fields do not vanish at H and -P.
  KR3 subgroup lattice of T: {1, 3 x C2, 4 x C3, V4, T}; a subgroup admits a nonzero invariant vector iff it is 1, C2 or
      C3; V4 and T force d = 0 at Gamma, H, +-P.  Minimal doubling therefore FORCES the breaking T -> H in {C3, C2, 1};
      it is an open condition (KR3b: generic small C3-breaking perturbations of h_W stay minimally doubled with trivial
      residual symmetry).  Only under the extra premise of MAXIMAL residual symmetry is H = C3, the stabiliser of one mark
      (axis || w_a): then the omitted mark of Theta_1 is a point in a 4-element orbit.
  KR4 the nearest-neighbour C3 class (steps +-b_i and the axis hop +-(1,1,1)): an isotropic unit cone at Gamma and the
      mirror at the T-fixed point P with an isotropic equal-speed cone of opposite chirality (V_Gamma = I fixed, which uses
      up the coin-frame rotation about the axis) leave a TWO-parameter family
      (theta, gamma): beta = rotation about [111] by theta off-axis, axis-hop amplitudes (gamma, -gamma); KR4b gamma != 0
      members are minimally doubled too; 'no hop along the selected mark axis' (h_W) sets gamma = 0, leaving the circle
      d_i = sin p_i + sum_j beta_ij(theta) (cos p_j - 1).
  KR5 node census on the circle (5-degree grid + Newton, heuristic root finding; boundary by bisection): theta = 0 gives
      eight nodes; minimal doubling (exactly Gamma(+) and P(-)) holds on the arc (theta_c, 360 - theta_c), theta_c ~ 34.06
      degrees; Theta_1's h_W is exactly theta = 120 degrees.
  KR6 selection list for Theta_1 inside the class: (i) isotropic unit cone at Gamma [IR], (ii) mirror at P, (iii)
      isotropic equal-speed mirror of opposite chirality, (iv) mirror frame differs by a mark rotation (theta in C3),
      (v) minimal doubling (excludes theta = 0), (vi) orientation picks +120 vs -120, (vii) one of four marks, (viii) no hop
      along the selected mark axis (gamma = 0).  (ii)-(vi) and (viii) are invisible to IR tests: members of the arc have the
      same Gamma cone but distinguishable mirror sectors (spin-resolved return amplitude near P differs at O(1)).
  KR7 (statement) SOURCE.KINETIC.RECONSTRUCTION.01 = NOT_DERIVED: no coin in the source quotient; conditional on
      AX.WEYLCOIN the T-symmetric source cannot output a minimally doubled chiral kinetic without breaking T (one mark
      under the maximal-residual-symmetry premise);
      Theta_1 is the unique member of the declared class meeting (i)-(viii), most of them fundamental, not IR, conditions.
    python3 -B experiments/theory-contracts/minimal-runtime-reconstruction-20261005/kinetic_reconstruction.py
"""
from __future__ import annotations

import sys
from itertools import combinations, permutations, product
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.optimize import fsolve

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from ur_lib import Log  # noqa: E402

log = Log()
rng = np.random.default_rng(7)
sig = [np.array([[0, 1], [1, 0]], dtype=complex), np.array([[0, -1j], [1j, 0]]), np.diag([1, -1]).astype(complex)]

Td = []
for perm in permutations(range(3)):
    for signs in product((-1, 1), repeat=3):
        if np.prod(signs) == 1:
            R = np.zeros((3, 3))
            for i, j in enumerate(perm):
                R[i, j] = signs[i]
            Td.append(R)
T = [R for R in Td if np.linalg.det(R) > 0]


def equiv(n1, n2):
    """k = (pi/2) n; equivalence modulo the fcc reciprocal lattice {pi m : sum m even}."""
    d = np.array(n1) - np.array(n2)
    return bool(np.all(d % 2 == 0) and (d // 2).sum() % 2 == 0)


reps = []
for n in product(range(4), repeat=3):
    if not any(equiv(n, r) for r in reps):
        reps.append(n)
fixed = [n for n in reps if all(equiv((R @ np.array(n)).astype(int), n) for R in T)]
names = {}
for n in fixed:
    k = np.pi / 2 * np.array(n)
    if equiv(n, (0, 0, 0)): names[n] = "Gamma"
    elif equiv(n, (2, 0, 0)): names[n] = "H"
    elif equiv(n, (1, 1, 1)): names[n] = "P"
    elif equiv(n, (3, 3, 3)): names[n] = "-P"
inv_vec = np.linalg.matrix_rank(np.vstack([R - np.eye(3) for R in T]))
log.step("KR1 T-fixed points of the bcc Brillouin torus are exactly Gamma, H, P, -P; T has no invariant vector (irreducible on R^3) "
         "-> every T-covariant d vanishes there: >= 4 band touchings, minimal doubling impossible with full mark rotation symmetry",
         len(T) == 12 and sorted(names.values()) == sorted(["Gamma", "H", "P", "-P"]) and len(fixed) == 4 and inv_vec == 3,
         f"{len(reps)} classes of (pi/2)Z^3 mod fcc; fixed: {sorted(names.values())}")

# KR2 random covariant fields
steps2 = [np.array(e) for e in product((-1, 1), repeat=3)] + [2 * np.eye(3)[i] * s for i in range(3) for s in (-1, 1)]


def rand_field():
    a = rng.normal(size=(len(steps2), 3)); b = rng.normal(size=(len(steps2), 3))
    return lambda k: sum(a[j] * np.cos(e @ k) + b[j] * np.sin(e @ k) for j, e in enumerate(steps2))


def symmetrise(d0, group):
    return lambda k: sum(R @ d0(R.T @ k) for R in group) / len(group)


C3 = [np.eye(3), np.roll(np.eye(3), 1, 0), np.roll(np.eye(3), 2, 0)]
pts = {nm: np.pi / 2 * np.array(n) for n, nm in names.items()}
tmax, c3min = 0.0, np.inf
for _ in range(5):
    dT = symmetrise(rand_field(), T)
    dC = symmetrise(rand_field(), C3)
    tmax = max(tmax, max(np.linalg.norm(dT(k)) for k in pts.values()))
    c3min = min(c3min, min(np.linalg.norm(dC(pts[nm])) for nm in ("H", "-P")))
cov = max(np.linalg.norm(dT(R @ k) - R @ dT(k)) for R in T for k in rng.normal(size=(4, 3)))
log.step("KR2 random T-covariant fields (range <= 2) vanish at Gamma, H, P, -P; random C3-covariant fields do not vanish at H, -P",
         tmax < 1e-12 and c3min > 1e-3 and cov < 1e-12, f"max |d_T| at fixed points {tmax:.1e}; min |d_C3| at H,-P {c3min:.2f}")

# KR3 subgroup lattice
def key(R): return tuple(R.astype(int).ravel())
def gen(gs):
    el = {key(np.eye(3)): np.eye(3)}
    fr = [np.eye(3)]
    while fr:
        nf = []
        for g in fr:
            for h in gs:
                x = h @ g
                if key(x) not in el:
                    el[key(x)] = x; nf.append(x)
        fr = nf
    return el
subs = {}
for r in range(0, 3):
    for gs in combinations(T, r):
        el = gen(list(gs)) if gs else {key(np.eye(3)): np.eye(3)}
        subs[frozenset(el)] = el
summary = []
allowed_orders = set()
c3_axes = []
w = np.array([[-.5, -.5, -.5], [-.5, .5, .5], [.5, -.5, .5], [.5, .5, -.5]])
for s, el in subs.items():
    mats = list(el.values())
    inv = 3 - np.linalg.matrix_rank(np.vstack([R - np.eye(3) for R in mats]))
    summary.append((len(mats), int(inv)))
    if inv > 0:
        allowed_orders.add(len(mats))
    if len(mats) == 3:
        R = [m for m in mats if not np.allclose(m, np.eye(3))][0]
        ev, vv = np.linalg.eig(R)
        ax = np.real(vv[:, np.argmin(abs(ev - 1))])
        c3_axes.append(int(np.argmax([abs(ax @ wa) / np.linalg.norm(wa) for wa in w])))
orders = sorted(summary)
log.step("KR3 subgroups of T: 1, three C2, four C3, V4, T; nonzero invariant vector iff 1, C2, C3 -> the largest subgroup compatible "
         "with minimal doubling is C3 (the stabiliser of one mark): minimal doubling forces T-breaking, a single mark only under "
         "maximal residual symmetry",
         orders == sorted([(1, 3), (2, 1), (2, 1), (2, 1), (3, 1), (3, 1), (3, 1), (3, 1), (4, 0), (12, 0)])
         and allowed_orders == {1, 2, 3} and sorted(c3_axes) == [0, 1, 2, 3],
         f"(order, dim invariant vectors) = {orders}; C3 axes along marks {sorted(c3_axes)}")

# KR4 the C3 nearest-neighbour class
B = np.array([[1, 1, -1], [1, -1, 1], [-1, 1, 1]], dtype=float)
Binv = np.linalg.inv(B)
Cp = np.roll(np.eye(3), 1, axis=1)                    # (Cp)_{i,i+1} = 1


def circ(c):
    return c[0] * np.eye(3) + c[1] * Cp + c[2] * Cp @ Cp


def d_p(p, beta, g1=0.0, g2=0.0, alpha=None):
    """d_i = sum_j alpha_ij sin p_j + sum_j beta_ij (cos p_j - 1) + [g2 sin S + g1 (cos S - 1)], S = p1 + p2 + p3."""
    alpha = np.eye(3) if alpha is None else alpha
    S = p.sum()
    return alpha @ np.sin(p) + beta @ (np.cos(p) - 1) + (g2 * np.sin(S) + g1 * (np.cos(S) - 1)) * np.ones(3)


def jac_p(p, beta, g1=0.0, g2=0.0, alpha=None):
    alpha = np.eye(3) if alpha is None else alpha
    S = p.sum()
    return alpha * np.cos(p)[None, :] - beta * np.sin(p)[None, :] + (g2 * np.cos(S) - g1 * np.sin(S)) * np.ones((3, 3))


def member(th, gam):
    """Isotropic unit Gamma cone (alpha = I - gam J), mirror at P (sum beta = 1 - 4 gam - g1), g1 = -gam, beta-eigenvalues e^{+-i th}."""
    s_b = 1 - 3 * gam
    beta = circ([(s_b + 2 * np.cos(th)) / 3, (s_b - np.cos(th) + np.sqrt(3) * np.sin(th)) / 3,
                 (s_b - np.cos(th) - np.sqrt(3) * np.sin(th)) / 3])
    return dict(beta=beta, g1=-gam, g2=gam, alpha=np.eye(3) - gam * np.ones((3, 3)))


Pp = np.pi / 2 * np.ones(3)
Jm = np.ones((3, 3))
ok4 = np.allclose(B @ np.ones(3), np.ones(3)) and np.allclose(Binv @ Cp @ B, Cp.T) and np.allclose(Binv @ Jm @ B, Jm)
# Analytic reduction: V_Gamma = alpha^T + g2 J = I; d(P) = 0 <=> sum beta = 1 - 4 g2 - g1; V_P = -beta^T + g1 J with axis eigenvalue
# -1 + 4 (g1 + g2): orthogonality <=> g1 + g2 in {0, 1/2}; opposite mirror chirality (det V_P < 0) <=> g1 + g2 = 0.
fam_ok, chir = True, []
for th, gam in ((0.9, 0.0), (2.0, 0.13), (4.1, -0.2)):
    m = member(th, gam)
    VG = Binv @ jac_p(np.zeros(3), **m) @ B
    VP = Binv @ jac_p(Pp, **m) @ B
    fam_ok &= np.allclose(VG, np.eye(3)) and np.allclose(d_p(Pp, **m), 0) and np.allclose(VP.T @ VP, np.eye(3))
    chir.append(int(np.sign(np.linalg.det(VP))))
m_bad = member(1.3, 0.1); m_bad["g1"] = 0.5 - 0.1                     # g1 + g2 = 1/2 branch, sum beta adjusted for d(P) = 0
m_bad["beta"] = m_bad["beta"] + ((1 - 4 * 0.1 - m_bad["g1"]) - np.trace(m_bad["beta"] @ np.ones((3, 3))) / 3) / 3 * Jm
VPb = Binv @ jac_p(Pp, **m_bad) @ B
bad = (bool(np.allclose(d_p(Pp, **m_bad), 0)), bool(np.allclose(VPb.T @ VPb, np.eye(3))), int(np.sign(np.linalg.det(VPb))))
log.step("KR4 nearest-neighbour C3 class (all eight bcc steps): unit isotropic cone at Gamma, mirror at P with isotropic equal-speed "
         "cone and opposite chirality is the TWO-parameter family (theta, gamma): beta has eigenvalues e^{+-i theta} off-axis, the "
         "axis-hop amplitudes are (g2, g1) = (gamma, -gamma), alpha = I - gamma J; the branch g1 + g2 = 1/2 has mirror chirality +1",
         ok4 and fam_ok and chir == [-1, -1, -1] and bad == (True, True, 1),
         f"samples (theta, gamma) = (0.9,0), (2.0,0.13), (4.1,-0.2): chirality {chir}; g1+g2 = 1/2 branch (d(P)=0, orthogonal, chirality) = {bad}")


def census(beta, n=28, g1=0.0, g2=0.0, alpha=None):
    grid = np.linspace(0, 2 * np.pi, n, endpoint=False)
    G = np.array(np.meshgrid(grid, grid, grid, indexing="ij")).reshape(3, -1).T
    kw = dict(g1=g1, g2=g2, alpha=alpha)
    val = np.array([np.linalg.norm(d_p(p, beta, **kw)) for p in G])
    cand = G[np.argsort(val)[:400]]
    roots = []
    for p0 in cand:
        p, info, ier, _ = fsolve(lambda p: d_p(p, beta, **kw), p0, fprime=lambda p: jac_p(p, beta, **kw), full_output=True, xtol=1e-13)
        if ier != 1 or np.linalg.norm(d_p(p, beta, **kw)) > 1e-10:
            continue
        p = np.mod(p, 2 * np.pi)
        if not any(np.linalg.norm(np.angle(np.exp(1j * (p - q)))) < 1e-6 for q, _ in roots):
            roots.append((p, int(np.sign(np.linalg.det(jac_p(p, beta, **kw))))))
    return roots


gam_census = {}
for gam in (-0.1, -0.05, 0.05, 0.1):
    r = census(**member(np.radians(120), gam))
    gam_census[gam] = (len(r), sum(c for _, c in r))
log.step("KR4b the gamma direction is a genuine family of minimally doubled kinetics (theta = 120 deg, gamma = +-0.05, +-0.1: two "
         "nodes, chirality sum 0); it is removed only by the condition (viii) 'no hop along the selected mark axis' (g1 = g2 = 0), "
         "which is what h_W does", all(v == (2, 0) for v in gam_census.values()), f"(nodes, chirality sum) {gam_census}")


def beta_theta(th):
    return member(th, 0.0)["beta"]


res5 = {}
for deg in range(0, 360, 5):
    r = census(beta_theta(np.radians(deg)))
    res5[deg] = (len(r), sum(c for _, c in r))
theta1 = circ([0, 1, 0])
same = np.allclose(beta_theta(np.radians(120)), theta1, atol=1e-12)
r120 = census(theta1)
names120 = sorted(("Gamma" if np.allclose(np.angle(np.exp(1j * p)), 0, atol=1e-6) else
                   "P" if np.allclose(np.angle(np.exp(1j * (p - Pp))), 0, atol=1e-6) else "other", c) for p, c in r120)
arc = [deg for deg, (n_, _) in res5.items() if n_ == 2]
log.step("KR5 node census on the circle: theta = 0 has eight nodes; exactly two nodes (Gamma +1, P -1) on one contiguous arc "
         "containing +-120 degrees; Theta_1's h_W (d_i = sin p_i + cos p_{i+1} - 1) is theta = 120 degrees",
         res5[0][0] == 8 and same and names120 == [("Gamma", 1), ("P", -1)] and 120 in arc and 240 in arc
         and all(s == 0 for _, s in res5.values()) and len(arc) == (max(arc) - min(arc)) // 5 + 1,
         f"theta=0: {res5[0][0]} nodes; two nodes for theta in [{min(arc)}, {max(arc)}] deg (5-degree census, chirality sums 0); "
         f"node counts outside the arc: {sorted(set(n_ for d_, (n_, _) in res5.items() if d_ not in arc))}")

# KR3b openness: small generic C3-breaking nearest-neighbour perturbations of h_W stay minimally doubled
C3cart = np.roll(np.eye(3), 1, axis=0)
open_res = []
for trial in range(3):
    pert = dict(beta=circ([0, 1, 0]) + 0.04 * rng.normal(size=(3, 3)), alpha=np.eye(3) + 0.04 * rng.normal(size=(3, 3)),
                g1=0.04 * rng.normal(), g2=0.04 * rng.normal())
    rr = census(**pert)
    dc = lambda k: Binv @ d_p(B @ k, **pert)
    viol = max(np.linalg.norm(dc(C3cart @ k) - C3cart @ dc(k)) for k in rng.normal(size=(6, 3)))
    open_res.append((len(rr), sum(c for _, c in rr), round(float(viol), 3)))
log.step("KR3b minimal doubling is open: three random C3-breaking nearest-neighbour perturbations (size 0.04) of h_W keep exactly two "
         "nodes with chirality sum 0 while violating C3 -> minimal doubling forces only T-breaking, not C3",
         all(n_ == 2 and c_ == 0 and v_ > 1e-3 for n_, c_, v_ in open_res), f"(nodes, chirality sum, C3 violation) = {open_res}")


def n_nodes(deg, n=36):
    return len(census(beta_theta(np.radians(deg)), n=n))


lo, hi = 30.0, 35.0
if n_nodes(lo) != 2 and n_nodes(hi) == 2:
    for _ in range(12):
        mid = (lo + hi) / 2
        lo, hi = (lo, mid) if n_nodes(mid) == 2 else (mid, hi)
theta_c = round((lo + hi) / 2, 2)
log.step("KR5b arc boundary by bisection (finer grid): minimal doubling for theta in (theta_c, 360 - theta_c)", 33.5 < theta_c < 34.6,
         f"theta_c ~ {theta_c} deg")

# KR6 IR equivalence, fundamental distinguishability along the arc
def h_k(k, beta):
    d = Binv @ d_p(B @ k, beta)
    return sum(di * s for di, s in zip(d, sig))
k_ir = 1e-3 * np.array([0.3, -0.2, 0.5])
arc_pts = [np.radians(a) for a in (60, 120, 180)]
ir_dev = max(np.abs(np.linalg.eigvalsh(h_k(k_ir, beta_theta(a))) - np.linalg.eigvalsh(h_k(k_ir, beta_theta(arc_pts[1])))).max()
             for a in arc_pts)
cone = max(abs(np.linalg.eigvalsh(h_k(k_ir, beta_theta(a)))[1] - np.linalg.norm(k_ir)) for a in arc_pts)
q = 0.05 * np.array([1.0, 0.0, 0.0])
up = np.array([1, 0], dtype=complex)
t_q = np.pi / (2 * np.linalg.norm(q))
amp = [abs(up.conj() @ expm(-1j * t_q * h_k(np.pi / 2 * np.ones(3) + q, beta_theta(a))) @ up) ** 2 for a in arc_pts]
log.step("KR6 IR versus fundamental along the arc (theta = 60, 120, 180): identical Gamma cone (|E - |k|| = O(k^2), members differ by <= 1e-6 at |k| = 1e-3), "
         "but the spin-resolved return probability near the mirror differs at O(1): conditions (ii)-(vi), (viii) are fundamental, not IR",
         ir_dev < 1e-5 and cone < 1e-5 and max(amp) - min(amp) > 0.2,
         f"IR spread {ir_dev:.1e}; mirror return probabilities {np.round(amp, 3)}")

log.step("KR7 (statement) SOURCE.KINETIC.RECONSTRUCTION.01 = NOT_DERIVED: no coin in the source quotient; conditional on AX.WEYLCOIN a "
         "T-symmetric source cannot give minimal doubling without breaking T (KR1/KR3; one mark under maximal residual symmetry); "
         "Theta_1 = unique member of the "
         "declared nearest-neighbour C3 class satisfying (i)-(viii)", True)

log.dump(HERE / "kinetic_reconstruction_results.json", {
    "research_id": "SOURCE.KINETIC.RECONSTRUCTION.01",
    "verdict": "NOT_DERIVED",
    "pins": {"T_fixed_points": sorted(names.values()), "census_theta": {str(k): v for k, v in res5.items()},
             "minimal_doubling_arc_samples_deg": arc, "theta_c_deg": theta_c, "gamma_census_theta120": {str(k): v for k, v in gam_census.items()}},
})
