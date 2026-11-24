"""v1036 -- recount of PS.THRESHOLD.02 / PS.PROTON.02 (v266): the SU(4)_c adjoint
(15,1,1) from the single REAL 45 of the E8 hull, shown side by side with the v266
count, at one and at two loops.  Correction proposal; it changes no ledger row.

Question.  v266 adds the E8-allowed (15,1,1) with Delta b4 = 4/3 (B_45 = (-9,-3,-7/3)
against minimal (-31/3,-3,-7/3)).  A real scalar in the 15 of SU(4) gives
Delta b4 = T(15)/6 = 2/3, a complex one 4/3.  The 248 contains exactly one (45,1)
(v247), and the 45 is real.  Is v266's survival verdict ("+(15,1,1) survives
Super-K, Hyper-K decides") an artefact of counting two real adjoints?

Method.  Group theory from scratch (Fractions): the general one- and two-loop gauge
coefficients (Machacek-Vaughn; Weyl fermions, complex scalars, real scalars at
half weight), calibrated on the SM b_i and b_ij.  Two-step chain
M_Z ->(SM)-> M_PS ->(PS)-> M_GUT with v266's tree matching a4=a3, a2L=a2,
a2R^-1=(5/3)a1^-1-(2/3)a3^-1 and v266's inputs alpha_i^-1(M_Z)=(59.01,29.59,8.47).
Two loops: SM gauge two-loop + top Yukawa (leading order), PS gauge two-loop; the
PS top Yukawa enters all three PS factors with the same coefficient and drops out
of the unification conditions (checked).  tau_p(p->e+pi0) with v266's formula
1e36 (M_GUT/1e16)^4 (alpha_U^-1/40)^2 yr and v266's x3 hadronic band, against
Super-K 2.4e34 yr and Hyper-K ~1.4e35 yr.

  [E] 1. CALIBRATION.  The general formula reproduces SM b_i=(41/10,-19/6,-7) and
        b_ij exactly, and v266's PS one-loop betas (-31/3,-3,-7/3) and (-9,-3,-7/3).
  [E] 2. COUNTING ERROR CONFIRMED.  Delta b4(real 15) = 2/3, Delta b4(complex 15) =
        4/3 = v266.  The PS content of the 248 holds exactly one (15,1,1), inside the
        real (45,1); v266's own check 5 says "only ONE 45".  v266 counts two.
  [E] 3. ONE LOOP SIDE BY SIDE.  v266 reproduced (minimal M_GUT=2.43e15,
        v266-(15,1,1) M_GUT=5.89e15); with the real 45: M_PS=4.11e13,
        M_GUT=3.61e15, alpha_U^-1=45.25, tau_p=2.17e34 yr.
  [X] 4. ONE LOOP, REAL 45: v266 CHECK 3 FAILS.  tau_p=2.2e34 yr < Super-K
        centrally; the x3 band [7.2e33, 6.5e34] straddles the bound -- not "proton
        safe now" by v266's own rule (tau/3 > Super-K), only band-marginal.
        Control: the same code with the complex count passes v266's rule
        (1.55e35 yr), so the verdict flips on the counting alone.
  [E] 5. TWO-LOOP CALIBRATION.  The 2L-SM + 1L-PS hybrid reproduces
        experiments/gauge-unification/results/pati_salam.json (M_GUT 1.29e15,
        2.80e15); RK4 and scipy DOP853 agree.
  [X] 6. TWO LOOPS KILL BOTH v266 CONTENTS.  Full two loops: minimal tau_p=2.5e32,
        v266-(15,1,1) 7.2e33, real 45 1.1e33 yr -- all below Super-K even at the
        top of the x3 band.  The v266 survivor fails at two loops independently of
        the counting; the repository's own hybrid already gave 7.7e33.
  [N] 7. E8 CONTENT SCAN (13440 contents, two loops).  With kappa=M_PS/M_s in
        [1,1.7] (PS.KAPPA.01): 11 contents Super-K-safe centrally, 101 in the x3
        band, tau_p = 8.7e33 ... 3.5e34 yr at kappa 1.00-1.43; the band top
        (1.06e35) lies under the Hyper-K reach 1.4e35 (4 band tops exceed 1e35).
        Every survivor needs the real (15,1,1) AND >= 4 of the six real (6,1,1)
        colour-triplet pairs light at M_PS; without the (15,1,1) the kappa-band
        maximum is 3.9e33 yr.  Minimal + one real (15,1,1) is not among them.
        (At one loop the same band holds 1331 / 2902 survivors: loop order decides.)
  [O] 8. RESIDUALS.  Tree matching (no MS-bar threshold constants), step thresholds,
        split multiplets in the scan (survival hypothesis), the x3 hadronic band,
        the gauging fork (QFT4D.FORK.01).

Verdict.  The counting error is real and decisive: v266's "+(15,1,1) survives
Super-K" holds only for two real adjoints, which the E8 hull does not supply, and
only at one loop.  With the single real 45 the content is band-marginal at one
loop and excluded at two loops.  The PS branch survives, if at all, only with other
E8 contents at tau_p ~ 1e34 yr, all within Hyper-K reach.

Python-only (fractions + numpy, scipy for the integrator cross-check).
"""
import itertools
import os
from fractions import Fraction as Fr

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve

from tfpt_constants import check, summary, reset, c3, Mbar

HERE = os.path.dirname(os.path.abspath(__file__))
PI = np.pi
M_Z = 91.1876
M_S = float(c3) ** 3.5 * float(Mbar)                    # scalaron scale ~3.06e13 GeV
AINV = np.array([59.01, 29.59, 8.47])                   # v266 inputs, GUT-normalised
DA3_ERR = 0.065                                          # alpha_s(M_Z)=0.1180(9) in alpha_3^-1
M_T, YT_MT = 172.69, 0.9369                              # top pole mass, y_t(M_t) (MS-bar NNLO)
TAU_SK = 2.4e34                                          # Super-K p->e+pi0 (90% CL), v266
TAU_HK = 1.4e35                                          # Hyper-K reach p->e+pi0, v266
HAD_BAND = 3.0                                           # v266's O(3) hadronic band
KAPPA_BAND = (1.0, 1.7)                                  # PS.KAPPA.01 / v266 check 1

# --------------------------------------------------------------------------- group theory
# irrep data (dim, Dynkin index T, quadratic Casimir C2)
SU4 = {"1": (1, Fr(0), Fr(0)), "4": (4, Fr(1, 2), Fr(15, 8)), "6": (6, Fr(1), Fr(5, 2)),
       "15": (15, Fr(4), Fr(4))}
SU3 = {"1": (1, Fr(0), Fr(0)), "3": (3, Fr(1, 2), Fr(4, 3))}
SU2 = {"1": (1, Fr(0), Fr(0)), "2": (2, Fr(1, 2), Fr(3, 4)), "3": (3, Fr(2), Fr(2))}


def u1(Y):
    """GUT-normalised hypercharge as a 'factor': dim 1, T = C2 = (3/5) Y^2."""
    q = Fr(3, 5) * Fr(Y) ** 2
    return (1, q, q)


def betas(c2_gauge, fields):
    """One- and two-loop gauge coefficients for a product group.

    fields: list of (kind, n, [(dim, T, C2) per factor]) with kind in
    'F' (Weyl fermion), 'C' (complex scalar), 'R' (real scalar, half weight).
    b_i   = -11/3 C2(G_i) + 2/3 sum_F S_i + 1/3 sum_C S_i + 1/6 sum_R S_i
    b_ij  = -34/3 C2(G_i)^2 d_ij + sum_F S_i (2 C2_j + 10/3 C2(G_i) d_ij)
            + sum_C S_i (4 C2_j + 2/3 C2(G_i) d_ij)      (real: half of complex)
    with S_i = T_i * prod_{k != i} dim_k.
    """
    m = len(c2_gauge)
    b = [Fr(-11, 3) * c2_gauge[i] for i in range(m)]
    B = [[(Fr(-34, 3) * c2_gauge[i] ** 2 if i == j else Fr(0)) for j in range(m)] for i in range(m)]
    for kind, n, reps in fields:
        w1 = {"F": Fr(2, 3), "C": Fr(1, 3), "R": Fr(1, 6)}[kind]
        eta = {"F": None, "C": Fr(1), "R": Fr(1, 2)}[kind]
        for i in range(m):
            S = reps[i][1]
            for k in range(m):
                if k != i:
                    S *= reps[k][0]
            b[i] += n * w1 * S
            for j in range(m):
                if kind == "F":
                    B[i][j] += n * S * (2 * reps[j][2] + (Fr(10, 3) * c2_gauge[i] if i == j else 0))
                else:
                    B[i][j] += n * eta * S * (4 * reps[j][2] + (Fr(2, 3) * c2_gauge[i] if i == j else 0))
    return b, B


# SM: factors (U(1)_Y, SU(2)_L, SU(3)_c)
SM_GAUGE = (Fr(0), Fr(2), Fr(3))
SM_GEN = [("F", 3, [u1(Fr(1, 6)), SU2["2"], SU3["3"]]),     # Q
          ("F", 3, [u1(Fr(-2, 3)), SU2["1"], SU3["3"]]),    # u^c
          ("F", 3, [u1(Fr(1, 3)), SU2["1"], SU3["3"]]),     # d^c
          ("F", 3, [u1(Fr(-1, 2)), SU2["2"], SU3["1"]]),    # L
          ("F", 3, [u1(Fr(1)), SU2["1"], SU3["1"]])]        # e^c
SM_HIGGS = [("C", 1, [u1(Fr(1, 2)), SU2["2"], SU3["1"]])]
CT_SM = np.array([17 / 10, 3 / 2, 2.0])                    # top-Yukawa 2L gauge coefficients

# PS: factors (SU(4)_c, SU(2)_L, SU(2)_R)
PS_GAUGE = (Fr(4), Fr(2), Fr(2))
PS_GEN = [("F", 3, [SU4["4"], SU2["2"], SU2["1"]]),           # (4,2,1)
          ("F", 3, [SU4["4"], SU2["1"], SU2["2"]])]           # (4bar,1,2)
# E8-hull scalar menu (v247): six real 10 = (6,1,1)+(1,2,2), four complex 16 =
# (4,2,1)+(4bar,1,2), one real 45 = (15,1,1)+(1,3,1)+(1,1,3)+(6,2,2)
MENU = {  # name: (kind, reps, max from one 248, min)
    "bd_r":   ("R", [SU4["1"], SU2["2"], SU2["2"]], 6, 1),
    "t6_r":   ("R", [SU4["6"], SU2["1"], SU2["1"]], 6, 0),
    "s412_c": ("C", [SU4["4"], SU2["1"], SU2["2"]], 4, 1),
    "s421_c": ("C", [SU4["4"], SU2["2"], SU2["1"]], 4, 0),
    "a15_r":  ("R", [SU4["15"], SU2["1"], SU2["1"]], 1, 0),
    "a31_r":  ("R", [SU4["1"], SU2["3"], SU2["1"]], 1, 0),
    "a13_r":  ("R", [SU4["1"], SU2["1"], SU2["3"]], 1, 0),
    "a622_r": ("R", [SU4["6"], SU2["2"], SU2["2"]], 1, 0),
}
NAMES = list(MENU)
# v266 minimal = complex bidoublet (two real) + complex (4bar,1,2)
MINIMAL = {"bd_r": 2, "s412_c": 1}
CONTENTS = {
    "minimal (v266)": MINIMAL,
    "+(15,1,1) complex (v266 = 2 real)": {**MINIMAL, "a15_r": 2},
    "+(15,1,1) real (E8: one 45)": {**MINIMAL, "a15_r": 1},
}


def ps_betas(content):
    fields = list(PS_GEN)
    for name, n in content.items():
        if n:
            kind, reps, _, _ = MENU[name]
            fields.append((kind, n, reps))
    return betas(PS_GAUGE, fields)


def to_np(b, B):
    return np.array([float(x) for x in b]), np.array([[float(x) for x in r] for r in B])


B1_SM, B2_SM = to_np(*betas(SM_GAUGE, SM_GEN + SM_HIGGS))


# --------------------------------------------------------------------------- 1 loop
def jac1(bp, b_sm=B1_SM):
    """One-loop Jacobian of (a4-a2L, a4-a2R) at M_GUT w.r.t. (L1, L2); bp shape (n,3)."""
    s3, s2, s2r = b_sm[2], b_sm[1], (5 / 3) * b_sm[0] - (2 / 3) * b_sm[2]
    J = np.empty((bp.shape[0], 2, 2))
    J[:, 0, 0] = -(s3 - s2) / (2 * PI)
    J[:, 1, 0] = -(s3 - s2r) / (2 * PI)
    J[:, 0, 1] = -(bp[:, 0] - bp[:, 1]) / (2 * PI)
    J[:, 1, 1] = -(bp[:, 0] - bp[:, 2]) / (2 * PI)
    return J


def solve_1loop(b_ps, ainv=AINV, b_sm=B1_SM):
    """Exact one-loop two-step solve; unknowns L1=ln(M_PS/M_Z), L2=ln(M_GUT/M_PS).
    alpha_i^-1(mu) = alpha_i^-1(M_Z) - b_i/(2 pi) ln(mu/M_Z) on each leg; at M_GUT
    a4 = a2L and a4 = a2R."""
    bp = np.atleast_2d(np.asarray(b_ps, float))
    a1, a2, a3 = ainv
    r0 = np.array([a3 - a2, a3 - ((5 / 3) * a1 - (2 / 3) * a3)])
    x = np.linalg.solve(jac1(bp, b_sm), np.broadcast_to(-r0, (bp.shape[0], 2))[..., None])[..., 0]
    L1, L2 = x[:, 0], x[:, 1]
    aU = a3 - b_sm[2] / (2 * PI) * L1 - bp[:, 0] / (2 * PI) * L2
    out = (M_Z * np.exp(L1), M_Z * np.exp(L1 + L2), aU, L1, L2)
    if np.ndim(b_ps) == 1:
        return tuple(float(v[0]) for v in out)
    return out


def tau_p(M_GUT, aU_inv):
    return 1e36 * (M_GUT / 1e16) ** 4 * (aU_inv / 40) ** 2


# --------------------------------------------------------------------------- 2 loops
def sm_rhs(y, B1, B2, yukawa):
    a = y[..., :3]
    al = 1.0 / a
    d = -B1 / (2 * PI) - al @ B2.T / (8 * PI ** 2)
    if not yukawa:
        return np.concatenate([d, np.zeros(y.shape[:-1] + (1,))], -1)
    yt = y[..., 3:4]
    g2 = 4 * PI * al
    d = d + CT_SM * yt ** 2 / (32 * PI ** 3)
    dyt = yt / (16 * PI ** 2) * (4.5 * yt ** 2 - 0.85 * g2[..., 0:1] - 2.25 * g2[..., 1:2] - 8 * g2[..., 2:3])
    return np.concatenate([d, dyt], -1)


def rk4(f, y, h, n):
    for _ in range(n):
        k1 = f(y)
        k2 = f(y + 0.5 * h * k1)
        k3 = f(y + 0.5 * h * k2)
        k4 = f(y + h * k3)
        y = y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    return y


class SMGrid:
    """Two-loop SM trajectory on a fine ln(mu) grid, linear interpolation."""

    def __init__(self, ainv=AINV, B2=B2_SM, yukawa=True, t_max=np.log(1e17 / M_Z), h=2e-3):
        self.yukawa = yukawa
        f = lambda y: sm_rhs(y, B1_SM, B2, yukawa)
        yt0 = YT_MT
        if yukawa:                                   # y_t(M_Z) from y_t(M_t) (fixed point)
            for _ in range(4):
                tt = np.log(M_T / M_Z)
                nn = 200
                yend = rk4(f, np.r_[ainv, yt0], tt / nn, nn)
                yt0 *= YT_MT / yend[3]
        self.yt_MZ = yt0
        n = int(np.ceil(t_max / h))
        self.t = np.linspace(0, t_max, n + 1)
        hh = self.t[1] - self.t[0]
        ys = np.empty((n + 1, 4))
        y = np.r_[ainv, yt0]
        ys[0] = y
        for k in range(n):
            y = rk4(f, y, hh, 1)
            ys[k + 1] = y
        self.y = ys

    def at(self, L1):
        return np.stack([np.interp(L1, self.t, self.y[:, i]) for i in range(4)], -1)


def ps_rhs(a, B1, B2):
    return -B1 / (2 * PI) - np.einsum("nij,nj->ni", B2, 1.0 / a) / (8 * PI ** 2)


def solve_2loop(B1p, B2p, grid, nsteps=240, iters=40):
    """Chord-Newton on (L1, L2) with the one-loop Jacobian; vectorised over contents."""
    B1p, B2p = np.atleast_2d(B1p), np.atleast_3d(B2p).reshape(-1, 3, 3)
    n = B1p.shape[0]
    _, _, _, L1, L2 = solve_1loop(B1p)
    Jinv = np.linalg.inv(jac1(B1p))
    for _ in range(iters):
        sm = grid.at(L1)
        a0 = np.stack([sm[:, 2], sm[:, 1], (5 / 3) * sm[:, 0] - (2 / 3) * sm[:, 2]], -1)
        h = (L2 / nsteps)[:, None]
        a = a0.copy()
        for _k in range(nsteps):
            k1 = ps_rhs(a, B1p, B2p)
            k2 = ps_rhs(a + 0.5 * h * k1, B1p, B2p)
            k3 = ps_rhs(a + 0.5 * h * k2, B1p, B2p)
            k4 = ps_rhs(a + h * k3, B1p, B2p)
            a = a + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        r = np.stack([a[:, 0] - a[:, 1], a[:, 0] - a[:, 2]], -1)
        dx = np.einsum("nij,nj->ni", Jinv, r)
        L1, L2 = L1 - dx[:, 0], L2 - dx[:, 1]
    res = np.abs(r).max(-1)
    return M_Z * np.exp(L1), M_Z * np.exp(L1 + L2), a[:, 0], res, grid.at(L1)[:, 3]


def solve_2loop_scipy(B1p, B2p, yt0, ainv=AINV):
    """Independent cross-check: DOP853 on both legs, fsolve on (L1, L2)."""
    f_sm = lambda t, y: sm_rhs(y, B1_SM, B2_SM, True)
    f_ps = lambda t, a: ps_rhs(a[None], B1p[None], B2p[None])[0]

    def top(x):
        L1, L2 = x
        s = solve_ivp(f_sm, (0, L1), np.r_[ainv, yt0], method="DOP853", rtol=1e-12, atol=1e-12).y[:, -1]
        a0 = np.array([s[2], s[1], (5 / 3) * s[0] - (2 / 3) * s[2]])
        return solve_ivp(f_ps, (0, L2), a0, method="DOP853", rtol=1e-12, atol=1e-12).y[:, -1]

    def resid(x):
        a = top(x)
        return [a[0] - a[1], a[0] - a[2]]

    _, _, _, L1, L2 = solve_1loop(B1p)
    x = fsolve(resid, [L1, L2], xtol=1e-13)
    return M_Z * np.exp(x[0]), M_Z * np.exp(x[0] + x[1]), top(x)[0]


def classify(t):
    if t * HAD_BAND < TAU_SK:
        return "EXCLUDED (even x3 band top < Super-K)"
    if t / HAD_BAND > TAU_SK:
        return "Super-K-safe (x3 band bottom > Super-K)"
    return "band-marginal (band straddles Super-K)"


def run():
    reset()
    print("v1036  PS.THRESHOLD.02/PS.PROTON.02 recount: real 45 vs v266's complex (15,1,1), 1 and 2 loops")

    # 1. calibration
    bsm, Bsm = betas(SM_GAUGE, SM_GEN + SM_HIGGS)
    want_b = [Fr(41, 10), Fr(-19, 6), Fr(-7)]
    want_B = [[Fr(199, 50), Fr(27, 10), Fr(44, 5)], [Fr(9, 10), Fr(35, 6), Fr(12)], [Fr(11, 10), Fr(9, 2), Fr(-26)]]
    b_min = ps_betas(CONTENTS["minimal (v266)"])[0]
    b_cx = ps_betas(CONTENTS["+(15,1,1) complex (v266 = 2 real)"])[0]
    b_re = ps_betas(CONTENTS["+(15,1,1) real (E8: one 45)"])[0]
    check("CALIBRATION [E]: general formula gives SM b_i=(41/10,-19/6,-7) and the SM b_ij "
          "[[199/50,27/10,44/5],[9/10,35/6,12],[11/10,9/2,-26]] exactly; PS minimal %s = v266 B_MIN "
          "(-31/3,-3,-7/3), two real (15,1,1) %s = v266 B_45 (-9,-3,-7/3)"
          % (tuple(str(x) for x in b_min), tuple(str(x) for x in b_cx)),
          bsm == want_b and Bsm == want_B
          and b_min == [Fr(-31, 3), Fr(-3), Fr(-7, 3)] and b_cx == [Fr(-9), Fr(-3), Fr(-7, 3)])

    # 2. counting: real vs complex 15; the 248 holds one (15,1,1), inside the real (45,1)
    d_real, d_cplx = b_re[0] - b_min[0], b_cx[0] - b_min[0]
    so10_to_ps = {"1": {(1, 1, 1): 1}, "10": {(6, 1, 1): 1, (1, 2, 2): 1},
                  "16": {(4, 2, 1): 1, (4, 1, 2): 1}, "16b": {(4, 2, 1): 1, (4, 1, 2): 1},
                  "45": {(15, 1, 1): 1, (1, 3, 1): 1, (1, 1, 3): 1, (6, 2, 2): 1}}
    e8_248 = [("45", 1), ("1", 15), ("10", 6), ("16", 4), ("16b", 4)]   # (SO(10) rep, SU(4)_fam dim)
    dim_ps = lambda k: k[0] * k[1] * k[2]
    total = sum(fam * sum(dim_ps(k) * m for k, m in so10_to_ps[r].items()) for r, fam in e8_248)
    n15 = sum(fam * so10_to_ps[r].get((15, 1, 1), 0) for r, fam in e8_248)
    check("COUNTING ERROR CONFIRMED [E]: Delta b4(real 15)=T(15)/6=%s, Delta b4(complex 15)=%s = v266's "
          "B_45-B_MIN; the PS content of the 248 (dim %d) holds exactly %d (15,1,1), inside the real (45,1) "
          "(v247; v266 check 5 itself says 'only ONE 45') -- v266 counts two real adjoints"
          % (d_real, d_cplx, total, n15),
          d_real == Fr(2, 3) and d_cplx == Fr(4, 3) and total == 248 and n15 == 1)

    # 3. one loop side by side
    one = {}
    for name, c in CONTENTS.items():
        b1, _ = to_np(*ps_betas(c))
        MPS, MG, aU, _, _ = solve_1loop(b1)
        one[name] = (float(MPS), float(MG), float(aU), float(tau_p(MG, aU)))
    print("\n  one loop (v266 method):")
    for name, (MPS, MG, aU, t) in one.items():
        print("    %-36s M_PS=%.3e (kappa %.2f)  M_GUT=%.3e  aU^-1=%.2f  tau_p=%.2e yr  [%.1e..%.1e]  %s"
              % (name, MPS, MPS / M_S, MG, aU, t, t / HAD_BAND, t * HAD_BAND, classify(t)))
    m1, c1, r1 = one["minimal (v266)"], one["+(15,1,1) complex (v266 = 2 real)"], one["+(15,1,1) real (E8: one 45)"]
    check("ONE LOOP SIDE BY SIDE [E]: v266 reproduced (minimal M_PS=%.2e, M_GUT=%.2e; complex 15 "
          "M_PS=%.2e, M_GUT=%.2e); real 45: M_PS=%.2e (kappa %.2f), M_GUT=%.2e, alpha_U^-1=%.2f, tau_p=%.2e yr"
          % (m1[0], m1[1], c1[0], c1[1], r1[0], r1[0] / M_S, r1[1], r1[2], r1[3]),
          abs(m1[1] / 2.4e15 - 1) < 0.03 and abs(c1[1] / 5.9e15 - 1) < 0.03
          and KAPPA_BAND[0] <= r1[0] / M_S <= KAPPA_BAND[1] and r1[0] < r1[1])

    # 4. one loop, real 45: v266's survival rule fails
    check("ONE LOOP REAL 45 -- v266 CHECK 3 FAILS [X]: tau_p=%.2e yr < Super-K %.1e centrally; band "
          "[%.1e, %.1e] straddles the bound, so v266's rule 'tau/3 > Super-K' is violated (v266 had "
          "%.2e with the complex count); band-marginal, not proton-safe"
          % (r1[3], TAU_SK, r1[3] / HAD_BAND, r1[3] * HAD_BAND, c1[3]),
          r1[3] < TAU_SK and r1[3] / HAD_BAND < TAU_SK < r1[3] * HAD_BAND and c1[3] / HAD_BAND > TAU_SK)
    # alpha_s sensitivity
    b1r = to_np(*ps_betas(CONTENTS["+(15,1,1) real (E8: one 45)"]))[0]
    t_lo = [tau_p(*solve_1loop(b1r, AINV + np.array([0, 0, s * DA3_ERR]))[1:3]) for s in (-1, 1)]
    print("    alpha_s error (+-%.3f in alpha_3^-1), real 45, 1L: tau_p = %.2e ... %.2e yr"
          % (DA3_ERR, min(t_lo), max(t_lo)))

    # 5. two-loop calibration: hybrid (2L SM gauge-only, 1L PS) with the experiment's inputs
    sw2, aem, als = 0.23122, 127.951, 0.1179
    ainv_exp = np.array([(3 / 5) * (1 - sw2) * aem, sw2 * aem, 1 / als])
    g_exp = SMGrid(ainv_exp, yukawa=False)
    zero = np.zeros((1, 3, 3))
    hyb = {}
    for name in ("minimal (v266)", "+(15,1,1) complex (v266 = 2 real)", "+(15,1,1) real (E8: one 45)"):
        b1, _ = to_np(*ps_betas(CONTENTS[name]))
        MPS, MG, aU, _, _ = solve_2loop(b1[None], zero, g_exp)
        hyb[name] = (MPS[0], MG[0], aU[0], tau_p(MG[0], aU[0]))
    exp_json = {"minimal (v266)": (3.7410e13, 1.2906e15), "+(15,1,1) complex (v266 = 2 real)": (3.6060e13, 2.7989e15)}
    dev = max(max(abs(hyb[k][0] / v[0] - 1), abs(hyb[k][1] / v[1] - 1)) for k, v in exp_json.items())
    # integrator cross-check (full 2L, real 45)
    grid = SMGrid()
    b1r, B2r = to_np(*ps_betas(CONTENTS["+(15,1,1) real (E8: one 45)"]))
    MPS_rk, MG_rk, aU_rk, _, _ = solve_2loop(b1r[None], B2r[None], grid)
    MPS_sp, MG_sp, aU_sp = solve_2loop_scipy(b1r, B2r, grid.yt_MZ)
    dint = max(abs(MPS_rk[0] / MPS_sp - 1), abs(MG_rk[0] / MG_sp - 1))
    print("\n  2L-SM + 1L-PS hybrid (experiment inputs): " + "; ".join(
        "%s M_PS=%.3e M_GUT=%.3e tau_p=%.2e" % (k.split(" (")[0], *hyb[k][:2], hyb[k][3]) for k in hyb))
    check("TWO-LOOP CALIBRATION [E]: hybrid reproduces experiments/gauge-unification/results/pati_salam.json "
          "(M_PS, M_GUT of minimal and complex-15 to %.2f%%); RK4 chord-Newton vs scipy DOP853+fsolve "
          "(full 2L, real 45) agree to %.1e" % (100 * dev, dint),
          dev < 0.01 and dint < 1e-6)

    # 6. full two loops: both v266 contents excluded
    two = {}
    for name, c in CONTENTS.items():
        b1, B2 = to_np(*ps_betas(c))
        MPS, MG, aU, res, yt = solve_2loop(b1[None], B2[None], grid)
        two[name] = (MPS[0], MG[0], aU[0], tau_p(MG[0], aU[0]), res[0], yt[0])
    print("\n  full two loops (SM 2L + y_t, PS 2L, tree matching):")
    for name, (MPS, MG, aU, t, res, _) in two.items():
        print("    %-36s M_PS=%.3e (kappa %.2f)  M_GUT=%.3e  aU^-1=%.2f  tau_p=%.2e yr  [%.1e..%.1e]  %s  (res %.0e)"
              % (name, MPS, MPS / M_S, MG, aU, t, t / HAD_BAND, t * HAD_BAND, classify(t), res))
    t2 = {k: v[3] for k, v in two.items()}
    check("TWO LOOPS KILL BOTH v266 CONTENTS [X]: minimal tau_p=%.2e, complex 15 %.2e, real 45 %.2e yr -- "
          "all below Super-K even at the top of the x3 band (the repository's own hybrid gave %.1e for "
          "the complex count)" % (t2["minimal (v266)"], t2["+(15,1,1) complex (v266 = 2 real)"],
                                  t2["+(15,1,1) real (E8: one 45)"], hyb["+(15,1,1) complex (v266 = 2 real)"][3]),
          all(t * HAD_BAND < TAU_SK for t in t2.values()) and max(v[4] for v in two.values()) < 1e-9)
    # alpha_s sensitivity at two loops
    t_a = []
    for s in (-1, 1):
        g_s = SMGrid(AINV + np.array([0, 0, s * DA3_ERR]))
        _, MG, aU, _, _ = solve_2loop(b1r[None], B2r[None], g_s)
        t_a.append(tau_p(MG[0], aU[0]))
    print("    alpha_s error, real 45, 2L: tau_p = %.2e ... %.2e yr (x3 band top %.1e < Super-K)"
          % (min(t_a), max(t_a), max(t_a) * HAD_BAND))

    # 6b. PS top Yukawa: universal coefficient (bidoublet Yukawa (4,2,1)(4bar,1,2)(1,2,2))
    pairs = 4 * 2 * 2                                    # SU(4) x SU(2)_L x SU(2)_R index pairs
    c_ps = [pairs * (2 * SU4["4"][2]) / 15, pairs * SU2["2"][2] / 3, pairs * SU2["2"][2] / 3]
    c_sm = [6 * (u1(Fr(1, 6))[2] + u1(Fr(-2, 3))[2]), 6 * SU2["2"][2] / 3, 6 * 2 * SU3["3"][2] / 8]
    rr = two["+(15,1,1) real (E8: one 45)"]
    dAU = float(c_ps[0]) * rr[5] ** 2 / (32 * PI ** 3) * np.log(rr[1] / rr[0])
    check("PS TOP YUKAWA DROPS OUT [E]: the same pair counting gives the SM coefficients %s = (17/10,3/2,2) "
          "and for the PS bidoublet Yukawa %s, equal in all three PS factors -> no shift of M_PS, M_GUT; "
          "alpha_U^-1 moves by <= %.4f" % (tuple(str(x) for x in c_sm), tuple(str(x) for x in c_ps), dAU),
          c_sm == [Fr(17, 10), Fr(3, 2), Fr(2)] and c_ps[0] == c_ps[1] == c_ps[2] and dAU < 0.01)

    # 7. E8 content scan at two loops
    rng = [range(MENU[k][3], MENU[k][2] + 1) for k in NAMES]
    cnt = np.array(list(itertools.product(*rng)), float)
    unit = {k: to_np(*betas(PS_GAUGE, [(MENU[k][0], 1, MENU[k][1])])) for k in NAMES}
    for k in NAMES:                                       # remove the gauge part counted per unit
        unit[k][1][np.diag_indices(3)] -= np.array([-34 / 3 * float(c) ** 2 for c in PS_GAUGE])
        unit[k][0][:] -= np.array([-11 / 3 * float(c) for c in PS_GAUGE])
    base1, base2 = to_np(*betas(PS_GAUGE, PS_GEN))
    B1s = base1 + sum(cnt[:, [i]] * unit[k][0] for i, k in enumerate(NAMES))
    B2s = base2 + sum(cnt[:, i, None, None] * unit[k][1] for i, k in enumerate(NAMES))
    target = {**{k: 0 for k in NAMES}, **CONTENTS["+(15,1,1) real (E8: one 45)"]}
    is_re45 = np.array([dict(zip(NAMES, r)) == target for r in cnt])
    j_re = int(np.where(is_re45)[0][0])
    col = {k: cnt[:, i] for i, k in enumerate(NAMES)}
    MPSs, MGs, aUs, ress, _ = solve_2loop(B1s, B2s, grid)
    tau = tau_p(MGs, aUs)
    kap = MPSs / M_S
    ok = (ress < 1e-8) & (MGs > MPSs)
    safe_c, safe_b = ok & (tau > TAU_SK), ok & (tau * HAD_BAND > TAU_SK)
    inb = ok & (kap >= KAPPA_BAND[0]) & (kap <= KAPPA_BAND[1])
    sel_c, sel_b = inb & safe_c, inb & safe_b
    no15 = inb & (col["a15_r"] == 0)
    MPS1, MG1, aU1, _, _ = solve_1loop(B1s)            # one-loop scan for contrast
    tau1, kap1 = tau_p(MG1, aU1), MPS1 / M_S
    inb1 = (kap1 >= KAPPA_BAND[0]) & (kap1 <= KAPPA_BAND[1])
    print("\n  E8 content scan, two loops: %d contents, %d converged; Super-K-safe central %d, x3 band %d"
          % (len(cnt), ok.sum(), safe_c.sum(), safe_b.sum()))
    print("    kappa in [%.1f, %.1f]: %d contents; Super-K-safe central %d, x3 band %d; band-safe tau_p = %.2e ... %.2e yr, "
          "kappa = %.2f ... %.2f, band top max %.2e (Hyper-K %.1e; band tops above 1e35: %d)"
          % (*KAPPA_BAND, inb.sum(), sel_c.sum(), sel_b.sum(), tau[sel_b].min(), tau[sel_b].max(),
             kap[sel_b].min(), kap[sel_b].max(), HAD_BAND * tau[sel_b].max(), TAU_HK,
             (sel_b & (HAD_BAND * tau > 1e35)).sum()))
    print("    every band-safe kappa-band content has a15_r=1: %s; t6_r >= %d; without the (15,1,1) the kappa-band "
          "maximum is tau_p=%.2e yr (band top %.2e)"
          % (bool((col["a15_r"][sel_b] == 1).all()), int(col["t6_r"][sel_b].min()), tau[no15].max(),
             HAD_BAND * tau[no15].max()))
    print("    one-loop scan for contrast: kappa-band Super-K-safe central %d, x3 band %d (loop order decides)"
          % ((inb1 & (tau1 > TAU_SK)).sum(), (inb1 & (HAD_BAND * tau1 > TAU_SK)).sum()))
    print("    largest M_GUT overall: %.2e GeV (kappa %.2f, tau_p %.1e yr); scan point = real-45 content: M_GUT=%.3e (direct %.3e)"
          % (MGs[ok].max(), kap[ok][np.argmax(MGs[ok])], tau[ok][np.argmax(MGs[ok])], MGs[j_re], rr[1]))
    for i in np.where(sel_c)[0][np.argsort(-tau[sel_c])][:4]:
        c = {k: int(v) for k, v in zip(NAMES, cnt[i]) if v}
        print("      central survivor: %s  kappa=%.2f  M_GUT=%.2e  tau_p=%.2e yr" % (c, kap[i], MGs[i], tau[i]))
    check("E8 CONTENT SCAN [N]: two loops, %d contents (all converged: %s); kappa in [1,1.7]: %d Super-K-safe "
          "centrally, %d in the x3 band, tau_p = %.1e ... %.1e yr at kappa <= %.2f; all of them need the real "
          "(15,1,1) AND >= %d real (6,1,1) at M_PS (without the 15: max %.1e yr); minimal + real 15 is not among "
          "them; the band top %.2e stays under the Hyper-K reach %.1e"
          % (len(cnt), bool(ok.all()), sel_c.sum(), sel_b.sum(), tau[sel_b].min(), tau[sel_b].max(),
             kap[sel_b].max(), int(col["t6_r"][sel_b].min()), tau[no15].max(), HAD_BAND * tau[sel_b].max(), TAU_HK),
          ok.all() and sel_b.sum() > 0 and HAD_BAND * tau[sel_b].max() < TAU_HK
          and (col["a15_r"][sel_b] == 1).all() and HAD_BAND * tau[no15].max() < TAU_SK
          and not sel_b[j_re] and abs(MGs[j_re] / rr[1] - 1) < 1e-6)

    # 8. residuals
    check("RESIDUALS [O]: tree matching at M_PS and M_GUT (no MS-bar threshold constants), step "
          "thresholds, split multiplets in the scan (survival hypothesis), the x3 hadronic band, and "
          "the gauging fork (QFT4D.FORK.01) stay open; none of them is a fit knob for the v266 verdict", True)

    return summary("v1036 PS recount: v266 counts (15,1,1) complex; real 45 -> band-marginal (1L), excluded (2L)")


if __name__ == "__main__":
    raise SystemExit(1 if run() else 0)
