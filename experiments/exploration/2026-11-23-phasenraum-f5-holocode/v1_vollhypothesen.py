"""V: Prüfstand für vollständige Hypothesen.

Jede Hypothese muss die ganze Lösung abdecken. Jede wird gegen denselben Katalog K1..K14 geprüft.
Ein Fehlschlag beendet die Prüfung nicht; am Ende steht eine Gesamttabelle.

Hypothesen
  H1  Zeitneutrale Fermionnaht: flache orientierte 2D-Naht mit einem Diracfeld, keine Zeitrichtung
      ausgezeichnet, Kragen der Breite l.
  H2  Turm der ganzzahligen Divisionsalgebren: Naht = Gaußsche Gerade Z[u] in den Hurwitz-Quaternionen,
      euklidische Raumzeit = Hurwitz-Gitter (größte assoziative Divisionsalgebra), Familien = die drei
      Gaußschen Geraden, inneres E8 = Hurwitz + Hurwitz, verklebt entlang (1+i).
  H3  Spektrales Produkt M4 x F: hyperkubischer 4-Torus, Dirac D = D_geom (x) 1 + gamma (x) D_F.
  H4  E8-GSO-Naht: 16 Majorana-Weyl-Fermionen auf der zeitneutralen Naht, Summe über Spinstrukturen.
  H5  Turm + GSO: H2, mit den Nahtfermionen aus H4; das innere E8 ist das Gitter der GSO-Naht.

Katalog
  K1  vier Marken mit zwei Blättern          K8  P2 (5 Plätze, 3+2) über P1 (P2.PARTITION.01)
  K2  quadratische Lage, tau = i               K9  drei Familien als gleiche Kopien
  K3  mu4 permutiert die Marken zyklisch       K10 Familiencharaktere {-1, i, -i}
  K4  Spinstruktur (A,A) für den Zustand       K11 E8 (240 Wurzeln) aus dem Prinzip
  K5  c3 = 1/(8 pi) ohne Wahl                  K12 3+1 aus dem Prinzip
  K6  Zustand und Zeit, Breite festgelegt      K13 Gravitation (Einstein mit 1/c3, TFPT bedingt)
  K7  flacher Sektor für das E8 Wort           K14 keine freien Parameter
  K15 Eichdynamik (Yang-Mills) in 4D           K16 Kontinuumslimes
  K17 Rahmungsanomalie der Naht behandelt      K18 Neutrinomischung gegen Daten
  K19 Existenz der Naht (Ursprung)
"""
import itertools
import json
from fractions import Fraction as Fr
from pathlib import Path

import numpy as np
from mpmath import mp, jtheta, exp, pi, qp, mpc

mp.dps = 25
facts = {}

# ---------------------------------------------------------------- 2D seam facts
# F1 time-neutral selection (orientation preserving exchange of the two cycles)
sol = []
for x in np.linspace(-0.5, 0.5, 101):
    for y in np.linspace(0.8, 2.0, 121):
        t = complex(x, y)
        if abs(t) < 1 - 1e-12:
            continue
        u = -1 / t
        if abs(t.imag - u.imag) < 1e-9 and abs(t.real - u.real) < 1e-9:
            sol.append((round(x, 6), round(y, 6)))
facts["F1_time_neutral_tau"] = sorted(set(sol))

# F2 centres c for which the quarter turn z -> i(z-c)+c permutes the four 2-torsion points cyclically
two_torsion = [complex(a / 2, b / 2) for a in (0, 1) for b in (0, 1)]


def red(z):
    return complex(round((z.real % 1) * 8) / 8 % 1, round((z.imag % 1) * 8) / 8 % 1)


cyc_centres, fix_centres = [], []
for a in range(4):
    for b in range(4):
        c = complex(a / 4, b / 4)
        perm = [two_torsion.index(red(1j * (p - c) + c)) if red(1j * (p - c) + c) in two_torsion else None
                for p in two_torsion]
        if None in perm:
            continue
        # cycle type
        seen, cycles = set(), []
        for s in range(4):
            if s not in seen:
                k, n = s, 0
                while k not in seen:
                    seen.add(k); k = perm[k]; n += 1
                cycles.append(n)
        (cyc_centres if sorted(cycles) == [4] else fix_centres).append((c, sorted(cycles)))
facts["F2_cyclic_centres"] = [str(c) for c, _ in cyc_centres]
facts["F2_other_centres"] = [(str(c), cy) for c, cy in fix_centres]

# F3 c3 per time-neutral orbit
facts["F3_c3"] = 2 / (4 * np.pi * 4)

# F4 collar width from time neutrality on thermal circle x seam
M = 30


def lev(lx, ly, n=40):
    a = np.arange(-M, M) + 0.5
    kx, ky = np.meshgrid(2 * np.pi * a / lx, 2 * np.pi * a / ly)
    return np.sort(np.sqrt(kx ** 2 + ky ** 2).ravel())[:n]


widths = [w for w in np.linspace(0.2, 1.2, 101) if np.max(np.abs(2 * w * lev(1, 1) - 1 * lev(2 * w, 1))) < 1e-9]
facts["F4_width"] = [round(w, 6) for w in widths]

# F5 zero modes of the seam Dirac per spin structure (square torus)
facts["F5_zero_modes_per_dirac"] = {"AA": 0, "AP": 0, "PA": 0, "PP": 2}

# F6 GSO E8 on the seam: 16 Majorana-Weyl
nmw = 16
facts["F6_E8_from_GSO"] = {"so16_bilinears": nmw * (nmw - 1) // 2, "ramond_chiral_ground_states": 2 ** (nmw // 2) // 2,
                           "total": nmw * (nmw - 1) // 2 + 2 ** (nmw // 2) // 2, "complex_zero_modes_PP": nmw // 2}


def ZE8(tau):
    nome = exp(1j * pi * tau)
    eta = exp(1j * pi * tau / 12) * qp(exp(2j * pi * tau))
    th = [jtheta(k, 0, nome) for k in (2, 3, 4)]
    return (th[0] ** 8 + th[1] ** 8 + th[2] ** 8) / (2 * eta ** 8)


z_i, z_s = ZE8(mpc(0, 1)), ZE8(-1 / mpc(0, 1))
t0 = mpc(0.13, 1.17)
facts["F6_ZE8_S_invariant"] = bool(abs(ZE8(t0) - ZE8(-1 / t0)) < 1e-15)
facts["F6_ZE8_at_i"] = float(z_i.real)

# ---------------------------------------------------------------- quaternion facts


def qmul(p, q):
    a1, b1, c1, d1 = p; a2, b2, c2, d2 = q
    return (a1 * a2 - b1 * b2 - c1 * c2 - d1 * d2, a1 * b2 + b1 * a2 + c1 * d2 - d1 * c2,
            a1 * c2 - b1 * d2 + c1 * a2 + d1 * b2, a1 * d2 + b1 * c2 - c1 * b2 + d1 * a2)


def qconj(q):
    return (q[0], -q[1], -q[2], -q[3])


def qnorm(q):
    return sum(x * x for x in q)


def hurwitz(q):
    return all(x.denominator == 1 for x in q) or all(x.denominator == 2 for x in q)


half = Fr(1, 2)
units = [q for q in itertools.product([Fr(v) for v in (-1, 0, 1)], repeat=4) if qnorm(q) == 1]
units += [q for q in itertools.product((half, -half), repeat=4)]
one, I, J, K = (Fr(1), Fr(0), Fr(0), Fr(0)), (Fr(0), Fr(1), Fr(0), Fr(0)), (Fr(0), Fr(0), Fr(1), Fr(0)), (Fr(0), Fr(0), Fr(0), Fr(1))
axes = [tuple(s * x for x in e) for e in (one, I, J, K) for s in (1, -1)]
facts["F7_units"] = len(units)
facts["F7_axis_orbit_left"] = len({qmul(u, one) for u in units} & set(axes)) == 8


def qorder(q):
    x, n = q, 1
    while x != one:
        x, n = qmul(x, q), n + 1
        if n > 12:
            return None
    return n


orders = [qorder(u) for u in units]
facts["F7_order_census"] = {k: orders.count(k) for k in sorted(set(orders))}
facts["F7_gaussian_lines"] = orders.count(4) // 2
facts["F7_eisenstein_lines"] = (orders.count(3) + orders.count(6)) // 4
sigma = (-half, half, half, half)
conj = lambda u, q: qmul(qmul(u, q), qconj(u))
img = [conj(sigma, x) for x in (I, J, K)]
facts["F7_sigma_cycles_ijk"] = sorted(img) == sorted([I, J, K]) and all(a != b for a, b in zip(img, (I, J, K)))
facts["F7_seam_quarter_on_lines"] = [conj(I, x) for x in (I, J, K)]

# F8 spin structures on R^4 / Hurwitz invariant under the units and invertible
Bm = np.array([[0.5, 0, 0, 0], [0.5, 1, 0, 0], [0.5, 0, 1, 0], [0.5, 0, 0, 1]])   # columns: (1+i+j+k)/2, i, j, k
Bd = np.linalg.inv(Bm).T


def lmat(u):
    return np.array([[float(x) for x in qmul(u, e)] for e in (one, I, J, K)]).T


def rmat(u):
    return np.array([[float(x) for x in qmul(e, u)] for e in (one, I, J, K)]).T


def spin_census(mats):
    good = []
    for e in itertools.product((0, 1), repeat=4):
        s = Bd @ (np.array(e) / 2)
        inv = all(np.allclose(c := np.linalg.solve(Bd, R @ s - s), np.round(c), atol=1e-9) for R in mats)
        if inv and e != (0, 0, 0, 0):
            good.append(e)
    return good


facts["F8_hurwitz_spin_left"] = spin_census([lmat(u) for u in units])
facts["F8_hurwitz_spin_left_right"] = spin_census([lmat(u) for u in units] + [rmat(u) for u in units])
# hypercubic Z^4 under signed permutations of determinant +1
Bd4 = np.eye(4)
hyper = []
for perm in itertools.permutations(range(4)):
    for signs in itertools.product((1, -1), repeat=4):
        R = np.zeros((4, 4))
        for r, c in enumerate(perm):
            R[r, c] = signs[r]
        if round(np.linalg.det(R)) == 1:
            hyper.append(R)
s = np.full(4, 0.5)
facts["F8_hypercubic_AAAA_invariant"] = all(np.allclose(R @ s - s, np.round(R @ s - s)) for R in hyper)

# F8b Lipschitz order Z<1,i,j,k> (units Q8): time neutral and invertible spin structure?
q8 = [tuple(sgn * x for x in e) for e in (one, I, J, K) for sgn in (1, -1)]
L4 = np.eye(4)
s4 = np.full(4, 0.5)
lip_mats = [lmat(u) for u in q8] + [rmat(u) for u in q8]
facts["F8b_lipschitz_units"] = len(q8)
facts["F8b_lipschitz_axis_transitive"] = len({qmul(u, one) for u in q8} & set(axes)) == 8
facts["F8b_lipschitz_AAAA_invariant"] = all(np.allclose(R @ s4 - s4, np.round(R @ s4 - s4)) for R in lip_mats)
Rsig = np.array([[float(x) for x in conj(sigma, e)] for e in (one, I, J, K)]).T
facts["F8b_sigma_preserves_lipschitz_and_AAAA"] = bool(np.allclose(Rsig, np.round(Rsig)) and
                                                     np.allclose(Rsig @ s4 - s4, np.round(Rsig @ s4 - s4)))
facts["F8b_lipschitz_gaussian_lines"] = sum(1 for u in q8 if qorder(u) == 4) // 2
# E8 from Z^8 = Lipschitz^2: D8 plus the half-integer glue with even coordinate sum
d8 = [v for v in itertools.product(range(-1, 2), repeat=8) if sum(x * x for x in v) == 2 and sum(v) % 2 == 0]
glue = [v for v in itertools.product((0.5, -0.5), repeat=8) if round(sum(v)) % 2 == 0]
facts["F9b_E8_from_Z8"] = {"D8_roots": len(d8), "glue": len(glue), "total": len(d8) + len(glue)}

# F9 E8 = {(x, y) in Hurwitz^2 : x = y mod (1+i)}
small = [q for q in itertools.product([Fr(v) for v in range(-2, 3)], repeat=4) if qnorm(q) <= 4]
small += [q for q in itertools.product([Fr(v, 2) for v in (-3, -1, 1, 3)], repeat=4) if qnorm(q) <= 4]
one_plus_i = (Fr(1), Fr(1), Fr(0), Fr(0))
inv_opi = (half, -half, Fr(0), Fr(0))                       # (1+i)^-1 = (1-i)/2


def congruent(x, y):
    d = tuple(a - b for a, b in zip(x, y))
    return hurwitz(qmul(inv_opi, d))


by_norm = {}
for q in small:
    by_norm.setdefault(qnorm(q), []).append(q)
count = {2: 0, 4: 0}
for target in (2, 4):
    for nx in range(0, target + 1):
        ny = target - nx
        for x in by_norm.get(nx, []):
            for y in by_norm.get(ny, []):
                if congruent(x, y):
                    count[target] += 1
facts["F9_E8_counts"] = count
facts["F9_hurwitz_norm_counts"] = {n: len(by_norm.get(n, [])) for n in range(5)}

# F10 associativity: quaternions yes, octonions (Cayley-Dickson) no


def omul(a, b):
    (p, q), (r, s) = a, b
    return (tuple(x - y for x, y in zip(qmul(p, r), qmul(qconj(s), q))),
            tuple(x + y for x, y in zip(qmul(s, p), qmul(q, qconj(r)))))


z4 = (Fr(0),) * 4
oi, oj, ol = (I, z4), (J, z4), (z4, one)
facts["F10_octonion_associator_nonzero"] = omul(omul(oi, oj), ol) != omul(oi, omul(oj, ol))
facts["F10_quaternion_associative"] = qmul(qmul(I, J), K) == qmul(I, qmul(J, K))

# F11 families as self-dual 2-forms: action of the seam quarter turn q -> z8 q z8 and of sigma on Lambda^2_+
z8 = (Fr(1), Fr(1), Fr(0), Fr(0))                       # 1 + i (normalised below)


def rot_from(fun):
    return np.array([fun(np.eye(4)[:, c]) for c in range(4)]).T


def qnum(v):
    return tuple(Fr(0) + 0 for _ in range(4))


def apply_q(a, v, b):
    """v in R^4 as quaternion, return a v b (numeric, a and b given as numpy quaternions)."""
    def mul(p, q):
        a1, b1, c1, d1 = p; a2, b2, c2, d2 = q
        return np.array([a1*a2-b1*b2-c1*c2-d1*d2, a1*b2+b1*a2+c1*d2-d1*c2, a1*c2-b1*d2+c1*a2+d1*b2, a1*d2+b1*c2-c1*b2+d1*a2])
    return mul(mul(a, v), b)


zeta8 = np.array([1, 1, 0, 0]) / np.sqrt(2)
omega = np.array([1, 1, 1, 1]) / 2
omega_bar = np.array([1, -1, -1, -1]) / 2
Rq = rot_from(lambda v: apply_q(zeta8, v, zeta8))          # seam quarter turn: rotation by 90 deg in the (1, i) plane
Rs = rot_from(lambda v: apply_q(omega, v, omega_bar))      # sigma: fixes 1, cycles i, j, k
pairs = list(itertools.combinations(range(4), 2))


def wedge2(R):
    W = np.zeros((6, 6))
    for c, (a, b) in enumerate(pairs):
        for r, (x, y) in enumerate(pairs):
            W[r, c] = R[x, a] * R[y, b] - R[y, a] * R[x, b]
    return W


# Hodge star on Lambda^2 R^4 (basis e01, e02, e03, e12, e13, e23)
star = np.zeros((6, 6))
for c, (a, b) in enumerate(pairs):
    rest = [x for x in range(4) if x not in (a, b)]
    perm = [a, b] + rest
    sgn = round(np.linalg.det(np.eye(4)[perm]))
    star[pairs.index(tuple(rest)), c] = sgn
wp, vp = np.linalg.eigh((star + star.T) / 2)
Vplus = vp[:, wp > 0]
def on_plus(R):
    W = wedge2(R)
    return Vplus.T @ W @ Vplus
evq = np.sort_complex(np.round(np.linalg.eigvals(on_plus(Rq)), 9))
evs = np.round(np.linalg.eigvals(on_plus(Rs)), 9)
facts["F11_quarter_on_selfdual"] = [complex(x) for x in evq]
facts["F11_quarter_on_selfdual_times_q_character"] = sorted([complex(-x) for x in evq], key=lambda c: (c.real, c.imag))
facts["F11_sigma_on_selfdual_order3"] = bool(np.allclose(np.linalg.matrix_power(on_plus(Rs), 3), np.eye(3)) and
                                            not np.allclose(on_plus(Rs), np.eye(3)))
facts["F11_quarter_fixes_normal_plane"] = bool(np.allclose(Rq[:, 2:], np.eye(4)[:, 2:]))


# ---------------------------------------------------------------- hardening of H6
# F12 basis independence of a chiral fermion seam: GSO sum of N Majorana-Weyl fermions must be S-invariant and
#     T-invariant up to a constant phase; only N = 0 mod 16 works, smallest N = 16 -> c = 8, lattice D8+ = E8


def f_gso(N, tau):
    nome = exp(1j * pi * tau)
    eta = exp(1j * pi * tau / 12) * qp(exp(2j * pi * tau))
    th2, th3, th4 = (jtheta(k, 0, nome) for k in (2, 3, 4))
    return (th3 ** (N // 2) + th4 ** (N // 2) + th2 ** (N // 2)) / (2 * eta ** (N // 2))


gso = {}
for N in (8, 16, 24, 32, 48):
    t1, t2 = mpc(0.11, 1.07), mpc(-0.23, 1.31)
    s_ok = all(abs(f_gso(N, t) / f_gso(N, -1 / t) - 1) < 1e-12 for t in (t1, t2))
    r1, r2 = f_gso(N, t1 + 1) / f_gso(N, t1), f_gso(N, t2 + 1) / f_gso(N, t2)
    t_ok = abs(abs(r1) - 1) < 1e-12 and abs(r1 - r2) < 1e-12
    gso[N] = {"S": bool(s_ok), "T_bis_auf_Phase": bool(t_ok), "Phase/2pi": float(mp.arg(r1) / (2 * pi)) if t_ok else None}
facts["F12_gso_modular"] = gso
facts["F12_smallest_N"] = min(N for N, v in gso.items() if v["S"] and v["T_bis_auf_Phase"])
facts["F12_central_charge"] = facts["F12_smallest_N"] / 2

# F13 the family group: <sigma, quarter turn> on Lambda^2_+ and its character against S4 (standard x sign)
gens3 = [on_plus(Rs), on_plus(Rq)]
grp = [np.eye(3)]
frontier = [np.eye(3)]
while frontier:
    nxt = []
    for g in frontier:
        for h in gens3:
            m = np.round(g @ h, 9)
            if not any(np.allclose(m, x) for x in grp):
                grp.append(m); nxt.append(m)
    frontier = nxt
facts["F13_family_group_order"] = len(grp)


def order_of(m):
    x, n = m, 1
    while not np.allclose(x, np.eye(len(m))):
        x, n = x @ m, n + 1
    return n


ours = sorted((order_of(g), round(float(np.trace(g)), 6)) for g in grp)
s4 = []
for perm in itertools.permutations(range(4)):
    P4 = np.zeros((4, 4))
    for a, b in enumerate(perm):
        P4[b, a] = 1
    std = np.trace(P4) - 1
    sign = round(np.linalg.det(P4))
    s4.append((order_of(P4), round(float(std * sign), 6)))
facts["F13_is_S4_standard_times_sign"] = ours == sorted(s4)
# H^1 of the four marks = standard representation: quarter turn as 4-cycle vs as (1,1,2)
four_cycle = np.roll(np.eye(4), 1, axis=0)
fix_swap = np.eye(4)[[0, 1, 3, 2]]
h1 = lambda P4: sorted(np.round(np.linalg.eigvals(P4), 9).tolist(), key=lambda c: (c.real, c.imag))
ev4 = [complex(x) for x in h1(four_cycle)]
ev4.remove(1 + 0j)
evf = [complex(x) for x in h1(fix_swap)]
evf.remove(1 + 0j)
lam_sign = sorted([-x for x in facts["F11_quarter_on_selfdual"]], key=lambda c: (c.real, c.imag))
facts["F13_H1_cyclic"] = sorted(ev4, key=lambda c: (c.real, c.imag))
facts["F13_H1_centre_at_mark"] = sorted(evf, key=lambda c: (c.real, c.imag))
facts["F13_selfdual_times_sign"] = lam_sign
facts["F13_centre_forced"] = (np.allclose(facts["F13_H1_cyclic"], lam_sign) and
                              not np.allclose(facts["F13_H1_centre_at_mark"], lam_sign))
facts["F13_quarter_turn_is_lattice_automorphism"] = bool(np.allclose(Rq, np.round(Rq)) and
                                                       np.allclose(Rq @ s4_ if (s4_ := np.full(4, .5)) is not None else 0,
                                                                   Rq @ s4_) and
                                                       np.allclose((Rq @ s4_ - s4_), np.round(Rq @ s4_ - s4_)))

# F14 chiral middle forms: three of them only in d = 4
from math import comb
facts["F14_chiral_middle_forms"] = {d: comb(d, d // 2) // 2 for d in range(2, 13, 2) if (d // 2) % 2 == 0}

# F15 inputs of TFPT's gravity derivation supplied by H6
facts["F15_gravity_inputs"] = {"c_seam": facts["F12_central_charge"], "c_TFPT_g_car_plus_N_fam": 5 + 3,
                               "mu4": 4, "Z2_times_chi": 2 * 2, "einstein_coeff_2pi_over_eta": 2 * np.pi * 4,
                               "one_over_c3": 1 / facts["F3_c3"]}

# ---------------------------------------------------------------- criteria
P, LES, F, O, IN = "PASS", "LESART", "FAIL", "OFFEN", "EINGABE"
cyc_ok = bool(facts["F2_cyclic_centres"])
seam_common = {
    "K2": (P, "orientierte Zeitneutralität erzwingt tau = i (F1)"),
    "K3": (LES if cyc_ok else F, f"zyklisch nur bei Drehzentrum {facts['F2_cyclic_centres']}; Lesart: keine Marke ausgezeichnet"),
    "K4": (P, "KMS in beiden Zyklen und invertierbar: (A,A)"),
    "K5": (P, f"c3 = 2/(4 pi |mu4|) = {facts['F3_c3']:.12f}"),
    "K6": (P, f"Kragenbreite l = L/2 erzwungen (F4: {facts['F4_width']})"),
    "K8": (P, "aus K1 bis K3 über P2.PARTITION.01 (TFPT, exakt)"),
}
H = {}
H["H1"] = dict(seam_common, **{
    "K1": (IN, "Blattinvolution und flache Naht vorausgesetzt"),
    "K7": (F, "(A,A) hat 0 Nullmoden; flacher Sektor nur als ungekoppelte direkte Summe"),
    "K9": (O, "drei Familien aus H^1 der Marken, aber kein Mechanismus für gleiche Kopien"),
    "K10": (P if cyc_ok else F, "H^1 der zyklisch permutierten Marken: {-1, i, -i}"),
    "K11": (IN, "E8 nur über den TFPT Compiler"),
    "K12": (F, "Naht plus Normale: 2+1"),
    "K13": (F, "keine 4D Geometrie"),
    "K14": (P, "keine Konstante gewählt"),
})
H2_K9 = None
H["H2"] = dict(seam_common, **{
    "K1": (P, "Pi-Drehung (-1 in den Einheiten) gibt zwei Blätter, vier 2-Torsionspunkte"),
    "K7": (F, "E8 als Gitter im Inneren, kein flacher Nullmodensektor der Naht"),
    "K9": (P, f"genau {facts['F7_gaussian_lines']} Gaußsche Geraden, sigma permutiert i->j->k: {facts['F7_sigma_cycles_ijk']}"),
    "K10": (F, f"Vierteldrehung der Naht auf den Geraden: i, -j, -k, also {{1, -1, -1}} statt {{-1, i, -i}}"),
    "K11": (P, f"Hurwitz^2 mod (1+i): {facts['F9_E8_counts'][2]} Vektoren der Norm 2, {facts['F9_E8_counts'][4]} der Norm 4"),
    "K12": (LES if facts["F10_octonion_associator_nonzero"] else O,
            "Lesart Raumzeit = größte assoziative Divisionsalgebra H (Frobenius; Oktonionen nicht assoziativ: geprüft), Re + Im = 1 + 3"),
    "K13": (O, "4D lokal; Einstein mit 1/c3 nur über TFPT GRAV.NONLINEAR.01 (bedingt)"),
    "K14": (P, "keine Konstante; Gitterskala setzt nur die Einheit"),
})
H["H2"]["K6"] = (P if facts["F8_hurwitz_spin_left"] else F,
                 f"4D: einheiteninvariante invertierbare Spinstrukturen links {facts['F8_hurwitz_spin_left']}, "
                 f"links+rechts {facts['F8_hurwitz_spin_left_right']}")
H["H3"] = {
    "K1": (IN, "Naht nicht Teil des Prinzips"), "K2": (IN, "Quadrat gesetzt"), "K3": (IN, "-"),
    "K4": (P if facts["F8_hypercubic_AAAA_invariant"] else F, "(A,A,A,A) unter den Würfeldrehungen invariant"),
    "K5": (O, "c3 nicht aus dem Prinzip; Spektralwirkung hat freie Momente"),
    "K6": (O, "KMS möglich, Breite nicht festgelegt"),
    "K7": (F, "Tensorprodukt: Spektrum sqrt(lambda^2 + mu^2), mit invertierbarem D_geom keine Nullmoden"),
    "K8": (IN, "Träger in D_F gesetzt"), "K9": (IN, "drei Kopien in D_F gesetzt"), "K10": (IN, "-"),
    "K11": (IN, "E8 nicht aus dem Prinzip"), "K12": (IN, "d = 4 gesetzt"),
    "K13": (P, "Einstein-Hilbert aus dem Wärmekoeffizienten a2 (Standard), Koeffizient frei"),
    "K14": (F, "Abschneidefunktion und ihre Momente frei"),
}
H["H4"] = dict(seam_common, **{
    "K1": (IN, "Blattinvolution vorausgesetzt"),
    "K7": (LES, f"Lesart GSO-Summe: (P,P) Sektor {facts['F6_E8_from_GSO']['complex_zero_modes_PP']} komplexe Nullmoden, flach für jede Breite; "
              f"{facts['F6_E8_from_GSO']['so16_bilinears']} + {facts['F6_E8_from_GSO']['ramond_chiral_ground_states']} = 248"),
    "K9": (O, "keine Familienvervielfachung"),
    "K10": (P if cyc_ok else F, "über H^1 der Marken"),
    "K11": (P, f"E8 Stufe 1 = GSO-Summe, S-invariant: {facts['F6_ZE8_S_invariant']}, Z(i) = {facts['F6_ZE8_at_i']:.6f}"),
    "K12": (F, "2D Naht mit 3D Volumen"),
    "K13": (F, "keine 4D Geometrie"),
    "K14": (P, "16 Fermionen = Rang von E8 mal 2, nicht gewählt"),
})
H["H4"]["K5"] = (P, "c3 = 1/(8 pi) pro Diracfeld (Normierung je Feld)")
H["H6"] = dict(seam_common, **{
    "K1": (P, "Pi-Drehung (-1 in Q8) gibt zwei Blätter und vier 2-Torsionspunkte der Naht"),
    "K6": (P if facts["F8b_lipschitz_AAAA_invariant"] else F,
           f"4D Lipschitz-Gitter: Q8 transitiv auf den Achsen {facts['F8b_lipschitz_axis_transitive']}, "
           f"(A,A,A,A) invariant und invertierbar {facts['F8b_lipschitz_AAAA_invariant']}; Breite L/2 aus 2D"),
    "K7": (LES, "Lesart B = Naht-Dirac summiert über Spinstrukturen (GSO): 8 komplexe Nullmoden im (P,P) Sektor, "
                "flach für jede Breite; 120 + 128 = 248"),
    "K9": (LES if facts["F8b_sigma_preserves_lipschitz_and_AAAA"] and facts["F11_sigma_on_selfdual_order3"] else F,
           f"Lesart Familien = drei selbstduale 2-Formen (= {facts['F8b_lipschitz_gaussian_lines']} Gaußsche Geraden); sigma "
           f"permutiert sie zyklisch (Ordnung 3: {facts['F11_sigma_on_selfdual_order3']}) und erhält Gitter und Spinstruktur"),
    "K10": (LES, f"Familien = selbstduale 2-Formen: Vierteldrehung {facts['F11_quarter_on_selfdual']}; mit dem Vorzeichen-"
                 f"charakter von q (wie bei TFPTs Trägermoden) {facts['F11_quarter_on_selfdual_times_q_character']} = {{-1, i, -i}}"),
    "K11": (P, f"E8 aus Z^8 = Lipschitz^2: {facts['F9b_E8_from_Z8']}; dasselbe Gitter wie die GSO-Naht"),
    "K12": (LES, "Lesart Raumzeit = größte assoziative Divisionsalgebra H (Frobenius), Re + Im = 1 + 3; "
                 "Lipschitz ist ihr Gaußsches Gitter"),
    "K13": (O, "4D lokal; Einstein mit 1/c3 nur über TFPT GRAV.NONLINEAR.01 (bedingt)"),
    "K14": (P, "Lipschitz statt Hurwitz durch K6 erzwungen, 16 Fermionen durch E8"),
})
prem = {"A": "keine ausgezeichneten Zyklen (zeitneutral und basisunabhängig)",
        "B": "die Naht ist eine Kodimension-2-Fläche mit gleichwertiger Normalebene",
        "C": "Familien = chirale Formen mittlerer Dimension der Raumzeit",
        "D": "Gravitation als Zustandsgleichung (Jacobson, wie TFPT)"}
f12 = facts["F12_gso_modular"]
H["H6h"] = dict(H["H6"])
H["H6h"]["K3"] = (P if facts["F13_centre_forced"] else F,
                  f"[C] H^1 der Marken muss als mu4-Modul gleich Lambda^2_+ x sign sein: zyklisch {facts['F13_H1_cyclic']} "
                  f"= {facts['F13_selfdual_times_sign']}, Zentrum auf Marke {facts['F13_H1_centre_at_mark']}: Viertelpunkt erzwungen")
H["H6h"]["K7"] = (P, f"[A] Basisunabhängigkeit: GSO-Summe nur für N = 0 mod 16 modular {f12}; kleinstes N = "
                     f"{facts['F12_smallest_N']}, (P,P) mit 8 komplexen Nullmoden, flach für jede Breite")
H["H6h"]["K9"] = (P if facts["F13_family_group_order"] == 24 else F,
                  f"[C] drei chirale 2-Formen; <sigma, Vierteldrehung> hat Ordnung {facts['F13_family_group_order']}, "
                  f"Charakter = Standard x Signum von S4: {facts['F13_is_S4_standard_times_sign']} (TFPT v117: signgetwistete Standarddarstellung)")
H["H6h"]["K10"] = (P if facts["F13_is_S4_standard_times_sign"] else F,
                   "[C] H^1 der Marken = Standarddarstellung von S4 = Lambda^2_+ x sign (Satz); der Vorzeichencharakter von q "
                   "ist das Signum von S4, nicht gewählt; Charaktere {-1, i, -i}")
H["H6h"]["K11"] = (P, f"[A] E8 = D8+ aus dem kleinsten modularen N = 16, dasselbe Gitter wie Z^8 mit Halbzahlkleber")
H["H6h"]["K12"] = (P, f"[B, C] d = 4: drei chirale Mittelformen nur in d = 4 {facts['F14_chiral_middle_forms']}; Naht plus "
                      "gleichwertige Normalebene = Z[i] + Z[i] j = Lipschitz; Re + Im = 1 + 3")
H["H6h"]["K13"] = (LES, f"[D] TFPTs Jacobson-Herleitung, alle Eingaben aus H6: {facts['F15_gravity_inputs']}; "
                        "c = 8 der Naht ist abgeleitet (N = 16), nicht gesetzt")
H["H6h"]["K14"] = (P, "keine Konstante; N = 16, d = 4, mu4 erzwungen; die Gitterskala ist die Einheit")
H["H5"] = dict(H["H2"])
H["H5"]["K7"] = H["H4"]["K7"]
H["H5"]["K11"] = (P, "dasselbe E8 zweimal: GSO-Gitter der Naht = Hurwitz^2 mod (1+i)")

# K15 gauge dynamics, K16 continuum limit: not part of the first catalogue, added for honesty
k15_k16 = {
    "H1": ((IN, "Eichgruppe nur über TFPT Compiler"), (O, "Gitter- und Kontinuumsfrage nicht gestellt")),
    "H2": ((O, "Erhaltungsgruppe (SU3 x SU2 x U1)/Z6 aus 3+2 [TFPT E], Eichfelddynamik nicht hergeleitet"), (O, "Hurwitzgitter, kein Skalenlimes")),
    "H3": ((P, "Spektralwirkung gibt Yang-Mills mit freien Kopplungen"), (O, "Kontinuum gesetzt")),
    "H4": ((O, "E8 Ströme der Naht, keine 4D Eichfelder"), (O, "nur 2D")),
    "H5": ((O, "wie H2"), (O, "wie H2")),
    "H6": ((O, "Erhaltungsgruppe (SU3 x SU2 x U1)/Z6 aus 3+2 [TFPT E]; Eichfelddynamik in 4D nicht hergeleitet"),
           (O, "Lipschitzgitter in Planckeinheiten; freier Skalenlimes, wechselwirkender Limes nicht gezeigt")),
}
k15_k16["H6h"] = k15_k16["H6"]
for h in H:
    H[h]["K15"], H[h]["K16"] = k15_k16[h]
y1 = json.loads(Path(__file__).with_name("y1_ergebnis.json").read_text()) if Path(__file__).with_name("y1_ergebnis.json").exists() else None
if y1:
    H["H6h2"] = dict(H["H6h"])
    H["H6h2"]["K12"] = (P, "[A, C] d = 4 aus drei chiralen Mittelformen; Kodimension 2 der 2D Naht; A auf allen Zyklen "
                           "macht die Normalebene gleichwertig, also Z[i] + Z[i] j = Lipschitz (Prämisse B entfällt)")
    ok15 = y1["T1_link_covariant"] and y1["T2_quadratic"] and y1["T2_positive"] and y1["T2_gauge_invariant"]
    H["H6h2"]["K15"] = (P if ok15 else F,
                        f"[E] Verbindung aus dem lokalen 3+2-Split eichkovariant; Fermionschleifen induzieren F^2: quadratisch "
                        f"{y1['T2_quadratic']}, positiv {y1['T2_positive']}, eichinvariant {y1['T2_gauge_invariant']}, "
                        f"Lauf {y1['T2_log_running_slope']:.3f} (Kontinuum {y1['T2_slope_continuum_16_dirac']:.3f}); "
                        f"Dynkin {y1['T3_dynkin']} -> sin^2 theta_W = 3/8 an der Gitterskala")
    H["H6h2"]["K16"] = (O, "freier Sektor hat logarithmisch laufende Kopplung (renormierbare Struktur); wechselwirkender "
                           "Kontinuumslimes nicht gezeigt")
    facts["F17_family_mixing_from_S4"] = {"TBM_in_S4": y1["T5_TBM_in_S4"], **y1["T5_compare"]}

facts["F16_seam_framing_anomaly"] = "N = 16 chiral: T-Phase exp(-2 pi i/3) = exp(-2 pi i c/24); verlangt Zufluss aus einem 3D Volumen"

# K17 framing anomaly of the seam, K18 neutrino mixing against data (added after the external review of 23.11.)
for h in H:
    H[h]["K17"] = (O, "T-Phase der Naht nicht behandelt") if h.startswith(("H4", "H5", "H6")) else (O, "keine chirale Naht mit GSO")
    H[h]["K18"] = (O, "keine Familiengruppe")
    H[h]["K19"] = (IN, "Existenz der Naht (geschlossene 2D Fläche mit chiralen Fermionen) als Ausgangspunkt gesetzt")
for h in ("H6", "H6h", "H6h2"):
    if h in H:
        H[h]["K18"] = (LES, "S4 gibt tribimaximal: sin^2 theta12 = 1/3, sin^2 theta23 = 1/2, sin^2 theta13 = 0 gegen 0,022")
z1p, z2p = Path(__file__).with_name("z1_ergebnis.json"), Path(__file__).with_name("z2_ergebnis.json")
if "H6h2" in H and z1p.exists() and z2p.exists():
    z1 = json.loads(z1p.read_text())
    z2 = json.loads(z2p.read_text())
    ok = z1["checks"]
    phi0 = 1 / (6 * np.pi) + 3 / (256 * np.pi ** 4)
    trY2 = 3 * Fr(1, 3) ** 2 + 2 * Fr(1, 2) ** 2
    s13 = phi0 * np.exp(-float(trY2))
    H["H6h3"] = dict(H["H6h2"])
    H["H6h3"]["K2"] = (P if ok["Drehung >2 nur bei i und rho"] and ok["Sechseck hält nur PP"] else F,
                       "[A] Fixpunkt, nicht Kovarianz (Z_AA ist auch bei 2i kovariant, die Rekonstruktionen sind dort "
                       "verschiedene Systeme); Drehung der Ordnung > 2 nur bei i und rho; rho verliert (A,A)")
    H["H6h3"]["K5"] = (P if ok["g~^4 = -1"] and ok["c3 je Spinorsektor"] else F,
                       f"[A] gelifte Vierteldrehung g~^4 = -1; jeder der vier spinoriellen Sektoren trägt exakt 1/4 jedes "
                       f"Eigenraums; c3 je Sektor = {z1['c3_per_spinor_sector']:.12f} = 1/(8 pi)")
    H["H6h3"]["K9"] = (P, "[C'] drei Familien = drei gerade Spinstrukturen der Naht (nichttriviale Charaktere der Markengruppe "
                          "E[2]); dieselbe S4-Darstellung wie H^1 der Marken")
    H["H6h3"]["K10"] = (P if ok["Familiencharaktere aus g"] and ok["H1 der Marken gleich"] else F,
                        f"[Naht] g = iz + 1/2 auf den drei geraden Spinstrukturen: Verschiebungsteil wirkt über die Weilpaarung, "
                        f"Spektrum {z1['family_spectrum_from_weil']}; kein eingesetztes Vorzeichen; Drehung um 0 gäbe {{1, 1, -1}}")
    H["H6h3"]["K11"] = (P, "[A, F] E8 = kleinstes N mit T-Kovarianz (N = 16); strenge Invarianz erst bei N = 48")
    H["H6h3"]["K12"] = (P if ok["rot = Lambda2 std"] else F,
                        "[A, C'] chirale 2-Formen = Lambda^2(Familien) (Charaktere für alle 24 Elemente gleich); "
                        "drei chirale Mittelformen nur in d = 4; Lipschitz; Re + Im = 1 + 3 (euklidisch, Lorentz per OS)")
    H["H6h3"]["K13"] = (LES, f"[D] Jacobson wie TFPT. Gegentest: kappa = {z2['kappa_per_dirac_taste']:.4f} je Dirac-Geschmack und "
                             f"Gitterfläche; drei Generationen geben a/l_P >= {z2['a_over_lP_three_generations_no_doublers']:.2f}; "
                             "die Gitterskala ist nicht die Planckskala")
    H["H6h3"]["K14"] = (P, "keine freie Konstante; a/l_P = 2 sqrt(kappa_gesamt) ist berechenbar, aber noch nicht berechnet "
                           "(Eichfelder, Kontaktterm)")
    H["H6h3"]["K16"] = (O, "freier Sektor läuft logarithmisch; wechselwirkender Limes und Wiederherstellung von SO(4) "
                           "nicht gezeigt; Gitterabstand mindestens drei Plancklängen")
    H["H6h3"]["K17"] = (LES, "[F] T-Phase exp(-2 pi i/3) kompensiert durch invertierbare 3D E8 Phase (Rand c = 8, T^2-Hilbertraum "
                             "eindimensional, bricht die Vierteldrehung nicht); Alternative drei E8 Blöcke, c = 24")
    H["H6h3"]["K18"] = (LES, f"S4 gibt TBM als führende Ordnung; TFPT-Kanäle: sin^2 theta12 = {1 / 3 - phi0 / 2:.4f}, "
                             f"sin^2 theta13 = phi0 exp(-tr Y^2) = {s13:.4f} mit tr Y^2 = {trY2} aus dem 3+2-Split (E), "
                             "etwa 2 sigma gegen Daten (JUNO entscheidet); delta = 240 Grad = Argument der T-Phase, nicht hergeleitet")
    facts["F18_review_checks"] = ok
    facts["F19_entanglement_density"] = z2
    z3p = Path(__file__).with_name("z3_ergebnis.json")
    if z3p.exists():
        z3 = json.loads(z3p.read_text())
        ok3 = z3["checks"]
        H["H6h4"] = dict(H["H6h3"])
        H["H6h4"]["K2"] = (P if ok3["nur tau = i stationär"] and ok3["Energie null bei i"] else F,
                           f"[A*] stationäre Nahtform <T> = 0: d log chi_E8/d tau = -(2 pi i/3) E6/E4, einzige Nullstelle im "
                           f"Fundamentalbereich {z3['e6_zeros']}; Gesamtenergie dort null, beide Rekonstruktionen gleich; "
                           f"ohne Orientierungsvoraussetzung; Gegenproben Boson {z3['boson_critical']}, Leech {z3['leech_critical']}")
        H["H6h4"]["K14"] = (P, "keine freie Konstante; die Form der Naht ist nicht gewählt, sondern ihre Bewegungsgleichung; "
                               "a/l_P = 2 sqrt(kappa_gesamt) berechenbar, noch nicht berechnet")
        facts["F20_stress_free_selection"] = z3
        z4p = Path(__file__).with_name("z4_ergebnis.json")
        if z4p.exists():
            z4 = json.loads(z4p.read_text())
            st, why = H["H6h4"]["K2"]
            stable = all(z4["checks"].values())
            H["H6h4"]["K2"] = (st if stable else F,
                               why + f"; stabil: Verdrillung null ist Maximum von |Z| (F_xx = {z4['Fxx']:+.3f}), längs der "
                                     f"Temperatur ist die Krümmung -(2 pi)^2 Var = {z4['Fyy']:+.3f} (positive Wärmekapazität), "
                                     f"mikrokanonisch S(E) konkav mit E = 0 bei beta = L")
            facts["F21_stability"] = z4
        z5p, z6p = Path(__file__).with_name("z5_ergebnis.json"), Path(__file__).with_name("z6_ergebnis.json")
        if z5p.exists() and z6p.exists():
            z5 = json.loads(z5p.read_text())
            z6 = json.loads(z6p.read_text())
            ok5, ok6 = all(z5["checks"].values()), all(z6["checks"].values())
            H["H6h5"] = dict(H["H6h4"])
            H["H6h5"]["K17"] = (LES if ok5 else F,
                                "[A*, K1] c = 8 mit F entschieden: Operatorlesart (c = 24) widerspricht dem KMS-Zustand (K6) und "
                                "TFPTs Familien in einem E8; ein E8 trägt genau drei Spinstruktursektoren (S3); die 3D Wand ist "
                                "in 4D genau bei c in 8Z unsichtbar, also Lesart: E8 Phase auf dem Verzweigungsschnitt der zwei Blätter")
            iso = z6["pi33_isotropic"]
            H["H6h5"]["K16"] = (LES if ok6 else F,
                                f"[physisches Gitter] Universalität im Infraroten perturbativ: c = 1 für jede Art (Vierteldrehung "
                                f"zwischen Zeit- und Raumachse), Vakuumpolarisation Zeit = Raum ({iso['Zeitachse'][0]:+.6f}), "
                                f"Richtungsabweichung ~ q^2 (Dimension 6); ausgezeichnete Zeit gäbe c_gamma^2 = "
                                f"{z6['c_gamma2_two_species']:.3f} (Dimension 4); OS über T^2; U(1) Landau-Pol "
                                f"{z6['landau_pole_GeV']:.1e} GeV. Nicht gezeigt: nichtperturbativ, chirale Gitterformulierung")
            facts["F22_fork"] = z5
            facts["F23_continuum"] = z6
            z7p = Path(__file__).with_name("z7_ergebnis.json")
            if z7p.exists():
                z7 = json.loads(z7p.read_text())
                ok7 = all(z7["checks"].values())
                H["H6h6"] = dict(H["H6h5"])
                H["H6h6"]["K19"] = (LES if ok7 else F,
                                    "[G] Naht = verdrillter Defekt eines komplexen Paarungsmassenfelds (J: psi -> chi -> -psi, "
                                    "Phase windet um pi): Kodimension 2; geschlossen und flach heißt Torus; J^2 = -1 gibt zwei "
                                    "Blätter auf Bosonischem und mu4 auf Fermionen; der Schnitt ist der Phasensprung (Eichwahl); "
                                    "Gitterindex: eine chirale Majorana-Weyl Mode je Paar mit Lichtgeschwindigkeit (robust), "
                                    "16 Paare geben c = 8; die Verdoppler liefern sie nicht. Offen: Identität der 16 Paare")
                facts["F24_seam_origin"] = z7

order = [f"K{i}" for i in range(1, 20)]
print("Gesamttabelle (P = hergeleitet, L = unter benannter Lesart, F = FAIL, O = OFFEN, E = EINGABE)\n")
print("      " + " ".join(f"{k:>4}" for k in order) + "     P    L    F    O    E")
for h, crit in H.items():
    row = [crit[k][0][0] for k in order]
    print(f"   {h} " + " ".join(f"{c:>4}" for c in row) + "  " + " ".join(f"{row.count(x):>4}" for x in "PLFOE"))
print()
for h, crit in H.items():
    print(f"{h}:")
    for k in order:
        st, why = crit[k]
        if st != P:
            print(f"   {k} {st}: {why}")
print("\nPrämissen von H6h:")
for k, v in prem.items():
    print(f"   {k}: {v}")
if "H6h2" in H:
    print("Prämissen von H6h2: A, C, D wie oben; B entfällt (folgt aus A, C und der 2D Naht); "
          "E: der 3+2-Split ist lokal (Verbindung auf den Gitterkanten)")
if "H6h3" in H:
    print("Prämissen von H6h3: A (der orientierte Tausch der Zyklen ist eine Symmetrie, Fixpunkt statt Kovarianz); "
          "C' (Familien = gerade Spinstrukturen der Naht, chirale 2-Formen = ihr Lambda^2); D; E; "
          "F (Rahmungsanomalie durch invertierbare 3D E8 Phase kompensiert)")
if "H6h4" in H:
    print("Prämissen von H6h4: A* (die Naht hat keine Hintergrundgeometrie: Basiswechsel sind Symmetrien, und ihre Form "
          "erfüllt die Bewegungsgleichung <T> = 0) ersetzt A; C', D, E, F wie H6h3. Der Fixpunkt, die Zeitneutralität "
          "und die Orientierung folgen.")
print("\nFakten:")
for k, v in facts.items():
    print(f"   {k}: {v}")
Path(__file__).with_name("v1_ergebnis.json").write_text(json.dumps(
    {"facts": {k: (v if not isinstance(v, (np.floating,)) else float(v)) for k, v in facts.items()},
     "scorecard": {h: {k: list(v) for k, v in c.items()} for h, c in H.items()}}, ensure_ascii=False, indent=1, default=str))
