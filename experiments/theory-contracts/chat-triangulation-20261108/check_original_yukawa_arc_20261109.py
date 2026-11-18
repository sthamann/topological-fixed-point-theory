"""Read the existing R032 and NS weak segment through both Yukawa ends."""
import contextlib
import io
import json
import runpy
import sys
from pathlib import Path

import numpy as np

repo = Path(sys.argv[1]).resolve()
base = repo / "experiments/theory-contracts/chat-triangulation-20261108"
sys.path.insert(0, str(repo / "verification"))
import v252_full_finite_triple as v

path = base / "check_r113_dressed_open_20261109.py"
args = sys.argv
try:
    sys.argv = [str(path), str(base)]
    with contextlib.redirect_stdout(io.StringIO()):
        source = runpy.run_path(str(path))
finally:
    sys.argv = args

nx = np.array([[0, 1], [1, 0]])
eye2 = np.eye(2)
q = source["quaternion"]
identity = np.eye(96)
gamma = v.grading()
swap = v.Kswap()
up_l = [v.IDX[("nu", "L", 0, k)] for k in range(3)]
down_l = [v.IDX[("e", "L", 0, k)] for k in range(3)]
nu_r = [v.IDX[("nu", "R", 0, k)] for k in range(3)]
e_r = [v.IDX[("e", "R", 0, k)] for k in range(3)]
left, right = up_l + down_l, nu_r + e_r
sq = lambda a: float(np.vdot(a, a).real)
zpair = lambda z: [float(z.real), float(z.imag)]


def lift(lam, weak, color=None):
    if color is None:
        color = np.exp(-1j*np.angle(lam)/3) * np.eye(3)
    assert abs(np.linalg.det(color) * lam - 1) < 2e-13
    a = v.rep(lam, weak, color)
    return a @ swap @ a.conj() @ swap


def decompose(u_h, hi, hj):
    b = q(hj).conj().T @ u_h @ q(hi)
    d = np.linalg.det(u_h)
    c, x = b[1, 1], b[0, 1]
    assert np.max(abs(b - np.array([[d * c.conjugate(), x],
                                     [-d * x.conjugate(), c]]))) < 2e-13
    assert abs(abs(c)**2 + abs(x)**2 - 1) < 2e-13
    return c, x


def use_both_yukawa_legs(hi, hj, lam, weak):
    # Both actual original columns, not hand-written replacement vertices.
    vi = source["inner_form"](hi)[np.ix_(left, right)]
    vj = source["inner_form"](hj)[np.ix_(left, right)]
    bare = source["inner_form"](np.array([0, 1]))[np.ix_(left, right)]
    ell = lift(lam, weak)
    ell_l = ell[np.ix_(left, left)]
    ell_r = ell[np.ix_(right, right)]
    b = np.linalg.solve(bare.conj().T, vj.conj().T @ ell_l @ vi) @ np.linalg.inv(bare)
    d = np.kron(np.diag([1, lam.conjugate()**2]), np.eye(3))
    assert np.max(abs(ell_r - d)) < 2e-13
    # After cancelling only the unchanged family matrices and the endpoint
    # right charge, the electron and neutrino readouts are conjugate.
    raw = d.conj().T @ b
    c, x = decompose(lam * weak, hi, hj)
    expected = np.array([[c.conjugate(), lam.conjugate()**2 * x],
                         [-lam**2 * x.conjugate(), c]])
    assert np.max(abs(raw - np.kron(expected, np.eye(3)))) < 2e-13
    return c, x


hi, hj = [np.array(h) / np.linalg.norm(h) for h in source["history"][:2]]
p, qj = [np.outer(h, h.conj()) for h in (hi, hj)]
kato = (qj @ p + (eye2-qj) @ (eye2-p)) / np.sqrt(np.trace(qj @ p).real)
phase = q(hj) @ q(hi).conj().T
rows = {}
for label, weak in (("native", eye2), ("ray", kato), ("vector", phase)):
    c, x = use_both_yukawa_legs(hi, hj, 1, weak)
    # Independent arbitrary changes of weak and hypercharge frame at both
    # ends. The transported electron overlap and the conversion probability
    # are invariant. The open conversion coefficient carries charge.
    ai, aj = np.exp(.27j), np.exp(-.43j)
    si, sj = q(np.array([.6, .8j])), q(np.array([.8j, .6]))
    di, dj = ai * si, aj * sj
    lam_new = aj * ai.conjugate()
    u_new = dj @ weak @ di.conj().T / lam_new
    gc, gx = use_both_yukawa_legs(di @ hi, dj @ hj, lam_new, u_new)
    assert abs(gc - c) + abs(abs(gx) - abs(x)) < 2e-13
    rows[label] = {"electron": zpair(c), "neutrino": zpair(c.conjugate()),
                   "conversion_probability": abs(x)**2}

# Actual constant-Higgs v252 source, NS double, previous Fold arc and width.
D = v.full_D(np.diag([1,2,3]), np.diag([4,5,6]), np.diag([2,3,4]),
             np.diag([1,3,5]), np.array([[2,1,0],[1,3,1],[0,1,4]]))
gen = v.rep(0, 1j * nx, np.zeros((3,3)))
Q = -1j * (gen + swap @ gen.conj() @ swap)
G = lift(1j, 1j * np.diag([1, -1]), 1j * np.eye(3))
assert np.max(abs(G @ D @ G - D)) == 0.
assert np.max(abs(G @ Q @ G + Q)) == 0.
plus, minus = (identity + G)/2, (identity - G)/2
s, ell = 1/8, np.pi/4
h0 = -1j * gamma @ D
fixed_h = np.array([0, 1])
mode_range = 384
tail_atol = float(np.sqrt(96) * np.exp(-s * (mode_range - 1 - np.linalg.norm(D, 2)))
                  / (np.pi * (1 - np.exp(-s))))


def original_cross(alpha, a, generator=Q, mass=h0):
    x = np.zeros_like(D)
    # This range is exactly paired under m -> -m-1 between quarter charts.
    for m in range(-mode_range, mode_range):
        p = m + alpha
        ev, u = np.linalg.eigh(p * gamma + mass + a * generator @ gamma)
        y = np.exp(-s * np.abs(ev))
        x += (np.exp(-1j * p * ell) / (2*np.pi)) * ((u * (y / (1 + y*y))) @ u.conj().T)
    square = generator @ generator
    w = identity - square + np.cos(ell*a)*square - 1j*np.sin(ell*a)*generator
    return w @ x


def probe_readouts(k):
    keep = (plus @ k @ plus)[np.ix_(up_l, up_l)]
    cross = (minus @ k @ plus)[np.ix_(down_l, up_l)]
    return sq(keep), sq(cross)


ns = {}
crosses = {}
for alpha in (.25, .75):
    by_strength = {}
    for a in (0., .5, 1.):
        weak = np.cos(ell*a) * eye2 - 1j*np.sin(ell*a) * nx
        c, x = use_both_yukawa_legs(fixed_h, fixed_h, 1, weak)
        k = original_cross(alpha, a)
        km = original_cross(alpha, -a)
        crosses[alpha,a] = (k,km)
        even, odd = (k+km)/2, (k-km)/2
        assert np.linalg.norm(G @ k @ G-km) < 1e-12
        assert np.linalg.norm(plus @ odd @ plus) < 1e-12
        assert np.linalg.norm(minus @ even @ plus) < 1e-12
        keep, cross = probe_readouts(k)
        assert abs(probe_readouts(km)[1] - cross) < 1e-17
        assert cross < 1e-30 if a == 0 else np.sqrt(cross) > 1e5 * tail_atol
        by_strength[str(a)] = {
            "wilson_electron": zpair(c), "wilson_conversion_amplitude": zpair(x),
            "wilson_conversion_probability": abs(x)**2,
            "yukawa_legs_stripped_neutrino_to_neutrino_norm_squared": keep,
            "yukawa_legs_stripped_neutrino_to_electron_norm_squared": cross,
            "conversion_fraction_among_these_two_observed_left_channels": cross/(keep+cross),
            "full_cross_matrix_norm": float(np.linalg.norm(k)),
            "positive_odd_all_96_norm_squared": sq(minus @ k @ plus),
            "opposite_a_paired_odd_matrix_norm": float(np.linalg.norm(minus @ even @ plus)),
        }
    ns[str(alpha)] = by_strength
quarter_conjugacy = max(float(np.linalg.norm(crosses[.75,a][1] - crosses[.25,a][0].conj()))
                        for a in (0., .5, 1.))
assert quarter_conjugacy < 1e-12

weak_z = np.diag([1, -1])
tz = v.rep(0, 1j * weak_z, np.zeros((3, 3)))
qz = -1j * (tz + swap @ tz.conj() @ swap)
assert np.linalg.norm(G @ qz @ G-qz) == 0.
reference = {
    "free_no_mass_weak_loop_response_change": {},
    "longitudinal_full_source_positive_odd_norm_squared": {},
    "transverse_full_source_change_at_trivial_weak_loop": {},
}
for alpha in (.25, .75):
    free_unshifted = original_cross(alpha, 0., mass=np.zeros_like(D))
    free_shifted = original_cross(alpha, 1., mass=np.zeros_like(D))
    free_change = float(np.linalg.norm(free_unshifted-free_shifted))
    assert free_change < 1e-12
    longitudinal = original_cross(alpha, 1., generator=qz)
    long_switch = sq(minus @ longitudinal @ plus)
    assert long_switch < 1e-28
    reference["free_no_mass_weak_loop_response_change"][str(alpha)] = free_change
    reference["longitudinal_full_source_positive_odd_norm_squared"][str(alpha)] = long_switch
    reference["transverse_full_source_change_at_trivial_weak_loop"][str(alpha)] = float(
        np.linalg.norm(crosses[alpha,1.][0]-crosses[alpha,0.][0]))
pairs = {}
averages = []
for label, axis in (("transverse", nx), ("longitudinal", weak_z)):
    links = [np.cos(ell)*eye2 + sign*1j*np.sin(ell)*axis for sign in (-1,1)]
    averages.append(sum(links)/2)
    amplitudes = np.array([decompose(link, fixed_h, fixed_h) for link in links])
    mean = amplitudes[:,0].mean()
    switch = np.mean(abs(amplitudes[:,1])**2)
    variance = np.mean(abs(amplitudes[:,0]-mean)**2)
    assert abs(1-abs(mean)**2-switch-variance) < 2e-13
    assert all(np.linalg.norm(np.linalg.matrix_power(link,8)-eye2) < 2e-13
               for link in links)
    pairs[label] = {"mean_electron": zpair(mean), "mean_conversion": float(switch),
                    "phase_variance": float(variance)}
assert np.linalg.norm(averages[0]-averages[1]) < 1e-14

print(json.dumps({
    "scope": "two given finite charts and original NS double with existing formal weak probe; physical QTheta, Padm and original positive measure not assigned",
    "R032_first_segment": rows,
    "NS_Fold_same_two_endpoints": ns,
    "NS_quarter_complement_conjugacy_error": quarter_conjugacy,
    "same_mean_Wilson_link_and_same_weak_loop_but_different_second_words": pairs,
    "same_arc_original_free_and_longitudinal_controls": reference,
    "absolute_cutoff": mode_range,
    "NS_matrix_frobenius_tail_bound": tail_atol
},indent=2))
