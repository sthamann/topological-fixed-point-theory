"""Tests on the given R032 exterior fiber. Set PYTHONPATH=verification."""
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_ivp
from scipy.linalg import expm, expm_frechet

base = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("quarter", base / "check_quarter_reference.py")
source = importlib.util.module_from_spec(spec)
spec.loader.exec_module(source)
native = source.native
dy = source.higgs_matrix
charged = np.flatnonzero(np.abs(source.mass_matrix).sum(axis=0) == 0)
active = charged[np.abs(dy[np.ix_(charged,charged)]).sum(axis=0) > 0]
assert len(active) == 84
gamma = source.grading[np.ix_(active, active)]
identity = np.eye(len(active))
charges = np.array([
    (1 if i < len(native.P) else -1) if native.P[i % len(native.P)][0] in ("u", "d") else 0
    for i in active
])
center = np.diag(np.exp(2j*np.pi*charges/3))
twists = [identity, center, center @ center]


def field(h2):
    q = np.diag([np.conj(h2), h2])
    projection = native.rep(1, np.zeros((2, 2)), np.zeros((3, 3)))
    rep = native.rep(0, q, np.zeros((3, 3)))
    linear = rep @ (dy @ projection - projection @ dy)
    linear = linear + linear.conj().T
    swap = native.Kswap()
    linear = linear + swap @ linear.conj() @ swap
    excluded = np.setdiff1d(np.arange(len(dy)), active)
    assert np.max(np.abs(linear[np.ix_(active, excluded)])) < 1e-14
    return -1j * gamma @ linear[np.ix_(active, active)]


mark = field(1)
assert np.max(np.abs(mark + 1j*gamma@dy[np.ix_(active, active)])) < 1e-14
positive = np.diag(gamma).real > 0
negative = ~positive
B = float(np.linalg.norm(mark[np.ix_(positive, negative)])**2)
colored = charges != 0
Bc = float(np.linalg.norm(mark[np.ix_(positive & colored, negative & colored)])**2)


def contract(transfer, derivative=None):
    determinants, variations = [], []
    for twist in twists:
        matrix = identity + twist @ transfer
        determinant = np.linalg.det(matrix)
        determinants.append(determinant)
        if derivative is not None:
            variations.append(determinant * np.trace(np.linalg.solve(matrix, twist @ derivative)))
    value = np.mean(determinants)
    assert value.real > 0 and abs(value.imag/value.real) < 1e-12
    if derivative is None:
        return np.log(value.real)
    variation = np.mean(variations)/value
    assert abs(variation.imag) < 1e-11
    return float(variation.real)


def sewn(segments, cuts=()):
    transfer = identity.copy()
    derivative = np.zeros_like(identity)
    start = 0
    for end, h, jet in segments:
        local_cuts = [cut for cut in cuts if start < cut < end]
        for stop in local_cuts + [end]:
            duration = stop - start
            if jet is None:
                step = expm(-duration*h)
                derivative = step @ derivative
            else:
                step, dstep = expm_frechet(-duration*h, -duration*jet)
                derivative = step @ derivative + dstep @ transfer
            transfer = step @ transfer
            start = stop
    return transfer, derivative


def direct(segments, strength=0):
    transfer = identity.astype(complex)
    start = 0
    for end, h, jet in segments:
        if jet is not None:
            h = h + strength*jet

        def ode(_, state):
            return (-h @ state.reshape(identity.shape)).ravel()

        solved = solve_ivp(ode, (start, end), transfer.ravel(),
                           rtol=3e-13, atol=3e-15, method="DOP853")
        assert solved.success
        transfer = solved.y[:, -1].reshape(identity.shape)
        start = end
    return transfer


def cut_test():
    p, T = .25, .8
    endpoints = [.13, .30, .50, .67, .8]
    p_gamma = p*gamma
    middle = field(.43)
    gaps = (.071, .247, .431, .723)
    q = np.diag([-1j, 1j])
    left = native.rep(1, q, np.eye(3))
    swap = native.Kswap()
    weak_frame = (left @ swap @ left.conj() @ swap)[np.ix_(active, active)]
    np.testing.assert_allclose(weak_frame @ field(.43) @ weak_frame.conj().T, field(.43j), atol=1e-14)
    np.testing.assert_allclose(weak_frame @ center, center @ weak_frame, atol=1e-14)
    np.testing.assert_allclose(weak_frame @ gamma, gamma @ weak_frame, atol=1e-14)
    result = []
    for outside in (field(.43), field(.43j)):
        hs = [p_gamma + outside, p_gamma, p_gamma + middle, p_gamma, p_gamma + outside]
        base_segments = [(end, h, None) for end, h in zip(endpoints, hs)]
        # A symmetric late window carries Y. The early U window is inside the first slab.
        segments = base_segments[:2] + [(.36, hs[2], None), (.44, hs[2], middle), (.50, hs[2], None)] + base_segments[3:]
        transfer, derivative = sewn(segments, gaps)
        expected = direct(segments)
        first = expm(-.13*hs[0])
        final = expm(-.13*hs[-1])
        inside = np.linalg.solve(final, transfer) @ np.linalg.inv(first)
        new_first = weak_frame @ first @ weak_frame.conj().T
        new_final = weak_frame @ final @ weak_frame.conj().T
        transported = new_final @ weak_frame @ inside @ weak_frame.conj().T @ new_first
        expected_frame = weak_frame @ transfer @ weak_frame.conj().T
        np.testing.assert_allclose(transported, expected_frame, atol=2e-12)
        np.testing.assert_allclose(contract(transported), contract(transfer), atol=2e-12)
        np.testing.assert_allclose(transfer, expected, rtol=1e-11, atol=2e-12)
        y_sewn = -contract(transfer, derivative)/.08

        def direct_mean(step):
            return -(contract(direct(segments, step)) - contract(direct(segments, -step)))/(.16*step)

        y_direct = (4*direct_mean(.001) - direct_mean(.002))/3
        assert abs(y_direct - y_sewn) < 1e-7
        # Separate scalar factors before the fermionic trace. The even T lift
        # has the identical unmarked generator; S removes the even field.
        durations = np.diff([0] + endpoints)
        log_sea = sum(dt*np.linalg.eigvalsh(h)[:len(active)//2].sum()
                      for dt, h in zip(durations, hs))
        log_A = contract(transfer) + log_sea
        ref_S_segments = [(end, p_gamma, None) for end in endpoints]
        ref_S_open = direct(ref_S_segments)
        np.testing.assert_allclose(ref_S_open, sewn(ref_S_segments, gaps)[0], atol=2e-12)
        ref_S = contract(ref_S_open) - T*p*len(active)/2
        ref_T_open = direct(base_segments)
        np.testing.assert_allclose(ref_T_open, transfer, atol=2e-12)
        ref_T = contract(ref_T_open) + log_sea
        assert abs(log_A - ref_T) < 1e-11
        U = .43**2
        y_with_early_U = -contract(transfer, U*derivative)/.08
        connected = y_with_early_U - U*y_sewn
        assert abs(connected) < 1e-13
        result.append(dict(early_U=U, log_fermion_A_over_S=log_A-ref_S,
                           log_fermion_A_over_T=log_A-ref_T,
                           late_Y_sewn=y_sewn, late_Y_direct=y_direct,
                           open_matrix_error=float(np.max(np.abs(transfer-expected))),
                           transport_error=float(np.max(np.abs(transported-expected_frame))),
                           fixed_history_connected_UY=connected))
    assert abs(result[0]["late_Y_sewn"]-result[1]["late_Y_sewn"]) > .1
    return result


def geometry_test(p):
    T = .8
    # Reflection-even disjoint profiles, including gaps where both vanish.
    ends = [.14, .27, .53, .66, .8]
    profiles = [(1, 0), (0, 0), (0, 1), (0, 0), (1, 0)]
    q, w = leggauss(28)

    def window(a, b):
        return (a+b)/2 + (b-a)*q/2, (b-a)*w/2

    t1, w1 = window(0, .14)
    t2, w2 = window(.66, .8)
    s, ws = window(.27, .53)
    distances = abs(np.r_[t1, t2][:, None] - s)
    weights = np.r_[w1, w2][:, None]*ws
    a = T - 2*distances
    pT = p*T
    denom = 1 + np.cosh(pT)
    integ = np.sum(weights*np.cosh(a*p))
    integ_prime = np.sum(weights*a*np.sinh(a*p))
    K = integ/denom
    K_prime = integ_prime/denom - K*T*np.sinh(pT)/denom
    s0 = 1/(2*denom)
    s_prime = -s0*T*np.tanh(pT/2)
    log_chi = 36*np.log1p(-3*s0)
    chi = np.exp(log_chi)
    alpha = 2*chi/(1+2*chi)
    alpha_prime = alpha*(1-alpha)*(-108*s_prime/(1-3*s0))
    correction = 3*Bc*alpha*s0/(1-3*s0)
    correction_prime = 3*Bc*(alpha_prime*s0/(1-3*s0)
                               + alpha*s_prime/(1-3*s0)**2)
    analytic_cross = K*(B+correction)
    analytic_radius = -p*(K_prime*(B+correction) + K*correction_prime)

    def geometry(x, z, radial_shift=0, insertion=False):
        momentum = p*np.exp(-radial_shift)
        parts = []
        for end, (f, g) in zip(ends, profiles):
            h = momentum*gamma + (x*f+z*g)*mark
            parts.append((end, h, -momentum*gamma if insertion else None))
        transfer, derivative = sewn(parts)
        return contract(transfer, derivative if insertion else None)

    def mixed(func, step):
        return sum(np.sign(x*z)*func(x,z) for x in (step,-step) for z in (step,-step))/(4*step**2)

    def mixed_extrap(func, step=.006):
        return (4*mixed(func, step/2) - mixed(func, step))/3

    observed_cross = mixed_extrap(geometry)
    inserted = mixed_extrap(lambda x,z: geometry(x,z,insertion=True))

    def direct_radius(step):
        plus = mixed_extrap(lambda x,z: geometry(x,z,step))
        minus = mixed_extrap(lambda x,z: geometry(x,z,-step))
        return (plus-minus)/(2*step)

    radius_fd = (4*direct_radius(.006) - direct_radius(.012))/3
    assert abs(analytic_cross-observed_cross) < 3e-7
    assert abs(analytic_radius-inserted) < 3e-7
    assert abs(analytic_radius-radius_fd) < 8e-6
    return dict(p=p, radius=1/(4*p), two_profile_cross=analytic_cross, radius_from_kernel=analytic_radius,
                radius_from_generator_insertion=inserted, radius_from_direct_perturbation=radius_fd,
                center_part_of_radius_score=-p*(K_prime*correction + K*correction_prime))


def infrared_test():
    assert np.sum(charges[negative]) % 3 == 0
    assert np.max(np.abs(center@mark-mark@center)) < 1e-13
    assert np.linalg.norm(mark[np.ix_(positive, positive)]) == 0
    assert np.linalg.norm(mark[np.ix_(negative, negative)]) == 0
    # All one particle-hole states produced by the centered mark are
    # neutral and have energy 2p. Exact vacuum spectral measure per r.
    np.testing.assert_allclose(B, np.linalg.norm(dy[np.ix_(active,active)])**2/2)
    q, w = leggauss(32)
    t = .07 + .02*q
    s = .33 + .04*q
    d = s[:,None] - t
    weights = np.outer(w/2,w/2)
    limit = B*np.sum(weights/d)
    rows = []
    for r in (4, 16, 64, 256, 1024, 4096):
        # Exact geometric sum of exp[-d(n+1/2)/r], including every mode.
        window_cross = B*np.sum(weights/(2*r*np.sinh(d/(2*r))))
        n = np.arange(r)
        static = B*np.sum(1/(n+.5))  # restricted to 0<E<=1
        eps = .01
        mass = B*np.count_nonzero((n+.5)/r <= eps)/r
        rows.append(dict(r=r, finite_windows=float(window_cross), static_low_energy=float(static),
                         low_energy_mass_over_eps=float(mass/eps)))
    slope = (rows[-1]["static_low_energy"]-rows[-2]["static_low_energy"])/np.log(4)
    assert abs(slope-B) < 1e-5
    assert abs(rows[-1]["finite_windows"]-limit) < 1e-6
    return dict(B=B, vacuum_center_trivial_on_marked_states=True,
                finite_windows_limit=float(limit), static_log_slope=float(slope), rows=rows)


print(json.dumps({"open_cut": cut_test(), "geometry": [geometry_test(p) for p in (.25, 5)], "infrared": infrared_test()}, indent=2))
