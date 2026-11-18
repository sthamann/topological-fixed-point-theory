"""Deterministic probes of three TFPT moving-cut hypotheses."""

import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
from scipy.integrate import quad
import sympy as sp


def load_original(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def collar_profile(w):
    r = (w - .5) / .25
    return np.exp(1 - 1 / (1 - r*r)) if abs(r) < 1 else 0.


def translation_test(pin):
    # Constant normal translations fail the Pin deck relation. Test a
    # deck-compatible flow instead: v(w)=sin(pi*w), zero at both ends.
    _, pi = sp.symbols("w"), sp.pi
    k = np.pi / 2
    modes = np.arange(-1, 2)
    sea = pin.old_covariance(modes, k)
    tangent = np.array([np.pi / 4, 0., np.pi / 4])
    leakage = np.vdot(tangent, sea @ tangent).real
    exact = pi**2 * (1 - 1 / sp.sqrt(5)) / 16
    assert abs(leakage - float(exact)) < 2e-15
    assert leakage > 0
    return {"same_vacuum_leakage_norm_squared": leakage,
            "exact": str(exact)}


def current_test(pin):
    # Unit current, q=i. Extend the original even collar bump with a
    # longitudinal harmonic that is not in the Möbius algebra.
    mean, _ = quad(collar_profile, .25, .75, epsabs=1e-13)
    alpha = .25
    x, y = 1/16, 9/16  # centers of the original disjoint arc probes
    s = x - y
    c = np.exp(2j*np.pi*(alpha-1)*s) / (1 - np.exp(-2j*np.pi*s))
    ks = 2*np.pi*(np.arange(-4, 4) + alpha)
    occ = (ks < 0).astype(float)
    dc = 0j
    nonzero_entries = 0
    for i, ki in enumerate(ks):
        for j, kj in enumerate(ks):
            if abs(i-j) == 2 and occ[i] != occ[j]:
                # Divided difference of the negative spectral projector
                # on the original reducing zero branch at epsilon=0.
                v = -mean*(ki+kj)/4
                dc += (occ[i]-occ[j])/(ki-kj)*v*np.exp(1j*(ki*x-kj*y))
                nonzero_entries += 1
    ward = -mean/(2*np.sqrt(2))
    spectral = -2*np.real(c.conjugate()*dc)
    assert nonzero_entries == 4
    assert abs(ward-spectral) < 2e-15 and ward != 0

    # The original normal-only bump vanishes near all boundary collars.
    # Its mixed, omitted response is visible already at order epsilon^2.
    m = np.arange(1, 65)
    hm = np.array([2*quad(lambda w: collar_profile(w)*np.cos(np.pi*j*w),
                          .25, .75, epsabs=1e-13)[0] for j in m])

    def positive_mode_gain(k):
        energy = np.hypot(k, np.pi*m)
        return float(np.sum(hm**2/16*((energy-k)/(energy+k))**2*(1-k/energy)))

    k = np.pi/2
    gain = positive_mode_gain(k)
    # All paired normal modes give a decreasing gain in |k|. For the
    # half-circle pair the current correction is an alternating sum of
    # two gains. Its first pair dominates its second.
    j = np.arange(512)
    pairs = np.array([positive_mode_gain(2*np.pi*(i+alpha))
                      + positive_mode_gain(2*np.pi*(i+1-alpha)) for i in j])
    response = float(np.sum((-1.)**j*pairs))
    assert np.all(np.diff(pairs) < 0) and response > 0
    direct = {}
    for eps in [.01, .005, .0025]:
        modes, d, _, _ = pin.dirac_block(64, eps, k)
        eigenvalues, eigenvectors = np.linalg.eigh(d)
        direct[str(eps)] = float(np.sum(abs(eigenvectors[64, eigenvalues < 0])**2)/eps**2)
    assert abs(direct["0.0025"]-gain) < 2e-6
    return {"unit_current_first_order": {"ward": ward, "spectral": spectral,
                                         "absolute_error": abs(ward-spectral)},
            "normal_bump_current_order_2_estimate": response,
            "normal_bump_cutoff_longitudinal_omission_bound":
                positive_mode_gain(2*np.pi*(512+alpha)) +
                positive_mode_gain(2*np.pi*(513-alpha)),
            "normal_bump_first_two_pair_difference": float(pairs[0]-pairs[1]),
            "normal_bump_positive_mode_order_2": gain,
            "normal_bump_direct_occupation_over_epsilon_squared": direct}


def refinement_test(pin, v17):
    k = np.pi/2
    times = np.linspace(0, 4, 401)
    low, high = [pin.run(n, .1, k, times) for n in (64, 128)]

    def curve_difference(key):
        return float(np.max(np.abs(np.asarray(low[key])-high[key])))

    # Exact spatial blocking on the original C6 Green component. Keep the
    # eliminated scalar determinant along with the effective kernel.
    a, b = sp.symbols("a b", positive=True)
    u = v17.U6_power(1)
    q = a*sp.eye(6) - b*(u+u.T)
    keep, out = [0, 2, 4], [1, 3, 5]
    s = q.extract(keep, keep) - q.extract(keep, out)*q.extract(out, out).inv()*q.extract(out, keep)
    expected = (a-b*b/a)*sp.eye(3) - (b*b/a)*sp.ones(3)
    assert (s-expected).applyfunc(sp.simplify) == sp.zeros(3)
    assert sp.factor(q.det()-a**3*s.det()) == 0
    assert sp.factor(q.inv()[0, 2]-s.inv()[0, 1]) == 0
    delta = ((794-7*sp.sqrt(9961))/2187)**sp.Rational(1, 6)
    rows = []
    for y in [sp.Integer(1), sp.Rational(1, 3), sp.Rational(2, 3)]:
        aa = y*y+delta*delta
        bb = y*delta
        rows.append({"y": str(y), "retained_G_02": float(s.inv()[0, 1].subs({a:aa,b:bb})),
                     "omitted_logdet_if_discarding_exterior": float(3*sp.log(aa)),
                     "omitted_delta_derivative": float(6*delta/aa)})
    return {"pin_modes": [low["mode_count"], high["mode_count"]],
            "pin_occupancy_curve_error": curve_difference("original_sea_zero_mode_occupancy"),
            "pin_survival_curve_error": curve_difference("zero_branch_survival"),
            "pin_peak_occupancy": max(high["original_sea_zero_mode_occupancy"]),
            "pin_memory_vs_full_error": [low["derived_memory_error"], high["derived_memory_error"]],
            "c6_three_sites_exact_schur_and_absolute_determinant": rows}


def main(root):
    pin = load_original(root / "experiments/theory-contracts/"
                        "original-admissible-source-closure-20261025/"
                        "geometric_source_return/source_geometric_machine.py", "pin_original")
    sys.path.insert(0, str(root / "verification"))
    v17 = load_original(root / "verification/v17_hexagonal_resolvent.py", "v17_original")
    print(json.dumps({"translation": translation_test(pin),
                      "current": current_test(pin),
                      "refinement": refinement_test(pin, v17)}, indent=2))


if __name__ == "__main__":
    main(Path(sys.argv[1]))
