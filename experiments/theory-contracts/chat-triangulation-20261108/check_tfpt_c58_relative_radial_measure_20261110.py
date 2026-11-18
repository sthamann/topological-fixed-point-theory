"""Finite relative subtraction and radial mass on the supplied C58 cap."""

import json
from pathlib import Path
import runpy
import sys

import numpy as np
from scipy.integrate import quad
from scipy.special import zeta


repo = Path(sys.argv[1]).resolve()
base = repo / "experiments/theory-contracts/chat-triangulation-20261108"
old = runpy.run_path(str(base / "check_tfpt_conditioned_r113_weights_20261110.py"))
source = old["source"]
majorana = source["majorana"]
higgs = source["inner_form"](np.array([0., 1.]))
ell = 1 / 8


def mass_squared(r):
    values = np.linalg.eigvalsh(majorana + r * higgs)
    assert np.max(abs(values + values[::-1])) < 1e-12
    return values[48:]**2


base_masses = mass_squared(0)
native_quadratic = np.trace(higgs @ higgs).real / 2
assert abs(np.trace(majorana @ higgs)) < 1e-14


def log_cap_finite_part(mass2, cutoff, include_tail=True):
    # All actual quarter momenta p=n+1/4, both signs.
    p = np.abs(np.arange(-cutoff, cutoff) + .25)[:, None]
    energy = np.sqrt(p*p + mass2)
    vacancy = mass2 / (2 * energy * (energy+p))
    excitation = np.log1p(vacancy * np.expm1(-2 * ell * energy))
    # Exact stable identity for E - p - m²/(2p).
    finite_bulk = -mass2**2 / (2 * p * (energy + p)**2)
    value = float(excitation.sum() + ell * finite_bulk.sum())
    if not include_tail:
        return value

    # Four terms of each convergent Taylor series at large p.
    def tail(order):
        return zeta(order, cutoff + .25) + zeta(order, cutoff + .75)

    m2, m4, m6, m8, m10 = [(mass2**k).sum() for k in range(1, 6)]
    value += (-tail(2)*m2/4 + 5*tail(4)*m4/32 - 11*tail(6)*m6/96 + 93*tail(8)*m8/1024
              + ell * (-tail(3)*m4/8 + tail(5)*m6/16 - 5*tail(7)*m8/128 + 7*tail(9)*m10/256))
    return value


def relative_log(radius, cutoff):
    masses = mass_squared(radius)
    assert abs(masses.sum() - base_masses.sum() - native_quadratic * radius**2) < 1e-10
    return log_cap_finite_part(masses, cutoff) - log_cap_finite_part(base_masses, cutoff)


rows = {str(r): relative_log(r, 1024) for r in (0., .1, .55, 1., 4., 8., 16.)}
assert rows["0.0"] == 0 and rows["16.0"] < rows["8.0"] < rows["4.0"] < -20
log_difference = max(abs(relative_log(r, 512) - v) for r, v in
                     ((float(k), value) for k,value in rows.items()))
assert log_difference < 1e-5


def integrated(c, cutoff):
    def density(r):
        return 2 * np.pi**2 * r**3 * np.exp(relative_log(r, cutoff) - c*r*r)
    mass, numerical_error = quad(density, 0, 8, epsabs=1e-12, epsrel=1e-10)
    return {"mass_on_radial_interval_0_8": mass, "quadrature_error_estimate": numerical_error,
            "cartesian_null_density_over_mass": 1/mass}


measure = {str(c): integrated(c, 512) for c in (-1., 0., 1.)}
comparison = integrated(0, 1024)["mass_on_radial_interval_0_8"]
assert abs(measure["0.0"]["mass_on_radial_interval_0_8"] - comparison) < 1e-10
assert measure["-1.0"]["mass_on_radial_interval_0_8"] > comparison > measure["1.0"]["mass_on_radial_interval_0_8"]

# Uniform tail estimate on r>=8. Keep only modes p<=m_min/4 from
# the largest mass; all discarded finite bulk and cap terms are <=0.
sigma = np.linalg.norm(higgs, 2)
m0 = np.linalg.norm(majorana, 2)
minimum_slope = sigma - m0/8
p_small = np.abs(np.arange(-512, 512) + .25)
p_small = p_small[p_small <= 2*minimum_slope]
decay = ell * minimum_slope**2 * np.sum(1/(4*p_small)) - 1
assert decay > 0
def base_tail(order):
    return zeta(order, 1024.25) + zeta(order, 1024.75)

# Use inequalities rather than the Taylor correction to bound the offset.
offset = (-log_cap_finite_part(base_masses, 1024, include_tail=False)
          + base_masses.sum() * base_tail(2) / 2
          + ell * (base_masses**2).sum() * base_tail(3) / 8)
tail_log_upper = (offset + np.log(np.pi**2) - 2*np.log(decay)
                  + np.log1p(64*decay) - 64*decay)
assert tail_log_upper < -186

print(json.dumps({"quarter": .25, "ell": ell, "exact_radial_quadratic_trace": native_quadratic,
                  "old_reference_higgs_zero_majorana_retained": True,
                  "finite_part_log_relative_to_r0": rows,
                  "maximum_log_512_1024_difference_on_test_radii": log_difference,
                  "conditional_counterterm_test": measure,
                  "log_tail_upper_bound_for_all_three_counterterms": float(tail_log_upper),
                  "mass_at_c0_cutoff_1024": comparison,
                  "scope": "Old constant-H C58 cap, declared one-log subtraction and d4H, no primitive original boson measure, loaded links or full P_adm."}, indent=2))
