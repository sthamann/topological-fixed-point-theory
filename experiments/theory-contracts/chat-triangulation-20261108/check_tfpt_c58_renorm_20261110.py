"""Compare finite sea prescriptions on the existing constant C58 source."""

import contextlib
import io
import json
import math
from pathlib import Path
import runpy
import sys

import numpy as np
from scipy.integrate import quad
from scipy.special import digamma


def compare(repo):
    path = (repo / 'experiments/theory-contracts/chat-triangulation-20261108'
            / 'check_tfpt_c58_relative_radial_measure_20261110.py')
    with contextlib.redirect_stdout(io.StringIO()):
        original = runpy.run_path(str(path))
    relative_log = original['relative_log']
    ell, coefficient = original['ell'], original['native_quadratic']
    gamma = float(np.euler_gamma)
    finite_quarter_sum = -float(digamma(.25) + digamma(.75)) / 2
    assert abs(finite_quarter_sum - (gamma + math.log(8))) < 1e-14

    schemes = {'modewise_subtracted': 0.,
               'actual_energy_Abel': finite_quarter_sum - 1 - gamma,
               'actual_energy_zeta': finite_quarter_sum - 1,
               'fixed_free_momentum_Abel': finite_quarter_sum - gamma,
               'fixed_free_momentum_zeta': finite_quarter_sum,
               'C59_AP_bulk_only_in_C58_cap_as_diagnostic': finite_quarter_sum - math.log(2)}
    largest_enhancement = ell * coefficient * max(schemes.values())
    # Concavity of sqrt(p²+x) bounds the sum of all 48 bulk terms by
    # 48 times the term at their mean m². Its ratio to m² decreases.
    momenta = np.abs(np.arange(-512, 512) + .25)
    mean_mass_squared_at_8 = (coefficient * 64 + original['base_masses'].sum()) / 48
    decrease = mean_mass_squared_at_8 / (2 * momenta *
                 (np.sqrt(momenta**2 + mean_mass_squared_at_8) + momenta)**2)
    decay = ell * coefficient * decrease.sum() - largest_enhancement
    assert decay > 0
    tail = (original['offset'] + math.log(math.pi**2) - 2 * math.log(decay)
            + math.log1p(64 * decay) - 64 * decay)
    assert tail < -100

    base_hessian = [2 * relative_log(step, 512) / step**2 for step in (.002, .001, .0005)]
    richardson = (4 * base_hessian[-1] - base_hessian[-2]) / 3
    assert abs(base_hessian[-1] - richardson) < 1e-3
    baseline = relative_log(.55, 512)
    result = {}
    for name, finite_constant in schemes.items():
        c = -ell * coefficient * finite_constant

        def density(r):
            return 2 * math.pi**2 * r**3 * math.exp(relative_log(r, 512) - c * r**2)

        mass, error = quad(density, 0, 8, epsabs=1e-13, epsrel=1e-10)
        moment, _ = quad(lambda r: density(r) * r**2, 0, 8, epsabs=1e-14, epsrel=1e-10)
        result[name] = {'c': c, 'quarter_sum_finite_constant': finite_constant,
                        'raw_log_at_r_0_55': baseline - c * .55**2,
                        'radial_mass_on_0_8': mass, 'quadrature_error': error,
                        'normalized_cartesian_null_density': 1 / mass,
                        'null_log_weight_radial_hessian': richardson - 2 * c,
                        'mean_r_squared': moment / mass}
    return {'scope': 'Same supplied C58 constant Higgs cap, old r=0 Majorana reference, '
            'declared d4H, circle radius and spectral scale 1; no physical P1 reference choice.',
            'native_A': coefficient, 'ell': ell,
            'base_radial_hessian_finite_difference_checks': base_hessian,
            'log_mass_tail_bound_r_ge_8_for_all_schemes': tail,
            'prescriptions': result}


if __name__ == '__main__':
    print(json.dumps(compare(Path(sys.argv[1]).resolve()), indent=2))
