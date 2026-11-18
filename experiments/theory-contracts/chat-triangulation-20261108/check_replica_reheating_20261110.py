"""R048 same-regulator conical/Einstein test and external-input reheating checks."""
import json
import sys
from pathlib import Path

import mpmath as mp
import sympy as sp

sys.path.insert(0, str(Path(sys.argv[1]) / "verification"))
from tfpt_constants import Mbar, c3, phi0
from v86_nstar_reheating import nstar_chain
from v246_unification_data_crosscheck import AINV, B1, M_Z
from v494_cosmo_running_mu_record import _observables, _phi_of_N, _phi_end


def replica():
    n, area, rs, xi, v, J, c, q = sp.symbols("n A Rsm xi v J C Delta")
    cone = (1/n-n)/12
    assert sp.diff(cone, n).subs(n, 1) == -sp.Rational(1, 6)
    integrated_r = n*rs + 4*sp.pi*(1-n)*area
    entropy = lambda s: sp.expand((n*sp.diff(s, n)-s).subs(n, 1))
    assert entropy(-c*integrated_r) == 4*sp.pi*c*area
    assert entropy(xi*v**2*integrated_r) == -4*sp.pi*xi*v**2*area
    C_scalar = (sp.Rational(1, 6)-xi)*J/(32*sp.pi**2)
    S_scalar = entropy(-C_scalar*integrated_r)
    entropy_canonical = J*area/(48*sp.pi)
    contact = -xi*J*area/(8*sp.pi)
    assert sp.expand(S_scalar - entropy_canonical - contact) == 0
    assert C_scalar.subs(xi, sp.Rational(1, 6)) == S_scalar.subs(xi, sp.Rational(1, 6)) == 0
    assert entropy(-q*integrated_r+q*integrated_r) == 0
    assert entropy(-q*integrated_r+q*rs) == 4*sp.pi*q*area - q*rs
    # The -q*rs term in a literal fixed denominator is a global replica offset.
    # Once subtracted by a metric-independent normalization (or by n copies),
    # it changes no geometric derivative and contributes no local area term.
    return {
        "local_area_per_EH_coefficient": str(4*sp.pi),
        "scalar_coefficient_of_minus_R": str(C_scalar),
        "scalar_ordinary_canonical_area": str(entropy_canonical),
        "scalar_Wald_contact_area": str(contact),
        "scalar_conformal_replica_area": str(S_scalar.subs(xi, sp.Rational(1, 6))),
        "classical_Higgs_vev_area": str(entropy(xi*v**2*integrated_r)),
        "identical_comoving_reference_area": str(entropy(-q*integrated_r+q*integrated_r)),
        "fixed_constant_reference_changes_no_local_area_coefficient": True
    }


def reheating():
    M = c3**(mp.mpf(7)/2)*Mbar
    mu = M/2
    b = [mp.mpf(41)/10, -mp.mpf(19)/6, -mp.mpf(7)]
    assert [float(bi) for bi in b] == list(B1)
    a0 = list(map(lambda x: mp.mpf(str(x)), AINV))
    alphas = [1/(ainv-bi/(2*mp.pi)*mp.log(mu/mp.mpf(str(M_Z))))
              for ainv,bi in zip(a0,b)]
    # Gorbunov-Tokareva (9), Jeong et al. (2.14).
    # The product b*alpha is invariant under changing U(1) normalization.
    assert mp.almosteq(b[0]*alphas[0], (mp.mpf(41)/6)*(mp.mpf(3)/5)*alphas[0])
    width_unit = M**3/(192*mp.pi*Mbar**2)
    widths = [width_unit*mult*(bi*ai)**2/(4*mp.pi**2)
              for mult, bi, ai in zip((1,3,8),b,alphas)]
    total = sum(widths)
    assert widths[2] > mp.mpf("0.9")*total

    def gamma_majorana(x):
        if x >= mp.mpf("0.5"):
            return mp.mpf("0")
        return 2*width_unit*x*x*(1-4*x*x)**(mp.mpf(3)/2)

    assert gamma_majorana(3) == 0
    x = sp.symbols('x', real=True, positive=True)
    shape = x*x*(1-4*x*x)**sp.Rational(3,2)
    assert sp.simplify(sp.diff(shape,x)-2*x*sp.sqrt(1-4*x*x)*(1-10*x*x)) == 0
    x1 = phi0**2/10  # v212 conditional light Majorana rung
    max_one = gamma_majorana(1/mp.sqrt(10))
    # v86 expresses widths in units 4*width_unit, keep all else fixed.
    def history(gamma):
        n, temp, width = nstar_chain(gamma/(4*width_unit), "0.05")
        assert mp.almosteq(gamma,width)
        ns_exact_v, r_exact_v, _, _ = _observables(_phi_of_N(float(n), _phi_end()))
        return {
            "Gamma_GeV": float(gamma),
            "T_GeV_H_equals_Gamma": float(temp),
            "N_star": float(n),
            "n_s_LO": float(1-2/n),
            "r_LO": float(12/n**2),
            "n_s_exact_potential_first_slow_roll": float(ns_exact_v),
            "r_exact_potential_first_slow_roll": float(r_exact_v)
        }

    return {
        "M_scal_GeV": float(M),
        "scale_M_over_2_GeV": float(mu),
        "v246_SM_alphas_GUT_normalized": list(map(float, alphas)),
        "Gamma_B_W_g_GeV": list(map(float, widths)),
        "Gamma_v212_light_candidate_GeV": float(gamma_majorana(x1)),
        "Gamma_v986_heavy_candidate_GeV": 0,
        "Gamma_upper_per_open_Majorana_GeV": float(max_one),
        "v86_as_published": history(16*width_unit),
        "minimal_SM_doublet_correct_width": history(4*width_unit),
        "conformal_gauge_only": history(total),
        "conformal_gauge_plus_v212_N1": history(total + gamma_majorana(x1)),
        "conformal_gauge_plus_two_open_Majoranas_max_production_only": history(total + 2*max_one)
    }


mp.mp.dps = 35
print(json.dumps({"replica": replica(), "reheating": reheating()}, indent=2))
