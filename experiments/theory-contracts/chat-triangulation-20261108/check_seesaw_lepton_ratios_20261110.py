"""Audit the frozen TFPT type-I texture as a charged-lepton mass correction.

The SM+N continuation is explicitly a diagnostic extrapolation, not the
physical Pati-Salam/GUT evolution. The source-to-MSbar boundary is assumed
for comparison. Running traces common to charged families cancel in shape.
"""
import json
import math
import runpy
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm
import yaml

repo = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(repo / "verification"))
import v420_fpole_qed_running as qed
import v466_seed_leptonmass_channel as seed
import v481_seesaw_carrier_ladder as sm
import v1001_nu_pentagon_chain as nu

source_file = repo / "experiments/nu-scalaron-falsification/hypotheses/nu_scalaron_v3.yaml"
with source_file.open() as stream:
    freeze = yaml.safe_load(stream)
scale = freeze["majorana_operator"]["eigenvalues_GeV"]
m2, m3 = float(scale["M2"]), float(scale["M3"])
m_scal = m3 / 3
mbar = sm.MBARF
loop = 1 / (16 * math.pi**2)
v = sm.V_EW
yt3, R3 = sm.run_sm_up(m3)
yt2, _ = sm.run_sm_up(m2)
y3 = yt3
y2 = y3 * freeze["dirac_texture"]["y2_over_y3"]
U = np.array((nu.v9_mixing_matrix() * nu.rotation13(nu.THETA_PENTAGON, nu.PHASE)).tolist(), complex)
assert np.allclose(U.conj().T @ U, np.eye(3), atol=1e-14)
assert math.isclose(abs(U[0,2])**2, freeze["observable_vector"]["mixing"]["sin2_theta13"])
assert math.isclose(abs(U[1,2])**2/(1-abs(U[0,2])**2), freeze["observable_vector"]["mixing"]["sin2_theta23"])
c3 = sm.C3
assert math.isclose(m_scal, c3**3.5 * mbar, rel_tol=1e-8)
P = [np.outer(col, col.conj()) for col in U.T]
coeff = np.array([1., -3., 2.])
source_masses = np.array([16/7, 4/3, 7/6]) * seed.PHI0**np.array([5, 3, 2])
pole = np.array([qed.ME, qed.MMU, qed.MTAU], float)
alpha = float(qed.ALPHA)
msbar = pole * (1 - alpha/math.pi*(1 + 3/4*np.log(float(qed.MZ)**2/pole**2)))


def shape(m):
    x, y = np.log(m[:2] / m[1:])
    return np.array([x, y, x-2*y])


def q(m):
    return float(m.sum()/np.sqrt(m).sum()**2)


def observable(m):
    return dict(log_ratios=shape(m).tolist(), K=float(np.exp(shape(m)[2])), Koide=q(m))


def gauge(t):
    ainv0 = np.array(sm.A_INV_MZ)
    b = np.array([41/10, -19/6, -7])
    return 4*math.pi/(ainv0-b*t/(2*math.pi))


def profiles(mu_hi, limit, yt0, y0, active):
    # Antusch D.33c,f, in the fixed orthogonal projectors, neglecting the
    # tiny charged-lepton and non-top quark matrices in these profiles.
    # I2,I3 integrate each retained squared Yukawa over logarithmic scale.
    def ode(t, state):
        top, yr2, yr3, *_ = state
        g1, g2, g3 = gauge(t)
        ys = np.array([yr2, yr3]) * np.array(active)
        trace = 3*top**2 + ys@ys
        dt = top * (1.5*top**2 + trace - 8*g3 - 2.25*g2 - .85*g1)
        dn = ys * (1.5*ys**2 + trace - 2.25*g2 - .45*g1)
        return np.r_[loop*np.r_[dt,dn], ys**2]
    lo, hi = (math.log(mu/sm.MZ) for mu in (limit, mu_hi))
    state = solve_ivp(ode,(lo,hi),np.r_[yt0,y0,0.,0.],rtol=2e-11,atol=3e-13).y[:,-1]
    return state


target_delta = shape(msbar)-shape(source_masses)
pole_delta = shape(pole)-shape(source_masses)
seed_fit_log = target_delta[0] / 2
tau_fit_log = seed_fit_log - target_delta[1]
tau_fit_y = math.sqrt(tau_fit_log/(1.5 * loop * math.log(mbar/m3)))
y_t_benchmark_mnu = 1e9*(y3*v)**2/(2*m3)*R3
aligned_inversion_mnu = 1e9*(tau_fit_y*v)**2/(2*m3)*R3

# The v2 source calibrates this frozen y2 at its own M2 threshold.
y2_at_m2 = y2
y2_at_m3 = profiles(m3, m2, yt2, [y2_at_m2, 0.], [1,0])[1]


def apply_response(k, m):
    return np.sort(np.linalg.svd(expm(k) @ np.diag(m), compute_uv=False))


def transfer(mu_hi):
    below = profiles(min(mu_hi,m3),m2,yt2,[y2_at_m2,0.],[1,0])
    integrals = below[3:5]
    n3 = 0.
    if mu_hi > m3:
        above = profiles(mu_hi,m3,yt3,[y2_at_m3,y3],[1,1])
        integrals = integrals + above[3:5]
        n3 = y3
    # Antusch D.33d (high to low). Zhang–Zhou 87 at each mu=Mi.
    # Replacing the linear matching matrix by exp shifts only higher orders.
    running = 1.5*loop*(integrals[0]*P[1] + integrals[1]*P[2])
    matching = 5/8*loop*(y2_at_m2**2*P[1] + n3**2*P[2])
    led = np.real(np.diag(running+matching))
    predicted = apply_response(running+matching,source_masses)
    return {"active_heaviest_M3": bool(mu_hi > m3),
            "integrated_nu_squared": integrals.tolist(),
            "leading_delta_log_ratios": [float(led[0]-led[1]),float(led[1]-led[2]),float(coeff@led)],
            "matrix_delta_log_ratios": (shape(predicted)-shape(source_masses)).tolist(),
            "finite_matching_delta_log_K": float(coeff@np.real(np.diag(matching))),
            "predicted_shape_with_neutrino_contribution": observable(predicted),
            "target_residual_after_native_neutrino": (shape(msbar)-shape(predicted)).tolist()}

out = {
  "scope": "native frozen TFPT texture plus explicitly assumed SM+N leading-order continuation; only charged-lepton flavor-split type-I contribution",
  "sources": ["verification/v420_fpole_qed_running.py",
              "verification/v466_seed_leptonmass_channel.py",
              "verification/v481_seesaw_carrier_ladder.py",
              "verification/v1001_nu_pentagon_chain.py",
              str(source_file.relative_to(repo)),
              "https://arxiv.org/abs/hep-ph/0501272 Eq D.33",
              "https://arxiv.org/abs/2107.12133 Eq 87"],
  "Mscal_GeV": m_scal, "heavy_M2_GeV": m2, "heavy_M3_GeV": m3,
  "Mbar_over_M3_log": math.log(mbar/m3), "y3": y3, "y2": y2, "m3_running_factor": R3,
  "v3_m3_eV_reproduced": y_t_benchmark_mnu, "pmns_column_3_abs_squared": np.abs(U[:,2])**2,
  "pmns_column_3_log_K_factor": float(coeff @ np.abs(U[:,2])**2),
  "dominant_pmns_aligned_positive_sign_requires_sin2_theta23_below":
      float((2-abs(U[0,2])**2)/(5*(1-abs(U[0,2])**2))),
  "source": observable(source_masses), "pole": observable(pole), "MSbar_MZ_QED_only": observable(msbar),
  "source_to_MSbar_delta_log_ratios": target_delta.tolist(),
  "source_to_pole_delta_log_ratios": pole_delta.tolist(),
  "source_over_MSbar_ratio_pct": (100*np.expm1(-target_delta[:2])).tolist(),
  "hypothetical_tau_only_decomposition_with_fitted_seed": {
    "seed_delta_log_phi": seed_fit_log, "tau_delta_log": tau_fit_log, "y3_fit": tau_fit_y,
    "m3_low_eV_with_same_kappa_running": aligned_inversion_mnu,
    "neutrino_sum_eV_with_v3_m2": aligned_inversion_mnu+freeze["observable_vector"]["masses"]["m2_eV"],
    "Cabibbo_from_native_phi": seed.OBS5[4]["fwd"](seed.PHI0),
    "Cabibbo_after_hypothetical_shared_seed": seed.OBS5[4]["fwd"](seed.PHI0*math.exp(seed_fit_log)),
    "Cabibbo_inverse_seed_from_existing_v466_observation": seed.OBS5[4]["inv"](seed.OBS5[4]["m"]),
    "native_phi": seed.PHI0
  },
  "shape_running_comparisons": {label:transfer(mu) for label,mu in (
      ("scalaron_scale_v163_comparison",m_scal),
      ("GUT_scale_v481_comparison",sm.M_GUT),
      ("reduced_Planck_scale_attachment_assumption",mbar))},
}
out["SM_no_N_tau_self_running_scale_over_Mbar"] = float(1.5*loop*(math.sqrt(2)*msbar[2]/(1000*v))**2 * math.log(mbar/sm.MZ))

# Independent earlier seam selector, not the frozen v3 texture. Hold the
# tree-level low-mass operator and heavy scales fixed under the CI freedom.
selector = runpy.run_path(str(repo/"experiments/theory-contracts/nu02_transport_selector.py"))
heavy = np.array([scale["M1"],m2,m3])
y0 = np.diag([0,y2,y3]) @ U.conj().T
inverse_mass = np.diag(1/heavy)
reference_kappa = y0.T @ inverse_mass @ y0
fixed_weights = loop*(1.5*np.log(mbar/heavy)+5/8)
ga, gb = [gauge(math.log(mu/sm.MZ)) for mu in (m3,mbar)]
negative_bound = np.array([.45,2.25]) @ np.maximum(ga[:2],gb[:2])
alternatives = {}
for z in (0.,selector["GAP"],-selector["GAP"]):
    Y = (np.sqrt(heavy)[:,None]*selector["R_of_z"](z)/np.sqrt(heavy)[None,:]) @ y0
    kappa_defect = float(np.linalg.norm(Y.T @ inverse_mass @ Y-reference_kappa)/np.linalg.norm(reference_kappa))
    assert kappa_defect < 1e-12
    high = float(np.linalg.svd(Y,compute_uv=False)[0])
    gram_weighted = Y.conj().T @ np.diag(fixed_weights) @ Y
    row = {"relative_tree_kappa_defect":kappa_defect, "largest_Ynu_singular_value":high,
           "fixed_leading_log_plus_matching_delta_log_K":float(coeff@np.diag(gram_weighted).real)}
    if 2.5*high**2 > negative_bound:
        # With charged Yukawas neglected, s=max singular value satisfies
        # ds/dlog(mu) >= loop*s*(2.5*s^2-negative_bound) above M3.
        # This comparison supplies an upper bound for the one-loop pole
        # and for reaching sqrt(4pi) inside the assumed SM+N continuation.
        bottom = 2.5-negative_bound/high**2
        reach = math.log((2.5-negative_bound/(4*math.pi))/bottom)/(2*loop*negative_bound)
        pole_log = math.log(2.5/bottom)/(2*loop*negative_bound)
        assert pole_log < math.log(mbar/m3)
        row.update(one_loop_perturbativity_scale_upper_bound_GeV=m3*math.exp(reach),
                   one_loop_pole_scale_upper_bound_GeV=m3*math.exp(pole_log))
    alternatives[str(z)] = row
out["NU_CI_02_seam_gap_hybrid_not_v3_texture"] = alternatives

assert out["source_to_MSbar_delta_log_ratios"][2] > .04
assert out["shape_running_comparisons"]["reduced_Planck_scale_attachment_assumption"]["matrix_delta_log_ratios"][2] < -.005
assert abs(out["shape_running_comparisons"]["scalaron_scale_v163_comparison"]["matrix_delta_log_ratios"][2]) < 2e-6
print(json.dumps(out, indent=2, default=lambda x: x.tolist()))
