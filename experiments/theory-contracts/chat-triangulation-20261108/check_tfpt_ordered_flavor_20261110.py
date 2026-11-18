"""Test the source-address clock and sheet action against the frozen seesaw trajectory."""
import contextlib
import io
import json
import math
from pathlib import Path
import runpy
import sys

import numpy as np

repo = Path(sys.argv[1]).resolve()
base = repo / "experiments/theory-contracts/chat-triangulation-20261108"
old_argv = sys.argv
try:
    driver = base / "check_seesaw_two_loop_control_20261110.py"
    sys.argv = [str(driver), str(repo), str(base/"sm_type_i_betas_generated_20261110.py")]
    with contextlib.redirect_stdout(io.StringIO()):
        scope = runpy.run_path(str(driver))
finally:
    sys.argv = old_argv

from v120_address_table import L_LEP

transfer = scope["transfer"]
env = transfer.__globals__
w, sm = scope["w"], scope["sm"]
source = w["source_masses"]
D = np.diag(source/source[-1])
identity = np.eye(3)
diag = lambda y: np.diag(y)
norm = lambda m: float(np.linalg.norm(m,2))
singular = lambda m: np.sort(np.linalg.svd(m,compute_uv=False))
output = []
evolve = env["evolve"]
original_mass_vector = w["source_masses"]


def capture(a,b,y,loops):
    result = evolve(a,b,y,loops)
    if b == sm.MZ:
        output.append(env["data"](result)[2][2])
    return result


env["evolve"] = capture
try:
    transfer(w["mbar"],2,True)
    final = output[-1]
    B = final @ np.linalg.inv(D)
    B /= abs(np.linalg.det(B))**(1/3)

    def replay_sheet(r):
        # Same dimensionless initial tau norm, modify only the charged
        # source shape before running the actual matrix beta functions.
        w["source_masses"] = original_mass_vector * np.array([1.,r,r])
        transfer(w["mbar"],2,True)
        return output[-1]
    target = w["shape"](w["msbar"])-w["shape"](original_mass_vector)
    omega = np.exp(2j*np.pi/3)
    xi = np.exp(1j*np.pi*np.array(list(L_LEP.values()))/3)
    assert np.max(np.abs(xi-np.array([omega,-omega,-1]))) < 1e-13
    # Unitary signed clock and unsigned A3 have fixed eigenvectors on the
    # three addresses. They need not commute with B to be mass-invisible.
    clock_tests = {}
    for label, T in [("signed_v120",diag(xi)),("unsigned_A3",diag([omega,omega,1]))]:
        forward, reverse = singular(T @ B @ D), singular(B @ T @ D)
        assert np.max(np.abs(forward-reverse)) < 1e-13
        clock_tests[label] = dict(commutator_norm=norm(T@B-B@T),
                                 max_order_delta_log_ratios=float(np.max(np.abs(w["shape"](forward)-w["shape"](reverse)))))

    # The algebraic signed Fuchs factor is a declared mass-interpretation
    # control, not a projected original self-energy or a new time choice.
    delta = ((794-7*math.sqrt(9961))/2187)**(1/6)
    a,b,c = 6*math.log(1.5),6*math.log(3),12*math.log(2)
    rows = []
    for y in (1.,1/3,2/3):
        R = y*delta/(y-delta)**2
        bn,cn = [a*math.log1p(n*R)/math.log1p(R) for n in (3,4)]
        for t in (.125,1.,4.):
            Fplus = math.expm1(-b*t)/math.expm1(-bn*t)
            Fminus = math.expm1(-(c-a)*t)/math.expm1(-(cn-a)*t)
            r = Fminus/Fplus
            assert r > 1
            F = diag([1.,r,r])
            before = w["shape"](singular(B @ F @ D))
            after = w["shape"](singular(F @ B @ D))
            joint_replay = w["shape"](singular(replay_sheet(r))) if t == 1. else None
            rows.append(dict(cusp=y,time=t,sheet_ratio=r,commutator_norm=norm(F@B-B@F),
                             frozen_before_shift=(before-w["shape"](singular(D))).tolist(),
                             frozen_after_shift=(after-w["shape"](singular(D))).tolist(),
                             before_minus_after=(before-after).tolist(),
                             nonlinear_two_loop_before_shift=None if joint_replay is None else
                             (joint_replay-w["shape"](singular(D))).tolist(),
                             nonlinear_before_minus_frozen_after=None if joint_replay is None else
                             (joint_replay-after).tolist()))
finally:
    env["evolve"] = evolve
    w["source_masses"] = original_mass_vector

eps = norm(B-identity)
rho = D[0,0]/D[1,1]
# If v is the lightest right singular vector of B F D, norm of
# the mu/tau coordinates is <= c_before/r by min-max and singular norm.
c_before = ((1+eps)/(1-eps))*rho
# For F B D use the lightest left singular vector u. Its mu/tau
# coordinates are <= c_after/r by applying D B^* F, bounding the
# complementary compression B^* by 1-eps, and testing B^-1 e1.
c_after = eps/(1-eps)+rho/(1-2*eps)
assert eps < .1 and max(c_before,c_after) < .1
upper = -1 + 4*max(c_before,c_after)**2
assert upper < -.95
base_shape = w["shape"](singular(B @ D))
assert base_shape[2] < w["shape"](singular(D))[2]
max_order = np.max(np.abs(np.array([row["before_minus_after"] for row in rows])),axis=0)
assert max_order[2] < .002

# Exact source quarter action. A positive threshold on the three retained
# lines can have different completions on the three complementary lines.
# Their two-ended readings differ by L^*L, already with one quarter link.
arc_path = base / "check_tfpt_mobius_full_sewing_20261109.py"
try:
    sys.argv = [str(arc_path),str(repo)]
    with contextlib.redirect_stdout(io.StringIO()):
        arc = runpy.run_path(str(arc_path))
finally:
    sys.argv = old_argv
sp = arc["sp"]
p6, q6, u6 = (arc[k] for k in ("P3","Q3","U6"))
L = sp.simplify(q6*u6*p6)
D6 = sp.simplify(p6*u6.H*q6*u6*p6)
assert D6 == sp.simplify(L.H*L) and D6.rank()==2
s = sp.symbols("s", positive=True)
assert sp.simplify(p6*u6.H*((p6+s*q6)-p6)*u6*p6-s*D6)==sp.zeros(6)
assert sp.simplify(p6*((p6+s*q6)-p6)*p6)==sp.zeros(6)

# Same original v3 mass/yukawa conventions as in the two-loop report.
Ynu = diag([0.,w["y2"],w["y3"]])@w["U"].conj().T
heavy = np.array([w["scale"]["M1"],w["m2"],w["m3"]])
gram = Ynu.conj().T @ Ynu
takagi = Ynu.T @ diag(1/heavy) @ Ynu
gram_diagonal = np.diag(gram).real
assert np.max(np.abs(takagi-takagi.T)) < 1e-27
# Pure right rotations of the charged Yukawa are mass-invisible.
right_clock = diag(xi)
assert np.max(np.abs(singular(B@D@right_clock)-singular(B@D))) < 1e-13
# A heavy right-neutrino relabelling carries the Majorana form with it.
# For comparison also keep M fixed, which is a physical misalignment,
# not merely the same field in a new basis.
heavy_cycle = np.array([[0.,1.,0.],[0.,0.,1.],[1.,0.,0.]])
Y_rot = heavy_cycle @ Ynu
M_rot = heavy_cycle @ diag(heavy) @ heavy_cycle.T
weight = diag(w["loop"]*(1.5*np.log(w["mbar"]/heavy)+5/8))
weighted_gram = Ynu.conj().T @ weight @ Ynu
assert norm(Y_rot.conj().T@Y_rot-gram) < 1e-15
assert norm(Y_rot.T@np.linalg.solve(M_rot,Y_rot)-takagi) < 1e-27
assert norm(Y_rot.conj().T@(heavy_cycle@weight@heavy_cycle.T)@Y_rot-weighted_gram) < 1e-15
assert norm(Y_rot.T@diag(1/heavy)@Y_rot-takagi) > norm(takagi)
shape_coefficients = [w["msbar"][0]/w["msbar"][1]/(w["seed"].PHI0**2),
                      w["msbar"][1]/w["msbar"][2]/w["seed"].PHI0]

print(json.dumps({
    "scope": "conditional three-family direct-sheet mass readout and actual frozen two-loop SM+N trajectory; separate exact six-line address algebra; no physical P_adm or PS/GUT mapping",
    "coordinate_order": ["log(e/mu)","log(mu/tau)","log(K)"],
    "required_delta": target.tolist(),
    "two_loop_no_Fuchs_delta": (base_shape-w["shape"](singular(D))).tolist(),
    "clock_phase_checks": clock_tests,
    "frozen_two_loop_transfer_spectral_norm_from_identity": eps,
    "d_e_over_d_mu": rho,
    "uniform_positive_sheet_logK_logr_derivative_upper_bound": upper,
    "order_difference_max_abs_nine_checks": max_order.tolist(),
    "mu_tau_filter": {
        "v3_sin_squared_theta23": float(abs(w["U"][1,2])**2/(1-abs(w["U"][0,2])**2)),
        "v3_Ynu_dagger_Ynu_diagonal": gram_diagonal.tolist(),
        "mu_minus_tau_diagonal": float(gram_diagonal[1]-gram_diagonal[2]),
        "actual_source_coefficients_e_mu_then_mu_tau": [12/7,8/7],
        "required_QED_matched_MZ_coefficients_e_mu_then_mu_tau": shape_coefficients,
        "required_delta_log_mu_tau": float(target[1]),
        "two_loop_delta_log_mu_tau": float(base_shape[1]-w["shape"](singular(D))[1]),
        "light_matrix_relative_change_with_heavy_rotation_but_fixed_M": norm(Y_rot.T@diag(1/heavy)@Y_rot-takagi)/norm(takagi),
        "heavy_basis_change_with_cotransformed_M_preserves_both_grams_and_Takagi": True,
    },
    "two_loop_ordered_sheet_controls": rows,
    "quarter_six_line_ambiguity": {
        "complement_probe_difference_formula": "s P U^* (I-P) U P",
        "rank": D6.rank(),
        "nonzero_eigenvalues": [str(value) for value in D6.eigenvals() if value!=0],
        "diagonal_on_actual_addresses": [str(value) for value in arc["trace_per_lepton"]],
    },
},indent=2,default=int))
