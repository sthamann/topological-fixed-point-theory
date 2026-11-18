"""Exact A3 family readouts and frozen, no-fit lepton pole comparisons."""

import json
import sys
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sp

repo = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(repo / "verification"))
from tfpt_constants import phi0, c3, Mbar
from v17_hexagonal_resolvent import U6_power
from v18_quark_yukawa import TAB
from v117_monodromy_weyl_a3 import M0_EXACT
from v120_address_table import L_LEP
from v252_full_finite_triple import IDX, full_D, quaternion, rep
from v420_fpole_qed_running import ME, MMU, MTAU, MZ

mp.mp.dps = 42
one = sp.eye(3)
half = sp.Rational(1, 2)
t = sp.symbols("t", real=True)
w = (-1+sp.I*sp.sqrt(3))/2
n = sp.Matrix([5, 3, 2])
contrast = sp.Matrix([1, -3, 2])
factor = mp.mpf(53)/54
u_corner = factor * phi0


def spectral_projector(unitary, eigenvalue, order):
    result = sum((sp.conjugate(eigenvalue)**j * unitary**j for j in range(order)),
                 sp.zeros(unitary.rows))
    return sp.simplify(result / order)


def source_data():
    rows = []
    c = []
    u6 = U6_power(1)
    m6 = sp.diag(M0_EXACT, -M0_EXACT)
    assert sp.factor(m6.charpoly().as_expr()) == sp.factor(u6.charpoly().as_expr())
    for label, (sheet, root) in zip(L_LEP, ((1, w), (-1, w), (-1, sp.Integer(1)))):
        winding, residue = divmod(L_LEP[label], 6)
        xi = sheet*root
        assert sp.simplify(sp.expand_complex(sp.exp(sp.I*sp.pi*residue/3))-xi) == 0
        p3 = spectral_projector(M0_EXACT, root, 3)
        p6 = spectral_projector(m6, xi, 6)
        assert p6 == (sp.diag(p3, sp.zeros(3)) if sheet == 1 else sp.diag(sp.zeros(3), p3))
        assert sp.simplify(p3*p3-p3) == sp.zeros(3) and sp.trace(p3) == 1
        squared_green = sp.simplify(1 / ((1-half*xi)*(1-half*sp.conjugate(xi))))
        coeff = 4**winding * squared_green if label != "tau" else sp.Rational(32,9)/(c[0]*c[1])
        c.append(coeff)
        rows.append({"family": label, "L": L_LEP[label], "sheet": sheet, "eigenvalue": xi,
                     "residue": residue, "winding": winding,
                     "squared_green": squared_green, "mass_coefficient": coeff})
    return c, rows


coefficients, addresses = source_data()
assert coefficients == [sp.Rational(16,7), sp.Rational(4,3), sp.Rational(7,6)]
assert sp.Abs(M0_EXACT[1,1]) == half
de, dm = 1+t+t*t, 1-t+t*t
dminus = sp.factor(sp.det(one-t*M0_EXACT))
dplus = sp.factor(sp.det(one+t*M0_EXACT))
assert dminus == sp.factor(1-t**3) and dplus == sp.factor(1+t**3)
assert sp.cancel(dminus/(1-t)) == de and sp.cancel(dplus/(1+t)) == dm
mode = sp.Matrix([4/de, 1/dm, sp.Rational(8,9)*de*dm])
det = sp.Matrix([2/dminus, sp.Rational(3,2)/dplus, sp.Rational(4,3)*dminus])
assert mode.subs(t,half) == det.subs(t,half) == sp.Matrix(coefficients)
slopes = [a.applyfunc(lambda x: sp.factor(sp.diff(x,t)/x).subs(t,half)) for a in (mode,det)]
assert slopes == [sp.Matrix([-sp.Rational(8,7),0,sp.Rational(8,7)]),
                  sp.Matrix([sp.Rational(6,7),-sp.Rational(2,3),-sp.Rational(6,7)])]
assert slopes[1]-slopes[0] == sp.Rational(4,3)*n-sp.Rational(14,3)*sp.ones(3,1)
assert contrast.dot(n) == contrast.dot(sp.ones(3,1)) == 0
assert all(contrast.dot(slope) == sp.Rational(8,7) for slope in slopes)
ratio_jacobians = [sp.Matrix([[-2, s[1]-s[0]], [-1, s[2]-s[1]]]) for s in slopes]
assert all(j.det() == -sp.Rational(8,7) for j in ratio_jacobians)


def k_from_coeff(a):
    return sp.factor(a[0]*a[2]**2/a[1]**3)


assert sp.cancel(k_from_coeff(det)/k_from_coeff(mode) -
                 (1-(2*t-1)**2*(t*t+2)/(3*dm**2))) == 0
assert sp.cancel(sp.prod(mode)) == sp.Rational(32,9)
assert sp.cancel(sp.prod(det)-4/(1+t**3)) == 0

mode_at = sp.lambdify(t, list(mode), "mpmath")
det_at = sp.lambdify(t, list(det), "mpmath")


def source_mass(c, u):
    return tuple(mp.mpf(ci)*u**int(ni) for ci, ni in zip(c,n))


def ratio(masses):
    e, mu, tau = masses
    return {"mu_over_e": mu/e, "tau_over_mu": tau/mu,
            "Q": (e+mu+tau)/(mp.sqrt(e)+mp.sqrt(mu)+mp.sqrt(tau))**2,
            "K": e*tau**2/mu**3}


def old_mobius(delta):
    a = (1+delta)/(1-delta)
    b = (mp.mpf(1)/3+delta)/(mp.mpf(1)/3-delta)
    return (1/(a*b)**2, mp.mpf(1), a*a)


def apply_lengths(lam, amplitudes):
    return tuple(lam**L_LEP[name]*amp for name,amp in zip(L_LEP,amplitudes))


# Show that the current real channels yield the old linear-fractional expression
# while the assigned electron and muon read the different complex pair.
y = sp.symbols("y", positive=True)
assert sp.cancel((1/(y-t))/(1/(y+t))-(y+t)/(y-t)) == 0
assert sp.cancel((y*y+y*t+t*t)/(y*y-y*t+t*t)-((y+t)/(y-t))**2) != 0

pole = ratio((ME,MMU,MTAU))
source = source_mass(coefficients,phi0)
corner = source_mass(coefficients,u_corner)
delta_ph = ((794-7*mp.sqrt(9961))/2187)**(mp.mpf(1)/6)
delta_old = mp.mpf(3)/5+phi0/6
lam = mp.sqrt(phi0*(1-phi0))
exact_amplitudes = tuple(m/lam**L_LEP[name] for name,m in zip(L_LEP,source))
v18_amplitudes = tuple(mp.mpf(TAB[name][1]) for name in L_LEP)
residual = [mp.log((ME,MMU,MTAU)[i]/corner[i]) for i in range(3)]
residual = [x-sum(residual)/3 for x in residual]
results = {
    "frozen_v420_poles": pole,
    "v20_v118_original_source": ratio(source),
    "v18_rounded_table_source": ratio(apply_lengths(lam,v18_amplitudes)),
    "v183_common_seed_lift": ratio(corner),
    "historical_mobius_unfitted": ratio(old_mobius(delta_old)),
    "v19_actual_bare_resolvent_at_delta_ph": ratio(tuple(lam**L_LEP[name]*delta_ph**row["residue"] /
                                                  (1-delta_ph**6) for name,row in zip(L_LEP,addresses))),
    "v183_on_length_amplitude_fixed_Lambda": ratio(apply_lengths(factor*lam,exact_amplitudes)),
    "v183_on_seed_amplitude_fixed_Lambda": ratio(apply_lengths(
        mp.sqrt(u_corner*(1-u_corner)),exact_amplitudes)),
}

for tag, delta in (("closed_other_branch",delta_ph),("historical_phase",delta_old)):
    for extension, fun in (("v20_mode",mode_at),("v118_det",det_at)):
        results[f"{extension}_replace_delta_with_{tag}_at_corner"] = ratio(source_mass(fun(delta),u_corner))

# Only to prove the single-coordinate failure, locate the unique t from mu/e
# on the real 0<t<1 branch. It is never used as a fit-free prediction.
de_dm_derivative = sp.factor(sp.diff(de/dm,t))
assert sp.cancel(de_dm_derivative-2*(1-t*t)/(dm*dm)) == 0
required_t = mp.findroot(lambda a: ratio(source_mass(mode_at(a),u_corner))["mu_over_e"]
                         - pole["mu_over_e"], (mp.mpf("0.40"),mp.mpf("0.50")))
mode_holdout = ratio(source_mass(mode_at(required_t),u_corner))
assert abs(mode_holdout["tau_over_mu"]-pole["tau_over_mu"]) > mp.mpf("0.5")
# Determinant lift: both mass ratios strictly decrease on 0<t<1.
det_ratio_grad = (sp.factor(sp.diff(det[1]/det[0],t)), sp.factor(sp.diff(det[2]/det[1],t)))
assert det_ratio_grad == (-9*t*t/(2*(t+1)**2*(t*t-t+1)**2), -sp.Rational(16,3)*t**5)

# Original finite v252 charged Yukawa consumer. Spectator choices do not
# enter this charged-lepton block, including the neutral Majorana coupling.
m = np.diag([float(mp.pi*x) for x in source])
zero = np.zeros((3,3),complex)
D = full_D(zero,zero,m,zero,np.eye(3))
left = [IDX[("e","L",0,g)] for g in range(3)]
right = [IDX[("e","R",0,g)] for g in range(3)]
assert np.array_equal(D[np.ix_(left,right)],m)
b = rep(1,np.zeros((2,2)),zero)
for h in (0.3,1.7):
    a = rep(0,quaternion(h,0),zero)
    fluctuation = a @ (D@b-b@D)
    assert np.array_equal(fluctuation[np.ix_(left,right)],h*m)
    assert np.allclose(np.linalg.svd(h*m,compute_uv=False) /
                       np.linalg.svd(m,compute_uv=False), h)

# Family-diagonal projection of the existing Hardy tensor in the same
# normalized Fuchs chart. Factor out kappa and P_gamma,+ as in the source.
cs = sp.Matrix([[0,0,1],[-1,0,0],[0,-1,0]])
b02 = sp.zeros(3)
b02[0,2] = 1
hardy = [sp.simplify(sp.trace(spectral_projector(cs,root,3)*b02)) for root in (w,w,sp.Integer(1))]
assert hardy == [w/3,w/3,sp.Rational(1,3)]
hardy_gram_y = mp.mpf(2)/3
raw_hardy_weights = [1/mp.sqrt(hardy_gram_y**2+delta_ph**2-2*hardy_gram_y*delta_ph*r)
                     for r in (-mp.mpf(1)/2,mp.mpf(1)/2,-mp.mpf(1))]
# Pure charged Yukawa part of v163 PyR@TE one-loop beta, frozen at source
# Yukawas for a leading-log scale diagnostic, not a full SM pole match.
scale_ratio = (c3**(mp.mpf(7)/2)*Mbar)/(MZ/1000)
beta_shift = tuple(-3*(mp.pi*mi)**2*mp.log(scale_ratio)/(32*mp.pi**2) for mi in source)
beta_ratio_shifts = [beta_shift[1]-beta_shift[0],beta_shift[2]-beta_shift[1]]
results["v163_charged_yukawa_cubic_leading_log_only"] = ratio(
    tuple(mi*mp.exp(shift) for mi,shift in zip(source,beta_shift)))

# Independent angular diagnostic from the older literature and the numerically
# identical seam-twist exponent. These are conditional substitutions, not an
# identification of an exponent with the mass azimuth.
theta = mp.mpf(2)/9
def masses_from_angle(q, azimuth):
    amplitude = mp.sqrt(2*(3*q-1))
    roots = [1+amplitude*mp.cos(azimuth+2*mp.pi*k/3) for k in (1,2,0)]
    assert all(x > 0 for x in roots)
    return tuple(x*x for x in roots)


def v82_q_after(steps):
    q = 3*ratio(source)["Q"]
    rho = ((q-2)/(5-q))*(mp.mpf(2)/3)**(6*steps)
    return (2+5*rho)/(3*(1+rho))


angle_candidates = {
    "exact_Koide": mp.mpf(2)/3,
    "v183_corner_Q": ratio(corner)["Q"],
    "v82_three_steps": v82_q_after(3),
    "v82_four_steps": v82_q_after(4),
}
angle_rows = {}
for name,q in angle_candidates.items():
    vals = ratio(masses_from_angle(q,theta))
    assert abs(vals["Q"]-q) < mp.mpf("1e-38")
    vals["mu_over_e_residual_ppm"] = 1e6*(vals["mu_over_e"]/pole["mu_over_e"]-1)
    vals["mu_over_e_residual_over_CODATA_2022_uncertainty"] = (
        (vals["mu_over_e"]-mp.mpf("206.7682827"))/mp.mpf("0.0000046"))
    angle_rows[name] = vals
assert all(abs(row["mu_over_e_residual_over_CODATA_2022_uncertainty"])>30
           for row in angle_rows.values())
pe, pm, pt = map(mp.sqrt, (ME,MMU,MTAU))
pole_azimuth = mp.atan2((pm-pe)/mp.sqrt(2),(2*pt-pe-pm)/mp.sqrt(6))

output = {
    "original_A3_responses": addresses,
    "readouts_no_fit": results,
    "delta_ph_other_branch": delta_ph,
    "historical_delta": delta_old,
    "v20_log_derivatives": list(slopes[0]),
    "v118_log_derivatives": list(slopes[1]),
    "common_log_K_derivative": contrast.dot(slopes[0]),
    "v20_t_inferred_from_one_ratio_diagnostic_only": required_t,
    "v20_other_ratio_holdout": mode_holdout,
    "needed_fixed_corner_logs_scale_removed": residual,
    "needed_log_ratio_shifts": [residual[1]-residual[0],residual[2]-residual[1]],
    "actual_v252_neutral_higgs_family_log_derivative": [1,1,1],
    "normalized_hardy_diagonal_coeffs_over_kappa": hardy,
    "bare_hardy_K_weights_at_original_y_two_thirds_zero_regulator": raw_hardy_weights,
    "v163_one_loop_cubic_log_ratio_shifts": beta_ratio_shifts,
    "pole_azimuth_rad": pole_azimuth,
    "conditional_theta_two_ninths": angle_rows,
}
print(json.dumps(output,indent=2,default=lambda x: mp.nstr(x,20) if isinstance(x,mp.mpf) else str(x)))
