"""Reproduce the original scalar corner transfer and test both mass shapes."""
import json
import mpmath as mp
import sympy as sp

import v183_koide_f_corner_transfer as corner
from v20_lepton_c_derivation import amp2, mu4, g_car, N_fam
from v420_fpole_qed_running import ME, MMU, MTAU

mp.mp.dps = 70

a, one = corner.A, corner.ONE
anchor_numerator = (a.T * corner.F * one)[0]
cubic_denominator = 2 * (one.T * corner.R * a)[0]
factor = anchor_numerator / cubic_denominator
assert (anchor_numerator, cubic_denominator, factor) == (53, 54, sp.Rational(53, 54))
assert corner._detB(corner.F) == 52
assert [(a.T * m * one)[0] for m in (corner.R,corner.Q,corner.K,corner.L)] == [32,21,35,56]

c_e, c_mu = mu4*amp2(2), amp2(5)
c_tau = sp.Rational(2**g_car, N_fam**2)/(c_e*c_mu)
powers = tuple(corner.K.row(2))
assert (c_e,c_mu,c_tau) == (sp.Rational(16,7),sp.Rational(4,3),sp.Rational(7,6))
assert powers == (5,3,2)
u = sp.symbols("u", positive=True)
mass_scale = sp.symbols("s", positive=True)
source_masses = [mass_scale * c * u**n for c,n in zip((c_e,c_mu,c_tau),powers)]
corner_masses = [m.subs(u,factor*u) for m in source_masses]


def K(m):
    e, mu, tau = m
    return e*tau**2/mu**3


assert sp.cancel(K(source_masses)) == sp.cancel(K(corner_masses)) == sp.Rational(21,16)
w = sp.Matrix([1,-3,2])
assert w.dot(one) == w.dot(sp.Matrix(powers)) == 0
Q_symbolic = sum(corner_masses)/sum(sp.sqrt(m) for m in corner_masses)**2
assert mass_scale not in sp.cancel(Q_symbolic).free_symbols

c3 = 1/(8*mp.pi)
u0 = 1/(6*mp.pi) + 48*c3**4
f = mp.mpf(53)/54


def masses(seed):
    return [mp.mpf(str(c))*seed**int(n) for c,n in zip((c_e,c_mu,c_tau),powers)]


def describe(m):
    e,mu,tau = m
    return dict(Q=float(sum(m)/sum(mp.sqrt(x) for x in m)**2), K=float(K(m)),
                e_mu=float(e/mu), mu_tau=float(mu/tau))


src, out, poles = masses(u0), masses(f*u0), [ME, MMU, MTAU]
input_points = [describe(m) for m in (src, out, poles)]
d_source,d_corner,d_pole = input_points
assert abs(d_corner["Q"] - 2/3) < 6e-7
assert abs(d_corner["K"]-d_pole["K"]) > .05
tau_unc = mp.mpf(".09")
lo, hi = describe([ME,MMU,MTAU-tau_unc]), describe([ME,MMU,MTAU+tau_unc])
assert lo["Q"] < d_corner["Q"] < hi["Q"]
assert lo["K"] > d_corner["K"]
u_from_poles = [mp.sqrt(7*ME/(12*MMU)), 7*MMU/(8*MTAU)]
assert abs(u_from_poles[0]-u_from_poles[1]) > mp.mpf(".001")
print(json.dumps(dict(operator={"a_F_one":int(anchor_numerator),
                          "two_one_R_a":int(cubic_denominator),"factor":str(factor),
                          "lepton_powers":list(map(int,powers))},
                      source=d_source, candidate=d_corner, archived_pole=d_pole,
                      pole_tau_band={k:[lo[k],hi[k]] for k in ("Q","K")},
                      candidate_relative_residuals={k:d_corner[k]/d_pole[k]-1
                                                    for k in ("e_mu","mu_tau","K")},
                      required_log_K_change=float(mp.log(K(poles)/K(out))),
                      seed0=float(u0), seed_candidate=float(f*u0),
                      implied_seeds=list(map(float,u_from_poles))),indent=2))
