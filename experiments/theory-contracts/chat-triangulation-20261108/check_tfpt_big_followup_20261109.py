"""Audit the existing TFPT lepton matching candidates on both mass ratios."""
import json
import mpmath as mp
import sympy as sp
from tfpt_constants import phi0
from v420_fpole_qed_running import ALPHA, ME, MMU, MTAU, MZ, koide

mp.mp.dps = 45
phi = mp.mpf(phi0)
source = [mp.mpf(16)/7*phi**5, mp.mpf(4)/3*phi**3, mp.mpf(7)/6*phi**2]
pole = [ME, MMU, MTAU]
a = ALPHA/mp.pi


def log_ratios(masses):
    e, mu, tau = masses
    return mp.log(e/mu), mp.log(mu/tau)


def checksum(masses):
    x, y = log_ratios(masses)
    return x-2*y


def to_msbar(m, scale):
    return m * (1 - a*(1+mp.mpf(3)/4*mp.log((scale/m)**2)))


def koide_shape(x, y):
    return koide(mp.exp(x+y), mp.exp(y), 1)


def shape_jacobian(masses):
    x, y = log_ratios(masses)
    dx = mp.diff(lambda t: koide_shape(t,y), x)
    dy = mp.diff(lambda t: koide_shape(x,t), y)
    return -2*dx-dy


source_log = checksum(source)
pole_log = checksum(pole)
needed = pole_log-source_log
one_loop = -mp.mpf(3)/2*a*source_log
literal_msbar = [to_msbar(m,MZ) for m in pole]
exponents = [b/a for a,b in zip(log_ratios(source),log_ratios(pole))]
assert needed > 0 > one_loop
assert koide(*source) < koide(*pole) < koide(*literal_msbar)
assert abs(exponents[0]-exponents[1]) > .006
assert all(abs(shape_jacobian(ms)) > .1 for ms in [source,pole])
# Direction remains the same under the one-loop power law of v371.
alpha_interval = [1/mp.mpf("137.036"), 1/mp.mpf("128")]
direction = [-mp.mpf(3)/2*(alpha/mp.pi)*source_log for alpha in alpha_interval]
assert all(z < 0 for z in direction)
# Linearized scalar curvature on a flat Euclidean chart.
k = sp.Matrix(sp.symbols("k0:4", real=True))
k2 = (k.T*k)[0]
transverse = sp.eye(4) - k*k.T/k2
curvature = k2*transverse
tt_part = transverse*curvature*transverse - transverse*sp.trace(transverse*curvature)/3
assert tt_part.applyfunc(sp.simplify) == sp.zeros(4)

print(json.dumps({k: float(v) for k,v in {
    "K_source":mp.exp(source_log), "K_pole_v420":mp.exp(pole_log),
    "required_delta_log_K":needed, "required_K_relative_change":mp.expm1(needed),
    "first_order_QED_delta_log_K":one_loop,
    "literal_v420_pole_minus_MSbar_delta_log_K":pole_log-checksum(literal_msbar),
    "Q_source":koide(*source), "Q_pole":koide(*pole),
    "Q_literal_v420_MSbar":koide(*literal_msbar),
    "exponent_from_e_mu":exponents[0], "exponent_from_mu_tau":exponents[1],
    "jacobian_at_source":shape_jacobian(source), "jacobian_at_pole":shape_jacobian(pole),
}.items()} | {"linear_scalar_curvature_TT_exactly_zero": True},indent=2))
