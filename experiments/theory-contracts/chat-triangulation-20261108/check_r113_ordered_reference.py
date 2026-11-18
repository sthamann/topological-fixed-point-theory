"""Checks on the original even R032 charged history and its two reference lifts."""
import importlib.util
import json
import sys

import numpy as np
from scipy.integrate import quad
from scipy.linalg import expm
from scipy.special import expit

spec = importlib.util.spec_from_file_location("r032", sys.argv[1])
source = importlib.util.module_from_spec(spec)
spec.loader.exec_module(source)
t = float(source.time)

# Positive pairs in the actual exterior coordinate chart, including its dual.
charged = np.where(np.abs(source.mass_matrix).sum(axis=0)==0)[0]
Y = source.higgs_matrix[np.ix_(charged,charged)]
G = source.grading[np.ix_(charged,charged)]
assert len(charged)==90
# Exclude the massless left-neutrino coordinates, which have no charged coupling.
connected = np.where(np.abs(Y).sum(axis=0)>0)[0]
Y = Y[np.ix_(connected,connected)]
G = G[np.ix_(connected,connected)]
assert len(connected)==84
assert np.linalg.norm(G@Y+Y@G)<1e-12
spectrum = np.linalg.eigvalsh(Y)
y = spectrum[spectrum>1e-9]
assert len(y)==42
y2 = np.sum(y*y)
assert abs(2*y2-np.linalg.norm(Y)**2)<1e-11


def sewn_vs_polar():
    momentum = 0.25
    amplitude = 0.5
    clock = 0.8
    h0 = momentum*G
    mark = -1j*G@Y
    h_plus = h0+amplitude*mark
    h_minus = h0-amplitude*mark
    k = expm(-clock*h_minus/4)@expm(-clock*h_plus/4)
    R = k.conj().T@k
    spectrum = np.linalg.eigvalsh(R)
    direct_log_eigs = np.sort(np.log(spectrum))
    e = np.hypot(momentum,amplitude*y)
    g = 2*np.arcsinh(momentum/e*np.sinh(clock*e/2))
    assert max(abs(direct_log_eigs-np.sort(np.r_[g,-g])))<2e-13
    raw = np.log1p(spectrum).sum()
    sea_sewn = -clock*e.sum()
    sea_polar = -g.sum()
    assert abs(raw+sea_sewn-(2*np.logaddexp(0,-g).sum()+sea_sewn-sea_polar))<1e-12
    assert sea_sewn-sea_polar < -0.01
    ref_s = expm(-clock*h0)
    # On this reflection-even history the original T coefficients are identical.
    t_ref_plus = h0+amplitude*mark
    t_ref_minus = h0-amplitude*mark
    k_t = expm(-clock*t_ref_minus/4)@expm(-clock*t_ref_plus/4)
    assert np.linalg.norm(k_t.conj().T@k_t-R)<1e-13
    s_log_ref = np.linalg.slogdet(np.eye(len(G))+ref_s)[1]-clock*momentum*len(y)
    log_sewn_S = raw+sea_sewn-s_log_ref
    return dict(p=momentum,time=clock,higgs_amplitude=amplitude,
                max_full_operator_log_spectrum_error=float(max(abs(direct_log_eigs-np.sort(np.r_[g,-g])))),
                log_sewn_minus_polar=float(sea_sewn-sea_polar),
                log_charged_sewn_quotient_S=float(log_sewn_S),log_charged_sewn_quotient_T=0.0)


def infrared():
    def polar_jet(p):
        z = t*p/2
        if np.isscalar(p) and p==0:
            return 0.0
        z_minus_tanh = np.where(z<0.01,z**3/3-2*z**5/15+17*z**7/315,z-np.tanh(z))
        return -2*expit(-t*p)*y2*z_minus_tanh/p**2

    def scalar_difference(p):
        return -y2*np.tanh(t*p/2)/p**2

    rows = []
    tail, _ = quad(lambda q: 2*scalar_difference(q),64,np.inf,epsabs=1e-11)
    for r in (32,64,128,256,512,1024):
        p = (2*np.arange(128*r)+1)/(4*r)
        polar = polar_jet(p).sum()/r
        sewn = polar+scalar_difference(p).sum()/r+tail
        rows.append([r,float(sewn),float(polar)])
    limit, err = quad(lambda q: 2*polar_jet(q),0,np.inf,epsabs=1e-10)
    measured = (rows[-1][1]-rows[-2][1])/np.log(2)
    expected = -t*y2
    assert abs(measured-expected)<0.002
    polar_extrapolated = (4*rows[-1][2]-rows[-2][2])/3
    assert abs(polar_extrapolated-limit)<1e-7
    return dict(squared_charged_pair_weight=float(y2),expected_sewn_slope=float(expected),
                measured_sewn_slope=float(measured),polar_density_limit=float(limit),
                extrapolated_polar_density=float(polar_extrapolated),
                polar_integral_error=err,rows_r_sewn_polar=rows)


if __name__=="__main__":
    print(json.dumps(dict(finite=sewn_vs_polar(),infrared=infrared()),indent=2))
