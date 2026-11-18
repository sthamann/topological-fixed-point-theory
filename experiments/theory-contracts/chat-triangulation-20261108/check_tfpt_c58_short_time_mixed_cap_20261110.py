"""Short-time mixed link jet on the actual old C58 Slater cap."""

import json
from pathlib import Path
import runpy
import sys

import numpy as np
from scipy.linalg import expm
from scipy.special import zeta


repo = Path(sys.argv[1]).resolve()
base = repo / "experiments/theory-contracts/chat-triangulation-20261108"
old = runpy.run_path(str(base / "check_tfpt_conditioned_r113_weights_20261110.py"))
native = old["native"]
gamma, swap = old["gamma"], old["swap"]
h, mr = old["clifford_higgs"], old["clifford_majorana"]
vl, vr = mr + .35*h, mr + .55*h
z2, z3 = np.zeros((2,2)), np.zeros((3,3))


def lift(p):
    return -1j*(p+swap@p.conj()@swap)


t3 = lift(native.rep(0, .5j*old["sigma_z"], z3))
y = lift(native.rep(3j, z2, -1j*np.eye(3)))/6
assert np.max(abs((t3+y)@vl-vl@(t3+y))) < 1e-13
assert np.max(abs(gamma@vl+vl@gamma)) < 1e-13
assert np.max(abs(t3@mr-mr@t3)) < 1e-13


def spectral_pair(sign):
    p = np.flatnonzero(sign*np.diag(gamma) < 0)
    q = np.flatnonzero(sign*np.diag(gamma) > 0)
    assert abs(t3.diagonal()[p].sum()) < 1e-12 and abs(y.diagonal()[p].sum()) < 1e-12
    br, bl = vr[np.ix_(p,q)], vl[np.ix_(p,q)]
    ur, sr, vhr = np.linalg.svd(br)
    ul, sl, vhl = np.linalg.svd(bl)
    d_t = t3.diagonal()[q,None]-t3.diagonal()[None,p]
    d_y = y.diagonal()[q,None]-y.diagonal()[None,p]
    coefficient = np.trace(br @ (-d_t*d_y*bl.conj().T)).real
    predicted = .35*.55*np.trace(h[np.ix_(p,q)]@h[np.ix_(q,p)]).real/4
    assert abs(coefficient-predicted) < 1e-12
    return (ur,sr,vhr,ul,sl,vhl,d_t,d_y,float(coefficient))


spectral = {sign: spectral_pair(sign) for sign in (-1, 1)}
coefficient = spectral[1][-1]
assert abs(coefficient-spectral[-1][-1]) < 1e-12


def contractions(p, ell, split, data):
    ur,sr,vhr,ul,sl,vhl,d_t,d_y,_ = data
    # Negative free subspace is the actual old Slater P0. After factoring
    # both diagonal caps, only the two off-diagonal return legs remain.
    def ratios(s, t):
        e = np.hypot(p[:,None], s)
        a = np.tanh(t*e)
        return s*a/(e+p[:,None]*a)
    right = -(ur*ratios(sr,ell*(1-split))[:,None,:])@vhr
    left = -(vhl.conj().T*ratios(sl,ell*split)[:,None,:])@ul.conj().T
    ident = np.eye(len(sr))
    mat = ident + right@left
    first = 1j*right@(d_t*left)
    second = 1j*right@(d_y*left)
    mixed = -right@(d_t*d_y*left)
    li,lj,lij = (np.linalg.solve(mat,x) for x in (first,second,mixed))
    jet = np.trace(lij-li@lj, axis1=-2,axis2=-1)
    # Distinguish the real log-modulus response from the open phase.
    # The central even combination removes the fundamental character.
    def angular(theta):
        return np.linalg.slogdet(ident+right@(np.exp(1j*theta*d_t)*left))[1]
    even = (angular(np.pi)+angular(3*np.pi)-angular(0)-angular(2*np.pi))/2
    return float(jet.real.sum()), float(even.sum()), float(jet.imag.sum())


def cap_read(ell, split, cutoff, include_tail):
    values = np.zeros(3)
    for sign, data in spectral.items():
        shift = .25 if sign==1 else .75
        modes = np.arange(cutoff)+shift
        for start in range(0,cutoff,96):
            values += contractions(modes[start:start+96],ell,split,data)
        if include_tail:
            # Leading large-p term. The cutoff comparisons verify the
            # actual size of neglected O(p^-4) and exponential terms.
            values[0] += coefficient/4*zeta(2,cutoff+shift)
    return {"mixed_log_jet":float(values[0]),"mixed_log_jet_over_ell":float(values[0]/ell),
            "even_combination":float(values[1]),"even_combination_over_ell":float(values[1]/ell),
            "phase_jet_before_tail":float(values[2]),"cutoff_each_sign":cutoff}


def boundary_coefficient(split):
    return -coefficient*(split*np.log(split)+(1-split)*np.log(1-split))


def dense_controls():
    controls = {}
    step, ell = .005, .125
    for momentum in (.25, -.75):
        sign = 1 if momentum > 0 else -1
        occupied = np.flatnonzero(momentum*np.diag(gamma) < 0)
        left = expm(-ell/2*(momentum*gamma+vl))[:,occupied]
        right = expm(-ell/2*(momentum*gamma+vr))[occupied,:]

        def log_modulus(theta, phi):
            link = np.exp(1j*(theta*np.diag(t3)+phi*np.diag(y)))
            return np.linalg.slogdet((right*link[None,:])@left)[1]

        finite_difference = sum(
            a*b*log_modulus(a*step,b*step) for a in (-1,1) for b in (-1,1)
        )/(4*step**2)
        factorized = contractions(np.array([abs(momentum)]),ell,.5,spectral[sign])[0]
        error = abs(finite_difference-factorized)
        assert error < 1e-7
        controls[str(momentum)] = {"full_matrix_finite_difference":float(finite_difference),
                                   "factorized":factorized,"absolute_error":float(error)}
    return controls


results = {}
for split in (.5, 1/3):
    prediction = boundary_coefficient(split)
    rows = {str(ell):cap_read(ell,split,max(128,int(np.ceil(16/ell))),True)
            for ell in (.125,.025,.01,.004)}
    check = cap_read(.004,split,8000,True)
    assert abs(check["mixed_log_jet_over_ell"]-rows["0.004"]["mixed_log_jet_over_ell"]) < 1e-6
    assert abs(check["mixed_log_jet_over_ell"]-prediction) < .04
    assert abs(rows["0.004"]["even_combination_over_ell"]) < abs(rows["0.125"]["even_combination_over_ell"])/20
    results[str(split)] = {"continuum_short_width_prediction":float(prediction),"rows":rows,"high_cutoff_control":check}

finite = {str(ell):cap_read(ell,.5,24,False) for ell in (.004,.001,.00025)}
assert finite["0.00025"]["mixed_log_jet_over_ell"] < finite["0.001"]["mixed_log_jet_over_ell"]/3
print(json.dumps({"trace_coefficient_per_sign":coefficient,"old_P0_and_free_S_reference":True,
                  "continuum_first":results,"fixed_24_modes_each_sign":finite,
                  "direct_full_matrix_controls":dense_controls(),
                  "scope":"C58 conditional old Slater cap with external native weak/physical-Y link. Width ell is collar width; no stochastic original H_B or physical OS map inferred."}, indent=2))
