"""Distinguish C58 angle jets, actual loop central term and conditional fixed points."""

import contextlib
from fractions import Fraction as F
import io
import json
from pathlib import Path
import runpy
import sys

import numpy as np
from scipy.special import zeta


repo = Path(sys.argv[1]).resolve()
base = repo / "experiments/theory-contracts/chat-triangulation-20261108"
with contextlib.redirect_stdout(io.StringIO()):
    old = runpy.run_path(str(base / "check_tfpt_c58_short_time_mixed_cap_20261110.py"))
    windows = runpy.run_path(str(base / "check_tfpt_c58_higgs_time_jet_20261110.py"))

t, y = (old[key].diagonal().real for key in ("t3","y"))
gamma = old["gamma"].diagonal().real.astype(int)
charges = {"Q":t+y,"T3":t,"Y":y}
exact = {key:[F(float(a)).limit_denominator(6) for a in values] for key,values in charges.items()}
pairs = (("Q","Q"),("Q","T3"),("T3","T3"),("Y","Y"),("Y","T3"))


def trace(a,b,weight):
    return sum(int(w)*x*z for w,x,z in zip(weight,exact[a],exact[b],strict=True))


quadratics = {}
for a,b in pairs:
    full = trace(a,b,np.ones(96))
    graded = trace(a,b,gamma)
    assert graded==0
    # First harmonic on the actual quarter sea. Restricted Fourier
    # multiplication has one crossing of the sea per grade for each unit shift.
    crosses = 0
    for i,g in enumerate(gamma):
        for mode in range(-2,3):
            if g*(mode+.25)>=0:
                continue
            crosses += exact[a][i]*exact[b][i]*(int(g*(mode-1+.25)>0)-int(g*(mode+1+.25)>0))
    assert crosses==-graded
    quadratics[f"{a},{b}"] = {"ungraded_full96":str(full),"ungraded_per_family16":str(full/6),
                              "signed_actual_96":str(graded),"first_harmonic_loop_cocycle":str(crosses)}

short = old["boundary_coefficient"](.5)
jets = {"Q,Q":0.,"Q,T3":0.,"T3,T3":float(-short),"Y,Y":float(-short),"Y,T3":float(short)}
assert np.max(abs((old["t3"]+old["y"])@old["vl"]-old["vl"]@(old["t3"]+old["y"])))<1e-14
assert np.max(abs((old["t3"]+old["y"])@old["vr"]-old["vr"]@(old["t3"]+old["y"])))<1e-14

# Bare det2 on the same relative occupied operator B(H) subtracts
# Tr(B-I). At H=0 the mixed Hessian is -Tr(B_f B_g), with B(0)=I.
# Use the native full 96 jet, at separated non-mirror windows.
def det2_separated(first,second,cutoff,mass):
    def one_mark(e,occupied,v,mark):
        a,b=mark
        _,dv = windows["transfer_derivatives"](e,v,b-a)
        return ((occupied.conj().T*np.exp(-(windows["width"]-b)*e))
                @ dv @ (np.exp(-a*e)[:,None]*occupied))

    total=0.
    for sign,shift in ((1,.25),(-1,.75)):
        grade = sign*gamma
        leg = mass[np.ix_(grade<0,grade>0)]
        vh = windows["higgs"][np.ix_(grade<0,grade>0)]
        d_square = leg@vh.conj().T + vh@leg.conj().T
        integrated_marks = (first[1]-first[0])*(second[1]-second[0])
        total -= integrated_marks*np.linalg.norm(d_square)**2*zeta(2,cutoff+shift)/4
        for p in sign*(np.arange(cutoff)+shift):
            for index,_ in windows["old"]["blocks"]:
                e,occupied,v = windows["geometry"](index,p,mass)
                cap = occupied.conj().T@(np.exp(-windows["width"]*e)[:,None]*occupied)
                le,lf = (np.linalg.solve(cap,one_mark(e,occupied,v,mark)) for mark in (first,second))
                total -= np.trace(le@lf).real
    return float(total)


det2_jet = {}
for cutoff in (256,512,1024):
    det2_jet[str(cutoff)] = det2_separated(windows["left"],windows["shifted"],cutoff,windows["majorana"])
assert abs(det2_jet["512"]-det2_jet["1024"])<1e-7
assert abs(det2_jet["1024"]-windows["cases"]["shifted"]["1024"][0])>.1
massless = det2_separated(windows["left"],windows["shifted"],16,0*windows["majorana"])
assert massless==0.

# The exponent read back from a model's own heat-character moment is
# the parameter inserted. U2 channels (j=1,q=0), (j=0,q=2).
heat_fixed_points = []
for t,u in ((.05,.03),(.122,.0732),(.25,.15)):
    moments = np.exp(-np.array([2*t,4*u]))
    recovered = -np.log(moments)/np.array([2,4])
    np.testing.assert_allclose(recovered,[t,u],atol=1e-15)
    heat_fixed_points.append({"input":[t,u],"one_iteration":recovered.tolist()})
# All these marginals may also be kept with the prior equal-Haar
# symmetric jump pair. In j=1,q=2 the difference is one eighteenth
# per unit Poisson rate. Neither marginal records this direction.
print(json.dumps({"central_terms":quadratics,
                  "external_cartan_cap_width_coefficients":jets,
                  "det2_separated_Higgs_jet":{"raw_C58_with_Majorana":windows["cases"]["shifted"]["1024"][0],
                     "det2_with_Majorana_by_cutoff":det2_jet,"det2_without_Majorana_exact":massless,
                     "raw_C58_without_Majorana":windows["cases"]["shifted"]["1024"][1]},
                  "conditional_marginal_heat_fixed_points":heat_fixed_points,
                  "same_marginal_jump_coefficient_difference_per_rate":"1/18",
                  "scope":"Old native circle sea, fixed external Cartan cap, and explicitly conditional heat moment readback. No physical C58 current central charge or R113 weight claimed."},indent=2))
