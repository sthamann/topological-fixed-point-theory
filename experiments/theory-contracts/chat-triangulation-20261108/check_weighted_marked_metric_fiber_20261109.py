"""Evaluate the known weighted marked R032 metric jets and cut data controls."""
import contextlib
import io
import json
import runpy
import sys
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.linalg import expm, solve

repo = Path(sys.argv[1]).resolve()
sys.path.insert(0,str(repo / "verification"))
base = repo / "experiments/theory-contracts/chat-triangulation-20261108"


def run_file(name, args):
    path = base / name
    old_args = sys.argv
    try:
        sys.argv = [str(path), *map(str,args)]
        with contextlib.redirect_stdout(io.StringIO()):
            return runpy.run_path(str(path))
    finally:
        sys.argv = old_args


source = run_file("check_original_whole_two_strips_20261109.py",[repo])
compose, strip, steps = (source[k] for k in ("compose","strip","strip_steps"))
I, zero = (source[k] for k in ("identity","zero"))
frames, links = (source[k] for k in ("frames","links"))
left, right = (source[k] for k in ("charged_left","charged_right"))


def readouts(jet,sea):
    r,r1,r2 = jet
    inverse = solve(I+r,I)
    inv_r1 = inverse@r1
    cov = I-inverse
    dcov = inv_r1@inverse
    ddcov = inverse@r2@inverse-2*inv_r1@inv_r1@inverse
    lx = source["traces"](jet,sea)
    return lx,np.array([cov,dcov,ddcov])


def selected(rows):
    return [[float(z.real),float(z.imag)] for z in rows]


branches = {}
for branch in ("A","S","T"):
    matrices,seas = steps(branch)
    whole,whole_sea = strip(matrices,seas,0,len(matrices))
    lx,C = readouts(whole,whole_sea)
    branches[branch]={"log_trace":lx,"covariance_jet":C,"splits":{}}
    for cut in (2,3):
        early,se = strip(matrices,seas,0,cut)
        late,sl = strip(matrices,seas,cut,len(matrices))
        glued = compose(late,compose((links[cut-1],zero,zero),early))
        lx_split,C_split = readouts(glued,se+sl)
        err = np.max(np.abs(C_split[:,left][:,:,right]-C[:,left][:,:,right]),axis=(1,2))
        assert max(err) < 2e-12
        branches[branch]["splits"][cut] = (lx_split,C_split,err)


def weighted_lr(logs,Cjet):
    values = np.array([frames[-1]@C@frames[0].conj().T for C in Cjet])
    c0,c1,c2 = values[:,left[0],right[0]]
    log_z,score,hessian = logs
    Z=np.exp(log_z)
    return (np.array([c0,c1,c2]),
            Z*np.array([c0, c1+score*c0, c2+2*score*c1+(hessian+score**2)*c0]))


A = branches["A"]
S = branches["S"]
T = branches["T"]
C,marked_S = weighted_lr(A["log_trace"]-S["log_trace"], A["covariance_jet"])
_,marked_T = weighted_lr(A["log_trace"]-T["log_trace"], A["covariance_jet"])
errors = {}
for cut in (2,3):
    la,ca,raw_error = A["splits"][cut]
    ls = S["splits"][cut][0]
    _,split_marked = weighted_lr(la-ls,ca)
    err = np.abs(split_marked-marked_S)
    assert max(err) < 2e-18
    errors[cut] = {"whole_charged_covariance_jet_max_errors": list(map(float,raw_error)),
                   "marked_A_over_S_jet_errors":list(map(float,err))}


def direct_metric(m):
    r,_ = source["direct_whole"]("A",m)
    c = solve(I+r,r)[left[0],right[0]]
    a_over_s = np.exp(source["direct_log"]("A",m)-source["direct_log"]("S",m))
    return a_over_s*c


step=.001
f2m,f1m,f0,f1p,f2p = (direct_metric(k*step) for k in (-2,-1,0,1,2))
direct = np.array([f0,(f2m-8*f1m+8*f1p-f2p)/(12*step),
                  (-f2m+16*f1m-30*f0+16*f1p-f2p)/(12*step**2)])
assert max(np.abs(direct-marked_S)) < 2e-13

# Passive test: move the intermediate weak basis without changing physics.
# Its link must acquire the opposite transport, including both derivatives.
real_action = source["source"]["real_action"]
generator = np.diag([1.,-1.])
J1 = (real_action(1j*generator)-real_action(-1j*generator))/2
J2 = (real_action(-np.eye(2))-I)/2
assert np.max(np.abs(J1.conj().T+J1)) < 1e-14
assert np.max(np.abs(J2-J1@J1)) < 1e-14
matrices,seas = steps("A")
whole,whole_sea = strip(matrices,seas,0,len(matrices))
frame_control = {}
for cut in (2,3):
    early,_ = strip(matrices,seas,0,cut)
    late,_ = strip(matrices,seas,cut,len(matrices))
    link = (links[cut-1],zero,zero)
    moved_early = compose((I,-J1,J2),early)
    moved_link = compose(link,(I,J1,J2))
    corrected = compose(late,compose(moved_link,moved_early))
    wrong = compose(late,compose(link,moved_early))
    jet_errors = [np.max(np.abs(a-b)) for a,b in zip(corrected,whole)]
    assert max(jet_errors) < 1e-12
    raw_error = frames[-1]@(wrong[1]-corrected[1])@frames[0].conj().T
    log_error = (source["traces"](wrong,whole_sea,False)
                 - source["traces"](corrected,whole_sea,False))
    frame_control[cut] = {
        "passively_moved_full_transfer_jet_errors": list(map(float,jet_errors)),
        "omitting_link_derivatives_charged_transfer_first_error":
            selected([raw_error[left[0],right[0]]])[0],
        "omitting_link_derivatives_closed_log_trace_errors": list(map(float,log_error)),
    }
assert abs(frame_control[3]["omitting_link_derivatives_closed_log_trace_errors"][1]) < 1e-12
assert abs(frame_control[3]["omitting_link_derivatives_closed_log_trace_errors"][2]) > 1

# The actual QFT002 color-center selector alone on the closed fixed history.
# The signed native quark/antiquark labels set its phase; lepton ends are neutral.
native = source["native"]
labels = native.P
center_charge = np.array([
    sign if label[0] in ("u", "d") else 0
    for sign in (-1, 1)
    for label in labels
])
assert len(center_charge) == len(I)
assert not np.any(center_charge[left+right])
omega = np.exp(2j*np.pi/3)
swap = source["source"]["swap"]
center_left = native.rep(1,np.eye(2),omega*np.eye(3))
center_native = center_left@swap@center_left.conj()@swap
center_phases = np.diag(center_native)
assert np.max(np.abs(center_native-np.diag(np.exp(2j*np.pi*center_charge/3)))) < 1e-14


def twisted_determinant_jet(matrix_jet, phase):
    r0,r1,r2 = matrix_jet
    x = I+phase*r0
    phase_det,log_det = np.linalg.slogdet(x)
    y = solve(x,phase*r1)
    return phase_det, log_det, np.trace(y), np.trace(solve(x,phase*r2)-y@y)


def color_trace_jet(matrix_jet, sea_jet):
    """One color-center insertion at the complete exterior trace."""
    assert max(np.max(np.abs((center_phases[:,None]-center_phases[None,:])*r))
               for r in matrix_jet) < 1e-12
    traces = np.zeros(3,dtype=complex)
    for power in range(3):
        phase = (center_phases**power)[:,None]
        phase_det,log_det,first_log,second_log = twisted_determinant_jet(matrix_jet,phase)
        traces += phase_det*np.exp(log_det)*np.array(
            [1, first_log, second_log+first_log**2])/3
    assert traces[0].real > 0
    assert np.max(np.abs(traces.imag)) < 1e-11*traces[0].real
    v,d1,d2 = traces.real
    return sea_jet + np.array([np.log(v),d1/v,d2/v-(d1/v)**2])


def tiny_center_correction(matrix_jet):
    """Stable log jet of 3 Z_center / Z_untwisted, before rounding to 1."""
    sig0,log0,d0,dd0 = twisted_determinant_jet(matrix_jet,1)
    sigz,logz,dz,ddz = twisted_determinant_jet(matrix_jet,center_phases[:,None])
    ratio = sigz/sig0*np.exp(logz-log0)
    first_log,second_log = dz-d0,ddz-dd0
    value,first,second = 2*np.real(ratio*np.array(
        [1,first_log,second_log+first_log**2]))
    return ratio, np.array([np.log1p(value),first/(1+value),
                           second/(1+value)-(first/(1+value))**2])


color_logs = {}
color_errors = {}
center_twist_ratios = {}
center_corrections = {}
for branch in ("A","S","T"):
    matrices,seas = steps(branch)
    whole,whole_sea = strip(matrices,seas,0,len(matrices))
    color_logs[branch] = color_trace_jet(whole,whole_sea)
    undressed = tuple(frames[-1]@jet@frames[0].conj().T for jet in whole)
    ratio,center_corrections[branch] = tiny_center_correction(undressed)
    center_twist_ratios[branch] = selected([ratio])[0]
    r = undressed[0]
    for power in (1,2):
        twisted_r = (center_phases**power)[:,None]*r
        twisted_c = solve(I+twisted_r,twisted_r)
        plain_c = solve(I+r,r)
        assert abs(twisted_c[left[0],right[0]]-plain_c[left[0],right[0]]) < 5e-13
    color_errors[branch] = {}
    for cut in (2,3):
        early,se = strip(matrices,seas,0,cut)
        late,sl = strip(matrices,seas,cut,len(matrices))
        glued = compose(late,compose((links[cut-1],zero,zero),early))
        error = np.max(np.abs(color_trace_jet(glued,se+sl)-color_logs[branch]))
        assert error < 1e-12
        color_errors[branch][cut] = float(error)
color_relative = color_logs["A"]-color_logs["S"]
color_relative_T = color_logs["A"]-color_logs["T"]
_,color_marked = weighted_lr(color_relative,A["covariance_jet"])


def direct_color_relative(m):
    log_values = {}
    for branch in ("A","S"):
        r,sea = source["direct_whole"](branch,m)
        log_values[branch] = color_trace_jet((r,zero,zero),np.array([sea,0.,0.]))[0]
    return log_values["A"]-log_values["S"]


f2m,f1m,f0,f1p,f2p = (direct_color_relative(k*step) for k in (-2,-1,0,1,2))
direct_color_jet = [f0,(f2m-8*f1m+8*f1p-f2p)/(12*step),
                    (-f2m+16*f1m-30*f0+16*f1p-f2p)/(12*step**2)]
assert np.max(np.abs(direct_color_jet-color_relative)) < 1e-6

tiny_center_relative = center_corrections["A"]-center_corrections["S"]


def direct_tiny_center(m):
    a,_ = source["direct_whole"]("A",m)
    s,_ = source["direct_whole"]("S",m)
    return tiny_center_correction((a,zero,zero))[1][0]-tiny_center_correction((s,zero,zero))[1][0]


f2m,f1m,f0,f1p,f2p = (direct_tiny_center(k*step) for k in (-2,-1,0,1,2))
tiny_center_check = [f0,(f2m-8*f1m+8*f1p-f2p)/(12*step),
                     (-f2m+16*f1m-30*f0+16*f1p-f2p)/(12*step**2)]
assert np.max(np.abs(tiny_center_check-tiny_center_relative)) < 1e-24

# Existing prior test supplies two reflected field histories with identical
# Higgs radii at each time and their retained color Z3 fiber (not full Padm).
windows=run_file("check_r113_two_kernel_windows.py",
                 [base/"check_r113_ordered_reference.py",base/"check_quarter_reference.py"])
paired=windows["paired"]
positions=paired.charged[paired.connected]
tags=paired.source.native.P
charge=np.array([(1 if i<len(tags) else -1) if tags[i%len(tags)][0] in ("u","d") else 0
                 for i in positions])
ident=np.eye(len(paired.G))
center=np.diag(np.exp(2j*np.pi*charge/3))
color_factors=(ident,center,center@center)
y=-1j*paired.G@paired.Y


def z3_log_trace(R):
    value=0j
    for W in color_factors:
        sign, logdet=np.linalg.slogdet(ident+W@R)
        value += sign*np.exp(logdet)/3
    assert abs(value.imag) < 1e-10*value.real and value.real>0
    return np.log(value.real)


def required_cut_factor(m,x_first,x_middle):
    p=.25*np.exp(-m)
    A1=expm(-.2*(p*paired.G+x_first*y))
    B1=expm(-.2*(p*paired.G+x_middle*y))
    # Two mirrored halves, both reusing original letters and time.
    R1=B1@A1
    R2=A1@B1
    R=R2@R1
    # This is the factor that would repair separately closed halves.
    # The additive sea factors cancel exactly in this diagnostic.
    return z3_log_trace(R)-z3_log_trace(R1)-z3_log_trace(R2)


radial_rows={}
for name,h in (("original_early_positive",.35),("same_norm_early_negative",-.35)):
    fm2,fm1,f0,fp1,fp2=(required_cut_factor(k*step,h,-.55) for k in (-2,-1,0,1,2))
    radial_rows[name]=[float(f0),float((fm2-8*fm1+8*fp1-fp2)/(12*step)),
                      float((-fm2+16*fm1-30*f0+16*fp1-fp2)/(12*step**2))]
gap=np.array(radial_rows["original_early_positive"])-np.array(radial_rows["same_norm_early_negative"])
assert abs(gap[0]) > 0.1 and abs(gap[1]) > .001

# Algebraic co-moving measure derivatives and mixed marked terms,
# without assigning values to the original unknown measure score.
m = sp.symbols("m")
sigma,kappa = sp.symbols("sigma kappa")
b0,b1,b2 = sp.symbols("b0 b1 b2")
density = sp.exp(sigma*m+kappa*m*m/2)
block = b0+b1*m+b2*m*m/2
assert sp.simplify(sp.diff(density*block,m,2).subs(m,0)-
                   (b2+2*sigma*b1+(sigma*sigma+kappa)*b0))==0

output={
 "metric_coordinate":"m=log radius, fixed duration and original holonomy",
 "same_fixed_original_R032_history": {
    "L0_R0_C_jet":selected(C),"relative_weighted_L0_R0_A_over_S_jet":selected(marked_S),
    "relative_weighted_L0_R0_A_over_T_jet":selected(marked_T),
    "relative_vacuum_A_over_S":float(np.exp(A["log_trace"][0]-S["log_trace"][0])),
    "relative_vacuum_log_A_over_S_jet":list(map(float,A["log_trace"]-S["log_trace"])),
    "relative_vacuum_log_A_over_T_jet":list(map(float,A["log_trace"]-T["log_trace"])),
    "direct_five_point_weighted_A_over_S_jet":selected(direct),
    "split_vs_whole":errors},
 "passive_intermediate_frame_control":frame_control,
 "original_color_center_factor_only": {
    "relative_vacuum_A_over_S":float(np.exp(color_relative[0])),
    "relative_vacuum_log_A_over_S_jet":list(map(float,color_relative)),
    "relative_vacuum_log_A_over_T_jet":list(map(float,color_relative_T)),
    "direct_five_point_relative_vacuum_log_A_over_S_jet":list(map(float,direct_color_jet)),
    "relative_weighted_L0_R0_A_over_S_jet":selected(color_marked),
    "twisted_exterior_over_untwisted_det":center_twist_ratios,
    "center_correction_to_untwisted_A_over_S_log_jet":list(map(float,tiny_center_relative)),
    "direct_center_correction_to_untwisted_A_over_S_log_jet":list(map(float,tiny_center_check)),
    "max_log_jet_errors_at_both_cuts":color_errors},
 "radial_cut_correction_counterexample":{
    "definition":"log[Z_full/(Z_half_1 Z_half_2)] in the prior conditional R032 Z3 branch",
    "same_radius_and_same_metric_distinct_Higgs_orientation":radial_rows,
    "differences_in_log_correction_and_first_two_metric_jets":list(map(float,gap))},
 "full_physical_original_fiber_evaluated":False}
print(json.dumps(output,indent=2))
