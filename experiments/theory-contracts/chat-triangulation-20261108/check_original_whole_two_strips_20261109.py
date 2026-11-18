"""Sew the unchanged R032 fixed history with the full exterior intermediate fiber."""
import contextlib
import io
import json
import runpy
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import eigh, expm, solve

repo = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(repo / "verification"))
base = repo / "experiments/theory-contracts/chat-triangulation-20261108"
source_path = base / "check_r113_dressed_open_20261109.py"
saved_args = sys.argv
try:
    sys.argv = [str(source_path), str(base)]
    with contextlib.redirect_stdout(io.StringIO()):
        source = runpy.run_path(str(source_path))
finally:
    sys.argv = saved_args

native = source["native"]
grading = source["grading"]
majorana = source["majorana"]
history = source["history"]
durations = source["durations"]
frames = source["frames"]
links = source["dressed_links"]
identity = np.eye(len(grading))
zero = np.zeros_like(identity)
charged_left = [native.IDX[("e", "L", 0, g)] for g in range(3)]
charged_right = [native.IDX[("e", "R", 0, g)] for g in range(3)]
charged = charged_left + charged_right


def generator(i, branch, radius_log):
    phi = source["inner_form"](history[i]) + majorana if branch != "S" else zero
    p = source["p"] * np.exp(-radius_log)
    return p * grading - 1j * grading @ phi, phi, p


def direct_whole(branch, radius_log):
    """Independent undressed direct ordered product via Hermitian spectra."""
    matrix = identity.copy()
    sea = 0.
    for i, width in enumerate(durations):
        h, _, _ = generator(i, branch, radius_log)
        values, vectors = eigh(h)
        heat = (vectors * np.exp(-width * values)) @ vectors.conj().T
        matrix = heat @ matrix
        sea += width * values[:len(values)//2].sum()
    return matrix, sea


def compose(left, right):
    """Multiplication of value, first derivative and second derivative jets."""
    l0, l1, l2 = left
    r0, r1, r2 = right
    return (l0 @ r0, l1 @ r0 + l0 @ r1, l2 @ r0 + 2*l1 @ r1 + l0 @ r2)


def strip_steps(branch):
    mats, seas = [], []
    for i, width in enumerate(durations):
        h, phi, p = generator(i, branch, 0.)
        frame = frames[i]
        h = frame.conj().T @ h @ frame
        a = -width*h
        da = width*p*grading
        block = np.block([[a, da, -da/2], [zero, a, da], [zero, zero, a]])
        exponential = expm(block)
        size = len(h)
        mats.append((exponential[:size,:size], exponential[:size,size:2*size],
                     2*exponential[:size,2*size:]))
        squares = np.maximum(np.linalg.eigvalsh(phi @ phi),0)
        energies = np.sqrt(p*p+squares)
        # E_sea is minus half the absolute spectrum on the actual symmetric fiber.
        sea = width/2 * np.array([-energies.sum(),
            (p*p/energies).sum(), -(2*p*p/energies-p**4/energies**3).sum()])
        seas.append(sea)
    return mats, seas


def strip(mats, seas, first, last):
    result = mats[first]
    for j in range(first+1,last):
        result = compose(mats[j], compose((links[j-1],zero,zero),result))
    return result, sum(seas[first:last])


def traces(jet, sea, require_positive=True):
    r0, r1, r2 = jet
    w = identity+r0
    sign, logdet = np.linalg.slogdet(w)
    if require_positive:
        assert abs(sign-1) < 3e-12
    f = solve(w,r1)
    return sea + np.array([logdet, np.trace(f).real,
                           np.trace(solve(w,r2)-f@f).real])


def direct_log(branch, radius_log):
    matrix, sea = direct_whole(branch, radius_log)
    return sea + np.linalg.slogdet(identity+matrix)[1]


def absmax(m):
    return float(np.max(np.abs(m)))


def complex_pair(z):
    return [float(z.real),float(z.imag)]


all_rows = {}
whole_traces = {}
two_by_branch = {}
early_and_late = {}
for branch in ("A","S","T"):
    mats,seas = strip_steps(branch)
    whole, whole_sea = strip(mats,seas,0,len(history))
    whole_trace = traces(whole,whole_sea)
    whole_traces[branch] = whole_trace
    rows = []
    direct, direct_sea = direct_whole(branch,0.)
    from_undressed = frames[-1] @ whole[0] @ frames[0].conj().T
    assert absmax(from_undressed-direct) < 5e-12
    assert abs(whole_sea[0]-direct_sea) < 5e-12
    for cut in (2,3):
        early, early_sea = strip(mats,seas,0,cut)
        late, late_sea = strip(mats,seas,cut,len(history))
        glued = compose(late,compose((links[cut-1],zero,zero),early))
        glued_sea = early_sea+late_sea
        glued_trace = traces(glued,glued_sea)
        delta = [absmax(x-y) for x,y in zip(glued,whole)]
        assert max(delta) < 1e-11
        assert absmax(glued_trace-whole_trace) < 1e-11
        rows.append({"cut_time":sum(durations[:cut]),"open_matrix_jet_errors":delta,
                     "sea_jet_error":absmax(glued_sea-whole_sea),
                     "exterior_trace_log_jet_errors":list(map(float,abs(glued_trace-whole_trace)))})
        if cut==3:
            two_by_branch[branch] = glued_trace
            early_and_late[branch] = (early,early_sea,late,late_sea,glued,glued_sea)
    all_rows[branch]=rows

relative_s = whole_traces["A"]-whole_traces["S"]
relative_t = whole_traces["A"]-whole_traces["T"]
assert absmax(relative_t) < 1e-12
assert absmax(relative_s-(two_by_branch["A"]-two_by_branch["S"])) < 1e-12
step = .001
def d(q):
    return direct_log("A",q)-direct_log("S",q)
fm2,fm1,f0,fp1,fp2 = [d(j*step) for j in (-2,-1,0,1,2)]
finite_diff = np.array([f0, (fm2-8*fm1+8*fp1-fp2)/(12*step),
                        (-fm2+16*fm1-30*f0+16*fp1-fp2)/(12*step**2)])
assert absmax(finite_diff-relative_s) < 4e-7

early,early_sea,late,late_sea,glued,glued_sea = early_and_late["A"]
r_dressed = glued[0]
r = frames[-1] @ r_dressed @ frames[0].conj().T
# Coherent open kernel coefficient and traced charged CAR value in the
# same exterior state, separate from the full physical OS projection.
open_lr = r[np.ix_(charged_left,charged_right)]
cov = solve(identity+r,r)
normalized_lr = cov[np.ix_(charged_left,charged_right)]
early_un = frames[2] @ early[0] @ frames[0].conj().T
late_un = frames[-1] @ late[0] @ frames[3].conj().T
assert absmax(late_un-early_un.conj().T) < 5e-12
drop = r[np.ix_(charged,charged)]-late_un[np.ix_(charged,charged)] @ early_un[np.ix_(charged,charged)]
assert min(np.linalg.eigvalsh(drop)) > -1e-10
assert np.linalg.norm(drop) > .01
assert absmax(open_lr-np.array(r)[np.ix_(charged_left,charged_right)]) < 1e-14
separately_closed = {}
for branch,(e,es,l,ls,_,_) in early_and_late.items():
    e_un = tuple(frames[2] @ x @ frames[0].conj().T for x in e)
    l_un = tuple(frames[-1] @ x @ frames[3].conj().T for x in l)
    assert max(absmax(x-y.conj().T) for x,y in zip(e_un,l_un)) < 5e-12
    # For mirrored halves the determinant phases and imaginary jets cancel.
    separately_closed[branch] = traces(e_un,es,False)+traces(l_un,ls,False)
bad = separately_closed["A"]-separately_closed["S"]
assert abs(bad[0]-relative_s[0]) > .05

output={
 "scope":"Unchanged R032 fixed reflected v252 diagnostic field history, full finite exterior fiber including sea, no dynamical Higgs integration or full physical P_adm.",
 "duration":sum(durations),
 "same_history_splits":all_rows,
 "closed_fermion_trace":{
   "log_A":float(whole_traces["A"][0]),"log_S":float(whole_traces["S"][0]),
   "A_over_S":float(np.exp(relative_s[0])),"log_A_over_S_metric_jet":list(map(float,relative_s)),
   "log_A_over_T_metric_jet":list(map(float,relative_t)),
   "independent_direct_finite_differences_S":list(map(float,finite_diff))},
 "open_fermion_return":{"sea_log_A":float(glued_sea[0]),"sea_factor_A":float(np.exp(glued_sea[0])),
                        "charged_L0_R0_coherent_coefficient_divided_by_sea":complex_pair(open_lr[0,0]),
                        "charged_L0_R2_coherent_coefficient_divided_by_sea":complex_pair(open_lr[0,2]),
                        "charged_L0_R0_coherent_coefficient_with_sea":complex_pair(np.exp(glued_sea[0])*open_lr[0,0]),
                        "charged_L0_R0_exterior_normalized_covariance":complex_pair(normalized_lr[0,0])},
 "premature_charged_only_cut":{"return_defect_trace":float(np.trace(drop).real),
                                "return_defect_operator_norm":float(np.linalg.norm(drop,2)),
                                "return_defect_L0_R0":complex_pair(drop[0,len(charged_left)])},
 "wrong_independent_closures":{"log_A_over_S_metric_jet":list(map(float,bad)),
                                "difference_in_relative_metric_jet":list(map(float,bad-relative_s))},
}
print(json.dumps(output,indent=2))
