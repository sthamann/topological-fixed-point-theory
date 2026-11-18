"""Window sum rule, quarter phase, thermal modular flow and native A2 clock."""

import ast
import contextlib
from itertools import combinations, product
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
    windows = runpy.run_path(str(base/"check_tfpt_c58_higgs_time_jet_20261110.py"))
    cover = runpy.run_path(str(base/"check_tfpt_pin_reflected_higgs_cut_20261110.py"))
h, m = (windows[k] for k in ("higgs","majorana"))
left, shifted, reflected = (windows[k] for k in ("left","shifted","shifted_reflected"))


def vacuum(radius, mass):
    p = np.abs(np.arange(-4096,4096)+.25)[:,None]/radius
    square = np.linalg.eigvalsh(mass)[48:]**2
    e = np.hypot(p,np.sqrt(square))
    vacancy = square/(2*e*(e+p))
    width = windows["width"]
    log_part = np.sum(np.log1p(vacancy*np.expm1(-2*width*e))
                      -width*square**2/(2*p*(e+p)**2))
    log_part -= square.sum()/4*radius**2*(zeta(2,4096.25)+zeta(2,4096.75))
    log_part -= width*(square**2).sum()/8*radius**3*(zeta(3,4096.25)+zeta(3,4096.75))
    return float(np.exp(log_part))


def window_reads(radius, amplitude, cutoff=768):
    values = np.zeros((2,2))
    mass = m+amplitude*h
    for sign, shift in ((1,.25),(-1,.75)):
        for p in sign*(np.arange(cutoff)+shift)/radius:
            for j,late in enumerate((shifted,reflected)):
                values[1,j] += windows["free_cross"](p,left,late)
            for indices,_ in windows["old"]["blocks"]:
                data = windows["geometry"](indices,p,mass)
                for j,late in enumerate((shifted,reflected)):
                    values[0,j] += windows["cross"](*data,left,late)
    a, free = values
    s = a[0]-(free[0]-free[1])/2
    t = a[0]-(a[0]+a[1])/2
    delta = s+t-a[0]
    assert abs(delta-((a[0]-free[0])-(a[1]-free[1]))/2)<1e-14
    v = vacuum(radius,mass)
    return {"A_S_T":[a[0],s,t],"S_plus_T_minus_A":delta,
            "modewise_S_vacuum_factor":v,"ratio":delta/v}


sum_rows = {f"R={radius},r={amplitude}":window_reads(radius,amplitude)
            for radius,amplitude in ((1.,0.),(1.,.35),(1.,.55),(.9,0.),(1.1,0.))}
assert abs(sum_rows["R=1.0,r=0.0"]["ratio"]-1)>.005

# Multicopy test adds the same native return probes, rather than replacing
# the observable with a product of returns.
re,im = cover["compress"](.1,64,.55,.2,.3)
phase_by_directions = {}
for count in (1,2,4,8,16):
    pair = (np.trace(re[:count,:count]).real,np.trace(im[:count,:count]).real)
    assert pair[1]>0
    phase_by_directions[str(count)] = float(np.degrees(np.arctan2(pair[1],pair[0])))
curve = {}
for u in (.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,.95):
    re_u,im_u = cover["compress"](.1,32,.55,u/3,2*u/3)
    curve[str(u)] = float(np.degrees(np.arctan2(np.trace(im_u).real,np.trace(re_u).real)))
curvature = curve["0.1"]-2*curve["0.2"]+curve["0.3"]
assert abs(curvature)>.05

energy,mark,_ = cover["spatial"](.1,12)
mass = float(np.linalg.eigvalsh(cover["phi_m"]+.55*cover["phi_h"])[-1])
k_values = np.linalg.eigvalsh(np.diag(energy)+mass*mark)
# Keep representative unsaturated levels for a direct odds comparison.
levels = k_values[abs(k_values)<10]
covariance = 1/(1+np.exp(2*levels))
modular = np.log1p(-covariance)-np.log(covariance)
assert max(abs(modular-2*levels))<1e-8
phi0 = 1/(6*np.pi)+3/(256*np.pi**4)
gap = cover["gaps"]["0.55"]["flat_spatial_Higgs_generator_gap"]

# Exact E8 roots in doubled coordinates, and the actual native permutation.
roots = []
for i,j in combinations(range(8),2):
    for a,b in product((-2,2),repeat=2):
        root = [0]*8
        root[i],root[j] = a,b
        roots.append(tuple(root))
roots += [x for x in product((-1,1),repeat=8) if x.count(-1)%2==0]
syntax = ast.parse((repo/"verification/v634_st31_structure.py").read_text())
idx = next(ast.literal_eval(node.value) for node in syntax.body if isinstance(node,ast.Assign)
           and any(isinstance(t,ast.Name) and t.id=="SIGMA_IDX" for t in node.targets))

def sigma(r):
    return tuple(r[i] for i in idx)

def dot(a,b):
    return sum(x*y for x,y in zip(a,b,strict=True))

def components(records):
    remaining = set(records)
    counts = []
    while remaining:
        root = remaining.pop()
        component, queue = {root},[root]
        while queue:
            here = queue.pop()
            following = {r for r in remaining if dot(r,here)!=0}
            remaining -= following
            component |= following
            queue += list(following)
        counts.append([len(component),int(np.linalg.matrix_rank(np.array(list(component))))])
    return sorted(counts)

fixed = [r for r in roots if sigma(r)==r]
beta = (0,0,0,2,0,-2,0,0)
orbit = [beta,sigma(beta),sigma(sigma(beta))]
assert all(dot(a,b)==-4 for a,b in combinations(orbit,2))
assert all(sum(x)==0 for x in zip(*orbit,strict=True))
centralizer_one = [r for r in roots if all(dot(r,a)==0 for a in orbit)]
assert components(fixed)==[[6,2],[6,2]] and components(centralizer_one)==[[72,6]]
d5 = [r for r in roots if all(x==0 for x in r[5:])]
assert len(d5)==40
d5_commuting = [r for r in d5 if r in centralizer_one]

# The UV coefficient of the original broken-angle collar response is
# the trace of the same native Higgs legs, times their two amplitudes.
gamma = windows["old"]["gamma"].diagonal().real
occupied = gamma<0
v_leg = h[np.ix_(occupied,~occupied)]
amplitude_pairs = ((.35,.55),(.175,.275),(.35,.35),(0.,.55))
stiffness = {str(pair):float(pair[0]*pair[1]*np.linalg.norm(v_leg)**2*np.log(2)/4)
             for pair in amplitude_pairs}
assert abs(stiffness[str((.35,.55))]-6.10840873983177)<1e-12

print(json.dumps({
    "C58_summation_check":sum_rows,
    "quarter_phase_by_number_of_native_directions":phase_by_directions,
    "quarter_phase_degrees_by_time_sum":curve,
    "nonzero_phase_second_difference_degrees":curvature,
    "modular":{"direct_max_error_to_2K":float(max(abs(modular-2*levels))),
               "direct_max_error_to_2piK":float(max(abs(modular-2*np.pi*levels))),
               "original_thermal_period":2,"spatial_gap":gap,
               "phi0_over_2":float(phi0/2),"difference":float(gap-phi0/2)},
    "native_sigma":{"fixed_root_components_size_rank":components(fixed),
                    "E6_orthogonal_to_one_A2_size_rank":components(centralizer_one),
                    "actual_D5_roots_commuting_with_that_A2":len(d5_commuting),
                    "glue_reflection_A2_orbit_doubled":orbit},
    "broken_angle_C58_Higgs_scaling":stiffness,
    "scope":"Old native C58 and static paired Pin checks, actual E8 sigma. No quantum boson measure, no physical modular inclusion or QTheta."
},indent=2))
