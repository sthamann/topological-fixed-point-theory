"""Radius sweep on the existing C58 windows and fixed Slater reference."""
from contextlib import redirect_stdout
from pathlib import Path
import io
import json
import runpy
import sys

import numpy as np
from scipy.integrate import quad
from scipy.special import zeta

repo = Path(sys.argv[1]).resolve()
base = repo/"experiments/theory-contracts/chat-triangulation-20261108"
with redirect_stdout(io.StringIO()):
    c58 = runpy.run_path(str(base/"check_tfpt_c58_higgs_time_jet_20261110.py"))

geometry, cross, free_cross = (c58[k] for k in ("geometry","cross","free_cross"))
higgs, mass = c58["higgs"],c58["majorana"]
early, late, reflected = c58["left"],c58["shifted"],c58["shifted_reflected"]
width = c58["width"]
blocks = [i for i,_ in c58["old"]["blocks"]]
mass2 = c58["radial"]["base_masses"]
strength = np.trace(higgs@higgs).real/2


def free_sum(radius, later):
    # Exact summed quarter lattice for s<t:
    # sum over n exp(-2|n+1/4|(t-s)/r)=1/(2 sinh((t-s)/(2r))).
    a,b = early
    c,d = later
    return strength*quad(lambda s:quad(lambda t:1/(2*np.sinh((t-s)/(2*radius))),c,d)[0],a,b)[0]


def values(radius, cutoff):
    momenta = (np.arange(-cutoff,cutoff)+.25)/radius
    windows = (late,reflected)
    jets = np.array([free_sum(radius,window) for window in windows])
    for p in momenta:
        jets -= [free_cross(p,early,window) for window in windows]
        for index in blocks:
            eig, occupied, v = geometry(index,p,mass)
            jets += [cross(eig,occupied,v,early,window) for window in windows]
    a, reverse = jets
    p = abs(momenta)[:,None]
    E = np.sqrt(p*p+mass2)
    log_empty = (np.log1p(mass2/(2*E*(E+p))*np.expm1(-2*width*E))
                 -width*mass2**2/(2*p*(E+p)**2)).sum()
    ref_s = (free_sum(radius,late)-free_sum(radius,reflected))/2
    ref_t = (a+reverse)/2
    return np.array([a,a-ref_s,a-ref_t,log_empty])


radii = (.5,1.,2.)
coarse = {r:values(r,512) for r in radii}
fine = {r:values(r,1024) for r in radii}
for r in radii:
    p = (np.arange(-1024,1024)+.25)/r
    np.testing.assert_allclose(sum(free_cross(k,early,late) for k in p),
                               free_sum(r,late),atol=2e-10)
    np.testing.assert_allclose(coarse[r][:3],fine[r][:3],atol=3e-9,rtol=0)
np.testing.assert_allclose(fine[1.][0],.9683881172,atol=1e-10)
assert max(fine[r][0]/r for r in radii)-min(fine[r][0]/r for r in radii)>1e-4

# The empty determinant converges with a power tail. Use the explicit
# C58 source expansion in the rescaled p=(m+1/4)/radius variable.
def empty_tail(radius,cutoff):
    z = lambda k:radius**k*(zeta(k,cutoff+.25)+zeta(k,cutoff+.75))
    a,b,c,d,e = [(mass2**i).sum() for i in range(1,6)]
    return (-z(2)*a/4 + 5*z(4)*b/32 - 11*z(6)*c/96 + 93*z(8)*d/1024
            +width*(-z(3)*b/8+z(5)*c/16-5*z(7)*d/128+7*z(9)*e/256))


for r in radii:
    coarse[r][3] += empty_tail(r,512)
    fine[r][3] += empty_tail(r,1024)
    np.testing.assert_allclose(coarse[r],fine[r],atol=3e-9,rtol=0)
np.testing.assert_allclose(fine[1.][3],-10.6101199202,atol=1e-9)
h = 1e-3
lower, upper = (values(float(np.exp(s)),512) for s in (-h,h))
lower[3] += empty_tail(np.exp(-h),512)
upper[3] += empty_tail(np.exp(h),512)
derivatives = (upper-lower)/(2*h)
np.testing.assert_allclose(derivatives,[.96856031,1.23099166,-.26240682,-12.9691542430],
                           atol=5e-6)
print(json.dumps({
    "scope":"Conditional C58 original fixed-background Slater cap; same width, marks, quarter, subtraction and S/T maps at all radii.",
    "radius_columns":["A","A_minus_S","A_minus_T","log_empty_A_over_S"],
    "radii":{str(r):list(map(float,fine[r])) for r in radii},
    "A_over_r":{str(r):float(fine[r][0]/r) for r in radii},
    "local_log_radius_derivatives_at_one":list(map(float,derivatives)),
    "local_log_slope_of_A":float(derivatives[0]/fine[1.][0]),
    "exact_summed_free_mode_kernel":"1/(2 sinh((t-s)/(2r))) for t>s",
    "physical_3plus1_spin2":None,
},indent=2))
