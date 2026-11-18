"""Original full Pin D-sea, mass rotation and its two-dimensional UV test."""

import contextlib
import io
import json
from pathlib import Path
import runpy
import sys

import numpy as np


repo = Path(sys.argv[1]).resolve()
base = repo/"experiments/theory-contracts/chat-triangulation-20261108"
with contextlib.redirect_stdout(io.StringIO()):
    source = runpy.run_path(str(base/"check_tfpt_conditioned_r113_weights_20261110.py"))
old = source["majorana"]
higgs = source["inner"](np.array([0.,1.]))
new = old+.55*higgs
old_mass,old_basis = np.linalg.eigh(old)
new_mass,new_basis = np.linalg.eigh(new)
cross = abs(old_basis.conj().T@new_basis)**2
mass_difference_squared = float(np.linalg.norm(new-old)**2)
prediction = mass_difference_squared/np.pi
mass_delta = abs(old_mass[:,None]-new_mass[None,:])
odd_error = (abs(old_mass[:,None])**3+abs(new_mass[None,:])**3)/2
even_error = abs(old_mass[:,None]**2-new_mass[None,:]**2)/2
error_4 = np.sum(cross*(even_error**2/2+mass_delta*odd_error))
error_6 = np.sum(cross*odd_error**2/2)


def one_mode(lam):
    e = np.hypot(lam,old_mass)
    f = np.hypot(lam,new_mass)
    # 1-dot(u,v) = |u-v|²/2, for u=(lambda,s)/hypot(lambda,s).
    # This avoids subtracting nearly equal numbers in the UV.
    dif1 = lam*(1/e[:,None]-1/f[None,:])
    dif2 = old_mass[:,None]/e[:,None]-new_mass[None,:]/f[None,:]
    return float(np.sum(cross*(dif1**2+dif2**2)/2))


def direct(lam):
    eye = np.eye(len(old))
    def P(mass):
        full = np.block([[lam*eye,mass],[mass,-lam*eye]])
        levels,vectors = np.linalg.eigh(full)
        occupied = vectors[:,levels<0]
        return occupied@occupied.conj().T
    return float(np.linalg.norm(P(old)-P(new))**2)


for lam in (np.pi/2,10.):
    np.testing.assert_allclose(one_mode(lam),direct(lam),atol=1e-13)

def lattice(limit, zero_only=False):
    n = np.arange(-int(limit/(2*np.pi))-2,int(limit/(2*np.pi))+3)
    m = np.arange(0,1 if zero_only else int(2*limit/np.pi)+1)
    alpha = np.where(m%2==0,.25,.75)
    k = 2*np.pi*(n[None,:]+alpha[:,None])
    e = np.hypot(k,np.pi*m[:,None]/2)
    weights = np.where(m[:,None]==0,1,2)+np.zeros_like(e,dtype=int)
    e = e.ravel()
    weights = weights.ravel()
    return e[e<=limit], weights[e<=limit]


# Exact low modes and explicit leading high-mode contribution with a bound.
# For A,B finite, |s/e-t/f| errors are O(lambda^-3) and lambda/e difference
# O(lambda^-2); its squared remainder gives O(lambda^-4) per geometric mode.
rows = {}
cutoffs = (32.,64.,128.,256.,512.,1024.)
low,low_weights = lattice(48.)
low_exact = sum(one_mode(x)*w for x,w in zip(low,low_weights,strict=True))
for cutoff in cutoffs:
    e,weights = lattice(cutoff)
    low_mask = e<=48.
    if cutoff<48.:
        hs = sum(one_mode(x)*w for x,w in zip(e,weights,strict=True))
    else:
        hs = low_exact+mass_difference_squared/2*np.sum(weights[~low_mask]/e[~low_mask]**2)
    line,line_w = lattice(cutoff,True)
    line_hs = sum(one_mode(x)*w for x,w in zip(line,line_w,strict=True))
    rows[str(cutoff)] = {"paired_projection_HS_square_full":float(hs),
                         "paired_projection_HS_square_normal_zero_only":float(line_hs),
                         "Dirac_modes_including_signs":int(weights.sum()),
                         "absolute_full_tail_approximation_error_bound":float(np.sum(
                             weights[~low_mask]*(error_4/e[~low_mask]**4+error_6/e[~low_mask]**6)))}
last = (rows["1024.0"]["paired_projection_HS_square_full"]
       -rows["512.0"]["paired_projection_HS_square_full"])/np.log(2)
assert abs(last/prediction-1)<.01
e,w = lattice(1024.)
shell = e>512.
slope_error = np.sum(w[shell]*(error_4/e[shell]**4+error_6/e[shell]**6))/np.log(2)
# A fixed high-mode term and a short exact shell control.
for lam in (64.,128.,256.):
    np.testing.assert_allclose(one_mode(lam),mass_difference_squared/(2*lam**2),rtol=.005)

k = np.pi/2
# For each actual mass, rotate sign(k)*tau3 into sign of the massive block.
theta = np.arctan(new_mass/k)
u = np.array([[np.cos(theta/2),-np.sin(theta/2)],
              [np.sin(theta/2),np.cos(theta/2)]]).transpose(2,0,1)
blocks = k*np.diag([1.,-1.])+new_mass[:,None,None]*np.array([[0.,1.],[1.,0.]])
diag = u.transpose(0,2,1)@blocks@u
expected = np.hypot(k,new_mass)[:,None,None]*np.diag([1.,-1.])
np.testing.assert_allclose(diag,expected,atol=1e-13)

print(json.dumps({"norm_difference_squared_native_mass_block":mass_difference_squared,
                  "full_flat_pin_logarithmic_coefficient":prediction,
                  "cutoff_rows":rows,"last_cutoff_doubling_log_slope":float(last),
                  "last_log_slope_approximation_error_bound":float(slope_error),
                  "single_mode_rotation_error":float(np.max(abs(diag-expected))),
                  "smallest_old_geometric_frequency":float(k),
                  "native_new_minmax_on_same_normal_zero_mode":list(map(float,
                       (min(np.hypot(k,new_mass)),max(np.hypot(k,new_mass))))),
                  "scope":"Given flat paired Pin Dirac, native M_R to M_R+.55H, pure spectral-sea comparison on all longitudinal and transverse modes. Not the thermal Green time, a physical q mapping, or a general F1 exclusion."},indent=2))
