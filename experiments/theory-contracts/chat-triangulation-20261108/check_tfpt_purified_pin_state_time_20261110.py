"""State-preserving purification and the source time/parity obstruction."""

import contextlib
import io
import json
from pathlib import Path
import runpy
import sys

import numpy as np
from numpy.polynomial.legendre import leggauss


repo = Path(sys.argv[1]).resolve()
base = repo/"experiments/theory-contracts/chat-triangulation-20261108"
with contextlib.redirect_stdout(io.StringIO()):
    native = runpy.run_path(str(base/"check_tfpt_conditioned_r113_weights_20261110.py"))
majorana = native["majorana"]
higgs = native["inner"](np.array([0.,1.]))
values,basis = np.linalg.eigh(majorana)
P_zero = basis[:,abs(values)<1e-12]@basis[:,abs(values)<1e-12].conj().T
assert np.linalg.matrix_rank(P_zero)==90
cross = (np.eye(96)-P_zero)@higgs@P_zero
assert np.linalg.norm(cross)>0.1

# Both actual scalar quarter spectra at the ORIGINAL angular-to-x conversion.
# For the massless normal generator K=H_b, E_m=m*pi/2.
records = {}
for numerator in (1,3):
    n = np.arange(-5,6)
    source = 2*np.pi*(n+numerator/4)
    normal_labels = 4*n+numerator
    target = normal_labels*np.pi/2
    assert np.array_equal(source,target)
    assert np.all(normal_labels%2==1)
    old_p = np.diag(source<0).astype(float)
    c_sqrt = np.diag(np.exp(-.5*np.logaddexp(0.,2*target)))
    complement_sqrt = np.diag(np.exp(-.5*np.logaddexp(0.,-2*target)))
    U = np.block([[c_sqrt,-complement_sqrt],[complement_sqrt,c_sqrt]])
    occ = U[:,:len(n)]
    covariance = occ@occ.T
    W = U@np.vstack([old_p,np.eye(len(n))-old_p])
    ambient = np.diag(np.r_[target,target])
    errs = {"CAR_isometry":np.max(abs(W.T@W-np.eye(len(n)))),
            "old_covariance_restored":np.max(abs(W.T@covariance@W-old_p)),
            "cover_covariance_restored":np.max(abs(covariance[:len(n),:len(n)]-c_sqrt@c_sqrt)),
            "projection":np.max(abs(covariance@covariance-covariance)),
            "time_intertwining":np.max(abs(ambient@W-W@np.diag(source)))}
    assert max(errs.values())<1e-12
    records[f"{numerator}/4"] = {"normal_modes":normal_labels.tolist(),
                                "normal_parity_relative_to_b":-1,
                                "direct_errors":{k:float(v) for k,v in errs.items()}}

# The canonical cross covariance F=(2 cosh K)^-1 is nonlocal on the
# same original spatial strip. Use two disjoint normalized test functions
# in the first spin component, in one massless native internal direction.
nodes,weights = leggauss(96)

def smeared_modes(m, start, end):
    w = (start+end+(end-start)*nodes)/2
    dx = weights*(end-start)/2
    bump = np.exp(1-1/(1-nodes*nodes))
    bump /= np.sqrt(np.sum(dx*bump**2))
    return np.exp(1j*(m*np.pi/2)[:,None]*w)@(dx*bump)/2

locality = {}
for limit in (16,32,64):
    normal = np.arange(-limit,limit+1)
    f = smeared_modes(normal,-.8,-.6)
    g = smeared_modes(normal,.3,.5)
    locality[str(limit)] = float(np.sum(f.conj()*g/(2*np.cosh(normal*np.pi/2))).real)
assert abs(locality["32"]-locality["64"])<1e-13 and locality["64"]>.01

# Local half-sided invariance test on the actual massless spin propagation.
# Away from the endpoint H_b=i*sigma3*d/dw shifts the two components
# in opposite directions. A bounded interval is left in both time signs.
interval=(.1,.4)
packet=(.375,.395)
for time_sign,spin_sign in ((1,-1),(-1,1)):
    # e^-itH f(w)=f(w+spin_sign*t).
    shift = -spin_sign*time_sign*.05
    image = tuple(x+shift for x in packet)
    assert image[0]>interval[1] and image[1]<1.

print(json.dumps({"minimal_gaussian_aux_dimension_per_n_strict_modes":"n",
                  "scalar_quarter_purified_maps":records,
                  "canonical_purification_disjoint_support_cross_covariance":locality,
                  "massless_bounded_interval_leakage_both_time_signs":1.,
                  "actual_Majorana_kernel_rank":int(np.linalg.matrix_rank(P_zero)),
                  "Higgs_out_of_Majorana_kernel_HS":float(np.linalg.norm(cross)),
                  "Higgs_out_of_Majorana_kernel_rank":int(np.linalg.matrix_rank(cross,tol=1e-10)),
                  "finite_Phi_at_0_55_kernel_rank":int(np.count_nonzero(abs(
                      np.linalg.eigvalsh(majorana+.55*higgs))<1e-12)),
                  "scope":"Constructive purification on the flat massless scalar spectral image with original 2pi x conversion. Opposite normal parity to tested Higgs return, not physical W or q mapping."},
                 indent=2))
