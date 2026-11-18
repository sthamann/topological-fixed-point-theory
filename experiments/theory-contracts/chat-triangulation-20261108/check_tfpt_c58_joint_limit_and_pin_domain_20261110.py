"""C58 Cartan crossover and the uncompressed original Pin domain."""

import contextlib
import io
import json
from pathlib import Path
import runpy
import sys

import numpy as np
from scipy.integrate import quad


repo = Path(sys.argv[1]).resolve()
contract = repo / "experiments/theory-contracts/chat-triangulation-20261108"
with contextlib.redirect_stdout(io.StringIO()):
    old = runpy.run_path(str(contract / "check_tfpt_c58_short_time_mixed_cap_20261110.py"))

strength = old["coefficient"]
crossover = {}
for split in (.5, 1/3):
    def integrand(u):
        return np.expm1(-2*split*u)*np.expm1(-2*(1-split)*u)/u**2

    rows = {}
    for y in (.1, 1., 4.):
        prediction = strength / 2 * quad(integrand, 0, y, epsabs=1e-12)[0]
        raw = {}
        for n in (64, 256, 1024):
            width = y/n
            raw[str(n)] = old["cap_read"](width,split,n,False)["mixed_log_jet_over_ell"]
        assert abs(raw["1024"]-prediction) < 2e-4
        rows[str(y)] = {"prediction":prediction,"observed":raw}
    crossover[str(split)] = rows

# In the same cap, T3+Y commutes exactly with BOTH Higgs/Majorana pieces.
# This checks the angle response in U and in X=U exp(-i phi T3).
angles = {}
for sign, data in old["spectral"].items():
    ur, sr, vhr, ul, sl, vhl, dt, dy, _ = data
    p = .25 if sign == 1 else .75
    def leg(s):
        e = np.hypot(p,s)
        a = np.tanh(.125/2*e)
        return s*a/(e+p*a)
    right = -(ur*leg(sr)[None,:])@vhr
    left = -(vhl.conj().T*leg(sl)[None,:])@ul.conj().T
    assert np.max(abs((dt+dy)*left)) < 1e-14
    mat = np.eye(len(sr))+right@left

    def log_hessian(first, second):
        d1 = np.linalg.solve(mat,1j*right@(first*left))
        d2 = np.linalg.solve(mat,1j*right@(second*left))
        d12 = np.linalg.solve(mat,-right@(first*second*left))
        return np.trace(d12-d1@d2).real

    angles[str(sign)] = {"raw_T3_Y":float(log_hessian(dt,dy)),
                         "raw_T3_T3":float(log_hessian(dt,dt)),
                         "same_cap_rotated_X_Y":float(log_hessian(dt,dt+dy))}
    assert abs(log_hessian(dt,dt+dy)) < 1e-14

pin = runpy.run_path(str(repo / ("experiments/theory-contracts/"
    "original-admissible-source-closure-20261025/geometric_source_return/source_geometric_machine.py")))
replay = runpy.run_path(str(contract / "test_original_two_face_transition_20261108.py"))
cut = {}
for epsilon in (0., .1):
    k = np.pi / 2
    modes, d, _, _ = pin["dirac_block"](64, epsilon, k)
    initial = np.flatnonzero(modes == 0)[0]
    others = modes != 0
    eig, u = np.linalg.eigh(d)
    z = 2j
    direct = np.sum(abs(u[initial,:])**2/(eig-z))
    compressed = 1/(d[initial,initial]-z)
    marked = replay["marked_green"](replay["transfer"](z,epsilon,k,[0.,.37,.71,1.],marked=True))
    cut[str(epsilon)] = {
        "norm_Q_D_P":float(np.linalg.norm(d[others,initial])),
        "true_resolvent": [float(direct.real),float(direct.imag)],
        "compressed_resolvent": [float(compressed.real),float(compressed.imag)],
        "compression_defect":float(abs(direct-compressed)),
        "marked_strip_error":float(abs(marked-direct)),
    }
    assert abs(marked-direct) < 3e-9
assert cut["0.0"]["compression_defect"] < 1e-12
assert cut["0.1"]["compression_defect"] > 1e-5

print(json.dumps({"crossover":crossover,"cartan_same_cap":angles,
                  "original_pin_normal_zero_branch":cut,
                  "scope":"C58 fixed external Cartan link and original geometry-only Pin cut; no full R113 or P_adm mapping."}, indent=2))
