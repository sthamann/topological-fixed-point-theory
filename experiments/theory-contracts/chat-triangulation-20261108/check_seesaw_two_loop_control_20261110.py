"""Two-loop SM+N continuation control using PyR@TE-generated dimensionless beta functions.

SM and sterile-neutrino matter content is continued above the TFPT PS/GUT
thresholds only for this diagnostic. Heavy thresholds and the v3 Dirac
entries at decoupling are supplied. Keep native v481's top-only quark input.
"""
import contextlib
import importlib.util
import io
import json
import math
import runpy
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

repo = Path(sys.argv[1]).resolve()
beta_path = Path(sys.argv[2]).resolve()
arguments = sys.argv
try:
    sys.argv = [str(Path(__file__).with_name("check_seesaw_lepton_ratios_20261110.py")), str(repo)]
    with contextlib.redirect_stdout(io.StringIO()):
        w = runpy.run_path(sys.argv[0])
finally:
    sys.argv = arguments

spec = importlib.util.spec_from_file_location("sm_type_i_rges", beta_path)
rge = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rge)

sm = w["sm"]
loop = w["loop"]
m2, m3 = w["m2"], w["m3"]
gauge0 = np.sqrt(4*math.pi/np.array(sm.A_INV_MZ))
ym = [np.zeros((3,3),complex) for _ in range(4)]
ym[0][-1,-1] = sm.YT_MZ
ym[2] = np.diag(w["msbar"]/(1000*sm.V_EW/math.sqrt(2)))
v2 = w["y2_at_m2"] * w["U"][:,1]
v3 = w["y3"] * w["U"][:,2]
t = lambda mu: math.log(mu/sm.MZ)


def flatten(gs, lam, mats):
    return np.r_[gs, lam, *[Y.flatten() for Y in mats]]


def data(y):
    return y[:3].real, y[3].real, y[4:].reshape((4,3,3))


def rhs(t_, y, loops):
    g, lam, mat = data(y)
    g1,g2,g3 = g
    u,d,e,n = map(np.asmatrix, mat)
    dg = np.zeros(3,complex)
    dlam = 0.
    db = np.zeros_like(mat)
    for j in range(1,loops+1):
        weight = loop**j
        dg += weight*np.array([rge.beta_g1(j,g1,g2,g3,u,d,e,n),
                               rge.beta_g2(j,g2,g1,g3,u,d,e,n),
                               rge.beta_g3(j,g3,g1,g2,u,d)])
        dlam += weight*rge.beta_lambda_(j,g1,g2,u,d,e,n,lam,g3)
        db += weight*np.array([rge.beta_Yu(j,g1,g2,g3,u,d,e,n,lam),
                               rge.beta_Yd(j,g1,g2,g3,u,d,e,n,lam),
                               rge.beta_Ye(j,g1,g2,u,d,e,n,g3,lam),
                               rge.beta_Ynu(j,g1,g2,u,d,e,n,g3,lam)])
    assert max(np.abs(dg.imag)) < 1e-12 and abs(dlam.imag) < 1e-12
    return flatten(dg.real,dlam.real,db)


def evolve(a,b,y,loops):
    sol = solve_ivp(lambda s,z: rhs(s,z,loops),(t(a),t(b)),y,rtol=4e-10,atol=3e-13)
    assert sol.success
    return sol.y[:,-1]


def match(y,i,vector,down):
    g, lam, mat = data(y)
    mat = mat.copy()
    if down:
        v = mat[3,:,i].copy()
    else:
        v = vector
    gram = np.outer(v,v.conj())
    scale = 1-loop*np.vdot(v,v).real/4
    charged = scale*np.eye(3)+(5/8)*loop*gram
    if down:
        mat[0:2] *= scale
        mat[2] = charged @ mat[2]
        mat[3,:,i] = 0.
    else:
        mat[0:2] /= scale
        mat[2] = np.linalg.solve(charged,mat[2])
        mat[3,:,i] = v
    return flatten(g,lam,mat)


def transfer(mu,loops,neutrinos,retain_uv_singular_vectors=False):
    y = flatten(gauge0,sm.LAM_MZ,ym)
    thresholds = [(m2,1,v2),(m3,2,v3)] if neutrinos else []
    thresholds = [item for item in thresholds if item[0] < mu]
    old = sm.MZ
    for edge,i,vec in thresholds:
        y = match(evolve(old,edge,y,loops),i,vec,False)
        old = edge
    y = evolve(old,mu,y,loops)
    g,lam,mat = data(y)
    mat = mat.copy()
    left, actual_UV_singular_values, right = np.linalg.svd(mat[2])
    source = w["source_masses"]*actual_UV_singular_values[0]/w["source_masses"][-1]
    mat[2] = (left @ np.diag(source[::-1]) @ right if retain_uv_singular_vectors
              else np.diag(source))
    y = flatten(g,lam,mat)
    old = mu
    for edge,i,vec in thresholds[::-1]:
        y = match(evolve(old,edge,y,loops),i,vec,True)
        old = edge
    down = evolve(old,sm.MZ,y,loops)
    me = np.sort(np.linalg.svd(data(down)[2][2],compute_uv=False))
    return (w["shape"](me)-w["shape"](source)).tolist(), w["observable"](me)


result = {}
for label,mu in (("scalaron_v163",w["m_scal"]),("GUT_v481",sm.M_GUT),("Mbar_attachment",w["mbar"])):
    result[label] = {}
    for j in (1,2):
        delta, shape = transfer(mu,j,True)
        base, _ = transfer(mu,j,False)
        result[label][j] = {
            "source_to_MZ_delta_log_ratios": delta,
            "SM_only_delta_log_ratios": base,
            "heavy_N_increment": (np.array(delta)-np.array(base)).tolist(),
            "predicted_K": shape["K"],
            "predicted_Koide": shape["Koide"],
        }
assert result["Mbar_attachment"][2]["heavy_N_increment"][2] < -.005
assert abs(result["Mbar_attachment"][2]["heavy_N_increment"][2]-result["Mbar_attachment"][1]["heavy_N_increment"][2]) < .001
retained_vectors, _ = transfer(w["mbar"],2,True,retain_uv_singular_vectors=True)
assert abs(retained_vectors[2]-result["Mbar_attachment"][2]["source_to_MZ_delta_log_ratios"][2]) < .001
print(json.dumps({
    "scope": "full two-loop dimensionless SM+N beta functions of explicitly specified continuation; finite one-loop dimension-four charged-fermion Yukawa matching; fixed M_i, no scalar/PS/GUT matching; not the TFPT UV matching calculation",
    "beta_origin": "locally installed PyR@TE 3, native models/SM.model plus 3 singlet nuR with Ynu Lbar Htilde nuR",
    "high_scale_charged_source_basis": "fixed flavor basis with the diagonal source spectrum; optional retain-upward-singular-vectors control",
    "Mbar_two_loop_with_upward_UV_singular_vectors_delta_log_ratios": retained_vectors,
    "outputs": result,
},indent=2))
