"""Exact common Spin4/R032 coefficient map, checked on the unchanged histories."""
import contextlib
import io
import json
from pathlib import Path
import runpy

import numpy as np
from scipy.linalg import eigh

base = Path(__file__).resolve().parent
with contextlib.redirect_stdout(io.StringIO()):
    spin = runpy.run_path(str(base / 'check_aps_higgs_normal_jet_20261110.py'))
    source = runpy.run_path(str(base / 'check_r113_mixed_metric_matter_20261110.py'))

gam = spin['gamma']
gamma5 = spin['gamma5']
I = source['I']
n = len(I)
grade = source['grade']
g = np.diag(grade).real
p = np.array([.15,.20,0.])  # Nonaxis direction at the actual original |p|=1/4.
abs_p = np.linalg.norm(p)
gn = gam[3]
normal = -1j * gn
clifford_mass = 1j * gn @ gamma5
A = 1j * gn @ sum(pj * gj for pj, gj in zip(p, gam)) / abs_p
Gamma = np.kron(gamma5, grade)
N = np.kron(normal, I)
free = abs_p * np.kron(A, I)


def maps():
    out = []
    for chirality in (1,-1):
        for sector in (1,-1):
            P = (np.eye(4)+chirality*gamma5)@(np.eye(4)+sector*A)/4
            assert abs(np.trace(P)-1) < 1e-13
            u = P[:, int(np.argmax(np.diag(P).real))]
            u /= np.linalg.norm(u)
            v = -1j*clifford_mass@u
            W = np.vstack([np.diag(np.where(g>0,u[k],v[k])) for k in range(4)])
            assert np.max(abs(Gamma@W-chirality*W)) < 1e-13
            assert np.max(abs(free@W-W@(sector*abs_p*grade))) < 1e-13
            out.append((chirality,sector,W))
    return out


maps_by_sector = maps()
W_full = np.hstack([W for _,_,W in maps_by_sector])
assert np.max(abs(W_full.conj().T@W_full-np.eye(4*n))) < 1e-13
K = source['original']['source']['native'].Kswap()
J_linear = np.kron(gam[0]@gam[2], K)
C_linear = Gamma @ J_linear
assert np.max(abs(J_linear@J_linear.conj()+np.eye(4*n))) < 1e-13
assert np.max(abs(C_linear@C_linear.conj()-np.eye(4*n))) < 1e-13
assert np.max(abs(C_linear@N.conj()+N@C_linear)) < 1e-13


def spectral(B, t):
    e,U = eigh(B)
    return (U*np.exp(-t*e))@U.conj().T, (U*(e>0))@U.conj().T


R4 = np.eye(4*n,dtype=complex)
R96 = [I.astype(complex) for _ in maps_by_sector]
errors = dict(operator=0., moving_projector=0., fixed_projector=0., real_pair=0.)
caps = []
for index, (field, width) in enumerate(zip(source['history'], source['durations'])):
    phi = source['inner_form'](field)+source['M']
    B = free + np.kron(clifford_mass, phi)
    B_mirror = -free + np.kron(clifford_mass, phi)
    errors['real_pair'] = max(errors['real_pair'],np.max(abs(C_linear@B.conj()-B_mirror@C_linear)))
    evolution, moving = spectral(B, width)
    R4 = evolution @ R4
    for k, (_,sector,W) in enumerate(maps_by_sector):
        h = sector*abs_p*grade-1j*grade@phi
        small_evolution, small_moving = spectral(h,width)
        R96[k] = small_evolution@R96[k]
        errors['operator'] = max(errors['operator'],np.max(abs(B@W-W@h)))
        errors['moving_projector'] = max(errors['moving_projector'],np.max(abs(moving@W-W@small_moving)))
        errors['fixed_projector'] = max(errors['fixed_projector'],np.max(abs((free@W/abs_p)-W@(sector*grade))))
    if index in (2,5):
        sign_plus, l_plus = np.linalg.slogdet(R96[0][np.ix_(g<0,g<0)])
        sign_minus, l_minus = np.linalg.slogdet(R96[1][np.ix_(g>0,g>0)])
        assert abs(sign_minus-sign_plus.conj()) < 1e-12 and abs(l_plus-l_minus) < 1e-12
        caps.append({'steps':index+1,'individual_phase':float(np.angle(sign_plus)),
                     'paired_S_relative_cap_log':float(l_plus+l_minus-n*abs_p*sum(source['durations'][:index+1]))})

errors['full_ordered_history'] = max(np.max(abs(R4@W-W@R))
    for R,(_,_,W) in zip(R96,maps_by_sector))
assert max(errors.values()) < 5e-12
print(json.dumps({'momentum':p.tolist(),'sectors':[list(x[:2]) for x in maps_by_sector],
    'sector_dimension':n,'errors':{k:float(v) for k,v in errors.items()},'fixed_caps':caps}, indent=2))
