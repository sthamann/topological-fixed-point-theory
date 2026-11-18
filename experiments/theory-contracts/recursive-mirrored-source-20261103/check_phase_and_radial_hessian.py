"""Original whole96 Higgs phase and radial spectral identities.

Uses the unchanged literal rng252 fixture. Its values and heat profiles
are diagnostic matrix probes, not selected physical parameters.
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import eigh, expm

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / 'verification'))
import v252_full_finite_triple as original

rng = np.random.default_rng(252)
def original_yukawa_fixture():
    return rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))
Yu, Yd, Ye, Yn = [original_yukawa_fixture() for _ in range(4)]
M0 = rng.normal(size=(3, 3))
M = M0 + M0.T
zero = np.zeros((3, 3))
DY = original.full_D(Yu, Yd, Ye, Yn, zero)
DM = original.full_D(zero, zero, zero, zero, M)

particle_T = np.array([
    (1 if name in ('nu', 'u') else -1) if chirality == 'L' else 0
    for name, chirality, color, generation in original.P
])
T = np.diag(np.r_[particle_T, -particle_T])
phase_derivative = 1j * (T @ DY - DY @ T)
charged_norm = np.trace(
    3 * Yu.conj().T @ Yu + 3 * Yd.conj().T @ Yd + Ye.conj().T @ Ye
).real
whole_norm = charged_norm + np.trace(Yn.conj().T @ Yn).real
K = original.Kswap()
phase = {
    'dimension': int(original.NH),
    'Majorana_commutator_norm': float(np.linalg.norm(T @ DM - DM @ T)),
    'phase_derivative_Hermiticity': float(np.linalg.norm(phase_derivative - phase_derivative.conj().T)),
    'phase_derivative_reality': float(np.linalg.norm(K @ np.conj(phase_derivative) @ K - phase_derivative)),
    'Tr_phase_derivative_squared': float(np.trace(phase_derivative @ phase_derivative).real),
    'Tr_DY_squared': float(np.trace(DY @ DY).real),
    'expected_4A': float(4 * whole_norm),
    'trace_residual': float(abs(np.trace(phase_derivative @ phase_derivative).real - 4 * whole_norm)),
}
for key in ('Majorana_commutator_norm', 'phase_derivative_Hermiticity', 'phase_derivative_reality', 'trace_residual'):
    assert phase[key] == 0

eigenvalues, U = eigh(DM)
E = U.conj().T @ DY @ U
radial = []
for heat_time in (.1, .4, 1.):
    fp = -2 * heat_time * eigenvalues * np.exp(-heat_time * eigenvalues ** 2)
    delta = eigenvalues[:, None] - eigenvalues[None, :]
    divided = np.empty_like(delta)
    distinct = abs(delta) > 1e-9
    np.divide(fp[:, None] - fp[None, :], delta, out=divided, where=distinct)
    fpp = (4 * heat_time ** 2 * eigenvalues ** 2 - 2 * heat_time) * np.exp(-heat_time * eigenvalues ** 2)
    divided[~distinct] = np.broadcast_to(fpp[:, None], delta.shape)[~distinct]
    whole = np.sum(divided * abs(E) ** 2).real
    closed = -8 * heat_time * (
        charged_norm + np.trace(expm(-heat_time * M @ M.conj().T) @ Yn.conj().T @ Yn).real
    )
    error = abs(whole - closed) / max(1, abs(whole))
    assert error < 1e-12
    radial.append({
        'diagnostic_heat_time': heat_time,
        'full96_radial_Hessian': float(whole),
        'closed_profile_formula': float(closed),
        'relative_error': float(error),
    })

# Full original Majorana compatible phase response, including its contact.
Phi = DY + DM
A = Phi @ Phi
B = 1j * (T @ Phi - Phi @ T)
L = 1j * (T @ A - A @ T)
Mcontact = -(T @ (T @ A - A @ T) - (T @ A - A @ T) @ T)
evals, U = eigh(A)
Le = U.conj().T @ L @ U
Be = U.conj().T @ B @ U
nodes, weights = np.polynomial.legendre.leggauss(96)
nodes, weights = (nodes + 1) / 2, weights / 2
phase_kernel = []
heat_time = .4
TL = []
TB = []
for a in nodes:
    w = np.exp(-heat_time * ((1-a) * evals[:, None] + a * evals[None, :]))
    TL.append(np.sum(abs(Le)**2 * w).real)
    TB.append(np.sum(abs(Be)**2 * w).real)
TL, TB = np.array(TL), np.array(TB)
contact = np.trace(Mcontact @ expm(-heat_time * A)).real
contact_residual = abs(contact - heat_time * np.dot(weights, TL))
assert contact_residual < 1e-10
for p_squared in (0., .1, 1., 10.):
    q = np.exp(-heat_time * nodes * (1-nodes) * p_squared)
    combined = np.dot(weights, (q-1)*TL + p_squared*q*TB) / (8*np.pi**2)
    uncombined = (np.dot(weights, q*(TL+p_squared*TB)) - contact/heat_time) / (8*np.pi**2)
    assert abs(combined-uncombined) < 1e-11
    phase_kernel.append({'diagnostic_p_squared': p_squared,
                         'whole96_phase_theta_squared_coefficient': float(combined),
                         'contact_cancellation_error': float(abs(combined-uncombined))})

# Complex symmetric Majorana check distinguishes MM* from the wrong M*M.
complex_majorana = []
for seed in (11, 42):
    r = np.random.default_rng(seed)
    mc = r.normal(size=(3,3)) + 1j*r.normal(size=(3,3))
    mc = mc + mc.T
    dc = original.full_D(zero,zero,zero,zero,mc)
    vals, uc = eigh(dc)
    ec = uc.conj().T @ DY @ uc
    fp = -2*heat_time*vals*np.exp(-heat_time*vals**2)
    dv = vals[:,None]-vals[None,:]
    dd = np.zeros_like(dv)
    mask = abs(dv)>1e-9
    np.divide(fp[:,None]-fp[None,:],dv,out=dd,where=mask)
    fpp = (4*heat_time**2*vals**2-2*heat_time)*np.exp(-heat_time*vals**2)
    dd[~mask]=np.broadcast_to(fpp[:,None],dv.shape)[~mask]
    direct = np.sum(dd*abs(ec)**2).real
    closed = -8*heat_time*(charged_norm+np.trace(expm(-heat_time*mc@mc.conj().T)@Yn.conj().T@Yn).real)
    wrong = -8*heat_time*(charged_norm+np.trace(expm(-heat_time*mc.conj().T@mc)@Yn.conj().T@Yn).real)
    assert abs(direct-closed)<1e-10
    assert abs(direct-wrong)>.1
    complex_majorana.append({'diagnostic_seed':seed,'direct_full96':float(direct),
                              'correct_MM_adjoint':float(closed),'wrong_adjointM_M':float(wrong)})

result = {
    'scope': 'Unchanged original finite96 fixture; exact Higgs phase identity and full-profile radial Hessian check only. No primitive Sadm identity or physical parameter/profile selection.',
    'phase': phase,
    'radial': radial,
    'whole96_phase_kernel': phase_kernel,
    'Majorana_contact_identity_error': float(contact_residual),
    'complex_Majorana_orientation': complex_majorana,
    'primitive_source_identification_proved': False,
}
print(json.dumps(result, indent=2))
if '--save' in sys.argv:
    (Path(__file__).parent / 'phase-and-radial-hessian-results.json').write_text(json.dumps(result, indent=2) + '\n')
