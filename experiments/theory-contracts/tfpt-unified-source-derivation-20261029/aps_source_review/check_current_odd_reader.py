"""Execute actual v252 -> two-face covariance -> original v258 -> odd reader.

Diagnostic integer couplings only; this is not a derivation of their values.
The original mu_geo and collar length remain separate symbols in the proof.
Numerical tests use ell=1/4, mu_geo=3/2 without calibrating their product.
"""
from pathlib import Path
import hashlib
import json
import sys
import numpy as np
import sympy as sy

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / 'verification'))
import v252_full_finite_triple as finite
import v258_dirac_covariance_induction as reader

ys = [np.diag(x).astype(complex) for x in
      ([1, 2, 3], [2, 3, 4], [3, 4, 5], [4, 5, 6])]
mr = np.array([[2, 1, 0], [1, 3, 1], [0, 1, 4]], complex)
phi0 = finite.full_D(*ys, mr)
dy = finite.full_D(*ys, np.zeros((3, 3)))
dmr = phi0 - dy
gamma = finite.grading()
K = finite.Kswap()
n = gamma.shape[0]
ident = np.eye(n)
Z = np.diag([1., -1.])
S = np.array([[0., 1.], [1., 0.]])
E = np.vstack(((ident + gamma)/2, 1j*(ident - gamma)/2))
U = (ident - 1j*gamma)/np.sqrt(2)

def exact(a):
    def cv(z):
        z = complex(z)
        if not z.real.is_integer() or not z.imag.is_integer():
            raise ValueError('not exact Gaussian integers')
        return sy.Integer(int(z.real)) + sy.I*sy.Integer(int(z.imag))
    return sy.SparseMatrix(sy.Matrix(a.tolist()).applyfunc(cv))

g, y, m, k = map(exact, (gamma, dy, dmr, K))
p, s, ell_s, mu_s = sy.symbols('p s ell mu', real=True)
phi = (1+s)*y + m
h = p*g - sy.I*g*phi
odd = (h-g*h*g)/2
unitary_numerator = sy.eye(n) - sy.I*g
checks = {
    'actual_96_coefficient': n == 96,
    'phi_odd_for_symbolic_s': g*phi+phi*g == sy.zeros(n),
    'odd_reader_removes_all_symbolic_momentum': odd == -sy.I*g*phi,
    'quarter_rotation_preserves_whole_phi':
        sy.expand(unitary_numerator*phi*unitary_numerator.H/2 + sy.I*g*phi) == sy.zeros(n),
    'charge_J_commutes_with_U_antilinearly':
        k*unitary_numerator.conjugate() == unitary_numerator*k,
    'reader_keeps_actual_yukawa_jet': odd.diff(s) == -sy.I*g*y,
    'reader_keeps_majorana': odd.subs(s,-1) == -sy.I*g*m,
}

ell, mu = 0.25, 1.5

def execute(sval, momentum):
    coeff = (1+sval)*dy + dmr
    B = np.kron(Z, momentum*ident) - np.kron(S, coeff)
    w, v = np.linalg.eigh(B)
    # Compute the actual full normal coefficient covariance before compression.
    a = (v * ((1+np.tanh(ell*w))/2)) @ v.conj().T
    x = (v * (1/(2*np.cosh(ell*w)))) @ v.conj().T
    G = np.block([[a,x],[x,np.eye(2*n)-a]])
    c22 = E.conj().T @ G[2*n:,2*n:] @ E
    ev = np.linalg.eigvalsh(c22)
    if not (ev.min()>1e-10 and ev.max()<1-1e-10):
        raise RuntimeError('v258 clipping would invalidate this diagnostic')
    induced = reader.induce(c22)
    extracted = mu*(induced-gamma@induced@gamma)/2
    target = 2*ell*mu*U@coeff@U.conj().T
    return extracted, {
        's':sval, 'p':momentum,
        'covariance_min':float(ev.min()),'covariance_max':float(ev.max()),
        'projector_error':float(np.linalg.norm(G@G-G)),
        'reader_error':float(np.linalg.norm(extracted-target)),
        'unprojected_reader_error':float(np.linalg.norm(induced-2*ell*(momentum*gamma-1j*gamma@coeff))),
    }

cases = []
for sval in (-0.2,0.,0.3):
    for momentum in (-0.75,0.25,0.75,1.25):
        _, row = execute(sval,momentum)
        cases.append(row)
step=1e-4
plus,_=execute(step,0.25)
minus,_=execute(-step,0.25)
jet_error=float(np.linalg.norm((plus-minus)/(2*step)-2*ell*mu*U@dy@U.conj().T))
checks['entire_actual_numeric_reader'] = max(r['reader_error'] for r in cases)<1e-10
checks['numeric_projector_before_reader'] = max(r['projector_error'] for r in cases)<1e-11
checks['actual_yukawa_jet_numeric'] = jet_error<1e-8
if not all(checks.values()):
    raise RuntimeError(checks)
result={
    'all_pass':True,'checks':checks,'cases':cases,'jet_error':jet_error,
    'ell_diagnostic':ell,'mu_geo_diagnostic':mu,'two_ell_mu':2*ell*mu,
    'formula':'mu_geo Pi_odd logodds(E* G22(B_s) E) = 2 ell mu_geo U Phi_s U*',
    'scope':'Exact symbolic coefficient proof and numerical actual-constructor/reader test. No P_F source selection, original mass or scale derivation, or physical time identification.',
    'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
       for p in (ROOT/'verification/v252_full_finite_triple.py',ROOT/'verification/v258_dirac_covariance_induction.py')},
}
out=Path(__file__).with_name('CURRENT-ODD-READER-CHECK.json')
out.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'all_pass':True,'checks':checks,'max_reader_error':max(r['reader_error'] for r in cases),'jet_error':jet_error,'output':str(out)}))
