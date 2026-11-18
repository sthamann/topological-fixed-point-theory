"""Exact actual v252 coefficient control; no physical mass selection."""
from pathlib import Path
import hashlib
import importlib.util
import json
import sys
import numpy as np
import sympy as s

ROOT=Path(__file__).resolve().parents[4]
SOURCE=ROOT/'verification/v252_full_finite_triple.py'
sys.path.insert(0,str(ROOT/'verification'))
spec=importlib.util.spec_from_file_location('actual_v252',SOURCE)
v=importlib.util.module_from_spec(spec)
spec.loader.exec_module(v)

def exact(arr):
    def entry(z):
        z=complex(z)
        if not(z.real.is_integer() and z.imag.is_integer()):
            raise ValueError('Constructor input is not a Gaussian integer')
        return s.Integer(int(z.real))+s.I*s.Integer(int(z.imag))
    return s.SparseMatrix(s.Matrix(arr.tolist()).applyfunc(entry))

ys=[np.diag(x).astype(complex) for x in
    ([1,2,3],[2,3,4],[3,4,5],[4,5,6])]
mr=np.array([[2,1,0],[1,3,1],[0,1,4]],complex)
phi=exact(v.full_D(*ys,mr))
g=exact(v.grading())
n=g.rows
one=s.eye(n)
Z=s.diag(1,-1)
N=s.Matrix([[0,-1],[1,0]])
S=N*Z
Gt=s.kronecker_product(Z,g)
Ic=s.I*Gt
C=s.kronecker_product(s.I*N,one)
plus,minus=(one+g)/2,(one-g)/2
E=plus.col_join(s.I*minus)
p=s.symbols('p',real=True)
B=s.kronecker_product(Z,p*one)-s.kronecker_product(S,phi)
h=p*g-s.I*g*phi
checks={
 'actual_dimension':n==96,
 'actual_phi_selfadjoint':phi.H==phi,
 'actual_phi_odd':phi*g+g*phi==s.zeros(n),
 'actual_total_B_even':Gt*B==B*Gt,
 'C_B_anticommutes':C*B==-B*C,
 'new_C_antilinear':C*Ic==-Ic*C,
 'E_isometry':E.H*E==one,
 'E_positive_total_chirality':Gt*E==E,
 'whole_original_mass_operator_intertwiner':B*E==E*h,
 'selfadjoint_h':h.H==h,
 'squared_h':h*h==p*p*one+phi*phi,
}
# The same original coefficient with arbitrary smooth scalar amplitude and
# connection, represented by independent value/first-derivative jets.
f,fp,a,ap=s.symbols('f fp a ap',real=True)
h1=-s.I*g
h0=a*g-s.I*g*f*phi
h0prime=ap*g-s.I*g*fp*phi
checks.update({
 'profile_second_derivative_coefficient':h1*h1==-one,
 'profile_first_derivative_coefficient':h1*h0+h0*h1==-2*s.I*a*one,
 'profile_zeroth_coefficient_retains_derivative':
     s.simplify(h1*h0prime+h0*h0-((a*a-s.I*ap)*one+f*f*phi*phi-fp*phi))==s.zeros(n),
})
if not all(checks.values()):
    raise RuntimeError(checks)
result={'all_pass':True,'checks':checks,'constructor_sha256':
    hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
    'scope':'Full actual finite coefficient and smooth classical profile jets; not physical quantization or bosonic state selection.'}
Path(__file__).with_name('finite-results.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result))
