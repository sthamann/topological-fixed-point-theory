"""Unwrap the original q/p word and test its native quarter cross response."""
from pathlib import Path
import sys
import sympy as s
import numpy as np
from scipy.linalg import solve

R=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
# Affine actions on the unwrapped argument, in units of pi.
def compose(f,g):
    return (f[0]*g[0],f[0]*g[1]+f[1])
def inv(f):
    return (f[0],-f[0]*f[1])
def word(*fs):
    out=(1,0)
    for f in fs:out=compose(out,f)
    return out
Q=(1,s.Rational(1,2)); S=(-1,0)
comm=word(Q,S,inv(Q),inv(S))
assert comm==(1,1)
assert word(inv(comm),comm)==(1,0) and word(comm,comm)==(1,2)
# The source-compatible sqrt(dz) frame carries a geometric spin sign.
assert s.I*(-s.I)==1 and s.I*s.I==-1
print('Unwrapped original q/p commutator: +pi; the two closures have winding 0 and 1, with geometric spin +1 and -1.')

# At x=2theta in the original NS/quarter dictionary, one z base turn
# shifts theta by pi. These are the known W0 and W1 multiplier exponents.
beta=[s.Rational(1,2),s.Rational(3,2)]
v0,v1=[s.simplify(s.exp(s.I*s.pi*b)) for b in beta]
assert (v0,v1)==(s.I,-s.I) and s.conjugate(v1)*v0==-1
print('Native longitudinal return on the double cover: V=T_pi, eigenvalues +i,-i. Opposite-quarter cross carries relative sign -1.')

# Re-run the same actual marked fixed history used in the original finite
# R032 quarter calculation. No new phase is supplied: V uses the modes.
f=R/'experiments/theory-contracts/chat-triangulation-20261108/check_native_quarter_winding_and_geometry_20261111.py'
previous=sys.argv
sys.argv=[str(f),str(R)]
ns={'__file__':str(f)}
exec(compile(f.read_text().split('baseline, _ = relative_cross()',1)[0],str(f),'exec'),ns)
sys.argv=previous
receipt,known=ns['relative_cross']()
a,b=[ns['arm'](x,False)[0] for x in (.25,.75)]
Id=ns['identity']; mark=ns['raising']
def test(A,B,M):
    ga=Id+A.conj().T@A; gb=Id+B.conj().T@B
    x=solve(gb,B.conj().T@M@A,assume_a='pos')
    return solve(ga,x.conj().T,assume_a='pos').conj().T
results=[]
for cuts in [(),(2,),(4,),(2,4)]:
    if cuts:
        A,B=ns['sew'](a,cuts),ns['sew'](b,cuts)
        MM=ns['frames'][-1].conj().T@mark@ns['frames'][-1]
    else:
        A,B=ns['direct'](a),ns['direct'](b)
        MM=mark
    x=test(A,B,MM)
    y=test(complex(v0)*A,complex(v1)*B,MM)
    co=test(complex(v0)*A,complex(v1)*B,-MM)
    assert np.array_equal(A.conj().T@A,(complex(v0)*A).conj().T@(complex(v0)*A))
    assert np.array_equal(B.conj().T@B,(complex(v1)*B).conj().T@(complex(v1)*B))
    errs=[ns['maxnorm'](y+x),ns['maxnorm'](co-x)]
    assert max(errs)<1e-12
    results.append((cuts,errs))
assert ns['maxnorm'](test(ns['direct'](a),ns['direct'](b),mark)[np.ix_(ns['u'],ns['d'])]-known)<1e-13
print('Original R032 direct and two-cut controls, errors for fixed/cotransported mark:',results)
print('Nonzero original open family block singular values:',receipt['uLdL_family_singular_values'])
print('Relative empty log weight unchanged:',receipt['relative_log_empty'])
print('With fixed weak marker the open block reverses sign. Transporting the marker too cancels the sign as a frame change.')

# Both blocks carry the same endpoint charge. Their mixed Hilbert-Schmidt
# contraction is invariant under the same unitary gauge frames on both arms.
ix=np.ix_(ns['u'],ns['d'])
A,B=ns['direct'](a),ns['direct'](b)
zero=test(A,B,mark)[ix]
one=test(complex(v0)*A,complex(v1)*B,mark)[ix]
norm=float(np.vdot(zero,zero).real)
interference=np.vdot(zero,one)
Gmix=np.array([[np.vdot(U,V) for V in (zero,one)] for U in (zero,one)])
assert interference==complex(-norm,0) and norm>0
assert np.max(np.abs(Gmix-norm*np.array([[1,-1],[-1,1]])))<1e-14
print('Neutral mixed contraction at fixed marker:',interference,'and self-norm:',norm)
print('Positive two-history Gram spectrum:',np.linalg.eigvalsh(Gmix).tolist())
