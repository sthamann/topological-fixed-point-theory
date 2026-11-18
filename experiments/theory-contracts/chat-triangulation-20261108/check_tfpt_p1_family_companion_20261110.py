"""Test an exterior family-root lift on the actual native paired Pin coefficients."""
from collections import defaultdict
import contextlib
from fractions import Fraction
import io
import json
from pathlib import Path
import runpy
import sys

import numpy as np
import sympy as sp

repo=Path(sys.argv[1]).resolve()
base=repo/'experiments/theory-contracts/chat-triangulation-20261108'
with contextlib.redirect_stdout(io.StringIO()):
    source=runpy.run_path(str(base/'check_tfpt_conditioned_r113_weights_20261110.py'))
rep=source['native']
tau2=np.array([[0,-1j],[1j,0]])
jpin=1j*tau2
H=.55
phi=source['majorana']+source['inner'](np.array([0,H]))
K=np.kron(-tau2,phi)
B=[np.kron(-tau2,source['inner'](np.array(h)))
   for h in ((0,1),(0,1j),(1,0),(1j,0))]
energies,V=np.linalg.eigh(K)
gaps=np.abs(energies[:,None]-energies)
same_energy=gaps<1e-10
distinct=float(np.min(gaps[~same_energy]))
assert distinct>1e-7
grading=np.kron(np.eye(2),source['gamma'])
charge_conj=np.kron(np.eye(2),source['swap'])
weak=[np.kron(np.eye(2),source['weak_lift'](g))
      for g in (np.diag([1j,-1j]),np.array([[0,1],[-1,0]]))]
color=np.kron(np.eye(2),np.diag([np.exp(2j*np.pi/3) if f in ('u','d') else 1
                                    for f,_,_,_ in rep.P]+[np.exp(-2j*np.pi/3) if f in ('u','d') else 1
                                    for f,_,_,_ in rep.P]))

def comm(a,b):
    return a@b-b@a

def hs(a):
    return float(np.linalg.norm(a))

def exterior_apply(matrix, state):
    out=defaultdict(complex)
    for indices,value in state.items():
        for slot,col in enumerate(indices):
            for row in np.flatnonzero(matrix[:,col]):
                target=list(indices)
                target[slot]=int(row)
                if len(set(target))<len(target):
                    continue
                sign=(-1)**sum(target[i]>target[j] for i in range(len(target)) for j in range(i+1,len(target)))
                out[tuple(sorted(target))]+=sign*matrix[row,col]*value
    return {k:v for k,v in out.items() if v}

def weight(state):
    return sum(abs(v)**2 for v in state.values())

roots={}
for i in range(3):
    for j in range(3):
        if i==j:
            continue
        e=np.zeros((rep.NP,rep.NP))
        for f,ch,c,g in rep.P:
            if g==j:
                e[rep.IDX[(f,ch,c,i)],rep.IDX[(f,ch,c,j)]]=1
        zero=np.zeros_like(e)
        # Complexified root of the real family action: particles E_ij,
        # antiparticles -E_ji. J_F exchanges the root with minus its adjoint.
        root=np.block([[e,zero],[zero,-e.T]])
        N=np.kron(jpin,root)
        A=N.conj().T
        assert hs(N@N)==0
        assert hs(charge_conj@N.conj()@charge_conj-A)==0
        assert hs(comm(N,grading))==0
        assert all(hs(comm(N,g))<1e-13 for g in (*weak,color))
        C=comm(N,A)
        assert hs(comm(C,N)-2*N)==0
        # Right-handed electron and its antiparticle: color, weak and
        # hypercharge cancel. Independent orthogonal modes carry two doublets.
        v0=tuple(sorted((rep.IDX[('e','R',0,j)],rep.NP+rep.IDX[('e','R',0,i)])))
        v={v0:1}
        w=exterior_apply(N,v)
        z=exterior_apply(N,w)
        assert [weight(x) for x in (v,w,z)]==[1,2,4]
        assert not exterior_apply(N,z)
        assert exterior_apply(A,w)=={v0:2}
        assert exterior_apply(A,z)=={k:2*c for k,c in w.items()}
        # Time average over this exact native Hermitian comparison K.
        # Distinct levels are numerically separated by >1e-7.
        stationary=V@((V.conj().T@N@V)*same_energy)@V.conj().T
        assert hs(comm(K,stationary))<1e-11
        roots[f'{i+1}<-{j+1}']={
           'exterior_0_1_2_power_norm2':[1,2,4],
           'pin_label_count_mod2_along_triplet':[[sum(k>=96 for k in t)%2 for t in x] for x in (v,w,z)],
           'actual_time_commutator_HS':hs(comm(K,N)),
           'higgs_commutator_HS':[hs(comm(b,N)) for b in B],
           'stationary_part_HS':hs(stationary),
           'stationary_part_normality_residual_HS':hs(comm(stationary,stationary.conj().T)),
           'stationary_part_square_HS':hs(stationary@stationary),
        }

# In each charged species, fixed color, particle/dual and tau2 block,
# K has a simple spectrum if the 3x3 Yukawa squared singular values
# are positive and distinct. This is an exact check on the binary inputs.
charged={}
phi_h=source['inner'](np.array([0.,1.]))
for f in ('e','u','d'):
    l=[rep.IDX[(f,'L',0,g)] for g in range(3)]
    r=[rep.IDX[(f,'R',0,g)] for g in range(3)]
    block=phi_h[np.ix_(l,r)]
    z=sp.Matrix([[sp.Rational(float(v.real))+sp.I*sp.Rational(float(v.imag))
                  for v in row] for row in block])
    x=sp.Poly((z*z.H).charpoly().as_expr(),extension=sp.I)
    assert x.degree()==3 and x.eval(0)!=0 and sp.gcd(x,x.diff()).degree()==0
    charged[f]={'exact_nonzero_simple_squared_singular_values':True,
                'exact_normality_of_stationary_one_body_root_on_charged_block':True}

old=json.loads((base/'tfpt-p1-local-spin-path-results-20261110.json').read_text())['old_p1_exact']
closed={}
for q,d in old.items():
    n=Fraction(d['reference_norm'])
    trip=Fraction(d['family_triplet_weight_zero_norm2'])
    w=Fraction(d['middle_W_norm2'])
    assert w==2*trip
    assert d['modes']['0']['early_vector_zero']
    closed[q]={'early_unweighted_closed_norm_over_q_reference':'0',
               'late_unweighted_closed_norm_over_q_reference':str(2*trip/n),
               'late_as_triplet_retention_of_charged_Bv':str(trip/Fraction(d['B_norm2']))}

print(json.dumps({
    'scope':'tested compound root, original native 192 coefficient chart on normal zero branch and H=.55; not the physical field representation',
    'minimum_resolved_energy_separation':distinct,
    'roots':roots,'old_closed_P1':closed,
    'charged_stationary_blocks':charged,
    'physical_P1_source':None,
},indent=2))
