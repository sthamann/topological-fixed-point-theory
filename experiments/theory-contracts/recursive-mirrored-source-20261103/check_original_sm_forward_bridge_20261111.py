"""Replay the original SM carrier and test the old covariance forward map."""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
import json
import math

import numpy as np
import sympy as s

from v258_dirac_covariance_induction import induce
from v724_phys_t3b_modular_flows import derive_internal_objects, Q


def carrier():
    y = [F(-1,3)]*3 + [F(1,2)]*2
    by_occupancy = Counter((sum(i<3 for i in bits),sum(i>=3 for i in bits),sum(y[i] for i in bits))
                           for degree in (0,2,4) for bits in combinations(range(5),degree))
    name = {(1,1):'Q',(2,0):'u^c',(2,2):'d^c',(3,1):'L',(0,2):'e^c',(0,0):'nu^c'}
    rows = {name[a,b]:(multiplicity,charge) for (a,b,charge),multiplicity in by_occupancy.items()}
    assert rows == {'Q':(6,F(1,6)),'u^c':(3,F(-2,3)),'d^c':(3,F(1,3)),
                    'L':(2,F(-1,2)),'e^c':(1,F(1)),'nu^c':(1,F(0))}
    linear,cubic=(sum(m*q**k for m,q in rows.values()) for k in (1,3))
    mixed_weak=3*rows['Q'][1]+rows['L'][1]
    mixed_color=2*rows['Q'][1]+rows['u^c'][1]+rows['d^c'][1]
    assert [linear,cubic,mixed_weak,mixed_color,2-1-1,(3+1)%2]==[0]*6
    trace=sum(m*q*q for m,q in rows.values())
    beta=(F(2,3)*F(3,5)*3*trace+F(1,3)*F(3,5)*2*F(1,2)**2,
          -F(22,3)+F(2,3)*3*4*F(1,2)+F(1,3)*F(1,2),
          -11+F(2,3)*3*4*F(1,2))
    assert beta==(F(41,10),F(-19,6),F(-7))
    return {'spinor':{k:[m,str(q)] for k,(m,q) in rows.items()},'anomalies':[0]*6,
            'beta':[str(v) for v in beta]}


def determinant_line():
    # External weak two-form contracts the already marked five-slot volume.
    matrix=s.Matrix(3,3,lambda a,j:s.LeviCivita(a,*list(combinations(range(3),2))[j],3,4))
    assert matrix.det() == 1 and matrix.rank()==3
    y=s.symbols('h0 h1',complex=True)
    h=s.Matrix(y)
    theta=s.Matrix([[0,1],[-1,0]])*h.conjugate()
    assert h.row_join(h).det()==0
    assert s.simplify(h.row_join(theta).det()+sum(v*s.conjugate(v) for v in y))==0
    return {'contraction_from_marked_volume':str(matrix), 'weak_pairing_h_h':0,
            'weak_pairing_h_conjugate':str(h.row_join(theta).det()),
            'hypercharge_color_determinant':-1,'hypercharge_weak_determinant':1,
            'hypercharge_higgs_and_adjoint_pair':0}


def original_forward_covariance():
    old=derive_internal_objects()
    cov=old.mother_covariances[math.inf]
    gamma=np.diag([1,1,1,-1,-1,-1])
    assert np.linalg.norm(gamma@cov@gamma-(np.eye(6)-cov))<1e-12
    logit=induce(cov)
    assert np.linalg.norm(logit[:3,:3])+np.linalg.norm(logit[3:,3:])<1e-11
    tangent=Q@np.diag([0,1,1])
    assert np.linalg.norm(logit[:3,3:]-tangent)<1e-11
    x=s.Symbol('x')
    exact_tangent=s.Matrix(tangent.astype(int))
    assert s.expand((exact_tangent.T*exact_tangent).charpoly(x).as_expr()-x*(x**2-10*x+5))==0
    left_null=s.Matrix([-2,1,0])/s.sqrt(5)
    right_null=s.eye(3)[:,0]
    assert left_null.T*exact_tangent==s.zeros(1,3) and exact_tangent*right_null==s.zeros(3,1)
    active_first=s.Matrix(Q.astype(int))*s.diag(1,0,0)
    assert (left_null.T*active_first*right_null)[0]==-3/s.sqrt(5)
    # For the v258 full logit on the half-occupied doublet, this is the
    # actual Fréchet coefficient relating a missing covariance to a mass.
    c=s.Symbol('c',positive=True)
    assert s.diff(s.log((1-c)/c),c).subs(c,s.Rational(1,2))==-4
    absolute=np.linalg.svd(logit[:3,3:],compute_uv=False)
    assert absolute[-1]<1e-11
    corner_gaps=[float(min(abs(np.linalg.eigvalsh(c)-.5))) for c in old.covariances[math.inf]]
    # A corner with no half-occupied state has no particle-hole
    # involution taking C to I-C on that odd-dimensional corner alone.
    assert sum(gap>1e-10 for gap in corner_gaps)>0
    return {'mother_covariance_spectrum':np.linalg.eigvalsh(cov).tolist(),
            'old_mass_free_logit_singulars':absolute.tolist(),'corner_half_occupancy_gaps':corner_gaps,
            'left_null':'(-2,1,0)/sqrt(5)','right_null':'(1,0,0)',
            'Q_first_column_null_overlap':'-3/sqrt(5)','logit_tangent_at_half_occupancy':-4}


def leakage_map():
    # Exact one-dimensional blocks of a pure projection. L2+C^2=C.
    c=s.Symbol('c',positive=True)
    g=s.log((1-c)/c)
    assert s.simplify((s.exp(g)+2+s.exp(-g))*c*(1-c))==1
    sizes=[(3,1),(48,16),(96,32)]
    # For any genuine pure completion of rank-three faithful C, one needs
    # three orthogonal complement channels. Demonstrate attainment of this
    # rank bound with a projection independent of old masses.
    c3=s.diag(s.Rational(1,4),s.Rational(2,5),s.Rational(1,2))
    b=(c3*(s.eye(3)-c3)).applyfunc(s.sqrt)
    two=c3.row_join(b).col_join(b.row_join(s.eye(3)-c3))
    assert two*two==two and b.rank()==3
    return {'one_family_complement':{'kept_complement_min_pure':[[n,m,n-m] for n,m in sizes]},
            'completion_identity':'|logit(c)| = 2 arcosh(1/(2 sigma(B))) for 0<c<1',
            'chirality_condition':'gamma C gamma = I-C is equivalent to an odd full logit on a faithful covariance'}


if __name__=='__main__':
    print(json.dumps({'original_carrier':carrier(),'line_test':determinant_line(),
                      'original_v724_v258':original_forward_covariance(),'pure_projection':leakage_map()},indent=2))
