#!/usr/bin/env python3
"""Exact conditional origin tests. No TFPT constants or dimensions are inputs.

Assumptions: isotropic rank-two cone R + R^d; an equivariant alternating
bilinear vector bracket; and, for the spacetime comparison only, Herm_2(C).
The first test selects d=3 inside this stated class, not physical spacetime.
"""
from itertools import combinations, product
import json
import sympy as s


def bracket_invariants(d):
    # c[i,j,k] is the coefficient of e_k in [e_i,e_j], i<j.
    variables = [(i,j,k) for i,j in combinations(range(d),2) for k in range(d)]
    # Exact pi rotations in every coordinate plane eliminate coefficients
    # with odd transformation sign. This is necessary for SO(d) invariance.
    candidates=[]
    for v in variables:
        if all(sum(x in (a,b) for x in v)%2 == 0
               for a,b in combinations(range(d),2)):
            candidates.append(v)
    if not candidates:
        return {'d':d,'alternating_coefficients':len(variables),
                'after_pi_rotations':0,'invariant_dimension':0}
    symbols=s.symbols('c:'+str(len(candidates)))
    lookup=dict(zip(candidates,symbols))
    def coefficient(i,j,k):
        if i==j: return s.Integer(0)
        return lookup.get((i,j,k),s.Integer(0)) if i<j else -lookup.get((j,i,k),s.Integer(0))
    def generator_action(a,b,j):
        return {b:1} if j==a else ({a:-1} if j==b else {})
    equations=[]
    for a,b in combinations(range(d),2):
        for i,j in combinations(range(d),2):
            for k in range(d):
                lhs=sum(coefficient(i,j,l)*generator_action(a,b,l).get(k,0) for l in range(d))
                rhs=sum(v*coefficient(l,j,k) for l,v in generator_action(a,b,i).items())
                rhs+=sum(v*coefficient(i,l,k) for l,v in generator_action(a,b,j).items())
                equations.append(s.expand(lhs-rhs))
    matrix,_=s.linear_eq_to_matrix(equations,symbols)
    basis=matrix.nullspace()
    return {'d':d,'alternating_coefficients':len(variables),
            'after_pi_rotations':len(candidates),'invariant_dimension':len(basis),
            'basis':[{str(k):str(v) for k,v in zip(candidates,w) if v} for w in basis]}


def cone_and_dynamics():
    t,x,y,z=s.symbols('t x y z',real=True)
    H=s.Matrix([[t+z,x-s.I*y],[x+s.I*y,t-z]])
    determinant=s.factor(H.det())
    assert s.expand(determinant-(t*t-x*x-y*y-z*z))==0
    # Rational determinant-one congruence: exact boost cosh(eta)=17/8,
    # sinh(eta)=15/8. The value is a verification sample, not a fit.
    B=s.diag(2,s.Rational(1,2))
    Hp=B*H*B.conjugate().T
    coordinates=[s.simplify(s.trace(Hp)/2),s.simplify((Hp[0,1]+Hp[1,0])/2),
                 s.simplify((Hp[1,0]-Hp[0,1])/(2*s.I)),s.simplify((Hp[0,0]-Hp[1,1])/2)]
    assert s.expand(Hp.det()-H.det())==0
    # Every full SO(d)-invariant bracket found in d=3 is proportional to
    # epsilon. Test Jacobi exactly, and the Pauli algebra correspondence.
    pauli=[s.Matrix([[0,1],[1,0]]),s.Matrix([[0,-s.I],[s.I,0]]),s.diag(1,-1)]
    for i,j in product(range(3),repeat=2):
        assert pauli[i]*pauli[j]-pauli[j]*pauli[i] == 2*s.I*sum((s.LeviCivita(i,j,k)*pauli[k] for k in range(3)),s.zeros(2))
    jacobi=[]
    for a,b,c,d in product(range(3),repeat=4):
        jacobi.append(sum(s.LeviCivita(b,c,k)*s.LeviCivita(a,k,d)+s.LeviCivita(c,a,k)*s.LeviCivita(b,k,d)+s.LeviCivita(a,b,k)*s.LeviCivita(c,k,d) for k in range(3)))
    assert all(v==0 for v in jacobi)
    return {'Herm2C_determinant':str(determinant),'boost_transformed_coordinates':list(map(str,coordinates)),
            'boost_determinant_residual':str(s.expand(Hp.det()-H.det())),
            'Pauli_and_Jacobi_verified':True,
            'rank_two_division_algebra_cones':{'R':'2+1','C':'3+1','H':'5+1','O_rank_two_only':'9+1'}}


def finite_boost_obstruction():
    # Verify the general trace proof on unrestricted Hermitian 2x2 symbols.
    # General proof: Tr(P[K,P])=Tr(KP^2)-Tr(KP^2)=0 for every finite size.
    k0,kx,ky,kz,p0,px,py,pz,c=s.symbols('k0 kx ky kz p0 px py pz c',real=True)
    K=s.Matrix([[k0+kz,kx-s.I*ky],[kx+s.I*ky,k0-kz]])
    P=s.Matrix([[p0+pz,px-s.I*py],[px+s.I*py,p0-pz]])
    A=s.I*(K*P-P*K)
    hs2=lambda M:s.simplify(s.trace(M.conjugate().T*M))
    cross=s.simplify(s.trace(P*A))
    residual=s.simplify(hs2(A-c*P)-hs2(A)-c*c*hs2(P))
    assert cross==0 and residual==0
    return {'trace_cross_term':str(cross),'squared_norm_identity_residual':str(residual),
            'theorem':'||i[K,P]-c P||_HS^2 = ||i[K,P]||_HS^2 + c^2 ||P||_HS^2',
            'relative_lower_bound_for_nonzero_P':'1, when divided by |c| ||P||_HS',
            'exact_finite_boost_for_c_nonzero':'P = 0 only',
            'scope':'finite Hermitian matrices, full Hilbert-Schmidt norm; not projected sectors or unbounded continuum limits'}


def local_tomography_baseline():
    # Standard n-level real/complex dimensions; quaternionic table is only
    # a naive same-field carrier comparison, NOT a valid compositing rule.
    dim={'R':lambda n:n*(n+1)//2,'C':lambda n:n*n,'H_naive':lambda n:n*(2*n-1)}
    return {field:{'one_bit':f(2),'product_effect_count':f(2)**2,'naive_four_level_dimension':f(4),
                   'matches_local_tomography_count':f(2)**2==f(4)} for field,f in dim.items()}


if __name__=='__main__':
    out={'assumptions':['homogeneous self-dual rank-two spin-factor cone','all orientation-preserving rotations act as symmetries',
                        'a nonzero alternating bilinear vector bracket is equivariant under those rotations'],
         'scan':[bracket_invariants(d) for d in range(2,9)],
         'cone':cone_and_dynamics(),'finite_boost':finite_boost_obstruction(),
         'tomography_comparison':local_tomography_baseline()}
    print(json.dumps(out,indent=2,ensure_ascii=False))
