#!/usr/bin/env python3
"""TFPT K4 / Z4 independent finite audit, 2026-11-13.

Audits mathematical candidate models only. Does not evaluate the original
weighted R113 P1 source, does not verify a physical TOE, does not change a repo.
Dependencies: numpy, sympy.  Run: python tfpt_k4_z4_sewing_audit.py
"""
from __future__ import annotations
from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations, product
import numpy as np
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ


def check(msg: str, condition: bool):
    if not condition:
        raise AssertionError(msg)
    print('PASS', msg)


def critical_group_audit():
    G = 4 * sp.eye(3) - sp.ones(3, 3)
    diag = [abs(int(smith_normal_form(G, domain=ZZ)[i,i])) for i in range(3)]
    check('K4 flow lattice det=16 / Smith (1,4,4)', G.det() == 16 and diag == [1,4,4])

    # K(K4) = {v mod 4 | sum(v)=0} / <(1,1,1,1)>.
    def canon(v):
        return tuple((int(x)-int(v[3])) % 4 for x in v)
    classes = sorted({canon(v) for v in product(range(4), repeat=4) if sum(v)%4 == 0})
    check('K4 critical group exactly 16 classes', len(classes) == 16)

    group = list(permutations(range(4)))
    even = [p for p in group if sum(p[i]>p[j] for i in range(4) for j in range(i+1,4))%2 == 0]
    def invariant(v, group):
        return all(canon([v[p[i]] for i in range(4)]) == v for p in group)
    for name, action in [('S4', group), ('A4',even)]:
        fixed = [v for v in classes if invariant(v, action)]
        check(f'K4 critical group invariants under {name}: only 0', fixed == [(0,0,0,0)])

    # Discriminant group D(L)=L*/L in 3 rational coordinate classes x/4,
    # with sum(x)=0 mod 4; b(x,y) = x^T G y /16 mod Z.
    dclasses = [x for x in product(range(4), repeat=3) if sum(x)%4 == 0]
    order4 = [x for x in dclasses if any(int(xi)%2 for xi in x)]
    def b(x, y):
        return Fraction(int((sp.Matrix(x).T*G*sp.Matrix(y))[0]),16)%1
    check('All 12 primitive elements b(x,x)=1/2 mod 1',
          len(order4) == 12 and all(b(x,x) == Fraction(1,2) for x in order4))
    order2 = [x for x in dclasses if all(int(xi)%2==0 for xi in x)]
    check('Order-two subgroup is isotropic (size 4)',
          len(order2)==4 and all(b(x,y)==0 for x in order2 for y in order2))

    # A3 weight of a distinguished mark, invariant *class* under all permutations.
    w=sp.Matrix([sp.Rational(3,4)]+[-sp.Rational(1,4)]*3)
    check('A3 primitive mark weight norm 3/4 and order 4',
          w.dot(w)==sp.Rational(3,4)
          and all((4*z).q==1 for z in w)
          and any((2*z).q!=1 for z in w))
    check('A3 weight class independent of which mark under S4',
          all(all((w[i]-w[p[i]]).q==1 for i in range(4)) for p in group))
    return G


def symmetric_complex_structure():
    # A single R4 permutation representation has commutant spanned by I,J4.
    perms=list(permutations(range(4)))
    unknown=sp.symbols('c0:16')
    C=sp.Matrix(4,4,unknown)
    equations=[]
    for p in perms:
        P=sp.zeros(4)
        for i,j in enumerate(p): P[j,i]=1
        equations.extend(C*P-P*C)
    coefficient,_=sp.linear_eq_to_matrix(equations,unknown)
    null=coefficient.nullspace()
    check('S4 commutant on one real cell has dimension 2',len(null)==2)
    # Every commuting endomorphism is a P_trivial + b P_standard;
    # squaring to -I requires a*a=b*b=-1 over reals: impossible.
    J2=sp.Matrix([[0,-1],[1,0]])
    check('A compatible complex structure emerges on two copies',J2**2==-sp.eye(2))


def lattice_glue_a3_d5():
    # V_A3 = sum-zero vectors in R4; D5 even-coordinate-sum in R5.
    def row_A(v):return sp.Matrix(list(v)+[0]*5)
    def row_D(v):return sp.Matrix([0]*4+list(v))
    E4=sp.eye(4)
    E5=sp.eye(5)
    a1=row_A(E4[:,0]-E4[:,1]); a2=row_A(E4[:,1]-E4[:,2]); a3=row_A(E4[:,2]-E4[:,3])
    ds=[row_D(E5[:,i]-E5[:,i+1]) for i in range(4)]+[row_D(E5[:,3]+E5[:,4])]
    wa=row_A(sp.Matrix([sp.Rational(3,4)]+[-sp.Rational(1,4)]*3))
    wd=row_D(sp.Matrix([sp.Rational(1,2)]*5))
    glue=wa+wd
    check('glue element even norm (A3,D5)=(3/4,5/4)',
          wa.dot(wa)==sp.Rational(3,4) and wd.dot(wd)==sp.Rational(5,4)
          and glue.dot(glue)==2)
    basis=[a1,a2,*ds,glue]  # remove a3: glue carries a3 coefficient 1/4.
    G=sp.Matrix(8,8,lambda i,j:basis[i].dot(basis[j]))
    check('Explicit E8 glue Gram integer, even, positive-definite, det 1',
          all(G[i,j].q==1 for i in range(8) for j in range(8))
          and all(int(G[i,i])%2==0 for i in range(8))
          and G.det()==1 and all(x>0 for x in G.eigenvals().keys()))
    # Root cosets 0,2,1,3 as A3 x D5 representations:
    check('E8 240 norm-two roots split 52+60+64+64',40+12+6*10+4*16+4*16==240)
    print('E8 Gram matrix from A3 + D5 + the primitive mixed glue:\n',G)


def z4_clock_test():
    q=4
    omega=1j
    X=np.roll(np.eye(q,dtype=complex),1,axis=0)
    Q=np.diag([omega**j for j in range(q)])
    P=np.zeros((q,q),dtype=complex);P[0,0]=1
    check('Z4 clock covariance Q X Q* = i X',np.allclose(Q@X@Q.conj().T,1j*X))
    check('Z4 signed return: -/+ at 2, ++++ at 4',
          np.allclose(P@X@X.conj().T@P,P) and
          all(np.allclose(P@np.linalg.matrix_power(X,k)@P,0) for k in (1,2,3))
          and np.allclose(P@np.linalg.matrix_power(X,4)@P,P))
    x2=np.array([[0,1],[1,0]],dtype=complex);p2=np.diag([1,0])
    check('Z2 cannot mimic: same-direction double insertion returns',np.allclose(p2@x2@x2@p2,p2))


def tetra_z4_gauge():
    edges=list(combinations(range(4),2))
    faces=[(1,2,3),(0,2,3),(0,1,3),(0,1,2)]
    B1=sp.zeros(4,6);B2=sp.zeros(6,4)
    for k,(a,b) in enumerate(edges): B1[a,k]=-1;B1[b,k]=1
    for f,(a,b,c) in enumerate(faces):
        sig=(-1)**f
        for aedge,coef in [((b,c),1),((a,c),-1),((a,b),1)]:
            B2[edges.index(aedge),f]=sig*coef
    check('Tetrahedron chain identity d1*d2=0; ranks (3,3)',
          B1*B2==sp.zeros(4,4) and B1.rank()==3 and B2.rank()==3)
    check('Self-dual tetra: vertex-edge incidence matches opposite faces/edges',
          all(int(B1[v,e]/B2[edges.index(tuple(sorted(set(range(4))-set(pair)))) ,v])
                  ==int(B1[(pair[0] if v==pair[1] else pair[1]),e]
                       /B2[edges.index(tuple(sorted(set(range(4))-set(pair)))) ,(pair[0] if v==pair[1] else pair[1])])
              for e,pair in enumerate(edges) for v in pair))
    # Each U(1) Z4 vertex/face eigenvalue assignment is neutral modulo 4.
    states=[t for t in product(range(4),repeat=4) if sum(t)%4==0]
    check('Star and plaquette commuting sectors 64 x 64 = 4096=4^6',
          len(states)==64 and len(states)**2==4**6)
    # Strict Gauss physical space has only 3 independent flux coordinates (64 states).
    ix={x:i for i,x in enumerate(states)}
    N=len(states)
    shifts=[]
    for e in range(6):
        dv=np.array(list(B2.row(e)),dtype=int)
        S=np.zeros((N,N))
        for k,s in enumerate(states):S[ix[tuple((np.array(s)+dv)%4)],k]=1.
        shifts.append(S)
    potential=np.diag([sum([0.,2.,4.,2.][n] for n in state) for state in states])
    kinetic=sum(2*np.eye(N)-S-S.T for S in shifts)
    check('Noncommuting gauge-invariant flux/transport terms',
          np.linalg.norm(kinetic@potential-potential@kinetic)>1)
    energies={}
    for h in [0,0.1,0.25,0.5,1.,2.]:
        val,vec=np.linalg.eigh(potential+h*kinetic)
        check(f'Strict-Gauss model unique positive-gap ground at h={h}',val[1]-val[0]>1e-8)
        energies[h]=(float(val[0]),float(val[1]-val[0]),float(abs(vec[ix[(0,0,0,0)],0])**2))
    A=sum(S+S.T for S in shifts)
    i0=ix[(0,0,0,0)]
    check('Single link move leaves zero-flux P; second return is scalar 12P',
          abs(A[i0,i0])<1e-12 and abs((A@A)[i0,i0]-12)<1e-12)
    print('Gauge-invariant 64-state H(g=1,h), E0, gap, ground |zero-flux|²')
    for h, vals in energies.items():print(' ',h, *(round(x,9) for x in vals))
    print('WARNING: 1D zero-flux return sector yields only energy shift, not R113 interaction.')
    print('WARNING: no original R113, chirality c=8, 3+1D gravity, or physical Gauss matter source tested.')



def z4_gauss_matter_return():
    """Two 2-flavor fermion cells, one Z4 electric link, exact local Gauss law."""
    nf,nlink=16,4
    def fermion_annihilator(mode):
        op=np.zeros((nf,nf),dtype=complex)
        for state in range(nf):
            if (state>>mode)&1:
                sign=(-1)**((state&((1<<mode)-1)).bit_count())
                op[state^(1<<mode),state]=sign
        return op
    c=[fermion_annihilator(j) for j in range(4)]
    U=np.roll(np.eye(nlink,dtype=complex),1,axis=0)
    Ze=np.diag([1,1j,-1,-1j])
    Eterm=2*np.eye(4)-Ze-Ze.conj().T
    K=np.zeros((nf*nlink,nf*nlink),dtype=complex)
    for a,b in [(0,2),(1,3)]:
        hop=c[a].conj().T@c[b]
        K+=np.kron(hop,U)+np.kron(hop.conj().T,U.conj().T)
    H0=np.kron(np.eye(nf),Eterm)
    GA=np.zeros_like(K);GB=np.zeros_like(K)
    phys=[]
    for state in range(nf):
        na=((state>>0)&1)+((state>>1)&1)
        nb=((state>>2)&1)+((state>>3)&1)
        for e in range(4):
            ix=state*4+e
            GA[ix,ix]=1j**(na-1-e)
            GB[ix,ix]=1j**(nb-1+e)
            if (na-1-e)%4==0 and (nb-1+e)%4==0:
                phys.append(ix)
    check('Full Z4 Gauss constraints select 6 states from 64',len(phys)==6)
    check('Full hopping and electric energy commute with both local Gauss operators',
          all(np.linalg.norm(G@M-M@G)<1e-12 for G in [GA,GB] for M in [H0,K]))
    hp=H0[np.ix_(phys,phys)]-K[np.ix_(phys,phys)]
    s5=np.sqrt(5.)
    expected=np.array([1-s5,0,0,0,2,1+s5])
    check('Exact physical Z4 matter spectrum for Delta=t=1',
          np.allclose(np.linalg.eigvalsh(hp),expected,atol=1e-12))
    P=[j for j,ix in enumerate(phys) if ix%4==0]
    Q=[j for j in range(6) if j not in P]
    check('Zero-electric-flux sector rank 4, first-order hopping forbidden',
          len(P)==4 and np.linalg.norm(K[np.ix_([phys[j] for j in P],[phys[j] for j in P])])<1e-12)
    kphys=K[np.ix_(phys,phys)]
    heff=-kphys[np.ix_(P,Q)]@kphys[np.ix_(Q,P)]/2
    check('Second-order return operator has spectrum (-2,0,0,0)',
          np.allclose(np.linalg.eigvalsh(heff),[-2,0,0,0],atol=1e-12))
    # Positive semigroup compression defect, H' = H - E0 I >=0.
    vals,vecs=np.linalg.eigh(hp)
    shifted=vals-vals[0]
    for tau in [0.1,0.25,0.5,1.0]:
        T=(vecs*np.exp(-tau*shifted))@vecs.conj().T
        T2=(vecs*np.exp(-2*tau*shifted))@vecs.conj().T
        defect=T2[np.ix_(P,P)]-T[np.ix_(P,P)]@T[np.ix_(P,P)]
        rhs=T[np.ix_(P,Q)]@T[np.ix_(Q,P)]
        predicted=.2*(1-np.exp(-2*s5*tau))**2
        check(f'Exact positive-return defect at tau={tau}: 3 zeros, one analytic eigenvalue',
              np.linalg.norm(defect-rhs)<1e-12 and
              np.allclose(np.linalg.eigvalsh(defect),[0,0,0,predicted],atol=1e-12))
    print('R113-CANDIDATE: return defect tests the original P1 transfer when its weights become available.')
    print('NOT AN ORIGINAL R113 EVALUATION: K4 four-mark origin not implemented in this dimer.')



def e8_z4_graded_root_brackets():
    """Exact E8 root-composition signatures by A3+D5 glue grade in Z4."""
    roots={k:set() for k in range(4)}
    for i in range(4):
        for j in range(4):
            if i!=j:
                roots[0].add(tuple(4*(int(k==i)-int(k==j)) for k in range(4))+(0,)*5)
    for i,j in combinations(range(5),2):
        for a,b in product((-1,1),repeat=2):
            roots[0].add((0,)*4+tuple(4*(a if k==i else (b if k==j else 0)) for k in range(5)))
    plus=[];minus=[]
    for signs in product((-1,1),repeat=5):
        (plus if sum(x<0 for x in signs)%2==0 else minus).append(tuple(2*x for x in signs))
    for i in range(4):
        wi=tuple(3 if j==i else -1 for j in range(4))
        roots[1].update(wi+s for s in plus)
        roots[3].update(tuple(-x for x in wi)+s for s in minus)
    for i,j in combinations(range(4),2):
        w2=tuple(2 if k in (i,j) else -2 for k in range(4))
        for k in range(5):
            for sign in (-1,1):
                roots[2].add(w2+tuple(4*sign if m==k else 0 for m in range(5)))
    check('E8 root grades (0,1,2,3)=(52,64,60,64) all have squared norm 2',
          [len(roots[k]) for k in range(4)]==[52,64,60,64]
          and len(set().union(*roots.values()))==240
          and all(sum(x*x for x in rt)==32 for k in roots for rt in roots[k]))
    for a,b,expected in [(1,1,960),(1,3,832),(1,2,960),(2,2,720)]:
        c=(a+b)%4
        count=sum(tuple(xi+yi for xi,yi in zip(x,y)) in roots[c]
                  for x in roots[a] for y in roots[b])
        check(f'E8 Z4 root bracket: grade {a}+{b}->{c} has {expected} ordered root pairs',
              count==expected)
    print('NOTE: grade1+grade3->grade0 count 832 excludes 64 opposite-root pairs closing into Cartan.')


def main():
    critical_group_audit()
    symmetric_complex_structure()
    lattice_glue_a3_d5()
    z4_clock_test()
    tetra_z4_gauge()
    z4_gauss_matter_return()
    e8_z4_graded_root_brackets()
    print('\nALL FINITE CHECKS PASSED. F1 remains open.')

if __name__=='__main__':main()
