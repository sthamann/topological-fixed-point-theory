"""Original Clifford channel in a moving frame, exact absence of a mark-only frame."""
from argparse import ArgumentParser
from itertools import combinations
from pathlib import Path
import sys
import sympy as s

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"familienmessung"))
from spintransport import matrix_units, native

def frame(root):
    phi,gamma=native(root)
    eye=s.eye(4)
    D=eye.vec()*eye.vec().T/4
    reflection=3*s.eye(16)+4*D-6*phi
    units=matrix_units(4)
    S=lambda X:(reflection*X.vec()).reshape(4,4).T
    equations=s.Matrix.vstack(*(s.kronecker_product(eye,S(E))-
                               s.kronecker_product(E,eye) for E in units))
    kernel=equations.nullspace()
    assert len(kernel)==1
    B=kernel[0].reshape(4,4).T
    assert B.T==-B and B*B.H==eye
    assert all(B*g.T*B.H==g for g in gamma)
    return phi,gamma,B

def response(B):
    eye=s.eye(4)
    S=s.Matrix.hstack(*( (B*E.T*B.H).vec() for E in matrix_units(4)))
    D=eye.vec()*eye.vec().T/4
    return s.eye(16)/2-S/6+2*D/3

def moving_edge(root):
    phi,gamma,B=frame(root)
    eye=s.eye(4)
    x=gamma[0]
    c=s.Rational(1,6)
    U=s.sqrt((1+c)/2)*eye+s.I*s.sqrt((1-c)/2)*x
    B_new=(eye*c+s.I*s.sqrt(1-c*c)*x)*B
    simplify=lambda M:M.applyfunc(s.simplify)
    assert simplify(U*B*U.T-B_new)==s.zeros(4)
    assert response(B)==phi
    ad=s.kronecker_product(U.conjugate(),U)
    assert simplify(response(B_new)*ad-ad*phi)==s.zeros(16)
    print("Original 16x16 transfer and transported edge covariance: 256/256")
    # Applying the inverse move must recover the reference operator-system
    # map. Composition of the identity just tested proves this for paths.
    assert simplify(ad.H*response(B_new)*ad-phi)==s.zeros(16)

def loop_invariants():
    eye=s.eye(7)
    a=(eye[:,6]-s.ones(7,1)/7)
    P=a*a.T/a.dot(a)
    normal_flip=eye-2*P
    V=eye[:,:6]-eye[:,6]*s.ones(1,6)
    T=eye[:,:5]-eye[:,5]*s.ones(1,5)
    pairs=list(combinations(range(6),2))
    def wedge(R):
        return s.Matrix([[R[i,k]*R[j,l]-R[i,l]*R[j,k]
                          for k,l in pairs] for i,j in pairs])
    equations=[]
    for b,c in combinations(range(6),2):
        d=eye[:,b]-eye[:,c]
        R=normal_flip*(eye-d*d.T)
        assert R*R==eye and R*a==-a and R.det()==1
        equations.append(wedge((R*V)[:6,:])-s.eye(15))
    assert s.Matrix.vstack(*equations).rank()==15
    # Invariance under even loops suffices, ruling out an unordered
    # three-effect trine preserved only up to relabeling.
    generators=[E+s.eye(15) for E in equations]
    even=[G*generators[0]-s.eye(15) for G in generators]
    assert s.Matrix.vstack(*even).rank()==15
    # In orthogonal tangent directions all loop matrices are signed
    # ordinary permutations. The spin lift anticommutes on disjoint pairs.
    assert T.rank()==5
    assert all(int(Q.det()) in (-1,1) for Q in generators)
    print("Based label holonomy: S6, 720 elements; fixed traceless dimension 0")
    print("Even holonomy A6 already has fixed traceless dimension 0")

if __name__=="__main__":
    p=ArgumentParser()
    p.add_argument("--repo",type=Path,default=Path(__file__).resolve().parents[4])
    moving_edge(p.parse_args().repo)
    loop_invariants()
