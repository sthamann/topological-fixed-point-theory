"""Coherent comparison of two primitive loops in the original four-spinor algebra."""
from argparse import ArgumentParser
from itertools import combinations
from pathlib import Path
import sys
import sympy as s
from rahmen import frame

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"familienmessung"))
from flaggen import transfers

def simplify(M):
    return M.applyfunc(lambda x:s.simplify(s.expand(x)))

def simplex():
    # Isometry from the six-label zero-sum space to tangent coordinates.
    tangent=s.Matrix([[s.S.One/s.sqrt(k*(k+1)) if j<k else
                       -k/s.sqrt(k*(k+1)) if j==k else s.S.Zero
                       for j in range(6)] for k in range(1,6)])
    d=s.sqrt(s.Rational(6,5))*tangent
    assert simplify(d.T*d)==(6*s.eye(6)-s.ones(6))/5
    assert simplify(d*d.T)==6*s.eye(5)/5
    return tangent,d

def main(root):
    phi,gammas,_=frame(root)
    T,d=simplex()
    directions=[d[:,j] for j in range(6)]
    I=s.eye(4)
    cliff=lambda v:sum((k*g for k,g in zip(v,gammas)),s.zeros(4))
    unitaries=[-cliff(v) for v in directions]
    for u in unitaries:
        assert u==u.H and simplify(u* u)==I
    probs={-1:s.Rational(3,5),1:s.Rational(2,5)}
    # One common terminal phase from identical positive-primitive convention.
    mixed=s.zeros(16)
    ad=lambda U:s.kronecker_product(U.conjugate(),U)
    for U in unitaries:
        mixed+=5*ad(U)/72
    assert simplify(7*s.eye(16)/12+mixed)==phi
    print("Six primitive axes + binary phase erasure + idle 1/6 reproduce native channel: 256/256")
    A=transfers()
    one=s.ones(21,1)
    base=7
    def prob(i,j):
        return (one.T*A[base]*A[i]*A[j]*A[base]*one)[0]/(21*8**4)
    mass=sum((2*prob(i,j) for i,j in combinations(range(1,7),2)),s.S.Zero)
    mean=s.zeros(5)
    sum_control=s.zeros(5)
    uniform=s.zeros(5)
    marginal=[s.S.Zero]*6
    for i,j in combinations(range(6),2):
        pair_weight=2*prob(i+1,j+1)/mass
        collinear=(i+1)^(j+1)==base
        assert pair_weight==(s.Rational(1,10) if collinear else s.Rational(7,120))
        assert prob(i+1,j+1)==prob(j+1,i+1)
        marginal[i]+=pair_weight/2
        marginal[j]+=pair_weight/2
        U,V=unitaries[i],unitaries[j]
        assert simplify(sum((ad((U+sign*V)/2) for sign in [-1,1]),s.zeros(16))
                        -(ad(U)+ad(V))/2)==s.zeros(16)
        if (i,j)==(0,1):
            out=[]
            for phase in (1,s.I,-1):
                lost=(U-phase*V)/2
                out.append(simplify(lost.H*lost).trace()/4)
            assert out==[s.Rational(3,5),s.Rational(1,2),s.Rational(2,5)]
            print("Antisymmetric port at relative phases 0, pi/2, pi in the tracial state:",out)
        for sign in (-1,1):
            K=(U+sign*V)/2
            assert simplify(K.H*K)==probs[sign]*I
            ray=simplify(d[:,i]+sign*d[:,j])
            proj=simplify(ray*ray.T/(ray.T*ray)[0])
            if sign==-1:
                mean+=pair_weight*proj
                uniform+=proj/15
            else:
                sum_control+=pair_weight*proj
    family=s.zeros(6)
    for i in range(1,7):
        j=i^base
        if i<j:
            v=s.eye(6)[:,i-1]-s.eye(6)[:,j-1]
            family+=v*v.T/2
    F=simplify(T*family*T.T)
    assert F.rank()==3 and simplify(F*F)==F
    assert simplify(mean-s.Rational(7,40)*s.eye(5)-F/24)==s.zeros(5)
    assert simplify(uniform)==s.eye(5)/5
    assert marginal==[s.Rational(1,6)]*6
    # The fixed original trine's two-axis contrast averages the squared
    # component of the direction along its second tangent vector to 1/5.
    midpoint=(1+s.Rational(1,5))/3
    assert midpoint==s.Rational(2,5)
    assert (1-s.Rational(1,36))*midpoint==s.Rational(7,18)
    print("Mean fixed-frame trine contrast at the second mark and half turn:",s.Rational(7,18),midpoint)
    print("Pair weights: 1/10 per collinear pair, 7/120 per other pair; branch marginal 1/6")
    print("Antisymmetric port norm squared",probs[-1],"and symmetric port",probs[1],"for 15/15 pairs")
    print("Antisymmetric axis response:",simplify(mean).eigenvals())
    print("Symmetric axis response:",simplify(sum_control).eigenvals())
    print("Uniform-pair antisymmetric control:",simplify(uniform).eigenvals())
    # Label selection conditional on a true change is unbiased in the
    # stationary two-word record; rows for a fixed previous mark.
    for u in range(1,8):
        assert len({(one.T*A[u]*A[v]*one)[0] for v in range(1,8) if u!=v})==1

if __name__=="__main__":
    p=ArgumentParser()
    p.add_argument("--repo",type=Path,default=Path(__file__).resolve().parents[4])
    main(p.parse_args().repo)
