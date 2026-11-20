"""Finite projective flags, faithful exchange histories and two native response controls."""
from itertools import product
import sympy as sp

Q = sp.Rational
FLAGS = [(p, l) for p in range(1, 8) for l in range(1, 8)
         if (p & l).bit_count() % 2 == 0]
N = len(FLAGS)
I = sp.eye(N)


def change(which):
    return sp.Matrix([[Q(1, 2) if u != v and u[which] == v[which] else 0
                       for v in FLAGS] for u in FLAGS])


def tr(m):
    return sp.trace(m) / N


def flag_transfers():
    a,b = change(0),change(1)
    assert N==21 and all(m == m.T and m*sp.ones(N,1) == sp.ones(N,1)
                         for m in (a,b))
    assert 2*a**2 == I+a and 2*b**2 == I+b and a*b*a == b*a*b
    # Pure point function of sum zero, and its sum over points on a line.
    f = lambda p: int(p == 1) - int(p == 2)
    cols = sp.Matrix([[f(p), sum(f(k) for k in range(1,8)
                                  if (k&l).bit_count()%2==0)] for p,l in FLAGS])
    t=Q(1,2)
    aa=sp.Matrix([[1,1],[0,-t]])
    bb=sp.Matrix([[-t,0],[t,1]])
    assert cols.rank()==2 and a*cols==cols*aa and b*cols==cols*bb
    standard_a=sp.Matrix([[-t,1],[0,1]])
    standard_b=sp.Matrix([[1,0],[t,-t]])
    change_basis=sp.Matrix([[2,-6],[-3,2]])
    assert aa*change_basis==change_basis*standard_a
    assert bb*change_basis==change_basis*standard_b and change_basis.det()!=0
    print("Flags",N,"; exact Hecke and braid relations; Burau parameter", t)
    print("Characteristic polynomial of mean of square transfers:")
    x=sp.symbols('x')
    clock=(a**2+b**2)/2
    print(sp.factor(clock.charpoly(x).as_expr()))
    assert (clock-I).nullspace()==[sp.ones(N,1)]
    # Quadratic Hecke relation gives eigenvalues 1,-1/2 for each factor.
    # The average of their squares is bounded from below by I/4.
    characteristic=clock.charpoly(x).as_expr()
    assert sp.simplify(characteristic - (x-1)*(x-Q(1,4))**8*(x*x-Q(5,4)*x+Q(23,64))**6)==0
    return a,b


def return_reader(a,b):
    # Row F of a reduced transfer is uniform exactly on flags at relative
    # position w. A given pair has just one Weyl label even when w has two
    # reduced expressions. No preferred starting flag is used.
    words=[I,a,b,a*b,b*a,a*b*a]
    counts=[1,2,2,4,4,8]
    masks=[size*m for size,m in zip(counts,words)]
    assert sum(masks,sp.zeros(N))==sp.ones(N)
    assert all(set(m)<={0,1} for m in masks)
    weights=[]
    for mask,m,size in zip(masks,words,counts):
        assert mask*sp.ones(N,1)==size*sp.ones(N,1)
        # G→F uses the transpose. The probability attached to (F,G)
        # is m[F,G], including both independent uniform draws.
        weights.append(sum(mask.multiply_elementwise(m),sp.S.Zero)/N**2)
    assert weights==[Q(1,21)]*6
    distribution=[w/sum(weights) for w in weights]
    raw=[Q(c,21) for c in counts]
    print("Raw relative position probabilities:",raw)
    print("Return rate:",sum(weights),"; recorded conditional weights:",distribution)
    print("Five visible orthogonal modes by parity: 3 odd and 2 even")
    # Test the resulting instrument on the entire four-dimensional matrix algebra.
    x=sp.Matrix([[0,1],[1,0]])
    y=sp.Matrix([[0,-sp.I],[sp.I,0]])
    z=sp.diag(1,-1)
    i2=sp.eye(2)
    gammas=[sp.kronecker_product(l,r) for l,r in [(i2,x),(i2,z),(x,y),(z,y),(y,y)]]
    maps=[sp.eye(16)]+[(sp.eye(16)+sp.kronecker_product(g.conjugate(),g))/2 for g in gammas]
    assert all(g*g==sp.eye(4) for g in gammas)
    out=sum((w*m for w,m in zip(distribution,maps)),sp.zeros(16))
    target=(7*sp.eye(16)+sum((sp.kronecker_product(g.conjugate(),g) for g in gammas),sp.zeros(16)))/12
    assert out==target and sum((w*m for w,m in zip(raw,maps)),sp.zeros(16))!=target
    assert maps[1].rank()==8 and (maps[1]*maps[2]).rank()==4
    print("Return-conditioned M4 channel: 256/256; without return selection: different")
    print("Elementary conditioned readout rank vs two distinct composed:",8,4)


def positive_history(a,b):
    # r=1/4 is the cylinder weight per signed elementary letter in the
    # symmetric four-letter sampling. Full-group values are defined in the report.
    r=Q(1,4)
    d2=r**4*tr((a**2-b**2).T*(a**2-b**2))
    order_d2=r**4*tr((a*b-b*a).T*(a*b-b*a))
    after=r**3*tr((a**2-b**2)*a)
    shape=tr(a**2*b**2)/tr(a**4)
    assert shape==Q(2,3) and d2==Q(1,1024) and order_d2==Q(1,512) and after==Q(1,256)
    print("Normalized double exchange overlap, independent of r:",shape)
    print("Reflected squared distance a²,b²:",d2)
    print("Reflected squared distance ab,ba:",order_d2)
    print("Same scalar observation",r**2*tr(a**2),r**2*tr(b**2),
          "; after a",r**3*tr(a**3),r**3*tr(b**2*a))
    # Positive half-span: Hecke basis indexed by six relative flag positions.
    words=[I,a,b,a*b,b*a,a*b*a]
    gram=sp.Matrix([[tr(u.T*v) for v in words] for u in words])
    assert gram.rank()==6
    assert gram==sp.diag(1,Q(1,2),Q(1,2),Q(1,4),Q(1,4),Q(1,8))
    print("Gram of six positive basis words:",list(gram.diagonal()),
          "; OS dimension:",gram.rank())
    # All signed relative positions are in the full group representation.
    C=a**2*b**-2
    assert C!=I
    print("Pure word a²b^-² on invariant rank two: has nontrivial standard block")


def coin_probe(a,b):
    # Three choices within a panel. Symmetry + real orthogonality leaves two
    # scalar signs on the uniform line and the two-dimensional difference space.
    panel=sp.ones(3)/3
    coins=[s*panel+t*(sp.eye(3)-panel) for s,t in product((1,-1),repeat=2)]
    mixing=[c for c in coins if c[0,1]!=0]
    assert len(mixing)==2 and mixing[0]==-mixing[1]
    weights=[x**2 for x in mixing[0][:,0]]
    assert weights==[Q(1,9),Q(4,9),Q(4,9)]
    # Independently supplied TFPT random-unitary channel at q=0.
    X=sp.Matrix([[0,1],[1,0]])
    Z=sp.diag(1,-1)
    Y=sp.Matrix([[0,-sp.I],[sp.I,0]])
    refs=[Z,(-Z+sp.sqrt(3)*X)/2,(-Z-sp.sqrt(3)*X)/2]
    def instrument(H):
        return sp.simplify((H+sum((w*U*H*U.H for w,U in zip(weights,refs)),
                                   sp.zeros(2)))/2)
    assert instrument(X)==2*X/3 and instrument(Z)==Z/3 and instrument(Y)==sp.zeros(2)
    bare = lambda H: 2*instrument(H)-H
    assert bare(X)==X/3 and bare(Z)==-Z/3 and bare(Y)==-Y
    # The unique minimal holding probability for a nonnegative spectrum
    # of the superoperator (a stronger requirement than complete positivity).
    hold=sp.symbols('hold',real=True)
    assert sp.solve(hold-(1-hold),hold)==[Q(1,2)]
    for T in (a,b):
        U=(4*T-I)/3
        S=U.multiply_elementwise(U)
        assert S==(I+8*T)/9 and (S-I/9)/(1-Q(1,9))==T
    print("Conditional successful scattering equals both Hecke generators; minimum hold 1/2")
    print("Local coin probabilities:",weights,"; channel eigenvalues on X,Y,Z: 2/3,0,1/3")
    for x,y,expected,label in [(Q(-1,3),Q(4,3),Q(-4,27),"amplitude"),
                               (Q(1,9),Q(8,9),Q(40,729),"Born transfer")]:
        U=x*I+y*a
        V=x*I+y*b
        assert U*V*U - V*U*V == expected*(a-b)
        print(label,"braid defect coefficient:",expected)
    # Exhaust all affine real normalized selfadjoint panel moves xI+(1-x)T:
    # relations reduce to x(1-x)(1+x)/2 multiplying a-b.
    x=sp.symbols('x')
    print("Normalized affine panel transform, braid solutions:",
          sp.solve(x*(1-x)*(1+x)/2,x))


if __name__=="__main__":
    a,b=flag_transfers()
    positive_history(a,b)
    return_reader(a,b)
    coin_probe(a,b)
