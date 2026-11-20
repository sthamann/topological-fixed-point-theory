"""Stationary two-result instrument from the two historical transfers."""
import sympy as sp


def simple(v):
    return sp.radsimp(sp.simplify(v))


def main():
    s=sp.sqrt
    X=sp.Matrix([[-1+s(2),s(7)],[s(7),-1-s(2)]])/4
    switch=sp.diag(1,-1)
    Y=switch*X*switch
    assert (X*Y*X-Y*X*Y).applyfunc(sp.simplify)==sp.zeros(2)
    L=lambda M: (X*M*X+Y*M*Y)/2
    basis=[sp.Matrix(2,2,[int(k==j) for k in range(4)]) for j in range(4)]
    Lmat=sp.Matrix.hstack(*(sp.Matrix(list(L(M))) for M in basis))
    lam=(3+s(57))/16
    assert Lmat.eigenvals()=={lam:1,(3-s(57))/16:1,sp.S(3)/8:1,-sp.S.Half:1}
    R=(s(57)+2*s(2))/7
    h=sp.diag(1,R)
    assert (L(h)-lam*h).applyfunc(sp.simplify)==sp.zeros(2)
    print('spectrum counted CP transfer:',Lmat.eigenvals())
    print('stationary density ratio after Perron normalization:',R**2)
    print('noncentral modular angular gap in dimensionless units: 2 log((sqrt(57)+2 sqrt(2))/7)')
    root=sp.diag(1,s(R))
    kraus=[root*g*root.inv()/s(2*lam) for g in (X,Y)]
    norm=lambda M: M.applyfunc(sp.simplify)
    assert norm(sum((k.T*k for k in kraus),sp.zeros(2))-sp.eye(2))==sp.zeros(2)
    assert norm(sum((k*h**2*k.T for k in kraus),sp.zeros(2))-h**2)==sp.zeros(2)
    # Born probabilities in chronological order, evaluated without square roots of h.
    def probability(word):
        M=sp.eye(2)
        for g in word:
            M={'a':X,'b':Y}[g]*M
        return simple(sp.trace(h*M*h*M.T)/((2*lam)**len(word)*sp.trace(h**2)))
    assert probability('a')==sp.Rational(1,2)
    p=probability('aabb')
    assert simple(p-probability('bbaa'))==0
    f=simple(probability('aabba')/p)
    g=simple(probability('bbaaa')/p)
    print('equal P(aabb)=P(bbaa)=',p)
    print('P(next a | aabb), P(next a | bbaa) =',f,g)
    print('difference:',simple(f-g),'approx',float(f-g))
    assert simple(f+g)==1 and simple(f-g)!=0
    # General q >= 1 normalized as p=1/q. In the bisector basis let
    # c=sqrt(p)/(1+p). The diagonal sector of L is a symmetric positive
    # 2x2 matrix. Its Perron eigenvector ratio is sqrt(1+b^2)+b.
    phex=7-4*s(3)
    b2=simple(4*(1-phex)**2*phex/(1+phex+phex**2)**2)
    assert b2==sp.Rational(16,75)
    Rhex=(s(273)+4*s(3))/15
    assert simple(Rhex-sp.sqrt(1+b2)-sp.sqrt(b2))==0
    print('hex tensor state: stationary spectral ratio =',Rhex**2)
    # Overall step normalization cancels in instrument probabilities.
    # This version of X,Y uses only square roots of rational integers.
    X=sp.Matrix([[1-2*s(3),s(15)],[s(15),-1-2*s(3)]])
    Y=switch*X*switch
    h=sp.diag(1,Rhex)
    lam=13+s(273)
    assert (L(h)-lam*h).applyfunc(sp.simplify)==sp.zeros(2)
    same=probability('aabb')
    assert simple(same-probability('bbaa'))==0
    f=simple(probability('aabba')/same)
    g=simple(probability('bbaaa')/same)
    assert simple(f+g)==1
    print('hex geometry P(next a | aabb), P(next a | bbaa):',float(f),float(g))
    print('hex geometric signed difference:',simple(f-g))


if __name__=='__main__':
    main()
