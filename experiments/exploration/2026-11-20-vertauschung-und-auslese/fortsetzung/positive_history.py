"""Finite matrices in a full, infinite-group mirror-positive braid kernel."""
from itertools import product
import sympy as sp

I=sp.eye(2)
P=sp.diag(1,0)
Q=sp.Matrix([[2,sp.sqrt(14)],[sp.sqrt(14),7]])/9
r=sp.Rational(1,6)
X=(-2*I+3*P)/2
Y=(-2*I+3*Q)/2


def inner(U,V):
    return sp.simplify(sp.trace(U.T*V)/2)


def main():
    assert Q*Q==Q and P*Q*P==sp.Rational(2,9)*P
    assert X*Y*X==Y*X*Y and X*Y!=Y*X
    assert X.eigenvals()==Y.eigenvals()=={sp.Rational(1,2):1,-1:1}
    assert (X*Y)**3==I/8
    print('X,Y eigenspectra',X.eigenvals(),'; (XY)^3 =',((X*Y)**3).tolist())
    print('sum of off-identity coefficient majorants <=',4*r/(1-2*r))
    products={'':I}
    for size in range(1,6):
        for letters in product('ab',repeat=size):
            word=''.join(letters)
            products[word]=products[word[:-1]]*({'a':X,'b':Y}[word[-1]])
    for u,v in [('aa','bb'),('aabb','bbaa')]:
        G=sp.Matrix([[inner(products[i],products[j]) for j in (u,v)] for i in (u,v)])
        contrast=sp.simplify((G[0,0]+G[1,1]-2*G[0,1])/(G[0,0]+G[1,1]))
        assert G.det()>0
        print('words',u,v,'Gram before common damping r^'+str(2*len(u)),G.tolist(),
              'normalized contrast',contrast)
    # Complete OS construction gives Gram matrices as vectors in real M_2
    # with the normalized Frobenius inner product. Compute exact rank for
    # low-degree words by taking the span of their four matrix entries.
    for n in range(4):
        vectors=[sp.Matrix(list(t)) for w,t in products.items() if len(w)<=n]
        print('positive histories through length',n,'OS rank',sp.Matrix.hstack(*vectors).rank())
    assert sp.Matrix.hstack(*(sp.Matrix(list(products[w])) for w in ('','a','b','ab'))).rank()==4
    # The reflected Hilbert space can be minimized by the exchange-symmetric
    # vector on the bisector of the two rank-one rays. No scalar transfer
    # representation can distinguish words differing only in order.
    left=sp.Matrix([1,0])
    right=sp.Matrix([sp.sqrt(2)/3,sp.sqrt(7)/3])
    v=left+right
    swap=2*(v*v.T)/(v.T*v)[0]-I
    assert (swap*P*swap-Q).applyfunc(sp.simplify)==sp.zeros(2)
    assert (swap*v-v).applyfunc(sp.simplify)==sp.zeros(2,1)
    assert sp.Matrix.hstack(v,X*v).rank()==2
    comm=(X**2*Y**2-Y**2*X**2)
    assert comm.T*comm==sp.Rational(7,128)*I
    print('exchange-symmetric pure boundary vector: minimal OS rank 2;'
          ' mixed-history squared distance before r^8: 7/128')
    # Square-step transfers are strictly positive. Original steps have a
    # negative eigenvalue and do not equal an exponential of a selfadjoint H.
    assert (X**2).is_positive_definite and (Y**2).is_positive_definite and X**2*Y**2!=Y**2*X**2

    # Hodge closure for two measured primitive directions, normalized trace on
    # their active span. The tensor-product norm fixes the scalar/quadrupole
    # mixture. On the *whole* symmetric power that is a different state.
    hexagon=sp.Matrix([1,sp.sqrt(3)])/2
    primitive=sp.Matrix([1,0])
    single=(primitive.T*hexagon)[0]**2
    assert single==sp.Rational(1,4)
    t=single**2
    root=7-4*sp.sqrt(3)  # root < 1 chosen to normalize the complementary mode
    assert sp.simplify(root/(1+root)**2-t)==0
    qform=(1-root**2)**2*(1-t)/(1+root**4)
    ordered=(1-root**2)**4*t*(1-t)/(2*root**4+(1-root**4)**2*t)
    assert sp.simplify(qform-sp.Rational(90,97))==0
    assert sp.simplify(ordered-sp.Rational(1080,1177))==0
    # Traceless reduction changes the relative weight and restores q=1.
    D=lambda v: v*v.T-I/2
    assert sp.trace(D(primitive)*D(hexagon))**2/(sp.trace(D(primitive)**2)
           *sp.trace(D(hexagon)**2))==sp.Rational(1,4)
    print('hex H1 / tensor-square / traceless pair overlaps:',single,t,single)
    print('tensor-square active span: q + 1/q = 14; contrasts aa/bb, mixed order',
          sp.simplify(qform),sp.simplify(ordered))


if __name__=='__main__':
    main()
