"""Operational readout and conditional projective control algebra."""
import sympy as sp


def simple(x):
    return sp.radsimp(sp.simplify(x))


def main():
    H=sp.diag(1,-1)
    a,b,d=sp.symbols('a b d',real=True)
    test=sp.Matrix([[a,b],[b,d]])
    vector=lambda m:sp.Matrix([m[0,0],m[0,1],m[1,0]])
    cols=[vector(M*H*M.inv()) for M in (sp.eye(2),test,test**2)]
    det=sp.factor(sp.Matrix.hstack(*cols).det())
    assert sp.factor(det-4*b**2*(a-d)*(a+d)/(a*d-b**2)**2)==0
    print('Lie-generator minor',det)
    s=sp.sqrt
    examples=[
        ('q=2', sp.Matrix([[-1+s(2),s(7)],[s(7),-1-s(2)]])/4,
         sp.diag(1,(s(57)+2*s(2))/7)),
        ('hex full quadratic',sp.Matrix([[1-2*s(3),s(15)],[s(15),-1-2*s(3)]]),
         sp.diag(1,(s(273)+4*s(3))/15))]
    for label,X,h in examples:
        Y=H*X*H
        assert simple(det.subs({a:X[0,0],b:X[0,1],d:X[1,1]}))!=0
        C=(X**2*Y**2-Y**2*X**2).applyfunc(sp.simplify)
        norm2=simple((C.T*C)[0,0])
        assert norm2>0 and (C**2+norm2*sp.eye(2)).applyfunc(simple)==sp.zeros(2)
        J=C/sp.sqrt(norm2)
        # Column-wise vectorization. Right J multiplication -> J^T tensor I.
        math_i=sp.kronecker_product(J.T,sp.eye(2))
        left_x=sp.kronecker_product(sp.eye(2),X)
        assert (math_i**2+sp.eye(4)).applyfunc(simple)==sp.zeros(4)
        assert math_i*left_x==left_x*math_i
        # The Born-effect operators, up to invertible common congruence by
        # h^(-1/2), are W^T h W. Include both double-record probes.
        raw=lambda M:M.T*h*M
        symvec=lambda M:sp.Matrix([M[0,0],M[1,1],M[0,1]])
        effects=(h,raw(X)-raw(Y),raw(X**2)+raw(Y**2))
        assert simple(sp.Matrix.hstack(*(symvec(e) for e in effects)).det())!=0
        assert simple((raw(X)-raw(Y))[0,1])!=0
        print(label,'raw Born effect rank 3; native imaginary direction on M2(R); Lie rank 3+3')
    # Determinant of Hermitian effect coordinate cone has Lorentz signature.
    t,x,y,z=sp.symbols('t x y z',real=True)
    Q=sp.Matrix([[t+z,x-sp.I*y],[x+sp.I*y,t-z]])
    assert sp.expand(Q.det())==t*t-x*x-y*y-z*z
    print('Hermitian cone determinant:',Q.det().expand())


if __name__=='__main__':
    main()
