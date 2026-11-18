"""Exact two-coefficient character-sector test for the original v117 ODE.

All computations are rational functions over Q(omega); no numeric ODE.
We use y=(h*p1/(a*b), h*z*p2/b, h*z*f), with a^2=8/9,b^2=5/9.
The third and second equations define p2,p1; the first is the residual.
"""
import sympy as S

w, om, c, d = S.symbols('w om c d')
mod = om**2+om+1

def red(expr):
    num, den = S.fraction(S.cancel(expr))
    # Work in Q(omega)(w): invert denominator modulo omega^2+omega+1.
    field = S.QQ.frac_field(w,c,d)
    np = S.Poly(num,om,domain=field).rem(S.Poly(mod,om,domain=field))
    dp = S.Poly(den,om,domain=field).rem(S.Poly(mod,om,domain=field))
    return S.cancel((np*S.invert(dp,S.Poly(mod,om,domain=field))).rem(
        S.Poly(mod,om,domain=field)).as_expr())

s=(1+w**3)/(1-w**3)  # z^2
l=3*w**2/(1-w**6)  # z'/z, derivative with respect to w
D=red(s*s-1)

def probe(H, exponent, label):
    # row3=h*z^exponent*f, row2=h*z^exponent*p2/b,
    # row1=h*z^(exponent-1)*p1/(a*b).
    k=red(S.diff(H,w)/(2*H))
    def ge(f,e):
        return red(e*f+(S.diff(f,w)+k*f)/l)
    f=c+d/H
    p2=red((D*ge(f,exponent)-s*s*f)/s)
    p1=red(D*ge(p2,exponent)-s*s*p2-S.Rational(5,9)*s*f)
    residual=red(ge(p1,exponent-1)-s*s/D*(2*p1+S.Rational(8,9)*p2))
    # c,d are constants in Q(omega), NOT restricted to rational/real values.
    # Do not split coefficients in omega before solving for them.
    eq=S.Poly(S.fraction(residual)[0],w).coeffs()
    candidate=None
    for equation in eq:
        ec=red(S.diff(equation,c))
        ed=red(S.diff(equation,d))
        if ed != 0:
            ratio=red(-ec/ed)
            if red(residual.subs({d:ratio*c})) == 0:
                candidate=ratio
            break
    sol=[(c,candidate*c)] if candidate is not None else [(S.Integer(0),S.Integer(0))]
    print(label,'solutions (c,d):',sol)
    for v in sol:
        if not any(v):
            continue
        sub={c:v[0],d:v[1]}
        print('p3=',red(f.subs(sub)))
        print('p2=',red(p2.subs(sub)))
        print('p1=',red(p1.subs(sub)))
        print('row1 residual=',red(residual.subs(sub)))
        if label.startswith('q '):
            ps=[red(p1.subs(sub).subs(c,1)),red(p2.subs(sub).subs(c,1)),red(f.subs(sub).subs(c,1))]
            p,q,r=ps
            direct=[red(D*ge(p,0)-s*s*(2*p+S.Rational(8,9)*q)),
                    red(D*ge(q,1)-p-s*s*q-S.Rational(5,9)*s*r),
                    red(D*ge(r,1)-s*q-s*s*r)]
            print('all original row residuals, scaled=',direct)
            P=S.Matrix([[red(t.subs(w,om**j*w)) for j in range(3)] for t in ps])
            detP=red(P.det())
            print('det of scaled column matrix P=',detP)
            print('det G /(z^4-1)=',red(s*detP/D)*S.Rational(27,10)/S.sqrt(2))
            kq=[-1,1,1]
            eps=[1,-1,-1]
            for j in range(3):
                aa=(-j)%3; pp=(-aa-1)%3; jj=(1-j)%3
                Hp=(w+om**pp)/(w+om**((pp+1)%3))
                fac=eps[aa]*om/Hp
                qrest=[]
                for row,t in enumerate(ps):
                    lhs=fac*t.subs(w,om**j/w)
                    rhs=kq[j]*t.subs(w,om**jj*w)
                    qrest.append(red(lhs-rhs if row<2 else lhs+rhs))
                print('Q column',j,'residuals=',qrest)
    return sol

probe((w+1)/(w+om),1,'q sector, row3=h*z*(c+d/H)')
probe((w-1)/(w-om),0,'existing p alternative, row3=h*(c+d/H)')
