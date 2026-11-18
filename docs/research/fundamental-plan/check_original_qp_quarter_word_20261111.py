"""Actual C3 quarter and hexagonal q/p side transport on original genus two."""
import sympy as s
w,e=s.symbols('w eta',nonzero=True)
rho=(1+s.I*s.sqrt(3))/2
I=s.I
def simple(z):return s.simplify(s.expand_complex(z)) if not z.free_symbols else s.simplify(s.expand(z))
def same(a,b):return all(simple(x-y)==0 for x,y in zip(a,b))
def compose(a,b):
    return tuple(simple(x.subs({w:b[0],e:b[1]},simultaneous=True)) for x in a)
identity=(w,e)
def power(a,n):
    r=identity
    for _ in range(n):r=compose(a,r)
    return r
Q=(1/w,-I*e/w**3)
R=(rho*w,e)
H=(w,-e)
z=lambda p:s.simplify(p[1]/(1-p[0]**3))
assert same(power(Q,2),H) and same(power(Q,4),identity)
assert same(power(R,6),identity)
assert s.cancel(z(Q)-I*z(identity))==0
assert s.cancel(s.cancel(z(R)*z(identity)-1).subs(e**2,1-w**6))==0
for n in (1,3,5):
    r=power(R,n)
    rel=compose(compose(compose(Q,r),power(Q,3)),power(R,6-n))
    predicted=compose(H,power(R,(-2*n)%6))
    assert same(rel,predicted) and s.cancel(z(rel)+z(identity))==0
    print(f'Original odd side step R^{n}: [Q,R^{n}] = H R^{(-2*n)%6}, z -> -z')
S=power(R,3)
closed=compose(compose(compose(compose(H,Q),S),power(Q,3)),S)
assert same(closed,identity)
print('For pure side exchange S=R^3, [Q,S]=H=Q^2. Inserting H closes the base; the spin/source lift of the corrected word is a separate test.')

# S0 spin line: Phi^2=(1+w^3) dw/eta, u=(1+w^3)/eta, Chi=Phi/u.
# Work with pullback operators on f(w,u) Phi, in this order.
u,a,b=s.symbols('u a b',nonzero=True)
def subst2(expr,p):return s.cancel(expr.subs({w:p[0],u:p[1]},simultaneous=True))
def lc(A,B):
    Ap,Am=A; Bp,Bm=B
    return (tuple(subst2(x,Ap) for x in Bp),s.cancel(Am*subst2(Bm,Ap)))
def eqspin(A,B):return all(s.cancel(x-y)==0 for x,y in zip(A[0],B[0])) and s.cancel(A[1]-B[1])==0
assert s.cancel((1+Q[0]**3)*s.diff(Q[0],w)/Q[1]/((1+w**3)/e)+I/w**2)==0
assert s.cancel((1+S[0]**3)*s.diff(S[0],w)/S[1]/((1+w**3)/e)+(1-w**3)/(1+w**3))==0
assert s.cancel((1+Q[0]**3)/Q[1]-I*(1+w**3)/e)==0
assert s.cancel(s.cancel(((1+S[0]**3)/S[1])*((1+w**3)/e)-1).subs(e**2,1-w**6))==0
lq=((1/w,I*u),a/w)
ls=((-w,1/u),b/u)
lH=lc(lq,lq)
# a^2=-i, b^2=-1 from the square. signs of a,b are immaterial.
lH=(lH[0],s.cancel(lH[1].subs(a**2,-I)))
assert lH==((w,-u),-I)
assert eqspin(lc(lq,ls),lc(lH,lc(ls,lq)))
assert eqspin(lc(lH,lH),((w,u),-1))
print('Actual S0 spin pullbacks: [Qhat,Shat]=Qhat^2 with the Q^2-induced half-turn.')
print('Correcting with (Qhat^2)^-1 yields +I. The other lift of the same base half-turn yields -I.')
