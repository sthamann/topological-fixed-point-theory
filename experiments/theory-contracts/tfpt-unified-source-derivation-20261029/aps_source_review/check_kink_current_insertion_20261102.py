"""One-mode independent check of the source Weyl and current order.

The one-quantum state below diagnoses the CCR identity; it is NOT an
identification of that state with any original TFPT quarter source state.
"""
import sympy as s
m=s.symbols('m',positive=True,real=True)
h,eps=s.symbols('h eps',real=True)
d,dc=s.symbols('d dc') # independent conjugate coordinates; dc=conj(d)

def F(x,xc):
    y=m*x*xc
    return s.exp(-y/2)*(1-y) # number state |1>

# Source current J(h)=sqrt(m)(h a* + h a) for real h, [a,a*]=1.
# W(eps h) W(d)=exp[-eps m h(d-dc)/2] W(d+eps h).
generated=s.diff(s.exp(-eps*m*h*(d-dc)/2)*F(d+eps*h,dc+eps*h),eps).subs(eps,0)/s.I
# Independent Fock calculation: alpha=i sqrt(m)d, <1|a*D(alpha)|1>
# = -alpha* exp(-|alpha|²/2), <1|a D(alpha)|1>
# = alpha (2-|alpha|²) exp(-|alpha|²/2).
operator=s.I*m*h*s.exp(-m*d*dc/2)*(dc+d*(2-m*d*dc))
assert s.simplify(generated-operator)==0
value=s.simplify(operator.subs({m:8,h:1,d:1+s.I,dc:1-s.I}))
assert s.simplify(value-(120-104*s.I)*s.exp(-8))==0
# Reversing current and Weyl word changes the answer by the CCR commutator.
reverse=s.diff(s.exp(+eps*m*h*(d-dc)/2)*F(d+eps*h,dc+eps*h),eps).subs(eps,0)/s.I
assert s.simplify(generated-reverse-s.I*m*h*(d-dc)*F(d,dc))==0
print('PASS one-mode Fock versus Weyl derivative; current order matters')
print('m=8, h=1, d=1+i:',value)
