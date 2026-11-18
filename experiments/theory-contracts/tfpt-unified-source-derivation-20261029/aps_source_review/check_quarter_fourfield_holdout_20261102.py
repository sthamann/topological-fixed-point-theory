#!/usr/bin/env python3
"""Exact CAR check on the existing R062 right-up colour triplet."""
import sympy as s

t, z = s.symbols('t z')
I = s.I
zero = s.zeros(8)
eye = s.eye(8)
ann = []
for mode in range(3):
    a = s.zeros(8)
    for state in range(8):
        if state & (1 << mode):
            sign = (-1) ** ((state & ((1 << mode)-1)).bit_count())
            a[state ^ (1 << mode), state] = sign
    ann.append(a)
for i in range(3):
    for j in range(3):
        assert ann[i]*ann[j].H + ann[j].H*ann[i] == (eye if i == j else zero)
N = sum((a.H*a for a in ann), zero)
O = 4*N*(N-eye)
Q = s.diag(*(z**state.bit_count() for state in range(8)))
rho = s.diag(*(t**state.bit_count()*(1-t)**(3-state.bit_count()) for state in range(8)))
assert Q.subs(z,I)*ann[0].H == I*ann[0].H*Q.subs(z,I)
direct = s.expand(s.trace(rho*Q*O))
assert s.simplify(direct-24*z**2*t**2*(1-t+z*t)) == 0
ordered = s.expand(s.trace(rho*ann[0]*Q*O*ann[0].H))
expected = 16*t*(1-t)**2*z**2 + 24*t**2*(1-t)*z**3
assert s.simplify(ordered-expected) == 0
ts = [s.Rational(i,100) for i in (52,54,56,58)]
ws = [[s.Rational(i,16) for i in weights] for weights in ((5,1,7,3),(3,7,1,5))]
moments = [[sum(w*x**j for w,x in zip(weights,ts)) for j in range(4)] for weights in ws]
assert moments[0][:3] == moments[1][:3] and moments[0][3] != moments[1][3]
values = [s.simplify(sum(w*ordered.subs({t:x,z:I}) for w,x in zip(weights,ts))) for weights in ws]
assert values[0] != values[1]
assert s.Rational(58,100) < (1+s.tanh(s.pi/2))/2
print('CAR, clock action, fourfield and ordered identities: PASS')
print('moments:', moments)
print('ordered quarter holdout:', values)
print('difference:', s.simplify(values[0]-values[1]))

# Same conditional source with all three coupled left partners retained.
# Bit pairs are (right,left); the displayed coherent 2-mode covariance has
# t=3/5, det(C)=1/5 and off-diagonal 1/5, compatible with p=1/4 and finite ell.
T, D = s.Rational(3,5), s.Rational(1,5)
rho2=s.Matrix([[D,0,0,0],[0,T-D,s.Rational(1,5),0],
               [0,s.Rational(1,5),1-T-D,0],[0,0,0,D]])
rho6=s.kronecker_product(rho2,rho2,rho2)
a=s.zeros(64)
for state in range(64):
    if state & 1:
        a[state ^ 1,state]=1
NR=s.diag(*(sum((state>>j)&1 for j in (0,2,4)) for state in range(64)))
O6=4*NR*(NR-s.eye(64))
Q6=s.diag(*(z**state.bit_count() for state in range(64)))
w0=D+z*(1-T-D)
w1=T-D+z*D
actual=s.expand(s.trace(rho6*a*Q6*O6*a.H))
assert s.simplify(actual-(16*z**2*w0**2*w1+24*z**3*w0*w1**2))==0
print('coupled six-mode complement, global quarter word: PASS', s.expand(actual.subs(z,I)))
