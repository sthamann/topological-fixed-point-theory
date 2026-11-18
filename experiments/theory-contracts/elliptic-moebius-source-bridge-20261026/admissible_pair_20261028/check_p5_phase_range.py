#!/usr/bin/env python3
"""Exact P5 interval on the previously defined SM phase family, not a P1 source."""
from pathlib import Path
import importlib.util
import json
import sympy as s

HERE = Path(__file__).resolve().parent
source = HERE.parent / 'hodge_complement_20261028' / 'hodge_complement_channel.py'
spec = importlib.util.spec_from_file_location('original_hodge', source)
h = importlib.util.module_from_spec(spec)
spec.loader.exec_module(h)
C = [s.Matrix(x.tolist()) for x in h.hodge_amplitudes()]
pairs = h.PAIRS

# Formal phases and their conjugates make the full result a polynomial identity.
a,b,c,A,B,Cc = s.symbols('a b c A B C')
z = [a if max(p)<3 else c if min(p)>=3 else b for p in pairs]
zc = [A if max(p)<3 else Cc if min(p)>=3 else B for p in pairs]
U, Uc = s.diag(*z),s.diag(*zc)
assert all((C[i].T*C[j]).trace()==6*(i==j) for i in range(5) for j in range(5))
# 1/18 = (1/6 for P5) * (1/3 for Kraus-normalized Choi).
response = s.expand(sum((C[i].T*C[j]*U).trace()*(C[i].T*C[j]*Uc).trace()
                        for i in range(5) for j in range(5))/18)
formula = ((a+4*b+c)*(A+4*B+Cc)/6 + (a+b)*(A+B))
assert s.expand(response-formula)==0
# Whole phase family preserves every individual event effect exactly.
for x in C:
    assert (Uc*x.T*x*U - x.T*x).subs({a*A:1,b*B:1,c*Cc:1}) == s.zeros(10)
# Actual normalization and all ten diagonal transition probabilities survive.
for v in [s.eye(10),s.diag(*[-1 if max(p)<3 else 1 for p in pairs]),
          s.diag(*[-1 if max(p)<3 or min(p)>=3 else 1 for p in pairs])]:
    assert sum((v.T*x.T*x*v for x in C),s.zeros(10))==3*s.eye(10)
    for x in C:
        for i in range(10):
            for j in range(10):
                assert (x*v)[i,j]**2 == x[i,j]**2
assert response.subs({a:1,b:1,c:1,A:1,B:1,Cc:1})==10
assert response.subs({a:-1,b:1,c:-1,A:-1,B:1,Cc:-1})==s.Rational(2,3)
assert response.subs({a:-1,b:1,c:1,A:-1,B:1,Cc:1})==s.Rational(8,3)
result = {
 'carrier': 'Lambda2(C5) in original PAIRS basis; original E8 Cm and P5',
 'family': 'Phi o Ad(a P_CC + b P_CW + c P_WW), |a|=|b|=|c|=1',
 'objective': 'Tr(P5 J) = |a+4b+c|^2/6 + |a+b|^2',
 'normalization': 'unnormalized Choi trace 10; divide objective by 10 for normalized probability',
 'exact_min': '2/3', 'min_phases': ['-1','1','-1'],
 'exact_max':'10', 'max_phases':['1','1','1'],
 'proof':'Gauge b=1. Put t=Re(a), s=|a+4|=sqrt(17+8t), 3<=s<=5. Minimizing over c gives (s-1)^2/6+(s^2-9)/4, strictly increasing for s>=3; hence 2/3 at s=3. Triangle bounds give 10; equal phases attain it. Continuous domain gives full interval.',
 'legacy_one_sided_control':'8/3 (normalized probability 4/15)',
 'coverage':'only this phase family, not all SM-covariant channels and not all processes',
 'source_bridge':'P5 embedding into local P1 A_adm and operational sewing rule remain unproved; none of these phase channels is thereby selected',
 'clock':'original U5 is scalar on Lambda2 C5 and cannot distinguish this family',
}
(HERE/'p5-phase-range-results.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
# Exact full-Petersen-population search, no SM covariance imposed.
# Every original Cm is a signed partial permutation with 6 nonzero entries.
# Flip exactly three; no cross-Cm overlap exists.
D=[]
for x in C:
    y=x.copy()
    nonzero=[(i,j) for i in range(10) for j in range(10) if x[i,j]!=0]
    assert len(nonzero)==6
    for i,j in nonzero[:3]: y[i,j] = -y[i,j]
    assert y.T*y == x.T*x
    D.append(y)
assert all((x.T*y).trace()==0 for x in C for y in D)
assert sum((y.T*y for y in D),s.zeros(10)) == 3*s.eye(10)
assert sum((y*y.T for y in D),s.zeros(10)) == 3*s.eye(10)
# Same full conditional populations, not only aggregated counts.
assert all(sum(y[i,j]**2 for y in D)==sum(x[i,j]**2 for x in C)
           for i in range(10) for j in range(10))
assert sum((x.T*y).trace()**2 for x in C for y in D)==0
# Positivity and trace 10 force 0 <= Tr(P5 J) <= 10 for any CPTP J.
# Convex mixtures between D and C achieve all values of that interval.
result['all_petersen_cptp']={
 'exact_min':'0','exact_max':'10',
 'construction':'Dm flips 3 of the 6 nonzero entries of each original Cm, with all other entries unchanged. J0 = (1/3) sum |Dm>><<Dm|.',
 'controls':'Exact orthogonality to all Cm, every individual Dm effect = Cm effect, TP, unital, every one of 100 populations unchanged.',
 'interval_proof':'0<=P5<=I and Tr(J)=10. Endpoints J0 and Jphi meet bounds; all convex mixtures preserve the full population and yield the interval.',
 'scope':'All CPTP channels with full fixed Petersen population; the exhibited zero-construction is not SM-covariant.'}
# Match Psi P5 value 2 with a different process in original SM phase family.
w=(1+2*s.sqrt(2)*s.I)/3
S=s.diag(*[-1 if max(p)<3 else w if min(p)>=3 else 1 for p in pairs])
Sc=s.conjugate(S)
assert s.simplify(sum((x.T*y*S).trace()*(x.T*y*Sc).trace() for x in C for y in C)/18)==2
# Swap-odd projector P45=(I-F)/2; obtain a second independent Choi response.
# P45 weight = (1/6) sum ||xS-(xS)^T||_HS^2 / 2.
anti=s.simplify(sum(sum(s.conjugate(q)*q for q in (x*S-(x*S).T)) for x in C)/12)
result['same_p5_as_psi']={'phases':['-1','1','(1+2i sqrt(2))/3'], 'P5': '2', 'P45':str(anti), 'psi_P45':'4 (prior exact six-mark original-source check)'}
assert anti != 4
(HERE/'p5-phase-range-results.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'all_petersen_cptp':result['all_petersen_cptp'],'same_p5_as_psi':result['same_p5_as_psi']},indent=2))
# Distinguish the unrestricted zero construction from SM-covariant channels.
# Under the actual real SU(3) generator E01-E10, dual acts by the same real
# skew matrix in this fixed chart. We use unnormalized M0=3 J0.
X=s.zeros(5); X[0,1]=1; X[1,0]=-1
r=s.Matrix(h.exterior_two_differential(__import__('numpy').array(X.tolist(),dtype=int)).tolist())
L=s.kronecker_product(r,s.eye(10))+s.kronecker_product(s.eye(10),r)
vec=lambda x:s.Matrix([x[i,j] for j in range(10) for i in range(10)])
M0=sum((vec(y)*vec(y).T for y in D),s.zeros(100))
comm=L*M0-M0*L
comm_norm2=sum(q*q for q in comm)
assert comm_norm2==576
result['all_petersen_cptp']['zero_candidate_sm_covariance']='fails: SU3 E01-E10 commutator of M0=3J0 has exact HS norm squared 576'
(HERE/'p5-phase-range-results.json').write_text(json.dumps(result,indent=2)+'\n')
print('Exact SU3 commutator norm^2 of M0=3J0:',comm_norm2)
