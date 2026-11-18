"""Test minimal central A3 winding as the native relative pair operation.

Conditional identification: P1 winding is NOT equated with A3 by this file.
"""
import itertools
import json
from pathlib import Path
import sys
import sympy as s
sys.path.insert(0,'/Users/stefanhamann/Projekte/tfpt-theoryv4')
from tfpt_explorer.neutral_source_response import _native_neutral_branch

n=_native_neutral_branch()
i=s.I
U=s.diag(1,i,-i)
M=s.Matrix([[0,-(1+i)/2,(1-i)/2],[-(1+i)/2,-i/2,-s.Rational(1,2)],[(1-i)/2,-s.Rational(1,2),i/2]])
U4=s.diag(1,U);M4=s.diag(1,M)
assert M4**3==s.eye(4)
x=s.symbols('x:16')
X=s.Matrix(4,4,x)
equations=list(X*U4-U4*X)+list(X*M4-M4*X)+[s.trace(X)]
solution=s.linsolve(equations,x)
row=next(iter(solution))
free=set().union(*(e.free_symbols for e in row))
assert len(free)==1
H=s.diag(-3,1,1,1)
assert X.subs(dict(zip(x,row)))==next(iter(free))*H
F=H/4
assert all(s.exp(2*s.pi*i*d)==i for d in F.diagonal())
assert s.trace(F)==0 and s.trace(F**2)==s.Rational(3,4)
# Proof over all integer lifts: eigenvalues =1/4+m_j, sum m_j=-1.
# Sum f_j²=sum m_j²-1/4 >=1-1/4. Equality requires exactly one m=-1,
# three m=0. Commuting with 1+irreducible3 chooses the fixed line.
assert s.Rational(1,4)*4-1==0
minima=[]
for shifts in itertools.product(range(-2,3),repeat=4):
    if sum(shifts)==-1:
        norm=sum((s.Rational(1,4)+m)**2 for m in shifts)
        if norm==s.Rational(3,4):minima.append(shifts)
assert len(minima)==4
assert all(tuple(sorted(m))==(-1,0,0,0) for m in minima)
weighted=s.zeros(8)
family=n['family']
v=-s.Matrix([0]*5+list(family[0]))/2
weights=[]
for r in n['spinor']:
    b=s.Matrix(r)/2
    a=family.index(tuple(r[5:]))
    w=(b.T*v)[0]
    assert w==F[a,a]
    weighted+=w*b*b.T
    weights.append(w)
assert weighted==n['R_matrix']
assert weighted==n['A48']-s.Rational(3,4)*n['A64']
assert s.trace(weighted)==0
assert s.trace(weighted**2)/2==48
# A sharper conditional route preserves the INTEGER P1 winding before
# taking the traceless A3 part. On the specified family C4, any commuting
# Hermitian one-parameter loop generator is diag(a,b,b,b). Closed loop
# means integer a,b; rank1 and determinant winding1 force a=1,b=0.
# This does NOT identify the original seam mode with this family C4.
P=s.diag(1,0,0,0)
assert s.trace(P)==1 and P.rank()==1
assert P*U4==U4*P and P*M4==M4*P
assert P-s.eye(4)/4==-F
theta=s.symbols('theta',real=True)
loop=s.diag(s.exp(2*s.pi*i*theta),1,1,1)
assert s.simplify(loop.det()-s.exp(2*s.pi*i*theta))==0
assert loop.subs(theta,1)==s.eye(4)
result={
 'relative_operation':'F=diag(-3,1,1,1)/4',
 'clock':'fixed line + actual original U/M three-dimensional representation',
 'traceless_clock_commutant_dimension':1,
 'central_holonomy':'exp(2πi F)=i I4',
 'minimal_coweight_norm_squared':'3/4',
 'infinite_integer_lift_proof':'sum m_j=-1 implies sum (m_j+1/4)^2=sum m_j^2-1/4 >=3/4; equality iff one m=-1, others0',
 'minimal_line_selection':'Clock fixes omitted line and acts irreducibly on retained3',
 'root_pair_weights':{'-3/4':weights.count(-s.Rational(3,4)),'1/4':weights.count(s.Rational(1,4))},
 'native_pair_response':'B48-(3/4)B64=R; Sugawara stress cancels automatically',
 'native_response_norm_squared':'48',
 'integer_winding_refinement':{
    'declared_class':'Hermitian constant generator of a closed one-parameter loop on the native family C4, commuting with U4/M4',
    'general_generator':'diag(a,b,b,b), a,b integer',
    'unit_winding_condition':'a+3*b=1',
    'rank_one_primitive_condition':'rank H=1 forces b=0, a=1',
    'unique_generator':'P_fixed=diag(1,0,0,0)',
    'loop':'diag(exp(2πit),1,1,1), determinant winding1',
    'traceless_part':'P_fixed-I4/4=-F',
    'relative_pair_response':'-R (orientation reversed gives R)',
    'first_added_identification':'the P1 primitive seam loop acts in this particular native family C4, through this constant one-parameter loop class; the source Hilbert space is larger',
 },
 'not_identified':['P1 geometric winding with this A3 central holonomy','Gaussian quarterturn with central iI','P1 primitive Hodge projection with relative pair observable','physical coupling normalization','physical alpha sewing time','3+1D source'],
 'verdict':'conditional exact origin candidate for R: minimal central winding + original clock; raw-seam link unproved'
}
Path(__file__).with_name('native_relative_winding_results.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(result,ensure_ascii=False,indent=2))
