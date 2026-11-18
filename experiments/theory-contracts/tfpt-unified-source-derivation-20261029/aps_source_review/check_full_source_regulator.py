"""Exact controls of full-source preservation; no new physical source model.

Reuse the existing original finite exterior algebra. Generic even densities
span every two-mode Grassmann functional, including its quartic interaction.
Bosonic Gaussian matrices test a universal regulator implication only.
"""
from pathlib import Path
import sys,json
import sympy as s
ROOT=Path(__file__).resolve().parents[4]
sys.path.insert(0,str(ROOT/'experiments/theory-contracts/seam-collar-multifield-20261016'))
from joint_source import add,clean,exponential,multiply,scale,variable
checks={}
def equal(p,q):return not clean({m:s.simplify(p.get(m,0)-q.get(m,0)) for m in set(p)|set(q)})
def derivative(p,j):
    return clean({m^(1<<j):(-1)**((m&((1<<j)-1)).bit_count())*v for m,v in p.items() if m&(1<<j)})
n=4;fields=[variable(n+j) for j in range(n)]
pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
aa=s.symbols('a0:6');bb=s.symbols('b0:6');t=s.symbols('t')
def quadratic(coeffs):
    out={}
    for a,(i,j) in zip(coeffs,pairs):out=add(out,scale(multiply(fields[i],fields[j]),a))
    return out
def differential(p,coeffs):
    out={}
    for a,(i,j) in zip(coeffs,pairs):out=add(out,scale(derivative(derivative(p,j),i),a))
    return out
def exp_differential(p,coeffs,sign=1):
    out=dict(p);term=dict(p)
    for k in range(1,n//2+1):
        term=scale(differential(term,coeffs),s.Rational(sign,k));out=add(out,term)
    return clean(out)
def integrate(p):
    mask=((1<<n)-1)<<n
    return clean({m^mask:v for m,v in p.items() if m&mask==mask})
def normalized(p):return scale(p,1/p[0])
rho={};cs=s.symbols('c0:8');even=[m for m in range(1<<n) if m.bit_count()%2==0]
for c,m in zip(cs,even):rho[m<<n]=c
source={}
for j in range(n):source=add(source,multiply(variable(j),fields[j]))
E=exponential(source,2*n)
z0=integrate(multiply(rho,E));V=quadratic(aa)
zr=integrate(multiply(multiply(rho,exponential(scale(V,-1),n)),E))
checks['generic_even_density_all8_coefficients_retained']=len(cs)==8 and all(any(p.has(c) for p in z0.values()) for c in cs)
checks['ordered_Grassmann_direct_integral_equals_source_differential']=equal(zr,exp_differential(z0,aa,-1))
checks['entire_quadratic_regulator_inverts_exactly']=equal(z0,exp_differential(zr,aa,1))
checks['two_regulator_stages_equal_sum']=equal(exp_differential(exp_differential(z0,aa,-1),bb,-1),exp_differential(z0,[a+b for a,b in zip(aa,bb)],-1))
zt=exp_differential(z0,[t*a for a in aa],-1)
checks['full_source_linear_RG_equation']=equal({m:s.diff(v,t) for m,v in zt.items()},scale(differential(zt,aa),-1))
zh=normalized(zt);Dh=differential(zh,aa)
rhs=add(scale(Dh,-1),scale(zh,Dh.get(0,0)))
checks['normalized_source_flow_keeps_vacuum_score']=equal({m:s.diff(v,t) for m,v in zh.items()},rhs)
checks['vacuum_log_flow_kept']=s.simplify(s.diff(zt[0],t)/zt[0]+Dh.get(0,0))==0
checks['normalization_requires_denominator']=equal(normalized(exp_differential(zr,aa,1)),normalized(z0))
def retained(p):return {m:v for m,v in p.items() if not m&12}
checks['discarding_exterior_source_before_unregulation_is_not_generic']=not equal(normalized(retained(z0)),normalized(retained(zr)))
# This is a test of an implication, not a TFPT measure or a new coupling.
C0=s.Matrix([[1,s.Rational(1,2)],[s.Rational(1,2),1]])
Rbad=s.ones(2);Rgood=s.Matrix([[2,-1],[-1,2]])
Cbad=(C0.inv()+Rbad).inv();Cgood=(C0.inv()+Rgood).inv()
P_adm=s.eye(2);theta=s.Matrix([[0,1],[1,0]])
checks['counterexample_starts_with_RP_gaussian']=C0.is_positive_definite and C0[0,1]>0
checks['counterexample_regulator_PSD']=Rbad.is_positive_semidefinite
checks['counterexample_structural_commutation']=Rbad*P_adm==P_adm*Rbad
checks['counterexample_reflection_covariance']=theta*Rbad*theta==Rbad
checks['counterexample_RP_cross_minus_one_over16']=Cbad[0,1]==-s.Rational(1,16)
checks['good_negative_cross_regulator_preserves_RP']=Cgood[0,1]>0
# Noncommuting same-side/cross matrices: necessity needs no joint diagonalization.
A=s.Matrix([[3,1],[1,2]]);B=s.diag(1,s.Rational(1,2));R=A.row_join(-B).col_join((-B).row_join(A))
cov=(s.eye(4)+R).inv()
checks['noncommuting_block_case']=A*B!=B*A
checks['good_full_block_PSD']=R.is_positive_semidefinite
checks['good_reflected_cross_covariance_PSD']=cov[:2,2:].is_positive_semidefinite
def mat(M):return [[str(x) for x in row] for row in M.tolist()]
checks={k:bool(v) for k,v in checks.items()}
out={'contract':'ORIGINAL-FULL-SOURCE-REGULATOR-20261030','all_passed':all(checks.values()),'checks':checks,
     'Grassmann_even_density_dimension':8,'Grassmann_quadratic_regulator_parameters':6,
     'generic_Z0_coefficients':{str(m):str(v) for m,v in z0.items()},
     'generic_regulated_scalar':str(zr.get(0,0)),
     'counterexample':{'C':mat(C0),'R':mat(Rbad),'P_adm':mat(P_adm),'regulated_C':mat(Cbad)},
     'good_regulator':{'R':mat(Rgood),'regulated_C':mat(Cgood)},
     'scope':'Generic finite exterior-algebra identities and finite tests of the analytic regulator cone. No original TFPT measure chosen; native68 states are not configuration fields. Whole RP/continuum/source origin not claimed.'}
Path(__file__).with_name('FULL-SOURCE-REGULATOR.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
if not out['all_passed']:raise RuntimeError('Full-source preservation check failed')
