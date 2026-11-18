"""Original real Higgs one-form and complete finite boson-source sewing.

Y coefficients label supplied matrices; none are fitted or interpreted as new
physical inputs. Boson moments remain those of the original unchosen source.
"""
from pathlib import Path
import sys, json
import sympy as s
import numpy as np
ROOT=Path(__file__).resolve().parents[4]
sys.path[:0]=[str(ROOT/'verification'),str(ROOT/'experiments/theory-contracts/seam-collar-multifield-20261016')]
import v252_full_finite_triple as original
from joint_source import add, scale, multiply, exponential, variable, clean
checks={}
def need(name, value):
    checks[name]=bool(value)
    if not value: raise RuntimeError(name)
def eq(a,b):
    return not clean({m:s.expand(a.get(m,0)-b.get(m,0)) for m in set(a)|set(b)})

# Four REAL coordinates of the original H=(h1,h2), no 1/sqrt(2) inserted.
Q=[np.array([[0,1],[-1,0]],complex),np.array([[0,1j],[1j,0]],complex),np.eye(2,dtype=complex),np.diag([-1j,1j])]
zero3=np.zeros((3,3),complex); zero2=np.zeros((2,2),complex)
bp=original.rep(1,zero2,zero3)
K=original.Kswap(); gamma=original.grading()
I=original.IDX
def expected(Qa,Ys):
    S=np.zeros((48,48),complex)
    for up,down,nc,Yup,Ydown in [('nu','e',1,Ys[3],Ys[2]),('u','d',3,Ys[0],Ys[1])]:
      for color in range(nc):
       for g in range(3):
        for k in range(3):
         for row,left in enumerate([up,down]):
          for col,right,Y in [(0,up,Yup),(1,down,Ydown)]:
           i=I[(left,'L',color,g)]; j=I[(right,'R',color,k)]
           S[i,j]=Qa[row,col]*Y[g,k];S[j,i]=np.conjugate(S[i,j])
    out=np.zeros((96,96),complex);out[:48,:48]=S;out[48:,48:]=np.conjugate(S)
    return out
count=0
for species in range(4):
 for g in range(3):
  for k in range(3):
   for phase in [1,1j]:
    Ys=[zero3.copy() for _ in range(4)];Ys[species][g,k]=phase
    D=original.full_D(*Ys,zero3)
    for Qa in Q:
     a=original.rep(0,Qa,zero3)
     om=a@(D@bp-bp@D)
     A=om+om.conjugate().T
     T=A+K@A.conjugate()@K
     if not np.array_equal(T,expected(Qa,Ys)):raise RuntimeError('actual oneform basis mismatch')
     if not np.array_equal(T,T.conjugate().T):raise RuntimeError('Hermitian')
     if np.any(T@gamma+gamma@T):raise RuntimeError('odd')
     if np.any(K@T.conjugate()@K-T):raise RuntimeError('KO reality')
     count+=1
need('all288_original_real_basis_images_equal_whole96_Higgs_derivatives',count==288)
need('all288_images_Hermitian_odd_and_JF_real',True)

# Generic quaternion factor, entire family matrices follow entrywise.
x=s.symbols('x1:5',real=True)
h1=x[0]+s.I*x[1];h2=x[2]+s.I*x[3]
qh=s.Matrix([[s.conjugate(h2),h1],[-s.conjugate(h1),h2]])
need('four_real_components_reconstruct_original_quaternion',qh==sum((x[i]*s.Matrix(Q[i]).applyfunc(s.nsimplify) for i in range(4)),s.zeros(2)))
need('quaternion_norm_keeps_whole_H_direction',s.simplify(qh.conjugate().T*qh-(sum(z*z for z in x))*s.eye(2))==s.zeros(2))

# Select four disjoint ORIGINAL directed bilinears (nu,e,u,d) solely to
# extract their coefficients. Other original legs are not declared absent.
y=s.symbols('Ynu Ye Yu Yd')
j={}
for k in range(4):j=add(j,scale(multiply(variable(2*k),variable(2*k+1)),y[k]))
H=s.symbols('h',real=True);mu=s.symbols('mu0:5');mom={0:s.Integer(1),1:mu[1],2:mu[2],3:mu[3],4:mu[4]}
def average(poly):
 out={}
 for mask,value in poly.items():
  p=s.Poly(s.expand(value),H)
  out[mask]=s.expand(sum(co*mom[power[0]] for power,co in p.terms()))
 return clean(out)
F=average(exponential(scale(j,-H),8))
u=add(F,{0:-s.Integer(1)})
logF={};power={0:s.Integer(1)}
for k in range(1,5):
 power=multiply(power,u);logF=add(logF,scale(power,s.Rational((-1)**(k+1),k)))
kappa=[None,mu[1],mu[2]-mu[1]**2,mu[3]-3*mu[2]*mu[1]+2*mu[1]**3,mu[4]-4*mu[3]*mu[1]-3*mu[2]**2+12*mu[2]*mu[1]**2-6*mu[1]**4]
expected_log={};power={0:s.Integer(1)}
for k in range(1,5):
 power=multiply(power,j);expected_log=add(expected_log,scale(power,(-1)**k*kappa[k]/s.factorial(k)))
need('all_moments_to8fermions_equal_full_cumulant_source',eq(logF,expected_log))
need('first_induced_nu_up_fourfermion_is_minus_covariance_times_given_Ys',s.expand(-logF[51]+kappa[2]*y[0]*y[2])==0)
need('sixfermion_retains_third_boson_cumulant',s.expand(-logF[63]-kappa[3]*y[0]*y[1]*y[2])==0)
need('eightfermion_retains_fourth_boson_cumulant',s.expand(-logF[255]+kappa[4]*s.prod(y))==0)

# Same original complex h2 direction: bar(nu_L)nu_R bar(u_R)u_L.
hu,hb=s.symbols('h2 h2bar');yn,yuc,G=s.symbols('Ynu conjugateYu G_h2_h2bar')
Bnu=scale(multiply(variable(0),variable(1)),yn)
Bu=scale(multiply(variable(2),variable(3)),yuc)
Fcomplex=exponential(scale(add(scale(Bnu,hb),scale(Bu,hu)),-1),4)
weight={m:s.expand(v).subs(hu*hb,G).subs({hu:0,hb:0}) for m,v in Fcomplex.items()}
need('complex_H2_cross_fourfermion_has_no_extra_factor_two',s.expand(weight[15]-G*yn*yuc)==0)

# Finite complete Pfaffian polynomial remains valid at singular K.
a,b,c,d=s.symbols('a b c d');M=s.Matrix([[a,b],[c,d]])
askew=s.zeros(4);askew[:2,2:]=-M.T;askew[2:,:2]=M
pf=askew[0,1]*askew[2,3]-askew[0,2]*askew[1,3]+askew[0,3]*askew[1,2]
need('old_charged_skew_pair_Pfaffian_equals_oriented_single_determinant',s.expand(pf+M.det())==0)
need('Pfaffian_not_squared_determinant',s.expand(pf**2-M.det())!=0)
need('singular_case_keeps_source_coefficients_without_inverse',pf.subs({a:1,b:0,c:0,d:0})==0 and s.diff(pf,d).subs({a:1,b:0,c:0,d:0})==-1)

# Logical control on the ORIGINAL proof step: commuting positive marginal
# reflected forms do not suffice for a joint reflected form.
from itertools import product
r=s.Rational(1,2);total=s.Integer(0);xx=s.Integer(0);yy=s.Integer(0);mixed=s.Integer(0);minimum=1
for xm,xp,ym,yp in product([-1,1],repeat=4):
 w=(1-r*xm*xp*ym*yp)/16;minimum=min(minimum,w)
 total+=w;xx+=w*xm*xp;yy+=w*ym*yp;mixed+=w*xm*xp*ym*yp
need('commuting_marginal_control_has_positive_normalized_measure',total==1 and minimum>0)
need('separate_marginal_reflected_forms_are_product_positive',xx==0 and yy==0)
need('mixed_reflected_form_is_negative_despite_commutation',mixed==-r)

# All-order even half-action factorization needed by the repaired RP proof.
A=add(scale(multiply(variable(0),variable(1)),s.Symbol('a')),scale(multiply(variable(2),variable(3)),s.Symbol('b')))
B=add(scale(multiply(variable(4),variable(5)),s.Symbol('c')),scale(multiply(variable(6),variable(7)),s.Symbol('d')))
need('even_half_actions_commute',eq(multiply(A,B),multiply(B,A)))
need('whole_even_interaction_exponential_sews_exactly',eq(exponential(scale(add(A,B),-1),8),multiply(exponential(scale(A,-1),8),exponential(scale(B,-1),8))))

result={'contract':'ORIGINAL-DYNAMICAL-HIGGS-SOURCE-20261031','all_passed':all(checks.values()),'checks':checks,'original_basis_images':count,'Yukawa_coordinates':'72 real coefficient basis directions of supplied four 3x3 matrices, no search or selection','first_common_cross_action_coefficient':'-G_h2_h2bar * Ynu[a,b] * conjugate(Yu[d,c])','counterexample_mixed_reflected_form':str(mixed),'scope':'Actual full96 original real Higgs consumer; complete induced finite source at given original boson moments. No physical boson measure chosen, joint original reflection not inferred from JF, no continuum/TOE promotion.'}
Path(__file__).with_name('DYNAMICAL-HIGGS-SOURCE.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
