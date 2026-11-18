"""R030/R031: local residue compatibility and conditional original96 consumer.
No new physical Y, profile, measure, spacetime, or ledger promotion.
"""
from pathlib import Path
import sys, json
import numpy as np
import sympy as s
ROOT=Path(__file__).resolve().parents[4]
sys.path.insert(0,str(ROOT/'verification'))
import v252_full_finite_triple as original
checks={}
def need(name, value):
    checks[name]=bool(value)
    if not value: raise AssertionError(name)
R,u,k,l,c,d,chi,A,B,C=s.symbols('R u k l c d chi A B C',real=True)
X=s.symbols('X',real=True)
# Local scalar-symbol heat identities, valid for each common linear trace.
E=-R/s.Integer(4)-X
a2=E+R/6
a4=E**2/2+R*E/6
need('same_linear_trace_mass_curvature_ratio_minus12',s.expand(a2).coeff(X)==-12*s.expand(a4).coeff(R).coeff(X))
printed2=-s.Rational(1,48); printedRu=s.Rational(1,96)
need('literal_printed_ratio_minus2_not_minus12',printed2/printedRu==-2 and printed2!=-12*printedRu)
need('printed_curvature_requires_sixfold_mass_coefficient',(-12*printedRu)/printed2==6)
need('mass_kept_requires_sixfold_smaller_curvature',(-printed2/12)/printedRu==s.Rational(1,6))
need('projector_quartic_to_mass_ratio_onehalf',s.simplify((k/2)/k)==s.Rational(1,2))
need('printed_quartic_to_mass_ratio_three',s.Rational(1,16)/(-printed2)==3)
xx=s.symbols('xx',real=True); hh=s.Function('h')(xx)
need('displayed_E_derivative_is_total_derivative_not_independent_kinetic',s.expand(s.diff(hh**2,xx,2)-2*s.diff(hh,xx)**2-2*hh*s.diff(hh,xx,2))==0)
need('common_weight_and_field_map_preserve_linear_ratio',s.simplify((-k*l)/((k*l)/12))==-12)
need('literal_ordinary48_rank2_printed_ED_is_raw_a2_over6',-s.Rational(48,12)/6==-s.Rational(16,24) and -s.Rational(2,1)/6==-s.Rational(16,48))
# Explicit Clifford trace independently checks kinetic factor and signs.
pa=[s.Matrix([[0,1],[1,0]]),s.Matrix([[0,-s.I],[s.I,0]]),s.diag(1,-1)]
gam=[s.kronecker_product(pa[0],p) for p in pa]+[s.kronecker_product(pa[1],s.eye(2))]
g5=gam[0]*gam[1]*gam[2]*gam[3]
need('euclidean_clifford_and_chirality',all(a*b+b*a==(2*s.eye(4) if i==j else s.zeros(4)) for i,a in enumerate(gam) for j,b in enumerate(gam)) and g5*g5==s.eye(4))
need('spin_trace_of_derivative_square_positive4',all(s.trace((s.I*a*g5)*(s.I*b*g5))==(4 if i==j else 0) for i,a in enumerate(gam) for j,b in enumerate(gam)))
need('spin_trace_of_derivative_term_zero',all(s.trace(s.I*a*g5)==0 for a in gam))
# Actual unmodified original constructor. Integer arithmetic stays exactly
# representable in complex128 for these controls; no rounded comparisons.
z3=np.zeros((3,3),complex);z2=np.zeros((2,2),complex)
Ys=[np.array([[1+1j,2,0],[0,1,1j],[1,-1j,2]],complex),np.array([[2,1j,1],[1,3,-1j],[0,1,1]],complex),np.array([[1,0,1j],[2,1,0],[0,-1j,2]],complex),np.array([[1,1j,0],[1,2,1],[0,1,3j]],complex)]
M=np.array([[2,1j,1],[1j,3,2],[1,2,1j]],complex)
DY=original.full_D(*Ys,z3); DM=original.full_D(z3,z3,z3,z3,M)
P=original.rep(1,z2,z3);K=original.Kswap()
def higgs(q):
    a=original.rep(0,q,z3);b=a@(DY@P-P@DY);b=b+b.conj().T
    return b+K@b.conj()@K+DM
def tr(a):return np.trace(a)
AA=tr(DY@DY);BB=tr(np.linalg.matrix_power(DY,4));CC=tr(DY@DY@DM@DM)
Y2=sum(m*tr(Y.conj().T@Y) for m,Y in zip([3,3,1,1],Ys))
Y4=sum(m*tr(np.linalg.matrix_power(Y.conj().T@Y,2)) for m,Y in zip([3,3,1,1],Ys))
T=tr(Ys[3].conj().T@Ys[3]@M@M.conj().T)
need('original96_A_equals4Y2',AA==4*Y2)
need('original96_B_equals4Y4',BB==4*Y4)
need('original96_C_equals2Majorana_trace',CC==2*T)
need('Majorana_C_is_positive_Frobenius_square',CC==tr((DY@DM).conj().T@(DY@DM)) and CC.real>=0 and CC.imag==0)
need('forbidden_original_M_Y_M_path_zero',np.array_equal(DM@DY@DM,np.zeros((96,96))))
for i,coords in enumerate([[1,2,2,4],[0,0,2,0],[2,0,0,1]]):
    x1,x2,x3,x4=coords;h1=x1+1j*x2;h2=x3+1j*x4
    q=np.array([[h2.conjugate(),h1],[-h1.conjugate(),h2]])
    F=higgs(q);uu=sum(x*x for x in coords)
    need('whole96_second_moment_'+str(i),tr(F@F)-tr(DM@DM)==AA*uu)
    need('whole96_fourth_moment_including_Majorana_'+str(i),tr(np.linalg.matrix_power(F,4))-tr(np.linalg.matrix_power(DM,4))==BB*uu**2+4*CC*uu)
Q=[np.array([[0,1],[-1,0]],complex),np.array([[0,1j],[1j,0]],complex),np.eye(2),np.diag([-1j,1j])]
L=[higgs(q)-DM for q in Q]
need('whole_four_real_derivative_Gram_A_identity',all(tr(a@b)==(AA if i==j else 0) for i,a in enumerate(L) for j,b in enumerate(L)))
# Conditional standard product Dirac, spin4 separately from internal96.
Z=d*A/(8*s.pi**2)
quad=-c*chi**2*A/(4*s.pi**2)+d*C/(2*s.pi**2)+d*A*R/(48*s.pi**2)
quartic=d*B/(8*s.pi**2)
need('canonical_mass_including_curvature_and_Majorana',s.simplify(-quad/Z-(2*c*chi**2/d-4*C/A-R/6))==0)
need('canonical_quartic_uses_fourth_moment_not_square_of_second',s.simplify(quartic/Z**2-8*s.pi**2*B/(d*A**2))==0)
need('printed_minus6_has_positive_local_quadratic_at_R0',s.expand(quad.subs({c:-6,d:1,R:0}))==3*A*chi**2/(2*s.pi**2)+C/(2*s.pi**2))
# Actual two-sheet dictionary must be supplied, not a freely frozen reference.
sigma,rr,ff,ss=s.symbols('sigma rr ff ss')
need('finite_part_scale_pole_mixing',s.expand(s.series(s.exp(2*sigma*ss)*(rr/ss+ff),ss,0,1).removeO()).coeff(ss,0)==ff+2*sigma*rr)
result={'contract':'ORIGINAL-HIGGS-LOCAL-HEAT-CONSUMER-20261031','all_passed':all(checks.values()),'checks':checks,'control_invariants':{'A':str(AA),'B':str(BB),'C':str(CC),'T':str(T)},'required_mass_to_Ru_ratio':-12,'printed_mass_to_Ru_ratio':-2,'mass_factor_if_Ru_kept':6,'scope':'Exact original96 moment/derivative identities; local standard same-trace residue compatibility only. Four-dimensional product Dirac is an explicit conditional consumer, not a derived whole original source, measure, continuation or physical vacuum. No new Yukawa fit or ledger promotion.'}
Path(__file__).with_name('ORIGINAL-HIGGS-HEAT.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
