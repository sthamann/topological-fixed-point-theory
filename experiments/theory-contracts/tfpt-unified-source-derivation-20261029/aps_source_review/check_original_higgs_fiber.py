"""Exact original full96 Higgs orbit and positive-z determinant reduction.
No physical Y, source measure, profile, or new model is selected.
"""
from pathlib import Path
import sys,json
import numpy as np
import sympy as s
ROOT=Path(__file__).resolve().parents[4]
sys.path.insert(0,str(ROOT/'verification'))
import v252_full_finite_triple as original
checks={}
def need(k,b):
 checks[k]=bool(b)
 if not b: raise AssertionError(k)
I=original.IDX; z3=np.zeros((3,3),complex); z2=np.zeros((2,2),complex)
Q=[np.array([[0,1],[-1,0]],complex),np.array([[0,1j],[1j,0]],complex),np.eye(2,dtype=complex),np.diag([-1j,1j])]
P=original.rep(1,z2,z3);K=original.Kswap()
def higgs(q,Ys,MR):
 D=original.full_D(*Ys,z3);a=original.rep(0,q,z3)
 b=a@(D@P-P@D);b=b+b.conj().T
 return b+K@b.conj()@K+original.full_D(z3,z3,z3,z3,MR)
def scaled_lift(q,r):
 V=np.eye(48,dtype=complex)*r
 for up,dn,nc in [('nu','e',1),('u','d',3)]:
  for color in range(nc):
   for gen in range(3):
    ix=[I[(up,'L',color,gen)],I[(dn,'L',color,gen)]]
    V[np.ix_(ix,ix)]=q
 W=np.zeros((96,96),complex);W[:48,:48]=V;W[48:,48:]=V.conj()
 return W
x=s.symbols('x1:5',real=True);h1=x[0]+s.I*x[1];h2=x[2]+s.I*x[3]
q=s.Matrix([[s.conjugate(h2),h1],[-s.conjugate(h1),h2]]);u=sum(a*a for a in x)
need('generic_quaternion_norm_and_determinant',s.simplify(q.conjugate().T*q-u*s.eye(2))==s.zeros(2) and s.expand(q.det()-u)==0)
# One exact integral quaternion, r=5; full class proof is block multiplication.
qn=np.array([[2-4j,1+2j],[-1+2j,2+4j]],complex);r=5;W=scaled_lift(qn,r)
need('actual_full96_scaled_lift_is_unitary',np.array_equal(W.conj().T@W,25*np.eye(96)))
count=0
for sp in range(4):
 for a in range(3):
  for b in range(3):
   for phase in [1,1j]:
    Ys=[z3.copy() for _ in range(4)];Ys[sp][a,b]=phase
    radial=original.full_D(*[r*Y for Y in Ys],z3)
    need_result=np.array_equal(W@radial@W.conj().T,25*higgs(qn,Ys,z3))
    if not need_result:raise AssertionError('Y orbit basis')
    count+=1
for a in range(3):
 for b in range(a,3):
  for phase in [1,1j]:
   M=z3.copy();M[a,b]=M[b,a]=phase
   D=original.full_D(z3,z3,z3,z3,M)
   if not np.array_equal(W@D@W.conj().T,25*D):raise AssertionError('Majorana orbit basis')
   count+=1
need('all84_original_real_coefficient_directions_transform_with_same_left_lift',count==84)
need('weak_lift_preserves_internal_grading_and_real_structure',np.array_equal(W@original.grading(),original.grading()@W) and np.array_equal(K@W.conj()@K,W))

# Dense Gaussian-integer coefficient control, independent original constructor.
Ys=[np.array([[1+1j,2,0],[0,1,1j],[1,-1j,2]],complex),np.array([[2,1j,1],[1,3,-1j],[0,1,1]],complex),np.array([[1,0,1j],[2,1,0],[0,-1j,2]],complex),np.array([[1,1j,0],[1,2,1],[0,1,3j]],complex)]
M=np.array([[2,1j,1],[1j,3,2],[1,2,1j]],complex)
rr=2;D=original.full_D(*[rr*Y for Y in Ys],M)
ixplus=[I[('nu','L',0,g)] for g in range(3)]+[48+I[('nu','R',0,g)] for g in range(3)]
ixminus=[I[('nu','R',0,g)] for g in range(3)]+[48+I[('nu','L',0,g)] for g in range(3)]
B=np.block([[rr*Ys[3],z3],[M.conj().T,rr*Ys[3].T]])
need('original_neutrino_chiral_block_is_exact6x6',np.array_equal(D[np.ix_(ixplus,ixminus)],B))
D2=D@D;det_checks=[]
def mat(a):return s.Matrix(a).applyfunc(s.nsimplify)
# Whole determinant is reconstructed from original disjoint invariant blocks;
# exact SymPy determinants on blocks, not floating det on a96 matrix.
nu=ixplus+ixminus
allslots=set(nu);direct_nu=mat(D2[np.ix_(nu,nu)])
for zz in [1,3]:
 left=direct_nu+zz*s.eye(12)
 right=mat(B.conj().T@B)+zz*s.eye(6)
 need('exact_neutrino_squared_determinant_z'+str(zz),left.det()==right.det()**2)
 det_checks.append(str(right.det()))
 Ymat=mat(Ys[3]);Mmat=mat(M);AA=zz*s.eye(3)+rr**2*s.conjugate(Ymat)*Ymat.T
 EE=zz*s.eye(3)+rr**2*Ymat.conjugate().T*Ymat+Mmat*Mmat.conjugate().T-rr**2*Mmat*Ymat.T*AA.inv()*s.conjugate(Ymat)*Mmat.conjugate().T
 need('exact_neutrino3x3_Schur_determinant_z'+str(zz),s.simplify(AA.det()*EE.det()-right.det())==0)
for sp,nc,j in [('u',3,0),('d',3,1),('e',1,2)]:
 for color in range(nc):
  ix=[I[(sp,ch,color,g)]+offset for offset in [0,48] for ch in ['L','R'] for g in range(3)]
  allslots.update(ix);block=mat(D2[np.ix_(ix,ix)])
  qsmall=mat(rr**2*Ys[j].conj().T@Ys[j])
  for zz in [1,3]:
   if (block+zz*s.eye(12)).det()!=(qsmall+zz*s.eye(3)).det()**4:raise AssertionError('charged det')
need('whole96_is_partitioned_with_color12_and_lepton4_multiplicities',allslots==set(range(96)))
need('all14_exact_charged_block_determinant_controls_passed',True)
# Coefficient of u=r² in log det(Phi(r)²+z), entire symmetric-MR class.
Ynu=mat(Ys[3]);Mm=mat(M);Yzero=original.full_D(*Ys,z3);Dzero=original.full_D(z3,z3,z3,z3,M)
Qzero=mat(Dzero@Yzero+Yzero@Dzero);YY=mat(Yzero@Yzero);DD=mat(Dzero@Dzero)
jet_controls=[]
for zz in [1,3]:
 R=(DD+zz*s.eye(96)).inv()
 exact=s.trace(R*YY)-s.trace(R*Qzero*R*Qzero)/2
 charged=4*sum(mult*s.trace(mat(Ys[j].conj().T@Ys[j]))/zz for j,mult in [(0,3),(1,3),(2,1)])
 neutral=4*s.trace((zz*s.eye(3)+Mm*Mm.conjugate().T).inv()*Ynu.conjugate().T*Ynu)
 need('whole96_first_mass_jet_reduces_to_positive3x3_resolvents_z'+str(zz),s.simplify(exact-charged-neutral)==0)
 jet_controls.append(str(s.simplify(exact)))
z=s.symbols('z',positive=True);a,b,c,d,m=s.symbols('a b c d m',real=True)
# 1-family algebra control of the all-family Schur identity.
Y=s.Matrix([[a+s.I*b]]);MM=s.Matrix([[m+s.I*d]]);tt=s.symbols('t',real=True);uu=tt**2
Bn=s.BlockMatrix([[tt*Y,s.zeros(1)],[MM.conjugate().T,tt*Y.T]]).as_explicit()
Kn=z*s.eye(2)+Bn.conjugate().T*Bn
need('generic_neutrino_polynomial_has_no_odd_radial_power',s.expand(Kn.det().subs(tt,-tt)-Kn.det())==0)
need('one_family_neutrino_determinant_closed_polynomial',s.expand(Kn.det()-((z+uu*(a*a+b*b))**2+z*(m*m+d*d)))==0)

# Fixed literal original bulk coefficients; this is a SOURCE MATCH test.
chi,w,Z=s.symbols('chi w Z',positive=True);U=s.symbols('u',nonnegative=True)
ED=-U/(48*s.pi**2);G4=U**2/(16*s.pi**2)
printed=s.expand(-6*chi**2*ED+G4)
target=-chi**2*w**2*U/(8*s.pi**2)+U**2/(16*s.pi**2)
need('printed_relative_bulk_has_positive_quadratic_and_quartic',s.diff(printed,U).subs(U,0)==chi**2/(8*s.pi**2) and s.diff(printed,U,2)==1/(8*s.pi**2))
need('literal_bulk_differs_from_stated_broken_higgs_potential',s.factor(printed-target)==chi**2*U*(1+w*w)/(8*s.pi**2))
need('literal_bulk_also_differs_from_minusV_Lagrangian',s.expand(-target-printed)==s.expand(chi**2*(w*w-1)*U/(8*s.pi**2)-U**2/(8*s.pi**2)))
need('relative_coefficient_sign_is_invariant_under_common_action_sign',s.sign(s.diff(printed,U).subs(U,0)*s.diff(printed,U,2))==1 and s.sign(s.diff(target,U).subs(U,0)*s.diff(target,U,2))==-1)
need('positive_field_rescaling_does_not_remove_sign_mismatch',s.diff(printed.subs(U,U/Z),U).subs(U,0)==chi**2/(8*s.pi**2*Z))
need('zero_higgs_bulk_tadpole_but_positive_radial_curvature',s.diff(printed.subs(U,tt*tt),tt).subs(tt,0)==0 and s.diff(printed.subs(U,tt*tt),tt,2).subs(tt,0)==chi**2/(4*s.pi**2))

# Real-coordinate derivative convention, no spurious factor1/2.
F=s.Function('F');f=F(u);grad=s.Matrix([s.diff(f,a) for a in x]);H=s.hessian(f,x)
# Substitute a separate scalar symbol to differentiate F without d/d(sum x²).
scalar=s.symbols('v');Fp=s.Subs(s.diff(F(scalar),scalar),scalar,u);Fpp=s.Subs(s.diff(F(scalar),scalar,2),scalar,u)
need('whole_four_real_Hessian_has_correct_coordinate_factors',s.simplify(H-2*Fp*s.eye(4)-4*Fpp*s.Matrix(x)*s.Matrix(x).T)==s.zeros(4))
# Independent all-mode-sum control; analytic convergence and identity in proof.
import math
NN=10000
for alpha in [0.25,0.75]:
 for vv in [1.,4.]:
  partial=math.fsum(1/((n+alpha)**2+vv) for n in range(-NN,NN+1))
  closed=math.pi/math.sqrt(vv)*math.sinh(2*math.pi*math.sqrt(vv))/(math.cosh(2*math.pi*math.sqrt(vv))-math.cos(2*math.pi*alpha))
  need('quarter_circle_closed_sum_with_positive_tail_'+str(alpha)+'_'+str(vv),0<=closed-partial<=2/(NN-1))
result={'contract':'ORIGINAL-HIGGS-FIBER-AND-RESIDUE-MATCH-20261031','all_passed':all(checks.values()),'checks':checks,'actual_original_coefficient_basis_images':count,'exact_neutrino_det_z1_z3':det_checks,'exact_full96_logdet_mass_jet_z1_z3':jet_controls,'literal_bulk_potential':str(printed),'target_potential_minus_literal':str(s.factor(target-printed)),'target_minusV_Lagrangian_minus_literal':str(s.expand(-target-printed)),'scope':'Whole constant original96 Higgs family, covariant internal lift with given Y and symmetric MR; conditional common-domain spectral consequence. Literal archived local residue match only; no replacement of full S_adm, no full Higgs covariance or physical instability verdict.'}
Path(__file__).with_name('ORIGINAL-HIGGS-FIBER.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
