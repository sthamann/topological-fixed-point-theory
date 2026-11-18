"""Scoped algebra controls at the existing full96 consumer and block reader.
No new physical projection, heavy sector, geometry or original trace is chosen.
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
def eq(a,b): return all(s.cancel(v)==0 for v in a-b)
I=s.I;eye=s.eye(2);zero=s.zeros(2)
pauli=[s.Matrix([[0,1],[1,0]]),s.Matrix([[0,-I],[I,0]]),s.diag(1,-1)]
gam=[s.kronecker_product(pauli[0],eye)]+[s.kronecker_product(pauli[1],x) for x in pauli]
g5=gam[0]*gam[1]*gam[2]*gam[3]
need('Euclidean_Clifford_and_chirality',all(eq(a*b+b*a,2*s.eye(4) if i==j else s.zeros(4)) for i,a in enumerate(gam) for j,b in enumerate(gam)) and eq(g5*g5,s.eye(4)))
odd=[I*g*g5 for g in gam];two=[(gam[i]*gam[j]-gam[j]*gam[i])/2 for i in range(4) for j in range(i+1,4)]
need('all_four_derivative_spin_components_orthogonal_to_all_six_curvatures',all(s.trace(a.H*b)==0 for a in odd for b in two))
need('all_four_derivative_components_have_exact_unit_projection',all(s.trace(a.H*b)==(4 if i==j else 0) for i,a in enumerate(odd) for j,b in enumerate(odd)))
need('derivative_spin_trace_and_mixed_scalar_trace_vanish',all(s.trace(a)==0 for a in odd))
need('derivative_spin_matrices_are_Hermitian',all(eq(a.H,a) for a in odd))
d=s.symbols('d0:4',real=True);qd=sum((x*a for x,a in zip(d,odd)),s.zeros(4))
need('whole_derivative_square_spin_trace',s.expand(s.trace(qd*qd)-4*sum(x*x for x in d))==0)
R,u=s.symbols('R u',real=True);E=-R/4-u
a2=s.expand(E+R/6);a4=s.expand((60*R*E+180*E*E)/360)
need('constant_field_heat_ratio_remains_minus12',a2.coeff(u)/a4.coeff(R,1).coeff(u,1)==-12)
# Real original96 constructor and oneform, diagnostic integer matrices only.
z3=np.zeros((3,3),complex);z2=np.zeros((2,2),complex)
ys=[np.array([[1+1j,2,0],[0,1,1j],[1,-1j,2]],complex),np.array([[2,1j,1],[1,3,-1j],[0,1,1]],complex),np.array([[1,0,1j],[2,1,0],[0,-1j,2]],complex),np.array([[1,1j,0],[1,2,1],[0,1,3j]],complex)]
mr=np.array([[2,1j,1],[1j,3,2],[1,2,1j]],complex)
DY=original.full_D(*ys,z3);DM=original.full_D(z3,z3,z3,z3,mr)
eR=original.rep(1,z2,z3);swap=original.Kswap()
def L(h):
    h1,h2=np.asarray(h,complex);q=np.array([[h2.conjugate(),h1],[-h1.conjugate(),h2]])
    a=original.rep(0,q,z3);w=a@(DY@eR-eR@DY);w=w+w.conj().T
    return w+swap@w.conj()@swap
h=np.array([1+2j,3-1j]);Y=L(h);Phi=Y+DM
P=(np.eye(96)+original.grading())/2;Q=np.eye(96)-P
A=4*(3*np.trace(ys[0].conj().T@ys[0]+ys[1].conj().T@ys[1])+np.trace(ys[2].conj().T@ys[2]+ys[3].conj().T@ys[3])).real
need('actual_original96_chirality_projection_annihilates_odd_vertex',np.array_equal(P@Phi@P,np.zeros((96,96))))
need('actual_original96_square_projection_defect_is_not_zero',np.array_equal(P@Phi@Phi@P,P@Phi@Q@Phi@P) and np.trace(P@Phi@Phi@P).real>0)
need('actual_original96_Higgs_defect_trace_is_Au_half',np.trace(P@Y@Y@P)==A*np.vdot(h,h)/2)
need('actual_original96_projection_defect_has_positive_Gram_form',np.array_equal(P@Phi@Q@Phi@P,(Q@Phi@P).conj().T@(Q@Phi@P)))
need('actual_original96_mass_square_commutes_with_finite_chirality',np.array_equal(P@Phi@Phi,Phi@Phi@P))
for k,v in enumerate([[1,0],[1j,0],[0,1],[0,1j]]):
    dv=L(v)
    need('actual_original96_derivative_moment_direction_'+str(k),np.trace(dv@dv)==A)
    need('actual_original96_weighted_heat_retains_half_kinetic_direction_'+str(k),np.array_equal(P@dv@P,np.zeros((96,96))) and 4*np.trace(P@dv@Q@dv@P)==2*A)
# Three-slot matrices below test identities only, not a TFPT normal model.
t,z=s.symbols('t z',real=True)
K0=s.Matrix([[2,1,1],[1,5,0],[1,0,7]])
V=s.Matrix([[1,1,-1],[1,2,1],[-1,1,-1]])
K=K0+t*V;M=K+2*s.eye(3)
AA=M[:1,:1];BB=M[:1,1:];CC=M[1:,:1];EE=M[1:,1:];RR=EE.inv();SS=AA-BB*RR*CC
invblock=SS.inv().row_join(-SS.inv()*BB*RR).col_join((-RR*CC*SS.inv()).row_join(RR+RR*CC*SS.inv()*BB*RR))
need('complete_block_inverse_including_all_exterior_entries',eq(M*invblock,s.eye(3)))
need('finite_determinant_includes_exterior_factor',s.cancel(M.det()-EE.det()*SS.det())==0)
eta=s.Matrix(s.symbols('e0:3'));bar=s.Matrix(1,3,s.symbols('b0:3'))
source=(bar[:,:1]-bar[:,1:]*RR*CC)*SS.inv()*(eta[:1,:]-BB*RR*eta[1:,:])+bar[:,1:]*RR*eta[1:,:]
need('complete_sources_have_two_shifted_legs_and_direct_exterior_term',s.cancel((bar*invblock*eta)[0]-source[0])==0)
# Spectral derivative uses one fixed t, while preserving all couplings.
Mz=K.subs(t,s.Rational(1,3))-z*s.eye(3)
Az,Bz,Cz,Ez=Mz[:1,:1],Mz[:1,1:],Mz[1:,:1],Mz[1:,1:]
Rz=Ez.inv();Sz=Az-Bz*Rz*Cz
need('Schur_spectral_derivative_has_exterior_weight',eq(Sz.diff(z),-s.eye(1)-Bz*Rz*Rz*Cz))
need('full_resolvent_trace_needs_log_derivative_and_exterior',s.cancel(s.trace(Mz.inv())-s.trace(Rz)+s.trace(Sz.inv()*Sz.diff(z)))==0)
need('naive_projected_resolvent_trace_is_not_whole_trace',s.cancel((s.trace(Mz.inv())-s.trace(Sz.inv())).subs(z,-2))!=0)
Mn=K.subs(t,s.Rational(1,3));om=s.Rational(3,2)
need('positive_squared_operator_resolvent_uses_both_signed_resolvents',eq((Mn*Mn+om**2*s.eye(3)).inv(),((Mn-I*om*s.eye(3)).inv()-(Mn+I*om*s.eye(3)).inv())/(2*I*om)))
# Explicit complete Higgs jets with nonzero background off-diagonal blocks.
A0,A1=AA.subs(t,0),AA.diff(t)
B0,B1=BB.subs(t,0),BB.diff(t);C0,C1=CC.subs(t,0),CC.diff(t)
E0,E1=EE.subs(t,0),EE.diff(t);R0=E0.inv()
rs=[(-R0*E1)**n*R0 for n in range(5)]
ss=[]
for n in range(5):
    val=A0 if n==0 else A1 if n==1 else s.zeros(1)
    val=val-sum((b*rs[n-i-j]*c for i,b in enumerate([B0,B1]) for j,c in enumerate([C0,C1]) if i+j<=n),s.zeros(1))
    ss.append(val)
    actual=SS.diff(t,n).subs(t,0)/s.factorial(n)
    need('complete_Schur_jet_'+str(n)+'_retains_background_and_exterior',eq(val,actual))
xx=[ss[0].inv()*sn for sn in ss]
l2=s.trace(xx[2]-xx[1]**2/2)-s.trace((R0*E1)**2)/2
l4=s.trace(xx[4]-xx[1]*xx[3]-xx[2]**2/2+xx[1]**2*xx[2]-xx[1]**4/4)-s.trace((R0*E1)**4)/4
G0=M.subs(t,0).inv()
need('full_logdet_H2_jet_equals_original_unsplit_matrix',s.cancel(l2+s.trace((G0*V)**2)/2)==0)
need('full_logdet_H4_jet_equals_original_unsplit_matrix',s.cancel(l4+s.trace((G0*V)**4)/4)==0)
# Noncommuting P block checks: matrix P dimension2 makes cyclic order matter.
M4=s.Matrix([[4,1,1,0],[1,6,0,1],[1,0,8,1],[0,1,1,10]])
V4=s.Matrix([[1,2,1,1],[2,-1,0,1],[1,0,2,1],[1,1,1,-2]])
T4=M4+t*V4;a,b,c,e=T4[:2,:2],T4[:2,2:],T4[2:,:2],T4[2:,2:];rr=e.subs(t,0).inv()
rn=[(-rr*e.diff(t))**n*rr for n in range(5)]
sn=[]
for n in range(5):
    val=a.subs(t,0) if n==0 else a.diff(t) if n==1 else s.zeros(2)
    val-=sum((bb*rn[n-i-j]*cc for i,bb in enumerate([b.subs(t,0),b.diff(t)]) for j,cc in enumerate([c.subs(t,0),c.diff(t)]) if i+j<=n),s.zeros(2))
    sn.append(val)
xn=[sn[0].inv()*v for v in sn]
need('noncommuting_Higgs_jet_blocks_are_actually_noncommuting',xn[1]*xn[2]!=xn[2]*xn[1])
L4=s.trace(xn[4]-xn[1]*xn[3]-xn[2]**2/2+xn[1]**2*xn[2]-xn[1]**4/4)-s.trace((rr*e.diff(t))**4)/4
need('noncommuting_fourth_jet_cyclic_ordering_is_correct',L4==-s.trace((M4.inv()*V4)**4)/4)
result={'contract':'ORIGINAL-HIGGS-SUPERCONNECTION-REDUCTION-20261031','all_passed':all(checks.values()),'checks':checks,'scope':'Existing conditional spin4/full96 original consumer; exact finite block identities and jets. Three/four-slot diagnostics are identity tests, not chosen physical normal operators. Actual normal projection, domain, weighted trace and infinite relative continuation remain original inputs; no new Lean proof or physical closure.'}
Path(__file__).with_name('ORIGINAL-HIGGS-SUPERCONNECTION-REDUCTION.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'all_passed':result['all_passed'],'count':len(checks),'failed':[k for k,v in checks.items() if not v]}))
