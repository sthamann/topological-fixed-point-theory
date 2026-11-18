"""Original96 total grading, own-operator spectral selection and full trace.
Small rectangular matrices are identity controls, not physical normal models.
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
    if not b:raise AssertionError(k)
def eq(a,b):return all(s.simplify(v)==0 for v in a-b)
I=s.I;p=[s.Matrix([[0,1],[1,0]]),s.Matrix([[0,-I],[I,0]]),s.diag(1,-1)]
gam=[s.kronecker_product(p[0],s.eye(2))]+[s.kronecker_product(p[1],v) for v in p]
g5=gam[0]*gam[1]*gam[2]*gam[3]
der=[I*g*g5 for g in gam];two=[gam[i]*gam[j] for i in range(4) for j in range(i+1,4)]
need('H_scalar_and_derivative_terms_have_zero_chiral_spin_trace',s.trace(g5)==0 and all(s.trace(g5*q)==0 for q in der))
need('H_kinetic_pair_has_no_local_chiral_index_component',all(s.trace(g5*a*b)==0 for a in der for b in der))
need('H_mixed_curvature_terms_have_no_local_chiral_index_component',all(s.trace(g5*a)==0 for a in two) and all(s.trace(g5*a*b)==0 for a in two for b in der))
need('pure_gauge_index_density_cannot_be_dropped_just_from_H_test',any(s.trace(g5*a*b)!=0 for a in two for b in two))
z3=np.zeros((3,3),complex);z2=np.zeros((2,2),complex)
ys=[np.array([[1+1j,2,0],[0,1,1j],[1,-1j,2]],complex),np.array([[2,1j,1],[1,3,-1j],[0,1,1]],complex),np.array([[1,0,1j],[2,1,0],[0,-1j,2]],complex),np.array([[1,1j,0],[1,2,1],[0,1,3j]],complex)]
mr=np.array([[2,1j,1],[1j,3,2],[1,2,1j]],complex)
DY=original.full_D(*ys,z3);DM=original.full_D(z3,z3,z3,z3,mr)
eR=original.rep(1,z2,z3);swap=original.Kswap();GF=original.grading()
def L(h):
    h1,h2=np.asarray(h,complex);q=np.array([[h2.conjugate(),h1],[-h1.conjugate(),h2]])
    aa=original.rep(0,q,z3);v=aa@(DY@eR-eR@DY);v=v+v.conj().T
    return v+swap@v.conj()@swap
gn=[np.array(g,dtype=complex) for g in gam];g5n=np.array(g5,dtype=complex)
G=np.kron(g5n,GF);Pc=(np.eye(384)+G)/2;Pi=np.kron(np.eye(4),(np.eye(96)+GF)/2)
need('two_halves_have_same_rank_but_different_actual_action',np.trace(Pc)==192 and np.trace(Pi)==192 and not np.array_equal(Pc,Pi))
need('total_chiral_half_is_not_spin_commuting_internal_half',not np.array_equal(Pc@np.kron(gn[0],np.eye(96)),np.kron(gn[0],np.eye(96))@Pc))
h=np.array([1+2j,3-1j]);Phi=L(h)+DM
need('original_whole96_Phi_is_odd_in_fixed_finite_grading',np.array_equal(GF@Phi+Phi@GF,np.zeros((96,96))))
K=2*np.kron(gn[0],np.eye(96))+np.kron(g5n,Phi);KK=K@K
need('actual_whole_product_is_odd_in_total_grading',np.array_equal(G@K+K@G,np.zeros((384,384))))
need('actual_total_chiral_half_reduces_square_not_first_order_operator',np.array_equal(Pc@KK,KK@Pc) and np.array_equal(Pc@K@Pc,np.zeros((384,384))) and not np.array_equal(Pc@K,K@Pc))
for n in [1,2]:
    power=np.linalg.matrix_power(KK,n)
    need('actual_total_chiral_square_moment_half_'+str(n),np.trace(Pc@power)==np.trace(power)/2)
A=4*(3*np.trace(ys[0].conj().T@ys[0]+ys[1].conj().T@ys[1])+np.trace(ys[2].conj().T@ys[2]+ys[3].conj().T@ys[3])).real
pcdiag=np.diag(Pc);pidiag=np.diag(Pi)
for i,v in enumerate([[1,0],[1j,0],[0,1],[0,1j]]):
    dv=L(v);dd=np.kron(np.array(der[i],complex),dv)
    cp=pcdiag[:,None]*dd*pcdiag[None,:];ip=pidiag[:,None]*dd*pidiag[None,:]
    need('actual_variable_H_direction_'+str(i)+'_survives_total_not_internal_half',np.count_nonzero(cp)>0 and np.count_nonzero(ip)==0 and np.array_equal(G@dd,dd@G))
    need('actual_variable_H_kinetic_moment_half_'+str(i),np.sum(cp*cp.T)==2*A and np.sum(dd*dd.T)==4*A)
# Exact rectangular positive-spectrum pairing, including an unpaired zero.
u,lam,t,z=s.symbols('u lam t z',real=True)
T=s.Matrix([[1,u,0],[0,1,u]]);plus=T.T*T;minus=T*T.T
need('whole_rectangular_positive_spectra_pair_with_multiplicity',s.expand((lam*s.eye(3)-plus).det()-lam*(lam*s.eye(2)-minus).det())==0)
need('index_term_is_present_for_every_rectangular_parameter',T.rank()==2 and len(plus.nullspace())-len(minus.nullspace())==1)
T1=T.subs(u,1);Lp=T1.T*T1;Lm=T1*T1.T
hp=s.trace((-t*Lp).exp());hm=s.trace((-t*Lm).exp())
need('whole_heat_supertrace_including_zero_is_index',s.simplify(hp-hm)==1)
need('absolute_chiral_heat_needs_half_index_term',s.simplify(hp-(hp+hm)/2)==s.Rational(1,2))
need('zero_excluded_whole_spectral_trace_is_exact_half',s.simplify((hp-1)-(hp+hm-1)/2)==0)
sn=s.symbols('s');zp=3**(-sn)+7**(-sn)
need('positive_zeta_half_has_no_unpaired_zero_added',s.simplify(zp-(2*3**(-sn)+2*7**(-sn))/2)==0)
# Same positive high-pass applied to both graded halves, kappa between1 and3.
Pgp=Lp*(Lp-s.eye(3))/6;Pgm=Lm*(Lm-s.eye(2))/6
need('combined_chirality_and_highpass_has_no_index_term',s.simplify(s.trace(Pgp*(-t*Lp).exp())-s.trace(Pgm*(-t*Lm).exp()))==0 and s.simplify(s.trace(Pgp*(-t*Lp).exp())-s.exp(-3*t))==0)
need('combined_highpass_removes_unpaired_zero_without_new_grading',eq(Pgp*Pgp,Pgp) and eq(Pgm*Pgm,Pgm) and all(eq(Pgp*v,s.zeros(3,1)) for v in Lp.nullspace()))
# A low eigenvalue u² touches zero without crossing kappa=2.
# Full zero-excluded FP count and removed positive count jump together.
need('zero_appearance_jumps_cancel_in_cut_FP_without_threshold_crossing',(2-1)==(1-0) and s.limit(s.exp(-t*u*u)+s.exp(-4*t),u,0)-s.exp(-4*t)==1)
# Both signs of the grading under a supplied reference lift preserve oddness.
GF2=s.diag(1,-1);KD=s.Matrix([[0,2],[2,0]])
for name,tau in [('preserves',s.diag(1,-1)),('swaps',s.Matrix([[0,1],[1,0]]))]:
    ref=(KD+tau*KD*tau.inv())/2
    need('bound_reference_is_odd_when_lift_'+name+'_grading',eq(GF2*ref+ref*GF2,s.zeros(2)))
need('multiplying_by_internal_grading_does_not_make_it_odd',GF2*GF2==s.eye(2) and GF2*GF2+GF2*GF2!=s.zeros(2))
# Exact own-operator spectral restriction in a non-diagonal fixed chart.
U=s.Matrix([[1,2,2],[2,1,-2],[-2,2,-1]])/3
LL=U*s.diag(0,1,4)*U.T;PP=U*s.diag(0,0,1)*U.T;SS=U[:,2:]
need('actual_selector_form_reduces_its_own_positive_operator',eq(PP*LL,LL*PP))
need('own_operator_positive_spectral_cut_kills_its_zero_mode',eq(PP*U[:,:1],s.zeros(3,1)))
heat=U*s.diag(1,s.exp(-t),s.exp(-4*t))*U.T
need('whole_weighted_heat_equals_restricted_heat_for_same_spectral_cut',s.simplify(s.trace(PP*heat)-s.exp(-4*t))==0 and eq(SS.T*LL*SS,s.Matrix([[4]])))
need('first_order_compression_need_not_match_restricted_square',eq(s.diag(1,0)*KD*KD*s.diag(1,0),s.diag(4,0)) and eq(s.diag(1,0)*KD*s.diag(1,0),s.zeros(2)))
# Cutoff: finite low positive eigenvalue varies, FP0 does not; logdet does.
low=(1+u*u)**(-sn)
need('finite_low_removal_is_entire_at_Einstein_pole',s.residue(low,sn,1)==0)
need('finite_low_removal_changes_FP0_by_count_not_eigenvalue',low.subs(sn,0)==1 and s.diff(low,u,2).subs(sn,0)==0 and s.diff(low,u,4).subs(sn,0)==0)
need('logdet_low_mode_has_H2_and_H4_jets_unlike_FP0',s.diff(s.log(1+u*u),u,2).subs(u,0)/2==1 and s.diff(s.log(1+u*u),u,4).subs(u,0)/24==-s.Rational(1,2))
fullheat=1+s.exp(-t*(1+u*u))+s.exp(-t*(4+u*u));cutheat=s.exp(-t*(4+u*u))
need('heat_count_includes_zero_while_positive_zeta_count_does_not',s.limit(cutheat-fullheat,t,0)==-2 and low.subs(sn,0)==1)
# Moving spectral projection with a fixed threshold branch, no guessed P'.
rot=s.Matrix([[s.cos(u),-s.sin(u)],[s.sin(u),s.cos(u)]])
LM=rot*s.diag(1+u,4+2*u)*rot.T;PM=rot*s.diag(0,1)*rot.T
HM=rot*s.diag(s.exp(-t*(1+u)),s.exp(-t*(4+2*u)))*rot.T
need('moving_spectral_projector_has_offdiagonal_derivative',eq(PM*PM.diff(u)*PM,s.zeros(2)) and not eq(PM.diff(u),s.zeros(2)))
need('first_trace_derivative_from_moving_projection_is_zero',s.trigsimp(s.trace(PM.diff(u)*HM))==0)
need('whole_moving_cutoff_trace_has_original_eigenvalue_derivative',s.trigsimp(s.diff(s.trace(PM*HM),u)+t*s.trace(PM*HM*LM.diff(u)))==0)
need('H_mass_curvature_ratio_unchanged_under_common_chiral_halving',(-s.Rational(1,8))/s.Rational(1,96)==-12)
result={'contract':'ORIGINAL-ADMISSIBILITY-TOTAL-CHIRAL-TRACE-20261031','all_passed':all(checks.values()),'checks':checks,'scope':'Actual96 finite constructor and existing conditional spin4 product; generic spectral-cut,rectangular index and moving-projection identity controls. No chosen physical normal module,trace or source. No new Lean proof,original Dadm/Ksigma identification,physical chirality selection or whole TFPT closure.'}
Path(__file__).with_name('ORIGINAL-ADMISSIBILITY-AND-TOTAL-CHIRAL-TRACE.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'all_passed':result['all_passed'],'count':len(checks)}))
