"""Original96 Higgs inner-fluctuation/scale map and finite source transport.
Standalone scoped controls; no ledger, new physical Y, or chosen Brel origin.
"""
from pathlib import Path
import sys,json
import numpy as np
import sympy as s
ROOT=Path(__file__).resolve().parents[4]
sys.path.insert(0,str(ROOT/'verification'))
import v252_full_finite_triple as original
checks={}
def need(name,b):
 checks[name]=bool(b)
 if not b:raise AssertionError(name)
zero3=np.zeros((3,3),complex);zero2=np.zeros((2,2),complex)
ys=[np.array([[1+1j,2,0],[0,1,1j],[1,-1j,2]],complex),np.array([[2,1j,1],[1,3,-1j],[0,1,1]],complex),np.array([[1,0,1j],[2,1,0],[0,-1j,2]],complex),np.array([[1,1j,0],[1,2,1],[0,1,3j]],complex)]
mr=np.array([[2,1j,1],[1j,3,2],[1,2,1j]],complex)
DY=original.full_D(*ys,zero3);DM=original.full_D(zero3,zero3,zero3,zero3,mr)
eR=original.rep(1,zero2,zero3);K=original.Kswap()
def L(h):
 h1,h2=np.asarray(h,dtype=complex);q=np.array([[h2.conjugate(),h1],[-h1.conjugate(),h2]])
 a=original.rep(0,q,zero3);one=a@(DY@eR-eR@DY);one=one+one.conj().T
 return one+K@one.conj()@K
hs=np.array([0,1],complex);h=np.array([1+2j,3-1j],complex)
def phi(h):return L(h)+DM
def change(h):return L(h-hs)
need('actual_unfluctuated_H_star_is_original_DY',np.array_equal(L(hs),DY))
need('actual_total_H_zero_keeps_Majorana_only',np.array_equal(phi([0,0]),DM))
need('original_oneform_change_is_total_minus_seed',np.array_equal(change(h),L(h)-DY))
need('seed_plus_change_is_whole_original_Higgs_operator',np.array_equal(DY+DM+change(h),phi(h)))
need('change_is_self_adjoint_and_contains_no_second_Majorana',np.array_equal(change(h),change(h).conj().T) and np.trace(change(h)@DM)==0)
for k,v in enumerate([np.array([1,0]),np.array([1j,0]),np.array([0,1]),np.array([0,1j])]):
 need('four_real_actual_H_directions_'+str(k),np.array_equal(change(h+v)-change(h),L(v)))
q=s.Rational(1,2);qf=float(q)
scaled=qf*phi(h);inner=qf*(DY+DM)+qf*change(h)
need('scaled_seed_and_scaled_inner_change_give_full_scaled_H_family',np.array_equal(inner,scaled))
need('adding_total_Higgs_again_double_counts_seed',not np.array_equal(qf*(DY+DM)+qf*phi(h),scaled))
outside=qf*(DY+DM)+change(h)
h_sigma=hs+(h-hs)/qf
need('outside_unscaled_change_is_exact_affine_pullback_of_same_family',np.array_equal(outside,qf*phi(h_sigma)))
need('outside_at_same_named_H_is_different_operator',not np.array_equal(outside,scaled))
need('affine_inverse_recovers_original_H',np.array_equal(hs+qf*(h_sigma-hs),h))
need('original_H_zero_is_not_zero_in_affine_coordinates',not np.array_equal(hs+(np.zeros(2)-hs)/qf,np.zeros(2)))
# Variable f, variable a,b and vector: exact differential commutator identity.
x=s.symbols('x',real=True);ff=s.exp(-s.Function('sigma')(x)/2)
bm=s.Matrix(2,2,[s.Function('b'+str(i))(x) for i in range(4)])
vv=s.Matrix([s.Function('v0')(x),s.Function('v1')(x)])
mass=s.Matrix([[1,s.I],[-s.I,2]])
def d(v):return -s.I*v.diff(x)+mass*v
def ds(v):return ff*d(ff*v)
comm=ds(bm*vv)-bm*ds(vv)-ff*(d(bm*(ff*vv))-bm*d(ff*vv))
need('variable_scale_inner_commutator_has_no_extra_df_term',all(s.simplify(z)==0 for z in comm))
# General four-real-field affine map, source and finite measure transport.
q0=s.symbols('q',positive=True);hv=s.Matrix(s.symbols('h0:4',real=True));hstar=s.Matrix([0,0,1,0]);j=s.Matrix(s.symbols('j0:4',real=True))
hh=hstar+(hv-hstar)/q0;back=hstar+q0*(hh-hstar)
need('whole_four_real_affine_inverse',s.simplify(back-hv)==s.zeros(4,1))
jac=hh.jacobian(hv).det()
need('four_real_jacobian_is_q_minus4',s.simplify(jac-q0**-4)==0)
need('source_is_qj_plus_required_affine_constant',s.simplify((j.T*hv)[0]-(q0*j.T*hh)[0]-(1-q0)*(j.T*hstar)[0])==0)
need('two_site_full_measure_jacobian_includes_both_sites',s.det(s.diag(q0*s.eye(4),s.Rational(3,2)*s.eye(4)))==q0**4*s.Rational(3,2)**4)
# A genuine convergent integral control. This Gaussian is ONLY a change-of-
# variables check, not a proposed TFPT boson action or physical measure.
z,t,js=s.symbols('z t js',real=True);a=s.symbols('a',positive=True)
hh1=1+(t-1)/q0;U=lambda v:a*v*v/2
transformed=s.expand(U(q0*z+1-q0)-q0*js*z)
az=s.diff(transformed,z,2);bz=s.diff(transformed,z).subs(z,0);cz=transformed.subs(z,0)
oldlog=js**2/(2*a);newlog=s.log(q0)+(1-q0)*js-s.log(az/a)/2+bz*bz/(2*az)-cz
need('convergent_integral_whole_source_and_Z0_transport',s.simplify(s.expand(newlog-oldlog))==0)
need('discarded_affine_source_factor_changes_mean',s.diff((1-q0)*js,js)==1-q0)
# Bare color3+weak2 module has no color-neutral weak-doublet endomorphism.
omega=(-1+s.sqrt(3)*s.I)/2
color=s.diag(omega,omega,omega,1,1);weak=s.diag(1,1,1,-1,-1)
cc=s.symbols('X0:25');X=s.Matrix(5,5,cc)
ec=list(color*X*color.inv()-X);ew=list(weak*X*weak+X)
lc=s.linear_eq_to_matrix(ec,cc)[0];lw=s.linear_eq_to_matrix(ew,cc)[0]
need('color_neutral_endomorphisms_are_diagonal_blocks_dim13',25-lc.rank()==13)
need('weak_center_odd_endomorphisms_are_cross_blocks_dim12',25-lw.rank()==12)
need('whole_bare_C5_Higgs_doublet_intertwiner_space_is_zero',s.linear_eq_to_matrix(ec+ew,cc)[0].rank()==25)
# A kernel sector and its operator are different objects, independent of
# any profile choice: algebraic PBP=0 follows BP=0 for a spectral kernel.
P=s.diag(1,1,0);bn=s.diag(0,0,7)
need('normal_zero_mode_projection_does_not_supply_a_Higgs_mass_vertex',bn*P==s.zeros(3) and P*bn*P==s.zeros(3))
ortho=s.Matrix([[1,2,2],[2,1,-2],[-2,2,-1]])/3
SS=ortho[:,:2];BB=ortho*s.diag(0,0,7)*ortho.T;fs=s.diag(s.Rational(1,2),s.Rational(1,3),s.Rational(1,5))
raw=fs.inv()*SS;GG=raw.T*raw;chol=GG.cholesky();Cnorm=chol.T.inv();TT=raw*Cnorm
need('whole_kernel_transports_by_inverse_f',fs*BB*fs*raw==s.zeros(3,2) and (fs*BB*fs).rank()==1)
need('kernel_Gram_is_compressed_q_inverse',GG==SS.T*fs.inv()**2*SS and GG.det()>0)
need('variable_normal_scale_Gram_need_not_be_scalar',GG[0,1]!=0 and GG[0,0]!=GG[1,1])
need('transported_mode_normalization_is_exact',s.simplify(TT.T*TT)==s.eye(2))
XX=s.Matrix([[1,2,0],[2,3,1],[0,1,4]])
need('whole_scaled_vertex_Gram_compression',s.simplify(TT.T*(fs*XX*fs)*TT-Cnorm.T*SS.T*XX*SS*Cnorm)==s.zeros(2))
need('whole_unscaled_vertex_Gram_compression',s.simplify(TT.T*XX*TT-Cnorm.T*SS.T*fs.inv()*XX*fs.inv()*SS*Cnorm)==s.zeros(2))
need('complex_two_mode_coefficient_Jacobian_is_inverse_det_Gram',s.simplify(GG.det()*Cnorm.det()**2)==1)
result={'contract':'ORIGINAL-HIGGS-INNER-SCALE-AND-SOURCE-20261031','all_passed':all(checks.values()),'checks':checks,'scope':'Generic central scaling and original96 oneform family; exact affine finite source transport. Bare C5 representation test excludes only direct H-linear endomorphism there. No Brel/index identification,chosen normal reduction,physical measure,or ledger promotion. Gaussian and three-slot matrices are identity diagnostics only.'}
Path(__file__).with_name('ORIGINAL-HIGGS-INNER-SCALE.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
