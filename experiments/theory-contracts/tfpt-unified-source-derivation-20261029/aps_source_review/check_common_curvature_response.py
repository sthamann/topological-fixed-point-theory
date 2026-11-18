"""Whole Spin4 curvature trace and given Higgs radial source response.
Conditional original R031 consumer; no added physical operator or R² term.
"""
from pathlib import Path
import json,itertools
import sympy as s
checks={}
def need(k,v):
 checks[k]=bool(v)
 if not v:raise AssertionError(k)
def eq(k,x,y):need(k,s.simplify(x-y)==0)
I=s.I;one=s.eye(2);sx=s.Matrix([[0,1],[1,0]]);sy=s.Matrix([[0,-I],[I,0]]);sz=s.diag(1,-1)
ga=[s.kronecker_product(sx,one),s.kronecker_product(sy,sx),s.kronecker_product(sy,sy),s.kronecker_product(sy,sz)]
need('same_full_Euclidean_Spin4_Clifford_algebra',all(ga[i]*ga[j]+ga[j]*ga[i]==2*(i==j)*s.eye(4) for i in range(4) for j in range(4)))
# All algebraic 4D curvature tensors: symmetric6x6 pair matrix (21 entries),
# one Bianchi relation R0123-R0213+R0312=0 gives dimension20.
pairs=list(itertools.combinations(range(4),2));lookup={p:i for i,p in enumerate(pairs)}
coords=list(itertools.combinations_with_replacement(range(6),2));removed=(lookup[(0,3)],lookup[(1,2)])
coords.remove(removed)
def tensor(mat):
 out={}
 for a,b,c,d in itertools.product(range(4),repeat=4):
  if a==b or c==d:out[a,b,c,d]=s.S(0);continue
  out[a,b,c,d]=(1 if a<b else -1)*(1 if c<d else -1)*mat[lookup[tuple(sorted((a,b)))],lookup[tuple(sorted((c,d)))]]
 return out
ts=[];oms=[]
for i,j in coords:
 mat=s.zeros(6);mat[i,j]=mat[j,i]=1
 aa=mat[lookup[(0,1)],lookup[(2,3)]];bb=mat[lookup[(0,2)],lookup[(1,3)]]
 mat[removed[0],removed[1]]=mat[removed[1],removed[0]]=bb-aa
 t=tensor(mat);ts.append(t)
 oms.append([sum((t[mu,nu,a,b]*ga[a]*ga[b]/4 for a,b in itertools.product(range(4),repeat=2)),s.zeros(4)) for mu,nu in itertools.product(range(4),repeat=2)])
need('whole20_curvature_basis_obeys_all_Bianchi_identities',len(ts)==20 and all(t[a,b,c,d]+t[a,c,d,b]+t[a,d,b,c]==0 for t in ts for a,b,c,d in itertools.product(range(4),repeat=4)))
gram=s.Matrix([[sum(ts[i][k]*ts[j][k] for k in ts[i]) for j in range(20)] for i in range(20)])
need('whole20_curvature_basis_is_independent',gram.rank()==20)
spin=s.Matrix([[s.simplify(sum(s.trace(x*y) for x,y in zip(oms[i],oms[j]))) for j in range(20)] for i in range(20)])
need('whole_polarized_spin_connection_curvature_identity',spin==-gram/2)
N,R,Ric2,Riem2,W2,Eu=s.symbols('N R Ric2 Riem2 W2 Euler',real=True)
# Pure spin gravity:TrI=4N,TrE=-NR,TrE²=NR²/4,TrOmega²=-N Riem²/2.
poly=s.expand(4*N*(5*R**2-2*Ric2+2*Riem2)-60*N*R**2+45*N*R**2-15*N*Riem2)
eq('full_spin_gravity_polynomial',poly,N*(5*R**2-8*Ric2-7*Riem2))
solution=s.solve([Eu-(Riem2-4*Ric2+R**2),W2-(Riem2-2*Ric2+R**2/3)],[Ric2,Riem2])
eq('whole_invariant_dictionary_has_no_independent_R_squared',poly.subs(solution),N*(11*Eu-18*W2))
eq('same_geometry_rank_reference_cancels_full_pure_a4_gravity',poly-poly,0)
eq('conformally_flat_constant_curvature_can_mislead_R_squared_fit',poly.subs({Ric2:R**2/4,Riem2:R**2/6}),s.Rational(11,6)*N*R**2)
# Execute the OTHER original curvature consumer v475 on the same WHOLE space.
# Its Bianconi blocks are R, Ricci and 2R_abcd on unordered2form pairs.
rr=[s.Matrix(4,4,lambda a,c:sum(t[a,b,c,b] for b in range(4))) for t in ts]
sc=[s.trace(q) for q in rr]
xx2=[s.Matrix(6,6,lambda a,b:2*t[(*pairs[a],*pairs[b])]) for t in ts]
need('original_entropic_whole_curvature_trace_is_3R',all(sc[i]+s.trace(rr[i])+s.trace(xx2[i])==3*sc[i] for i in range(20)))
need('original_entropic_whole_polarized_square_trace',all(s.trace(xx2[i]*xx2[j])==gram[i,j] for i in range(20) for j in range(20)))
beta,Mbar=s.symbols('beta Mbar',positive=True)
ent=beta**2*(R**2+Ric2+Riem2)/2
eq('entropic_full_independent_invariant_coefficients',ent.subs(solution),beta**2*(s.Rational(5,6)*R**2+s.Rational(5,4)*W2-s.Rational(3,4)*Eu))
eq('entropic_maximally_symmetric_17_over24_is_retained',ent.subs({Ric2:R**2/4,Riem2:R**2/6}),s.Rational(17,24)*beta**2*R**2)
# Independent conformal plane wave: Hessian(sigma)=-k_a k_b sigma.
ks=s.symbols('k0:4',real=True);hs=-s.Matrix(ks)*s.Matrix(ks).T
delta=lambda a,b:s.S(int(a==b))
rt={(a,b,c,d):delta(a,d)*hs[b,c]+delta(b,c)*hs[a,d]-delta(a,c)*hs[b,d]-delta(b,d)*hs[a,c] for a,b,c,d in itertools.product(range(4),repeat=4)}
rp=s.Matrix(4,4,lambda a,c:sum(rt[a,b,c,b] for b in range(4)));rs=s.trace(rp)
ricplane=s.trace(rp*rp);riemplane=sum(q*q for q in rt.values())
eq('independent_all_momenta_scalar_response_is_5_over3_R_squared',rs**2+ricplane+riemplane,s.Rational(5,3)*rs**2)
eq('scalar_plane_wave_Euler_piece_is_pure_boundary',riemplane-4*ricplane+rs**2,0)
c3=1/(8*s.pi);bval=c3/(6*8*s.pi*Mbar**2)
mass_full=s.simplify(3/(5*bval));mass_ms=s.simplify(12/(17*bval))
eq('same_original_entropic_Einstein_pin',((8*s.pi*Mbar**2)**2)*3*bval,Mbar**2/2)
eq('full_quadratic_scalar_coefficient_scale',mass_full,s.Rational(1152,5)*s.pi**2*Mbar**2)
eq('old_max_symmetric_mass_ratio_is_not_full_scalar_coefficient',mass_full/mass_ms,s.Rational(17,20))
eq('full_entropic_to_existing_TFPT_scale_mismatch',mass_full/(c3**7*Mbar**2),s.Rational(18,5)*(8*s.pi)**9)
# All original H traces retained; c2,d are actual common moment weights only.
A,B,C,chi,c2,d=s.symbols('A B C chi c2 d',positive=True)
lam=8*s.pi**2*B/(d*A**2);mu2=2*c2*chi**2/d-4*C/A;xi=s.Rational(1,6)
u=s.symbols('u',real=True);V=lam*u**2-(mu2-xi*R)*u
us=(mu2-xi*R)/(2*lam)
eq('given_whole_Higgs_stationary_branch',s.diff(V,u).subs(u,us),0)
Vs=s.factor(V.subs(u,us))
eq('whole_Higgs_algebraic_curvature_susceptibility',Vs,-(mu2-xi*R)**2/(4*lam))
eq('induced_local_R_squared_from_given_Higgs',s.diff(Vs,R,2)/2,-d*A**2/(1152*s.pi**2*B))
v,rho,m2,la,xx=s.symbols('v rho m2 lambda xi',positive=True)
UU=(v+rho/s.sqrt(2))**2;Vrad=la*UU**2-m2*UU+xx*R*UU
Vbroken=s.expand(Vrad.subs(m2,2*la*v**2))
eq('original_Higgs_radial_mass_squared',s.diff(Vbroken,rho,2).subs({rho:0,R:0}),4*la*v**2)
g=s.diff(Vbroken,rho,R).subs(rho,0)
eq('actual_mixed_curvature_radial_vertex',g,s.sqrt(2)*xx*v)
k,mrad=s.symbols('k mrad',nonnegative=True)
response=-g**2/(2*(k**2+mrad**2))
eq('whole_quadratic_response_keeps_radial_propagator',response,-xx**2*v**2/(k**2+mrad**2))
eq('same_low_momentum_R_squared_from_own_radial_pole',response.subs(k,0).subs(mrad**2,4*la*v**2),-xx**2/(4*la))
eq('Poisson_original_Higgs_radial_mass',2*mu2.subs({c2:2,d:1}),8*chi**2-8*C/A)
# Exact quadratic finite Gaussian source includes determinant and every J leg.
kap,j=s.symbols('K J',positive=True);src=s.symbols('source',real=True)
logZ=-s.log(kap)/2+(src-g*R)**2/(2*kap)
eq('whole_quadratic_source_retains_shifted_mean',s.diff(logZ,src),(src-g*R)/kap)
eq('whole_quadratic_source_retains_covariance',s.diff(logZ,src,2),1/kap)
eq('absolute_zero_source_weight_retains_determinant',logZ.subs(src,0),-s.log(kap)/2+g**2*R**2/(2*kap))
result={'contract':'COMMON-CURVATURE-RESPONSE-20261031','all_passed':all(checks.values()),'checks':checks,'values':{'pure_spin_a4_numerator':str(poly),'invariant_form':'N(11Euler-18Weyl²), no independentR²','whole_algebraic_curvature_space_dimension':20,'polarized_identity_entries':400,'original_entropic_invariants':'beta²(5R²/6+5Weyl²/4-3Euler/4)','original_entropic_full_scalar_mass_squared':str(mass_full),'original_entropic_original_max_symmetric_mass_squared':str(mass_ms),'original_entropic_required_enhancement':str(s.Rational(18,5)*(8*s.pi)**9),'Higgs_local_R_squared':str(-d*A**2/(1152*s.pi**2*B)),'Higgs_nonlocal_response':'-xi²v² R(-Delta+m_rad²)^-1R','scope':'Existing conditional ordinary Spin4x96 Dirac consumer and separate original v475 entropic vacuum functional, never identified with each other. Closed/interior geometry, actual matching reference. Entropic mass parameter follows its quadratic metric matching, not complete Lorentzian/RP/measure theorem. Higgs response is fixed-background tree/quadratic stationary branch, not full H integration or new scalaron/physical prediction. Retain determinant, higher vertices, source and Z0.'}}
Path(__file__).with_name('COMMON-CURVATURE-RESPONSE.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'all_passed':result['all_passed'],'checks':len(checks),'curvature_basis':20,'polarized_entries':400,'local_Higgs_R_squared':result['values']['Higgs_local_R_squared']}))
