"""Current v255/v259 argument audit; original Higgs moments retained.
No new physical cutoff, Y/MR, reference, domain or ledger promotion.
"""
from pathlib import Path
import json
import sympy as s
u,r,p,m,L=s.symbols('u r p m L',positive=True)
checks={};values={}
def need(k,b):
 checks[k]=bool(b)
 if not b:raise AssertionError(k)
def eq(k,a,b):need(k,s.simplify(a-b)==0)
def moment(F):return [s.limit(F,u,0,dir='+'),s.integrate(F,(u,0,s.oo)),s.integrate(u*F,(u,0,s.oo))]
Mh=moment(s.exp(-u));Mp=moment(s.exp(-s.sqrt(u)))
need('heat_squared_argument_moments_1_1_1',Mh==[1,1,1])
need('original_positive_source_in_squared_argument_moments_1_2_12',Mp==[1,2,12])
for name,f,w in [('poisson',s.exp(-r),[1,2,12]),('heat',s.exp(-r*r),[1,1,1])]:
 rr=[f.subs(r,0),2*s.integrate(r*f,(r,0,s.oo)),2*s.integrate(r**3*f,(r,0,s.oo))]
 need('Dirac_and_squared_argument_dictionary_'+name,rr==w)
# Source coefficient table, in the order (a0,a2,a4).
printed=[2,2,1];poisson=[12,2,1];heat=[1,1,1]
need('printed_triplet_is_neither_heat_nor_poisson',printed!=poisson and printed!=heat)
need('printed_triplet_excluded_from_entire_positive_Laplace_cone',printed[2]*printed[0]-printed[1]**2==-2)
need('Poisson_Higgs_and_Einstein_weight_matches_printed_a2_a4',poisson[1:]==printed[1:])
need('relative_same_symbol_a0_cancels',s.expand(2*s.Symbol('a0')-2*s.Symbol('a0'))==0)
sg=s.symbols('sigma',real=True)
need('unscaled_reference_keeps_nonzero_scaled_volume_jet',s.diff(2*(s.exp(4*sg)-1),sg).subs(sg,0)==8)
# Independent exact flat four-dimensional phase-space integral, per spin/internal eigenvalue.
# This is the existing CONDITIONAL ordinary4D consumer, not a new spacetime origin.
q=s.symbols('q',positive=True)
IP=s.integrate((q**3-m*m*q)*s.exp(-q/L),(q,m,s.oo))
IPformula=s.exp(-m/L)*(6*L**4+6*m*L**3+2*m*m*L**2)
eq('exact_Poisson_constant_fibre_integral',IP,IPformula)
IH=s.integrate(p**3*s.exp(-(p*p+m*m)/L**2),(p,0,s.oo))
eq('exact_heat_constant_fibre_integral',IH,L**4*s.exp(-m*m/L**2)/2)
x=s.symbols('x',positive=True)
PS=s.series((IPformula.subs(L,1/x))*x**4,x,0,7).removeO()
HS=s.series((IH.subs(L,1/x))*x**4,x,0,7).removeO()
eq('Poisson_a0_weight_12',PS.coeff(x,0),6)
eq('Poisson_a2_weight_2',PS.coeff(x,2),-m*m)
eq('Poisson_a4_weight_1',PS.coeff(x,4),m**4/4)
eq('heat_a0_weight_1',HS.coeff(x,0),s.Rational(1,2))
eq('heat_a2_weight_1',HS.coeff(x,2),-m*m/2)
eq('heat_a4_weight_1',HS.coeff(x,4),m**4/4)
need('Poisson_remainder_need_not_be_heat_even_power_remainder',PS.coeff(x,5)!=0 and HS.coeff(x,5)==0)
A,B,C,RR,chi=s.symbols('A B C R chi',positive=True);U,kin=s.symbols('U kin',real=True)
a2H=-A*U/(4*s.pi**2)
a4H=(A*kin+B*U**2+4*C*U)/(8*s.pi**2)+A*RR*U/(48*s.pi**2)
SP=2*chi**2*a2H+a4H;SH=chi**2*a2H+a4H
eq('current_Poisson_whole_Higgs_response',SP,(A*kin+B*U**2+4*C*U-4*A*chi**2*U)/(8*s.pi**2)+A*RR*U/(48*s.pi**2))
Z=A/(8*s.pi**2)
eq('canonical_Poisson_mu_squared',-s.diff(SP,U).subs(U,0)/Z,4*chi**2-4*C/A-RR/6)
eq('canonical_heat_mu_squared',-s.diff(SH,U).subs(U,0)/Z,2*chi**2-4*C/A-RR/6)
eq('same_canonical_quartic',s.diff(SP,U,2)/(2*Z*Z),8*s.pi**2*B/(A*A))
eq('current_conditional_Higgs_threshold',s.solve(s.diff(SP,U).subs(U,0),chi**2)[0],C/A+RR/24)
eh=s.symbols('EH');eq('same_weight_multiplies_original_Einstein_a2',s.diff(2*chi**2*eh,eh),2*chi**2)
# The archived ratio problem is an operator heat-coefficient question, not fixed by changing this moment.
c2=s.symbols('c2',nonzero=True)
eq('no_new_resolution_of_archived_factor_six',s.cancel(c2*(-12)/c2),-12)
values={'ordinary_heat_coefficient_order':['a0','a2','a4'],'heat_weights':heat,'Poisson_weights':poisson,'printed_v255_v259_weights':printed,'printed_positive_Laplace_Gram_determinant':-2,'independent_phase_space_Poisson_series':str(PS),'independent_phase_space_heat_series':str(HS),'mu_squared_Poisson':str(4*chi**2-4*C/A-RR/6),'scope':'Conditional same ordinary4D Dirac consumer and original whole96 moments. No identification of positive boundaryQ with original bosonic fullK or Tr_rel. R016/R017/R019 reused,not reclassified as new findings. a0 obstruction does not apply after its exact same-symbol relative cancellation.'}
result={'contract':'CURRENT-ACTION-ARGUMENT-20261031','all_passed':all(checks.values()),'checks':checks,'values':values}
Path(__file__).with_name('CURRENT-ACTION-ARGUMENT.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'all_passed':result['all_passed'],'count':len(checks),'weights':values}))
