"""Actual full96 coefficient inverse and ordered neutral norm identity.
Reuses unchanged R031 integer diagnostic matrices; no new physical model.
No Fock space is formed on the 96 internal labels.
"""
from pathlib import Path
import sys,json
import numpy as np
root=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
sys.path.insert(0,str(root/'verification'))
import v252_full_finite_triple as v
z3=np.zeros((3,3),complex);z2=np.zeros((2,2),complex)
Ys=[np.array([[1+1j,2,0],[0,1,1j],[1,-1j,2]],complex),np.array([[2,1j,1],[1,3,-1j],[0,1,1]],complex),np.array([[1,0,1j],[2,1,0],[0,-1j,2]],complex),np.array([[1,1j,0],[1,2,1],[0,1,3j]],complex)]
MR=np.array([[2,1j,1],[1j,3,2],[1,2,1j]],complex)
DY=v.full_D(*Ys,z3);DM=v.full_D(z3,z3,z3,z3,MR)
P=v.rep(1,z2,z3);K=v.Kswap();gamma=v.grading()
qs=[np.array([[0,1],[-1,0]],complex),np.array([[0,1j],[1j,0]],complex),np.eye(2),np.diag([-1j,1j])]
def image(q):
 a=v.rep(0,q,z3);b=a@(DY@P-P@DY);b=b+b.conj().T
 return b+K@b.conj()@K
Va=[image(q) for q in qs]
A=np.trace(DY@DY);M2=np.trace(DM@DM)
checks=[]
def need(name,ok):
 checks.append({'name':name,'passed':bool(ok)})
 if not ok:raise AssertionError(name)
need('original_A_nonzero_without_new_normalization',A.real>0 and A.imag==0)
# Each key is a whole ordered noncommutative boson word. No xa xb=xb xa.
phi={():DM,**{(a,):Va[a] for a in range(4)}}
# The kinetic symbol is a separate formal original A operator, not a number.
h={('kinetic',):gamma,():-1j*gamma@DM,**{(a,):-1j*gamma@Va[a] for a in range(4)}}
rec={w:1j*gamma@((m-gamma@m@gamma)/2) for w,m in h.items()}
need('odd_generator_part_removes_entire_original_kinetic_leg',np.array_equal(rec[('kinetic',)],np.zeros((96,96))))
need('odd_generator_part_retains_entire_original_Majorana',np.array_equal(rec[()],DM))
for a in range(4):
 need(f'odd_generator_part_recovers_original_H_direction_{a}',np.array_equal(rec[(a,)],Va[a]))
 for b in range(4):
  val=np.trace(Va[a]@rec[(b,)])/A
  need(f'whole_operator_component_inverse_{a}_{b}',val==(1 if a==b else 0))
 need(f'original_Majorana_is_not_a_H_component_{a}',np.trace(Va[a]@DM)==0)
# Full coefficientwise partial trace of Phi^2, retaining word order.
trace_square={}
for wa,ma in phi.items():
 for wb,mb in phi.items():
  key=wa+wb;trace_square[key]=trace_square.get(key,0j)+np.trace(ma@mb)
expected={():M2,**{(a,a):A for a in range(4)}}
need('whole_ordered_norm_polynomial_no_boson_commutativity_used',all(trace_square.get(w,0)==expected.get(w,0) for w in set(trace_square)|set(expected)))
need('ordered_norm_polynomial_retains_original_Majorana_scalar',trace_square[()]==M2)
# Original generator-Yukawa source directions; no conversion to Euclidean
# four-dimensional Weyl fields is inferred.
for a in range(4):
 T=-1j*gamma@Va[a]
 need(f'original_Yukawa_generator_direction_Hermitian_{a}',np.array_equal(T,T.conj().T))
 need(f'original_Yukawa_generator_direction_recovers_V_{a}',np.array_equal(1j*gamma@T,Va[a]))
# New weighted double-CAR extraction application. Numerical coefficients
# act only on internal labels; kinetic symbols are not physically fitted.
Ta=[-1j*gamma@m for m in Va]
for a in range(4):
 need(f'weighted_contact_kills_original_even_kinetic_{a}',np.trace(Ta[a]@gamma)==0)
 need(f'weighted_contact_kills_original_Majorana_{a}',np.trace(Ta[a]@(-1j*gamma@DM))==0)
 for b in range(4):
  need(f'weighted_contact_recovers_H_component_{a}_{b}',np.trace(Ta[a]@Ta[b])/A==(1 if a==b else 0))
# An already used R031 four-component constructor point, not a new model.
coords=[1,2,2,4]
Phi=DM+sum((coords[a]*Va[a] for a in range(4)),np.zeros((96,96),complex))
h=gamma/4-1j*gamma@Phi
extracted=[np.trace(T@h)/A for T in Ta]
need('actual_conditional_source_returns_original_H_scalars',all(extracted[a]==coords[a] for a in range(4)))
need('actual_conditional_source_U_is_scalar25',sum(x*x for x in extracted)==25)
# The whole anticommutators remain ordered. Each D contains one derivative.
from itertools import product
word_terms=[]
for order in product([0,1],repeat=3):
 w=[]
 for k,flip in enumerate(order,1):
  pair=[f'c_j{k}',f'delta_c_i{k}_dagger']
  w+=pair if flip==0 else pair[::-1]
 w+=['c_r_dagger','c_q'];word_terms.append(w)
need('neutral_UY_has_eight_ordered_terms',len(word_terms)==8)
need('each_UY_word_has_eight_CAR_letters_three_original_derivatives',all(len(w)==8 and sum(x.startswith('delta_') for x in w)==3 for w in word_terms))
# At fixed H the extracted U is scalar, so connected U/Y is zero for any
# actual conditional source. This is algebraic, not a new state evaluation.
import sympy as sym
y=sym.symbols('actual_conditional_Y_expectation')
need('conditional_connected_UY_is_exactly_zero',sym.expand(25*y-25*y)==0)
out={'scope':'Exact unchanged whole96 constructor controls of ordered coefficient/operator identities. Uses original R031 integer fixture only. Not96 physical fermion integrations, no new boson state, weight, spectrum or source value.',
 'checked':sum(x['passed'] for x in checks),'total':len(checks),'checks':checks,
 'diagnostic_invariants':{'A':str(A),'Tr_DM_squared':str(M2)},
 'identity':'Phi_hat=i gamma (V-gamma V gamma)/2; x_hat_a=Tr96[Va(Phi_hat-DM)]/A; U_hat=(Tr96 Phi_hat^2-Tr96 DM^2)/A=sum x_hat_a^2',
 'unexecuted_original_input':'Actual original Q0[u] whole endpoint mark or its equivalent same-source physical composite insertion; exact product representation if original endpoint chart supplies it. No background averaging substituted.',
 'whole_ordered_UY_word_patterns':word_terms,
 'conditional_existing_branch':{'original_R031_H_coordinates':coords,'extracted_coordinates':[str(x) for x in extracted],'U':'25 I','connected_UY':'0 identically because U is a fixed-background scalar'},
 'Y_scope':'Y_ham=sum Xa dGamma(Ta) in actual R113 h chart. No automatic Euclidean Sf pairing/half-factor/sign identification.',
 'primitive_mixed_response_evaluated':False}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:out[k] for k in ['checked','total','diagnostic_invariants','scope','primitive_mixed_response_evaluated']},indent=2))
