"""Source-derived, exact reachable process for original return and native J.

No new Hamiltonian or physical history is chosen. The alphabet is explicitly
the already existing single-beta response L and original Cartan quarter J.
"""
from pathlib import Path
import sys,json,time
import sympy as s
ROOT=Path(__file__).resolve().parents[4]
sys.path.insert(0,str(ROOT/'experiments/theory-contracts/joint-native-algebra-selection-20261016'))
import check_source_return_memory as src
m=src.m;Q=m.Q
start=time.time()
b=[Q.zero,Q.zero,Q(0,-1),Q.zero];c=[Q.zero,Q.zero,Q.zero,Q(4,-4)]
source=src.full_source(b,c);neut=src.basis(4);inverse={}
for mon in neut:
 for n,value in source({mon:Q.one}).items():
  inverse.setdefault(n,[]).append((mon,m.conj(value)*Q(m.weight(n))/m.weight(mon)))
def adj(f):
 out={}
 for n,value in f.items():
  for mon,z in inverse.get(n,[]):out[mon]=out.get(mon,Q.zero)+z*value
 return {mon:z for mon,z in out.items() if z}
def L(f):return adj(source(f))
def clock(f,inverse=False):
 return {mon:a*Q(0,1)**((-1 if inverse else 1)*sum(1 if k%8<4 else -1 for k in mon)) for mon,a in f.items()}
ops=[L,clock,lambda f:clock(f,True)]
polys,nn=m.code();vectors=[];norms=[];provenance=[]
def add_vector(v,label):
 for u,n in zip(vectors,norms):
  z=m.inner(u,v)/n
  if z:v=m.add((Q.one,v),(-z,u))
 if not v:return False
 norm=m.inner(v,v);m.need(norm!=Q.zero,'Nonzero exact positive Fock norm')
 vectors.append(v);norms.append(norm);provenance.append(label)
 return True
for j,p in enumerate(polys):add_vector(p,['code',j])
# Source-derived first exterior layer, retained for the nested-elimination test.
for j,p in enumerate(polys):add_vector(L(p),['first_return',j])
first_layer=len(vectors)
index=0
while index<len(vectors):
 for a,op in enumerate(ops):add_vector(op(vectors[index]),[a,index])
 index+=1
 if index%10==0:print('processed',index,'actual reachable',len(vectors),'seconds',round(time.time()-start,2),flush=True)
d=len(vectors);G=s.diag(*[Q.to_sympy(n) for n in norms]);matrices=[];checks={}
for a,op in enumerate(ops):
 cols=[]
 for j,v in enumerate(vectors):
  out=op(v);coords=[m.inner(u,out)/n for u,n in zip(vectors,norms)]
  remainder=m.add((Q.one,out),*[(-z,u) for z,u in zip(coords,vectors) if z])
  m.need(not remainder,'Exact original-vector invariant closure')
  cols.append(s.Matrix([Q.to_sympy(x) for x in coords]))
 matrices.append(s.Matrix.hstack(*cols))
A,C,Ci=matrices
checks['actual_full726_source_used']=len(neut)==726
checks['actual_reachable_Gram_strictly_positive']=all(Q.to_sympy(n)>0 for n in norms)
checks['first_source_layer_has_dimension10']=first_layer==10
checks['original_L_selfadjoint_in_actual_metric']=A.H*G==G*A
checks['native_J_unitary_in_actual_metric']=C.H*G*C==G
checks['native_J_fourth_power_identity']=C**4==s.eye(d)
checks['native_J_adjoint_kept']=Ci*C==s.eye(d) and C.H*G==G*Ci
checks['original_code_native_J_identity']=C[:,:5]==s.eye(d)[:,:5]
checks['all_original_generators_exactly_closed']=True
def schur(M,k):
 return M[:k,:k]-M[:k,k:]*M[k:,k:].inv(method='DM')*M[k:,:k] if k<len(M) else M
# A negative spectral point gives invertible compressions since L>=0.
M=s.eye(d)+A
direct=schur(M,5);stage1=schur(M,first_layer);nested=schur(stage1,5)
checks['two_stage_equals_direct_full_source_elimination']=nested==direct
checks['coarse_inverse_is_actual_full_source_resolvent']=direct.inv(method='DM')==M.inv(method='DM')[:5,:5]
checks['memory_term_is_nonzero']=direct!=M[:5,:5]
# Clock-bearing ordered histories, not just scalar repetitions of L.
word_defect=(A*C*A-A*A)[:5,:5]
word_commutator=A*C-C*A
checks['clock_and_return_act_noncommutatively']=word_commutator!=s.zeros(d)
checks['intermediate_native_clock_is_detected']=word_defect!=s.zeros(5)
actual_moment=(A*A)[:5,:5]
naive_moment=A[:5,:5]**2
checks['source_memory_witness_reused']=s.simplify(actual_moment[0,0]-naive_moment[0,0])==s.Rational(203312,9)
payload=dict(contract='ORIGINAL-RECURSIVE-WORD-SOURCE-20261030',all_passed=all(checks.values()),checks=checks,original_input_dimension=726,reachable_dimension=d,first_source_layer_dimension=first_layer,generator_alphabet=['original single-beta L=Sdagger S','native Cartan J','native Cartan Jdagger'],metric_diagonal=[str(x) for x in G.diagonal()],operators={label:src.mat(M) for label,M in zip(['L','J','Jdagger'],matrices)},original_vector_basis=[{'provenance':p,'polynomial':[[list(mon),str(Q.to_sympy(v))] for mon,v in sorted(f.items())]} for p,f in zip(provenance,vectors)],clock_ordered_word_defect=src.mat(word_defect),direct_coarse_I_plus_L_schur=src.mat(direct),elapsed_seconds=time.time()-start,scope='Exact complete reachable/observable word process for the specified original L,J,Jdagger alphabet in the full neutral grade4 block. General all-word and nested-composition proofs are separate. This is not all E8 current operations, a physical time law, or a full TFPT origin/space-time theorem.')
Path(__file__).with_name('RECURSIVE-WORD-SOURCE.json').write_text(json.dumps(payload,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({k:v for k,v in payload.items() if k not in ['operators','original_vector_basis','metric_diagonal','clock_ordered_word_defect','direct_coarse_I_plus_L_schur']},ensure_ascii=False,indent=2),flush=True)
if not payload['all_passed']:raise RuntimeError('Exact source check failed')
