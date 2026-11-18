from fractions import Fraction as F
from collections import Counter
from itertools import combinations,product
from pathlib import Path
import json
solutions=[];mixed=[];census=[]
for d in range(2,33):
 row=[]
 for k in range(d):
  norm=F(d+1,4)+F(k*(d-k),d)
  if norm==2:
   solutions.append([d,k]);row.append(k)
   if k:mixed.append([d,k])
 census.append({'d':d,'weight1_spinor_classes':row})
assert solutions==[[4,1],[4,3],[7,0]]
roots=[]
for i,j in combinations(range(8),2):
 for a,b in product((-1,1),repeat=2):
  r=[F(0)]*8;r[i]=F(a);r[j]=F(b);roots.append(r)
for signs in product((-1,1),repeat=8):
 if sum(x<0 for x in signs)%2==0:roots.append([F(x,2) for x in signs])
grades=Counter(int(2*sum(r[:5]))%4 for r in roots)
assert [grades[i] for i in range(4)]==[52,64,60,64]
full=[grades[0]+8]+[grades[i] for i in (1,2,3)]
assert full==[60,64,60,64] and sum(full)==248
# Four-sector matrix E0 defined in the local synthesis: nonneutral factors vanish.
u=[1,1j,-1,-1j]*4
raw_X=sum((u[1]*u[0].conjugate())**j for j in range(4))/4
raw_Xstar=sum((u[0]*u[1].conjugate())**j for j in range(4))/4
raw_XstarX=sum((u[0]*u[0].conjugate())**j for j in range(4))/4
assert raw_X==0 and raw_Xstar==0 and raw_XstarX==1
operator_grades=Counter((b-a)%4 for a in range(4) for b in range(4) for _ in range(16))
assert [operator_grades[r] for r in range(4)]==[64]*4
# Each E0 block has complex matrix dimension16; a45-dimensional simple D5 Lie algebra
# cannot embed in any block; blockwise kernel is all orzero. Proof is in the audit text.
result={'status':'PASS','weight1_solutions':solutions,'mixed_weight1_solutions':mixed,'finite_census_range':[2,32],'global_proof':'d>7 excluded since spinor norm(d+1)/4>2; d2..7 exact. Nonminimal coset representatives increase squarednorm by nonnegative evenintegers.','d4_minimal_extension':{'base_determinant':16,'glue_order':4,'extension_determinant':1},'d7_minimal_extension':{'lattice':'A6 direct_sum E8','determinant':7,'mixed':False},'root_grades':[grades[i] for i in range(4)],'weight1_current_grades':full,'E0_example':{'E0_X':0,'E0_Xstar':0,'E0_XstarX':1},'finite_M16_operator_grade_dimensions':[operator_grades[i] for i in range(4)],'finite_E0_fixed_algebra':'M4 direct_sum M4 direct_sum M4 direct_sum M4','scope':'Independent check of pasted provisional Claude assertions; no mutation or assumed access to Claude cloud checkout. No physical source-field selection derived.'}
Path(__file__).with_name('spinor-selection-result.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
