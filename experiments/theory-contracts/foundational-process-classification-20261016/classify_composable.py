"""Finite classical process classification. No TFPT imports or target fit.

All single-copy Boolean process functions and explicit two-copy counterexamples.
Positive composability follows by triangular recursion, not finite sampling.
"""
from itertools import product,permutations
import json,time
from pathlib import Path
start=time.time();N=3
bits=list(product(range(2),repeat=N));states=tuple(range(8))
# Little-endian world encoding; a local bit function is its two output bits.
contexts=[]
for fs in product(range(4),repeat=N):
 contexts.append(tuple(sum(((fs[k]>>((x>>k)&1))&1)<<k for k in range(N)) for x in states))
valid=[]
for cs in product(range(16),repeat=N):
 w=[]
 for o in states:
  inp=0
  for k,c in enumerate(cs):
   others=[j for j in range(N) if j!=k]
   idx=sum(((o>>j)&1)<<i for i,j in enumerate(others))
   inp|=((c>>idx)&1)<<k
  w.append(inp)
 if all(sum(w[out[x]]==x for x in states)==1 for out in contexts):
  deps=[]
  for k in range(N):
   deps.append({j for j in range(N) if j!=k and any(((w[o]^w[o^(1<<j)])>>k)&1 for o in states)})
  orders=[pi for pi in permutations(range(N)) if all(deps[k]<=set(pi[:i]) for i,k in enumerate(pi))]
  valid.append({'components':cs,'map':w,'orders':orders})
print('single_copy',len(valid),'fixed_order',sum(bool(x['orders']) for x in valid),'no_first',sum(all(c not in (0,15) for c in x['components']) for x in valid),flush=True)
# Each local input/output is two bits, encoded copy0+2copy1.
xs=list(product(range(4),repeat=N));enc=lambda x:sum(v<<(2*k) for k,v in enumerate(x))
x_by_enc=sorted(xs,key=enc)
ops=list(permutations(range(4)))+[tuple([k]*4) for k in range(4)]
# Prefer a small reusable witness bank; exhaustive search only if needed.
bank=[];failed=[];witnesses=[]
for num,v in enumerate(x for x in valid if not x['orders']):
 w=v['map'];ww=[]
 for x in x_by_enc:
  a=w[sum((z&1)<<k for k,z in enumerate(x))];b=w[sum((z>>1)<<k for k,z in enumerate(x))]
  ww.append(sum((((a>>k)&1)+2*((b>>k)&1))<<(2*k) for k in range(N)))
 def fixed(fs):
  return [x for x,xx in enumerate(x_by_enc) if ww[enc(tuple(fs[k][xx[k]] for k in range(N)))]==x]
 found=None
 for fs in bank:
  fixeds=fixed(fs)
  if len(fixeds)!=1:found=(fs,fixeds);break
 if found is None:
  for fs in product(ops,repeat=N):
   fixeds=fixed(fs)
   if len(fixeds)!=1:
    found=(fs,fixeds);bank.append(fs);break
 if found is None:failed.append(v)
 else:witnesses.append({'components':v['components'],'map':w,'local_joint_maps':found[0],'fixed_inputs':found[1]})
assert len(valid)==744 and sum(bool(x['orders']) for x in valid)==488
assert len(witnesses)==256 and not failed
assert all(all(sorted(f)==list(range(4)) for f in v['local_joint_maps']) for v in witnesses)
print('non_fixed_order_with_failure',len(witnesses),'not_separated',len(failed),'witness_bank',len(bank),'seconds',time.time()-start,flush=True)
result={'single_copy_valid':len(valid),'fixed_total_order':sum(bool(x['orders']) for x in valid),'no_constant_component':sum(all(c not in (0,15) for c in x['components']) for x in valid),'non_fixed_total_order':sum(not x['orders'] for x in valid),'non_fixed_order_two_copy_failures':len(witnesses),'not_separated':len(failed),'witness_operation_set':'24 reversible maps plus 4 constant maps on two local bits','valid_processes':valid,'two_copy_failure_witnesses':witnesses,'witness_bank':bank,'scope':'Binary deterministic three-party one-shot functions only. Fixed-order positives follow for all joint operations/copies by triangular recursion; negatives have explicit two-copy witnesses. No general-N or physical quantum/TFPT theorem.'}
result['all_recorded_two_copy_witnesses_are_reversible']=True
result['fixed_point_count_histogram']={str(k):sum(len(x['fixed_inputs'])==k for x in witnesses) for k in (0,1,2)}
open(Path(__file__).with_name('composition_classification.json'),'w').write(json.dumps(result,indent=2)+'\n')
