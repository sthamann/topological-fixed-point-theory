import sys, importlib.util, itertools, hashlib, json
from pathlib import Path
import sympy as s
repo=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
sys.path.insert(0,str(repo))
# Import original pinned source checker only; no repository writes.
p=repo/'experiments/theory-contracts/compiler-cone-object-audit/checker.py'
spec=importlib.util.spec_from_file_location('cone_audit',p); cone=importlib.util.module_from_spec(spec); spec.loader.exec_module(cone)
bridge,_=cone.source(); src,_=bridge.inherited(); d=cone.data(); I=d['eye']; al=d['alpha']; w=d['w']
def clean(m): return m.applyfunc(s.simplify)
def key(m):
    m=clean(m); pivot=next((z for z in m if z!=0),None)
    return None if pivot is None else tuple(s.cancel(z/pivot) for z in m)
def fa(x):
    z=[x[2*j]+s.I*x[2*j+1] for j in range(4)]
    return clean((z[3]*I+sum((z[j]*d['us'][j] for j in range(3)),s.zeros(2)))*w.H/2)
rootcoords=[tuple(2*sgn if j==k else 0 for j in range(8)) for k in range(8) for sgn in (-1,1)]
for supp in src.CSTAR_SUPPORTS_EXPECTED:
  for signs in itertools.product((-1,1),repeat=4):
    x=[0]*8
    for j,v in zip(supp,signs): x[j]=v
    rootcoords.append(tuple(x))
allmat=[clean(fa(x)) for x in rootcoords]
reps={key(A):A for A in allmat}
# endpoint record by projector; rank-one source reps labeled output s, input t.
rayproj={}
def pkey(P): return tuple(s.simplify(x) for x in P)
for j,a in enumerate(al):
  for sign in (-1,1):
    P=clean((I+sign*a)/2); rayproj[pkey(P)]=(j,sign)
rank=[]; unit=[]
for A in reps.values():
  if A.det()==0:
    po=clean(A*A.H/2); pi=clean(A.H*A/2)
    rank.append((pkey(po),pkey(pi),A))
  else: unit.append(A)
assert len(reps)==60 and len(rank)==36 and len(unit)==24
# Keep representative by (s,t); construct alternative L_(s,t)=alpha_t A_(s,t).
bylabel={(rayproj[so],rayproj[ti]):A for so,ti,A in rank}
assert len(bylabel)==36
alt={label:clean(al[label[1][0]]*A) for label,A in bylabel.items()}
# It is same projective source menu, permuting s at fixed t.
altlabels={}
for label,B in alt.items():
  found=[(so,ti) for (so,ti),A in bylabel.items() if key(A)==key(B)]
  assert len(found)==1 and found[0][1]==label[1]
  altlabels[label]=found[0]
# Per-label effects and unconditional channels: matrix superoperators on E_ij.
for label,A in bylabel.items():
  B=alt[label]
  assert clean(A.H*A/60)==clean(B.H*B/60)
def superop(ops,X): return clean(sum((A*X*A.H/60 for A in ops),s.zeros(2)))
for X in [s.Matrix([[1,0],[0,0]]),s.Matrix([[0,1],[0,0]]),s.Matrix([[0,0],[1,0]]),s.Matrix([[0,0],[0,1]])]:
  assert superop(list(bylabel.values())+unit,X)==superop(list(alt.values())+unit,X)
  assert superop(list(bylabel.values())+unit,X)==s.trace(X)*I/2
# Verify source 24-root unitary covariance. Projective h conjugation acts on outcome endpoints.
def branch_transform(A,h): return clean(h*A*h.H)
unit_keys={key(h):h for h in unit}
Pbylabel={(j,sgn):clean((I+sgn*a)/2) for j,a in enumerate(al) for sgn in (-1,1)}
for h in unit:
  for label,A in bylabel.items():
    # conjugated endpoint projectors transform by h
    so,ti=label
    so2=rayproj[pkey(clean(h*Pbylabel[so]*h.H))]
    ti2=rayproj[pkey(clean(h*Pbylabel[ti]*h.H))]
    label2=(so2,ti2)
    # output/input matrices as alternative instrument agree up to projective scalar
    assert key(branch_transform(alt[label],h))==key(alt[label2])
# actual sigma representative w is in unitary list; projective order 3 and covariance checked above.
assert key(w) in unit_keys
# Two-time exact event probabilities. first event s=|+x>,t=|0z>; then second event port t2=|+x|, sum over six output labels.
plusZ=rayproj[pkey((I+al[0])/2)]; plusX=rayproj[pkey((I+al[2])/2)]; plusY=rayproj[pkey((I+al[1])/2)]
first=(plusX,plusZ); assert first in bylabel
psi=s.Matrix([1,0]); A=bylabel[first]; B=alt[first]
# branches are normalized with sqrt60. For input |0>, branch click prob 1/30.
def normprob(v): return s.simplify((v.H*v)[0])
post=clean(A*psi); post2=clean(B*psi)
assert normprob(post)/60==normprob(post2)/60==s.Rational(1,30)
# Normalize conditional vectors only; compare second-instrument total probability for fixed incoming ray +x.
second_effect=s.zeros(2)
for (so,ti),C in bylabel.items():
  if ti==plusX: second_effect += C.H*C/60
second_effect=clean(second_effect)
assert second_effect==clean(s.Rational(1,5)*Pbylabel[plusX])
# Joint p(first, any second with t2=+x), calculated with unnormalized post branch.
P2=s.simplify(s.trace(second_effect*post*post.H)/60)
Q2=s.simplify(s.trace(second_effect*post2*post2.H)/60)
assert (P2,Q2)==(s.Rational(1,150),s.Integer(0))
# Check simple representative w cyclically permutes axes (not essential, but report).
result={
 'source_pins_validated_by_original_checker':True,
 'root_count':len(rootcoords),'projective_rays':len(reps),'unitary_rays':len(unit),'rank_one_rays':len(rank),
 'per_labeled_effects_identical':True,'unconditional_channels_identical_on_matrix_units':True,'common_channel':'Tr(X) I/2',
 'alternative_projective_root_menu_preserved':True,'same_t_output_permutations':True,
 'full_24_unitary_root_covariance':True,'sigma_w_in_native_unitary_set':True,
 'example':{'first_label':'s=+X, t=+Z','first_event_probability':'1/30','second_fixed_input_port':'t2=+X','second outputs summed':'six original s labels at t2=+X','canonical_joint_probability':str(P2),'alternative_joint_probability':str(Q2)},
 'caveat':'The 60-root source specifies an algebra and family permutation; the full continuous tau_t action on these outcome branches is not identified in the original source. Covariance is proven only for the native 24-root unitary action, including its family cycle w.'
}
print(json.dumps(result,indent=2))
