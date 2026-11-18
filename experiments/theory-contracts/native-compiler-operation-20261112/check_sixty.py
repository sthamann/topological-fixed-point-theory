from pathlib import Path
import importlib.util,itertools as it,json
import numpy as np
import sympy as s
p=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4/experiments/theory-contracts/compiler-correlated-event-clock-20260919/source_channel.py')
spec=importlib.util.spec_from_file_location('native_source',p);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
rays=m.source_rays();bell,six,smaps=m.reflection_actions(rays);outer,pentads=m.outer_mark_action(six)
# exact Gaussian integer arithmetic happens in source ray/dictionary operations
Ps=[np.outer(z,z.conj())/4 for z in rays]
classes={t:[i for i,p in enumerate(six) if p==t] for t in sorted(set(six))}
assert all(np.array_equal(sum((Ps[i] for i in ids),np.zeros((4,4))),np.eye(4)) for ids in classes.values())
assert all(all(np.vdot(rays[i],rays[j])==0 for i,j in it.combinations(ids,2)) for ids in classes.values())
rows=[];returns=[];phase_hist=[]
S60=s.Matrix(m.superoperator_integer(bell).real.astype(int));C60=S60/60
# Checks of projectors prove complete spectral decomposition, not sampled eigenvalues.
def projections(C,vals):
 out=[]
 for a in vals:
  Q=s.eye(C.rows)
  for b in vals:
   if a!=b:Q=Q*(C-b*s.eye(C.rows))/(a-b)
  out.append(Q)
 assert sum(out,s.zeros(C.rows))==s.eye(C.rows)
 assert all(C*Q==v*Q for v,Q in zip(vals,out))
 return out
vals60=[s.Integer(1),s.Rational(1,3),s.Rational(1,5),s.Rational(1,15)]
P60=projections(C60,vals60)
# All six marked subsets: common one-copy map and concrete two-copy return.
for q in range(6):
 idx=[i for i,p in enumerate(outer) if p[q]!=q];bg=[i for i in range(60) if i not in idx]
 groups=[]
 for selected in [idx,bg,list(range(60))]:
  R=[np.eye(4)-2*Ps[i] for i in selected]
  one=sum((np.kron(r.conj(),r) for r in R),np.zeros((16,16)))/len(R)
  target=np.eye(16)/5+np.outer(np.eye(4).reshape(16),np.eye(4).reshape(16))/5
  assert np.max(np.abs(one-target))<1e-14
  # Repeat exactly on units multiplied by5*len; binary fractions are exact.
  numer=sum((np.kron((2*r).conj(),2*r) for r in R),np.zeros((16,16)))
  assert np.array_equal(5*numer,4*len(R)*np.eye(16)+4*len(R)*np.outer(np.eye(4).reshape(16),np.eye(4).reshape(16)))
  groups.append(sum(abs((np.eye(4)-2*Ps[i])[0,0])**4 for i in selected)/len(selected))
 returns.append(groups)
 # For selected set group into five orthogonal classes and count oriented lex triples.
 selected_classes=[ids for ids in classes.values() if set(ids).issubset(idx)]
 assert len(selected_classes)==5
 assert all(abs(np.vdot(rays[a],rays[b]))**2==4 for c,d in it.combinations(selected_classes,2) for a in c for b in d)
 hist={}
 for c,d,e in it.combinations(selected_classes,3):
  for a,b,k in it.product(c,d,e):
   z=np.vdot(rays[a],rays[b])*np.vdot(rays[b],rays[k])*np.vdot(rays[k],rays[a])
   key=(int(z.real),int(z.imag));hist[key]=hist.get(key,0)+1
 assert sorted(hist.values())==[160,160,320]
 phase_hist.append({str(k):v for k,v in hist.items()})
 if q==0:
  C20=s.Matrix(m.superoperator_integer([bell[i] for i in idx]).real.astype(int))/20
  C40=(3*C60-C20)/2
  assert C20*C60==C60*C20
  vals20=[s.Integer(1),s.Rational(3,5),s.Integer(0),s.Rational(1,5),-s.Rational(1,5)]
  P20=projections(C20,vals20)
  for a,Q in zip(vals60,P60):
   for b,R in zip(vals20,P20):
    dim=s.trace(Q*R)
    if dim:rows.append([int(dim),str(a),str(b),str((3*a-b)/2)])
expected=sorted([[1,'1','1','1'],[5,'1/3','3/5','1/5'],[4,'1/3','0','1/2'],[45,'1/5','1/5','1/5'],[30,'1/15','1/5','0'],[15,'1/15','-1/5','1/5']])
assert sorted(rows)==expected
t=s.symbols('t', real=True)
family=[[dim,str(s.factor(t*s.sympify(b)+(1-t)*s.sympify(c)))] for dim,a,b,c in rows]
family_returns={str(t0):str(s.Rational(13,40)-3*t0/40) for t0 in [s.Rational(1,4),s.Rational(1,3),s.Rational(2,5)]}
assert list(family_returns.values())==['49/160','3/10','59/200']
res={'positive_HS_mixture_interval':'0 <= t <= 1/2 (t is selected20 prior weight)', 'mixture_spectrum':family, 'three_positive_full_support_returns':family_returns, 'native_source_checks':len(m.checks),'orthogonal_contexts':15,'complete_MUB_subsets':6,'joint_spectrum':rows,'returns_20_40_60':returns,'bargmann_integer_histograms':phase_hist,'phase_convention':'one orientation per triple of ordered basis classes, lexicographic order; conjugate on reversal','origin_selected':False}
Path(__file__).with_name('sixty-results.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res,indent=2))
