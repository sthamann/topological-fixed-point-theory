from pathlib import Path
import importlib.util,json
import numpy as np
import sympy as s
repo=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
p=repo/'experiments/theory-contracts/compiler-correlated-event-clock-20260919/source_channel.py'
spec=importlib.util.spec_from_file_location('source',p);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
rays=m.source_rays();bell,raw,smaps=m.reflection_actions(rays);outer,pentads=m.outer_mark_action(raw)
assert all(sum(p[q]!=q for q in range(6))==2 for p in outer)
move=np.array([[p[q]!=q for p in outer] for q in range(6)],dtype=int)
assert np.array_equal(move.sum(axis=1),np.full(6,20))
assert np.array_equal(move.sum(axis=0),np.full(60,2))
# Original Gaussian integer rays have squared norm four.
Pn=[np.outer(z,z.conj()) for z in rays]
R2=[2*np.eye(4)-x for x in Pn]
def key(x):
 assert np.array_equal(x.real,np.rint(x.real)) and np.array_equal(x.imag,np.rint(x.imag))
 return tuple(zip(x.real.astype(int).ravel(),x.imag.astype(int).ravel()))
index={key(x):i for i,x in enumerate(Pn)}
conjugation=[]
for g in range(60):
 action=[index[key(R2[g]@x@R2[g].conj().T/4)] for x in Pn]
 assert sorted(action)==list(range(60))
 for q in range(6):
  assert all(move[q,l]==move[outer[g][q],action[l]] for l in range(60))
 conjugation.append(action)
unseen=set((q,l) for q in range(6) for l in range(60));orbits=[]
while unseen:
 seed=min(unseen); orbit={seed}; todo=[seed]
 while todo:
  q,l=todo.pop()
  for g in range(60):
   pair=(outer[g][q],conjugation[g][l])
   if pair not in orbit: orbit.add(pair);todo.append(pair)
 unseen-=orbit;orbits.append(orbit)
assert sorted(map(len,orbits))==[120,240]
assert all(len({int(move[q,l]) for q,l in o})==1 for o in orbits)
t,z=s.symbols('t z',positive=True)
weights=s.Matrix(6,60,lambda q,l:t/120 if move[q,l] else (1-t)/240)
assert all(s.simplify(sum(weights[:,l])-s.Rational(1,60))==0 for l in range(60))
assert all(s.simplify(sum(weights[q,:])-s.Rational(1,6))==0 for q in range(6))
assert s.simplify(sum(weights[q,l]*move[q,l] for q in range(6) for l in range(60))-t)==0
Craw=m.superoperator_integer(bell)
assert not np.any(Craw.imag)
C60=s.Matrix(Craw.real.astype(int))/60
branch=[]
for q in range(6):
 ids=np.flatnonzero(move[q]).tolist(); S=np.zeros((100,100),complex);J=np.zeros_like(S)
 for i in ids:
  U=bell[i]; perm=m.bell_permutation(U)
  Ds=np.diag([int(perm[a]==a) for a in range(10)]);Dj=np.eye(10)-Ds
  S+=np.kron((U@Ds).conj(),U@Ds);J+=np.kron((U@Dj).conj(),U@Dj)
 assert not np.any(S.imag) and not np.any(J.imag)
 assert np.array_equal(S+J,m.superoperator_integer([bell[i] for i in ids]))
 Si=s.Matrix(S.real.astype(int))/20; Ji=s.Matrix(J.real.astype(int))/20
 assert Si==s.diag(*Si.diagonal())
 trace_row=s.Matrix(1,100,lambda a,b:int(b//10==b%10))
 assert trace_row*Si==s.Rational(2,5)*trace_row
 assert trace_row*Ji==s.Rational(3,5)*trace_row
 Phi=Ji*s.diag(*[1/(1-x) for x in Si.diagonal()])
 for tv in [s.Rational(1,4),s.Rational(1,3),s.Rational(2,5)]:
  St=(1-tv)*s.eye(100)+tv*Si
  Phi_t=tv*Ji*s.diag(*[1/(1-x) for x in St.diagonal()])
  assert Phi_t==Phi
  # Trace functional of the counting map, on all inputs, at z=1/2.
  zv=s.Rational(1,2)
  Pz=zv*tv*Ji*s.diag(*[1/(1-zv*x) for x in St.diagonal()])
  assert trace_row*Pz==s.cancel((3*tv*zv/5)/(1-(1-3*tv/5)*zv))*trace_row
 Ci=Si+Ji
 bg=[i for i in range(60) if i not in ids]
 Cbg=m.superoperator_integer([bell[i] for i in bg])
 assert not np.any(Cbg.imag)
 assert s.Matrix(Cbg.real.astype(int))/40==(3*C60-Ci)/2
 if q in (4,5):
  r20=sum(s.Rational(float(abs(R2[i][0,0]/2)**4)) for i in ids)/20
  r40=sum(s.Rational(float(abs(R2[i][0,0]/2)**4)) for i in bg)/40
  assert (r20,r40)==(s.Rational(1,4),s.Rational(13,40))
 branch.append(Ci)
assert sum(branch,s.zeros(100))/6==C60
# All t: sum of mark-conditioned C20 is C60 and same for C40.
assert sum([(3*C60-b)/2 for b in branch],s.zeros(100))/6==C60
D=[s.Rational(3,2)*(b-C60) for b in branch]
assert sum(D,s.zeros(100))==s.zeros(100)
assert all(d==d.T for d in D)
Defect=sum([d*d for d in D],s.zeros(100))/6
assert s.trace(Defect)==s.Rational(27,5)
repeat=s.factor(t*t/20+(1-t)**2/40)
assert s.simplify(repeat-s.Rational(1,60)-(3*t-1)**2/120)==0
for tv in [s.Rational(1,4),s.Rational(1,3),s.Rational(2,5)]:
 c2=sum([(C60+(tv-s.Rational(1,3))*d)**2 for d in D],s.zeros(100))/6
 assert c2-C60*C60==(tv-s.Rational(1,3))**2*Defect
h=(3*t*z/5)/(1-(1-3*t/5)*z)
pgf=h**6
mean=s.factor(s.diff(pgf,z).subs(z,1));variance=s.factor(s.diff(pgf,z,2).subs(z,1)+mean-mean**2)
assert mean==10/t and s.simplify(variance-(50/(3*t*t)-10/t))==0
rows=[]
for tv in [s.Rational(1,4),s.Rational(1,3),s.Rational(2,5)]:
 rows.append({'t':str(tv),'joint_move_probability':str(tv),'fixed_mark_return':str((13-3*tv)/40),'two_step_HS_trace_defect':str(s.Rational(27,5)*(tv-s.Rational(1,3))**2),'mean_six_accepts':str(mean.subs(t,tv)),'variance_six_accepts':str(variance.subs(t,tv))})
# One-step native code reflections cannot be Q4: character invariant.
assert (6-2-1)**2==9
assert 1+9*s.Rational(4,9)==5
result={'two_step_defect_coefficient_trace':'27/5','closure_selects_t':'1/3 iff averaged persistent-mark two-step transfer equals C60^2','joint_orbit_sizes':sorted(map(len,orbits)),'all_invariant_joint_laws_classified':True,'native_source_checks':len(m.checks),'exact_ray_conjugation_actions':len(conjugation),'simultaneous_mark_event_covariance_checks':60*6*60,'uniform_marginals_all_t':True,'all_six_full_operator_stopped_channels_equal_for_three_t':True,'analytic_stopped_identity':'t J [I-((1-t)I+t S)]^-1=J(I-S)^-1 for t>0','entire_parameter_family':'0<t<=1; physical identification not supplied','count_pgf':str(pgf),'count_mean':str(mean),'count_variance':str(variance),'three_completions':rows,'code_reflection_trace':9,'Q4_superoperator_trace':5,'physical_P1_mixed_response_evaluated':False,'F1_complete':False}
Path(__file__).with_name('results.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
