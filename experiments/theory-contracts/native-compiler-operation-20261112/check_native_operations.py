from pathlib import Path
import contextlib,io,json
import sympy as s
base=Path(__file__).parent
prior=base.parent/'compiler-attachment-audit-20261112/audit.py'
ns={'__file__':str(prior)}
with contextlib.redirect_stdout(io.StringIO()):exec(prior.read_text().split("(base/'results.json').write_text",1)[0],ns)
L,G,P,synt,N=[ns[x] for x in ['L','G','P','synt','N']]
K=ns['K'];ident=s.kronecker_product(P,P)
# K acts on V tensor V in ambient R6 tensor R6.
raw=ns['ns'];RR=[(L*r*G.inv()*L.T).applyfunc(s.simplify) for r in raw['R']]
PM=[];perms=[]
for matching in synt:
 perm=list(range(6))
 for a,b in matching:perm[a]=b;perm[b]=a
 Q=s.Matrix(6,6,lambda i,j:int(i==perm[j]));perms.append(Q)
 PM.append(P*(s.eye(6)+Q)*P/2)
matching=[next(j for j,y in enumerate(PM) if r==y) for r in RR]
assert sorted(matching)==list(range(15))
assert sum(PM,s.zeros(6))==6*P
F=ns['Fs'];pi=[x/2 for x in F]
Born=s.Matrix(15,15,lambda a,m:s.trace(pi[a]*PM[m]))
assert Born==(s.ones(15)-N)/2
assert all(sum(Born[:,j])==6 for j in range(15)) # sum pure root projectors = 3P; each R has trace2
assert all(sum(Born[i,:])==6 for i in range(15))
Qpair=sum((s.kronecker_product(r,r) for r in PM),s.zeros(36))
assert Qpair==s.Rational(3,2)*ident+s.Rational(9,2)*K
Plus=Qpair/6
Minus=sum((s.kronecker_product(P-r,P-r) for r in PM),s.zeros(36))/9
Ref=sum((s.kronecker_product(P-2*r,P-2*r) for r in PM),s.zeros(36))/15
assert Plus==ident/4+3*K/4
assert Minus==ident/2+K/2
assert Ref==-ident/5+6*K/5
for r,t in zip(PM,perms):assert P-2*r==-P*t*P
J=s.Rational(2,5)*Plus;Idle=s.Rational(3,5)*Minus
# Evaluate on all actual spectral sectors; ambient complement ignored.
k=s.symbols('k');stop=(1+3*k)/(7-3*k)
assert [s.factor(stop.subs(k,v)) for v in [1,s.Rational(4,9),0]]==[1,s.Rational(7,17),s.Rational(1,7)]
# New pasted text's partial isometry from Duads to internal Petersen edges.
A=ns['A'];M=ns['M'];edges=ns['edges']
J15=s.ones(15);E1=J15/15;E9=(A+3*s.eye(15)-9*E1)/4;Q=E1+E9
indices=[i for i,e in enumerate(edges) if 5 not in e]
Restr=s.Matrix(10,15,lambda i,j:int(j==indices[i]))
B=Restr*Q;H=B*B.T
AP=s.Matrix(10,10,lambda i,j:A[indices[i],indices[j]])
P1=s.ones(10)/10;P5=(AP+2*s.eye(10)-5*P1)/3;P4=s.eye(10)-P1-P5
# Polar normalizer exactly determined on Petersen irreducibles.
hs=[s.trace(H*pr)/s.trace(pr) for pr in [P1,P4,P5]]
Hinvroot=sum((pr/s.sqrt(h) for pr,h in zip([P1,P4,P5],hs)),s.zeros(10))
FF=(Hinvroot*B).applyfunc(s.simplify)
assert (FF*FF.T-s.eye(10)).applyfunc(s.simplify)==s.zeros(10)
assert (FF.T*FF-Q).applyfunc(s.simplify)==s.zeros(15)
assert (FF*M*FF.T-(AP/3)**2-P5/3).applyfunc(s.simplify)==s.zeros(10)
assert s.trace(P1+P4)==5 # pasted 'six' is an arithmetic error
res={'original_checks':len(raw['checks']),'rank2_plane_to_syntheme':matching,'born_overlap':'(1-N)/2','native_POVM_transition':'(1-N)/12','plus_channel':'id/4+3K/4','minus_channel':'id/2+K/2','reflection_channel':'-id/5+6K/5','first_plus_stopped_spectrum':{'1':1,'7/17':9,'1/7':15},'partial_isometry_gram_eigenvalues':list(map(str,hs)),'partial_isometry_exact':True,'petersen_matching_dimensions':5,'petersen_difference_rank':5,'physical_source_derived':False,'F1_complete':False}
(base/'results.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res,indent=2))
