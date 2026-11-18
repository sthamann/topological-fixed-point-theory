"""Actual v252 real-bimodule Y insertion, full Dirac contacts and seam target type.
No masses, physical coupling values, measure or scalar response are fitted.
"""
import ast,itertools,json,hashlib,sys
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import numpy as np
import sympy as sp
REPO=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4');OUT=Path('work/y-caller');OUT.mkdir(exist_ok=True,parents=True)
path=REPO/'verification/v252_full_finite_triple.py';t=ast.parse(path.read_text())
# Original constructor definitions and their literal basis data only, no run().
names={'NG','GEN_FIELDS','P','IDX','NP','NH'};funcs={'build_index','dirac_block','full_D','grading','Kswap','rep'};nodes=[]
for n in t.body:
    if isinstance(n,ast.FunctionDef) and n.name in funcs:nodes.append(n)
    elif isinstance(n,ast.Assign):
        targets={z.id for target in n.targets for z in ast.walk(target) if isinstance(z,ast.Name)}
        if targets & names:nodes.append(n)
env={'np':np};exec(compile(ast.Module(body=nodes,type_ignores=[]),str(path),'exec'),env)
P,NP,NH,rep,K=(env[k] for k in ['P','NP','NH','rep','Kswap']);K=K()
checks=Counter()
def ck(v,label):
    if not bool(v):raise RuntimeError(label)
    checks[label]+=1
# Infinitesimal unimodular real-bimodule lift: lambda=i/2, m=-i I3/6.
# Scale by6 before representation, so every source matrix operation is integral.
X=rep(3j,np.zeros((2,2),complex),-1j*np.eye(3))
Y6=(X+K@X.conj()@K)/1j
X3=rep(0j,np.diag([3j,-3j]),np.zeros((3,3),complex))
T36=(X3+K@X3.conj()@K)/1j;Q6=Y6+T36
ck(np.array_equal(Y6,np.diag(np.diag(Y6))),'actual real-bimodule hypercharge diagonal')
yp=[]
for name,ch,c,g in P:
    yp.append(1 if ch=='L' and name in ['u','d'] else -3 if ch=='L' else {'u':4,'d':-2,'e':-6,'nu':0}[name])
ck(np.array_equal(np.diag(Y6),yp+[-v for v in yp]),'real lift produces original physical L/R charges and their conjugates')
gam=env['grading']();ix=np.flatnonzero(np.diag(gam)==1)
# Native even exterior16 is the all-left-handed chart, including charge conjugates of R.
native=[]
for bits in itertools.product([0,1],repeat=5):
    if sum(bits)%2:continue
    y=-2*sum(bits[:3])+3*sum(bits[3:]);t3=3*(bits[3]-bits[4]);native.append((y,t3))
actual=Counter((int(Y6[i,i].real),int(T36[i,i].real)) for i in ix)
ck(actual==Counter({k:3*v for k,v in Counter(native).items()}),'native all-left48 maps to gamma-plus full96, not particle48')
ck(Counter(int(x) for x in yp)!=Counter(y for y,t in native for _ in range(3)),'particle48 is a different chiral charge chart')

# Exhaust all real and imaginary Yukawa basis directions and symmetric Majorana directions.
Z=np.zeros((3,3),complex);basis=[]
for f in range(4):
    for i,j in itertools.product(range(3),repeat=2):
        for phase in [1,1j]:
            args=[Z.copy() for _ in range(5)];args[f][i,j]=phase;basis.append(('Yukawa',args))
for i in range(3):
    for j in range(i,3):
        for phase in [1,1j]:
            args=[Z.copy() for _ in range(5)];args[4][i,j]=phase;args[4][j,i]=phase;basis.append(('Majorana',args))
def comm(a,b):return a@b-b@a
nonzero=Counter()
for kind,args in basis:
    D=env['full_D'](*args);CY=comm(Y6,D);CQ=comm(Q6,D)
    ck(np.count_nonzero(CQ)==0,'physical EM direction commutes with whole neutral-H Dirac basis')
    hphi=-1j*gam@D
    ck(np.array_equal(comm(Y6,hphi),-1j*gam@CY),'original collar h retains actual hypercharge Higgs commutator')
    if kind=='Yukawa':
        ck(np.count_nonzero(CY)>0,'hypercharge variation has actual Higgs Yukawa contact')
        ck(np.array_equal(comm(Y6,CY),9*D),'all Yukawa double hypercharge commutators equal9 D')
    else:ck(np.count_nonzero(CY)==0,'whole Majorana block neutral under Y and EM')
    nonzero[kind]+=1
moments={str(n):{'native_chiral48':str(sum(F(int(Y6[i,i].real),6)**n for i in ix)),'full96':str(sum(F(int(v.real),6)**n for v in np.diag(Y6)))} for n in [0,2,4,5]}

# Original minimal seam response: delta is an INPUT amplitude, not an output of Y norm.
a,d=sp.symbols('a d',positive=True);q=d*sp.exp(-2*a);W=-sp.Rational(5,4)*sp.log(1-q)
fourth=sp.factor(sp.diff(W,a,4));expected=20*q*(1+4*q+q*q)/(1-q)**4
ck(sp.simplify(fourth-expected)==0,'actual supplied seam functional fourth derivative')
result={'claim':'SOURCE.Y.DIRAC.CALLER.CHECK.20261106','original_v252_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'real_bimodule_Y6_particle_diagonal':yp,'native48_target':'gamma_F=+1 subspace of full96, physical L particles plus conjugates of physical R; opposite gamma block via J','native_joint_Y_T3_weight_match':True,'full_D_basis_directions':dict(nonzero),'Dirac_contacts':{'Qem':'[Qem,D_F]=0 for every original Yukawa and symmetric MR coefficient on the neutral-H chart','hypercharge':'[Y,[Y,D_F]]=D_Yukawa/4, [Y,D_Majorana]=0; all72 real Yukawa and12 real Majorana basis directions checked','scope':'Original finite internal coefficient, not full spacetime/gauge/geometry integration; Higgs and connection variations required for local gauge Ward identity'},'charge_moments':moments,'minimal_seam_response':{'definition':'W_Sigma(a)=-(5/4) log(1-delta exp(-2a))','fourth_derivative':'20 q(1+4q+q^2)/(1-q)^4, q=delta exp(-2a)','delta_input':'delta=48 c3^4 is supplied before differentiation; replacing it by a variable d leaves the algebra valid and does not select its value','target_type':'Zero-source seam defect amplitude; not the coefficient of the fourth Y source variation'},'joint_physical_response':None,'delta_48_derived':False,'fundamental_solution':False,'checks':dict(checks)}
result['original_neutral_H_chart']='R078 q(H)=[[bar h2,h1],[-bar h1,h2]]; original v252 diagonal mass edges correspond to q(H)=I, H=(0,1).'
result['first_time_intertwining_test']='Native [H8,JY(0)]=0, but original fixed-H collar [Y,h]=-i gammaF[Y,Phi] is nonzero for Yukawa directions. A simultaneous rigid H8->frozen h, JY0->fermionic Y96 map fails this test. Dynamic original bosonic/Higgs charge and HB are retained; no full-source no-go.'
result['finite_affine_mode_trace_obstruction']='Full affine Y modes cannot all be represented on96: [JY(1),JY(-1)]=(5/6)I would have trace80 rather than0. Only the zero-mode internal action is being represented; original unbounded spacetime factors remain.'
(OUT/'dirac-y-contact-results.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
