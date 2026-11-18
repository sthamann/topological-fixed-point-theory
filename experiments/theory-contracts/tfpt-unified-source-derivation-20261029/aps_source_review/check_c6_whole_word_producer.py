"""Whole original C6 component; no replacement physical stochastic source."""
import importlib.util
import itertools
import json
from pathlib import Path
import sys
import mpmath as mp
import sympy as s

repo=Path(__file__).resolve().parents[4]
sys.path.insert(0,str(repo/'verification'))
spec=importlib.util.spec_from_file_location('native_c6',repo/'verification/v17_hexagonal_resolvent.py')
native=importlib.util.module_from_spec(spec);spec.loader.exec_module(native)
U=native.U6_power(1);I=s.eye(6);C=(U+U.T)/2
checks=[]
def need(ok,name):
    if not ok:raise AssertionError(name)
    checks.append(name)
need(U**6==I,'unchanged original U6')
need(4*C**4-5*C**2+I==s.zeros(6),'original cosine polynomial')
half=s.Rational(1,2)
for n in range(11):
    need(sum(s.rf(half,k)*s.rf(half,n-k)/(s.factorial(k)*s.factorial(n-k)) for k in range(n+1))==1,f'all splits of fixed length{n}')
for n in range(1,9):
    counts=[0]*6
    for word in itertools.product([-1,1],repeat=n):counts[sum(word)%6]+=1
    need(list((C**n)[0,:])==[s.Rational(c,2**n) for c in counts],f'full native endpoint sum length{n}')
    need((C**n).trace()==s.Rational(6*counts[0],2**n),f'all based closed words length{n}')
def edges(word):
    v=0;path=[]
    for d in word:
        w=(v+d)%6;path.append((v,w));v=w
    return path
need(sum((a,b)==(1,2) for a,b in edges([1,1,-1,-1]))==1,'ordered word visits marked edge')
need(sum((a,b)==(1,2) for a,b in edges([1,-1,1,-1]))==0,'same endpoint different ordered word')
y,d,e,nu=s.symbols('y d e nu',positive=True)
a=y*y+d*d+e*e;k=(y-d)**2+e*e;q=2*y*d/a
mean=nu*q/(1-q)
score_abs=-nu*s.diff(a,d)/a+mean*(1/d-s.diff(a,d)/a)
need(s.simplify(score_abs+nu*s.diff(k,d)/k)==0,'absolute row derivative retains amplitude')
norm_score=s.diff(q,d)*(mean/q-nu/(1-q))
need(s.simplify(norm_score)==0,'normalized source score mean0')
need(s.simplify(s.diff(q,d)*(s.Symbol('n')/q-nu/(1-q))-nu*s.diff(k,d)/k-(-2*nu*d/a+s.Symbol('n')*(1/d-2*d/a)))==0,'wordwise normalized and absolute derivative relation')
# Real original D tangent and a nonuniform edge control, without physical J claim.
D=y*I-d*U;Q=D.T*D+e*e*I
V=s.zeros(6);V[0,1]=-1
L=D.T*V+V.T*D
expected=s.zeros(6);expected[0,1]=-y;expected[1,0]=-y;expected[1,1]=2*d
need(L==expected,'native edge Gram first derivative ordered correctly')
z=s.symbols('z')
need(s.diff(-s.log(1+z),z).subs(z,0)==-1,'negative log determinant series sign')

mp.mp.dps=55
delta=((mp.mpf(794)-7*mp.sqrt(9961))/2187)**(mp.mpf(1)/6)
cos=[mp.cos(mp.pi*j/3) for j in range(6)]
records=[]
for cusp in [mp.mpf(1),mp.mpf(1)/3,mp.mpf(2)/3]:
    av=cusp*cusp+delta*delta;bv=cusp*delta;kv=(cusp-delta)**2;qv=2*bv/av
    # Full positive series via six endpoint probabilities, no 2^n enumeration.
    row=[mp.mpf(1),mp.mpf(0),mp.mpf(0),mp.mpf(0),mp.mpf(0),mp.mpf(0)]
    accum=[mp.mpf(0)]*6;deriv=[mp.mpf(0)]*6
    mass=mp.mpf(0);coef=mp.sqrt(1-qv);n=0
    while True:
        score=-delta/av+n*(1/delta-2*delta/av)
        for j in range(6):
            accum[j]+=coef*row[j]
            deriv[j]+=coef*row[j]*score
        mass+=coef
        tail=qv**(n+1)/mp.sqrt(1-qv)
        if tail<mp.mpf('1e-40'):break
        n+=1;coef*=qv*(n-mp.mpf('.5'))/n
        row=[(row[(j-1)%6]+row[(j+1)%6])/2 for j in range(6)]
    eig=[1/mp.sqrt(av-2*bv*c) for c in cos]
    krow=[sum(eig[j]*mp.cos(mp.pi*j*r/3) for j in range(6))/6 for r in range(6)]
    error=max(abs(accum[j]/mp.sqrt(kv)-krow[j]) for j in range(6))
    need(error<=tail/mp.sqrt(kv)+mp.mpf('1e-48'),f'original cusp{cusp} whole half-word series with tail')
    need(abs(mass-1)<=tail+mp.mpf('1e-48'),f'original cusp{cusp} normalized length mass')
    drow=[-sum((delta-cusp*cos[j])/(av-2*bv*cos[j])**mp.mpf('1.5')*mp.cos(mp.pi*j*r/3) for j in range(6))/6 for r in range(6)]
    dtail=(mp.sqrt(1-qv)/mp.sqrt(kv)*qv**(n+1)*(abs(delta/av)/(1-qv)+abs(1/delta-2*delta/av)*((n+1)-n*qv)/(1-qv)**2))
    need(max(abs(deriv[j]/mp.sqrt(kv)-drow[j]) for j in range(6))<=dtail+mp.mpf('1e-45'),'whole absolute delta word derivative matches source operator')
    # Independent original cubic Green formula on same native cosine matrix.
    Cmp=mp.matrix([[mp.mpf(str(C[i,j])) for j in range(6)] for i in range(6)])
    Imp=mp.eye(6)
    f=lambda c:1/mp.sqrt(av-2*bv*c)
    E1=(f(1)+f(-1))/2;Eh=(f(mp.mpf('.5'))+f(mp.mpf('-.5')))/2
    O1=(f(1)-f(-1))/2;Oh=(f(mp.mpf('.5'))-f(mp.mpf('-.5')))/2
    oldsqrt=(4*Eh-E1)/3*Imp+(8*Oh-O1)/3*Cmp+4*(E1-Eh)/3*Cmp**2+4*(O1-2*Oh)/3*Cmp**3
    need(max(abs(oldsqrt[0,r]-krow[r]) for r in range(6))<mp.mpf('1e-48'),'old whole square-root interpolation independently agrees')
    poly=((av**3-5*av*bv*bv)*Imp+2*bv*(av*av-5*bv*bv)*Cmp+4*av*bv*bv*Cmp**2+8*bv**3*Cmp**3)/((av*av-bv*bv)*(av*av-4*bv*bv))
    need(max(abs(poly[0,r]-sum(1/(av-2*bv*cos[j])*mp.cos(mp.pi*j*r/3) for j in range(6))/6) for r in range(6))<mp.mpf('1e-46'),'old cubic Green agrees at actual cusp')
    # Square literal endpoint row -> Green, retaining original amplitude.
    need(max(abs(sum(krow[t]*krow[(r-t)%6] for t in range(6))-poly[0,r]) for r in range(6))<mp.mpf('1e-46'),'two original square-root factors sew to Green')
    loop=mp.fsum(qv**m/m*(2+4*mp.power(2,-m)) for m in range(2,n+1,2))
    F=-mp.fsum(mp.log(av-2*bv*c) for c in cos)
    looptail=6*qv**(n+1)/((n+1)*(1-qv))
    need(abs(-6*mp.log(av)+loop-F)<=looptail+mp.mpf('1e-47'),'same original determinant closes all loops with1/n')
    records.append({'y':str(cusp),'delta_original':str(delta),'q':str(qv),'kappa':str(kv),'mean_length_half':str(qv/(2*(1-qv))),'terms':n+1,'tail_bound':str(tail),'absolute_K_error':str(error),'half_length_first_120':[float(mp.sqrt(1-qv)*mp.rf(mp.mpf('.5'),j)/mp.factorial(j)*qv**j) for j in range(120)]})
out={'result':'PASS','checks':len(checks),'checks_detail':checks,'native_builder':'verification/v17_hexagonal_resolvent.py:U6_power','original_cusps':records,'scope':'Canonical whole-word expansion of bare original component; not primitive physical probability law','complete_TFPT_solution':False,'ledger_changed':False}
Path(__file__).with_name('C6-WHOLE-WORD-PRODUCER.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'result':out['result'],'checks':len(checks),'cusps':[{k:v for k,v in r.items() if k not in ['half_length_first_120','delta_original']} for r in records],'complete_TFPT_solution':False},indent=2))
