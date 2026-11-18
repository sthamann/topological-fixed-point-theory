"""Execute old native v990 caller; preserve its parameters and actual local builders."""
from pathlib import Path
import sys,io,contextlib,json
import numpy as np
repo=Path(__file__).resolve().parents[4]
sys.path.insert(0,str(repo/'verification'))
import v990_wj_transduction_mechanism as v
locals990={}
def capture(frame,event,arg):
    if event=='return' and frame.f_code.co_name=='run' and frame.f_code.co_filename==v.__file__:
        locals990.update(frame.f_locals)
old=sys.getprofile();log=io.StringIO()
try:
    sys.setprofile(capture)
    with contextlib.redirect_stdout(log):v.run()
finally:sys.setprofile(old)
if '[FAIL]' in log.getvalue() or '13 passed, 0 failed' not in log.getvalue():
    raise RuntimeError('Original v990 checks failed')
H=locals990['H_W0'];O=locals990['O_Z'];beta=v.BETA
E,V=np.linalg.eigh(H);Oe=V.conj().T@O@V
Es=E-E.min();rawp=np.exp(-beta*Es);p=rawp/rawp.sum()
mean=np.dot(p,np.diag(Oe));equal=np.dot(p,np.diag(Oe@Oe))-mean**2
kubo=v.kubo_connected(Oe,Oe,p,E)
def connected(tau):
    return np.sum(p[:,None]*np.exp(tau*(E[:,None]-E[None,:]))*Oe*Oe.T)-mean**2
def full_word(times,ops):
    if len(times)!=len(ops) or any(times[j]>times[j+1] for j in range(len(times)-1)):
        raise ValueError('ordered times/fields required')
    if not times:return 1.+0j
    if times[0]<0 or times[-1]>beta:raise ValueError('original thermal interval')
    out=np.diag(np.exp(-(beta-times[-1])*Es)).astype(complex)
    for j in range(len(times)-1,-1,-1):
        dt=times[j]-(times[j-1] if j else 0)
        out=out@(V.conj().T@ops[j]@V)@np.diag(np.exp(-dt*Es))
    return np.trace(out)/rawp.sum()
checks=[]
def need(ok,name):
    if not ok:raise RuntimeError(name)
    checks.append(name)
Ps=[v.w_ferm(i,4,6) for i in range(4)]
need(max(np.linalg.norm(a@b.conj().T+b.conj().T@a-(np.eye(96) if i==j else 0)) for i,a in enumerate(Ps) for j,b in enumerate(Ps))<1e-12,'actual four-mode CAR representation')
vac=np.eye(96,dtype=complex)
for a in Ps:vac=vac@(np.eye(96)-a.conj().T@a)
need(np.linalg.matrix_rank(vac,tol=1e-10)==6 and np.linalg.norm(vac@vac-vac)<1e-12,'native CAR vacuum multiplicity6')
nodes,weights=np.polynomial.legendre.leggauss(64)
integrated=sum(w*connected(beta*(x+1)/2)/2 for x,w in zip(nodes,weights))
need(abs(integrated-kubo)<1e-12,'original Kubo equals full time integral')
need(abs(equal-kubo)>1e-6,'instantaneous covariance is not original Hessian')
need(max(abs(full_word([0,t],[O,O])-mean**2-connected(t)) for t in [0,beta/3,beta/2,beta])<1e-12,'actual full timed trace matches original spectral reader')
need(abs(connected(beta/3)-connected(2*beta/3))<1e-12,'whole original thermal two-point symmetry')
data={'result':'PASS','native_v990_checks':13,'additional_controls':len(checks),'checks':checks,'space':'C6_link tensor Fock(C4),not internal96 one-particle SM','common_CAR_vacuum_rank':6,'beta':beta,'mean_OZ':float(mean.real),'equal_time_connected':float(equal.real),'integrated_Kubo_connected':float(kubo.real),'gap_equal_vs_integrated':float((equal-kubo).real),'direct_integral_error':float(abs(integrated-kubo)),'midtime_connected':float(connected(beta/2).real),'noncommuting_H_OZ_HS':float(np.linalg.norm(H@O-O@H)),'full_four_word_at_original_beta_quarters':float(full_word([0,beta/4,beta/2,3*beta/4],[O]*4).real),'time_curve':[[float(t/beta),float(connected(t).real)] for t in np.linspace(0,beta,81)],'complete_TFPT_solution':False,'ledger_changed':False}
Path(__file__).with_name('RECOVERED-NATIVE-SOURCE.json').write_text(json.dumps(data,indent=2)+'\n')
Path(__file__).with_name('RECOVERED-NATIVE-SOURCE-v990-run.txt').write_text(log.getvalue())
print(json.dumps({k:z for k,z in data.items() if k!='time_curve'},indent=2))
