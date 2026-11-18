"""Evolve the ACTUAL conformally varied comparison Dirac operator.

Probe metric g_eps=exp(2 eps h(w))(dx^2+dw^2). After unitary density transport,
D_eps=exp(-eps h/2) D_0 exp(-eps h/2), on the original fixed Pin domain.
The profile and eps are controlled inputs, not a law selected by P1.
No projection/reset occurs between steps. Finite normal-mode Galerkin runs
are compared with an analytic all-mode first-response formula.
"""
import argparse
import json
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent


def require(condition,message):
    if not condition:
        raise RuntimeError(message)


def profile_values(w,profile):
    if profile=='quadratic':return w*w,2*w
    # A smooth even bump, zero near w=0 and both original physical boundaries.
    # Its location/width are diagnostic inputs, not TFPT predictions.
    r=(w-.5)/.25
    f=np.zeros_like(w);derivative=np.zeros_like(w)
    inside=abs(r)<1
    rr=r[inside]
    f[inside]=np.exp(1-1/(1-rr*rr))
    derivative[inside]=f[inside]*(-2*rr/(1-rr*rr)**2)/.25
    return f,derivative


def dirac_block(nmax, epsilon, k, b=1,profile='collar_bump'):
    # In the original Pin sector containing the normal zero branch, the
    # allowed transverse H indices are even: E_m=pi*m, J=b.
    modes=np.arange(-nmax,nmax+1)
    nodes,weights=np.polynomial.legendre.leggauss(max(512,8*nmax))
    w=(nodes+1)/2; weights=weights/2
    orders=np.arange(-2*nmax,2*nmax+1)
    f,_=profile_values(w,profile)
    moments=np.cos(np.pi*orders[:,None]*w[None,:])@(weights*np.exp(-epsilon*f))
    hankel=moments[(modes[:,None]+modes[None,:])+2*nmax]
    d=b*(k+0.5j*np.pi*(modes[:,None]-modes[None,:]))*hankel
    require(np.max(abs(d-d.conj().T))<1e-12,'Hermitian original Dirac block')
    return modes,d,w,weights


def old_covariance(modes,k,b=1):
    nmax=len(modes)//2
    c=np.zeros((len(modes),len(modes)),complex)
    c[nmax,nmax]=float(b*k<0)
    for m in range(1,nmax+1):
        indices=[nmax+m,nmax-m]
        d=b*np.array([[0,k+1j*np.pi*m],[k-1j*np.pi*m,0]])
        c[np.ix_(indices,indices)]=(np.eye(2)-d/np.sqrt(k*k+np.pi*np.pi*m*m))/2
    require(np.linalg.norm(c@c-c)<1e-12,'Original pure D sea')
    return c


def run(nmax,epsilon,k,times,profile='collar_bump'):
    modes,d,w,weights=dirac_block(nmax,epsilon,k,profile=profile)
    c0=old_covariance(modes,k)
    eigenvalues,q=np.linalg.eigh(d)
    initial=np.zeros(len(modes),complex);initial[nmax]=1
    coefficients=q.conj().T@initial
    # e^(itD) on test functions gives the usual e^(-itE) annihilator phase.
    vectors=q@(np.exp(1j*eigenvalues[:,None]*times[None,:])*coefficients[:,None])
    norm=np.sum(abs(vectors)**2,axis=0)
    occupancy=np.sum(vectors.conj()*(c0@vectors),axis=0).real
    require(np.max(abs(norm-1))<1e-11,'Unitary source recursion')
    require(occupancy.min()>-1e-12 and occupancy.max()<1+1e-12,'CAR covariance positivity')
    dt=times[1]-times[0]
    step=(q*np.exp(1j*dt*eigenvalues)[None,:])@q.conj().T
    repeated=initial.copy()
    for _ in range(len(times)-1):repeated=step@repeated
    recurrence_error=float(np.linalg.norm(repeated-vectors[:,-1]))
    require(recurrence_error<1e-10,'Repeated source steps versus direct evolution')
    # Eliminate every exterior mode algebraically, without projecting/resetting
    # the physical recursion. For unitary U, C=P U Q need not equal B^*.
    outside=np.arange(len(modes))!=nmax
    K=step[nmax,nmax]
    B=step[outside,nmax]
    C=step[nmax,outside]
    E=step[np.ix_(outside,outside)]
    memory=[];wrong_memory=[];v=B.copy()
    for _ in range(len(times)-1):
        memory.append(C@v)
        wrong_memory.append(B.conj()@v)
        v=E@v
    def reduced(kernel):
        x=[1+0j]
        for n in range(len(times)-1):
            x.append(K*x[-1]+sum(kernel[n-1-j]*x[j] for j in range(n)))
        return np.asarray(x)
    reduced_response=reduced(memory)
    memory_error=float(np.max(abs(reduced_response-vectors[nmax])))
    wrong_memory_error=float(np.max(abs(reduced(wrong_memory)-vectors[nmax])))
    require(memory_error<1e-10,'Derived full source memory versus full evolution')
    if abs(epsilon)>=.01:
        require(wrong_memory_error>1e-7,'Incorrect adjoint memory control did not fire')
    column=d[:,nmax].copy();column[nmax]=0
    measured_return=float(np.linalg.norm(column)**2)
    f,df=profile_values(w,profile)
    mean=float(weights@np.exp(-epsilon*f))
    mean2=float(weights@np.exp(-2*epsilon*f))
    normal2=float(weights@(df*df*np.exp(-2*epsilon*f)))
    exact_return=k*k*(mean2-mean*mean)+epsilon*epsilon*normal2/4
    weak_coefficient=k*k*(weights@(f*f)-(weights@f)**2)+(weights@(df*df))/4
    return dict(nmax=nmax,mode_count=len(modes),epsilon=epsilon,k=k,
                return_form_finite=measured_return,return_form_all_mode=exact_return,
                return_form_missing_tail=exact_return-measured_return,
                weak_geometry_return_coefficient=float(weak_coefficient),
                recurrence_error=recurrence_error,norm_error=float(np.max(abs(norm-1))),
                derived_memory_error=memory_error,incorrect_adjoint_memory_error=wrong_memory_error,
                memory_kernel=[[float(z.real),float(z.imag)] for z in memory],
                times=times.tolist(),zero_branch_survival=(abs(vectors[nmax])**2).tolist(),
                original_sea_zero_mode_occupancy=occupancy.tolist(),
                final_occupancy=float(occupancy[-1]))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--epsilon',type=float,default=.1)
    parser.add_argument('--profile',choices=['collar_bump','quadratic'],default='collar_bump')
    parser.add_argument('--output',type=Path,default=HERE/'runs/source_geometric_machine.json')
    args=parser.parse_args()
    require(abs(args.epsilon)<=.5,'Use a bounded geometric probe |epsilon|<=0.5')
    times=np.linspace(0,4,401)
    k=np.pi/2  # Actual b=+1,c=+i old quarter-circle lowest positive momentum.
    runs=[run(n,args.epsilon,k,times,args.profile) for n in [16,32,64,128]]
    convergence=[]
    for a,b in zip(runs,runs[1:]):
        convergence.append(dict(n_from=a['nmax'],n_to=b['nmax'],
            max_survival_difference=float(np.max(abs(np.array(a['zero_branch_survival'])-b['zero_branch_survival']))),
            max_occupancy_difference=float(np.max(abs(np.array(a['original_sea_zero_mode_occupancy'])-b['original_sea_zero_mode_occupancy'])))))
    result=dict(source='Actual fixed Pin Dirac conformal metric variation; same initial pure D sea',
                metric_profile=args.profile+'; explicitly chosen diagnostic profile',
                native_input='b=+1,c=+i,q=i,k=pi/2; all compatible transverse modes retained up to stated Galerkin bound',
                analytic_zero_geometry_return_coefficient='k^2 Var(f)+<f_prime^2>/4',
                runs=runs,convergence=convergence,
                scope='All-mode first-response form proved analytically; finite-mode unitary CAR evolution is numerically tested. No primitive history selection or P1 physical map claimed.',
                physical_gates_closed=[],complete_TFPT_solution=False)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(runs=[{key:value for key,value in x.items() if key not in ['times','zero_branch_survival','original_sea_zero_mode_occupancy','memory_kernel']} for x in runs],convergence=convergence),indent=2))


if __name__=='__main__':main()
