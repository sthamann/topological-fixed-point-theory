"""Deterministic integration of the existing finite R048 radial reweighting.

Requires the actual original finite-regulator quadratic precision, mean,
cell volumes and bound original lambda/vev. No fallback model/default exists.
precision is the matrix of the original quadratic FORM in field coordinates,
including original cell factors, rather than a bare differential expression.
The Gaussian baseline Z_quad[q] stays with its caller; this returns its
dimensionless multiplicative reweighting, and deterministic analytic error.
The integration domain is the supplied R048 Gaussian R^N domain. The local
radial gauge chart does not establish a global Higgs domain/Jacobian. Actual
domain indicators or nonpolynomial angular/ghost/Jacobian factors need their
own retained integration/error control; they are not silently included here.

Floating point evaluation is not interval arithmetic. The analytic truncation
and cell error bound is rigorous; hardware rounding needs separate enclosure
before the computed floating point endpoint is called a certified interval.
"""
from pathlib import Path
import argparse, itertools, json, math
import numpy as np
from scipy.special import ndtr

def add_term(poly,power,value):
    poly[power]=poly.get(power,0)+value
    if poly[power]==0:del poly[power]

def shift_to_x(poly,a,n):
    """Exact algebraic recipe eta=x-3a; floats are NOT interval endpoints."""
    out={}
    for power,coef in poly.items():
        for ks in itertools.product(*(range(p+1) for p in power)):
            c=coef
            for p,k in zip(power,ks):c*=math.comb(p,k)*(-3*a)**(p-k)
            add_term(out,tuple(ks),c)
    return out

def derivative(poly,j):
    out={}
    for power,coef in poly.items():
        if power[j]:
            new=list(power);new[j]-=1
            add_term(out,tuple(new),coef*power[j])
    return out

def polynomial_envelope(poly,c):
    def mk(k,j):return 1 if k==0 else (k/(2*math.e*c[j]))**(k/2)
    return math.fsum(abs(coef)*math.prod(mk(k,j)for j,k in enumerate(power))for power,coef in poly.items())

def insertion_bounds(poly_eta,a,beta,c,largest_cov_eigenvalue):
    n=len(beta);poly=shift_to_x(poly_eta,a,n)
    B=polynomial_envelope(poly,c);D=[]
    for j in range(n):
        dp=derivative(poly,j)
        for power,coef in poly.items():
            # d/dx_j exp(-W)=-beta_j*(x_j^3-6a*x_j^2+9a^2*x_j)*exp(-W).
            for degree,factor in ((3,1),(2,-6*a),(1,9*a*a)):
                new=list(power);new[j]+=degree
                add_term(dp,tuple(new),-coef*beta[j]*factor)
        D.append(polynomial_envelope(dp,c))
    return B,math.sqrt(largest_cov_eigenvalue)*float(np.linalg.norm(D))

def original_insertions(data,n,vev):
    one={(0,)*n:1.0};out={'Z':one}
    # Actual neutral field polynomial of the given R048 coordinate chart.
    for i in data.get('u_indices',list(range(n))):
        if not (0<=i<n):raise ValueError('Requested u source is not an original cell coordinate.')
        p={(0,)*n:vev*vev};power=[0]*n;power[i]=1;p[tuple(power)]=math.sqrt(2)*vev
        power[i]=2;p[tuple(power)]=0.5;out[f'u[{i}]']=p
    for name,terms in data.get('original_polynomial_insertions',{}).items():
        if name in out:raise ValueError('Insertion name collides with an original neutral source.')
        p={}
        for term in terms:
            powers=tuple(term['powers']);coeff=term['coefficient']
            if len(powers)!=n or any(not isinstance(k,int)or k<0 for k in powers):
                raise ValueError('Original insertion powers must be nonnegative integral multiindices.')
            if isinstance(coeff,list):coeff=complex(*coeff)
            add_term(p,powers,coeff)
        out[name]=p
    requested=data.get('normalize_insertions')
    if requested is not None:
        if not isinstance(requested,list) or any(not isinstance(name,str) for name in requested):
            raise ValueError('normalize_insertions must list supplied original insertion names.')
        missing=set(requested)-set(out)
        if missing:
            raise ValueError('Requested normalized insertion has no supplied original polynomial: '+', '.join(sorted(missing)))
    denominator=data.get('original_denominator_polynomial','Z')
    if denominator not in out:
        raise ValueError('Original denominator polynomial was not provided.')
    return out

def poly_value(poly,eta):
    return sum(coef*math.prod(float(x)**k for x,k in zip(eta,power)) for power,coef in poly.items())

def plan_original_integral(data):
    precision=np.asarray(data['original_quadratic_precision'],dtype=float)
    mean=np.asarray(data['original_quadratic_mean'],dtype=float)
    volumes=np.asarray(data['original_cell_volumes'],dtype=float)
    lam=float(data['original_lambda']);vev=float(data['original_vev'])
    eps=float(data['absolute_error_for_bounded_expectation'])
    if not (0<eps<1 and lam>0 and vev>0 and math.isfinite(lam) and math.isfinite(vev)):
        raise ValueError('Requires the original positive broken branch and 0<error<1.')
    n=len(mean)
    if n<1 or precision.shape!=(n,n) or volumes.shape!=(n,):
        raise ValueError('Original coordinate dimensions disagree.')
    if not np.all(np.isfinite(precision)) or not np.all(np.isfinite(mean)) or not np.all(np.isfinite(volumes)) or not np.all(volumes>0):
        raise ValueError('Original inputs must be finite with positive cell volumes.')
    if not np.allclose(precision,precision.T,rtol=1e-13,atol=1e-14):
        raise ValueError('Supply the original symmetric quadratic FORM, with cell factors.')
    eig=np.linalg.eigvalsh(precision)
    if eig[0]<=0:raise ValueError('The existing Gaussian caller requires positive precision.')
    cov=np.linalg.inv(precision);L=np.linalg.cholesky(cov)
    a=math.sqrt(2)*vev;beta=lam*volumes;c=beta*a*a/2
    insertions=original_insertions(data,n,vev)
    bounds={name:insertion_bounds(poly,a,beta,c,1/float(eig[0]))for name,poly in insertions.items()}
    Bmax=max(1,max(B for B,Lp in bounds.values()))
    Lmax=max(Lp for B,Lp in bounds.values())
    T=math.sqrt(2*math.log(4*n*Bmax/eps))
    requested_h=eps/(Lmax*math.sqrt(n))
    cells=max(1,math.ceil(2*T/requested_h));h=2*T/cells
    C=float(np.sum(27*beta*a**4/4))
    diag=np.diag(cov)
    EV=float(np.sum(beta*(a*(mean**3+3*mean*diag)+(mean**4+6*mean**2*diag+3*diag**2)/4)))
    return dict(n=n,cells=cells,h=h,T=T,C=C,EV=EV,insertion_bounds=bounds,
                truncation_bound=2*n*math.exp(-T*T/2),
                covariance=cov,L=L,mean=mean,beta=beta,a=a,insertions=insertions)

def integrate_original_integral(data, execute=False):
    p=plan_original_integral(data)
    out={k:v for k,v in p.items() if k not in ('covariance','L','mean','beta','insertions')}
    out['evaluations']=p['cells']**p['n']
    errors={name:Lp*p['h']*math.sqrt(p['n'])/2+B*p['truncation_bound']for name,(B,Lp)in p['insertion_bounds'].items()}
    out['analytic_absolute_errors']=errors
    out['log_lower_bound_original_reweighting']=-p['EV']
    out['log_upper_bound_original_reweighting']=p['C']
    out['baseline']='Multiply by the same original unnormalized Z_quad[q], retaining its reference/eta.'
    out['rounding_enclosed']=False
    out['scope']='Existing finite fixed-background R048 radial interaction only; S_extra and other fields retained outside.'
    out['domain_scope']='Given unrestricted R048 Gaussian domain only. Local gauge chart does not establish global Higgs domain, angles or Jacobian. Extra domain factors require their actual retained control.'
    if not execute:return out
    maxeval=int(data['algorithm_max_evaluations'])
    if out['evaluations']>maxeval:
        raise ValueError(f"Original input needs {out['evaluations']} evaluations for this bound; exceeds supplied algorithm limit {maxeval}.")
    edges=np.linspace(-p['T'],p['T'],p['cells']+1)
    mids=(edges[1:]+edges[:-1])/2
    probs=ndtr(edges[1:])-ndtr(edges[:-1])
    terms={name:[]for name in p['insertions']}
    for indices in itertools.product(range(p['cells']),repeat=p['n']):
        z=mids[list(indices)];eta=p['mean']+p['L']@z
        W=np.sum(p['beta']/4*(eta+3*p['a'])**2*((eta-p['a'])**2+2*p['a']**2))
        weighted=float(np.prod(probs[list(indices)])*math.exp(-float(W)))
        for name,poly in p['insertions'].items():terms[name].append(poly_value(poly,eta)*weighted)
    estimates={name:complex(math.fsum(complex(x).real for x in ts),math.fsum(complex(x).imag for x in ts))for name,ts in terms.items()}
    Q=estimates['Z'].real;err=errors['Z']
    low=max(0,Q-err);high=min(1,Q+err)
    out.update(bounded_expectation_midpoint=Q,
               analytic_interval_before_rounding=[low,high],
               log_original_reweighting_estimate=p['C']+math.log(Q) if Q>0 else None,
               log_original_reweighting_interval_before_rounding=[p['C']+math.log(low) if low>0 else None,p['C']+math.log(high)])
    denomname=data.get('original_denominator_polynomial','Z')
    if denomname not in estimates:raise ValueError('Original denominator polynomial was not provided.')
    qz=estimates[denomname];ez=errors[denomname];lower=abs(qz)-ez
    out['denominator_disc_excludes_zero']=bool(lower>0)
    out['unnormalized_original_insertions']={}
    out['normalized_original_insertions']={}
    normalize=set(data.get('normalize_insertions',[name for name in estimates if name=='Z' or name.startswith('u[')]))
    for name,qp in estimates.items():
        entry={'unnormalized_bounded_estimate':[qp.real,qp.imag],'analytic_absolute_error':errors[name]}
        out['unnormalized_original_insertions'][name]=entry
        if lower>0 and name in normalize:
            ratio=qp/qz
            ratioerror=errors[name]/lower+abs(qp)*ez/(abs(qz)*lower)
            out['normalized_original_insertions'][name]={**entry,'normalized_estimate':[ratio.real,ratio.imag],'analytic_ratio_error':ratioerror}
    return out

def integrate_until_original_tolerance(data):
    """Refine analytic quadrature error, keeping ALL physical inputs fixed.
    This remains a floating point implementation, not an interval certificate.
    A complex original denominator is permitted only when its error disc
    excludes zero. No positivity of a det/Pf-weighted source is fabricated.
    """
    atol=float(data['normalized_absolute_tolerance'])
    rtol=float(data['normalized_relative_tolerance'])
    if atol<0 or rtol<0 or atol+rtol==0:raise ValueError('Supply a positive algorithmic target tolerance.')
    current=dict(data);history=[]
    for rung in range(int(data['algorithm_max_refinements'])):
        out=integrate_original_integral(current,execute=True)
        good=out['denominator_disc_excludes_zero']
        if good:
            for name,item in out['normalized_original_insertions'].items():
                estimate=abs(complex(*item['normalized_estimate']))
                good &= item['analytic_ratio_error'] <= (atol+rtol*estimate)/(1+rtol)
        history.append({'rung':rung,'evaluations':out['evaluations'],'analytic_tolerance_met':bool(good),
                        'denominator_disc_excludes_zero':out['denominator_disc_excludes_zero']})
        if good:
            out['refinement_history']=history
            out['analytic_tolerance_met_before_rounding']=True
            return out
        current['absolute_error_for_bounded_expectation']/=2
    raise ValueError('The supplied algorithm limit was reached without the original denominator/tolerance condition; no value certified.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('actual_original_input',type=Path)
    ap.add_argument('--execute',action='store_true')
    ap.add_argument('--refine',action='store_true')
    args=ap.parse_args()
    data=json.loads(args.actual_original_input.read_text())
    out=integrate_until_original_tolerance(data) if args.refine else integrate_original_integral(data,args.execute)
    print(json.dumps(out,indent=2))
