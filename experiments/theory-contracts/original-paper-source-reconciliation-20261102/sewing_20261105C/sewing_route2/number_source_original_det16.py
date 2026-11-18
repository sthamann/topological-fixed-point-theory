"""Whole nonlinear original R104 N source, finite exact reduction, native check.

All sources are unnormalized and use the original H=2(I-P_phi), not a new
Hamiltonian or a selected physical temperature. Profiles are in written
left-to-right operator order. Complex h is an analytic source insertion.
"""
from pathlib import Path
import json
import math
import sys

import numpy as np
from scipy import linalg,sparse
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import expm_multiply

sys.path.insert(0,str(Path(__file__).parent))
import twined_original_det16 as source
original=source.original
NATIVE_COMPONENT_INVENTORY=set()


def corner_hamiltonian():
    p=original.phi
    return np.array([[1,-np.exp(-1j*p)],[-np.exp(1j*p),1]],dtype=complex)


def corner_source(profile):
    """Factors in listed left-to-right order; reverse chronological input."""
    result=np.eye(2,dtype=complex)
    n=np.diag([0.,16.])
    h2=corner_hamiltonian()
    for dt,h in profile:result=result@linalg.expm(-dt*(h2-h*n))
    return result


def reduced(profile,r):
    beta=sum(dt for dt,h in profile)
    chi=sum(dt*h for dt,h in profile)
    det_r=np.linalg.det(r)
    remaining=math.exp(-2*beta)*(np.linalg.det(np.eye(16)+np.exp(chi)*r)-1-np.exp(16*chi)*det_r)
    corner=np.trace(np.diag([1,det_r])@corner_source(profile))
    return complex(remaining+corner)


def center_reduced(profile):
    beta=sum(dt for dt,h in profile)
    chi=sum(dt*h for dt,h in profile)
    return complex(math.exp(-2*beta)*(1820*np.exp(4*chi)+12870*np.exp(8*chi)+1820*np.exp(12*chi))+np.trace(corner_source(profile)))


def constant_center(beta,h):
    chi=beta*h
    return complex(math.exp(-2*beta)*(1820*np.exp(4*chi)+12870*np.exp(8*chi)+1820*np.exp(12*chi))
                   +2*np.exp(-beta*(1-8*h))*np.cosh(beta*np.sqrt(1+64*complex(h)**2)))


def native_full_source(profile,H,N):
    # FULL ORIGINAL native sparse input. Its actual sparse graph determines the
    # components, without using the displayed H2 or any character formula.
    # Generic sparse Pade expm would solve 65536 sparse right-hand sides.
    result=original.I.astype(complex)
    for dt,h in profile:
        A=(-dt*(H-h*N)).tocsr()
        graph=A.copy();graph.setdiag(0);graph.eliminate_zeros()
        graph.data=np.ones(graph.nnz,dtype=np.int8)
        count,labels=connected_components(graph,directed=False)
        sizes=np.bincount(labels)
        nonsingle=np.flatnonzero(sizes>1)
        NATIVE_COMPONENT_INVENTORY.add((int(count),tuple(int(sizes[x]) for x in nonsingle)))
        diagonal=np.exp(A.diagonal())
        correction_rows=[];correction_cols=[];correction_values=[]
        for component in nonsingle:
            indices=np.flatnonzero(labels==component)
            block=linalg.expm(A[np.ix_(indices,indices)].toarray())-np.diag(diagonal[indices])
            for i,row in enumerate(indices):
                for j,col in enumerate(indices):
                    correction_rows.append(row);correction_cols.append(col);correction_values.append(block[i,j])
        factor=sparse.diags(diagonal,format='csr')+sparse.coo_matrix((correction_values,(correction_rows,correction_cols)),shape=A.shape).tocsr()
        result=result@factor
    return result.tocsr()


def zcomplex(z):return [float(complex(z).real),float(complex(z).imag)]


def main():
    original.initialize_native_controls()
    occupancy=np.fromiter((i.bit_count() for i in range(original.D)),dtype=np.int8,count=original.D)
    N=sparse.diags(occupancy.astype(float),format='csr')
    H=2*(original.I-original.P)
    profiles=[('zero',[(1.,0.)]),('real_positive',[(1.,.1)]),
              ('real_negative',[(1.,-.2)]),('complex',[(1.,.1+.07j)]),
              ('zero_discriminant',[(1.,1j/8)]),
              ('nonconstant_real',[(.2,.1),(.3,-.2),(.5,.3)]),
              ('nonconstant_complex',[(.2,.1+.07j),(.3,-.05+.02j),(.5,.12-.09j)])]
    twist_diagonals={label:q**occupancy for label,q in source.TWISTS.items()}
    center_diagonal=(occupancy%4==0).astype(float)
    generators,pairs=original.native.positive_weyl_generators()
    G=generators[pairs.index((0,1))]
    angle=.2718281828459045 # same original v1002 diagnostic gauge insertion
    r_cartan=linalg.expm(1j*angle*G)
    lift_cartan=np.exp(1j*angle*np.fromiter((sum(G[k,k] for k in range(16) if state>>k&1) for state in range(original.D)),dtype=complex,count=original.D))
    checks={};rows=[];maximum=0.;maxrel=0.;comparisons=0
    checks['original_H_corner_matches_separate_H2']=bool(np.max(np.abs(H[np.ix_([0,original.D-1],[0,original.D-1])].toarray()-corner_hamiltonian()))<1e-14)
    checks['actual_Cartan_gauge_unitary_det_one']=bool(np.max(np.abs(r_cartan.conj().T@r_cartan-np.eye(16)))<1e-14 and abs(np.linalg.det(r_cartan)-1)<1e-14)
    for name,profile in profiles:
        native=native_full_source(profile,H,N)
        diagonal=native.diagonal()
        comparisons_this=[]
        for label,q in source.TWISTS.items():
            direct=np.dot(twist_diagonals[label],diagonal)
            value=reduced(profile,q*np.eye(16,dtype=complex))
            error=float(abs(direct-value));comparisons+=1
            maximum=max(maximum,error);maxrel=max(maxrel,error/max(1.,abs(direct)))
            checks[name+'_'+label]=error<2e-10*max(1.,abs(direct))
            comparisons_this.append({'twist':label,'native':zcomplex(direct),'reduced':zcomplex(value),'absolute_error':error})
        direct=np.dot(center_diagonal,diagonal)
        value=center_reduced(profile)
        error=float(abs(direct-value));comparisons+=1
        maximum=max(maximum,error);maxrel=max(maxrel,error/max(1.,abs(direct)))
        checks[name+'_P4']=error<2e-10*max(1.,abs(direct))
        comparisons_this.append({'twist':'P4','native':zcomplex(direct),'reduced':zcomplex(value),'absolute_error':error})
        checks[name+'_center_is_actual_group_average']=abs(value-sum(reduced(profile,q*np.eye(16)) for q in source.TWISTS.values())/4)<2e-10*max(1.,abs(value))
        direct_cartan=np.dot(lift_cartan,diagonal)
        value_cartan=reduced(profile,r_cartan)
        error=float(abs(direct_cartan-value_cartan));comparisons+=1
        maximum=max(maximum,error);maxrel=max(maxrel,error/max(1.,abs(direct_cartan)))
        checks[name+'_actual_noncentral_Cartan']=error<2e-10*max(1.,abs(direct_cartan))
        comparisons_this.append({'twist':'actual noncentral Cartan','native':zcomplex(direct_cartan),'reduced':zcomplex(value_cartan),'absolute_error':error})
        if len(profile)==1:
            beta,h=profile[0]
            checks[name+'_constant_source_cosh_formula']=abs(value-constant_center(beta,h))<2e-10*max(1.,abs(value))
        checks[name+'_full_corner_operator_equality']=bool(np.max(np.abs(native[np.ix_([0,original.D-1],[0,original.D-1])].toarray()-corner_source(profile)))<2e-11*max(1.,np.max(np.abs(corner_source(profile)))))
        rows.append({'profile':name,'segments':[{'dt':dt,'h':zcomplex(h)} for dt,h in profile],'beta':sum(dt for dt,h in profile),'chi':zcomplex(sum(dt*h for dt,h in profile)),'comparisons':comparisons_this})
    # Exact equal-time and two-time source derivatives as the kick-source
    # coefficients, without a finite Taylor truncation in the h-family.
    s,t=.3,.7;a,b=math.exp(-2*s),math.exp(-2*t)
    E_s=a*original.I-math.expm1(-2*s)*original.P
    E_t=b*original.I-math.expm1(-2*t)*original.P
    native_N2=np.dot(center_diagonal,(E_s@N@E_t@N).diagonal())
    phi2=128*((1-a)*b+a*(1-b))+64*(1-a)*(1-b)
    expected_N2=1115136*a*b+phi2
    checks['actual_native_N2_source_coefficient']=abs(native_N2-expected_N2)<2e-10*max(1.,abs(native_N2))
    checks['N2_kick_corner_coefficient']=abs(np.trace(corner_source([(s,0.)])@np.diag([0.,16.])@corner_source([(t,0.)])@np.diag([0.,16.]))-(256*a*b+phi2))<1e-12
    checks['zero_source_original_center_partition']=abs(center_reduced([(1.,0.)])-(1+16511*math.exp(-2)))<1e-10
    checks['zero_width_original_center_rank']=abs(center_reduced([(0.,0.)])-16512)<1e-10
    # Equal-time Tr(P4 e^-betaH N²), a distinct limit from time integration.
    eq=np.dot(center_diagonal,(native_full_source([(1.,0.)],H,N)@(N@N)).diagonal())
    checks['actual_equal_time_N_squared']=abs(eq-(1115136*math.exp(-2)+128*(1-math.exp(-2))))<2e-10*max(1.,abs(eq))
    # A second native method: sparse exponential ACTION on the full 65536D
    # operator, with no component split, for a nonconstant complex profile.
    profile=profiles[-1][1]
    basis=np.zeros((original.D,2),dtype=complex);basis[[0,original.D-1],[0,1]]=1
    action=basis.copy()
    for dt,h in reversed(profile):
        action=expm_multiply((-dt*(H-h*N)).tocsr(),action)
    product=native_full_source(profile,H,N)
    checks['whole_native_sparse_exponential_action_matches_component_product']=bool(np.max(np.abs(action-product@basis))<2e-11*max(1.,np.max(np.abs(action))))
    derivative_errors=[float(abs(native_N2-expected_N2)),float(abs(eq-(1115136*math.exp(-2)+128*(1-math.exp(-2)))))]
    profile_maximum=maximum
    maximum=max(maximum,*derivative_errors)
    checks={k:bool(v) for k,v in checks.items()};assert all(checks.values()),{k:v for k,v in checks.items() if not v}
    report={'status':'pass','scope':'Whole finite original R104 nonlinear gauge-invariant N source, original H=2(I-P_phi), canonical Gamma(r) boundary; no primitive S_adm/v1013 gauge-map or physical beta/h calibration.',
            'whole_profile_formula':'exp(-2beta)[det(I+exp(chi)r)-1-exp(16chi)det(r)] + Tr2[diag(1,det r) product_j exp(-dt_j(H2-h_jN2))]',
            'profile_order':'Listed left-to-right OPERATOR product order. For source samples listed chronologically early-to-late, reverse the list for conventional time ordering (later times on the left).',
            'center_profile_formula':'exp(-2beta)[1820exp(4chi)+12870exp(8chi)+1820exp(12chi)]+Tr2 product_j exp(-dt_j(H2-h_jN2))',
            'constant_center_corner_formula':'2 exp[-beta(1-8h)] cosh[beta sqrt(1+64h²)] (entire even dependence, including discriminant zero)',
            'native_full_dimension':original.D,'native_full_sparse_exponential_comparisons':comparisons,'checks_passed':sum(checks.values()),'checks_total':len(checks),
            'checked':sum(checks.values()),'total':len(checks),'native_comparisons':comparisons,
            'native_exponential_method':'Entire actual sparse matrix; graph-derived connected components, plus independent full sparse exponential action on nonconstant complex profile. No displayed H2 or character formula used for native trace.',
            'native_discovered_component_inventories':[{'component_count':n,'nontrivial_sizes':list(v)} for n,v in sorted(NATIVE_COMPONENT_INVENTORY)],
            'max_absolute_error':maximum,'max_scaled_error':maxrel,'checks':checks,'profiles':rows,
            'profile_max_absolute_error':profile_maximum,'derivative_native_comparisons':2,'derivative_max_absolute_error':max(derivative_errors),
            'N2_source_coefficients':{'s':s,'t':t,'native':zcomplex(native_N2),'closed':zcomplex(expected_N2),'equal_time_native':zcomplex(eq)},
            'antecedents_retained':['Actual original full rest/gauge representation and physical field map','Actual bound reference denominator, never reset to1','Same thermal/Euclidean time and source ordering']}
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ['profiles','checks']},indent=2))


if __name__=='__main__':main()
