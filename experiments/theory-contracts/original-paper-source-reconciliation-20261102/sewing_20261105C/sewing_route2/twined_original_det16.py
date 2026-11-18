"""Actual R104 source with its original Spin10 central boundary insertion.

Uses the original CAR ordering, rank-one projector and time units. No physical
S_adm or v1013 gauge-map identity is assumed. This is unnormalized throughout.
"""
from pathlib import Path
from functools import lru_cache
import importlib.util
import hashlib
import json
import math
import sys

import numpy as np
from scipy import sparse

ROOT = Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
SOURCE = ROOT/'experiments/theory-contracts/original-local-source-variable-background-20261102/check_det16_source.py'
sys.path.insert(0,str(SOURCE.parent))
spec = importlib.util.spec_from_file_location('original_r104',SOURCE)
original = importlib.util.module_from_spec(spec)
spec.loader.exec_module(original)
TWISTS = {'1':1+0j,'i':1j,'-1':-1+0j,'-i':-1j}


def gate(letter):
    return isinstance(letter,tuple) and len(letter)==2 and letter[0]=='U'


def apply_gated_word(letters,state):
    amplitude=1+0j
    for letter in reversed(letters):
        if gate(letter):
            amplitude*=letter[1]**state.bit_count()
        else:
            state,sign=original.apply_word([letter],state)
            if state is None:return None,0j
            amplitude*=sign
    return state,amplitude


def phi_coefficient(letters):
    if not any(gate(letter) for letter in letters):
        return original.word_phi(letters)
    value=0j
    for state,ket in ((0,1),(original.D-1,np.exp(1j*original.phi))):
        end,amp=apply_gated_word(letters,state)
        if end==0:value+=amp*ket/2
        elif end==original.D-1:value+=amp*ket*np.exp(-1j*original.phi)/2
    return value


def fock_trace(letters):
    if not any(gate(letter) for letter in letters):
        return original.word_fock_trace(letters)
    modes=sorted({letter[0] for letter in letters if not gate(letter)})
    bits={}
    total_q=1+0j
    for letter in reversed(letters):
        if gate(letter):total_q*=letter[1]
        else:
            mode,dagger=letter
            if mode not in bits:bits[mode]=0 if dagger else 1
    state=sum(bit<<mode for mode,bit in bits.items())
    end,amp=apply_gated_word(letters,state)
    # No division by a determinant, character or (1+q), including q=-1.
    return 0j if end!=state else amp*(1+total_q)**(16-len(modes))


def generalized_renewal(letters,spans,*,phi_eval=None,trace_eval=None):
    """Original cyclic rank-one algorithm on CAR letters and actual U_q."""
    phi_eval=phi_coefficient if phi_eval is None else phi_eval
    trace_eval=fock_trace if trace_eval is None else trace_eval
    k=len(letters)
    if len(spans)!=k+1 or any(not math.isfinite(t) or t<0 for t in spans):
        raise ValueError('one nonnegative original interval around every letter')
    beta=sum(spans)
    if not k:return original.D*math.exp(-2*beta)-math.expm1(-2*beta)
    delta=[spans[0]+spans[-1]]+spans[1:-1]
    doubled=letters+letters
    @lru_cache(maxsize=None)
    def edge(i,j):
        b=-math.expm1(-2*delta[i%k])
        a=math.exp(-2*sum(delta[r%k] for r in range(i+1,j)))
        return b*a*phi_eval(doubled[i:j])
    result=math.exp(-2*beta)*trace_eval(letters)
    for anchor in range(k):
        dp={anchor:1+0j}
        for j in range(anchor+1,k):
            dp[j]=sum(dp[h]*edge(h,j) for h in range(anchor,j))
        result+=sum(dp[h]*edge(h,anchor+k) for h in range(anchor,k))
    return result


def twined_word(letters,spans,q):
    if q not in TWISTS.values():raise ValueError('actual q=1,i,-1,-i required')
    # Tr(U E0 A0 ... E_k) = Tr(E0 A0 ... E_k U E_0).
    return generalized_renewal(list(letters)+[('U',q)],list(spans)+[0.0])


def zcomplex(z):return [float(complex(z).real),float(complex(z).imag)]


def main():
    original.initialize_native_controls()
    occupancy=np.fromiter((n.bit_count() for n in range(original.D)),dtype=np.int8,count=original.D)
    U={name:sparse.diags(q**occupancy,format='csr') for name,q in TWISTS.items()}
    cases=[
        ('empty',[]),
        ('number2',[(0,True),(0,False)]),
        ('density4',[(0,True),(0,False),(1,True),(1,False)]),
        ('crossed4',[(0,True),(1,False),(1,True),(0,False)]),
        ('density8',[(i,d) for i in range(4) for d in (True,False)]),
        ('omega16',[(i,False) for i in range(16)]),
        ('omega_dagger16',[(i,True) for i in reversed(range(16))]),
        ('density32',[(i,d) for i in range(16) for d in (True,False)]),
    ]
    rows=[];maximum=0.;checks={};count=0
    for beta in (0.0,1.0,2*math.pi):
        a=math.exp(-2*beta);b=-math.expm1(-2*beta)
        for name,letters in cases:
            k=len(letters)
            spans=[beta*(i+1)/sum(range(1,k+2)) for i in range(k+1)]
            matrices=[original.cre[i] if dagger else original.ann[i] for i,dagger in letters]
            for label,q in TWISTS.items():
                direct=original.word_direct(matrices+[U[label]],spans+[0.0])
                fast=twined_word(letters,spans,q)
                error=abs(direct-fast);maximum=max(maximum,error);count+=1
                if error>2e-10*max(1.,abs(direct)):
                    raise AssertionError((beta,name,label,direct,fast,error))
                rows.append({'beta':beta,'word':name,'length':k,'q':label,'native_sparse':zcomplex(direct),'renewal':zcomplex(fast),'absolute_error':error})
                if name=='empty':
                    expected=a*(1+q)**16+b
                    checks[f'Z0_{beta}_{label}']=abs(fast-expected)<1e-11*max(1.,abs(expected))
                if name=='omega16':
                    delta=[spans[0]+spans[-1]]+spans[1:-1]
                    expected=np.exp(1j*original.phi)/2*sum(
                        (-1)**j*q**(16-j)*(-math.expm1(-2*t))*math.exp(-2*(beta-t))
                        for j,t in enumerate(delta))
                    checks[f'whole_det_coherence_{beta}_{label}']=abs(fast-expected)<2e-13
                if name=='number2' and beta==0:
                    checks[f'zero_time_number_character_{label}']=abs(fast-q*(1+q)**15)<1e-11
                if name=='density32' and beta==0 and label=='-1':
                    checks['zero_Z_twined_sector_retains_nonzero_marked_coefficient']=abs(fast-1)<1e-14 and abs(twined_word([],[0.],q))<1e-14
    # Ordering probe with the actual U_i between a creation and annihilation.
    letters=[(0,False),('U',1j),(0,True)]
    direct=original.trace(original.ann[0]@U['i']@original.cre[0])
    checks['internal_twist_uses_intermediate_occupancy']=abs(fock_trace(letters)-direct)<1e-11
    checks['discarding_internal_twist_fails']=abs(fock_trace(letters)-original.word_fock_trace([(0,False),(0,True)]))>1
    # Independent actual center projector, retained only at the original boundary.
    P4=sparse.diags((occupancy%4==0).astype(float),format='csr')
    checks['actual_center_projector_rank']=abs(original.trace(P4)-16512)<1e-12
    for beta in (0.,1.,2*math.pi):
        a=math.exp(-2*beta);b=-math.expm1(-2*beta)
        checks[f'center_selected_Z0_{beta}']=abs(sum(twined_word([],[beta],q) for q in TWISTS.values())/4-(16512*a+b))<1e-10
        for name,letters in cases:
            spans=[beta*(i+1)/sum(range(1,len(letters)+2)) for i in range(len(letters)+1)]
            matrices=[original.cre[i] if dagger else original.ann[i] for i,dagger in letters]
            direct=original.word_direct(matrices+[P4],spans+[0.])
            fast=sum(twined_word(letters,spans,q) for q in TWISTS.values())/4
            error=abs(direct-fast);maximum=max(maximum,error);count+=1
            checks[f'original_boundary_center_projection_{beta}_{name}']=error<2e-10*max(1.,abs(direct))
    # Independent prefix-phase formulas supplied by two separate reviewers.
    for beta in (1.,2*math.pi):
        omega_letters=[(i,False) for i in range(16)]
        uniform=[beta/32]+[beta/16]*15+[beta/32]
        for label,q in TWISTS.items():
            direct=original.word_direct([original.ann[i] for i in range(16)]+[U[label]],uniform+[0.])
            fast=twined_word(omega_letters,uniform,q)
            error=abs(direct-fast);maximum=max(maximum,error);count+=1
            checks[f'uniform_omega16_native_{beta}_{label}']=error<2e-13
            checks[f'uniform_omega16_phase_formula_{beta}_{label}']=abs(fast-np.exp(1j*original.phi)/2*(-math.expm1(-beta/8))*math.exp(-15*beta/8)*sum((-1)**j*q**(16-j) for j in range(16)))<2e-13
        checks[f'uniform_omega16_untwined_zero_{beta}']=abs(twined_word(omega_letters,uniform,1+0j))<2e-13
    s0,t0=.3,.7
    aa,at=math.exp(-2*s0),math.exp(-2*t0)
    bs,bt=-math.expm1(-2*s0),-math.expm1(-2*t0)
    for label,q in TWISTS.items():
        fast=twined_word([(0,True),(0,False)],[s0,t0,0.],q)
        expected=aa*at*q*(1+q)**15+(bs*at+q*aa*bt)/2
        checks[f'independent_twined_two_time_CAR_formula_{label}']=abs(fast-expected)<2e-11
        letters=[(0,True),(0,False),(1,True),(1,False)]
        fast=twined_word(letters,[s0,0.,t0,0.,0.],q)
        expected=aa*at*q*q*(1+q)**14+(bs*at+aa*bt)/2+bs*bt/4
        checks[f'independent_twined_two_time_density_formula_{label}']=abs(fast-expected)<2e-11
    number=sparse.diags(occupancy.astype(float),format='csr')
    number_rows=[]
    for beta in (1.,2*math.pi):
        ss,tt=.3*beta,.7*beta
        aa,at=math.exp(-2*ss),math.exp(-2*tt)
        bb,bt=-math.expm1(-2*ss),-math.expm1(-2*tt)
        all_twined=[]
        for label,q in TWISTS.items():
            chi2=16*q*(1+q)**15+240*q*q*(1+q)**14
            expected=aa*at*chi2+128*(bb*at+aa*bt)+64*bb*bt
            direct=original.word_direct([number,number,U[label]],[ss,tt,0.,0.])
            error=abs(direct-expected);maximum=max(maximum,error);count+=1
            checks[f'gauge_invariant_number_word_{beta}_{label}']=error<2e-10*max(1.,abs(direct))
            all_twined.append(direct)
            number_rows.append({'beta':beta,'s':ss,'t':tt,'q':label,'native':zcomplex(direct),'closed':zcomplex(expected)})
        direct4=original.word_direct([number,number,P4],[ss,tt,0.,0.])
        zz=1+16511*aa*at
        expected4=1115136*aa*at+128*(bb*at+aa*bt)+64*bb*bt
        connected=(64*(aa+at)+58240*aa*at)/zz
        count+=1;maximum=max(maximum,abs(direct4-expected4))
        checks[f'gauge_invariant_number_center_average_{beta}']=abs(direct4-sum(all_twined)/4)<2e-10*max(1.,abs(direct4))
        checks[f'gauge_invariant_number_center_covariance_{beta}']=abs(direct4/zz-64-connected)<2e-11
        number_rows.append({'beta':beta,'s':ss,'t':tt,'q':'center average','unnormalized':zcomplex(direct4),'Z4':zz,'mean_N':8,'connected':connected})
    # Removing U_q is exactly the old evaluator, not a new thermal source.
    for name,letters in cases:
        spans=[1/(len(letters)+1)]*(len(letters)+1)
        checks['untwined_original_'+name]=abs(twined_word(letters,spans,1+0j)-original.word_renewal(letters,spans))<2e-11
    checks={key:bool(value) for key,value in checks.items()}
    assert all(checks.values()),{k:v for k,v in checks.items() if not v}
    report={
        'status':'pass','scope':'Original R104 H=2(I-P_phi), phi=v1002.PHI, all CAR grades, original Spin10 central U_q=q^N. No primitive S_adm or v1013 gauge embedding claimed.',
        'source':str(SOURCE.relative_to(ROOT)),'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        'native_comparisons':count,'max_absolute_error':maximum,'checks':checks,'checks_passed':sum(checks.values()),'checks_total':len(checks),
        'whole_word_theorem':'The unchanged rank-one cyclic renewal is valid for arbitrary inserted bounded operators; exact bit action of U_q supplies every ordered segment and spectator character without inverses.',
        'use_in_existing_group_average':'Replace z_m(g_x) by the computed unnormalized twined ordered coefficient for each original mirror source word, retaining the entire original rest trace and the original bound reference denominator. Applies only where the actual U_m is this supplied canonical representation.',
        'zero_source_formula':'exp(-2 beta)*(1+q)^16+1-exp(-2 beta)',
        'omega16_multitime_formula':'exp(i*phi)/2 * sum_(j=0..15) (-1)^j q^(16-j) (1-exp(-2 delta_j)) exp(-2(beta-delta_j)), delta_0=t_0+t_16 and delta_j=t_j for j>0',
        'center_selected_source':'P4=(1/4)sum_q U_q has rank16512, Z4=16512 exp(-2 beta)+1-exp(-2 beta). Whole marked coefficients are the same group average, P4 remains at the original boundary. This is not the full v1013/Spin10/P1 Gauss map.',
        'gauge_invariant_number_source':'N=sum n_i commutes with every canonical Gamma(r). For original P4, connected C_N(s,t)=[64(exp(-2s)+exp(-2t))+58240exp(-2(s+t))]/[1+16511exp(-2(s+t))].',
        'number_source_values':number_rows,
        'values':rows,
    }
    out=Path(__file__).with_suffix('.json');out.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='values'},indent=2))


if __name__=='__main__':main()
