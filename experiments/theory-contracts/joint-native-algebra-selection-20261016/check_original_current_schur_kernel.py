"""Actual original CAR packet -> Hall/Schur representation, complete current kernel.

Same source states and maps, no new clock, coupling, mode-pair selection or
identification of folded local fields with point-local cover fields.
"""
import json,math
from functools import lru_cache
from pathlib import Path
import sympy as s
import check_mixed_return as m
import check_original_density_current as car
from check_fold_quartic_transport import cl,series_response
Q=m.Q

@lru_cache(None)
def maya_partition(bits):
    occupied=[r for r in range(car.HIGH,car.LOW-1,-1) if (bits>>(r-car.LOW))&1]
    lam=tuple(r+i for i,r in enumerate(occupied,1))
    m.need(all(x>=0 for x in lam) and all(a>=b for a,b in zip(lam,lam[1:])), 'Actual charge-zero Maya partition')
    return tuple(x for x in lam if x)

def to_schur(v):
    out={}
    for bits,c in v.items():
        lam=maya_partition(bits)
        m.need(sum(lam)==car.energy(bits),'Original source energy is Hall weighted degree')
        out[lam]=out.get(lam,Q.zero)+c
    return car.clean(out)

@lru_cache(None)
def ribbons(lam):
    """All removable Young rim ribbons, checked by partition/skew shape."""
    if not lam:return {}
    row,col=0,lam[0];rim=[]
    while col>0:
        rim.append((row,col))
        if row+1<len(lam) and lam[row+1]>=col:row+=1
        else:col-=1
    out={}
    for start in range(len(rim)):
        removed={}
        for stop in range(start,len(rim)):
            r,c=rim[stop];removed.setdefault(r,set()).add(c)
            mu=tuple(n-len(removed.get(i,())) for i,n in enumerate(lam))
            if any(a<b for a,b in zip(mu,mu[1:])):continue
            if any(cells!=set(range(mu[i]+1,lam[i]+1)) for i,cells in removed.items()):continue
            target=tuple(x for x in mu if x);degree=stop-start+1
            sign=(-1)**(len(removed)-1)
            out.setdefault(degree,[]).append((target,sign))
    return out

def lower(v,mode):
    out={}
    for lam,c in v.items():
        for mu,sign in ribbons(lam).get(mode,[]):out[mu]=out.get(mu,Q.zero)+sign*c
    return car.clean(out)

@lru_cache(None)
def add_ribbons(lam,mode):
    """Power-sum multiplication via infinite Maya beads, no finite-window cutoff."""
    count=len(lam)+mode+2
    occupied=[(lam[i-1] if i<=len(lam) else 0)-i for i in range(1,count+1)]
    occ=set(occupied);out=[]
    for old in occupied:
        new=old+mode
        if new in occ:continue
        sign=(-1)**sum(old<x<new for x in occupied)
        target=sorted(occ-{old}|{new},reverse=True)
        mu=tuple(v+i for i,v in enumerate(target,1))
        m.need(all(x>=0 for x in mu),'Infinite Maya tail is fully retained')
        mu=tuple(x for x in mu if x)
        m.need(sum(mu)==sum(lam)+mode,'Power sum raises exactly source degree')
        out.append((mu,sign))
    return out

def raise_mode(v,mode):
    out={}
    for lam,c in v.items():
        for mu,sign in add_ribbons(lam,mode):out[mu]=out.get(mu,Q.zero)+sign*c
    return car.clean(out)

def stress(v,mode,alpha):
    out={};degree=max((sum(lam) for lam in v),default=0);delta=Q.from_sympy(alpha-s.Rational(1,2))
    if mode>0:
        for n in range(1,max(1,degree-mode+1)):car.add(out,raise_mode(lower(v,mode+n),n),Q.one)
        for n in range(1,mode):car.add(out,lower(lower(v,n),mode-n),Q.from_sympy(s.Rational(1,2)))
        car.add(out,lower(v,mode),delta)
    elif mode<0:
        k=-mode
        for n in range(1,degree+1):car.add(out,raise_mode(lower(v,n),n+k),Q.one)
        for n in range(1,k):car.add(out,raise_mode(raise_mode(v,n),k-n),Q.from_sympy(s.Rational(1,2)))
        car.add(out,raise_mode(v,k),delta)
    else:raise ValueError('Use weighted degree for L0 including ground term')
    return car.clean(out)

def source_packets():
    ld=m._lattice_data();T=s.Matrix.hstack(*map(s.Matrix,m._TARGET_SIMPLE))
    O=cl(T*ld['basis'].inv()/s.sqrt(2));U=s.zeros(8,4)
    for a in range(4):U[2*a,a]=1/s.sqrt(2);U[2*a+1,a]=-s.I/s.sqrt(2)
    M=s.Matrix(json.loads((m.R/'experiments/theory-contracts/v16-operator-pullback-20261014/fold-results.json').read_text())['explicit_Weyl_matrix'])
    V=cl(M.T*O*U);vq=[[Q.from_sympy(V[k,a]) for k in range(8)] for a in range(4)]
    ps,norms=m.code()
    @lru_cache(None)
    def word_state(word):
        if not word:return {car.VAC:Q.one}
        n,a=divmod(word[-1],8);out={}
        for k,c in enumerate(vq[a]):
            if c:car.add(out,car.density(word_state(word[:-1]),8*(n+1),k),c)
        return car.clean(out)
    def polynomial_state(poly):
        out={}
        for word,c in poly.items():car.add(out,word_state(word),c)
        return car.clean(out)
    result={}
    for alpha in (s.Rational(1,4),s.Rational(3,4)):
        delta=s.Matrix([(k+alpha)/8-s.Rational(1,2) for k in range(8)]);b=cl(V.T*delta)
        states=[polynomial_state(series_response(f,[Q.from_sympy(x) for x in b],4)) for f in ps]
        result[alpha]=states
    return result,norms

def main():
    packets,norms=source_packets();old=json.loads(Path(__file__).with_name('original_density_current_results.json').read_text());rows={}
    for alpha,states in packets.items():
        ps=[to_schur(v) for v in states];scale=s.diag(*[s.sqrt(n) for n in norms])
        gram=s.Matrix(5,5,lambda i,j:Q.to_sympy(car.inner(ps[i],ps[j])))
        m.need(cl(scale.inv()*gram*scale.inv())==s.Matrix(old['cases'][str(alpha)]['full_input_gram']), 'Hall orthonormal Schur metric reproduces actual full state Gram')
        m.need(all(car.charges(w)==(0,)*8 for v in states for w in v),'Actual native packet has empty8core, not arbitrary source energy32 states')
        for mode in range(1,33):
            for v,p in zip(states,ps):
                m.need(lower(p,mode)==to_schur(car.density(v,-mode)), 'Independent ribbon formula equals full CAR current vector for every allowed lowering mode')
        print(json.dumps({'alpha':str(alpha),'all32_CAR_ribbon_current_maps':'PASS'}),flush=True)
        for mode in (1,2):
            for v,p in zip(states,ps):
                original=car.density(v,-mode,weights=lambda r:Q.from_sympy(r-s.Rational(mode,2)+alpha))
                m.need(stress(p,mode,alpha)==to_schur(original),'Original geometric stress is the fixed Hall Sugawara operator, full vectors')
        for p in ps:
            lhs=lower(raise_mode(p,1),1);car.add(lhs,raise_mode(lower(p,1),1),-Q.one);car.add(lhs,p,-Q.one)
            m.need(not car.clean(lhs),'Full Hall creation/annihilation CCR, no finite-sector closure')
        # Entire normal-ordered first-mode generating matrix, no perturbative cutoff.
        lowered=ps;coefficients=[]
        for k in range(33):
            H=s.Matrix(5,5,lambda i,j:Q.to_sympy(car.inner(lowered[i],lowered[j])))
            coefficients.append({'power':k,'raw_gram':m.mat(H),'factorial_denominator':str(math.factorial(k)**2)})
            lowered=[lower(v,1) for v in lowered]
        m.need(all(not v for v in lowered),'Normal ordered first-mode jet terminates by source energy, not chosen truncation')
        rows[str(alpha)]={'weighted_degree':32,'norms':norms,'schur_support_sizes':[len(v) for v in ps],
            'all32_current_vector_maps_match_CAR':True,'original_stress_modes_1_2_match_CAR':True,
            'normal_ordered_first_mode_jet':coefficients,
            'schur_coefficients':[[{'partition':list(lam),'coefficient':str(Q.to_sympy(c))} for lam,c in sorted(v.items())] for v in ps]}
        print(json.dumps({'alpha':str(alpha),'full_Hall_kernel_checks':'PASS','support_sizes':[len(v) for v in ps]}),flush=True)
    result={'status':'PASS_EXACT_ORIGINAL_CURRENT_HALL_STATE_AND_ORDERED_KERNEL','date':'2026-10-17','cases':rows,
        'Hall_dictionary':'Actual source Maya positions r_i=lambda_i-i; Schur states orthonormal; J_-m=p_m, J_m=m partial_p_m',
        'all_profile_normal_kernel':'K_ij(u,v)=<P_i(p),exp(sum_m u_m p_m) P_j(p_m+m v_m)>Hall; finite weighted u/v degrees <=32, with every state coefficient retained',
        'Weyl_kernel':'exp(-sum_m>0 m|phi_m|^2/2) K_ij(i phi_m,i conjugate(phi_m)); Gaussian retains ALL modes; finite polynomial factor uses only m<=32',
        'scope':'Original single-circle Q0 CAR state/current representation, fixed prior NS/native dictionary, no new physical source/ref/clock/coupling. Hall STATE/current intertwiner is NOT a point-local field map: native folded quartics remain supported on eight original preimage arcs, not reidentified with source-local weight32 fields. All smooth current Weyl histories have original BCH phases; finite stress-mode words follow Sugawara. Arbitrary exponentiated geometry or unsmeared coincident products need separate analytic controls. No all-E8/P1 or physical3+2 claim.'}
    Path(__file__).with_name('original_current_hall_kernel_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'])

if __name__=='__main__':main()
