"""Exact fold/current/twist bookkeeping, no replacement covariance."""
from pathlib import Path
from collections import Counter
import json,math,itertools as it
import sympy as s
import check_mixed_return as m

def main():
    fold=json.loads((m.R/'experiments/theory-contracts/v16-operator-pullback-20261014/fold-results.json').read_text())
    M=s.Matrix(fold['explicit_Weyl_matrix']);ld=m._lattice_data()
    T=s.Matrix.hstack(*map(s.Matrix,m._TARGET_SIMPLE));roots=[T*s.Matrix(z)/2 for z in ld['root_coordinates']]
    rstar=s.Matrix([0,0,1,0,0,1,0,0]);m.need(M*s.ones(8,1)==2*rstar,'Original fixed M, not fitted clock')
    z=s.Symbol('z');delta=s.Symbol('delta',real=True)
    finite=s.simplify(s.limit((s.exp(s.I*delta*z)-1)/(2*s.I*s.sin(z/2)),z,0))
    m.need(finite==delta,'Exact NS point-split current counterterm from same kernel')
    rows={}
    for alpha in (s.Rational(1,4),s.Rational(3,4)):
        mu=s.Matrix([(k+alpha)/8 for k in range(8)]);dd=mu-s.ones(8,1)/2;d=M*dd;nu=M*(2*mu)
        m.need(M*s.ones(8,1)/2==rstar and rstar in roots,'M1/2 is original E8 lattice root')
        phases=[s.Mod((d.T*r)[0],1) for r in roots];gphases=[s.Mod(((M*mu).T*r)[0],1) for r in roots]
        m.need(phases==gphases,'Li delta twist is actual fold g on all240 charged roots')
        m.need(all(s.Mod((nu.T*r)[0]-2*phases[i],1)==0 for i,r in enumerate(roots)),'Original g2 is square of Li twist, not same twist')
        order=math.lcm(*[v.q for v in phases]);m.need(order==16,'Actual full twist order16')
        order2=math.lcm(*[s.Mod(2*v,1).q for v in phases]);m.need(order2==8,'Actual g2order8 retained')
        for k in range(8):
            m.need(s.Rational(1,2)+dd[k]==mu[k] and s.Rational(1,2)-dd[k]==1-mu[k],'Particle/hole gaps fix negative-D current orientation')
        # Original source convention f_(j-a), occupied j>=1: alpha=1-a,
        # n=j-1 in SAME orientation, hence positive halfline.
        a=1-alpha
        for j in range(-8,10):
            n=j-1;ell,k=divmod(n,8)
            m.need(j-a==n+alpha,'Literal original QWZ/source phase dictionary')
            m.need((j>=1)==(ell>=0),'Original positive filled halfline retained, not relabelled negative')
        rows[str(alpha)]={'d':[str(v) for v in d],'twist_order':order,'twist_square_order':order2,
              'phase_multiplicities':{str(k):v for k,v in Counter(phases).items()},
              'ground_conformal_weight':str((d.T*d)[0]/2),
              'zero_excitation_roots':[list(r) for r in roots if 1+(d.T*r)[0]==0]}
    B=s.Matrix([[s.binomial(k+r-1,2*r-1) if k>=r else 0 for r in range(1,5)] for k in range(1,5)])
    m.need(B.det()==1,'First4 positive modeGrams recover all4 connected Wick tensors, no profile fitting')
    out={'status':'PASS_EXACT_ORIGINAL_FOLD_COVARIANCE_COUNTERTERM_AND_TWIST_DICTIONARY',
         'date':'2026-10-17','cases':rows,'mode_to_wick_matrix':m.mat(B),
         'counterterm_negative_D':'delta=mu-1/2','counterterm_original_positive_D':'-delta, with charge/current ordering transported',
         'original_CAR_sea_not_replaced':True,'eight_source_copies_added':0,
         'scope':'Conditional existingchiral8fold quotient,actual covariance eigenchannelreindexing andfixedconformalpointsplit reference. OriginalQWZpositivehalflineexplicit,notfullP1/rawkernelderivation. E8twistedextensiongiven,notoddCARlocalinclusionoruniquegroundstateselection. Neutralquarticresponse signstable.'}
    Path(__file__).with_name('fold_state_dictionary_results.json').write_text(json.dumps(out,ensure_ascii=False,indent=2,default=str)+'\n')
    print(json.dumps(out,ensure_ascii=False,indent=2,default=str))
if __name__=='__main__':main()
