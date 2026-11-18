"""All-position shifted quartic Wick kernel from original fold delta.
No mode cutoff, selected q/flag, optimization or physical instrument.
"""
import json,itertools as it,math
from pathlib import Path
from collections import Counter
import sympy as s
import check_mixed_return as m
from check_fold_quartic_transport import cl,zero

def main():
    ld=m._lattice_data();T=s.Matrix.hstack(*map(s.Matrix,m._TARGET_SIMPLE));O=T*ld['basis'].inv()/s.sqrt(2)
    U=s.zeros(8,4)
    for j in range(4):U[2*j,j]=1/s.sqrt(2);U[2*j+1,j]=-s.I/s.sqrt(2)
    M=s.Matrix(json.loads((m.R/'experiments/theory-contracts/v16-operator-pullback-20261014/fold-results.json').read_text())['explicit_Weyl_matrix'])
    ps,ns=m.code();x=s.symbols('x0:4')
    polys=[sum(m.Q.to_sympy(c)*s.prod(x[j] for j in mon) for mon,c in f.items())/s.sqrt(ns[i]) for i,f in enumerate(ps)]
    raw=s.Matrix([[2,2,-1,-1,-1,-1],[0,0,-1,-1,1,1],[0,0,-1,1,-1,1],[0,0,1,-1,-1,1],[1,-1,0,0,0,0]])
    N=cl(s.diag(2,s.sqrt(12),s.sqrt(12),s.sqrt(12),s.sqrt(24))*raw/s.sqrt(48));rows={}
    for alpha in (s.Rational(1,4),s.Rational(3,4)):
        delta=s.Matrix([(k+alpha)/8-s.Rational(1,2) for k in range(8)])
        b=cl(U.T*O.T*M*delta);sub=dict(zip(x,b));Ds=[]
        for r in range(5):
            monos=list(it.combinations_with_replacement(range(4),r))
            E=s.Matrix([[s.simplify((s.diff(f,*[x[j] for j in mon]) if mon else f).subs(sub,simultaneous=True))/s.sqrt(math.prod(math.factorial(n) for n in Counter(mon).values())) for f in polys] for mon in monos])
            Ds.append(cl(E.H*E))
        m.need(Ds[4]==s.eye(5),'Universal r4 tensor retains original quartic normalization')
        m.need(zero(Ds[3]-(b.H*b)[0]*s.eye(5)),'Original erasure20 forces scalar r3 response')
        ker=Ds[1].nullspace();m.need(len(ker)==1,'Original first new response has rank4, unique kernel direction')
        k=cl(ker[0]);matches=[q for q in range(6) if s.Matrix.hstack(k,N[:,q]).rank()==1]
        means=s.Matrix([s.simplify(f.subs(sub,simultaneous=True)) for f in polys])
        meanmatches=[q for q in range(6) if s.Matrix.hstack(means,N[:,q]).rank()==1]
        # Compare entire state-dependent response tensors with every actual flag.
        commutes=[];r1commutes=[]
        # Carry raw code metric to Gaussian rationals for cheap exact zero tests.
        scale=s.diag(*[s.sqrt(n) for n in ns]);rawDs=[cl(scale.inv()*D*scale) for D in Ds]
        for q in range(6):
            for pair in it.combinations([j for j in range(6) if j!=q],2):
                vv=[N[:,j]-(s.sqrt(6)-1)*N[:,q]/5 for j in pair]
                P=cl(sum([v*v.T for v in vv],s.zeros(5)));Pr=cl(scale.inv()*P*scale)
                z1=zero(rawDs[1]*Pr-Pr*rawDs[1])
                if z1:r1commutes.append([q,list(pair)])
                if z1 and all(zero(D*Pr-Pr*D) for D in rawDs):commutes.append([q,list(pair)])
        # Full local kernel fixes every positive-frequency mode Gram.
        B4=cl(4*Ds[1]+10*Ds[2]+6*Ds[3]+Ds[4])
        full=json.loads(Path(__file__).with_name('fold_quartic_transport_results.json').read_text())['cases'][str(alpha)]
        m.need(zero(B4-s.Matrix(full['full_level4_response_gram'])),'Independent all-position Wick kernel reproduces entire direct Fock mode4 Gram')
        # b->-b flips D_r amplitudes by (-1)^(4-r), hence leaves all D_r.
        # conjugation transports all tensors and real candidate flags together.
        rows[str(alpha)]={'positive_shift':[str(v) for v in b],'onepoint':[str(v) for v in means],
          'onepoint_simplex_matches':meanmatches,'wick_tensors':[m.mat(D) for D in Ds],
          'r1_rank':Ds[1].rank(),'r1_kernel':[str(v) for v in k],
          'r1_kernel_simplex_matches':matches,'r1_commuting_flags':r1commutes,
          'joint_commuting_flags':commutes,'mode4_from_wick_tensors':m.mat(B4),
          'scalar_r3':str(s.simplify((b.H*b)[0]))}
        print(json.dumps({'alpha':str(alpha),'mean_matches':meanmatches,'kernel_matches':matches,'r1_commutes':r1commutes,'joint_commutes':commutes}),flush=True)
    out={'status':'PASS_EXACT_ALL_POSITION_FOLD_QUARTIC_KERNEL','date':'2026-10-17','cases':rows,
         'all_mode_formula':'B_k=sum_r1to4 binomial(k+r-1,2r-1) D_r for k>=1, binomial zero for k<r; B0=D0',
         'kernel':'sum_r0to4 D_r / ((z-w)^(2r) z^(4-r) w^(4-r))',
         'scope':'Givenconformalpointsplit/currentdictionary,singleexistingnegative-Dchiralorientationfoldsea,existingM/O. OriginalQWZpositivehalflinerequirescurrent/particle-holeorientationtransport;notautomaticstate-adaptedcurrentidentification,nophysicalmodefilter/profile/instrument,newcoupling,allP1selectionno-go. Bothquarterphases;orientation/chargeinversion b->-b keepsDtensors.'}
    Path(__file__).with_name('fold_response_tensors_results.json').write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
if __name__=='__main__':main()
