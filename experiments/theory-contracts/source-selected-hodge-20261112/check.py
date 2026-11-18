from pathlib import Path
import contextlib, io, runpy, json
import numpy as np
import sympy as s
ROOT=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
BASE=Path(__file__).with_name('common_construct.py')
with contextlib.redirect_stdout(io.StringIO()):
    b=runpy.run_path(str(BASE))
C,C5,C6=b['C'],b['C5'],b['C6']
Ps,Pa=b['Psym'],b['Panti']
PH=C5/2
F=s.zeros(100)
for a in range(10):
    for c in range(10):F[a+10*c,c+10*a]=1
P45=(s.eye(100)-F)/2
P50=(s.eye(100)+F)/2-PH
assert PH*PH==PH and s.trace(PH)==5
assert PH*Ps==PH and PH*Pa==s.zeros(100)
assert PH*C==C*PH==s.Rational(2,15)*PH
branch=PH*C*PH
assert branch==C5/15
rem=C-branch
assert rem==s.Rational(2,15)*(Ps-PH)+s.Rational(8,15)*Pa
partial=lambda M:s.Matrix(10,10,lambda a,c:sum(M[j+10*a,j+10*c] for j in range(10)))
assert partial(branch)==s.eye(10)/15
assert partial(rem)==s.Rational(14,15)*s.eye(10)
for P,r in [(PH,5),(P45,45),(P50,50)]:
    assert P*P==P and s.trace(P)==r and partial(P)==s.Rational(r,10)*s.eye(10)
weights=lambda M:[s.trace(P*M)/10 for P in (PH,P45,P50)]
assert weights(C)==[s.Rational(1,15),s.Rational(4,5),s.Rational(2,15)]
assert weights(C6)==[s.Rational(1,3),0,s.Rational(2,3)]
# All original twenty immediate jumps and all stay ranks, in actual D chart.
D=b['D']; ranksets={'stay':set(),'jump':set()};jump_skew=0
for k in b['selected']:
    U=b['bell'][k];mask=np.diag(np.diag(U)!=0).astype(int)
    for label,K in [('stay',U@mask),('jump',U@(np.eye(10)-mask))]:
        M=s.Matrix(D@K@D.conj().T).applyfunc(s.nsimplify)
        ranksets[label].add(M.rank())
        if label=='jump':
            assert M.T==-M
            assert PH*b['vec'](M)==s.zeros(100,1)
            jump_skew+=1
assert ranksets=={'stay':{4},'jump':{6}}
# Nonzero epsilon contraction has rank6 for every complex vector, by GL5 orbit;
# exact basis and a complex isotropic vector are independent diagnostic examples.
Cs=[s.Matrix(H) for H in b['Cs']]
sample_ranks=[]
for v in ([1,0,0,0,0],[1,s.I,0,0,0],[1,2,s.I,3,-2]):
    M=sum((v[i]*Cs[i] for i in range(5)),s.zeros(10))
    sample_ranks.append(M.rank())
assert sample_ranks==[6,6,6]
# A full carrier SU5-covariant submap is a positive linear combination of
# PH,P45,P50. Forbidden Petersen diagonal entries force coefficients45,50=0.
for d,P in [(5,PH),(45,P45),(50,P50)]:
    J=s.Rational(10,d)*P
    # source input duad01, same01/intersect02/disjoint23
    duads=b['duads'];a=duads.index((0,1))
    obs=[J[z+10*a,z+10*a] for z in [a,duads.index((0,2)),duads.index((2,3))]]
    expected={5:[0,0,s.Rational(1,3)],45:[0,s.Rational(1,9),s.Rational(1,9)],50:[s.Rational(1,5),s.Rational(1,10),s.Rational(1,15)]}[d]
    assert obs==expected
# Resolve the original stopped histories by exact number n of genuine stays.
S20=s.Matrix(b['S']).applyfunc(s.nsimplify)/20
J20=s.Matrix(b['J']).applyfunc(s.nsimplify)/20
for a in range(10):
    assert S20[a+10*a,a+10*a]==s.Rational(2,5)
for a,c in b['edges']:
    assert S20[a+10*c,a+10*c]==0
    assert S20[c+10*a,c+10*a]==0
for count in (1,2,3):
    Jn=b['W']*b['choi'](J20*S20**count)*b['W'].conjugate().T
    expected=s.Rational(1,5)*s.Rational(2,5)**count*(Ps+Pa)
    assert Jn==expected
    assert partial(Jn)==s.Rational(3,5)*s.Rational(2,5)**count*s.eye(10)
    assert PH*Jn*PH==s.Rational(1,10)*s.Rational(2,5)**count*C5
# Twirl is a comparison, never substituted for the original channel.
result={'verdict':'established_scoped; no physical P1 selection',
 'chart':'actual common G240/A5 D chart; dual-SU5 extension declared',
 'maximal_dual_su5_covariant_subchannel':'Phi5/15',
 'canonical_Hodge_branch_probability':'1/15',
 'remainder_probability':'14/15',
 'classical_history_postprocessing_maximum':'0',
 'first_activation':{'zero_stays_probability':'0','n_stays_joint_success':'(1/10)*(2/5)^n for n>=1','success_given_positive_stay_count':'1/6','success_total':'(2/5)*(1/6)=1/15','conditional_mean_stays':'5/3','scope':'source-event count, not physical time; coherence only among histories of the same count suffices'},
 'stay_ranks':sorted(ranksets['stay']),'jump_ranks':sorted(ranksets['jump']),
 'skew_immediate_jumps':jump_skew,
 'all_nonzero_epsilon_matrices_rank':6,
 'rank_proof':'GL5 equivariance of wedge-by-v, rank choose(4,2)=6; examples diagnostic only',
 'twirl_weights_5_45_50':{'native_stopped':list(map(str,weights(C))),'six_matching':list(map(str,weights(C6)))},
 'physical_status':'Coherent spectral environment projection is a selected operation. Original physical source equality, event field map, state/time and weight identification not supplied by this finite theorem.'}
Path(__file__).with_name('selection_results.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
