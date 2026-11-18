"""Exact source-sector test of the existing Ward30 response, not a new model."""
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import sys
import sympy as sp

HERE=Path(__file__).resolve().parent
REPO=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
sys.path.insert(0,str(REPO))
from tfpt_explorer.current_block_geometry import _actual_current_source
ROOT=REPO/'experiments/theory-contracts'
FOLD=ROOT/'v16-operator-pullback-20261014/fold-results.json'
PHASE=ROOT/'elliptic-moebius-source-bridge-20261026/hodge_complement_20261028/e8-clifford-response-results.json'
CURRENT=REPO/'tfpt_explorer/current_block_geometry.py'
CHEVALLEY=REPO/'verification/v498_celestial_wp5b_singular_vector.py'
PINS={FOLD:'253f370a9c84b70c28169ed0aad83d1866245804227e44aa0189e823382f5639',
      PHASE:'67fe9e3f06140c5c23452e9a7bfaefc02441ff1e7b55708e55409cd0a398e7ee',
      CURRENT:'c9efcac4ba33d4ee3dc20f8c3c69a0f6ff6b2084893a5edf4aaae6c204d9056b',
      CHEVALLEY:'9f5a2d62523283595b754fea6f70ddc6501a251daf86752cbc67ded252b96180'}
checks=[]
def need(v,label):
    if not v:raise RuntimeError(label)
    checks.append(label)
def key(r):return ','.join(map(str,r))
def ladder(mask,i,create):
    if ((mask>>i)&1)==create:return None
    return mask^(1<<i),(-1)**((mask&((1<<i)-1)).bit_count())
def hist(xs):return {str(k):v for k,v in sorted(Counter(xs).items())}
for p,h in PINS.items():need(hashlib.sha256(p.read_bytes()).hexdigest()==h,'source pin '+p.name)
M=[[F(x) for x in row] for row in json.loads(FOLD.read_text())['explicit_Weyl_matrix']]
def pull(r):return tuple(sum(M[j][i]*F(r[j],2) for j in range(8)) for i in range(8))
def grade(r):
    p=pull(r)
    if all(z.denominator==1 for z in p):return 0
    need(all(z.denominator==2 for z in p),'source spinor has uniform half integral charge')
    return 1
source=_actual_current_source();alg=source['chevalley'];roots=tuple(sorted(alg.roots));R=set(roots)
half=[r for r in roots if all(abs(z)==1 for z in r)]
masks={r:sum((z>0)<<i for i,z in enumerate(r)) for r in half}
hin=[r for r in half if (masks[r]&31).bit_count()==2]
hout=[r for r in half if (masks[r]&31).bit_count()==3]
creator=[r for r in roots if sum(z!=0 for z in r[:5])==1 and sum(z!=0 for z in r[5:])==1 and sum(r[:5])==2]
ns=[r for r in creator if grade(r)==0];spin=[r for r in creator if grade(r)==1]
need(len(hin)==len(hout)==40 and len(creator)==30,'actual typed Ward source dimensions')
need(Counter(grade(r) for r in hin)==Counter({0:20,1:20}),'input40 raw NS20 plus spinor20')
need(Counter(grade(r) for r in hout)==Counter({0:20,1:20}),'output40 raw NS20 plus spinor20')
need((len(ns),len(spin))==(14,16),'creation30 raw NS14 plus spinor16')
w=tuple(row[0] for row in M)
need(w==(F(1,2),F(1,2),0,0,0,F(1,2),F(1,2),0),'raw parity cocharacter in the target chart')
for r in hin+hout:
    mask=masks[r];pc=((mask&3).bit_count())%2;pf=(((mask>>5)&3).bit_count())%2
    need(grade(r)==(pc+pf)%2,'source grading equals carrier parity times family parity')
need(hist((mask&3).bit_count()%2 for mask in range(32) if mask.bit_count()==2)=={'0':4,'1':6},'carrier grading has dimensions4 plus6')

ph=json.loads(PHASE.read_text());gains={r:Counter() for r in hin};trans=Counter();columns=[]
for b in creator:
    m=next(i for i in range(5) if b[i]);j=next(i for i in range(5,8) if b[i])-5
    for a in hin:
        raw=alg.bracket(alg.ridx[b],alg.ridx[a]);ca=ladder(masks[a]&31,m,True);fa=ladder(masks[a]>>5,j,b[j+5]>0)
        if ca is None or fa is None:
            need(not raw,'actual bracket zeros equal tensor zeros');continue
        g=tuple(x+y for x,y in zip(a,b));need(g in hout,'actual bracket has target output')
        need(set(raw)=={alg.ridx[g]},'one signed output per actual root column')
        c=raw[alg.ridx[g]];phase=int(ph['state_phases'][key(a)])*int(ph['state_phases'][key(g)])*int(ph['root_phases'][key(b)])
        need(c==phase*ca[1]*fa[1],'actual signed bracket equals pinned joint CAR tensor action')
        need(grade(g)==(grade(a)+grade(b))%2,'each Kraus has exact source sector covariance')
        gains[a][grade(b)]+=1;trans[(grade(a),grade(g),grade(b))]+=1;columns.append((a,b,g,c))
need(len(columns)==360,'full30 has360 signed nonzero columns')
need(trans==Counter({(0,0,0):84,(1,1,0):84,(0,1,1):96,(1,0,1):96}),'complete source sector transition census')
gain_hist=Counter((grade(a),gains[a][0],gains[a][1]) for a in hin)
need(gain_hist==Counter({(q,5,4):6 for q in (0,1)})+Counter({(q,4,5):12 for q in (0,1)})+Counter({(q,3,6):2 for q in (0,1)}),'uniform14 restriction has nonconstant effects3,4,5')

# Exact classification of diagonal positive weights on the same NS14 roots.
# This is a comparison family, not a claim that P1 chooses those weights.
A=sp.Matrix([[int(tuple(x+y for x,y in zip(a,b)) in R) for b in ns] for a in hin])
need(A.rank()==12,'NS14 trace-preservation equations have exact rank12')
weights=[]
for b in ns:
    m=next(i for i in range(5) if b[i]);weights.append(F(1,6) if m<2 else F(1,3))
need(A*sp.Matrix([sp.Rational(z.numerator,z.denominator) for z in weights])==sp.ones(40,1),'explicit positive NS14 weights give exact TP')
for m in range(5):
    for fm in (0,3,5,6):
        effect=F(0)
        for b,weight in zip(ns,weights):
            if b[m]!=2:continue
            j=next(i for i in range(5,8) if b[i])-5
            if ladder(fm,j,b[j+5]>0) is not None:effect+=weight
        need(effect==F(1,3),'operator-valued family gain is I4/3 for every carrier Kraus')
need(sum(trans[k] for k in trans if k[2]==0)==168,'NS14 witness uses168 actual signed columns')

# The sector label is part of the joint object and can distinguish responses
# whose full carrier marginal is the same Hodge channel.
beta=(0,2,0,0,0,2,0,0);color=(0,-2,2,0,0,0,0,0);gamma=tuple(x+y for x,y in zip(beta,color))
need(beta in ns and gamma in spin and bool(alg.bracket(alg.ridx[color],alg.ridx[beta])),'physical color action takes NS14 support outside itself')
need(grade(color)==1,'this physical color current is spinorial in the raw source frame')
allowed_marginal=0
for a in hin:
    if grade(a):continue
    for b in hin:
        if grade(b) or masks[a]>>5!=masks[b]>>5:continue
        need(((masks[a]&3).bit_count()-(masks[b]&3).bit_count())%2==0,'NS-only inputs have carrier block diagonal reduced states')
allowed_marginal=4*4+6*6
need(allowed_marginal==52,'NS-only carrier observable algebra has52 not100 directions')

result={'verdict':'PASS_SCOPED_SOURCE_SECTOR_RESPONSE_TEST',
 'source_sector_dimensions':{'input40':[20,20],'output40':[20,20],'creation30':[14,16]},
 'sector_transition_columns':{','.join(map(str,k)):v for k,v in sorted(trans.items())},
 'uniform_NS14_input_effect_histogram':{','.join(map(str,k)):v for k,v in sorted(gain_hist.items())},
 'actual_signed_columns':360,
 'raw_source_grading_target_cocharacter':list(map(str,w)),
 'grading_factorization':'(-1)^(N_C,first2) tensor (-1)^(N_F,first2)',
 'carrier_parity_dimensions':[4,6],'NS_only_carrier_matrix_directions':allowed_marginal,
 'NS14_alternative':{'effect_matrix_rank':12,'affine_weight_parameters':2,
   'positive_TP_witness':'1/6 on eight roots with carrier m=1,2; 1/3 on six roots with carrier m=3,4,5',
   'carrier_marginal':'exact full Hodge Phi on the stipulated joint40 input, by operator-valued family gain I4/3',
   'preserves_raw_source_sector':True,'full_physical_P2_color_covariance':False,
   'P1_selection_derived':False},
 'scope':'Actual given E8 current response and same pinned Fold/CAR chart. Uniform30 and normalized NS14 have identical full carrier marginal but different joint source-sector transport. NS14 witness is not the unchanged30 process, not a P1-derived weighting, and not an embedding of all40 inputs into the raw NS vacuum Hilbert space. NS-only carrier reduction has a4+6 grading. No general no-go for DHR observable maps or other P1 sources.',
 'sources':{str(p.relative_to(REPO)):h for p,h in PINS.items()},'checks':len(checks)}
(HERE/'source-sector-response-results.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
