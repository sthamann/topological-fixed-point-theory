"""Actual compact positive two-current Gram, compared to native quartic projector.
No arbitrary projection: build all symmetric two-current states in E8 level1.
"""
import runpy,json,itertools
from contextlib import redirect_stdout
from io import StringIO
from collections import defaultdict
from fractions import Fraction as F
import numpy as np
from scipy.sparse import csr_matrix,diags
with redirect_stdout(StringIO()):s=runpy.run_path(str(__import__('pathlib').Path(__file__).with_name('e7_quartic_projector.py')))
ch,af,G,K,pairs,ck=(s[k] for k in ['ch','env','G','K','pairs','ck']);af=af['af'];vac={():F(1)}
weights=defaultdict(list);states=[]
for p,(i,j) in enumerate(pairs):
    weights[tuple(a+b for a,b in zip(ch.roots[G[i]],ch.roots[G[j]]))].append(p)
    state={}
    for a,b in [(i,j),(j,i)]:
        for word,c in af.act(G[a],-1,af.act(G[b],-1,vac)).items():state[word]=state.get(word,F(0))+c/2
    states.append({w:c for w,c in state.items() if c})
rows=[];cols=[];vals=[];count=0
mult=[1 if i==j else 2 for i,j in pairs]
for block in weights.values():
    for p,q in itertools.product(block,repeat=2):
        i,j=pairs[p];inner=F(0)
        for a,b in [(i,j),(j,i)]:
            st=af.act(ch.opp[G[a]],1,states[q]);st=af.act(ch.opp[G[b]],1,st)
            inner+=st.get((),F(0))/2
        inner*=mult[p]*mult[q]
        ck(inner.denominator==1,'actual symmetric pair compact Gram integral')
        if inner:rows.append(p);cols.append(q);vals.append(int(inner))
        count+=1
Gram=csr_matrix((vals,(rows,cols)),shape=K.shape,dtype=np.int64)
W=diags(mult,dtype=np.int64,format='csr');Pnum=K+2*s['I']
target=W@Pnum
rr,cc=Gram.nonzero();ratios={F(int(Gram[i,j]),int(target[i,j])) for i,j in zip(rr,cc) if target[i,j]}
ck(len(ratios)==1,'one source-derived compact Gram proportionality')
factor=next(iter(ratios));delta=Gram*factor.denominator-target*factor.numerator;delta.eliminate_zeros()
ck(delta.nnz==0,'whole1596 source Gram equals compact metric times quartic projector')
sym=K.T@W-W@K;sym.eliminate_zeros();ck(sym.nnz==0,'quartic channel operator self-adjoint in inherited compact Sym2 metric')
result={'claim':'MATH.ESERIES.ACTUAL.TWO.CURRENT.GRAM.20261106','state_definition':'S_ab=(J_a(-1)J_b(-1)+J_b(-1)J_a(-1))Omega/2; coordinates use2*S_ab for a<b and S_aa on diagonal','inherited_Sym2_metric':'W diagonal1 for a=b and2 for a<b','actual_Gram_identity':f'Gram={factor}*W*(K+2I)={factor*24}*W*P133','compared_weight_allowed_entries':count,'Gram_nonzero_entries':Gram.nnz,'Gram_rank':133,'Gram_nullity':1463,'scope':'Actual positive compact source Gram in original level1 E8 vacuum. Null1463 subspace vanishes in the unitary integrable quotient;133 survives. This is a derived internal two-current state channel, not the primitive physical admissibility projector or physical fermion spectrum.','checks':dict(s['env']['checks'])}

# Formal level variation diagnoses why the original k=1 channel is special.
# A four-current vacuum Gram has degree at most2 in k, by affine commutation.
# All k=0 pairings vanish; two further levels determine the full polynomial.
af2=type(af)(ch,2);vals2=[];rows2=[];cols2=[]
states2=[]
for i,j in pairs:
    state={}
    for a,b in [(i,j),(j,i)]:
        for word,c in af2.act(G[a],-1,af2.act(G[b],-1,vac)).items():state[word]=state.get(word,F(0))+c/2
    states2.append({w:c for w,c in state.items() if c})
for block in weights.values():
    for p,q in itertools.product(block,repeat=2):
        i,j=pairs[p];inner=F(0)
        for a,b in [(i,j),(j,i)]:
            st=af2.act(ch.opp[G[a]],1,states2[q]);st=af2.act(ch.opp[G[b]],1,st)
            inner+=st.get((),F(0))/2
        inner*=mult[p]*mult[q]
        ck(inner.denominator==1,'actual level2 pair Gram integral')
        if inner:rows2.append(p);cols2.append(q);vals2.append(int(inner))
Gram2=csr_matrix((vals2,(rows2,cols2)),shape=K.shape,dtype=np.int64)
diff=Gram2-W@(6*s['I']+K);diff.eliminate_zeros()
ck(diff.nnz==0,'whole1596 level2 Gram fixes full level polynomial')
result['formal_level_polynomial']='Gram_k=W*(2*k*(k-1)*I+(k/2)*(K+2I))=W*(2*k*(k+5)*P133+2*k*(k-1)*P1463)'
result['formal_level_scope']='Same original current brackets with symbolic level diagnostic. Degree at most2 follows affine commutation; full matrix verified at k=1 and2, with zero constant term. Does not change original level1 or derive its physical selection.'
result['positive_level_consequence']='k>0 and positivity of these actual two-current norms imply k>=1; k=1 is the rank133 boundary. This necessary condition alone does not prove unitarity at arbitrary k or select k=1 without its existing premise.'
result['checks']=dict(s['env']['checks'])
s['OUT'].joinpath('two-current-gram-results.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
