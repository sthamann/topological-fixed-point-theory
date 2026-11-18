"""Exact P1 closed return from the original Cartan oscillator modes."""
from contextlib import redirect_stdout
from itertools import product
import io
import json
from pathlib import Path
import runpy
import sys

repo=Path(sys.argv[1]).resolve()
base=repo/'experiments/theory-contracts/chat-triangulation-20261108'
with redirect_stdout(io.StringIO()):
    s=runpy.run_path(str(base/'check_tfpt_family_middle_ordered_cut_20261110.py'))
sp,Q,q,r,d1,V=(s[k] for k in ('sp','Q','q','r','d1','V'))
old,vertex,middle,inner,combine=(s[k] for k in ('old','vertex','middle','inner','combine'))
axes=s['axes']

def mul(a,f):
    return {word:f*c for word,c in a.items() if f*c}

def minus(a,b):
    return combine([*a.items(),*((w,-c) for w,c in b.items())])

def replace_axes(p,columns):
    terms=[]
    for mon,c in p.items():
        factors=[[(degree,i,x) for i,x in columns[j].items()] for degree,j in mon]
        for branch in product(*factors):
            word=tuple(sorted((degree,i) for degree,i,_ in branch))
            a=c
            for _,_,x in branch:
                a*=x
            terms.append((word,a))
    return combine(terms)

overlap=[Q.from_sympy(x) for x in (axes.T*q).applyfunc(sp.simplify)]
reflect=[combine([(j,Q.one),(4,-overlap[j])]) for j in range(6)]
for j in range(6):
    coeff=sp.zeros(6,1)
    for i,v in reflect[j].items():
        coeff[i]=Q.to_sympy(v)
    assert (axes*coeff-(sp.eye(8)-q*q.T)*axes[:,j]).applyfunc(sp.simplify)==sp.zeros(8,1)

def qcount(p):
    terms=[]
    # Q_osc = sum_n J^q_{-n} J^q_n/(2n); each annihilator
    # contributes n times the bilinear pairing. This avoids an extra cutoff.
    for mon,c in p.items():
        for pos,(degree,j) in enumerate(mon):
            word=tuple(sorted(mon[:pos]+mon[pos+1:]+((degree,4),)))
            terms.append((word,c*overlap[j]*Q.from_sympy(sp.Rational(1,2))))
    return combine(terms)

def norm(a):
    return sp.simplify(inner(a,a))

def even_minus_odd_polynomial(p):
    k1=qcount(p)
    k2=minus(qcount(k1),k1)
    k3=minus(qcount(k2),mul(k2,Q.from_sympy(sp.Integer(2))))
    return combine([*p.items(),*mul(k1,Q.from_sympy(sp.Integer(-2))).items(),
                    *mul(k2,Q.from_sympy(sp.Integer(2))).items(),*mul(k3,Q.from_sympy(sp.Rational(-4,3))).items()])

out={}
for a in map(sp.Rational,('1/4','3/4')):
    drift=sp.Matrix([j+a-4 for j in range(8)])
    c=[Q.from_sympy(z) for z in old['cl'](V.T*drift/8)]
    records=[]
    for i,poly in enumerate(s['polys']):
        v=old['vertex_zero'](old['series_response'](poly,c,4))
        B=vertex(v,d1,q)
        assert all(s['grade'](mon)==3 for mon in B)
        flipped=replace_axes(B,reflect)
        assert even_minus_odd_polynomial(B)==flipped
        D=minus(B,flipped)
        late0=vertex(B,q,r)
        twice=vertex(late0,-q,r+q)
        assert twice in (D,mul(D,-Q.one))
        assert sp.simplify(inner(B,D))==norm(late0)
        early0=vertex(vertex(v,q,q),d1,2*q)
        assert bool(early0)==(i==0)
        late1=vertex(B,q,r,1)
        lower1=vertex(B,-q,r,-1)
        assert sp.simplify(norm(late1)-norm(lower1)-norm(B))==0
        assert norm(late1)!=norm(late0)
        early1=vertex(vertex(v,q,q,1),d1,2*q)
        assert early1
        records.append({'packet':i,'early_zero_mode_norm2':str(norm(early0)),
                        'late_zero_mode_norm2':str(norm(late0)),
                        'early_local_mode_norm2':str(norm(early1)),
                        'late_local_mode_norm2':str(norm(late1)),
                        'opposite_local_mode_norm2':str(norm(lower1)),
                        'local_Ward_difference_equals_unmarked_norm2':True,
                        'same_neutral_parity_insufficient_for_local_mode':True})
    out[str(a)]=records
print(json.dumps({'actual_old_source_all_five_packets':out,'physical_values':None},indent=2))
