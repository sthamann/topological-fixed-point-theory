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
    v=old['vertex_zero'](old['series_response'](s['polys'][1],c,4))
    B=vertex(v,d1,q)
    W=middle(v)
    assert all(s['grade'](mon)==3 and len(mon)<=3 for mon in B)
    flipped=replace_axes(B,reflect)
    assert replace_axes(flipped,reflect)==B
    assert even_minus_odd_polynomial(B)==flipped
    odd=mul(minus(B,flipped),Q.from_sympy(sp.Rational(1,2)))
    assert replace_axes(odd,reflect)==mul(odd,-Q.one)
    W2=vertex(vertex(B,q,r),-q,r+q)
    D=minus(B,flipped)
    assert W2 in (D,mul(D,-Q.one)) # original adjoint convention fixed by positivity
    assert sp.simplify(inner(B,D))==norm(W)==2*norm(odd)
    assert not vertex(v,q,q)
    E=32+drift.dot(r)
    t,delta=s['t'],s['delta_t']
    return_value=sp.simplify(norm(W)/norm(v))
    reflect_value=sp.simplify(inner(B,flipped)/norm(v))
    projections=[]
    for k in range(4):
        b=B
        for l in range(4):
            if k!=l:
                b=mul(minus(qcount(b),mul(b,Q.from_sympy(sp.Integer(l)))),Q.from_sympy(sp.Rational(1,k-l)))
        assert qcount(b)==mul(b,Q.from_sympy(sp.Integer(k)))
        projections.append(b)
    assert combine([(w,c) for p in projections for w,c in p.items()])==B
    assert all(inner(x,y)==0 for i,x in enumerate(projections) for y in projections[i+1:])
    assert norm(projections[1])+norm(projections[3])==norm(odd)
    out[str(a)]={
        'unmarked_charged_B_norm2_over_q_reference':str(sp.simplify(norm(B)/norm(v))),
        'with_closed_Cartan_parity_same_reference':str(reflect_value),
        'late_difference_unweighted_exact':str(return_value),
        'early_null_mode_closed_return':0,
        'late_reflected_old_time_same_reference':float((return_value*sp.exp(-2*E*(t+delta))).evalf(20)),
        'polynomial_degree':3,
        'Cartan_count_k0_through_k3_weight_over_input':[str(norm(p)/norm(v)) for p in projections],
    }
print(json.dumps({'old_P1_same_reference':out,'physical_weighted_value':None},indent=2))
