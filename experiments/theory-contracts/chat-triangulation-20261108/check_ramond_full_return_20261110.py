#!/usr/bin/env python3
"""Full same-state Ramond return in the existing E8 parent, separate from P_adm."""
from collections import Counter
from functools import lru_cache
from itertools import combinations, permutations, product
from math import factorial
from pathlib import Path
import json
import sys

import sympy as sp

ROOT = Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path("/Users/stefanhamann/Projekte/tfpt-theoryv4")
LEG = ROOT / "experiments/theory-contracts/joint-native-algebra-selection-20261016"
sys.path[:0] = [str(LEG),str(ROOT)]
import check_mixed_return as m
from check_fold_quartic_transport import cl, series_response

M = sp.Matrix(json.loads((ROOT / "experiments/theory-contracts/v16-operator-pullback-20261014/fold-results.json").read_text())["explicit_Weyl_matrix"])
q = M.T * (sp.eye(8)[:,7]-sp.eye(8)[:,6])
ld = m._lattice_data()
T = sp.Matrix.hstack(*map(sp.Matrix,m._TARGET_SIMPLE))
O = cl(T*ld["basis"].inv()/sp.sqrt(2))
U = sp.zeros(8,4)
for a in range(4):
    U[2*a,a]=1/sp.sqrt(2)
    U[2*a+1,a]=-sp.I/sp.sqrt(2)
V = cl(M.T*O*U)
c = cl(V.T*q)
Q = m.Q
zq = [Q.from_sympy(x) for x in c]

def exp_y(weight):
    for n1,n2,n3 in product(range(4),range(2),range(2)):
        if n1+2*n2+3*n3!=weight:
            continue
        yield ((1,4),)*n1+((2,4),)*n2+((3,4),)*n3, Q.from_sympy(
            sp.Rational(1,factorial(n1)*factorial(n2)*factorial(n3)*2**n2*3**n3)
        )

def vertex_zero(response):
    result={}
    for word, coeff in response.items():
        factors = [(i//8+1,i%8) for i in word]
        for n in range(1,len(factors)+1):
            for deleted in combinations(range(len(factors)),n):
                weight=sum(factors[k][0] for k in deleted)
                other=tuple(factors[k] for k in range(len(factors)) if k not in deleted)
                a=coeff
                for k in deleted:
                    a=a*(-zq[factors[k][1]])
                for y,b in exp_y(weight-1):
                    key=tuple(sorted(other+y))
                    result[key]=result.get(key,Q.zero)+a*b
    assert all(sum(n for n,a in k)==3 for k in result)
    return {k:v for k,v in result.items() if v}

inner_metric = sp.eye(5)
for i in range(4):
    inner_metric[i,4]=sp.conjugate(c[i])
    inner_metric[4,i]=c[i]
inner_metric[4,4]=2
Pswap = sp.eye(8)[[4,5,6,7,0,1,2,3], :]
assert Pswap*q==q
swap_metric = inner_metric.copy()
swap_metric[:4,:4]=cl(V.H*Pswap*V)
assert cl(swap_metric.H)==swap_metric
inner_q = [[[Q.from_sympy(metric[a,b]) for b in range(5)] for a in range(5)]
           for metric in (inner_metric,swap_metric)]

@lru_cache(None)
def contractions(x,y,swap=False):
    if Counter(n for n,a in x)!=Counter(n for n,a in y):
        return Q.zero
    total=Q.zero
    for perm in permutations(y):
        p=Q.one
        for (n,a),(m,b) in zip(x,perm):
            if n!=m:
                p=Q.zero
                break
            p=p*Q.from_sympy(sp.Integer(n))*inner_q[swap][a][b]
        total=total+p
    return total

def inner(x,y,swap=False):
    result=Q.zero
    for k,v in x.items():
        for l,w in y.items():
            result=result+m.conj(v)*w*contractions(k,l,swap)
    return Q.to_sympy(result)

def input_word_map(response):
    return {tuple((i//8+1,i%8) for i in word):coeff for word,coeff in response.items()}

def projected_response(x):
    # Replace h_a by its q-perpendicular component; all y factors then vanish.
    return {k:v for k,v in x.items() if all(n==1 and a!=4 for n,a in k)}

def projected_metric_check(left,right):
    from itertools import permutations
    mat=cl(V.H*(sp.eye(8)-q*q.T/2)*V)
    g=0
    for x,v in left.items():
        for y,w in right.items():
            g+=sp.conjugate(Q.to_sympy(v))*Q.to_sympy(w)*sum(
                sp.prod(mat[a,b] for (_,a),(_,b) in zip(x,perm)) for perm in permutations(y)
            )
    return sp.simplify(g)

def encode(mat):
    return [[str(v) for v in row] for row in mat.tolist()]

def spectral_gap(alpha):
    w=(1,-1,0,-1,-1,-1,0,-1)
    d=[(j+alpha)/8-sp.Rational(1,2) for j in range(8)]
    candidates=[]
    for signs in product((-1,1),repeat=8):
        p=[sp.Rational(k,2) for k in signs]
        if signs.count(-1)%2==0 and sum(p)==0 and sum(x*y for x,y in zip(p,w))%3==0:
            candidates.append((1+sum(x*y for x,y in zip(d,p)),p))
    minimum=min(x for x,_ in candidates)
    # Any |p_i|>=3/2 costs at least 1-max|δ| above unrestricted min 0.
    assert sum(map(abs,d))==2 and minimum < 1-max(map(abs,d))
    return {"K_minimum":str(minimum),"H_minimum":str(8*minimum),
            "minimizers":[list(map(str,p)) for k,p in candidates if k==minimum],
            "higher_coordinate_exclusion_bound":str(1-max(map(abs,d)))}

def main():
    f,norms = m.code()
    scale=sp.diag(*[sp.sqrt(n) for n in norms])
    out={}
    old=json.loads((LEG/"original_density_current_results.json").read_text())["cases"]
    for alpha in (sp.Rational(1,4),sp.Rational(3,4)):
        delta=sp.Matrix([(j+alpha)/8-sp.Rational(1,2) for j in range(8)])
        assert delta.dot(q)==0
        b=[Q.from_sympy(x) for x in cl(V.T*delta)]
        prepared=[series_response(x,b,4) for x in f]
        image=list(map(vertex_zero,prepared))
        G=cl(scale.inv()*sp.Matrix(5,5,lambda i,j:inner(image[i],image[j]))*scale.inv())
        Gi=cl(scale.inv()*sp.Matrix(5,5,lambda i,j:Q.to_sympy(m.inner(prepared[i],prepared[j])))*scale.inv())
        assert cl(Gi-sp.Matrix(old[str(alpha)]["full_input_gram"]))==sp.zeros(5)
        for difference in (G-Gi/5,sp.Rational(17,20)*Gi-G):
            assert all(sp.simplify(difference[:k,:k].det())>0 for k in range(1,6))
        Gi_swap=cl(scale.inv()*sp.Matrix(5,5,lambda i,j:inner(input_word_map(prepared[i]),input_word_map(prepared[j]),True))*scale.inv())
        G_swap=cl(scale.inv()*sp.Matrix(5,5,lambda i,j:inner(image[i],image[j],True))*scale.inv())
        assert G_swap.is_hermitian and Gi_swap.is_hermitian
        in_plus=cl((Gi+Gi_swap)/2)
        in_minus=cl((Gi-Gi_swap)/2)
        out_plus=cl((G+G_swap)/2)
        out_minus=cl((G-G_swap)/2)
        assert [in_plus.rank(),in_minus.rank()]==[5,5]
        assert [out_plus.rank(),out_minus.rank()]==[5,4]
        assert out_minus[0,0]==0
        projected=list(map(projected_response,image))
        Gp=cl(scale.inv()*sp.Matrix(5,5,lambda i,j:projected_metric_check(projected[i],projected[j]))*scale.inv())
        assert sp.factor(Gp.charpoly().as_expr())==sp.factor((sp.Symbol("lambda")-1)*
            (sp.Symbol("lambda")-sp.Rational(1,2))**3*(sp.Symbol("lambda")-sp.Rational(1,8)))
        assert G.is_hermitian and all(sp.simplify((G-Gp)[:k,:k].det())>=0 for k in range(1,6))
        out[str(alpha)] = {"full_return_Gram":encode(G), "input_Gram":encode(Gi),
            "full_return_polynomial":str(sp.factor(G.charpoly().as_expr())),
            "exact_normalized_relative_return_bounds":["1/5","17/20"],
            "relative_return_polynomial":str(sp.factor((G-sp.Symbol("lambda")*Gi).det())),
            "individual_normalized_return_coefficients":[str(sp.simplify(G[i,i]/Gi[i,i])) for i in range(5)],
            "actual_center_Q0_s_sector_K_and_H_gap":spectral_gap(alpha),
            "old_half_swap_native_lift_controls":{
                "input_plus_Gram":encode(in_plus), "input_minus_Gram":encode(in_minus),
                "input_plus_and_minus_ranks":[in_plus.rank(),in_minus.rank()],
                "input_plus_fraction":[str(sp.simplify(in_plus[i,i]/Gi[i,i])) for i in range(5)],
                "out_plus_Gram_kappa_plus":encode(out_plus), "out_plus_Gram_kappa_minus":encode(out_minus),
                "output_plus_ranks_kappa_plus_minus":[out_plus.rank(),out_minus.rank()],
                "same_plus_projection_at_input_and_output_ranks_kappa_plus_minus":[out_plus.rank(),0],
                "output_plus_fractions_kappa_plus":[str(sp.simplify(out_plus[i,i]/G[i,i])) for i in range(5)],
                "output_plus_fractions_kappa_minus":[str(sp.simplify(out_minus[i,i]/G[i,i])) for i in range(5)]}}
    Path(__file__).with_name("ramond-full-return-results-20261110.json").write_text(json.dumps(out,indent=2)+"\n")
    for alpha,x in out.items():
        print(alpha,x["individual_normalized_return_coefficients"],x["actual_center_Q0_s_sector_K_and_H_gap"])
        p=x["old_half_swap_native_lift_controls"]
        print(p["input_plus_and_minus_ranks"],p["output_plus_ranks_kappa_plus_minus"],p["input_plus_fraction"],p["output_plus_fractions_kappa_plus"])

if __name__=="__main__":
    main()
