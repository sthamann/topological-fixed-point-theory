"""Full quadratic CAR data of the native five inputs and five q outputs."""

from functools import lru_cache
from itertools import product
import json
from pathlib import Path
import runpy
import sys

import numpy as np
import sympy as sp


repo = Path(sys.argv[1]).resolve()
base = repo/"experiments/theory-contracts/chat-triangulation-20261108"
old = runpy.run_path(str(base/"check_ramond_full_return_20261110.py"))
M, q, V, m, Q, cl = (old[k] for k in ("M","q","V","m","Q","cl"))
basis = sp.Matrix.hstack(V,q)
vectors = [[Q.from_sympy(basis[k,a]) for k in range(8)] for a in range(5)]
levels = list(range(-5,6))
width = len(levels)


def bit(channel, level):
    return channel*width+level-levels[0]


def sea(occupy_zero):
    return sum(1<<bit(i,n) for i in range(8) for n in levels
               if n<0 or (n==0 and occupy_zero[i]))


@lru_cache(None)
def step(state, order, channel):
    result = []
    for n in levels:
        if n+order not in levels:
            continue
        b,c = bit(channel,n),bit(channel,n+order)
        if state>>b&1 and not state>>c&1:
            between = ((1<<c)-1)^((1<<(b+1))-1)
            sign = (-1)**((state&between).bit_count())
            result.append((state^(1<<b)^(1<<c),sign))
    return result


@lru_cache(None)
def word_state(vac, word):
    if not word:
        return {vac:Q.one}
    order,axis = word[-1]
    ret = {}
    for state, coefficient in word_state(vac,word[:-1]).items():
        for channel,factor in enumerate(vectors[axis]):
            if not factor:
                continue
            for out, sign in step(state,order,channel):
                ret[out] = ret.get(out,Q.zero)+coefficient*factor*sign
    return {state:a for state,a in ret.items() if a}


def as_car(vac, polynomial):
    result = {}
    for word,value in polynomial.items():
        for state, coeff in word_state(vac,word).items():
            result[state] = result.get(state,Q.zero)+coeff*value
    return {state:a for state,a in result.items() if a}


def inner(a,b):
    return sum((m.conj(value)*b.get(key,Q.zero) for key,value in a.items()),Q.zero)


def diagonal_covariance(state, vac, grade):
    weights = [(s,m.conj(a)*a) for s,a in state.items()]
    norm = sum((a for _,a in weights),Q.zero)
    changed = {k:Q.zero for k in range(8*width)}
    for bits,weight in weights:
        difference = [(i,n,int((bits>>bit(i,n))&1)-int((vac>>bit(i,n))&1))
                      for i in range(8) for n in levels
                      if ((bits^vac)>>bit(i,n))&1]
        assert sum(n*sign for i,n,sign in difference)==grade
        assert all(sum(sign for i,n,sign in difference if i==j)==0 for j in range(8))
        assert all(n not in (levels[0],levels[-1]) for i,n,sign in difference)
        for i,n,sign in difference:
            changed[bit(i,n)] += weight*sign
    values = [Q.to_sympy(changed[k]/norm)+int(vac>>k&1) for k in range(8*width)]
    assert all(x.is_Rational and 0<=x<=1 for x in values)
    # Every c_i*c_j off diagonal changes a conserved individual charge or
    # the exact oscillator grade. Every cc lowers total charge by two.
    return {"car_basis_terms":len(state),"trace_N":str(sum(values)),
            "exact_empty":values.count(0),"exact_filled":values.count(1),
            "N_diagonal_by_channel_then_mode":[str(x) for x in values],
            "F_identically_zero":True,"N_off_diagonal_zero":True}


assert all(abs(x)==sp.Rational(1,2) for x in q)
vac_ns = sea([False]*8)
vac_r = sea([x>0 for x in q])
original = json.loads((base/"ramond-full-return-results-20261110.json").read_text())
polynomials,norms = m.code()
scale = sp.diag(*[1/sp.sqrt(n) for n in norms])
cases = {}
for alpha in (sp.Rational(1,4),sp.Rational(3,4)):
    delta = sp.Matrix([(j+alpha)/8-sp.Rational(1,2) for j in range(8)])
    drift = [Q.from_sympy(x) for x in cl(V.T*delta)]
    prepare = [old["series_response"](f,drift,4) for f in polynomials]
    ns = [as_car(vac_ns,old["input_word_map"](f)) for f in prepare]
    ramond = [as_car(vac_r,old["vertex_zero"](f)) for f in prepare]
    for key, states in (("input_Gram",ns),("full_return_Gram",ramond)):
        gram = cl(scale*sp.Matrix(5,5,lambda i,j:Q.to_sympy(inner(states[i],states[j])))*scale)
        assert gram==sp.Matrix(original[str(alpha)][key])
    cases[str(alpha)] = {"input":[diagonal_covariance(x,vac_ns,4) for x in ns],
                         "q_output":[diagonal_covariance(x,vac_r,3) for x in ramond]}

# The four original label projections resolve the regular two-involution
# action. The actual linear spin maps anticommute, since sigma3 sigma1=-sigma1 sigma3.
labels = list(product((1,-1),repeat=2))  # c/i, b
hmove = sp.zeros(4)
rmove = sp.zeros(4)
for i,(c,b) in enumerate(labels):
    hmove[labels.index((-c,-b)),i]=1
    rmove[labels.index((-c,b)),i]=b
assert hmove**2==rmove**2==sp.eye(4)
assert hmove*rmove*hmove*rmove==-sp.eye(4)
# Enumerate the label-resolved operator basis using the words I,H,R,HR.
units = set()
for i,j in product(range(4),repeat=2):
    if any(word[i,j]!=0 for word in (sp.eye(4),hmove,rmove,hmove*rmove)):
        units.add((i,j))
assert len(units)==16
# Without resolved labels the involutions have a two-dimensional Pauli model.
x = sp.Matrix([[0,1],[1,0]])
z = sp.diag(1,-1)
assert x*z*x*z==-sp.eye(2)

print(json.dumps({"covariance_levels":levels,"R_charge":list(map(str,q)),
                  "common_unexcited_tail":"filled below -5 and empty above +5",
                  "cases":cases,"input_and_output_full_grams_match_exactly":True,
                  "four_pin_labels":{"minimum_with_distinct_labels":4,
                     "resolved_operator_algebra_dimension":16,
                     "linear_spin_move_commutator_loop":-1,
                     "minimum_for_unresolved_abstract_involutions":2},
                  "scope":"Actual native charged Cartan CAR packet states at specified grades. No W mapping into the original physical cover, no QTheta."},
                 indent=2))
