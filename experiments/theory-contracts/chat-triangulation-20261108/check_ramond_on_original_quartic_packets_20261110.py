#!/usr/bin/env python3
"""An exact projected charged response on the five original grade-32 packets."""
from pathlib import Path
import json
from math import factorial
from itertools import permutations
from collections import Counter
import sys

import sympy as sp

ROOT = Path(sys.argv[1]).resolve()
LEG = ROOT / "experiments/theory-contracts/joint-native-algebra-selection-20261016"
sys.path[:0] = [str(LEG),str(ROOT)]
import check_mixed_return as m
from check_fold_quartic_transport import cl

FOLD = ROOT / "experiments/theory-contracts/v16-operator-pullback-20261014/fold-results.json"
M = sp.Matrix(json.loads(FOLD.read_text())["explicit_Weyl_matrix"])
Q = m.Q
E = sp.eye(8)
q = M.T*(E[:,7]-E[:,6])
ld = m._lattice_data()
T = sp.Matrix.hstack(*map(sp.Matrix,m._TARGET_SIMPLE))
O = cl(T*ld["basis"].inv()/sp.sqrt(2))
U = sp.zeros(8,4)
for a in range(4):
    U[2*a,a]=1/sp.sqrt(2)
    U[2*a+1,a]=-sp.I/sp.sqrt(2)
V = cl(M.T*O*U)
c = cl(V.T*q)
assert cl(V.H*V)==sp.eye(4)
assert cl(V.T*V)==sp.zeros(4)
assert cl(c.H*c)==sp.eye(1)
vp = cl((sp.eye(8)-q*q.T/2)*V)
assert cl(vp.H*vp)==cl(sp.eye(4)-sp.conjugate(c)*c.T/2)
assert cl(vp.H*vp).det().simplify()==sp.Rational(1,2)

# f_i is the unshifted quartic and is exactly the all-mode-1 component
# of each full source R_{b,4} f_i in both quarters.
quartics, norms = m.code()
ders = []
for f in quartics:
    result = {}
    for mon,coeff in f.items():
        for i,a in enumerate(mon):
            sub = mon[:i]+mon[i+1:]
            result[sub] = result.get(sub,0)+Q.to_sympy(coeff)*c[a]
    ders.append({k:sp.simplify(v) for k,v in result.items() if sp.simplify(v)!=0})
keys = sorted(set().union(*ders))
matrix = sp.Matrix(len(keys),len(ders),lambda i,j:ders[j].get(keys[i],0))
scale = sp.diag(*[sp.sqrt(n) for n in norms])
matrix = cl(matrix*scale.inv())
assert matrix.rank()==5
gram = cl(matrix.H * sp.diag(*[
    sp.prod(factorial(x) for x in Counter(k).values()) for k in keys
]) * matrix)
assert gram==sp.eye(5)
projected_metric = cl(vp.H*vp)
perms = list(permutations(range(3)))
hall_projected = sp.Matrix(len(keys),len(keys),lambda i,j:sp.simplify(sum(
    sp.prod(projected_metric[keys[i][a],keys[j][b]] for a,b in enumerate(perm)) for perm in perms
)))
projected_gram = cl(matrix.H*hall_projected*matrix)
assert projected_gram.rank()==5
# Cubic projected mode-1 metric induced by vp has all eigenvalues >= (1/2)^3.
out = {
    "ramond_q":list(map(str,q)),
    "contraction_Vt_q":list(map(str,c)),
    "number_of_original_packets":len(quartics),
    "rank_of_projected_ramond_mode_response":matrix.rank(),
    "response_norm_lower_bound":"norm^2(J_q,0 P) >= v* Gram v/8 in the original native-normalized basis",
    "Gram": [[str(x) for x in row] for row in gram.tolist()],
    "Gram_characteristic_polynomial":str(sp.factor(gram.charpoly().as_expr())),
    "exact_projected_response_Gram": [[str(x) for x in row] for row in projected_gram.tolist()],
    "exact_projected_response_Gram_characteristic_polynomial":str(sp.factor(projected_gram.charpoly().as_expr())),
    "common_original_K":4,
    "applies_to_both_original_quarters":True,
    "scope":"same original grade-32 inputs, given E8 extension, full physical P_adm and local P1 q remain uncomputed",
}
print(json.dumps(out, indent=2))
