#!/usr/bin/env python3
"""Exact joint Hodge-span TP identity and a sub-TP correlation counterexample.

Uses the existing integer epsilon amplitudes.  This proves a conditional finite
operator statement; it does not construct the original TFPT event map.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path
import sympy as s


SOURCE = Path(__file__).resolve().parents[1] / "hodge_complement_20261028" / "hodge_complement_channel.py"
spec = importlib.util.spec_from_file_location("tfpt_hodge_existing", SOURCE)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)
C = [s.Matrix(x.tolist()) for x in mod.hodge_amplitudes()]
I = s.eye(10)
checks = []

def check(condition, label):
    if not condition:
        raise RuntimeError(label)
    checks.append(label)

# The identity decomposes uniquely in the actual 25-product operator basis.
products = [C[n].H*C[m] for m in range(5) for n in range(5)]
basis = s.Matrix.hstack(*(p.reshape(100,1) for p in products))
check(basis.rank() == 25, "25 products are linearly independent")
target = s.Matrix([s.Rational(int(m == n), 3) for m in range(5) for n in range(5)])
check(basis*target == I.reshape(100,1), "unique TP coefficient matrix is I5/3")
check(sum((x.H*x for x in C), s.zeros(10)) == 3*I, "original amplitude normalization")
check(all((C[m].H*C[n]).trace() == 6*int(m == n) for m in range(5) for n in range(5)), "five orthogonal epsilon amplitudes")

# Symbolic, operator-valued coefficients on the family space, represented by
# an arbitrary 2x2 matrix per index pair.  This confirms the extraction used
# in the analytic proof; the proof itself holds in any family dimension.
g = s.symbols("g0:100")
G = [s.Matrix(2,2,g[4*j:4*j+4]) for j in range(25)]
joint = sum((s.kronecker_product(products[j], G[j]) for j in range(25)), s.zeros(20))
for u in range(2):
    for v in range(2):
        block = joint.extract([2*i+u for i in range(10)], [2*i+v for i in range(10)])
        coeff = s.Matrix([G[j][u,v] for j in range(25)])
        check(block.reshape(100,1) == basis*coeff, f"operator-valued coefficient extraction {u}{v}")

# A and B have identical carrier marginal but different environment correlation.
def ket(label):
    v=s.zeros(10,1); v[mod.PAIR_INDEX[label],0]=1; return v
e0=s.Matrix([1,0]); e1=s.Matrix([0,1]); p0=e0*e0.H; p1=e1*e1.H
v12=ket((1,2)); v02=ket((0,2)); v34=ket((3,4))
def proj(v): return v*v.H
def prodstate(v,e): return proj(s.kronecker_product(v,e))
A=(prodstate(v12,e0)+prodstate(v02,e1))/2
B=(prodstate(v12,e1)+prodstate(v02,e0))/2
def partial_family(M):
    return s.Matrix(10,10,lambda i,j:sum(M[2*i+k,2*j+k] for k in range(2)))
L=s.kronecker_product(C[0],p0)+s.kronecker_product(C[1],p1)
defect=s.eye(20)-L.H*L
check(partial_family(A)==partial_family(B), "same carrier marginal")
check(defect==defect.H and defect*defect==defect, "sub-TP defect is a positive projection")
check(partial_family(L*A*L.H)==proj(v34), "first correlated state outputs 34")
check(partial_family(L*B*L.H)==s.zeros(10), "second correlated state has zero outcome probability")
check((L*A*L.H).trace()==1 and (L*B*L.H).trace()==0, "outcome probabilities differ 1 versus 0")
print(f"PASS {len(checks)} exact checks")
for item in checks: print("  " + item)
