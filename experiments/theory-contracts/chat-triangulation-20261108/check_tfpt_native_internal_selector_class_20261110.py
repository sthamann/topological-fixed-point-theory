"""Test determinant candidates on the actual finite native field algebra."""

import contextlib
import io
from itertools import combinations, product
import json
from pathlib import Path
import runpy
import sys

import numpy as np


repo = Path(sys.argv[1]).resolve()
base = repo/"experiments/theory-contracts/chat-triangulation-20261108"
with contextlib.redirect_stdout(io.StringIO()):
    source = runpy.run_path(str(base/"check_tfpt_conditioned_r113_weights_20261110.py"))
native, gamma, swap, dm = (source[k] for k in ("native","gamma","swap","majorana"))
lepton = np.array([a[0] in ("nu","e") for a in native.P]*2)
left = np.array([a[1]=="L" for a in native.P]*2)
labels = np.arange(96)
eye = np.eye(96)
q_lepton = np.diag(1-2*lepton.astype(int))
lr = np.diag(2*left.astype(int)-1)
c2 = np.zeros((2,2),complex)
c3 = np.zeros((3,3),complex)


def lift(matrix):
    return -1j*(matrix+swap@matrix.conj()@swap)


weak = [lift(native.rep(0,1j*w/2,c3)) for w in
        (np.array([[0,1],[1,0]]),np.array([[0,-1j],[1j,0]]),np.diag([1,-1]))]
color_diagonal = [lift(native.rep(0,c2,1j*np.diag(d))) for d in ([1,-1,0],[1,1,-2])]
hyper = lift(native.rep(3j,c2,-1j*np.eye(3)))/6
color_roots = []
for a,b in ((0,1),(1,2)):
    block = c3.copy()
    block[a,b] = block[b,a] = 1j
    color_roots.append(lift(native.rep(0,c2,block)))
    block = c3.copy()
    block[a,b],block[b,a] = 1,-1
    color_roots.append(lift(native.rep(0,c2,block)))
higgs = [source["inner"](np.array(v)) for v in ([1,0],[1j,0],[0,1],[0,1j])]
diagonals = [gamma,weak[2],hyper,*color_diagonal]
assert all(np.max(abs(op-np.diag(op.diagonal())))==0 for op in diagonals)
eigen = np.stack([op.diagonal().real for op in diagonals],axis=1)

prime = 65521
basis = []
for i in labels:
    z = np.zeros((96,96),complex)
    z[i,i] = 1
    basis.append(z)
for i,j in combinations(labels,2):
    if not np.array_equal(eigen[i],eigen[j]):
        continue
    for q in (1,1j):
        z = np.zeros((96,96),complex)
        z[i,j],z[j,i] = q,np.conj(q)
        basis.append(z)
projected = {}
for z in basis:
    real = z+swap@z.conj()@swap
    k = np.flatnonzero(real)
    if not len(k):
        continue
    first = real.ravel()[k[0]]
    sign = 1 if (first.real or first.imag)>0 else -1
    key = tuple((int(i),complex(sign*real.ravel()[i])) for i in k)
    projected[key] = sign*real
basis = np.stack(list(projected.values()))
# These Hermitian J-compatible basis vectors have disjoint real/imag
# supports and are an exact rational basis of the retained candidate space.
support = np.concatenate([basis.real.reshape(len(basis),-1)!=0,
                          basis.imag.reshape(len(basis),-1)!=0],axis=1)
assert np.max(support.sum(axis=0))==1
constraints = weak[:2]+color_roots


def residue(value):
    n,d = float(value).as_integer_ratio()
    return n*pow(d,-1,prime)%prime


reduce_mod = np.vectorize(residue,otypes=[np.int64])
br,bi = reduce_mod(basis.real),reduce_mod(basis.imag)


def product_mod(a,b,c,d):
    return (a@c-b@d)%prime,(a@d+b@c)%prime


def blocks(op):
    a,b=reduce_mod(op.real),reduce_mod(op.imag)
    left_pair = product_mod(br,bi,a,b)
    right_pair = product_mod(a,b,br,bi)
    def flatten(pair):
        return np.concatenate(pair,axis=1).reshape(len(basis),-1).T
    return flatten(left_pair),flatten(right_pair)


def add_gram(normal, equations):
    rows = equations[np.any(equations,axis=1)]
    normal += rows.T@rows
    normal %= prime


def rank_mod(matrix):
    matrix = matrix.copy()
    pivot = 0
    for col in range(matrix.shape[1]):
        candidates = np.flatnonzero(matrix[pivot:,col])
        if not len(candidates):
            continue
        i = pivot+candidates[0]
        matrix[[pivot,i]]=matrix[[i,pivot]]
        matrix[pivot] = matrix[pivot]*pow(int(matrix[pivot,col]),-1,prime)%prime
        matrix[pivot+1:] = (matrix[pivot+1:]
                             -matrix[pivot+1:,col,None]*matrix[pivot])%prime
        pivot += 1
        if pivot==len(matrix):
            break
    return pivot


k = len(basis)
gauge = np.zeros((k,k),dtype=np.int64)
for op in constraints:
    l,r = blocks(op)
    add_gram(gauge,(l-r)%prime)
higgs_blocks = [blocks(op) for op in higgs]
dm_blocks = blocks(dm)
results = {}
for sign in (1,-1):
    one = gauge.copy()
    paired = np.block([[gauge,np.zeros_like(gauge)],[np.zeros_like(gauge),gauge]])
    for (l,r),eta in zip(higgs_blocks+[dm_blocks],[sign]*4+[1],strict=True):
        add_gram(one,(l-eta*r)%prime)
        add_gram(paired,np.column_stack([-eta*r,l])%prime)
        add_gram(paired,np.column_stack([l,-eta*r])%prime)
    dimension_one = k-rank_mod(one)
    dimension_pair = 2*k-rank_mod(paired)
    assert dimension_one==2 and dimension_pair==3
    witness = eye if sign==1 else lr
    quark_lr = np.diag(np.where(lepton,0,lr.diagonal()))
    lepton_proj = np.diag(lepton.astype(int))
    quark_proj = eye-lepton_proj
    choices = [(witness@quark_proj,witness@quark_proj),
               (witness@lepton_proj,witness@lepton_proj),
               (witness@quark_lr,-witness@quark_lr)]
    for first,second in choices:
        for value in (first,second):
            assert np.array_equal(value,swap@value.conj()@swap)
            for op in diagonals+constraints:
                assert np.array_equal(value@op,op@value)
        for op,eta in zip(higgs+[dm],[sign]*4+[1],strict=True):
            assert np.array_equal(second@op,eta*op@first)
            assert np.array_equal(first@op,eta*op@second)
    involutions = []
    for q_left,q_right,l_sign in product((1,-1),repeat=3):
        a = (q_left+q_right)/2
        b = (q_left-q_right)/2
        first = a*choices[0][0]+l_sign*choices[1][0]+b*choices[2][0]
        second = a*choices[0][1]+l_sign*choices[1][1]+b*choices[2][1]
        assert np.array_equal(first@first,eye) and np.array_equal(second@second,eye)
        involutions.append([int(np.count_nonzero(value.diagonal()==1)) for value in (first,second)])
    results[str(sign)] = {"common_internal_dimension_exact":dimension_one,
                          "label_dependent_internal_dimension_exact":dimension_pair,
                          "number_of_label_dependent_involutions":len(involutions),
                          "positive_ranks_by_label":sorted(involutions),
                          "modular_rank_certificate_prime":prime,
                          "Hermitian_J_compatible_input_dimension_per_label":k}

# Both I and the quark/lepton sign commute with the actual entire fixed
# field and preserve the weak determinant on each native doublet (sign squared).
# They give different selections on the unnormalized rank-96 Higgs span.
proj = (eye+q_lepton)/2
assert np.linalg.matrix_rank(proj)==72
assert np.linalg.matrix_rank(proj@higgs[-2])==72

print(json.dumps({"native_internal_candidate_relations":results,
                  "exact_even_candidates":{"I_positive_rank":96,"quark_lepton_sign_positive_rank":72,
                    "both_preserve_native_doublet_determinant":True},
                  "scope":"Label-preserving constant real-compatible Hermitian internal matrices on v252 and reflection-related Pin labels, actual gauge/gradation and Higgs/Majorana relations; no representation of physical QTheta or Dadm."},
                 indent=2))
