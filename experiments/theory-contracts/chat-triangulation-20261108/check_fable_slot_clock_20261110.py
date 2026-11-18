#!/usr/bin/env python3
"""Exact audit of the slot reconstruction, from native v634/v17/v987 inputs."""
import ast
from collections import Counter, deque
import runpy
from itertools import combinations, product
import json
from pathlib import Path
import sys

import numpy as np
import sympy as sp
from sympy.combinatorics import Permutation, PermutationGroup

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "verification"))
import v17_hexagonal_resolvent as v17
from v987_os_dilation_package import B_SYM

N = 8
I = np.eye(N, dtype=np.int64)
e = I.T
roots = []
for i, j in combinations(range(N), 2):
    for a, b in product((-2, 2), repeat=2):
        roots.append(a * e[i] + b * e[j])
roots += [r for r in product((-1, 1), repeat=N) if sum(r) % 4 == 0]
R = np.array(sorted(map(tuple, roots)))
root_idx = {tuple(r): i for i, r in enumerate(R)}
ell = e[:5].sum(axis=0)
y6 = np.array([-2, -2, -2, 3, 3, 0, 0, 0])
weak = 2 * (e[3] - e[4])
beta = e[3] - e[5]
mirror = I - np.outer(beta, beta)
v634 = ast.parse((ROOT / "verification/v634_st31_structure.py").read_text())
native_idx = next(ast.literal_eval(n.value) for n in v634.body if isinstance(n, ast.Assign)
                  and any(getattr(t, "id", "") == "SIGMA_IDX" for t in n.targets))
j_vec = next(n for n in v634.body if isinstance(n, ast.FunctionDef) and n.name == "J_vec")
ns = {"np": np}
exec(compile(ast.Module(body=[j_vec], type_ignores=[]), "v634:J_vec", "exec"), ns)
J0 = np.column_stack([ns["J_vec"](v) for v in e])
sigma = I[native_idx, :]
c = J0 @ sigma

# Use an integral E8 simple-root basis. A mark is its eight residues mod 4.
A = np.column_stack([np.array([1,-1,-1,-1,-1,-1,-1,1]), 2*(e[0]+e[1]),
                     *(2*(e[i+1]-e[i]) for i in range(6))])
Ainv = np.array(4 * sp.Matrix(A).inv(), dtype=np.int64)
assert np.array_equal(Ainv @ A, 4*I)
assert abs(int(sp.Matrix(A).det())) == 2**8
X = R @ Ainv.T // 4
assert np.array_equal(X @ A.T, R)
simple = [I - np.outer(e[i], (A.T @ A)[i]//4) for i in range(N)]
assert all(np.array_equal(s@s, I) for s in simple)
mark = tuple(ell @ A % 4)
target = tuple(ell @ mirror @ A % 4)
G = A.T @ A
Gi = np.array(4 * sp.Matrix(G).inv(), dtype=np.int64)
def inverse(g):
    return Gi @ g.T @ G // 4

def key(g):
    return g.astype(np.int8).tobytes()

def on_roots(g):
    return [root_idx[tuple(r)] for r in (X @ g.T) @ A.T]

def orbit(gens, start, conjugate=False):
    seen = {(key(start) if conjugate else start): (start if conjugate else I)}
    todo = deque(seen)
    while todo:
        k = todo.popleft()
        g = seen[k]
        for s in gens:
            if conjugate:
                h = s @ g @ s
                z = key(h)
            else:
                h = g @ s
                z = tuple(np.array(k) @ s % 4)
            if z not in seen:
                seen[z] = h
                todo.append(z)
    return seen

def grade_census(g2):
    b = R.copy()
    vals = []
    for _ in range(12):
        n = b @ g2.T
        assert not (n % 2).any()
        b = n // 2
        vals.append(int(np.count_nonzero(((b-R) @ ell) % 4)))
    assert np.array_equal(b,R)
    return tuple(vals)

def cycle_list(points, f):
    unseen = set(points)
    cycles = []
    while unseen:
        a = min(unseen)
        cy = []
        while a in unseen:
            unseen.remove(a)
            cy.append(a)
            a = f(a)
        assert a == cy[0]
        cycles.append(cy)
    return cycles

def mixed_and_clock():
    assert np.array_equal(c @ c @ c, -J0) and np.array_equal(np.linalg.matrix_power(c,4),sigma)
    assert np.array_equal(np.linalg.matrix_power(c,9),J0)
    C = grade_census(2*c)
    assert C == (112,)*5 + (128,) + (112,)*5 + (0,)
    print("native powers 1..12:", C)
    assert int((c@weak)@y6)==-2
    ring = {tuple(r) for r in R if max(abs(r))==2 and
            (sum(abs(r[6:])) in (0,4))}
    assert len(ring)==64 and {tuple(c@r) for r in ring}==ring
    print("D6+D2 roots / split symmetric difference:",len(ring),48)
    print("native characteristic polynomial:", sp.factor(sp.Matrix(c).charpoly().as_expr()))
    sets = [set(np.argmax(abs(np.linalg.matrix_power(c,k)[:,5:]),axis=0)) for k in range(6)]
    assert [next(iter(s-{6,7}))+1 for s in sets] == [6,1,4,5,2,3]
    split = [{tuple(r) for r in R if max(abs(r))==2 and
              (set(np.flatnonzero(r)) <= s or not set(np.flatnonzero(r)) & s)} for s in sets]
    assert set(len(a^b) for a,b in combinations(split,2)) == {48}
    markings = [(np.linalg.matrix_power(c,k) @ ell) % 4 for k in range(12)]
    classes = [(R @ m) % 4 for m in markings]
    assert {np.count_nonzero(classes[i]!=classes[j]) for i,j in combinations(range(6),2)} == {112}
    assert all(np.count_nonzero(classes[k]!=classes[k+6])==128 for k in range(6))
    common = set.intersection(*split)
    assert len(common)==4  # Only D2 is in every D5+A3 subsystem.
    ring_axis = [tuple(r) for r in R if max(abs(r))==2 and sum(abs(r[6:]))==2]
    images = lambda r: tuple(c @ np.array(r))
    mixed = cycle_list(ring_axis,images)
    assert sorted(map(len,mixed)) == [12]*4
    # Import exact cover successor from v623, without running unrelated kernel checks.
    t = ast.parse((ROOT / "verification/v623_covered_seam.py").read_text())
    succ = next(n for n in t.body if isinstance(n,ast.FunctionDef) and n.name=="successor")
    cns = {"N":16,"MARK_BONDS":(0,4,8,12)}
    exec(compile(ast.Module(body=[succ],type_ignores=[]),"v623:successor","exec"),cns)
    cover = cycle_list(product(range(16),range(3)),lambda p:cns["successor"](*p,(1,)*4))[0]
    lift = {p:cover[(k+4)%48] for k,p in enumerate(cover)}
    cover_orbs = cycle_list(cover,lift.__getitem__)
    bij = {a:b for ac,bc in zip(mixed,cover_orbs) for a,b in zip(ac,bc)}
    assert len(cover)==len(bij)==len(set(bij.values()))==48
    assert all(bij[images(a)]==lift[b] for a,b in bij.items())
    print("mixed roots and cover cycles:", list(map(len,mixed)), "explicit C12 bijection")
    fold = ROOT/"experiments/theory-contracts/v16-operator-pullback-20261014/fold-results.json"
    M = sp.Matrix(json.loads(fold.read_text())["explicit_Weyl_matrix"])
    assert M.T*M==sp.eye(8)
    folded = M.T*sp.Matrix(sigma)*M
    print("source sigma(e1+e3):", list(folded * sp.Matrix(e[0]+e[2])))

def census():
    marks = orbit(simple,mark)
    assert target in marks
    stab_order = 696729600//len(marks)
    assert stab_order*len(marks)==696729600
    print("mark orbit / stabilizer:",len(marks),stab_order,flush=True)
    gens = []
    permgroup = PermutationGroup([Permutation(list(range(240)))])
    for k,g in marks.items():
        if permgroup.order()==stab_order:
            break
        for s in simple:
            z = tuple(np.array(k) @ s % 4)
            h = g@s@inverse(marks[z])
            assert tuple(np.array(mark) @ h % 4)==mark
            hperm = Permutation(on_roots(h))
            if not permgroup.contains(hperm):
                gens.append(h)
                permgroup = PermutationGroup([*permgroup.generators,hperm])
    assert permgroup.order()==stab_order
    seen = {key(I):I}
    q = deque(seen)
    while q:
        h = seen[q.popleft()]
        for gen in gens:
            g = h@gen
            k = key(g)
            if k not in seen:
                seen[k] = g
                q.append(k)
    assert len(seen)==stab_order
    candidates = []
    for h in seen.values():
        g = h @ marks[target]
        if np.array_equal(g@g@g,I) and not np.array_equal(g,I):
            g2 = A@g@Ainv//2  # doubled operator in orthogonal coordinates
            assert np.array_equal(2*A@g,g2@A)
            assert not (((R@g2.T//2)@ell-R@mirror.T@ell)%4).any()
            candidates.append(g2)
    Jb = Ainv@J0@A//4
    Js = np.array([A@j@Ainv//2 for j in orbit(simple,Jb,True).values()])
    assert len(Js)==15120
    counts = Counter()
    commuting_j0 = []
    for g in candidates:
        mask = (Js @ g == g @ Js).all(axis=(1,2))
        counts[(int((8+np.trace(g))//3), int(mask.sum()))] += 1
        if np.array_equal(g@J0,J0@g):
            dY = (R@g.T//2-R)@y6
            keep = (((R@g.T//2-R)@ell)%4)==2
            commuting_j0.append((g, int((g@weak//2)@y6), len(set(dY[keep]))))
    assert counts == {(6,0):12,(4,36):192}
    print("order-three (fixed dimension, commuting J) census:",dict(counts))
    print("J0 candidates (12Y of e4-e5 image, distinct DeltaY on glue changes):",
          [(q,n) for _,q,n in commuting_j0])
    chosen = [g for g,q,n in commuting_j0 if q==4 and n==8]
    assert len(chosen)==1 and np.array_equal(chosen[0],2*sigma)
    weak_pairs = []
    for g,_,_ in commuting_j0:
        pairs=[]
        for i,j in combinations(range(5),2):
            y = np.array([-2]*5+[0]*3)
            y[[i,j]]=3
            w=2*(e[i]-e[j])
            RR=R@g.T//2
            if abs(int((g@w//2)@y))==4 and len(set(((RR-R)@y)[((RR-R)@ell)%4==2]))==8:
                pairs.append((i+1,j+1))
        weak_pairs.append((bool(np.array_equal(g,2*sigma)),pairs))
    print("weak pairs with either orientation for each J0 candidate:",weak_pairs)
    Js_sigma = Js[(Js@(2*sigma)==(2*sigma)@Js).all(axis=(1,2))]
    patterns = Counter()
    for j in Js_sigma:
        monomial = bool(np.all(np.count_nonzero(j,axis=0)==1))
        patterns[(monomial,grade_census(j@sigma))] += 1
    assert patterns == {
        (False,(184,112)*2+(184,128)+(184,112)*2+(184,0)):24,
        (True,(112,)*5+(128,)+(112,)*5+(0,)):12,
    }
    for (monomial,powers), num in sorted(patterns.items()):
        print("native sigma alternative J:",monomial,num,powers)

def transfer():
    t = sp.symbols("t",real=True)
    u = v17.U6_power(1)
    eye = sp.eye(6)
    K = (sp.Rational(1,2)+t)*eye + (sp.Rational(2,9)-t)*(u+u.T) + t*(u**2+u.T**2) + (sp.Rational(1,18)-t)*u**3
    v = sp.zeros(6,3)
    for j,pair in enumerate(((1,2),(4,5),(0,3))):
        for k in pair:
            v[k,j]=1/sp.sqrt(2)
    assert sp.simplify(K*v-v*B_SYM)==sp.zeros(6,3)
    alt = sp.Matrix([(-1)**j for j in range(6)])
    assert sp.simplify(K*alt-6*t*alt)==sp.zeros(6,1)
    a = sp.symbols("a",positive=True)
    green = (eye+a*(2*eye-u-u.T)).inv()
    assert sp.simplify(green*alt-alt/(1+4*a))==sp.zeros(6,1)
    print("ring char poly:",sp.factor(K.charpoly().as_expr()), "native alt Green: 1/(1+4r)")
    # Native Clifford channel from Kernessenz p.11: 7/12 id + 1/12 sum Ad Gamma.
    Q = sp.Matrix([[0,1],[1,0]])
    Y = sp.Matrix([[0,-sp.I],[sp.I,0]])
    Z = sp.diag(1,-1)
    kp=sp.kronecker_product
    gam=[kp(Q,sp.eye(2)),kp(Y,sp.eye(2)),*(kp(Z,a) for a in (Q,Y,Z))]
    blocks=[(sp.eye(4),1),*((g,sp.Rational(1,3)) for g in gam),
            *((g*h,sp.Rational(2,3)) for g,h in combinations(gam,2))]
    assert sp.Matrix.hstack(*(m.reshape(16,1) for m,_ in blocks)).rank()==16
    for m,lam in blocks:
        assert 7*m+sum((g*m*g for g in gam),sp.zeros(4))==12*lam*m
    print("native Clifford spectrum:",Counter(lam for _,lam in blocks))

def charges96():
    # Extract the actual v252 real-bimodule zero mode as in the original Y audit.
    path = ROOT/"verification/v252_full_finite_triple.py"
    t = ast.parse(path.read_text())
    names={"NG","GEN_FIELDS","P","IDX","NP","NH"}
    funcs={"build_index","Kswap","rep"}
    nodes=[n for n in t.body if
           (isinstance(n,ast.FunctionDef) and n.name in funcs) or
           (isinstance(n,ast.Assign) and any(z.id in names for t in n.targets
                    for z in ast.walk(t) if isinstance(z,ast.Name)))]
    env={"np":np}
    exec(compile(ast.Module(body=nodes,type_ignores=[]),str(path),"exec"),env)
    K=env["Kswap"]()
    X=env["rep"](3j,np.zeros((2,2),complex),-1j*np.eye(3))
    Y=(X+K@X.conj()@K)/1j
    assert np.array_equal(Y,np.diag(np.diag(Y))) and np.array_equal(Y,-K@Y@K)
    actual=sum(int(x.real)**5 for x in np.diag(Y))/sp.Integer(6)**5
    def moment(removed,parity):
        coeff=np.delete(y6[:6],removed)
        return sum(int(coeff@signs)**5 for signs in product((-1,1),repeat=5)
                   if signs.count(-1)%2==parity)/sp.Integer(12)**5
    frames=[(moment(k,0),moment(k,1)) for k in range(6)]
    assert actual==0 and frames[:5]==[(0,0)]*5 and set(frames[5])=={-sp.Rational(5,9),sp.Rational(5,9)}
    print("fifth hypercharge moments: native 96",actual,"half-spinors on each fixed-Y frame",frames)

def packets():
    # Use the unchanged five inputs and full q-current consumer.
    p = runpy.run_path(str(HERE/"check_ramond_full_return_20261110.py"))
    M,q,V,cl = (p[x] for x in ("M","q","V","cl"))
    sig = M.T*sp.Matrix(sigma)*M
    J = M.T*sp.Matrix(J0)*M
    assert sig*q==q and q.dot(J*q)==0 and J*q!=q
    f,_ = p["m"].code()
    Q = p["Q"]
    qmat = sp.Matrix.hstack(V,q)
    met = cl(qmat.H*sig*qmat)
    inner_q = p["inner"].__globals__["inner_q"]
    inner_q[1] = [[Q.from_sympy(met[a,b]) for b in range(5)] for a in range(5)]
    p["contractions"].cache_clear()
    for alpha in map(sp.Rational, ("1/4","3/4")):
        delta=sp.Matrix([(j+alpha)/8-sp.Rational(1,2) for j in range(8)])
        b=[Q.from_sympy(x) for x in cl(V.T*delta)]
        fields=[p["vertex_zero"](p["series_response"](x,b,4)) for x in f]
        # Phase of a Weyl lift on the fixed q momentum is common and cannot change rank.
        G=cl(sp.Matrix(5,5,lambda i,j:p["inner"](fields[i],fields[j],False)))
        C=cl(sp.Matrix(5,5,lambda i,j:p["inner"](fields[i],fields[j],True)))
        rank = G.rank(),C.rank()
        print("q-packet Gram / native sigma compression rank",alpha,rank,
              "zero sigma diagonals",[i+1 for i in range(5) if C[i,i]==0],flush=True)
        assert rank==(5,5) and all(C[i,i]!=0 for i in range(5))
    print("native J fixed-q compression rank: 0 (different momentum)",flush=True)


if __name__ == "__main__":
    mixed_and_clock()
    transfer()
    census()
    charges96()
    packets()
