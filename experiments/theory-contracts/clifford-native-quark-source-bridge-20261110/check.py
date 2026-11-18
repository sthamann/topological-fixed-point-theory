"""Exact finite checks and a replay against the existing TFPT quark builder."""

from contextlib import redirect_stdout
from itertools import combinations, combinations_with_replacement, product
from pathlib import Path
import io
import json
from math import factorial, prod
import runpy
import sys

import numpy as np
from scipy.linalg import polar
import sympy as s

repo = Path(sys.argv[1]).resolve()
contract = repo / "experiments/theory-contracts/chat-triangulation-20261108"
I = s.I
X = s.Matrix([[0, 1], [1, 0]])
Y = s.Matrix([[0, -I], [I, 0]])
Z = s.diag(1, -1)
eye = s.eye(4)
tensor = s.kronecker_product
g = [tensor(X, X), tensor(X, Y), tensor(X, Z), tensor(Y, s.eye(2)), tensor(Z, s.eye(2))]

def phi(a):
    return (7*a + sum((b*a*b for b in g), s.zeros(4))) / 12

assert all(a*b+b*a == 2*int(j == k)*eye for j, a in enumerate(g) for k, b in enumerate(g))
assert phi(eye) == eye
assert all(phi(a) == a/3 for a in g)
assert all(phi(a*b) == 2*a*b/3 for a, b in combinations(g, 2))
p = (eye-g[-1])/2
q = eye-p
for i, j in product((2, 3), repeat=2):
    unit = s.zeros(4)
    unit[i, j] = 1
    assert phi(unit) == 2*unit/3 + s.trace(unit)*q/6
    assert p*phi(phi(unit))*p == 4*unit/9 + s.trace(unit)*p/18
assert g[0]*g[1]*g[2]*g[3]*g[4] in (eye, -eye)
u, e = s.symbols("u e", positive=True)
density = 6*u/(1+u)**4
assert s.integrate(density, (u, 0, s.oo)) == 1
assert s.integrate(u*density, (u, 0, s.oo)) == 2
assert s.cancel(s.integrate(density, (u, 0, e)) - e**2*(3+e)/(1+e)**3) == 0
assert s.cancel(density - density.subs(u, 1/u)/u**2) == 0

roots = []
for i, j in combinations(range(8), 2):
    for a, b in product((-2, 2), repeat=2):
        roots.append(tuple(a if k == i else b if k == j else 0 for k in range(8)))
roots.extend(v for v in product((-1, 1), repeat=8) if v.count(-1) % 2 == 0)
assert len(set(roots)) == 240
for k, target in [(2, 60), (4, 36), (6, 30)]:
    for indices in combinations_with_replacement(range(8), k):
        powers = [indices.count(i) for i in range(8)]
        actual = sum(prod(v[i] for i in indices) for v in roots)
        expected = 0
        if all(n % 2 == 0 for n in powers):
            expected = target * 2**k * prod(factorial(n) for n in powers) * factorial(k//2)
            expected //= factorial(k) * prod(factorial(n//2) for n in powers)
        assert actual == expected
m8_axis = sum(s.Rational(v[0], 2)**8 for v in roots)
m8_diag = sum(s.Rational(v[0]+v[1], 2)**8/16 for v in roots)
assert (m8_diag-m8_axis)/s.factorial(8) == s.Rational(1, 3840)

def glue(v):
    half = v[0] % 2
    return (half + sum(x-half for x in v[:5])) % 4

def rotate(v):
    return tuple(-v[i+1] if i % 2 == 0 else v[i-1] for i in range(8))

def state(v):
    a, b = glue(v), glue(rotate(v))
    assert (b-a) % 2 == 0
    return a, ((b-a)//2) % 2

states = set()
common_roots = []
for v in roots:
    j = rotate(v)
    assert rotate(j) == tuple(-x for x in v)
    a, bit = state(v)
    assert state(j) == ((a+2*bit) % 4, (bit+a) % 2)
    states.add((a, bit))
    if glue(v) == 0 and glue(j) == 0:
        common_roots.append(v)
assert len(states) == 8
assert s.Matrix(common_roots).rank() == 6
axis_56 = [tuple(4*int(i == j) for i in range(8)) for j in (4, 5)]
assert all(glue(v) == glue(rotate(v)) == 0 for v in axis_56)
assert s.Matrix(common_roots + axis_56).rank() == 8

x, a = s.symbols("x a", nonzero=True)
f = x**5 + a*x**4
assert s.factor(s.diff(f, x)) == x**3*(5*x+4*a)
assert f.subs(x, -4*a/5) == 256*a**5/3125
assert s.diff(f, x, 2).subs(x, -4*a/5) != 0

def chain(n):
    return s.eye(n) - s.Matrix(n, n, lambda i, j: int(j == i+1))

seifert = tensor(chain(4), chain(2))
gram = seifert + seifert.T
assert (gram.det(), gram[:6, :6].det(), gram[6:, 6:].det(), gram[:6, 6:].rank()) == (1, 3, 3, 2)
ux, uz = (s.eye(2)-I*X)/s.sqrt(2), (s.eye(2)-I*Z)/s.sqrt(2)
comm = (ux*uz*ux.H*uz.H).applyfunc(s.simplify)
assert s.trace(comm) == 1 and (comm**3).applyfunc(s.simplify) == -s.eye(2)
assert s.gcd(s.cyclotomic_poly(30, x), s.cyclotomic_poly(12, x)*s.cyclotomic_poly(3, x)) == 1
assert s.totient(s.ilcm(4,30)) == 2*s.totient(30)
assert min((15*n+m,n,m) for n in range(1,40) for m in range(16) if 3*n%9 == (15*n+m)%16 == 0) == (48,3,3)

sys.path.insert(0, str(repo / "verification"))
with redirect_stdout(io.StringIO()):
    sys.argv = ["check_r113_dressed_open_20261109.py", str(contract)]
    native = runpy.run_path(str(contract / sys.argv[0]))
rep, grading, majorana, inner = (native[key] for key in ("native","grading","majorana","inner_form"))
modes = [rep.IDX[(f,side,0,j)] for side in ("L","R") for f in ("u","d") for j in range(3)]
block = np.ix_(modes, modes)
g0 = grading[block]
assert np.all(majorana[block] == 0)
assert np.array_equal(g0, np.diag([g0[0,0]]*6+[-g0[0,0]]*6))

def exact(matrix):
    return s.Matrix([[s.Rational(float(v.real))+I*s.Rational(float(v.imag)) for v in row] for row in matrix])

def qh(h):
    a,b = h
    return np.array([[np.conj(b),a],[-np.conj(a),b]])

inner0 = inner(np.array([0.,1.]))[block]
up, down = inner0[:3,6:9], inner0[3:6,9:12]
yu, yd = exact(up), exact(down)
a, b = yu*yu.H, yd*yd.H
assert yu.det() != 0 and yd.det() != 0
assert s.expand((a*b-b*a).det()) != 0
singular = np.r_[np.linalg.svd(up, compute_uv=False), np.linalg.svd(down, compute_uv=False)]
Uu, Du = polar(up)
Ud, Dd = polar(down)
zero = np.zeros((3,3))
U = np.block([[Uu,zero],[zero,Ud]])
D = np.block([[Du,zero],[zero,Dd]])
fields = np.array([[1,0],[1j,0],[0,1],[0,1j],[.22+.17j,.4-.10j]], complex)
for h in fields:
    B = np.kron(qh(h), np.eye(3)) @ U @ D
    K = np.block([[np.zeros((6,6)),B],[B.conj().T,np.zeros((6,6))]])
    np.testing.assert_allclose(inner(h)[block], K, atol=2e-14)
    for p5 in (-.75,.25):
        eig = np.sqrt(p5**2 + np.linalg.norm(h)**2*singular**2)
        np.testing.assert_allclose(np.linalg.eigvalsh(p5*g0-1j*g0@K), np.sort(np.r_[-eig,eig]), atol=2e-14)
for unit_h in fields[:4]:
    Q = qh(unit_h)
    np.testing.assert_allclose(Q.conj().T@Q, np.eye(2), atol=1e-14)
    np.testing.assert_allclose(np.linalg.det(Q), 1, atol=1e-14)
print(json.dumps({"checks":["Clifford eigenspaces and four logical matrix units","two-step return","normalized radial integrals","E8 moments through order six and eighth-order contrast","eight glue states and rank six versus eight","E6/A2 distinguished-cycle intersection","pulse and anomaly arithmetic","exact native Yukawa nonsingular commutator","native field factorization and signed auxiliary spectral gap"],"m8":[str(m8_axis),str(m8_diag)],"native_Y_singular_values":singular.tolist(),"cross_intersection_rank":2},indent=2))
