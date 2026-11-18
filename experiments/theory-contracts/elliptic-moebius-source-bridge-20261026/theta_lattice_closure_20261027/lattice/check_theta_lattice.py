"""Exact theta-frame lattice and retained spin-clock check.

Gaussian root coordinates z represent physical vectors z/sqrt(2).
No E8 roots are supplied to the construction; the standard shell is compared
only after the affine-support generators and their lattice are computed.
"""
import itertools
import json
from collections import Counter
from pathlib import Path

import sympy as s
from sympy.matrices.normalforms import hermite_normal_form


def require(test, message):
    if not test:
        raise RuntimeError(message)


def real_tuple(v):
    return tuple(int(t) for z in v for t in (s.re(z), s.im(z)))


def complex_vector(v):
    return s.Matrix([v[2*j] + s.I*v[2*j+1] for j in range(4)])


def integral(M):
    return all(t.is_Integer for t in M)


units = [1, s.I, -1, -s.I]
pi = 1 + s.I
roots = set()
counts = {}
for d in range(3):
    block = set()
    if d == 0:
        for j, u in itertools.product(range(4), units):
            v = [0]*4
            v[j] = u*pi**2
            block.add(real_tuple(v))
    elif d == 1:
        for a, b in itertools.combinations(range(4), 2):
            for u, t in itertools.product(units, repeat=2):
                v = [0]*4
                v[a], v[b] = u*pi, u*t*pi
                block.add(real_tuple(v))
    else:
        for u, a, b, sign in itertools.product(units, units, units, [1, -1]):
            block.add(real_tuple([u, u*a, u*b, u*sign*a*b]))
    counts[str(d)] = len(block)
    roots |= block
require(counts == {'0': 16, '1': 96, '2': 128}, 'support counts')
require(len(roots) == 240, 'root shell cardinality')
require(all(sum(x*x for x in v) == 4 for v in roots), 'root norms')

# Integer basis in z-coordinates; the physical real Gram is B^T B/2.
B = hermite_normal_form(s.Matrix.hstack(*(s.Matrix(v) for v in sorted(roots))))
require(B.shape == (8, 8), 'rank eight')
Gram = B.T*B/2
require(integral(Gram), 'integral Gram')
require(all(Gram[j, j] % 2 == 0 for j in range(8)), 'even Gram')
require(abs(B.det()) == 16 and Gram.det() == 1, 'unimodularity')

# Check the explicit congruence description on Z^8 / 2Z^8.
def admitted(v):
    residues = [(v[2*j] + v[2*j+1]) % 2 for j in range(4)]
    return len(set(residues)) == 1 and sum(v[0::2]) % 2 == 0 and sum(v[1::2]) % 2 == 0

residue_count = sum(admitted(v) for v in itertools.product(range(2), repeat=8))
require(residue_count == 16, 'congruence index')
require(all(admitted(tuple(B[:, j])) for j in range(8)), 'basis congruences')

# Compare actual roots with the original standard E8 shell after an explicit
# orthogonal change of real frame. Outputs below use doubled E8 coordinates.
mapped = {tuple(t for j in range(4) for t in (v[2*j]+v[2*j+1], v[2*j]-v[2*j+1])) for v in roots}
standard = set()
for j, k in itertools.combinations(range(8), 2):
    for a, b in itertools.product([-2, 2], repeat=2):
        v = [0]*8
        v[j], v[k] = a, b
        standard.add(tuple(v))
for v in itertools.product([-1, 1], repeat=8):
    if sum(x < 0 for x in v) % 2 == 0:
        standard.add(v)
require(mapped == standard, 'literal standard E8 root equality')

def real_matrix(M):
    out = s.zeros(8)
    for j in range(4):
        for k in range(4):
            a, b = s.re(M[j, k]), s.im(M[j, k])
            out[2*j, 2*k] = a
            out[2*j, 2*k+1] = -b
            out[2*j+1, 2*k] = b
            out[2*j+1, 2*k+1] = a
    return out

I2, X, Z = s.eye(2), s.Matrix([[0, 1], [1, 0]]), s.diag(1, -1)
S = s.diag(1, s.I)
CNOT = s.zeros(4)
for x, y in itertools.product(range(2), repeat=2):
    CNOT[2*x+(x ^ y), 2*x+y] = 1
# zeta8^{-1}H = (1-i)/2 [[1,1],[1,-1]]; all gates are Gaussian rational.
Htilde = (1-s.I)*s.Matrix([[1, 1], [1, -1]])/2
gates = {'X1':s.kronecker_product(X,I2),'X2':s.kronecker_product(I2,X),
         'Z1':s.kronecker_product(Z,I2),'Z2':s.kronecker_product(I2,Z),
         'S1':s.kronecker_product(S,I2),'S2':s.kronecker_product(I2,S),
         'CNOT':CNOT,'H1_tilde':s.kronecker_product(Htilde,I2),
         'H2_tilde':s.kronecker_product(I2,Htilde)}
for name, g in gates.items():
    A = s.simplify(B.inv()*real_matrix(g)*B)
    require(integral(A) and abs(A.det()) == 1, 'lattice gate ' + name)
    require(A.T*Gram*A == Gram, 'metric gate ' + name)

# Frozen source-compatible representative: sigma=(0 1 2), fixes 3;
# U is the signed permutation in PAULI-SPIN-CENTRAL-LIFT equation (13).
sigma = s.zeros(4)
for j,k in enumerate([1,2,0,3]):
    sigma[k,j] = 1
U = s.Matrix([[1,0,0,0],[0,0,0,-1],[0,0,1,0],[0,1,0,0]])
V = s.kronecker_product(X,I2)*U
require(V*V == s.eye(4), 'rational spin representative square')
require(V*sigma*V == sigma**-1, 'retained cubic conjugation')
for name,g in [('sigma',sigma),('V',V),('J',s.I*s.eye(4))]:
    A = s.simplify(B.inv()*real_matrix(g)*B)
    require(integral(A) and abs(A.det())==1, 'retained rational action ' + name)

zeta8 = (1+s.I)/s.sqrt(2)
alpha = zeta8*V
require(s.simplify(alpha**2) == s.I*s.eye(4), 'alpha square J')
require(s.simplify(alpha**4) == -s.eye(4), 'alpha fourth power')
require(s.simplify(alpha**8) == s.eye(4), 'alpha eighth power')
require(s.simplify(alpha*sigma*alpha.inv()) == sigma**-1, 'alpha conjugation')
Aalpha = s.simplify(B.inv()*real_matrix(alpha)*B)
require(not integral(Aalpha), 'alpha must not preserve this lattice')
# Its two branch transitions are V and iV, both exact lattice maps.
for g in [V, s.I*V]:
    require(integral(s.simplify(B.inv()*real_matrix(g)*B)), 'two-lattice transition')

c = s.I*sigma
require(c**12 == s.eye(4) and all(c**n != s.eye(4) for n in [1,2,3,4,6]), 'C12 order')
require(c**4 == sigma and c**9 == s.I*s.eye(4), 'C12 component identities')

out = {'verdict':'EXACT_E8_LATTICE_AND_TWO_BRANCH_SPIN_ACTION',
       'input':'level-two Schrodinger frame, affine support phases; no E8 input in generator construction',
       'roots_by_support_dimension':counts,'roots':len(roots),'complex_root_lines':len(roots)//4,
       'integer_basis_determinant':int(abs(B.det())), 'physical_gram_determinant':int(Gram.det()),
       'physical_gram':[[int(t) for t in Gram.row(j)] for j in range(8)],
       'congruence_residues_mod_two':residue_count,'standard_e8_shell_equal':True,
       'compensated_clifford_generators':list(gates),'all_gates_integral_isometries':True,
       'retained_alpha_order':8,'retained_alpha_square':'iI',
       'alpha_preserves_single_constructed_lattice':False,
       'alpha_exchanges_Lambda_and_zeta8_Lambda':True,
       'native_cyclic_clock_order':12,'native_clock_identities':['c^4=sigma','c^9=iI'],
       'scope':'abstract integral construction in a theta frame; no physical P1 realization, no uniqueness among every Heisenberg-stable lattice'}
print(json.dumps(out,ensure_ascii=False,sort_keys=True,indent=2))
