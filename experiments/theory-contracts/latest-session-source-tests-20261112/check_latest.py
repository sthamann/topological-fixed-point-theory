"""Exact checks of three source identifications, not a proposed P1 model."""
from pathlib import Path
import itertools as it
import importlib.util
import json
import hashlib
import numpy as np
import sympy as s

HERE = Path(__file__).resolve().parent
ROOT = Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
SOURCE = ROOT/'experiments/theory-contracts/compiler-correlated-event-clock-20260919/source_channel.py'
spec = importlib.util.spec_from_file_location('source_channel', SOURCE)
original = importlib.util.module_from_spec(spec)
spec.loader.exec_module(original)

# The actual length-four sl2 string of E8_1, derived in PROOF.txt.
# Nonorthonormal basis {v,e v,e^2 v,e^3 v} has the Gram matrix below.
e = s.zeros(4)
for k in range(3): e[k+1,k] = 1
gram = s.diag(1,3,12,36)
f = gram.inv()*e.T*gram
h = s.diag(3,1,-1,-3)
assert f*e-e*f == h
assert h*e-e*h == -2*e
v = s.eye(4)[:,0]
assert ((e**2*v).T*gram*(e**2*v))[0] == 12
assert e**4 == s.zeros(4)

# Independent grade-one check: the original E8 root lattice, norm(alpha)^2=2.
roots = []
for i,j in it.combinations(range(8),2):
    for a,b in it.product([-1,1],repeat=2):
        r = [0]*8; r[i]=a; r[j]=b; roots.append(tuple(r))
for signs in it.product([-1,1],repeat=8):
    if sum(x<0 for x in signs)%2==0:
        roots.append(tuple(s.Rational(x,2) for x in signs))
charges = [r[0]+r[5] for r in roots]
assert len(roots)==240 and sum(c*c for c in charges)==120
counts = {int(c):charges.count(c) for c in set(charges)}
assert counts == {-2:1,-1:56,0:126,1:56,2:1}
grade_one_trace = 12*counts[-2]+4*counts[-1]
assert grade_one_trace == 236
q,VH = s.symbols('q VH')
X = q/(1-q)*(VH-1)-2*q*q/(1-q)**2
A = (4*X+2*VH)/(1-q*q)
assert s.simplify(A-(2*VH/(1-q)**2-4*q/(1-q)**3))==0
E4 = 1+sum(240*s.divisor_sigma(n,3)*q**n for n in range(1,6))
thermal_word = q*s.diff(E4,q)/(E4*(1-q)**2)-4*q/(1-q)**3
assert thermal_word.series(q,0,2).removeO()==236*q

# The duad-syntheme incidence is the original Doily combinatorics.
duads = list(it.combinations(range(6),2))
def matchings(xs):
    if not xs:
        yield ()
        return
    a = xs[0]
    for b in xs[1:]:
        for tail in matchings([x for x in xs if x not in (a,b)]):
            yield tuple(sorted(((a,b),)+tail))
synthemes = list(matchings(list(range(6))))
N = s.Matrix([[int(d in t) for t in synthemes] for d in duads])
B = N.T*N
assert B.eigenvals() == {s.Integer(9):1,s.Integer(4):9,s.Integer(0):5}
P = (13*B-B*B)/36
assert P*P == P and s.trace(P)==10
inner = inner_twist = 0
transposition_char = None
for g5 in it.permutations(range(5)):
    g = g5+(5,)
    perm = [synthemes.index(tuple(sorted(tuple(sorted((g[a],g[b]))) for a,b in t))) for t in synthemes]
    ch_virtual = sum(P[j,perm[j]] for j in range(15))
    fixed = sum(g5[i]==i for i in range(5))
    fixed_square = sum(g5[g5[i]]==i for i in range(5))
    ch_wedge = (fixed*fixed-fixed_square)//2
    ch_standard = fixed-1
    inversions = sum(g5[i]>g5[j] for i in range(5) for j in range(i+1,5))
    inner += ch_virtual*ch_wedge
    inner_twist += (-1)**inversions*ch_virtual*ch_wedge
    if g5 == (1,0,2,3,4): transposition_char = [int(ch_virtual),ch_wedge]
assert inner/120 == 1 and inner_twist == 0
assert transposition_char == [4,2]

# Original C4 matrices, not a fitted four-dimensional representation.
paulis = original.PAULIS
omega = np.eye(4).ravel()/2
Pomega = np.outer(omega,omega)
avg = sum(np.kron(U,U.conj()) for U in paulis)/16
assert np.array_equal(avg,Pomega)
rays = original.source_rays()
reflections = [np.eye(4)-np.outer(z,z.conj())/2 for z in rays]
for U in paulis+reflections:
    V = np.kron(U,U.conj())
    assert np.array_equal(V@omega,omega)
    assert np.array_equal(V@Pomega@V.conj().T,Pomega)
# Exact rational p avoids floating tolerance; denominators powers of two.
for p in [s.Rational(1,16),s.Rational(17,32),s.Integer(1)]:
    rho = p*s.Matrix(Pomega)+(1-p)*(s.eye(16)-s.Matrix(Pomega))/15
    rho = rho.applyfunc(s.nsimplify)
    marginal = s.Matrix(4,4,lambda a,b:sum(rho[4*a+j,4*b+j] for j in range(4)))
    assert marginal == s.eye(4)/4 and s.trace(rho)==1

# Same discrete quarter and half returns do not specify a continuous loop.
windings = [8*(1+4*m) for m in range(3)]
assert windings == [8,40,72]

result = {
 'verdict':'exact scoped source-identification checks passed; physical P1 closure not established',
 'native_E8':{'norm_e2_v_squared':12,'vacuum_e2_zero':'state-specific','gibbs_lower_bound':'12*exp(-beta)/Z_beta','exact_gibbs_word':'q*E4_prime(q)/(E4(q)*(1-q)^2)-4*q/(1-q)^3','grade_one_trace':grade_one_trace,'root_charge_counts':counts,'scope':'canonical E8 lattice vacuum module Gibbs trace, no chemical potential; not identified with P1'},
 'doily':{'rank':10,'spectrum':{'9':1,'4':9,'0':5},'transposition_characters':transposition_char,'S5_Hom_dimension':1,'maximum_intertwiner_rank':4,'one_sided_sign_twist_Hom_dimension':0},
 'conjugate_C4':{'pauli_fixed_dimension':1,'representation':'1+15','invariant_density_parameter':'0 <= p <= 1','both_marginals':'I4/4','native_Pauli_trivial_V5_absent':True},
 'quarter_loop':{'same_discrete_quarter':True,'windings':windings,'scope':'does not assume a primitive continuous source action'},
 'source_channel_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
 'physical_source_derived':False,'Z0_computed':False,'dynamics_uniquely_derived':False}
(HERE/'results.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
