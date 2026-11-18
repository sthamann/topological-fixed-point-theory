"""Exact tests on the original neutral circle, native E8 roots and Q4."""
from collections import Counter
from itertools import combinations, permutations, product
from pathlib import Path
import runpy
from contextlib import redirect_stdout
import io

import sympy as s
from sympy import QQ
from sympy.polys.matrices import DomainMatrix

repo = Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
base = repo / 'docs/research/fundamental-plan'

# Four quarter-conjugates of the original geometric parabolic group.
theta = s.symbols('theta', real=True)
vectors = [s.expand_trig(1-s.cos(theta-k*s.pi/2)) for k in range(4)]
assert s.simplify(vectors[0]+vectors[2]-2) == 0
assert s.simplify(vectors[1]+vectors[3]-2) == 0
bracket = vectors[0]*s.diff(vectors[2], theta)-vectors[2]*s.diff(vectors[0], theta)
assert s.simplify(bracket+2*s.sin(theta)) == 0
x = s.symbols('x')
a, b, c, d, e, f = s.symbols('a b c d e f')
u, v = a+b*x+c*x*x, d+e*x+f*x*x
assert s.Poly(u*s.diff(v,x)-v*s.diff(u,x),x).all_coeffs() == [
    b*f-c*e, 2*(a*f-c*d), a*e-b*d]
print('Original quarter support: pointwise rank 1 everywhere; opposite cuts do not commute.')
print('The commuting Möbius coefficients have all 2x2 minors zero.')

# Ten actual Q4 effect rays in the already verified orthonormal slot coordinates.
slots = list(combinations(range(5),2))
centre = (4+s.sqrt(6))/10
rays = [s.sqrt(s.Rational(2,3))*s.Matrix(
    [int(i in pair)-centre for i in range(5)]) for pair in slots]
enc = s.Matrix.hstack(*(s.kronecker_product(v,v) for v in rays))
K = (enc*enc.T/2).applyfunc(s.expand)
Y = s.diag(*([s.Rational(-1,3)]*3+[s.Rational(1,2)]*2))
charges = [int(6*(Y[i,i]-Y[j,j])/5) for i,j in product(range(5),repeat=2)]
h = s.eye(25)-K
neutral = [i for i,q in enumerate(charges) if q == 0]
charged = [i for i,q in enumerate(charges) if q != 0]
assert Counter(charges) == {-1:6,0:13,1:6}
coupling = h.extract(neutral,charged)
leakage = (coupling*coupling.T).applyfunc(s.expand)
assert DomainMatrix.from_Matrix(leakage).convert_to(QQ.algebraic_field(s.sqrt(6))).rank() == 4
assert s.expand(leakage.trace()-(19912-2432*s.sqrt(6))/140625) == 0
reflect = s.eye(5)-2*s.ones(5)/5
pair_reflect = s.kronecker_product(reflect,reflect)
opposite_h = (pair_reflect*h*pair_reflect).applyfunc(s.expand)
opposite_coupling = opposite_h.extract(neutral,charged)
assert s.expand(s.trace(opposite_coupling*opposite_coupling.T)-
                (19912+2432*s.sqrt(6))/140625) == 0
assert reflect*Y != Y*reflect
pairs = list(product(range(25),range(-1,2)))
totals = [charges[i]+n for i,n in pairs]
lift = s.Matrix(75,75,lambda a,b: h[pairs[a][0],pairs[b][0]] if totals[a]==totals[b] else 0)
E = s.zeros(75,25)
for i,q in enumerate(charges):
    E[pairs.index((i,-q)),i] = 1
assert E.T*E == s.eye(25) and lift*E == E*h
assert lift*s.diag(*totals) == s.diag(*totals)*lift
print('Original Q4 slots: exact 75-dimensional neutral intertwiner and rank-four intermediate-charge term.')

adj = s.Matrix(10,10,lambda i,j:int(not (set(slots[i]) & set(slots[j]))))
ones = s.ones(10)
T = adj/3
P5 = -(adj+2*s.eye(10))*(adj-3*s.eye(10))/6
negative = (adj-s.eye(10))*(adj-3*s.eye(10))/15
absolute = T+s.Rational(4,3)*negative
readout = (enc.T*enc/2).applyfunc(s.expand)
assert readout == s.Rational(4,9)*s.eye(10)+ones/18
assert T**2+P5/3 == readout
assert all(absolute[i,j] == -s.Rational(1,45) for i,j in product(range(10),repeat=2) if adj[i,j])
print('Petersen correction equals original Q4 readout; absolute spectral lift has negative neighbor entries.')

# Extract the native E8 root coordinates used in the existing compiler check.
with redirect_stdout(io.StringIO()):
    old = runpy.run_path(str(base/'check_reflection_origin_eml_20261111.py'))
root_set = {tuple(s.Rational(i,4) for i in raw) for raw,_ in old['roots']}
assert len(root_set) == 240
def dot(v,w):
    return sum(i*j for i,j in zip(v,w))

def subtract(v,w):
    return tuple(i-j for i,j in zip(v,w))

def negate(v):
    return tuple(-i for i in v)

positive = {r for r in root_set if next(i for i in r if i)>0}
sums = {tuple(i+j for i,j in zip(r,t)) for r in positive for t in positive}
simple = sorted(positive-sums)
basis = s.Matrix.hstack(*(s.Matrix(r) for r in simple))
gram = basis.T*basis
assert len(simple)==8 and gram.det()==1
dual = gram.inv()*basis.T
coord = {r:tuple(dual*s.Matrix(r)) for r in root_set}
assert all(basis*s.Matrix(coord[r])==s.Matrix(r) and all(n.is_Integer for n in coord[r]) for r in root_set)
form = s.eye(8)+s.Matrix(8,8,lambda i,j:gram[i,j] if i>j else 0)
assert form+form.T == gram
def epsilon(r,t):
    power = (s.Matrix(1,8,coord[r])*form*s.Matrix(coord[t]))[0]
    return (-1)**power

# This A2 is an explicit choice in the original moving A4, not a selected physical embedding.
first = tuple(s.eye(9)[:,0]-s.eye(9)[:,1])
second = tuple(s.eye(9)[:,1]-s.eye(9)[:,2])
assert first in root_set and second in root_set and dot(first,second)==-1
native_c5 = s.diag(old['C'],s.eye(4))
assert native_c5**5 == s.eye(9)
fixed_roots = [r for r in root_set if native_c5*s.Matrix(r)==s.Matrix(r)]
assert len(fixed_roots)==20
fixed_span = s.Matrix.hstack(*(s.Matrix(r) for r in fixed_roots))
moving_span = s.Matrix.hstack(*(s.eye(9)[:,i]-s.eye(9)[:,i+1] for i in range(4)))
assert fixed_span.rank()==moving_span.rank()==4
assert fixed_span.T*moving_span == s.zeros(20,4)
# An A2 root plane preserved by either order-five map must be fixed pointwise:
# its root automorphism group is dihedral of order 12. The other side is the
# orthogonal A4; the two fixed root spans therefore have zero intersection.
assert native_c5*s.Matrix(first) != s.Matrix(first)
print('No A2 root subsystem is fixed by both original side clocks, or by C_x and H.')

# The two A4 simple bases provide the quarter J(x,y)=(-y,x).
def a4_basis(roots):
    positive = {r for r in roots if next(x for x in r if x)>0}
    sums = {tuple(x+y for x,y in zip(r,t)) for r in positive for t in positive}
    for order in permutations(sorted(positive-sums)):
        trial = s.Matrix.hstack(*(s.Matrix(r) for r in order))
        if trial.T*trial==moving_span.T*moving_span:
            return trial
    raise AssertionError('A4 simple basis missing')
second_a4 = a4_basis(fixed_roots)
a4_dual = (moving_span.T*moving_span).inv()
# The sign of the chosen second basis changes J to -J, with the same mod-2 result.
J = second_a4*a4_dual*moving_span.T-moving_span*a4_dual*second_a4.T
assert {tuple(J*s.Matrix(r)) for r in root_set} == root_set
J_basis = dual*J*basis
assert J_basis**2 == -s.eye(8)
fixed_mod2 = [v for v in product(range(2),repeat=8)
              if all(x%2==0 for x in (J_basis-s.eye(8))*s.Matrix(v))]
assert len(fixed_mod2)==16
for v in fixed_mod2:
    if any(v):
        functional = s.Matrix(1,8,v)*gram
        assert (functional*s.Matrix(v))[0]%4==0
        assert sum(dot(functional,w)%2==0 for w in coord.values())==112
print('Original lattice side quarter: 15 invariant mod-2 character kernels, each with 112 D8 roots.')

weights = sorted(r for r in root_set if (dot(r,first),dot(r,second))==(1,0))
orthogonal = {r for r in root_set if dot(r,first)==dot(r,second)==0}
a2_roots = {first,second,negate(first),negate(second),
            tuple(i+j for i,j in zip(first,second)),tuple(-i-j for i,j in zip(first,second))}
assert a2_roots <= root_set
middle = [subtract(r,first) for r in weights]
low = [subtract(r,second) for r in middle]
mixed = set(weights+middle+low)
assert len(orthogonal)==72 and len(weights)==27 and len(mixed)==81
assert root_set == orthogonal|a2_roots|mixed|{negate(r) for r in mixed}
transport_middle = [epsilon(negate(first),r) for r in weights]
transport_low = [q*epsilon(negate(second),r) for q,r in zip(transport_middle,middle)]
def cubic_sign(i,j,k):
    r,t,w = weights[i],middle[j],low[k]
    if any(a+b+c for a,b,c in zip(r,t,w)):
        return 0
    assert dot(r,t)==-1
    return epsilon(r,t)*transport_middle[j]*transport_low[k]

triads = {}
for indices in combinations(range(27),3):
    signs = {cubic_sign(*p) for p in permutations(indices)}
    assert len(signs)==1
    sign = signs.pop()
    if sign:
        triads[indices] = sign
assert len(triads)==45
assert all(cubic_sign(i,i,k)==cubic_sign(i,k,i)==cubic_sign(k,i,i)==0
           for i,k in product(range(27),repeat=2))
print('Original E8 chosen A2: 72 + 6 + 81 + 81 roots; native-lattice bracket yields 45 symmetric cubic terms.')

# Check equivalence to the actual binary geometry asserted in the attachment.
def q(v):
    a,b = v>>3,v&7
    return ((a&b).bit_count()+(a&1)+(b&1))%2

points = [v for v in range(1,64) if q(v)==0]
binary_lines = {tuple(sorted((i,j,i^j))) for i,j in combinations(points,2) if q(i^j)==0}
assert (len(points),63-len(points),len(binary_lines))==(27,36,45)
def graph(nodes,lines):
    return {i:{j for line in lines if i in line for j in line if j!=i} for i in nodes}
binary_graph, root_graph = graph(points,binary_lines),graph(range(27),triads)
def match(mapping):
    if len(mapping)==27:
        return mapping
    used = set(mapping.values())
    options = {
        i:[j for j in root_graph if j not in used and
           all((m in binary_graph[i])==(n in root_graph[j]) for m,n in mapping.items())]
        for i in binary_graph if i not in mapping
    }
    i = min(options,key=lambda i:len(options[i]))
    for j in options[i]:
        found = match({**mapping,i:j})
        if found is not None:
            return found
    return None

iso = match({})
assert iso is not None
assert {tuple(sorted(iso[i] for i in line)) for line in binary_lines}==set(triads)
p = points[0]
collinear = binary_graph[p]
assert (len(collinear),len(points)-1-len(collinear))==(10,16)
assert Counter(sum(i in collinear for i in line) for line in binary_lines)=={1:40,2:5}
print('All binary Picard lines match all native E8 cubic triads via an explicit isomorphism.')
print('Binary labels -> chosen native roots (indices in sorted highest-weight list):',sorted(iso.items()))

# A three-line against three-line product with equal incidence and opposite sign.
# This remains invariant under arbitrary nonzero rescalings of the weight lines.
seen = {}
witness = None
for lines in combinations(triads,3):
    support = tuple(sorted(i for line in lines for i in line))
    sign = triads[lines[0]]*triads[lines[1]]*triads[lines[2]]
    if support in seen and seen[support][0] != sign:
        witness = (seen[support][1],lines)
        break
    seen[support] = (sign,lines)
assert witness is not None
inverse = {v:k for k,v in iso.items()}
print('Opposite sign products at identical label multiplicities:',
      [([tuple(inverse[i] for i in line) for line in terms],
        [triads[t] for t in terms]) for terms in witness])
