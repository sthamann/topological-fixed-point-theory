"""Direct integer-charge CAR endpoint and exact arc support in the actual fold."""
from itertools import combinations
from pathlib import Path
from contextlib import redirect_stdout
import io, runpy
import json,sys
import sympy as s

repo=Path(sys.argv[1])
M=s.Matrix(json.loads((repo/"experiments/theory-contracts/v16-operator-pullback-20261014/fold-results.json")
                       .read_text())["explicit_Weyl_matrix"])
d=s.Matrix([0,1,1,-1,0,0,1,0])
# n = -1, 0 suffice for the four exact ground-mode occupations.
filled=sum(1<<(2*k) for k in range(8))
bits=filled
sign=1
word=[(1,0,True),(2,0,True),(3,-1,False),(6,0,True)]
for k,n,create in reversed(word):
    j=2*k+n+1
    assert bool(bits>>j&1)!=create
    sign*=(-1)**((bits&((1<<j)-1)).bit_count())
    bits^=1<<j
assert sign**2==1
assert s.Matrix([sum(((bits>>(2*k+j))&1)-((filled>>(2*k+j))&1) for j in range(2))
                 for k in range(8)])==d
assert list(M*d)==[0,0,0,0,1,1,-1,1]

# Strip the quarter row phases from the original unitary; they do not affect support.
root=(1+s.I)/s.sqrt(2)
columns=[1,2,6]
rows=[[s.expand_complex(root**((j*k)%8)) for k in columns] for j in range(8)]
all_minors={triple:s.simplify(s.det(s.Matrix([rows[j] for j in triple]))) for triple in combinations(range(8),3)}
support={triple for triple,val in all_minors.items() if val!=0}
counts={j:sum(j in triple for triple in support) for j in range(8)}
norm=s.simplify(sum(s.conjugate(v)*v for v in all_minors.values())/8**3)
assert norm==1 and all(counts.values())
# With annihilation on a different arc, the commutator with each individual arc occupation is nonzero.
assert all(any(j in triple and k not in triple for triple in support for k in range(8)) for j in range(8))
energies={}
for alpha in (s.Rational(1,4),s.Rational(3,4)):
    E=sum((8*n+k+alpha)*(1 if create else -1) for k,n,create in word)
    assert E==14+2*alpha
    energies[str(alpha)]=str(E)
with redirect_stdout(io.StringIO()):
    old = runpy.run_path(str(repo / "experiments/theory-contracts/chat-triangulation-20261108"
                             / "check_tfpt_family_middle_ordered_cut_20261110.py"))
q, r, b, K = (old[key] for key in ("q", "r", "d1", "Q"))
raise_mode, inner = old["vertex"], old["inner"]
ground_b = {(): K.one}
middle = raise_mode(ground_b, q, b)
ground_d = raise_mode(middle, q, r)
# The old coefficient engine omits the fixed lattice cocycle on each step.
assert middle == {((1, 4),): K.one} and ground_d == {(): -2 * K.one}
assert [inner(v, v) for v in (ground_b, middle, ground_d)] == [1, 2, 4]
assert M * q == s.eye(8)[:, 7] - s.eye(8)[:, 6]
for quarter, energy in energies.items():
    alpha = s.Rational(quarter)
    drift = s.Matrix([k + alpha - 4 for k in range(8)])
    assert drift.dot(q) == 0
    assert all(16 + drift.dot(p) == s.Rational(energy) for p in (b, r, d))

print(json.dumps({"original_CAR_ground_norm":1,"H":energies,"nonzero_three_creator_arc_triples":len(support),
                  "all_possible_triples":len(all_minors),"triples_containing_each_arc":counts,
                  "terms_with_four_distinct_arcs":5*len(support),
                  "terms_with_three_distinct_arcs":3*len(support),
                  "family_triplet_squared_norms_before_normalization":[1,2,4]},indent=2))
