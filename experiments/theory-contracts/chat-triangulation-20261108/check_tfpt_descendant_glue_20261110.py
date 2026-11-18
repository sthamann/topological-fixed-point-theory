"""Exact joint q,r sectors in the original E8 lattice and folded coordinates."""
from collections import Counter, defaultdict
from itertools import combinations, product
import json
from pathlib import Path
import sys
import sympy as S

repo = Path(sys.argv[1])
M = S.Matrix(json.loads((repo / "experiments/theory-contracts/v16-operator-pullback-20261014/fold-results.json").read_text())["explicit_Weyl_matrix"])
unit = S.eye(8)
q = unit[:,7]-unit[:,6]
r = unit[:,4]+unit[:,5]
roots = [a*unit[:,i]+b*unit[:,j] for i,j in combinations(range(8),2) for a,b in product((-1,1), repeat=2)]
roots += [S.Matrix(p)/2 for p in product((-1,1),repeat=8) if p.count(-1)%2==0]
key = lambda v: tuple(v)
root_set = {key(v) for v in roots}
assert {key(M*v) for v in roots}==root_set
assert q.dot(r)==0 and q.dot(q)==r.dot(r)==2 and key(q) in root_set and key(r) in root_set
assert {key(M.T*q),key(M.T*r)} <= root_set
assert key(M.T*(r-q)) == (1,0,0,0,1,-1,0,1)
transverse = [v for v in roots if v.dot(q)==v.dot(r)==0]
assert len(transverse)==60
# The simple roots from an explicitly regular chamber have a D6 Cartan matrix.
chamber = S.Matrix([100**i for i in range(8)])
positive = [v for v in transverse if chamber.dot(v)>0]
positive_set = {key(v) for v in positive}
simple = [v for v in positive if all(key(v-w) not in positive_set for w in positive)]
assert len(simple)==6
G = S.Matrix.hstack(*simple).T*S.Matrix.hstack(*simple)
assert sorted(sum(G[i,j] == -1 for j in range(6)) for i in range(6))==[1,1,1,2,2,3]
assert G.det()==4
basis = S.Matrix.hstack(q,r,*simple)
assert abs(basis.det())==4
# Four cosets of the direct-sum lattice inside the unit-covolume E8 lattice.
# Pairing with q,r modulo 2 distinguishes all four: use representatives
# that are themselves E8 roots. Each group is a full quotient class.
projections = defaultdict(set)
counts = Counter()
reps = {}
for v in roots:
    a,b = v.dot(q),v.dot(r)
    label = f"{int(a)%2}{int(b)%2}"
    counts[label]+=1
    if label!="00":
        c = v-(a*q+b*r)/2
        projections[label].add(key(c))
        reps[label]=v
assert counts=={"00":64,"01":64,"10":64,"11":48}
details={}
for label in ("01","10","11"):
    comps=projections[label]
    norm2={S.Matrix(c).dot(S.Matrix(c)) for c in comps}
    assert len(norm2)==1
    details[label]={"root_count":counts[label],"transverse_minimal_vectors":len(comps),
                    "transverse_squared_length":str(norm2.pop())}
    # Verify that changing the representative within the residue class
    # changes the complement only by its integral root lattice.
    assert all(all(t.is_Integer for t in basis.inv()*(S.Matrix(c)-S.Matrix(next(iter(comps))))) for c in comps)
assert [details[k]["transverse_minimal_vectors"] for k in ("01","10","11")]==[32,32,12]
assert [details[k]["transverse_squared_length"] for k in ("01","10","11")]==["3/2","3/2","1"]
assert all((p.dot(q)%2,p.dot(r)%2)==(0,0) for p in (q,r-q,r,r+q,2*q))
# All neutral and charged P1 words under study shift momenta by integers
# in the q,r root sublattice. Hence they preserve these two residues.
assert all(all(t.is_Integer for t in basis.inv()*p) for p in (q,r-q,r,r+q,2*q))
print(json.dumps({"direct_sublattice":"A1(q) + A1(r) + D6", "index":4,
                  "root_count":dict(sorted(counts.items())),
                  "nontrivial_residue_sectors":details,"old_P1_residue":"00",
                  "physical_gluing_computed":False},indent=2))
