"""Compare the pinned native pair tensor with the original E8 root bracket."""
import json
from pathlib import Path
import runpy
import sys

repo = Path(sys.argv[1]).resolve()
sys.path.insert(0,str(repo))
from tfpt_explorer.current_block_geometry import _affine_module

affine = _affine_module()
c = affine.Chevalley(affine.e8_std())
w = runpy.run_path(str(repo / "experiments/theory-contracts/universalraum-native-operations-ground-response-20260915/common.py"))
W = w["load_tensor"]()
fweights,bweights = (list(map(tuple,m.tolist())) for m in w["weights"]())
bindex = {v:a for a,v in enumerate(bweights)}
assert set(fweights + bweights) <= c.ridx.keys()
_, support = w["channel_support"](W)
equations = []
for col,(i,j) in enumerate(w["PAIRS"]):
    a = bindex.get(tuple(x+y for x,y in zip(fweights[i],fweights[j],strict=True)))
    bracket = c.bracket(c.ridx[fweights[i]],c.ridx[fweights[j]])
    if a is None:
        assert not bracket and (i,j) not in support
        continue
    coefficient = bracket.get(c.ridx[bweights[a]],0)
    assert len(bracket)==1 and coefficient in (-1,1)
    assert support[(i,j)] == (a,int(W[a,col]))
    rhs = int(W[a,col] != coefficient)
    equations.append(((1<<i)|(1<<j)|(1<<(64+a)),rhs))

basis={}
for row,rhs in equations:
    while row:
        pivot = row.bit_length()-1
        if pivot not in basis:
            basis[pivot] = (row,rhs)
            break
        row ^= basis[pivot][0]
        rhs ^= basis[pivot][1]
    assert row or rhs==0
solution=0
for pivot,(row,rhs) in sorted(basis.items()):
    if (row & solution).bit_count()%2 != rhs:
        solution |= 1<<pivot
assert len(equations)==480
assert all((row&solution).bit_count()%2==rhs for row,rhs in equations)
print(json.dumps({"columns_checked":len(w["PAIRS"]),"nonzero_brackets":len(equations),
                  "phase_equation_rank":len(basis),"allowed_phase_choices":2**(124-len(basis)),
                  "fermion_signs":[(-1)**(solution>>i&1) for i in range(64)],
                  "boson_signs":[(-1)**(solution>>(64+a)&1) for a in range(60)],
                  "scope":"finite antisymmetric root-bracket tensor, not dynamical CAR fields"},indent=2))
