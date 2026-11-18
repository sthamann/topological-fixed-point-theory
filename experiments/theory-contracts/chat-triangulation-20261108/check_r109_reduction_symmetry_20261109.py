"""Test the actual R109 W bank before invariant environment integration."""
from pathlib import Path
import json
import runpy
import sys

import numpy as np
from scipy.sparse import csr_matrix

repo = Path(sys.argv[1]).resolve()
native = runpy.run_path(str(repo / "experiments/theory-contracts/universalraum-keyD-native-instruments-20260915/native_source.py"))
w = csr_matrix(native["W"])
assert (w @ w.T).toarray().tolist() == (8*np.eye(60,dtype=int)).tolist()
by_sector = {"Spin10":0,"SU4":0}
for i,x in enumerate(native["LIE1"]):
    lift = native["exterior_square"](x)
    image = w @ lift
    eight_b = image @ w.T
    defect = 8*image - eight_b @ w
    defect.eliminate_zeros()
    assert defect.nnz == 0
    anti = eight_b + eight_b.getH()
    anti.eliminate_zeros()
    assert anti.nnz == 0
    by_sector["Spin10" if i<45 else "SU4"] += 1
print(json.dumps({"scope":"same native R109 64 CAR/60 boson W bank; checks Lie action, not a primitive R113 identification",
                  "W_source_sha256":native["PINS"]["W_sha256"],
                  "nonzero_vertices":int(w.nnz),
                  "exact_Lie_equivariance_generators":by_sector,
                  "relation":"W L_X = B_X W, B_X = W L_X W^*/8, B_X^* = -B_X",
                  "invariant_environment_reduction_preserves_Spin10":True,
                  "source_context":"consequence for full invariant trace or invariant rest cap; actual R113 rest cap not assigned"},indent=2))
