"""Test W covariance for the actual q-root with unequal spinor weights."""
from collections import defaultdict
import json
from pathlib import Path
import runpy
import sys
import numpy as np
import sympy as S
from scipy.sparse import csr_matrix

repo = Path(sys.argv[1])
native = runpy.run_path(str(repo / "experiments/theory-contracts/universalraum-native-operations-ground-response-20260915/common.py"))
W = csr_matrix(native["load_tensor"]())
fermion_weights, boson_weights = native["weights"]()
q = np.array([0,0,0,0,0,0,-2,2])
family_pairs = [(i,j) for i in range(4) for j in range(4)
                if np.array_equal(fermion_weights[i]-fermion_weights[j],q)]
assert len(family_pairs)==1
i,j = family_pairs[0]
equations = defaultdict(lambda: [0]*16)
for s in range(16):
    X = np.zeros((64,64),dtype=np.int64)
    X[4*s+i,4*s+j] = 1
    a = W @ native["exterior_square_derivation"](X)
    d = (8*a - ((a@W.T)@W)).tocoo()
    for row,col,val in zip(d.row,d.col,d.data,strict=True):
        if val:
            equations[int(row),int(col)][s] = int(val)
A = S.Matrix([v for v in equations.values() if any(v)])
assert A.cols==16 and A*S.ones(16,1)==S.zeros(A.rows,1)
basis = S.Matrix([list(v) for v in A.T.columnspace()])
assert basis.shape==(15,16) and basis.nullspace()==[S.ones(16,1)]
N = np.zeros((64,64),dtype=np.int64)
N[4*np.arange(16)+i,4*np.arange(16)+j] = 1
H = N @ N.T - N.T @ N
assert np.array_equal(H@N-N@H,2*N)
assert np.trace(N.T@N)==16
Nb = native["boson_generator"](W.toarray(),N).real.astype(np.int64)
assert np.count_nonzero(Nb)==20
print(json.dumps({"native_family_pair_zero_based":[i,j],"spinor_copies":16,
                  "covariance_constraint_rank":basis.rows,
                  "independent_integer_constraints":[[int(x) for x in row] for row in basis.tolist()],
                  "covariant_diagonal_weights":[1]*16,
                  "trace_N_adjoint_N":int(np.trace(N.T@N)),
                  "root_normalized_chiral_complex_CAR_level_if_promoted":16,
                  "old_E8_root_level":1,"physical_anomaly_computed":False},indent=2))
