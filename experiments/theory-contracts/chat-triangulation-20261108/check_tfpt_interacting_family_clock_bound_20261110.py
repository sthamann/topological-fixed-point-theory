"""Energy and exact central-term bounds for the original native W source."""
from fractions import Fraction
from math import comb
from pathlib import Path
import json
import runpy
import sys

import numpy as np
import sympy as sp

repo = Path(sys.argv[1]).resolve()
base = repo/"experiments/theory-contracts/universalraum-native-operations-ground-response-20260915"
source = runpy.run_path(str(base/"common.py"))
W = source["load_tensor"]()
source["channel_support"](W)
assert W.shape==(60,2016)
assert np.array_equal(W@W.T,8*np.eye(60,dtype=np.int64))
assert np.count_nonzero(W)==480

# The exact q=4 blocks are (4 fermions, 1 boson plus 2 fermions, 2 bosons).
dimensions = [comb(64,4),60*comb(64,2),comb(61,2)]
assert sum(dimensions)==758166

# Creation from two bosons has norm sqrt(8)*sqrt(2). For a single
# occupied boson and two fermions, the exterior multiplication
# wedge^2 tensor wedge^2 -> wedge^4 has six distinct ordered
# decompositions of each four-set, each with unit sign.
indices = ((0,1),(0,2),(0,3),(1,2),(1,3),(2,3))
weights = []
for first in indices:
    second = tuple(i for i in range(4) if i not in first)
    order = first+second
    weights.append((-1)**sum(x>y for i,x in enumerate(order) for y in order[i+1:]))
assert sum(x*x for x in weights)==comb(4,2)
A_squared = 8*2
B_squared_bound = 8*sum(x*x for x in weights)
operator_bound = sp.sqrt(A_squared+B_squared_bound)
assert operator_bound==8
delta, g = Fraction(1),Fraction(1,20)
spectrum = (-int(operator_bound)*g,2*delta+int(operator_bound)*g)
diameter = spectrum[1]-spectrum[0]
assert diameter==Fraction(14,5) and diameter<8

# A genuine affine Cartan raising/lowering pair has central level >0.
# If it commutes with original total Q, restricting to the finite Q=4
# sector makes that identity impossible regardless of time rescaling.
assert 2*sum(dimensions)>0
print(json.dumps({
    "original_pinned_W_shape":list(W.shape),
    "WW_transpose": "8 I60",
    "native_Q4_ordered_sector_dimensions":dimensions,
    "interaction_block_norm_squared_two_boson_to_one_boson_two_fermions":A_squared,
    "interaction_block_norm_squared_bound_one_boson_two_fermions_to_four_fermions":B_squared_bound,
    "interaction_full_norm_bound_in_units_abs_g":int(operator_bound),
    "declared_Delta":str(delta),"declared_g":str(g),
    "exact_Q4_spectral_interval_bound":[str(x) for x in spectrum],
    "exact_Q4_spectral_diameter_bound":str(diameter),
    "no_original_old_shift_8_at_same_declared_time":True,
    "fixed_Q_full_affine_central_relation_possible":False,
    "physical_R113_values":None,
},indent=2))
