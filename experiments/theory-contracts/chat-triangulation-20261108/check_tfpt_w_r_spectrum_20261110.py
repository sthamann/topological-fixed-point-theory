"""Reconstruct original W dynamics in the r and 2r Cartan sectors."""
from fractions import Fraction
from itertools import combinations
import json
from math import factorial
from pathlib import Path
import runpy
import sys

import numpy as np
import sympy as sp

repo = Path(sys.argv[1]).resolve()
contracts = repo / "experiments/theory-contracts"
common = runpy.run_path(str(contracts / "universalraum-native-operations-ground-response-20260915/common.py"))
tensor = common["load_tensor"]()
support, lookup = common["channel_support"](tensor)
fermion_weights, boson_weights = common["weights"]()
r = np.array([0, 0, 0, 0, 2, 2, 0, 0])
q = np.array([0, 0, 0, 0, 0, 0, -2, 2])
assert r @ r == 8 and r @ q == 0
row, = np.flatnonzero(np.all(boson_weights == r, axis=1))
modes = np.flatnonzero(fermion_weights @ r == 4).tolist()
pairs = support[row]
assert len(modes) == 16 and len(pairs) == 8
assert sorted(j for ij, _ in pairs for j in ij) == modes
assert set(np.unique(fermion_weights @ r)) == {-4, 0, 4}
assert np.flatnonzero(boson_weights @ r == 8).tolist() == [int(row)]
assert np.max(boson_weights @ r) == 8
for b, terms in enumerate(support):
    for (i, j), coefficient in terms:
        assert abs(coefficient) == 1
        assert np.array_equal(fermion_weights[i] + fermion_weights[j], boson_weights[b])

g, d, E = sp.symbols("g Delta E", real=True)
ann, create = common["ann"], common["create"]

def original_sector(k):
    # Unnormalized b occupation basis, squared norm b!. The r-dot bound
    # certifies that no other boson row or fermion mode can occur.
    labels = []
    for b in range(k, -1, -1):
        for slots in combinations(modes, 2 * (k - b)):
            w = np.sum(fermion_weights[list(slots)], axis=0)
            if np.array_equal(w, (k - b) * r):
                labels.append((b, sum(1 << j for j in slots)))
    index = {key: i for i, key in enumerate(labels)}
    interaction = np.zeros((len(labels), len(labels)), dtype=int)
    for j, (b, mask) in enumerate(labels):
        for (a, c), w in pairs:
            if (u := ann(mask, a)) and (v := ann(u[0], c)):
                interaction[index[b + 1, v[0]], j] += w * u[1] * v[1]
            if b and (u := create(mask, c)) and (v := create(u[0], a)):
                interaction[index[b - 1, v[0]], j] += b * w * u[1] * v[1]
    mass = np.diag([factorial(b) for b, _ in labels])
    assert np.array_equal(mass @ interaction, interaction.T @ mass)
    return labels, mass, interaction

labels2, mass2, inter2 = original_sector(1)
labels4, mass4, inter4 = original_sector(2)

def bright_basis(labels, k):
    index = {key: i for i, key in enumerate(labels)}
    bright = np.zeros((len(labels), k+1), dtype=int)
    bright[index[(k, 0)], 0] = 1
    for order in range(1,k+1):
        for group in combinations(pairs, order):
            mask, sign = 0, 1
            for (i, j), w in group:
                mask, s = create(mask, j)
                mask, t = create(mask, i)
                sign *= w*s*t
            bright[index[(k-order, mask)], order] = sign
    return bright

F2, F4 = bright_basis(labels2, 1), bright_basis(labels4, 2)
A2 = np.array([[0,8],[1,0]])
A4 = np.array([[0,8,0],[2,0,7],[0,2,0]])
assert np.array_equal(inter2 @ F2, F2 @ A2)
assert np.array_equal(inter4 @ F4, F4 @ A4)
assert (F2.T @ mass2 @ F2).diagonal().tolist() == [1,8]
assert (F4.T @ mass4 @ F4).diagonal().tolist() == [2,8,28]
assert len(labels2) == 9
# Number of Q4 states that contain zero, one or two complete pair channels.
occupancy_census = {}
for b, mask in labels4:
    filled = sum((mask >> i & 1) and (mask >> j & 1) for (i,j), _ in pairs)
    occupancy_census[b, filled] = occupancy_census.get((b,filled),0) + 1
    if b==0 and filled==0:
        assert not np.any(inter4[:,labels4.index((b,mask))])
assert occupancy_census == {(2,0):1,(1,1):8,(0,2):28,(0,0):occupancy_census.get((0,0),0)}
matched_two = [i for i,(b,mask) in enumerate(labels4)
               if b == 0 and sum((mask >> a & 1) and (mask >> c & 1) for (a,c),_ in pairs) == 2]
matched_one = [i for i,(b,_) in enumerate(labels4) if b == 1]
incidence = inter4[np.ix_(matched_one, matched_two)]
gauge_one = F4[matched_one, 1]
gauge_two = F4[matched_two, 2]
assert len(matched_two) == 28 and len(matched_one) == 8
incidence = gauge_one[:,None]*incidence*gauge_two[None,:]
assert set(incidence.ravel()) == {0,1}
assert np.array_equal(incidence @ incidence.T, 6*np.eye(8,dtype=int)+np.ones((8,8),dtype=int))
assert occupancy_census[(0,0)] == 12
H2 = sp.diag(d,0)+g*sp.Matrix(A2)
H4 = sp.diag(2*d,d,0)+g*sp.Matrix(A4)
e1_poly = E**2-d*E-8*g**2
assert sp.expand(sp.rem((H4-2*E*sp.eye(3)).det(),e1_poly,E)-4*g**2*(d-E)) == 0

# Check the full occupation matrix, not just the analytic bright subspace.
Delta, coupling = 1.,1/20
energy_plus = (Delta+np.sqrt(Delta**2+32*coupling**2))/2
norms = np.sqrt(mass4.diagonal())
raw = Delta*np.diag([b for b,_ in labels4])+coupling*inter4
full = norms[:,None]*raw/norms[None,:]
spectrum = np.linalg.eigvalsh(full)
top_bright = np.max(np.linalg.eigvalsh(np.array(
    [[2*Delta,4*coupling,0],[4*coupling,Delta,np.sqrt(14)*coupling],
     [0,np.sqrt(14)*coupling,0]])))
assert abs(spectrum[-1] - top_bright) < 1e-13
assert spectrum[-1] < 2*energy_plus

# Positive exact-energy Q2 state has no twice-applied, same-frequency image.
# Naked positive-energy dressing (linear combination of P^* and b^*) still
# has a strictly positive Q4 component before spectral filtering.
a, b = sp.symbols("a b", nonnegative=True)
unfiltered_second_norm2 = 112*a**4+32*a**2*b**2+2*b**4
assert sp.expand(2*(8*a**2+b**2)**2-16*a**4-unfiltered_second_norm2) == 0

print(json.dumps({
    "q_root_doubled":q.tolist(),"r_root_doubled":r.tolist(),"native_r_boson_row":int(row),
    "original_r_pair_modes":pairs,"q2_dimension":len(labels2),"q4_dimension":len(labels4),
    "q4_census_by_boson_number_and_complete_pairs":{str(k):v for k,v in occupancy_census.items()},
    "q2_upper_energy":energy_plus,
    "q4_spectral_decomposition":"3 bright eigenvalues + 7 copies of each root of E^2-Delta*E-6*g^2 + 32 zeros","q4_full_maximum":spectrum[-1],
    "twice_q2_upper_energy_minus_q4_max":2*energy_plus-spectrum[-1],
    "symbolic_det_H4bright_minus_2E_on_shell":"4*g^2*(Delta-E)",
    "symbolic_low_order_Q4_norm2":str(unfiltered_second_norm2),
    "physical_weight":None,"physical_APS_mapping":None,
},indent=2))

# Read the old level-one vertex in the same published fold chart.
from contextlib import redirect_stdout
import io
with redirect_stdout(io.StringIO()):
    old = runpy.run_path(str(contracts/"chat-triangulation-20261108/check_tfpt_family_middle_ordered_cut_20261110.py"))
s, K, old_r, old_q, vertex = (old[name] for name in ("sp","Q","r","q","vertex"))
assert tuple(int(v) for v in 2*old["M"]*old_r) == tuple(r)
assert tuple(int(v) for v in 2*old["M"]*old_q) == tuple(q)
vacuum = {(): K.one}
zero = s.zeros(8,1)
first = vertex(vacuum,old_r,zero,1)
second_equal_mode = vertex(first,old_r,old_r,1)
second_third_mode = vertex(first,old_r,old_r,3)
assert first == vacuum and not second_equal_mode
assert s.simplify(old["inner"](second_third_mode,second_third_mode)) == 1
assert (2*old_r).dot(2*old_r)/2 - 2*old_r.dot(old_r)/2 == 2
print("Old folded current: first norm² 1, repeated -1 norm² 0, -3 after -1 norm² 1.")
old_path = vacuum
for k in range(1,9):
    old_path = vertex(old_path,old_r,(k-1)*old_r,2*k-1)
    assert s.simplify(old["inner"](old_path,old_path)) == 1
    assert (k*old_r).dot(k*old_r)/2 == k*k
print("The old charge-k*r word is nonzero and has L0=k² through k=8 (exact).")
# Work in the occupation basis (b^*)^n instead of normalized b states.
# Fermion signs, boson creation and integer annihilation coefficients give
# integral matrices for 20H on every finite Q block when Delta=1 and g=1/20.
for alpha, trace_value in [(Fraction(1,4),Fraction(1160,13)),(Fraction(3,4),Fraction(248,3))]:
    first_old = 8 + 2*alpha - 2
    second_old = 32 + 2*(2*alpha-2)
    calibration_ratio = second_old/first_old
    assert 20*calibration_ratio == trace_value
    assert trace_value.denominator != 1
    print(f"Quarter {alpha}: second/first={calibration_ratio}, trace of scaled target={trace_value}, not integral.")
