"""Use the original R113 Higgs coefficients and labelled Calderon graph."""

import contextlib
from fractions import Fraction
import io
import json
from pathlib import Path
import runpy
import sys

import numpy as np
from scipy.linalg import block_diag
from scipy.special import expit
import sympy as sp

repo = Path(sys.argv[1]).resolve()
base = repo / "experiments/theory-contracts/chat-triangulation-20261108"
sys.path.insert(0, str(repo / "verification"))
entry = "check_r113_dressed_open_20261109.py"
sys.argv = [entry, str(base)]
with contextlib.redirect_stdout(io.StringIO()):
    original = runpy.run_path(str(base / entry))
rep = original["native"]
gamma = original["grading"]
majorana = original["majorana"]
T0 = -1j * gamma @ original["inner_form"](np.array([1, 0]))
T1 = -1j * gamma @ original["inner_form"](np.array([1j, 0]))
T = (T0 - 1j*T1)/2
assert np.count_nonzero(gamma@T + T@gamma) == 0
assert np.count_nonzero(T@T) == 0
n = len(T)
grade = gamma.diagonal().real.astype(int)
root = np.array([[0, 1], [0, 0]])
left = rep.rep(1, root, np.zeros((3,3))) - rep.rep(1, np.zeros((2,2)), np.zeros((3,3)))
weak_root = left - rep.Kswap() @ left.T @ rep.Kswap()
# The conjugate action of an infinitesimal complex root is the negative transpose.
assert np.trace(gamma @ weak_root.conj().T @ weak_root) == 0

h_charge = {"nu": {"L": -3, "R": 0}, "e": {"L": -3, "R": -6},
            "u": {"L": 1, "R": 4}, "d": {"L": 1, "R": -2}}
labels = rep.P
y6_p = np.array([h_charge[f][side] for f, side, _, _ in labels])
y6 = np.diag(np.r_[y6_p, -y6_p])
assert np.count_nonzero(y6@T - T@y6 - 3*T) == 0
weak_weights = {"nu": 1, "e": -1, "u": 1, "d": -1}
weak_p = np.array([weak_weights[f] if side == "L" else 0 for f, side, _, _ in labels])
weak_twice = np.diag(np.r_[weak_p, -weak_p])
assert np.count_nonzero(weak_twice@T-T@weak_twice-T) == 0

def exact(matrix):
    return sp.Matrix([[sp.Rational(float(v.real))+sp.I*sp.Rational(float(v.imag))
                       for v in row] for row in matrix])

def connected_blocks(M):
    edges = np.argwhere(M != 0)
    row_to_cols = {i:set() for i in np.unique(edges[:,0])}
    col_to_rows = {i:set() for i in np.unique(edges[:,1])}
    for i,j in edges:
        row_to_cols[i].add(j)
        col_to_rows[j].add(i)
    unseen = set(row_to_cols)
    while unseen:
        rows, cols = {unseen.pop()}, set()
        while True:
            next_cols = set.union(*(row_to_cols[i] for i in rows))
            next_rows = set.union(*(col_to_rows[j] for j in next_cols))
            if next_cols == cols and next_rows == rows:
                break
            rows, cols = next_rows, next_cols
        unseen.difference_update(rows)
        yield sorted(rows), sorted(cols)

yu, yd, ye, ynu = map(exact, original["qref"].yukawa)
for upper, lower in ((yu, yd), (ynu, ye)):
    a, b = upper@upper.H, lower@lower.H
    assert sp.expand((a@b-b@a).det()) != 0
assert np.count_nonzero(majorana[:rep.NP,rep.NP:]) > 0

ranks = {-1: 0, 1: 0}
for rows, cols in connected_blocks(T):
    assert len(rows) == len(cols) == 3
    assert exact(T[np.ix_(rows,cols)]).det() != 0
    assert len(set(grade[cols])) == 1
    ranks[grade[cols[0]]] += len(cols)
assert ranks == {-1:24, 1:24}
u,s,vh = np.linalg.svd(T, full_matrices=False)
polar = u[:, :48] @ vh[:48, :]
np.testing.assert_allclose(polar@polar.conj().T@polar, polar, atol=2e-15)
assert s[0] - s[47] > 1

def absolute_square(v):
    return Fraction(float(v.real))**2 + Fraction(float(v.imag))**2

signed_norm = sum(grade[j].item()*absolute_square(T[i,j])
                  for i,j in zip(*np.nonzero(T), strict=True))
column_norms = {sum(absolute_square(v) for v in T[:,j])
                for j in np.flatnonzero(np.any(T != 0,axis=0))}
assert len(column_norms) > 1
norm_sq = sum(absolute_square(v) for v in T.ravel())
Z = np.zeros_like(T)
S = block_diag(gamma,-gamma)
diagonal = block_diag(T,Z)
crossed = np.block([[Z,T],[Z,Z]])
root_polar = np.block([[Z,polar],[Z,Z]])
assert np.count_nonzero(S@diagonal-diagonal@S) != 0
assert np.count_nonzero(S@crossed-crossed@S) == 0
assert np.count_nonzero(crossed@crossed) == 0

# Actual weak charge generators on the particle/antiparticle blocks.
# These are a separate test, not an identification of a weak root with r.
q_particle = np.array([f in ("u","d") for f, _, _, _ in labels])
masks = [np.r_[q_particle, np.zeros(rep.NP,bool)],
         np.r_[np.zeros(rep.NP,bool), q_particle], np.tile(~q_particle,2)]
retained_ranks = []
active = np.any(T != 0,axis=0)
for mask in masks:
    assert np.count_nonzero(T[np.ix_(mask,~mask)]) == 0
    assert np.count_nonzero(majorana[np.ix_(mask,~mask)]) == 0
    counts = [int(np.count_nonzero(mask & active & (grade==g))) for g in (-1,1)]
    retained_ranks.append(counts)
assert retained_ranks == [[9,9],[9,9],[6,6]]
expected_weak = [9, -9, 0]
for mask, level in zip(masks, expected_weak, strict=True):
    assert np.trace(gamma @ (weak_root[:,mask] @ weak_root[:,mask].conj().T)) == level

higgs = original["inner_form"](np.array([.22+.17j, .4-.10j]))
V = -1j*gamma@(majorana+higgs)
def graph(p):
    eig, Q = np.linalg.eigh(p*gamma+V)
    c = expit(-4*np.pi*eig)
    x = np.sqrt(expit(-4*np.pi*eig)*expit(4*np.pi*eig))
    C = (Q*c)@Q.conj().T
    X = (Q*x)@Q.conj().T
    return np.block([[np.eye(n)-C,X],[X,C]])

rows = {}
for alpha in (.25,.75):
    graphs = {j:graph(j+alpha) for j in range(-256,258)}
    result = {}
    for cutoff in (16,64,256):
        sums = {}
        for key,A in [("one_face",diagonal),("crossed",crossed)]:
            value = sum(float(np.linalg.norm(graphs[j+1]@A-A@graphs[j])**2)
                        for j in range(-cutoff,cutoff))
            sums[key] = value
        sums["one_face_per_mode"] = sums["one_face"]/(2*cutoff)
        for key,A in [("raw",crossed),("polar",root_polar)]:
            sums[key+"_central_slope"] = float(np.trace(
                (graphs[cutoff]-graphs[-cutoff]) @ (A @ A.conj().T)).real)
        result[str(cutoff)] = sums
    assert abs(result["256"]["raw_central_slope"]+float(signed_norm)) < .01
    assert abs(result["256"]["polar_central_slope"]) < .01
    assert abs(result["256"]["one_face_per_mode"]-float(norm_sq)) < 1
    rows[str(alpha)] = result
print(json.dumps({
    "grade_ranks": ranks,
    "hypercharge_difference": "1/2",
    "weak_cartan_difference": "1/2",
    "distinct_exact_nonzero_column_norms_squared": len(column_norms),
    "one_face_HS_squared_per_mode_limit": float(norm_sq),
    "raw_crossed_cocycle_per_mode": -float(signed_norm),
    "exact_raw_crossed_cocycle_per_mode": str(-signed_norm),
    "polar_crossed_cocycle_per_mode": 0,
    "crossed_rank_by_grade_on_quark_particle_quark_antiparticle_leptons": retained_ranks,
    "nonzero_singular_value_interval": [float(s[47]), float(s[0])],
    "weak_levels_quark_particle_quark_antiparticle_leptons": expected_weak,
    "numerical_controls": rows,
}, indent=2))
