"""Dressed open R032 operator and original S/T reference check."""
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import expm, expm_frechet

base = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("quarter", base / "check_quarter_reference.py")
qref = importlib.util.module_from_spec(spec)
spec.loader.exec_module(qref)
native = qref.native
grading = qref.grading
majorana = qref.mass_matrix
yukawa = qref.higgs_matrix
swap = native.Kswap()
right_projection = native.rep(1, np.zeros((2, 2)), np.zeros((3, 3)))
identity = np.eye(len(grading))


def quaternion(h):
    a, b = h
    return np.array([[b.conjugate(), a], [-a.conjugate(), b]])


def real_action(g):
    left = native.rep(1, g, np.eye(3))
    return left @ swap @ left.conj() @ swap


def inner_form(h):
    rep = native.rep(0, quaternion(h), np.zeros((3, 3)))
    part = rep @ (yukawa @ right_projection - right_projection @ yukawa)
    part = part + part.conj().T
    return part + swap @ part.conj() @ swap


def ordered(blocks, links, jets):
    transfer, response = blocks[0], jets[0]
    for block, link, jet in zip(blocks[1:], links, jets[1:]):
        response = jet @ link @ transfer + block @ link @ response
        transfer = block @ link @ transfer
    return transfer, response


def maximum(a, b):
    return float(np.max(np.abs(a-b)))


p = .25
half_higgs = [(0.22+0.17j, 0.4-0.10j), (-0.23j, .18+0.12j), (.31-.27j, .16j)]
history = half_higgs + half_higgs[::-1]
durations = [.09, .14, .17, .17, .14, .09]
norms = [np.linalg.norm(h) for h in history]
frames = [real_action(quaternion(h)/r) for h, r in zip(history, norms)]
links = [identity] * (len(history)-1)
dressed_links = [v.conj().T @ u for u, v in zip(frames, frames[1:])]
assert np.linalg.norm(dressed_links[0]@dressed_links[1]-dressed_links[1]@dressed_links[0]) > .1
assert max(maximum(W.conj().T@W, identity) for W in frames) < 1e-14
assert max(maximum(W.conj().T@majorana@W, majorana) for W in frames) < 1e-14


def wing(branch, dressed, transport):
    blocks, jets, sea = [], [], 0.
    for i, (h, r, width) in enumerate(zip(history, norms, durations)):
        coord = (0j, complex(r)) if dressed else h
        linear = inner_form(coord) if branch != "S" else np.zeros_like(yukawa)
        bare = majorana if branch != "S" else np.zeros_like(majorana)
        mark = -1j*grading@linear
        generator = p*grading - 1j*grading@bare + mark
        if i in (2, 3) and branch == "A":
            block, jet = expm_frechet(-width*generator, -width*mark)
        else:
            block, jet = expm(-width*generator), np.zeros_like(generator)
        blocks.append(block)
        jets.append(jet)
        sea += width*np.linalg.eigvalsh(generator)[:len(generator)//2].sum()
    transfer, derivative = ordered(blocks, transport, jets)
    return transfer, derivative, sea


old = {j: wing(j, False, links) for j in "AST"}
new = {j: wing(j, True, dressed_links) for j in "AST"}
residuals = {}
for j in old:
    matrix, mark, sea = old[j]
    stripped, dmark, dsea = new[j]
    residuals[j] = [maximum(frames[-1]@stripped@frames[0].conj().T, matrix),
                    maximum(frames[-1]@dmark@frames[0].conj().T, mark),
                    abs(sea-dsea)]
assert max(v for row in residuals.values() for v in row) < 1e-11
assert maximum(old["A"][0], old["T"][0]) < 1e-12
assert abs(old["A"][2]-old["T"][2]) < 1e-12

# Deliberately omit the transformed connections, retaining the same U path.
wrong = wing("A", True, links)
bad = maximum(frames[-1]@wrong[0]@frames[0].conj().T, old["A"][0])
assert bad > .1
assert maximum(old["S"][0], new["S"][0]) < 1e-12
# Recheck the new smooth static IR estimator on the recorded one-domain
# R032 vacuum measure. No full admissible selector is supplied by this branch.
archive = json.loads((base / "r113-temporal-cross.json").read_text())
weight = archive["densities"]["squared_charged_weight"]
radius = 2**18
energy = (np.arange(radius) + .5) / radius
z = .004
slope = weight*np.sum(z/(energy+z)**2)/radius
slope_limit_at_z = weight/(1+z)
assert abs(slope - slope_limit_at_z) < 2e-5
print(json.dumps(dict(full96_matrix_mark_sea_errors=residuals,
                      drop_dressed_links_error=bad,
                      three_norms_squared=[r*r for r in norms[:3]],
                      r032_vacuum_IR_slope_at_z=float(slope),
                      r032_vacuum_IR_limit_at_fixed_z=float(slope_limit_at_z),
                      open_boson_integral_evaluated=False), indent=2))
