"""Read both local q/r orderings through the original neutral CAR currents."""

from contextlib import redirect_stdout
from functools import cache
from itertools import combinations, product
import io
import json
from pathlib import Path
import runpy
import sys

import sympy as s

repo = Path(sys.argv[1]).resolve()
base = repo / "experiments/theory-contracts/chat-triangulation-20261108"
with redirect_stdout(io.StringIO()):
    direct = runpy.run_path(str(base / "check_tfpt_family_middle_ordered_cut_20261110.py"))
    reader = runpy.run_path(str(base / "check_ramond_closed_car_return_20261110.py"))
    cocycle = runpy.run_path(str(repo / "experiments/theory-contracts/"
                                 "canonical-seam-fermionization-20261025/"
                                 "check_same_source_affine_defect.py"))["eps_ext"]

q, r, b = (direct[k] for k in ("q", "r", "d1"))
assert q.dot(b) == -2 and b == r-q
axes = s.Matrix.hstack(direct["V"], q, r)
Q = direct["Q"]
vertex = direct["vertex"]
pairing = direct["pairing"]
old = direct["old"]
kernel = reader["kernel"]
kernel.__globals__["zero_power"] = (0,0,0,0)
initial = q
edges = list(combinations(range(4),2))
order = {"early": ((q,b),(0,-2)), "late": ((b,q),(-2,0))}

def charge_label(p):
    half = int(not p[0].is_Integer)
    return tuple(int(v-s.Rational(half,2)) for v in p), half

def cocycle_phase(charges, adjoints):
    phase = 1+0j
    p = initial
    for a in charges:
        phase *= cocycle(charge_label(a), charge_label(p))
        p += a
    assert p == initial
    for a in adjoints:
        phase /= cocycle(charge_label(-a), charge_label(a))
    return s.Integer(int(phase.real)) + s.I * int(phase.imag)

def radial_coefficient(charges, target):
    powers = [int(charges[i].dot(charges[j])) for i,j in edges]
    base_power = [int(a.dot(initial)) for a in charges]
    for (_,j),d in zip(edges,powers,strict=True):
        base_power[j] += d

    @cache
    def binomials(d, maximum):
        return [Q.from_sympy((-1)**k*s.binomial(d,k)) for k in range(maximum+1)]

    @cache
    def read(normal):
        remaining = [t-a-p for t,a,p in zip(target,base_power,normal,strict=True)]
        assert sum(remaining) == 0
        prefixes = [sum(remaining[:i+1]) for i in range(3)]
        if min(prefixes) < 0:
            return Q.zero
        maxima = [min(prefixes[i:j]) for i,j in edges]
        coefficients = [binomials(d,v) for d,v in zip(powers,maxima,strict=True)]
        # Flows from earlier to later radial insertions; a positive flow
        # raises the early power and lowers the late power.
        value = Q.zero
        for k01 in range(maxima[0]+1):
            for k02 in range(maxima[1]+1):
                k03 = remaining[0]-k01-k02
                if not 0 <= k03 <= maxima[2]:
                    continue
                for k12 in range(maxima[3]+1):
                    k13 = remaining[1]+k01-k12
                    k23 = remaining[2]+k02+k12
                    if not (0 <= k13 <= maxima[4] and 0 <= k23 <= maxima[5]):
                        continue
                    term = Q.one
                    for seq,k in zip(coefficients,(k01,k02,k03,k12,k13,k23),strict=True):
                        term *= seq[k]
                    value += term
        return value
    return read

def neutral_pair(left,right):
    before, target = order[right]
    reverse = order[left][0][::-1]
    charges = before + tuple(-a for a in reverse)
    assert sum(charges,s.zeros(8,1)) == s.zeros(8,1)
    phase = cocycle_phase(charges,reverse)
    read = kernel(axes,axes,charges,radial_coefficient(charges,target+(-2,-2)))
    return phase, read

def first_cut_parity(order_name):
    intermediate = initial + order[order_name][0][0]
    characters = [s.exp(2*s.pi*s.I*x) for x in intermediate]
    assert len(set(characters)) == 1
    return characters[0]

assert first_cut_parity("early") == 1
assert first_cut_parity("late") == -1
forms = {pair:neutral_pair(*pair) for pair in product(order,repeat=2)}
result = {}
for alpha in (s.Rational(1,4), s.Rational(3,4)):
    drift = s.Matrix([j+alpha-4 for j in range(8)])
    displacement = [Q.from_sympy(z) for z in old["cl"](direct["V"].T*drift/8)]
    packets = [old["vertex_zero"](old["series_response"](p,displacement,4))
               for p in direct["polys"]]
    evaluated = {
        "early":[vertex(vertex(v,q,q,1),b,2*q) for v in packets],
        "late":[vertex(vertex(v,b,q),q,r,1) for v in packets],
    }
    row = {}
    for (left,right),(phase,reader_fn) in forms.items():
        neutral = s.Matrix(5,5,lambda i,j:phase*reader_fn(packets[i],packets[j]))
        reference = pairing(evaluated[left],evaluated[right])
        assert neutral == reference, (alpha,left,right,phase,neutral-reference)
        row[left+"_"+right] = {"matched_entries":25, "phase":str(phase),
                               "right_first_cut_character": str(first_cut_parity(right)),
                               "first_two":direct["encode"](neutral[:2,:2])}
    result[str(alpha)] = row
print(json.dumps(result,indent=2))
