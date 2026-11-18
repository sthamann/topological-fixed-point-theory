"""Same-energy Ramond input shell of the original Fold counter-mode probe."""
from contextlib import redirect_stdout
from collections import defaultdict
from itertools import product
from pathlib import Path
import io
import json
import runpy
import sys

import numpy as np

repo = Path(sys.argv[1]).resolve()
base = repo / "experiments/theory-contracts/chat-triangulation-20261108"
with redirect_stdout(io.StringIO()):
    rectangle = runpy.run_path(str(base / "check_tfpt_higgs_rectangle_countermodes_20261111.py"))
    cover = runpy.run_path(str(base / "check_tfpt_cover_selector_compiler_20261110.py"))

sp, q, Q, vertex, phase, signed = (
    cover[k] for k in ("sp", "q", "Q", "vertex", "phase", "signed")
)
target = tuple(signed["target"])
gauge_rows = signed["rows"]
q2 = tuple(map(int, 2*q))


def shell(alpha, sign):
    target_energy = int(8 * (32 + 2 * sign * alpha))
    def cost(j, x):
        return 8*x*x + int(4*(j + alpha - 4))*x
    minima = [min(cost(j, x) for x in (-1, 1)) for j in range(8)]
    assert min(cost(j, x) + sum(minima) - minima[j]
               for j, x in product(range(8), (-7, 7))) > target_energy
    assert all(abs(4*(j + alpha - 4)) < 16 for j in range(8))
    choices = [[(x, cost(j, x)) for x in (-5, -3, -1, 1, 3, 5)] for j in range(8)]
    found = []

    def visit(xs, accumulated, charges):
        j = len(xs)
        if j == 8:
            extra = target_energy - accumulated
            if tuple(charges) == target and (sum(xs) - sum(q2)) % 4 == 0 and extra % 64 == 0:
                p = sp.Matrix(xs) / 2
                found.append((p, extra // 64))
            return
        for x, energy in choices[j]:
            if accumulated + energy + sum(minima[j+1:]) <= target_energy:
                visit(xs+[x], accumulated+energy,
                      [charge+row[j]*x for charge,row in zip(charges,gauge_rows,strict=True)])

    visit([], 0, [0]*4)
    # Independent meet-in-the-middle check of the exact allowed momenta.
    right = defaultdict(list)
    for xs in product(*[[x for x, _ in choices[j]] for j in range(4, 8)]):
        charge = tuple(sum(row[j+4]*x for j, x in enumerate(xs)) for row in gauge_rows)
        e = sum(cost(j+4, x) for j, x in enumerate(xs))
        right[charge].append((xs, e))
    independently_found = set()
    for xs in product(*[[x for x, _ in choices[j]] for j in range(4)]):
        charge = tuple(t-sum(row[j]*x for j, x in enumerate(xs))
                       for t, row in zip(target, gauge_rows, strict=True))
        e = sum(cost(j, x) for j, x in enumerate(xs))
        for ys, f in right[charge]:
            extra = target_energy-e-f
            if extra >= 0 and extra % 64 == 0 and (sum(xs+ys)-sum(q2)) % 4 == 0:
                independently_found.add((tuple(sp.Matrix(xs+ys)/2), extra//64))
    assert independently_found == {(tuple(p), n) for p, n in found}
    return found


def partitions(n, floor=1):
    if n == 0:
        return [()]
    return [(k,) + part for k in range(floor, n+1) for part in partitions(n-k, k)]


def lowered_columns(candidates, spectator_grade):
    incoming = [(p, word) for p, n in candidates
                for word in partitions(n-spectator_grade)]
    columns = []
    for p, word in incoming:
        image = {}
        for root, sign in ((q, 1), (-q, phase(q,-q))):
            coefficients = vertex({tuple((n,4) for n in word): Q.one}, root, p, -1)
            for monomial, value in coefficients.items():
                image[(tuple(p+root), monomial)] = Q.to_sympy(value)*sign*phase(root,p)
        columns.append(image)
    rows = sorted(set().union(*columns))
    matrix = sp.Matrix(len(rows), len(columns), lambda i,j: columns[j].get(rows[i],0))
    return incoming, matrix


output = {}
for quarter, rows in signed["controls"].items():
    alpha = sp.Rational(quarter)
    output[quarter] = {}
    for sign in (-1, 1):
        momenta = [q+sp.Matrix([sp.sympify(x) for x in row["charge_old"]])
                   for row in rows if row["sign"] == sign]
        candidates = shell(alpha, sign)
        assert all(any(p==m and n==0 for p,n in candidates) for m in momenta)
        odd = [(p,n) for p,n in candidates if q.dot(p)%2==1]
        assert len(candidates)==len(odd)+2
        rank_details = []
        for spectator_grade, multiplicity in ((0,1),(1,7),(2,35)):
            incoming, matrix = lowered_columns(odd, spectator_grade)
            rank_details.append({"spectator_grade": spectator_grade,
                                 "spectator_multiplicity": multiplicity,
                                 "input_q_dimension": len(incoming), "lower_rank": matrix.rank()})
        total_odd_dim = sum(x["spectator_multiplicity"]*x["input_q_dimension"] for x in rank_details)
        total_odd_rank = sum(x["spectator_multiplicity"]*x["lower_rank"] for x in rank_details)
        assert total_odd_dim == sum((1,8,44)[n] for _,n in odd)
        assert total_odd_rank > 0
        output[quarter][str(sign)] = {"full_old_shell_dimension": 2+total_odd_dim,
                                      "extra_q_odd_dimension": total_odd_dim,
                                      "extra_q_odd_lower_mode_rank": total_odd_rank,
                                      "extra_q_odd_rank_details": rank_details,
                                      "momentum_grade_channels": [
            {"momentum": [str(x) for x in p], "oscillator_grade": n,
             "old_cartans_q_and_sum": [str(q.dot(p)), str(sum(p))],
             "original_u_or_v": any(p==m and n==0 for m in momenta)}
            for p,n in candidates
        ]}

def c58_anchor(momentum):
    transfer, word, gamma = (rectangle[k] for k in ("transfer", "word", "gamma"))
    native = rectangle["native"]
    conjugation = native["swap"]
    np.testing.assert_allclose(conjugation @ gamma.conj() @ conjugation.conj().T,
                               -gamma, atol=1e-14)
    widths = rectangle["durations"]
    field = rectangle["width"]
    blocks = []
    for h, d in zip((0, field, field*(1+1j), field*1j, 0), widths, strict=True):
        phi = native["majorana"] + native["inner_form"](np.array([0,h]))
        np.testing.assert_allclose(conjugation @ phi.conj() @ conjugation.conj().T,
                                   phi, atol=1e-14)
        block, _ = transfer(momentum,h,d)
        opposite, _ = transfer(-momentum,h,d)
        np.testing.assert_allclose(conjugation @ block.conj() @ conjugation.conj().T,
                                   opposite, atol=1e-13)
        blocks.append(block)
    empty, _ = transfer(momentum,0,sum(widths))
    p0 = momentum*np.diag(gamma).real<0
    assert p0.sum()==48
    np.testing.assert_allclose(conjugation @ np.diag(p0) @ conjugation.conj().T,
                               np.diag(~p0), atol=1e-14)
    halves = [word(blocks),word(blocks[::-1]),empty]
    occupied_columns = [block[:,p0] for block in halves]
    determinants = [np.linalg.slogdet(y.conj().T @ y) for y in occupied_columns]
    assert all(abs(sign-1)<1e-13 for sign,_ in determinants)
    logs = [float(value) for _,value in determinants]
    independent = [float(2*np.log(np.linalg.svd(y, compute_uv=False)).sum()) for y in occupied_columns]
    np.testing.assert_allclose(logs, independent, atol=3e-13)
    return {"log_reflected_quotients_same_empty": [x-logs[2] for x in logs[:2]],
            "reflected_quotients_same_empty": [float(np.exp(x-logs[2])) for x in logs[:2]],
            "log_order_contrast": logs[0]-logs[1]}


print(json.dumps({
    "old_Ramond_energy_shells": output,
    "original_C58_occupied_boundary_at_finite_momentum": {
        str(p):c58_anchor(p) for p in (.25, .75, -.25, -.75)
    },
    "scope": "Old Ramond energy shells and given native C58 occupied boundary, separately. No common physical packet representation."
}, indent=2))
