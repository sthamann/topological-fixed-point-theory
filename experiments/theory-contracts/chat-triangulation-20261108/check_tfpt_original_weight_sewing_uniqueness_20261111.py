"""Check two positive time sewings of the same original Higgs sign fibre."""
import json
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

import numpy as np

repo = Path(sys.argv[1]).resolve()
package = repo / (
    "experiments/theory-contracts/original-paper-source-reconciliation-20261102"
    "/contribution_20261106Q"
)
sys.path.insert(0, str(package / "runtime"))
from assemble import Assembler

native_path = package / "domain/native-diagnostic-v1.json"
coarse = Assembler(native_path)
background = {"genes": [2, 2, 2]}
old = coarse.assemble(background)["h_N"]
domain = json.loads(native_path.read_text())
domain["modes"] = list(range(-7, 7))
domain["field_proposal"]["slots"].append(
    {"kind": "cos", "harmonic": 8, "component": 0, "unit": "13/100"}
)
with TemporaryDirectory() as directory:
    path = Path(directory) / "same-q-refinement.json"
    path.write_text(json.dumps(domain))
    fine = Assembler(path)
    fields = {s: fine.assemble({"genes": background["genes"]+[s+1]}) for s in (-1, 1)}

modes = domain["modes"]
indices = np.concatenate([96*modes.index(n) + np.arange(96) for n in coarse.domain["modes"]])
for field in fields.values():
    assert np.array_equal(field["h_N"][np.ix_(indices, indices)], old)


def real_derivative(matrix):
    return (matrix - fine.K @ matrix @ fine.K)/1j


rep = fine.original.rep
six_y = real_derivative(rep(3j, np.zeros((2, 2)), -1j*np.eye(3)))
two_t3 = real_derivative(rep(0, 1j*np.diag([1, -1]), np.zeros((3, 3))))
electromagnetic = np.exp(1j*np.pi*np.diag(six_y/6 + two_t3/2).real)
half_turn = np.concatenate([(-1.0)**n*electromagnetic for n in modes])
plus, minus = (fields[s]["h_N"] for s in (1, -1))
np.testing.assert_allclose(half_turn[:, None]*plus*half_turn.conj()[None, :], minus,
                           rtol=0, atol=1e-13)

amp = F(13, 100)**2
for sign, field in fields.items():
    h = {n: v[[0, 2]]+1j*v[[1, 3]] for n, v in field["higgs_fourier"].items()}
    marks = [2*sum(np.vdot(v, h[n+k]) for n, v in h.items() if n+k in h) for k in (7, 16)]
    np.testing.assert_allclose(marks, [float(sign*amp), float(amp/2)], rtol=0, atol=1e-15)

signs = (-1, 1)


def transition(a, b, correlation):
    return (1+a*b*correlation)/2


def path_mass(word, correlation):
    weight = F(1, 2)
    for a, b in zip(word, word[1:]):
        weight *= transition(a, b, correlation)
    return weight


def check(correlation, half):
    assert 0 < correlation < 1 and half**2 == correlation
    for a, c in product(signs, repeat=2):
        assert sum(transition(a,b,half)*transition(b,c,half) for b in signs) == transition(a,c,correlation)
        assert sum(path_mass((a,b,c), correlation) for b in signs) == path_mass((a,c),correlation**2)
    words = list(product(signs, repeat=4))
    distribution = {w: path_mass(w, correlation) for w in words}
    assert sum(distribution.values()) == 1 and all(p > 0 for p in distribution.values())
    assert all(distribution[w] == distribution[w[::-1]] == distribution[tuple(-s for s in w)] for w in words)
    assert all(sum(p for w,p in distribution.items() if w[i]==1) == F(1,2) for i in range(4))
    response = [sum(amp**2*w[0]*w[i]*p for w,p in distribution.items()) for i in (1, 2)]
    assert response == [amp**2*correlation, amp**2*correlation**2]

    # An arbitrary complex two-site half-word basis is spanned by its
    # four sign indicators. The bond reflection matrix is an exact sum
    # of two positive rank-one forms, independently checked entrywise.
    halfwords = list(product(signs, repeat=2))
    for left, right in product(halfwords, repeat=2):
        lp, rp = transition(*left,correlation), transition(*right,correlation)
        gram = distribution[left[::-1]+right]
        constant_part = lp*rp/4
        odd_part = correlation*(left[0]*lp)*(right[0]*rp)/4
        assert gram == constant_part + odd_part

    # Conditioning on the complete sign at the middle reproduces the
    # reflected two-step moment. Keeping only the common coarse label
    # makes the two conditional means vanish.
    complete = sum((amp*correlation*s)**2/2 for s in signs)
    collapsed = (sum(amp*correlation*s/2 for s in signs))**2
    assert complete == response[1] and collapsed == 0
    return {
        "odd_transfer_eigenvalue_at_tau": str(correlation),
        "same_sign_probability": str(transition(1,1,correlation)),
        "one_time_U7_mean": "0", "one_time_U16_mean": str(amp/2),
        "U7_joint_moment_tau_and_2tau": [str(x) for x in response],
        "U7_joint_moment_tau_and_2tau_decimal": [float(x) for x in response],
        "collapsed_cut_U7_joint_moment_2tau": str(collapsed),
    }


print(json.dumps({
    "original_refined_signs_share_coarse_matrix_and_full_single_slice_spectrum": True,
    "sign_marks_U7": [str(-amp), str(amp)],
    "A": check(F(1,4), F(1,2)),
    "B": check(F(1,16), F(1,4)),
    "scope": "Two projective, reflection-positive time extensions of the checked scalar original sign fibre. No identification with full P_adm or the original joint time action.",
}, indent=2))
