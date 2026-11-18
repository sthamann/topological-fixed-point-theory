#!/usr/bin/env python3
"""Exact collective-order test on the existing native four-register source.
No physical Hamiltonian, bath, intervention or additional register is introduced.
"""
from __future__ import annotations
import argparse
import hashlib
import itertools as it
import json
import math
from pathlib import Path
import sys
import numpy as np

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))
from tfpt_explorer.process import _code_basis
from tfpt_explorer.cartan_source import _roots, _rays


def need(condition, message):
    if not bool(condition):
        raise RuntimeError(message)


def sym_scaled(n):
    words = list(it.product(range(4), repeat=n))
    index = {w: i for i, w in enumerate(words)}
    scale = math.factorial(n)
    out = np.zeros((4**n, 4**n), dtype=np.int64)
    groups = {}
    for w in words:
        groups.setdefault(tuple(sorted(w)), []).append(index[w])
    for ids in groups.values():
        need(scale % len(ids) == 0, "Symmetric orbit denominator")
        out[np.ix_(ids, ids)] = scale // len(ids)
    return out


def reduce_four(matrix, keep):
    env = tuple(j for j in range(4) if j not in keep)
    perm = keep + env
    tensor = matrix.reshape((4,)*8).transpose(perm + tuple(j+4 for j in perm))
    d, e = 4**len(keep), 4**len(env)
    return np.einsum("aibi->ab", tensor.reshape(d, e, d, e))


def moment_numerator(rays, n):
    # Actual unit vector is (a_0+i b_0,...,a_3+i b_3)/2.
    # Return Gaussian-integer numerator, denominator 60*4**n.
    words = list(it.product(range(4), repeat=n))
    real = np.zeros((4**n, 4**n), dtype=np.int64)
    imag = np.zeros_like(real)
    for root in rays:
        z = [(root[2*j], root[2*j+1]) for j in range(4)]
        vr, vi = [], []
        for word in words:
            a, b = 1, 0
            for j in word:
                a, b = a*z[j][0]-b*z[j][1], a*z[j][1]+b*z[j][0]
            vr.append(a)
            vi.append(b)
        vr, vi = np.array(vr, dtype=np.int64), np.array(vi, dtype=np.int64)
        real += np.outer(vr, vr) + np.outer(vi, vi)
        imag += np.outer(vi, vr) - np.outer(vr, vi)
    return real, imag


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=Path(__file__).with_name("results.json"))
    args = parser.parse_args()
    _, supports = _code_basis()
    need([len(s) for s in supports] == [4,12,12,12,24], "Literal original supports")
    words = list(it.product(range(4), repeat=4))
    index = {w:i for i,w in enumerate(words)}
    pc24 = np.zeros((256,256), dtype=np.int64)
    for support in supports:
        ids = [index[w] for w in support]
        pc24[np.ix_(ids, ids)] = 24 // len(support)
    ps24 = sym_scaled(4)
    identity = np.eye(256, dtype=np.int64)
    need(np.array_equal(pc24@pc24, 24*pc24), "Code projector")
    need(np.array_equal(ps24@ps24, 24*ps24), "Symmetric projector")
    need(np.array_equal(pc24@ps24, 24*pc24), "Actual code is symmetric")
    need(int(np.trace(pc24)) == 120 and int(np.trace(ps24)) == 840, "Exact ranks 5 and 35")
    # rho_C - rho_S = (7*pc24-ps24)/840.
    difference = 7*pc24-ps24
    reductions = {}
    for n in (1,2,3):
        subsets = list(it.combinations(range(4), n))
        need(all(not np.any(reduce_four(difference, k)) for k in subsets), "All proper marginals equal")
        reductions[str(n)] = {"subsets":len(subsets), "exact_equal":True}
    roots = _roots()
    rays = _rays(roots)
    need(len(roots)==240 and len(rays)==60, "Original roots and native J lines")
    low = {}
    for n in (1,2,3):
        real, imag = moment_numerator(rays, n)
        sym = sym_scaled(n)
        dim = math.comb(n+3,3)
        need(not np.any(imag), "Lower moment real")
        need(np.array_equal(math.factorial(n)*dim*real, 60*4**n*sym), "Exact uniform lower moment")
        low[str(n)] = {"symmetric_dimension":dim, "uniform_exact":True}
    real, imag = moment_numerator(rays, 4)
    need(not np.any(imag), "Fourth moment imaginary part vanishes")
    need(np.array_equal(real, 16*(pc24+ps24)), "40 M4 = Pi_C + Pi_sym4")
    # Recover the original code projector FROM the source moment, no target fit.
    need(np.array_equal(real-16*ps24, 16*pc24), "Literal source reconstruction of code projector")
    h48 = 48*identity-pc24-ps24
    need(not np.any(h48@pc24), "Parent annihilates code")
    need(np.array_equal(h48@(ps24-pc24), 24*(ps24-pc24)), "Symmetric complement energy 1/2")
    need(np.array_equal(h48@(24*identity-ps24), 48*(24*identity-ps24)), "Nonsymmetric energy 1")
    g48 = pc24+ps24
    need(np.array_equal(g48@g48, 24*(3*pc24+ps24)), "Square of algebraic filter")
    mutant = pc24.copy()
    ids = [index[w] for w in supports[0]]
    mutant[np.ix_(ids,ids)] = 0
    mutant[0,0] = 24
    need(np.any(reduce_four(7*mutant-ps24,(0,1,2))), "Wrong-amplitude negative control fires")
    wrong_source = real.copy()
    wrong_source[0,0] += 1
    need(np.any(wrong_source-16*(pc24+ps24)), "Wrong-source moment control fires")
    files = ["tfpt_explorer/process.py", "tfpt_explorer/cartan_source.py",
             "experiments/theory-contracts/joint-native-algebra-selection-20261016/SOURCE-SELECTION.md"]
    result = {
        "verdict":"EXACT_SCOPED_COLLECTIVE_ORDER_FOUR; PHYSICAL_ORIGIN_OPEN",
        "arithmetic":"Gaussian integer tensor products and exact int64 matrix identities; no tolerances",
        "support_sizes":[len(s) for s in supports],
        "code_rank":5, "symmetric_rank":35, "proper_marginals":reductions,
        "root_count":240, "native_complex_lines":60, "lower_moments":low,
        "fourth_moment_identity":"40 M4 = Pi_C + Pi_sym4",
        "source_reconstructs_literal_code_projector":True,
        "parent_spectrum":{"0":5,"1/2":30,"1":221},
        "parent_gap":"1/2", "global_trace_distance":"6/7",
        "no_k_local_exact_ground_selection":"For k<=3, ground(C) implies Sym^4(C^4) is ground as well; arbitrary Hermitian terms, frustration allowed, no ancillas",
        "algebraic_filter":"(20 M4)^n = Pi_C + 2^(-n)(Pi_sym4-Pi_C), n>=1; no physical execution asserted",
        "negative_controls":{"wrong_amplitude":True,"wrong_source_moment":True},
        "source_sha256":{f:hashlib.sha256((REPO/f).read_bytes()).hexdigest() for f in files},
        "not_derived":["Physical choice of E8/J or interaction", "thermalization or physical dissipative law", "3+2 event or clock selection", "many-cell or continuum origin", "48*c3^4 seed coefficient"]
    }
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True)+"\n")
    print(json.dumps({"code_rank":5,"symmetric_rank":35,"proper_marginals_exact":True,"lower_moments_exact":True,"fourth_moment_exact":True,"parent_spectrum":result["parent_spectrum"],"negative_controls":result["negative_controls"]},sort_keys=True))


if __name__ == "__main__":
    main()
