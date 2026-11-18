"""Construct source-compatible crossed-lattice lifts; no fitted H, state or time.

The first search is deliberately finite. Its algebraic solutions are conditional
defect extensions, never a selection of the original physical seam from P1/P2.
"""
from __future__ import annotations
import itertools
import json
from pathlib import Path
import numpy as np

SCHEMA = "tfpt.crossed_defect.v1"
ENUMERATION_SIZE = 4**9
SOURCE_FILES = ["domain_spec.json"]
PERM = [4, 5, 6, 7, 0, 1, 2, 3]
PAIRS = list(itertools.combinations(range(8), 2))
_CACHE = {}


def canonicalize(g):
    # X -> exp(i theta)X is a choice of section. U_s carries its inverse phase,
    # so zeta is gauge here and cannot generate novel physical candidates.
    return {"schema": SCHEMA, "chi_exp": [int(x) % 4 for x in g["chi_exp"]],
            "zeta_exp": 0, "kappa_exp": int(g["kappa_exp"]) % 4}


def sample(rng, n):
    return [{"schema": SCHEMA, "chi_exp": rng.integers(0, 4, 8).tolist(),
             "zeta_exp": 0, "kappa_exp": int(rng.integers(4))} for _ in range(n)]


def mutate(rng, parent):
    g = canonicalize(parent)
    # Structural block-copy mutations can discover a common character;
    # the common value is inherited, never hard-coded to the known solution.
    if rng.random() < .3:
        i = int(rng.integers(8))
        for j in rng.choice(8, size=int(rng.integers(1, 9)), replace=False):
            g["chi_exp"][int(j)] = g["chi_exp"][i]
    else:
        for j in rng.choice(9, size=int(rng.integers(1, 4)), replace=False):
            if j == 8:
                g["kappa_exp"] = int(rng.integers(4))
            else:
                g["chi_exp"][int(j)] = int(rng.integers(4))
    return g


def crossover(rng, a, b):
    return canonicalize({"chi_exp": [a["chi_exp"][j] if rng.random() < .5 else b["chi_exp"][j]
                                     for j in range(8)],
                         "kappa_exp": a["kappa_exp"] if rng.random() < .5 else b["kappa_exp"]})


def valid(g):
    return (set(g) == {"schema", "chi_exp", "zeta_exp", "kappa_exp"} and
            g["schema"] == SCHEMA and len(g["chi_exp"]) == 8 and
            all(type(x) is int and 0 <= x < 4 for x in g["chi_exp"] + [g["zeta_exp"], g["kappa_exp"]]))


def proof_candidate(g):
    return canonicalize(g)


def _normal_product(n, c, m, d, chi, zeta):
    """Exact source normal form E_n X^c E_m X^d, independent CPU oracle."""
    bilinear = lambda u, v: sum(u[i] * v[j] for i in range(8) for j in range(i))
    nm = n + m
    exponent = 2 * bilinear(n, m) + c * int(chi @ m) + c*d*zeta + 2*c*d*bilinear(nm, np.ones(8, dtype=int))
    return nm + c*d, (c+d) % 2, exponent


def _associator_coefficients():
    # Associator exponent is affine in each charge after source bilinear
    # cancellation. The fixed certificate proves that cancellation universally.
    # These 200 CUDA probes are a search filter, not the universal proof.
    charge_probes = [np.zeros(24, dtype=int)] + list(np.eye(24, dtype=int))
    rows, constants = [], []
    for bits in itertools.product(range(2), repeat=3):
        c, d, f = bits
        for probe in charge_probes:
            n, m, l = probe[:8], probe[8:16], probe[16:]
            def exponent(gene):
                chi, zeta = gene[:8], gene[8]
                nm, cd, a = _normal_product(n, c, m, d, chi, zeta)
                _, _, b = _normal_product(nm, cd, l, f, chi, zeta)
                ml, df, e = _normal_product(m, d, l, f, chi, zeta)
                _, _, h = _normal_product(n, c, ml, df, chi, zeta)
                return a + b - e - h
            zero = exponent(np.zeros(10, dtype=int))
            rows.append([exponent(e) - zero for e in np.eye(10, dtype=int)])
            constants.append(zero)
    return np.asarray(rows, dtype=np.int64), np.asarray(constants, dtype=np.int64)


def _time_fixtures():
    a = np.array([1, 1, 1, 1, 0, 0, 0, 0], dtype=np.int64)
    probes = [np.zeros(8, dtype=np.int64)] + list(np.eye(8, dtype=np.int64))
    rows = []
    for alpha_n in (1, 3):
        delta32 = np.arange(8, dtype=np.int64)*4 + alpha_n - 16
        for c in (0, 1):
            for n in probes:
                p2 = 2*n+c
                image2 = p2[PERM]+2*a
                energy64 = lambda t: 8*int(t@t)+int(delta32@t)
                ward = energy64(image2)-energy64(p2)-(16*int(p2.sum())+48+8*alpha_n)
                charge = int(image2.sum()-p2.sum())-8
                # 16 times the fine-clock difference, and 8 times its
                # quarter character after removing an integer multiple of 8.
                q = int(p2.sum())//2
                clock8 = (int(np.arange(8) @ (image2-p2))//2 + alpha_n) % 8
                target8 = (6+alpha_n+4*(q % 2)) % 8
                rows.append([ward, charge, (clock8-target8) % 8])
    return np.asarray(rows, dtype=np.int64)


def evaluate(candidates, device):
    import torch
    if device not in _CACHE:
        coeff, const = _associator_coefficients()
        _CACHE[device] = (torch.tensor(coeff, device=device, dtype=torch.int64),
                          torch.tensor(const, device=device, dtype=torch.int64),
                          torch.tensor(_time_fixtures(), device=device, dtype=torch.int64))
    coeff, const, time_fixed = _CACHE[device]
    genes = torch.tensor([g["chi_exp"]+[g["zeta_exp"],g["kappa_exp"]] for g in candidates],
                         device=device, dtype=torch.int64)
    chi = genes[:, :8]
    square = (2*chi-2).remainder(4)
    q_fixed = chi.sum(dim=1, keepdim=True).remainder(4)
    difference = torch.stack([(chi[:,i]-chi[:,j]).remainder(4) for i,j in PAIRS], dim=1)
    pairs = torch.stack([(chi[:,i]+chi[:,j]-2).remainder(4) for i,j in PAIRS], dim=1)
    p = (chi-chi[:,PERM]).remainder(4)
    kappa = (2*genes[:,9:10]).remainder(4)
    composition = (-chi[:,4:].sum(dim=1,keepdim=True)).remainder(4)
    assoc = ((genes[:,None,:]*coeff[None,:,:]).sum(dim=2)+const).remainder(4)
    # Exact roots-of-unity phases and source-scaled integer Ward probes.
    groups = [torch.cat([square,q_fixed,assoc],dim=1),torch.cat([difference,pairs,p,kappa],dim=1),composition]
    scores = torch.stack([(r != 0).sum(dim=1) for r in groups],dim=1)
    fixed_count = int((time_fixed != 0).sum().item())
    finite = (scores.sum(dim=1)==0) & (fixed_count==0)
    score_cpu, feasible_cpu = scores.cpu().tolist(), finite.cpu().tolist()
    result = []
    for g, score, feasible in zip(candidates, score_cpu, feasible_cpu):
        uniform = len(set(g["chi_exp"])) == 1
        orientation = str(g["chi_exp"][0]) if uniform else "mixed"
        result.append({"objectives": score+[fixed_count], "feasible": bool(feasible),
                       "descriptor": [orientation,g["kappa_exp"],sum(x%2 for x in g["chi_exp"])],
                       "details": {"tensor_device": str(genes.device),"dtype": str(genes.dtype),
                                   "associator_probes":int(assoc.shape[1]),"time_probes":int(time_fixed.shape[0]),
                                   "source_action":"U_s R_c = W4; U_s=i^(-zeta) X",
                                   "source_selection_proved":False,"physical_solution":False,
                                   "verdict_ceiling":"conditional_crossed_defect_algebra",
                                   "additional_premises":["given eight-channel CAR source", "adjoined half-lattice sector",
                                                          "same quadratic time on extension", "dual-Z2-invariant state extension"],
                                   "original_target":"derive source, marking, state and time jointly from P1/P2",
                                   "unresolved_origin_obligation":True}})
    return result


def exhaustive_generator_classification():
    # This algebraic finite classifier is independent of the genetic run.
    # C4 forces all entries equal; C2 makes this common entry odd; C5 makes
    # kappa even. The two rephasing-invariant orientations remain distinct.
    return [{"schema":SCHEMA,"chi_exp":[c]*8,"zeta_exp":0,"kappa_exp":k}
            for c in (1,3) for k in (0,2)]


def enumerate_candidates(shard=0, shards=1):
    """Disjoint base-four enumeration of every canonical genotype exactly once."""
    if not 0 <= shard < shards:
        raise ValueError("Require 0 <= shard < shards")
    for index in range(shard, 4**9, shards):
        value = index
        genes = []
        for _ in range(9):
            genes.append(value % 4)
            value //= 4
        yield {"schema": SCHEMA, "chi_exp": genes[:8],
               "zeta_exp": 0, "kappa_exp": genes[8]}


if __name__ == "__main__":
    print(json.dumps(evaluate(exhaustive_generator_classification(), "cpu"), indent=2))
