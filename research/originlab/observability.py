"""Exact generated word spaces on declared observed and extended alphabets."""
from itertools import combinations
from pathlib import Path

import sympy as sp

from reconstruction import matrix, require
from verify import digest

Q4_SOURCE = Path("experiments/theory-contracts/compiler-attachment-audit-20261112/PROOF.txt")
Q4_SHA = "0faee2e5bfec552b4815dfc8e1e526355fcdaaddc551e9fafa26c72712501681"


def word_space(operations):
    dim = next(iter(operations.values())).rows
    seen,words,ranks = {},[],[]
    layer = [((),sp.eye(dim))]
    while layer:
        following = []
        for word,op in layer:
            vector = list(op)
            for index,basis in sorted(seen.items()):
                coefficient = vector[index]
                if coefficient:
                    vector = [a-coefficient*b for a,b in zip(vector,basis)]
            pivot = next((i for i,x in enumerate(vector) if x!=0),None)
            if pivot is None:
                continue
            seen[pivot] = [x/vector[pivot] for x in vector]
            words.append((word,op))
            following.extend(((name,*word),letter*op) for name,letter in sorted(operations.items()))
        if not ranks or len(words)!=ranks[-1]:
            ranks.append(len(words))
        layer = following
    return ranks,words


def compare_alphabets(data):
    raw = data["operations"]
    require(isinstance(raw,dict) and raw and all(isinstance(k,str) and k for k in raw),
            "named word operations are required")
    ops = {k:matrix(v) for k,v in raw.items()}
    dim = next(iter(ops.values())).rows
    require(all(m.shape==(dim,dim) for m in ops.values()),"all marks must use one square carrier")
    base = data["accessible"]
    require(isinstance(base,list) and len(set(base))==len(base) and set(base)<=set(ops),
            "accessible operation labels must be distinct and declared")
    ranks,limited = word_space({k:ops[k] for k in base}) if base else ([1],[((),sp.eye(dim))])
    total,full = word_space(ops)
    hidden = None
    if total[-1]!=ranks[-1]:
        null = sp.Matrix([list(m) for _,m in limited]).nullspace()
        for word,op in full:
            if any((sp.Matrix([list(op)])*v)[0]!=0 for v in null):
                hidden = list(word)
                break
        require(hidden is not None,"missing exact extension separator")
    return {"status":"EXACT_ORDERED_ALPHABET_RANK","carrier_dimension":dim,
            "accessible_ranks_by_depth":ranks,"extended_ranks_by_depth":total,
            "accessible_stable_rank":ranks[-1],"extended_stable_rank":total[-1],
            "new_algebraic_directions":total[-1]-ranks[-1],"first_extension_word":hidden,
            "scope":"linear operator span in the declared carrier; state support and physical interventions are separate"}


def compiler_control(root):
    path = root/Q4_SOURCE
    require(digest(path)==Q4_SHA,"the original Q4 comparison frame changed")
    neutral = {}
    for i,j in combinations(range(5),2):
        sign = sp.Matrix([1 if k in (i,j,5) else -1 for k in range(6)])
        neutral[f"P{i}{j}"] = sign*sign.T/6
    reference = sp.ones(6,1)/6
    direction = sp.eye(6)[:,0]-sp.eye(6)[:,1]
    bridge = reference*direction.T
    ops = {**neutral,"X":bridge,"X_adjoint":bridge.H}
    result = compare_alphabets({"operations":{k:[[str(x) for x in r] for r in m.tolist()] for k,m in ops.items()},
                                "accessible":list(neutral)})
    require(result["accessible_ranks_by_depth"]==[1,11,26] and result["extended_ranks_by_depth"]==[1,13,34,36],
            "the declared Bell reference comparison did not replay")
    return {**result,"source_sha256":Q4_SHA,
            "reference_role":"explicit algebraic reference bridge, physical admissibility open"}
