"""Construct and distinguish finite ordered return generators."""
from itertools import product
from collections import Counter

import sympy as sp
import z3

from language import positive
from reconstruction import equivalent, require, scalar
from sewing_dynamics import HADAMARD, clifford_cells


def exact_values(result,variables):
    model = result.model()
    return [sp.Rational(model.eval(x).as_fraction()) for x in variables]


def canonical_weights(amplitudes):
    return list(HADAMARD.T*sp.Matrix(amplitudes)/4)


def opposite_signs(signs):
    return len(signs)==4 and all(type(s) is int and s in (-1,1) for s in signs) and sp.prod(signs)==-1


def squared_response(signs,weights,gap):
    return sum(s*w for s,w in zip(signs,weights))**2/gap


def sector_words():
    a,b,pa,pb = clifford_cells()
    eye = sp.eye(16)
    selector = (eye+pa)*(eye+pb)/4
    basis = sp.Matrix.hstack(*selector.columnspace())
    q_basis = sp.Matrix.hstack(*(eye-selector).columnspace())
    links = [sp.I*x*y for x,y in zip(a,b)]
    return basis,q_basis,links


def even_projectors(basis,links):
    projectors = []
    for signs in HADAMARD.tolist():
        projector = sp.eye(4)
        for i in range(1,4):
            pair = (basis.H*links[0]*links[i]*basis).applyfunc(sp.expand)
            projector = projector*(sp.eye(4)+signs[i]*pair)/2
        projectors.append(projector.applyfunc(sp.expand))
    return projectors


def verify_words(basis,links,weights,gap,observations):
    # The inverse solver only sees signed amplitudes. Independently replay
    # the CAR matrices and test each physically labelled pair insertion.
    first = sum((w*m for w,m in zip(weights,links)),sp.zeros(16))
    restricted = (basis.H*first*first*basis).applyfunc(sp.expand)/gap
    squares = HADAMARD*sp.Matrix(weights)
    for index,projector in enumerate(even_projectors(basis,links)):
        require(projector.trace()==1 and restricted*projector==observations[index]*projector and
                squares[index]**2/gap==observations[index],"inverted return failed independent CAR replay")
    return first,restricted


def find_route(basis,q_basis,links,weights,gap,observations):
    # Start with block words in V and a selected boundary. Products are
    # composed before compression. Their intermediate projection is derived
    # from the defect, not imposed at each event.
    first,expected = verify_words(basis,links,weights,gap,observations)
    direct = basis.H*first*basis
    boundary = basis*basis.H
    complement = q_basis*q_basis.H
    require(boundary+complement==sp.eye(16) and not any(boundary*complement),"incomplete block decomposition")
    counterexamples = []
    blocks = [("PVP",direct)]
    for name,inner in (("P",basis),("Q",q_basis)):
        entrance = inner.H*first*basis
        propagator = sp.eye(inner.cols)/gap
        label = "PVPVP" if name=="P" else "PVQ(QH0Q)^-1QVP"
        blocks.append((label,entrance.H*propagator*entrance))
    for name,candidate in blocks:
        difference = (candidate-expected).applyfunc(sp.expand)
        if not any(difference):
            return name,counterexamples
        row,col = next((i,j) for i,j in product(range(4),repeat=2) if difference[i,j])
        counterexamples.append({"candidate":name,"word":["V","V"],"entry":[row,col],
                                "expected":str(expected[row,col]),"obtained":str(candidate[row,col]),
                                "scope":"declared P, return scale and reached complement"})
    raise ArithmeticError("no block word reproduced the positive return")


def equivalence_probe(return_matrix,projectors):
    # Each block supplies a normalized trace for one labelled initial
    # sector. The argument is the independently replayed return matrix.
    row = sp.Matrix([[int(i==j) for i in range(4) for j in range(4)]])
    state = sum(projectors,sp.zeros(4)).reshape(16,1)/4
    alphabet = {f"sector_{i}":sp.kronecker_product(p,sp.eye(4)) for i,p in enumerate(projectors)}
    return state,row,{**alphabet,"return_generator":sp.kronecker_product(return_matrix,sp.eye(4))}


def completion_search(data):
    gap = scalar(data.get("gap",8))
    basis,q_basis,links = sector_words()
    require(gap.is_Rational and gap>0,"positive rational complement gap required")
    require(data.get("input_role")=="declared_finite_oracle","a declared oracle role is required")
    readout = data["projected_return_observations"]
    require(len(readout)==4,"one squared response per labelled even sector is required")
    observations = list(map(scalar,readout))
    require(all(y.is_Rational and y>=0 for y in observations) and any(observations),
            "four nonnegative rational responses, including a nonzero one, are required")
    variables = z3.Reals("x0 x1 x2 x3")
    # This rational subproblem has exactly two signs for each nonzero square.
    # QF_NRA certifies every assignment; a nonsquare is outside this contract.
    square_roots = [sp.sqrt(gap*y) for y in observations]
    if not all(x.is_Rational for x in square_roots):
        return {"status":"OUTSIDE_RATIONAL_WEIGHT_CLASS","physical_solution":False}
    solver = z3.SolverFor("QF_NRA")
    solver.add(*(x*x==z3.RealVal(str(gap*y)) for x,y in zip(variables,observations)))
    # Fix the simultaneous sign, an actual conjugation by one local parity.
    first_nonzero = next(i for i,y in enumerate(observations) if y)
    solver.add(variables[first_nonzero]>0)
    weights = []
    status = solver.check()
    while status==z3.sat:
        amps = exact_values(solver,variables)
        weights.append(canonical_weights(amps))
        solver.add(z3.Or(*(x!=z3.RealVal(str(a)) for x,a in zip(variables,amps))))
        status = solver.check()
    require(status==z3.unsat and len(weights)==2**(sum(bool(y) for y in observations)-1),
            "sector enumeration failed its analytic census")
    weights.sort(key=lambda w:tuple(w))
    if data.get("require_positive_weights"):
        weights = [w for w in weights if all(x>0 for x in w)]
    used = set()
    for sample in data.get("opposite_parity_observations",[]):
        signs,value = sample["signs"],scalar(sample["value"])
        require(opposite_signs(signs) and value.is_Rational and value>=0,
                "a measured opposite parity response must be a nonnegative rational")
        used.add(tuple(s*signs[0] for s in signs))
        weights = [w for w in weights if squared_response(signs,w,gap)==value]
    if not weights:
        return {"status":"NO_CANDIDATE_IN_DECLARED_CLASS","physical_solution":False}
    symbol,counterexamples = find_route(basis,q_basis,links,weights[0],gap,observations)

    # Keep the clock/reference mode fixed. Compare the unordered absolute
    # weights too, so even an arbitrary permutation/sign change cannot
    # explain a reported pair with different full spectra.
    signatures = {tuple(sorted(abs(x) for x in w)):w for w in weights}
    split = list(signatures.values())
    witnesses = [split[0]] if len(split)==1 else [split[0],split[-1]]
    verified = []
    for witness in witnesses:
        _,actual = verify_words(basis,links,witness,gap,observations)
        verified.append(actual)
    match = len(witnesses)==2 and verified[0]==verified[1]
    projectors = even_projectors(basis,links)
    oracle_equivalence = equivalent(*(equivalence_probe(r,projectors) for r in verified)) if match else None
    if oracle_equivalence is not None:
        oracle_equivalence["operator_identity_in_independent_replay"] = True
    held_out = []
    for s in product((1,-1),repeat=3):
        signs = (1,*s)
        if not opposite_signs(signs) or signs in used:
            continue
        counts = Counter(squared_response(signs,w,gap) for w in weights)
        if len(counts)>1:
            vals = [squared_response(signs,w,gap) for w in witnesses]
            separated = (len(weights)**2-sum(c**2 for c in counts.values()))//2
            held_out.append({"kind":"opposite_local_parity","signs":list(signs),"answers":list(map(str,vals)),
                             "separation":str(max(vals)-min(vals)), "pairs_separated":separated,
                             "answer_sizes":{str(k):v for k,v in sorted(counts.items())}})
    held_out.sort(key=lambda x:(-x["pairs_separated"],-sp.Rational(x["separation"])))
    remaining = [] if len(witnesses)<2 else [
        {"kind":"full_generator_power","power":2*k,
         "answers":[str(sum(sum(a*b for a,b in zip(s,w))**(2*k)
                             for s in product((1,-1),repeat=4))/16) for w in witnesses]}
        for k in range(1,5)
    ]
    next_time = next((x for x in remaining if len(set(x["answers"]))>1),None)
    return {"status":"UNDERDETERMINED_ON_DECLARED_INTERFACE" if match else "FINITE_CANDIDATES_ENUMERATED",
            "input_role":"declared_finite_oracle","gap":str(gap),
            "contract":{"kinematics":"crossed CAR with fixed labels and local parity",
                        "state":"normalized readout on each labelled even joint sector",
                        "dynamics":"diagonal Majorana coupling and fixed gap on the reached complement",
                        "observations":"four squared amplitudes, absolute energy reference fixed"},
            "learned_operation":symbol,"counterexamples":counterexamples,
            "number_of_candidates_up_to_overall_sign":len(weights),
            "candidates":[list(map(str,w)) for w in weights],
            "witnesses":[list(map(str,w)) for w in witnesses],"same_interface":oracle_equivalence,
            "next_discriminating_observation":held_out[0] if held_out else None,
            "first_full_trace_power_separator":next_time,
            "all_discriminating_observations":held_out,"physical_solution":False}


def held_out_search(data):
    """Replay both answers of a computed experiment before a real answer exists."""
    gap = scalar(data["gap"])
    signs = data["signs"]
    require(opposite_signs(signs),
            "an opposite local parity observation is required")
    weights = [list(map(scalar,w)) for w in data["witnesses"]]
    require(len(weights)>=2 and all(len(w)==4 for w in weights),"two or more witness generators required")
    require(gap>0,"positive gap required")
    # Independent block for the opposite local parities.
    a,b,pa,pb = clifford_cells()
    eye = sp.eye(16)
    p = (eye+pa)*(eye-pb)/4
    modes = [sp.I*x*y for x,y in zip(a,b)]
    from sewing_dynamics import joint_projector
    pair = joint_projector(modes,signs)+joint_projector(modes,[-s for s in signs])
    sector = (pair*p).applyfunc(sp.expand)
    require(sector.trace()==1,"held-out probe selected the wrong sector")
    branches = {}
    for witness in weights:
        generator = sum((w*m for w,m in zip(witness,modes)),sp.zeros(16))
        ret = (sector*generator*generator*sector).applyfunc(sp.expand)/gap
        answer = sp.expand(ret.trace())
        require(positive(ret),"held-out return is not positive")
        branches.setdefault(str(answer),[]).append(list(map(str,witness)))
    return {"status":"HELD_OUT_PROBE_EXECUTED","input_role":"declared_finite_oracle",
            "hypothetical_branches":branches,"separates_witnesses":len(branches)>1,
            "remaining_candidate_counts":{v:len(w) for v,w in branches.items()},
            "original_observation":None}
