"""Finite causal process tomography and operational completeness."""
from itertools import combinations, product
from math import prod

import sympy as sp

from language import positive
from reconstruction import matrix, require, scalar
from routes import partial_trace


def causal_residual(operator, dims, normalization):
    # Slots are ordered [O_k,I_k,...,O_1,I_1]. Recover each shorter
    # comb by partial trace, then verify Tr_Ok C_k=I_Ik tensor C_(k-1).
    residual = []
    for n in range(len(dims),0,-2):
        axes = dims[-n:]
        left = partial_trace(operator,axes,list(range(1,n)))
        operator = partial_trace(left,axes[1:],list(range(1,n-1)))/axes[1]
        residual.extend(left-sp.kronecker_product(sp.eye(axes[1]),operator))
    residual.append(operator[0,0]-normalization)
    return residual


def causal_process(operator, dims):
    return positive(operator) and not any(causal_residual(operator,dims,1))


def hermitian_basis(d):
    for i in range(d):
        yield sp.eye(d)[:,i]*sp.eye(d)[i,:]
    for i,j in combinations(range(d),2):
        e = sp.eye(d)[:,i]*sp.eye(d)[j,:]
        yield e+e.T
        yield sp.I*(e-e.T)


def effect_basis(dims, slots):
    require(len(slots)==len(dims)//2 and all(isinstance(slot,dict) and slot for slot in slots),
            "give each time slot a nonempty dictionary of labelled interventions")
    parsed = []
    for axes,slot in zip(zip(dims[::2],dims[1::2]),slots):
        size = prod(axes)
        choices = []
        for label,data in sorted(slot.items()):
            if isinstance(data,dict):
                operation = matrix(data["effect"])
                state = matrix(data["input_state"])
            else:
                require(axes[1]==1,"interventions with an input leg need a declared input state")
                operation,state = matrix(data),sp.eye(1)
            require(operation.shape==(size,size) and positive(operation),
                    "supplied interventions must be positive effects on the declared slot")
            require(state.shape==(axes[1],axes[1]) and positive(state) and state.trace()==1,
                    "input state for the intervention must be a density operator")
            # A normalized one-slot tester sums to I_output tensor rho_input.
            # Local trace-nonincreasing channel Choi operators alone need not
            # be bounded testers when linked with a causal comb.
            require(positive(sp.kronecker_product(sp.eye(axes[0]),state)-operation),
                    "intervention effect exceeds its normalization operator")
            choices.append((label,operation))
        parsed.append(choices)
    return [(tuple(label for label,_ in word),sp.kronecker_product(*(op.T for _,op in word)))
            for word in product(*parsed)]


def solve(data):
    dims = data["dimensions"]
    require(isinstance(dims,list) and dims and len(dims)%2==0 and
            all(type(d) is int and d>=1 for d in dims) and prod(dims)<=16,
            "declare positive [output,input,...] dimensions of total size <=16")
    d = prod(dims)
    basis = list(hermitian_basis(d))
    accessible = dict(effect_basis(dims,data["interventions"]))
    # Use the supplied matrix as a causal witness if present. Alternative
    # intervention sets are compared on exactly the same affine family.
    reference = matrix(data["process"]) if "process" in data else None
    if reference is not None:
        require(reference.shape==(d,d) and causal_process(reference,dims),
                "supplied process fails positivity, normalization or ordered causality")
    holdout = {}
    if "extension_interventions" in data:
        holdout = dict(effect_basis(dims,data["extension_interventions"]))

    def row(effect):
        return [(effect*b).trace() for b in basis]

    measured = {}
    for q in data.get("observations",[]):
        word = tuple(q["word"])
        require(word in accessible and word not in measured,"observations need distinct accessible words")
        value = scalar(q["value"])
        require(value.is_Rational and 0<=value<=1,"observed probability must lie in [0,1]")
        measured[word] = value
    if reference is not None:
        replay = {word:(effect*reference).trace() for word,effect in accessible.items()}
        require(all(v.is_Rational and 0<=v<=1 for v in replay.values()),
                "supplied process and intervention convention do not give probabilities")
        require(all(replay[w]==v for w,v in measured.items()),"process does not reproduce supplied probabilities")
        if data.get("observe_all_accessible",False):
            measured = replay
    elif data.get("observe_all_accessible",False):
        raise ValueError("observe_all_accessible needs a causal reference process")

    columns = [causal_residual(b,dims,0) for b in basis]
    constraints = sp.Matrix(columns).T
    rhs = sp.zeros(constraints.rows,1)
    rhs[-1] = 1
    channel_rank = constraints.rank()
    measured_rows = sp.Matrix(len(measured),len(basis),
                              [v for word in measured for v in row(accessible[word])])
    design = constraints.col_join(measured_rows)
    target = rhs.col_join(sp.Matrix(len(measured),1,list(measured.values())))
    if design.rank()!=design.row_join(target).rank():
        return {"status":"REFUTED_SCOPED","reason":"probabilities violate the common causal linear identities"}
    hidden = design.nullspace()
    observable = {word:row(effect) for word,effect in accessible.items()}
    extension = {word:row(effect) for word,effect in holdout.items()}

    def separates(r):
        return any((sp.Matrix([r])*vec)[0]!=0 for vec in hidden)

    next_accessible = next((list(word) for word,r in observable.items() if separates(r)),None)
    next_extension = next((list(word) for word,r in extension.items() if separates(r)),None)
    reconstructed = reference
    if not hidden and reference is None:
        coefficients = design.gauss_jordan_solve(target)[0]
        reconstructed = sum((c*b for c,b in zip(coefficients,basis)),sp.zeros(d))
        if not causal_process(reconstructed,dims):
            return {"status":"REFUTED_SCOPED","reason":"unique linear candidate is not a positive causal process"}
    witness = None
    if reference is not None and reference.det()>0 and hidden:
        probe = extension.get(tuple(next_extension)) if next_extension is not None else None
        vector = next((v for v in hidden if probe is None or (sp.Matrix([probe])*v)[0]!=0))
        perturbation = sum((c*b for c,b in zip(vector,basis)),sp.zeros(d))
        step = sp.S.Half
        while not (positive(reference+step*perturbation) and positive(reference-step*perturbation)):
            step /= 2
        if probe is not None:
            witness = [str((holdout[tuple(next_extension)]*(reference+sign*step*perturbation)).trace())
                       for sign in (1,-1)]
    return {"status":"EXACT_FINITE_PROCESS_TOMOGRAPHY",
            "causal_affine_dimension":len(basis)-channel_rank,
            "linear_family_dimension":len(hidden),"full_linear_reconstruction":not hidden,
            "full_reconstruction_on_accessible_operations":next_accessible is None,
            "full_reconstruction_on_extension":next_extension is None if holdout else None,
            "causal_process_exists":True if reconstructed is not None else None,
            "rank_of_accessible_constraints":constraints.col_join(
                sp.Matrix(len(observable),len(basis),[x for r in observable.values() for x in r])).rank()-channel_rank,
            "observations":len(measured),"next_accessible_word":next_accessible,
            "next_extension_word":next_extension,"positive_pair_extension_responses":witness,
            "reconstructed_process":[list(map(str,r)) for r in reconstructed.tolist()] if not hidden else None,
            "scope":"supplied ordered finite interventions; physical admissibility and the origin of time require the source"}


def controls():
    # A two-slot causal comb with only prepare/measure interventions on
    # two emitted qubits, so all trace conditions reduce to a unit trace.
    z = {"z0":[[1,0],[0,0]],"z1":[[0,0],[0,1]]}
    tomographic = {**z,"x":[["1/2","1/2"],["1/2","1/2"]],
                   "y":[["1/2",[0,"-1/2"]],[[0,"1/2"],"1/2"]]}
    process = [["1/4" if i==j else 0 for j in range(4)] for i in range(4)]
    base = {"dimensions":[2,1,2,1],"process":process,"observe_all_accessible":True}
    restricted = solve({**base,"interventions":[z,z],"extension_interventions":[tomographic,tomographic]})
    full = solve({**base,"interventions":[tomographic,tomographic]})
    return {"status":"EXACT_TOMOGRAPHY_CONTROLS","restricted":restricted,"complete":full}
