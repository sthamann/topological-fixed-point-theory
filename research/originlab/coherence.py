"""Sparse CAR and exact finite memory tests for coherent boundary data."""
from functools import lru_cache
from itertools import combinations, product
from math import comb, lcm
from pathlib import Path

import sympy as sp

from reconstruction import PAULI, equivalent, matrix, reachable, require, scalar
from verify import digest

R104_SOURCES = (
    Path("experiments/theory-contracts/original-local-source-variable-background-20261102/"
         "DET16-WHOLE-SOURCE-REDUCTION.txt"),
    Path("verification/v1002_spin10_mirror_projector.py"),
    Path("experiments/theory-contracts/original-paper-source-reconciliation-20261102/"
         "sewing_20261105C/sewing_route1/r104-whole-number-source-review.txt"),
)


def car_letter(state, mode, creation):
    occupied = bool(state & (1 << mode))
    if occupied == creation:
        return None, 0
    sign = -1 if (state & ((1 << mode)-1)).bit_count() % 2 else 1
    return state ^ (1 << mode), sign


def car_word(vector, word):
    for mode, creation in reversed(word):
        mapped = {}
        for state, value in vector.items():
            target, sign = car_letter(state, mode, creation)
            if sign:
                mapped[target] = sign*value
        vector = mapped
    return vector


def inner(left, right):
    return sp.simplify(sum(sp.conjugate(value)*right.get(key, 0)
                           for key, value in left.items()))


def thermal(vector, ground, a):
    amplitude = inner(ground, vector)
    return {key:sp.simplify(a*vector.get(key,0)+(1-a)*amplitude*ground.get(key,0))
            for key in set(vector)|set(ground)}


def exact_small_ranks(n):
    # Repeated actual CAR application on sparse vectors. Common normalization
    # can be omitted. The phase i does not affect the support proof.
    frontier = [{0:1, (1<<n)-1:sp.I}]
    for depth in range(max(2,(n+1)//2)+1):
        signatures = {tuple(sorted(v.items())):v for v in frontier if v}
        singletons = {next(iter(v)) for v in signatures.values() if len(v)==1}
        # At depth zero the cap vector is one direction. From depth one
        # the only nonsingleton is that original cap vector.
        dim = len(singletons)+int(not {0,(1<<n)-1}<=singletons)
        yield dim
        frontier = [*signatures.values(),*(car_word(v,[(mode,creation)])
                    for v in signatures.values() for mode,creation in product(range(n),(False,True)))]


def word_ranks(n):
    depths = max(2,(n+1)//2)
    # At k>=2 the endpoint components are separated by n_j and 1-n_j.
    # A basis vector is then reachable iff its Hamming distance to one
    # endpoint is <= k. At k=1 the cap is still one endpoint direction.
    return [1, 1+min(2*n,(1<<n)-2), *[
        sum(comb(n,j) for j in range(n+1) if min(j,n-j)<=k)
        for k in range(2,depths+1)
    ]]


@lru_cache(maxsize=4)
def operator_cut_rank(left, right):
    # I, |0><0|, |F><F|, |0><F|, |F><0| on each side. Use their
    # Hilbert-Schmidt Gram rather than allocating a Fock-space matrix.
    def gram(n):
        return sp.Matrix([[2**n,1,1,0,0],[1,1,0,0,0],[1,0,1,0,0],
                          [0,0,0,1,0],[0,0,0,0,1]])
    # All five coefficients are strictly positive for 0 < beta < infinity.
    # Both Gram matrices are invertible when both sides have >=2 modes.
    # If exactly one side is 1-mode the other still gives rank 4. If both
    # sides are 1-mode the Gram matrices coincide and a positive diagonal
    # coefficient preserves their rank.
    coefficient = sp.eye(5)
    return (gram(left)*coefficient*gram(right)).rank()


@lru_cache(maxsize=1)
def proper_support_control():
    cap = sp.Matrix([1,*([0]*14),sp.I])/sp.sqrt(2)
    orthogonal = sp.conjugate(cap)
    checked = 0
    for support in range(4):
        for sites in combinations(range(4),support):
            for axes in product((1,2,3),repeat=support):
                word = dict(zip(sites,axes))
                op = sp.kronecker_product(*(PAULI[word.get(i,0)] for i in range(4)))
                require((cap.H*op*cap)[0]==(orthogonal.H*op*orthogonal)[0],
                        "a proper-support operator distinguished the two phase caps")
                checked += 1
    return checked


def r104(root, modes=16, phase_pi="1/7", modulus=4):
    require(type(modes) is int and 2<=modes<=64,"coherent occupation family requires 2 to 64 modes")
    require(type(modulus) is int and modulus>0,"positive number-sector modulus required")
    phase = scalar(phase_pi)
    require(phase.is_Rational,"phase must be a real rational multiple of pi")
    phase %= 2
    phase_text = f"exp(i*pi*({phase}))"
    z = sp.exp(sp.I*sp.pi*phase)
    full = (1 << modes)-1
    ground = {0:1/sp.sqrt(2),full:z/sp.sqrt(2)}
    p0 = sp.diag(1,0)
    p1 = sp.diag(0,1)
    cap = sp.diag(1,z)
    # Local products kill all mixed occupation strings; one nonempty cut
    # already has a Schmidt matrix with its two nonzero diagonal entries.
    require(p0*p1==sp.zeros(2) and (cap*p0**modes).trace()==1 and
            (cap*p1**modes).trace()==z,"coherent two-channel cap did not replay")
    spatial = sp.kronecker_product(p0,p0)+sp.kronecker_product(p1,p1)
    k = (modes+1)//2
    create = [(i,True) for i in range(k)]
    destroy = [(i,False) for i in reversed(range(k,modes))]
    overlap = inner(car_word(ground,create),car_word(ground,destroy))
    determinant = [(i,False) for i in reversed(range(modes))]
    split_sign = (-1)**(k*(modes-k))
    require(sp.simplify(sp.expand_complex(overlap-split_sign*z/2))==0 and
            sp.simplify(sp.expand_complex(inner(ground,car_word(ground,determinant))-z/2))==0,
            "CAR orientation or fixed-reference response failed")
    # A fixed determinant insertion sees the phase, even with time steps.
    time_word = inner(ground, thermal(car_word(thermal(ground,ground,sp.Rational(2,3)),determinant),
                                      ground,sp.Rational(1,2)))
    require(sp.simplify(sp.expand_complex(time_word-z/2))==0,"ordered fixed-reference witness failed")
    # On each nontrivial proper occupation subset the complement is
    # orthogonal between the two cap branches. Petz reduces to the two
    # matching classical labels on the B support.
    ghz2 = sp.Matrix([[sp.Rational(1,2),sp.conjugate(z)/2],[z/2,sp.Rational(1,2)]])
    mixture = sp.eye(2)/2
    diff = ghz2-mixture
    require((diff*diff).applyfunc(sp.simplify)==sp.eye(2)/4 and diff.trace()==0,
            "wrong trace distance for the reconstructed mixture")
    overlap_probability = (ghz2*mixture).trace()
    psd = 1<<modes
    selected = sum(comb(modes,j) for j in range(0,modes+1,modulus))
    cap_is_selected = modes % modulus == 0
    beta = sp.Symbol("beta",nonnegative=True)
    partition = selected*sp.exp(-2*beta)+(1 if cap_is_selected else sp.Rational(1,2))*(1-sp.exp(-2*beta))
    ranks = word_ranks(modes)
    original = (modes,phase,modulus)==(16,sp.Rational(1,7),4)
    refs = {str(path):digest(root/path) for path in R104_SOURCES}
    require(not original or list(refs.values())==[
        "fc43fcfef453bb278b8f5c28123078331a36d1e4350672b45c47a5d99c992897",
        "84f7cd7987887ee071d7e05872285ab274ae329f6bcfa693be69bb52add19a4c",
        "32e9d5d630ebb061cd4fa41aec9be6e9d094364df27e20b59264751a84df4c22",
    ],"original R104 sources changed; review the literal DET factor again")
    result = {"status":"EXACT_R104_FACTOR" if original
              else "EXACT_CONDITIONAL_VARIANT",
              "scope":"separated DET factor on occupation modes with the stated Hamiltonian",
              "original_sources":refs,
              "modes":modes,"phase":phase_text,"fock_dimension":psd,
              "ground_preserves_fermion_parity":modes % 2 == 0,
              "markov":{"all_nonempty_mode_tripartitions":"REFUTED_SCOPED" if modes>=3 else "NO_TRIPARTITION",
                        "conditional_information_bits":1 if modes>=3 else None,
                        "petz_overlap_probability":str(overlap_probability) if modes>=3 else None,
                        "petz_trace_distance":"1/2" if modes>=3 else None,
                        "proof":"on any proper subset the complement kills the off-diagonal cap term"},
              "cap":{"state_bond_dimension":2,
                     "operator_schmidt_ranks_positive_finite_time":
                         [operator_cut_rank(min(c,2),min(modes-c,2)) for c in range(1,modes)],
                     "spatial_transfer_spectrum":sorted(map(str,sp.flatten(
                         [[e]*m for e,m in spatial.eigenvals().items()]))),
                     "temporal_gap":"2","temporal_ground_rank":1,
                     "scope":"mode bipartitions, not a spatial locality proof"},
              "coherent_gram":{"ranks":ranks,"first_full_depth":ranks.index(psd),
                  "first_phase_sensitive_gram_depth":k,
                  "minimum_equal_time_phase_degree":modes,
                  "signed_split_overlap":str(overlap),
                  "determinant_word":[[i+1,"annihilate"] for i,_ in determinant],
                  "fixed_reference_timed_value":str(time_word)},
              "number_sector":{"modulus":modulus,"dimension":selected,
                               "commutes_with_hamiltonian":cap_is_selected,
                               "absolute_filtered_trace":str(sp.collect(sp.expand(partition),sp.exp(-2*beta))),
                               "absolute_unfiltered_trace":str(1+(psd-1)*sp.exp(-2*beta))},
              "frame_equivalence":{"all_number_preserving_ordered_words":"SAME_AFTER_FRAME_TRANSPORT",
                "proof":"U=exp(i theta N/n) conjugates state and heat, and fixes each number-preserving insertion",
                "fixed_determinant_reference":"EQUAL_TO_PHASE_ZERO" if phase==0 else "DISTINGUISHES_FROM_PHASE_ZERO",
                "original_p1_reference_bound":False},
              "locality":{"vacuum_to_full_hamiltonian":str(-z),
                          "minimum_CAR_degree_for_that_entry":modes},
              "physical_original_values":dict.fromkeys(("Z0","ZU","ZY","ZUY"))}
    for small in (2,3,4,6):
        require(list(exact_small_ranks(small))==word_ranks(small),"independent sparse CAR rank replay failed")
    # Every proper-support operator has a zero cross term because the
    # missing complement is orthogonal. Equal energies at the spectral
    # minimum put both orthogonal cap states in the ground eigenspace.
    result["locality"]["indistinguishable_orthogonal_phase_caps_for_every_proper_support"] = True
    result["locality"]["unique_ground_impossible_for_sum_of_proper_support_terms"] = True
    result["locality"]["four_mode_proper_support_pauli_checks"] = proper_support_control()
    return result


def det_word(modes, phase_pi, word, intervals, number_modulus=None):
    """Absolute thermal CAR words, with an optional single outer number selector."""
    require(type(modes) is int and 2<=modes<=64,"2 to 64 occupation modes required")
    phase = scalar(phase_pi)
    require(phase.is_Rational,"phase must be a real rational multiple of pi")
    require(isinstance(word,list) and all(isinstance(item,list) and len(item)==2 and
            type(item[0]) is int and 1<=item[0]<=modes and item[1] in ("create","annihilate")
            for item in word),"invalid CAR word; use 1-based mode and create/annihilate")
    require(isinstance(intervals,list) and len(intervals)==len(word)+1,
            "one Euclidean interval before, between and after all CAR insertions required")
    times = [scalar(v) for v in intervals]
    require(all(t.is_Rational and t>=0 for t in times),"nonnegative rational intervals required")
    require(number_modulus is None or type(number_modulus) is int and number_modulus>0,
            "invalid optional number modulus")
    letters = [(i-1,op=="create") for i,op in word]
    if number_modulus is not None:
        letters.insert(0,None)
        times.insert(0,sp.S.Zero)
    denominator = lcm(*(int(t.q) for t in times))
    q = sp.Symbol("q")
    z = sp.Symbol("z",nonzero=True)
    full = (1 << modes)-1

    def apply(initial, block):
        state, sign = initial, 1
        for letter in reversed(block):
            if letter is None:
                if state.bit_count() % number_modulus:
                    return None,0
            else:
                state, factor = car_letter(state,*letter)
                if not factor:
                    return None,0
                sign *= factor
        return state,sign

    def cap_value(block):
        value = sp.S.Zero
        for initial,amplitude in ((0,sp.S.One),(full,z)):
            state,sign = apply(initial,block)
            if state in (0,full):
                value += sign*amplitude*(1 if state==0 else 1/z)/2
        return value

    def fock_trace():
        forced = {}
        for letter in reversed(letters):
            if letter is not None:
                mode,create = letter
                forced.setdefault(mode,not create)
        initial = sum(1<<mode for mode,occupied in forced.items() if occupied)
        state, sign, selected_occupancy = initial, 1, None
        for letter in reversed(letters):
            if letter is None:
                selected_occupancy = state.bit_count()
            else:
                state,factor = car_letter(state,*letter)
                if not factor:
                    return 0
                sign *= factor
        if state!=initial:
            return 0
        free = modes-len(forced)
        if selected_occupancy is None:
            return sign*(1 << free)
        return sign*sum(comb(free,n) for n in range((-selected_occupancy) % number_modulus,
                                                    free+1,number_modulus))

    count = len(letters)
    beta = sum(times)
    if not count:
        value = 1+full*q**int(beta*denominator)
    else:
        delta = [times[0]+times[-1],*times[1:-1]]
        exponents = [int(t*denominator) for t in delta]
        cyclic = letters+letters
        value = q**sum(exponents)*fock_trace()

        @lru_cache(maxsize=None)
        def edge(i,j):
            intervening = sum(exponents[k % count] for k in range(i+1,j))
            return (1-q**exponents[i % count])*q**intervening*cap_value(cyclic[i:j])

        # For each lowest projector cut, sum over the increasing paths of
        # later cuts and close the last path around the trace.
        for i in range(count):
            paths = {i:sp.S.One}
            for j in range(i+1,count):
                paths[j] = sp.expand(sum(paths[h]*edge(h,j) for h in range(i,j)))
            value += sum(paths[h]*edge(h,i+count) for h in range(i,count))
    value = sp.expand(value)
    degree = sum(1 if op=="create" else -1 for _,op in word)
    actual = value.subs({q:sp.exp(-sp.Rational(2,denominator)),z:sp.exp(sp.I*sp.pi*phase)})
    return {"status":"EXACT_ABSOLUTE_DET_WORD","modes":modes,"phase_pi":str(phase),
            "time_intervals":intervals,"word":word,"number_modulus":number_modulus,
            "number_charge":degree,"phase_dependence":z in value.free_symbols,
            "phase_exponents":sorted({int(term.as_powers_dict().get(z,0)) for term in sp.Add.make_args(value)}) if value else [],
            "polynomial_in_q_z":str(value),"q_definition":f"exp(-2/{denominator})",
            "value":str(actual),"scope":"separated DET factor with supplied Euclidean intervals"}


def trace_model(data):
    cap = matrix(data["cap"])
    ops = {k:matrix(v) for k,v in data["operations"].items()}
    require(cap.rows==cap.cols and ops and all(isinstance(k,str) and k for k in ops) and
            all(m.shape==cap.shape for m in ops.values()),
            "cap and marked operations must share one square memory")
    dim = cap.rows
    # vec_row(A X)=(A tensor I)vec_row(X); trace(BX)=vec_row(B^T)^T vec_row(X).
    return (sp.eye(dim).reshape(dim**2,1),cap.T.reshape(1,dim**2),
            {k:sp.kronecker_product(a,sp.eye(dim)) for k,a in ops.items()})


def word_value(model, word):
    initial, row, operations = model
    require(isinstance(word,list) and all(isinstance(k,str) and k in operations for k in word),
            "word has an undeclared operation")
    for label in reversed(word):
        initial = operations[label]*initial
    return (row*initial)[0]


def compare_memory(left, right):
    require(isinstance(left.get("interface"),dict) and bool(left["interface"]) and
            left["interface"]==right.get("interface"),"same labelled state/time/reference interface required")
    a,b = trace_model(left),trace_model(right)
    result = equivalent(a,b)
    result["scope"] = "all absolute words, including the empty word, on the supplied finite memory alphabet"
    result["original_source_identity_proved"] = False
    if "word" in result:
        result["shortest_word_length"] = len(result["word"])
    return result


def cap_constraints(data):
    """Solve for all caps obeying the supplied exact linear word responses."""
    ops = {k:matrix(v) for k,v in data["operations"].items()}
    require(ops,"at least one labelled operation required")
    d = next(iter(ops.values())).rows
    require(all(x.shape==(d,d) for x in ops.values()),"cap algebra must be square")
    initial,_,letters = trace_model({"cap":[[int(i==j) for j in range(d)] for i in range(d)],
                                      "operations":data["operations"]})

    def row(word):
        vec = initial
        for label in reversed(word):
            require(label in letters,"word has undeclared label")
            vec = letters[label]*vec
        return vec.T

    observations = data["observations"]
    design = sp.Matrix.vstack(sp.zeros(0,d*d),*(row(q["word"]) for q in observations))
    values = sp.Matrix(len(observations),1,[scalar(q["value"]) for q in observations])
    if design.row_join(values).rank()!=design.rank():
        return {"status":"REFUTED_SCOPED","reason":"inconsistent coherent word data"}
    particular = design.gauss_jordan_solve(values)[0] if observations else sp.zeros(d*d,1)
    particular = particular.subs(dict.fromkeys(particular.free_symbols,0))
    kernel = design.nullspace()
    basis = reachable(initial,letters)
    invisible = all((vector.T*direction)[0]==0 for _,vector in basis for direction in kernel)
    next_word = next((list(word) for word,vector in basis
                      if any((vector.T*direction)[0]!=0 for direction in kernel)),None)
    return {"status":"ALL_CAPS_OPERATIONALLY_EQUIVALENT" if invisible else "UNDERDETERMINED",
            "ambient_linear_dimension":d*d,"affine_cap_dimension":len(kernel),
            "all_word_rank":len(basis),"next_observable":next_word,
            "one_affine_solution":[[str(v) for v in line] for line in particular.reshape(d,d).T.tolist()],
            "scope":"complex-linear cap constraints on supplied operations; positivity, physical admissibility and source identity are separate"}


def controls(root, modes=16, phase_pi="1/7", modulus=4):
    result = r104(root,modes,phase_pi,modulus)
    word = result["coherent_gram"]["determinant_word"]
    result["absolute_word_controls"] = {
        "empty":det_word(modes,phase_pi,[],[1],modulus),
        "determinant":det_word(modes,phase_pi,word,[1]+[0]*modes,modulus),
        "neutral_with_intermediate_time":det_word(modes,phase_pi,
            [[1,"create"],[1,"annihilate"]],["1/5","2/5","3/5"],modulus),
    }
    return result
