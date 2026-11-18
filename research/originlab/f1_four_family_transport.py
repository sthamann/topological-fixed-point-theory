"""Exact ordered return through the fourth family of the D8 half-spinor."""

from itertools import combinations, product

import sympy as sp

from reconstruction import require


def letter(state, mode, create):
    occupied = (state >> mode) & 1
    if occupied == create:
        return None
    phase = (-1)**((state & ((1 << mode)-1)).bit_count())
    return state ^ (1 << mode), phase


def action_of(states, i, j, a, b):
    action = {}
    for s in states:
        first = letter(s,j,b)
        if first is None:
            continue
        last = letter(first[0],i,a)
        if last is not None:
            action[s] = (last[0],first[1]*last[1])
    return action


def matrix_of(coordinates, action):
    size = len(coordinates)
    return sp.MutableSparseMatrix(size,size,{(coordinates[t],coordinates[s]):sign
                                            for s,(t,sign) in action.items()})


def gauge_neutral_control(states, coordinates, kept, outside, complement):
    # The 3+2 number-preserving gauge algebra retains both singlet
    # occupation levels. The family root acts on the other three modes.
    seed = next(s for s in outside if s & 31 == 0)
    other = seed ^ 24
    active = {seed,complement[seed]}
    inactive = {other,complement[other]}
    root = action_of(states,5,6,1-((seed>>5)&1),1-((seed>>6)&1))
    E = matrix_of(coordinates,{s:v for s,v in root.items() if s & 31 in (0,31)})
    A = E+E.H
    require(A**3==A and sum(A.row(coordinates[s]).norm() for s in outside)==2,
            "gauge-neutral family return must use two omitted states")
    gauge = [matrix_of(coordinates,action_of(states,i,j,1,0))
             for block in (range(3),range(3,5)) for i,j in product(block,repeat=2) if i!=j]
    identity = sp.eye(128)
    zero = sp.zeros(128)
    require(all(A*g==g*A for g in gauge),"the family insertion must be a gauge singlet")
    P = sp.diag(*(int(s in kept) for s in states))
    t1,t2 = sp.Rational(1,2),sp.Rational(1,4)
    numerators, spectra, refs = [], [], []
    for high in (active,inactive):
        diagonal = [sp.S.One if s in kept else t1 if s in high else t2 for s in states]
        T = sp.diag(*diagonal)
        require(all(T*g-g*T==zero for g in gauge)
                and all(diagonal[coordinates[s]]==diagonal[coordinates[complement[s]]] for s in states),
                "the comparison transfer must preserve SM blocks and the reflection pair")
        trace = lambda op:sp.trace(P*op*P)
        numerators.append([trace(T**3),trace(T**2*(identity+A)*T),
                           trace(T*A*T**2),trace(T*A*T*(identity+A)*T)])
        spectra.append(sorted(diagonal))
        refs.append(sp.trace(T**3))
    require(numerators==[[96,96,0,1],[96,96,0,sp.Rational(1,2)]]
            and spectra[0]==spectra[1] and refs[0]==refs[1]==sp.Rational(3095,32),
            "the gauge-neutral full-word alternative must have a different connected value")
    return {"connected_UY":["1/96","1/192"],"common_reflected_trace_reference":"3095/32",
            "symmetry":"number-preserving SU(3) x SU(2) x U(1), occupation-complement paired transfer",
            "same_gauge_resolved_spectrum_claimed":False}


def odd_root_intermediate():
    # Root coordinates are doubled. D8 roots and the even-parity
    # half-spinor are then the two eigenspaces of the Z2 grading.
    integer = {tuple(a if k==i else b if k==j else 0 for k in range(8))
               for i,j in combinations(range(8),2) for a,b in product((-2,2),repeat=2)}
    spinor = {s for s in product((-1,1),repeat=8) if s.count(-1)%2==0}
    counts = {sum(tuple(x+y for x,y in zip(a,b)) in integer for b in spinor) for a in spinor}
    require(counts=={28},"spinorial brackets must have a uniform even-root census")
    outgoing = counts.pop()
    return {"odd_root_operators":len(spinor),"even_root_targets_per_odd_operator":outgoing,
            "additional_cartan_image_per_odd_operator":1,"conditional_odd_projection_is_zero":True,
            "scope":"root-support audit of adjoint E8; includes the Cartan image of the opposite root"}


def single_window_cartan_square(states, actions, keep):
    def second_moment(section):
        weights = [[(1 if (s>>i)&1 else -1) for i in range(5)] for s in section]
        return sp.Matrix(weights).T*sp.Matrix(weights)

    subset = set(keep)
    require(second_moment(sorted(subset)) == len(subset)*sp.eye(5),
            "the retained charge covariance must be isotropic")
    for projection in (False,True):
        for action in actions.values():
            pairs = [(s,t) for s,(t,_) in action.items()
                     if s in subset and (not projection or t in subset)]
            for section in zip(*pairs):
                require(second_moment(section) == len(section)*sp.eye(5),
                        "positive one-window Cartan square must have zero connected response")
    t3_counts = {value:sum((((s>>3)&1)-((s>>4)&1))**2==value for s in subset)
                 for value in (0,1)}
    return {"mixed_roots_checked":len(actions),"positive_same_window_square_connected":"0",
            "checked_at_both_insertions_and_with_intermediate_family_projection":True,
            "T3_squared_multiplicities":{"0":t3_counts[0],"1/4":t3_counts[1]},
            "scope":"all real D5 Cartan linear combinations, unit or outer-only family restriction, uniform trace"}


def controls():
    states = tuple(s for s in range(256) if s.bit_count()%2 == 0)
    coordinates = {s:i for i,s in enumerate(states)}
    family = lambda s: s >> 5
    weight = lambda s: min(family(s),7-family(s))
    complement = {s:s ^ 255 for s in states}
    roots = tuple(product(range(5),range(5,8),(0,1),(0,1)))
    actions = {}
    for i,j,a,b in roots:
        action = action_of(states,i,j,a,b)
        require(len(action)==32 and not set(action).intersection(t for t,_ in action.values()),
                "each mixed root must act as a nilpotent signed partial isometry")
        actions[i,j,a,b] = action

    family_reports = []
    for omitted in range(4):
        outside = [s for s in states if weight(s)==omitted]
        kept = set(states)-set(outside)
        require(len(kept)==96 and len(outside)==32,"wrong family restriction")
        inner = sum(s in kept and t in kept for action in actions.values() for s,(t,_) in action.items())
        leaks = sum(s in kept and t not in kept for action in actions.values() for s,(t,_) in action.items())
        require((inner,leaks)==(960,480),"mixed root census changed")
        arrival = [{t for s,(t,_) in action.items() if s in kept and t not in kept}
                   for action in actions.values()]
        active = [a for a in arrival if a]
        require(len(active)==60 and all(len(a)==8 for a in active),"every mixed root must leak")
        one = sp.Matrix([[int(s in a) for s in outside] for a in active])
        supports = sorted(set(tuple(sorted(a)) for a in active))
        two = sp.Matrix([[int(s in a and s in b) for s in outside]
                         for a,b in combinations(map(set,supports),2)])
        ranks = one.rank(),one.col_join(two).rank()
        require(ranks==(12,32),"diagonal fourth-family probe ranks changed")
        reflected_axes = [s for s in outside if s<complement[s]]
        folding = sp.Matrix([[int(s in (r,complement[r])) for r in reflected_axes] for s in outside])
        reflected_ranks = (one*folding).rank(),(one.col_join(two)*folding).rank()
        require(reflected_ranks==(6,16),"reflected diagonal probe ranks changed")
        parity_profile = sp.Matrix([(-1)**((s&3).bit_count()) for s in outside])
        paired_responses = two*parity_profile
        require(sum(parity_profile)==0 and one*parity_profile==sp.zeros(one.rows,1)
                and max(abs(v) for v in paired_responses)==4,
                "an isospectral reflection-paired perturbation must survive every single return")

        root = roots[0]
        data = actions[root]
        E = matrix_of(coordinates,data)
        A = E+E.H
        U,Y = sp.eye(128)+A,A
        require(A**3==A and U*U==2*U-(sp.eye(128)-A**2),
                "the early comparison insertion is positive on each root block")
        P = sp.diag(*(int(s in kept) for s in states))
        Q = sp.eye(128)-P
        used = {t for s,(t,_) in data.items() if s in kept and t not in kept}
        used |= {s for s,(t,_) in data.items() if s not in kept and t in kept}
        require(len(used)==16 and {complement[s] for s in used}==used,
                "the excursion support must be reflection paired")
        numerators, references, reduced, spectra = [], [], [], []
        for high,low in ((sp.Rational(1,2),sp.Rational(1,4)),(sp.Rational(1,4),sp.Rational(1,2))):
            diagonal = [sp.S.One if s in kept else high if s in used else low for s in states]
            require(all(diagonal[coordinates[s]]==diagonal[coordinates[complement[s]]] for s in states),
                    "reflection must preserve the transfer")
            T = sp.diag(*diagonal)
            trace = lambda op: sp.trace(P*op*P)
            numerators.append([trace(T**3),trace(T**2*U*T),trace(T*Y*T**2),trace(T*Y*T*U*T)])
            references.append(sp.trace(T**3))
            reduced.append(trace(T*Y*T*P*U*T))
            spectra.append(sorted(diagonal))
            require(numerators[-1][3]-reduced[-1]==trace(T*Y*T*Q*U*T),"missing exact excursion")
        require(spectra[0]==spectra[1] and references==[sp.Rational(393,4)]*2
                and numerators==[[96,96,0,40],[96,96,0,36]] and reduced==[32,32],
                "equal marginals, spectra and reflected reference must leave the mixed value ambiguous")
        family_reports.append({"omitted_complement_pair":[omitted,7-omitted],"inner_root_entries":inner,
                               "leaking_root_entries":leaks,"single_return_probe_rank":ranks[0],
                               "paired_overlap_probe_rank":ranks[1],"hidden_single_return_dimensions":32-ranks[0],
                               "reflection_paired_probe_ranks":list(reflected_ranks),
                               "isospectral_pair_contrast_in_paired_probe":"1",
                               "raw_ordered_numerators":[list(map(str,row)) for row in numerators],
                               "common_reflected_trace_reference":str(references[0]),
                               "connected_UY":["5/12","3/8"],"prematurely_compressed_connected_UY":"1/3",
                               "cartan_square":single_window_cartan_square(states,actions,kept),
                               "gauge_neutral_counterexample":gauge_neutral_control(
                                   states,coordinates,kept,outside,complement)})
    return {"status":"EXACT_FOUR_FAMILY_ORDERED_TRANSPORT","full_dimension":128,
            "projected_dimension":96,"mixed_root_operators":len(actions),"families":family_reports,
            "odd_e8_intermediate":odd_root_intermediate(),
            "full_cartan_square":single_window_cartan_square(states,actions,states),
            "scope":"D8 half-spinor in the E8 root chart; declared positive comparison transfers and marks",
            "native_r113_mark_identification_proved":False,"physical_solution":False}
