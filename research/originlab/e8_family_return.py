"""Compose the signed E8 currents before selecting three families."""

from collections import Counter
from fractions import Fraction

from sympy import Matrix, symbols

from affine_operators import (LatticeFock, MARKED_GAMMA, MARKED_P1, MARKED_P2,
                              ROOTS, add, dot, integral, lattice_frame, neg, require)
from relative_history import connected_history


def family_label(root):
    """Pair a fundamental family weight with its conjugate in the other chirality."""
    weight = root[5:]
    return neg(weight) if root[:5].count(-1) % 2 else weight


def mixed_histories(unweighted, attenuated):
    """Classify all couplings between equiprobable scalar and transfer marginals."""
    require(unweighted > attenuated > 0, "need two distinguishable positive reflected responses")
    wide,narrow = str(unweighted),str(attenuated)
    report = connected_history({
        "histories":["u0_unweighted","u0_attenuated","u2_unweighted","u2_attenuated"],
        "mark_coefficients":dict(Z0=[1,1,1,1],ZU=[0,0,2,2],ZY=[wide,narrow,wide,narrow],
                                 ZUY=[0,0,str(2*unweighted),str(2*attenuated)]),
        "relative_constraints":{
            "uniform_early_marker":{"coefficients":[0,0,2,2],"value":1},
            "uniform_transfer_choice":{"coefficients":[1,0,1,0],"value":"1/2"},
        },
    })
    limit = (unweighted-attenuated)/2
    require(report["vertices"]==2 and sorted(Fraction(row["connected_UY"]) for row in report["witnesses"])
            == [-limit,limit], "mixed history bounds do not exhaust the declared couplings")
    return {**report,"connected_interval":[str(-limit),str(limit)],
            "fixed_ZY":str((unweighted+attenuated)/2)}


def common_middle_comparison(matter, family_labels, groups, steps, hypercharge):
    """Test a positive middle operator shared by the signed currents."""
    charges = (Fraction(-1, 6),) * 3 + (Fraction(1, 4),) * 2
    mean_u = sum(charge**2 for charge in charges)
    targets = tuple((mean_u - 2*charge**2)/4 for charge in charges)
    dual = (Fraction(2, 11),) * 3 + (Fraction(5, 22),) * 2

    def orbit(root):
        signs = tuple(map(sum, (root[:3], root[3:5])))
        return min(signs, neg(signs))

    orbits = Counter(map(orbit, matter))
    flipped = {root: tuple((hypercharge[root] - 2*charge*root[i])**2
                           for i, charge in enumerate(charges)) for root in matter}
    require(all(sum(squares) == hypercharge[root]**2 + 4*mean_u
                for root, squares in flipped.items()), "five-coordinate square identity failed")
    certificates = [sum(weight*(square/4-target)
                        for weight, square, target in zip(dual, squares, targets))
                    for squares in flipped.values()]
    sharp_gap = Fraction(7, 528)
    require(sum(dual) == 1 and min(certificates) == sharp_gap,
            "common positive middle operator must have a signed residual of at least 7/528")

    # At fixed family coordinate, each arrival has exactly one predecessor
    # per carrier coordinate. Its family label is independent of that coordinate.
    for j in range(5, 8):
        masks = []
        for i in range(5):
            predecessors = {add(root, gamma): root for gamma in groups[i, j] for root in steps[gamma]}
            require(len(predecessors) == len(matter)
                    and all(hypercharge[root]**2 == flipped[target][i]
                            for target, root in predecessors.items()), "incorrect matter return")
            masks.append({target: family_label(root) for target, root in predecessors.items()})
        require(all(mask == masks[0] for mask in masks),
                "outer three-family selection must give a common arrival mask at fixed family coordinate")
        require(set(masks[0].values()) == set(family_labels), "family selection changed")

    def profile(root, key):
        return Fraction(len(matter), orbits[key]) if orbit(root) == key else Fraction(0)

    orbit_keys = sorted(orbits)
    equations = ({}, {})

    def observed(weight, omitted=None, *, record_equations=False):
        retained = [r for r in matter if family_label(r) != omitted]
        require(sum(hypercharge[r]**2 for r in retained) / len(retained) == mean_u,
                "the unweighted outer Higgs mark must preserve its reference value")
        retained_set = set(retained)
        error = Fraction(0)
        for (i, _), quartet in groups.items():
            for gamma in quartet:
                pairs = [(root, add(root, gamma)) for root in steps[gamma] & retained_set]
                y = sum(weight(i, target) for _, target in pairs)/len(retained)
                uy = sum(hypercharge[root]**2*weight(i, target) for root, target in pairs)/len(retained)
                require(y == Fraction(1, 4), "positive comparison must preserve every unmarked return")
                error = max(error, abs(uy-targets[i]))
                if record_equations:
                    unmarked, marked = Counter(), Counter()
                    for root, target in pairs:
                        key = orbit(target)
                        unmarked[key] += Fraction(1, len(retained))
                        marked[key] += hypercharge[root]**2 / len(retained)
                    equations[i >= 3][tuple(unmarked[key] for key in orbit_keys)] = y
                    equations[i >= 3][tuple(marked[key] for key in orbit_keys)] = uy
        return error

    optimum = observed(lambda _, root: Fraction(43, 88)*profile(root, (-3, -2))
                       + Fraction(45, 88)*profile(root, (-1, -2)))
    require(optimum == sharp_gap, "positive witness must attain the dual bound on all 60 currents")

    # Use one positive gauge-invariant diagonal block for each type of current.
    colour = dict(zip(orbit_keys, (Fraction(n, 67) for n in (167,127,39,39,63,39))))
    weak = dict(zip(orbit_keys, (Fraction(n, 13) for n in (4,4,4,28,4,28))))
    require(all(min(block.values()) > 0 for block in (colour,weak)),
            "each middle block must be strictly positive on every destination")

    def graded_middle(i, root):
        return (colour if i < 3 else weak)[orbit(root)]

    selections = (None,*sorted(family_labels))
    require(all(observed(graded_middle, omitted, record_equations=True) == 0 for omitted in selections),
            "the same two middle blocks must reproduce all split marks for every outer selection")
    # The scalar generators considered here depend only on these occupancy
    # orbits. Positive time rescaling preserves their common ordering.
    w = Matrix(symbols("w:6"))
    order_bounds = (Fraction(125, 114), Fraction(16, 9))
    ordered_differences = (w[1]-w[3], w[3]-w[1])
    residuals = (Fraction(7, 228)*w[0] + Fraction(55, 228)*w[2], (w[0]+w[1])/9)
    minimizers = ((Fraction(131, 76), Fraction(235, 228)),
                  (Fraction(16, 45), Fraction(104, 45)))
    for rows, difference, bound, residual, minimizer in zip(
            equations, ordered_differences, order_bounds, residuals, minimizers):
        coefficients = Matrix(list(rows))
        values = Matrix(list(rows.values()))
        solution, free = coefficients.gauss_jordan_solve(values)
        substitution = dict(zip(w, solution))
        require(len(free) == 2 and (difference-bound-residual).subs(substitution).expand() == 0,
                "opposite occupancy order must follow from every fitted positive block")
        minimizing = solution.subs(dict(zip(free, minimizer)))
        require(min(minimizing) >= 0 and difference.subs(dict(zip(w, minimizing))) == bound,
                "nonnegative occupancy weights must attain the order bound")

    return {"common_middle_minimax_gap_full_trace":str(sharp_gap),
            "common_middle_lower_bound_after_any_outer_family_omission":str(sharp_gap),
            "dual_weights_by_carrier_coordinate":list(map(str, dual)),
            "sharp_positive_semidefinite_witness_orbit_masses":["43/88", "45/88"],
            "carrier_dependent_positive_definite_witness": {
                "occupancy_orbits":[list(key) for key in sorted(orbits)],
                "colour_diagonal":list(map(str, colour.values())),
                "weak_diagonal":list(map(str, weak.values())),
                "currents_per_selection":sum(map(len, groups.values())),
                "outer_selections":len(selections),
            },
            "shared_scalar_generator_with_channel_time_rescaling": {
                "status":"EXCLUDED_ON_OCCUPANCY_ORBITS",
                "contrasting_orbits":[list(orbit_keys[i]) for i in (1, 3)],
                "colour_first_minus_second_minimum":str(order_bounds[0]),
                "weak_second_minus_first_minimum":str(order_bounds[1]),
            },
            "scope":"U=hypercharge squared, unweighted outer trace, each Y=1/4, "
                    "reflected returns through a common positive middle operator"}


def scalar_marker_comparison(fock,frame,matter,labels):
    """Compare one positive early insertion with two signed Cartan factors."""
    charges = {root:integral(frame,root) for root in matter}
    states = {root:fock.vacuum(charge) for root,charge in charges.items()}
    hypercharge = {root:Fraction(sum(a*b for a,b in zip((-2,-2,-2,3,3,0,0,0),root)),12)
                   for root in matter}
    groups = {}
    for root in ROOTS:
        support = tuple(i for i,value in enumerate(root) if value)
        if not root[0]%2 and len(support)==2 and support[0]<5<=support[1]:
            groups.setdefault(support,[]).append(root)
    require(len(groups)==15 and all(len(rows)==4 for rows in groups.values()),
            "the 60 signed bifundamental currents must form 15 quartets")

    steps = {}
    for quartet in groups.values():
        for gamma in quartet:
            operator = integral(frame,gamma)
            returns = set()
            for root,state in states.items():
                image = fock.vertex(operator,0,state)
                norm = fock.inner(image,image)
                target = add(root,gamma)
                require((norm==1 and target in states) if image else target not in states,
                        "bifundamental vertex norm must equal its matter root incidence")
                if image:
                    returns.add(root)
            steps[gamma] = returns
        require(all(sum(root in steps[gamma] for gamma in quartet)==1 for root in matter),
                "a signed quartet must partition the matter endpoints")

    readouts = []
    for omitted in (None,*sorted(labels)):
        retained = {root for root in matter if family_label(root)!=omitted}
        dimension = len(retained)
        u = sum(hypercharge[root]**2 for root in retained)/dimension
        differences = Counter()
        attenuation_floors = {"colour":set(),"weak":set()}
        for gamma,domain in steps.items():
            hits = retained & domain
            y = Fraction(len(hits),dimension)
            single = sum(hypercharge[root]**2 for root in hits)/dimension
            split = sum(hypercharge[root]*hypercharge[add(root,gamma)] for root in hits)/dimension
            weak = next(i for i in range(5) if gamma[i])>=3
            expected = Fraction(-1,32) if weak else Fraction(-1,72)
            require(u==Fraction(5,24) and y==Fraction(1,4) and single-u*y==0 and split-u*y==expected,
                    "the scalar square and split Cartan marker must separate on every signed current")
            differences[str(split-u*y)] += 1
            # For 0 <= Y_eff <= Y, the trace of the centered insertion is
            # minimized by deleting every positive centered diagonal entry.
            floor = -sum(max(hypercharge[root]**2-u,0) for root in hits)/dimension
            attenuation_floors["weak" if weak else "colour"].add(floor)
        readouts.append({"omitted_family":list(omitted) if omitted else None,"matter_dimension":dimension,
                         "ZU":str(u),"ZY":"1/4","ZUY_with_single_square":"5/96",
                         "connected_with_single_square":"0","split_connected_counts":dict(differences),
                         "attenuation_only_connected_floors":{
                             key:list(map(str,sorted(values))) for key,values in attenuation_floors.items()
                         }})

    require(readouts[0]["attenuation_only_connected_floors"]["weak"]==["-19/768"]
            and all(Fraction(value)>Fraction(-1,32) for row in readouts
                    for value in row["attenuation_only_connected_floors"]["weak"]),
            "a positive contraction of the fixed return must miss the weak split value")

    return {"status":"EXCLUDED_SINGLE_MARKER_MATCH_ON_FIXED_RETURN_MENU",
            "vertex_norms_checked":len(matter)*len(steps),"signed_quartets":len(groups),
            "sum_of_returns_per_quartet":"identity on matter endpoints",
            "sum_of_connected_responses_per_quartet_for_any_common_early_marker":"0",
            "optimal_worst_absolute_residual_to_split_cartan":"1/32","readouts":readouts,
            "native_probe":{"word":["U-minus-mean","sum of four signed returns at fixed carrier and family coordinate"],
                            "split_target_by_carrier":{"colour":"-1/18","weak":"-1/8"},
                            "normalization":"shared selected empty response"},
            "weak_split_match_by_attenuation_only":False,
            "positive_middle_comparison":common_middle_comparison(matter, labels, groups, steps, hypercharge),
            "scope":"fixed unweighted returns and one early insertion common to the four signs"}


def controls():
    """Exact level-one grade-one calculation with a declared finite reference."""
    frame,_ = lattice_frame()
    fock = LatticeFock([[dot(a,b) for b in frame] for a in frame])
    left,middle,right = (integral(frame,r) for r in (neg(MARKED_P1),MARKED_GAMMA,MARKED_P2))
    matter = [r for r in ROOTS if r[0] % 2]
    even = {integral(frame,r) for r in ROOTS if not r[0] % 2}
    matter_charges = {integral(frame,r) for r in matter}
    labels = Counter(family_label(r) for r in matter)
    require(len(labels)==4 and set(labels.values())=={32} and len(even)==112,
            "the D8 grade-one sector must contain four matter families and 112 even roots")
    # Degree one also contains eight zero-charge oscillator states (Cartan).
    cartan = {next(iter(state)) for state in fock.basis((0,)*8,1)}
    require(len(cartan)==8, "grade-one Cartan basis changed")
    transfer_factors = (Fraction(1), Fraction(1,2))
    returns = {}
    dropped_cartan = []
    counts = Counter()
    for root in matter:
        state = fock.vacuum(integral(frame,root))
        early = fock.vertex(right,0,state)
        later = fock.vertex(middle,0,early)
        require(all((key[0] in even and not key[1]) or key in cartan
                    for step in (early,later) for key in step), "word left the even grade-one sector")
        require(all(key[0] not in matter_charges for key in early),
                "the first odd insertion must vanish under an intervening matter projection")
        counts["after_first"] += bool(early)
        counts["after_middle"] += bool(later)
        closed = fock.vertex(left,0,later)
        coefficient = closed.get(next(iter(state)),Fraction(0))
        require(closed==fock.combine((coefficient,state)), "neutral word must be diagonal on the matter root")
        values = []
        for factor in transfer_factors:
            weighted = fock.vertex(middle,0,fock.combine((factor,early)))
            weighted = fock.vertex(left,0,fock.combine((factor,weighted)))
            values.append(fock.inner(weighted,weighted))
        require(values==[coefficient**2,coefficient**2/16],
                "two positive even-sector transfers must scale the reflected response by 1/16")
        returns[root] = (coefficient,values)
        root_only = fock.vertex(middle,0,{key:v for key,v in early.items() if key not in cartan})
        root_only = fock.vertex(left,0,{key:v for key,v in root_only.items() if key not in cartan})
        if root_only!=closed:
            require(closed and not root_only, "unexpected root-only truncation residual")
            dropped_cartan.append({"start_root":list(root),
                                   "cartan_after_step":1 if any(key in cartan for key in early) else 2})

    def readout(retained):
        reference = sum(fock.inner(fock.vacuum(integral(frame,r)), fock.vacuum(integral(frame,r)))
                        for r in retained)
        word = sum(returns[r][0] for r in retained)
        bare,attenuated = (sum(returns[r][1][i] for r in retained)/reference for i in (0,1))
        return {"matter_dimension":len(retained), "word_rank":sum(bool(returns[r][0]) for r in retained),
                "signed_trace":str(word),"reflected_empty_trace":str(reference),
                "reflected_return":str(bare),"attenuated_reflected_return":str(attenuated),
                "reflected_mixed_histories":mixed_histories(bare,attenuated)}

    selected = [{"omitted_family":list(label),
                 **readout([r for r in matter if family_label(r)!=label])} for label in sorted(labels)]
    full = readout(matter)
    require(full["word_rank"]==14 and full["signed_trace"]=="-14" and len(dropped_cartan)==2
            and sorted(row["word_rank"] for row in selected)==[9,9,12,12],
            "marked family return census changed")
    return {"status":"EXACT_E8_FAMILY_RETURN", "ambient_dimension":248,
            "word":"J^(-p1)_0 J^gamma_0 J^p2_0", "ordered_prefix_support":dict(counts),
            "full_matter":full,"three_family_readouts":selected,
            "premature_matter_projection_rank":0,"root_only_word_rank":12,
            "returns_requiring_Cartan":dropped_cartan,
            "reflected_reference":"trace of the identity on the retained matter endpoints",
            "transfer_control":"equal history masses; even-sector transfer eigenvalue 1 or 1/2 per gap; identity on matter",
            "scalar_marker_comparison":scalar_marker_comparison(fock,frame,matter,labels),
            "scope":"level-one grade one with finite reflected trace; Y is the squared return in this control",
            "physical_source_identity_proved":False}
