"""Colour Ward relations and the neutral route to the marked affine word."""

from collections import Counter
from functools import lru_cache
from itertools import combinations, combinations_with_replacement
from math import prod

import sympy as sp

from affine_operators import (LatticeFock, MARKED_GAMMA, MARKED_P1, MARKED_P2, ROOTS, ROOT_SET,
                              add, charge_signature, closure, components, dot, integral, lattice_frame,
                              marked_word, neg, simple_roots, unit)
from reconstruction import require


SECTORS = tuple(tuple(sorted(root for root in ROOTS if charge_signature(root)[1]==c)) for c in range(3))
COLOUR_SIMPLE = (add(unit(0),neg(unit(1))),add(unit(1),neg(unit(2))))
TRIPLET_WEIGHTS = {(1,0),(-1,1),(0,-1)}
COLOUR_WEIGHTS = ({(0,0)},TRIPLET_WEIGHTS,{(-a,-b) for a,b in TRIPLET_WEIGHTS})
E6 = tuple(root for root in SECTORS[0] if all(dot(root,c)==0 for c in COLOUR_SIMPLE))
HYPERCHARGE = sp.Matrix([-2,-2,-2,3,3,0,0,0])/6
NEUTRAL_RAISE = add(unit(3),unit(5))


@lru_cache(maxsize=1)
def _lattice():
    frame,_ = lattice_frame()
    fock = LatticeFock([[dot(a,b) for b in frame] for a in frame])
    return frame,fock,{r:integral(frame,r) for r in ROOTS}


def _phase(root):
    return 1 if root > neg(root) else -1


def _bracket():
    _,fock,coords = _lattice()
    require(all(fock.sign(coords[r],coords[neg(r)])==-1 for r in ROOTS),
            "opposite root phases must give positive adjoints")

    @lru_cache(maxsize=None)
    def bracket(a,b):
        total = add(a,b)
        if total not in ROOT_SET:
            return None,0
        return total,_phase(a)*_phase(b)*_phase(total)*fock.sign(coords[a],coords[b])
    return bracket


def _component(graph, seed):
    seen = {seed}
    pending = [seed]
    while pending:
        fresh = graph[pending.pop()]-seen
        pending.extend(fresh)
        seen.update(fresh)
    return seen


def _component_sizes(graph):
    sizes,remaining = [],set(graph)
    while remaining:
        component = _component(graph,min(remaining))
        sizes.append(len(component))
        remaining -= component
    return sorted(sizes)


def _is_pair_root(root):
    return all(x%2==0 for x in root)


def _ward_class(neutral, charged, bracket):
    pairs = list(combinations(charged,2))
    nodes = {ab for ab in pairs if bracket(*ab)[0] is not None}
    graph = {ab:set() for ab in nodes}
    orbit = {r:set() for r in charged}
    rows = 0
    for h in neutral:
        moves = {r:bracket(h,r) for r in charged}
        for r,(s,_) in moves.items():
            if s is not None:
                orbit[r].add(s)
        for a,b in pairs:
            terms = []
            t,s = bracket(a,b)
            if t is not None:
                _,v = bracket(h,t)
                if v:
                    terms.append(((a,b),s*v))
            t,s = moves[a]
            if t is not None:
                _,v = bracket(t,b)
                if v:
                    terms.append((tuple(sorted((t,b))),-s*v))
            t,s = moves[b]
            if t is not None:
                _,v = bracket(a,t)
                if v:
                    terms.append((tuple(sorted((a,t))),-s*v))
            if terms:
                require(len(terms)==2 and sum(value for _,value in terms)==0,
                        "neutral Ward identity must equate two normalized coefficients")
                i,j = (ab for ab,_ in terms)
                require(i!=j,"the Ward row must constrain distinct coefficients")
                graph[i].add(j)
                graph[j].add(i)
                rows += 1
    a,b = min(nodes)
    c = add(a,b)
    target,s = bracket(neg(b),c)
    families = _component_sizes(graph)
    require(target==a and bracket(a,b)[1]*s==-1 and dot(a,b)==-1 and
            bracket(neg(b),a)[0] is None and
            all(bracket(neg(x),neg(y))[1]==-bracket(x,y)[1] for x,y in nodes),
            "adjoints and Jacobi must fix the cubic scale modulus")
    triangles = {tuple(sorted((x,y,neg(add(x,y))))) for x,y in nodes}
    for x,y,z in triangles:
        require(bracket(x,y)[1]==bracket(y,z)[1]==bracket(z,x)[1],
                "the reference cubic must be cyclic with respect to the invariant form")
        for i,j in (((x,y),(y,z)),((y,z),(x,z))):
            graph[i].add(j)
            graph[j].add(i)
    return {"unordered_bracket_pairs":len(nodes),"oriented_cubic_entries":2*len(nodes),
            "unordered_cubic_triangles":len(triangles),"ward_rows_checked":len(neutral)*len(pairs),
            "nonzero_ward_rows":rows,"relation_rank":len(nodes)-len(families),
            "relative_complex_parameters":len(families),"ward_component_sizes":families,
            "cyclic_ward_component_sizes":_component_sizes(graph),
            "neutral_root_orbit_sizes":_component_sizes(orbit),
            "jacobi_witness":[list(a),list(b),list(neg(b))]}


def _root_strings(bracket):
    result = []
    lowering = neg(MARKED_P1)
    for i,roots in enumerate(SECTORS):
        weights = Counter(dot(MARKED_P1,r) for r in roots)
        if i==0:
            weights[0] += 8
        triplets,doublets = weights[2],weights[1]
        singlets = weights[0]-triplets
        transitions = Counter((_is_pair_root(r),_is_pair_root(s)) for r in roots
                              if (s:=bracket(lowering,r)[0]) is not None)
        require(weights[-1]==doublets and weights[-2]==triplets and
                sum(transitions.values())==doublets and
                transitions[True,True]==transitions[False,False]==0,
                "the marked lowering strings changed")
        result.append({"colour":i,"singlets":singlets,"doublets":doublets,"triplets":triplets,
                       "lowering_rank":doublets+2*triplets,"lowering_square_rank":triplets,
                       "D8_to_spinor":transitions[True,False],"spinor_to_D8":transitions[False,True],
                       "casimir_trace":str(sp.Rational(3,2)*doublets+6*triplets)})
    return result


@lru_cache(maxsize=1)
def _e6_bounds():
    simple = simple_roots(E6)
    gram = sp.Matrix([[dot(a,b) for b in simple] for a in simple])
    inverse = gram.inv()
    positives = []
    for root in E6:
        coords = tuple(inverse*sp.Matrix([dot(r,root) for r in simple]))
        require(all(x.is_Integer for x in coords),"E6 roots must have integral simple coordinates")
        if all(x>=0 for x in coords):
            positives.append(coords)
    require(len(positives)==36 and len(simple)==6 and MARKED_P1 in E6,"wrong E6 complement")

    def dimension(weight):
        return int(prod((sum(root)+sum(x*y for x,y in zip(root,weight)))/sum(root) for root in positives))

    fundamental = list(zip(*sp.eye(6).tolist()))
    dims = [dimension(w) for w in fundamental]
    y = (sp.Matrix(simple).T/2).gauss_jordan_solve(HYPERCHARGE)[0]
    require(all((6*x).is_Integer for x in y),"the E6 hypercharge centre needs order three")
    grades = [int(6*x)%3 for x in y]
    small = [w for w,d in zip(fundamental,dims) if d<78]
    require(sorted(dims)==[27,27,78,351,351,2925] and
            [g for g,d in zip(grades,dims) if d<78]==[1,2] and
            all(dimension(add(a,b))>=78 for a,b in combinations_with_replacement(small,2)),
            "the smallest nontrivial centre-neutral E6 representation must have dimension 78")
    return dims,grades


def _colour_weak_profiles(by_weight):
    profiles = [Counter() for _ in COLOUR_WEIGHTS]
    for (y,w,_,_),n in by_weight.items():
        profiles[y%3][y,w] += n
    for profile,weights in zip(profiles,COLOUR_WEIGHTS):
        for (y,w),n in profile.items():
            require({(a,b) for q,t,a,b in by_weight if (q,t)==(y,w)}==weights and
                    all(by_weight[y,w,a,b]*len(weights)==n for a,b in weights),
                    "a charge eigenspace does not have the fixed physical colour weights")
    return profiles


def _e6_weight_profiles():
    # The two minuscule modules are the coloured parts of the adjoint
    # E8 representation with the physical colour factor removed.
    by_weight = Counter((int(6*charge_signature(r)[2]),(r[3]-r[4])//2,
                         *(dot(c,r) for c in COLOUR_SIMPLE)) for r in (*E6,*SECTORS[1],*SECTORS[2]))
    by_weight[0,0,0,0] += 6
    profiles = [Counter({weight:n//len(colour) for weight,n in p.items()})
                for p,colour in zip(_colour_weak_profiles(by_weight),COLOUR_WEIGHTS)]
    require([sum(p.values()) for p in profiles]==[78,27,27],"wrong E6 charge profiles")
    return profiles


def _charge_counts(weights):
    counts = Counter()
    for (y,_),n in weights.items():
        counts[y] += n
    return dict(sorted(counts.items()))


def _weak_highest_weights(weights):
    # For a compact SU(2) action each spin contributes one state per weight.
    return Counter({(y,w):n-weights[y,w+2] for (y,w),n in weights.items()
                    if w>=0 and n!=weights[y,w+2]})


def _exterior_traces(charges, profiles):
    traces = [sum(y*n for y,n in p.items())//6 for p in profiles]
    require(traces[0]==0 and traces[1]==-traces[2],"native centre traces must be conjugates")
    neutral,positive,negative = [sum(p.values()) for p in profiles]
    require(positive==negative,"native centre dimensions must be conjugates")
    t = sp.Symbol("t")
    # Differentiate the exterior character in Y, then project to colour one
    # at the two primitive roots of unity. Colour two has the opposite trace.
    trace_poly = traces[1]*sp.Poly(t,t)*sp.Poly(1+t,t)**neutral*sp.Poly(1-t+t**2,t)**(positive-1)
    wedge_two = [0,0,0]
    for a,b in combinations(charges,2):
        wedge_two[(a+b)%3] += a+b
    two = int(trace_poly.nth(2))
    require(wedge_two==[0,6*two,-6*two],"exterior charge projection failed its pair census")
    total = int(trace_poly.eval(1))
    even = (total+int(trace_poly.eval(-1)))//2
    return {"native_Y_traces_by_colour":traces,
            "native_exterior_two_Y_traces_by_colour":[0,two,-two],
            "native_full_exterior_Y_traces_by_colour":[0,total,-total],
            "native_even_exterior_Y_traces_by_colour":[0,even,-even],
            "native_odd_exterior_Y_traces_by_colour":[0,total-even,even-total],
            "exterior_degrees_with_zero_colour_trace":[n for n in range(len(charges)+1) if not trace_poly.nth(n)],
            "exterior_trace_generating_function_colour_one":f"{traces[1]}*t*(1+t)^{neutral}*(1-t+t^2)^{positive-1}"}


def native_carrier_check(gauge_weights):
    """Audit joint 6Y, 2T3 and the two colour Cartan weights on the native fibre."""
    by_weight = Counter(tuple(map(int,row)) for row in gauge_weights)
    weak_profiles = _colour_weak_profiles(by_weight)
    profiles = [_charge_counts(p) for p in weak_profiles]
    require(profiles==[{-6:3,-3:6,0:6,3:6,6:3},{-2:9,1:18,4:9},{-4:9,-1:18,2:9}],
            "the native full96 hypercharge census changed")
    dims,grades = _e6_bounds()
    target = _e6_weight_profiles()
    multiplicities = [Counter({weight:n//len(colour) for weight,n in p.items()})
                      for p,colour in zip(weak_profiles,COLOUR_WEIGHTS)]
    target_highest = list(map(_weak_highest_weights,target))
    native_highest = list(map(_weak_highest_weights,multiplicities))
    for old,new in zip(native_highest,target_highest):
        require(all(n>0 for n in (*old.values(),*new.values())) and old<=new,
                "the candidate E6 modules must carry all native weak representations at each hypercharge")
    # Weyl's formula is monotone in the dominant weight coordinates.
    # The singlet colour contains nonzero Y and requires a centre-neutral
    # E6 irrep. Each triplet colour requires a nontrivial E6 irrep.
    minimum = min(d for d,g in zip(dims,grades) if g==0)+6*min(dims)
    require(minimum==sum(sum(p.values())*c for p,c in zip(target,(1,3,3))),
            "the candidate modules must attain the lower dimension bound")
    charges = [y for p in profiles for y,n in p.items() for _ in range(n)]
    return {"native_one_particle_dimensions":[24,36,36],"native_colour_multiplicity_dimensions":[24,12,12],
            **_exterior_traces(charges,profiles),
            "neutral_E6_on_unchanged_Y_excluded":True,
            "missing_terminal_6Y_in_two_neutral_current_word":5,
            "minimum_charge_preserving_neutral_completion_dimension":minimum,
            "completion_profile_by_6Y":list(map(_charge_counts,target)),
            "completion_weak_highest_weights_per_colour_basis":[
                [{"6Y":y,"2T3_highest":w,"multiplicity":n} for (y,w),n in sorted(p.items())]
                for p in target_highest],
            "completion_modules":["E6 adjoint times colour singlet","E6 27 times colour 3","E6 dual 27 times dual colour 3"],
            "scope":"finite E6 and physical colour SU(3), with pinned native SU(2), Y and colour actions; time, original weighted source and field identifications are open"}


def _marked_reordering(bracket, size):
    frame,fock,coords = _lattice()
    a,b,c = neg(MARKED_P1),MARKED_GAMMA,MARKED_P2
    aa,bb,cc = (coords[r] for r in (a,b,c))
    d = coords[NEUTRAL_RAISE]
    empty = fock.vacuum()
    early = fock.vertex(cc,-1,empty)
    interim = fock.vertex(aa,0,early)
    terminal = add(a,c)
    stop = add(terminal,NEUTRAL_RAISE)
    second = fock.vertex(d,0,interim)
    require(terminal==neg(b) and all(charge_signature(r)[1]==0 for r in (a,NEUTRAL_RAISE)) and
            charge_signature(terminal)[1:]==(2,sp.Rational(1,3)) and
            charge_signature(stop)[1:]==(2,sp.Rational(5,6)) and
            fock.inner(second,second)==1 and not fock.vertex(d,0,early),
            "the two neutral current steps must reach the missing native terminal charge")
    for m in range(1,size+1):
        state = fock.vertex(cc,-m,empty)
        ordered = fock.vertex(d,0,fock.vertex(aa,0,state))
        require(fock.inner(ordered,ordered)==m and not fock.vertex(d,0,state),
                "the ordered Ward return must preserve the affine mode norm")
    length = marked_word(fock,frame,size)
    coefficient = bracket(a,b)[1]
    signs = _phase(a)*_phase(b)*_phase(c)
    for (m,n),value in length.items():
        state = fock.vertex(cc,-n,empty)
        swapped = fock.vertex(bb,m,fock.vertex(aa,n-m,state)).get(((0,)*8,()),0)
        require(signs*value==coefficient*min(m,n) and swapped==-value,
                "both orderings must obey their signed affine Ward coefficients")
    return {"tested_mode_pairs":len(length),"marked_coefficient":coefficient,
            "identity":"W(m,n) = -<J^(-gamma)_(-m) Omega, J^(-p1)_(n-m) J^p2_(-n) Omega>",
            "early_to_terminal_colour":[2,2,2],"early_to_terminal_Y":["1/3","1/3","5/6"],
            "two_neutral_roots_doubled":[list(a),list(NEUTRAL_RAISE)],
            "terminal_doubled_root":list(stop),"reverse_two_current_word_norm_squared":0,
            "forward_two_current_word_norm_squared_at_mode_m":"m"}


def _positive_mode_test():
    _,fock,coords = _lattice()
    empty = fock.vacuum()
    neutral = fock.vertex(coords[MARKED_P1],-1,empty)
    lowering = fock.vertex(coords[neg(MARKED_P1)],1,neutral)
    require(_phase(neg(MARKED_P1))*_phase(MARKED_P1)==-1 and
            lowering=={state:-value for state,value in empty.items()},
            "the positive neutral mode must return the level-one vacuum")
    return {"neutral_current_positive_mode_vacuum_coefficient":1,
            "charged_grade_one_annihilated_by_neutral_positive_modes":True,
            "charged_reason":"neutral modes preserve colour and lower affine grade; the only grade-zero state is the neutral vacuum"}


def folded_bridge(fold_cases):
    result = []
    terminal = add(neg(MARKED_GAMMA),NEUTRAL_RAISE)
    for case in fold_cases:
        d = sp.Matrix(list(map(sp.Rational,case["folded_cartan"])))
        early,horizontal_late,end = (d.dot(sp.Matrix(r)/2) for r in
                                     (MARKED_P2,neg(MARKED_GAMMA),terminal))
        old_late,old_early = map(sp.Rational,case["end_exponents_in_K_time"])
        require(early==old_early and horizontal_late==early-old_late and end==0,
                "the original folded clock must transport the two charge shifts")
        difference = old_late-horizontal_late
        result.append({"quarter":case["quarter"],"two_neutral_word_K_shifts":list(map(str,(early,horizontal_late,end))),
                       "marked_over_negative_horizontal_left_exponent":str(difference)})
    require([row["marked_over_negative_horizontal_left_exponent"] for row in result]==["-3/2","-5/4"],
            "the fixed folded cyclic Ward exponent changed")
    return {"cases":result,"relation":"W(x,y)=-x^difference H(x,y)",
            "affine_two_neutral_word_normalized_terminal_heat_at_mode_m":"exp(-t*m)",
            "scope":"fixed folded affine chart, identical signed radial phases and affine vacuum; original weighted source must bind both preparations and the neutral modes"}


def controls():
    from marked_source import controls as marked_source_controls

    neutral = SECTORS[0]
    pair_counts = [sum(map(_is_pair_root,roots)) for roots in SECTORS]
    spinor_counts = [len(roots)-n for roots,n in zip(SECTORS,pair_counts)]
    d8 = [r for r in neutral if _is_pair_root(r)]
    require([len(r) for r in SECTORS]==[78,81,81] and components(neutral)==["A2","E6"]
            and components(d8)==["A2","D5"],"colour-neutral root systems changed")
    su3_c = sp.Rational(8,1+3)
    free_c = sp.Rational(45,1+8)+su3_c+1
    extended_c = sp.Rational(78,1+12)+su3_c
    require(free_c==extended_c==8,"the neutral Sugawara central charges must agree")
    require(closure([*d8,MARKED_P1])==set(neutral) and closure([*neutral,MARKED_P2])==ROOT_SET,
            "the two missing source steps must close the neutral and charged sectors in order")
    bracket = _bracket()
    dims,grades = _e6_bounds()
    d8_ward = [_ward_class(d8,coloured,bracket) for coloured in SECTORS[1:]]
    charged_ward = [_ward_class(neutral,coloured,bracket) for coloured in SECTORS[1:]]
    require(all(r["ward_component_sizes"]==[30,60,240,480] and
                r["cyclic_ward_component_sizes"]==[90,720] for r in d8_ward) and
            all(r["ward_component_sizes"]==[810] and r["neutral_root_orbit_sizes"]==[81]
                for r in charged_ward),"neutral spinor completion must connect all charged coefficients")
    return {"status":"EXACT_AFFINE_COLOUR_WARD_CHECKS","grade_one_dimensions":[86,81,81],
            "D8_pair_roots_by_colour":pair_counts,"spinor_roots_by_colour":spinor_counts,
            "positive_modes":_positive_mode_test(),
            "neutral_level_one_central_charges":{"D5_A2_U1":int(free_c),"E6_A2":int(extended_c)},
            "neutral_root_counts_along_extension":[len(d8),len(neutral),len(ROOTS)],
            "D8_neutral_ward":d8_ward,"charged_ward":charged_ward,"cubic_scale_modulus_squared":1,
            "phase_freedom":"conjugate rephasing of the two charged sectors changes the cubic scale by a cube",
            "neutral_marked_lowering":_root_strings(bracket),"marked_ordering":_marked_reordering(bracket,4),
            "E6_fundamental_dimensions":dims,"E6_fundamental_centre_grades":grades,
            "E6_profiles_by_6Y":list(map(_charge_counts,_e6_weight_profiles())),
            "marked_source":marked_source_controls(
                tuple(r["neutral_root_orbit_sizes"][0] for r in charged_ward), COLOUR_SIMPLE, bracket),
            "physical_source_identity_proved":False,"original_values":None}
