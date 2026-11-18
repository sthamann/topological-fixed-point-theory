"""Replay the weak word in the original coupled W Hamiltonian."""

import ast
from collections import Counter, defaultdict, deque
from hashlib import sha256
from itertools import combinations, product
from math import factorial, gcd, prod
from pathlib import Path

import sympy as sp

from affine_operators import LatticeFock, ROOT_SET, add, charge_signature, dot, integral, lattice_frame, neg
from affine_time import radial_word
from reconstruction import require
from verify import digest
from weak_source import triangle_components


CONTRACT = Path("experiments/theory-contracts/universalraum-keyD-native-instruments-20260915")
PAIR_DRESSING = CONTRACT/"check_gap_closure.py"
PINNED_INPUTS = {
    CONTRACT/"native_source.py":"380577f85d2afa8b50f91a0991769eaf8267ee5c09f327d87f17e1b4c0af1672",
    Path("experiments/theory-contracts/source-native-retarded-selection-20260922/PROOF.md"):
        "e68bd415f03b83bfdd17d819bea40e1432a0cb70e354c0c51ee1c2ad3595ec1c",
}
INPUTS = (*PINNED_INPUTS,PAIR_DRESSING)
DRESSED_WORD_TIME_RANK = 4


def _load_original(root):
    source = (root/CONTRACT/"native_source.py").resolve()
    tree = ast.parse(source.read_text(),filename=str(source))
    # Execute the original W, J and sparse C3 constructors and the original
    # CAR helpers. The unrelated dense C3 archival hash would exceed 1 GB.
    last = next(i for i,node in enumerate(tree.body) if isinstance(node,ast.Assign) and
                any(isinstance(t,ast.Name) and t.id=="_GAMMA" for t in node.targets))
    program = ast.Module(body=tree.body[:last],type_ignores=[])
    namespace = {"__name__":"originlab_native_weak_w","__file__":str(source)}
    exec(compile(program,str(source),"exec"),namespace)
    return namespace


def _pair_creators(w):
    from numpy import flatnonzero

    return tuple(tuple((w["PAIRS"][column],int(row[column])) for column in flatnonzero(row))
                 for row in w["W"])


def _boson_creator(state, label):
    return {(fermions,tuple(sorted((*bosons,label)))):value
            for (fermions,bosons),value in state.items()}


def _original_pair_dressing(root, w):
    source = root/PAIR_DRESSING
    text = source.read_text()
    value, = [node.value for node in ast.parse(text).body if isinstance(node,ast.Assign)
              and any(isinstance(name,ast.Name) and name.id=="bplus" for name in node.targets)]
    expression = ast.get_source_segment(text,value)
    require(sha256(expression.encode()).hexdigest()=="7050274c48b145d1e87b31e16b859f18992c0d3daae8a479d33f7f5301a35090",
            "review the changed original neutral-pair operator")
    dressing = eval(compile(ast.Expression(value),str(source),"eval"),
                    {"bar":w["BAR"],"signs":w["ETA"]})
    require(len(dressing)==30 and all(not annihilators and tag is None for _,annihilators,tag in dressing),
            "original pair dressing must have thirty purely bosonic creator pairs")
    return [(bosons,int(coefficient)) for (bosons,_,_),coefficient in dressing.items()]


def _neutral_pair(state, dressing):
    out = defaultdict(int)
    for (fermions,bosons),value in state.items():
        for pair,coefficient in dressing:
            out[fermions,tuple(sorted((*bosons,*pair)))] += coefficient*value
    return out


def _lower_boson_raise_pair(state, creators, action):
    """T_- on monomials of boson creators and normalized CAR occupation kets."""
    out = defaultdict(int)
    for (mask,bosons),value in state.items():
        for label,count in Counter(bosons).items():
            index = bosons.index(label)
            reduced = bosons[:index]+bosons[index+1:]
            for pair,coefficient in creators[label]:
                result = action(mask,pair,dagger=True)
                if result is not None:
                    target,sign = result
                    out[target,reduced] += count*coefficient*sign*value
    return out


def _inner(left, right):
    # Bosonic creation monomials here are unnormalized. Their exact squared
    # norms are products of factorials. All used coefficients are real.
    if len(left)>len(right):
        left,right = right,left
    return sum(coeff*right[key]*prod(factorial(n) for n in Counter(key[1]).values())
               for key,coeff in left.items() if coeff and key in right)


def _raise_boson_lower_pair(state, pairs, action):
    out = defaultdict(int)
    for (mask,bosons),value in state.items():
        occupied = []
        remaining = mask
        while remaining:
            lowest = remaining & -remaining
            occupied.append(lowest.bit_length()-1)
            remaining ^= lowest
        for pair in combinations(occupied,2):
            vertex = pairs.get(pair)
            if vertex is not None:
                label,coefficient = vertex
                target,sign = action(mask,pair)
                out[target,tuple(sorted((*bosons,label)))] += coefficient*sign*value
    return out


def _primitive(vector):
    nonzero = {key:value for key,value in vector.items() if value}
    factor = gcd(*nonzero.values())
    return {key:value//factor for key,value in nonzero.items()} if factor>1 else nonzero


def _reachable(seed, lower, raise_, dimension):
    # Eliminate occupation amplitudes over Q separately at each boson number.
    # The queue closes only after applying both cubic terms to every basis vector.
    basis, pivots = [], defaultdict(list)
    pending = deque([seed])
    while pending:
        vector = _primitive(pending.popleft())
        if not vector:
            continue
        degree = len(next(iter(vector))[1])
        for pivot,other in pivots[degree]:
            if pivot not in vector:
                continue
            a,b = other[pivot],vector[pivot]
            common = gcd(a,b)
            residue = defaultdict(int,{key:value*(a//common) for key,value in vector.items()})
            for key,value in other.items():
                residue[key] -= value*(b//common)
            vector = _primitive(residue)
        if not vector:
            continue
        require(len(basis)<dimension,"original W cyclic space exceeded its verified rank")
        pivots[degree].append((min(vector),vector))
        lowered = lower(vector)
        basis.append((degree,vector,lowered))
        pending.extend((lowered,raise_(vector)))
    require(len(basis)==dimension and basis[0][1]==seed,"original W cyclic rank or seed changed")
    return basis


def _time_blocks(basis):
    degrees = [n for n,_,_ in basis]
    metric = sp.Matrix([[_inner(v,u) if n==m else 0 for m,u,_ in basis] for n,v,_ in basis])
    lower_forms = sp.Matrix([[_inner(v,lower_u) if n==m-1 else 0 for m,_,lower_u in basis]
                             for n,v,_ in basis])
    inverse = metric.inv()
    lower = inverse*lower_forms
    return degrees,metric,lower,inverse*lower.T*metric


def _time_reduction(incoming, cap, middle, w, creators, coefficient):
    def lower(vector):
        return _lower_boson_raise_pair(vector,creators,w["pair_action"])

    def raise_(vector):
        return _raise_boson_lower_pair(vector,w["_lookup"],w["pair_action"])

    right,left = (_reachable(seed,lower,raise_,size)
                  for seed,size in ((incoming,DRESSED_WORD_TIME_RANK),(cap,13)))
    raised = [(n+1,_boson_creator(v,middle)) for n,v,_ in right]
    middle_form = sp.Matrix([[_inner(v,u) if n==m else 0 for m,u in raised] for n,v,_ in left])
    delta,coupling = sp.symbols("Delta abs_g",real=True)

    def rows(matrix):
        return [[str(value) for value in row] for row in matrix.tolist()]

    def block(basis):
        degrees,metric,down,up = _time_blocks(basis)
        generator = delta*sp.diag(*degrees)+coupling*(up+down)
        powers = [sp.eye(len(basis))[:,0]]
        for _ in range(3):
            powers.append((generator*powers[-1]).applyfunc(sp.expand))
        return sp.Matrix.hstack(*powers),{
            "boson_numbers":degrees,"gram":rows(metric),"lowering_matrix":rows(down),
            "raising_matrix":rows(up),
        }

    rp,right_block = block(right)
    lp,left_block = block(left)
    jet = (lp.T*middle_form*rp).applyfunc(sp.expand)
    det = sp.factor(jet.det(method="domain-ge"))
    require(jet[0,0]==31*coefficient and
            sp.expand(jet[:2,:2].det())==14198*coupling**2 and det==22096281600*coupling**12,
            "original W full two-time rank certificate changed")
    norm = _inner(incoming,incoming)*_inner(cap,cap)
    # N_b phase rotations remove arg(g). The phase of b† cancels the
    # difference between the seed boson numbers three and two.
    return {"right":right_block,"left":left_block,"middle_form":rows(middle_form),
            "convention":"T_plus is the Gram adjoint of T_minus. H=Delta*N_b+abs_g*(T_plus+T_minus). "
                         "The word contracts the two heat-evolved coordinate e_0 vectors through the middle form.",
            "representative_bracket_sign":coefficient,
            "unnormalized_moment_minor_orders_zero_to_three":str(det),
            "normalized_moment_minor_orders_zero_to_three":str(sp.factor(det/norm**2)),
            "mixed_time_rank_at_nonzero_g":DRESSED_WORD_TIME_RANK,"mixed_time_rank_at_zero_g":1}


def _affine_rank_check():
    # Divide the signed original word by nonvanishing row and column
    # factors. The remaining Cauchy kernel has full rank for every size
    # on distinct positive widths. Verify a short exact grid as a control.
    widths = [sp.Rational(1,i) for i in (2,3,4,5,6)]
    checks = []
    for n in (2,3,4,5):
        signed = sp.Matrix([[radial_word(x,y) for y in widths[:n]] for x in widths[:n]])
        divided = sp.Matrix([[1/(1+x*y) for y in widths[:n]] for x in widths[:n]])
        require(divided.det()!=0 and sp.cancel(signed.det(),extension=sp.I)!=0,
                "signed affine Cauchy kernel lost rank")
        checks.append(n)
    return checks


def _sign_solution(equations):
    pivots = {}
    for row,rhs in equations:
        while row:
            pivot = row.bit_length()-1
            if pivot not in pivots:
                pivots[pivot] = (row,rhs)
                break
            existing,wanted = pivots[pivot]
            row ^= existing
            rhs ^= wanted
        else:
            if rhs:
                return None
    solution = 0
    for pivot,(row,rhs) in sorted(pivots.items()):
        if (row&solution).bit_count()%2 != rhs:
            solution |= 1 << pivot
    require(all((row&solution).bit_count()%2==rhs for row,rhs in equations),"sign solution failed replay")
    return {"rank":len(pivots),"sign_flip_mask_hex":hex(solution),"sign_flips":solution.bit_count()}


def _signature_checks(w,f,b,paths):
    frame,_ = lattice_frame()
    fock = LatticeFock([[dot(a,c) for c in frame] for a in frame])
    fw,bw = (tuple(integral(frame,root) for root in weights) for weights in (f,b))
    nf,nb = len(f),len(b)
    triples = [(int(k),int(col)//nf,int(col)%nf,int(w["J"][k,col])) for k,col in zip(*w["J"].nonzero())]
    pairs = [(int(label),*w["PAIRS"][column],int(w["W"][label,column]))
             for label,column in zip(*w["W"].nonzero())]

    def constraint(i,j,a,old,new):
        return (1 << i)^(1 << j)^(1 << (nf+a)),int(old!=new)

    def opposite_constraint(a,target):
        return ((1 << (nf+a))^(1 << (nf+w["BAR"][a])),int(w["ETA"][a]!=target))

    signs = [constraint(k,r,a,value,fock.sign(bw[a],fw[r])*fock.sign(fw[k],neg(fw[k])))
             for k,a,r,value in triples]
    signs += [constraint(i,j,a,value,fock.sign(fw[i],fw[j])) for a,i,j,value in pairs]
    signs += [opposite_constraint(a,fock.sign(bw[a],neg(bw[a]))) for a in range(nb)]
    common = _sign_solution(signs)
    require(common is not None and common["rank"]==116,"one shared root sign chart must match W and J")
    first,sign = signs[0]
    require(_sign_solution([(first,sign^1),*signs[1:]]) is None,"one changed cubic coefficient must fail the shared sign chart")

    # Lift the actual colour, weak-slot and last-slot Weyl transpositions
    # to signed one-particle permutations of the original W generators.
    generators = []
    for first,last in ((0,1),(1,2),(3,4),(6,7)):
        def permute(v):
            return tuple(v[last] if i==first else v[first] if i==last else value for i,value in enumerate(v))

        a_map,r_map = ([lookup[permute(weight)] for weight in lookup] for lookup in (b,f))
        require(all(a_map[w["BAR"][a]]==w["BAR"][a_map[a]] for a in range(nb)),"root partners must transport with the symmetry")
        equation = []
        for k,a,r,value in triples:
            new = int(w["J"][r_map[k],nf*a_map[a]+r_map[r]])
            require(new!=0,"a J word lost its support under a coordinate transposition")
            equation.append(constraint(k,r,a,value,new))
        for a,i,j,value in pairs:
            mapped = (r_map[i],r_map[j])
            new = int(w["W"][a_map[a],w["PAIR_INDEX"][tuple(sorted(mapped))]])
            require(new!=0,"a W pair lost its support under a coordinate transposition")
            equation.append(constraint(i,j,a,value,new if mapped[0]<mapped[1] else -new))
        equation.extend(opposite_constraint(a,w["ETA"][a_map[a]]) for a in range(nb))
        solution = _sign_solution(equation)
        require(solution is not None,"W/J/partner signs must admit a common signed symmetry")
        generators.append(((first+1,last+1),a_map,r_map,solution))
    all_paths = {tuple(row["indices_cap_boson_fermion_zero_based"]) for terms in paths.values() for row in terms}
    orbit = {min(all_paths)}
    for _ in range(len(all_paths)):
        orbit |= {(r_map[k],a_map[a],r_map[r]) for k,a,r in orbit for _,a_map,r_map,_ in generators}
        if len(orbit)==len(all_paths):
            break
    require(orbit==all_paths,"the signed W symmetries must connect all twelve weak components")
    return {"common_phase_chart":{**common,"constraints":len(signs),"free_signs":nf+nb-common["rank"],
                                  "single_coefficient_mutant_rejected":True},
            "time_and_pair_dressing_preserving_symmetries":[{"coordinate_transposition":list(coords),**solution}
                                                           for coords,_,_,solution in generators],
            "weak_components_related_by_actual_W_symmetries":len(orbit),
            "all_time_lower_to_upper_ratio":"1 wherever the dressed word is nonzero"}


def replay(root):
    import numpy as np

    refs = {str(path):digest(root/path) for path in INPUTS}
    require(all(refs[str(path)]==sha for path,sha in PINNED_INPUTS.items()),"review the changed original W source")
    w = _load_original(root)
    dressing = _original_pair_dressing(root,w)
    f, b = ({tuple(map(int,weight)):i for i,weight in enumerate(w[key])} for key in ("FW","BW"))
    j = w["J"]
    require(len(f)==64 and len(b)==60 and {charge_signature(r)[0] for r in f}=={3} and
            {charge_signature(r)[0] for r in b}=={2},"original W Cartan grading changed")
    count = 0
    for (p,f_index),(r,b_index) in product(f.items(),b.items()):
        total = add(r,p)
        actual = np.flatnonzero(j[:,len(f)*b_index+f_index]).tolist()
        expected = [f[neg(total)]] if total in ROOT_SET else []
        require(actual==expected,"the native W cubic and affine roots disagree on bracket support")
        count += len(actual)
    require(count==960,"original W charge action count changed")
    creators = _pair_creators(w)
    rows = {}
    time_reduction = None
    for path,words in triangle_components().items():
        if path not in ("weak_spinor_minus","weak_spinor_plus"):
            continue
        rows[path] = []
        for left,middle,right in words:
            k,a,r = f[left],b[middle],f[right]
            coefficient = int(j[k,len(f)*a+r])
            cap = {(1 << (column%len(f)),(column//len(f),)):int(j[k,column])
                   for column in map(int,np.flatnonzero(j[k]))}
            incoming = {(1 << r,()):1}
            require(abs(coefficient)==1 and _inner(cap,cap)==15 and
                    _inner(cap,_boson_creator(incoming,a))==coefficient,
                    "original charge-three W endpoint did not recover the bare word")
            require(not any(_lower_boson_raise_pair(cap,creators,w["pair_action"]).values()),
                    "the actual W cap must be in ker C3 adjoint")

            incoming = _neutral_pair(incoming,dressing)
            cap = _neutral_pair(cap,dressing)
            leaked_in = _lower_boson_raise_pair(incoming,creators,w["pair_action"])
            leaked_out = _lower_boson_raise_pair(cap,creators,w["pair_action"])
            a0 = _inner(cap,_boson_creator(incoming,a))
            mixed = _inner(leaked_out,_boson_creator(leaked_in,a))
            norms = [_inner(v,v) for v in (incoming,cap,leaked_in,leaked_out)]
            require((a0,mixed,norms)==(31*coefficient,458*coefficient,[30,465,465,6870]),
                    "the original W dressed short-time response changed")
            if time_reduction is None:
                time_reduction = {"indices_cap_boson_fermion_zero_based":[k,a,r],
                                  **_time_reduction(incoming,cap,a,w,creators,coefficient)}
            rows[path].append({"indices_cap_boson_fermion_zero_based":[k,a,r],"native_bracket_sign":coefficient,
                               "dressed_zero_time_word":a0,"mixed_transfer_leakage":mixed,
                               "mixed_jet_determinant_over_abs_g_squared":a0*mixed,
                               "endpoint_normalized_mixed_jet_determinant_over_abs_g_squared":
                                   str(sp.Rational(a0*mixed,norms[0]*norms[1]))})
    phase = _signature_checks(w,f,b,rows)
    return {"status":"EXACT_NATIVE_W_WEAK_WORDS_WITH_SCOPED_CLOCK_AND_PARITY_LIMITS",
            "source_refs":refs,"constructor_checks":len(w["checks"]),
            "cross_chart_bracket_support_pairs":len(f)*len(b),"nonzero_bracket_pairs":count,
            "bare_word":"native_bracket_sign * exp(-Delta * late_time) / sqrt(15)",
            "bare_affine_grade_one_time_match":False,
            "paired_boson_dressing_endpoint_charges":[5,7],"dressed_components":rows,
            "dressing":"original B_+ with opposite-root bosons; it carries zero Cartan weight and native number charge +4",
            "dressed_word_mixed_jet_scoped":"at zero interval times, all original W coupling terms included",
            "dressed_word_exact_time_reduction":time_reduction,
            "phase_and_symmetry_replay":phase,
            "all_time_normalized_early_first_probe":{"weak_spinor_minus":"-3*I/5","weak_spinor_plus":"3*I/5"},
            "prescribed_cartan_extension":{"commutes_with_H_by_original_vertices":480,
                "rule":"add the literal folded Q_d to H, assigning the original d charge to each boson and fermion",
                "neutral_pair_carries_zero_Q_d":True,
                "restores_folded_lower_to_upper_ratio":"sqrt(left_radial_factor * right_radial_factor)",
                "scope":"comparison with a prescribed clock extension; it retains finite rank in each component"},
            "original_w_fock_glue_rule":"grade modulo 4 is -N_f + 2*N_b, so its half character is exactly CAR parity",
            "grade_preserving_even_spinor_field_possible":False,
            "coordinate_empty_energy":"0 in Q=0",
            "known_Q2_energy_below_empty_at_nonzero_g":"(Delta - sqrt(Delta**2 + 32*abs(g)**2))/2",
            "all_fixed_finite_weighted_charge_sectors_have_finite_two_time_rank":True,
            "signed_affine_cauchy_ranks_checked":_affine_rank_check(),
            "scope":"fixed finite native W channels with the coordinate empty state and the given H; the sign comparison covers supported coefficient tensors, and physical R113 source/reference and arbitrary correlated weights are not inferred",
            "physical_source_identity_proved":False,"original_r113_values":None}
