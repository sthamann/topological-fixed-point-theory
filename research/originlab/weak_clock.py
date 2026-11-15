"""Resolve the weak doublet weights in the original folded E8 clock chart."""

from itertools import permutations, product

import sympy as sp

from affine_operators import (LatticeFock, MARKED_GAMMA, MARKED_P2, ROOTS, add, charge_signature, dot,
                              integral, lattice_frame, neg)
from affine_time import radial_word
from folded_time import marked_clock
from language import positive
from native_weak_w import DRESSED_WORD_TIME_RANK
from reconstruction import require, scalar
from weak_source import CHANNELS, COWEIGHT, WEAK_ROOT, WEAK_TORUS_ENTRY, axis, triangle_components


def _phase_sectors():
    z = WEAK_TORUS_ENTRY
    return {"spinor":-sp.I, **{
        f"{path}_{side}":sp.expand(phase*value)
        for path,phase in (("d8",1),("weak_spinor_minus",-sp.I),("weak_spinor_plus",sp.I))
        for side,value in (("upper",z),("lower",sp.conjugate(z)))
    }}


def _phase_matrix():
    sectors = _phase_sectors()
    require(len(set(sectors.values()))==7,"the joint action needs seven distinct eigenvalues")
    return sectors,sp.Matrix([[sp.expand(value**k) for value in sectors.values()] for k in range(8)])


def _root_transport(words):
    """Use two commuting compact D8 root actions to reach the minus path."""
    frame,_ = lattice_frame()
    fock = LatticeFock([[dot(a,b) for b in frame] for a in frame])
    original = {middle:right for _,middle,right in words["spinor"]}
    lifts,early_maps = [],set()
    for left,middle,right in words["weak_spinor_minus"]:
        colour = next(i for i in range(3) if middle[i])
        other_colours = [i for i in range(3) if i!=colour]
        removed_weak_slot = next(i for i in (3,4) if right[i]<0)
        a = add(axis(other_colours[0]),axis(other_colours[1]))
        b = add(axis(removed_weak_slot),axis(6))
        p = original[middle]
        old_late = add(middle,p)
        require(add(add(p,neg(a)),neg(b))==right and add(middle,right)==neg(left) and
                dot(a,b)==0 and all(dot(r,middle)==dot(r,COWEIGHT)==0 and
                                    dot(r,p)==dot(r,old_late)==1 for r in (a,b)),
                "the root rotations must preserve the middle field and signed winding")
        aa,bb,pp,rp,gg,pl = (integral(frame,v) for v in (a,b,p,right,middle,old_late))

        def transport(v):
            return fock.vertex(neg(aa),0,fock.vertex(neg(bb),0,v))

        def adjoint(v):
            return fock.vertex(bb,0,fock.vertex(aa,0,v))

        phase = fock.sign(neg(bb),pp)*fock.sign(neg(aa),integral(frame,add(p,neg(b))))
        for other in original.values():
            if other!=p:
                value = fock.vertex(integral(frame,other),-1,fock.vacuum())
                require(not transport(value),"the two root operations must select their own colour weight")
        for n,m in product(range(1,5),repeat=2):
            early = fock.vertex(pp,-n,fock.vacuum())
            transformed = transport(early)
            require(transformed==fock.combine((phase,fock.vertex(rp,-n,fock.vacuum()))) and
                    adjoint(transformed)==early,
                    "compact D8 root rotations did not preserve the marked descendant norm")
            late = fock.vertex(gg,n-m,early)
            transformed_late = transport(late)
            require(transformed_late==fock.vertex(gg,n-m,transformed) and
                    adjoint(transformed_late)==late,
                    "two-step transport must preserve the affine middle action and adjoint")
            outgoing = fock.vertex(pl,-m,fock.vacuum())
            response = fock.inner(outgoing,late)
            require(abs(response)==min(m,n) and fock.inner(transport(outgoing),transformed_late)==response,
                    "the full signed middle response must survive the root lift")
        lifts.append({"colour_slot":colour+1,"removed_weak_slot":removed_weak_slot+1,
                      "lowering_roots_doubled":[list(a),list(b)],"relative_phase":phase})
        early_maps.add((p,right))
        charge = neg(add(a,b))
        _,triality,hypercharge = charge_signature(charge)
        require(dot(charge,charge)==4 and triality==1 and hypercharge==sp.Rational(1,6) and
                abs(dot(charge,WEAK_ROOT))==1,"wrong charges for the composite source insertion")
    marked,weak = ({pair[k] for pair in early_maps} for k in (0,1))
    require(len(marked)==3 and len(weak)==2 and early_maps==set(product(marked,weak)),
            "the six lowering operations must span all maps between the early carriers")
    # The simple colour and weak root strings identify the two carriers.
    # Zero mode commutators on root currents give the same action at every grade.
    colour = [add(axis(i),neg(axis(i+1))) for i in (0,1)]
    require(all(sorted(dot(r,c) for r in marked)==[-1,0,1] for c in colour) and
            all(charge_signature(r)==(1,2,sp.Rational(1,3)) for r in marked) and
            all(charge_signature(r)==(3,0,sp.S.Half) for r in weak) and
            all(dot(r,c)==0 for r,c in product(weak,colour)) and
            all(dot(r,WEAK_ROOT)==0 for r in marked) and
            sorted(dot(r,WEAK_ROOT) for r in weak)==[-1,1],
            "mark and half spinor must carry the original colour dual triplet and weak doublet")
    return lifts,{"marked_carrier":"(bar(3),1)_(1/3)","weak_carrier":"(1,2)_(1/2)",
                  "early_ranks_per_grade":[len(marked),len(weak)],"map_space_rank":len(early_maps),
                  "insertion_multiplet":"(3,2)_(1/6)","cartan_weight_squared":4,
                  "embedding_requires_auxiliary_multiplet":"(bar(3),2)_(-1/6) for an equivariant embedding into the weak carrier tensor an auxiliary carrier",
                  "scope":"on the early affine cyclic carriers; a physical P1 operation is required to supply this action"}


def _time_stabilizer(background):
    # Fixing e8 restricts every E8 lattice isometry to the D7 roots in its
    # orthogonal hyperplane. Their signed permutations must preserve the
    # parity of the spinor coset. Fixing e1+e6 then leaves this D5 and A1.
    free = (1,2,3,4,6)
    stabilizer, orbit = [], set()
    for swap,permutation,signs in product((False,True),permutations(free),product((-1,1),repeat=5)):
        if signs.count(-1)%2:
            continue
        action = [(i,1) for i in range(8)]
        action[0],action[5] = ((5,1),(0,1)) if swap else ((0,1),(5,1))
        for i,j,sign in zip(free,permutation,signs):
            action[i] = (j,sign)
        value = tuple(sign*MARKED_P2[j] for j,sign in action)
        orbit.add(value)
        if all(background[i]==sign*background[j] for i,(j,sign) in enumerate(action)):
            stabilizer.append(value)
    return orbit,stabilizer


def _clock_exponents(background, words):
    def eta(root):
        return sum(d*r/2 for d,r in zip(background,root))

    rows = {}
    for path,tensors in words.items():
        pieces = [(dot(right,WEAK_ROOT),eta(right),eta(neg(left))) for left,_,right in tensors]
        if path == "spinor":
            require(len({left for _,_,left in pieces})==1,"marked late exponent must be fixed")
        else:
            upper = sorted((right,left) for weight,right,left in pieces if weight==1)
            lower = sorted((right,left) for weight,right,left in pieces if weight==-1)
            require([(right+sp.S.Half,left+sp.S.Half) for right,left in upper]==lower,
                    "the two weak weights must differ by one half in both intervals")
            require(len(set(right for right,_ in upper))==1,"early colour singlet must have one right exponent per weak weight")
        rows[path] = pieces
    return rows


def _known_weight_rank(matrix):
    # Retain the calibrated weak lower to upper coefficient s. On the real
    # folded time chart s=sqrt(x*y). The four effective columns are still
    # independent after replacing the old equal weights by (1,s).
    s = sp.Symbol("s",real=True)
    columns = [matrix[:4,0],*[matrix[:4,i]+s*matrix[:4,i+1] for i in (1,3,5)]]
    determinant = sp.expand(sp.Matrix.hstack(*columns).det())
    imaginary = sp.factor(sp.im(determinant))
    expected = sp.Rational(32,15625)*(s+1)*(481*s*s-616*s+5751)
    require(sp.expand(imaginary-expected)==0 and 616**2-4*481*5751<0,
            "calibrated four-sector matrix lost its rank certificate")
    return str(imaginary)


def _noncommuting_heat_check(cases):
    # Two original early heat factors at grade one, with a horizontal
    # SU(3) and SU(2) rotation on their respective carriers.
    bases = sp.Rational(2,3),sp.Rational(4,5)
    rotation = sp.Matrix([[3,-4],[4,3]])/5

    def check(key, link, upper):
        weights = [[base**(16*(1+sp.Rational(t))) for t in case[key]]
                   for case,base in zip(cases,bases)]
        a,b = (sp.diag(*diagonal) for diagonal in weights)
        bound = (max if upper else min)
        edge = (bound(weights[0])*bound(weights[1]))**2
        b = link*b*link.T
        mixed = b*a
        difference = edge*sp.eye(a.rows)-mixed.T*mixed
        require(a*b!=b*a and positive(difference if upper else -difference),
                "the folded norm bound failed its noncommuting transport control")
        return edge

    ratio = check("weak_spinor_early_tilts",rotation,True)/check("marked_early_tilts",sp.diag(rotation,1),False)
    require(ratio==sp.Rational(1024,18225),"unexpected same-history folded suppression")
    return str(ratio)


def controls(root):
    folded = marked_clock(root)
    words = triangle_components()
    sectors,matrix = _phase_matrix()
    lifts,colour_lift = _root_transport(words)
    marked_residues,weak_residues = set(),set()
    cases = []
    endpoints = {right for _,middle,right in words["weak_spinor_minus"] if middle==MARKED_GAMMA}
    e8_perp = [root for root in ROOTS if root[7]==0]
    require(len(e8_perp)==84 and all(all(value%2==0 for value in r) for r in e8_perp),
            "the fixed winding hyperplane must carry D7")
    x,y = sp.symbols("x y",positive=True)
    for case in folded["cases"]:
        quarter = sp.Rational(case["quarter"])
        background = tuple(map(sp.Rational,case["folded_cartan"]))
        exponents = _clock_exponents(background,words)
        original = sorted({right for _,right,_ in exponents["spinor"]})
        weak = sorted({right for _,right,_ in exponents["weak_spinor_minus"]})
        original_late = {left for _,_,left in exponents["spinor"]}
        weak_late = {left for _,_,left in exponents["weak_spinor_minus"]}
        require(exponents["weak_spinor_minus"]==exponents["weak_spinor_plus"],
                "swapping e7/e8 should preserve the fixed fold and reverse winding")
        marked_residues.update(value%1 for value in original)
        weak_residues.update(value%1 for value in weak)
        gap = weak[0]-original[-1]
        require(gap==(1-quarter)/4 and gap>0 and min(weak_late)-max(original_late)==gap,
                "the original and weak bands must separate in both intervals")
        orbit,stabilizer = _time_stabilizer(background)
        require(len(orbit)==16 and endpoints<=orbit and len(stabilizer)==4 and endpoints.isdisjoint(stabilizer),
                "unexpected same-middle, same-winding folded clock stabilizer")
        cases.append({"quarter":str(quarter),"marked_early_tilts":list(map(str,original)),
                      "weak_spinor_early_tilts":list(map(str,weak)),"early_gap_K":str(gap),
                      "marked_late_tilts":list(map(str,sorted(original_late))),
                      "weak_spinor_late_tilts":list(map(str,sorted(weak_late))),"late_gap_K":str(gap),
                      "early_gap_H8":str(8*gap),"fixed_winding_middle_orbit_size":len(orbit),
                      "clock_preserving_lattice_maps":len(stabilizer),
                      "clock_preserving_maps_to_weak_endpoints":0,
                      "fixed_colour_co_transport":[_co_transport(background,lift) for lift in lifts
                                                    if lift["colour_slot"]==1],
                      "normalized_affine_fold_at_x_y_quarter": {
                          path:str(sp.simplify(radial_word(sp.Rational(1,4),sp.Rational(1,4))*
                                               sum(sp.Rational(1,4)**(r+l) for _,r,l in rows)/len(rows)))
                          for path,rows in exponents.items()},
                      "affine_signed_fold_factor_times_radial_word":{
                          path:str(sum(x**l*y**r for _,r,l in rows)/len(rows))
                          for path,rows in exponents.items()},
                      "interval_tilts_by_component":{path:[[w,str(r),str(l)] for w,r,l in rows]
                                                     for path,rows in exponents.items()}})
    require(marked_residues.isdisjoint(weak_residues),
            "marked and weak spectra must remain disjoint after integer shifts on either fixed chart")
    target = matrix[:7,:]
    outside = sp.ones(8,1)   # Character 1 is outside the seven measured sectors.
    residual = sp.simplify(outside[7]-(matrix[7,:]*target.inv()*outside[:7,:])[0])
    require(residual==sp.Rational(144,125)*(1+sp.I),"seventh-power outside-sector control failed")
    return {"status":"EXACT_ORIGINAL_WEAK_FOLDED_CLOCK_COMPARISON","source_refs":folded["source_refs"],
            "cases":cases,"root_lifts":lifts,"root_lift_mode_pairs_count":16*len(lifts),
            "joint_colour_lift":colour_lift,
            "marked_early_residues":list(map(str,sorted(marked_residues))),
            "weak_early_residues":list(map(str,sorted(weak_residues))),
            "original_reflected_leg_norm_bound":"exp(-2 sum_j t_j (1-q_j)/4) times the marked norm squared on either interval",
            "norm_scope":"equal incident norms at the same grade, products of common folded heat steps and horizontal unitary gauge links in each leg irrep, before branch-dependent selection",
            "noncommuting_two_quarter_norm_ratio_bound":_noncommuting_heat_check(cases),
            "weak_lower_to_upper_literal_fold":"sqrt(left_radial_factor * right_radial_factor)",
            "fold_factor_rank_on_two_distinct_positive_widths_per_axis":{"marked":1,"weak":2},
            "fold_factor_rectangular_ratio_on_each_axis_(1/4,1/16)":{"marked":"1","weak":"85/81"},
            "weak_normalized_first_joint_probe_at_x_y_quarter":"12/25 - 3*I/5 for the minus-winding path",
            "phase_sectors":{name:str(value) for name,value in sectors.items()},
            "calibrated_four_column_imaginary_determinant_times_(1+s)^3":_known_weight_rank(matrix),
            "independent_power_seven_outside_character_residual":str(residual),
            "time_rank_test_widths":[str(sp.Rational(i,i+1)**16) for i in range(4,9)],
            "physical_source_identity_proved":False,"original_r113_values":None}


def _co_transport(background,lift):
    # The two reflections commute with the eighth-axis mirror. Transport
    # d on the first sheet and its mirror on the second sheet together.
    directions = [sp.Matrix(root)/2 for root in lift["lowering_roots_doubled"]]
    d = sp.Matrix(background)
    shift = -sum((d.dot(root)*root for root in directions),sp.zeros(8,1))
    mirror = sp.diag(1,1,1,1,1,1,1,-1)
    transport = (sp.eye(8)-directions[0]*directions[0].T)*(sp.eye(8)-directions[1]*directions[1].T)
    require(transport*d==d+shift and transport*mirror==mirror*transport and shift!=sp.zeros(8,1),
            "the reflected branch must co-transport its folded Cartan")
    return {"removed_weak_slot":lift["removed_weak_slot"],"folded_cartan_shift":list(map(str,shift)),
            "folded_cartan_shift_norm_squared":str(shift.dot(shift)),
            "ray_energy_increase_in_fixed_chart":str(-sum(d.dot(root) for root in directions)),
            "covariance_of_original_eighth_axis_mirror":True}


def checked_sectors(data, contract):
    """Recover resolved amplitudes from a fixed inserted action on supplied whole forms."""
    sectors,matrix = _phase_matrix()
    inverse = matrix[:7,:].inv()
    samples = data["samples"]
    require(isinstance(samples,list) and samples,"supply a whole time rectangle and a holdout")
    seen,training,holdouts = set(),set(),set()
    rows,amplitudes = [],{}
    for sample in samples:
        x,y = (scalar(sample[k]) for k in ("left_radial_factor","right_radial_factor"))
        holdout = sample["held_out"]
        require(type(holdout) is bool and x.is_Rational and y.is_Rational and 0<x<1 and 0<y<1 and (x,y) not in seen,
                "need distinct rational time pairs and a Boolean holdout label")
        seen.add((x,y))
        (holdouts if holdout else training).add((x,y))
        forms = sample["power_forms"]
        require(isinstance(forms,list) and len(forms)==8,"supply whole absolute word forms for powers zero to seven")
        values = sp.Matrix([contract(form) for form in forms])
        weights = (inverse*values[:7,:]).applyfunc(sp.simplify)
        recovered = dict(zip(sectors,weights))
        amplitudes[x,y] = recovered
        expected = sp.sqrt(x*y)
        balance = {}
        for name in CHANNELS:
            if name=="spinor":
                continue
            upper,lower = (recovered[f"{name}_{side}"] for side in ("upper","lower"))
            informative = upper!=0 or lower!=0
            balance[name] = {"measured_ratio":str(sp.simplify(lower/upper)) if upper else None,
                             "matches_literal_fold":sp.simplify(lower-expected*upper)==0 if informative else None,
                             "matches_equal_weights":lower==upper if informative else None}
        prediction = sp.simplify((matrix[7,:]*weights)[0])
        rows.append({"radial_factors":[str(x),str(y)],"held_out":holdout,
                     "absolute_sector_amplitudes":{k:str(v) for k,v in recovered.items()},
                     "absolute_uninserted_word":str(values[0]),
                     "power_seven_absolute":str(values[7]),"power_seven_prediction":str(prediction),
                     "power_seven_residual":str(sp.simplify(values[7]-prediction)),
                     "literal_fold_lower_to_upper":str(expected),"weak_weight_comparison":balance})
    left,right = (sorted({pair[k] for pair in training}) for k in (0,1))
    require(len(left)>=2 and len(right)>=2 and training==set(product(left,right)) and
            any(x not in left and y not in right for x,y in holdouts),
            "supply a calibration rectangle and a holdout with both widths unseen")
    if any(row["power_seven_residual"]!="0" for row in rows):
        status = "SUPPLIED_SEVEN_SECTOR_MODEL_REFUTED"
    elif all(all(value=="0" for value in row["absolute_sector_amplitudes"].values()) for row in rows):
        status = "INCONCLUSIVE_ZERO_RESPONSE"
    else:
        status = "SUPPLIED_SEVEN_SECTOR_MODEL_REPLAYED"
    rank_report = {"status":"UNVERIFIED_SECTOR_AMPLITUDES"}
    if status!="SUPPLIED_SEVEN_SECTOR_MODEL_REFUTED":
        ranks = {name:sp.Matrix([[amplitudes[x,y][name] for y in right] for x in left]).to_DM().rank()
                 for name in sectors}
        separating = [name for name,rank in ranks.items()
                      if name.startswith("weak_spinor_") and rank>DRESSED_WORD_TIME_RANK]
        rank_report = {
            "status": "FIXED_NATIVE_W_DRESSING_EXCLUDED" if separating else "RANK_DOES_NOT_SEPARATE",
            "calibration_shape":[len(left),len(right)],"complex_sector_ranks":ranks,
            "fixed_native_w_rank_bound":DRESSED_WORD_TIME_RANK,"separating_sectors":separating,
            "scope":"one fixed B_+ dressed W word with its original Hamiltonian, allowing separable clock factors",
        }
    return {"status":status,"samples":rows,"resolved_time_rank":rank_report,
            "physical_source_identity_proved":False,
            "scope":"finite powers of the same early insertion on the supplied whole source; clock ratios require their actual time map"}
