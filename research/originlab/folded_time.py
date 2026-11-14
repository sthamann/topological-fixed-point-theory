"""Compare the marked affine word in the original radial and quarter times."""

from itertools import permutations, product
from pathlib import Path
import json

import sympy as sp

from affine_time import marked_charges, radial_grid, radial_word
from reconstruction import require
from source_return import REGIONAL
from verify import digest


FOLD = Path("experiments/theory-contracts/v16-operator-pullback-20261014/fold-results.json")
TIME = Path("experiments/theory-contracts/joint-source-clock-field-closure-20261102/"
            "SOURCE-TIME-AND-AFFINE-CLOCK-TEST.txt")
WINDING = Path("docs/research/fundamental-plan/TFPT-E8-P1-LADUNGSWINDUNG-20261112.txt")
MIRROR = Path("docs/research/fundamental-plan/TFPT-HALBTAKT-UND-QUARTERZEIT-20261111.txt")
FULL_MIRROR = Path("docs/research/fundamental-plan/TFPT-GELADENE-NAHT-VERKLEBUNG-20261111.txt")
STABILIZER = Path("docs/research/fundamental-plan/TFPT-MARKIERTER-STABILISATOR-20261112.txt")
SOURCES = {
    FOLD: "253f370a9c84b70c28169ed0aad83d1866245804227e44aa0189e823382f5639",
    TIME: "1707f253a9f8794f63fede4444e2ea0ff673a84530eccabf0b4f9baf7b0624db",
    WINDING: "2a68509aa343f1705d21c4839c22f1ebecb858843f9080aec6b14429393bd473",
    MIRROR: "ddca975ad7d873330d6682c1907a8617536b2116161df4ec5bdbdcbaa2518c90",
    FULL_MIRROR: "0fbc0d768bcf78be1f25e77f78782edc351c1e24822ee8f85fe9255b4d6aae67",
    STABILIZER: "5db3a43389d9da66a8490a10aa8598b461d4ec48350ff3ab1435cf179f896cc3",
    REGIONAL: "6888980fcb80f292c4da087a5bc3f486395097be914794cb22b6ce876f40c64d",
}
HYPERCHARGE = sp.Matrix([-2,-2,-2,3,3,0,0,0])/6
QUARTERS = (sp.Rational(1,4),sp.Rational(3,4))


def folded_grid(left, right, eta_left, eta_right):
    """The original two transfers on the given E8 lattice comparison."""
    return sp.diag(*(x**eta_left for x in left))*radial_grid(left,right)*sp.diag(*(y**eta_right for y in right))


def mixed_ratio(matrix):
    if not (matrix[0,1] and matrix[1,0]):
        return None
    return sp.cancel(matrix[0,0]*matrix[1,1]/(matrix[0,1]*matrix[1,0]),extension=sp.I)


def visible_defect(matrix):
    """First-row and first-column calibration predicts all held-out entries."""
    return [(m,n) for m in range(1,matrix.rows) for n in range(1,matrix.cols)
            if sp.cancel(matrix[m,n]*matrix[0,0]-matrix[m,0]*matrix[0,n],extension=sp.I) != 0]


def positive_gauge_average(background):
    # All gauge rotations preserve the projection onto Y. The Weyl average
    # removes precisely the noncentral components of this Cartan element.
    frames = []
    for color,weak in product(permutations(range(3)),permutations((3,4))):
        frames.append(sp.Matrix([*(background[i] for i in (*color,*weak)),*background[5:]]))
    average = sum(frames,sp.zeros(8,1))/len(frames)
    overlap = background.dot(HYPERCHARGE)
    require(overlap<0 and all(v.dot(HYPERCHARGE)==overlap for v in frames),
            "fixed gauge orbit must preserve the nonzero hypercharge component")
    require(average[0]==average[1]==average[2] and average[3]==average[4],
            "Weyl average failed to remove the noncentral gauge components")
    return {"hypercharge_inner_product":str(overlap),"coefficient_of_Y":str(overlap/HYPERCHARGE.dot(HYPERCHARGE)),
            "fixed_cartan_weyl_average":list(map(str,average)),
            "minimum_norm_squared_of_fixed_quarter_gauge_averages":str(average.dot(average))}


def coset_bottom(background,shift,parity):
    # With |d_i|<=1/2, moving any coordinate back by two toward zero
    # preserves the coset parity and strictly lowers energy outside these sets.
    require(max(map(abs,background))<=sp.S.Half,"mirror minimum needs bounded original Cartan shifts")
    values = (-2,0,2) if shift==0 else (-3,-1,1,3)
    weights = [int(16*d) for d in background]
    require(list(background)==[sp.Rational(n,16) for n in weights],"unexpected mirror energy denominator")
    positions = [v for v in product(values,repeat=8) if (sum(v)-8*shift)%4==2*parity]
    scores = [sum(4*v*v+d*v for v,d in zip(row,weights)) for row in positions]
    minimum = min(scores)
    return {"K_energy":str(sp.Rational(minimum,32)),"ground_momenta_doubled":
            [list(row) for row,score in zip(positions,scores) if score==minimum]}


def mirrored_isometry(background,first,second,original_backgrounds):
    # The actual eighth-axis reflection maps the whole even lattice to the
    # other even lattice. It agrees with p -> p-e8 only on the pinned legs.
    e8 = sp.eye(8)[:,7]
    mirror = sp.eye(8)-2*e8*e8.T
    reflected = mirror*background
    shift = reflected-background
    require(mirror.T*mirror==sp.eye(8) and mirror*HYPERCHARGE==HYPERCHARGE and shift==e8/2,
            "wrong same-clock reflected Cartan")
    source_legs = (first,second)
    target_legs = tuple(mirror*p for p in source_legs)
    require(not any(all(p.dot(reflected-d).is_Integer for p in target_legs) for d in original_backgrounds),
            "reflected Cartan is equal to an old quarter modulus in the fixed charge chart")
    fixed = [8*(p.dot(p)/2+background.dot(p)) for p in target_legs]
    transported = [8*(p.dot(p)/2+reflected.dot(p)) for p in target_legs]
    require([x-y for x,y in zip(fixed,transported)]==[2,2] and
            transported==[8*(p.dot(p)/2+background.dot(p)) for p in source_legs],
            "full mirror relative energies changed")
    # Compare K on the source after pulling back the unchanged mirror
    # background. On the fixed marked spinor lines its extra charge is 1/4.
    phase_on_spinor = sp.simplify(sp.exp(2*sp.pi*sp.I*shift.dot(first)))
    require(phase_on_spinor==sp.I and shift.dot(2*e8)==1,"wrong pulled-back quarter character")
    x,y = sp.symbols("x y",positive=True)
    fixed_factor = x**shift.dot(first)*y**shift.dot(second)
    require(fixed_factor==(x*y)**sp.Rational(1,4),"wrong same-background full mirror time shift")
    return {"target_empty_impulse":"0","target_lattice":"r(E8) = D8 union the opposite spinor coset",
            "fixed_background_marked_energies_in_H8_time":list(map(str,fixed)),
            "transported_background_marked_energies_in_H8_time":list(map(str,transported)),
            "required_increment_to_cartan":list(map(str,shift)),"fixed_Y_preserved":True,
            "restores_all_relative_energies_with_original_lattice_isometry":True,
            "marked_vacuum_preserved_by_transport":True,
            "untouched_background_relative_to_transported_word_factor":"(x*y)**(1/4)",
            "relative_character_of_pulled_back_fixed_generator_at_K_angle_2pi":str(phase_on_spinor),
            "physical_reflection_identity_proved":False}


def marked_clock(root):
    refs = {str(path):digest(root/path) for path in SOURCES}
    require(all(refs[str(path)]==sha for path,sha in SOURCES.items()),"review changed original fold, time or winding")
    frame = sp.Matrix(json.loads((root/FOLD).read_bytes())["explicit_Weyl_matrix"])
    require(frame.T*frame==sp.eye(8),"the original fold frame is no longer orthogonal")
    left,middle,right = (sp.Matrix([sp.Rational(v) for v in row["root"]]) for row in marked_charges())
    p1,p2 = -left,right
    require(p1-p2==middle,"wrong serial transfer charges")
    require(tuple(map(tuple,(frame.T*p1,frame.T*p2)))==((1,1,0,0,0,0,0,0),(0,1,-1,0,0,0,0,0)),
            "wrong two original folded end lines")
    e = sp.eye(8)
    twist = -e[:,7]
    p1_reflected = p1+twist
    p2_reflected = p2+twist
    require(p1_reflected.dot(p1_reflected)==p2_reflected.dot(p2_reflected)==2 and
            p1_reflected-p2_reflected==middle and twist.dot(HYPERCHARGE)==0,
            "wrong marked reflected roots or neutral reference shift")
    x,y = sp.symbols("x y",positive=True)
    mixed_derivative = sp.cancel(x*y*sp.diff(sp.diff(radial_word(x,y),x)/radial_word(x,y),y))
    require(sp.cancel(mixed_derivative+x*y/(1+x*y)**2)==0,"marked affine mixed derivative changed")
    probe = (sp.S.Half,sp.Rational(1,4))
    affine_ratio = mixed_ratio(radial_grid(probe,probe))
    require(affine_ratio==sp.Rational(81,85),"wrong normalized affine rectangular response")
    cases = []
    original_backgrounds = [frame*sp.Matrix([(k+quarter)/8-sp.S.Half for k in range(8)]) for quarter in QUARTERS]
    require(all(d.dot(d)==sp.Rational(85,128) for d in original_backgrounds),
            "quarter background norm changed")
    for quarter,background in zip(QUARTERS,original_backgrounds):
        eta_left,eta_right = background.dot(p1),background.dot(p2)
        gauge_shifts = [background[0]-background[1],background[1]-background[2],background[3]-background[4]]
        # Pin one original nonzero Cartan component by a gauge-invariant scalar.
        charge_overlap = background.dot(HYPERCHARGE)
        require(eta_right==-sp.Rational(1,8) and eta_left==quarter/4-sp.Rational(7,8) and
                charge_overlap==-sp.Rational(1,12)-quarter/12,"wrong original folded time marks")
        # The original fold has the empty impulse and one charged zero line.
        # Check both whole cosets with the same exact bounded minimization.
        zero_line = e[:,0]+e[:,3]
        ground = [coset_bottom(background,shift,0) for shift in (0,1)]
        require([s["K_energy"] for s in ground]==["0","1/8"] and
                ground[0]["ground_momenta_doubled"]==[[0]*8,list(map(int,2*zero_line))],
                "unexpected original folded vacuum and charged zero line")
        # Keep the actual mirror shift and test the whole translated module
        # in the same fixed background. In the original D8 chart the two
        # vector ground momenta are e1 and e4, with energy zero.
        reflection_gaps = [(p.dot(p)-twist.dot(twist))/2+background.dot(p-twist)
                           for p in (e[:,0],p1_reflected,p2_reflected)]
        require(reflection_gaps==[-sp.Rational(3,4),quarter/4-sp.Rational(3,8),sp.Rational(3,8)],
                "wrong fixed-background mirror ground or marked energies")
        mirror_bottom = {"vector":coset_bottom(background,0,1),"opposite_spinor":coset_bottom(background,1,1)}
        require(mirror_bottom["vector"]["K_energy"]=="0" and mirror_bottom["opposite_spinor"]["K_energy"]=="1/8",
                "unexpected bottom of the two original mirror cosets")
        weak = (e[:,3]-e[:,4])/2
        require([v.dot(HYPERCHARGE) for v in (twist,p1_reflected)]==[0,0] and
                [v.dot(weak) for v in (twist,p1_reflected)]==[0,0] and
                all(len({v[j] for j in range(3)})==1 for v in (twist,p1_reflected)),
                "the mirror witness must have the same Y, color and weak Cartan charges")
        # The reflection comparison keeps the two marked charge rays. The
        # equality of gaps is necessary for an antiunitary fixing time units.
        original_gaps = [1+eta_left,1+eta_right]
        require([a-b for a,b in zip(original_gaps,reflection_gaps[1:])]==[sp.S.Half,sp.S.Half],
                "unexpected locked marked-mirror energy difference")
        variables = sp.Matrix(sp.symbols("p0:8",real=True))
        paired_backgrounds = []
        for sign in (1,-1):
            # Translation to the twisted coset, keeping the clock on every
            # momentum, fixes the other background by the linear coefficient.
            target_background = sign*background-twist
            moved = sign*variables+twist
            difference = (moved.dot(moved)-twist.dot(twist)-variables.dot(variables))/2
            difference += target_background.dot(moved-twist)-background.dot(variables)
            norm = target_background.dot(target_background)
            require(sp.expand(difference)==0 and norm!=background.dot(background),
                    "wrong required reflected co-transport of the original clock")
            paired_backgrounds.append({"charge_orientation":sign,"cartan":list(map(str,target_background)),
                                       "norm_squared":str(norm)})
        # A postulated time dilation calibrates the leading slopes only if
        # independent leg scales are allowed. The rectangular invariant tests both.
        rate_left,rate_right = (16*(1+eta) for eta in (eta_left,eta_right))
        width_scales = ([r**rate_left for r in probe],[r**rate_right for r in probe])
        calibrated = mixed_ratio(radial_grid(*width_scales))
        require(calibrated!=affine_ratio,"independent leading-rate calibration unexpectedly matched the full grid")
        # The exponents have denominator 16. Use exact rational radial factors.
        exact_probe = tuple(r**16 for r in probe)
        require(mixed_ratio(folded_grid(exact_probe,exact_probe,eta_left,eta_right))==
                mixed_ratio(radial_grid(exact_probe,exact_probe)),"time twist is not a separate end factor")
        cases.append({"quarter":str(quarter),"folded_cartan":list(map(str,background)),
                      "end_exponents_in_K_time":list(map(str,(eta_left,eta_right))),
                      "marked_energies_in_H8_time":list(map(str,(8*(1+eta_left),8*(1+eta_right)))),
                      "left_to_right_primary_rate_ratio":str((1+eta_left)/(1+eta_right)),
                      "untransported_sm_simple_root_shifts":list(map(str,gauge_shifts)),
                      "inner_product_with_original_hypercharge":str(charge_overlap),
                      "positive_gauge_average":positive_gauge_average(background),
                      "original_full_fold_ground_momenta_doubled":ground[0]["ground_momenta_doubled"],
                      "twisted_mirror_gaps_relative_to_minus_e8_reference_in_H8_time":
                          dict(zip(("e1_control","reflected_p1","reflected_p2"),map(str,(8*g for g in reflection_gaps)))),
                      "fixed_background_translated_lattice_bottom_by_coset":mirror_bottom,
                      "clock_preserving_mirror_backgrounds":paired_backgrounds,
                      "whole_lattice_mirror":mirrored_isometry(background,p1,p2,original_backgrounds),
                      "mixed_ratio_with_independent_leading_rate_scales":str(calibrated)})
    # The common contribution is invariant under each original gauge conjugation.
    # This consequence is linear and extends to noncommuting positive segments
    # by taking their common-support infinitesimal tangent.
    return {"status":"EXACT_ORIGINAL_FOLDED_RADIAL_COMPARISON","source_refs":refs,"cases":cases,
            "clock_choice_on_fixed_word":{"one_common_time_reparametrization":False,
                "independent_radial_leg_rescaling_reproduces_whole_word":False,
                "linear_cartan_twist_factor":"x**eta_left * y**eta_right",
                "normalization_invariant_mixed_response":"-x*y/(1+x*y)**2",
                "symmetric_exact_reduced_radial_probe_ratio":str(affine_ratio)},
            "gauge_frame_repair":{"all_positive_original_quarter_gauge_orbit_sums_have_nonzero_hypercharge_projection":True,
                "scope":"positive short-time segments without extra generators, fixed physical hypercharge and unchanged radial part",
                "opposite_background_repair_if_both_are_originally_admitted":"exp(-t*(L0+D)/2) exp(-t*(L0-D)/2) = exp(-t*L0)"},
            "mirror_fixed_background":{"scope":"same folded clock on v+E8, v=-e8 as the selected stationary reference",
                "reference_is_ground_state":False,"original_selector_or_reflected_time_identified":False},
            "physical_r113_endpoint_values":None}


def observed_two_time_grid(observed,radial,coordinates,times,root):
    report = {"mixed_time_ratios": {"measured":None,"radial":str(mixed_ratio(radial))},
              "separable_deviation_mismatch_indices":None,"folded_clock_comparison":None}
    ratio = mixed_ratio(observed)
    if ratio is not None:
        report["mixed_time_ratios"]["measured"] = str(ratio)
    # Normalize by the nonzero analytic affine grid. A zero matrix or a
    # missing pivot cannot support a relative calibration.
    factors = sp.Matrix(observed.rows,observed.cols,
                        lambda m,n:sp.cancel(observed[m,n]/radial[m,n],extension=sp.I))
    if factors[0,0]:
        report["separable_deviation_mismatch_indices"] = visible_defect(factors)
    if root is None:
        report["folded_clock_comparison"] = {"status":"ORIGINAL_FOLD_SOURCES_NOT_AVAILABLE"}
        return report
    # A common exponential time map is needed to interpret the fold rate.
    # x_i**t_j = x_j**t_i on both axes, with one shared rate and fixed origin.
    indexed = [(t,x) for axis,scales in zip(times,coordinates) for t,x in zip(axis,scales)]
    t0,x0 = indexed[0]
    def same_rate(claim):
        t,x = claim
        r = t/t0
        if r.is_Integer:
            return x==x0**r
        # A reduced rational can be a rational q-th power only when
        # both numerator and denominator have exact q-th roots.
        top,top_exact = sp.integer_nthroot(x0.p,r.q)
        bottom,bottom_exact = sp.integer_nthroot(x0.q,r.q)
        return top_exact and bottom_exact and x==sp.Rational(top,bottom)**r.p
    if not all(same_rate(entry) for entry in indexed[1:]):
        report["folded_clock_comparison"] = {"status":"SUPPLIED_RADIAL_TIME_MAP_NOT_A_COMMON_EXPONENTIAL"}
        return report
    folded = marked_clock(root)
    comparisons,targets = [],[]
    for row in folded["cases"]:
        etas = tuple(sp.Rational(x) for x in row["end_exponents_in_K_time"])
        target = folded_grid(*coordinates,*etas)
        targets.append(target)
        difference = [(m,n) for m in range(target.rows) for n in range(target.cols)
                      if sp.simplify(observed[m,n]-target[m,n])!=0]
        relative_mismatches = None
        if observed[0,0]:
            relative_mismatches = [(m,n) for m in range(target.rows) for n in range(target.cols)
                                   if sp.simplify(observed[m,n]*target[0,0]-target[m,n]*observed[0,0])!=0]
        comparisons.append({"quarter":row["quarter"],"matches_up_to_one_global_normalization":
                            None if relative_mismatches is None else not relative_mismatches,
                            "relative_mismatch_indices":relative_mismatches,
                            "absolute_matches_supplied_grid":not difference,
                            "absolute_mismatch_indices":difference})
    calibration = sp.Matrix([[target[i,0] for target in targets] for i in (0,1)])
    coefficients = calibration.inv()*observed.extract((0,1),(0,))
    coefficients = coefficients.applyfunc(sp.simplify)
    prediction = coefficients[0]*targets[0]+coefficients[1]*targets[1]
    mixture_mismatches = [(m,n) for m in range(observed.rows) for n in range(observed.cols)
                          if (m,n) not in ((0,0),(1,0)) and sp.simplify(observed[m,n]-prediction[m,n])!=0]
    mixture = {"calibration_indices":[[0,0],[1,0]],"coefficients":list(map(str,coefficients)),
               "positive_normalized_mixture":not mixture_mismatches and
                    all(c.is_nonnegative is True for c in coefficients) and sum(coefficients)==1,
               "held_out_mismatch_indices":mixture_mismatches,"additional_left_time_held_out":observed.rows>2}
    report["folded_clock_comparison"] = {"status":"EXACT_COMPARISON_ON_SUPPLIED_GRID","quarters":comparisons,
                                        "stationary_mixture_of_original_quarters":mixture}
    return report
