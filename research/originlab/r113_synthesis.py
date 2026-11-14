"""Exact R113 synthesis for bounded history and charged spectral contracts."""
from itertools import combinations, product
from math import comb
from pathlib import Path
import subprocess
import sys

import sympy as sp

from reconstruction import matrix, require, scalar
from return_moments import reconstruct_return, rows
from source import MARKS, NATIVE
from verify import digest

MARKED_SOURCE = Path("docs/research/fundamental-plan/p1-original-marken-und-metrik-20261112/")
NATIVE_GRADE_CHECK = MARKED_SOURCE / "check_tfpt_native_middle_grade.py"
NATIVE_GRADE_SHA = "7c69f98aaa4ad5fb2a25d08ad9344b232933024f85d408421e6ecffc8c6e485c"
NATIVE_DOMAIN_SHA = "7fbe04de691cfedbde4eec25735727e86daf06aaedba9ecd8b314a4e0cc2d9bf"


def history_range(data):
    """Vertices of a rational probability polytope on declared histories."""
    histories = data["histories"]
    require(isinstance(histories,list) and histories and all(isinstance(h,str) and h for h in histories)
            and len(histories)==len(set(histories)),"use distinct named histories")
    n = len(histories)
    constraints = data["constraints"]
    require(isinstance(constraints,dict) and "Z0" in constraints,"empty response constraint is required")
    reference = scalar(data["reference"])
    require(reference.is_positive,"the same positive absolute reference is required")
    coefficients = matrix([[1]*n,*[item["coefficients"] for item in constraints.values()]])
    target = sp.Matrix([1,*[reference*scalar(item["value"]) for item in constraints.values()]])
    mixed = matrix([data["mixed_coefficients"]])
    require(coefficients.cols==mixed.cols==n and all(v.is_Rational for v in (*coefficients,*target,*mixed)),
            "use real rational constraints and one coefficient per history")
    reduced, pivots = coefficients.row_join(target).rref()
    if pivots and pivots[-1]==n:
        return {"status":"EXCLUDED_IN_DECLARED_HISTORY_CLASS","reason":"inconsistent linear constraints"}
    rank = len(pivots)
    if comb(n,rank)>4096:
        return {"status":"BOUND_EXCEEDED","basis_count":comb(n,rank)}
    rhs = reduced[:rank,n]
    extrema = None
    count = 0
    for support in combinations(range(n),rank):
        basis = reduced[:rank,list(support)]
        if basis.det()==0:
            continue
        solution = basis.inv()*rhs
        if any(v<0 for v in solution):
            continue
        weights = sp.zeros(n,1)
        for i,value in zip(support,solution):
            weights[i] = value
        count += 1
        value = (mixed*weights)[0]/reference
        if extrema is None:
            extrema = [(value,weights),(value,weights)]
        if value<extrema[0][0]:
            extrema[0] = (value,weights)
        if value>extrema[1][0]:
            extrema[1] = (value,weights)
    if extrema is None:
        return {"status":"EXCLUDED_IN_DECLARED_HISTORY_CLASS","reason":"no nonnegative normalized history weights"}
    low,high = (item[0] for item in extrema)
    return {"status":"MIXED_RESPONSE_FIXED" if low==high else "MIXED_RESPONSE_UNDERDETERMINED",
            "history_labels":histories,"retained_constraints":list(constraints),
            "feasible_bases":count,"mixed_response_interval":[str(low),str(high)],
            "witnesses":[{"value":str(v),"weights":list(map(str,w))} for v,w in extrema],
            "scope":"all nonnegative probability weights on these fixed histories and over the declared reference"}


def spectral_pair(data):
    """Exhaust real insertion signs on a fixed simple charged transfer spectrum."""
    levels = tuple(map(scalar,data["spectrum"]))
    n = len(levels)
    require(1<=n<=10 and len(set(levels))==n and all(v.is_Rational and 0<v<=1 for v in levels),
            "one to ten distinct rational transfer levels in (0,1] required")
    diagonals = matrix(data["diagonal_observations"])
    require(diagonals.shape==(n,2) and all(v.is_Rational for v in diagonals),
            "supply n exact real observations for each of the two individual marks")
    vandermonde = sp.Matrix([[value**k for value in levels] for k in range(n)])
    weights = vandermonde.inv()*diagonals
    if any(v<0 for v in weights):
        return {"status":"EXCLUDED_IN_DECLARED_SPECTRAL_CLASS","reason":"negative single-mark spectral weight"}
    observed = {int(k):scalar(v) for k,v in data["mixed_observations"].items()}
    require(len(observed)==len(data["mixed_observations"]) and
            all(0<=k<=64 and v.is_Rational for k,v in observed.items()),"invalid known mixed observation")
    amplitudes = [sp.sqrt(weights[i,0]*weights[i,1]) for i in range(n)]
    transfer = sp.diag(*levels)
    first = sp.Matrix([sp.sqrt(v) for v in weights[:,0]])
    candidates = []
    for signs in product(*[(-1,1) if a else (1,) for a in amplitudes]):
        overlaps = [sign*a for sign,a in zip(signs,amplitudes)]
        if all(sp.simplify(sum(a*v**k for a,v in zip(overlaps,levels))-target)==0
               for k,target in observed.items()):
            candidates.append((overlaps,sp.Matrix([sign*sp.sqrt(v) for sign,v in zip(signs,weights[:,1])])))
    if not candidates:
        return {"status":"EXCLUDED_IN_DECLARED_SPECTRAL_CLASS","reason":"no real phase assignment matches"}
    report = {"status":"MIXED_RESPONSE_FIXED","compatible_mixed_spectra":len(candidates),
              "spectrum":list(map(str,levels)),"single_mark_spectral_weights":rows(weights),
              "known_mixed_steps":sorted(observed),"absolute_reference":"1 (declared control)",
              "scope":"fixed simple spectrum, real marks, common neutral vacuum and charged intermediate sector"}
    if len(candidates)==1:
        return {**report,"mixed_spectral_weights":list(map(str,candidates[0][0]))}
    # Distinct weights on n fixed levels cannot agree on n consecutive moments.
    separation = next((k for k in range(n) if len({
        sp.simplify(sum(a*v**k for a,v in zip(overlaps,levels))) for overlaps,_ in candidates})>1),None)
    require(separation is not None,"distinct mixed spectra must separate by the Vandermonde bound")
    witnesses = []
    # A single cyclic insertion can need n powers before the next level
    # first shows a stable Hankel rank.
    steps = 2*(n+1)
    for overlaps,second in candidates:
        value = sp.simplify(sum(a*v**separation for a,v in zip(overlaps,levels)))
        if witnesses and value==witnesses[0][0]:
            continue
        marks = first.row_join(second)
        moments = [(marks.H*transfer**k*marks).applyfunc(sp.simplify) for k in range(steps+2)]
        replay = reconstruct_return(moments[:steps],moments[steps:])
        require(replay["status"]=="EXACT_ACCESSIBLE_POSITIVE_TRANSFER","spectral witness did not replay")
        witnesses.append((value,{"mixed_value":str(value),"mark_vectors":rows(marks),
                                 "moment_zero":rows(moments[0]),"positive_transfer_replay":replay}))
        if len(witnesses)==2:
            break
    return {**report,"status":"MIXED_RESPONSE_UNDERDETERMINED",
            "next_observation":{"word":["Y_adjoint",f"T^{separation}","U"],"step":separation,
                                "word_convention":"U and Y prepare same-charge vectors; Y_adjoint returns to neutral",
                                "witness_values":[str(v) for v,_ in witnesses]},
            "witnesses":[w for _,w in witnesses]}


def control():
    return spectral_pair({"spectrum":["1/2","1/4","1/8"],
                          "diagonal_observations":[[3,6],["7/8","19/8"],["21/64","69/64"]],
                          "mixed_observations":{"0":2}})


def from_input(data):
    from complex_moments import analyze
    from source_jets import gaussian_jet

    require(data.get("input_role")=="declared_finite_oracle","declare the synthesis input as a finite oracle")
    routes = {"fixed_histories":history_range,"real_charged_spectrum":spectral_pair,
              "complex_charged_spectrum":analyze,"gaussian_source_jet":gaussian_jet}
    kind = data.get("kind")
    require(kind in routes,"unknown R113 synthesis class")
    return {"input_role":data["input_role"],**routes[kind](data)}


def native_synthesis(root, domain, endpoint):
    """Eliminate native grade-zero words before attempting endpoint integration."""
    from affine_time import marked_charges
    from original_paths import original_word
    from root_doublet import radial_symmetry_check
    from complex_moments import control as complex_control

    if domain.resolve()!=(root/NATIVE).resolve():
        return {"status":"CUSTOM_DOMAIN_REQUIRES_OWN_SOURCE_RULES","control":control()}
    word = original_word(root)
    refs = word["source_refs"]
    domain_sha = digest(domain)
    grade_sha = digest(root/NATIVE_GRADE_CHECK)
    require((domain_sha,grade_sha)==(NATIVE_DOMAIN_SHA,NATIVE_GRADE_SHA),
            "review changed native domain or center-grading checker")
    checked = subprocess.run([sys.executable,"-B",str(root/NATIVE_GRADE_CHECK),str(root)],
                             cwd=root,capture_output=True,text=True,timeout=30)
    require(checked.returncode==0 and domain_sha in checked.stdout,
            f"native center-grading replay failed: {checked.stderr.strip()}")
    leg_data = marked_charges()
    grades = [leg["color_triality"] for leg in reversed(leg_data)]
    path = [0]
    for grade in grades:
        path.append((path[-1]+grade)%3)
    require(path[-1]==0 and path[1]!=0 and path[1]!=path[2],"changed original marked colour path")
    return {"status":"ORIGINAL_CHARGED_PRODUCER_REQUIRED" if endpoint["status"]=="BLOCKED" else
            "ORIGINAL_ENDPOINT_REQUIRES_SOURCE_REVIEW",
            "domain_sha256":domain_sha,"source_refs":{**refs,str(NATIVE_GRADE_CHECK):grade_sha},
            "original_values":endpoint.get("original_values",dict.fromkeys(MARKS)),
            "native_endpoint":endpoint.get("native_endpoint_readout"),
            "native_grade_replay":checked.stdout.splitlines(),
            "scoped_exclusion":{"operator_class":"algebra generated by the reviewed native centre-preserving Q matrices",
                                "closure_rule":"products, sums, adjoints and commutators preserve the centre blocks",
                                "ordered_color_trialities":grades,"required_sector_path":path,
                                "required_blocks":[{"from":start,"to":stop,"operator":leg["field"],
                                                    "hypercharge":leg["Y"]}
                                                   for start,stop,leg in zip(path,path[1:],reversed(leg_data))],
                                "middle_required_block":{"from":path[1],"to":path[2],"operator":leg_data[1]["field"],
                                                         "hypercharge":leg_data[1]["Y"]}},
            "conditional_original_affine_word":word["vacuum_word"],
            "original_root_string_check":radial_symmetry_check(root),
            "weighted_root_doublet":endpoint.get("marked_root_doublet"),
            "next_original_observation":"Identify the two distinct full Cartan weights in colour 2 and the neutral zero-mode link, or bind the charged horizontal doublet directly on the original reflected quotient. Then test the whole ordered word on a held-out two-time rectangle with the same absolute reference.",
            "control":control(),"complex_control":complex_control(),"physical_source_derived":False}
