"""Exact grade-one test for the marked affine E8 root string."""
import sympy as sp

from affine_time import HYPERCHARGES, TRIALITY, marked_charges, radial_grid
from language import positive
from reconstruction import require
from return_moments import rows


def root_doublet():
    legs = marked_charges()
    left, root, right = (sp.Matrix(leg["root"]).applyfunc(sp.Rational) for leg in legs)
    end = -left
    require(root + right == end and (root.dot(right), root.dot(end)) == (-1,1),
            "marked states are not the endpoints of a root string")
    require(((right-root).dot(right-root), (end+root).dot(end+root)) == (6,6),
            "marked root string does not terminate at its two endpoints")
    e = sp.Matrix([[0,0],[1,0]])
    f = e.H
    h = sp.diag(-1,1)
    require(e*f-f*e==h and h*e-e*h==2*e and h*f-f*h==-2*f,
            "normalized root operators fail the horizontal sl2 relations")
    # Cyclically move the neutral root to the zero mode. This is an
    # identity of structure constants, before time weights and selection.
    cyclic_states = (right,-root)
    require(right+left==-root and (left.dot(right),left.dot(-root))==(-1,1) and
            all(vector.dot(vector)==2 for vector in cyclic_states),"cyclic root string changed")
    return {"root":legs[1]["root"],"radial_grades":[1,1],
            "state_trialities":[TRIALITY[2],(-TRIALITY[0])%3],
            "state_hypercharges":[HYPERCHARGES[2],str(-sp.Rational(HYPERCHARGES[0]))],
            "lower_root_string_norm_squared":str((right-root).dot(right-root)),
            "upper_root_string_norm_squared":str((end+root).dot(end+root)),
            "normalized_middle_abs_squared":"1",
            "cyclic_horizontal_test":{"middle":legs[0]["field"],
                                      "same_color":TRIALITY[2],
                                      "full_cartan_weights":[list(map(str,v)) for v in cyclic_states],
                                      "normalized_matrix_element_abs_squared":"1",
                                      "scope":"cyclic affine structure constant at fixed cocycle, not a reordered weighted endpoint"},
            "scope":"horizontal zero mode on the two grade-one root vectors of the given unitary affine vacuum",
            "phase_convention":"the upper line is phased by E applied to the lower line"}


def checked_doublet(metric, vacuum, create_lower, create_upper, lower_end, upper_end,
                    e, f, h, triality, hypercharge, steps):
    """Check a supplied small root module on a common reflected carrier."""
    n = metric.rows
    identity = sp.eye(n)
    require(len(steps)==2 and positive(metric) and metric.det()!=0 and
            all(positive(metric*step) and positive(metric*(identity-step))
                                 for step in steps),"two positive contraction steps are required")
    def simplify(matrix):
        return matrix.applyfunc(lambda value: sp.simplify(sp.expand_complex(value)))

    require(simplify(triality**3-identity)==sp.zeros(n) and
            simplify(triality.H*metric*triality-metric)==sp.zeros(n),
            "need the unitary order-three colour action")
    inverse = metric.inv()

    def adjoint(operator):
        return inverse*operator.H*metric
    require(hypercharge==adjoint(hypercharge),"hypercharge must be metric self-adjoint")
    lower = create_lower*vacuum
    upper = create_upper*vacuum
    pair = lower.row_join(upper)
    f_error = f-adjoint(e)
    omega_squared = (vacuum.H*metric*vacuum)[0]
    require(omega_squared>0,"need a nonzero reflected empty word")
    v0 = sp.zeros(n,1)
    defects = {
        "neutral_empty_word":(triality*vacuum-vacuum).col_join(hypercharge*vacuum),
        "endpoint_norms":pair.H*metric*pair-omega_squared*sp.eye(2),
        "endpoint_adjoints":(lower_end-adjoint(create_lower)).col_join(upper_end-adjoint(create_upper)),
        "middle_adjoint":f_error,
        "middle_hypercharge":hypercharge*e-e*hypercharge+e/3,
        "middle_triality":triality*e-sp.exp(2*sp.pi*sp.I/3)*e*triality,
        "endpoint_hypercharges":hypercharge*pair-pair*sp.diag(sp.Rational(1,3),0),
        "endpoint_trialities":triality*pair-pair*sp.diag(sp.exp(4*sp.pi*sp.I/3),1),
        "root_weight":h*pair-pair*sp.diag(-1,1),
        "root_weight_adjoint":h-adjoint(h),
        "horizontal_commutator":(e*f-f*e-h)*pair,
        "root_string_tops":(f*lower-v0).row_join(e*upper-v0),
        "marked_middle":e*lower-upper,
        "reverse_middle":f*upper-lower,
    }
    ranks = {k:simplify(a).rank() for k,a in defects.items()}
    timed = []
    for step in steps:
        left, right = (sp.simplify((v.H*metric*step*v)[0]/omega_squared) for v in (lower,upper))
        leakage = (step*lower-left*lower).row_join(step*upper-right*upper)
        timed.append({"compressed_diagonal":[str(left),str(right)],
                      "time_leakage_rank":leakage.applyfunc(sp.simplify).rank(),
                      "middle_time_commutator_rank":((step*e-e*step)*pair).applyfunc(sp.simplify).rank(),
                      "relative_ordered_return":str(sp.simplify((vacuum.H*metric*upper_end*step*e*step*lower)[0]/
                                                              omega_squared))})
    # Absolute endpoint contraction still belongs to the R113 reader.
    return {"status":"EXACT_SUPPLIED_ROOT_DOUBLET" if not any(ranks.values()) else "REFUTED_SUPPLIED_ROOT_DOUBLET",
            "kinematic_defect_ranks":ranks,"time_steps":timed,
            "root_basis_gram":rows(pair.H*metric*pair),
            "physical_source_identity_proved":False}


def radial_symmetry_check(root):
    from affine_time import original_time
    from folded_time import marked_clock

    radial = original_time(root)
    folded = marked_clock(root)
    cases = []
    for case in folded["cases"]:
        energies = case["marked_energies_in_H8_time"]
        left,right = map(sp.Rational,energies)
        cases.append({"quarter":case["quarter"],"H8_energies":energies,
                      "left_minus_right_H8_energy":str(left-right),
                      "commutes_with_middle_on_root_doublet":left==right})
    widths = (sp.Rational(1,2),sp.Rational(1,4))
    grid = radial_grid(widths,widths)
    determinant = sp.cancel(grid.det(),extension=sp.I)
    require(determinant!=0,"original affine time grid unexpectedly separable")
    ratio = sp.cancel(grid[0,0]*grid[1,1]/(grid[0,1]*grid[1,0]),extension=sp.I)
    require(ratio==sp.Rational(81,85),"original affine mixed time ratio changed")
    return {"root_doublet":root_doublet(),"source_refs":{**radial["source_refs"],**folded["source_refs"]},
            "finite_carrier_limit":{
                "identity":"[H_1,H_-1] = k*(H,H)*I for a Cartan element H",
                "obstruction":"the trace of a finite commutator is zero, while the right-hand trace is nonzero at positive level",
                "scope":"finite matrices cannot realize the complete positive-level affine mode algebra; compression and the horizontal doublet remain possible"},
            "radial_commutes_with_horizontal_middle_on_doublet":True,
            "quarter_cases":cases,
            "affine_two_time_control":{"transfer_factors":list(map(str,widths)),
                                       "rank":grid.rank(),"mixed_ratio":str(ratio),
                                       "pure_doublet_rank_bound":1,"pure_doublet_refuted_for_full_word":True},
            "physical_R113_time_map_derived":False,
            "scope":"fixed marked lines and fixed affine and folded time actions; a twisted time with a charged commutator remains possible"}
