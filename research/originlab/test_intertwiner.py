from intertwiner import solve
from source_jets import gaussian_jet


def test_rectangular_inclusion_and_no_map_certificate():
    source = [[1,0],[0,1]]
    base = {"domain_metric":source,"target_metric":[[1,0,0],[0,1,0],[0,0,1]],
            "relations":{"charge":{"domain":[[1,0],[0,2]],"target":[[2,0,0],[0,1,0],[0,0,3]]}}}
    positive = solve(base)
    assert positive["status"]=="EXACT_RECTANGULAR_ISOMETRY"
    assert positive["linear_solution_dimension"]==2
    base["relations"]["charge"]["target"] = [[3,0,0],[0,3,0],[0,0,3]]
    impossible = solve(base)
    assert impossible["status"]=="EXCLUDED_RECTANGULAR_INTERTWINER"
    assert len(impossible["conflicting_relations"])==1


def test_gaussian_contact_and_singular_reference():
    data = {"input_role":"declared_finite_oracle","kinetic":[[2]],"U":[[3]],"Y":[[5]],"UY":[[7]],
            "scalar_weight_jet":["11","13","17","19"],"reference":"2"}
    result = gaussian_jet(data)
    assert result["unnormalized_numerators"]=={"Z0":"22","ZU":"59","ZY":"89","ZUY":"231"}
    assert result["kinetic_mixed_contact_numerator"]=="77"
    assert result["connected_UY"]=="-169/484"
    data["kinetic"]=[[0]]
    singular = gaussian_jet(data)
    assert singular["base_kinetic_singular"] and singular["connected_UY"] is None
    assert singular["unnormalized_numerators"]["ZUY"]=="193"
