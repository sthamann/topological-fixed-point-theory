import sympy as sp
from f1_family_index import controls, index_data, projectors, finite_toeplitz, paired_completion

def test_compensating_winding():
    q,p=projectors();z=sp.symbols("z",nonzero=True)
    assert (q+z*p).det().factor()==z**3
    assert (z**-2*q+z*p).det().factor()==z
    assert index_data(-2,1)["total_index"]==-1

def test_conjugate_cancel():
    for k0,k3 in [(0,1),(-2,1),(1,0),(3,-2)]:
        assert index_data(k0,k3)["total_index"]+index_data(-k0,-k3)["total_index"]==0

def test_finite_boundary_unitary():
    for n in (3,5):
        t=finite_toeplitz(n,-2,1)
        assert 4*n-t.rank()==4*n-t.T.rank()==5
        u=paired_completion(n)
        assert u.T*u==sp.eye(8*n) and u*u.T==sp.eye(8*n)

def test_original_not_proved():
    result=controls()
    assert result["status"]=="EXACT_FINITE_FAMILY_INDEX_OBSTRUCTION"
    assert result["original_r113_values"] is None
    assert result["physical_source_identity_proved"] is False
    assert result["unit_determinant_winding_with_three_relative"]==dict(k0=-2,k3=1)


def test_nontrivial_winding_is_not_an_L0_internal_symmetry():
    from f1_family_index import radial_time_covariance
    data=radial_time_covariance()
    assert [row["commutes_with_L0"] for row in data]==[True,False,False,False]


def test_same_four_marks_and_determinant_cannot_fix_relative_index():
    from f1_family_index import finite_mark_blindness
    evidence=finite_mark_blindness()
    assert evidence["same_four_mark_matrices"]
    assert evidence["same_full_determinant_function"]=="z"
    assert evidence["same_positive_reflected_gram"]=="I4"
    assert evidence["separated_relative_indices"]==[0,-3]


def test_holomorphic_one_sheet_forces_no_three_family_winding():
    from f1_family_index import analytic_one_sheet_obstruction
    receipt=analytic_one_sheet_obstruction()
    assert receipt["allowed_zero_counts"]==dict(common=1,relative=0)
    assert receipt["relative_family_index_from_one_sheet"]==0
    assert receipt["mixed_family_symbol_requires_reference_pole_order_at_least"]==2
    assert receipt["original_physical_analyticity_proved"] is False


def test_four_marks_and_determinant_remain_equal_during_SU4_unwinding():
    from f1_family_index import based_su4_unwinding
    receipt=based_su4_unwinding()
    assert receipt["same_four_mark_values_throughout"]
    assert receipt["same_full_determinant_throughout"]
    assert receipt["destroys_fixed_S4_block_equivariance"]
    assert receipt["physical_mixing_admissibility_proved"] is False


def test_winding_sampling_requires_original_derivative_bound():
    from f1_family_index import certified_winding_sampling
    receipt=certified_winding_sampling()
    assert not receipt["source_derivative_bound_supplied"]
    assert not receipt["original_winding_certified"]
    assert [x["certified"] for x in receipt["test_examples"]]==[True,True,False]
