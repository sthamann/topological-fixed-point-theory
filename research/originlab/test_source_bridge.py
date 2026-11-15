from pathlib import Path
import os

from flint import arb
import numpy as np

from source_bridge import neutral_return_budget, regional_transport, replay
from spatial_memory import spatial_cut


def test_native_regional_span_has_only_one_commuting_direction():
    regional = regional_transport()
    assert regional["span_rank"] == 3
    assert regional["bracket_wedge_determinant"] == "2"
    assert regional["max_independent_commuting_rank"] == 1
    assert regional["witness"]["value_at_half"] == "-1/2"


def test_neutral_rank_shortfall():
    assert neutral_return_budget(1, 1)["status"] == "NECESSARY_RANK_OK"
    assert neutral_return_budget(3, 1)["min_additional_neutral_rank"] == 2


def test_spatial_return_applies_the_adjoint_in_opposite_directions():
    edge = np.array([[1, 2j], [0, 1]], dtype=complex)
    diagonal = {n: np.array([[n, 1], [1, -n]], dtype=complex) for n in range(-3, 3)}
    modes = tuple(diagonal)
    full = np.block([[diagonal[n] if n == m else edge if m == n + 1 else
                      edge.conj().T if m == n - 1 else np.zeros_like(edge)
                      for m in modes] for n in modes])
    z = 2j
    direct = np.linalg.inv(z * np.eye(len(full)) - full)[4:8, 4:8]
    result, _, _ = spatial_cut(z, diagonal.__getitem__, edge, 2)
    np.testing.assert_allclose(result, direct, atol=1e-14)
    assert neutral_return_budget(3, 1)["status"] == "EXCLUDED_SCOPED"


def test_original_fixed_history_keeps_the_scalar_separate_from_the_trace():
    root = Path(os.environ.get("TFPT_REPO_ROOT", Path(__file__).resolve().parents[2]))
    result = replay(root)
    fibre = result["charged_fixed_history"]
    assert fibre["pairs"] == 42
    for key, published_rounding in (
        ("log_charged_trace_over_S", "[-19.7198830602 +/- 1e-10]"),
        ("log_sewn_over_polar_scalar", "[-19.5442075498 +/- 1e-10]"),
    ):
        assert arb(fibre[key]).overlaps(arb(published_rounding))
    assert arb(fibre["log_unmarked_scalar_over_S"]) < arb(fibre["log_charged_trace_over_S"])
    assert fibre["physical_R113"] == dict.fromkeys(("Z0", "F0", "CUY"))
    profiles = result["pointwise_affine_profiles"]
    assert profiles["two_evaluation_sectors"]["charged_rank_outside_unit_return"] == 0
    assert profiles["transfer_mixing_evaluation_sectors"]["noncommuting_profile_indices"] == [1]
    rejected = profiles["scalar_average_with_extra_charged_path"]
    assert rejected["profile_product_mismatches"] == []
    assert rejected["profile_return_mismatches"] == [[1, 1]]
    assert rejected["charged_rank_outside_unit_return"] == 1
    spatial = result["native_spatial_memory"]
    assert spatial["exact_edge_rank"] == 96
    assert [row["krylov_rank_through_depth"] for row in spatial["depths"][:3]] == [192, 384, 576]
    assert spatial["depths"][1]["first_different_return_power"] == 4
    assert not spatial["finite_static_exact_realization_exists"]
    assert not spatial["absolute_weighted_R113_produced"]
    assert not result["physical_solution"]
