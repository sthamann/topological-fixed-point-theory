from pathlib import Path
import os

from flint import arb

from source_bridge import neutral_return_budget, regional_transport, replay


def test_native_regional_span_has_only_one_commuting_direction():
    regional = regional_transport()
    assert regional["span_rank"] == 3
    assert regional["bracket_wedge_determinant"] == "2"
    assert regional["max_independent_commuting_rank"] == 1
    assert regional["witness"]["value_at_half"] == "-1/2"


def test_neutral_rank_shortfall():
    assert neutral_return_budget(1, 1)["status"] == "NECESSARY_RANK_OK"
    assert neutral_return_budget(3, 1)["min_additional_neutral_rank"] == 2
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
    assert not result["physical_solution"]
