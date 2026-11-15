import unittest
from pathlib import Path
import os

from native_weak_bridge import NATIVE, replay as native_weak_bridge
from r113_synthesis import control, from_input, history_range, spectral_pair


class R113SynthesisTests(unittest.TestCase):
    def test_native_weak_graft_fails_even_with_free_coloured_matrices(self):
        root = Path(os.environ.get("TFPT_REPO_ROOT",Path(__file__).resolve().parents[2]))
        report = native_weak_bridge(root,root/NATIVE)
        self.assertEqual(report["rank_per_colour_and_weak_weight"],[[[3]*5]*2]*3)
        self.assertEqual(report["native_root_string"]["family_block_allowed_parameters"],0)
        self.assertEqual([row["native_polar_rank"] for row in report["polar_car_parity"]],[48,48])
        self.assertTrue(all(row["even_car_allowed_odd_weight_space"]==0 for row in report["polar_car_parity"]))
        self.assertEqual(report["det16_spectator"]["maximum_time_commuting_odd_operator_rank"],65534)
        self.assertEqual(report["det16_spectator"]["odd_ground_rank"],0)
        self.assertEqual(report["det16_spectator"]["odd_ground_return"],"exp(-2*t) * ||B phi||^2")
        self.assertTrue(all(value is None for value in report["original_values"].values()))

    def test_synthesizes_and_separates_two_positive_transfers(self):
        result = control()
        self.assertEqual(result["compatible_mixed_spectra"],2)
        self.assertEqual(result["next_observation"]["step"],1)
        self.assertEqual(set(result["next_observation"]["witness_values"]),{"7/8","9/8"})
        for witness in result["witnesses"]:
            self.assertEqual(witness["moment_zero"],[["3","2"],["2","6"]])
            replay = witness["positive_transfer_replay"]
            self.assertEqual(replay["status"],"EXACT_ACCESSIBLE_POSITIVE_TRANSFER")
            self.assertEqual(replay["rank"],3)
            self.assertEqual(replay["matched_held_out_moments"],2)

    def test_feedback_resolves_phase_in_declared_class(self):
        data = {"kind":"real_charged_spectrum","input_role":"declared_finite_oracle",
                "spectrum":["1/2","1/4","1/8"],
                "diagonal_observations":[[3,6],["7/8","19/8"],["21/64","69/64"]],
                "mixed_observations":{"0":2,"1":"9/8"}}
        result = from_input(data)
        self.assertEqual(result["status"],"MIXED_RESPONSE_FIXED")
        self.assertEqual(result["mixed_spectral_weights"],["2","1","-1"])
        data["mixed_observations"]["1"] = 1
        self.assertEqual(spectral_pair(data)["status"],"EXCLUDED_IN_DECLARED_SPECTRAL_CLASS")

    def test_collinear_marks_have_enough_powers_for_hankel_replay(self):
        result = spectral_pair({"spectrum":["1/2","1/4","1/8"],
                                "diagonal_observations":[[3,3],["7/8","7/8"],["21/64","21/64"]],
                                "mixed_observations":{}})
        self.assertEqual(result["status"],"MIXED_RESPONSE_UNDERDETERMINED")
        self.assertEqual(result["next_observation"]["step"],0)

    def test_fixed_history_interval_preserves_reference(self):
        data = {"histories":["00","01","10","11"],"reference":2,
                "constraints":{"Z0":{"coefficients":[2,2,2,2],"value":1},
                               "ZU":{"coefficients":[0,0,2,2],"value":"1/2"},
                               "ZY":{"coefficients":[0,2,0,2],"value":"1/2"}},
                "mixed_coefficients":[0,0,0,"4/3"]}
        result = history_range(data)
        self.assertEqual(result["mixed_response_interval"],["0","1/3"])
        self.assertEqual({tuple(w["weights"]) for w in result["witnesses"]},
                         {("0","1/2","1/2","0"),("1/2","0","0","1/2")})
        result = history_range({**data,"mixed_coefficients":[0,0,2,2]})
        self.assertEqual(result["status"],"MIXED_RESPONSE_FIXED")
        self.assertEqual(result["mixed_response_interval"],["1/2","1/2"])
        data["constraints"]["ZU"]["value"] = 2
        self.assertEqual(history_range(data)["status"],"EXCLUDED_IN_DECLARED_HISTORY_CLASS")
