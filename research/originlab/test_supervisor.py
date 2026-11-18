import json
import tempfile
import unittest
from pathlib import Path

from autonomy import HERE, cache_key, completeness, enqueue, result_rows, setup, worker, work_order


class SupervisorTests(unittest.TestCase):
    def test_unchanged_inputs_reuse_a_result_in_a_later_round(self):
        spec = {"lane":"regulator","grid":24,"cutoff":4,"time":"1/10"}
        with tempfile.TemporaryDirectory() as tmp:
            db = Path(tmp)/"queue.sqlite3"
            setup(db)
            enqueue(db,"first",[(1,spec)])
            self.assertEqual(worker(db,"first",HERE.parents[1],90,"basis"),0)
            enqueue(db,"second",[(1,spec)])
            self.assertEqual(worker(db,"second",HERE.parents[1],90,"basis"),1)
            enqueue(db,"third",[(1,spec)])
            self.assertEqual(worker(db,"third",HERE.parents[1],90,"changed source"),0)
            first,second,third = (result_rows(db,phase)[0] for phase in ("first","second","third"))
            self.assertEqual(first["result"],second["result"])
            self.assertEqual(first["result"],third["result"])
            self.assertEqual(second["queue_status"],"done")

    def test_proof_and_model_proposal_jobs_always_run(self):
        self.assertIsNone(cache_key("basis",json.dumps({"lane":"f1_solution","stage":3})))
        self.assertIsNone(cache_key("basis",json.dumps({"lane":"adaptive_dsl","sol_timeout":45})))
        self.assertIsNone(cache_key(None,json.dumps({"lane":"rewrite"})))
        self.assertNotEqual(cache_key("a",json.dumps({"lane":"rewrite"})),
                            cache_key("b",json.dumps({"lane":"rewrite"})))

    def test_work_order_names_the_inputs_no_stage_can_produce(self):
        blocked = {"status":"BLOCKED","missing":["original_blocked_fiber_packet"],
                   "first_operation":"evaluate_original_full_blocked_endpoint_source"}
        order = work_order(blocked,{"formal_acceptance_registered":False},
                           {"matching":[],"rejected":[]},[],backend="elan")
        self.assertFalse(order["f1_proof_reachable"])
        self.assertEqual([item["id"] for item in order["items"]],
                         ["native_acceptance_contract","original_source_packet"])
        self.assertEqual(len(order["items"][0]["open_gate_statements"]),8)
        self.assertEqual(order["items"][1]["blocks"],["G3","G4"])

    def test_a_registered_contract_and_backend_make_the_proof_reachable(self):
        failed = [{"spec":{"lane":"scenario"},"result":{"error_type":"TimeoutExpired"}}]
        order = work_order({"status":"EXACT_SUPPLIED_CONTRACTION"},{"status":"F1_CONTRACT_OPEN"},
                           {"matching":[],"rejected":[]},failed,backend="elan")
        self.assertTrue(order["f1_proof_reachable"])
        self.assertEqual([item["id"] for item in order["items"]],["failed_jobs"])
        missing = work_order({},{"status":"F1_CONTRACT_OPEN"},{"matching":[],"rejected":[]},[],backend="")
        self.assertFalse(missing["f1_proof_reachable"])

    def test_completeness_is_computed_from_the_registry_and_the_job_grammar(self):
        from source import NATIVE
        root = HERE.parents[1]
        audit = completeness(root,root/NATIVE,None,{"matching":[],"rejected":[]},
                             {"formal_acceptance_registered":False})
        self.assertFalse(audit["complete"])
        self.assertEqual(audit["gates_without_acceptance_predicate"],[f"G{n}" for n in range(8)])
        self.assertEqual(audit["source_reviewed_packets"],0)
        self.assertEqual(audit["lanes_still_growing_with_stage"],["coherent_origin","f1_solution"])


if __name__ == "__main__":
    unittest.main()
