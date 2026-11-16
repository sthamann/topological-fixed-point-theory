"""Arithmetic tests for the reflected pair probe at the R113 consumer."""

import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import sympy as sp

from source import COMPONENTS, evaluate
from verify import digest


class ReflectedPairCapTests(unittest.TestCase):
    def setUp(self):
        folder = TemporaryDirectory()
        self.addCleanup(folder.cleanup)
        self.root = Path(folder.name)
        evidence = self.root/"evidence.txt"
        evidence.write_text("arithmetic control only\n")
        ref = {"path":evidence.name,"sha256":digest(evidence)}
        self.native = self.root/"native.json"
        self.base = {"domain_id":"paired-fixture","source_originals":[ref]}
        self.native.write_text(json.dumps(self.base))
        self.anchor = {"path":self.native.name,"sha256":digest(self.native)}
        self.attachment = self.root/"source.json"
        self.embedding = sp.Matrix([[1,0,0,1],[0,1,0,0],[0,0,1,0]])
        self.metric = sp.diag(2,2,3)
        diagonal = sp.diag(3,5,7)
        selector = sp.Matrix([[1,1,0],[1,1,0],[0,0,0]])/2
        self.probe = {
            "source_id":f"paired-fixture:{self.anchor['sha256']}","time_id":"two-windows",
            "operator_dictionary_binding":"fixture","two_window_order_binding":"fixture",
            "absolute_reference_binding":"fixture","paired_support_binding":"fixture",
            "ordered_word":["A","B","A","B"],
            "paired_projector_form":self.form(sp.diag(1,0,0)),
            "ordered_abab_form":self.form(sp.diag(-3,5,7)),
        }
        self.bundle = {
            "schema":"tfpt.original-blocked-source.v1","domain_id":self.base["domain_id"],
            "native_inputs_sha256":self.anchor["sha256"],"input_role":"arithmetic_control_only",
            "matrix_convention":"sesquilinear_forms_in_original_word_basis",
            "ensemble_scope":"whole_original_ensemble",
            "representation":"original_exact_blocked_endpoint_source",
            "source_refs":[ref],"component_bindings":dict.fromkeys(COMPONENTS,"fixture"),
            "retained_words":[[],["pair"],["other"],["dependent"]],
            "time":{"kind":"original_euclidean_ordered","id":"two-windows",
                    "T":2,"t_U":"1/2","t_Y":"3/2"},
            "gram_form":self.form(sp.eye(3)),"physical_selector_form":self.form(selector),
            "boundary_pairing_form":self.form(sp.diag(1,0,0)),
            "scalar_factor":[2,0],"reference_scalar":[5,0],
            "ordered_source_forms":{
                "Z0":self.form(diagonal),"ZU":self.form(sp.eye(3)),
                "ZY":self.form(sp.diag(2,3,0)),
                "ZUY":self.form(sp.Matrix([[0,1,0],[1,0,0],[0,0,0]])),
            },
            "reflected_pair_cap":self.probe,
        }

    def form(self, operator):
        matrix = self.embedding.H*self.metric*operator*self.embedding
        return [[[str(sp.re(value)),str(sp.im(value))] for value in row] for row in matrix.tolist()]

    def read(self):
        data = dict(self.base,original_domain_binding=self.anchor,whole_source_bundle=self.bundle)
        self.attachment.write_text(json.dumps(data))
        return evaluate(self.root,self.attachment)

    def test_original_ordered_word_and_selection_with_metric_and_null_quotient(self):
        result = self.read()
        probe = result["reflected_pair_cap"]
        self.assertEqual(result["quotient_rank"],3)
        self.assertEqual(probe["status"],"EXACT_PAIRED_WORD_MATCH")
        self.assertEqual(probe["bare_numerator"],"4")
        self.assertEqual(result["values"]["Z0"],"4/5")
        self.assertEqual(probe["absolute_reference"],"5")
        self.assertEqual(probe["quotient_basis_indices"],[0,1,2])
        self.assertEqual(probe["paired_numerator"],"3/2")
        self.assertEqual(probe["paired_ratio"],"3/8")
        self.assertEqual(probe["normalized_ABAB"],"1/4")
        self.assertEqual(probe["predicted_ordered_numerator"],"1")
        self.assertEqual(probe["measured_ordered_numerator"],"1")
        self.assertEqual(probe["measured_ordered_value_over_reference"],"1/5")
        self.assertEqual(probe["measured_operator"],probe["predicted_operator"])
        self.assertEqual(probe["restricted_first_responses"],
                         {"ZY":{"rank":1,"selected_numerator":"1"},
                          "ZUY":{"rank":0,"selected_numerator":"0"}})
        self.assertFalse(probe["physical_source_identity_proved"])

    def test_same_contraction_hides_preselector_matrix_residual(self):
        self.probe["ordered_abab_form"] = self.form(sp.diag(-3,5,8))
        probe = self.read()["reflected_pair_cap"]
        self.assertEqual(probe["status"],"PAIRED_WORD_MATRIX_MISMATCH")
        self.assertEqual(probe["contracted_residual"],"0")
        self.assertEqual(probe["selected_matrix_residual_rank"],0)
        self.assertEqual(probe["matrix_residual"][2][2],"1")

    def test_same_history_squared_defect_reads_before_selector(self):
        data = {"original_mixed_moment_binding":"fixture", "form":self.form(sp.diag(0,0,3))}
        self.probe["same_history_defect"] = data
        result = self.read()["reflected_pair_cap"]["same_history_defect"]
        self.assertEqual(result,{"form_rank":1,"selected_rank":0,"numerator":"0",
                                 "derived_from_fiber_transfers":False,
                                 "value_over_reference":"0","nonzero_hidden_by_selector":True,
                                 "nonzero_selected_hidden_by_contraction":False})
        data["form"] = self.form(sp.diag(2,0,0))
        result = self.read()["reflected_pair_cap"]["same_history_defect"]
        self.assertEqual(result["selected_rank"],1)
        self.assertEqual(result["value_over_reference"],"1/5")
        data["form"] = self.form(sp.diag(0,0,-1))
        with self.assertRaisesRegex(ValueError,"must be Hermitian PSD"):
            self.read()

    def test_same_history_squared_defect_from_mixed_moment_blocks(self):
        projection = self.form(sp.diag(1,1,0))
        data = {"original_mixed_moment_binding":"fixture",
                "moments":{"input_norm":self.form(sp.eye(3)),
                           "transfer_overlap":projection,"transfer_norm":projection}}
        self.probe["same_history_defect"] = data
        result = self.read()["reflected_pair_cap"]["same_history_defect"]
        self.assertEqual((result["form_rank"],result["selected_rank"]),(1,0))
        data["moments"]["transfer_norm"] = self.form(sp.zeros(3))
        with self.assertRaisesRegex(ValueError,"lack a common positive Gram"):
            self.read()

    def test_nonreducing_bare_transfer_does_not_generate_a_prediction(self):
        self.bundle["ordered_source_forms"]["Z0"] = self.form(sp.Matrix([[3,1,0],[1,5,0],[0,0,7]]))
        probe = self.read()["reflected_pair_cap"]
        self.assertEqual(probe["status"],"PAIR_CAP_DOES_NOT_REDUCE_BARE_TRANSFER")
        self.assertEqual(probe["bare_commutator_rank"],2)
        self.assertEqual(probe["measured_ordered_numerator"],"1")
        self.assertIsNone(probe["predicted_ordered_numerator"])
        self.assertIsNone(probe["paired_ratio"])

    def test_refuses_bad_binding_and_null_relations(self):
        self.probe["time_id"] = "another-time"
        with self.assertRaisesRegex(ValueError,"same original source"):
            self.read()
        self.probe["time_id"] = "two-windows"
        self.probe["paired_projector_form"][3][3] = ["42",0]
        with self.assertRaisesRegex(ValueError,"does not descend"):
            self.read()
        self.probe["paired_projector_form"] = self.form(sp.diag(2,0,0))
        with self.assertRaisesRegex(ValueError,"metric-orthogonal paired projector"):
            self.read()

    def test_missing_probe_remains_blocked_for_native_and_sparse_inputs(self):
        self.bundle.pop("reflected_pair_cap")
        for result in (evaluate(self.root,self.native),self.read()):
            probe = result["reflected_pair_cap"]
            self.assertEqual(probe["status"],"BLOCKED")
            self.assertIsNone(probe["normalized_ABAB"])
            self.assertIsNone(probe["predicted_ordered_numerator"])
