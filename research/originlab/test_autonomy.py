import unittest

from bridges import exact_tetrahedron, tensor_control, run as bridges
from regulator import certify
from scenarios import run
from stitch import derive, frontier
from synthesis import equivalent_examples, synthesize


class AutonomousResearchTests(unittest.TestCase):
    def test_cylinder_and_moebius_distinguish_equal_betti(self):
        cylinder = run("cylinder",3,"glued_boundary")
        moebius = run("moebius",3,"glued_boundary")
        self.assertEqual(cylinder["invariants"]["betti_Q"],moebius["invariants"]["betti_Q"])
        self.assertEqual(cylinder["invariants"]["boundary_components"],2)
        self.assertEqual(moebius["invariants"]["boundary_components"],1)
        self.assertFalse(moebius["invariants"]["orientable"])

    def test_faces_and_holonomy(self):
        torus = run("torus",3,"face_gluing")
        klein = run("klein",3,"oriented_local_system")
        self.assertEqual(torus["invariants"]["betti_Q"],[1,2,1])
        self.assertEqual(klein["invariants"]["betti_F2"],[1,2,1])
        self.assertEqual(klein["invariants"]["betti_Q"],[1,1,0])
        self.assertGreater(klein["twist_support"],0)
        self.assertEqual(torus["response"]["fixed_modes"],2)
        self.assertGreater(run("circle",3,"hidden_history")["response"]["omitted_path_rank"],0)

    def test_real_backends(self):
        self.assertEqual(equivalent_examples()["status"],"EXACT_TYPED_IDENTITIES")
        synthesis = synthesize(limit=4)
        self.assertEqual(len({tuple(r["truth_table"]) for r in synthesis["rules"]}),4)
        self.assertGreater(synthesis["questions"][0]["pairs_separated"],0)
        self.assertEqual(certify(32,4,"1/10")["status"],"CERTIFIED_ENCLOSURE")

    def test_composition_keeps_source_identity(self):
        rule = [("join",("charge","regions"),"shared")]
        self.assertNotIn("shared",derive([("a","charge","x"),("b","regions","y")],rule)["a"])
        self.assertIn("shared",derive([("a","charge","x"),("a","regions","y")],rule)["a"])
        result = frontier()
        self.assertEqual(len(result["common_origin_routes"]),2)
        self.assertIn("original_weighted_form",result["shared_source_obligations"])

    def test_uploaded_bridges(self):
        self.assertEqual(exact_tetrahedron()["commutator_norm"],"sqrt(2)/12")
        self.assertTrue(tensor_control()["operator_positive_Pauli"]["positive_block"])
        self.assertFalse(tensor_control()["operator_positive_Pauli"]["commutative_finite"])
        self.assertEqual(bridges(8)["charged_finite_control"]["relative_phase"],"I")


if __name__ == "__main__":
    unittest.main()
