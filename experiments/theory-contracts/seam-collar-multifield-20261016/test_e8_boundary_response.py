import json
from pathlib import Path
import unittest

import sympy as sp

import e8_boundary_response as e8


class BoundaryDictionaryTests(unittest.TestCase):
    def test_published_probe_replays_with_its_own_charges(self):
        report = json.loads(Path(__file__).with_name("results").joinpath(
            "e8-boundary-response.json").read_text())
        probe = report["nonabelian_probe"]
        charges = tuple(tuple(sp.Rational(x) for x in a) for a in probe["charges"])
        self.assertEqual(tuple(map(sum, zip(*charges))), e8.ZERO)
        self.assertEqual([e8.grade(a) for a in charges], probe["grades"])
        self.assertEqual(str(e8.vertex_response(charges, (3, 2, 1))),
                         probe["three_point_at_3_2_1"])

    def test_elementary_car_vector_is_not_a_local_e8_vacuum_field(self):
        vector = (sp.S.One,) + (sp.S.Zero,) * 7
        with self.assertRaisesRegex(ValueError, "outside"):
            e8.vertex_response((vector, e8.neg(vector)), (2, 1))

    def test_coincident_sources_need_contact_data(self):
        a = (e8.HALF,) * 8
        with self.assertRaisesRegex(ValueError, "contact"):
            e8.current_response((a, e8.neg(a)), (2, 1), ((a, 2),))

    def test_two_cartan_insertions_have_charge_and_level_terms(self):
        a = (e8.HALF,) * 8
        response = e8.adjoint_sign(a) * e8.current_response(
            (a, e8.neg(a)), (3, 1), ((a, 5), (a, 4)))
        # Root pair=1/4. Means: mu(5)=1/2, mu(4)=4/3.
        # Level-one Cartan contraction a.a/(5-4)^2=2.
        self.assertEqual(response, sp.Rational(2, 3))


if __name__ == "__main__":
    unittest.main()
