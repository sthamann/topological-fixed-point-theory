"""Regression tests for the proof grammar and generative source adapter."""

from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import sympy as sp

from proof_search import GaussianProducer, LeanKernel, affine_forecast, terms_at_cost


class ProofSearchTests(unittest.TestCase):
    def test_generated_forms_preserve_order_and_cutoff(self):
        producer = GaussianProducer(1)
        windows = tuple(map(sp.Rational,("1/2","2/3","3/5")))
        u,y = ((1,-1),),((1,1),)
        p=(2,1)
        forms=producer.forms(p,p,windows,u,y)
        self.assertEqual(forms,{"Q0":sp.Rational(4,625),"QU":0,"QY":0,"QUY":sp.Rational(8,625)})
        self.assertEqual(producer.forms((*p,0,0),(*p,0,0),windows,u,y),forms)
        self.assertEqual(producer.forms(p,p,windows,y,u)["QUY"],sp.Rational(12,625))
        self.assertEqual(producer.gram(p,(1,1)),0)

    def test_grammar_reaches_nested_constructors(self):
        constructors={"Model.seed":0,"Model.map":1}
        self.assertIn("(Model.map (Model.map Model.seed))",list(terms_at_cost(constructors,3)))
        self.assertNotIn("(Model.map (Model.map Model.seed))",list(terms_at_cost(constructors,2)))

    def test_changed_goal_hash_rejected_before_elaboration(self):
        with TemporaryDirectory() as temp:
            goal=Path(temp)/"FixedGoal.lean"
            goal.write_text("changed")
            kernel=LeanKernel(Path(temp),"unused",1)
            with self.assertRaisesRegex(ValueError,"digest mismatch"):
                kernel.prepare(((goal,"0"*64),))

    def test_affine_producer_generates_existing_charged_word(self):
        result=affine_forecast()
        self.assertEqual(result["empty_charged_form"],"6*u**6*v**6")
        self.assertEqual(result["ordered_charged_form"],"4*u**6*v**4")
        self.assertIsNone(result["absolute_r113_reference"])
