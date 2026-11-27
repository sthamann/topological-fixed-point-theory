"""Tests for the ledger knowledge map."""
import unittest

import ledger_map as lm


class LedgerMap(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rows = lm.load()

    def test_classes(self):
        self.assertEqual(lm.classify("Identity [I]"), "exakt")
        self.assertEqual(lm.classify("Identity (exact) [E]/[C]"), lm.PARTIAL)
        self.assertEqual(lm.classify("Open (research contract)"), "offen")
        self.assertEqual(lm.classify("Physical (conjecture)"), "bedingt")
        self.assertEqual(lm.classify("Kill-test [X]"), "Kill-Test")
        self.assertEqual(lm.classify("Axiom (declared input)"), "gesetzt")
        self.assertEqual(lm.classify("Firewall (standing discipline rule)"), "Regel oder Synthese")

    def test_every_claim_has_a_class_and_an_area(self):
        self.assertGreater(len(self.rows), 700)
        self.assertTrue(all(r["class"] in lm.CLASSES for r in self.rows))
        exact = sum(1 for r in self.rows if r["class"] in ("exakt", lm.PARTIAL))
        self.assertGreater(exact, len(self.rows) // 2)

    def test_recorded_branches_are_found(self):
        ids = {cid for cid, *_ in lm.find("twistor", self.rows)}
        self.assertTrue({"WOIT.OS.TWISTOR.01", "WOIT.BETA3.DUALITY.01"} <= ids)
        self.assertTrue(lm.find("birefringence", self.rows))

    def test_map_is_current(self):
        self.assertEqual(lm.TARGET.read_text(), lm.render(self.rows))

    def test_atlas_summary_is_current_and_complete(self):
        summary = lm.atlas_summary(self.rows)
        self.assertEqual(lm.ATLAS_SUMMARY.read_text(), lm.summary_text(self.rows))
        self.assertEqual(sum(area["claims"] for area in summary["areas"]), len(self.rows))
        self.assertEqual(len(summary["areas"]), len(lm.AREAS))
        twistor = next(area for area in summary["areas"] if area["title"].startswith("Twistor"))
        self.assertIn("WOIT", twistor["blocks"])


if __name__ == "__main__":
    unittest.main()
