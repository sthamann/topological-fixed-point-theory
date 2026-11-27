"""Tests for the theory-map reference checker."""
import unittest

import map_refs


class MapReferences(unittest.TestCase):
    def test_excerpt_is_found_at_the_cited_lines_or_once_elsewhere(self):
        lines = ["header", "AX.P1.01,c3 = 1/(8 pi) boundary seam kernel", "AX.P2.01,g_car = 5 carrier", "tail"]
        self.assertEqual(map_refs.locate(lines, "AX.P1.01,c3 = 1/(8 pi) boundary …", 2, 2), (2, 2))
        self.assertEqual(map_refs.locate(lines, "AX.P2.01,g_car = 5 carrier", 2, 2), (3, 3))
        self.assertEqual(map_refs.locate(lines, "header AX.P1.01,c3", 1, 2), (1, 2))

    def test_a_vanished_or_ambiguous_excerpt_is_not_refreshed(self):
        lines = ["same text", "same text", "other"]
        self.assertIsNone(map_refs.locate(lines, "text that is gone", 1, 1))
        self.assertIsNone(map_refs.locate(lines, "same text", 3, 3))

    def test_installed_map_is_current(self):
        report = map_refs.audit()
        self.assertGreater(report["checked"], 100)
        self.assertEqual(report["needs_review"], [])
        self.assertEqual(report["refreshed"], [], "run python3 semantic_explorer/map_refs.py --refresh")


if __name__ == "__main__":
    unittest.main()
