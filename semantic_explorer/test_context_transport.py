import json
from pathlib import Path
import tempfile
import unittest

from semantic_explorer.build_context import atomic


class ContextTransportTests(unittest.TestCase):
    def test_huge_exact_contract_integers_remain_lossless(self):
        huge = 10 ** 879 + 37
        original = {"attrs": {"exact": [huge, -huge], "normal": 3,
                              "flag": True, "optional": None}}
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "context.json"
            atomic(path, original)
            parsed = json.loads(path.read_text())
        self.assertEqual(parsed["attrs"]["exact"], [str(huge), str(-huge)])
        self.assertEqual(original["attrs"]["exact"], [huge, -huge])
        self.assertEqual(parsed["attrs"]["normal"], 3)
        self.assertIs(parsed["attrs"]["flag"], True)
        self.assertIsNone(parsed["attrs"]["optional"])

    def test_rust_integer_boundaries_keep_native_numbers(self):
        values = [-(1 << 63), (1 << 64) - 1,
                  -(1 << 63) - 1, 1 << 64]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "context.json"
            atomic(path, values)
            parsed = json.loads(path.read_text())
        self.assertEqual(parsed, values[:2] + [str(x) for x in values[2:]])


if __name__ == "__main__":
    unittest.main()
