"""Focused regression for explicit metadata-only publication omissions."""
import csv
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch


spec = importlib.util.spec_from_file_location(
    "research_inventory_export", Path(__file__).with_name("build-research-inventory.py")
)
inventory = importlib.util.module_from_spec(spec)
spec.loader.exec_module(inventory)


class ResearchInventoryTests(unittest.TestCase):
    def test_shared_status_has_live_public_evidence(self):
        root = Path(__file__).resolve().parents[1]
        areas = json.loads((root / "lib/atlas-status.json").read_text())
        graph = json.loads((root / "lib/semantic-atlas.json").read_text())
        foundation = json.loads((root / "lib/foundation-index.json").read_text())
        evidence = {r["id"] for r in graph["findings"] + graph["patterns"]}
        evidence.update("claim:" + claim["id"] for claim in foundation["claims"])
        evidence.add("verification/script_registry.csv")
        self.assertEqual(
            {a["id"] for a in areas},
            {"original-source", "compiler", "matter-alpha", "clocks-quantum",
             "gravity", "verification", "completion"},
        )
        for area in areas:
            self.assertIn(area["evidence"], evidence, area["id"])
            self.assertTrue(area["scope_kind"], area["id"])
            self.assertTrue(area["related_evidence"], area["id"])
            self.assertTrue(set(area["related_evidence"]) <= evidence, area["id"])
            for field in ("title", "result", "scope"):
                self.assertTrue(all(area[field].get(lang) for lang in ("de", "en")), area["id"])
            self.assertTrue(area["href"].startswith("/") and not area["href"].startswith("//"))


    def test_resolved_tasks_and_registered_modules_keep_their_evidence(self):
        root = Path(__file__).resolve().parents[1]
        graph = json.loads((root / "lib/semantic-atlas.json").read_text())
        foundation = json.loads((root / "lib/foundation-index.json").read_text())
        findings = {finding["id"] for finding in graph["findings"] if finding["sourceRefIds"]}
        resolved = [task for task in graph["tasks"] if task["status"] == "resolved"]
        self.assertTrue(resolved)
        for task in resolved:
            self.assertTrue(task["resolutionIds"], task["id"])
            self.assertTrue(set(task["resolutionIds"]) <= findings, task["id"])
        for module in foundation["modules"]:
            self.assertTrue(module["synopsis"], module["id"])
            self.assertTrue(module["snapshotSourceHref"], module["id"])
            self.assertEqual(module["reviewState"], "registered-source-not-executed")
        self.assertIn("around 2", inventory.module_synopsis("The estimates straddle around 2.")[0])

    def test_reported_wolfram_runs_match_the_ledger(self):
        root = Path(__file__).resolve().parents[2]
        constants = (root / "website/lib/verification-runs.ts").read_text()
        with (root / "verification/status_ledger.csv").open(newline="") as file:
            claims = {row["claim_id"]: row for row in csv.DictReader(file)}
        for constant, claim in (("WOLFRAM_BASE", "GATE.WOLFRAM.01"),
                                ("WOLFRAM_EXTENSION", "GATE.WOLFRAM.02")):
            count = constants.split(f'export const {constant} = "')[1].split('"')[0]
            self.assertIn(count, claims[claim]["canonical_status"], claim)

    def test_withheld_source_keeps_metadata_and_reason_but_cannot_publish_payload(self):
        ref = {
            "id": "source:db19a580b8ce8c1fa806be0b",
            "path": "experiments/tfpt-discovery/ccc_kms_bridge_closure_probe.py",
            "sha256": "a" * 64,
            "lineStart": 1,
            "lineEnd": 2,
            "href": None,
            "freshness": "not-published",
            "publicationState": "withheld-root",
            "historicalSource": False,
        }
        omissions = [{"sourceId": ref["id"], "path": ref["path"],
                      "reason": "outside approved scientific publication roots"}]
        filename = "sources/" + ref["id"].split(":", 1)[1] + ".txt"
        with tempfile.TemporaryDirectory() as directory, patch.object(inventory, "PUBLIC", Path(directory)):
            inventory.validate_source_ref(ref, omissions, [])
            for payload in ({"href": "/research/leak.txt"}, {"remoteHref": "https://example.com/leak"},
                            {"excerptSha256": "b" * 64}, {"text": "unpublished source bytes"}):
                with self.subTest(payload=payload), self.assertRaises(AssertionError):
                    inventory.validate_source_ref({**ref, **payload}, omissions, [])
            with self.assertRaises(AssertionError):
                inventory.validate_source_ref(ref, [], [])
            with self.assertRaises(AssertionError):
                inventory.validate_source_ref(ref, omissions, [{"path": filename, "sha256": "b" * 64}])
            excerpt = Path(directory) / filename
            excerpt.parent.mkdir()
            excerpt.write_text("unpublished source bytes")
            with self.assertRaises(AssertionError):
                inventory.validate_source_ref(ref, omissions, [])


if __name__ == "__main__":
    unittest.main()
