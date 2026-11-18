import json
import tempfile
import unittest
from pathlib import Path
from semantic_explorer import source_inventory

class PlanInventoryTests(unittest.TestCase):
    def test_live_plan_sources_are_readable_current_and_not_research_claims(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory).resolve()
            base=root/"docs/research/fundamental-plan"
            base.mkdir(parents=True)
            (base/"PLAN.md").write_text("# Current plan\nOriginal context.\n")
            progress=base/"progress.json"
            progress.write_text(json.dumps({"revision":71,"next_task":"F1","continuation":"Actual current continuation"}))
            inventory=source_inventory.build(root,set())
            plan=next(n for n in inventory["nodes"]if n["id"]=="plan:fundamental")
            self.assertEqual(plan["attrs"]["revision"],71)
            self.assertFalse(plan["attrs"]["proof_coverage"])
            originals={s["path"]:s for s in inventory["sources"]}
            self.assertEqual(set(originals),{"docs/research/fundamental-plan/PLAN.md","docs/research/fundamental-plan/progress.json"})
            old=originals["docs/research/fundamental-plan/progress.json"]
            progress.write_text(json.dumps({"revision":72,"next_task":"F1","continuation":"Changed actual instruction"}))
            refreshed=source_inventory.build(root,set())
            new=next(s for s in refreshed["sources"]if s["path"]==old["path"])
            self.assertNotEqual(old["sha256"],new["sha256"])
            self.assertNotEqual(old["id"],new["id"])
            self.assertFalse(any(n["type"] in {"claim","research_finding"}for n in inventory["nodes"]))

if __name__=="__main__":
    unittest.main()
