"""CPU regressions for persistence, genetic composition, dedup and diversity."""
import contextlib
import io
import json
from pathlib import Path
import sqlite3
import tempfile
import unittest

import numpy as np

import gpu_search as engine


FIXTURE = '''
def sample(rng,n):
    return [[int(rng.integers(-1000,1001)),int(rng.integers(-1000,1001))] for _ in range(n)]
def canonicalize(g): return [int(x) for x in g]
def valid(g): return len(g)==2 and all(abs(x)<=2000 for x in g)
def mutate(rng,g):
    g[int(rng.integers(2))]+=int(rng.choice([-7,-3,-1,1,3,7]))
    return g
def crossover(rng,a,b): return [a[0],b[1]]
def evaluate(gs,device):
    assert device=='cpu'
    return [dict(objectives=[float((g[0]-7)**2+(g[1]+3)**2),float(abs(g[0])+abs(g[1]))],
        feasible=True,descriptor=[g[0]%4,g[1]%4],details=dict(construction=g)) for g in gs]
'''


class SearchTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.plugin_path = self.root / "fixture.py"
        self.plugin_path.write_text(FIXTURE)

    def tearDown(self):
        self.tmp.cleanup()

    def args(self, name="run", generations=3, island=0):
        return engine.parser().parse_args([
            "--problem", str(self.plugin_path), "--output", str(self.root / name),
            "--device", "cpu", "--population", "12", "--batch", "12",
            "--migration-every", "2", "--seed", str(31 + island),
            "--island", str(island), "--generations", str(generations)])

    def run_search(self, args):
        with contextlib.redirect_stdout(io.StringIO()):
            return engine.run(args)

    def state(self, name, island=0):
        return json.loads((self.root / name / f"island-{island}.json").read_text())

    def test_resume_matches_uninterrupted(self):
        self.run_search(self.args("full", 6))
        self.run_search(self.args("resumed", 2))
        result = self.run_search(self.args("resumed", 6))
        self.assertEqual(result["total_evaluated_here"], 72)
        full, resumed = self.state("full"), self.state("resumed")
        for field in ("population", "rng", "python_rng", "generation"):
            self.assertEqual(full[field], resumed[field])

    def test_source_change_rejects_resume(self):
        self.run_search(self.args())
        self.plugin_path.write_text(FIXTURE + "\n# objective changed\n")
        with self.assertRaisesRegex(ValueError, "source or search configuration"):
            self.run_search(self.args(generations=4))

    def test_stop_sentinel_preserves_then_resumes(self):
        args = self.args()
        directory = Path(args.output)
        directory.mkdir()
        (directory / "STOP").touch()
        self.assertEqual(self.run_search(args)["reason"], "stop_sentinel")
        self.assertEqual(self.state("run")["total_evaluated_here"], 0)
        (directory / "STOP").unlink()
        self.assertEqual(self.run_search(args)["total_evaluated_here"], 36)

    def test_exhaustion_does_not_loop_or_repeat_evaluations(self):
        self.plugin_path.write_text(FIXTURE + '''
def sample(rng,n): return [[0,0] for _ in range(n)]
def mutate(rng,g): return g
def crossover(rng,a,b): return a
''')
        result = self.run_search(self.args(generations=100))
        self.assertEqual(result["reason"], "no_novel_candidates")
        self.assertEqual(result["generation"], 1)
        self.assertEqual(result["total_evaluated_here"], 1)

    def test_shared_islands_deduplicate_and_migrate(self):
        self.run_search(self.args("shared", 2, island=0))
        self.run_search(self.args("shared", 2, island=1))
        db = sqlite3.connect(self.root / "shared" / "archive.sqlite3")
        count, distinct = db.execute("SELECT count(*),count(DISTINCT key) FROM candidates").fetchone()
        self.assertEqual(count, distinct)
        self.assertEqual(count, 48)
        self.assertEqual(db.execute("SELECT count(DISTINCT owner) FROM candidates").fetchone()[0], 2)
        owner0 = {r[0] for r in db.execute("SELECT key FROM candidates WHERE owner='0'")}
        self.assertTrue(owner0 & {r["key"] for r in self.state("shared", 1)["population"]})
        db.close()

    def test_selection_retains_feasible_diverse_constructions(self):
        rows = [{"key": str(i), "genotype": [i], "result": dict(
            objectives=[float(i), float(30-i)], feasible=i != 0,
            descriptor=[i % 5])} for i in range(30)]
        chosen = engine.select(rows, 10, np.random.default_rng(8))
        self.assertTrue(all(r["result"]["feasible"] for r in chosen))
        self.assertGreaterEqual(len({tuple(r["result"]["descriptor"]) for r in chosen}), 4)

    def test_reject_invalid_objective(self):
        with self.assertRaisesRegex(ValueError, "Nonfinite"):
            engine.validate_result(dict(objectives=[float("nan")], feasible=True, descriptor=[]))
        with self.assertRaisesRegex(ValueError, "explicit boolean"):
            engine.validate_result(dict(objectives=[1.0], feasible="yes", descriptor=[]))

    def test_blocked_pareto_matches_scalar_definition(self):
        rng = np.random.default_rng(73)
        scores = rng.integers(0, 10, size=(513, 4)).astype(float)
        eligible = rng.random(513) > .2
        expected = [i for i in range(len(scores)) if eligible[i] and not any(
            eligible[j] and engine.dominates(scores[j], scores[i]) for j in range(len(scores)))]
        self.assertEqual(engine.pareto_indices(scores, eligible), expected)

    def test_two_exhaustive_shards_complete_deduplicated_union(self):
        self.plugin_path.write_text(FIXTURE + '''
ENUMERATION_SIZE=8
def sample(rng,n): return [[int(rng.integers(8)),0] for _ in range(n)]
def enumerate_candidates(shard,shards):
    for x in range(shard,ENUMERATION_SIZE,shards): yield [x,0]
''')
        for island in (0,1):
            args=self.args("finite",1,island)
            args.population=args.batch=2
            args.exhaustive_after=True
            args.exhaustive_shards=2
            result=self.run_search(args)
            self.assertEqual(result["classification"]["classification_count"],4)
            self.assertEqual(self.state("finite",island)["generation"],1)
            self.assertEqual(self.state("finite",island)["total_evaluated_here"],2)
        union=json.loads((self.root/"finite/classification-union.json").read_text())
        self.assertTrue(union["complete"])
        self.assertEqual(union["classified_unique"],8)
        self.assertEqual(union["feasible_unique"],8)
        self.assertEqual(len(union["survivors"]),8)
        status=json.loads((self.root/"finite/island-1-status.json").read_text())
        self.assertEqual(status["reason"],"classification_complete")
        rerun=self.run_search(args)
        self.assertEqual(rerun["classification"]["classification_new_evaluations"],0)
        repeated=json.loads((self.root/"finite/classification-union.json").read_text())
        self.assertEqual(union["canonical_key_digest"],repeated["canonical_key_digest"])
        self.assertEqual(union["classification_result_digest"],repeated["classification_result_digest"])

    def test_overlapping_shards_are_rejected(self):
        self.plugin_path.write_text(FIXTURE + '''
ENUMERATION_SIZE=8
def enumerate_candidates(shard,shards):
    for x in range(8): yield [x,0]
''')
        for island in (0,1):
            args=self.args("overlap",0,island)
            args.exhaustive_after=True
            args.exhaustive_shards=2
            if island==0:
                self.run_search(args)
            else:
                with self.assertRaisesRegex(ValueError,"shards overlap"):
                    self.run_search(args)


if __name__ == "__main__":
    unittest.main()
