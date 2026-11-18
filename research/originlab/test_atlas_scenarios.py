import unittest

from atlas_scenarios import assess, contexts, open_inputs, scenarios, spanning
from autonomy import HERE, execute


def text(value):
    return {"de": value, "en": value}


def finding(name, status, premises, path="a.py"):
    return {"id": name, "status": status, "premise_ids": premises, "problem": "problem:x",
            "sources": [{"path": path, "sha256": "0"}]}


def task(name, status, premises, depends=()):
    return {"id": name, "status": status, "kind": "construction", "problem": "problem:x",
            "title": text(name), "success_criterion": text("c"), "next_action": text("n"),
            "required_premise_ids": premises, "depends_on": list(depends)}


def strategy(name, outcome, premises):
    return {"id": name, "outcome": outcome, "premise_ids": premises, "problem": "problem:x",
            "changes_required": text("change"), "finding_ids": []}


KB = {
    "revision": "r", "premises": dict.fromkeys(["a", "b", "c", "d", "e"]),
    "topics": {"problem:x": ["topic:origin"]},
    "findings": [finding("f1", "established_scoped", ["a", "b"]),
                 finding("f2", "excluded_scoped", ["c"]),
                 finding("f3", "established_scoped", ["d"], "b.py"),
                 finding("f4", "conditional", ["a", "c"], "c.py")],
    "tasks": [task("t1", "open", ["a", "c"]), task("t2", "open", ["d"], ["t1"]),
              task("t3", "resolved", ["a"]), task("t4", "open", ["e"])],
    "strategies": [strategy("s-dead", "excluded_scoped", ["a"]), strategy("s-reuse", "reuse_result", ["a"])],
    "connections": [],
}


# Unrelated single-premise findings keep the shared premise below the method-rule share.
PAD = [f"p{n}" for n in range(6)]
KB["premises"].update(dict.fromkeys(PAD))
KB["findings"] += [finding(f"pad-{name}", "established_scoped", [name], f"{name}.py") for name in PAD]


class AtlasScenarioTests(unittest.TestCase):
    def test_joint_use_forms_contexts(self):
        universal, label = contexts(KB)
        self.assertEqual(universal, [])
        self.assertEqual(len({label["a"], label["b"], label["c"]}), 1)
        self.assertEqual(len(set(label.values())), 3+len(PAD))

    def test_open_inputs_carry_closed_and_reusable_routes(self):
        universal, label = contexts(KB)
        rows = open_inputs(KB, universal, label)
        self.assertEqual([row["task"] for row in rows], ["t1", "t4", "t2"])
        self.assertEqual(rows[2]["open_prerequisites"], ["t1"])
        self.assertEqual(rows[0]["closed_routes"][0]["strategy"], "s-dead")
        self.assertEqual(rows[0]["reusable_routes"], ["s-reuse"])

    def test_scenarios_count_identifications_and_unlinked_contexts(self):
        result = scenarios(KB)
        self.assertEqual(result["minimum_identifications"], 2)
        self.assertEqual(len(result["recorded_links"]), 1)
        first = result["scenarios"][0]
        self.assertEqual(first["assumptions"], 2)
        self.assertEqual(first["contexts_without_recorded_link"], [["e"]])
        self.assertEqual(first["all_needed_contexts_joined"]["open_tasks_posed_inside"], 3)
        self.assertFalse(result["physical_solution"])

    def test_spanning_reports_groups_without_a_link(self):
        chosen, groups = spanning({"x", "y", "z"}, [{"contexts": ["x", "y"], "evidence": 1}])
        self.assertEqual(len(chosen), 1)
        self.assertEqual(groups, [["x", "y"], ["z"]])

    def test_a_closed_route_needs_a_changed_premise(self):
        base = {"id": "s", "statement": "one source", "adopts": ["a", "c"], "resolves": ["t1"]}
        self.assertEqual(assess(KB, base)["status"], "ADMISSIBLE_FOR_CHECKS")
        repeat = assess(KB, {**base, "routes": ["s-dead"]})
        self.assertEqual(repeat["status"], "REFUTED_BY_RECORDED_EXCLUSION")
        self.assertEqual(assess(KB, {**base, "routes": ["s-dead"], "changed_premises": ["a"]})["status"],
                         "REFUTED_BY_RECORDED_EXCLUSION")
        self.assertEqual(assess(KB, {**base, "routes": ["s-dead"], "changed_premises": ["d"]})["status"],
                         "ADMISSIBLE_FOR_CHECKS")
        self.assertEqual(assess(KB, {**base, "resolves": ["t2"]})["status"], "INCOMPLETE_PREMISES")
        self.assertEqual(assess(KB, {**base, "adopts": ["zz"]})["status"], "REJECTED_UNKNOWN_RECORD")
        with self.assertRaises(ValueError):
            assess(KB, {**base, "physical_solution": True})

    def test_live_atlas_and_a_supplied_scenario_with_an_executed_check(self):
        root = HERE.parents[1]
        live = execute({"lane": "atlas_scenarios"}, root)
        self.assertEqual(live["status"], "ATLAS_SCENARIOS_ENUMERATED")
        self.assertEqual(live["minimum_identifications"], len(live["needed_contexts"])-1)
        self.assertEqual(len(live["open_inputs"]),
                         live["scenarios"][0]["all_needed_contexts_joined"]["open_tasks"])
        row = next(row for row in live["open_inputs"] if row["contexts"])
        scenario = {"id": "control", "statement": "control", "adopts": row["required_premises"],
                    "resolves": [row["task"]]}
        check = {"job": {"lane": "custom_complex", "simplices": [[0, 1], [1, 2], [0, 2]]}}
        observed = execute(check["job"], root)["status"]
        passed = execute({"lane": "atlas_scenario_case",
                          "scenario": {**scenario, "checks": [{**check, "expect": observed}]}}, root)
        self.assertEqual(passed["status"], "DECLARED_PREDICTIONS_HOLD")
        failed = execute({"lane": "atlas_scenario_case",
                          "scenario": {**scenario, "checks": [{**check, "expect": "OTHER"}]}}, root)
        self.assertEqual(failed["status"], "DECLARED_PREDICTION_FAILED")
        self.assertFalse(passed["physical_solution"])
        replay = execute({"lane": "atlas_scenario_case", "scenario": {**scenario, "checks": [
            {"job": "rewrite", "expect": "EXACT_TYPED_IDENTITIES"}]}}, root)
        self.assertEqual(replay["status"], "ADMISSIBLE_FOR_CHECKS")
        broken = execute({"lane": "atlas_scenario_case", "scenario": {**scenario, "checks": [
            {"job": {"lane": "no_such_job"}, "expect": "X"}]}}, root)
        self.assertEqual(broken["status"], "REJECTED_MALFORMED_CHECK")
        invented = execute({"lane": "atlas_scenario_case", "scenario": {**scenario, "checks": [
            {"job": {"lane": "r113_orientation", "late": True, "time": "dtn"}, "expect": "X"}]}}, root)
        self.assertEqual(invented["status"], "REJECTED_MALFORMED_CHECK")
        self.assertIn("r113_orientation", invented["error"])

    def test_model_answer_is_parsed_as_untrusted_data(self):
        from scenario_model import parse
        text = ('Here: [{"id":"join-a","statement":"s","adopts":["a"],'
                '"new_premises":[{"id":"hypothesis:rule-one","statement":"r","rule":"f(x)=x",'
                '"predicts":["A1","A5"]}]},'
                '{"id":"one-anchor","statement":"s","adopts":["a"],"new_premises":[{"id":"hypothesis:rule-two",'
                '"statement":"r","rule":"f","predicts":["A1"]}]},'
                '{"id":"BAD ID","statement":"s","adopts":["a"]},'
                '{"id":"extra","statement":"s","adopts":["a"],"physical_solution":true}]')
        accepted, rejected = parse(text)
        self.assertEqual([row["id"] for row in accepted], ["model-join-a"])
        self.assertEqual(len(rejected), 3)
        with self.assertRaises(ValueError):
            parse("no array")

    def test_proposal_uses_the_live_record_and_reports_an_absent_model(self):
        from scenario_model import propose
        root = HERE.parents[1]
        seen = {}

        def fake(model, prompt, timeout, budget):
            seen["prompt"] = prompt
            return '[{"id":"one","statement":"s","adopts":[]}]', 0.0

        result = propose(root, "m", 5, 1.0, [{"id": "old", "status": "REFUTED_BY_CHECK"}],
                         {"custom_complex": {"required": ["simplices"], "optional": []}}, {}, ask=fake)
        self.assertEqual(result["status"], "PROPOSED")
        self.assertEqual(result["scenarios"][0]["id"], "model-one")
        self.assertIn("earlier_scenarios", seen["prompt"])
        self.assertLess(result["prompt_chars"], 120000)

        def absent(*_):
            raise FileNotFoundError("claude")

        self.assertEqual(propose(root, "m", 5, 1.0, [], {}, {}, ask=absent)["status"], "MODEL_UNAVAILABLE")

    def test_a_model_written_check_runs_isolated_and_reports_its_verdict(self):
        import os
        from pathlib import Path
        from scenario_agent import investigate, run_script
        root = HERE.parents[1]
        probe = (
            "import json, os, urllib.request\n"
            "def tried(action):\n"
            "    try:\n"
            "        action(); return 'allowed'\n"
            "    except Exception:\n"
            "        return 'denied'\n"
            "net = tried(lambda: urllib.request.urlopen('https://example.com', timeout=3))\n"
            "home = tried(lambda: os.listdir(%r))\n"
            "write = tried(lambda: open(%r, 'w'))\n"
            "print(json.dumps({'verdict': 'PREDICTION_FAILS', 'quantity': 'rank', 'derived': net+home, "
            "'target': write, 'reason': 'probe'}))\n") % (str(Path.home()), str(root/"originlab-probe"))
        result = run_script(root, probe)
        self.assertEqual(result["status"], "REFUTED_BY_COMPUTATION")
        self.assertEqual(result["computed"]["derived"], "denieddenied")
        self.assertEqual(result["computed"]["target"], "denied")
        reads = ("import json, os, subprocess, sys\n"
                 "repo = os.environ['ORIGINLAB_REPO']\n"
                 "head = open(os.path.join(repo, 'AGENTS.md')).read(1)\n"
                 "child = subprocess.run([sys.executable, '-c', 'print(6*7)'], capture_output=True, text=True)\n"
                 "print(json.dumps({'verdict': 'PREDICTION_MATCHES', 'quantity': 'q', 'derived': "
                 "child.stdout.strip(), 'target': '42', 'reason': str(len(head))}))\n")
        ran = run_script(root, reads)
        self.assertEqual(ran["status"], "SURVIVES_COMPUTED_PREDICTION")
        self.assertEqual(ran["computed"]["derived"], "42")
        self.assertFalse((root/"originlab-probe").exists())
        self.assertEqual(run_script(root, "raise SystemExit(2)")["status"], "CHECK_FAILED_TO_RUN")
        self.assertEqual(run_script(root, "print('x')")["status"], "CHECK_OUTPUT_MALFORMED")
        live = execute({"lane": "atlas_scenarios"}, root)
        link = live["recorded_links"][0]["contexts"]
        scenario = {"id": "control", "statement": "control", "adopts": link, "identifies": [link]}
        seen = {}

        def fake(root, model, prompt, timeout, budget):
            seen["prompt"] = prompt
            script = ("import json\nprint(json.dumps({'verdict':'PREDICTION_MATCHES','quantity':'q',"
                      "'derived':1,'target':1,'reason':'r'}))")
            header = json.dumps({"rationale": "r {braces}", "sources_used": [], "target": {"value": 1}})
            return "The spectrum {1, (2/3)^6} is tested.\n"+header+"\n```python\n"+script+"\n```\n", 0.0

        import json
        grade = lambda name: (lambda *_: {"grade": name, "reason": "r", "cost_usd": 0.0})
        verdict = investigate(root, scenario, "m", 5, 1.0, ask=fake, referee=grade("INDEPENDENT_AND_SPECIFIC"))
        self.assertEqual(verdict["status"], "SURVIVES_COMPUTED_PREDICTION")
        weak = investigate(root, scenario, "m", 5, 1.0, ask=fake, referee=grade("FITTED_TO_TARGET"))
        self.assertEqual(weak["status"], "TEST_NOT_DISCRIMINATING")
        self.assertEqual(weak["computed_status"], "SURVIVES_COMPUTED_PREDICTION")
        rule = {"id": "hypothesis:new-rule", "statement": "a rule the record lacks",
                "rule": "w(x) = x", "predicts": ["A1", "A2"]}
        fresh = {"id": "new", "statement": "s", "adopts": link[:1], "new_premises": [rule],
                 "identifies": [["hypothesis:new-rule", link[0]]]}
        self.assertEqual(execute({"lane": "atlas_scenario_case", "scenario": fresh}, root)
                         ["new_hypotheses"]["hypothesis:new-rule"]["predicts"], ["A1", "A2"])
        single = investigate(root, fresh, "m", 5, 1.0, ask=fake, referee=grade("INDEPENDENT_AND_SPECIFIC"))
        self.assertEqual(single["status"], "TOO_FEW_INDEPENDENT_TARGETS")

        def two(root, model, prompt, timeout, budget):
            script = ("import json\nprint(json.dumps({'verdict':'PREDICTION_MATCHES','quantity':'q',"
                      "'derived':1,'target':1,'targets':[{'group':'A1','derived':1,'target':1},"
                      "{'group':'A2','derived':2,'target':2}],'reason':'r'}))")
            return "{\"rationale\": \"r\"}\n```python\n"+script+"\n```\n", 0.0

        both = investigate(root, fresh, "m", 5, 1.0, ask=two, referee=grade("INDEPENDENT_AND_SPECIFIC"))
        self.assertEqual(both["status"], "SURVIVES_COMPUTED_PREDICTION")
        self.assertEqual(investigate(root, fresh, "m", 5, 1.0, ask=two, referee=grade("RULE_NOT_USED"))
                         ["status"], "TEST_NOT_DISCRIMINATING")
        self.assertIn("a rule the record lacks", seen["prompt"])
        self.assertIn("The self-reference", seen["prompt"])
        self.assertEqual(verdict["input_role"], "unreviewed_model_written_check")
        self.assertFalse(verdict["physical_solution"])
        self.assertIn(link[0], seen["prompt"])
    def test_a_research_thread_revises_until_its_rule_survives(self):
        import json
        from scenario_agent import ANCHOR_SOURCES, graded, research
        root = HERE.parents[1]
        self.assertTrue(all((root/path).is_file() for paths in ANCHOR_SOURCES.values() for path in paths))
        self.assertEqual(graded("REFUTED_BY_COMPUTATION", "FITTED_TO_TARGET"), "HYPOTHESIS_REJECTED")
        self.assertEqual(graded("SURVIVES_COMPUTED_PREDICTION", "FITTED_TO_TARGET"), "TEST_NOT_DISCRIMINATING")
        self.assertEqual(graded("REFUTED_BY_COMPUTATION", "INDEPENDENT_AND_SPECIFIC"), "REFUTED_BY_COMPUTATION")
        prompts = []

        def answer(verdict, groups):
            targets = [{"group": g, "derived": 1, "target": 1} for g in groups]
            script = "import json\nprint(json.dumps(%r))" % {"verdict": verdict, "quantity": "q",
                                                              "targets": targets, "reason": "r"}
            script = script.replace("'", '"')
            header = json.dumps({"id": "x", "statement": "s", "rule": "f", "rationale": "r"})
            return header+"\n```python\n"+script+"\n```\n"

        def ask(root, model, prompt, timeout, budget):
            prompts.append(prompt)
            return answer("PREDICTION_FAILS" if len(prompts) == 1 else "PREDICTION_MATCHES", ["A1", "A2"]), 0.5

        grades = iter(["FITTED_TO_TARGET", "INDEPENDENT_AND_SPECIFIC"])
        referee = lambda *_: {"grade": next(grades), "reason": "r", "cost_usd": 0.1}
        result = research(root, ["A1", "A2"], "m", 5, 1.0, ask=ask, referee=referee)
        self.assertEqual([row["status"] for row in result["attempts"]],
                         ["HYPOTHESIS_REJECTED", "SURVIVES_COMPUTED_PREDICTION"])
        self.assertEqual(result["id"], "research-a1-a2")
        self.assertAlmostEqual(result["cost_usd"], 1.2)
        self.assertIn("previous_rule", prompts[1])
        self.assertIn("e8_family_return.py", prompts[0])
        self.assertFalse(result["physical_solution"])
        from scenario_agent import run_with_repairs
        fixed = ("import json\nprint(json.dumps({'verdict':'PREDICTION_FAILS','quantity':'q','reason':'r'}))")
        outcome, script, cost = run_with_repairs(root, "assert 1 == 2", "m", 1.0,
                                                 ask=lambda *a: ("```python\n"+fixed+"\n```", 0.2))
        self.assertEqual(outcome["status"], "REFUTED_BY_COMPUTATION")
        self.assertTrue(outcome["repaired"])
        self.assertEqual((script.strip(), cost), (fixed, 0.2))
        none = research(root, ["A1", "A2"], "m", 5, 1.0,
                        ask=lambda *a: (answer("NO_TESTABLE_CONSEQUENCE", []), 0.0), referee=referee)
        self.assertEqual(len(none["attempts"]), 1)
        self.assertEqual(none["status"], "NO_TESTABLE_CONSEQUENCE_FOUND")

    def test_every_anchor_group_is_in_the_brief(self):
        from scenario_model import ANCHORS, THEORY
        text = THEORY.read_text(encoding="utf-8")
        self.assertTrue(all(f"- {name}," in text for name in ANCHORS))


if __name__ == "__main__":
    unittest.main()
