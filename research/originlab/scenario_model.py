"""Ask a language model for solution scenarios over the compiled Atlas record.

The model receives the open inputs, contexts, closed routes and earlier
verdicts. Its answer is untrusted data: a list of scenarios that name Atlas
records and typed jobs. Nothing in the answer is executed as code.
"""

from pathlib import Path
from tempfile import TemporaryDirectory
import json
import re
import subprocess

FIELDS = {"id", "statement", "adopts", "identifies", "resolves", "routes", "changed_premises",
          "new_premises"}
THEORY = Path(__file__).with_name("THEORY-BRIEF.md")
ANCHORS = tuple(f"A{n}" for n in range(1, 12))
MAX_SCENARIOS = 6
RULES = """You propose solution scenarios for an open research programme. Return only a JSON array
of at most 6 objects. Use no tools. Each object has:
  "id": short lowercase slug,
  "statement": one or two sentences stating which recorded contexts describe the same source and why
               that would resolve the named tasks,
  "adopts": premise ids the scenario assumes together,
  "identifies": optional list of [premise id, premise id] pairs read as the same object,
  "resolves": open task ids the scenario addresses; every required premise of such a task must lie
              in a context of an adopted or identified premise,
  "routes": optional strategy ids the scenario uses,
  "changed_premises": premise ids that a closed route did not have; required to reuse a closed route,
  "new_premises": optional list of {"id": "hypothesis:<slug>", "statement": text, "rule": text,
                  "predicts": [anchor group ids]} for a rule the record does not contain.
                  "statement" names the object the rule acts on and what it fixes that the declared
                  inputs leave open. "rule" is the explicit mathematical definition: an operator,
                  map, weight or constraint with its formula, concrete enough to program. "predicts"
                  names at least two anchor groups (A1 to A11) from different constructions whose
                  values this one rule yields. A sentence that restates the open problem, or a rule
                  whose only consequence is one value that an existing script already produces, is
                  rejected.
Use only record ids that appear below, plus your own hypothesis ids. At least two of your scenarios
must introduce a new premise: the reviewed record says the declared inputs leave the four endpoint
marks undetermined, so a missing rule is expected. Read the theory brief first. Prefer scenarios that join contexts without a recorded link, that
address several open tasks through one shared input, and that differ from the earlier scenarios.
Study the earlier scenarios: each lists its rule and why its test ended as it did. Do not repeat
a rule that failed for vagueness or inconsistency; state a sharper one or a different one. A rule
must define every symbol it uses, act on a concrete space of stated dimension, and its formula must
give the claimed values when evaluated literally.
Each admissible scenario with an "identifies" pair or a new premise is tested by a computation:
a consequence of the scenario is derived and compared with an exact recorded value that was not
used as input. Propose only what yields such a consequence. A pair may join a hypothesis id with a
record premise.
State no physical derivation as established."""


def clip(text, size):
    return text if len(text) <= size else text[:size-1]+"…"


def digest(row):
    """What the proposer needs from an earlier verdict: the rule and why its test ended as it did."""
    test = row.get("computed_comparison") or {}
    return {"id": row["id"], "statement": clip(row.get("statement", ""), 300), "status": row.get("status"),
            "rules": [{"id": new["id"], "rule": clip(new.get("rule", ""), 400), "predicts": new.get("predicts")}
                      for new in row.get("new_premises", [])],
            "test": {"status": test.get("status"), "reason": clip(str((test.get("computed") or {}).get("reason")), 500),
                     "referee": test.get("referee")} if test else None}


def brief(kb, atlas, history, catalogue, known):
    """Bounded text of the record a proposal needs."""
    title = lambda p: clip(kb["premises"][p]["title"]["de"], 110)
    needed = {p for row in atlas["open_inputs"] for p in row["required_premises"]}
    closed, tasks = {}, []
    for row in atlas["open_inputs"]:
        for route in row["closed_routes"][:6]:
            closed[route["strategy"]] = clip(route["changes_required"], 220)
        tasks.append({"task": row["task"], "title": clip(row["title"], 120),
                      "criterion": clip(row["success_criterion"], 260),
                      "requires": row["required_premises"], "contexts": row["contexts"],
                      "open_prerequisites": row["open_prerequisites"],
                      "closed_routes": [r["strategy"] for r in row["closed_routes"][:6]],
                      "reusable_routes": row["reusable_routes"][:6]})
    first = atlas["scenarios"][0] if atlas["scenarios"] else {}
    data = {
        "record": atlas["record"], "minimum_identifications": atlas["minimum_identifications"],
        "needed_contexts": atlas["needed_contexts"],
        "contexts_without_recorded_link": first.get("contexts_without_recorded_link", []),
        "recorded_links": atlas["recorded_links"][:20],
        "premises": {p: title(p) for p in sorted(needed | {c["context"] for c in atlas["needed_contexts"]})},
        "open_tasks": tasks, "closed_routes_and_required_change": closed,
        "earlier_scenarios": [digest(row) for row in history[-60:]],
    }
    return (RULES+"\n\n# Theory brief\n\n"+THEORY.read_text(encoding="utf-8")+"\n\n# Record\n\n"
            +json.dumps(data, ensure_ascii=False, separators=(",", ":")))


def parse(text):
    match = re.search(r"\[.*\]", text, re.S)
    if match is None:
        raise ValueError("the model returned no JSON array")
    rows = json.loads(match[0])
    if not isinstance(rows, list):
        raise ValueError("the model must return a JSON array")
    accepted, rejected = [], []
    for row in rows[:MAX_SCENARIOS]:
        ok = (isinstance(row, dict) and {"id", "statement", "adopts"} <= set(row) <= FIELDS
              and isinstance(row["id"], str) and re.fullmatch(r"[a-z0-9][a-z0-9-]{2,60}", row["id"])
              and isinstance(row["statement"], str) and len(row["statement"]) <= 1200
              and all(isinstance(row.get(key, []), list)
                      for key in ("adopts", "identifies", "resolves", "routes", "changed_premises",
                                  "new_premises"))
              and all(isinstance(new, dict) and set(new) == {"id", "statement", "rule", "predicts"}
                      and isinstance(new["id"], str) and re.fullmatch(r"hypothesis:[a-z0-9-]{3,60}", new["id"])
                      and all(isinstance(new[key], str) and 0 < len(new[key]) <= 2000
                              for key in ("statement", "rule"))
                      and isinstance(new["predicts"], list) and len(set(new["predicts"])) >= 2
                      and set(new["predicts"]) <= set(ANCHORS)
                      for new in row.get("new_premises", [])))
        (accepted if ok else rejected).append(row)
    # The prefix marks the origin in every later report.
    return [{**row, "id": "model-"+row["id"]} for row in accepted], rejected


def ask(model, prompt, timeout, budget):
    """Run the model without tools, settings or session storage in an empty directory."""
    with TemporaryDirectory(prefix="originlab-scenarios-") as folder:
        command = ["claude", "-p", "--tools", "", "--strict-mcp-config", "--setting-sources", "",
                   "--no-session-persistence", "--output-format", "json", "--model", model,
                   "--max-budget-usd", str(budget)]
        run = subprocess.run(command, input=prompt, text=True, capture_output=True, cwd=folder,
                             timeout=timeout, check=True)
    answer = json.loads(run.stdout)
    if answer.get("is_error"):
        raise ValueError(f"model call failed: {answer.get('subtype')}")
    return answer["result"], answer.get("total_cost_usd")


def propose(root, model, timeout, budget, history, catalogue, known, ask=ask):
    from atlas_scenarios import PROGRAM, knowledge, scenarios
    if not (root/PROGRAM).is_file():
        return {"status": "ATLAS_UNAVAILABLE", "scenarios": [], "physical_solution": False}
    kb = knowledge(root)
    prompt = brief(kb, scenarios(kb), history, catalogue, known)
    try:
        text, cost = ask(model, prompt, timeout, budget)
        accepted, rejected = parse(text)
    except (OSError, ValueError, KeyError, subprocess.SubprocessError) as error:
        return {"status": "MODEL_UNAVAILABLE", "model": model, "error_type": type(error).__name__,
                "error": clip(str(error), 400), "scenarios": [], "physical_solution": False}
    return {"status": "PROPOSED", "model": model, "atlas_revision": kb["revision"],
            "prompt_chars": len(prompt), "cost_usd": cost, "scenarios": accepted,
            "malformed": len(rejected), "input_role": "untrusted_model_proposal",
            "physical_solution": False}
