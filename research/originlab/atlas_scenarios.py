"""Compile the reviewed Atlas program into contexts, open inputs and solution scenarios.

Every reviewed finding holds under its recorded premises. Premises used
together in one finding belong to one context. A scenario states which
contexts describe the same source and reports what follows from the record.
"""

from collections import Counter, defaultdict
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import argparse
import json

PROGRAM = Path("semantic_explorer/research/program.json")
DEAD = ("excluded_scoped", "insufficient_for_goal")
EXECUTABLE = (".py", ".lean")
SCOPE = ("bookkeeping over reviewed Atlas records; a shared context or source file is recorded "
         "evidence for a link and proves no identification")


def knowledge(root):
    raw = (root/PROGRAM).read_bytes()
    data = json.loads(raw)
    rows = {kind: [] for kind in ("findings", "tasks", "strategies", "connections")}
    for problem in data["problems"]:
        for kind, target in rows.items():
            target.extend({**row, "problem": problem["id"]} for row in problem.get(kind, []))
    return {"revision": data["revision"], "sha256": sha256(raw).hexdigest(),
            "premises": {row["id"]: row for row in data["premises"]},
            "topics": {problem["id"]: problem.get("topics", []) for problem in data["problems"]},
            **rows}


def contexts(kb, universal_share=0.3):
    """Label every premise with its connected component of joint use."""
    adjacent = defaultdict(set)
    nodes = set()
    for finding in kb["findings"]:
        premises = set(finding["premise_ids"])
        nodes |= premises
        for left, right in combinations(premises, 2):
            adjacent[left].add(right)
            adjacent[right].add(left)
    for row in (*kb["tasks"], *kb["strategies"]):
        nodes |= set(row.get("required_premise_ids", row.get("premise_ids", [])))
    # A premise attached to a large share of all others is a method rule.
    universal = {node for node in nodes if len(adjacent[node]) > universal_share*len(nodes)}
    label = {}
    for start in sorted(nodes-universal):
        if start in label:
            continue
        stack, members = [start], set()
        while stack:
            node = stack.pop()
            if node not in members:
                members.add(node)
                stack.extend(adjacent[node]-universal-members)
        label.update(dict.fromkeys(members, start))
    return sorted(universal), label


def context_records(kb, universal, label):
    """Findings, problems and source files recorded inside each context."""
    record = defaultdict(lambda: {"findings": Counter(), "problems": set(), "sources": set(), "premises": 0})
    for premise, context in label.items():
        record[context]["premises"] += 1
    for finding in kb["findings"]:
        for context in {label[p] for p in finding["premise_ids"] if p not in universal}:
            record[context]["findings"][finding["status"]] += 1
            record[context]["problems"].add(finding["problem"])
            record[context]["sources"].update(source["path"] for source in finding["sources"])
    return record


def open_inputs(kb, universal, label):
    """Open tasks without an open prerequisite, with the routes already decided."""
    open_ids = {task["id"] for task in kb["tasks"] if task["status"] == "open"}
    routes = defaultdict(list)
    for strategy in kb["strategies"]:
        routes[strategy["problem"]].append(strategy)
    shared = defaultdict(set)
    for link in kb["connections"]:
        if link["kind"] == "shares_unresolved_input":
            shared[link["problem"]].add(link["target_problem_id"])
            shared[link["target_problem_id"]].add(link["problem"])
    rows = []
    for task in kb["tasks"]:
        if task["id"] not in open_ids:
            continue
        required = [p for p in task["required_premise_ids"] if p not in universal]
        rows.append({
            "task": task["id"], "problem": task["problem"], "kind": task["kind"],
            "title": task["title"]["de"], "success_criterion": task["success_criterion"]["de"],
            "next_action": task["next_action"]["de"], "required_premises": task["required_premise_ids"],
            "contexts": sorted({label[p] for p in required}),
            "open_prerequisites": sorted(set(task["depends_on"]) & open_ids),
            "closed_routes": [{"strategy": s["id"], "outcome": s["outcome"],
                               "changes_required": s["changes_required"]["de"]}
                              for s in routes[task["problem"]] if s["outcome"] in DEAD],
            "reusable_routes": [s["id"] for s in routes[task["problem"]] if s["outcome"] == "reuse_result"],
            "shares_unresolved_input_with": sorted(shared[task["problem"]]),
        })
    return sorted(rows, key=lambda row: (bool(row["open_prerequisites"]), row["task"]))


def links(kb, record, needed):
    """Recorded evidence that two needed contexts concern the same object."""
    connected = Counter()
    for link in kb["connections"]:
        connected[frozenset((link["problem"], link["target_problem_id"]))] += 1
    rows = []
    for left, right in combinations(sorted(needed), 2):
        sources = record[left]["sources"] & record[right]["sources"]
        problems = record[left]["problems"] & record[right]["problems"]
        crossing = sum(connected[frozenset((a, b))] for a in record[left]["problems"]
                       for b in record[right]["problems"] if a != b)
        if sources or problems or crossing:
            rows.append({"contexts": [left, right], "shared_sources": sorted(sources),
                         "shared_problems": sorted(problems), "problem_connections": crossing,
                         "evidence": len(sources)+len(problems)+crossing})
    return sorted(rows, key=lambda row: (-row["evidence"], row["contexts"]))


def spanning(needed, candidates, skip=()):
    """Join the needed contexts along the strongest recorded links."""
    parent = {context: context for context in needed}

    def find(node):
        while parent[node] != node:
            parent[node] = parent[parent[node]]
            node = parent[node]
        return node

    chosen = []
    for row in candidates:
        if tuple(row["contexts"]) in skip:
            continue
        left, right = map(find, row["contexts"])
        if left != right:
            parent[left] = right
            chosen.append(row)
    groups = defaultdict(list)
    for context in sorted(needed):
        groups[find(context)].append(context)
    return chosen, sorted(groups.values(), key=lambda group: (-len(group), group))


def consequences(kb, universal, label, environment):
    """What the record yields once the listed contexts are read as one source."""
    inside = lambda premises: all(p in universal or label[p] in environment for p in premises)
    status = Counter(f["status"] for f in kb["findings"] if inside(f["premise_ids"]))
    topics = Counter(topic for f in kb["findings"] if inside(f["premise_ids"])
                     for topic in kb["topics"][f["problem"]])
    open_tasks = [task for task in kb["tasks"] if task["status"] == "open"]
    return {"findings": dict(sorted(status.items())), "findings_by_topic": dict(sorted(topics.items())),
            "closed_routes_inside": sum(1 for s in kb["strategies"]
                                        if s["outcome"] in DEAD and inside(s["premise_ids"])),
            "open_tasks_posed_inside": sum(1 for task in open_tasks if inside(task["required_premise_ids"])),
            "open_tasks": len(open_tasks)}


def scenarios(kb, limit=8):
    """Enumerate the ways to join all contexts the open tasks need."""
    universal, label = contexts(kb)
    record = context_records(kb, universal, label)
    inputs = open_inputs(kb, universal, label)
    needed = {context for row in inputs for context in row["contexts"]}
    candidates = links(kb, record, needed)
    base, _ = spanning(needed, candidates)
    # Each alternative drops one chosen link and lets the next best one reconnect.
    variants = [()] + [(tuple(row["contexts"]),) for row in base]
    seen, rows = set(), []
    for skip in variants:
        chosen, groups = spanning(needed, candidates, skip)
        key = tuple(sorted(tuple(row["contexts"]) for row in chosen))
        if key in seen:
            continue
        seen.add(key)
        main = set(groups[0])
        rows.append({
            "identifications": [{"contexts": row["contexts"], "evidence": row["evidence"],
                                 "replayable_sources": [p for p in row["shared_sources"]
                                                        if p.endswith(EXECUTABLE)]} for row in chosen],
            "contexts_without_recorded_link": [group for group in groups[1:]],
            "assumptions": len(chosen)+len(groups)-1,
            "largest_joined_environment": consequences(kb, universal, label, main),
            "all_needed_contexts_joined": consequences(kb, universal, label, needed),
            "status": "COMPLETE_UNDER_ASSUMED_IDENTIFICATIONS",
        })
    rows.sort(key=lambda row: (len(row["contexts_without_recorded_link"]),
                               -sum(item["evidence"] for item in row["identifications"])))
    sizes = Counter(label.values())
    return {
        "status": "ATLAS_SCENARIOS_ENUMERATED", "atlas_revision": kb["revision"], "scope": SCOPE,
        "physical_solution": False,
        "record": {"premises": len(kb["premises"]), "findings": len(kb["findings"]),
                   "contexts": len(sizes), "largest_context_premises": max(sizes.values()),
                   "method_premises": universal},
        "open_inputs": inputs,
        "needed_contexts": [{"context": context, "premises": record[context]["premises"],
                             "findings": dict(record[context]["findings"]),
                             "open_tasks": sum(context in row["contexts"] for row in inputs)}
                            for context in sorted(needed, key=lambda c: -record[c]["premises"])],
        "recorded_links": [{k: row[k] for k in ("contexts", "evidence", "problem_connections")}
                           | {"shared_sources": len(row["shared_sources"])} for row in candidates],
        "minimum_identifications": len(needed)-1,
        "scenarios": rows[:limit],
    }


def assess(kb, proposal):
    """Check a supplied scenario against the record before any job runs."""
    universal, label = contexts(kb)
    allowed = {"id", "statement", "adopts", "identifies", "resolves", "routes", "changed_premises",
               "checks", "new_premises"}
    if not isinstance(proposal, dict) or not {"id", "statement", "adopts"} <= set(proposal) <= allowed:
        raise ValueError("a scenario needs id, statement and adopts, and only registered fields")
    # A rule absent from the record enters as a declared hypothesis with its own context.
    hypotheses = {row["id"]: {key: row.get(key) for key in ("statement", "rule", "predicts")}
                  for row in proposal.get("new_premises", [])}
    if any(not name.startswith("hypothesis:") or name in kb["premises"] for name in hypotheses):
        raise ValueError("a new premise needs an unused id that starts with hypothesis:")
    label = {**label, **{name: name for name in hypotheses}}
    premises = [*proposal["adopts"], *(p for pair in proposal.get("identifies", []) for p in pair),
                *proposal.get("changed_premises", []), *hypotheses]
    unknown = sorted(set(premises)-set(kb["premises"])-set(hypotheses))
    tasks = {task["id"]: task for task in kb["tasks"]}
    strategies = {strategy["id"]: strategy for strategy in kb["strategies"]}
    unknown += sorted(set(proposal.get("resolves", []))-set(tasks))
    unknown += sorted(set(proposal.get("routes", []))-set(strategies))
    if unknown:
        return {"id": proposal["id"], "status": "REJECTED_UNKNOWN_RECORD", "unknown": unknown}
    environment = {label[p] for p in premises if p not in universal}
    changed = set(proposal.get("changed_premises", [])) | set(hypotheses)
    # A closed route may be reused only with a premise its exclusion did not have.
    repeated = [{"strategy": name, "outcome": strategies[name]["outcome"],
                 "changes_required": strategies[name]["changes_required"]["de"]}
                for name in proposal.get("routes", [])
                if strategies[name]["outcome"] in DEAD and not changed-set(strategies[name]["premise_ids"])]
    unmet = {task: sorted(p for p in tasks[task]["required_premise_ids"]
                          if p not in universal and label[p] not in environment)
             for task in proposal.get("resolves", [])}
    unmet = {task: missing for task, missing in unmet.items() if missing}
    open_tasks = {task["id"] for task in kb["tasks"] if task["status"] == "open"}
    status = ("REFUTED_BY_RECORDED_EXCLUSION" if repeated else
              "INCOMPLETE_PREMISES" if unmet else "ADMISSIBLE_FOR_CHECKS")
    return {"id": proposal["id"], "status": status, "repeated_closed_routes": repeated,
            "tasks_with_missing_premises": unmet,
            "open_tasks_addressed": len(set(proposal.get("resolves", [])) & open_tasks),
            "open_tasks": len(open_tasks), "new_hypotheses": hypotheses,
            "consequences": consequences(kb, universal, label, environment-set(hypotheses)),
            "input_role": "candidate_requires_source_review", "physical_solution": False}


def evidence_integrity(root, kb):
    """Compare the recorded source digests of all findings with the working tree."""
    digests, state = {}, Counter()
    stale = set()
    for finding in kb["findings"]:
        for source in finding["sources"]:
            path = source["path"]
            if path not in digests:
                file = root/path
                digests[path] = sha256(file.read_bytes()).hexdigest() if file.is_file() else None
            kind = ("missing" if digests[path] is None else
                    "current" if digests[path] == source.get("sha256") else "changed")
            state[kind] += 1
            if kind != "current":
                stale.add(finding["id"])
    return {"source_references": dict(sorted(state.items())), "distinct_files": len(digests),
            "findings_with_changed_or_missing_source": len(stale)}


def run(root, limit=8):
    kb = knowledge(root)
    return {**scenarios(kb, limit), "evidence_integrity": evidence_integrity(root, kb)}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--limit", type=int, default=8)
    parser.add_argument("--assess", type=Path, help="JSON list of supplied scenarios")
    args = parser.parse_args()
    if args.assess:
        kb = knowledge(args.repo)
        print(json.dumps([assess(kb, row) for row in json.loads(args.assess.read_text())],
                         ensure_ascii=False, indent=1))
    else:
        print(json.dumps(run(args.repo, args.limit), ensure_ascii=False, indent=1))
