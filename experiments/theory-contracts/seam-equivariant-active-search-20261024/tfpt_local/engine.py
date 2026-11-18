"""Small real proposal -> exact witness -> repair path, with resumable files."""
from __future__ import annotations

import json
from pathlib import Path
import random
import re
import time
import urllib.request

from .mathcore import Source, bits, digest, exact_gf2

SYSTEM = (
    "Solve an exact binary character correction problem from an E8 lattice field lift. "
    "Each equation is XOR over the listed variable indices. Variables a0..a7 are binary integers. "
    "Return JSON only: {\"bits\":[a0,a1,a2,a3,a4,a5,a6,a7]}. "
    "The array order is a0 first, NOT a binary numeral. Satisfy ALL equations. "
    "If impossible return {\"status\":\"inconsistent\",\"rows\":[equation_indices]} where XORing "
    "those distinct equations yields 0=1. Do not claim a physics derivation."
)
PROPOSAL_SCHEMA = {
    "type": "object", "properties": {"bits": {"type": "array", "items": {"type": "integer", "enum": [0, 1]},
                                                "minItems": 8, "maxItems": 8}},
    "required": ["bits"], "additionalProperties": False,
}
PROGRAM_SCHEMA = {"anyOf": [PROPOSAL_SCHEMA, {
    "type": "object", "properties": {"operation": {"type": "string", "enum": ["gf2_elimination"]}},
    "required": ["operation"], "additionalProperties": False}]}
TOOL_SYSTEM = SYSTEM + (
    ' You also have the exact symbolic operation gf2_elimination. For a binary linear system, '
    'select that operation rather than doing arithmetic yourself: return {"operation":"gf2_elimination"}. '
    'It acts ONLY on the supplied original equations. You cannot change the task or checker.'
)


def prompt(task):
    # Retain all constraints for the certificate, but avoid redundant prompt tokens.
    unique = {}
    for i, (r, b) in enumerate(zip(task["rows"], task["rhs"])):
        if r or b:
            unique.setdefault((r, b), i)
    equations = [{"equation_index": i, "xor_indices": [j for j in range(8) if r >> j & 1], "equals": b}
                 for (r, b), i in unique.items()]
    return ("Find the character correction. All additions are modulo 2.\n"
            + json.dumps({"equations": equations}, separators=(",", ":")))


def parse_response(raw):
    try:
        result = json.loads(raw)
    except json.JSONDecodeError:
        stripped = re.sub(r"^```(?:json)?\s*|\s*```$", "", raw.strip())
        result = json.loads(stripped)
    if not isinstance(result, dict):
        raise ValueError("Response is not a JSON object")
    return result


def save_json(path, data):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n")
    temp.replace(path)


def make_dataset(source, directory):
    """64 distinct equations in ONE source family, with explicit extra marks.

    The extra free-character marks are mathematical task inputs, not physics.
    Splits measure unseen instances only; no claim of held-out family transfer.
    """
    pool = []
    seen = set()
    for b in [0, 4, 16, 20]:
        for m in range(16):
            fixed = {str(j): (m >> i) & 1 for i, j in enumerate([0, 5, 6, 7])}
            task = source.task(b, fixed)
            if task["equation_key"] in seen:
                raise ArithmeticError("Duplicate mathematical instance")
            seen.add(task["equation_key"])
            solved = exact_gf2(task["rows"], task["rhs"])
            if solved["status"] != "solved" or solved["solution_count"] != 1:
                raise ArithmeticError("Expected one exactly verified marked lift")
            answer = {"bits": bits(solved["solution"])}
            if source.check(task, answer)["status"] != "verified":
                raise ArithmeticError("Invalid training target")
            pool.append({"task": task, "answer": answer})
    random.Random(20261024).shuffle(pool)
    splits = {"train": pool[:36], "valid": pool[36:48], "test": pool[48:]}
    directory = Path(directory); directory.mkdir(parents=True, exist_ok=True)
    for name, items in splits.items():
        (directory/f"{name}.tasks.jsonl").write_text("".join(json.dumps(x)+"\n" for x in items))
        chats = [{"messages": [{"role": "system", "content": SYSTEM}, {"role": "user", "content": prompt(x['task'])},
                               {"role": "assistant", "content": json.dumps(x['answer'], separators=(',', ':'))}]} for x in items]
        (directory/f"{name}.jsonl").write_text("".join(json.dumps(x)+"\n" for x in chats))
    save_json(directory/"manifest.json", {"source_id": source.identity, "counts": {k: len(v) for k, v in splits.items()},
        "unique_equations": len(seen), "split": "disjoint equation instances within one fixed E8 reflection family",
        "family_generalization_measured": False, "additional_marks": "a0,a5,a6,a7 explicitly specified as mathematical inputs",
        "labels": "exact Gaussian elimination plus direct original lattice action check"})
    return splits


def read_tasks(path):
    return [json.loads(line)["task"] for line in Path(path).read_text().splitlines() if line.strip()]


class Ollama:
    def __init__(self, model="qwen3.6:35b", endpoint="http://127.0.0.1:11434", schema=True, think=False, num_predict=256, tools=False, response_schema=None):
        self.model, self.endpoint, self.schema = model, endpoint.rstrip("/"), schema
        self.think, self.num_predict = think, num_predict
        self.tools = tools
        self.response_schema = response_schema
        inventory = json.load(urllib.request.urlopen(self.endpoint + "/api/tags", timeout=10))
        row = next((x for x in inventory["models"] if x["name"] == model), None)
        if row is None:
            raise ValueError("Requested local model is not installed")
        self.identity = {"provider": "ollama", "model": model, "digest": row["digest"], "schema": schema,
                         "think": think, "num_predict": num_predict, "temperature": 0, "seed": 20261024, "tools": tools,
                         "response_schema": response_schema}

    def __call__(self, messages):
        body = {"model": self.model, "messages": messages, "stream": False, "think": self.think,
                "format": self.response_schema if self.response_schema is not None else ((PROGRAM_SCHEMA if self.tools else PROPOSAL_SCHEMA) if self.schema else "json"),
                "options": {"temperature": 0, "seed": 20261024, "num_predict": self.num_predict}}
        request = urllib.request.Request(self.endpoint+"/api/chat", data=json.dumps(body).encode(),
                                         headers={"Content-Type": "application/json"})
        started = time.monotonic()
        response = json.load(urllib.request.urlopen(request, timeout=120))
        return response["message"]["content"], {"elapsed_seconds": time.monotonic()-started,
            "prompt_tokens": response.get("prompt_eval_count"), "generated_tokens": response.get("eval_count"),
            "done_reason": response.get("done_reason"), "model": response.get("model")}


def run_task(source, task, backend, outdir, max_attempts=3, initial=None, allow_tools=False):
    system = TOOL_SYSTEM if allow_tools else SYSTEM
    identity = {"task": task["id"], "source": source.identity, "backend": backend.identity,
                "max_attempts": max_attempts, "initial": initial,
                "engine": digest(Path(__file__).read_text()),
                "mathematical_checker": digest(Path(__file__).with_name('mathcore.py').read_text()),
                "system": system, "allow_tools": allow_tools}
    key = digest(identity)
    path = Path(outdir)/f"{key}.json"
    if path.exists():
        result = json.loads(path.read_text())
        if result["identity"] != identity:
            raise ValueError("Stale run identity")
        for attempt in result.get('attempts', []):
            if attempt.get('candidate') is not None:
                if source.check(task, attempt['candidate']) != attempt['report']:
                    raise ValueError("Persisted mathematical certificate was modified")
        if result.get('attempts') and result['status'] != result['attempts'][-1]['report']['status']:
            raise ValueError("Persisted verdict was modified")
        return result, True
    messages = [{"role": "system", "content": system}, {"role": "user", "content": prompt(task)}]
    attempts = []
    for index in range(max_attempts):
        try:
            if index == 0 and initial is not None:
                raw, metrics = json.dumps(initial), {"elapsed_seconds": 0, "origin": "declared_test_input"}
            else:
                raw, metrics = backend(messages)
        except Exception as e:
            result = {"identity": identity, "task_id": task['id'], "attempts": attempts,
                      "status": "backend_error", "error": f"{type(e).__name__}: {e}"}
            save_json(path, result)
            return result, False
        try:
            proposal = parse_response(raw)
            candidate = proposal
            tool_trace = None
            if allow_tools and proposal == {"operation": "gf2_elimination"}:
                solved = exact_gf2(task["rows"], task["rhs"])
                if solved["status"] == "solved":
                    candidate = {"bits": bits(solved["solution"])}
                else:
                    candidate = {"status": "inconsistent", "rows": [i for i in range(len(task["rows"]))
                        if solved["row_combination"] >> i & 1]}
                tool_trace = {"operation": "gf2_elimination", "result": solved,
                              "provenance": "host_exact_solver_not_neural_arithmetic"}
            report = source.check(task, candidate)
        except (ValueError, KeyError, TypeError) as e:
            candidate, proposal, tool_trace = None, None, None
            report = {"status": "invalid", "error": str(e)}
        attempts.append({"raw": raw, "proposal": proposal, "candidate": candidate, "tool_trace": tool_trace,
                         "report": report, "metrics": metrics})
        messages.append({"role": "assistant", "content": raw})
        if report["status"].startswith("verified"):
            break
        messages.append({"role": "user", "content": "Exact checker feedback: " + json.dumps(report)
                         + ". Revise the eight bits and satisfy ALL original equations."})
    result = {"identity": identity, "task_id": task["id"], "attempts": attempts,
              "status": attempts[-1]["report"]["status"], "messages": messages,
              "weight_learning": bool(backend.identity.get('adapter_sha256')), "feedback_used": len(attempts) > 1,
              "host_solver_used": any(x.get('tool_trace') for x in attempts),
              "physical_completion": False}
    save_json(path, result)
    return result, False


def benchmark(source, tasks, backend, outdir, max_attempts=3, allow_tools=False):
    started = time.monotonic(); results = []; resumed = 0
    for i, task in enumerate(tasks):
        result, cached = run_task(source, task, backend, outdir, max_attempts, allow_tools=allow_tools)
        resumed += int(cached); results.append(result)
        print(json.dumps({"task": i, "status": result["status"], "attempts": len(result["attempts"]), "cached": cached}), flush=True)
    live = [x for r in results for x in r["attempts"] if x["metrics"].get("origin") != "declared_test_input"]
    summary = {"backend": backend.identity, "tasks": len(tasks), "verified_first_attempt": sum(
        bool(r['attempts']) and r['attempts'][0]['report']['status'].startswith('verified') for r in results),
        "verified_after_feedback": sum(r["status"].startswith("verified") for r in results),
        "model_calls_recorded": len(live), "resumed_tasks": resumed,
        "sum_recorded_generation_seconds": sum(x["metrics"]["elapsed_seconds"] for x in live),
        "this_invocation_seconds": time.monotonic()-started,
        "results": [r["task_id"] for r in results], "task_set_digest": digest([t["id"] for t in tasks]),
        "weight_learning": bool(backend.identity.get('adapter_sha256')), "scope": "finite marked reflection instances, one mathematical family"}
    summary['host_solver_tasks'] = sum(bool(r.get('host_solver_used')) for r in results)
    save_json(Path(outdir)/"summary.json", summary)
    return summary
