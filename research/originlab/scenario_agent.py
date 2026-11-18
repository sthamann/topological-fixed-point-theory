"""Test one scenario by a model-written prediction computed from the original sources.

A model with read-only file tools inspects the sources of the contexts a
scenario names, derives one exact consequence of the scenario and compares
it with a recorded value that the derivation did not use. The script may
read the repository and run its existing checks. It runs without network,
without writes outside a scratch directory and without access to the rest
of the home directory. A second model call referees whether the comparison
can decide anything. Script, verdict and grade are unreviewed.
"""

from hashlib import sha256
from pathlib import Path
from tempfile import TemporaryDirectory
import json
import re
import subprocess
import sys

THEORY = Path(__file__).with_name("THEORY-BRIEF.md")
VERDICTS = {"PREDICTION_MATCHES": "SURVIVES_COMPUTED_PREDICTION",
            "PREDICTION_FAILS": "REFUTED_BY_COMPUTATION",
            "NO_TESTABLE_CONSEQUENCE": "NO_TESTABLE_CONSEQUENCE_FOUND"}
DECISIVE = ("SURVIVES_COMPUTED_PREDICTION", "REFUTED_BY_COMPUTATION")
# Where each anchor group of the theory brief is constructed.
ANCHOR_SOURCES = {
    "A1": ["experiments/theory-contracts/origin-diagnosis-new-hypotheses-20261116/check_word.py",
           "experiments/theory-contracts/origin-diagnosis-new-hypotheses-20261116/WORD.md"],
    "A2": ["research/originlab/e8_family_return.py",
           "docs/research/fundamental-plan/evidence/E8-FAMILY-RETURN-SCALAR-20261116.txt"],
    "A3": ["research/originlab/f1_four_family_transport.py",
           "docs/research/fundamental-plan/evidence/R113-FOUR-FAMILY-TRANSPORT-20261116.txt"],
    "A4": ["docs/research/fundamental-plan/evidence/E8-POSITIVE-MIDDLE-BLOCKS-20261116.txt"],
    "A5": ["docs/research/fundamental-plan/TFPT-E8-P1-LADUNGSWINDUNG-20261112.txt"],
    "A6": ["docs/research/fundamental-plan/TFPT-ZWEI-FUENFERTAKTE-AM-ORIGINALRAND-20261112.txt",
           "docs/research/fundamental-plan/check_tfpt_five_clock_source_20261112.py"],
    "A7": ["docs/research/fundamental-plan/p1-original-marken-und-metrik-20261112/"
           "TFPT-P1-ORIGINAL-MARKEN-UND-METRIK-20261112.txt"],
    "A8": ["verification/v4_flavor_matrix.py"],
    "A9": ["research/originlab/coherent_origin.py"],
    "A10": ["research/originlab/affine_operators.py"],
    "A11": ["verification/v3_em_alpha.py", "verification/v56_unique_attractor.py"],
}
RESEARCH = """You look for one mechanism behind two recorded results of a research programme. The
two anchor groups named below were obtained by different constructions. Use only the Read, Grep and
Glob tools on this repository. Read both constructions closely: the space each value lives on, the
operators, the factors and the normalization.

Then decide whether one operator, space or structure produces values of both groups with the same
factors on the same spaces the record uses. Matching final numbers by other factors is no
mechanism. Finding none is a valid and useful result.

Answer with a JSON object followed by one fenced python block:
  {"id": short lowercase slug,
   "statement": what object the common rule acts on and what it fixes,
   "rule": the explicit mathematical definition, every symbol defined, concrete enough to program,
   "rationale": three to eight sentences comparing the two mechanisms factor by factor,
   "sources_used": [repository paths with line ranges],
   "target": {"value": recorded values, "source": paths and lines}}
```python
<script>
```
The script follows this contract. It uses the standard library, sympy, numpy and scipy with exact
arithmetic where the claim is exact. It may read files below the directory in ORIGINLAB_REPO and
run existing repository checks with subprocess and sys.executable to read off recorded objects. It
writes only below TMPDIR, uses no network, finishes within 240 seconds and asserts its
intermediate facts. It programs the rule once, derives the values of both groups from that single
implementation, compares factors and spaces with the recorded constructions, and prints as its last
line one JSON object
  {"verdict": "PREDICTION_MATCHES" | "PREDICTION_FAILS" | "NO_TESTABLE_CONSEQUENCE",
   "quantity": name, "derived": value, "target": value, "alternative": value an alternative gives,
   "targets": [{"group": anchor group id, "quantity": name, "derived": value, "target": value}],
   "reason": one sentence}.
Use NO_TESTABLE_CONSEQUENCE when you conclude that no common mechanism exists; then the reason
states the structural obstruction, and "rule" may say that none was found."""
GRADES = ("INDEPENDENT_AND_SPECIFIC", "SAME_CONSTRUCTION", "GENERIC_PROPERTY", "FITTED_TO_TARGET",
          "RULE_NOT_USED", "SCRIPT_UNSOUND")
RULES = """You test one scenario of a research programme by an exact computation. The scenario
equates recorded contexts or adds a new rule (a hypothesis). Use only the Read, Grep and Glob tools
on this repository. Read the cited sources and, where useful, the existing verification scripts.

Derive one consequence that holds if the scenario is true and could fail if it is false, and
compare it with an exact recorded value: a number from the theory brief or from a source of the
other side. The target value must not enter the derivation. Good consequences are specific: a
rank, multiplicity, exact rational coefficient, charge spectrum, sign, cross ratio or the
conjugacy of two explicit operators. A property that every object of the class has is no test,
and neither is one formula written in two notations.

Equal final numbers are not enough. Read in the sources how each target value arises: its
factors, the space it is a trace or rank on, its normalization. Your derivation must produce the
same factors on the same space. If the record obtains 1/96 as a unit response over a 96-dimensional
reference and the rule obtains it as (2/3)/64, the rule fails even though the numbers agree. State
this mechanism comparison in the rationale and assert it in the script where it can be computed.

If the scenario adds a hypothesis, the test is stricter. Program the stated rule once, as one
function or operator. From that single implementation derive values of at least two anchor groups
named in "predicts", which come from different constructions. For each target the computation
must differ from the existing construction of that group: do not call, import or re-implement the
script that produced the target, and do not branch on the group to select a separate formula.
Existing scripts may be run only to read off target values. If the rule as stated is too vague to
program, or yields fewer than two such values, answer NO_TESTABLE_CONSEQUENCE.

Answer with a JSON object followed by one fenced python block:
  {"rationale": two to six sentences: the consequence, why it follows from the scenario, where the
                target value is recorded, and which alternative would give a different value,
   "sources_used": [repository paths with line ranges],
   "target": {"value": the recorded value, "source": path and lines}}
```python
<script>
```
The script uses the standard library, sympy, numpy and scipy with exact arithmetic where the
claim is exact. It may read files below the directory in the environment variable ORIGINLAB_REPO
and may run existing repository checks with subprocess and sys.executable. It writes only below
TMPDIR, uses no network, finishes within 240 seconds, asserts its intermediate facts, and prints
as its last line one JSON object
  {"verdict": "PREDICTION_MATCHES" | "PREDICTION_FAILS" | "NO_TESTABLE_CONSEQUENCE",
   "quantity": name, "derived": value, "target": value, "alternative": value an alternative gives,
   "targets": [{"group": anchor group id or null, "quantity": name, "derived": value,
                "target": value}],
   "reason": one sentence}.
List every compared value in "targets". PREDICTION_MATCHES requires all of them to agree.
Choose NO_TESTABLE_CONSEQUENCE when the scenario fixes nothing computable; do not invent a
consequence and do not copy the target into the derivation. A verdict is evidence within the
stated settings and establishes no physical identification."""
REFEREE = """You referee a computed test in a research programme. A scenario equates recorded
contexts or adds a hypothesis. A script derived one consequence and compared it with a recorded
target. Decide whether the outcome bears on the scenario. Answer with one JSON object and nothing
else: {"grade": one of INDEPENDENT_AND_SPECIFIC, SAME_CONSTRUCTION, GENERIC_PROPERTY,
FITTED_TO_TARGET, RULE_NOT_USED, SCRIPT_UNSOUND, "reason": one or two sentences}.
SAME_CONSTRUCTION: both sides are one construction or formula, so agreement was forced; this
includes re-running or re-implementing the computation that produced the target. RULE_NOT_USED: the
printed values do not depend on the hypothesis, or separate formulas per target replace one rule.
GENERIC_PROPERTY: every object of the class has the property. FITTED_TO_TARGET: the target value
or a free parameter tuned to it enters the derivation, or the number agrees while the record
obtains it from different factors or on a different space than the rule does. SCRIPT_UNSOUND: a value is assigned without
being computed, or the assertions do not establish the printed values, or the derivation does not
use the scenario. INDEPENDENT_AND_SPECIFIC only if none of these applies. Be strict."""


def clip(text, size):
    return text if len(text) <= size else text[:size-1]+"…"


def dossier(kb, scenario, limit=40):
    """Premises, hypotheses and findings of the contexts the scenario names."""
    hypotheses = {row["id"]: {key: row.get(key) for key in ("statement", "rule", "predicts")}
                  for row in scenario.get("new_premises", [])}
    named = sorted({*scenario["adopts"], *(p for pair in scenario.get("identifies", []) for p in pair)}
                   - set(hypotheses))
    premises = {p: {"title": kb["premises"][p]["title"]["de"],
                    "statement": clip(kb["premises"][p]["statement"]["de"], 700)} for p in named}
    rows = [f for f in kb["findings"] if set(f["premise_ids"]) & set(named)]
    # Findings with the fewest premises state the setting of one side most directly.
    rows.sort(key=lambda f: (len(f["premise_ids"]), f["id"]))
    per_side, findings = {}, []
    for finding in rows:
        side = tuple(sorted(set(finding["premise_ids"]) & set(named)))
        if per_side.get(side, 0) >= max(4, limit//max(1, len(named))) or len(findings) >= limit:
            continue
        per_side[side] = per_side.get(side, 0)+1
        findings.append({"id": finding["id"], "status": finding["status"], "premises": finding["premise_ids"],
                         "title": finding["title"]["de"], "statement": clip(finding["statement"]["de"], 600),
                         "scope": clip(finding["scope"]["de"], 300),
                         "sources": [{k: s.get(k) for k in ("path", "line_start", "line_end")}
                                     for s in finding["sources"][:4]]})
    return {"scenario": scenario["statement"], "identifies": scenario.get("identifies", []),
            "new_hypotheses": hypotheses, "premises": premises, "findings": findings}


def claude(arguments, prompt, cwd, timeout, budget, model):
    command = ["claude", "-p", *arguments, "--strict-mcp-config", "--setting-sources", "",
               "--no-session-persistence", "--output-format", "json", "--model", model,
               "--max-budget-usd", str(budget)]
    run = subprocess.run(command, input=prompt, text=True, capture_output=True, cwd=cwd,
                         timeout=timeout, check=True)
    answer = json.loads(run.stdout)
    if answer.get("is_error"):
        raise ValueError(f"model call failed: {answer.get('subtype')}")
    return answer["result"], answer.get("total_cost_usd")


def ask(root, model, prompt, timeout, budget):
    tools = ["--tools", "Read,Grep,Glob", "--allowedTools", "Read", "Grep", "Glob"]
    return claude(tools, prompt, root, timeout, budget, model)


def read_answer(text):
    """Split the model answer into its JSON header and the fenced script."""
    blocks = re.findall(r"```python\n(.*?)```", text, re.S)
    if not blocks:
        raise ValueError("the model returned no fenced script")
    # Prose may contain braces, for example a set of eigenvalues, before the header.
    header = None
    for match in re.finditer(r"\{", text):
        try:
            candidate, _ = json.JSONDecoder().raw_decode(text[match.start():])
        except ValueError:
            continue
        if isinstance(candidate, dict) and "rationale" in candidate:
            header = candidate
            break
    script = blocks[-1]
    if header is None or len(script) > 80000:
        raise ValueError("the model returned no usable header and script")
    return header, script


def plain(model, prompt, timeout, budget):
    """Model call without tools in an empty directory."""
    with TemporaryDirectory(prefix="originlab-referee-") as folder:
        return claude(["--tools", ""], prompt, folder, timeout, budget, model)


def referee(model, scenario, header, script, computed, budget, ask=plain):
    """A tool-free model call grades whether the test can decide anything."""
    case = {"scenario": scenario["statement"], "new_hypotheses": scenario.get("new_premises", []),
            "rationale": header.get("rationale"), "target": header.get("target"),
            "sources_used": header.get("sources_used"), "computed": computed, "script": script}
    try:
        text, cost = ask(model, REFEREE+"\n\n"+json.dumps(case, ensure_ascii=False), 300, budget)
        grade = json.loads(re.search(r"\{.*\}", text, re.S)[0])
        if grade.get("grade") not in GRADES:
            raise ValueError("unknown grade")
    except (OSError, ValueError, KeyError, TypeError, subprocess.SubprocessError) as error:
        return {"grade": "REFEREE_UNAVAILABLE", "reason": clip(str(error), 300), "cost_usd": None}
    return {"grade": grade["grade"], "reason": clip(str(grade.get("reason")), 600), "cost_usd": cost}


def profile(root, work):
    """Seatbelt rules: no network, writes only in the scratch directory, no other home files."""
    readable = {str(root), sys.prefix, sys.base_prefix}
    return "\n".join([
        "(version 1)", "(allow default)", "(deny network*)", "(deny file-write*)",
        f'(allow file-write* (subpath "{work}") (literal "/dev/null") (literal "/dev/dtracehelper"))',
        f'(deny file-read* (subpath "{Path.home()}"))',
        *(f'(allow file-read* (subpath "{path}"))' for path in sorted(readable)),
        '(allow file-read-metadata (subpath "/Users"))'])


def run_script(root, script, seconds=300):
    if sys.platform != "darwin":
        return {"status": "SANDBOX_UNAVAILABLE", "reason": "the isolation profile requires macOS sandbox-exec"}
    root = root.resolve()
    with TemporaryDirectory(prefix="originlab-check-") as folder:
        work = Path(folder).resolve()
        (work/"check.py").write_text(script, encoding="utf-8")
        (work/"profile.sb").write_text(profile(root, work))
        # Repository checks started by the script inherit the same isolation.
        env = {"PATH": "/usr/bin:/bin", "TMPDIR": str(work), "HOME": str(work), "PYTHONHASHSEED": "0",
               "PYTHONDONTWRITEBYTECODE": "1", "ORIGINLAB_REPO": str(root)}
        try:
            run = subprocess.run(["sandbox-exec", "-f", str(work/"profile.sb"), sys.executable, "-B",
                                  str(work/"check.py")], cwd=work, env=env, text=True,
                                 capture_output=True, timeout=seconds)
        except subprocess.TimeoutExpired:
            return {"status": "CHECK_TIMED_OUT"}
    lines = [line for line in run.stdout.strip().splitlines() if line.strip()]
    if run.returncode or not lines:
        return {"status": "CHECK_FAILED_TO_RUN", "returncode": run.returncode,
                "stderr": run.stderr.strip()[-1500:]}
    try:
        verdict = json.loads(lines[-1])
        status = VERDICTS[verdict["verdict"]]
    except (ValueError, KeyError, TypeError):
        return {"status": "CHECK_OUTPUT_MALFORMED", "last_line": clip(lines[-1], 400)}
    return {"status": status, "computed": {key: verdict.get(key) for key in
                                           ("verdict", "quantity", "derived", "target", "alternative", "targets", "reason")}}


REPAIR = """A Python check for a research programme failed to run. Return only the corrected
script in one fenced python block. Repair programming errors and brittle steps, for example an
assertion on the literal text of a file. Never delete or weaken an assertion that compares a
computed value with an expected one: if such a comparison fails, that is the finding. Then keep the
computation, drop the assert, and let the final JSON line report PREDICTION_FAILS with that
quantity, its derived and target values. Keep the output contract: the last printed line is one
JSON object with verdict, quantity, derived, target, alternative, targets and reason."""


def run_with_repairs(root, script, model, budget, repairs=3, ask=None, run_script=None):
    """Run a check; a tool-free model call repairs a script that stops before its verdict."""
    ask, run_script = ask or plain, run_script or globals()["run_script"]
    spent, outcome = 0.0, run_script(root, script)
    for _ in range(repairs):
        if outcome["status"] not in ("CHECK_FAILED_TO_RUN", "CHECK_OUTPUT_MALFORMED", "CHECK_TIMED_OUT"):
            break
        report = json.dumps({key: outcome.get(key) for key in ("status", "stderr", "last_line")})
        try:
            text, cost = ask(model, REPAIR+"\n\n# Failure\n\n"+report+"\n\n# Script\n\n```python\n"
                             +script+"\n```", 300, min(budget, 1.0))
            blocks = re.findall(r"```python\n(.*?)```", text, re.S)
            if not blocks:
                break
        except (OSError, ValueError, KeyError, TypeError, subprocess.SubprocessError):
            break
        spent += cost or 0
        script = blocks[-1]
        outcome = run_script(root, script)
        outcome["repaired"] = True
    return outcome, script, spent


def graded(status, grade):
    """Only an independent and specific prediction counts for or against a scenario."""
    if grade == "INDEPENDENT_AND_SPECIFIC":
        return status
    # A failed prediction of a rule that was fitted or never used leaves no rule to keep.
    if status == "REFUTED_BY_COMPUTATION" and grade in ("FITTED_TO_TARGET", "RULE_NOT_USED"):
        return "HYPOTHESIS_REJECTED"
    return "TEST_NOT_DISCRIMINATING"


def research(root, groups, model, timeout, budget, rounds=3, ask=ask, run_script=run_script, referee=referee):
    """Search for a common mechanism of two anchor groups, revising after each computed outcome."""
    name = "research-"+"-".join(group.lower() for group in groups)
    base = (RESEARCH+"\n\n# Theory brief\n\n"+THEORY.read_text(encoding="utf-8")+"\n\n# Task\n\n"
            +json.dumps({"anchor_groups": {group: ANCHOR_SOURCES[group] for group in groups}}))
    prompt, spent, attempts, result = base, 0.0, [], None
    for attempt in range(1, rounds+1):
        try:
            text, cost = ask(root, model, prompt, timeout, budget)
            header, script = read_answer(text)
        except (OSError, ValueError, KeyError, TypeError, subprocess.SubprocessError) as error:
            attempts.append({"attempt": attempt, "status": "AGENT_UNAVAILABLE", "error": clip(str(error), 300)})
            break
        spent += cost or 0
        rule = {"id": "hypothesis:"+name, "statement": str(header.get("statement")),
                "rule": str(header.get("rule")), "predicts": list(groups)}
        scenario = {"id": name, "statement": rule["statement"], "adopts": [], "new_premises": [rule]}
        outcome, script, repair_cost = run_with_repairs(root, script, model, budget, run_script=run_script)
        spent += repair_cost
        if outcome["status"] in DECISIVE:
            targets = outcome["computed"].get("targets")
            reached = {row.get("group") for row in targets if isinstance(row, dict)} if isinstance(targets, list) else set()
            if not set(groups) <= reached:
                outcome["computed_status"], outcome["status"] = outcome["status"], "TOO_FEW_INDEPENDENT_TARGETS"
        if outcome["status"] in DECISIVE:
            review = referee(model, scenario, header, script, outcome["computed"], min(budget, 1.0))
            spent += review.get("cost_usd") or 0
            outcome.update(referee=review, computed_status=outcome["status"],
                           status=graded(outcome["status"], review["grade"]))
            # Two groups that one existing construction already produces are one result, not two.
            if review["grade"] == "SAME_CONSTRUCTION" and outcome["computed_status"] == "SURVIVES_COMPUTED_PREDICTION":
                outcome["status"] = "GROUPS_SHARE_ONE_CONSTRUCTION"
        result = {**outcome, "scenario": scenario, "rationale": clip(str(header.get("rationale")), 2500),
                  "target": header.get("target"), "sources_used": header.get("sources_used"),
                  "script": script, "script_sha256": sha256(script.encode()).hexdigest()}
        attempts.append({"attempt": attempt, "status": outcome["status"],
                         "computed_status": outcome.get("computed_status"),
                         "reason": (outcome.get("computed") or {}).get("reason"),
                         "referee": outcome.get("referee")})
        # A surviving rule or a reasoned absence of one ends the thread.
        if outcome["status"] in ("SURVIVES_COMPUTED_PREDICTION", "NO_TESTABLE_CONSEQUENCE_FOUND",
                                 "GROUPS_SHARE_ONE_CONSTRUCTION"):
            break
        feedback = {"previous_rule": rule, "previous_script": script,
                    "outcome": {key: outcome.get(key) for key in
                                ("status", "computed_status", "computed", "referee", "stderr", "last_line")}}
        prompt = (base+"\n\n# Your previous attempt and its outcome\n\n"+json.dumps(feedback, ensure_ascii=False)
                  +"\n\nRevise: repair the error, state a different rule that meets the objection, or "
                   "conclude NO_TESTABLE_CONSEQUENCE with the obstruction.")
    if result is None:
        return {"id": name, "status": "AGENT_UNAVAILABLE", "attempts": attempts, "cost_usd": spent,
                "input_role": "unreviewed_model_written_check", "physical_solution": False, "model": model}
    return {"id": name, **result, "attempts": attempts, "cost_usd": spent, "groups": list(groups),
            "input_role": "unreviewed_model_written_check", "physical_solution": False, "model": model}


def investigate(root, scenario, model, timeout, budget, ask=ask, run_script=run_script, referee=referee):
    from atlas_scenarios import knowledge
    kb = knowledge(root)
    base = {"id": scenario["id"], "input_role": "unreviewed_model_written_check",
            "physical_solution": False, "model": model}
    prompt = (RULES+"\n\n# Theory brief\n\n"+THEORY.read_text(encoding="utf-8")+"\n\n# Scenario\n\n"
              +json.dumps(dossier(kb, scenario), ensure_ascii=False, separators=(",", ":")))
    try:
        text, cost = ask(root, model, prompt, timeout, budget)
        header, script = read_answer(text)
    except (OSError, ValueError, KeyError, TypeError, subprocess.SubprocessError) as error:
        return {**base, "status": "AGENT_UNAVAILABLE", "error_type": type(error).__name__,
                "error": clip(str(error), 400)}
    outcome, script, repair_cost = run_with_repairs(root, script, model, budget, run_script=run_script)
    cost = (cost or 0)+repair_cost
    if scenario.get("new_premises") and outcome["status"] in DECISIVE:
        # One rule has to reach values of two different constructions.
        claimed = {group for row in scenario["new_premises"] for group in row.get("predicts", [])}
        targets = outcome["computed"].get("targets")
        reached = {row.get("group") for row in targets if isinstance(row, dict)} if isinstance(targets, list) else set()
        if len(reached & claimed) < 2:
            outcome["computed_status"] = outcome["status"]
            outcome["status"] = "TOO_FEW_INDEPENDENT_TARGETS"
    if outcome["status"] in DECISIVE:
        review = referee(model, scenario, header, script, outcome["computed"], min(budget, 1.0))
        outcome["referee"] = review
        outcome["computed_status"] = outcome["status"]
        outcome["status"] = graded(outcome["status"], review["grade"])
    return {**base, **outcome, "cost_usd": cost, "rationale": clip(str(header.get("rationale")), 1800),
            "target": header.get("target"), "sources_used": header.get("sources_used"),
            "script": script, "script_sha256": sha256(script.encode()).hexdigest()}
