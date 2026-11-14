"""Run the independent OriginLab research campaigns."""
from pathlib import Path
from datetime import datetime, timezone
from importlib.metadata import version
from tempfile import NamedTemporaryFile
import fcntl
import re
from concurrent.futures import ProcessPoolExecutor
import argparse
import json
import os
import platform
import subprocess
import sys
import time

from language import Term
from search import proposal_prompt, proposals, run
from source import NATIVE, evaluate_original, source_input_key
from backends import inverse_rank, lean_check
from verify import car_replay, clock_test, digest, eigenvalue_test, r102

HERE = Path(__file__).resolve().parent
CAMPAIGNS = ("r113", "r113_synthesis", "source_bridge", "marked_search", "double_cell", "moving_carrier", "affine_polarization", "affine_operators", "affine_car", "r102", "bell", "clock", "identifiability", "four_mark_selection",
             "spectral_source_controls", "empirical_source_audit", "p0", "f1")
FINITE_CONTROLS = ("double_cell", "moving_carrier", "affine_polarization", "affine_operators", "affine_car")


def atomic_json(path, payload):
    with NamedTemporaryFile(mode="w",dir=path.parent,delete=False,encoding="utf-8") as f:
        temp = Path(f.name)
        try:
            json.dump(payload,f,indent=2,sort_keys=True,allow_nan=False)
            f.write("\n")
            f.flush()
            os.fsync(f.fileno())
            os.replace(temp,path)
        finally:
            temp.unlink(missing_ok=True)


def atlas_inputs(data):
    ids = {"problem:f1-session-synthesis-20261113",
           "problem:originlab-double-cell-parity-20261114",
           "problem:automatic-origin-response-decision",
           "problem:original-joint-short-time-kernel-20261105"}
    return {
        problem["id"]: {
            "findings":[{"id":f["id"],"status":f["status"],"premise_ids":f.get("premise_ids",[]),
                         "sources":[s for s in f.get("sources",[]) if isinstance(s,dict)]}
                        for f in problem.get("findings",[])],
            "open_tasks":[t["id"] for t in problem.get("tasks",[]) if t.get("status") == "open"],
        }
        for problem in data["problems"] if problem["id"] in ids
    }


def current_context(root):
    from hashlib import sha256
    progress = root/"docs/research/fundamental-plan/progress.json"
    plan = json.loads(progress.read_text())
    atlas = root/"semantic_explorer/research/program.json"
    raw = atlas.read_bytes()
    data = json.loads(raw)
    return dict(plan_revision=plan["revision"],next_task=plan["next_task"],
                plan_recorded_atlas_revision=plan.get("current_atlas_revision"),
                atlas_program_sha256=sha256(raw).hexdigest(),
                atlas_revision=data["revision"],
                atlas_inputs=atlas_inputs(data),
                git_head=subprocess.check_output(["git","-C",str(root),"rev-parse","HEAD"],
                                                 text=True).strip(),
                python=platform.python_version(),sympy=version("sympy"))


def campaign_context(args):
    from source_bridge import INPUTS

    if args.finite_controls_only:
        return dict(plan_revision=None,atlas_revision=None,atlas_provenance="not_used_for_finite_controls",
                    git_head=subprocess.check_output(["git","-C",str(args.repo),"rev-parse","HEAD"],
                                                     text=True).strip(),
                    python=platform.python_version(),sympy=version("sympy"))
    context = current_context(args.repo)
    domain = args.domain or args.repo/NATIVE
    context["domain_sha256"] = digest(domain)
    context["source_input_key"] = source_input_key(domain,args.repo)
    if "source_bridge" in args.campaign:
        context["source_bridge_inputs"] = {
            str(path): digest(args.repo/path) if (args.repo/path).is_file() else None
            for path in (NATIVE, *INPUTS)
        }
    context["extra_inputs"] = {str(p.resolve()):digest(p)
                               for p in (args.proposal_file,args.finite_observations,args.r113_synthesis_input,
                                         args.marked_search_input) if p}
    return context


def finite_observations(path):
    if path is None:
        return []
    raw = json.loads(path.read_text())
    if not isinstance(raw,list):
        raise ValueError("finite observations must be a JSON list")
    for item in raw:
        if (not isinstance(item,dict) or set(item) != {"dimension","word","value","evidence_id"}
                or type(item["dimension"]) is not int or not 1 <= item["dimension"] <= 12
                or not isinstance(item["word"],list) or not 1 <= len(item["word"]) <= 12
                or any(type(t) is not int or t not in (0,1) for t in item["word"])
                or not isinstance(item["value"],str)
                or re.fullmatch(r"-?[0-9]+(?:/[1-9][0-9]*)?",item["value"]) is None
                or not isinstance(item["evidence_id"],str) or not item["evidence_id"].strip()):
            raise ValueError("invalid finite observation, expected dimension, word, rational value and evidence_id")
    return raw


def execute(campaign, root, args, suggestion_rows, retained_rows, prior_questions):
    if campaign == "benchmarks":
        from benchmarks import benchmark_suite
        return benchmark_suite(root,args.lean_toolchain)
    if campaign == "r113":
        return evaluate_original(root,args.domain)
    if campaign == "r113_synthesis":
        from r113_synthesis import from_input, native_synthesis
        if args.r113_synthesis_input:
            return from_input(json.loads(args.r113_synthesis_input.read_text()))
        domain = args.domain or root/NATIVE
        return native_synthesis(root,domain,evaluate_original(root,domain))
    if campaign == "source_bridge":
        from source_bridge import replay
        return replay(root)
    if campaign == "marked_search":
        from marked_search import run as marked_run
        return marked_run(json.loads(args.marked_search_input.read_text()) if args.marked_search_input else None)
    if campaign == "double_cell":
        from double_cell import replay
        return replay()
    if campaign == "moving_carrier":
        from moving_carrier import replay
        return replay()
    if campaign == "affine_polarization":
        from affine_polarization import controls
        return controls()
    if campaign == "affine_operators":
        from affine_operators import controls
        return controls()
    if campaign == "affine_car":
        from affine_car import controls
        return controls()
    if campaign == "r102":
        result = r102(3,args.rank_limit)
        result["inverse_smt"] = inverse_rank(3,args.rank_limit)
        result["lean"] = lean_check(HERE/"R102.lean",args.lean_toolchain)
        if result["inverse_smt"]["status"] == "EXHAUSTED" and (
            result["inverse_smt"]["solutions"] != sorted({x["n"] for x in result["solutions"]})
        ):
            raise ArithmeticError("inverse synthesis and exact scan disagree")
        return result
    if campaign == "bell":
        return car_replay(root,args.word_depth)
    if campaign == "clock":
        return clock_test()
    if campaign == "identifiability":
        return eigenvalue_test()
    if campaign in ("four_mark_selection","spectral_source_controls","empirical_source_audit"):
        from autonomy import execute as execute_job
        return execute_job({"lane":campaign},root)
    return run(args.size,args.seed,args.mutations,min(args.word_depth,6),
               tuple(map(Term.parse,suggestion_rows)),tuple(map(Term.parse,retained_rows)),prior_questions,
               finite_observations(args.finite_observations))


def one_round(args, directory, index, retained, prior_questions, manifest):
    from evidence_gates import assess

    directory.mkdir(parents=True,exist_ok=False)
    proposal_rows, model_status = [], {}
    if "p0" in args.campaign:
        if args.proposal_file:
            raw = json.loads(args.proposal_file.read_text())
            if not isinstance(raw,list) or len(raw) > 128:
                raise ValueError("proposal file must be a JSON array of up to 128 terms")
            proposal_rows.extend(Term.parse(x).data() for x in raw)
            model_status["proposal_file"] = {"sha256":digest(args.proposal_file),"count":len(raw)}
        for model in args.model:
            try:
                rows = [x.data() for x in proposals(model,proposal_prompt(retained),args.model_timeout)]
                proposal_rows.extend(rows)
                model_status[model] = dict(status="proposed",count=len(rows))
            except (OSError,ValueError,subprocess.SubprocessError) as exc:
                model_status[model] = dict(status="unavailable",error=str(exc))
        atomic_json(directory/"proposals.json",dict(status=model_status,terms=proposal_rows))
    outcomes = {}

    def record(campaign, future):
        try:
            result = future()
        except Exception as error:
            result = dict(status="ERROR",error=f"{type(error).__name__}: {error}")
        outcomes[campaign] = result
        atomic_json(directory/f"{campaign}.json",result)
        print(f"{campaign}: {result['status']}",flush=True)

    # Start with the original R113 path before running independent search lanes.
    independent = [c for c in args.campaign if c != "f1"]
    if "r113" in independent:
        independent.remove("r113")
        record("r113",lambda:execute("r113",args.repo,args,[],[],[]))
    if args.jobs == 1:
        for c in independent:
            record(c,lambda:execute(c,args.repo,args,proposal_rows,retained,prior_questions))
    else:
        with ProcessPoolExecutor(max_workers=args.jobs) as pool:
            jobs = {c:pool.submit(execute,c,args.repo,args,proposal_rows,retained,prior_questions) for c in independent}
            for c,future in jobs.items():
                record(c,future.result)
    if "f1" in args.campaign:
        def f1_frontier():
            source_result = outcomes.get("r113")
            bell = outcomes.get("bell",{}).get("discriminator",{})
            return dict(status="OPEN",source_status=(source_result or {}).get("status","NOT_RUN"),
                        atlas_problem="problem:f1-session-synthesis-20261113",
                        next_original_inputs=[
                            "absolute weighted endpoint Gram and bound empty reference",
                            "same-source four oriented operations and charged primitive loop dictionary",
                            "original charged ordered field-to-word map",
                            "time and regional forms between distinct local embeddings",
                        ],
                        finite_separator=bell.get("questions",[])[:1],
                        r113_synthesis=outcomes.get("r113_synthesis",{}).get("status","NOT_RUN"),
                        original_source_bridge_status=outcomes.get("source_bridge",{}).get("status","NOT_RUN"),
                        r113_next_original_observation=outcomes.get("r113_synthesis",{}).get("next_original_observation"),
                        finite_marked_search_status=outcomes.get("marked_search",{}).get("status","NOT_RUN"),
                        double_cell_control_status=outcomes.get("double_cell",{}).get("status","NOT_RUN"),
                        moving_carrier_control_status=outcomes.get("moving_carrier",{}).get("status","NOT_RUN"),
                        affine_polarization_status=(source_result or {}).get("affine_polarization",{}).get("status","NOT_RUN"),
                        affine_operator_source_status=outcomes.get("affine_operators",{}).get("status","NOT_RUN"),
                        free_car_doublet_source_status=outcomes.get("affine_car",{}).get("status","NOT_RUN"),
                        mark_alphabet_completion=outcomes.get("affine_operators",{}).get("next_original_observation"),
                        physical_evidence=assess(source_result,outcomes.get("r113_synthesis")),
                        separator_requires_original_field_identification=True,
                        energy_moments=(source_result or {}).get("energy_moments"),
                        regional_motion=((source_result or {}).get("physical_reconstruction",{})
                                         .get("regional",{}).get("regional_motion")),
                        source_obligations=(source_result or {}).get("obligations",{}))
        record("f1",f1_frontier)
    p0 = outcomes.get("p0",{})
    gates = assess(outcomes.get("r113"),outcomes.get("r113_synthesis"))["gates"]
    next_action = ("Review the exact finite controls and their original-source follow-up."
                   if args.finite_controls_only else
                   "Produce the original marked endpoint packet or equivalent weighted "
                   "form before claiming physical transfer or regional time.")
    summary = dict(round=index,created_at=datetime.now(timezone.utc).isoformat(),
                   context=manifest,gates=gates,physical_solution=False,
                   outcomes={k:v["status"] for k,v in outcomes.items()},
                   pending_experiments=p0.get("questions",[])[:3],
                   next_action=next_action,
                   survivor_count=len(p0.get("survivors",[])))
    atomic_json(directory/"summary.json",summary)
    return summary,p0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo",type=Path,default=HERE.parents[1])
    parser.add_argument("--out",type=Path,required=True)
    parser.add_argument("--domain",type=Path)
    parser.add_argument("--r113-synthesis-input",type=Path,help="Use a declared finite oracle in the R113 synthesis campaign")
    parser.add_argument("--marked-search-input",type=Path,help="Replay finite phase, clock and cell observations")
    parser.add_argument("--campaign",action="append",choices=(*CAMPAIGNS,"benchmarks"))
    parser.add_argument("--finite-controls-only",action="store_true",
                        help="Run public finite controls without a local Atlas or original source")
    parser.add_argument("--rank-limit",type=int,default=96)
    parser.add_argument("--word-depth",type=int,default=5)
    parser.add_argument("--size",type=int,default=4)
    parser.add_argument("--mutations",type=int,default=64)
    parser.add_argument("--seed",type=int,default=731)
    parser.add_argument("--lean-toolchain")
    parser.add_argument("--proposal-file",type=Path)
    parser.add_argument("--finite-observations",type=Path)
    parser.add_argument("--model",action="append",default=[])
    parser.add_argument("--model-timeout",type=int,default=120)
    parser.add_argument("--jobs",type=int,default=2)
    parser.add_argument("--rounds",type=int,default=1,help="New rounds this invocation")
    parser.add_argument("--watch-seconds",type=int,default=0,
                        help="Wait for source changes between rounds; 0 runs immediately")
    parser.add_argument("--resume",action="store_true")
    args = parser.parse_args()
    args.repo = args.repo.resolve()
    args.domain = args.domain.resolve() if args.domain else None
    defaults = FINITE_CONTROLS if args.finite_controls_only else CAMPAIGNS
    args.campaign = list(dict.fromkeys(args.campaign or defaults))
    if args.finite_controls_only and (set(args.campaign)-set(FINITE_CONTROLS) or args.watch_seconds):
        parser.error("the Atlas-free mode only supports finite controls and runs immediately")
    if args.r113_synthesis_input and "r113_synthesis" not in args.campaign:
        parser.error("--r113-synthesis-input requires the r113_synthesis campaign")
    if args.marked_search_input and "marked_search" not in args.campaign:
        parser.error("--marked-search-input requires the marked_search campaign")
    if not (1 <= args.size <= 6 and 0 <= args.mutations <= 10000 and
            1 <= args.jobs <= 16 and 1 <= args.word_depth <= 12 and
            3 <= args.rank_limit <= 1000 and args.rounds >= 1 and
            args.watch_seconds >= 0 and args.model_timeout > 0):
        parser.error("search parameters outside supported bounds")
    args.out.mkdir(parents=True,exist_ok=True)
    lock = (args.out/".writer").open("a")
    try:
        fcntl.flock(lock,fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        lock.close()
        parser.error("a campaign writer already owns this output")
    try:
        state_path = args.out/"checkpoint.json"
        state = json.loads(state_path.read_text()) if args.resume and state_path.exists() else None
        if state_path.exists() and not args.resume:
            parser.error("output already contains a campaign, use --resume or a fresh directory")
        if args.resume and state is None:
            parser.error("no checkpoint to resume")
        prior = state["round"] if state else 0
        retained = state["survivors"] if state else []
        prior_questions = state.get("questions",[]) if state else []
        configuration = {k:(str(v.resolve()) if isinstance(v,Path) else v)
                         for k,v in vars(args).items() if k not in {"resume","rounds","out","jobs","watch_seconds"}}
        if not args.finite_controls_only:
            configuration.pop("finite_controls_only")
        engine = {p.name:digest(p) for p in [*HERE.glob("*.py"), HERE/"gates.json",HERE/"R102.lean"]}
        if state and (configuration != state["configuration"] or engine != state["engine"]):
            parser.error("configuration or engine changed, start a fresh run")
        had_error = False
        for i in range(prior+1,prior+args.rounds+1):
            context = campaign_context(args)
            if args.watch_seconds and state is not None:
                while state["context"] == context:
                    time.sleep(args.watch_seconds)
                    context = campaign_context(args)
            target = args.out/f"round-{i:04d}"
            if target.exists():
                # Recover an interrupted round without discarding its artifacts.
                archive = target.with_name(target.name+f".interrupted-{time.time_ns()}")
                target.rename(archive)
            args.seed = configuration["seed"]+i-1
            atomic_manifest = dict(context,seed=args.seed,configuration=configuration,engine=engine)
            summary,p0 = one_round(args,target,i,retained,prior_questions,atomic_manifest)
            retained = [x["term"] for x in p0.get("survivors",[])[:128]]
            prior_questions = p0.get("questions",[])[:3]
            had_error |= any(v in {"ERROR","BENCHMARKS_INCOMPLETE"} for v in summary["outcomes"].values())
            state = dict(round=i,survivors=retained,questions=prior_questions,configuration=configuration,
                         engine=engine,context=context)
            atomic_json(state_path,state)
            print(f"saved {target/'summary.json'}",flush=True)
        return 1 if had_error else 0
    finally:
        lock.close()


if __name__ == "__main__":
    sys.exit(main())
