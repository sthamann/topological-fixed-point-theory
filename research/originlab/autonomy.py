"""Durable autonomous scenario search. Workers run typed jobs from one queue."""
from concurrent.futures import ThreadPoolExecutor
from contextlib import closing
from hashlib import sha256
from importlib.metadata import PackageNotFoundError, version
from itertools import product
from collections import defaultdict
from pathlib import Path
from datetime import datetime,timezone
import argparse
import fcntl
import json
import shutil
import sqlite3
import subprocess
import sys
import time

from main import atomic_json, current_context, finite_observations
from scenarios import DYNAMICS, TOPOLOGIES, discriminate
from source import NATIVE, missing_source_result, source_input_key
from verify import digest

HERE = Path(__file__).resolve().parent
ENGINE = {p.name:digest(p) for p in sorted(HERE.glob("*.py"))}
ENGINE.update({p.name:digest(p) for p in [*HERE.glob("*.lean"), HERE/"gates.json",
                                       HERE/"f1_native_contract.json",HERE/"THEORY-BRIEF.md"]})

CASE_FIELDS = {
    "complete_mark_graphs":(set(),{"maximum"}),
    "clifford_sewing":({"input"},set()),
    "finite_gauge_link":({"order"},set()),
    "link_sewing":({"order"},{"graphs"}),
    "return_generator_search":({"input"},set()),
    "return_moment_reconstruction":({"input"},set()),
    "four_mark_refinement":({"input"},set()),
    "clifford_candidate":({"input"},set()),
    "primitive_charge_compatibility":({"input"},set()),
    "spectral_source_candidate":({"input"},set()),
    "custom_complex":({"simplices"},set()),
    "route_candidate":({"candidate"},set()),
    "compare_reconstructions":({"left","right"},set()),
    "compare_memories":({"left","right"},set()),
    "cap_constraints":({"input"},set()),
    "det_word":({"modes","phase_pi","word","intervals"},{"number_modulus"}),
    "r104_variant":(set(),{"modes","phase_pi","modulus"}),
    "r113_orientation":(set(),{"early_radii","late","momentum","time"}),
    "r113_synthesis_case":({"input"},set()),
    "r113_pair_search":(set(),set()),
    "p1_field_search":(set(),{"input"}),
    "q0_rule_search":(set(),{"depth"}),
    "p1_origin_hypotheses":(set(),{"hypotheses","max_cone_order"}),
    "coherent_origin":(set(),{"stage"}),
    "proof_search":(set(),{"first_cost","shells","toolchain"}),
    "f1_solution":(set(),{"contract","sha256","stage","seconds","toolchain"}),
    "marked_search":(set(),{"input"}),
    "rectangular_intertwiner":({"input"},set()),
    "relational_transition":({"input"},set()),
    "transport_network":({"input"},set()),
    "transported_caps":({"input"},set()),
    "charge_visibility":({"input"},set()),
    "flux_inference":({"input"},set()),
    "oriented_shifts":({"input"},set()),
    "stationary_symmetry":({"input"},set()),
    "process_tomography":({"input"},set()),
    "observable_alphabets":({"input"},set()),
    "gaussian_source":({"input"},set()),
    "joint_source":({"input"},set()),
    "atlas_scenario_case":({"scenario"},set()),
    "origin_dsl":({"size","seed","profile"},{"observations","suggested","mutations","depth"}),
}

SOURCE_FILES = (
    NATIVE,
    Path("verification/predictions_frozen.json"),
    Path("verification/v382_alpha_quillen_exact.py"),
    Path("experiments/theory-contracts/clifford-native-quark-source-bridge-20261110/PROOF.txt"),
    Path("experiments/theory-contracts/original-paper-source-reconciliation-20261102/"
         "contribution_20261106Q/runtime/source_adapter.py"),
)


def packed(data):
    return json.dumps(data,sort_keys=True,separators=(",",":"),allow_nan=False)


def installed_versions(packages):
    versions = {}
    for package in packages:
        try:
            versions[package] = version(package)
        except PackageNotFoundError:
            versions[package] = None
    return versions


def token(data):
    return sha256(packed(data).encode()).hexdigest()


def connection(path):
    db = sqlite3.connect(path,timeout=20,isolation_level=None)
    db.execute("PRAGMA journal_mode=WAL")
    db.execute("PRAGMA synchronous=FULL")
    return db


def setup(path):
    with closing(connection(path)) as db:
        db.execute("""CREATE TABLE IF NOT EXISTS jobs (
            id TEXT PRIMARY KEY, phase TEXT NOT NULL, spec TEXT NOT NULL, priority INTEGER NOT NULL,
            status TEXT NOT NULL DEFAULT 'queued', attempts INTEGER NOT NULL DEFAULT 0,
            result TEXT)""")
        db.execute("CREATE TABLE IF NOT EXISTS cache (key TEXT PRIMARY KEY, result TEXT NOT NULL)")
        # The supervisor holds the only writer lease for this campaign.
        db.execute("UPDATE jobs SET status='queued' WHERE status='running'")


def enqueue(path,phase,specs):
    with closing(connection(path)) as db:
        db.executemany("INSERT OR IGNORE INTO jobs(id,phase,spec,priority) VALUES(?,?,?,?)",
                       [(token([phase,spec]),phase,packed(spec),priority) for priority,spec in specs])


def claim(db,phase):
    return db.execute("""UPDATE jobs SET status='running', attempts=attempts+1 WHERE id=(
        SELECT id FROM jobs WHERE phase=? AND status='queued' ORDER BY priority DESC,id LIMIT 1
        ) RETURNING id,spec,attempts""",(phase,)).fetchone()


# Proof lanes depend on the Lean backend and its time budget.
UNCACHED = {"f1_solution","proof_search","lean","atlas_scenarios","atlas_scenario_case",
            "scenario_proposals","scenario_agent","scenario_research"}


def cache_key(basis,spec):
    """Identify a job by its inputs, independent of epoch, stage and round."""
    data = json.loads(spec)
    if basis is None or data["lane"] in UNCACHED or data.get("sol_timeout"):
        return None
    return token([basis,data])


def worker(path,phase,root,timeout,basis=None):
    reused = 0
    with closing(connection(path)) as db:
        while (row := claim(db,phase)):
            job,spec,attempt = row
            key = cache_key(basis,spec)
            hit = key and db.execute("SELECT result FROM cache WHERE key=?",(key,)).fetchone()
            if hit:
                db.execute("UPDATE jobs SET status='done',result=? WHERE id=?",(hit[0],job))
                reused += 1
                continue
            try:
                limit = max(timeout,json.loads(spec).get("model_timeout",0)+60)
                child = subprocess.run([sys.executable,"-B",__file__,"--job","-","--repo",str(root)],
                                       input=spec,text=True,capture_output=True,timeout=limit,check=True)
                result = json.loads(child.stdout)
                status = "done"
            except (subprocess.SubprocessError,OSError,ValueError) as error:
                stderr = getattr(error,"stderr",None)
                result = {"status":"ERROR","error_type":type(error).__name__,
                          "error":str(error) if isinstance(error,ValueError) else None,
                          "detail":str(stderr)[-3000:]}
                status = "queued" if attempt < 2 else "failed"
            db.execute("UPDATE jobs SET status=?,result=? WHERE id=?",(status,packed(result),job))
            if key and status=="done" and "INCONCLUSIVE" not in str(result.get("status")):
                db.execute("INSERT OR REPLACE INTO cache(key,result) VALUES(?,?)",(key,packed(result)))
    return reused


def result_rows(path,phase):
    with closing(connection(path)) as db:
        return [{"spec":json.loads(spec),"result":json.loads(data),"queue_status":status}
                for spec,data,status in db.execute(
                    "SELECT spec,result,status FROM jobs WHERE phase=? AND result IS NOT NULL ORDER BY priority DESC,id",(phase,))]


def registered_origin_campaigns():
    from main import FINITE_CONTROLS

    already_dispatched = {"p1_origin_hypotheses","f1_solution","f1_native_pair_selection",
                          "f1_reflected_cap","f1_binary_clifford_origin","f1_spherical_pairing_selection"}
    return tuple(name for name in FINITE_CONTROLS if name.startswith(("p1_","f1_"))
                 and name not in already_dispatched)


def execute(spec,root):
    lane = spec["lane"]
    if lane=="registered_origin_campaign":
        from main import execute as run_campaign
        campaign = spec["campaign"]
        if campaign not in registered_origin_campaigns():
            raise ValueError("unknown registered origin campaign")
        return run_campaign(campaign,root,None,[],[],[])
    if lane=="atlas_scenarios":
        from atlas_scenarios import PROGRAM, run
        return run(root) if (root/PROGRAM).is_file() else {"status":"ATLAS_UNAVAILABLE","physical_solution":False}
    if lane=="scenario_proposals":
        from scenario_model import propose
        catalogue = {name:{"required":sorted(need),"optional":sorted(extra)}
                     for name,(need,extra) in CASE_FIELDS.items() if name!="atlas_scenario_case"}
        return propose(root,spec["model"],spec["model_timeout"],spec["budget"],spec["history"],
                       catalogue,spec["known"])
    if lane=="scenario_research":
        from scenario_agent import research
        return research(root,spec["groups"],spec["model"],spec["attempt_timeout"],spec["budget"],spec["rounds"])
    if lane=="scenario_agent":
        from scenario_agent import investigate
        return investigate(root,spec["scenario"],spec["model"],spec["model_timeout"],spec["budget"])
    if lane=="atlas_scenario_case":
        from atlas_scenarios import assess, knowledge
        verdict = assess(knowledge(root),spec["scenario"])
        checks = []
        if verdict["status"]=="ADMISSIBLE_FOR_CHECKS":
            free = {job["lane"] for _,job in specifications(3,Path("."),None,[]) if set(job)=={"lane"}}
            for check in spec["scenario"].get("checks",[]):
                job = check.get("job") if isinstance(check,dict) else None
                job = {"lane":job} if isinstance(job,str) else job
                try:
                    if set(check)!={"job","expect"} or job.get("lane")=="atlas_scenario_case":
                        raise ValueError("a scenario check names one typed job and its expected status")
                    if job!={"lane":job["lane"]} or job["lane"] not in free:
                        valid_case(job)
                except (ValueError,AttributeError,KeyError,TypeError) as error:
                    return {**verdict,"status":"REJECTED_MALFORMED_CHECK","error":str(error)[:300],"checks":checks}
                try:
                    observed = execute(job,root).get("status")
                except Exception as error:  # untrusted job arguments reach many lane validators
                    return {**verdict,"status":"REJECTED_MALFORMED_CHECK","checks":checks,
                            "error":f"{job['lane']}: {type(error).__name__}: {str(error)[:300]}"}
                checks.append({"job":job,"expected_status":check["expect"],"status":observed,
                               "passed":observed==check["expect"],"discriminating":len(job)>1})
            # A typed job does not depend on the scenario. Its outcome is a
            # declared prediction about the record and decides no identification.
            if any(not c["passed"] for c in checks):
                verdict["status"] = "DECLARED_PREDICTION_FAILED"
            elif any(c["discriminating"] for c in checks):
                verdict["status"] = "DECLARED_PREDICTIONS_HOLD"
        return {**verdict,"checks":checks}
    if lane=="discovered_original_packet":
        from original_discovery import replay
        return replay(root,Path(spec["domain"]),spec["entry"])
    if lane=="discovered_q0_comparison":
        from q0_rule_search import native_comparison, search
        return native_comparison(search(spec["depth"]),spec["endpoint"])
    if lane=="complete_mark_graphs":
        from sewing_dynamics import complete_graphs
        return {"status":"EXACT_GRAPH_RANK_SELECTION",**complete_graphs(spec.get("maximum",7))}
    if lane=="clifford_sewing":
        from sewing_dynamics import two_cells
        return two_cells(spec["input"])
    if lane=="finite_gauge_link":
        from sewing_dynamics import finite_gauge_link
        return finite_gauge_link(spec["order"])
    if lane=="link_sewing":
        from sewing_dynamics import link_sewing
        return {"status":"EXACT_FINITE_LINK_SEWING",**link_sewing(spec["order"],spec.get("graphs"))}
    if lane in ("return_generator_search","return_generator_followup"):
        from return_search import completion_search, held_out_search
        return (completion_search if lane=="return_generator_search" else held_out_search)(spec["input"])
    if lane=="return_moment_reconstruction":
        from return_moments import from_input
        return from_input(spec["input"])
    if lane=="q0_rule_search":
        from q0_rule_search import run
        return run(root, Path(spec["domain"]) if "domain" in spec else None, spec.get("depth", 3))
    if lane=="r113_pair_search":
        from r113_pair_search import search
        return search(root)
    if lane=="p1_origin_hypotheses":
        from p1_origin_hypotheses import search
        return search(spec.get("hypotheses"), spec.get("max_cone_order", 32))
    if lane=="coherent_origin":
        from coherent_origin import search
        return search(spec.get("stage",2))
    if lane=="proof_search":
        from proof_search import DEFAULT_TOOLCHAIN, run
        return run(spec.get("first_cost",2),spec.get("shells",2),
                   spec.get("toolchain") or DEFAULT_TOOLCHAIN)
    if lane=="f1_solution":
        from f1_solution import search
        from proof_search import DEFAULT_TOOLCHAIN
        return search(Path(spec["contract"]) if spec.get("contract") else None,spec.get("sha256"),
                      spec.get("stage",2),spec.get("toolchain") or DEFAULT_TOOLCHAIN,spec.get("seconds",20),
                      root=root,domain=Path(spec["domain"]) if "domain" in spec else None,
                      q0_depth=min(spec.get("stage",2)-1,4))
    if lane=="original_affine_return":
        from source_return import affine_return
        return affine_return(root)
    if lane=="original_affine_time":
        from affine_time import original_time
        return original_time(root)
    if lane=="original_folded_time":
        from folded_time import marked_clock
        return marked_clock(root)
    if lane=="original_return_transport":
        from source import evaluate_original
        from source_return import affine_transport_constraint
        return affine_transport_constraint(root,evaluate_original(root,Path(spec["domain"])))
    if lane=="four_mark_selection":
        from selection import clifford_controls, original_quark_control, selection_controls
        return {**selection_controls(),"clifford":clifford_controls(),
                "native_fixed_quark_block":original_quark_control(root)}
    if lane in ("four_mark_refinement","clifford_candidate","primitive_charge_compatibility"):
        from selection import charge_compatibility, matrix_candidate, refinement_search
        checks = {"four_mark_refinement":refinement_search,"clifford_candidate":matrix_candidate,
                  "primitive_charge_compatibility":charge_compatibility}
        return checks[lane](spec["input"])
    if lane=="spectral_source_candidate":
        from spectral_candidate import finite_candidate
        return finite_candidate(root,spec["input"])
    if lane=="spectral_source_controls":
        from spectral_candidate import spectral_controls
        return spectral_controls(root)
    if lane=="empirical_source_audit":
        from empirical import audit
        return audit(root)
    if lane in ("r104_variant","det_word","compare_memories","cap_constraints"):
        import coherence
        if lane=="r104_variant":
            return coherence.controls(root,**{k:v for k,v in spec.items() if k!="lane"})
        if lane=="det_word":
            return coherence.det_word(**{k:v for k,v in spec.items() if k!="lane"})
        if lane=="compare_memories":
            return coherence.compare_memory(spec["left"],spec["right"])
        return coherence.cap_constraints(spec["input"])
    if lane=="r113_orientation":
        from r113 import witness
        return witness(root,**{k:v for k,v in spec.items() if k!="lane"})
    if lane=="r113_synthesis_case":
        from r113_synthesis import from_input
        return from_input(spec["input"])
    if lane=="marked_search":
        from marked_search import run as marked_run
        return marked_run(spec.get("input"))
    if lane=="rectangular_intertwiner":
        from intertwiner import solve
        return solve(spec["input"])
    if lane=="r113_synthesis":
        from r113_synthesis import native_synthesis
        from source import evaluate_original
        domain = Path(spec["domain"])
        return native_synthesis(root,domain,evaluate_original(root,domain))
    if lane=="source_bridge":
        from source_bridge import replay
        return replay(root)
    if lane=="weak_source":
        from weak_source import affine_controls
        return affine_controls()
    if lane=="colour_ward":
        from colour_ward import controls
        return controls()
    if lane=="weak_folded_clock":
        from weak_clock import controls
        return controls(root)
    if lane=="p1_field_search":
        from p1_field_search import controls, from_input
        return from_input(spec["input"]) if "input" in spec else controls(root)
    if lane=="weak_return_rank":
        from weak_return_rank import controls
        from weak_return_completion import controls as completion
        return {**controls(),"completion_orbits":completion()}
    if lane=="four_mark_geometry":
        from four_mark_geometry import controls
        return controls()
    if lane=="original_regional_word":
        from original_paths import original_word
        return original_word(root)
    if lane in ("process_tomography","tomography_controls"):
        import tomography
        return tomography.controls() if lane=="tomography_controls" else tomography.solve(spec["input"])
    if lane in ("observable_alphabets","compiler_alphabets"):
        from observability import compare_alphabets, compiler_control
        return compiler_control(root) if lane=="compiler_alphabets" else compare_alphabets(spec["input"])
    if lane=="joint_source":
        from source import joint_source
        return joint_source(spec["input"])
    if lane in ("relational_transition","transport_network","transported_caps","charge_visibility",
                "flux_inference","oriented_shifts",
                "stationary_symmetry","gaussian_source","relational_controls"):
        import relational
        handlers = {"relational_transition":relational.relative_transition,
                    "transport_network":relational.network,"transported_caps":relational.transported_caps,
                    "charge_visibility":relational.charge_visibility,
                    "flux_inference":relational.flux_inference,
                    "oriented_shifts":relational.oriented_shifts,
                    "stationary_symmetry":relational.covariant_stationary,"gaussian_source":relational.gaussian}
        return relational.controls() if lane=="relational_controls" else handlers[lane](spec["input"])
    if lane=="r102":
        from backends import inverse_rank
        from verify import r102
        return {**r102(3,64),"inverse_smt":inverse_rank(3,64)}
    if lane == "six_routes":
        from routes import run
        return run()
    if lane == "route_candidate":
        from routes import examine
        return examine(spec["candidate"])
    if lane == "custom_complex":
        from space import custom_complex
        return custom_complex(spec["simplices"])
    if lane == "compare_reconstructions":
        from reconstruction import compare_documents
        return compare_documents(spec["left"],spec["right"])
    if lane == "reconstruction":
        from reconstruction import controls
        return controls()
    if lane == "complexes":
        from space import finite_space
        return finite_space(spec["vertices"],spec["dimension"])
    if lane == "rational_equivalence":
        from space import finite_process_space
        return finite_process_space(spec["dimension"],spec["height"])
    if lane == "scenario":
        from scenarios import run
        return run(spec["topology"],spec["resolution"],spec["dynamics"],spec["metric"])
    if lane == "r113":
        from source import evaluate_original
        return evaluate_original(root,Path(spec["domain"]))
    if lane == "stitch":
        from stitch import frontier
        return {"status":"CONDITIONAL_ROUTES",**frontier()}
    if lane == "bridges":
        from bridges import run
        return run()
    if lane == "lean":
        from backends import lean_check
        return lean_check(HERE/"R102.lean",spec["toolchain"])
    if lane == "rewrite":
        from synthesis import equivalent_examples
        return equivalent_examples()
    if lane == "synthesis":
        from synthesis import synthesize
        return synthesize(limit=spec["limit"],split_depth=spec["split_depth"])
    if lane in ("origin_dsl","adaptive_dsl"):
        from search import profile_job, proposal_prompt, proposals
        sol = None
        if spec.get("sol_timeout"):
            model = "gpt-6.1-sol"
            try:
                generated = [t.data() for t in proposals(model,proposal_prompt([]),spec["sol_timeout"])]
                sol = {"status":"PROPOSED","model":model,"terms":generated}
                spec = {**spec,"suggested":[*spec.get("suggested",[]),*generated]}
            except (OSError,ValueError,subprocess.SubprocessError) as error:
                sol = {"status":"UNAVAILABLE","model":model,"error_type":type(error).__name__}
        result = profile_job(spec,spec.get("feedback"))
        return {**result,"excluded_count":len(result["excluded"]),"sol_proposals":sol}
    if lane == "regulator":
        from regulator import certify
        return certify(spec["grid"],spec["cutoff"],spec["time"])
    raise ValueError(f"unknown lane: {lane}")


def specifications(resolution,domain,toolchain,focus,seed_research=True,packets=()):
    yield 100,{"lane":"r113","domain":str(domain)}
    for entry in packets:
        yield 100,{"lane":"discovered_original_packet","domain":str(domain),"entry":entry}
    yield 100,{"lane":"r113_synthesis","domain":str(domain)}
    yield 99,{"lane":"atlas_scenarios"}
    yield 99,{"lane":"q0_rule_search","domain":str(domain),"depth":min(resolution-1,4)}
    yield 98,{"lane":"p1_origin_hypotheses"}
    yield 98,{"lane":"coherent_origin","stage":resolution-1}
    yield 100,{"lane":"source_bridge"}
    yield 99,{"lane":"colour_ward"}
    yield 99,{"lane":"weak_source"}
    yield 99,{"lane":"weak_folded_clock"}
    yield 99,{"lane":"p1_field_search"}
    yield 99,{"lane":"weak_return_rank"}
    yield 98,{"lane":"four_mark_geometry"}
    yield 100,{"lane":"original_affine_time"}
    yield 100,{"lane":"original_folded_time"}
    yield 97,{"lane":"f1_solution","domain":str(domain),"stage":resolution-1,"toolchain":toolchain}
    if seed_research:
        yield 97,{"lane":"proof_search","toolchain":toolchain}
        for campaign in registered_origin_campaigns():
            yield 89,{"lane":"registered_origin_campaign","campaign":campaign}
        yield 98,{"lane":"marked_search"}
        yield 100,{"lane":"original_affine_return"}
        yield 99,{"lane":"complete_mark_graphs"}
        for weight,parity in ((0,[1,1]),(1,[1,1]),(2,[1,1]),(4,[1,1]),(1,[1,-1])):
            yield 95,{"lane":"clifford_sewing","input":{"reference_weight":weight,"parity_signs":parity}}
        for order in (2,4):
            yield 95,{"lane":"finite_gauge_link","order":order}
            yield 95,{"lane":"link_sewing","order":order}
        yield 98,{"lane":"return_generator_search","input":{
            "input_role":"declared_finite_oracle","gap":8,"require_positive_weights":True,
            "projected_return_observations":["81/2",8,2,"1/2"]}}
        yield 98,{"lane":"return_moment_reconstruction","input":{
            "input_role":"declared_finite_oracle","moments":[[[x]] for x in (1,"1/2","5/16","7/32","41/256","61/512")],
            "held_out":[[["365/4096"]],[["547/8192"]]]}}
    yield 100,{"lane":"original_regional_word"}
    yield 100,{"lane":"four_mark_selection"}
    yield 99,{"lane":"spectral_source_controls"}
    yield 98,{"lane":"empirical_source_audit"}
    yield 99,{"lane":"r113_orientation"}
    yield 98,{"lane":"r104_variant"}
    yield 97,{"lane":"six_routes"}
    yield 96,{"lane":"relational_controls"}
    yield 95,{"lane":"reconstruction"}
    yield 94,{"lane":"tomography_controls"}
    yield 93,{"lane":"compiler_alphabets"}
    yield 92,{"lane":"r102"}
    yield 90,{"lane":"bridges"}
    yield 87,{"lane":"complexes","vertices":min(resolution+1,5),"dimension":min(3,resolution)}
    yield 86,{"lane":"rational_equivalence","dimension":2,"height":1 if resolution==3 else 2}
    yield 85,{"lane":"stitch"}
    yield 80,{"lane":"synthesis","limit":64,"split_depth":1}
    yield 70,{"lane":"rewrite"}
    if toolchain:
        yield 65,{"lane":"lean","toolchain":toolchain}
    for grid in (24,48,96):
        yield 60,{"lane":"regulator","grid":grid,"cutoff":4,"time":"1/10"}
    pool = focus or product(TOPOLOGIES,DYNAMICS,("unit","alternating"))
    for topology,dynamics,metric in pool:
        yield 40,{"lane":"scenario","topology":topology,"dynamics":dynamics,
                  "metric":metric,"resolution":min(resolution,10)}


def valid_case(row):
    lane = row.get("lane") if isinstance(row,dict) else None
    if not isinstance(lane,str) or lane not in CASE_FIELDS:
        raise ValueError(f"unknown custom job: {lane!r}")
    required,optional = CASE_FIELDS[lane]
    fields = set(row)-{"lane"}
    if not required<=fields<=required|optional:
        raise ValueError(f"invalid fields for {lane}: {sorted(fields)}")
    if lane=="origin_dsl" and (type(row["size"]) is not int or not 1<=row["size"]<=8 or
            type(row.get("depth",5)) is not int or not 1<=row.get("depth",5)<=8):
        raise ValueError("OriginDSL size and test depth must be between 1 and 8")
    return row


def extra_specs(path):
    if path is None:
        return []
    rows = json.loads(path.read_text())
    if not isinstance(rows,list) or any(not isinstance(row,dict) for row in rows):
        raise ValueError("case manifest must be a list of job objects")
    return list({packed(valid_case(row)):row for row in rows}.values())


def watched_sources(root, domain):
    from affine_time import SOURCES
    from folded_time import SOURCES as FOLDED_SOURCES
    from source_return import ORIGINAL_FILES
    from r113_synthesis import NATIVE_GRADE_CHECK
    from native_weak_bridge import INPUTS as WEAK_INPUTS
    from native_weak_w import INPUTS as W_INPUTS
    from p1_field_search import TANGENT
    from source_bridge import INPUTS
    from f1_reflected_cap import ASSEMBLER
    from r113_radial_order import RADIAL_DATA
    paths = {root/path for path in SOURCE_FILES}
    paths.update(root/path for path in (*ORIGINAL_FILES,*SOURCES,*FOLDED_SOURCES,NATIVE_GRADE_CHECK,NATIVE,
                                       *INPUTS,*WEAK_INPUTS,*W_INPUTS,TANGENT,ASSEMBLER,RADIAL_DATA))
    data = json.loads(domain.read_bytes()) if domain.is_file() else {}
    for record in data.get("source_originals",[]):
        paths.add((root/record["path"]).resolve())
    return {str(path):digest(path) if path.is_file() else None for path in sorted(paths)}


def refinement_focus(rows):
    groups = defaultdict(list)
    for row in rows:
        key = tuple(row["response"]["return_moments"])
        groups[key].append(row)
    # Advance the topology/transfer combinations whose four return moments
    # were indistinguishable, prioritizing distinct mechanism or topology.
    selected = set()
    for group in groups.values():
        if len({(r["topology"],r["dynamics"]) for r in group})>1:
            selected.update((r["topology"],r["dynamics"],r["metric"]) for r in group)
    return sorted(selected)


def dsl_spec(args, resolution, previous, observations, proposals):
    return {"lane":"adaptive_dsl","size":min(resolution+1,args.max_dsl_size),"seed":730+resolution,
            "observations":observations,"suggested":proposals,"feedback":previous,
            "sol_timeout":args.sol_timeout if args.sol_proposals else None}


def research_followups(rows,domain):
    jobs = []
    depth = next((row["spec"]["depth"] for row in rows if row["spec"]["lane"]=="q0_rule_search"),3)
    for row in rows:
        spec,result = row["spec"],row["result"]
        if row["queue_status"]!="done":
            continue
        if spec["lane"]=="return_generator_search" and (probe := result.get("next_discriminating_observation")):
            jobs.append((92,{"lane":"return_generator_followup","input":{
                "gap":result["gap"],"signs":probe["signs"],"witnesses":result["candidates"]}}))
        if spec["lane"]=="discovered_original_packet" and result.get("status","").startswith("EXACT_"):
            jobs.append((97,{"lane":"discovered_q0_comparison","depth":depth,"endpoint":{
                "input_role":result.get("input_role"),"status":result["status"],
                "discovered_packet":result["discovered_packet"],"original_values":result.get("original_values"),
                "q0_rule_probe":result.get("q0_rule_probe")}}))
    original = next((r["result"] for r in rows if r["spec"]["lane"]=="r113"),{})
    anchor = next((r["result"] for r in rows if r["spec"]["lane"]=="original_affine_return"),{})
    if original.get("status")=="BLOCKED" and anchor.get("status")=="EXACT_ORIGINAL_AFFINE_RETURN":
        jobs.append((100,{"lane":"original_return_transport","domain":str(domain)}))
    return jobs


def discovered_source_evidence(original,replays):
    if original.get("status")!="BLOCKED":
        return original,None
    candidates = {row["result"]["discovered_packet"]["sha256"]:row["result"] for row in replays
                  if row["queue_status"]=="done" and row["result"].get("status") in
                  ("EXACT_SUPPLIED_CONTRACTION","EXACT_RELATIVE_CONTRACTION","EXACT_RADIAL_MOMENT_BOUND")}
    if not candidates:
        return original,None
    if len(candidates)>1:
        return original,{"status":"MULTIPLE_SOURCE_PACKETS_REQUIRE_REVIEW","sha256":sorted(candidates)}
    replay, = candidates.values()
    return replay,{"status":"REPLAYED_SOURCE_REQUIRES_REVIEW","packet":replay["discovered_packet"]}


def work_order(original,complete,discovered,errors,backend=None):
    """Name the inputs that a higher search stage cannot produce."""
    from f1_solution import native_registry
    from original_discovery import PACKET_SCHEMA, SOURCE_DIRS
    backend = shutil.which("elan") if backend is None else backend
    items = []
    if complete.get("formal_acceptance_registered") is False:
        items.append({
            "id":"native_acceptance_contract","blocks":sorted(native_registry()["gates"]),
            "reason":"No pinned native F1 contract is registered. Every stage returns the search bound.",
            "open_gate_statements":{gate:row["open"] for gate,row in native_registry()["gates"].items()},
            "pickup":"Register {file, sha256} as 'contract' in f1_native_contract.json. "
                     "The engine digest changes and the supervisor opens a new epoch."})
    if original.get("status")=="BLOCKED":
        items.append({
            "id":"original_source_packet","blocks":["G3","G4"],
            "reason":original.get("reason") or "The native consumer received no original source operation.",
            "first_operation":original.get("first_operation"),"missing":original.get("missing",[]),
            "obligations":original.get("obligations",{}),"rejected_packets":discovered.get("rejected",[]),
            "pickup":f"A {PACKET_SCHEMA} packet with input_role source_reviewed_original, bound to the "
                     f"native domain digest, under {', '.join(SOURCE_DIRS)} is discovered and replayed."})
    if not backend:
        items.append({"id":"lean_backend","blocks":sorted(native_registry()["gates"]),
                      "reason":"elan is not on PATH, so no proof certificate can be checked.",
                      "pickup":"Install elan. The pinned toolchain is fetched on first use."})
    if errors:
        items.append({"id":"failed_jobs","blocks":[],
                      "reason":"These jobs failed twice and contribute no evidence.",
                      "jobs":[{"lane":row["spec"]["lane"],"error_type":row["result"].get("error_type")}
                              for row in errors]})
    proof_inputs = {"native_acceptance_contract","lean_backend"}
    return {"f1_proof_reachable":not any(item["id"] in proof_inputs for item in items),
            "scope":"inputs outside the implemented search grammars; a listed item is a research task, "
                    "and its presence alone establishes no physical claim",
            "items":items}


def completeness(root,domain,toolchain,discovered,complete,far=64):
    """Compute whether this configuration can decide every gate and still grow."""
    from f1_solution import native_registry
    registry = native_registry()
    contract = registry["contract"]
    predicates = dict.fromkeys(registry["gates"])
    if contract is not None and complete.get("same_model_gates"):
        predicates.update(complete["same_model_gates"])
    gates = {gate:{"acceptance_predicate":predicates[gate],
                   "scoped_evidence":{path:(root/path).is_file() for path in row["evidence"]}}
             for gate,row in registry["gates"].items()}
    # A lane whose job is identical at two far stages has exhausted its grammar.
    stage = lambda n: {packed(spec) for _,spec in specifications(n,domain,toolchain,[],seed_research=False)}
    growing = sorted({json.loads(spec)["lane"] for spec in stage(far+1)-stage(far)})
    undecidable = sorted(gate for gate,row in gates.items() if row["acceptance_predicate"] is None)
    packets = len(discovered.get("matching",[]))
    return {"complete":not undecidable and bool(growing) and (packets>0 or contract is not None),
            "gates":gates,"gates_without_acceptance_predicate":undecidable,
            "source_reviewed_packets":packets,
            "lanes_still_growing_with_stage":growing,
            "scope":"software reachability of the registered gates; it establishes no physical claim"}


def scenario_verdicts(rows):
    """Compact record of tested scenarios for the next proposal."""
    verdicts = []
    computed = {row["spec"]["scenario"]["id"]:row["result"] for row in rows
                if row["spec"]["lane"]=="scenario_agent" and row["queue_status"]=="done"}
    for row in rows:
        if row["spec"]["lane"]!="atlas_scenario_case":
            continue
        scenario,result = row["spec"]["scenario"],row["result"]
        verdicts.append({"id":scenario["id"],"statement":scenario["statement"][:400],"scenario":scenario,
                         "status":result.get("status"),
                         "repeated_closed_routes":[r["strategy"] for r in result.get("repeated_closed_routes",[])],
                         "tasks_with_missing_premises":result.get("tasks_with_missing_premises",{}),
                         "error":result.get("error"),
                         "new_premises":scenario.get("new_premises",[]),
                         "computed_comparison":{k:computed[scenario["id"]].get(k) for k in
                                                ("status","computed_status","computed","referee","rationale","target",
                                                 "script_sha256","cost_usd")}
                                               if scenario["id"] in computed else None,
                         "failed_checks":[c["job"]["lane"] for c in result.get("checks",[]) if not c["passed"]]})
    return verdicts


def search_round(args,db,out,epoch,resolution,context,domain,focus,extras,dsl,history=(None,0),
                 scenario_history=(),spent_before=0.0):
    round_key = token([epoch,resolution,context["atlas_program_sha256"],context["git_head"],source_input_key(domain,args.repo),ENGINE,
                       context["watched_sources"],context["discovered_original_sources"],
                       focus,args.lean_toolchain,args.timeout,extras,dsl])
    round_dir = out/round_key
    round_dir.mkdir(exist_ok=True)
    specs = list(specifications(resolution,domain,args.lean_toolchain,focus,
                                seed_research=resolution==args.start_resolution,
                                packets=context["discovered_original_sources"]["matching"]))
    enqueue(db,round_key,specs)
    enqueue(db,round_key,[(75,spec) for spec in extras])
    enqueue(db,round_key,[(55,dsl)])
    basis = token([context["atlas_program_sha256"],context["git_head"],source_input_key(domain,args.repo),
                   ENGINE,context["watched_sources"],context["discovered_original_sources"]])
    reused = 0
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        def run_queued():
            nonlocal reused
            timeout = args.timeout * max(1,resolution-args.start_resolution+1) if args.until_proved else args.timeout
            futures = [pool.submit(worker,db,round_key,args.repo,timeout,basis) for _ in range(args.workers)]
            reused += sum(future.result() for future in futures)
            return result_rows(db,round_key)

        rows = run_queued()
        followups = research_followups(rows,domain)
        if followups:
            enqueue(db,round_key,followups)
            rows = run_queued()
        if args.scenario_model and not args.scenario_proposals:
            agents = []
            from itertools import combinations
            order = ("A1","A7","A5","A8","A11","A2","A9","A10","A6","A3","A4")
            covered = {row["id"] for row in scenario_history}
            pairs = [pair for pair in combinations(order,2)
                     if "research-"+"-".join(group.lower() for group in pair) not in covered]
            threads = [(94,{"lane":"scenario_research","groups":list(pair),"model":args.scenario_model,
                            "attempt_timeout":args.scenario_agent_timeout,"rounds":args.scenario_research_rounds,
                            "model_timeout":args.scenario_research_rounds*(args.scenario_agent_timeout+400),
                            "budget":args.scenario_agent_budget_usd})
                       for pair in pairs[:args.scenario_research_max]]
            if threads:
                enqueue(db,round_key,threads)
                rows = run_queued()
        elif args.scenario_model:
            # Parameter-free jobs of this stage give the model observed outcomes to predict against.
            known = {row["spec"]["lane"]:row["result"].get("status") for row in rows
                     if set(row["spec"])=={"lane"} and row["queue_status"]=="done"}
            ask = {"lane":"scenario_proposals","model":args.scenario_model,
                   "model_timeout":args.scenario_timeout,"budget":args.scenario_budget_usd,
                   "history":list(scenario_history),"known":known}
            enqueue(db,round_key,[(96,ask)])
            rows = run_queued()
            proposed = next(row["result"] for row in rows if row["spec"]==ask)
            cases = []
            for scenario in proposed.get("scenarios",[]):
                try:
                    cases.append((96,valid_case({"lane":"atlas_scenario_case","scenario":scenario})))
                except ValueError:
                    continue
            followups = [*followups,(96,ask),*cases]
            if cases:
                enqueue(db,round_key,cases)
                rows = run_queued()
            # Scenarios consistent with the record get a computed comparison of their sources.
            admissible = [row["spec"]["scenario"] for row in rows
                          if row["spec"]["lane"]=="atlas_scenario_case" and row["queue_status"]=="done"
                          and row["result"].get("status")=="ADMISSIBLE_FOR_CHECKS"
                          and (row["spec"]["scenario"].get("identifies")
                               or row["spec"]["scenario"].get("new_premises"))]
            # Scenarios that add a rule come first: the record cannot produce those itself.
            admissible.sort(key=lambda scenario: not scenario.get("new_premises"))
            current = {scenario["id"] for scenario in admissible}
            backlog = [row["scenario"] for row in scenario_history
                       if row.get("status")=="ADMISSIBLE_FOR_CHECKS"
                       and (row.get("computed_comparison") or {}).get("status") in (None,"AGENT_UNAVAILABLE")
                       and row.get("scenario",{}).get("new_premises") and row["id"] not in current]
            admissible = [*backlog,*admissible]
            agents = [(95,{"lane":"scenario_agent","scenario":scenario,"model":args.scenario_model,
                           "model_timeout":args.scenario_agent_timeout,"budget":args.scenario_agent_budget_usd})
                      for scenario in admissible[:args.scenario_agent_max]]
            # Research threads look for one mechanism behind two anchor groups that no thread covered yet.
            from itertools import combinations
            # Groups on different spaces come first: the E8 matter roots, the 96-dimensional triple,
            # the conformal seam, the flavour lattice and the fixed point.
            order = ("A1","A7","A5","A8","A11","A2","A9","A10","A6","A3","A4")
            covered = {row["id"] for row in scenario_history}
            pairs = [pair for pair in combinations(order,2)
                     if "research-"+"-".join(group.lower() for group in pair) not in covered]
            threads = [(94,{"lane":"scenario_research","groups":list(pair),"model":args.scenario_model,
                            "attempt_timeout":args.scenario_agent_timeout,"rounds":args.scenario_research_rounds,
                            "model_timeout":args.scenario_research_rounds*(args.scenario_agent_timeout+400),
                            "budget":args.scenario_agent_budget_usd})
                       for pair in pairs[:args.scenario_research_max]]
            if agents or threads:
                enqueue(db,round_key,[*agents,*threads])
                rows = run_queued()
    scenarios = [row["result"] for row in rows if row["spec"]["lane"]=="scenario" and row["queue_status"]=="done"]
    baseline = {packed(spec) for _,spec in specs+followups}
    baseline.add(packed(dsl))
    results = {row["spec"]["lane"]:row["result"] for row in rows if packed(row["spec"]) in baseline}
    original = results.get("r113",{})
    bridge = results.get("bridges",{})
    dsl_result = results.get("adaptive_dsl",{})
    questions = discriminate(scenarios)
    first_search = results.get("return_generator_search",{})
    prior_research = dsl["feedback"].get("directed_research_round")
    new_research = {
        "round":round_key,"synthesis_status":first_search.get("status"),
        "synthesis_witnesses":first_search.get("witnesses"),
        "moment_reconstruction":results.get("return_moment_reconstruction"),
        "followups":[r for r in rows if r["spec"]["lane"] in ("return_generator_followup","original_return_transport")],
    } if resolution==args.start_resolution else None
    from space import coverage
    from evidence_gates import assess
    accepted = [r["result"] for r in rows if r["spec"]["lane"]=="f1_solution"
                and r["queue_status"]=="done" and r["result"].get("status") in
                ("F1_PROVED","F1_CONTRACT_PROVED")]
    complete = (max(accepted,key=lambda row:row["status"]=="F1_PROVED")
                if accepted else results.get("f1_solution",{}))
    original_replays = [row for row in rows if row["spec"]["lane"]=="discovered_original_packet"]
    r113_evidence,source_choice = discovered_source_evidence(original,original_replays)
    physical = assess(r113_evidence,results.get("r113_synthesis"),complete)
    frontier = {
        "state":("F1_PROVED" if physical["physical_solution"] else
                 "SOURCE_REVIEW_REQUIRED" if source_choice else
                 {"BLOCKED":"SOURCE_REQUIRED","EXACT_SUPPLIED_CONTRACTION":"OPEN",
                 "ORIGINAL_NATIVE_VALUES_REQUIRE_REVIEW":"SOURCE_REVIEW_REQUIRED"}.get(original.get("status"),"ERROR")),
        "atlas_revision":context["atlas_revision"],
        "discovered_original_sources":context["discovered_original_sources"],
        "registered_origin_campaigns":[{"campaign":row["spec"]["campaign"], "result":row["result"],
                                        "queue_status":row["queue_status"]}
                                       for row in rows if row["spec"]["lane"]=="registered_origin_campaign"],
        "discovered_original_replays":original_replays,
        "discovered_source_selection":source_choice,
        "discovered_q0_comparisons":[row for row in rows if row["spec"]["lane"]=="discovered_q0_comparison"],
        "physical_solution":physical["physical_solution"],
        "physical_evidence":physical,
        "gates":{gate:{"status":"OPEN","required":items}
                 for gate,items in json.loads((HERE/"gates.json").read_text()).items()},
        "primary_original_test":original.get("obligations",{}),
        "native_endpoint_readout":original.get("native_endpoint_readout"),
        "bridge_edges":bridge.get("edges",[]),
        "same_source_search":results.get("stitch",{}),
        "coverage":coverage(),
        "custom_cases":[r for r in rows if r["spec"] in extras],
        "reconstruction":results.get("reconstruction",{}),
        "six_routes":{k:v for k,v in results.get("six_routes",{}).items() if k!="variations"},
        "six_routes_original":original.get("six_routes"),
        "original_reconstruction":original.get("physical_reconstruction"),
        "original_energy_moments":original.get("energy_moments"),
        "weighted_root_doublet":original.get("marked_root_doublet"),
        "original_r113_synthesis":results.get("r113_synthesis"),
        "original_source_bridge":results.get("source_bridge"),
        "weighted_weak_gauge_responses":original.get("weak_gauge_responses"),
        "weighted_weak_clock_sectors":original.get("weak_clock_sectors"),
        "original_weak_folded_clock":results.get("weak_folded_clock"),
        "weighted_p1_charged_jet":original.get("p1_charged_jet"),
        "p1_field_search":results.get("p1_field_search"),
        "q0_rule_search":results.get("q0_rule_search"),
        "p1_origin_hypotheses":results.get("p1_origin_hypotheses"),
        "formal_proof_search":results.get("proof_search"),
        "complete_formal_contract":complete or None,
        "weighted_weak_parity_bridge":original.get("weak_parity_bridge"),
        "weighted_weak_return_grade_two":original.get("weak_return_grade_two"),
        "affine_weak_source":results.get("weak_source"),
        "affine_colour_ward":results.get("colour_ward"),
        "affine_weak_return_rank":results.get("weak_return_rank"),
        "four_mark_geometry":results.get("four_mark_geometry"),
        "finite_marked_search":results.get("marked_search"),
        "generated_complexes":{k:v for k,v in results.get("complexes",{}).items()
                               if k!="representatives_by_homology"},
        "operational_classes":len(results.get("rational_equivalence",{}).get("equivalence_classes",[])),
        "scenario_observables":questions,
        "cellular_questions":results.get("synthesis",{}).get("questions",[]),
        "origin_dsl_questions":dsl_result.get("questions",[]),
        "origin_dsl_feedback":dsl_result.get("feedback",{}),
        "coherent_det_factor":results.get("r104_variant",{}),
        "relational_controls":results.get("relational_controls",{}),
        "process_tomography":results.get("tomography_controls",{}),
        "compiler_observability":results.get("compiler_alphabets",{}),
        "original_relational_paths":original.get("relational_original"),
        "original_regional_word":results.get("original_regional_word",{}),
        "original_process_tomography":original.get("process_tensor"),
        "original_joint_source_checks":original.get("joint_source_checks"),
        "original_ordered_return":original.get("ordered_return"),
        "original_return_moment_reconstruction":original.get("ordered_return_moment_reconstruction"),
        "original_affine_anchor":original.get("affine_return_anchor"),
        "original_affine_time":results.get("original_affine_time"),
        "original_folded_time":results.get("original_folded_time"),
        "weighted_marked_affine_time_grid":original.get("marked_affine_time_grid"),
        "directed_research":new_research or prior_research,
        "sewing_candidates":[r for r in rows if r["spec"]["lane"] in
             ("complete_mark_graphs","clifford_sewing","finite_gauge_link","link_sewing",
              "return_generator_search","return_moment_reconstruction","original_affine_return")],
        "four_mark_original":original.get("four_mark_original"),
        "four_mark_selection":results.get("four_mark_selection",{}),
        "spectral_source_controls":results.get("spectral_source_controls",{}),
        "empirical_source_audit":results.get("empirical_source_audit",{}),
        "open_source_tests":[
            {"source":"original_native_r113",
             "test":"construct the weighted regulated configuration integral and common empty reference",
             "observed":original.get("status")},
            {"source":"original_marked_e8_middle_in_r113",
             "test":"evaluate the two independent intervals of the three signed insertions in the same weighted form",
             "observed":(original.get("marked_affine_time_grid") or {}).get("status","MISSING_WEIGHTED_GRID")},
            {"source":"original_folded_e8_mirror_in_r113",
             "test":"compare the transported full lattice mirror and the translated half-defect, including their respective bound empty lines on the same weighted form",
             "known_obstruction":"the transported lattice mirror keeps the clock only with r(d); the fixed-background half-defect puts one marked line below its proposed reference",
             "observed":original.get("status")},
            {"source":"original_four_mark_operations",
             "test":"square of one physically marked face, then all six cyclic anticommutators and charged weights",
             "observed":original.get("four_mark_original",{}).get("status","MISSING_INPUT")},
        ],
        "r113_orientation_witness":results.get("r113_orientation",{}),
        "charged_original_observation":original.get("charged_bridge"),
        "refinement":{"resolution":resolution,"next_resolution":resolution+1 if resolution<args.max_resolution else None},
    }
    atlas = results.get("atlas_scenarios") or {}
    frontier["atlas_scenarios"] = atlas
    frontier["supplied_scenarios"] = [row for row in rows if row["spec"]["lane"]=="atlas_scenario_case"]
    frontier["scenario_proposals"] = results.get("scenario_proposals")
    frontier["scenario_comparisons"] = [row for row in rows
                                        if row["spec"]["lane"] in ("scenario_agent","scenario_research")]
    for row in frontier["scenario_comparisons"]:
        if row["result"].get("script"):
            folder = out/"scenario-checks"
            folder.mkdir(exist_ok=True)
            (folder/f"{row['result']['id']}.py").write_text(row["result"]["script"],encoding="utf-8")
    tested = {row["id"]:row for row in (*scenario_history,*scenario_verdicts(rows))}
    for row in rows:
        result = row["result"]
        if row["spec"]["lane"]=="scenario_research" and row["queue_status"]=="done":
            scenario = result.get("scenario") or {"id":result["id"],"statement":"","new_premises":[]}
            tested[result["id"]] = {
                "id":result["id"],"statement":scenario["statement"][:400],"scenario":scenario,
                "status":"RESEARCH_THREAD","new_premises":scenario.get("new_premises",[]),
                "attempts":result.get("attempts"),
                "computed_comparison":{key:result.get(key) for key in (
                    "status","computed_status","computed","referee","rationale","target","script_sha256","cost_usd")}}
    # A scenario from the backlog was assessed in an earlier stage; attach its new test.
    for row in rows:
        name = row["spec"].get("scenario",{}).get("id")
        if row["spec"]["lane"]=="scenario_agent" and row["queue_status"]=="done" and name in tested:
            tested[name] = {**tested[name],"computed_comparison":{
                key:row["result"].get(key) for key in ("status","computed_status","computed","referee",
                                                       "rationale","target","script_sha256","cost_usd")}}
    calls = [frontier["scenario_proposals"] or {},*(row["result"] for row in frontier["scenario_comparisons"])]
    spent = spent_before+sum((call.get("cost_usd") or 0)+((call.get("referee") or {}).get("cost_usd") or 0)
                             for call in calls)
    errors = [row for row in rows if row["queue_status"]=="failed"]
    order = work_order(original,complete,context["discovered_original_sources"],errors)
    # A stage makes progress only if a gate, source or acceptance state moves.
    signature = token([frontier["state"],{gate:row["status"] for gate,row in physical["gates"].items()},
                       original.get("status"),original.get("missing"),complete.get("status"),
                       source_choice,[item["id"] for item in order["items"]],
                       sorted(packed(row["spec"]["lane"]) for row in errors)])
    stalled = history[1]+1 if history[0]==signature else 0
    order.update(stalled_rounds=stalled,round=round_key,resolution=resolution,
                 completeness=completeness(args.repo,domain,args.lean_toolchain,
                                           context["discovered_original_sources"],complete))
    frontier["work_order"] = order
    report = {"epoch":epoch,"atlas_revision":context["atlas_revision"],
              "git_head":context["git_head"],"domain_sha256":digest(domain) if domain.is_file() else None,
              "engine_sha256":token(ENGINE),
              "created_at":datetime.now(timezone.utc).isoformat(),
              "watched_sources_sha256":token(context["watched_sources"]),
              "original_status":original.get("status","ERROR"),"scenarios":len(scenarios),
              "discovered_source_selection":source_choice,
              "p1_field_status":results.get("p1_field_search",{}).get("status"),
              "binary_rules_found":len(results.get("synthesis",{}).get("rules",[])),
              "r104_status":results.get("r104_variant",{}).get("status"),
              "r113_orientation_status":results.get("r113_orientation",{}).get("status"),
              "relational_status":results.get("relational_controls",{}).get("status"),
              "tomography_status":results.get("tomography_controls",{}).get("status"),
              "compiler_observability_status":results.get("compiler_alphabets",{}).get("status"),
              "original_regional_word_status":results.get("original_regional_word",{}).get("status"),
              "four_mark_selection_status":results.get("four_mark_selection",{}).get("status"),
              "spectral_source_status":results.get("spectral_source_controls",{}).get("status"),
              "empirical_source_status":results.get("empirical_source_audit",{}).get("status"),
              "directed_research_round":(new_research or prior_research or {}).get("round"),
              "directed_research_status":(new_research or prior_research or {}).get("synthesis_status"),
              "sol_proposal_status":(dsl_result.get("sol_proposals") or {}).get("status","DISABLED"),
              "origin_dsl_feedback":dsl_result.get("feedback",{}),
              "custom_cases":len(extras),
              "generated_complexes":results.get("complexes",{}).get("count"),
              "reconstruction_status":results.get("reconstruction",{}).get("status"),
              "reused_results":reused,"stalled_rounds":stalled,
              "f1_proof_reachable":order["f1_proof_reachable"],
              "search_complete":order["completeness"]["complete"],
              "gates_without_acceptance_predicate":order["completeness"]["gates_without_acceptance_predicate"],
              "lanes_still_growing_with_stage":order["completeness"]["lanes_still_growing_with_stage"],
              "work_order":[item["id"] for item in order["items"]],
              "atlas_scenarios":{
                  "status":atlas.get("status"),"contexts_needed":len(atlas.get("needed_contexts",[])),
                  "minimum_identifications":atlas.get("minimum_identifications"),
                  "identifications_without_recorded_link":len(
                      (atlas.get("scenarios") or [{}])[0].get("contexts_without_recorded_link",[])),
                  "open_inputs":len(atlas.get("open_inputs",[])),
                  "findings_with_changed_source":(atlas.get("evidence_integrity") or {}).get(
                      "findings_with_changed_or_missing_source"),
                  "model":(frontier["scenario_proposals"] or {}).get("status","DISABLED"),
                  "model_cost_usd":(frontier["scenario_proposals"] or {}).get("cost_usd"),
                  "spent_usd":round(spent,2),
                  "supplied":{row["spec"]["scenario"]["id"]:row["result"].get("status")
                              for row in frontier["supplied_scenarios"]},
                  "computed_comparisons":{row["result"].get("id"):row["result"].get("status")
                                          for row in frontier["scenario_comparisons"]}},
              "failures":len(errors),"best_readout":questions[0] if questions else None,
              "physical_solution":physical["physical_solution"],"round":round_key,"resolution":resolution,
              "complete_contract_status":complete.get("status","NOT_RUN"),
              "tool_versions":installed_versions(("sympy","cvc5","egglog","python-flint"))}
    for name,data in (("results",rows),("frontier",frontier),("summary",report)):
        atomic_json(round_dir/f"{name}.json",data)
    atomic_json(round_dir/"engine.json",ENGINE)
    atomic_json(round_dir/"sources.json",context["watched_sources"])
    atomic_json(out/"latest.json",report)
    atomic_json(out/"work-order.json",order)
    atomic_json(round_dir/"review-outbox.json",
                {"status":"SOURCE_REVIEW_REQUIRED","atlas_revision":context["atlas_revision"],
                 "candidate_artifacts":{k:digest(round_dir/f"{k}.json")
                                        for k in ("results","frontier","summary","engine","sources")},
                 "physical_status_changes":[]})
    print(packed(report),flush=True)
    archived = {(item["test_scope"],item["semantic_key"]):item
                for item in (*dsl["feedback"].get("excluded",[]),*dsl_result.get("excluded",[]))
                if "semantic_key" in item}
    feedback = {"retained":[item["term"] for item in dsl_result.get("survivors",[])],
                "questions":dsl_result.get("questions",[]),"excluded":list(archived.values()),
                "directed_research_round":new_research or prior_research,
                "physical_solution":physical["physical_solution"]}
    atomic_json(out/"scenarios.json",list(tested.values()))
    return round_key,refinement_focus(scenarios),feedback,{"signature":signature,"stalled":stalled,
                                                         "scenario_history":list(tested.values())[-200:],
                                                         "scenario_spent_usd":spent}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo",type=Path,default=HERE.parents[1])
    parser.add_argument("--out",type=Path)
    parser.add_argument("--domain",type=Path)
    parser.add_argument("--cases",type=Path,help="Additional typed jobs from a JSON case manifest")
    parser.add_argument("--finite-observations",type=Path)
    parser.add_argument("--proposal-file",type=Path)
    parser.add_argument("--sol-proposals",action="store_true",
                        help="Send the finite default grammar to OpenAI Sol for candidate proposals")
    parser.add_argument("--sol-timeout",type=int,default=45)
    parser.add_argument("--max-dsl-size",type=int,default=5)
    parser.add_argument("--workers",type=int,default=4)
    parser.add_argument("--start-resolution",type=int,default=3)
    parser.add_argument("--max-resolution",type=int,default=4)
    parser.add_argument("--timeout",type=int,default=90)
    parser.add_argument("--lean-toolchain")
    parser.add_argument("--watch-seconds",type=int,default=0,
                        help="Monitor original source revisions after bounded search is exhausted")
    parser.add_argument("--until-proved",action="store_true",
                        help="Keep increasing proof and model stage until native F1 is proved or interrupted")
    parser.add_argument("--scenario-model",
                        help="Claude model alias or id that proposes Atlas scenarios each stage, e.g. opus")
    parser.add_argument("--scenario-timeout",type=int,default=600)
    parser.add_argument("--scenario-budget-usd",type=float,default=2.0,
                        help="Spending limit of one proposal call")
    parser.add_argument("--scenario-agent-max",type=int,default=2,
                        help="Scenarios per stage whose sources a model compares by a sandboxed computation")
    parser.add_argument("--scenario-agent-timeout",type=int,default=900)
    parser.add_argument("--scenario-agent-budget-usd",type=float,default=3.0)
    parser.add_argument("--scenario-proposals",action=argparse.BooleanOptionalAction,default=True,
                        help="Let the model propose scenarios; --no-scenario-proposals runs research threads only")
    parser.add_argument("--scenario-research-max",type=int,default=0,
                        help="Research threads per stage, each seeking one mechanism behind two anchor groups")
    parser.add_argument("--scenario-research-rounds",type=int,default=3,
                        help="Attempts per research thread, each revising after its computed outcome")
    parser.add_argument("--scenario-total-budget-usd",type=float,default=20.0,
                        help="The campaign stops with status 4 once model calls have cost this much")
    parser.add_argument("--stall-rounds",type=int,default=3,
                        help="With --until-proved, stop raising the stage after this many stages without a "
                             "change of gate, source or acceptance state; 0 never stops")
    parser.add_argument("--job",help=argparse.SUPPRESS)
    args = parser.parse_args()
    args.repo = args.repo.resolve()
    if args.job:
        spec = json.load(sys.stdin) if args.job=="-" else json.loads(args.job)
        try:
            result = execute(spec,args.repo)
        except FileNotFoundError as error:
            result = missing_source_result(args.repo,error)
            if result is None:
                raise
        print(packed(result))
        return
    if not args.out or not (3<=args.start_resolution<=args.max_resolution<=10 and
                            1<=args.workers<=32 and 1<=args.max_dsl_size<=8 and
                            args.timeout>0 and args.watch_seconds>=0 and args.sol_timeout>0 and
                            args.stall_rounds>=0 and args.scenario_timeout>0 and
                            args.scenario_budget_usd>0 and args.scenario_agent_max>=0 and
                            args.scenario_research_max>=0 and 1<=args.scenario_research_rounds<=6 and
                            args.scenario_agent_timeout>0 and args.scenario_agent_budget_usd>0 and
                            (not args.sol_proposals or args.sol_timeout<args.timeout)):
        parser.error("output and valid bounded worker/resolution parameters required")
    domain = (args.domain or args.repo/NATIVE).resolve()
    if args.domain and not domain.is_file():
        parser.error("the supplied original domain does not exist")
    out = args.out.resolve()
    out.mkdir(parents=True,exist_ok=True)
    with (out/".supervisor").open("a") as lock:
        try:
            fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:
            parser.error("campaign already has a supervisor")
        db = out/"queue.sqlite3"
        setup(db)
        cursor = out/"campaign.json"
        previous = json.loads(cursor.read_text()) if cursor.exists() else None
        epoch = 0 if previous is None else previous["epoch"]
        scenario_archive = []
        spent_archive = 0.0
        from original_discovery import inventory
        while True:
            context = current_context(args.repo)
            context["watched_sources"] = watched_sources(args.repo,domain)
            context["discovered_original_sources"] = inventory(args.repo,domain)
            extras = extra_specs(args.cases)
            observations = finite_observations(args.finite_observations)
            proposals = json.loads(args.proposal_file.read_text()) if args.proposal_file else []
            if not isinstance(proposals,list):
                parser.error("proposal file must contain a JSON array of OriginDSL terms")
            from language import Term
            proposals = [Term.parse(term).data() for term in proposals]
            context_key = token([context["atlas_program_sha256"],source_input_key(domain,args.repo),extras,
                                 context["plan_revision"],ENGINE,
                                 context["watched_sources"],context["discovered_original_sources"],
                                 context["git_head"],str(args.repo),str(domain),
                                 args.start_resolution,args.lean_toolchain,
                                 args.max_dsl_size,observations,proposals,args.sol_proposals,args.sol_timeout,args.timeout,
                                 args.scenario_model,args.scenario_proposals])
            if previous and previous["context"] != context_key:
                epoch += 1
                archives = previous.get("feedback",{}).get("excluded",[]) if previous.get("engine")==ENGINE else []
                # Tested scenarios stay known across epochs; the Atlas revision may have moved.
                scenario_archive = previous.get("scenario_history",[])
                spent_archive = previous.get("scenario_spent_usd",0.0)
                previous = None
            else:
                archives = []
            if previous and previous.get("feedback",{}).get("physical_solution"):
                return
            if (args.scenario_model and previous
                    and previous.get("scenario_spent_usd",0.0)>=args.scenario_total_budget_usd):
                print(f"scenario budget spent: {previous['scenario_spent_usd']:.2f} USD; see {out/'scenarios.json'}",
                      file=sys.stderr)
                sys.exit(4)
            # A proposing model changes the search every stage, so only its budget ends the run.
            if (args.until_proved and args.stall_rounds and not args.scenario_model and previous
                    and previous.get("stalled",0)>=args.stall_rounds):
                # Further stages repeat the same blocked state. Wait for a
                # changed source, contract or engine, which opens a new epoch.
                if not args.watch_seconds:
                    print(f"stalled after {previous['stalled']} unchanged stages; see {out/'work-order.json'}",
                          file=sys.stderr)
                    sys.exit(3)
                time.sleep(args.watch_seconds)
                continue
            start = previous["resolution"]+1 if previous else args.start_resolution
            for resolution in range(start,start+1 if args.until_proved else args.max_resolution+1):
                focus = previous.get("focus",[]) if previous and resolution>args.start_resolution+1 else []
                feedback_input = previous.get("feedback",{}) if previous else {"excluded":archives}
                dsl = dsl_spec(args,resolution,feedback_input,observations,proposals)
                history = (previous.get("signature"),previous.get("stalled",0)) if previous else (None,0)
                round_key,next_focus,feedback,progress = search_round(
                    args,db,out,epoch,resolution,context,domain,focus,extras,dsl,history,
                    previous.get("scenario_history",[]) if previous else scenario_archive,
                    previous.get("scenario_spent_usd",0.0) if previous else spent_archive)
                previous = {"context":context_key,"epoch":epoch,"resolution":resolution,
                            "round":round_key,"focus":next_focus,"feedback":feedback,"engine":ENGINE,
                            **progress}
                atomic_json(cursor,previous)
                if feedback.get("physical_solution"):
                    return
            if args.until_proved:
                continue
            if not args.watch_seconds:
                break
            time.sleep(args.watch_seconds)


if __name__ == "__main__":
    main()
