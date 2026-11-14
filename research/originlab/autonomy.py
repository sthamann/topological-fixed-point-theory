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
import sqlite3
import subprocess
import sys
import time

from main import atomic_json, current_context, finite_observations
from scenarios import DYNAMICS, TOPOLOGIES, discriminate
from source import NATIVE, source_input_key
from verify import digest

HERE = Path(__file__).resolve().parent
ENGINE = {p.name:digest(p) for p in sorted(HERE.glob("*.py"))}
ENGINE.update({name:digest(HERE/name) for name in ("gates.json","R102.lean")})

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


def worker(path,phase,root,timeout):
    with closing(connection(path)) as db:
        while (row := claim(db,phase)):
            job,spec,attempt = row
            try:
                child = subprocess.run([sys.executable,"-B",__file__,"--job","-","--repo",str(root)],
                                       input=spec,text=True,capture_output=True,timeout=timeout,check=True)
                result = json.loads(child.stdout)
                status = "done"
            except (subprocess.SubprocessError,OSError,ValueError) as error:
                stderr = getattr(error,"stderr",None)
                result = {"status":"ERROR","error_type":type(error).__name__,
                          "error":str(error) if isinstance(error,ValueError) else None,
                          "detail":str(stderr)[-3000:]}
                status = "queued" if attempt < 2 else "failed"
            db.execute("UPDATE jobs SET status=?,result=? WHERE id=?",(status,packed(result),job))


def result_rows(path,phase):
    with closing(connection(path)) as db:
        return [{"spec":json.loads(spec),"result":json.loads(data),"queue_status":status}
                for spec,data,status in db.execute(
                    "SELECT spec,result,status FROM jobs WHERE phase=? AND result IS NOT NULL ORDER BY priority DESC,id",(phase,))]


def execute(spec,root):
    lane = spec["lane"]
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


def specifications(resolution,domain,toolchain,focus,seed_research=True):
    yield 100,{"lane":"r113","domain":str(domain)}
    yield 100,{"lane":"r113_synthesis","domain":str(domain)}
    yield 100,{"lane":"source_bridge"}
    yield 100,{"lane":"original_affine_time"}
    yield 100,{"lane":"original_folded_time"}
    if seed_research:
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
    yield 80,{"lane":"synthesis","limit":min(32,8*resolution),"split_depth":1}
    yield 70,{"lane":"rewrite"}
    if toolchain:
        yield 65,{"lane":"lean","toolchain":toolchain}
    for grid in (24,48,96):
        yield 60,{"lane":"regulator","grid":grid,"cutoff":4,"time":"1/10"}
    pool = focus or product(TOPOLOGIES,DYNAMICS,("unit","alternating"))
    for topology,dynamics,metric in pool:
        yield 40,{"lane":"scenario","topology":topology,"dynamics":dynamics,
                  "metric":metric,"resolution":resolution}


def extra_specs(path):
    if path is None:
        return []
    rows = json.loads(path.read_text())
    if not isinstance(rows,list) or any(not isinstance(row,dict) for row in rows):
        raise ValueError("case manifest must be a list of job objects")
    for row in rows:
        lane = row.get("lane")
        if not isinstance(lane,str) or lane not in CASE_FIELDS:
            raise ValueError(f"unknown custom job: {lane!r}")
        required,optional = CASE_FIELDS[lane]
        fields = set(row)-{"lane"}
        if not required<=fields<=required|optional:
            raise ValueError(f"invalid fields for {lane}: {sorted(fields)}")
        if lane=="origin_dsl" and (type(row["size"]) is not int or not 1<=row["size"]<=8 or
                type(row.get("depth",5)) is not int or not 1<=row.get("depth",5)<=8):
            raise ValueError("OriginDSL size and test depth must be between 1 and 8")
    return list({packed(row):row for row in rows}.values())


def watched_sources(root, domain):
    from affine_time import SOURCES
    from folded_time import SOURCES as FOLDED_SOURCES
    from source_return import ORIGINAL_FILES
    from r113_synthesis import NATIVE_GRADE_CHECK
    from source_bridge import INPUTS
    paths = {root/path for path in SOURCE_FILES}
    paths.update(root/path for path in (*ORIGINAL_FILES,*SOURCES,*FOLDED_SOURCES,NATIVE_GRADE_CHECK,NATIVE,*INPUTS))
    for record in json.loads(domain.read_bytes()).get("source_originals",[]):
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
    for row in rows:
        spec,result = row["spec"],row["result"]
        if row["queue_status"]!="done":
            continue
        if spec["lane"]=="return_generator_search" and (probe := result.get("next_discriminating_observation")):
            jobs.append((92,{"lane":"return_generator_followup","input":{
                "gap":result["gap"],"signs":probe["signs"],"witnesses":result["candidates"]}}))
    original = next((r["result"] for r in rows if r["spec"]["lane"]=="r113"),{})
    anchor = next((r["result"] for r in rows if r["spec"]["lane"]=="original_affine_return"),{})
    if original.get("status")=="BLOCKED" and anchor.get("status")=="EXACT_ORIGINAL_AFFINE_RETURN":
        jobs.append((100,{"lane":"original_return_transport","domain":str(domain)}))
    return jobs


def search_round(args,db,out,epoch,resolution,context,domain,focus,extras,dsl):
    round_key = token([epoch,resolution,context["atlas_program_sha256"],context["git_head"],source_input_key(domain,args.repo),ENGINE,
                       context["watched_sources"],focus,args.lean_toolchain,extras,dsl])
    round_dir = out/round_key
    round_dir.mkdir(exist_ok=True)
    specs = list(specifications(resolution,domain,args.lean_toolchain,focus,
                                seed_research=resolution==args.start_resolution))
    enqueue(db,round_key,specs)
    enqueue(db,round_key,[(75,spec) for spec in extras])
    enqueue(db,round_key,[(55,dsl)])
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        def run_queued():
            futures = [pool.submit(worker,db,round_key,args.repo,args.timeout) for _ in range(args.workers)]
            for future in futures:
                future.result()
            return result_rows(db,round_key)

        rows = run_queued()
        followups = research_followups(rows,domain)
        if followups:
            enqueue(db,round_key,followups)
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
    frontier = {
        "state":{"BLOCKED":"SOURCE_REQUIRED","EXACT_SUPPLIED_CONTRACTION":"OPEN",
                 "ORIGINAL_NATIVE_VALUES_REQUIRE_REVIEW":"SOURCE_REVIEW_REQUIRED"}.get(original.get("status"),"ERROR"),
        "atlas_revision":context["atlas_revision"],
        "physical_solution":False,
        "physical_evidence":assess(original,results.get("r113_synthesis")),
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
    errors = [row for row in rows if row["queue_status"]=="failed"]
    report = {"epoch":epoch,"atlas_revision":context["atlas_revision"],
              "git_head":context["git_head"],"domain_sha256":digest(domain),"engine_sha256":token(ENGINE),
              "created_at":datetime.now(timezone.utc).isoformat(),
              "watched_sources_sha256":token(context["watched_sources"]),
              "original_status":original.get("status","ERROR"),"scenarios":len(scenarios),
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
              "failures":len(errors),"best_readout":questions[0] if questions else None,
              "physical_solution":False,"round":round_key,"resolution":resolution,
              "tool_versions":installed_versions(("sympy","cvc5","egglog","python-flint"))}
    for name,data in (("results",rows),("frontier",frontier),("summary",report)):
        atomic_json(round_dir/f"{name}.json",data)
    atomic_json(round_dir/"engine.json",ENGINE)
    atomic_json(round_dir/"sources.json",context["watched_sources"])
    atomic_json(out/"latest.json",report)
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
                "directed_research_round":new_research or prior_research}
    return round_key,refinement_focus(scenarios),feedback


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
    parser.add_argument("--job",help=argparse.SUPPRESS)
    args = parser.parse_args()
    args.repo = args.repo.resolve()
    if args.job:
        spec = json.load(sys.stdin) if args.job=="-" else json.loads(args.job)
        print(packed(execute(spec,args.repo)))
        return
    if not args.out or not (3<=args.start_resolution<=args.max_resolution<=10 and
                            1<=args.workers<=32 and 1<=args.max_dsl_size<=8 and
                            args.timeout>0 and args.watch_seconds>=0 and args.sol_timeout>0 and
                            (not args.sol_proposals or args.sol_timeout<args.timeout)):
        parser.error("output and valid bounded worker/resolution parameters required")
    domain = (args.domain or args.repo/NATIVE).resolve()
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
        while True:
            context = current_context(args.repo)
            context["watched_sources"] = watched_sources(args.repo,domain)
            extras = extra_specs(args.cases)
            observations = finite_observations(args.finite_observations)
            proposals = json.loads(args.proposal_file.read_text()) if args.proposal_file else []
            if not isinstance(proposals,list):
                parser.error("proposal file must contain a JSON array of OriginDSL terms")
            from language import Term
            proposals = [Term.parse(term).data() for term in proposals]
            context_key = token([context["atlas_program_sha256"],source_input_key(domain,args.repo),extras,
                                 context["plan_revision"],ENGINE,
                                 context["watched_sources"],
                                 context["git_head"],str(args.repo),str(domain),
                                 args.start_resolution,args.max_resolution,args.lean_toolchain,
                                 args.max_dsl_size,observations,proposals,args.sol_proposals,args.sol_timeout])
            if previous and previous["context"] != context_key:
                epoch += 1
                archives = previous.get("feedback",{}).get("excluded",[]) if previous.get("engine")==ENGINE else []
                previous = None
            else:
                archives = []
            start = previous["resolution"]+1 if previous else args.start_resolution
            for resolution in range(start,args.max_resolution+1):
                focus = previous.get("focus",[]) if previous and resolution>args.start_resolution+1 else []
                feedback_input = previous.get("feedback",{}) if previous else {"excluded":archives}
                dsl = dsl_spec(args,resolution,feedback_input,observations,proposals)
                round_key,next_focus,feedback = search_round(args,db,out,epoch,resolution,context,domain,focus,extras,dsl)
                previous = {"context":context_key,"epoch":epoch,"resolution":resolution,
                            "round":round_key,"focus":next_focus,"feedback":feedback,"engine":ENGINE}
                atomic_json(cursor,previous)
            if not args.watch_seconds:
                break
            time.sleep(args.watch_seconds)


if __name__ == "__main__":
    main()
