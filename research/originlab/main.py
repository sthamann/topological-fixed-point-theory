"""Run the independent OriginLab research campaigns."""
from pathlib import Path
from datetime import datetime, timezone
from importlib.metadata import version
from itertools import count
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
from source import NATIVE, evaluate_original, missing_source_result, source_input_key
from backends import inverse_rank, lean_check
from verify import car_replay, clock_test, digest, eigenvalue_test, r102

HERE = Path(__file__).resolve().parent
CAMPAIGNS = ("r113", "r113_synthesis", "source_bridge", "marked_search", "double_cell", "moving_carrier", "affine_polarization", "affine_operators", "affine_car", "colour_ward", "weak_source", "weak_folded_clock", "p1_field_search", "weak_return_rank", "four_mark_geometry", "f1_family_index", "r102", "bell", "clock", "identifiability", "four_mark_selection",
             "spectral_source_controls", "empirical_source_audit", "composed_search", "seam_cohomology", "p1_source_falsifiers", "p1_charged_modulus", "p1_four_current_gate", "p1_e8_geometric_lift", "p1_seam_dirac", "p1_spinor_cover", "p1_geometric_gram", "p1_minimal_spinor_marks", "p1_principal_su3", "p1_principal_su3_ward", "p1_moving_divisor_origin", "p1_operator_closure", "p1_modular_source_no_go", "p1_modulus_score", "p1_original_j_witness", "p1_elliptic_determinant_source", "f1_source_theorem_audit", "f1_functorial_source", "f1_hidden_history_dilation", "f1_covariant_quantum_semigroup", "f1_spin2_stieltjes_gate", "f1_single_source_clifford_bridge", "f1_ramond_e8_grade_bridge", "f1_qutrit_car_e8_bridge", "f1_e8_gravitational_inflow", "f1_p1_chiral_index_gate", "f1_single_original_source_gate", "f1_pregeometric_history_compiler", "f1_compiler_kms_no_go", "f1_history_gram_universal_no_go", "f1_oriented_seam_source", "f1_p1_physical_necessity_routes", "f1_seam_projective_consistency", "f1_p1_cocycle_boundary_gate", "f1_klein_twisted_source", "f1_klein_phase_space_double", "f1_origin_scenarios", "p0", "f1")
FINITE_CONTROLS = ("double_cell", "moving_carrier", "affine_polarization", "affine_operators", "affine_car", "colour_ward", "weak_source", "weak_return_rank", "four_mark_geometry", "f1_family_index", "seam_cohomology", "p1_source_falsifiers", "p1_charged_modulus", "p1_four_current_gate", "p1_e8_geometric_lift", "p1_seam_dirac", "p1_spinor_cover", "p1_geometric_gram", "p1_minimal_spinor_marks", "p1_principal_su3", "p1_principal_su3_ward", "p1_moving_divisor_origin", "p1_operator_closure", "p1_modular_source_no_go", "p1_modulus_score", "p1_original_j_witness", "p1_elliptic_determinant_source", "f1_source_theorem_audit", "f1_functorial_source", "f1_hidden_history_dilation", "f1_covariant_quantum_semigroup", "f1_spin2_stieltjes_gate", "f1_single_source_clifford_bridge", "f1_ramond_e8_grade_bridge", "f1_qutrit_car_e8_bridge", "f1_e8_gravitational_inflow", "f1_p1_chiral_index_gate", "f1_single_original_source_gate", "f1_pregeometric_history_compiler", "f1_compiler_kms_no_go", "f1_history_gram_universal_no_go", "f1_oriented_seam_source", "f1_p1_physical_necessity_routes", "f1_seam_projective_consistency", "f1_p1_cocycle_boundary_gate", "f1_klein_twisted_source", "f1_klein_phase_space_double", "f1_origin_scenarios")
FINITE_CONTROLS += ("f1_prephysical_interpretation",)
CAMPAIGNS += ("f1_prephysical_interpretation",)
FINITE_CONTROLS += ("q0_rule_search", "p1_origin_hypotheses")
CAMPAIGNS += ("q0_rule_search", "p1_origin_hypotheses")
CAMPAIGNS += ("proof_search",)
FINITE_CONTROLS += ("f1_spherical_pairing_selection",)
CAMPAIGNS += ("f1_spherical_pairing_selection",)
FINITE_CONTROLS += ("f1_binary_clifford_origin",)
CAMPAIGNS += ("f1_binary_clifford_origin",)
CAMPAIGNS += ("f1_solution",)
FINITE_CONTROLS += ("f1_native_pair_selection",)
CAMPAIGNS += ("f1_native_pair_selection",)
FINITE_CONTROLS += ("f1_reflected_cap",)
CAMPAIGNS += ("f1_reflected_cap",)
FINITE_CONTROLS += ("r113_pair_search",)
CAMPAIGNS += ("r113_pair_search",)
FINITE_CONTROLS += ("f1_two_arc_purity",)
CAMPAIGNS += ("f1_two_arc_purity",)


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
           "problem:originlab-affine-neutral-return-20261114",
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
    plan = json.loads(progress.read_text()) if progress.is_file() else {}
    atlas = root/"semantic_explorer/research/program.json"
    raw = atlas.read_bytes() if atlas.is_file() else None
    data = json.loads(raw) if raw is not None else {"problems":[]}
    return dict(plan_revision=plan.get("revision"),next_task=plan.get("next_task"),
                plan_recorded_atlas_revision=plan.get("current_atlas_revision"),
                atlas_program_sha256=sha256(raw).hexdigest() if raw is not None else None,
                atlas_revision=data.get("revision"),
                atlas_provenance="canonical_local" if raw is not None else "local_atlas_unavailable",
                atlas_inputs=atlas_inputs(data),
                git_head=subprocess.check_output(["git","-C",str(root),"rev-parse","HEAD"],
                                                 text=True).strip(),
                python=platform.python_version(),sympy=version("sympy"))


def campaign_context(args):
    from source_bridge import INPUTS

    research_options = {}
    if "q0_rule_search" in args.campaign or ("f1_solution" in args.campaign and args.f1_contract is None):
        research_options["q0_depth"] = args.q0_depth
    if args.p1_hypothesis_input:
        research_options["p1_hypothesis_input_sha256"] = digest(args.p1_hypothesis_input)
    if args.proof_contract:
        from proof_search import read_contract
        _, modules = read_contract(args.proof_contract)
        research_options["proof_contract"] = {
            "manifest":digest(args.proof_contract),
            "modules":{file.name:digest(file) for file,_ in modules},
        }
    if "f1_solution" in args.campaign:
        from f1_solution import contract_fingerprints
        research_options["f1_solution_contract"] = contract_fingerprints(args.f1_contract,args.f1_contract_sha256)
        if args.f1_contract is None:
            from r113_radial_order import RADIAL_DATA
            radial = args.repo/RADIAL_DATA
            research_options["r113_radial_order_data_sha256"] = digest(radial) if radial.is_file() else None
    if "f1_native_pair_selection" in args.campaign or ("f1_solution" in args.campaign and args.f1_contract is None):
        from f1_native_pair_selection import RING_SOURCE
        source = args.repo/RING_SOURCE
        research_options["native_pair_ring_sha256"] = digest(source) if source.is_file() else None
    if {"f1_reflected_cap", "r113_pair_search"}.intersection(args.campaign) or ("f1_solution" in args.campaign and args.f1_contract is None):
        from f1_reflected_cap import ASSEMBLER
        research_options["reflected_cap_inputs"] = {
            str(path): digest(args.repo/path) if (args.repo/path).is_file() else None
            for path in (ASSEMBLER, NATIVE)
        }
    if args.finite_controls_only:
        return dict(plan_revision=None,atlas_revision=None,atlas_provenance="not_used_for_finite_controls",
                    git_head=subprocess.check_output(["git","-C",str(args.repo),"rev-parse","HEAD"],
                                                     text=True).strip(),
                    python=platform.python_version(),sympy=version("sympy"), **research_options)
    context = current_context(args.repo)
    context.update(research_options)
    domain = args.domain or args.repo/NATIVE
    context["domain_sha256"] = digest(domain) if domain.is_file() else None
    context["source_input_key"] = source_input_key(domain,args.repo)
    if "source_bridge" in args.campaign:
        context["source_bridge_inputs"] = {
            str(path): digest(args.repo/path) if (args.repo/path).is_file() else None
            for path in (NATIVE, *INPUTS)
        }
    if "r113_synthesis" in args.campaign and domain.resolve()==(args.repo/NATIVE).resolve():
        from native_weak_bridge import INPUTS as weak_inputs
        from native_weak_w import INPUTS as w_inputs

        context["native_weak_source_inputs"] = {
            str(path):digest(args.repo/path) if (args.repo/path).is_file() else None
            for path in (*weak_inputs,*w_inputs)
        }
    if "weak_folded_clock" in args.campaign:
        from folded_time import SOURCES as folded_sources

        context["weak_folded_clock_inputs"] = {
            str(path):digest(args.repo/path) if (args.repo/path).is_file() else None for path in folded_sources
        }
    if "p1_field_search" in args.campaign and not args.p1_field_input:
        from p1_field_search import INPUTS as p1_sources

        context["p1_field_sources"] = {
            str(path):digest(args.repo/path) if (args.repo/path).is_file() else None for path in p1_sources
        }
    if args.composed_square_resolution is not None:
        context["composed_square_resolution"] = args.composed_square_resolution
    context["extra_inputs"] = {str(p.resolve()):digest(p)
                               for p in (args.proposal_file,args.finite_observations,args.r113_synthesis_input,
                                         args.marked_search_input,args.p1_field_input,
                                         getattr(args,"composed_search_input",None)) if p}
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
    if campaign == "composed_search":
        from composed_search import run as search_composed, square_control
        return search_composed(json.loads(args.composed_search_input.read_text()) if args.composed_search_input
                               else square_control(args.composed_square_resolution or 8))
    if campaign == "seam_cohomology":
        from seam_cohomology import controls
        return controls()
    if campaign == "p1_source_falsifiers":
        from p1_source_falsifiers import controls
        return controls()
    if campaign == "p1_charged_modulus":
        from p1_charged_modulus import controls
        return controls()
    if campaign == "p1_four_current_gate":
        from p1_four_current_gate import controls
        return controls()
    if campaign == "p1_e8_geometric_lift":
        from p1_e8_geometric_lift import controls
        return controls()
    if campaign == "p1_seam_dirac":
        from p1_seam_dirac import controls
        return controls()
    if campaign == "p1_spinor_cover":
        from p1_spinor_cover import controls
        return controls()
    if campaign == "p1_geometric_gram":
        from p1_geometric_gram import controls
        return controls()
    if campaign == "p1_minimal_spinor_marks":
        from p1_minimal_spinor_marks import controls
        return controls()
    if campaign == "p1_principal_su3":
        from p1_principal_su3 import controls
        return controls()
    if campaign == "p1_principal_su3_ward":
        from p1_principal_su3_ward import controls
        return controls()
    if campaign == "p1_moving_divisor_origin":
        from p1_moving_divisor_origin import controls
        return controls()
    if campaign == "p1_operator_closure":
        from p1_operator_closure import controls
        return controls()
    if campaign == "p1_modular_source_no_go":
        from p1_modular_source_no_go import controls
        return controls()
    if campaign == "p1_modulus_score":
        from p1_modulus_score import controls
        return controls()
    if campaign == "p1_original_j_witness":
        from p1_original_j_witness import controls
        return controls()
    if campaign == "p1_elliptic_determinant_source":
        from p1_elliptic_determinant_source import controls
        return controls()
    if campaign == "f1_source_theorem_audit":
        from f1_source_theorem_audit import controls
        return controls()
    if campaign == "f1_functorial_source":
        from f1_functorial_source import controls
        return controls()
    if campaign == "f1_hidden_history_dilation":
        from f1_hidden_history_dilation import controls
        return controls()
    if campaign == "f1_covariant_quantum_semigroup":
        from f1_covariant_quantum_semigroup import controls
        return controls()
    if campaign == "f1_spin2_stieltjes_gate":
        from f1_spin2_stieltjes_gate import controls
        return controls()
    if campaign == "f1_single_source_clifford_bridge":
        from f1_single_source_clifford_bridge import controls
        return controls()
    if campaign == "f1_ramond_e8_grade_bridge":
        from f1_ramond_e8_grade_bridge import controls
        return controls()
    if campaign == "f1_qutrit_car_e8_bridge":
        from f1_qutrit_car_e8_bridge import controls
        return controls()
    if campaign == "f1_e8_gravitational_inflow":
        from f1_e8_gravitational_inflow import controls
        return controls()
    if campaign == "f1_p1_chiral_index_gate":
        from f1_p1_chiral_index_gate import controls
        return controls()
    if campaign == "f1_single_original_source_gate":
        from f1_single_original_source_gate import controls
        return controls()
    if campaign == "f1_pregeometric_history_compiler":
        from f1_pregeometric_history_compiler import controls
        return controls()
    if campaign == "f1_compiler_kms_no_go":
        from f1_compiler_kms_no_go import controls
        return controls()
    if campaign == "f1_history_gram_universal_no_go":
        from f1_history_gram_universal_no_go import controls
        return controls()
    if campaign == "f1_oriented_seam_source":
        from f1_oriented_seam_source import controls
        return controls()
    if campaign == "f1_p1_physical_necessity_routes":
        from f1_p1_physical_necessity_routes import controls
        return controls()
    if campaign == "f1_seam_projective_consistency":
        from f1_seam_projective_consistency import controls
        return controls()
    if campaign == "f1_p1_cocycle_boundary_gate":
        from f1_p1_cocycle_boundary_gate import controls
        return controls()
    if campaign == "f1_klein_twisted_source":
        from f1_klein_twisted_source import controls
        return controls()
    if campaign == "f1_origin_scenarios":
        from f1_origin_scenarios import controls
        return controls()
    if campaign == "f1_klein_phase_space_double":
        from f1_klein_phase_space_double import controls
        return controls()
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
    if campaign == "weak_source":
        from weak_source import affine_controls
        return affine_controls()
    if campaign == "p1_field_search":
        from p1_field_search import controls, from_input
        return from_input(json.loads(args.p1_field_input.read_text())) if args.p1_field_input else controls(root)
    if campaign in ("weak_return_rank","weak_folded_clock","four_mark_geometry","colour_ward"):
        from autonomy import execute as execute_job
        return execute_job({"lane":campaign},root)
    if campaign == "f1_family_index":
        from f1_family_index import controls
        return controls()
    if campaign == "f1_prephysical_interpretation":
        from f1_prephysical_interpretation import controls
        return controls()
    if campaign == "f1_spherical_pairing_selection":
        from f1_spherical_pairing_selection import controls
        return controls()
    if campaign == "f1_binary_clifford_origin":
        from f1_binary_clifford_origin import controls
        return controls()
    if campaign == "f1_native_pair_selection":
        from f1_native_pair_selection import controls
        return controls(root)
    if campaign == "f1_reflected_cap":
        from f1_reflected_cap import controls
        return controls(root)
    if campaign == "f1_two_arc_purity":
        from f1_two_arc_purity import controls
        return controls()
    if campaign == "r113_pair_search":
        from r113_pair_search import search
        return search(root)
    if campaign == "q0_rule_search":
        from q0_rule_search import run as q0_search
        source_root = None if args.finite_controls_only else root
        return q0_search(source_root, args.domain, args.q0_depth)
    if campaign == "p1_origin_hypotheses":
        from p1_origin_hypotheses import search
        options = json.loads(args.p1_hypothesis_input.read_text()) if args.p1_hypothesis_input else {}
        return search(**options)
    if campaign == "proof_search":
        from proof_search import DEFAULT_TOOLCHAIN, run as proof_control, search_contract
        options = {"first_cost":getattr(args,"proof_cost",args.proof_first_cost),
                   "shells":args.proof_shells,"toolchain":args.lean_toolchain or DEFAULT_TOOLCHAIN}
        return search_contract(args.proof_contract,**options) if args.proof_contract else proof_control(**options)
    if campaign == "f1_solution":
        from f1_solution import search
        from proof_search import DEFAULT_TOOLCHAIN
        hypotheses = json.loads(args.p1_hypothesis_input.read_text()) if args.p1_hypothesis_input else None
        return search(args.f1_contract,args.f1_contract_sha256,args.f1_stage,
                      args.lean_toolchain or DEFAULT_TOOLCHAIN,args.f1_check_seconds, root=root,
                      domain=args.domain,q0_depth=args.q0_depth,p1_hypotheses=hypotheses)
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
            result = missing_source_result(args.repo,error) if isinstance(error,FileNotFoundError) else None
            if result is None:
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
            complete = outcomes.get("f1_solution",{})
            source_result = outcomes.get("r113") or complete.get("original_source")
            bell = outcomes.get("bell",{}).get("discriminator",{})
            proven = complete.get("status") in ("F1_PROVED","F1_CONTRACT_PROVED")
            evidence = assess(source_result,outcomes.get("r113_synthesis"),complete)
            return dict(status=complete["status"] if proven else "OPEN",
                        physical_solution=evidence["physical_solution"],
                        complete_contract=complete or None,
                        source_status=(source_result or {}).get("status","NOT_RUN"),
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
                        native_weak_higgs_bridge=outcomes.get("r113_synthesis",{}).get("native_weak_higgs_bridge"),
                        native_w_weak_words=outcomes.get("r113_synthesis",{}).get("native_w_weak_words"),
                        finite_marked_search_status=outcomes.get("marked_search",{}).get("status","NOT_RUN"),
                        finite_mark_transfer=outcomes.get("marked_search",{}).get("transfer"),
                        harmonic_seam_search=outcomes.get("marked_search",{}).get("harmonic_seam"),
                        composed_search=outcomes.get("composed_search"),
                        origin_scenarios=outcomes.get("f1_origin_scenarios"),
                        original_composed_search=(source_result or {}).get("composed_search"),
                        double_cell_control_status=outcomes.get("double_cell",{}).get("status","NOT_RUN"),
                        moving_carrier_control_status=outcomes.get("moving_carrier",{}).get("status","NOT_RUN"),
                        affine_polarization_status=(source_result or {}).get("affine_polarization",{}).get("status","NOT_RUN"),
                        affine_profile_returns=(source_result or {}).get("affine_profile_returns"),
                        weak_split_returns=(source_result or {}).get("weak_split_returns"),
                        weak_gauge_responses=(source_result or {}).get("weak_gauge_responses"),
                        weak_clock_sectors=(source_result or {}).get("weak_clock_sectors"),
                        p1_charged_jet=(source_result or {}).get("p1_charged_jet"),
                        p1_field_search=outcomes.get("p1_field_search"),
                        weak_folded_clock_status=outcomes.get("weak_folded_clock",{}).get("status","NOT_RUN"),
                        weak_parity_bridge=(source_result or {}).get("weak_parity_bridge"),
                        weak_return_grade_two=(source_result or {}).get("weak_return_grade_two"),
                        weak_source_status=outcomes.get("weak_source",{}).get("status","NOT_RUN"),
                        fixed_middle_search=outcomes.get("weak_source",{}).get("fixed_middle_search"),
                        weak_source_next_observation=outcomes.get("weak_source",{}).get("next_original_observation"),
                        joint_weak_return_support=outcomes.get("weak_return_rank",{}).get("completion_orbits",{}).get("support_by_grade"),
                        four_mark_geometry=outcomes.get("four_mark_geometry"),
                        affine_colour_ward=outcomes.get("colour_ward"),
                        q0_rule_search=outcomes.get("q0_rule_search"),
                        formal_proof_search=outcomes.get("proof_search"),
                        p1_origin_hypothesis_status=outcomes.get("p1_origin_hypotheses",{}).get("status","NOT_RUN"),
                        affine_operator_source_status=outcomes.get("affine_operators",{}).get("status","NOT_RUN"),
                        mark_alphabet_completion=outcomes.get("affine_operators",{}).get("next_original_observation"),
                        physical_evidence=evidence,
                        separator_requires_original_field_identification=True,
                        energy_moments=(source_result or {}).get("energy_moments"),
                        regional_motion=((source_result or {}).get("physical_reconstruction",{})
                                         .get("regional",{}).get("regional_motion")),
                        source_obligations=(source_result or {}).get("obligations",{}))
        record("f1",f1_frontier)
    p0 = outcomes.get("p0",{})
    complete = outcomes.get("f1_solution",{})
    evidence = assess(outcomes.get("r113") or complete.get("original_source"),
                      outcomes.get("r113_synthesis"),complete)
    next_action = ("Review the exact finite controls and their original-source follow-up."
                   if args.finite_controls_only else
                   "Produce the original marked endpoint packet or equivalent weighted "
                   "form before claiming physical transfer or regional time.")
    if complete.get("status") in ("F1_PROVED","F1_CONTRACT_PROVED"):
        next_action = "Replay the common-model proof against the pinned contract."
    summary = dict(round=index,created_at=datetime.now(timezone.utc).isoformat(),
                   context=manifest,gates=evidence["gates"],physical_solution=evidence["physical_solution"],
                   complete_contract_status=complete.get("status","NOT_RUN"),
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
    parser.add_argument("--composed-search-input",type=Path,help="Search an exact composition of finite operator fragments")
    parser.add_argument("--composed-square-resolution",type=int,help="Set the rational grid denominator in the built-in composed search")
    parser.add_argument("--p1-field-input",type=Path,help="Resolve a supplied two-face charged jet")
    parser.add_argument("--q0-depth",type=int,default=3,choices=(1,2,3,4))
    parser.add_argument("--p1-hypothesis-input",type=Path,help="Select registered hypotheses and cone bound from JSON")
    parser.add_argument("--proof-contract",type=Path,help="Trusted pinned Lean goal and term grammar")
    parser.add_argument("--proof-first-cost",type=int,default=2)
    parser.add_argument("--proof-shells",type=int,default=2)
    parser.add_argument("--f1-contract",type=Path,help="Complete formal goal and its eight gate predicates")
    parser.add_argument("--f1-contract-sha256",help="Externally fixed SHA256 of the complete contract")
    parser.add_argument("--f1-check-seconds",type=int,default=20)
    parser.add_argument("--until-proof",action="store_true",help="Advance the common model search until a proof is accepted")
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
    parser.add_argument("--rounds",type=int,default=1,
                        help="New rounds this invocation, including native rounds with --until-proof")
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
    if args.composed_search_input and "composed_search" not in args.campaign:
        parser.error("--composed-search-input requires the composed_search campaign")
    if args.composed_square_resolution is not None and (
            args.composed_square_resolution<1 or "composed_search" not in args.campaign or args.composed_search_input):
        parser.error("--composed-square-resolution requires a positive denominator and the built-in composed search")
    if args.p1_field_input and "p1_field_search" not in args.campaign:
        parser.error("--p1-field-input requires the p1_field_search campaign")
    if args.p1_hypothesis_input and not (
        "p1_origin_hypotheses" in args.campaign or ("f1_solution" in args.campaign and not args.f1_contract)
    ):
        parser.error("--p1-hypothesis-input requires p1_origin_hypotheses or native f1_solution")
    if args.proof_contract and "proof_search" not in args.campaign:
        parser.error("--proof-contract requires the proof_search campaign")
    if args.proof_first_cost < 2 or args.proof_shells < 1:
        parser.error("proof search starts at cost two and requires at least one shell")
    if (args.f1_contract is None) != (args.f1_contract_sha256 is None):
        parser.error("--f1-contract and --f1-contract-sha256 must be supplied together")
    if args.f1_contract and "f1_solution" not in args.campaign:
        parser.error("--f1-contract requires the f1_solution campaign")
    if args.until_proof and (args.watch_seconds or "f1_solution" not in args.campaign):
        parser.error("--until-proof requires f1_solution and is incompatible with --watch-seconds")
    if args.f1_check_seconds < 1:
        parser.error("--f1-check-seconds must be positive")
    if args.domain and not args.domain.is_file():
        parser.error("the supplied original domain does not exist")
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
                         for k,v in vars(args).items()
                         if k not in {"resume","rounds","out","jobs","watch_seconds","until_proof"}}
        if not args.finite_controls_only:
            configuration.pop("finite_controls_only")
        engine = {p.name:digest(p) for p in [*HERE.glob("*.py"), *HERE.glob("*.lean"),
                                            HERE/"gates.json", HERE/"f1_native_contract.json"]}
        if state and (configuration != state["configuration"] or engine != state["engine"]):
            parser.error("configuration or engine changed, start a fresh run")
        if args.until_proof and state and state.get("complete_contract_status") in {
            "F1_PROVED","F1_CONTRACT_PROVED",
        }:
            print(f"proof already accepted in round {prior}",flush=True)
            return 0
        had_error = False
        rounds = count(prior+1) if args.until_proof else range(prior+1,prior+args.rounds+1)
        for i in rounds:
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
            args.proof_cost = args.proof_first_cost + (i-1)*args.proof_shells
            args.f1_stage = i+1
            atomic_manifest = dict(context,seed=args.seed,configuration=configuration,engine=engine)
            summary,p0 = one_round(args,target,i,retained,prior_questions,atomic_manifest)
            retained = [x["term"] for x in p0.get("survivors",[])[:128]]
            prior_questions = p0.get("questions",[])[:3]
            had_error |= any(v in {"ERROR","BENCHMARKS_INCOMPLETE"} for v in summary["outcomes"].values())
            state = dict(round=i,survivors=retained,questions=prior_questions,configuration=configuration,
                         engine=engine,context=context,
                         complete_contract_status=summary["complete_contract_status"])
            atomic_json(state_path,state)
            print(f"saved {target/'summary.json'}",flush=True)
            verdict = summary["complete_contract_status"]
            if args.until_proof:
                if verdict in {"F1_PROVED","F1_CONTRACT_PROVED","PROOF_BACKEND_INCONCLUSIVE","ERROR"}:
                    break
                if verdict=="NATIVE_SEARCH_BOUND_REACHED" and i>=prior+args.rounds:
                    break
        return 1 if had_error else 0
    finally:
        lock.close()


if __name__ == "__main__":
    sys.exit(main())
