"""Executable acceptance benchmarks for OriginLab, with scoped results."""
from pathlib import Path
from tempfile import TemporaryDirectory
import json
from shutil import copyfile
import subprocess
import sys

from affine_time import marked_charges
from backends import inverse_rank, lean_check
from coherence import r104
from marked_search import cell_search
from root_doublet import radial_symmetry_check
from selection import tetrahedron
from source import MARKS, NATIVE, evaluate_original
from verify import car_replay, digest, r102

OPERATOR_ENGINE = Path("experiments/theory-contracts/seam-constructive-search-20261015/operator_synthesis")
COLOR_DIR = Path("experiments/theory-contracts/tfpt-unified-source-derivation-20261029/aps_source_review")


def cegis(root):
    engine = root/OPERATOR_ENGINE/"clifford_search.py"
    with TemporaryDirectory() as directory:
        subprocess.run([sys.executable,"-B",str(engine),"--run-dir",directory,"--seed","7"],
                       cwd=root,check=True,capture_output=True,timeout=60)
        receipt = json.loads((Path(directory)/"results.json").read_text())
    return receipt


def benchmark_suite(root,toolchain=None):
    results = {}

    def record(identifier, calculation):
        try:
            passed,details = calculation()
            results[identifier] = {"detail":details,"status":"PASS" if passed else "FAIL"}
        except (ValueError,ArithmeticError,OSError,subprocess.SubprocessError) as error:
            results[identifier] = {"status":"ERROR","error":f"{type(error).__name__}: {error}"}

    def b01():
        r = cegis(root)
        trace = r["counterexample_trace"]
        return (r["candidate_count"]==2048 and len(r["survivors"])==1 and
                r["survivors"][0]["word"]==[0,2,4,6,8] and r["survivors"][0]["antilinear"] and bool(trace),
                {"generated":r["candidate_count"],"survivors":r["survivors"],
                 "learned_counterexamples":trace,"engine_sha256":r["engine_sha256"]})

    def b02():
        r = car_replay(root,5)
        q = r["discriminator"]
        return r["shared_moments"]==2245 and q["shortest_length"]==5 and bool(q["questions"]),r

    def b03():
        r = r102(3,96)
        inverse = inverse_rank(3,96)
        proof = lean_check(Path(__file__).with_name("R102.lean"),toolchain)
        values = sorted({row["n"] for row in r["solutions"]})
        passed = (values==[4] and
                  (inverse["status"]=="UNAVAILABLE" or (inverse["status"]=="EXHAUSTED" and inverse["solutions"]==values)) and
                  (not toolchain or proof["status"]=="KERNEL_CHECKED"))
        return passed,{"selected_ranks":values,"inverse":inverse,"lean":proof,
                       "scope":"declared integral discriminant class"}

    def b04():
        r = r104(root)
        return (r["markov"]["all_nonempty_mode_tripartitions"]=="REFUTED_SCOPED" and
                r["cap"]["state_bond_dimension"]==2 and r["coherent_gram"]["first_phase_sensitive_gram_depth"]==8),r

    def b05():
        r = radial_symmetry_check(root)
        passed = (r["affine_two_time_control"]["mixed_ratio"]=="81/85" and
                  all(not row["commutes_with_middle_on_root_doublet"] for row in r["quarter_cases"]))
        return passed,r

    def b06():
        checker = root/COLOR_DIR/"check_original_color_sectors.py"
        with TemporaryDirectory() as directory:
            scratch = Path(directory)
            copy = scratch/COLOR_DIR/checker.name
            copy.parent.mkdir(parents=True)
            copyfile(checker,copy)
            (scratch/"verification").symlink_to(root/"verification",target_is_directory=True)
            payload = subprocess.run([sys.executable,"-B",str(copy)],cwd=scratch,
                                     check=True,capture_output=True,timeout=60,text=True)
        receipt = json.loads(payload.stdout.splitlines()[-1])
        grades = marked_charges()
        return (receipt["verdict"]=="SCOPED_ORIGINAL_COLOR_CENTER_FIXED_SOURCE_AND_CHARGED_FIELD_WORDS" and
                [r["color_triality"] for r in grades]==[0,1,2],
                {"original_checker":str(checker.relative_to(root)),"checker_sha256":digest(checker),
                 "receipt":receipt,"marked_trialities":[r["color_triality"] for r in grades]})

    def b07():
        r = evaluate_original(root,root/NATIVE)
        native = r.get("native_endpoint_readout",{})
        return (r["status"]=="BLOCKED" and r.get("original_values")==dict.fromkeys(MARKS) and
                native.get("whole_status")=="original_full_source_operation_not_supplied"),{
                    "physical_status":r["status"],"detail":r}

    def b08():
        r = cell_search({})
        return (r["status"]=="EXACT_FINITE_SEWING_SEARCH" and r["shared"]["g_identified"] is False and
                r["next_probe"]["separated_pairs"]==6),r

    def b09():
        from evidence_gates import assess
        example = {"status":"EXACT_SUPPLIED_CONTRACTION","input_role":"arithmetic_control_only",
                   "values":dict.fromkeys(MARKS,"1"),"physical_source_identity_proved":False}
        r = assess(example)
        return (r["physical_status"]=="OPEN" and all(g["status"]!="PASS" for g in r["gates"].values())),r

    for number,check in enumerate((b01,b02,b03,b04,b05,b06,b07,b08,b09),1):
        record(f"B{number:02}",check)
    return {"status":"BENCHMARKS_PASSED" if all(r["status"]=="PASS" for r in results.values())
            else "BENCHMARKS_INCOMPLETE","benchmarks":results,"physical_solution":False}
