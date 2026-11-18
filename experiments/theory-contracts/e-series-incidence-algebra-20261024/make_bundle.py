#!/usr/bin/env python3
"""Export the original files and the new incidence proof without rewriting them."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import zipfile

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]

RUNNER = '''#!/usr/bin/env python3
import argparse
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
parser = argparse.ArgumentParser()
parser.add_argument("--optimized", action="store_true")
args = parser.parse_args()
manifest = json.loads((ROOT/"manifest.json").read_text())
for relative, digest in manifest.items():
    path = ROOT/relative
    if not path.is_file() or sha256(path.read_bytes()).hexdigest() != digest:
        raise RuntimeError("Missing or changed bundled source: " + relative)
tests = [
    ("compiler", "experiments/theory-contracts/e-series-joint-source-20261023/compiler_closure.py"),
    ("clock-cover", "experiments/theory-contracts/e-series-incidence-algebra-20261024/clock_cover_audit.py"),
    ("selection", "experiments/theory-contracts/e-series-joint-source-20261023/stabilizer_selection.py"),
    ("relative", "experiments/theory-contracts/e-series-joint-source-20261023/relative_transport.py"),
    ("incidence", "experiments/theory-contracts/e-series-incidence-algebra-20261024/checker.py"),
]
records = []
for mode in (["normal", "optimized"] if args.optimized else ["normal"]):
    out = ROOT/"check-results"/mode
    out.mkdir(parents=True, exist_ok=True)
    for name, relative in tests:
        cmd = [sys.executable, "-B"] + (["-OO"] if mode == "optimized" else []) + [str(ROOT/relative)]
        if name == "incidence":
            cmd += ["--output-dir", str(out)]
        run = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=False, timeout=300)
        (out/(name+".stdout")).write_bytes(run.stdout)
        if run.returncode:
            sys.stderr.buffer.write(run.stderr)
            raise RuntimeError(name + " failed in " + mode)
        if run.stderr:
            sys.stderr.buffer.write(run.stderr)
        records.append({"mode": mode, "checker": name, "exit_code": run.returncode,
                        "stdout_sha256": sha256(run.stdout).hexdigest()})
        print(mode + " " + name + ": PASS", flush=True)
if args.optimized:
    for filename in [name+".stdout" for name, _ in tests] + ["results.json", "geometry.json"]:
        normal = (ROOT/"check-results/normal"/filename).read_bytes()
        optimized = (ROOT/"check-results/optimized"/filename).read_bytes()
        if normal != optimized:
            raise RuntimeError("Normal/optimized output differs: " + filename)
    print("All five checker outputs and both geometry files are byte-identical.", flush=True)
receipt = {"verdict": "PASS", "manifest_files_checked": len(manifest),
           "normal_and_optimized": args.optimized, "runs": records}
(ROOT/"check-results/receipt.json").write_text(json.dumps(receipt, indent=2)+"\\n")
'''


def make_bundle(destination):
    paths = set()
    for folder in ("e-series-joint-source-20261023", "e-series-compiler-closure-20261023",
                   "e-series-incidence-algebra-20261024"):
        paths.update(p for p in (HERE.parent/folder).iterdir()
                     if p.is_file() and p.suffix in (".py", ".txt", ".json"))
    for relative in (
        "verification/v1_e8_glue.py", "verification/v634_st31_structure.py",
        "tfpt_explorer/process.py",
        "experiments/theory-contracts/seam-joint-ward-state-20261014/operator_transport_check.py",
        "experiments/theory-contracts/marked-seam-monodromy-20261023/THEORY.txt",
        "experiments/theory-contracts/marked-seam-monodromy-20261023/checker.py",
        "experiments/theory-contracts/global-g31-lattice-lift-20261023/THEOREMS.txt",
    ):
        paths.add(REPO/relative)
    contents = {str(path.relative_to(REPO)): path.read_bytes() for path in sorted(paths)}
    contents["run_checks.py"] = RUNNER.encode()
    contents["READ-ME.txt"] = (HERE/"READ-ME.txt").read_bytes()
    contents["requirements.txt"] = b"numpy\nsympy\n"
    manifest = {name: sha256(data).hexdigest() for name, data in sorted(contents.items())}
    contents["manifest.json"] = (json.dumps(manifest, indent=2)+"\n").encode()
    destination.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(destination, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for name, data in sorted(contents.items()):
            archive.writestr("TFPT-E6-E7-E8-Pruefpaket/"+name, data)
    with zipfile.ZipFile(destination) as archive:
        if archive.testzip() is not None:
            raise RuntimeError("Archive CRC check failed")
        for name, data in contents.items():
            if archive.read("TFPT-E6-E7-E8-Pruefpaket/"+name) != data:
                raise RuntimeError("Archive source differs: "+name)
    return {"archive": str(destination), "files": len(contents),
            "bytes": destination.stat().st_size, "sha256": sha256(destination.read_bytes()).hexdigest()}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()
    print(json.dumps(make_bundle(args.destination.resolve()), indent=2))
