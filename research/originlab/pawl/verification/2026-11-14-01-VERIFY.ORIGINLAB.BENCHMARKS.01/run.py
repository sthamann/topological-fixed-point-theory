"""Frozen execution gate. Physical F1 must remain unpromoted in control campaigns."""
import json
import os
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

here = Path(__file__).resolve().parent
repo = Path(os.environ.get("TFPT_REPO_ROOT",here.parents[4]))
entry = Path(os.environ.get("ORIGINLAB_ENTRY",here.parents[2]/"main.py"))
with TemporaryDirectory() as target:
    command = [sys.executable,"-B",str(entry),"--repo",str(repo),
               "--out",target,"--campaign","benchmarks","--jobs","1"]
    if toolchain := os.environ.get("ORIGINLAB_LEAN_TOOLCHAIN"):
        command.extend(("--lean-toolchain",toolchain))
    subprocess.run(command,check=True)
    receipt = json.loads((Path(target)/"round-0001/benchmarks.json").read_text())
    assert set(receipt["benchmarks"])=={f"B{i:02}" for i in range(1,10)}
    assert all(row["status"]=="PASS" for row in receipt["benchmarks"].values())
    assert receipt["physical_solution"] is False
    assert receipt["benchmarks"]["B07"]["detail"]["physical_status"]=="BLOCKED"
    assert receipt["benchmarks"]["B09"]["detail"]["physical_status"]=="OPEN"
    Path(__file__).with_name("last_run.json").write_text(json.dumps(receipt,indent=2)+"\n")
