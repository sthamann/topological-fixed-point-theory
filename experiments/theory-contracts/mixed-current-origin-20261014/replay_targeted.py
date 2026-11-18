"""Run the unchanged mixed-current checker with pinned historical inputs.

This is a targeted historical replay, not ``historical_overlay``: the latter
also checks the active ledger against its manifest and currently rejects an
unrelated ledger edit. No checker pins, active files or receipts are changed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "verification"))
from historical_source_replay import checked, reader_only_adapter, reconstruct_tex, source_paths


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True, type=Path,
                        help="new result JSON outside the repository")
    args = parser.parse_args()
    output = args.out.resolve()
    if output.is_relative_to(ROOT) or output.exists():
        parser.error("--out must be a new path outside the repository")
    manifest_path = ROOT / "verification/historical_replay/manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    replacements = {}
    for record in manifest["inputs"]:
        kind = record["kind"]
        if kind == "reader_only_python_adapter":
            original = checked(ROOT, record["snapshot"], record["historical_sha256"])
            current = checked(ROOT, record["path"], record["current_sha256"])
            reader_only_adapter(original, current)
        elif kind == "exact_modular_tex":
            checked(ROOT, record["path"], record["current_sha256"])
            original = reconstruct_tex(ROOT, record)
        elif kind == "relocated_neutral_helper":
            original = checked(ROOT, record["current_path"], record["historical_sha256"])
        elif kind == "frozen_external_historical_input":
            original = checked(ROOT, record["snapshot"], record["historical_sha256"])
        elif kind == "retained_ledger_projection":
            continue
        else:
            raise ValueError(f"Unknown manifest kind: {kind}")
        if sha(original) != record["historical_sha256"]:
            raise ValueError(f"Historical input changed: {record['path']}")
        replacements[record["path"]] = original

    checker = "experiments/theory-contracts/mixed-current-origin-20261014/checker.py"
    with tempfile.TemporaryDirectory(prefix="tfpt-mixed-targeted-") as temporary:
        overlay = Path(temporary)
        for relative in source_paths(ROOT):
            target = overlay / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(ROOT / relative, target)
        for relative, data in replacements.items():
            target = overlay / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(data)
        checker_output = overlay / "target-result.json"
        env = dict(os.environ)
        env["PYTHONPATH"] = str(overlay / "verification")
        subprocess.run([sys.executable, "-B", str(overlay / checker), "--repo",
                        str(overlay), "--out", str(checker_output)],
                       env=env, check=True)
        report = {
            "scope": "targeted_historical_replay_not_full_historical_overlay",
            "ledger_projection_checked": False,
            "current_source_pin_pass": False,
            "active_files_or_pins_modified": False,
            "manifest_sha256": sha(manifest_path.read_bytes()),
            "checker_sha256": sha((ROOT / checker).read_bytes()),
            "historical_inputs": {name: sha(data) for name, data in replacements.items()},
            "result": json.loads(checker_output.read_text(encoding="utf-8")),
        }
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open("x", encoding="utf-8") as handle:
            json.dump(report, handle, indent=2, ensure_ascii=False)
            handle.write("\n")
    print(f"Targeted historical replay report: {output}")


if __name__ == "__main__":
    main()
