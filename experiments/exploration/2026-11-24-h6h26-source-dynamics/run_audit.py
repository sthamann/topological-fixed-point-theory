"""Reproduce the H6h26/H6h27 follow-up within the existing OriginLab lane.

Use --lane to also run the historical finite composed universe. The default
certificate run is fast and does not overwrite any historical result.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "research" / "originlab"))

from process_source_dynamics import source_dynamics_audit


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--lane", action="store_true")
    parser.add_argument("--output", type=Path, default=HERE / "certificate_results.json")
    args = parser.parse_args()
    result = source_dynamics_audit()
    if args.lane:
        from process_rules import run
        result = run(certified=True)
    result["reproduction"] = {
        "python": sys.version.split()[0],
        "file_hashes": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                        for p in sorted((ROOT / "research" / "originlab").glob("process_*.py"))},
        "historical_sources": {
            name: hashlib.sha256((ROOT / "experiments" / "exploration" /
                "2026-11-23-phasenraum-f5-holocode" / name).read_bytes()).hexdigest()
            for name in ("z48_zuordnung_quelle.py", "z49_flavourgruppen.py",
                         "z50_metrischer_stern.py", "z51_kante_quartik.py")},
    }
    args.output.write_text(json.dumps(result, indent=2, ensure_ascii=False,
                                      sort_keys=True, allow_nan=False)+"\n")
    audit = result["source_dynamics"] if args.lane else result
    print(json.dumps({"all_checks_pass": audit["all_checks_pass"],
                      "physical_solution": audit["physical_solution"],
                      "checks": audit["checks"], "output": str(args.output)},
                     ensure_ascii=False, sort_keys=True))
    if not audit["all_checks_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
