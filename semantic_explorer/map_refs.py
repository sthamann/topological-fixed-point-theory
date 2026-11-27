#!/usr/bin/env python3
"""Keep the source references of the authored theory map current without changing what they cite.

A reference stores the hash of the whole file.  Every later commit to that file makes the reference stale,
although the cited lines may still say the same thing; the verification ledger changes almost daily, so map
entries that cite it went stale and the server test for current references stayed red.

This tool re-verifies instead of re-reviewing.  A stale reference is refreshed only when its recorded excerpt
still stands at the cited lines, or is found exactly once elsewhere in the file (the lines then move with it).
A reference whose excerpt is gone is reported and left untouched: that one needs a human or agent review.

    python3 map_refs.py --check      exit 1 when any local reference is stale
    python3 map_refs.py --refresh    refresh the verifiable ones, report the rest
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import re
import sys

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
MAP = HERE / "theory_map.json"


def norm(text):
    return re.sub(r"\s+", " ", text).strip().rstrip("…").rstrip(".").strip()


def references(value):
    if isinstance(value, dict):
        if isinstance(value.get("path"), str) and isinstance(value.get("sha256"), str) \
                and value.get("source_location_state") != "moved_to_external_repository":
            yield value
        for child in value.values():
            yield from references(child)
    elif isinstance(value, list):
        for child in value:
            yield from references(child)


def locate(lines, excerpt, start, end):
    """Return the 1-based line range at which the excerpt stands, or None."""
    wanted = norm(excerpt)
    if not wanted:
        return (start, end) if 0 < start <= end <= len(lines) else None
    if 0 < start <= end <= len(lines) and wanted in norm(" ".join(lines[start - 1:end])):
        return start, end
    span = max(end - start, 0)
    head = wanted[:80]
    found = [i + 1 for i, line in enumerate(lines) if head in norm(line)
             and wanted in norm(" ".join(lines[i:i + span + 1]))]
    return (found[0], found[0] + span) if len(found) == 1 else None


def audit(refresh=False):
    text = MAP.read_text(encoding="utf-8")
    data = json.loads(text)
    report = {"checked": 0, "current": 0, "refreshed": [], "needs_review": []}
    cache = {}
    for ref in references(data):
        report["checked"] += 1
        path = REPO / ref["path"]
        if not path.is_file():
            report["needs_review"].append({"path": ref["path"], "issue": "missing file"})
            continue
        if ref["path"] not in cache:
            raw = path.read_bytes()
            cache[ref["path"]] = (sha256(raw).hexdigest(), raw.decode("utf-8", "replace").splitlines())
        digest, lines = cache[ref["path"]]
        start, end = ref.get("line_start"), ref.get("line_end")
        if ref["sha256"] == digest and isinstance(start, int) and isinstance(end, int) and 0 < start <= end <= len(lines):
            report["current"] += 1
            continue
        place = locate(lines, str(ref.get("excerpt") or ""), start or 0, end or 0)
        if place is None or not ref.get("excerpt"):
            report["needs_review"].append({"path": ref["path"], "line_start": start, "line_end": end,
                                           "issue": "excerpt no longer stands at or near the cited lines"})
            continue
        report["refreshed"].append({"path": ref["path"], "lines": list(place)})
        if refresh:
            ref["sha256"], ref["line_start"], ref["line_end"] = digest, place[0], place[1]
    if refresh and report["refreshed"]:
        MAP.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--check", action="store_true")
    group.add_argument("--refresh", action="store_true")
    args = parser.parse_args(argv)
    report = audit(refresh=args.refresh)
    stale = len(report["needs_review"]) + (0 if args.refresh else len(report["refreshed"]))
    print(json.dumps({"theory_map_references": report, "stale": stale,
                      "notice": "A refreshed hash confirms that the cited text still stands; it is no review of its content."},
                     ensure_ascii=False))
    return 1 if stale else 0


if __name__ == "__main__":
    sys.exit(main())
