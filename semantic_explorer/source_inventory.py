"""Read-only, source-native inventory for the TFPT Atlas builders.

This module replaces Atlas' dependence on the generated legacy theory graph.
It reads only stable ledgers/catalogs and a small allowlist of experiment
records; it never executes experiment code or generates a second graph file.
Relationships are emitted only from explicit IDs present in those sources.
"""
from __future__ import annotations

import hashlib
import json
import re
import argparse
import csv
import sys
from pathlib import Path
from typing import Any

CLAIM_TOKEN_RE = re.compile(r"\b([A-Z][A-Z0-9_-]*(?:\.[A-Z0-9_-]+){1,5})\b")
VERDICT_RE = re.compile(r"(?i)\b(?:verdict|status)\s*[:=]\s*[`* ]*([A-Z][A-Z0-9_ -]{1,50})")
DATE_SUFFIX_RE = re.compile(r"-(\d{8})$")
HEADING_RE = re.compile(r"^#{1,3}\s+(.+)$", re.M)
SECTION_RE = re.compile(r"^(#{1,6})\s+(.+?)\s*$", re.M)
GATE_ID_RE = re.compile(r"^T[1-8]$")
GATE_RANGE_OPEN_RE = re.compile(r"T1\s*[–—-]\s*T8[^.]*?(bleiben? offen|remains? open|stay open)", re.I)
GATE_SINGLE_OPEN_RE = re.compile(r"\bT([1-8])\b")
FINDING_TRIGGER_RE = re.compile(r"widerlegt|refuted|no[- ]?go|gekillt|killed|widerspruch|contradiction|bestätigt|confirmed|bewiesen|proved|offen", re.I)
FINDING_REFUTES_RE = re.compile(r"widerlegt|refuted|widerspruch|contradiction", re.I)
FINDING_NOGO_RE = re.compile(r"no[- ]?go|gekillt|killed", re.I)
FINDING_SUPPORTS_RE = re.compile(r"bestätigt|confirmed|bewiesen|proved", re.I)
FINDING_UNPROVED_RE = re.compile(
    r"\bnicht(?:\s+(?:noch|bisher|damit|dadurch|daraus|hier|vollständig|direkt|streng|rigoros|automatisch|bereits|formal)){0,4}"
    r"(?:\s+aus\s+[^\s,;:.!?]+(?:\s+[^\s,;:.!?]+){0,5})?\s+(?:bewiesen|hergeleitet)\b"
    r"|\bnot(?:\s+(?:yet|been|thereby|thus|here|fully|directly|rigorously|automatically|currently|already|formally)){0,4}"
    r"(?:\s+from\s+[^\s,;:.!?]+(?:\s+[^\s,;:.!?]+){0,5})?\s+(?:proved|proven|derived)\b", re.I)
FINDING_RESOLVED_RE = re.compile(
    r"(?:\b(?:zuvor|vorher|vormals|damals|previously|formerly)\b.{0,80}"
    r"\b(?:offen(?:gelassen\w*)?|open|unresolved)\w*.{0,140}"
    r"\b(?:konkret\s+)?(?:ausgewertet|beantwortet|geklärt|gelöst|geschlossen|erledigt|resolved|answered|settled)\b"
    r"|\b(?:offen(?:e|en|er|es)?|open|unresolved)\w*.{0,100}"
    r"\b(?:wird|wurde|ist|sind|has been|have been|was|were)\b.{0,80}"
    r"\b(?:beantwortet|geklärt|gelöst|geschlossen|erledigt|resolved|answered|settled)\b)", re.I)
FINDING_OPEN_RE = re.compile(
    r"\b(?:offene\s+(?:herkunfts)?frage|offene\s+lücke|open\s+(?:question|issue|gap)|unresolved\s+(?:question|issue|gap))\b"
    r"|\b(?:bleibt|bleiben|verbleibt|verbleiben|remains?|still\s+is|weiterhin)\b.{0,100}"
    r"\b(?:offen|ungeklärt|unbeantwortet|unresolved|not\s+derived)\b"
    r"|\b(?:offen|ungeklärt|unbeantwortet|unresolved|not\s+derived)\b.{0,100}"
    r"\b(?:noch|weiterhin|remains?|remains?\s+open)\b", re.I)
NEXTXT_LINE_RE = re.compile(
    r"^(?:#\s*(\d{4}-\d{2}-\d{2})(?:\s*\([^)]*\))?\s*[—–-]?\s+|"
    r"\[(\d{4}-\d{2}-\d{2})\]\s*)(.*\S)\s*$"
)
VERDICT_PATTERNS = [
    re.compile(r"\*\*Verdict:\s*`([^`]+)`\.\*\*"),
    re.compile(r"Verdict:\s*\*\*`([^`]+)`\*\*"),
    re.compile(r"Contract-Verdict nach pr.registriertem Enum:\s*\n+\s*`?([A-Z][A-Z0-9_]*)`?"),
    re.compile(r"Inhaltliches(?:\s+pr.registriertes)?\s+Verdict:\s*\*\*`([^`]+)`\*\*"),
    re.compile(r"\*\*Verdict:\*\*\s*`([^`]+)`"),
    re.compile(r"\*\*Verdict:\s*`([^`]+)`\*\*"),
    re.compile(r"(?<![A-Za-z])Verdict:\s*`([^`]+)`"),
    re.compile(r"Success verdict:\s*`([^`]+)`"),
    re.compile(r"(?<![A-Za-z])[Vv]erdict:\s*\*\*`([^`]+)`\*\*"),
    re.compile(r"\*\*Verdict:\s*([A-ZÄÖÜ][A-Za-zÄÖÜäöü0-9_ +.-]{1,60}?)\s*\.\*\*"),
    re.compile(r"\*\*[A-Za-z-]*[Vv]erdict:\s*([A-ZÄÖÜ]{2,}(?:\s+[A-ZÄÖÜ]+){0,3})\s*\.?\*\*"),
    re.compile(r"(?<![A-Za-z])Verdicts?:\s*\*\*[^\n`]*`([A-Z0-9_]+)`"),
    re.compile(r"Record verdict:\s*`?([A-Z][A-Z0-9_]*)`?"),
    re.compile(r"(?<![A-Za-z])[Vv]erdict\s*[:=]\s*\*\*([A-Z][A-Z0-9_ -]{1,50})\*\*"),
]
TEXT_FILES = ("PROOF.txt", "PROOF.md", "RESULTS.md", "README.md")
JSON_FILES = ("contract_index.json", "validation.json", "result.json", "results.json")
SKIP_DIRS = {"__pycache__", ".git", "node_modules", ".lake", ".venv", "venv", "data", "logs", "outputs", "_archive"}
MAX_DECLARED_ARTIFACT_BYTES = 262_144
MAX_DECLARED_SOURCE_BYTES = 1_048_576


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _source(root: Path, path: Path, sources: dict[str, dict[str, Any]]) -> str | None:
    try:
        resolved = path.resolve()
        resolved.relative_to(root.resolve())
        data = resolved.read_bytes()
    except (OSError, ValueError):
        return None
    rel = resolved.relative_to(root.resolve()).as_posix()
    digest = _sha(data)
    sid = "source:" + hashlib.sha256(f"{rel}|{digest}".encode()).hexdigest()[:20]
    text = data.decode("utf-8", errors="replace")
    sources.setdefault(sid, {
        "id": sid, "path": rel, "sha256": digest,
        "line_start": 1, "line_end": len(text.splitlines()),
        "line_count": len(text.splitlines()), "language": "mixed",
        "kind": resolved.suffix.lstrip(".").lower() or "text", "bytes": len(data),
    })
    return sid


def _node(node_id: str, kind: str, label: str, attrs: dict[str, Any], summary: str = "") -> dict[str, Any]:
    return {"id": node_id, "type": kind, "label": label, "summary": summary, "attrs": attrs}


def _edge(edges: list[dict[str, Any]], src: str, dst: str, kind: str, **attrs: Any) -> None:
    if src == dst:
        return
    edges.append({"src": src, "dst": dst, "type": kind, "attrs": attrs})


def _slugify(text: Any, width: int = 60) -> str:
    slug = re.sub(r"[^a-z0-9]+", "-", str(text or "").lower()).strip("-")
    return slug[:width].strip("-") or "row"


def _read_json(path: Path) -> dict[str, Any] | None:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    return value if isinstance(value, dict) else None


def _read_json_trailing_stripped(path: Path) -> dict[str, Any] | None:
    """Accept run artefacts with a non-JSON progress line after their object."""
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
        idx = text.rfind("}")
        if idx >= 0:
            text = text[:idx + 1]
        value = json.loads(text)
    except (OSError, ValueError):
        return None
    return value if isinstance(value, dict) else None


def _clip(value: Any, width: int) -> str:
    if value is None:
        return ""
    text = " ".join(str(value).split())
    return text if len(text) <= width else text[:width - 3] + "..."


def _harvest_json_attrs(rec: dict[str, Any], data: dict[str, Any]) -> None:
    """Keep the source parser's independent run-artefact field merge."""
    if rec["exact_checks"] is None:
        exact = data.get("exact_checks", data.get("check_count"))
        if exact is not None:
            rec["exact_checks"] = exact
    if rec["promotion"] is None and data.get("promotion") is not None:
        rec["promotion"] = bool(data.get("promotion"))
    for token in data.get("closed_gate_ids") or []:
        token = str(token)
        if GATE_ID_RE.match(token):
            if token not in rec["closes"]:
                rec["closes"].append(token)
        else:
            rec["claims_tested_raw"].append(token)
    missing = data.get("missing") or []
    if missing and not rec["missing"]:
        rec["missing"] = [_clip(item, 80) for item in list(missing)[:10]]
    if rec["t1_t8_closed"] is None and "T1_T8_closed" in data:
        rec["t1_t8_closed"] = data.get("T1_T8_closed")
    not_derived = data.get("not_derived") or []
    if not_derived and not rec["not_derived"]:
        rec["not_derived"] = [_clip(item, 80) for item in list(not_derived)[:10]]


def _contract_json_candidates(folder: Path) -> list[Path]:
    out: list[Path] = []
    seen: set[str] = set()
    for pattern in ("*normal*.json", "*optimized*.json", "replay.json"):
        for path in sorted(folder.glob(pattern)):
            if path.name not in seen:
                seen.add(path.name)
                out.append(path)
    results = folder / "results" / "results.json"
    if results.is_file():
        out.append(results)
    for path in sorted(folder.glob("*.json")):
        if path.name not in seen:
            seen.add(path.name)
            out.append(path)
    return out


def _find_verdict_in_text(text: str | None) -> str | None:
    if not text:
        return None
    for pattern in VERDICT_PATTERNS:
        match = pattern.search(text)
        if match:
            return match.group(1).strip().rstrip(".").strip()
    return None


def _verdict_class(raw: Any) -> str:
    if not raw:
        return "unknown"
    normalized = re.sub(r"[^a-z0-9]+", "_", str(raw).lower()).strip("_")
    if not normalized or normalized.upper() in {"UNKNOWN", "UNBEKANNT", "KEIN_ERGEBNIS"}:
        return "unknown"
    upper = normalized.upper()
    negated = any(term in normalized for term in ("not_confirmed", "not_refuted", "unrefuted", "nicht_widerlegt"))
    kill_passed = "kill_test" in normalized and any(term in normalized for term in ("passed", "survived"))
    if not negated and not kill_passed:
        if "REFUT" in upper or "widerlegt" in normalized:
            return "refuted"
        if any(term in upper for term in ("KERNEL_VIOLATION", "KILLED", "NO_GO")) or "gekillt" in normalized:
            return "killed"
    if "CANDIDATE" in upper:
        return "candidate"
    if any(term in upper for term in ("PARTIAL", "PARTIELL", "MARGINAL")) or "conditional" in normalized:
        return "partial"
    if any(term in normalized for term in ("bleibt_offen", "bleiben_offen", "remains_open")):
        return "open"
    if "FAIL" in upper or upper.startswith("BLOCKED") or upper in {"NEIN", "NO", "NOT_VIABLE", "NOTVIABLE"}:
        return "fail"
    if any(term in upper for term in ("PASS", "HOLDS", "CERTIFIED", "BESTANDEN", "PROVED")) or "best_tigt" in normalized or "confirmed" in normalized or ("byteidentisch" in normalized and "pass" in normalized) or upper in {"CONSISTENT", "CLEAN", "EXAKT_ERHALTEN", "JA", "YES"}:
        return "pass"
    return "other"


def _extract_gate_mentions(*texts: str | None) -> list[str]:
    opened: set[str] = set()
    for text in texts:
        if not text:
            continue
        for sentence in re.split(r"(?<=[.!?])\s+", " ".join(text.splitlines())):
            if GATE_RANGE_OPEN_RE.search(sentence):
                opened.update(f"T{i}" for i in range(1, 9))
            elif re.search(r"offen|open", sentence, re.I):
                opened.update(f"T{gate}" for gate in GATE_SINGLE_OPEN_RE.findall(sentence))
    return sorted(opened)


def _finding_kind(sentence: str) -> str:
    if FINDING_REFUTES_RE.search(sentence):
        return "refutes"
    if FINDING_NOGO_RE.search(sentence):
        return "no_go"
    # Negation belongs to the proof predicate, not to any nearby sentence word.
    # In particular, "not only proved" and "does not refute" do not match.
    if FINDING_UNPROVED_RE.search(sentence):
        return "open"
    if FINDING_SUPPORTS_RE.search(sentence):
        return "supports"
    # A past open item that the source says was evaluated or answered is a
    # completed finding. Do not let the word "open" override its tense/status.
    if FINDING_RESOLVED_RE.search(sentence):
        return "supports"
    # Only affirmative unresolved language becomes an open finding. Other
    # sentences mentioning openness remain descriptive; contract_index fields
    # such as leaves_open are the explicit status authority.
    if FINDING_OPEN_RE.search(sentence):
        return "open"
    return "other"


def _declared_contract_artifacts(root: Path, folder: Path, indexed: dict[str, Any] | None) -> list[tuple[str, Path]]:
    """Resolve only explicitly declared proof/checker files inside this repository.

    A declaration is provenance, not authority to read arbitrary paths: absolute
    paths, traversal, unsupported file types, files outside the repository and
    oversized artifacts are excluded.
    """
    if not isinstance(indexed, dict):
        return []
    specs = {"proof": {".md", ".txt"}, "checker": {".py"}}
    found: list[tuple[str, Path]] = []
    repo_root = root.resolve()
    folder = folder.resolve()
    for field, suffixes in specs.items():
        raw = indexed.get(field)
        if not isinstance(raw, str) or not raw.strip():
            continue
        declared = Path(raw.strip())
        if declared.is_absolute() or ".." in declared.parts:
            continue
        candidates = (folder / declared, repo_root / declared, repo_root / "verification" / declared)
        for candidate in candidates:
            try:
                resolved = candidate.resolve()
                resolved.relative_to(repo_root)
                size = resolved.stat().st_size
            except (OSError, ValueError):
                continue
            if resolved.is_file() and resolved.suffix.lower() in suffixes and size <= MAX_DECLARED_ARTIFACT_BYTES:
                found.append((field, resolved))
                break
    return found


def _declared_source_manifest(root: Path, folder: Path, indexed: dict[str, Any] | None) -> tuple[Path | None, list[dict[str, Any]]]:
    """Read a declared source manifest and resolve its bounded repository files.

    Supported shapes are ``sources`` or ``files`` arrays of objects with a
    repository-relative ``path`` and optional SHA-256. Paths outside the repo,
    traversal, oversized files, and non-file entries are ignored. Hash
    mismatches remain visible in metadata while the current source is still
    registered with its actual content hash.
    """
    if not isinstance(indexed, dict):
        return None, []
    raw = indexed.get("source_manifest")
    if not isinstance(raw, str) or not raw.strip():
        return None, []
    repo_root = root.resolve()
    folder = folder.resolve()
    declared = Path(raw.strip())
    if ".." in declared.parts:
        return None, []
    candidates = (folder / declared, repo_root / declared)
    manifest_path = None
    for candidate in candidates:
        try:
            resolved = candidate.resolve()
            resolved.relative_to(repo_root)
            if resolved.is_file() and resolved.stat().st_size <= MAX_DECLARED_ARTIFACT_BYTES:
                manifest_path = resolved
                break
        except (OSError, ValueError):
            continue
    if manifest_path is None:
        return None, []
    payload = _read_json(manifest_path)
    if not isinstance(payload, dict):
        return manifest_path, []
    entries = payload.get("sources")
    if not isinstance(entries, list):
        entries = payload.get("files")
    if not isinstance(entries, list):
        return manifest_path, []
    resolved_entries: list[dict[str, Any]] = []
    for entry in entries:
        if not isinstance(entry, dict) or not isinstance(entry.get("path"), str):
            continue
        raw_path = entry["path"].strip()
        if not raw_path:
            continue
        rel_path = Path(raw_path)
        if ".." in rel_path.parts:
            continue
        path_candidates = [rel_path] if rel_path.is_absolute() else [repo_root / rel_path, manifest_path.parent / rel_path]
        source_path = None
        for candidate in path_candidates:
            try:
                resolved = candidate.resolve()
                resolved.relative_to(repo_root)
                if resolved.is_file() and resolved.stat().st_size <= MAX_DECLARED_SOURCE_BYTES:
                    source_path = resolved
                    break
            except (OSError, ValueError):
                continue
        if source_path is None:
            continue
        expected = entry.get("sha256")
        if not isinstance(expected, str) or not re.fullmatch(r"[0-9a-fA-F]{64}", expected):
            expected = None
        actual = _sha(source_path.read_bytes())
        resolved_entries.append({
            "path": source_path.relative_to(repo_root).as_posix(),
            "expected_sha256": expected.lower() if expected else None,
            "actual_sha256": actual,
            "hash_matches": (actual == expected.lower()) if expected else None,
        })
    return manifest_path, resolved_entries


def _extract_findings(readme: str | None, results: str | None, cap: int = 10) -> list[dict[str, Any]]:
    findings: list[dict[str, Any]] = []
    seen: set[str] = set()
    for where, text in (("readme", readme), ("results", results)):
        if not text:
            continue
        section = None
        for block in re.split(r"\n\s*\n", text):
            block = block.strip()
            if not block:
                continue
            heading = SECTION_RE.match(block)
            if heading:
                section = _clip(re.sub(r"[#*`]", "", heading.group(2)), 120)
                continue
            if block.startswith("#") or (section and "verdict-enum" in section.lower()):
                continue
            paragraph = " ".join(block.splitlines())
            for sentence in re.split(r"(?<=[.!?])\s+", paragraph):
                if not (FINDING_TRIGGER_RE.search(sentence) or FINDING_UNPROVED_RE.search(sentence)):
                    continue
                clean = re.sub(r"\*\*|__|`", "", sentence).strip()
                clean = " ".join(clean.split())
                if len(clean) < 25 or re.match(r"^-?\s*[A-Z][A-Z0-9_]{3,}\s*—", clean):
                    continue
                key = clean[:120]
                if key in seen:
                    continue
                seen.add(key)
                findings.append({"kind": _finding_kind(clean), "text": _clip(clean, 200), "section": section, "where": where})
                if len(findings) >= cap:
                    return findings
    return findings


def _apply_findings_overrides(root: Path, contract_name: str, findings: list[dict[str, Any]]) -> list[dict[str, Any]]:
    path = root / "verification" / "findings_overrides.json"
    payload = _read_json(path)
    if not payload:
        return findings
    overrides = payload.get("overrides") if isinstance(payload.get("overrides"), dict) else payload
    kept: list[dict[str, Any]] = []
    for finding in findings:
        key = f"{contract_name}#{hashlib.sha1(finding['text'].encode()).hexdigest()[:12]}"
        override = overrides.get(key) if isinstance(overrides, dict) else None
        if not isinstance(override, dict):
            kept.append(finding)
            continue
        if override.get("action") == "drop":
            continue
        if override.get("action") == "correct" and override.get("kind") in {"refutes", "no_go", "supports", "open"}:
            finding = {**finding, "kind": override["kind"], "override": "llm"}
        kept.append(finding)
    return kept


def _contract_nodes(root: Path, valid_claims: set[str], sources: dict[str, dict[str, Any]], edges: list[dict[str, Any]]) -> list[dict[str, Any]]:
    result: list[dict[str, Any]] = []
    contracts = root / "experiments" / "theory-contracts"
    if not contracts.is_dir():
        return result
    override_source = root / "verification" / "findings_overrides.json"
    override_source_id = _source(root, override_source, sources) if override_source.is_file() else None
    next_txt_path = root / "experiments" / "next.txt"
    next_txt_source_id = _source(root, next_txt_path, sources) if next_txt_path.is_file() else None
    contract_names = {
        path.name for path in contracts.iterdir()
        if path.is_dir() and not path.name.startswith(".") and path.name != "__pycache__"
    }
    contract_names_ordered = sorted(contract_names, key=lambda value: (-len(value), value))
    slug_pattern = re.compile(
        r"(?<![a-z0-9-])(?:" + "|".join(re.escape(name) for name in contract_names_ordered) + r")(?![a-z0-9-])",
        re.I,
    ) if contract_names_ordered else None
    next_txt_pairs: set[tuple[str, str]] = set()
    if next_txt_path.is_file() and slug_pattern is not None:
        try:
            next_lines = next_txt_path.read_text(encoding="utf-8", errors="replace").splitlines()
        except OSError:
            next_lines = []
        for line in next_lines:
            if not NEXTXT_LINE_RE.match(line):
                continue
            slugs = list(dict.fromkeys(slug_pattern.findall(line)))
            for index, left in enumerate(slugs):
                for right in slugs[index + 1:]:
                    next_txt_pairs.add((left, right))
    for folder in sorted(contracts.iterdir(), key=lambda p: p.name):
        if not folder.is_dir() or folder.name.startswith(".") or folder.name in SKIP_DIRS:
            continue
        indexed_path = folder / "contract_index.json"
        validation_path = folder / "validation.json"
        indexed = _read_json(indexed_path)
        validation = _read_json(validation_path)
        declared_artifacts = _declared_contract_artifacts(root, folder, indexed)
        manifest_path, manifest_entries = _declared_source_manifest(root, folder, indexed)
        # Mirror the legacy source parser's precedence: contract_index owns
        # declared fields, validation/run artefacts still fill independent
        # evidence fields even when an index supplies the verdict.
        record: dict[str, Any] = {
            "name": folder.name, "date": None, "question": None,
            "verdict_raw": None, "verdict_class": "unknown", "checker": None,
            "parse_source": "none", "claims_tested": [], "claims_tested_raw": [],
            "unresolved_claim_tokens": [], "gates": [], "kills": [], "kills_free": [],
            "exact_checks": None, "promotion": None, "missing": [],
            "t1_t8_closed": None, "not_derived": [], "closes": [],
            "leaves_open": [], "prose_gate_mentions": [], "topics": [], "scope": None,
            "unresolved_topic_ids": [], "findings": [], "summary_de": None,
            "relations_llm": [],
        }
        date_match = DATE_SUFFIX_RE.search(folder.name)
        if date_match:
            record["date"] = date_match.group(1)
        readme_path = folder / "README.md"
        results_path = folder / "RESULTS.md"
        readme_text = None
        results_text = None
        try:
            if readme_path.is_file():
                readme_text = readme_path.read_text(encoding="utf-8", errors="replace")
                heading = HEADING_RE.search(readme_text)
                if heading:
                    record["question"] = _clip(re.sub(r"[#*`]", "", heading.group(1)), 200)
        except OSError:
            pass
        try:
            if results_path.is_file():
                results_text = results_path.read_text(encoding="utf-8", errors="replace")
        except OSError:
            pass
        checker_files = sorted(path.name for path in folder.glob("*.py") if path.name != "__init__.py")
        if checker_files:
            record["checker"] = checker_files[0]
        refs: list[str] = []
        source_paths = [indexed_path, validation_path, *[folder / name for name in TEXT_FILES],
                        *[path for _field, path in declared_artifacts],
                        *([manifest_path] if manifest_path else []),
                        *[root / entry["path"] for entry in manifest_entries]]
        source_id_by_path: dict[str, str] = {}
        for path in source_paths:
            if path.is_file():
                source_id = _source(root, path, sources)
                if source_id and source_id not in refs:
                    refs.append(source_id)
                if source_id:
                    source_id_by_path[path.resolve().relative_to(root.resolve()).as_posix()] = source_id
        # (a) contract_index is the declared record and relation authority.
        if indexed is not None:
            raw = indexed.get("verdict")
            record["verdict_raw"] = str(raw) if isinstance(raw, (str, int, float)) else None
            if record["verdict_raw"]:
                record["parse_source"] = "contract_index_llm" if indexed.get("source") == "llm_backfill" else "contract_index"
            record["question"] = _clip(indexed.get("question") or record["question"] or "", 200) or None
            record["date"] = str(indexed.get("date") or record["date"] or "") or None
            record["checker"] = indexed.get("checker") or record["checker"]
            record["summary_de"] = _clip(indexed.get("summary_de"), 400) or None
            record["scope"] = _clip(indexed.get("scope"), 1200) or None
            declared_topics = indexed.get("topics")
            if isinstance(declared_topics, list):
                known_topics = {
                    "topic:origin", "topic:e8", "topic:clocks", "topic:flavor",
                    "topic:alpha", "topic:quantum", "topic:gravity", "topic:cosmology",
                    "topic:rh", "topic:validation", "topic:universalraum",
                    "topic:metrology", "topic:spacetime", "topic:process",
                }
                for item in declared_topics:
                    topic_id = str(item)
                    if topic_id in known_topics and topic_id not in record["topics"]:
                        record["topics"].append(topic_id)
                    elif topic_id not in record["unresolved_topic_ids"]:
                        record["unresolved_topic_ids"].append(topic_id)
            record["claims_tested_raw"].extend(str(token) for token in (indexed.get("claims_tested") or []))
            record["gates"] = list(indexed.get("gates") or [])
            record["kills"].extend(str(item) for item in (indexed.get("kills") or []))
            for relation in indexed.get("relations") or []:
                if isinstance(relation, dict) and relation.get("target"):
                    record["relations_llm"].append({"type": relation.get("type"), "target": str(relation["target"])})
            _harvest_json_attrs(record, indexed)
        # (b) validation contributes independent fields and only supplies a
        # verdict when the contract index did not declare one.
        if validation is not None:
            _harvest_json_attrs(record, validation)
            for token in validation.get("claims_tested") or []:
                record["claims_tested_raw"].append(str(token))
            if record["verdict_raw"] is None:
                raw = validation.get("verdict")
                if not isinstance(raw, (str, int, float)):
                    raw = validation.get("status")
                if isinstance(raw, (str, int, float)) and str(raw).strip():
                    record["verdict_raw"] = str(raw)
                    record["parse_source"] = "validation"
        # (c) prose result, then (d) prioritized result artefacts. JSON attrs
        # are merged whether or not they supply the eventual verdict.
        if record["verdict_raw"] is None and results_text:
            raw = _find_verdict_in_text(results_text)
            if raw:
                record["verdict_raw"] = raw
                record["parse_source"] = "results_md"
        for candidate in _contract_json_candidates(folder):
            source_id = _source(root, candidate, sources)
            if source_id and source_id not in refs:
                refs.append(source_id)
            payload = _read_json_trailing_stripped(candidate)
            if payload is None:
                continue
            is_named = candidate.name != "results.json" or candidate.parent.name == "results"
            if not is_named and not any(key in payload for key in ("status", "verdict")):
                continue
            _harvest_json_attrs(record, payload)
            if record["verdict_raw"] is None:
                raw = payload.get("verdict")
                if not isinstance(raw, (str, int, float)):
                    raw = payload.get("decision")
                if not isinstance(raw, (str, int, float)):
                    raw = payload.get("status")
                if isinstance(raw, (str, int, float)) and str(raw).strip():
                    record["verdict_raw"] = str(raw)
                    record["parse_source"] = "run_json"
        # (e) and (f): keep the same last-resort text verdict behavior.
        if record["verdict_raw"] is None and readme_text:
            raw = _find_verdict_in_text(readme_text)
            if raw:
                record["verdict_raw"] = raw
                record["parse_source"] = "readme_md"
        record["verdict_class"] = _verdict_class(record["verdict_raw"])
        if record["question"] is None:
            record["question"] = folder.name
        declared_open = indexed.get("leaves_open") if indexed else None
        declared_open = declared_open if isinstance(declared_open, list) else []
        # Only structured contract metadata can declare a gate as open. Prose
        # mentions stay available for navigation, but a sentence such as
        # “T1–T8 remain open” cannot create eight status edges when the index
        # explicitly declares gates=[] or omits per-gate scope.
        record["leaves_open"] = list(dict.fromkeys(str(item) for item in declared_open))
        record["prose_gate_mentions"] = _extract_gate_mentions(readme_text, results_text)
        record["findings"] = _extract_findings(readme_text, results_text, cap=10)
        record["findings"] = _apply_findings_overrides(root, folder.name, record["findings"])

        claims: list[str] = []
        raw_claims: list[str] = []
        for token in record.pop("claims_tested_raw"):
            claim = token if token in valid_claims else None
            if claim is None:
                matches = [candidate for candidate in valid_claims if candidate.endswith("." + token)]
                if len(matches) == 1:
                    claim = matches[0]
            if claim and claim not in claims:
                claims.append(claim)
                raw_claims.append(token)
        record["unresolved_claim_tokens"] = []
        for text in (readme_text, results_text):
            if not text:
                continue
            for token in CLAIM_TOKEN_RE.findall(text):
                claim = token if token in valid_claims else None
                matches = [] if claim else [candidate for candidate in valid_claims if candidate.endswith("." + token)]
                if claim is None and len(matches) == 1:
                    claim = matches[0]
                elif claim is None and len(matches) > 1 and token not in record["unresolved_claim_tokens"]:
                    record["unresolved_claim_tokens"].append(token)
                if claim and claim not in claims and len(claims) < 20:
                    claims.append(claim)
        record["claims_tested"] = claims
        # Explicit claim IDs and explicit gate mentions are source-backed
        # relationships. The graph uses goals for the old gate:T* targets.
        for claim in claims:
            _edge(edges, f"contract:{folder.name}", f"claim:{claim}", "tests", source="source_contract_claim_reference", reviewed=False)
        for gate in record["leaves_open"]:
            goal = gate if gate.startswith("goal:") else f"goal:{gate}"
            if goal.rsplit(":", 1)[-1] in {f"T{i}" for i in range(1, 9)}:
                _edge(edges, f"contract:{folder.name}", goal, "leaves_open", source="source_gate_annotation", reviewed=False)
        for gate in record["closes"]:
            goal = gate if gate.startswith("goal:") else f"goal:{gate}"
            if goal.rsplit(":", 1)[-1] in {f"T{i}" for i in range(1, 9)}:
                _edge(edges, f"contract:{folder.name}", goal, "closes", source="explicit_closed_gate_ids", reviewed=False)
        for relation in record["relations_llm"]:
            target = relation.get("target")
            kind = relation.get("type")
            if isinstance(target, str) and target and kind in {"references", "followup_of", "contradicts", "supports", "uses", "extends", "refines"}:
                target_name = target.removeprefix("contract:")
                if target_name in contract_names:
                    _edge(edges, f"contract:{folder.name}", f"contract:{target_name}", kind, source="contract_index_relation", reviewed=False)
        for topic_id in record["topics"]:
            _edge(edges, f"contract:{folder.name}", topic_id, "about", source="contract_index_topics", reviewed=False)
        # Match only exact existing contract slugs in README/RESULTS, mirroring
        # the legacy enumerator. These links communicate citation/reference,
        # not derivation or proof.
        for where, text in (("readme", readme_text), ("results", results_text)):
            if not text or slug_pattern is None:
                continue
            for target in sorted(set(slug_pattern.findall(text))):
                if target.lower() != folder.name.lower():
                    _edge(edges, f"contract:{folder.name}", f"contract:{target}", "references", source=where, reviewed=False)
        # Current next.txt entries can explicitly co-mention multiple known
        # contract names. Retain that as a navigation reference just like the
        # legacy source enumerator did; it is not a causal relation.
        for left, right in next_txt_pairs:
            if folder.name == left:
                _edge(edges, f"contract:{left}", f"contract:{right}", "references", source="next_txt", reviewed=False)
        # Old records split textual negative descriptions from exact ledger
        # kills. Keep the same distinction and emit kill edges only for exact
        # claim IDs.
        for killed in record["kills"]:
            if killed in valid_claims:
                _edge(edges, f"contract:{folder.name}", f"claim:{killed}", "kills", source="contract_index_kills", reviewed=False)
        attrs = {
            **record,
            "original_contract": indexed or validation or {},
            "kills_free": [_clip(item, 120) for item in record["kills"] if item not in valid_claims],
            "claims_tested": claims,
            "source_refs": refs,
            "source_manifest_entries": [
                {**entry, "source_ref": source_id_by_path.get(entry["path"])}
                for entry in manifest_entries
            ],
            "inventory_source": "source_native_contract_directory",
            "attempt_kind": "contract",
            "reviewed": False,
            "notice": "Source-inventoried contract record; status and claims remain as stated in source.",
        }
        if override_source_id:
            attrs["source_refs"].append(override_source_id)
        if next_txt_source_id and any(left == folder.name for left, _right in next_txt_pairs):
            attrs["source_refs"].append(next_txt_source_id)
        node_id = f"contract:{folder.name}"
        label = str(attrs.get("question") or folder.name)
        result.append(_node(node_id, "contract", label, attrs, str(attrs.get("summary_de") or attrs.get("verdict_raw") or "")))
    return result


def _experiment_nodes(root: Path, valid_claims: set[str], sources: dict[str, dict[str, Any]], edges: list[dict[str, Any]]) -> list[dict[str, Any]]:
    result: list[dict[str, Any]] = []
    exp_root = root / "experiments"
    if not exp_root.is_dir():
        return result
    for folder in sorted(exp_root.iterdir(), key=lambda p: p.name):
        if not folder.is_dir() or folder.name.startswith(".") or folder.name in SKIP_DIRS or folder.name == "theory-contracts":
            continue
        readme = folder / "README.md"
        results_file = folder / "results" / "results.json"
        preregs = sorted((folder / "hypotheses").glob("*.yaml")) if (folder / "hypotheses").is_dir() else []
        preregs += sorted(folder.glob("prereg*.yaml"))
        attrs: dict[str, Any] = {
            "path": folder.relative_to(root).as_posix(),
            "has_prereg": bool(preregs),
            "has_results_dir": (folder / "results").is_dir(),
            "firewall": "experiments",
            "reviewed": False,
        }
        for path in (readme, results_file):
            if path.is_file():
                sid = _source(root, path, sources)
                if sid:
                    attrs.setdefault("source_refs", []).append(sid)
        if preregs:
            try:
                raw = preregs[0].read_bytes()
                attrs["prereg_sha256"] = _sha(raw)
            except OSError:
                pass
        if results_file.is_file():
            payload = _read_json(results_file)
            if payload:
                status = payload.get("verdict", payload.get("status"))
                if status is None and isinstance(payload.get("overall"), dict):
                    status = payload["overall"].get("verdict")
                if status is not None:
                    attrs["verdict_raw"] = str(status)
        node_id = f"experiment:{folder.name}"
        result.append(_node(node_id, "experiment", folder.name, attrs, "Source-inventoried experiment directory; no result is independently validated."))
        texts: list[str] = []
        for path in (readme, results_file):
            try:
                if path.is_file():
                    texts.append(path.read_text(encoding="utf-8", errors="replace"))
            except OSError:
                pass
        for claim in sorted(set(token for token in CLAIM_TOKEN_RE.findall("\n".join(texts)) if token in valid_claims)):
            _edge(edges, node_id, f"claim:{claim}", "tests", source="exact_claim_id_in_experiment_record")
    return result


def _rh_nodes(root: Path, valid_claims: set[str], sources: dict[str, dict[str, Any]], edges: list[dict[str, Any]]) -> list[dict[str, Any]]:
    path = root / "rh" / "catalog" / "rh_semantic_catalog.json"
    catalog = _read_json(path)
    if not catalog:
        return []
    source_id = _source(root, path, sources)
    result: list[dict[str, Any]] = []
    for record in catalog.get("records", []):
        if not isinstance(record, dict):
            continue
        relpath = record.get("path")
        if not isinstance(relpath, str) or not relpath:
            continue
        node_id = f"rh:{relpath}"
        attrs = {key: record.get(key) for key in ("outcome", "failure_class", "family", "mechanism", "question", "round") if key in record}
        attrs.update({"firewall": "rh-catalog", "source_refs": [source_id] if source_id else [], "reviewed": False})
        result.append(_node(node_id, "rh", relpath, attrs, str(record.get("mechanism") or "")))
        for claim in record.get("ledger_ids", []):
            if isinstance(claim, str) and claim in valid_claims:
                _edge(edges, node_id, f"claim:{claim}", "references", source="rh_semantic_catalog")
                if record.get("outcome") == "KILLED":
                    _edge(edges, node_id, f"claim:{claim}", "kills", source="rh_semantic_catalog")
    return result


def _scorecard_nodes(root: Path, valid_claims: set[str], sources: dict[str, dict[str, Any]], edges: list[dict[str, Any]]) -> list[dict[str, Any]]:
    path = root / "experiments" / "evidence_scorecard.json"
    payload = _read_json(path)
    if not payload:
        return []
    source_id = _source(root, path, sources)
    experiment_names = [
        directory.name for directory in sorted((root / "experiments").iterdir(), key=lambda p: p.name)
        if directory.is_dir() and directory.name != "theory-contracts" and len(directory.name) >= 5
    ] if (root / "experiments").is_dir() else []
    result = []
    for row in payload.get("rows", []):
        if not isinstance(row, dict):
            continue
        domain = str(row.get("domain") or "?").strip()
        observable = str(row.get("observable") or "").strip()
        slug = _slugify(observable)
        node_id = f"scorecard:{domain}:{slug}"
        attrs = {key: row.get(key) for key in ("domain", "observable", "status", "stage", "evidence_class", "pull_sigma", "kill_condition", "source") if key in row}
        attrs.update({"source_refs": [source_id] if source_id else [], "firewall": "experiments", "reviewed": False})
        result.append(_node(node_id, "scorecard", f"{domain}: {observable[:60]}", attrs))
        for claim in sorted(set(token for token in CLAIM_TOKEN_RE.findall(json.dumps(row, ensure_ascii=False)) if token in valid_claims)):
            _edge(edges, node_id, f"claim:{claim}", "tests", source="scorecard_exact_claim_token")
        source_text = str(row.get("source") or "").lower().replace("_", "-")
        for experiment in experiment_names:
            if experiment in source_text:
                _edge(edges, node_id, f"experiment:{experiment}", "sourced_from", source="scorecard_source_field")
    return result


def _failure_nodes(root: Path, known_ids: set[str], valid_claims: set[str], edges: list[dict[str, Any]]) -> list[dict[str, Any]]:
    path = root / "verification" / "graph_v2" / "failure_scopes.json"
    payload = _read_json(path)
    if not payload:
        return []
    result = []
    for record in payload.get("records", []):
        if not isinstance(record, dict) or not isinstance(record.get("id"), str):
            continue
        node_id = f"failure:{record['id']}"
        attrs = {key: value for key, value in record.items() if key != "id"}
        attrs["reviewed"] = bool(record.get("provenance", {}).get("reviewed"))
        attrs["inventory_source"] = "verification/graph_v2/failure_scopes.json"
        result.append(_node(node_id, "failure", record["id"], attrs))
        for evidence in record.get("evidence", []):
            raw_target = str(evidence)
            if raw_target.startswith(("contract:", "claim:", "rh:")):
                target = raw_target
            elif raw_target in valid_claims:
                target = f"claim:{raw_target}"
            else:
                target = f"contract:{raw_target}"
            if target in known_ids:
                _edge(edges, node_id, target, "excludes", source="failure_scopes", reviewed=False)
    return result


def _plan_nodes(root: Path, sources: dict[str, dict[str, Any]], edges: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Keep the live operational plan in the sole Atlas without promoting proof status."""
    base=root/"docs/research/fundamental-plan"
    if not base.is_dir():
        return []
    progress=_read_json(base/"progress.json")
    if not isinstance(progress,dict):
        progress={}
    attrs={"path":"docs/research/fundamental-plan/progress.json",
           "revision":progress.get("revision"),"updated_at":progress.get("updated_at"),
           "next_task":progress.get("next_task"),"continuation":str(progress.get("continuation",""))[:4000],
           "proof_coverage":False,"authority":"current user research plan; navigation only"}
    nodes=[_node("plan:fundamental","research_plan","Current fundamental research plan",attrs,
                 "Live plan and continuation; no mathematical or physical status promotion.")]
    for path in sorted(base.iterdir()):
        if not path.is_file() or path.suffix.lower()not in{".md",".json",".txt"}:
            continue
        sid=_source(root,path,sources)
        if sid is None:
            continue
        rel=path.relative_to(root).as_posix()
        nid="research_plan_document:"+rel
        nodes.append(_node(nid,"research_plan_document",path.name,
                           {"path":rel,"source_refs":[sid],"proof_coverage":False,
                            "scope":"local research operation; original scientific sources remain authoritative"},
                           "Canonical plan source, available through hash-checked source reads."))
        _edge(edges,"plan:fundamental",nid,"has_plan_document",proof_coverage=False)
    return nodes


def build(root: Path, valid_claims: set[str]) -> dict[str, Any]:
    """Return a deterministic in-memory source inventory; perform no writes."""
    root = root.resolve()
    sources: dict[str, dict[str, Any]] = {}
    edges: list[dict[str, Any]] = []
    nodes = []
    nodes.extend(_plan_nodes(root, sources, edges))
    nodes.extend(_contract_nodes(root, valid_claims, sources, edges))
    nodes.extend(_experiment_nodes(root, valid_claims, sources, edges))
    nodes.extend(_rh_nodes(root, valid_claims, sources, edges))
    nodes.extend(_scorecard_nodes(root, valid_claims, sources, edges))
    known_ids = {node["id"] for node in nodes} | {f"claim:{claim}" for claim in valid_claims}
    nodes.extend(_failure_nodes(root, known_ids, valid_claims, edges))
    # Eliminate duplicate source-backed relations and keep output stable.
    unique_edges = {}
    for edge in edges:
        key = (edge["src"], edge["dst"], edge["type"])
        unique_edges.setdefault(key, edge)
    return {
        "nodes": sorted(nodes, key=lambda node: node["id"]),
        "edges": sorted(unique_edges.values(), key=lambda edge: (edge["src"], edge["type"], edge["dst"])),
        "sources": sorted(sources.values(), key=lambda source: source["path"]),
    }


def check_against_corpus(root: Path, inventory: dict[str, Any], corpus_path: Path | None = None) -> dict[str, Any]:
    """Check imported source hashes and record IDs against the last Atlas corpus."""
    corpus_path = corpus_path or root / "semantic_explorer" / "data" / "corpus.json"
    try:
        corpus = json.loads(corpus_path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        return {"ok": False, "error": f"cannot read Atlas corpus: {type(exc).__name__}", "stale_sources": [], "missing_records": []}
    old_sources = {source.get("path"): source.get("sha256") for source in corpus.get("sources", []) if isinstance(source, dict)}
    stale_sources = sorted(
        source["path"] for source in inventory["sources"]
        if old_sources.get(source["path"]) != source["sha256"]
    )
    old_records = {node.get("id") for node in corpus.get("nodes", []) if isinstance(node, dict)}
    missing_records = sorted(node["id"] for node in inventory["nodes"] if node["id"] not in old_records)
    ledger_path = root / "verification" / "status_ledger.csv"
    graph_path = root / "semantic_explorer" / "data" / "knowledge.json"
    stale_claims: list[str] = []
    ledger_check_available = ledger_path.is_file() and graph_path.is_file()
    if ledger_check_available:
        try:
            with ledger_path.open(encoding="utf-8-sig", newline="") as handle:
                current_ledger = {row.get("claim_id", ""): row for row in csv.DictReader(handle) if row.get("claim_id")}
            graph = json.loads(graph_path.read_text(encoding="utf-8"))
            graph_claims = {
                node.get("id", "").removeprefix("claim:"): node.get("attrs", {}).get("ledger")
                for node in graph.get("nodes", [])
                if isinstance(node, dict) and node.get("type") == "claim"
            }
            stale_claims = sorted(
                claim_id for claim_id in (set(current_ledger) | set(graph_claims))
                if current_ledger.get(claim_id) != graph_claims.get(claim_id)
            )
        except (OSError, ValueError):
            ledger_check_available = False
    return {
        "ok": not stale_sources and not missing_records and (not ledger_check_available or not stale_claims),
        "stale_sources": stale_sources,
        "missing_records": missing_records,
        "stale_claims": stale_claims,
        "ledger_check_available": ledger_check_available,
        "checked_sources": len(inventory["sources"]),
        "checked_records": len(inventory["nodes"]),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Check source inventory freshness against the saved Atlas corpus.")
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--check", action="store_true", help="check current source hashes and inventory IDs; do not write files")
    args = parser.parse_args()
    root = args.repo.resolve()
    ledger = root / "verification" / "status_ledger.csv"
    try:
        with ledger.open(encoding="utf-8-sig", newline="") as handle:
            valid_claims = {row["claim_id"].strip() for row in csv.DictReader(handle) if row.get("claim_id")}
    except OSError as exc:
        print(json.dumps({"ok": False, "error": f"cannot read ledger: {type(exc).__name__}"}))
        return 1
    inventory = build(root, valid_claims)
    if args.check:
        result = check_against_corpus(root, inventory)
        if not result.get("ledger_check_available"):
            result["ok"] = False
        summary = {key: value for key, value in result.items() if key not in ("stale_sources", "missing_records", "stale_claims")}
        summary["stale_source_count"] = len(result.get("stale_sources", []))
        summary["missing_record_count"] = len(result.get("missing_records", []))
        summary["stale_claim_count"] = len(result.get("stale_claims", []))
        summary["stale_source_samples"] = result.get("stale_sources", [])[:10]
        summary["missing_record_samples"] = result.get("missing_records", [])[:10]
        summary["stale_claim_samples"] = result.get("stale_claims", [])[:10]
        print(json.dumps(summary, ensure_ascii=False, sort_keys=True))
        return 0 if result["ok"] else 1
    parser.error("--check is required; this module does not write graph files")
    return 2


if __name__ == "__main__":
    sys.exit(main())
