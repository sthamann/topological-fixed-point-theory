#!/usr/bin/env python3
"""Build and validate a source-grounded second-pass semantic review packet.

This utility does not call a model and does not change graph annotations. Its
output is a proposal packet for a separate reviewer; source excerpts are only
included after repository-boundary and SHA-256 checks.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

try:
    import extract
except ModuleNotFoundError:  # package import for tests
    from . import extract


PROMPT = """Review the supplied machine-generated semantic annotation against only the attached original-source excerpts and static observations. Return JSON matching the schema exactly. Treat every finding as a proposal, never as an authoritative correction or claim verdict. Do not infer that a script was executed or that a mathematical/physical claim is proven.

Prioritize concrete source disagreements and omitted dependencies. In particular, distinguish a callable argument passed at a call site from a value the callable actually reads; check the static callable_audit. Distinguish a selected/set input from an internally computed helper, while recognizing that a computed expression can still encode a selected model choice. Do not call syntactic constants proof gaps by themselves.

Each issue must identify one exact annotation field path and include at least one short verbatim quote from an attached source excerpt. The quoted text must be copied exactly, with a valid source_id. If a source is unavailable, stale, or insufficient, report that in coverage and do not assert a contradiction from it. A source citation establishes what that source says, not that the source's scientific claim is true. Keep a clean distinction between code-level facts, source/ledger statements, and unresolved interpretation.

Do not rewrite the annotation, change ledger status, promote or demote a claim, or imply human review. Use `no_findings` only when the selected fields and attached evidence were checked and no source-supported issue was found; otherwise use `incomplete` when coverage is insufficient."""


RESPONSE_SCHEMA: dict[str, Any] = {
    "type": "object",
    "additionalProperties": False,
    "required": ["schema_version", "claim_id", "review_state", "reviewed", "issues", "coverage", "notice"],
    "properties": {
        "schema_version": {"type": "integer", "const": 1},
        "claim_id": {"type": "string"},
        "review_state": {"type": "string", "enum": ["reviewed_with_findings", "no_findings", "incomplete"]},
        "reviewed": {"type": "boolean", "const": False},
        "issues": {
            "type": "array",
            "items": {
                "type": "object",
                "additionalProperties": False,
                "required": ["field_path", "code", "severity", "finding", "suggested_correction", "evidence"],
                "properties": {
                    "field_path": {"type": "string"},
                    "code": {"type": "string", "enum": ["source_contradiction", "missing_dependency", "scope_mismatch", "unsupported_interpretation", "stale_reference", "other"]},
                    "severity": {"type": "string", "enum": ["critical", "high", "medium", "low"]},
                    "finding": {"type": "string"},
                    "suggested_correction": {"type": "string"},
                    "evidence": {
                        "type": "array",
                        "minItems": 1,
                        "items": {
                            "type": "object",
                            "additionalProperties": False,
                            "required": ["source_id", "relation", "quote"],
                            "properties": {
                                "source_id": {"type": "string"},
                                "relation": {"type": "string", "enum": ["contradicts", "supports", "limits"]},
                                "quote": {"type": "string", "minLength": 1},
                            },
                        },
                    },
                },
            },
        },
        "coverage": {
            "type": "object",
            "additionalProperties": False,
            "required": ["fields_checked", "source_ids_checked", "unavailable_sources", "limitations"],
            "properties": {
                "fields_checked": {"type": "array", "items": {"type": "string"}},
                "source_ids_checked": {"type": "array", "items": {"type": "string"}},
                "unavailable_sources": {"type": "array", "items": {"type": "string"}},
                "limitations": {"type": "array", "items": {"type": "string"}},
            },
        },
        "notice": {"type": "string", "const": "Unreviewed second-pass proposal; not an authoritative correction or proof verdict."},
    },
}


REVIEW_FIELDS = (
    "title_de",
    "overview_de",
    "interpretation_de",
    "input_audit",
    "candidate_review",
    "script_analysis",
    "results",
    "proof_gaps",
    "derivations",
    "relationships",
    "unknowns",
)


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def walk_refs(value: Any) -> set[str]:
    refs: set[str] = set()
    if isinstance(value, dict):
        source_refs = value.get("source_refs")
        if isinstance(source_refs, list):
            refs.update(x for x in source_refs if isinstance(x, str))
        for child in value.values():
            refs.update(walk_refs(child))
    elif isinstance(value, list):
        for child in value:
            refs.update(walk_refs(child))
    return refs


def canonical_path(repo_root: Path, relative: str) -> Path:
    candidate = (repo_root / relative).resolve(strict=True)
    candidate.relative_to(repo_root.resolve(strict=True))
    return candidate


def excerpt_source(source: dict[str, Any], repo_root: Path, max_lines: int) -> tuple[dict[str, Any] | None, str | None]:
    sid = source.get("id", "")
    relative = source.get("path", "")
    if not isinstance(relative, str) or not relative:
        return None, "source has no repository-relative path"
    try:
        path = canonical_path(repo_root, relative)
        raw = path.read_bytes()
    except (OSError, ValueError) as exc:
        return None, f"source unavailable or outside repository ({type(exc).__name__})"

    expected_hash = source.get("sha256")
    if not isinstance(expected_hash, str) or not expected_hash:
        return None, "source record has no SHA-256 for freshness verification"
    if sha256(raw) != expected_hash:
        return None, "source SHA-256 differs from graph record"

    start = max(1, int(source.get("line_start", 1) or 1))
    end = max(start, int(source.get("line_end", start) or start))
    try:
        if isinstance(source.get("text"), str):
            text = source["text"]
        else:
            text = raw.decode("utf-8")
        lines = text.splitlines()
    except UnicodeDecodeError:
        return None, "binary source lacks a stored text excerpt"

    # Stored excerpts already correspond to the source record's line range.
    # Otherwise, source ranges index the complete UTF-8 file.
    if isinstance(source.get("text"), str):
        selected = lines[:max_lines]
        shown_start = start
        truncated = len(lines) > max_lines
    else:
        selected = lines[start - 1 : end]
        shown_start = start
        truncated = len(selected) > max_lines
        selected = selected[:max_lines]
    if not selected:
        return None, "source range has no readable lines"
    record = {
        "source_id": sid,
        "path": relative,
        "sha256": expected_hash,
        "sha256_verified": True,
        "line_start": shown_start,
        "line_end": shown_start + len(selected) - 1,
        "page": source.get("page"),
        "truncated": truncated,
        "text": "\n".join(selected),
    }
    return record, None


def annotation_payload(attrs: dict[str, Any]) -> tuple[dict[str, Any] | None, bool]:
    """Read either a direct extraction or the semantic payload of a priority overlay."""
    candidates = []
    direct = attrs.get("semantic")
    if isinstance(direct, dict):
        candidates.append((direct, attrs.get("enrichment_state")))
    priority = attrs.get("priority_annotation")
    if isinstance(priority, dict) and isinstance(priority.get("semantic"), dict):
        candidates.append((priority["semantic"], priority.get("enrichment_state")))
    populated = [
        (semantic, state) for semantic, state in candidates
        if any(bool(semantic.get(field)) for field in REVIEW_FIELDS)
    ]
    if not populated:
        return None, False
    # A merged context graph keeps base and priority annotations separate. The
    # priority annotation is the intended overlay when it is itself populated.
    chosen = populated[-1]
    return chosen[0], chosen[1] == "enriched" if chosen[1] is not None else True


def review_graph_paths(graph_path: Path) -> list[Path]:
    """Return the requested graph plus sibling base/priority graphs for fallback."""
    requested = graph_path.resolve()
    candidates = [requested]
    if requested.name == "context.json":
        candidates.extend([requested.with_name("knowledge.json"), requested.parent / "priority" / "knowledge.json"])
    elif requested.parent.name == "priority":
        candidates.append(requested.parent.parent / "knowledge.json")
    else:
        candidates.append(requested.parent / "priority" / "knowledge.json")
    return list(dict.fromkeys(path for path in candidates if path.is_file()))


def location_document_sources(ledger: dict[str, Any] | None, repo_root: Path, claim: str) -> list[dict[str, Any]]:
    """Resolve explicit TeX location tags to narrow, hash-checked source records."""
    if not isinstance(ledger, dict):
        return []
    locations = [x.strip() for x in str(ledger.get("location", "")).split(";") if x.strip()]
    if not locations:
        return []
    sources: list[dict[str, Any]] = []
    for tag in locations:
        # Tags such as tfpt_1/tfpt_2 name the primary theory volumes. Resolve
        # only unambiguous root TeX files; never walk archived or generated trees.
        matches = sorted(p for p in repo_root.glob(tag + "*.tex") if p.is_file())
        if len(matches) != 1:
            continue
        path = matches[0]
        try:
            raw = path.read_bytes()
            text = raw.decode("utf-8")
        except (OSError, UnicodeDecodeError):
            continue
        lines = text.splitlines()
        ranges = extract.location_claim_ranges(ledger, tag, lines)
        for start, end in ranges:
            rel = str(path.resolve().relative_to(repo_root.resolve()))
            digest = sha256(raw)
            sid = "source:" + sha256(f"{rel}:{digest}:{start}:{end}".encode())[:20]
            sources.append({"id": sid, "path": rel, "sha256": digest,
                "line_start": start, "line_end": end,
                "text": "\n".join(lines[start - 1:end])})
    return sources


def build_packet(graph_path: Path, repo_root: Path, claim_id: str, max_sources: int = 32, max_lines: int = 160) -> dict[str, Any]:
    full_id = claim_id if claim_id.startswith("claim:") else f"claim:{claim_id}"
    graphs = [json.loads(path.read_text(encoding="utf-8")) for path in review_graph_paths(graph_path)]
    if not graphs:
        raise ValueError(f"review graph does not exist: {graph_path}")
    indexed_graphs = []
    for graph in graphs:
        nodes = {node["id"]: node for node in graph.get("nodes", []) if isinstance(node, dict) and isinstance(node.get("id"), str)}
        claim_node = nodes.get(full_id)
        if not claim_node:
            continue
        attrs = claim_node.get("attrs", {})
        semantic, current = annotation_payload(attrs if isinstance(attrs, dict) else {})
        if semantic is None or not current:
            continue
        priority = attrs.get("priority_annotation", {}) if isinstance(attrs, dict) else {}
        provenance = (priority.get("provenance", {}) if isinstance(priority, dict)
                      and semantic is priority.get("semantic") else attrs.get("provenance", {}))
        stamp = str(provenance.get("created_at", "") if isinstance(provenance, dict) else "")
        indexed_graphs.append((stamp, graph, nodes, claim_node, semantic))
    if not indexed_graphs:
        present = any(any(n.get("id") == full_id for n in g.get("nodes", []) if isinstance(n, dict)) for g in graphs)
        if present:
            raise ValueError(f"claim has no source-current enriched semantic annotation: {full_id}")
        raise ValueError(f"claim node not found in requested/base/priority graphs: {full_id}")
    _, graph, nodes, claim, semantic = max(indexed_graphs, key=lambda item: item[0])
    attrs = claim.get("attrs", {})
    selected_semantic = {key: semantic[key] for key in REVIEW_FIELDS if key in semantic}
    if not any(bool(value) for value in selected_semantic.values()):
        raise ValueError(f"claim semantic annotation is empty: {full_id}")

    refs = sorted(walk_refs(semantic))
    source_index = {source.get("id"): source for candidate_graph in graphs for source in candidate_graph.get("sources", []) if isinstance(source, dict)}
    unavailable: list[dict[str, str]] = []
    excerpts: list[dict[str, Any]] = []
    for sid in refs[:max_sources]:
        source = source_index.get(sid)
        if source is None:
            unavailable.append({"source_id": sid, "reason": "source ID is absent from graph"})
            continue
        excerpt, problem = excerpt_source(source, repo_root, max_lines)
        if problem:
            unavailable.append({"source_id": sid, "reason": problem})
        elif excerpt:
            excerpts.append(excerpt)
    if len(refs) > max_sources:
        unavailable.extend({"source_id": sid, "reason": "source omitted by packet limit"} for sid in refs[max_sources:])

    script_ids = {
        item.get("script_id")
        for group in semantic.get("script_analysis", [])
        if isinstance(group, dict)
        for item in [group]
        if isinstance(item.get("script_id"), str)
    }
    for edge in graph.get("edges", []):
        if edge.get("src") == full_id and edge.get("type") in {"uses_script", "checks", "verified_by", "has_script", "implemented_check"}:
            script_ids.add(edge.get("dst"))
        elif edge.get("dst") == full_id and edge.get("type") in {"used_by", "checks", "verified_by", "has_script", "implemented_check"}:
            script_ids.add(edge.get("src"))
    static_scripts = []
    for sid in sorted(x for x in script_ids if isinstance(x, str)):
        script_node = nodes.get(sid)
        if not script_node or script_node.get("type") != "script":
            continue
        static = script_node.get("attrs", {}).get("static_inputs", {})
        static_scripts.append({
            "script_id": sid,
            "path": script_node.get("attrs", {}).get("path"),
            "callable_audit": static.get("callable_audit", []),
            "syntactic_input_candidates": static.get("syntactic_input_candidates", []),
            "checks": static.get("checks", []),
            "notice": static.get("notice", "Syntactic candidates are not proof gaps."),
        })

    for source in location_document_sources(attrs.get("ledger"), repo_root, full_id):
        excerpt, problem = excerpt_source(source, repo_root, max_lines)
        if problem:
            unavailable.append({"source_id": source["id"], "reason": problem})
        elif excerpt and all(item.get("source_id") != excerpt["source_id"] for item in excerpts):
            excerpts.append(excerpt)
    if not excerpts:
        raise ValueError(f"claim has no readable, hash-verified source excerpts: {full_id}")

    full_annotation_sha256 = sha256(json.dumps(semantic, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode())
    semantic_fingerprint_sha256 = sha256("\n".join(str(semantic.get(k, "")) for k in ("title_de", "overview_de", "interpretation_de")).encode())
    return {
        "schema_version": 1,
        "claim_id": full_id,
        "reviewed": False,
        "proposal_only": True,
        "ledger": attrs.get("ledger") if isinstance(attrs, dict) else None,
        "semantic": selected_semantic,
        "source_annotation_sha256": full_annotation_sha256,
        "source_semantic_sha256": semantic_fingerprint_sha256,
        "static_script_observations": static_scripts,
        "source_excerpts": excerpts,
        "unavailable_sources": unavailable,
        "review_prompt": PROMPT,
        "response_schema": RESPONSE_SCHEMA,
        "notice": "Second-pass review is a proposal. It cannot overwrite ledger status, semantic annotations, or source files.",
    }


def validate_response(response: dict[str, Any], packet: dict[str, Any]) -> list[str]:
    errors: list[str] = []
    if not isinstance(response, dict):
        return ["response must be a JSON object"]
    top_fields = {"schema_version", "claim_id", "review_state", "reviewed", "issues", "coverage", "notice"}
    if set(response) != top_fields:
        errors.append("response must contain exactly the schema's top-level fields")
    if isinstance(response.get("schema_version"), bool) or response.get("schema_version") != 1:
        errors.append("schema_version must equal 1")
    if response.get("claim_id") != packet.get("claim_id"):
        errors.append("claim_id does not match packet")
    if response.get("reviewed") is not False:
        errors.append("reviewed must remain false")
    if not isinstance(response.get("claim_id"), str) or not response["claim_id"]:
        errors.append("claim_id must be a nonempty string")
    if not isinstance(response.get("review_state"), str) or response.get("review_state") not in {"reviewed_with_findings", "no_findings", "incomplete"}:
        errors.append("review_state is not allowed")
    if response.get("notice") != RESPONSE_SCHEMA["properties"]["notice"]["const"]:
        errors.append("notice must preserve the proposal-only wording")
    coverage = response.get("coverage")
    if not isinstance(coverage, dict) or set(coverage) != {"fields_checked", "source_ids_checked", "unavailable_sources", "limitations"}:
        errors.append("coverage must contain exactly the schema's required fields")
        coverage = None
    elif any(not isinstance(coverage.get(field), list) or any(not isinstance(x, str) for x in coverage[field]) for field in coverage):
        errors.append("coverage fields must be arrays of strings")
        coverage = None
    if coverage is not None:
        expected_fields = sorted(packet.get("semantic", {}))
        allowed_context_fields = {"ledger", "static_script_observations"}
        def selected_field(value: str) -> str | None:
            # Accept a precise path beneath the selected annotation root, in
            # dotted or JSON-pointer form, while rejecting unrelated fields.
            normalized = value.replace("/", ".").replace("[", ".").replace("]", "")
            parts = [part for part in normalized.split(".") if part]
            matches = [field for field in expected_fields if field in parts]
            return matches[-1] if len(matches) == 1 else None
        expected_sources = sorted(s.get("source_id") for s in packet.get("source_excerpts", []) if isinstance(s, dict) and isinstance(s.get("source_id"), str))
        expected_unavailable = sorted(s.get("source_id") for s in packet.get("unavailable_sources", []) if isinstance(s, dict) and isinstance(s.get("source_id"), str))
        checked_fields = coverage.get("fields_checked", [])
        checked_sources = coverage.get("source_ids_checked", [])
        canonical_checked = [selected_field(x) for x in checked_fields]
        if len(set(checked_fields)) != len(checked_fields) or any(x is None and not any(raw == f or raw.startswith(f + ".") or raw.startswith(f + "[") for f in allowed_context_fields) for raw, x in zip(checked_fields, canonical_checked)):
            errors.append("coverage.fields_checked contains an unselected field or duplicate: " + repr(checked_fields))
        if len(set(checked_sources)) != len(checked_sources) or not set(checked_sources).issubset(expected_sources):
            errors.append("coverage.source_ids_checked must name attached source excerpts without duplicates")
        if response.get("review_state") == "no_findings" and (any(x is None and not any(raw == f or raw.startswith(f + ".") or raw.startswith(f + "[") for f in allowed_context_fields) for raw, x in zip(checked_fields, canonical_checked)) or sorted({x for x in canonical_checked if x is not None}) != expected_fields):
            errors.append("no_findings requires coverage of every selected annotation field")
        if response.get("review_state") == "no_findings" and sorted(checked_sources) != expected_sources:
            errors.append("no_findings requires coverage of every attached source excerpt")
        if sorted(set(coverage.get("unavailable_sources", []))) != expected_unavailable:
            errors.append("coverage.unavailable_sources must match packet source limitations")
        if response.get("review_state") == "no_findings" and (expected_unavailable or not expected_sources):
            errors.append("no_findings is not allowed when evidence is unavailable or no sources are attached")
    valid_sources = {source["source_id"]: source["text"] for source in packet.get("source_excerpts", [])}
    issue_fields = {"field_path", "code", "severity", "finding", "suggested_correction", "evidence"}
    evidence_fields = {"source_id", "relation", "quote"}
    issue_codes = {"source_contradiction", "missing_dependency", "scope_mismatch", "unsupported_interpretation", "stale_reference", "other"}
    severities = {"critical", "high", "medium", "low"}
    issues = response.get("issues")
    if not isinstance(issues, list):
        errors.append("issues must be an array")
        issues = []
    def valid_field_path(value: Any) -> bool:
        if not isinstance(value, str):
            return False
        normalized = value.replace("/", ".").replace("[", ".").replace("]", "")
        parts = [part for part in normalized.split(".") if part]
        matches = [field for field in packet.get("semantic", {}) if field in parts]
        return len(matches) == 1
    for index, issue in enumerate(issues):
        if not isinstance(issue, dict) or set(issue) != issue_fields:
            errors.append(f"issues[{index}] does not match the strict issue schema")
            continue
        if not isinstance(issue.get("code"), str) or issue.get("code") not in issue_codes or not isinstance(issue.get("severity"), str) or issue.get("severity") not in severities:
            errors.append(f"issues[{index}] has an invalid code or severity")
        if any(not isinstance(issue.get(key), str) or not issue[key] for key in ["field_path", "finding", "suggested_correction"]):
            errors.append(f"issues[{index}] has an empty or non-string required text field")
        if not valid_field_path(issue.get("field_path")):
            errors.append(f"issues[{index}] field_path does not identify a selected annotation field: {issue.get('field_path')!r}")
        if not isinstance(issue.get("evidence"), list) or not issue["evidence"]:
            errors.append(f"issues[{index}] requires at least one evidence quote")
            continue
        for evidence in issue.get("evidence", []):
            if not isinstance(evidence, dict) or set(evidence) != evidence_fields:
                errors.append(f"issues[{index}] evidence does not match the strict evidence schema")
                continue
            if not isinstance(evidence.get("source_id"), str) or not evidence["source_id"] or not isinstance(evidence.get("quote"), str) or not evidence["quote"]:
                errors.append(f"issues[{index}] evidence source_id and quote must be nonempty strings")
                continue
            if not isinstance(evidence.get("relation"), str) or evidence.get("relation") not in {"contradicts", "supports", "limits"}:
                errors.append(f"issues[{index}] has an invalid evidence relation")
            sid, quote = evidence["source_id"], evidence["quote"]
            source_text = valid_sources.get(sid)
            if source_text is None:
                errors.append(f"issues[{index}] cites a source not present in the packet: {sid}")
            elif not quote or quote not in source_text:
                errors.append(f"issues[{index}] contains a quote not present verbatim in its cited excerpt")
    if response.get("review_state") == "no_findings" and response.get("issues"):
        errors.append("no_findings must have an empty issues array")
    if response.get("review_state") == "reviewed_with_findings" and not response.get("issues"):
        errors.append("reviewed_with_findings requires at least one issue")
    cited: set[str] = set()
    for issue in issues:
        if not isinstance(issue, dict) or not isinstance(issue.get("evidence"), list):
            continue
        for evidence in issue["evidence"]:
            if isinstance(evidence, dict) and isinstance(evidence.get("source_id"), str):
                cited.add(evidence["source_id"])
    if not cited.issubset(valid_sources):
        errors.append("one or more issue sources were not attached to the packet")
    return errors


def parse_model_response(response: dict[str, Any]) -> dict[str, Any]:
    if not isinstance(response, dict) or response.get("status") != "completed":
        raise ValueError("model response did not complete")
    texts = [c.get("text") for out in response.get("output", []) if isinstance(out, dict)
             for c in out.get("content", []) if isinstance(c, dict) and c.get("type") == "output_text"]
    if not texts or any(not isinstance(t, str) for t in texts):
        raise ValueError("model response contains no structured output")
    try:
        parsed = json.loads("".join(texts))
    except json.JSONDecodeError:
        raise ValueError("model response is not valid JSON") from None
    if not isinstance(parsed, dict):
        raise ValueError("model response must be a JSON object")
    return parsed


def run_review(packet: dict[str, Any], model: str, key_file: str | None = None) -> dict[str, Any]:
    semantic = packet.get("semantic")
    if not isinstance(semantic, dict) or not any(bool(semantic.get(field)) for field in REVIEW_FIELDS):
        raise ValueError("review preflight failed: packet has no populated semantic annotation")
    excerpts = packet.get("source_excerpts")
    if packet.get("unavailable_sources"):
        raise ValueError("review preflight failed: one or more referenced sources are unavailable")
    if not isinstance(excerpts, list) or not excerpts or not all(
        isinstance(item, dict) and isinstance(item.get("source_id"), str)
        and isinstance(item.get("path"), str) and isinstance(item.get("sha256"), str)
        and item.get("sha256_verified") is True
        and isinstance(item.get("text"), str) and item["text"].strip()
        for item in excerpts
    ):
        raise ValueError("review preflight failed: packet has no readable source excerpts")
    key = extract.resolve_key(key_file)
    request = {
        "model": model,
        "store": False,
        "instructions": PROMPT,
        "input": json.dumps(packet, ensure_ascii=False),
        "text": {"format": {"type": "json_schema", "name": "semantic_second_pass_review", "strict": True, "schema": RESPONSE_SCHEMA}},
        "max_output_tokens": 12000,
    }
    parsed = None
    errors: list[str] = []
    for attempt in range(2):
        response = extract.api("responses", request, key)
        parsed = parse_model_response(response)
        errors = validate_response(parsed, packet)
        if not errors:
            break
        if attempt == 0:
            request["input"] = json.dumps({
                "source_packet": packet,
                "invalid_response": parsed,
                "validation_errors": errors,
                "repair_instruction": "Return a corrected response conforming exactly to the schema. Copy every evidence quote verbatim from its cited excerpt; remove any issue you cannot support with an exact quote. Coverage paths must name selected annotation fields or attached static context. Do not add evidence or conclusions not present in the packet.",
            }, ensure_ascii=False)
    if errors:
        raise ValueError("model response validation failed: " + "; ".join(errors))
    semantic = packet.get("semantic", {})
    source_semantic_sha256 = packet.get("source_semantic_sha256") or hashlib.sha256("\n".join(str(semantic.get(k, "")) for k in ("title_de", "overview_de", "interpretation_de")).encode()).hexdigest()
    source_annotation_sha256 = packet.get("source_annotation_sha256") or hashlib.sha256(json.dumps(semantic, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode()).hexdigest()
    source_hashes = [
        {"source_id": source["source_id"], "path": source["path"], "sha256": source["sha256"]}
        for source in packet.get("source_excerpts", [])
        if isinstance(source, dict) and isinstance(source.get("source_id"), str) and isinstance(source.get("path"), str) and isinstance(source.get("sha256"), str)
    ]
    source_records = [
        {key: source[key] for key in ("source_id", "path", "sha256", "line_start", "line_end", "page") if key in source}
        for source in packet.get("source_excerpts", [])
        if isinstance(source, dict) and isinstance(source.get("source_id"), str)
        and isinstance(source.get("path"), str) and isinstance(source.get("sha256"), str)
    ]
    return {
        "proposal_only": True,
        "review_model": model,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "source_annotation_sha256": source_annotation_sha256,
        "source_semantic_sha256": source_semantic_sha256,
        "source_hashes": source_hashes,
        "sources": source_records,
        "review": parsed,
    }


def atomic_write_new(path: Path, value: dict[str, Any]) -> None:
    path = path.resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise FileExistsError("proposal output already exists; choose a new path")
    fd, temp_name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as stream:
            json.dump(value, stream, ensure_ascii=False, indent=2)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        # Refuse replacing a path created concurrently after the initial check.
        os.link(temp_name, path)
    finally:
        try:
            os.unlink(temp_name)
        except FileNotFoundError:
            pass


def main() -> int:
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--graph", type=Path, default=here / "data/knowledge.json")
    parser.add_argument("--repo-root", type=Path, default=here.parent)
    parser.add_argument("--claim", default="EM.FP.01")
    parser.add_argument("--max-sources", type=int, default=32)
    parser.add_argument("--max-lines", type=int, default=160)
    parser.add_argument("--out", type=Path, help="Write packet JSON here; stdout if omitted")
    parser.add_argument("--run-review", action="store_true", help="Call the model for a validated, proposal-only second pass")
    parser.add_argument("--model", default="gpt-5.6-terra", help="Model for --run-review (default: gpt-5.6-terra)")
    parser.add_argument("--key-file", help="Key file passed to extract.resolve_key; the key is never printed or saved")
    parser.add_argument("--validate-response", type=Path, help="Validate model response JSON against the generated packet")
    parser.add_argument("--packet", type=Path, help="Packet JSON for --validate-response")
    args = parser.parse_args()
    try:
        if args.validate_response:
            if not args.packet:
                parser.error("--packet is required with --validate-response")
            packet = json.loads(args.packet.read_text(encoding="utf-8"))
            response = json.loads(args.validate_response.read_text(encoding="utf-8"))
            errors = validate_response(response, packet)
            if errors:
                print("\n".join(errors), file=sys.stderr)
                return 1
            print("response quotes and source IDs validated")
            return 0
        packet = build_packet(args.graph, args.repo_root, args.claim, args.max_sources, args.max_lines)
        if args.run_review:
            if not args.out:
                parser.error("--out is required with --run-review for a separate proposal file")
            proposal = run_review(packet, args.model, args.key_file)
            atomic_write_new(args.out, proposal)
            print(f"validated proposal saved: {args.out}")
            return 0
        encoded = json.dumps(packet, ensure_ascii=False, indent=2)
        if args.out:
            args.out.parent.mkdir(parents=True, exist_ok=True)
            args.out.write_text(encoded + "\n", encoding="utf-8")
        else:
            print(encoded)
        return 0
    except (OSError, ValueError, json.JSONDecodeError, RuntimeError) as exc:
        print(f"review packet error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
