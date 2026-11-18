"""Read the canonical Atlas snapshot and project it for legacy consumers.

This is an adapter only: the Atlas remains the sole graph store.  The projected
items are short-lived, read-only view objects and are never persisted.
"""
from __future__ import annotations

import json
from pathlib import Path
import subprocess
from collections import Counter, defaultdict
from typing import Any


ATTEMPT_KIND_MAP = {
    "contract": "contract",
    "experiment": "experiment",
    "rh": "rh",
    "scorecard": "scorecard",
    "failure": "failure",
}


def load_atlas_export(repo_root: Path | None = None, *, timeout: int = 180) -> dict[str, Any]:
    """Return the canonical Atlas graph via its read-only `run.sh export` route.

    There is deliberately no fallback to a generated or cached legacy graph.
    """
    root = (repo_root or Path(__file__).resolve().parents[1]).resolve()
    command = [str(root / "semantic_explorer" / "run.sh"), "export"]
    try:
        result = subprocess.run(
            command, cwd=root, check=False, capture_output=True, text=True,
            encoding="utf-8", errors="replace", timeout=timeout,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise RuntimeError(f"TFPT Knowledge Atlas export failed: {exc}") from exc
    if result.returncode != 0:
        detail = result.stderr.strip()[-2000:]
        raise RuntimeError(f"TFPT Knowledge Atlas export exited {result.returncode}: {detail}")
    try:
        graph = json.loads(result.stdout)
    except json.JSONDecodeError as exc:
        raise RuntimeError(f"TFPT Knowledge Atlas export returned invalid JSON: {exc}") from exc
    if not isinstance(graph, dict) or not isinstance(graph.get("nodes"), list) or not isinstance(graph.get("edges"), list):
        raise RuntimeError("TFPT Knowledge Atlas export is missing its nodes/edges arrays")
    return graph


def _source_index(graph: dict[str, Any]) -> dict[str, dict[str, Any]]:
    sources = graph.get("sources")
    if isinstance(sources, dict):
        return {str(key): value for key, value in sources.items() if isinstance(value, dict)}
    if isinstance(sources, list):
        return {str(item["id"]): item for item in sources if isinstance(item, dict) and item.get("id")}
    return {}


def _resolved_sources(refs: Any, sources: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    if not isinstance(refs, list):
        refs = [refs] if refs else []
    result: list[dict[str, Any]] = []
    seen: set[tuple[str, int | None, int | None]] = set()
    for ref in refs:
        source = ref if isinstance(ref, dict) else sources.get(str(ref))
        if not isinstance(source, dict) or not source.get("path"):
            continue
        path = str(source["path"])
        if path.startswith("/") or ".." in Path(path).parts:
            continue
        item: dict[str, Any] = {"path": path}
        if source.get("line_start") is not None:
            item["line"] = source["line_start"]
            item["line_start"] = source["line_start"]
        if source.get("line_end") is not None:
            item["line_end"] = source["line_end"]
        if source.get("sha256"):
            item["sha256"] = source["sha256"]
        key = (path, item.get("line_start"), item.get("line_end"))
        if key not in seen:
            seen.add(key)
            result.append(item)
    return result


def _short_id(node: dict[str, Any], prefix: str) -> str:
    node_id = str(node.get("id") or "")
    if node_id.startswith(prefix):
        return node_id[len(prefix):]
    return node_id.split(":", 1)[-1]


def _node_path(node: dict[str, Any], sources: list[dict[str, Any]], root: Path) -> str:
    attrs = node.get("attrs") if isinstance(node.get("attrs"), dict) else {}
    path = attrs.get("path") or node.get("path")
    if path:
        candidate = Path(str(path))
        try:
            if candidate.is_absolute():
                path = str(candidate.resolve().relative_to(root))
            else:
                path = str(candidate)
            if ".." not in Path(path).parts:
                return path
        except (OSError, ValueError):
            pass
    return str(sources[0]["path"]) if sources else "semantic_explorer/data/knowledge.json"


def project_evidence_graph(graph: dict[str, Any], repo_root: Path | None = None) -> dict[str, Any]:
    """Project Atlas nodes into the compact item types consumed by TFPT Explorer.

    This view intentionally contains claims, registered scripts, attempts, and
    one inventory record per Lean source file. The complete Atlas graph remains
    available separately to relation lookups.
    """
    root = (repo_root or Path(__file__).resolve().parents[1]).resolve()
    source_index = _source_index(graph)
    source_nodes: dict[str, dict[str, Any]] = {}
    source_edges: list[dict[str, Any]] = []
    lean_by_path: dict[str, list[dict[str, Any]]] = {}
    node_by_id = {
        str(n.get("id") or ""): n for n in graph.get("nodes", []) if isinstance(n, dict)
    }

    for node in graph.get("nodes", []):
        if not isinstance(node, dict):
            continue
        node_id = str(node.get("id") or "")
        node_type = str(node.get("type") or "")
        attrs = node.get("attrs") if isinstance(node.get("attrs"), dict) else {}
        refs = _resolved_sources(attrs.get("source_refs") or attrs.get("sources"), source_index)
        path = _node_path(node, refs, root)

        if node_type == "lean_theorem" or node_type == "lean_module":
            lean_path = str(attrs.get("path") or path)
            if lean_path.startswith("/") or ".." in Path(lean_path).parts:
                continue
            lean_by_path.setdefault(lean_path, []).append(node)
            continue

        projected_type = node_type
        if node_type == "attempt":
            projected_type = ATTEMPT_KIND_MAP.get(str(attrs.get("attempt_kind") or ""), "")
            if not projected_type:
                continue
        elif node_type not in {"claim", "script"}:
            continue

        if node_type == "claim":
            ledger = attrs.get("ledger") if isinstance(attrs.get("ledger"), dict) else {}
            # Keep the complete Atlas attributes while flattening the original
            # ledger row for existing read-only consumers.
            attrs = {**attrs, **ledger}
            label = str(ledger.get("claim_id") or _short_id(node, "claim:") or node.get("label") or node_id)
            item_id = f"claim:{label}"
            path = "verification/status_ledger.csv"
        elif node_type == "script":
            raw = str(attrs.get("path") or node.get("label") or _short_id(node, "script:"))
            label = Path(raw).name.removesuffix(".py")
            item_id = f"script:{label}"
            path = raw if raw and not Path(raw).is_absolute() and ".." not in Path(raw).parts else path
        else:
            kind = projected_type
            node_id = str(node.get("id") or "")
            if node_id.startswith(f"{kind}:"):
                suffix = node_id[len(kind) + 1:]
            elif node_type == "attempt" and node_id.startswith("attempt:"):
                suffix = node_id[len("attempt:"):]
                if suffix.startswith(f"{kind}:"):
                    suffix = suffix[len(kind) + 1:]
            else:
                suffix = node_id.split(":", 1)[-1]
            if kind == "contract":
                # Contract IDs are stable handles; their human question is
                # descriptive metadata, never the identifier shown as title.
                label = suffix
                item_id = f"contract:{suffix}"
            elif kind in {"experiment", "rh", "scorecard", "failure"}:
                label = str(node.get("label") or attrs.get("question") or suffix)
                item_id = f"{kind}:{suffix}"
                if kind == "rh" and suffix and not Path(suffix).is_absolute() and ".." not in Path(suffix).parts:
                    candidate = root / suffix
                    if candidate.is_file():
                        # RH probe records are keyed by their exact source path;
                        # keep that identity so fresh audit overlays still join.
                        path = suffix
            else:
                label = str(node.get("label") or suffix)
                item_id = f"{kind}:{suffix}"

        item_attrs = attrs
        status = (attrs.get("status_raw") or attrs.get("verdict_raw") or attrs.get("status")
                  or attrs.get("outcome") or attrs.get("marker") or attrs.get("verdict_class") or "unclassified")
        verdict = attrs.get("verdict_raw") or attrs.get("verdict_class") or attrs.get("outcome") or "unclassified"
        title = label if projected_type in {"claim", "script", "contract", "failure"} else (
            attrs.get("title") or attrs.get("question") or node.get("label") or label
        )
        description = (node.get("summary") or attrs.get("summary_de") or attrs.get("what_web")
                       or attrs.get("statement") or attrs.get("mechanism") or attrs.get("canonical_status") or "")
        scope = attrs.get("scope") or attrs.get("canonical_status") or attrs.get("firewall") or ""
        if not refs:
            refs = [{"path": path}]
        item = {
            "id": item_id,
            "type": projected_type,
            "label": label,
            "title": str(title),
            "status": str(status),
            "verdict": str(verdict),
            "scope": str(scope),
            "description": str(description),
            "path": path,
            "claims": [label] if projected_type == "claim" else [],
            "scripts": [label] if projected_type == "script" else [],
            "sources": refs,
            "stage_ids": [],
            "last_run": None,
            "original_status": str(status),
            "status_origin": "atlas_ledger" if projected_type == "claim" else "atlas_context",
            "execution_status": "not_run_in_explorer" if projected_type in {"contract", "experiment"} else None,
            "attributes": item_attrs,
        }
        source_nodes[node_id] = item

    # Atlas also registers Lean source files that contain no indexed theorem or
    # module declaration. Preserve the complete source-file inventory rather
    # than silently dropping those files from the Explorer's file view.
    for source in source_index.values():
        path = str(source.get("path") or "")
        if not path.endswith(".lean") or path.startswith("/") or ".." in Path(path).parts:
            continue
        lean_by_path.setdefault(path, [])

    # Lean is a static declaration inventory, grouped by source file to match
    # the Explorer's file-oriented evidence cards. No kernel status is inferred.
    for path, declarations in sorted(lean_by_path.items()):
        ids = [str(n.get("id") or "") for n in declarations]
        source_by_path: dict[str, dict[str, Any]] = {}
        for declaration in declarations:
            attrs = declaration.get("attrs") or {}
            for source in _resolved_sources(attrs.get("source_refs"), source_index):
                source_by_path.setdefault(source["path"], {"path": source["path"], **({"sha256": source["sha256"]} if source.get("sha256") else {})})
        source_refs = list(source_by_path.values()) or [{"path": path}]
        source_nodes[f"lean:{path}"] = {
            "id": f"lean:{path}", "type": "lean", "label": path, "title": path,
            "status": "statically inventoried; kernel/build not checked",
            "verdict": "static_inventory_only",
            "scope": "Declaration inventory only; no Lean build or kernel check is implied.",
            "description": f"{len(declarations)} statically catalogued Lean declarations.",
            "path": path, "claims": [], "scripts": [], "sources": source_refs,
            "stage_ids": [], "last_run": None,
            "original_status": "static_inventory_only", "status_origin": "atlas_lean_inventory",
            "execution_status": "not_run_in_explorer",
            "attributes": {"declaration_count": len(declarations), "declaration_ids": ids,
                           "inventory_method": "static", "kernel_checked": False},
        }

    all_nodes = list(node_by_id.values())
    type_by_id = {node_id: str(n.get("type") or "") for node_id, n in node_by_id.items()}
    claims_by_node: dict[str, set[str]] = defaultdict(set)
    scripts_by_node: dict[str, set[str]] = defaultdict(set)
    for edge in graph.get("edges", []):
        if not isinstance(edge, dict):
            continue
        src, dst = str(edge.get("src") or ""), str(edge.get("dst") or "")
        for neighbor, related in ((dst, src), (src, dst)):
            if type_by_id.get(related) == "claim":
                claims_by_node[neighbor].add(_short_id(node_by_id.get(related, {}), "claim:"))
            if type_by_id.get(related) == "script":
                script_node = node_by_id.get(related, {})
                script_attrs = script_node.get("attrs") if isinstance(script_node.get("attrs"), dict) else {}
                scripts_by_node[neighbor].add(Path(str(script_attrs.get("path") or script_node.get("label") or related)).name.removesuffix(".py"))
    for original_id, item in source_nodes.items():
        if item["type"] != "script":
            item["claims"] = sorted(claims_by_node.get(original_id, set()))
        if item["type"] != "claim":
            item["scripts"] = sorted(scripts_by_node.get(original_id, set()))
        item["stage_ids"] = []

    # Preserve edge directions and only endpoint IDs; relation lookups use the
    # full canonical snapshot held by the server.
    for edge in graph.get("edges", []):
        if not isinstance(edge, dict):
            continue
        src, dst = str(edge.get("src") or ""), str(edge.get("dst") or "")
        if src in source_nodes and dst in source_nodes:
            source_edges.append({"src": src, "dst": dst, "type": edge.get("type"), "attrs": edge.get("attrs", {})})

    type_counts = Counter(item["type"] for item in source_nodes.values())
    attempt_counts = Counter(item["type"] for item in source_nodes.values() if item["type"] in ATTEMPT_KIND_MAP.values())
    meta = dict(graph.get("meta") or {})
    meta["projection_counts"] = dict(type_counts)
    meta["attempt_counts"] = dict(attempt_counts)
    meta["atlas_node_count"] = len(all_nodes)
    meta["atlas_edge_count"] = len(graph.get("edges", []))
    meta["lean_source_files"] = len(lean_by_path)
    return {"meta": meta, "nodes": list(source_nodes.values()), "edges": source_edges,
            "sources": graph.get("sources") or []}


__all__ = ["load_atlas_export", "project_evidence_graph"]
