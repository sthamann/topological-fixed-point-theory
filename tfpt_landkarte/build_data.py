#!/usr/bin/env python3
"""Exportiert den lebenden Atlasstand als Browserdaten für die TFPT Landkarte.

Liest nur: semantic_explorer/research/program.json, semantic_explorer/theory_map.json
und website/public/verification/predictions_frozen.json. Schreibt data/atlas.js.
Es werden keine Statuswerte verändert und keine Prüfungen ausgeführt.
"""
from __future__ import annotations

import datetime as _dt
import json
import pathlib

ROOT = pathlib.Path(__file__).resolve().parent.parent
OUT = pathlib.Path(__file__).resolve().parent / "data" / "atlas.js"
STATEMENT_LIMIT = 520


def de(value) -> str:
    if isinstance(value, dict):
        return (value.get("de") or value.get("en") or "").strip()
    return (value or "").strip()


def clip(text: str, limit: int = STATEMENT_LIMIT) -> str:
    text = " ".join(text.split())
    return text if len(text) <= limit else text[: limit - 1].rstrip() + "…"


def source_paths(sources) -> list[str]:
    paths = []
    for src in sources or []:
        if isinstance(src, dict):
            path = src.get("path") or src.get("source_path") or src.get("id") or ""
            if path and src.get("line_start"):
                path = f"{path}:{src['line_start']}"
        else:
            path = str(src)
        if path and path not in paths:
            paths.append(path)
    return paths[:4]


def main() -> None:
    program = json.loads((ROOT / "semantic_explorer/research/program.json").read_text())
    theory_map = json.loads((ROOT / "semantic_explorer/theory_map.json").read_text())
    frozen = json.loads(
        (ROOT / "website/public/verification/predictions_frozen.json").read_text()
    )

    problems = []
    for problem in program["problems"]:
        findings = [
            {
                "t": de(f.get("title")),
                "s": f.get("status") or "",
                "x": clip(de(f.get("statement"))),
                "g": clip(de(f.get("scope")), 320),
                "p": source_paths(f.get("sources")),
            }
            for f in problem.get("findings", [])
        ]
        strategies = [
            {
                "t": de(s.get("title")),
                "o": s.get("outcome") or "",
                "x": clip(de(s.get("explanation"))),
                "c": clip(de(s.get("changes_required")), 320),
            }
            for s in problem.get("strategies", [])
        ]
        tasks = [
            {
                "t": de(t.get("title")),
                "s": t.get("status") or "",
                "n": clip(de(t.get("next_action")), 320),
            }
            for t in problem.get("tasks", [])
        ]
        problems.append(
            {
                "id": problem["id"],
                "title": de(problem.get("title")),
                "q": clip(de(problem.get("question")), 600),
                "topics": [t.split(":", 1)[-1] for t in problem.get("topics", [])],
                "f": findings,
                "s": strategies,
                "k": tasks,
            }
        )

    patterns = [
        {
            "id": p["id"],
            "label": p.get("label_de", ""),
            "kind": p.get("relation_type", ""),
            "summary": p.get("summary_de", ""),
            "conditions": p.get("conditions", [])[:4],
            "counter": p.get("counterevidence", [])[:4],
        }
        for p in theory_map.get("patterns", [])
    ]

    map_nodes = [
        {
            "id": n["id"],
            "label": n.get("label_de", ""),
            "stage": n.get("stage", ""),
            "state": n.get("state", ""),
            "summary": clip(n.get("summary_de", ""), 600),
            "scope": clip(n.get("scope_note", ""), 400),
        }
        for n in theory_map.get("nodes", [])
    ]

    payload = {
        "generated": _dt.datetime.now().strftime("%Y-%m-%d %H:%M"),
        "revision": program.get("revision", ""),
        "problems": problems,
        "patterns": patterns,
        "mapNodes": map_nodes,
        "mapEdges": theory_map.get("edges", []),
        "frozen": {
            "freeze_date": frozen.get("freeze_date"),
            "predictions": frozen.get("predictions", []),
            "assigned": frozen.get("assigned_texture_values", []),
            "bands": frozen.get("conditional_bands", []),
        },
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(
        "window.ATLAS = " + json.dumps(payload, ensure_ascii=False, separators=(",", ":")) + ";\n"
    )
    findings = sum(len(p["f"]) for p in problems)
    strategies = sum(len(p["s"]) for p in problems)
    print(
        f"atlas.js: {len(problems)} Probleme, {findings} Befunde, {strategies} Strategien, "
        f"{OUT.stat().st_size / 1e6:.2f} MB, Revision {payload['revision']}"
    )


if __name__ == "__main__":
    main()
