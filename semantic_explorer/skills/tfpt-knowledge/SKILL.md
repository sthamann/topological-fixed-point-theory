---
name: tfpt-knowledge
description: Query the TFPT Knowledge Atlas for papers, exact sources, Python/Lean/Wolfram evidence and recorded research decisions. Use it to check ideas and to continue named tasks; in Exploration mode it is a check after the idea phase, not the starting point.
---

The `tfpt-knowledge` MCP server is the sole theory atlas. The ledger
(`verification/status_ledger.csv`) and original sources own claim status; Atlas
links, graph paths and AI annotations (`reviewed:false`) are navigation, not proof.

## Mode

Read the repository `AGENTS.md` and pick the mode first.

- **Exploration:** write your own ideas following `docs/research/EXPLORATION.md`
  before any Atlas call. Then use the Atlas to check each idea: has it been tried
  (`search_knowledge`, `list_attempts`), which exclusion applies and under which
  exact premises. A scoped exclusion only blocks the idea if its premises really
  match; a different premise makes the idea legitimate.
- **Continuation:** `graph_status`, `theory_overview`, then
  `research_state(query, lang, max_chars=12000)` for the named task, together with
  the top of `PLAN.md` and the compact `progress.json`.

## Retrieval

- Papers: `list_papers`/`search_papers` → `paper_outline` → `paper_section` →
  `paper_evidence` → `read_source` (explicit `line_start`/`line_end`, ≤200 lines
  per call). Never pass a section ID to `paper_outline`.
- Prior work: `research_context`, `research_brief`, `compare_aspects` (2–6 IDs),
  `solution_chain`, `list_attempts`; open `get_dossier` / `audit_inputs` only for
  decisive items. For `packet_state="partial"` follow `budget.continuation`.
- Before repeating a recorded strategy: `assess_research_plan` with exact
  `problem_id`, `strategy_id`, `premise_ids`, and name the changed premise.
- Missing reviewed coverage is not evidence that a question is open; a missing
  match is not evidence that it was never tried.

## Discipline

Equal dimensions, spectra or group orders do not supply a physical map; check
group actions, charges, grading, adjoints, domains and normalization. Keep
conjunctive premises together; shared seeds give correlated, not independent
evidence. A working auxiliary model is not an origin derivation. Lean inventory
is static: require the statement, axioms and a kernel/build receipt. A
reproduction of a measured target without TFPT inputs or fit may be reported as
**Kandidat** (see `AGENTS.md`).

## Recording

After a source-supported result or correction: load
`research_state(query=problem ID, editable=true, max_chars=60000)`, update the
complete record and call `record_research_update(expected_revision=…)`. Reload on
conflict; never hand-edit generated graph data. Exploration results stay in
`experiments/exploration/` until a Continuation session decides to adopt them.

## Interface checks

An old MCP process can expose obsolete tools: reconnect
(`semantic_explorer/run.sh mcp`). Web explorer: `semantic_explorer/run.sh serve`
(http://127.0.0.1:8731). Freshness: `python3 semantic_explorer/refresh.py --check`.
Detailed historic routing: `references/` (load only the section you need).
