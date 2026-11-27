---
name: tfpt-knowledge
description: Query the TFPT Knowledge Atlas for papers, exact sources, Python/Lean/Wolfram evidence and recorded research decisions. Use it to check ideas and to continue named tasks; in Exploration mode it supplies details and checks ideas; it does not set the next step.
---

The `tfpt-knowledge` MCP server is the sole theory atlas. The ledger
(`verification/status_ledger.csv`) and original sources own claim status; Atlas
links, graph paths and AI annotations (`reviewed:false`) are navigation, not proof.

## Mode

Read the repository `AGENTS.md` and pick the mode first.

- **Exploration:** start from `docs/research/TFPT-STAND.md` and
  `docs/research/EXPLORATION.md`. Use the Atlas for details on what works and for
  checking each idea, not to continue a contract chain: has it been tried
  (`search_knowledge`, `list_attempts`), which exclusion applies and under which
  exact premises. A scoped exclusion only blocks the idea if its premises really
  match; a different premise makes the idea legitimate.
- **Original source routes:** many records end with "the original R113 call stays
  blocked". That is scoped to the original source of those contracts. It is never
  a gate for another origin model or for F1, and never the default next step
  (user decision of 2026-11-23, `AGENTS.md`).
- **Before calling anything new or open:** `theory_overview` carries
  `ledger_areas` (every active ledger claim counted by area and class; the full
  list is `docs/research/LEDGER-KARTE.md`). `search_knowledge` takes German or
  English questions, drops function words, weighs rare concepts higher and lists
  ledger claims, reviewed findings and theory-map aspects (`map:…`) before
  sub-records and archived paper passages. Check both before the verdict.
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

For an existing problem, call
`research_state(query=problem ID, editable=true, record_ids=[])` for its revision.
To correct existing records, pass their exact IDs in `record_ids`. After reviewing
the original sources, call `record_research_update(expected_revision=…, patch=…, reason=…)`
with the `problem_id` and only the new or changed bilingual records. Use `problem`
for creation or explicit full replacement. Reload the affected records on a
revision conflict. Never hand-edit generated graph data. Exploration results stay in
`experiments/exploration/` until a Continuation session decides to adopt them.

## Interface checks

An old MCP process can expose obsolete tools: reconnect
(`semantic_explorer/run.sh mcp`). Web explorer: `semantic_explorer/run.sh serve`
(http://127.0.0.1:8731). Freshness: `python3 semantic_explorer/refresh.py --check`.
Detailed historic routing: `references/` (load only the section you need).
