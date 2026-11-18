---
name: tfpt-knowledge
description: Query the single current TFPT Knowledge Atlas for papers, exact sources, Python/Lean/Wolfram evidence and source-reviewed research decisions. Start fundamental research from the live shared plan, preserve established results and scoped exclusions, and trace actual prerequisites before proposing a new derivation.
---

Use the `tfpt-knowledge` MCP server for theory records and the separate
codebase-memory MCP graph for code discovery. The Atlas is a navigation and
decision tool; original sources and `verification/status_ledger.csv` own proof
scope and claim status. No import, graph edge or AI summary upgrades a verdict.

## Start from the current task and live state

Reconstruct the user's requested outcome. Read the required whole-system entry
map and applicable project rules once, preserving the connected compiler,
clocks, states, E8, flavor, alpha and physical completion goals. Current user
instructions take precedence over dated instructions and handoffs. For an Atlas
maintenance or audit request, repair and verify the requested system; do not
start a separate scientific continuation merely because a plan exists.

For fundamental research the canonical shared entry is the repository's
[PLAN.md](../../../docs/research/fundamental-plan/PLAN.md),
[progress.json](../../../docs/research/fundamental-plan/progress.json), newest
[RESULTS.md](../../../docs/research/fundamental-plan/RESULTS.md) entries and
[SEARCH-METHODS.md](../../../docs/research/fundamental-plan/SEARCH-METHODS.md).
Read `continuation_basis`, completed subtasks, `next_task`, the applicable task's
`next_action`, dependencies and active owner in the live progress record. Use
[SOURCE-DICTIONARY.md](../../../docs/research/fundamental-plan/SOURCE-DICTIONARY.md)
for the actual original operator/version route. The 30 October directive and
later user corrections supersede historical next actions; no task or result
number in this skill is a permanent priority. Read those files directly when
local access is available. With MCP-only access, locate their registered source
handles through the Atlas, read them, and check freshness. If current plan
access is missing, report that retrieval gap rather than treating an old ready
task as today's priority. Preserve existing owners and concurrent work.

Call `graph_status` and `theory_overview`, then
`research_state(query=the actual scientific question or exact topic/problem/goal ID,
lang="de" or "en", max_chars=24000)`. Reconcile its current-plan context and
original sources with the live plan before selecting a task. A historical
problem's unresolved task can be valid scoped evidence while being superseded
as an operative priority. Do not reopen a resolved theorem because a broad
bounded query did not show it. Check `source_state`, `coverage`,
`topic_coverage`, `problem_index` and `budget`; missing reviewed coverage is not
absence of previous research. Use `search_knowledge` and exact IDs to resolve
omissions before claiming a new gap.

For original-source, state/time, Higgs or 3+2 questions, use the compact
[research continuity router](references/source-bridge-continuity.md). Its dated
archives preserve earlier findings and exclusions, never a competing current
work queue. Load only the relevant historical section and its original source.
Do not read the full historical handoff into every new agent's context.

For historical dynamics, recovered correspondence or mixed-response work, use
[recovered historical sources](references/historical-dynamics-source-recovery.md).
It preserves version-specific curvature and kinetic conventions and separates
private provenance from anonymized scientific evidence.

## Retrieve complete evidence with bounded calls

For manuscript questions require `list_papers`, `search_papers`, `paper_outline`,
`paper_section`, `paper_evidence` and `paper_coverage`. Follow
[paper navigation](references/papers-navigation.md) for exact arguments and
scope. `search_papers` can return a `section:` ID: open it with `paper_section`,
then copy one of its exact `document_ids` to `paper_outline` when parent context
is needed. `list_papers` returns document IDs directly. Never pass a section ID
to `paper_outline` or invent an ID from a heading.

Use `paper_evidence` and bounded inclusive `read_source` ranges for the actual
statement, context and cited Python/Lean/Wolfram scope. For a load-bearing
passage request explicit `line_start` and `line_end` and continue until all
required lines are read, at most 200 per call. Inspect any preview, returned
extent, continuation and freshness metadata; an excerpt or bounded response is
not the full registered source. A source hash error calls for source refresh,
not a scientific reversal. A stored excerpt, including extracted PDF text,
retains its actual coverage limitation.

`summary_state: complete` describes current full-body chunk coverage with
validated source-grounded AI summaries/quotes, not scientific review. Also
inspect `context_state`: available source text can lack original prerequisites.
AI paper annotations remain `reviewed:false`; quotation matching does not prove
semantic entailment. Keep pending, error, stale and incomplete records visible.
List/search default to public scope. Use `scope: local_private` or `scope: all`
explicitly for authorized local research; retain exact-ID scope labels and
never copy private sources or local overlays automatically into public output.

For compact prior-attempt context use `research_context(id, max_chars=12000)`
and `research_brief`, after the durable `research_state` route. Inspect their
omission/coverage fields. Use `compare_aspects` on 2–6 exact IDs for shared
premises and divergent assumptions, and `solution_chain`, `topic_overview` or
`list_attempts` for their respective typed relations. Open only decisive
`get_dossier`, `audit_inputs` and source references. `get_dossier` is compact
by default; request `expanded=true` only for deliberately bounded detail. Do
not substitute a whole-graph dump for missing compact tools.

For `packet_state="partial"`, follow the complete `budget.continuation`,
including `detail_offset` and `expected_revision`. Keep whole findings,
conjunctive premise IDs and `active_task_handles`; a handle is navigation, not
the complete statement. Never mix pages across revisions. On old clients an
empty bounded `problems` array can be budget exhaustion. Reconnect to a current
catalog. If one editable record exceeds the tool limit, read only that complete
problem from `semantic_explorer/research/program.json` when local access is
available, retaining its latest MCP revision and every field. Without the
complete record preserve the existing state and report the retrieval limit.

## Preserve scope and advance the actual derivation

Before repeating a recorded strategy call `assess_research_plan` with its
exact `problem_id`, `strategy_id` and proposed `premise_ids`. Keep conjunctive
premises together and identify the changed premise that makes a retry useful.
Text similarity only suggests candidates. A ready task is an actionable
investigation, not proof of success or current global priority.

Compare actual maps, products, group actions, charges, grading, domains,
adjoints, state, locality, time and normalization before transferring a result.
Equal dimensions, spectra or group orders do not supply the physical map.
Distinguish internal labels, spacetime dimensions, field components and gauge
or global symmetries. Name the first newly chosen geometry, coupling, state or
readout and work on its origin when that is the user's missing derivation.
Preserve the actual source's full charge/mode sectors, exterior paths and
memory; an observer/input slice need not be closed after every operation.
Separate a supplied operator or response from an executable physical
intervention and its selection rule.

`required_together` identifies one derivation's joint premises. `references`,
`mentions`, `implemented_check`, `formalization_candidate`, `related_structure`,
`shares_seed` and `consistency_check` do not prove the target. Shared seeds give
correlated readouts. Audit script inputs and later linked proofs before calling
an old local gap globally unresolved. Lean inventory is static: read the exact
statement and axiom dependencies and require a matching kernel/build receipt
before claiming formal validation. Script presence implies no execution.
Preserve mathematical consequences, conditional physical identifications,
prospective independent tests and full physical completion as distinct scopes.
A working auxiliary model is not the requested origin derivation.

## Record a source-supported update and verify the interface

After a meaningful source-supported result, failure or correction, load
`research_state(query=exact problem ID, editable=true, max_chars=60000)`. Update
the complete bilingual record, related findings, resolved gaps/tasks,
prerequisites and source references together. Use
`record_research_update(expected_revision=returned revision, problem=complete
updated record, premises=new definitions if needed, reason=source-backed
explanation)`. Reload and reconcile on conflict; never overwrite a partial
packet or hand-edit generated graph data. Retain established subresults and
scoped negatives without promoting a parent physical goal. Reconcile the live
plan's results, progress and next action when that research program is active.
Scientific status belongs in the Atlas and original records. Personal memory
changes require explicit user authorization. Do not promise background updates
without an active session or requested automation.

Check the actual connected catalog for `research_state`, `assess_research_plan`,
`record_research_update`, the compact context tools and paper tools. An old MCP
process can expose obsolete capabilities. Reconnect the client; a fresh
`semantic_explorer/run.sh mcp` stdio connection can be used without stopping
other clients. Package versions and a current web page do not prove a current
client catalog. The local web explorer runs at http://127.0.0.1:8731 through
`semantic_explorer/run.sh serve`.

Use `python3 semantic_explorer/refresh.py --check` to check source/ledger
freshness without AI calls, and the documented refresh route for maintenance.
`semantic_explorer/run.sh export` emits the same composed graph as the service.
The retired graph/primer/query CLI is no second current research interface.
Structural checks, successful tool calls and traversal completeness do not
certify every mathematical or physical claim.
