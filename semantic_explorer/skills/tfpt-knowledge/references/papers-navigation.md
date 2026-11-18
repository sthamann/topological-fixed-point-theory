# Paper navigation in the shared TFPT Atlas

Use the single `tfpt-knowledge` service. The source registry owns exact
`document:` and `section:` identities; the local paper overlay adds searchable,
source-grounded AI descriptions. It neither replaces original sources nor
updates the claim ledger or curated research decisions.

## Tools and bounded arguments

| Tool | Arguments | Purpose |
|---|---|---|
| `paper_coverage` | none | Inspect actual imported/expected sections, summary states, source freshness and public/private counts. |
| `list_papers` | `query?`, `scope?`, `offset=0`, `limit=20`, `max_chars=16000` | Find document views and their exact IDs. |
| `search_papers` | `query`, `document_id?`, `evidence_kind?`, `scope?`, `offset=0`, `limit=20`, `max_chars=16000` | Search source bodies and available semantic descriptions; constrain a known document or evidence kind. |
| `paper_outline` | `document_id`, `offset=0`, `limit=30`, `max_chars=16000` | Read ordered sections and copy the precise target section ID. |
| `paper_section` | `section_id`, `max_chars=16000` | Read one section's source locations, semantic summary state, assumptions/results/limits and source references. |
| `paper_evidence` | `id`, `evidence_type?`, `scope?`, `offset=0`, `limit=20`, `max_chars=16000` | Follow typed relations from a document, section, claim, Python script or Lean declaration. Returns source handles for `read_source`. |
| `read_source` | `source_id`, `line_start?`, `line_end?` | Read the registered original source with hash validation; inclusive line ranges, at most 200 lines per call. |

Evidence filters match stored relation/evidence type strings; use a value
returned by the service rather than guessing an enum. Paged tools expose
`total`, `offset`, `returned`, `omitted` and `next_offset`, plus a character
budget with `actual_chars`, `max_chars` and `truncated`. Follow that continuation;
an empty bounded response does not prove absence. `paper_section` has no offset:
use its source references with bounded `read_source` ranges when context is cut.

## Practical next-session route

For a fundamental research continuation, first read the live
`docs/research/fundamental-plan/PLAN.md`, `progress.json`, newest `RESULTS.md`
entries and `SEARCH-METHODS.md`, including established `continuation_basis`
results. Reconcile the actual question and current owner with `research_state`.
The paper route below retrieves evidence; it does not select a new research
priority or supersede the live plan.

1. Call `graph_status` and `theory_overview`, then `paper_coverage`. Distinguish
   inventory completeness, completed source-grounded AI coverage, current source
   hashes and curated review. Vector availability is another separate measure.
2. For example, search `carrier_rank_pascal_unique` with `search_papers`, a
   small limit and `scope: public`. Search can return exact `section:` IDs;
   open such a hit with `paper_section(section_id=the returned ID)`. Copy one
   of that section's `document_ids` into `paper_outline(document_id=...)` for
   parent context. `list_papers` supplies document IDs directly. Never pass a
   section ID to `paper_outline` or invent an ID from a filename or heading.
3. Open the exact section with `paper_section`. Read the full original section,
   its parent context and referenced prerequisites. A summary or isolated
   theorem name can omit a hypothesis or a prior definition.
4. Use `paper_evidence` and `read_source` to inspect actual source wording and
   evidence. For a load-bearing passage request explicit `line_start` and
   `line_end`, at most 200 inclusive lines per call, and continue until the
   required source range and prerequisites are read. Inspect any preview,
   returned extent, continuation and freshness metadata. A stored excerpt or
   bounded response is not the complete original; extracted PDF text retains
   its coverage limits. Preserve the returned relation type and match basis: explicit
   citation, lexical topic suggestion and same-identifier join are different.
   A graph edge is not a derivation. Read Lean statements and axiom dependencies;
   do not imply a fresh kernel run or script execution.
5. Call `research_state` for the exact topic/problem or actual scientific
   question. Use the normal `assess_research_plan`/versioned update workflow
   before repeating a strategy or recording a source-reviewed result. Do not
   manufacture a new open problem merely because a reviewed match is missing.

## Meaning of the paper records

For sections, `summary_state` is `pending`, `complete`, `error`, `stale` or
`source_error`; an aggregate document can also be `partial`.
`summary_state_as_recorded` retains the imported state.
Only source-current, validated summaries/quotes covering every bounded chunk
of the complete source body count as `complete`. Read `source_current` and
`source_error`; source changes require re-extraction, not a scientific reversal.
Also read `context_state` (`current`, `stale` or `incomplete`),
`context_complete`, `context_current` and `context_error`. `source_current`
checks the physical section; document prerequisites can still be missing or
stale. An available-body summary may be complete while original document
context is incomplete; `semantic_complete` then remains false. Literal quote
span validation does not establish that the generated statement follows from
that quotation. Read the original premises and qualification yourself.
AI summaries, concepts, assumptions, results, limits and quotations remain
`reviewed:false`. They help retrieve and understand a source; they are not
curated research review, proof certification or an upgrade of the ledger.

`paper_coverage` reports `expected_sections`, `imported_sections`, `complete`,
`error`, `stale`, `pending`, `source_error`, `public_sections`, `private_sections`,
`documents`, `physical_sources`, `current_chunks`, `total_chunks` and
`semantic_complete`. Context counts distinguish `context_stale_sections` and
`context_incomplete_sections`; paper embedding coverage is separate.
Report actual counts and remaining states. Do not infer
complete AI coverage from a complete inventory or from several good examples.
Even `semantic_complete: true` does not mean every statement has been reviewed
or proved. Do not use whole-graph vector coverage as a scientific completeness score.

Publication scope is `public` or `local_private`; content scope distinguishes
active, archive, exploration, historical, support and private material. List and
search default to `public`; `local_private` (alias `private`) explicitly selects
private records, `all` includes both. Content selectors `active`, `archive`,
`exploration`, `historical` and `support` select public material. Exact private
IDs can be read directly and remain visibly labelled. Historical results keep
their dates and status; private results must not be published automatically.

The public importer does not consult private records/cache. Explicit local
import may include local/private inputs in ignored `data/papers.json` and its
local cache; this is not a new public corpus. No private paper overlay belongs
in website downloads or a public content-registry regeneration.

## Maintain source and semantic coverage

From the repository root, `python3 semantic_explorer/papers.py --check` is a
read-only coverage/source/cache check. Use `--offline` to regenerate navigation
from exact current sources and valid cached descriptions/vectors without model
calls. `python3 semantic_explorer/papers.py --include-local --embeddings`
resumes missing or invalid source-grounded descriptions and then their vectors.
`--section section:ID` limits a deliberate extraction; complete global coverage
still requires the remaining sections. `--public-only` explicitly excludes
private inputs. An ordinary refresh preserves the previously selected private
scope and valid vectors and does not race an active writer.

Inspect `paper_coverage` after maintenance. Source changes, incomplete includes,
invalid generated prose and missing vectors are different states. Check exact
source references and the original ledger rather than repairing a scientific
status to satisfy an inventory. Generated descriptions can use plain-text math
names and ASCII transliteration; quoted original source remains exact UTF-8.

## Missing tools or stale clients

Check the actual connected MCP catalog. A fresh `semantic_explorer/run.sh mcp`
session rebuilds a stale service binary without terminating existing clients;
reconnect old clients to expose new tools. The HTTP counterparts are
`/api/papers`, `/api/papers/outline`, `/api/papers/section`, `/api/papers/search`,
`/api/papers/evidence` and `/api/papers/coverage` on the local-only service.
If the paper layer is unavailable, use bounded public registry queries and
original sources and state the retrieval limitation. Never substitute a full
private graph dump or claim an unavailable semantic review.

## Paper versus Atlas status review (2026-10-18)

Use `research_state` for source-reviewed scope; a claim annotation's statement
that a source was not supplied refers to its original extraction packet, not
to absence from the entire repository or to an unresolved mathematical result.
`get_dossier` retains `cache_freshness`. `research_brief` carries annotation
freshness on result/gap samples and does not report complete gap coverage when
cited annotation sources are stale or unverified. Even current annotation
sources do not establish complete scientific coverage.

The targeted paper/source comparison corrected `FORM.GLUE.01`,
`FORM.LADDER.01` and `AX.P2.01`. General AI extraction provenance remains
`reviewed:false`; named source corrections do not imply a fresh Lean run.
Find the durable reviewed results through these exact records:

- `problem:carrier-pascal-glue-formal-scope`: unbounded natural-number Pascal
  uniqueness under the supplied exhaustion equation; precise finite glue
  element/orbit declarations. Do not merge the separate bounded PascalLadder
  result with the universal CarrierRankUniqueness theorem.
- `problem:physical-time-projector-transfer`: additionally retains
  `finding:shared-register-gkls-embedding-exact-q-interval` and
  `finding:marked-m4-channel-all-powers-homogeneous-gksl-excluded`. These refer
  to two different fixed channel families. Same finite system and
  time-homogeneous GKLS/GKSL are essential premises; larger environments,
  time-dependent execution, GNS Hilbert transfer and physical source time
  are separate questions.

Alpha global-root uniqueness, scoped flavor ratios, current-response
construction, positive-word-kernel GNS reconstruction and native carry/memory
were already represented. The comparison did not establish physical source
selection, a chiral interacting 3+1D theory, or complete gravitational origin.
It was a targeted comparison with counterchecks, not an exhaustive independent
proof review of every paper or every legacy annotation.
