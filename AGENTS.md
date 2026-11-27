# Working in TFPT

This repository owns physical TFPT, its compiler mathematics and its physical
interfaces. Prime-number, RH, Weil-window and factorisation research belongs to
[Prime Front](https://github.com/sthamann/prime-front); do not reintroduce it here.
Current explicit user instructions always take precedence over this file.

## Pick the mode first — user decision of 2026-11-20

Every research session runs in exactly one of two modes. Decide from the user's
words; when unclear, ask once.

**Exploration** — the user asks for new ideas, a fresh approach, creativity,
"völlig neu", "anders", "Neuansatz", a broad search, a complete solution, or
asks whether the research has run aground ("verrannt", "festgefahren",
"Grundannahmen prüfen").
Goal: a complete TFPT solution (one principle that yields the inputs P1 and P2,
reproduces the existing results and supplies state, time, 3+1 and gravity).
Read [`docs/research/TFPT-STAND.md`](docs/research/TFPT-STAND.md) (what works,
patterns, what is missing) and follow
[`docs/research/EXPLORATION.md`](docs/research/EXPLORATION.md). Atlas and the
catalogue `docs/catalog/TFPT_Pruefkatalog_2026-11-20.xlsx` are available for
details at any time. TFPT-STAND is a selection: before calling anything new or
open, look it up in [`docs/research/LEDGER-KARTE.md`](docs/research/LEDGER-KARTE.md)
(all active ledger claims by area, generated) or with
`research/originlab/ledger_map.py --find 'pattern'`; the Atlas gives the same
view in `theory_overview.ledger_areas` and `search_knowledge` answers German
and English questions with ledger claims first. Do not read PLAN.md, RESULTS.md, progress.json or contract
chains beforehand, and do not start from an existing TFPT object: those objects
are what a new idea must produce.

**Continuation** — the user asks to continue the plan, a named task, contract or
result. Read `docs/research/fundamental-plan/PLAN.md` (top section), the compact
`progress.json` and the newest `RESULTS.md` entries, then the live Atlas. Continue
the named task. Historic per-session records are in `progress-history.json`; read
them only for a specific ID.

Neither mode turns an automatic "Nächster Schritt" from the last contract into an
obligation. Before a new contract in either mode, check in one line whether it is
just another variant of the last three contracts; if so, pick something else.

**Review of foundations** (part of Exploration): when the user asks whether the
work has run aground, list which recorded identifications, gates or search
defaults a result depends on, and check each one against its source instead of
assuming it. Recorded identifications can be wrong; a new principle may
reclassify them, provided it says which one and why.

## Original source routes are not gates — user decision of 2026-11-23

R113 and the other original source contracts of November 2026 test one specific
claim: that the original source of those contracts yields its marked endpoint
values. They are scoped tests of that route, not prerequisites of a fundamental
solution, of F1, or of any other origin model.

- Never cite missing R113 values (Z0, ZU, ZY, ZUY, the blocked fibre packet, the
  absolute empty reference) as the reason a different candidate stays open, and
  never propose "unblock R113" as the default next step. Name the concrete map the
  candidate lacks instead.
- R113 is required only when a result claims to be that original source or claims
  its specific physical identification.
- A new principle must reproduce the established results (`TFPT-STAND.md`) or
  reclassify their scope explicitly; it need not reuse the R113 computation.
- OriginLab reports the R113 call as the `original_source` route; its block does
  not block other routes. New origin models are checked through their own
  candidate contract (`research/originlab/F1-CANDIDATE-CONTRACT.md`). The empty
  native contract in `f1_native_contract.json` is a missing acceptance
  configuration, not a research result.

## Scientific discipline (both modes)

- Distinguish declared input, exact consequence, numerical observation,
  conditional physical identification and open derivation. Ledger
  (`verification/status_ledger.csv`) and original sources own claim status.
- Matching numbers or dimensions do not establish an operator dictionary.
  Name the first newly chosen geometry, state, coupling or selection rule.
- Smallest decisive test first; new models need a question, success criterion
  and stopping rule. Preserve scoped negative results.
- A result may be marked **Kandidat**: it reproduces a target from the list
  without TFPT inputs or fitted constants. That is a positive outcome worth
  reporting, not a proof and not a ledger promotion.
- AI annotations, graph edges and Atlas links are navigation, not proof.

## Atlas and tools

After a ledger change run `python3 semantic_explorer/refresh.py` (offline, no
model calls): it regenerates the ledger map, the Atlas ledger summary and the
source hashes of the theory map; `refresh.py --check` and
`research/originlab/test_ledger_map.py` fail while any of them is stale.

`tfpt-knowledge` (`.mcp.json`) is the sole theory-atlas interface:
`graph_status`, `theory_overview`, `research_state`, `assess_research_plan`,
`record_research_update` (with `expected_revision`). For an existing entry,
use `research_state(editable=true, record_ids=[])` on its exact problem ID for
the revision, or request the exact records to edit. Send only changed records
in `record_research_update(patch=…)`; reserve the full `problem` payload for
creation or explicit full replacement. Code discovery:
`codebase-memory-mcp`, then targeted text search. Content registry:
`python3 verification/make_content_registry.py --query "topic"`, details in
`docs/registry/INDEX.md`. Search programs live in `research/originlab`
(`.venv/bin/python`); extend them instead of writing parallel tools.
Catalog of all checks and experiments: `docs/catalog/TFPT_Pruefkatalog_2026-11-20.xlsx`.

## Recording results

Exploration: one folder under `experiments/exploration/<date>-<slug>/` with
`IDEAS.md`, scripts and `VERDICT.md` (Kandidat / verworfen / offen + why).
Continuation: contract under `experiments/theory-contracts/`, entry in
`RESULTS.md`, update the compact `progress.json` (`next_task`, `continuation`,
task statuses). Persist source-reviewed changes through the Atlas update.
Never store per-session `latest_*` records in `progress.json`; append them to
`progress-history.json`.

## Source ownership and editing

- `documents/` holds editable paper modules; `tex-artefacts/` shared style and
  fragments. `docs/registry/`, `verification.tex`, document maps, script index,
  website metadata and Atlas data are generated — regenerate, never hand-edit.
- Exploratory scripts stay under `experiments/` until an explicit reviewed promotion.
- `_archive/` is historical; do not rewrite old verdicts.
- When moving TeX, preserve order, macros, labels and references; use the shared
  include reader. See `docs/LATEX_MAINTENANCE.md` for build, registry, website and
  manifest regeneration (`bash build.sh audit`, manifests `--check` last).
- For an executable or scientific change run the affected verification and its
  negative control before editing paper, ledger or status. Report exactly what ran.

## Working tree and delivery

Preserve unrelated modified and untracked research; never reset or broadly add
another session's files. Report what changed, what passed and what is unverified.
Do not ask again for an already authorized reversible step.
