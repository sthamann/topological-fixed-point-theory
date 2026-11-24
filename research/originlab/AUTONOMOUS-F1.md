# Autonomous F1 search

From a clean checkout, install `research/originlab/requirements.txt`. The
optional Lean proof lanes use the pinned toolchain in `proof_search.py`.

```sh
python3 -B research/originlab/autonomy.py --out /tmp/originlab-f1 --until-proved
```

The same command and output path resume the search after interruption. Use
`--start-resolution 3 --max-resolution 4` without `--until-proved` for
a finite run. Run again with the same path and a higher bound to extend it.

The supervisor discovers standalone and inline R113 packets under the
research source directories. An exact domain identity and digest, the
source-reviewed input designation, and the existing endpoint checks are
required for replay. The search records rejected inputs separately. Replay
uses a temporary domain binding and leaves the original source untouched.
Successful replays automatically enter the affine Q0 comparison. Their
source and dependency digests participate in checkpoint invalidation.

Every resolution schedules F1 synthesis, even with an empty native proof
registry. An uninterrupted `--until-proved` run increases the proof and
coherent model stage, refreshes the source inventory after each stage, and
allocates more worker time at higher stages. The finite topology searches
retain their own supported resolution limits. SQLite stores jobs and
results. The supervisor recovers interrupted jobs and keeps previous
evidence by epoch. `latest.json`, the stage's `frontier.json` and its review
outbox record the outcome and identify the evidence used.

A job with unchanged engine, sources and specification reuses its stored
result in later stages and epochs. Proof jobs and model proposals always
run. `latest.json` reports the number of reused results.

Each stage records whether a gate, the source state or the acceptance
state moved. `work-order.json` lists the inputs that no higher stage can
produce, what each one blocks, and how the supervisor picks it up. After
`--stall-rounds` unchanged stages (default 3) an `--until-proved` run
exits with status 3. With `--watch-seconds N` it waits instead and opens a
new epoch as soon as a source, packet, contract or engine file changes:

```sh
python3 -B research/originlab/autonomy.py --out /tmp/originlab-f1 --until-proved --watch-seconds 300
```

`f1_proof_reachable` is false while `f1_native_contract.json` registers no
contract or the Lean backend is absent. In that configuration no stage
returns `F1_PROVED`. The work order marks the missing contract as
`missing_acceptance_configuration`: a configuration gap, not a research
result. A new origin model is checked through its own candidate contract
and reaches `CANDIDATE_PROVED`; see `F1-CANDIDATE-CONTRACT.md`.

## Routes and the original R113 source

The original R113 endpoint call is one route (`original_source`). Since
23 November 2026 it runs after the source-free generators
(`genesis_foreign`, `genesis_rules`, `process_rules`, `structure_scenarios`,
`atlas_scenarios`) in both `main.py` and the supervisor. The frontier
`state` follows the acceptance evidence (`OPEN`, `CANDIDATE_PROVED`,
`F1_PROVED`); the R113 status is reported separately as
`original_source_route_state`. A missing R113 packet blocks G3 of that
route only; G4 records the dependency under `route_dependencies` and keeps
its own status. Its open inputs are listed as `original_route_open_inputs`
and are not the next step of the search. `main.py` fingerprints every
Atlas problem in its run context, not only the original-source problems.

`work-order.json` also carries a computed completeness audit: for each
gate the registered acceptance predicate and the presence of its scoped
evidence files, the number of source-reviewed packets, and the lanes
whose jobs still change at far stages. `search_complete` in `latest.json`
is true only when every gate has a predicate, a source or contract is
bound, and at least one lane still grows.

## Scenarios from the whole Atlas

Every stage runs `atlas_scenarios.py` on the reviewed Atlas program. It
groups premises that findings use together into contexts, lists each open
task with its required premises and the routes already closed or reusable,
and enumerates the ways to join all contexts the open tasks need. Each
scenario names its identifications, the recorded evidence for each, the
contexts without any recorded link, and the findings available once the
contexts are read as one source. It also compares all recorded source
digests with the working tree.

```sh
python3 -B research/originlab/atlas_scenarios.py > /tmp/atlas-scenarios.json
```

A supplied scenario enters the queue through `--cases` as an
`atlas_scenario_case` job. It names the premises it adopts, the tasks it
resolves and the Atlas strategies it uses. A closed strategy is accepted
only with a changed premise that its exclusion did not have. Contexts and
shared sources are recorded evidence for a link; they prove no
identification.

## Generated structure scenarios

Every stage also runs `structure_scenarios.py`. It needs no model: it
generates objects with computed features, derives the fingerprint from the
established Atlas findings, ranks objects, pairs and triples against random
fingerprints, and queues its best scenarios as `atlas_scenario_case` jobs
with a `hypothesis:structure-…` premise. `latest.json` reports the best
object, the number of significant layers, the unexplained numbers and the
withheld p-values under `structure_scenarios`.

## Model proposals and computed tests

```sh
python3 -B research/originlab/autonomy.py --out /tmp/originlab-f1 --until-proved \
  --scenario-model opus --scenario-agent-max 3
```

With `--scenario-model` each stage runs three model steps. All prompts
carry `THEORY-BRIEF.md`. The brief opens with the task in plain physical
terms and treats the R113 route as one route among others; a scenario
needs the R113 marks only if it claims to be that original source.

1. Proposal. A tool-free call receives the compiled record and earlier
   verdicts and returns up to six scenarios. A scenario may add a rule
   absent from the record as `new_premises` with an id `hypothesis:…`.
2. Test. For up to `--scenario-agent-max` admissible scenarios, those with
   a new rule first, a call with read-only file tools derives one
   consequence and writes a script that compares it with a recorded exact
   value. The script may read the repository and run its checks. It runs
   under `sandbox-exec` without network, with writes confined to a scratch
   directory and without access to other home files.
3. Referee. A tool-free call grades the test. Only the grade
   `INDEPENDENT_AND_SPECIFIC` lets the outcome count as
   `SURVIVES_COMPUTED_PREDICTION` or `REFUTED_BY_COMPUTATION`. Otherwise
   the scenario is reported as `TEST_NOT_DISCRIMINATING`.

4. Research threads. With `--scenario-research-max N`, up to N calls per
   stage take two anchor groups of the brief with their source files and
   look for one operator that yields both with the recorded factors on the
   recorded spaces. A thread revises after each computed outcome, up to
   `--scenario-research-rounds` times, and may conclude that no common
   mechanism exists. `--no-scenario-proposals` runs these threads alone.

A script that stops before its verdict is repaired by a tool-free call,
up to three times. A failed comparison of a computed with an expected
value must be reported as a failed prediction. With a model the campaign
ends only on success or when `--scenario-total-budget-usd` is spent,
checked between stages; an interrupted stage is not booked.

`scenarios.json` keeps every tested scenario with its verdict, and
`scenario-checks/` keeps the scripts. Budgets per call are set with
`--scenario-budget-usd` and `--scenario-agent-budget-usd`. Proposals,
scripts, verdicts and grades are unreviewed model output. They change no
gate and no Atlas record.

Physical acceptance requires a checked native proof covering the eight
gates on the same producer. Source replay and comparison candidates retain
their scope in the report. The current source inventory and the available
formal predicates limit what the program can prove; unbounded stage
enumeration covers the implemented grammars, whose completeness for the
physical problem is open.
