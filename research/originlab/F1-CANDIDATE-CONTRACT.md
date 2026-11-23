# Candidate F1 contract for a new origin model

Added 23 November 2026. `f1_native_contract.json` registers no native contract
(`contract_status: NOT_CONFIGURED`). That is a missing acceptance
configuration, not a research result. A new origin model does not need it,
and it does not need the original R113 endpoint values: it is checked
through its own candidate contract.

## What a candidate contract is

Two JSON files and one or more Lean modules in one folder, for example
`experiments/exploration/<date>-<slug>/contract/`:

`lean-contract.json`

```json
{"modules": [{"file": "Candidate.lean", "sha256": "<sha256 of the file>"}],
 "goal": "Candidate.Goal",
 "producers": {"Candidate.model": 0},
 "proofs": {"Candidate.complete": 0}}
```

`contract.json`

```json
{"lean_contract": "lean-contract.json",
 "lean_contract_sha256": "<sha256 of lean-contract.json>",
 "gates": {"G0": "Candidate.G0", "G1": "Candidate.G1", "G2": "Candidate.G2",
           "G3": "Candidate.G3", "G4": "Candidate.G4", "G5": "Candidate.G5",
           "G6": "Candidate.G6", "G7": "Candidate.G7"}}
```

The Lean module defines `Goal` and one predicate per gate on the model type,
and a theorem whose statement is the conjunction of `Goal` and all eight
gate predicates on the same producer. `gates.json` lists what each gate
asks; the predicates must state that content for the candidate's own source,
state and time. G3 asks for the candidate's own source identity and
responses, not for the R113 marks. The fixture in `test_f1_solution.py`
shows the minimal syntax.

## Running it

```sh
.venv/bin/python -B research/originlab/main.py --out /tmp/candidate \
  --campaign f1_solution --campaign f1 \
  --f1-contract path/to/contract.json --f1-contract-sha256 <sha256 of contract.json>
```

A checked proof returns `F1_CONTRACT_PROVED`. `evidence_gates.assess` then
reports `physical_status: CANDIDATE_PROVED` with `CANDIDATE_PASS` on all
eight gates. `physical_solution` stays false: a contract written by the same
session that wrote the model can encode the wrong physics.

## From candidate to native

The user or a reviewing session reads the Lean predicates against
`gates.json` and the cited sources. If they state the physical content of
each gate, register the contract in `f1_native_contract.json` as
`"contract": {"file": ..., "sha256": ...}` and set `contract_status` to
`REGISTERED`. The next run returns `F1_PROVED` for that contract.

## Gate statuses

`PASS`, `CANDIDATE_PASS`, `PARTIAL`, `BLOCKED` (a named prerequisite within
the same route is missing), `NOT_TESTED`, `NOT_APPLICABLE`, `REFUTED`.
The missing R113 packet sets G3 of the `original_source` route to `BLOCKED`
and records that route's dependency under `G4.route_dependencies`. The G4
status itself stays `NOT_TESTED` for every other model.
