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

Physical acceptance requires a checked native proof covering the eight
gates on the same producer. Source replay and comparison candidates retain
their scope in the report. The current source inventory and the available
formal predicates limit what the program can prove; unbounded stage
enumeration covers the implemented grammars, whose completeness for the
physical problem is open.
